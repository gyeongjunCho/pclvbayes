args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
repo <- normalizePath(file.path(script_dir, "../.."))
source(file.path(script_dir, "profile_helpers.R"))
source(file.path(repo, "benchmarks/mtist/mtist_adapter.R"))
config_path <- Sys.getenv("PCLV_PROFILE_CONFIG", file.path(script_dir, "configs/profile.R"))
config <- source(config_path)$value
result_dir <- file.path(script_dir, "results", format(Sys.time(), "%Y%m%d-%H%M%S"))
dir.create(result_dir, recursive = TRUE)

profile_library <- tempfile("pclv-profile-library-")
dir.create(profile_library)
on.exit(unlink(profile_library, recursive = TRUE), add = TRUE)
install_timing <- time_expression("setup", "isolated_package_install", {
  withr::with_libpaths(profile_library, action = "prefix",
    devtools::install(repo, quiet = TRUE, upgrade = FALSE,
                      dependencies = FALSE, keep_source = TRUE))
})
.libPaths(c(profile_library, .libPaths()))
library(pclvbayes)
study <- load_mtist_study(config$dataset_id, Sys.getenv("MTIST_ROOT", "~/mtist"))
model_compile <- time_expression("setup", "canonical_model_load_or_compile", pclvbayes:::get_pclv_model(quiet = TRUE))

meta <- study$sample_metadata
meta$Sample <- rownames(meta)
mat <- t(study$relative)
component <- list()
component[[1L]] <- time_expression("component", "global_cv_spline_smoothing", {
  for (z in seq_len(config$component_repetitions))
    sm <- pclvbayes:::.precompute_spline_smoothed(mat, meta, study$taxa, eps = 1e-6, min_unique_times = 3L)
  sm
})
sm <- component[[1L]]$value
component[[2L]] <- time_expression("component", "canonical_pair_preprocessing", {
  for (z in seq_len(config$component_repetitions))
    pair <- pclvbayes:::.make_pair_inputs_glv(
      sm_mat = sm, meta_df = meta, j = study$taxa[[2L]], i = study$taxa[[1L]],
      min_pairs = 4L, zero_mode_alr = "minpos_time", minpos_alpha = .5,
      minpos_base = "ij", eps_fixed = 1e-6, lib_eps_c = .65,
      rest_floor_frac = 1, alr_cap = 12, smooth_scale = "logra",
      alr_spline_cv = TRUE, nz_partner_min_frac = .15)
  pair
})
pair <- component[[2L]]$value
component[[3L]] <- time_expression("component", "kfold_split_index_construction", {
  for (z in seq_len(config$component_repetitions * 20L))
    splits <- pclvbayes:::.make_repkfold_splits(pair$subject, config$kfold_K, config$kfold_R, config$seed)
  splits
})
stan_probe <- list(N = nrow(pair), y = pair$y, xi = pair$xi, xj = pair$xj,
                   S = length(unique(pair$subject)), sid = as.integer(factor(pair$subject)),
                   prev = integer(nrow(pair)), dt = numeric(nrow(pair)))
json_path <- tempfile("pclv-profile-", fileext = ".json")
on.exit(unlink(json_path), add = TRUE)
component[[4L]] <- time_expression("component", "cmdstan_json_serialization", {
  for (z in seq_len(config$component_repetitions * 20L)) cmdstanr::write_stan_json(stan_probe, json_path)
  invisible(json_path)
})
component[[1L]]$timing$calls <- config$component_repetitions
component[[2L]]$timing$calls <- config$component_repetitions
component[[3L]]$timing$calls <- config$component_repetitions * 20L
component[[4L]]$timing$calls <- config$component_repetitions * 20L

# A direct one-direction fold probe measures real fold preparation, sampling,
# Kalman scoring, and aggregation even when a deliberately short main fit is
# not diagnostically eligible to enter its own K-fold stage.
ordered <- pair[order(pair$subject, pair$time), , drop = FALSE]
prev <- integer(nrow(ordered)); dt <- numeric(nrow(ordered))
for (ix in split(seq_len(nrow(ordered)), ordered$subject)) {
  if (length(ix) > 1L) {
    prev[ix[-1L]] <- ix[-length(ix)]
    dt[ix[-1L]] <- diff(ordered$time[ix])
  }
}
stan_base <- list(N = nrow(ordered), y = ordered$y, xi = ordered$xi, xj = ordered$xj,
                  S = length(unique(ordered$subject)), sid = as.integer(factor(ordered$subject)),
                  prev = prev, dt = dt)
sample_args <- list(chains = config$chains, parallel_chains = config$chains,
                    iter_warmup = config$iter_warmup, iter_sampling = config$iter_sampling,
                    seed = config$seed, init = 0.2, adapt_delta = 0.98,
                    max_treedepth = 14L, metric = "diag_e", refresh = 0L)
kfold_probe <- time_expression("kfold_probe", "fold_preparation_sampling_scoring_aggregation", {
  pclvbayes:::.repkfold_eval(model_compile$value, stan_base, sample_args, ordered,
                             K = config$kfold_K, R = config$kfold_R, seed = config$seed,
                             silent_sampler = TRUE, max_retries = config$max_retries,
                             n_workers_kfold = config$kfold_workers, min_pairs = 4L,
                             freeze_retry_hypers = TRUE)
}, scope = "combined fold preparation, CmdStan execution, Kalman scoring, and aggregation")

run_fit <- function(label, outer_workers) {
  rprof <- file.path(result_dir, paste0(label, ".Rprof"))
  gc(); before_mem <- platform_memory(); before_tmp <- directory_snapshot(tempdir())
  Rprof(rprof, interval = 0.01, memory.profiling = TRUE)
  on.exit(Rprof(NULL), add = TRUE)
  captured_warnings <- character()
  timed <- withCallingHandlers(time_expression(label, "end_to_end_fit", fit_pclv_bayes(
    study$physeq, "subject", "time", taxa_vec = study$taxa,
    chains = config$chains, iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling, seed = config$seed,
      n_workers_outer = outer_workers, n_workers_kfold = config$kfold_workers,
    kfold_K = config$kfold_K, kfold_R = config$kfold_R
  ), scope = "includes preprocessing, CmdStan startup/sampling, K-fold, and assembly"),
    warning = function(w) captured_warnings <<- c(captured_warnings, conditionMessage(w)))
  Rprof(NULL)
  after_mem <- platform_memory(); after_tmp <- directory_snapshot(tempdir())
  list(fit = timed$value, timing = timed$timing, before_mem = before_mem,
       after_mem = after_mem, before_tmp = before_tmp, after_tmp = after_tmp,
       hotspots = summarise_rprof(rprof), rprof = rprof,
       warnings = unique(captured_warnings))
}

sequential <- run_fit("sequential", config$sequential_outer_workers)
parallel <- run_fit("parallel", config$parallel_outer_workers)
equivalent <- isTRUE(all.equal(scientific_result_signature(sequential$fit),
                               scientific_result_signature(parallel$fit), tolerance = 0))

timings <- do.call(rbind, c(lapply(component, `[[`, "timing"),
                           list(install_timing$timing, model_compile$timing, kfold_probe$timing,
                                sequential$timing, parallel$timing)))
utils::write.csv(timings, file.path(result_dir, "phase_timings.csv"), row.names = FALSE)
utils::write.csv(rbind(transform(sequential$hotspots, run = "sequential"),
                       transform(parallel$hotspots, run = "parallel")),
                 file.path(result_dir, "rprof_hotspots.csv"), row.names = FALSE)
utils::write.csv(object_sizes(list(study = study, smoothed_matrix = sm, pair_input = pair,
                                   stan_data = stan_probe, sequential_fit = sequential$fit,
                                   parallel_fit = parallel$fit)),
                 file.path(result_dir, "object_sizes.csv"), row.names = FALSE)
resource <- rbind(cbind(run = "sequential_before", sequential$before_mem, sequential$before_tmp[-1]),
                  cbind(run = "sequential_after", sequential$after_mem, sequential$after_tmp[-1]),
                  cbind(run = "parallel_before", parallel$before_mem, parallel$before_tmp[-1]),
                  cbind(run = "parallel_after", parallel$after_mem, parallel$after_tmp[-1]))
utils::write.csv(resource, file.path(result_dir, "resources.csv"), row.names = FALSE)
failures <- function(x, run) {
  raw <- x$raw
  data.frame(run = run,
             direction = c(paste0(raw$j, "->", raw$i), paste0(raw$i, "->", raw$j)),
             ok = c(raw$direction_ok_ij, raw$direction_ok_ji),
             diagnostic_class = c(raw$diagnostic_class_ij, raw$diagnostic_class_ji))
}
utils::write.csv(rbind(failures(sequential$fit, "sequential"), failures(parallel$fit, "parallel")),
                 file.path(result_dir, "run_status.csv"), row.names = FALSE)
writeLines(c("[sequential]", sequential$warnings, "[parallel]", parallel$warnings),
           file.path(result_dir, "warnings.txt"))
cpu_lines <- tryCatch(readLines("/proc/cpuinfo", warn = FALSE), error = function(e) character())
cpu_model <- sub("^[^:]+: *", "", grep("^model name", cpu_lines, value = TRUE)[1L])
mem_lines <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character())
total_ram_kb <- suppressWarnings(as.numeric(sub("^MemTotal:\\s*([0-9]+).*", "\\1",
                                                grep("^MemTotal:", mem_lines, value = TRUE)[1L])))
metadata <- list(repository_sha = system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
                 timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), config = config,
                 dataset_id = config$dataset_id, result_equivalent = equivalent,
                 system_time_utility = system_time_available(), R = R.version.string,
                 cmdstan = as.character(cmdstanr::cmdstan_version()), session = capture.output(sessionInfo()))
metadata$operating_system <- paste(Sys.info()[c("sysname", "release")], collapse = " ")
metadata$cpu_model <- cpu_model
metadata$logical_cores <- parallel::detectCores(logical = TRUE)
metadata$total_ram_kb <- total_ram_kb
metadata$peak_total_process_tree_rss_kb <- NA_real_
metadata$peak_temporary_disk_bytes <- NA_real_
jsonlite::write_json(metadata, file.path(result_dir, "metadata.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
saveRDS(list(timings = timings, resource = resource, equivalent = equivalent), file.path(result_dir, "profile.rds"))
cat("PROFILE_RESULT_DIR=", result_dir, "\n", sep = "")
print(timings); cat("sequential_parallel_equivalent=", equivalent, "\n", sep = "")
