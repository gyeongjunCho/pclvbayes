args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
repo <- normalizePath(file.path(script_dir, "../.."))
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "ten_species_helpers.R"))
source(file.path(script_dir, "v021_truth_isolation.R"))
`%||%` <- function(x, fallback) if (is.null(x)) fallback else x
config <- source(file.path(script_dir, "configs", "ten_species_confirmation.R"))$value
result_dir <- file.path(script_dir, "results", config$result_label)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
stage_b_path <- file.path(script_dir, "results", config$stage_b_result_label, "direction_results.tsv")
stage_b_config <- jsonlite::read_json(file.path(dirname(stage_b_path), "config.json"), simplifyVector = TRUE)
stage_b <- utils::read.delim(stage_b_path, check.names = FALSE, stringsAsFactors = FALSE)
evaluation_targets <- select_confirmation_directions(stage_b)
fit_targets <- lapply(seq_len(nrow(evaluation_targets)), function(i)
  build_v021_confirmation_target(evaluation_targets[i, , drop = FALSE]))
targets <- do.call(rbind, fit_targets)
rm(stage_b)
expected <- data.frame(
  source = c("species_1", "species_2", "species_3", "species_7", "species_8", "species_5"),
  target = c("species_7", "species_9", "species_7", "species_4", "species_4", "species_7"))
if (nrow(targets) != 6L || !setequal(paste(targets$source, targets$target),
                                     paste(expected$source, expected$target)))
  stop("Retained Stage B converged directions do not match the reviewed six.")
if (length(unique(targets$task_id)) != 6L)
  stop("Confirmation requires six distinct unordered pairs.")

stage_b_chain_iterations <- 90 * stage_b_config$chains *
  (stage_b_config$iter_warmup + stage_b_config$iter_sampling) +
  6 * stage_b_config$kfold_K * stage_b_config$chains *
  (stage_b_config$iter_warmup + stage_b_config$iter_sampling)
long_chain_iterations <- 6 * config$chains * (config$iter_warmup + config$iter_sampling)
projected_seconds <- stage_b_config$elapsed_wall_seconds *
  long_chain_iterations / stage_b_chain_iterations
conservative_seconds <- 3 * projected_seconds
projected_temp_bytes <- 85 * 1024^2
cat("TARGETS\n")
print(targets[c("target", "source", "task_id", "direction_index", "seed")], row.names = FALSE)
cat(sprintf("distinct_pairs=6 outer_workers=%d chains_per_direction=%d maximum_simultaneous_chains=%d\n",
            config$n_workers_outer, config$chains,
            config$n_workers_outer * config$chains))
cat(sprintf("projection_point_seconds=%.3f projection_conservative_seconds=%.3f projected_temp_mib=%.1f\n",
            projected_seconds, conservative_seconds, projected_temp_bytes / 1024^2))
if (conservative_seconds > config$maximum_projected_seconds ||
    projected_temp_bytes > config$maximum_projected_temp_bytes)
  stop("Projected resource ceiling exceeded.")

benchmark_library <- tempfile("pclv-confirm-library-")
dir.create(benchmark_library)
on.exit(unlink(benchmark_library, recursive = TRUE), add = TRUE)
withr::with_libpaths(benchmark_library, action = "prefix",
  devtools::install(repo, quiet = TRUE, upgrade = FALSE,
                    dependencies = FALSE, keep_source = TRUE))
.libPaths(c(benchmark_library, .libPaths()))
library(pclvbayes)
root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
loaded_study <- load_mtist_study(config$dataset_id, root)
inference_study <- list(
  physeq = .v021_copy(loaded_study$physeq),
  taxa = as.character(.v021_copy(loaded_study$taxa))
)
rm(loaded_study)

fit_formals <- formals(fit_pclv_bayes)
control_names <- setdiff(names(fit_formals), c("physeq", "subject_col", "time_col", "taxa_vec"))
evaluation <- new.env(parent = environment(fit_pclv_bayes))
controls <- list()
for (name in control_names) {
  value <- eval(fit_formals[[name]], evaluation)
  controls[[name]] <- value
  assign(name, value, evaluation)
}
overrides <- list(
  chains = config$chains, iter_warmup = config$iter_warmup,
  iter_sampling = config$iter_sampling, seed = config$seed,
  progress = "none",
  n_workers_outer = config$n_workers_outer, n_workers_kfold = 1L,
  kfold_K = 5L, kfold_R = 1L
)
controls[names(overrides)] <- overrides
validated <- pclvbayes:::.validate_fit_pclv_inputs(
  inference_study$physeq, "subject", "time", inference_study$taxa, controls)
runtime <- pclvbayes:::.prepare_fit_runtime(validated)
if (inherits(runtime, "pclv_failure")) stop(runtime$reason)
ctx <- runtime$ctx
ctx$n_workers_kfold_eff <- 1L
exe <- ctx$mod_exe_file
exe_before <- file.info(exe)
approved_runtime_context <- build_v021_runtime_context(ctx)

inference_config <- controls[intersect(names(controls), v021_inference_config_fields)]
inference_specs <- lapply(fit_targets, function(target)
  build_v021_inference_spec(inference_study, config = inference_config, task = target))
fit_jobs <- lapply(inference_specs, make_v021_confirmation_fit_closure,
  runtime_context = approved_runtime_context,
  fit_direction = pclvbayes:::.fit_direction_main_posterior)

snapshot <- function(path) {
  entries <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, include.dirs = TRUE, no.. = TRUE) else character()
  info <- if (length(entries)) file.info(entries) else data.frame(size = numeric(), isdir = logical())
  c(files = sum(!info$isdir), directories = sum(info$isdir),
    bytes = sum(info$size[!info$isdir], na.rm = TRUE))
}
proc_value <- function(name) {
  line <- grep(paste0("^", name, ":"), readLines(sprintf("/proc/%s/status", Sys.getpid())),
               value = TRUE)
  as.numeric(sub("^[^:]+:\\s*([0-9]+).*", "\\1", line[[1L]]))
}
owned_root <- file.path(result_dir, "cmdstan-owned")
dir.create(owned_root, showWarnings = FALSE)
old_options <- options(glvpair.output_root = owned_root)
on.exit(options(old_options), add = TRUE)
owned_before <- snapshot(owned_root)
cache_root <- tools::R_user_dir("glvpair", which = "cache")
cache_before <- snapshot(cache_root)
rss_before <- proc_value("VmRSS")
started <- Sys.time()

future::plan(future::multisession, workers = config$n_workers_outer)
on.exit(future::plan(future::sequential), add = TRUE)
results <- furrr::future_map(fit_jobs, function(job) job(),
  .options = furrr::furrr_options(seed = FALSE, packages = "pclvbayes"))
future::plan(future::sequential)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
rss_after <- proc_value("VmRSS")
hwm <- proc_value("VmHWM")
owned_after <- snapshot(owned_root)
cache_after <- snapshot(cache_root)
exe_after <- file.info(exe)

collapse_named <- function(x) {
  x <- unlist(x, use.names = TRUE)
  if (!length(x)) return(NA_character_)
  paste(paste0(names(x), "=", format(as.numeric(x), digits = 10)), collapse = ";")
}
collapse_reasons <- function(x) paste(unlist(x, use.names = FALSE), collapse = ";")
residual_text <- function(x, parameter) {
  d <- x$chain_residual_summary[[1L]]
  d <- d[d$parameter == parameter, , drop = FALSE]
  if (!nrow(d)) return(NA_character_)
  paste(apply(d[c("chain", "mean", "median", "q05", "q95")], 1, paste, collapse = ":"), collapse = ";")
}
rows <- lapply(seq_along(results), function(i) {
  b <- evaluation_targets[i, , drop = FALSE]
  x <- results[[i]]
  if (inherits(x, "pclv_failure")) x <- pclvbayes:::.failed_direction_result(x)
  long_sign <- if (is.finite(x$a_mean)) sign(x$a_mean) else NA_real_
  stage_sign <- sign(b$posterior_mean[[1L]])
  data.frame(
    task_id = b$task_id, direction_index = b$direction_index,
    target = b$target, source = b$source, seed = b$seed,
    truth_coefficient = b$truth_coefficient, truth_sign = b$truth_sign,
    stage_b_class = b$diagnostic_class,
    stage_b_interaction_identifiable = b$interaction_identifiable,
    stage_b_residual_identifiable = b$residual_identifiable,
    stage_b_mean = b$posterior_mean, stage_b_median = b$posterior_median,
    stage_b_sign_probability = b$posterior_sign_probability,
    stage_b_lfsr = b$lfsr, stage_b_sign = stage_sign,
    long_class = x$diagnostic_class,
    long_interaction_identifiable = x$interaction_identifiable,
    long_residual_identifiable = x$residual_identifiable,
    long_mean = x$a_mean, long_median = x$a_median %||% NA_real_,
    long_positive_probability = x$positive_sign_probability %||% NA_real_,
    long_negative_probability = x$negative_sign_probability %||% NA_real_,
    long_lfsr = x$lfsr %||% NA_real_, long_sign = long_sign,
    sign_correct = is.finite(long_sign) && long_sign == b$truth_sign,
    chain_means = collapse_named(x$chain_aij_means),
    chain_medians = collapse_named(x$chain_aij_medians),
    chain_positive_probabilities = collapse_named(x$chain_aij_positive_probabilities),
    chain_negative_probabilities = collapse_named(x$chain_aij_negative_probabilities),
    chain_sign_agreement = x$chain_sign_agreement,
    rhat = x$diag$worst_rhat, ess_bulk = x$diag$min_ess_bulk,
    ess_tail = x$diag$min_ess_tail, divergences = x$diag$n_divergent,
    treedepth_hits = x$diag$n_treedepth_hit,
    chain_ebfmi = collapse_named(x$diag$ebfmi_chain),
    sigma_summaries = residual_text(x, "sigma"),
    sd_ou_summaries = residual_text(x, "sd_ou"),
    phi_summaries = residual_text(x, "phi"),
    lambda_summaries = residual_text(x, "lambda"),
    nu_summaries = residual_text(x, "nu"),
    residual_regime_disagreement = x$residual_regime_disagreement,
    diagnostic_reasons = collapse_reasons(x$indeterminate_reason),
    outcome = confirmation_outcome("converged", stage_sign, x$diagnostic_class, long_sign),
    kfold_attempted = FALSE, elpd_calculated = FALSE, stacking_calculated = FALSE,
    stringsAsFactors = FALSE)
})
comparison <- do.call(rbind, rows)
utils::write.table(comparison, file.path(result_dir, "confirmation_comparison.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
metrics <- data.frame(
  metric = c("prespecified", "stage_b_correct", "retained_stage_b_sign", "changed_sign",
             "still_converged", "downgraded", "long_correct_all_six",
             "long_correct_still_converged"),
  numerator = c(6L, 3L, sum(comparison$long_sign == comparison$stage_b_sign, na.rm = TRUE),
    sum(comparison$long_sign != comparison$stage_b_sign, na.rm = TRUE),
    sum(comparison$long_class == "converged"), sum(comparison$long_class != "converged"),
    sum(comparison$sign_correct), sum(comparison$sign_correct[comparison$long_class == "converged"])),
  denominator = c(6L, 6L, 6L, 6L, 6L, 6L, 6L, sum(comparison$long_class == "converged")))
utils::write.table(metrics, file.path(result_dir, "confirmation_metrics.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
environment <- list(
  elapsed_wall_seconds = elapsed, projected_seconds = projected_seconds,
  conservative_projected_seconds = conservative_seconds,
  projected_temp_bytes = projected_temp_bytes,
  parent_rss_before_kb = rss_before, parent_rss_after_kb = rss_after,
  parent_hwm_kb = hwm, result_object_bytes = as.numeric(object.size(results)),
  owned_before = as.list(owned_before), owned_after = as.list(owned_after),
  cache_before = as.list(cache_before), cache_after = as.list(cache_after),
  executable = exe, executable_mtime_before = as.numeric(exe_before$mtime),
  executable_mtime_after = as.numeric(exe_after$mtime),
  executable_reused = identical(as.numeric(exe_before$mtime), as.numeric(exe_after$mtime)),
  recompilation_occurred = !identical(as.numeric(exe_before$mtime), as.numeric(exe_after$mtime)),
  maximum_simultaneous_chains = config$n_workers_outer * config$chains,
  kfold_invocations = 0L, elpd_calculations = 0L, stacking_calculations = 0L)
jsonlite::write_json(environment, file.path(result_dir, "environment.json"),
                     pretty = TRUE, auto_unbox = TRUE, na = "null")
jsonlite::write_json(config, file.path(result_dir, "config.json"),
                     pretty = TRUE, auto_unbox = TRUE)
cat(sprintf("RESULT_DIR=%s\nactual_seconds=%.3f\n", result_dir, elapsed))
print(comparison[c("source", "target", "stage_b_class", "long_class", "outcome",
                   "stage_b_sign", "long_sign", "truth_sign", "sign_correct")], row.names = FALSE)
print(metrics, row.names = FALSE)
print(environment)
