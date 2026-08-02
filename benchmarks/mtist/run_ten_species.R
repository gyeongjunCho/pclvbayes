args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
repo <- normalizePath(file.path(script_dir, "../.."))
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "score_mtist.R"))
source(file.path(script_dir, "ten_species_helpers.R"))
config_path <- Sys.getenv("PCLV_MTIST_10_CONFIG",
  file.path(script_dir, "configs", "ten_species_smoke.R"))
config <- source(config_path)$value
root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
selected <- select_mtist_10species(root)
if (!identical(as.integer(config$dataset_id), as.integer(selected$did[[1L]])))
  stop("Configuration does not match deterministic 10-species selection.")
study <- load_mtist_study(config$dataset_id, root)
if (length(study$taxa) != 10L) stop("Ten-species benchmark requires exactly 10 taxa.")

result_dir <- file.path(script_dir, "results", config$result_label)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(result_dir, ".monitor-running"))
cat("pcLVbayes 10-species MTIST benchmark\n")
cat("config:", normalizePath(config_path), " dataset:", config$dataset_id, "\n")

benchmark_library <- tempfile("pclv-ten-species-library-")
dir.create(benchmark_library)
on.exit(unlink(benchmark_library, recursive = TRUE), add = TRUE)
withr::with_libpaths(benchmark_library, action = "prefix",
  devtools::install(repo, quiet = TRUE, upgrade = FALSE,
                    dependencies = FALSE, keep_source = TRUE))
.libPaths(c(benchmark_library, .libPaths()))
library(pclvbayes)
model <- pclvbayes:::get_pclv_model(quiet = TRUE)
exe <- model$exe_file()
exe_info_before <- file.info(exe)
owned_root <- file.path(result_dir, "cmdstan-owned")
dir.create(owned_root, showWarnings = FALSE)
old_options <- options(glvpair.output_root = owned_root)
on.exit(options(old_options), add = TRUE)

dir_snapshot <- function(path) {
  entries <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                        include.dirs = TRUE, no.. = TRUE)
  info <- if (length(entries)) file.info(entries) else data.frame(size = numeric(), isdir = logical())
  c(files = sum(!info$isdir), directories = sum(info$isdir),
    bytes = sum(info$size[!info$isdir], na.rm = TRUE))
}
proc_rss <- function(pid) {
  path <- sprintf("/proc/%s/status", pid)
  if (!file.exists(path)) return(NA_real_)
  line <- grep("^VmRSS:", readLines(path, warn = FALSE), value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(sub("^VmRSS:\\s*([0-9]+).*", "\\1", line[[1L]]))
}
baseline <- dir_snapshot(owned_root)
parent_rss_before <- proc_rss(Sys.getpid())
started <- Sys.time()
fit <- fit_pclv_bayes(
  study$physeq, "subject", "time", taxa_vec = study$taxa,
  chains = config$chains, iter_warmup = config$iter_warmup,
  iter_sampling = config$iter_sampling, seed = config$seed,
  n_workers_outer = config$n_workers_outer, n_workers_kfold = config$n_workers_kfold,
  kfold_K = config$kfold_K, kfold_R = config$kfold_R,
 )
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
parent_rss_after <- proc_rss(Sys.getpid())
status_lines <- readLines(sprintf("/proc/%s/status", Sys.getpid()), warn = FALSE)
parent_hwm <- suppressWarnings(as.numeric(sub("^VmHWM:\\s*([0-9]+).*", "\\1", grep("^VmHWM:", status_lines, value = TRUE)[1L])))
resources <- list(peak_process_tree_rss_kb = NA_real_, peak_descendants = NA_integer_,
                  peak_owned_disk_bytes = NA_real_, peak_owned_files = NA_integer_)
final_owned <- dir_snapshot(owned_root)
if (inherits(fit, "pclv_failure") || is.null(fit$raw)) stop("Top-level benchmark fit failed.")

directions <- expand_mtist_directions(fit, study$truth, config$dataset_id,
                                      config$seed, config$conservative_threshold)
if (nrow(fit$raw) != 45L || nrow(directions) != 90L) stop("Pair/direction count invariant failed.")
matrices <- build_mtist_benchmark_matrices(directions, study$taxa,
                                           config$conservative_threshold)
diagonal <- aggregate_mtist_diagonal(directions, study$taxa)
for (r in seq_len(nrow(diagonal))) {
  tx <- diagonal$taxon[[r]]
  matrices$posterior[tx, tx] <- diagonal$aggregated_value[[r]]
  matrices$conservative[tx, tx] <- diagonal$aggregated_value[[r]]
}
validate_mtist_prediction(matrices$posterior, study$taxa)
validate_mtist_prediction(matrices$conservative, study$taxa)
scores <- data.frame(metric = c("posterior_es_with_diagonal", "posterior_es_without_diagonal",
  "conservative_es_with_diagonal", "conservative_es_without_diagonal"), value = c(
  official_mtist_es(study$truth, matrices$posterior, root, FALSE),
  official_mtist_es(study$truth, matrices$posterior, root, TRUE),
  official_mtist_es(study$truth, matrices$conservative, root, FALSE),
  official_mtist_es(study$truth, matrices$conservative, root, TRUE)))
coverage <- coverage_sign_metrics(directions)
stages <- stage_counts(directions)
asymmetry <- direction_asymmetry(directions)
classes <- as.data.frame(table(directions$diagnostic_class), stringsAsFactors = FALSE)
names(classes) <- c("diagnostic_class", "count")

write_matrix <- function(x, file) utils::write.table(x, file.path(result_dir, file), sep = ",",
  row.names = FALSE, col.names = FALSE, quote = FALSE, na = "NA")
write_matrix(matrices$posterior, "posterior_mean_matrix.csv")
write_matrix(matrices$conservative, "conservative_matrix.csv")
write_matrix(study$truth, "truth_matrix.csv")
write_tsv <- function(x, file) utils::write.table(x, file.path(result_dir, file), sep = "\t",
  row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(directions, "direction_results.tsv"); write_tsv(matrices$masks, "direction_masks.tsv")
write_tsv(diagonal, "diagonal_summary.tsv"); write_tsv(scores, "official_scores.tsv")
write_tsv(coverage, "coverage_metrics.tsv"); write_tsv(stages, "stage_counts.tsv")
write_tsv(classes, "diagnostic_classes.tsv"); write_tsv(asymmetry, "direction_asymmetry.tsv")

exe_info_after <- file.info(exe)
config$selected_by <- "noise=0.01, even sampling, >=10 series, >=15 timepoints; lowest dataset ID"
config$elapsed_wall_seconds <- elapsed
config$taxon_order <- study$taxa
config$truth_id <- study$record$ground_truth[[1L]]
config$pair_tasks <- 45L; config$directed_interactions <- 90L
jsonlite::write_json(config, file.path(result_dir, "config.json"), pretty = TRUE, auto_unbox = TRUE)
environment <- list(
  pcLVbayes_sha = system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
  mtist_sha = mtist_git_sha(root), mtist_root = root,
  start_timestamp = format(started, tz = "UTC", usetz = TRUE),
  end_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R = R.version.string, cmdstan = as.character(cmdstanr::cmdstan_version()),
  os = paste(Sys.info()[c("sysname", "release")], collapse = " "),
  logical_cores = parallel::detectCores(logical = TRUE),
  executable = exe, executable_mtime_before = as.numeric(exe_info_before$mtime),
  executable_mtime_after = as.numeric(exe_info_after$mtime),
  executable_reused = identical(as.numeric(exe_info_before$mtime), as.numeric(exe_info_after$mtime)),
  result_object_bytes = as.numeric(object.size(fit)),
  owned_baseline = as.list(baseline), owned_final = as.list(final_owned),
  parent_rss_before_kb = parent_rss_before, parent_rss_after_kb = parent_rss_after,
  parent_hwm_kb = parent_hwm,
  peak_process_tree_rss_kb = resources$peak_process_tree_rss_kb %||% NA_real_,
  peak_descendants = resources$peak_descendants %||% NA_integer_,
  peak_owned_disk_bytes = resources$peak_owned_disk_bytes %||% NA_real_,
  peak_owned_files = resources$peak_owned_files %||% NA_integer_)
cpu <- tryCatch(readLines("/proc/cpuinfo", warn = FALSE), error = function(e) character())
environment$cpu_model <- sub("^[^:]+: *", "", grep("^model name", cpu, value = TRUE)[1L] %||% NA_character_)
mem <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character())
environment$total_ram_kb <- suppressWarnings(as.numeric(sub("^MemTotal:\\s*([0-9]+).*", "\\1", grep("^MemTotal:", mem, value = TRUE)[1L])))
jsonlite::write_json(environment, file.path(result_dir, "environment.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
console_summary <- c(
  "pcLVbayes 10-species MTIST benchmark",
  paste0("result_dir=", result_dir), paste0("elapsed_seconds=", elapsed),
  capture.output(print(classes)), capture.output(print(stages)),
  capture.output(print(scores)), capture.output(print(coverage)))
writeLines(console_summary, file.path(result_dir, "console.log"))
cat("RESULT_DIR=", result_dir, "\n", sep = "")
cat("elapsed_seconds=", elapsed, " pairs=", nrow(fit$raw), " directions=", nrow(directions), "\n", sep = "")
print(classes); print(stages); print(scores); print(coverage)
