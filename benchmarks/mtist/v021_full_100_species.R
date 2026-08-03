# Benchmark-only full 100-species execution contract for ROADMAP V021-06.

v021_full_schema <- "v021_full_100_species_v1"
v021_full_status_schema <- "v021_full_100_species_status_v1"
v021_full_selection_rule <- paste(
  "MTIST 100-species metadata; noise=0.01; even sampling; at least 10 series;",
  "at least 15 timepoints; lowest dataset ID; no truth fields inspected"
)

build_v021_full_config <- function(output_root) {
  if (!is.character(output_root) || length(output_root) != 1L || is.na(output_root) ||
      !nzchar(output_root)) stop("output_root must be one path.")
  list(
    benchmark_schema = v021_full_schema,
    dataset_id = "1",
    seed = 20260802L,
    taxa_count = 100L,
    unordered_pair_count = 4950L,
    directed_task_count = 9900L,
    selection_rule = v021_full_selection_rule,
    chains = 4L,
    iter_warmup = 2000L,
    iter_sampling = 2000L,
    nominal_retained_draws = 8000L,
    maximum_simultaneous_fits = 3L,
    run_kfold = FALSE,
    use_pathfinder = FALSE,
    output_root = normalizePath(output_root, mustWork = FALSE)
  )
}

validate_v021_full_config <- function(config) {
  fields <- c(
    "benchmark_schema", "dataset_id", "seed", "taxa_count",
    "unordered_pair_count", "directed_task_count", "selection_rule", "chains",
    "iter_warmup", "iter_sampling", "nominal_retained_draws",
    "maximum_simultaneous_fits", "run_kfold", "use_pathfinder", "output_root"
  )
  if (!is.list(config) || !identical(names(config), fields) ||
      !identical(config$benchmark_schema, v021_full_schema) ||
      !identical(config$dataset_id, "1") || !identical(config$seed, 20260802L) ||
      !identical(config$taxa_count, 100L) ||
      !identical(config$unordered_pair_count, 4950L) ||
      !identical(config$directed_task_count, 9900L) ||
      !identical(config$selection_rule, v021_full_selection_rule) ||
      !identical(config$chains, 4L) || !identical(config$iter_warmup, 2000L) ||
      !identical(config$iter_sampling, 2000L) ||
      !identical(config$nominal_retained_draws, 8000L) ||
      !identical(config$maximum_simultaneous_fits, 3L) ||
      !identical(config$run_kfold, FALSE) || !identical(config$use_pathfinder, FALSE))
    stop("Invalid V021-06 full benchmark configuration.")
  invisible(TRUE)
}

select_v021_100_species_dataset <- function(metadata, dataset_exists) {
  required <- c("did", "n_species", "noise", "n_timeseries", "n_timepoints",
                "sampling_scheme")
  if (!is.data.frame(metadata) || !all(required %in% names(metadata)))
    stop("MTIST metadata lacks required truth-free selection fields.")
  forbidden <- names(metadata)[.v021_truth_name(names(metadata))]
  if (length(forbidden)) metadata <- metadata[setdiff(names(metadata), forbidden)]
  eligible <- metadata$n_species == 100L & metadata$noise == 0.01 &
    metadata$sampling_scheme == "even" & metadata$n_timeseries >= 10L &
    metadata$n_timepoints >= 15L & vapply(metadata$did, dataset_exists, logical(1))
  if (!any(eligible)) stop("No eligible local 100-species observational dataset exists.")
  metadata[which(eligible)[order(metadata$did[eligible])][[1L]], required, drop = FALSE]
}

build_v021_full_task_table <- function(taxa, seed, dataset_id = "1") {
  if (!is.character(taxa) || length(taxa) != 100L || anyNA(taxa) ||
      any(!nzchar(taxa)) || anyDuplicated(taxa)) stop("Exactly 100 unique taxa are required.")
  index <- ten_species_direction_index(taxa, seed)
  if (nrow(index) != 9900L || length(unique(index$task_id)) != 4950L ||
      !identical(index$direction_index, seq_len(9900L)))
    stop("Canonical 100-species direction enumeration failed.")
  index$dataset_id <- as.character(dataset_id)
  index <- index[c("dataset_id", "task_id", "direction_index", "target", "source", "seed")]
  rownames(index) <- NULL
  index
}

v021_full_paths <- function(config, timestamp = format(Sys.time(), "%Y%m%d-%H%M%S")) {
  validate_v021_full_config(config)
  root <- config$output_root
  list(
    output_root = root,
    checkpoint_root = file.path(root, "checkpoints"),
    manifest = file.path(root, "manifest.rds"),
    monitor_root = file.path(root, "monitor"),
    status = file.path(root, "status.rds"),
    summary = file.path(root, "summary.rds"),
    pid_file = file.path(root, "v021_full_100_species.pid"),
    log_file = file.path(root, paste0("v021_full_100_species_", timestamp, ".log"))
  )
}

build_v021_full_launch_command <- function(runner, config_path, log_file, pid_file) {
  values <- c(runner, config_path, log_file, pid_file)
  if (anyNA(values) || any(!nzchar(values))) stop("Launch paths must be non-empty.")
  sprintf("PCLV_V021_FULL_CONFIG=%s nohup Rscript %s > %s 2>&1 < /dev/null & echo $! > %s",
          shQuote(config_path), shQuote(runner), shQuote(log_file), shQuote(pid_file))
}

write_v021_full_status <- function(path, state, config, completed = 0L,
                                    running = 0L, failed = 0L, skipped = 0L,
                                    peak_chains = 0L, peak_slots = 0L,
                                    reason = NA_character_) {
  record <- list(
    status_schema = v021_full_status_schema,
    state = state,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    parent_pid = as.integer(Sys.getpid()),
    output_root = config$output_root,
    completed = as.integer(completed), running = as.integer(running),
    failed = as.integer(failed), skipped = as.integer(skipped),
    peak_chains = as.integer(peak_chains), peak_slots = as.integer(peak_slots),
    reason = reason
  )
  .v021_atomic_save_rds(record, path)
  invisible(record)
}

validate_v021_full_restart_states <- function(manifest) {
  plan <- plan_v021_checkpoint_restart(manifest)
  terminal <- plan$manifest$execution_state %in% c("completed", "failed", "skipped")
  if (any(plan$resume_indices %in% which(terminal)))
    stop("Terminal tasks cannot enter the V021-06 resume plan.")
  plan
}
