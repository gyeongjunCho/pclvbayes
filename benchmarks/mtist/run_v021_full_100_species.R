args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg))
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
repo <- normalizePath(file.path(script_dir, "../.."))
pkgload::load_all(repo, quiet = TRUE)
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "ten_species_helpers.R"))
source(file.path(script_dir, "v021_truth_isolation.R"))
source(file.path(script_dir, "v021_resource_policy.R"))
source(file.path(script_dir, "v021_checkpoint_manifest.R"))
source(file.path(script_dir, "v021_diagnostic_features.R"))
source(file.path(script_dir, "v021_four_chain_preflight.R"))
source(file.path(script_dir, "v021_full_100_species.R"))

config_path <- Sys.getenv(
  "PCLV_V021_FULL_CONFIG",
  file.path(script_dir, "configs", "v021_full_100_species.R"))
config <- source(config_path)$value
validate_v021_full_config(config)
paths <- v021_full_paths(config)

existing <- if (dir.exists(paths$output_root))
  list.files(paths$output_root, all.files = TRUE, no.. = TRUE) else character()
recognized_launch <- grepl("^v021_full_100_species_[0-9]{8}-[0-9]{6}\\.log$|^v021_full_100_species\\.pid$",
                           existing)
if (length(existing) && !all(recognized_launch) &&
    !all(c("config.rds", "manifest.rds") %in% existing))
  stop("V021-06 output root contains unrecognized or conflicting content.")
dir.create(paths$output_root, recursive = TRUE, showWarnings = FALSE)
for (path in c(paths$checkpoint_root, paths$monitor_root,
               file.path(paths$output_root, "inference_artifacts"),
               file.path(paths$output_root, "feature_records"),
               file.path(paths$output_root, "worker_registries"),
               file.path(paths$output_root, "cmdstan-owned")))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
validate_v021_resource_policy(policy)
derivation <- derive_safe_outer_concurrency(
  policy, build_v021_operation_spec("main_fit", 3L))
if (!identical(policy$policy_schema, "v021_resource_policy_v2") ||
    !identical(derivation$projected_active_cmdstan_chains, 12L) ||
    !identical(derivation$projected_active_cmdstan_processes, 12L))
  stop("V021-06 requires the verified policy-v2 12-chain/12-slot contract.")
thread_environment <- v021_single_thread_environment()
validate_v021_single_thread_environment(thread_environment)

launch_record <- list(
  launch_schema = "v021_full_100_species_launch_v1",
  parent_pid = as.integer(Sys.getpid()),
  start_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  output_root = paths$output_root,
  pid_file = paths$pid_file,
  manifest_path = paths$manifest,
  checkpoint_root = paths$checkpoint_root,
  config_path = normalizePath(config_path, mustWork = TRUE),
  resource_policy_schema = policy$policy_schema,
  logical_host_threads = policy$logical_host_threads,
  reserved_host_threads = policy$reserved_host_threads,
  maximum_active_cmdstan_chains = policy$maximum_active_cmdstan_chains,
  maximum_cmdstan_process_slots = policy$maximum_cmdstan_process_slots,
  resume_command = sprintf("PCLV_V021_FULL_CONFIG=%s Rscript %s",
    shQuote(normalizePath(config_path, mustWork = TRUE)),
    shQuote(normalizePath(file.path(script_dir, "run_v021_full_100_species.R"),
                          mustWork = TRUE)))
)
.v021_atomic_save_rds(launch_record, file.path(paths$output_root, "launch_record.rds"))
.v021_atomic_save_rds(policy, file.path(paths$output_root, "resource_policy.rds"))
.v021_atomic_save_rds(thread_environment, file.path(paths$output_root, "thread_environment.rds"))
write_v021_full_status(paths$status, "initializing", config)
.v021_atomic_save_rds(list(state = "initializing", timestamp = launch_record$start_timestamp),
                       file.path(paths$monitor_root, "latest.rds"))

root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
metadata <- .read_indexed_csv(.mtist_paths(root)$metadata)
selected_dataset <- select_v021_100_species_dataset(
  metadata, function(did) file.exists(file.path(.mtist_paths(root)$datasets,
                                                paste0("dataset_", did, ".csv"))))
if (!identical(as.character(selected_dataset$did[[1L]]), config$dataset_id))
  stop("Configured dataset does not match deterministic V021-06 selection.")
observations <- load_v021_preflight_observations(root, as.integer(config$dataset_id))
if (length(observations$taxa) != config$taxa_count)
  stop("V021-06 observational input does not contain exactly 100 taxa.")
tasks <- build_v021_full_task_table(observations$taxa, config$seed, config$dataset_id)
manifest_input <- tasks[c("task_id", "direction_index", "target", "source", "seed")]
expected_manifest <- build_v021_checkpoint_manifest(
  manifest_input, config$chains, paths$checkpoint_root)
tasks$pair_id <- expected_manifest$pair_id
tasks <- tasks[c("dataset_id", "pair_id", "task_id", "direction_index",
                 "target", "source", "seed")]

if (file.exists(paths$manifest)) {
  saved_config <- readRDS(file.path(paths$output_root, "config.rds"))
  validate_v021_full_config(saved_config)
  if (!identical(saved_config, config)) stop("V021-06 restart configuration changed.")
  manifest <- read_v021_checkpoint_manifest(paths$manifest)
  identity <- c("task_id", "pair_id", "direction_id", "direction_index",
                "target", "source", "seed", "chain_seeds", "output_location")
  if (!identical(manifest[identity], expected_manifest[identity]))
    stop("V021-06 restart manifest identity changed.")
  restart <- validate_v021_full_restart_states(manifest)
  manifest <- restart$manifest
  write_v021_manifest_atomic(manifest, paths$manifest)
} else {
  manifest <- expected_manifest
  write_v021_manifest_atomic(manifest, paths$manifest)
  saveRDS(config, file.path(paths$output_root, "config.rds"))
}
.v021_atomic_save_rds(tasks, file.path(paths$output_root, "task_table.rds"))

fit_formals <- formals(fit_pclv_bayes)
control_names <- setdiff(names(fit_formals), c("physeq", "subject_col", "time_col", "taxa_vec"))
evaluation <- new.env(parent = environment(fit_pclv_bayes))
controls <- list()
for (name in control_names) {
  controls[[name]] <- eval(fit_formals[[name]], evaluation)
  assign(name, controls[[name]], evaluation)
}
controls[c(
  "eps", "zero_mode_alr", "minpos_alpha", "minpos_base", "eps_fixed",
  "lib_eps_c", "rest_floor_frac", "smooth_scale", "alr_spline_df",
  "alr_spline_spar", "alr_spline_cv", "metric", "quiet", "progress_every",
  "silent_sampler", "max_retries", "use_pathfinder_init", "pf_num_paths",
  "pf_draws", "pf_history_size", "pf_max_lbfgs_iters", "pf_psis_resample"
)] <- list(
  1e-6, "minpos_time", 0.5, "ij", 1e-6, 0.65, 1.0, "logra", NULL,
  NULL, TRUE, "diag_e", FALSE, 1L, FALSE, 3L, TRUE, 8L, 1000L, 50L,
  200L, TRUE
)
controls[c("chains", "iter_warmup", "iter_sampling", "seed", "progress",
           "n_workers_outer", "n_workers_kfold", "kfold_K", "kfold_R")] <- list(
  config$chains, config$iter_warmup, config$iter_sampling, config$seed,
  "none", 3L, 1L, 5L, 1L)
preexisting_executable <- normalizePath(file.path(repo, "inst", "stan", "pclv"),
                                        mustWork = TRUE)
preexisting_info <- unclass(file.info(preexisting_executable)[c("size", "mtime")])
validated <- pclvbayes:::.validate_fit_pclv_inputs(
  observations$physeq, "subject", "time", observations$taxa, controls)
runtime <- pclvbayes:::.prepare_fit_runtime(validated)
if (inherits(runtime, "pclv_failure")) stop(runtime$reason)
runtime$ctx$n_workers_kfold_eff <- 1L
runtime$ctx$max_retries <- 0L
runtime$ctx$use_pathfinder_init <- FALSE
executable <- normalizePath(runtime$ctx$mod_exe_file, mustWork = TRUE)
executable_info <- unclass(file.info(executable)[c("size", "mtime")])
if (!identical(executable, preexisting_executable) ||
    !identical(executable_info, preexisting_info))
  stop("V021-06 executable reuse was not established before worker launch.")
launch_record$executable_path <- executable
launch_record$worker_compilation_count <- 0L
.v021_atomic_save_rds(launch_record, file.path(paths$output_root, "launch_record.rds"))
approved_runtime <- build_v021_runtime_context(runtime$ctx)
inference_config <- controls[intersect(names(controls), v021_inference_config_fields)]
study <- list(physeq = observations$physeq, taxa = observations$taxa)

peak_chains <- 0L
peak_slots <- 0L
batch_number <- 0L

status_counts <- function() {
  states <- table(factor(manifest$execution_state, levels = v021_checkpoint_states))
  list(completed = unname(states[["completed"]]), running = unname(states[["running"]]),
       failed = unname(states[["failed"]]), skipped = unname(states[["skipped"]]))
}

update_status <- function(state = "running", reason = NA_character_) {
  counts <- status_counts()
  write_v021_full_status(paths$status, state, config, counts$completed, counts$running,
    counts$failed, counts$skipped, peak_chains, peak_slots, reason)
}

mark_running <- function(indices) {
  for (i in indices) {
    record <- .v021_record_from_row(manifest[i, , drop = FALSE], "running")
    write_v021_checkpoint_atomic(record)
    manifest <<- .v021_apply_checkpoint(manifest, i, record)
    write_v021_manifest_atomic(manifest, paths$manifest)
  }
  update_status()
}

store_outcome <- function(i, result, elapsed) {
  if (inherits(result, "try-error") || inherits(result, "pclv_failure")) {
    reason <- if (inherits(result, "pclv_failure"))
      paste(result$stage %||% "fit", result$reason %||% "unknown", sep = ":")
    else as.character(result)
    checkpoint <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "failed", elapsed_time = elapsed,
      failure_reason = reason,
      retry_history = if (inherits(result, "pclv_failure"))
        result$details$retry_history %||% list() else list(),
      pathfinder_used = FALSE)
  } else {
    result$.predictive_context <- NULL
    retained <- v021_preflight_retained_record(tasks[i, , drop = FALSE], result)
    inference_results <- setNames(vector("list", length(v021_inference_result_fields)),
                                  v021_inference_result_fields)
    inference_results$posterior_coefficients <- retained$posterior_mean
    inference_results$posterior_summaries <- retained[c("posterior_mean", "posterior_median",
                                                        "posterior_sd")]
    inference_results$posterior_intervals <- c(lower = retained$posterior_interval_lower,
                                               upper = retained$posterior_interval_upper)
    inference_results$posterior_sign_probabilities <- retained$posterior_sign_probability
    inference_results$significance_decisions <- NA
    inference_results$diagnostic_classes <- retained$diagnostic_class
    inference_results$bayesian_eligibility <- retained$bayesian_eligible
    inference_results$kfold_results <- list(state = "not_executed")
    inference_results$elpd_results <- list(state = "not_executed", value = NA_real_)
    inference_results$stacking_results <- list(state = "not_executed", value = NA_real_)
    inference_results$matrices <- list()
    inference_results$masks <- list()
    inference_results$feature_inputs <- retained
    inference_results$stan_data <- list(N = retained$n_pairs)
    inference_results$execution <- list(state = "completed", elapsed_seconds = elapsed)
    inference_results$status <- list(original_reporting_state = retained$diagnostic_class)
    inference_results$unavailable_values <- list(kfold = NA_real_, elpd = NA_real_,
                                                 stacking = NA_real_)
    inference_results$zero_placeholders <- list()
    artifact <- finalize_v021_inference_artifact(list(
      artifact_schema = v021_retained_fixture_schema, artifact_state = "retained",
      execution_state = "completed", inference_results = inference_results))
    feature <- build_v021_diagnostic_feature_record(
      retained, observations$design,
      list(chains = config$chains, iter_sampling = config$iter_sampling))
    .v021_atomic_save_rds(artifact, file.path(paths$output_root, "inference_artifacts",
      paste0(manifest$direction_id[[i]], ".rds")))
    .v021_atomic_save_rds(feature, file.path(paths$output_root, "feature_records",
      paste0(manifest$direction_id[[i]], ".rds")))
    checkpoint <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "completed", elapsed_time = elapsed,
      retry_history = result$retry_history[[1L]] %||% list(), pathfinder_used = FALSE,
      completion_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
      result = artifact)
  }
  write_v021_checkpoint_atomic(checkpoint)
  manifest <<- .v021_apply_checkpoint(manifest, i, checkpoint)
  write_v021_manifest_atomic(manifest, paths$manifest)
}

run_batch <- function(indices) {
  validate_v021_preflight_launch_capacity(policy, length(indices))
  batch_number <<- batch_number + 1L
  batch_id <- sprintf("batch-%05d", batch_number)
  mark_running(indices)
  specs <- lapply(indices, function(i)
    build_v021_inference_spec(study, inference_config, tasks[i, , drop = FALSE]))
  jobs <- lapply(specs, make_v021_confirmation_fit_closure,
    runtime_context = approved_runtime,
    fit_direction = pclvbayes:::.fit_direction_main_posterior)
  started <- proc.time()[["elapsed"]]
  stop_file <- tempfile("v021-full-monitor-stop-")
  start_file <- tempfile("v021-full-fit-start-")
  monitor_start_file <- tempfile("v021-full-monitor-start-")
  registry_file <- tempfile("v021-full-registry-", fileext = ".rds")
  final_file <- file.path(paths$monitor_root, paste0(batch_id, ".rds"))
  latest_file <- file.path(paths$monitor_root, "latest.rds")
  parent_pid <- Sys.getpid()
  fit_jobs <- lapply(seq_along(indices), function(j) parallel::mcparallel({
    while (!file.exists(start_file)) Sys.sleep(0.01)
    with_v021_single_thread_environment(jobs[[j]])
  }, silent = TRUE))
  fit_pids <- vapply(fit_jobs, function(job) as.integer(job$pid), integer(1))
  fit_registry <- do.call(rbind, lapply(seq_along(fit_pids), function(j)
    register_v021_preflight_worker(
      fit_pids[[j]], parent_pid, "direction_fit_worker", batch_id,
      manifest$direction_id[[indices[[j]]]])))
  stop_owned <- function(snapshot) {
    owned <- fit_pids
    repeat {
      children <- snapshot$pid[snapshot$ppid %in% owned]
      expanded <- unique(c(owned, children))
      if (identical(sort(expanded), sort(owned))) break
      owned <- expanded
    }
    for (pid in rev(owned)) system2("kill", c("-INT", as.character(pid)))
  }
  monitor_job <- parallel::mcparallel({
    while (!file.exists(monitor_start_file)) Sys.sleep(0.01)
    registry <- readRDS(registry_file)
    snapshots <- list()
    repeat {
      snap <- tryCatch(capture_v021_preflight_process_snapshot(parent_pid, executable),
                       error = identity)
      if (inherits(snap, "error")) {
        payload <- list(error = conditionMessage(snap), snapshots = snapshots)
        .v021_atomic_save_rds(payload, latest_file)
        .v021_atomic_save_rds(payload, final_file)
        if (length(snapshots)) stop_owned(tail(snapshots, 1L)[[1L]])
        break
      }
      snapshots[[length(snapshots) + 1L]] <- snap
      monitor <- monitor_v021_process_snapshots(
        list(snap), policy, parent_pid, executable, worker_registry = registry)
      payload <- list(snapshot = snap, latest_monitor = monitor, batch_id = batch_id)
      .v021_atomic_save_rds(payload, latest_file)
      if (!identical(monitor$monitoring_state, "verified") ||
          !identical(monitor$compliance_status, "compliant")) {
        payload$snapshots <- snapshots
        payload$error <- monitor$reason %||% "resource_ceiling_exceeded"
        .v021_atomic_save_rds(payload, latest_file)
        .v021_atomic_save_rds(payload, final_file)
        stop_owned(snap)
        break
      }
      if (file.exists(stop_file)) {
        payload$snapshots <- snapshots
        .v021_atomic_save_rds(payload, final_file)
        break
      }
      Sys.sleep(1)
    }
    TRUE
  }, silent = TRUE)
  monitor_registry <- register_v021_preflight_worker(
    as.integer(monitor_job$pid), parent_pid, "resource_monitor", batch_id,
    paste0(batch_id, "-monitor"))
  registry <- rbind(fit_registry, monitor_registry)
  validate_v021_worker_registry(registry)
  .v021_atomic_save_rds(registry, registry_file)
  .v021_atomic_save_rds(registry, file.path(paths$output_root, "worker_registries",
                                            paste0(batch_id, ".rds")))
  file.create(monitor_start_file)
  file.create(start_file)
  results <- unname(parallel::mccollect(fit_jobs))
  file.create(stop_file)
  parallel::mccollect(monitor_job)
  payload <- readRDS(final_file)
  if (!is.null(payload$error)) stop("Live process monitor failed: ", payload$error)
  batch_monitor <- monitor_v021_process_snapshots(
    payload$snapshots, policy, parent_pid, executable, worker_registry = registry)
  validate_v021_preflight_monitor(batch_monitor, policy)
  peak_chains <<- max(peak_chains, batch_monitor$observed_peak_active_cmdstan_chains)
  peak_slots <<- max(peak_slots, batch_monitor$observed_peak_active_cmdstan_processes)
  elapsed <- (proc.time()[["elapsed"]] - started) / length(indices)
  for (j in seq_along(indices)) store_outcome(indices[[j]], results[[j]], elapsed)
  update_status()
  invisible(TRUE)
}

old_options <- options(glvpair.output_root = file.path(paths$output_root, "cmdstan-owned"))
on.exit(options(old_options), add = TRUE)
update_status("running")

execution_error <- NULL
tryCatch(with_v021_single_thread_environment(function() {
  repeat {
    restart <- validate_v021_full_restart_states(manifest)
    manifest <<- restart$manifest
    write_v021_manifest_atomic(manifest, paths$manifest)
    if (!length(restart$resume_indices)) break
    run_batch(head(restart$resume_indices, config$maximum_simultaneous_fits))
  }
}), error = function(e) execution_error <<- e,
interrupt = function(e) execution_error <<- e)

if (!is.null(execution_error)) {
  update_status("failed", conditionMessage(execution_error))
  .v021_atomic_save_rds(list(
    failure_schema = "v021_full_100_species_failure_v1",
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    reason = conditionMessage(execution_error), manifest_path = paths$manifest,
    monitor_path = file.path(paths$monitor_root, "latest.rds")),
    file.path(paths$output_root, "failure_payload.rds"))
  stop(conditionMessage(execution_error))
}

executable_after <- unclass(file.info(executable)[c("size", "mtime")])
summary <- list(
  summary_schema = "v021_full_100_species_summary_v1",
  state = "completed",
  end_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  task_count = nrow(manifest), pair_count = length(unique(manifest$task_id)),
  execution_states = as.list(table(manifest$execution_state)),
  peak_active_chains = peak_chains, peak_cmdstan_slots = peak_slots,
  executable_path = executable,
  executable_reused = identical(preexisting_info, executable_after),
  worker_compilation_count = 0L,
  truth_absent = TRUE,
  manifest_path = paths$manifest,
  checkpoint_root = paths$checkpoint_root)
.v021_atomic_save_rds(summary, paths$summary)
update_status("completed")
cat("V021-06 completed.\n")
