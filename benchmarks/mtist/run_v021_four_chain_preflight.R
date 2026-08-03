args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
repo <- normalizePath(file.path(script_dir, "../.."))
pkgload::load_all(repo, quiet = TRUE)
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "ten_species_helpers.R"))
source(file.path(script_dir, "v021_truth_isolation.R"))
source(file.path(script_dir, "v021_resource_policy.R"))
source(file.path(script_dir, "v021_checkpoint_manifest.R"))
source(file.path(script_dir, "v021_diagnostic_features.R"))
source(file.path(script_dir, "v021_four_chain_preflight.R"))

config <- source(file.path(script_dir, "configs", "v021_four_chain_preflight.R"))$value
validate_v021_four_chain_preflight_config(config)
if (dir.exists(config$output_root) && length(list.files(config$output_root, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite an existing non-empty preflight output root: ", config$output_root)
dir.create(config$output_root, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(config$output_root, "checkpoints")
artifact_dir <- file.path(config$output_root, "inference_artifacts")
feature_dir <- file.path(config$output_root, "feature_records")
dir.create(checkpoint_dir); dir.create(artifact_dir); dir.create(feature_dir)
monitor_dir <- file.path(config$output_root, "monitor")
dir.create(monitor_dir)
manifest_path <- file.path(config$output_root, "manifest.rds")

policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
operation <- build_v021_operation_spec("main_fit", 3L)
derivation <- derive_safe_outer_concurrency(policy, operation)
thread_environment <- v021_single_thread_environment()
validate_v021_single_thread_environment(thread_environment)

root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
observations <- load_v021_preflight_observations(root, as.integer(config$dataset_id))
index <- ten_species_direction_index(observations$taxa, config$seed)
selected <- select_v021_preflight_tasks(index, config$selected_direction_count)
manifest_input <- selected[c("task_id", "direction_index", "target", "source", "seed")]
manifest <- build_v021_checkpoint_manifest(manifest_input, config$chains, checkpoint_dir)
selected$pair_id <- manifest$pair_id
selected <- selected[c("dataset_id", "pair_id", "task_id", "direction_index",
                       "target", "source", "seed")]

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
preexisting_executable_info <- unclass(file.info(preexisting_executable)[c("size", "mtime")])
validated <- pclvbayes:::.validate_fit_pclv_inputs(
  observations$physeq, "subject", "time", observations$taxa, controls)
runtime <- pclvbayes:::.prepare_fit_runtime(validated)
if (inherits(runtime, "pclv_failure")) stop(runtime$reason)
runtime$ctx$n_workers_kfold_eff <- 1L
runtime$ctx$max_retries <- 0L
runtime$ctx$use_pathfinder_init <- FALSE
executable <- normalizePath(runtime$ctx$mod_exe_file, mustWork = TRUE)
executable_before <- unclass(file.info(executable)[c("size", "mtime")])
compilation_count <- if (identical(executable, preexisting_executable) &&
                         identical(executable_before, preexisting_executable_info)) 0L else 1L
approved_runtime <- build_v021_runtime_context(runtime$ctx)
inference_config <- controls[intersect(names(controls), v021_inference_config_fields)]
study <- list(physeq = observations$physeq, taxa = observations$taxa)
specs <- lapply(seq_len(nrow(selected)), function(i)
  build_v021_inference_spec(study, inference_config, selected[i, , drop = FALSE]))
jobs <- lapply(specs, make_v021_confirmation_fit_closure,
  runtime_context = approved_runtime,
  fit_direction = pclvbayes:::.fit_direction_main_posterior)
rm(study)

write_v021_manifest_atomic(manifest, manifest_path)
saveRDS(config, file.path(config$output_root, "config.rds"))
saveRDS(selected, file.path(config$output_root, "prospective_selection.rds"))
saveRDS(policy, file.path(config$output_root, "resource_policy.rds"))
saveRDS(thread_environment, file.path(config$output_root, "thread_environment.rds"))

cat("V021-05 BOUNDED PREFLIGHT\n")
print(selected, row.names = FALSE)
cat("selection_rule:", config$selection_rule, "\n")
cat("output_root:", config$output_root, "\n")
cat(sprintf("chains=%d warmup=%d retained=%d nominal_draws=%d max_simultaneous_fits=%d projected_chains=%d\n",
  config$chains, config$iter_warmup, config$iter_sampling, config$nominal_retained_draws,
  config$maximum_simultaneous_fits, derivation$projected_active_cmdstan_chains))
print(thread_environment)
cat("executable:", executable, "\nmanifest:", manifest_path,
    "\ncheckpoints:", checkpoint_dir, "\n")

all_snapshots <- list()
all_worker_registry <- NULL
execution_order <- integer()
preflight_root_pid <- Sys.getpid()

mark_running <- function(indices) {
  for (i in indices) {
    running <- .v021_record_from_row(manifest[i, , drop = FALSE], "running")
    write_v021_checkpoint_atomic(running)
    manifest <<- .v021_apply_checkpoint(manifest, i, running)
    write_v021_manifest_atomic(manifest, manifest_path)
  }
}

store_outcome <- function(i, result, elapsed) {
  execution_order <<- c(execution_order, manifest$direction_index[[i]])
  if (inherits(result, "pclv_failure")) {
    checkpoint <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "failed", elapsed_time = elapsed,
      failure_reason = paste(result$stage %||% "fit", result$reason %||% "unknown", sep = ":"),
      retry_history = result$details$retry_history %||% list(), pathfinder_used = FALSE)
  } else {
    result$.predictive_context <- NULL
    retained <- v021_preflight_retained_record(selected[i, , drop = FALSE], result)
    inference_results <- setNames(vector("list", length(v021_inference_result_fields)),
                                  v021_inference_result_fields)
    inference_results$posterior_coefficients <- retained$posterior_mean
    inference_results$posterior_summaries <- retained[c("posterior_mean", "posterior_median", "posterior_sd")]
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
    inference_results$unavailable_values <- list(kfold = NA_real_, elpd = NA_real_, stacking = NA_real_)
    inference_results$zero_placeholders <- list()
    artifact <- finalize_v021_inference_artifact(list(
      artifact_schema = v021_retained_fixture_schema, artifact_state = "retained",
      execution_state = "completed", inference_results = inference_results))
    feature <- build_v021_diagnostic_feature_record(
      retained, observations$design,
      list(chains = config$chains, iter_sampling = config$iter_sampling))
    saveRDS(artifact, file.path(artifact_dir, paste0(manifest$direction_id[[i]], ".rds")))
    saveRDS(feature, file.path(feature_dir, paste0(manifest$direction_id[[i]], ".rds")))
    checkpoint <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "completed", elapsed_time = elapsed,
      retry_history = result$retry_history[[1L]] %||% list(), pathfinder_used = FALSE,
      completion_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), result = artifact)
  }
  write_v021_checkpoint_atomic(checkpoint)
  manifest <<- .v021_apply_checkpoint(manifest, i, checkpoint)
  write_v021_manifest_atomic(manifest, manifest_path)
}

run_batch <- function(indices, controlled_interrupt = FALSE) {
  validate_v021_preflight_launch_capacity(policy, length(indices))
  mark_running(indices)
  started <- proc.time()[["elapsed"]]
  stop_file <- tempfile("v021-monitor-stop-")
  batch_id <- sprintf("batch-%02d", length(all_snapshots) + 1L)
  snapshot_file <- file.path(monitor_dir, paste0(batch_id, "-final.rds"))
  latest_file <- file.path(monitor_dir, "latest.rds")
  interrupt_file <- file.path(monitor_dir, paste0(batch_id, "-interrupt.rds"))
  start_file <- tempfile("v021-fit-start-")
  monitor_start_file <- tempfile("v021-monitor-start-")
  registry_file <- tempfile("v021-worker-registry-", fileext = ".rds")
  fit_jobs <- lapply(indices, function(i) parallel::mcparallel({
    while (!file.exists(start_file)) Sys.sleep(0.01)
    with_v021_single_thread_environment(jobs[[i]])
  }, silent = TRUE))
  fit_worker_pids <- vapply(fit_jobs, function(job) as.integer(job$pid), integer(1))
  fit_registry <- do.call(rbind, lapply(seq_along(fit_worker_pids), function(j)
    register_v021_preflight_worker(
      fit_worker_pids[[j]], preflight_root_pid, "direction_fit_worker", batch_id,
      manifest$direction_id[[indices[[j]]]])))
  stop_fit_processes <- function(snapshot) {
    owned <- fit_worker_pids
    repeat {
      children <- snapshot$pid[snapshot$ppid %in% owned]
      expanded <- unique(c(owned, children))
      if (identical(sort(expanded), sort(owned))) break
      owned <- expanded
    }
    for (pid in rev(owned)) system2("kill", c("-INT", as.character(pid)))
  }
  controlled_worker_pid <- if (controlled_interrupt) fit_jobs[[1L]]$pid else NA_integer_
  monitor_job <- parallel::mcparallel({
    while (!file.exists(monitor_start_file)) Sys.sleep(0.01)
    batch_worker_registry <- readRDS(registry_file)
    snapshots <- list()
    interruption_sent <- FALSE
    repeat {
      snap <- tryCatch(capture_v021_preflight_process_snapshot(
        root_pid = preflight_root_pid, known_model_executables = executable), error = identity)
      if (inherits(snap, "error")) {
        payload <- list(error = conditionMessage(snap), snapshots = snapshots)
        .v021_atomic_save_rds(payload, latest_file)
        .v021_atomic_save_rds(payload, snapshot_file)
        if (length(snapshots)) stop_fit_processes(tail(snapshots, 1L)[[1L]])
        break
      }
      snapshots[[length(snapshots) + 1L]] <- snap
      classified <- monitor_v021_process_snapshots(
        list(snap), policy, preflight_root_pid, executable,
        worker_registry = batch_worker_registry)
      payload <- list(snapshots = snapshots, latest_monitor = classified)
      .v021_atomic_save_rds(payload, latest_file)
      if (!identical(classified$monitoring_state, "verified") ||
          !identical(classified$compliance_status, "compliant")) {
        payload$error <- classified$reason %||% "resource_ceiling_exceeded"
        .v021_atomic_save_rds(payload, latest_file)
        .v021_atomic_save_rds(payload, snapshot_file)
        stop_fit_processes(snap)
        break
      }
      live_chains <- classified$records$pid[
        classified$records$classification == "cmdstan_chain"]
      if (controlled_interrupt && !interruption_sent && length(live_chains)) {
        for (pid in live_chains) system2("kill", c("-INT", as.character(pid)))
        system2("kill", c("-INT", as.character(controlled_worker_pid)))
        event <- list(timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
                      worker_pid = controlled_worker_pid,
                      chain_pids = as.integer(live_chains),
                      reason = "controlled_preflight_interrupt_during_sampling")
        .v021_atomic_save_rds(event, interrupt_file)
        interruption_sent <- TRUE
      }
      if (file.exists(stop_file)) {
        .v021_atomic_save_rds(payload, snapshot_file)
        break
      }
      Sys.sleep(0.1)
    }
    TRUE
  }, silent = TRUE)
  monitor_registry <- register_v021_preflight_worker(
    as.integer(monitor_job$pid), preflight_root_pid, "resource_monitor", batch_id,
    paste0(batch_id, "-monitor"))
  batch_worker_registry <- rbind(fit_registry, monitor_registry)
  validate_v021_worker_registry(batch_worker_registry)
  combined_registry <- if (is.null(all_worker_registry)) batch_worker_registry else
    rbind(all_worker_registry, batch_worker_registry)
  validate_v021_worker_registry(combined_registry)
  all_worker_registry <<- combined_registry
  .v021_atomic_save_rds(batch_worker_registry, registry_file)
  .v021_atomic_save_rds(all_worker_registry,
                        file.path(config$output_root, "worker_registry.rds"))
  file.create(monitor_start_file)
  file.create(start_file)
  results <- parallel::mccollect(fit_jobs)
  file.create(stop_file)
  parallel::mccollect(monitor_job)
  monitor_payload <- readRDS(snapshot_file)
  if (!is.null(monitor_payload$error)) stop("Live process monitor failed: ", monitor_payload$error)
  all_snapshots <<- c(all_snapshots, monitor_payload$snapshots)
  elapsed <- (proc.time()[["elapsed"]] - started) / length(indices)
  if (controlled_interrupt) {
    if (!file.exists(interrupt_file)) stop("Controlled interruption was not observed during sampling.")
    i <- indices[[1L]]
    event <- readRDS(interrupt_file)
    incomplete <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "incomplete", elapsed_time = elapsed,
      interrupted = TRUE, interruption_reason = event$reason,
      failure_reason = event$reason)
    write_v021_checkpoint_atomic(incomplete)
    manifest <<- .v021_apply_checkpoint(manifest, i, incomplete)
    write_v021_manifest_atomic(manifest, manifest_path)
    return(event)
  }
  for (j in seq_along(indices)) store_outcome(indices[[j]], results[[j]], elapsed)
  invisible(TRUE)
}

old_options <- options(glvpair.output_root = file.path(config$output_root, "cmdstan-owned"))
on.exit(options(old_options), add = TRUE)
execution_error <- NULL
tryCatch(with_v021_single_thread_environment(function() {
  run_batch(c(1L, 2L, 3L))

  # Interrupt a live sampling task, then reconcile and resume only that task.
  i <- config$controlled_interrupt_direction
  interruption_event <- run_batch(i, controlled_interrupt = TRUE)
  completed_before <- vapply(manifest$output_location[manifest$execution_state == "completed"],
                              tools::md5sum, character(1))
  restart <- plan_v021_checkpoint_restart(manifest)
  if (!identical(restart$resume_indices, i)) stop("Restart did not select only the incomplete task.")
  manifest <<- restart$manifest
  run_batch(i)
  completed_after <- vapply(names(completed_before), tools::md5sum, character(1))
  restart_record <<- list(
    controlled_interruption = TRUE, resumed_direction_index = selected$direction_index[[i]],
    interruption_event = interruption_event,
    completed_tasks_not_rerun = identical(completed_before, completed_after),
    incomplete_task_resumed = selected$direction_index[[i]] %in% execution_order,
    manifest_reconciled = identical(readRDS(manifest_path), manifest),
    execution_order = execution_order)
}), error = function(e) execution_error <<- e,
interrupt = function(e) execution_error <<- e)

if (!is.null(execution_error)) {
  orphan_audit <- tryCatch({
    snapshot <- capture_v021_preflight_process_snapshot(Sys.getpid(), executable)
    classified <- classify_v021_process_snapshot(
      snapshot, Sys.getpid(), executable, worker_registry = all_worker_registry)
    live_classes <- c("cmdstan_chain", "pathfinder_process", "cmdstan_diagnostic",
                      "unknown_potential_cmdstan")
    list(passed = !any(classified$classification %in% live_classes),
         active_preflight_cmdstan_processes = sum(classified$classification %in% live_classes))
  }, error = function(e) list(passed = FALSE,
                              active_preflight_cmdstan_processes = NA_integer_,
                              reason = conditionMessage(e)))
  saveRDS(list(
    failure_schema = "v021_four_chain_preflight_failure_v1",
    state = "failed", stage = "execution", reason = conditionMessage(execution_error),
    manifest = manifest, monitor_latest = file.path(monitor_dir, "latest.rds"),
    orphan_process_check = orphan_audit),
    file.path(config$output_root, "failure_payload.rds"))
  stop(conditionMessage(execution_error))
}

monitor <- monitor_v021_process_snapshots(
  all_snapshots, policy, Sys.getpid(), executable,
  worker_registry = all_worker_registry)
validate_v021_preflight_monitor(monitor, policy)
saveRDS(monitor$records, file.path(config$output_root, "process_tree_snapshots.rds"))
saveRDS(restart_record, file.path(config$output_root, "restart_record.rds"))
executable_after <- unclass(file.info(executable)[c("size", "mtime")])
executable_record <- list(
  executable_path = executable, compiled_before_workers = TRUE,
  compilation_count = compilation_count,
  executable_before = executable_before, executable_after = executable_after,
  reuse_verified = identical(executable_before, executable_after),
  worker_side_compilation = FALSE, worker_compilation_count = 0L)
saveRDS(executable_record, file.path(config$output_root, "executable_reuse.rds"))
features <- lapply(list.files(feature_dir, full.names = TRUE), readRDS)
invisible(lapply(features, validate_v021_diagnostic_feature_record))
orphan_snapshot <- capture_v021_preflight_process_snapshot(Sys.getpid(), executable)
orphan_classified <- classify_v021_process_snapshot(
  orphan_snapshot, Sys.getpid(), executable, worker_registry = all_worker_registry)
orphan_process_check <- list(
  passed = !any(orphan_classified$classification %in%
                  c("cmdstan_chain", "pathfinder_process", "cmdstan_diagnostic",
                    "unknown_potential_cmdstan")),
  active_preflight_cmdstan_processes = sum(orphan_classified$classification %in%
    c("cmdstan_chain", "pathfinder_process", "cmdstan_diagnostic",
      "unknown_potential_cmdstan")))
saveRDS(orphan_process_check, file.path(config$output_root, "orphan_process_check.rds"))
summary <- build_v021_preflight_summary(config, selected, policy, monitor, manifest,
                                        executable_record, restart_record, features, TRUE,
                                        orphan_process_check)
validate_v021_preflight_summary(summary)
saveRDS(summary, file.path(config$output_root, "summary.rds"))
saveRDS(list(manifest = manifest, failures = manifest[manifest$execution_state == "failed", ]),
        file.path(config$output_root, "failure_records.rds"))
if (!identical(summary$state, "passed")) stop("V021-05 preflight failed: ",
                                               paste(summary$failure_reasons, collapse = "; "))
cat("V021-05 preflight passed. Peak chains:", summary$observed_peak_active_chains,
    "peak CmdStan process slots:", summary$observed_peak_cmdstan_process_slots, "\n")
