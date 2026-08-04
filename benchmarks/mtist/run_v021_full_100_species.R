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
source(file.path(script_dir, "v021_robustness_features.R"))
source(file.path(script_dir, "v021_four_chain_preflight.R"))
source(file.path(script_dir, "v021_full_100_species.R"))

run_v021_full_100_species <- function() {
config_path <- Sys.getenv(
  "PCLV_V021_FULL_CONFIG",
  file.path(script_dir, "configs", "v021_full_100_species.R"))
config <- source(config_path)$value
validate_v021_full_config(config)
if (!identical(v021_current_cpu_affinity(), config$cpu_affinity))
  stop("V021-06 controller CPU affinity must be ", config$cpu_affinity, ".")

# Optional foreground-only debug mode. When enabled, both variables are required
# and the debug output root must differ from the configured production root.
debug_directions_raw <- trimws(Sys.getenv("PCLV_V021_DEBUG_DIRECTIONS", ""))
debug_output_root_raw <- trimws(Sys.getenv("PCLV_V021_DEBUG_OUTPUT_ROOT", ""))
debug_mode <- nzchar(debug_directions_raw)
production_output_root <- normalizePath(config$output_root, mustWork = FALSE)

if (debug_mode) {
  if (!nzchar(debug_output_root_raw))
    stop("PCLV_V021_DEBUG_OUTPUT_ROOT is required in debug mode.")

  debug_output_root <- normalizePath(debug_output_root_raw, mustWork = FALSE)
  if (identical(debug_output_root, production_output_root))
    stop("Debug output root must differ from the production output root.")

  pieces <- trimws(strsplit(debug_directions_raw, ",", fixed = TRUE)[[1L]])
  if (!length(pieces) || any(!nzchar(pieces)) ||
      any(!grepl("^[0-9]+$", pieces)))
    stop("PCLV_V021_DEBUG_DIRECTIONS must contain comma-separated positive integers.")

  debug_directions <- suppressWarnings(as.integer(pieces))
  if (anyNA(debug_directions) || any(debug_directions < 1L) ||
      anyDuplicated(debug_directions))
    stop("Debug direction indices must be unique positive integers.")

  config$output_root <- debug_output_root
}

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
               file.path(paths$output_root, "sampler_scale_audits"),
               file.path(paths$output_root, "robustness_feature_records"),
               file.path(paths$output_root, "worker_registries"),
               file.path(paths$output_root, "batch_ownership"),
               file.path(paths$output_root, "failure_traces"),
               file.path(paths$output_root, "cleanup_audits"),
               file.path(paths$output_root, "cmdstan-owned")))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

policy <- build_v021_full_resource_policy(config)
validate_v021_resource_policy(policy)
derivation <- derive_safe_outer_concurrency(
  policy,
  build_v021_operation_spec("main_fit", config$maximum_simultaneous_fits))
kfold_derivation <- derive_safe_outer_concurrency(
  policy,
  build_v021_operation_spec(
    "kfold_fit", config$maximum_simultaneous_kfold_fits))
if (!identical(policy$policy_schema, "v021_resource_policy_v5") ||
    !identical(policy$main_chains, 4L) ||
    !identical(policy$main_parallel_chains, 1L) ||
    !identical(policy$kfold_parallel_chains, 1L) ||
    !identical(policy$proposed_outer_concurrency, 12L) ||
    !identical(policy$maximum_concurrent_kfold_fits, 12L) ||
    !identical(derivation$per_fit_simultaneous_chain_slots, 1L) ||
    !identical(derivation$projected_active_cmdstan_chains, 12L) ||
    !identical(derivation$projected_active_cmdstan_processes, 12L) ||
    !identical(kfold_derivation$per_fit_simultaneous_chain_slots, 1L) ||
    !identical(kfold_derivation$projected_active_cmdstan_chains, 12L) ||
    !identical(kfold_derivation$projected_active_cmdstan_processes, 12L))
  stop("V021-06 requires the verified policy-v5 main-and-K-fold 12x1 contract.")
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
  cpu_affinity = config$cpu_affinity,
  resume_command = if (debug_mode) sprintf(
    "PCLV_V021_FULL_CONFIG=%s PCLV_V021_DEBUG_DIRECTIONS=%s PCLV_V021_DEBUG_OUTPUT_ROOT=%s Rscript %s",
    shQuote(normalizePath(config_path, mustWork = TRUE)),
    shQuote(paste(debug_directions, collapse = ",")),
    shQuote(paths$output_root),
    shQuote(normalizePath(file.path(script_dir, "run_v021_full_100_species.R"),
                          mustWork = TRUE))
  ) else sprintf(
    "PCLV_V021_FULL_CONFIG=%s Rscript %s",
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

git_commit <- trimws(system2(
  "git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)[[1L]])
prepared_execution <- prepare_v021_full_execution(
  config, observations$taxa,
  provenance = list(code_commit = git_commit,
                    benchmark_schema = config$benchmark_schema),
  initialize = FALSE)

if (identical(tolower(Sys.getenv("PCLV_V021_MANIFEST_DRY_RUN", "false")), "true")) {
  resource_dry_run <- prepare_v021_resource_dry_run(prepared_execution, policy)
  integration_dry_run <- run_v021_full_dry_run_audit(
    prepared_execution, policy, paths$output_root, config$maximum_task_attempts)
  .v021_atomic_save_rds(
    integration_dry_run, file.path(paths$output_root, "dry_run_audit.rds"))
  cat(sprintf(
    paste0("V021-03/V021-02 dry run: manifest=%s runnable=%d completed=%d ",
           "waves=%d maximum_chains=%d sampling_launched=false\n"),
    prepared_execution$manifest$manifest_hash,
    length(prepared_execution$plan$runnable_indices),
    length(prepared_execution$plan$completed_indices),
    length(resource_dry_run$waves), resource_dry_run$maximum_active_cmdstan_chains))
  quit(save = "no", status = 0L)
}

execution_paths <- initialize_v021_execution_root(
  prepared_execution$manifest, paths$output_root)
execution_plan <- plan_v021_execution_resume(
  prepared_execution$manifest, paths$output_root, config$maximum_task_attempts,
  persist_reconciliation = TRUE)
execution_status <- execution_plan$status
attempt_ledger <- execution_plan$attempt_ledger

controller_ownership <- acquire_v021_controller_ownership(
  paths$output_root, prepared_execution$manifest$manifest_hash,
  prepared_execution$manifest$configuration_hash)
on.exit(release_v021_controller_ownership(paths$output_root,
                                           controller_ownership), add = TRUE)
reservation_path <- file.path(paths$output_root, "chain_reservations.rds")
active_reservations <- if (file.exists(reservation_path))
  reconcile_v021_dead_controller_reservations(
    read_v021_reservations(reservation_path, policy), policy,
    controller_alive = FALSE) else new_v021_reservations()
write_v021_reservations_atomic(active_reservations, reservation_path, policy)

if (debug_mode) {
  missing_directions <- setdiff(debug_directions, tasks$direction_index)
  if (length(missing_directions))
    stop("Unknown debug direction indices: ",
         paste(missing_directions, collapse = ", "))

  # match() preserves the exact user-requested order and canonical identities.
  tasks <- tasks[
    match(debug_directions, tasks$direction_index),
    ,
    drop = FALSE
  ]
  if (!identical(tasks$direction_index, debug_directions))
    stop("Debug direction selection did not preserve canonical indices.")

  message(
    "V021-06 DEBUG MODE: directions ",
    paste(tasks$direction_index, collapse = ", "),
    "; output root = ",
    paths$output_root
  )
}

manifest_input <- tasks[c("task_id", "direction_index", "target", "source", "seed")]
expected_manifest <- build_v021_checkpoint_manifest(
  manifest_input, config$chains, paths$checkpoint_root,
  allow_incomplete_pairs = debug_mode)
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
  "none", config$maximum_simultaneous_fits, 1L, 5L, 1L)
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
launch_record$scheduler <- list(
  strategy = "chain_slot_rolling_window",
  chains_per_direction = as.integer(config$chains),
  parallel_chains_per_direction = as.integer(config$main_parallel_chains),
  maximum_active_main_fits = as.integer(config$maximum_simultaneous_fits),
  predictive_strategy = "kfold_global_fold_slot_rolling",
  predictive_scheduling_unit = "direction_repetition_fold",
  kfold_parallel_chains_per_fold =
    as.integer(config$kfold_parallel_chains),
  maximum_active_kfold_fits =
    as.integer(config$maximum_simultaneous_kfold_fits),
  main_kfold_overlap = FALSE,
  window_size = as.integer(config$rolling_window_size),
  polling_seconds = 0.05)
.v021_atomic_save_rds(launch_record, file.path(paths$output_root, "launch_record.rds"))
approved_runtime <- build_v021_runtime_context(runtime$ctx)
inference_config <- controls[intersect(names(controls), v021_inference_config_fields)]
study <- list(physeq = observations$physeq, taxa = observations$taxa)

durable_peaks <- recover_v021_monitor_peaks(paths$output_root, policy, executable)
peak_chains <- durable_peaks$chains
peak_slots <- durable_peaks$processes
batch_number <- 0L
rolling_window_size <- as.integer(config$rolling_window_size)
trace_root <- file.path(paths$output_root, "failure_traces")
parent_trace <- new_v021_failure_trace_context(trace_root, "parent")

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
    started_attempt <- start_v021_task_attempt(
      execution_status, attempt_ledger, prepared_execution$manifest,
      tasks$direction_index[[i]], format(Sys.time(), tz = "UTC", usetz = TRUE),
      worker_provenance = "v021_full_controller")
    execution_status <<- started_attempt$status
    attempt_ledger <<- started_attempt$ledger
    write_v021_task_status_atomic(
      execution_status, execution_paths$status, prepared_execution$manifest)
    write_v021_attempt_ledger_atomic(
      attempt_ledger, execution_paths$attempts, prepared_execution$manifest)
    record <- .v021_record_from_row(manifest[i, , drop = FALSE], "running")
    write_v021_checkpoint_atomic(record)
    manifest <<- .v021_apply_checkpoint(manifest, i, record)
    write_v021_manifest_atomic(manifest, paths$manifest)
  }
  update_status()
}

store_outcome <- function(i, result, elapsed) {
  set_v021_failure_trace_phase(parent_trace, "fit_return_handling")
  if (inherits(result, "try-error") || inherits(result, "pclv_failure") ||
      inherits(result, "v021_traced_child_error")) {
    reason <- if (inherits(result, "pclv_failure"))
      paste(result$stage %||% "fit", result$reason %||% "unknown", sep = ":")
    else if (inherits(result, "v021_traced_child_error")) result$condition_message
    else as.character(result)
    checkpoint <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "failed", elapsed_time = elapsed,
      failure_reason = reason,
      retry_history = if (inherits(result, "pclv_failure"))
        result$details$retry_history %||% list() else list(),
      pathfinder_used = FALSE)
    finished_attempt <- finish_v021_task_attempt(
      execution_status, attempt_ledger, prepared_execution$manifest,
      tasks$direction_index[[i]], "failed",
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      terminal_reason = reason, failure_class = class(result)[[1L]])
  } else {
    result$.predictive_context <- NULL
    set_v021_failure_trace_phase(parent_trace, "scientific_state_classification")
    retained <- v021_preflight_retained_record(tasks[i, , drop = FALSE], result)
    set_v021_failure_trace_phase(parent_trace, "scientific_state_classification",
                                 completed = TRUE)
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
    inference_results$kfold_results <- if (retained$kfold_attempted)
      list(state = if (retained$kfold_completed) "completed" else "failed",
           folds_ok = retained$kfold_folds_ok, folds_failed = retained$kfold_folds_failed)
    else list(state = "not_eligible", n_successful_test_observations = 0L)
    inference_results$elpd_results <- list(
      state = if (retained$elpd_available) "available" else retained$aggregate_elpd_missing_state,
      value = retained$aggregate_elpd,
      method = "student-t-scale-mixture-kalman-ou-q16")
    inference_results$matrices <- list()
    inference_results$masks <- list()
    inference_results$feature_inputs <- retained
    inference_results$stan_data <- list(N = retained$n_pairs)
    inference_results$execution <- list(state = "completed", elapsed_seconds = elapsed)
    inference_results$status <- list(original_reporting_state = retained$diagnostic_class)
    inference_results$unavailable_values <- list(kfold = NA_real_, elpd = NA_real_)
    inference_results$zero_placeholders <- list()
    set_v021_failure_trace_phase(parent_trace, "artifact_finalization")
    artifact <- finalize_v021_inference_artifact(list(
      artifact_schema = v021_retained_fixture_schema, artifact_state = "retained",
      execution_state = "completed", inference_results = inference_results))
    set_v021_failure_trace_phase(parent_trace, "artifact_finalization", completed = TRUE)
    set_v021_failure_trace_phase(parent_trace, "feature_generation")
    feature <- build_v021_diagnostic_feature_record(
      retained, observations$design,
      list(chains = config$chains, iter_sampling = config$iter_sampling))
    .v021_atomic_save_rds(artifact, file.path(paths$output_root, "inference_artifacts",
      paste0(manifest$direction_id[[i]], ".rds")))
    .v021_atomic_save_rds(feature, file.path(paths$output_root, "feature_records",
      paste0(manifest$direction_id[[i]], ".rds")))
    .v021_atomic_save_rds(result$.accepted_scale_audit,
      file.path(paths$output_root, "sampler_scale_audits",
                paste0(manifest$direction_id[[i]], ".rds")))
    observed_support <- calculate_v021_observed_support_features(
      v021_full_observed_support_input(approved_runtime, tasks[i, , drop = FALSE]))
    posterior_stability <- extract_v021_posterior_stability_features(list(
      finalized = TRUE, terminal_state = "completed",
      posterior_mean = retained$posterior_mean,
      posterior_median = retained$posterior_median,
      posterior_sd = retained$posterior_sd,
      interval_lower = retained$posterior_interval_lower,
      interval_upper = retained$posterior_interval_upper,
      psp = retained$p_sign2, lfsr = retained$lfsr,
      rhat_max = retained$rhat, bulk_ess_min = retained$ess_bulk,
      tail_ess_min = retained$ess_tail,
      divergence_count = retained$divergences,
      treedepth_hit_count = retained$treedepth_hits,
      ebfmi_min = retained$ebfmi_min,
      retry_count = max(length(result$retry_history[[1L]]) - 1L, 0L),
      diagnostic_class = retained$diagnostic_class,
      predictive_eligible = retained$bayesian_eligible,
      elpd = NA_real_, successful_test_observation_count = 0L,
      folds_attempted = 0L, folds_failed = 0L))
    robustness <- build_v021_robustness_feature_artifact(
      prepared_execution$manifest, tasks$direction_index[[i]], observed_support,
      posterior_stability, code_provenance = git_commit)
    write_v021_robustness_feature_artifact_atomic(
      robustness, file.path(paths$output_root, "robustness_feature_records",
                            paste0(manifest$direction_id[[i]], ".rds")),
      prepared_execution$manifest)
    set_v021_failure_trace_phase(parent_trace, "feature_generation", completed = TRUE)
    checkpoint <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "completed", elapsed_time = elapsed,
      retry_history = result$retry_history[[1L]] %||% list(), pathfinder_used = FALSE,
      completion_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
      result = artifact)
    completion <- build_v021_completion_artifact(
      prepared_execution$manifest, tasks$direction_index[[i]], paths$output_root,
      posterior_summary = list(mean = retained$posterior_mean,
                               median = retained$posterior_median,
                               sd = retained$posterior_sd,
                               lower = retained$posterior_interval_lower,
                               upper = retained$posterior_interval_upper),
      diagnostics = list(class = retained$diagnostic_class,
                         rhat = retained$rhat, bulk_ess = retained$ess_bulk,
                         tail_ess = retained$ess_tail,
                         divergences = retained$divergences,
                         treedepth_hits = retained$treedepth_hits),
      psp_lfsr = list(PSP = retained$p_sign2, LFSR = retained$lfsr),
      predictive_eligibility = retained$bayesian_eligible,
      subject_elpd = list(state = if (retained$kfold_attempted)
        if (retained$kfold_completed) "completed" else "failed" else "not_applicable",
        n_successful_test_observations =
          v021_preflight_successful_test_observation_count(result$kfold_success_total),
        method = "student-t-scale-mixture-kalman-ou-q16"),
      seed_split_provenance = list(direction_seed = retained$seed,
                                   kfold_seed = config$kfold_seed))
    completion_path <- file.path(execution_paths$artifact_root,
                                 paste0(manifest$direction_id[[i]], ".rds"))
    write_v021_completion_artifact_atomic(
      completion, completion_path, prepared_execution$manifest, paths$output_root)
    finished_attempt <- finish_v021_task_attempt(
      execution_status, attempt_ledger, prepared_execution$manifest,
      tasks$direction_index[[i]], "completed",
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      artifact_path = completion_path, result_root = paths$output_root)
  }
  execution_status <<- finished_attempt$status
  attempt_ledger <<- finished_attempt$ledger
  write_v021_task_status_atomic(
    execution_status, execution_paths$status, prepared_execution$manifest)
  write_v021_attempt_ledger_atomic(
    attempt_ledger, execution_paths$attempts, prepared_execution$manifest)
  set_v021_failure_trace_phase(parent_trace, "checkpoint_write")
  write_v021_checkpoint_atomic(checkpoint)
  set_v021_failure_trace_phase(parent_trace, "checkpoint_write", completed = TRUE)
  manifest <<- .v021_apply_checkpoint(manifest, i, checkpoint)
  set_v021_failure_trace_phase(parent_trace, "manifest_update")
  write_v021_manifest_atomic(manifest, paths$manifest)
  set_v021_failure_trace_phase(parent_trace, "manifest_update", completed = TRUE)
}

run_batch <- function(indices) {
  if (!length(indices)) return(invisible(TRUE))

  set_v021_failure_trace_phase(parent_trace, "batch_launch")
  validate_v021_controller_ownership(
    controller_ownership, paths$output_root,
    prepared_execution$manifest$manifest_hash,
    prepared_execution$manifest$configuration_hash)

  maximum_active <- as.integer(config$maximum_simultaneous_fits)
  validate_v021_preflight_launch_capacity(
    policy, min(length(indices), maximum_active))

  batch_number <<- batch_number + 1L
  batch_id <- sprintf("batch-%05d", batch_number)
  parent_trace$batch_id <- batch_id
  parent_trace$task_ids <- as.integer(manifest$task_id[indices])
  parent_trace$direction_ids <- as.character(manifest$direction_id[indices])

  specs <- lapply(indices, function(i)
    build_v021_inference_spec(study, inference_config, tasks[i, , drop = FALSE]))
  jobs <- lapply(specs, make_v021_confirmation_fit_closure,
    runtime_context = approved_runtime,
    fit_direction = pclvbayes:::.fit_direction_main_posterior)
  worker_traces <- lapply(seq_along(indices), function(j)
    new_v021_failure_trace_context(
      trace_root, "outer_worker", batch_id,
      manifest$task_id[[indices[[j]]]], manifest$direction_id[[indices[[j]]]]))

  results <- vector("list", length(indices))
  elapsed_by_position <- rep(NA_real_, length(indices))
  stop_file <- tempfile("v021-full-monitor-stop-")
  monitor_start_file <- tempfile("v021-full-monitor-start-")
  registry_file <- tempfile("v021-full-registry-", fileext = ".rds")
  final_file <- file.path(paths$monitor_root, paste0(batch_id, ".rds"))
  latest_file <- file.path(paths$monitor_root, "latest.rds")
  parent_pid <- Sys.getpid()

  fit_job_tracker <- new_v021_child_job_tracker()
  monitor_job_tracker <- new_v021_child_job_tracker()
  predictive_job_tracker <- new_v021_child_job_tracker()
  predictive_monitor_job_tracker <- new_v021_child_job_tracker()
  active_meta <- new.env(hash = TRUE, parent = emptyenv())

  collect_ready_jobs <- function(tracker) {
    if (!is.environment(tracker) || is.null(tracker$jobs))
      stop("Invalid V021 child-job tracker.")

    tracked <- tracker$jobs
    if (!length(tracked)) return(list())

    collected <- parallel::mccollect(tracked, wait = FALSE)
    if (is.null(collected) || !length(collected))
      return(list())
    if (!is.list(collected) || is.null(names(collected)) ||
        anyNA(names(collected)) || any(!nzchar(names(collected))))
      stop("Ready V021 child results must be a PID-named list.")

    tracked_pids <- vapply(
      tracked, function(job) as.integer(job$pid), integer(1))
    collected_pids <- suppressWarnings(as.integer(names(collected)))
    if (anyNA(collected_pids) || anyDuplicated(collected_pids) ||
        any(!collected_pids %in% tracked_pids))
      stop("Ready V021 child results did not match tracked child PIDs.")

    tracker$jobs <- tracked[!tracked_pids %in% collected_pids]
    collected
  }

  reserve_and_mark <- function(j) {
    i <- indices[[j]]
    attempt_number <- v021_task_attempt_number(
      execution_status, tasks$direction_index[[i]], next_attempt = TRUE)
    active_reservations <<- reserve_v021_capacity(
      active_reservations, policy, manifest$direction_id[[i]],
      attempt_number, "main_fit")
    write_v021_reservations_atomic(active_reservations, reservation_path, policy)
    mark_running(i)
  }

  spawn_fit_job <- function(j) {
    gate_file <- tempfile("v021-full-fit-start-")
    job <- parallel::mcparallel({
      while (!file.exists(gate_file)) Sys.sleep(0.01)
      set_v021_failure_trace_phase(worker_traces[[j]], "worker_execution")
      run_v021_traced_child(function() {
        value <- withr::with_options(
          list(pclvbayes.parallel_chains_override =
                 as.integer(config$main_parallel_chains)),
          with_v021_single_thread_environment(jobs[[j]])
        )
        set_v021_failure_trace_phase(
          worker_traces[[j]], "worker_execution", completed = TRUE)
        value
      }, worker_traces[[j]])
    }, silent = TRUE)

    pid <- as.integer(job$pid)
    pid_key <- as.character(pid)
    fit_job_tracker$jobs[[pid_key]] <- job
    active_meta[[pid_key]] <- list(
      position = as.integer(j),
      start_elapsed = proc.time()[["elapsed"]])

    list(
      gate_file = gate_file,
      registry = register_v021_preflight_worker(
        pid, parent_pid, "direction_fit_worker", batch_id,
        manifest$direction_id[[indices[[j]]]])
    )
  }

  initial_count <- min(length(indices), maximum_active)
  initial_spawned <- vector("list", initial_count)
  for (j in seq_len(initial_count)) {
    reserve_and_mark(j)
    initial_spawned[[j]] <- spawn_fit_job(j)
  }
  fit_registry <- do.call(
    rbind, lapply(initial_spawned, function(value) value$registry))

  cleanup_path <- file.path(
    paths$output_root, "cleanup_audits", paste0(batch_id, ".rds"))
  ownership <- new_v021_batch_ownership(
    batch_id, parent_pid, fit_registry, executable, paths$output_root, cleanup_path)
  parent_trace$cleanup_audit_path <- cleanup_path
  batch_exit_condition <- NULL

  registry <- fit_registry
  persist_registry <- function() {
    validate_v021_worker_registry(registry)
    updated_ownership <- ownership
    updated_ownership$worker_registry <- registry
    ownership <<- updated_ownership
    .v021_atomic_save_rds(registry, registry_file)
    .v021_atomic_save_rds(
      registry,
      file.path(
        paths$output_root, "worker_registries", paste0(batch_id, ".rds")))
    .v021_atomic_save_rds(list(
      ownership_schema = "v021_full_batch_ownership_v1",
      batch_id = batch_id,
      parent_pid = as.integer(parent_pid),
      known_model_executable = executable,
      output_root = paths$output_root,
      worker_registry = registry,
      registration_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)),
      file.path(
        paths$output_root, "batch_ownership", paste0(batch_id, ".rds")))
    invisible(TRUE)
  }

  runtime_inventory_reader <- function() {
    snapshot <- capture_v021_preflight_process_snapshot(parent_pid, executable)
    classify_v021_process_snapshot(
      snapshot, parent_pid, executable, worker_registry = ownership$worker_registry)
  }
  runtime_alive_reader <- function(pids, start_times) {
    if (!length(pids)) return(integer())
    keep <- vapply(seq_along(pids), function(k) {
      stat <- tryCatch(
        readLines(
          file.path("/proc", pids[[k]], "stat"), warn = FALSE, n = 1L),
        error = function(e) character())
      if (!length(stat)) return(FALSE)
      parsed <- tryCatch(
        .v021_parse_linux_stat(stat, pids[[k]]), error = identity)
      !inherits(parsed, "error") &&
        identical(parsed$start_time, start_times[[k]]) &&
        !identical(parsed$process_state, "Z")
    }, logical(1))
    as.integer(pids[keep])
  }
  runtime_signaler <- function(pids, signal) {
    option <- switch(
      signal,
      SIGINT = "-INT",
      SIGTERM = "-TERM",
      stop("Unsupported cleanup signal."))
    invisible(lapply(as.integer(pids), function(pid)
      system2(
        "kill", c(option, as.character(pid)),
        stdout = FALSE, stderr = FALSE)))
  }
  runtime_reaper <- function()
    reap_v021_uncollected_jobs(
      fit_job_tracker, monitor_job_tracker,
      predictive_job_tracker, predictive_monitor_job_tracker)

  on.exit({
    cleanup_result <- tryCatch(
      cleanup_v021_owned_batch(
        ownership, batch_exit_condition, "parent", parent_trace$phase,
        runtime_inventory_reader, runtime_signaler, runtime_alive_reader,
        runtime_reaper),
      error = identity)
    if (inherits(cleanup_result, "error"))
      message("V021 batch cleanup failure: ", conditionMessage(cleanup_result))

    release_result <- tryCatch({
      owned_ids <- manifest$direction_id[indices]
      active <- active_reservations$state == "reserved" &
        active_reservations$task_identity %in% owned_ids
      active_reservations$state[active] <<- "worker_terminated"
      write_v021_reservations_atomic(
        active_reservations, reservation_path, policy)
      TRUE
    }, error = identity)
    if (inherits(release_result, "error"))
      message(
        "V021 batch reservation cleanup failure: ",
        conditionMessage(release_result))
  }, add = TRUE)

  stop_owned <- function(snapshot, registry_snapshot = NULL) {
    registry_now <- registry_snapshot
    if (is.null(registry_now))
      registry_now <- tryCatch(
        readRDS(registry_file), error = function(e) registry)
    fit_rows <- registry_now[
      registry_now$worker_role == "direction_fit_worker", , drop = FALSE]
    if (!nrow(fit_rows)) return(invisible(TRUE))

    matched <- match(snapshot$pid, fit_rows$pid)
    same_identity <- !is.na(matched) &
      as.character(snapshot$start_time) ==
        as.character(fit_rows$start_time[matched])
    owned <- as.integer(snapshot$pid[same_identity])

    repeat {
      children <- snapshot$pid[snapshot$ppid %in% owned]
      expanded <- unique(c(owned, children))
      if (identical(sort(expanded), sort(owned))) break
      owned <- expanded
    }
    for (pid in rev(owned))
      system2(
        "kill", c("-INT", as.character(pid)),
        stdout = FALSE, stderr = FALSE)
    invisible(TRUE)
  }

  monitor_trace <- new_v021_failure_trace_context(
    trace_root, "monitor_child", batch_id,
    manifest$task_id[indices], manifest$direction_id[indices])
  monitor_job <- parallel::mcparallel(run_v021_traced_child(function() {
    while (!file.exists(monitor_start_file)) Sys.sleep(0.01)
    snapshots <- list()
    last_registry <- readRDS(registry_file)

    repeat {
      set_v021_failure_trace_phase(monitor_trace, "monitor_polling")
      registry_read <- tryCatch(readRDS(registry_file), error = identity)
      if (!inherits(registry_read, "error"))
        last_registry <- registry_read
      snap <- tryCatch(
        capture_v021_preflight_process_snapshot(parent_pid, executable),
        error = identity)

      if (inherits(registry_read, "error") || inherits(snap, "error")) {
        reason <- if (inherits(registry_read, "error"))
          conditionMessage(registry_read) else conditionMessage(snap)
        payload <- list(error = reason, snapshots = snapshots)
        .v021_atomic_save_rds(payload, latest_file)
        .v021_atomic_save_rds(payload, final_file)
        if (!inherits(snap, "error"))
          stop_owned(snap, last_registry)
        else if (length(snapshots))
          stop_owned(tail(snapshots, 1L)[[1L]], last_registry)
        break
      }

      registry_now <- registry_read
      snapshots[[length(snapshots) + 1L]] <- snap
      monitor <- monitor_v021_process_snapshots(
        list(snap), policy, parent_pid, executable,
        worker_registry = registry_now)
      payload <- list(
        snapshot = snap, latest_monitor = monitor, batch_id = batch_id)
      .v021_atomic_save_rds(payload, latest_file)

      if (!identical(monitor$monitoring_state, "verified") ||
          !identical(monitor$compliance_status, "compliant")) {
        payload$snapshots <- snapshots
        payload$error <- monitor$reason %||% "resource_ceiling_exceeded"
        .v021_atomic_save_rds(payload, latest_file)
        .v021_atomic_save_rds(payload, final_file)
        stop_owned(snap, registry_now)
        break
      }

      set_v021_failure_trace_phase(
        monitor_trace, "monitor_polling", completed = TRUE)
      if (file.exists(stop_file)) {
        set_v021_failure_trace_phase(monitor_trace, "monitor_shutdown")
        payload$snapshots <- snapshots
        .v021_atomic_save_rds(payload, final_file)
        set_v021_failure_trace_phase(
          monitor_trace, "monitor_shutdown", completed = TRUE)
        break
      }
      Sys.sleep(1)
    }
    TRUE
  }, monitor_trace), silent = TRUE)
  monitor_job_tracker$jobs <- setNames(
    list(monitor_job), as.character(monitor_job$pid))
  monitor_registry <- register_v021_preflight_worker(
    as.integer(monitor_job$pid), parent_pid, "resource_monitor", batch_id,
    paste0(batch_id, "-monitor"))
  registry <- rbind(registry, monitor_registry)
  persist_registry()

  launch_registered_fit <- function(j) {
    reserve_and_mark(j)
    spawned <- spawn_fit_job(j)
    registry <<- rbind(registry, spawned$registry)
    persist_registry()
    if (!file.create(spawned$gate_file))
      stop("Failed to release a registered rolling fit worker.")
    invisible(TRUE)
  }

  record_batch_condition <- function(condition) {
    batch_exit_condition <<- condition
    if (file.exists(latest_file))
      parent_trace$monitor_state <- tryCatch(
        readRDS(latest_file),
        error = function(e) list(read_error = conditionMessage(e)))
  }

  tryCatch({
    if (!file.create(monitor_start_file))
      stop("Failed to start the rolling resource monitor.")
    for (spawned in initial_spawned)
      if (!file.create(spawned$gate_file))
        stop("Failed to release an initial rolling fit worker.")
    set_v021_failure_trace_phase(parent_trace, "batch_launch", completed = TRUE)

    set_v021_failure_trace_phase(parent_trace, "fit_return_handling")
    next_position <- initial_count + 1L

    repeat {
      while (length(fit_job_tracker$jobs) < maximum_active &&
             next_position <= length(indices)) {
        launch_registered_fit(next_position)
        next_position <- next_position + 1L
      }

      if (!length(fit_job_tracker$jobs) &&
          next_position > length(indices))
        break

      collected <- collect_ready_jobs(fit_job_tracker)
      if (!length(collected)) {
        Sys.sleep(0.05)
        next
      }

      collected_pids <- names(collected)
      for (k in seq_along(collected)) {
        pid_key <- collected_pids[[k]]
        meta <- active_meta[[pid_key]]
        if (is.null(meta))
          stop("Collected rolling child was not present in active metadata.")

        j <- as.integer(meta$position)
        i <- indices[[j]]
        results[j] <- list(collected[[k]])
        elapsed_by_position[[j]] <-
          proc.time()[["elapsed"]] - as.numeric(meta$start_elapsed)
        rm(list = pid_key, envir = active_meta)

        attempt_number <- v021_task_attempt_number(
          execution_status, tasks$direction_index[[i]])
        active_reservations <<- release_v021_capacity(
          active_reservations, policy, manifest$direction_id[[i]],
          attempt_number, "main_fit", worker_terminated = TRUE)
        write_v021_reservations_atomic(
          active_reservations, reservation_path, policy)
      }
    }

    set_v021_failure_trace_phase(
      parent_trace, "fit_return_handling", completed = TRUE)

    # Predictive evaluation is still controller-side and sequential. Stop the
    # monitor first so K-fold forks cannot inherit the monitor pipe.
    if (!file.create(stop_file))
      stop("Failed to request rolling monitor shutdown.")
    set_v021_failure_trace_phase(parent_trace, "monitor_shutdown")
    monitor_collected <- unname(
      collect_v021_tracked_jobs(monitor_job_tracker))
    if (length(monitor_collected) != 1L)
      stop("Monitor child did not return exactly one collected result.")
    monitor_result <- monitor_collected[[1L]]
    if (inherits(monitor_result, "v021_traced_child_error"))
      stop(monitor_result$condition_message)

    wait_v021_collected_child_exit(
      as.integer(monitor_job$pid),
      expected_start_time = as.character(monitor_registry$start_time[[1L]]))
    set_v021_failure_trace_phase(
      parent_trace, "monitor_shutdown", completed = TRUE)

    payload <- readRDS(final_file)
    parent_trace$monitor_state <-
      payload$latest_monitor %||% payload$error %||% NULL
    if (!is.null(payload$error))
      stop("Live process monitor failed: ", payload$error)

    batch_monitor <- monitor_v021_process_snapshots(
      payload$snapshots, policy, parent_pid, executable,
      worker_registry = registry)
    validate_v021_preflight_monitor(batch_monitor, policy)
    peak_chains <<- max(
      peak_chains, batch_monitor$observed_peak_active_cmdstan_chains)
    peak_slots <<- max(
      peak_slots, batch_monitor$observed_peak_active_cmdstan_processes)

    run_predictive_rolling <- function() {
      valid_positions <- which(!vapply(
        results,
        inherits,
        logical(1),
        what = c("try-error", "pclv_failure", "v021_traced_child_error")
      ))
      if (!length(valid_positions))
        return(list(results = results, elapsed = rep(0, length(indices))))

      direction_plans <- lapply(valid_positions, function(j) {
        i <- indices[[j]]
        prepare_v021_predictive_fold_plan(
          results[[j]], j, manifest$direction_id[[i]])
      })
      queued <- build_v021_predictive_fold_queue(direction_plans)
      direction_plans <- queued$plans
      fold_tasks <- queued$tasks
      predictive_results <- results
      predictive_elapsed <- rep(0, length(indices))
      for (plan in direction_plans)
        predictive_results[[plan$position]] <- plan$result

      if (!length(fold_tasks))
        return(list(results = predictive_results,
                    elapsed = predictive_elapsed))

      maximum_predictive <- min(
        length(fold_tasks),
        as.integer(config$maximum_simultaneous_kfold_fits)
      )
      derive_safe_outer_concurrency(
        policy,
        build_v021_operation_spec("kfold_fit", maximum_predictive)
      )

      predictive_batch_id <- paste0(batch_id, "-kfold")
      predictive_stop_file <- tempfile("v021-full-kfold-monitor-stop-")
      predictive_monitor_start_file <- tempfile(
        "v021-full-kfold-monitor-start-")
      predictive_registry_file <- tempfile(
        "v021-full-kfold-registry-", fileext = ".rds")
      predictive_final_file <- file.path(
        paths$monitor_root, paste0(predictive_batch_id, ".rds"))
      predictive_latest_file <- file.path(paths$monitor_root, "latest.rds")
      fold_results <- vector("list", length(fold_tasks))
      predictive_meta <- new.env(hash = TRUE, parent = emptyenv())
      predictive_registry <- NULL
      predictive_ownership <- NULL
      predictive_exit_condition <- NULL
      predictive_cleanup_installed <- FALSE

      # Fail closed even if setup fails before durable predictive ownership is
      # installed. Every fold child is initially held behind a private gate.
      on.exit({
        if (!predictive_cleanup_installed) {
          jobs <- predictive_job_tracker$jobs
          if (length(jobs)) {
            pids <- vapply(
              jobs, function(job) as.integer(job$pid), integer(1))
            try(runtime_signaler(pids, "SIGINT"), silent = TRUE)
            Sys.sleep(0.1)
            try(runtime_signaler(pids, "SIGTERM"), silent = TRUE)
          }
          try(reap_v021_uncollected_jobs(
            predictive_job_tracker, predictive_monitor_job_tracker),
            silent = TRUE)
          try({
            fold_ids <- vapply(
              fold_tasks, `[[`, character(1), "fold_identity")
            active <- active_reservations$state == "reserved" &
              active_reservations$job_type == "kfold_fit" &
              active_reservations$task_identity %in% fold_ids
            active_reservations$state[active] <<- "worker_terminated"
            write_v021_reservations_atomic(
              active_reservations, reservation_path, policy)
          }, silent = TRUE)
        }
      }, add = TRUE)

      predictive_traces <- lapply(seq_along(fold_tasks), function(q) {
        task <- fold_tasks[[q]]
        i <- indices[[task$position]]
        new_v021_failure_trace_context(
          trace_root, "predictive_fold_worker", predictive_batch_id,
          manifest$task_id[[i]], task$fold_identity)
      })

      reserve_predictive <- function(q) {
        task <- fold_tasks[[q]]
        i <- indices[[task$position]]
        attempt_number <- v021_task_attempt_number(
          execution_status, tasks$direction_index[[i]])
        active_reservations <<- reserve_v021_capacity(
          active_reservations, policy, task$fold_identity,
          attempt_number, "kfold_fit")
        write_v021_reservations_atomic(
          active_reservations, reservation_path, policy)
      }

      release_predictive <- function(q) {
        task <- fold_tasks[[q]]
        i <- indices[[task$position]]
        attempt_number <- v021_task_attempt_number(
          execution_status, tasks$direction_index[[i]])
        active_reservations <<- release_v021_capacity(
          active_reservations, policy, task$fold_identity,
          attempt_number, "kfold_fit", worker_terminated = TRUE)
        write_v021_reservations_atomic(
          active_reservations, reservation_path, policy)
      }

      spawn_predictive_job <- function(q) {
        task <- fold_tasks[[q]]
        gate_file <- tempfile("v021-full-kfold-fold-start-")
        trace <- predictive_traces[[q]]
        job <- parallel::mcparallel({
          while (!file.exists(gate_file)) Sys.sleep(0.01)
          set_v021_failure_trace_phase(trace, "predictive_fold_execution")
          run_v021_traced_child(function() {
            value <- withr::with_options(
              list(pclvbayes.parallel_chains_override =
                     as.integer(config$kfold_parallel_chains)),
              with_v021_single_thread_environment(function() {
                run_v021_predictive_fold_task(task)
              })
            )
            set_v021_failure_trace_phase(
              trace, "predictive_fold_execution", completed = TRUE)
            value
          }, trace)
        }, silent = TRUE)

        pid <- as.integer(job$pid)
        pid_key <- as.character(pid)
        predictive_job_tracker$jobs[[pid_key]] <- job
        predictive_meta[[pid_key]] <- list(
          task_index = as.integer(q),
          direction_position = as.integer(task$position),
          start_elapsed = proc.time()[["elapsed"]])
        list(
          gate_file = gate_file,
          registry = register_v021_preflight_worker(
            pid, parent_pid, "predictive_fold_worker", predictive_batch_id,
            task$fold_identity)
        )
      }

      initial_count <- min(length(fold_tasks), maximum_predictive)
      initial_indices <- seq_len(initial_count)
      initial_spawned <- vector("list", initial_count)
      for (k in seq_along(initial_indices)) {
        q <- initial_indices[[k]]
        reserve_predictive(q)
        initial_spawned[[k]] <- spawn_predictive_job(q)
      }
      predictive_fold_registry <- do.call(
        rbind, lapply(initial_spawned, function(value) value$registry))
      predictive_registry <- predictive_fold_registry

      predictive_cleanup_path <- file.path(
        paths$output_root, "cleanup_audits",
        paste0(predictive_batch_id, ".rds"))
      predictive_ownership <- new_v021_batch_ownership(
        predictive_batch_id, parent_pid, predictive_fold_registry,
        executable, paths$output_root, predictive_cleanup_path)

      persist_predictive_registry <- function() {
        validate_v021_worker_registry(predictive_registry)
        updated <- predictive_ownership
        updated$worker_registry <- predictive_registry
        predictive_ownership <<- updated
        .v021_atomic_save_rds(
          predictive_registry, predictive_registry_file)
        .v021_atomic_save_rds(
          predictive_registry,
          file.path(
            paths$output_root, "worker_registries",
            paste0(predictive_batch_id, ".rds")))
        .v021_atomic_save_rds(list(
          ownership_schema = "v021_full_batch_ownership_v1",
          batch_id = predictive_batch_id,
          parent_pid = as.integer(parent_pid),
          known_model_executable = executable,
          output_root = paths$output_root,
          worker_registry = predictive_registry,
          registration_timestamp =
            format(Sys.time(), tz = "UTC", usetz = TRUE)),
          file.path(
            paths$output_root, "batch_ownership",
            paste0(predictive_batch_id, ".rds")))
        invisible(TRUE)
      }

      predictive_inventory_reader <- function() {
        snapshot <- capture_v021_preflight_process_snapshot(
          parent_pid, executable)
        classify_v021_process_snapshot(
          snapshot, parent_pid, executable,
          worker_registry = predictive_ownership$worker_registry)
      }
      predictive_reaper <- function()
        reap_v021_uncollected_jobs(
          predictive_job_tracker, predictive_monitor_job_tracker)

      on.exit({
        cleanup_result <- tryCatch(
          cleanup_v021_owned_batch(
            predictive_ownership, predictive_exit_condition,
            "parent", "predictive_return_handling",
            predictive_inventory_reader, runtime_signaler,
            runtime_alive_reader, predictive_reaper),
          error = identity)
        if (inherits(cleanup_result, "error"))
          message(
            "V021 predictive cleanup failure: ",
            conditionMessage(cleanup_result))

        release_result <- tryCatch({
          fold_ids <- vapply(
            fold_tasks, `[[`, character(1), "fold_identity")
          active <- active_reservations$state == "reserved" &
            active_reservations$job_type == "kfold_fit" &
            active_reservations$task_identity %in% fold_ids
          active_reservations$state[active] <<- "worker_terminated"
          write_v021_reservations_atomic(
            active_reservations, reservation_path, policy)
          TRUE
        }, error = identity)
        if (inherits(release_result, "error"))
          message(
            "V021 predictive reservation cleanup failure: ",
            conditionMessage(release_result))
      }, add = TRUE)
      predictive_cleanup_installed <- TRUE

      stop_owned_predictive <- function(snapshot, registry_snapshot = NULL) {
        registry_now <- registry_snapshot
        if (is.null(registry_now))
          registry_now <- tryCatch(
            readRDS(predictive_registry_file),
            error = function(e) predictive_registry)
        fit_rows <- registry_now[
          registry_now$worker_role == "predictive_fold_worker",
          , drop = FALSE]
        if (!nrow(fit_rows)) return(invisible(TRUE))

        matched <- match(snapshot$pid, fit_rows$pid)
        same_identity <- !is.na(matched) &
          as.character(snapshot$start_time) ==
            as.character(fit_rows$start_time[matched])
        owned <- as.integer(snapshot$pid[same_identity])
        repeat {
          children <- snapshot$pid[snapshot$ppid %in% owned]
          expanded <- unique(c(owned, children))
          if (identical(sort(expanded), sort(owned))) break
          owned <- expanded
        }
        for (pid in rev(owned))
          system2(
            "kill", c("-INT", as.character(pid)),
            stdout = FALSE, stderr = FALSE)
        invisible(TRUE)
      }

      predictive_monitor_trace <- new_v021_failure_trace_context(
        trace_root, "predictive_monitor_child", predictive_batch_id,
        manifest$task_id[indices[valid_positions]],
        manifest$direction_id[indices[valid_positions]])
      predictive_monitor_job <- parallel::mcparallel(
        run_v021_traced_child(function() {
          while (!file.exists(predictive_monitor_start_file))
            Sys.sleep(0.01)
          snapshots <- list()
          last_registry <- readRDS(predictive_registry_file)

          repeat {
            set_v021_failure_trace_phase(
              predictive_monitor_trace, "predictive_monitor_polling")
            registry_read <- tryCatch(
              readRDS(predictive_registry_file), error = identity)
            if (!inherits(registry_read, "error"))
              last_registry <- registry_read
            snap <- tryCatch(
              capture_v021_preflight_process_snapshot(
                parent_pid, executable),
              error = identity)

            if (inherits(registry_read, "error") ||
                inherits(snap, "error")) {
              reason <- if (inherits(registry_read, "error"))
                conditionMessage(registry_read) else conditionMessage(snap)
              payload <- list(error = reason, snapshots = snapshots)
              .v021_atomic_save_rds(payload, predictive_latest_file)
              .v021_atomic_save_rds(payload, predictive_final_file)
              if (!inherits(snap, "error"))
                stop_owned_predictive(snap, last_registry)
              else if (length(snapshots))
                stop_owned_predictive(
                  tail(snapshots, 1L)[[1L]], last_registry)
              break
            }

            registry_now <- registry_read
            snapshots[[length(snapshots) + 1L]] <- snap
            monitor <- monitor_v021_process_snapshots(
              list(snap), policy, parent_pid, executable,
              worker_registry = registry_now)
            payload <- list(
              snapshot = snap, latest_monitor = monitor,
              batch_id = predictive_batch_id)
            .v021_atomic_save_rds(payload, predictive_latest_file)

            if (!identical(monitor$monitoring_state, "verified") ||
                !identical(monitor$compliance_status, "compliant")) {
              payload$snapshots <- snapshots
              payload$error <- monitor$reason %||%
                "resource_ceiling_exceeded"
              .v021_atomic_save_rds(payload, predictive_latest_file)
              .v021_atomic_save_rds(payload, predictive_final_file)
              stop_owned_predictive(snap, registry_now)
              break
            }

            set_v021_failure_trace_phase(
              predictive_monitor_trace, "predictive_monitor_polling",
              completed = TRUE)
            if (file.exists(predictive_stop_file)) {
              set_v021_failure_trace_phase(
                predictive_monitor_trace, "predictive_monitor_shutdown")
              payload$snapshots <- snapshots
              .v021_atomic_save_rds(payload, predictive_final_file)
              set_v021_failure_trace_phase(
                predictive_monitor_trace, "predictive_monitor_shutdown",
                completed = TRUE)
              break
            }
            Sys.sleep(1)
          }
          TRUE
        }, predictive_monitor_trace),
        silent = TRUE)
      predictive_monitor_job_tracker$jobs <- setNames(
        list(predictive_monitor_job),
        as.character(predictive_monitor_job$pid))
      predictive_monitor_registry <- register_v021_preflight_worker(
        as.integer(predictive_monitor_job$pid), parent_pid,
        "predictive_resource_monitor", predictive_batch_id,
        paste0(predictive_batch_id, "-monitor"))
      predictive_registry <- rbind(
        predictive_registry, predictive_monitor_registry)
      persist_predictive_registry()

      launch_registered_predictive <- function(q) {
        reserve_predictive(q)
        spawned <- spawn_predictive_job(q)
        predictive_registry <<- rbind(
          predictive_registry, spawned$registry)
        persist_predictive_registry()
        if (!file.create(spawned$gate_file))
          stop("Failed to release a registered predictive fold worker.")
        invisible(TRUE)
      }

      if (!file.create(predictive_monitor_start_file))
        stop("Failed to start the predictive resource monitor.")
      for (spawned in initial_spawned)
        if (!file.create(spawned$gate_file))
          stop("Failed to release an initial predictive fold worker.")

      next_task <- initial_count + 1L
      repeat {
        while (length(predictive_job_tracker$jobs) < maximum_predictive &&
               next_task <= length(fold_tasks)) {
          launch_registered_predictive(next_task)
          next_task <- next_task + 1L
        }

        if (!length(predictive_job_tracker$jobs) &&
            next_task > length(fold_tasks))
          break

        collected <- collect_ready_jobs(predictive_job_tracker)
        if (!length(collected)) {
          Sys.sleep(0.05)
          next
        }

        for (pid_key in names(collected)) {
          meta <- predictive_meta[[pid_key]]
          if (is.null(meta))
            stop(
              "Collected predictive fold child was absent from active metadata.")
          q <- as.integer(meta$task_index)
          fold_results[[q]] <- collected[[pid_key]]
          elapsed <- proc.time()[["elapsed"]] -
            as.numeric(meta$start_elapsed)
          predictive_elapsed[[as.integer(meta$direction_position)]] <-
            predictive_elapsed[[as.integer(meta$direction_position)]] + elapsed
          rm(list = pid_key, envir = predictive_meta)
          release_predictive(q)
        }
      }

      if (!file.create(predictive_stop_file))
        stop("Failed to request predictive monitor shutdown.")
      predictive_monitor_collected <- unname(
        collect_v021_tracked_jobs(predictive_monitor_job_tracker))
      if (length(predictive_monitor_collected) != 1L)
        stop(
          "Predictive monitor child did not return exactly one result.")
      predictive_monitor_result <- predictive_monitor_collected[[1L]]
      if (inherits(
          predictive_monitor_result, "v021_traced_child_error"))
        stop(predictive_monitor_result$condition_message)
      wait_v021_collected_child_exit(
        as.integer(predictive_monitor_job$pid),
        expected_start_time = as.character(
          predictive_monitor_registry$start_time[[1L]]))

      predictive_payload <- readRDS(predictive_final_file)
      parent_trace$monitor_state <-
        predictive_payload$latest_monitor %||%
          predictive_payload$error %||% NULL
      if (!is.null(predictive_payload$error))
        stop(
          "Live predictive process monitor failed: ",
          predictive_payload$error)
      predictive_monitor <- monitor_v021_process_snapshots(
        predictive_payload$snapshots, policy, parent_pid, executable,
        worker_registry = predictive_registry)
      validate_v021_preflight_monitor(predictive_monitor, policy)
      peak_chains <<- max(
        peak_chains,
        predictive_monitor$observed_peak_active_cmdstan_chains)
      peak_slots <<- max(
        peak_slots,
        predictive_monitor$observed_peak_active_cmdstan_processes)

      for (plan in direction_plans) {
        if (isTRUE(plan$eligible)) {
          predictive_results[[plan$position]] <-
            finalize_v021_predictive_direction(
              plan, fold_results[plan$queue_indices])
        }
      }
      list(results = predictive_results, elapsed = predictive_elapsed)
    }

    # The main monitor is fully collected before predictive workers are forked.
    # K-fold fold fits then run as one global 12x1 rolling phase across
    # direction, repetition, and fold. Durable
    # artifacts remain finalized in canonical direction order.
    predictive <- run_predictive_rolling()
    results <- predictive$results
    predictive_elapsed <- predictive$elapsed

    for (j in seq_along(indices)) {
      i <- indices[[j]]
      result <- results[[j]]
      elapsed <- elapsed_by_position[[j]]
      if (!is.finite(elapsed)) elapsed <- 0
      if (is.finite(predictive_elapsed[[j]]))
        elapsed <- elapsed + predictive_elapsed[[j]]

      store_outcome(i, result, elapsed)
      results[j] <- list(NULL)
      update_status()
    }
  }, error = function(condition)
       record_and_resignal_v021_condition(condition, record_batch_condition),
     interrupt = function(condition)
       record_and_resignal_v021_condition(condition, record_batch_condition))

  invisible(TRUE)
}

old_options <- options(glvpair.output_root = file.path(paths$output_root, "cmdstan-owned"))
on.exit(options(old_options), add = TRUE)
update_status("running")

execution_error <- NULL
parent_condition_handler <- function(condition) {
  if (is.null(parent_trace$last_trace)) capture_v021_failure_trace(condition, parent_trace)
}
record_execution_condition <- function(condition) execution_error <<- condition
execution_error <- tryCatch(
  withCallingHandlers(
    tryCatch(with_v021_single_thread_environment(function() {
      repeat {
        restart <- validate_v021_full_restart_states(manifest)
        manifest <<- restart$manifest
        write_v021_manifest_atomic(manifest, paths$manifest)
        execution_plan <- plan_v021_execution_resume(
          prepared_execution$manifest, paths$output_root,
          config$maximum_task_attempts, persist_reconciliation = TRUE)
        execution_status <<- execution_plan$status
        attempt_ledger <<- execution_plan$attempt_ledger
        if (!length(execution_plan$runnable_indices)) break
        runnable_ordinals <- prepared_execution$manifest$tasks$task_ordinal[
          execution_plan$runnable_indices]
        if (debug_mode)
          runnable_ordinals <- runnable_ordinals[
            runnable_ordinals %in% tasks$direction_index]
        if (!length(runnable_ordinals)) break
        runnable_indices <- match(runnable_ordinals, manifest$direction_index)
        if (anyNA(runnable_indices)) stop("V021-03 runnable task identity mismatch.")
        run_batch(head(runnable_indices, rolling_window_size))
      }
    }), error = function(condition)
         record_and_resignal_v021_condition(condition, record_execution_condition),
       interrupt = function(condition)
         record_and_resignal_v021_condition(condition, record_execution_condition)),
    error = parent_condition_handler, interrupt = parent_condition_handler),
  error = identity, interrupt = identity)
if (!inherits(execution_error, "condition")) execution_error <- NULL

if (!is.null(execution_error)) {
  set_v021_failure_trace_phase(parent_trace, "summary_failure_persistence")
  update_status("failed", conditionMessage(execution_error))
  trace_paths <- list.files(trace_root, pattern = "^v021_full_failure_trace_v1-.*\\.rds$",
                            full.names = TRUE)
  failure_persistence <- persist_v021_failure_payload(list(
    failure_schema = "v021_full_100_species_failure_v1",
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    reason = conditionMessage(execution_error), manifest_path = paths$manifest,
    monitor_path = file.path(paths$monitor_root, "latest.rds"),
    cleanup_audit_path = parent_trace$cleanup_audit_path,
    parent_trace = parent_trace$last_trace,
    child_trace_paths = trace_paths),
    file.path(paths$output_root, "failure_payload.rds"), trace_root)
  if (!isTRUE(failure_persistence$persisted))
    message("V021 failure-payload persistence failure: ",
            failure_persistence$persistence_error)
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
set_v021_failure_trace_phase(parent_trace, "summary_failure_persistence")
.v021_atomic_save_rds(summary, paths$summary)
update_status("completed")
cat("V021-06 completed.\n")
}

run_v021_full_100_species()
