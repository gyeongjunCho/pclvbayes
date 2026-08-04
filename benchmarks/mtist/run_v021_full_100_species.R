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

policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
validate_v021_resource_policy(policy)
derivation <- derive_safe_outer_concurrency(
  policy, build_v021_operation_spec("main_fit", 3L))
if (!identical(policy$policy_schema, "v021_resource_policy_v3") ||
    !identical(derivation$projected_active_cmdstan_chains, 12L) ||
    !identical(derivation$projected_active_cmdstan_processes, 12L))
  stop("V021-06 requires the verified policy-v3 12-chain/12-slot contract.")
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

durable_peaks <- recover_v021_monitor_peaks(paths$output_root, policy, executable)
peak_chains <- durable_peaks$chains
peak_slots <- durable_peaks$processes
batch_number <- 0L
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
  set_v021_failure_trace_phase(parent_trace, "batch_launch")
  validate_v021_controller_ownership(
    controller_ownership, paths$output_root,
    prepared_execution$manifest$manifest_hash,
    prepared_execution$manifest$configuration_hash)
  validate_v021_preflight_launch_capacity(policy, length(indices))
  batch_number <<- batch_number + 1L
  batch_id <- sprintf("batch-%05d", batch_number)
  parent_trace$batch_id <- batch_id
  parent_trace$task_ids <- as.integer(manifest$task_id[indices])
  parent_trace$direction_ids <- as.character(manifest$direction_id[indices])
  for (i in indices) {
    attempt_number <- v021_task_attempt_number(
      execution_status, tasks$direction_index[[i]], next_attempt = TRUE)
    active_reservations <<- reserve_v021_capacity(
      active_reservations, policy, manifest$direction_id[[i]],
      attempt_number, "main_fit")
  }
  write_v021_reservations_atomic(active_reservations, reservation_path, policy)
  mark_running(indices)
  set_v021_failure_trace_phase(parent_trace, "batch_launch", completed = TRUE)
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
  worker_traces <- lapply(seq_along(indices), function(j)
    new_v021_failure_trace_context(
      trace_root, "outer_worker", batch_id,
      manifest$task_id[[indices[[j]]]], manifest$direction_id[[indices[[j]]]]))
  fit_jobs <- lapply(seq_along(indices), function(j) parallel::mcparallel({
    while (!file.exists(start_file)) Sys.sleep(0.01)
    set_v021_failure_trace_phase(worker_traces[[j]], "worker_execution")
    run_v021_traced_child(function() {
      value <- with_v021_single_thread_environment(jobs[[j]])
      set_v021_failure_trace_phase(worker_traces[[j]], "worker_execution",
                                   completed = TRUE)
      value
    }, worker_traces[[j]])
  }, silent = TRUE))
  fit_pids <- vapply(fit_jobs, function(job) as.integer(job$pid), integer(1))
  fit_registry <- do.call(rbind, lapply(seq_along(fit_pids), function(j)
    register_v021_preflight_worker(
      fit_pids[[j]], parent_pid, "direction_fit_worker", batch_id,
      manifest$direction_id[[indices[[j]]]])))
  cleanup_path <- file.path(paths$output_root, "cleanup_audits",
                            paste0(batch_id, ".rds"))
  ownership <- new_v021_batch_ownership(
    batch_id, parent_pid, fit_registry, executable, paths$output_root, cleanup_path)
  parent_trace$cleanup_audit_path <- cleanup_path
  batch_exit_condition <- NULL
  runtime_inventory_reader <- function() {
    snapshot <- capture_v021_preflight_process_snapshot(parent_pid, executable)
    classify_v021_process_snapshot(
      snapshot, parent_pid, executable, worker_registry = ownership$worker_registry)
  }
  runtime_alive_reader <- function(pids, start_times) {
    if (!length(pids)) return(integer())
    keep <- vapply(seq_along(pids), function(k) {
      stat <- tryCatch(readLines(file.path("/proc", pids[[k]], "stat"),
                                 warn = FALSE, n = 1L),
                       error = function(e) character())
      if (!length(stat)) return(FALSE)
      parsed <- tryCatch(.v021_parse_linux_stat(stat, pids[[k]]), error = identity)
      !inherits(parsed, "error") && identical(parsed$start_time, start_times[[k]]) &&
        !identical(parsed$process_state, "Z")
    }, logical(1))
    as.integer(pids[keep])
  }
  runtime_signaler <- function(pids, signal) {
    option <- switch(signal, SIGINT = "-INT", SIGTERM = "-TERM",
                     stop("Unsupported cleanup signal."))
    invisible(lapply(as.integer(pids), function(pid)
      system2("kill", c(option, as.character(pid)), stdout = FALSE, stderr = FALSE)))
  }
  fit_job_tracker <- new_v021_child_job_tracker(fit_jobs)
  monitor_job_tracker <- new_v021_child_job_tracker()
  runtime_reaper <- function()
    reap_v021_uncollected_jobs(fit_job_tracker, monitor_job_tracker)
  on.exit({
    cleanup_result <- tryCatch(cleanup_v021_owned_batch(
      ownership, batch_exit_condition, "parent", parent_trace$phase,
      runtime_inventory_reader, runtime_signaler, runtime_alive_reader,
      runtime_reaper), error = identity)
    if (inherits(cleanup_result, "error"))
      message("V021 batch cleanup failure: ", conditionMessage(cleanup_result))
    release_result <- tryCatch({
      owned_ids <- manifest$direction_id[indices]
      active <- active_reservations$state == "reserved" &
        active_reservations$task_identity %in% owned_ids
      active_reservations$state[active] <<- "worker_terminated"
      write_v021_reservations_atomic(active_reservations, reservation_path, policy)
      TRUE
    }, error = identity)
    if (inherits(release_result, "error"))
      message("V021 batch reservation cleanup failure: ", conditionMessage(release_result))
  }, add = TRUE)
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
  monitor_trace <- new_v021_failure_trace_context(
    trace_root, "monitor_child", batch_id,
    manifest$task_id[indices], manifest$direction_id[indices])
  monitor_job <- parallel::mcparallel(run_v021_traced_child(function() {
    while (!file.exists(monitor_start_file)) Sys.sleep(0.01)
    registry <- readRDS(registry_file)
    snapshots <- list()
    repeat {
      set_v021_failure_trace_phase(monitor_trace, "monitor_polling")
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
      set_v021_failure_trace_phase(monitor_trace, "monitor_polling", completed = TRUE)
      if (file.exists(stop_file)) {
        set_v021_failure_trace_phase(monitor_trace, "monitor_shutdown")
        payload$snapshots <- snapshots
        .v021_atomic_save_rds(payload, final_file)
        set_v021_failure_trace_phase(monitor_trace, "monitor_shutdown",
                                     completed = TRUE)
        break
      }
      Sys.sleep(1)
    }
    TRUE
  }, monitor_trace), silent = TRUE)
  monitor_job_tracker$jobs <- list(monitor_job)
  monitor_registry <- register_v021_preflight_worker(
    as.integer(monitor_job$pid), parent_pid, "resource_monitor", batch_id,
    paste0(batch_id, "-monitor"))
  registry <- rbind(fit_registry, monitor_registry)
  validate_v021_worker_registry(registry)
  ownership$worker_registry <- registry
  .v021_atomic_save_rds(registry, registry_file)
  .v021_atomic_save_rds(registry, file.path(paths$output_root, "worker_registries",
                                            paste0(batch_id, ".rds")))
  .v021_atomic_save_rds(list(
    ownership_schema = "v021_full_batch_ownership_v1",
    batch_id = batch_id,
    parent_pid = as.integer(parent_pid),
    known_model_executable = executable,
    output_root = paths$output_root,
    worker_registry = registry,
    registration_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    file.path(paths$output_root, "batch_ownership", paste0(batch_id, ".rds")))
  record_batch_condition <- function(condition) {
    batch_exit_condition <<- condition
    if (file.exists(latest_file))
      parent_trace$monitor_state <- tryCatch(readRDS(latest_file), error = function(e)
        list(read_error = conditionMessage(e)))
  }
  tryCatch({
    file.create(monitor_start_file)
    file.create(start_file)
  set_v021_failure_trace_phase(parent_trace, "fit_return_handling")
  results <- unname(collect_v021_tracked_jobs(fit_job_tracker))
  for (i in indices) {
    attempt_number <- v021_task_attempt_number(
      execution_status, tasks$direction_index[[i]])
    active_reservations <<- release_v021_capacity(
      active_reservations, policy, manifest$direction_id[[i]],
      attempt_number, "main_fit", worker_terminated = TRUE)
  }
  write_v021_reservations_atomic(active_reservations, reservation_path, policy)
  set_v021_failure_trace_phase(parent_trace, "fit_return_handling", completed = TRUE)
  # Close the concurrent-main-fit monitor before spawning any controller-side
  # K-fold process. Forked K-fold utilities must not inherit the monitor pipe.
  file.create(stop_file)
  set_v021_failure_trace_phase(parent_trace, "monitor_shutdown")
  monitor_collected <- unname(
    collect_v021_tracked_jobs(monitor_job_tracker)
  )
  if (length(monitor_collected) != 1L)
    stop("Monitor child did not return exactly one collected result.")
  monitor_result <- monitor_collected[[1L]]
  if (inherits(monitor_result, "v021_traced_child_error"))
    stop(monitor_result$condition_message)

  # The scoped collection above is the single authoritative reap operation.
  # Do not mutate parallel's private child registry or collect unrelated jobs.
  wait_v021_collected_child_exit(
    as.integer(monitor_job$pid),
    expected_start_time = as.character(monitor_registry$start_time[[1L]])
  )
  set_v021_failure_trace_phase(parent_trace, "monitor_shutdown", completed = TRUE)
  payload <- readRDS(final_file)
  parent_trace$monitor_state <- payload$latest_monitor %||% payload$error %||% NULL
  if (!is.null(payload$error)) stop("Live process monitor failed: ", payload$error)
  batch_monitor <- monitor_v021_process_snapshots(
    payload$snapshots, policy, parent_pid, executable, worker_registry = registry)
  validate_v021_preflight_monitor(batch_monitor, policy)
  peak_chains <<- max(peak_chains, batch_monitor$observed_peak_active_cmdstan_chains)
  peak_slots <<- max(peak_slots, batch_monitor$observed_peak_active_cmdstan_processes)
  # Main fits finish as one canonical batch. Predictive evaluation then follows
  # direction order, one eligible direction at a time, under the shared budget.
  evaluate_predictive <- function(i, result) {
    direction_id <- manifest$direction_id[[i]]
    attempt_number <- v021_task_attempt_number(
      execution_status, tasks$direction_index[[i]])
    evaluated <- with_v021_predictive_reservation(
      result, active_reservations, policy, direction_id, attempt_number,
      evaluator = pclvbayes:::.add_predictive_evaluation,
      persist = function(value) {
        active_reservations <<- value
        write_v021_reservations_atomic(value, reservation_path, policy)
      })
    active_reservations <<- evaluated$reservations
    evaluated$result
  }
  for (j in seq_along(indices)) {
    if (inherits(results[[j]], c("try-error", "pclv_failure", "v021_traced_child_error")))
      next
    i <- indices[[j]]
    results[[j]] <- evaluate_predictive(i, results[[j]])
  }
  elapsed <- (proc.time()[["elapsed"]] - started) / length(indices)
  for (j in seq_along(indices)) store_outcome(indices[[j]], results[[j]], elapsed)
  update_status()
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
        run_batch(head(runnable_indices, config$maximum_simultaneous_fits))
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

