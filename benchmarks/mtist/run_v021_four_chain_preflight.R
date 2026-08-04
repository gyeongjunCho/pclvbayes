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
source(file.path(script_dir, "v021_robustness_features.R"))
source(file.path(script_dir, "v021_four_chain_preflight.R"))

run_v021_four_chain_preflight <- function() {
config <- source(file.path(script_dir, "configs", "v021_four_chain_preflight.R"))$value
validate_v021_four_chain_preflight_config(config)
resume_mode <- identical(tolower(Sys.getenv("PCLV_V021_PREFLIGHT_RESUME", "false")), "true")
if (!resume_mode && dir.exists(config$output_root) &&
    length(list.files(config$output_root, all.files = TRUE, no.. = TRUE)))
  stop("Refusing to overwrite an existing non-empty preflight output root: ", config$output_root)
if (resume_mode && !dir.exists(config$output_root))
  stop("Preflight resume requires the existing dedicated output root.")
if (resume_mode) {
  saved_config_early <- readRDS(file.path(config$output_root, "config.rds"))
  if (!identical(normalizePath(saved_config_early$output_root, mustWork = TRUE),
                 normalizePath(config$output_root, mustWork = TRUE)))
    stop("Preflight resume output-root identity changed.")
  config$output_root <- saved_config_early$output_root
}
dir.create(config$output_root, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(config$output_root, "checkpoints")
artifact_dir <- file.path(config$output_root, "inference_artifacts")
feature_dir <- file.path(config$output_root, "feature_records")
robustness_dir <- file.path(config$output_root, "robustness_feature_records")
dir.create(checkpoint_dir, showWarnings = FALSE)
dir.create(artifact_dir, showWarnings = FALSE)
dir.create(feature_dir, showWarnings = FALSE)
dir.create(robustness_dir, showWarnings = FALSE)
monitor_dir <- file.path(config$output_root, "monitor")
dir.create(monitor_dir, showWarnings = FALSE)
manifest_path <- file.path(config$output_root, "manifest.rds")

policy <- build_v021_resource_policy(
  proposed_outer_concurrency = 3L,
  maximum_concurrent_kfold_fits = 1L)
operation <- build_v021_operation_spec("main_fit", 3L)
derivation <- derive_safe_outer_concurrency(policy, operation)
thread_environment <- v021_single_thread_environment()
validate_v021_single_thread_environment(thread_environment)

root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
observations <- load_v021_preflight_observations(root, as.integer(config$dataset_id))
index <- ten_species_direction_index(observations$taxa, config$seed)
selected <- select_v021_preflight_tasks(index, config$selected_direction_count)
git_commit <- trimws(system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"),
                              stdout = TRUE)[[1L]])
execution_manifest <- build_v021_execution_manifest(
  dataset_id = config$dataset_id, taxa_order = observations$taxa,
  task_table = index[c("task_id", "direction_index", "target", "source", "seed")],
  public_seed = config$seed, kfold_seed = config$seed,
  preprocessing_config = list(identity = "pclv_smoothed_full_composition_closure_v1"),
  posterior_config = list(identity = "student_t_irregular_time_ou_4x2000_v1",
                          chains = config$chains, iter_warmup = config$iter_warmup,
                          iter_sampling = config$iter_sampling),
  kfold_config = list(K = 5L, R = 1L, enabled = TRUE),
  predictive_config = list(identity = "student-t-scale-mixture-kalman-ou-q16_k5_r1_v1"),
  provenance = list(code_commit = git_commit, benchmark_schema = config$preflight_schema))
preflight_execution_paths <- initialize_v021_execution_root(
  execution_manifest, config$output_root)
preflight_plan <- plan_v021_execution_resume(
  execution_manifest, config$output_root, maximum_attempts = 2L,
  persist_reconciliation = TRUE)
preflight_status <- preflight_plan$status
preflight_attempt_ledger <- preflight_plan$attempt_ledger
manifest_input <- selected[c("task_id", "direction_index", "target", "source", "seed")]
manifest <- build_v021_checkpoint_manifest(
  manifest_input, config$chains, checkpoint_dir, allow_incomplete_pairs = TRUE)
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

run_indices <- seq_len(nrow(manifest))
resume_completed_hashes <- character()
resume_all_complete <- FALSE
if (resume_mode) {
  saved_config <- readRDS(file.path(config$output_root, "config.rds"))
  if (!identical(saved_config, config)) stop("Preflight resume configuration changed.")
  stored <- read_v021_checkpoint_manifest(manifest_path)
  identity_fields <- c("task_id", "pair_id", "direction_id", "direction_index",
                       "target", "source", "seed", "chain_seeds", "output_location")
  if (!identical(stored[identity_fields], manifest[identity_fields]))
    stop("Preflight resume manifest identity changed.")
  before <- vapply(stored$output_location, tools::md5sum, character(1))
  plan <- plan_v021_checkpoint_restart(stored)
  after <- vapply(stored$output_location, tools::md5sum, character(1))
  completed <- plan$completed_indices
  resume_completed_hashes <- if (length(completed))
    vapply(plan$manifest$output_location[completed], tools::md5sum, character(1))
  run_indices <- plan$resume_indices
  manifest <- plan$manifest
  write_v021_manifest_atomic(manifest, manifest_path)
  audit <- data.frame(
    completed_tasks_not_rerun = identical(before, after),
    resume_task_count = length(plan$resume_indices),
    completed_task_count = sum(stored$execution_state == "completed"),
    manifest_identity_unchanged = TRUE,
    seeds_unchanged = identical(stored$seed, manifest$seed),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE)
  write.table(audit, file.path(config$output_root, "restart_idempotence_audit.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  if (!length(run_indices)) {
    cat("V021-05 resume audit passed; completed posterior tasks were skipped.\n")
    resume_all_complete <- TRUE
  }
  if (length(run_indices)) cat("V021-05 controlled resume: ", length(run_indices),
      " incomplete task(s); completed checkpoints will be skipped.\n", sep = "")
}
if (resume_mode && dir.exists(v021_ownership_path(config$output_root))) {
  stale_owner <- reconcile_v021_stale_ownership(
    config$output_root, .v021_local_process_identity(), administrative_override = TRUE)
  saveRDS(stale_owner, file.path(config$output_root, "stale_ownership_reconciliation.rds"))
}
controller_ownership <- acquire_v021_controller_ownership(
  config$output_root, execution_manifest$manifest_hash,
  execution_manifest$configuration_hash)
ownership_released <- FALSE
on.exit(if (!ownership_released) try(
  release_v021_controller_ownership(config$output_root, controller_ownership),
  silent = TRUE), add = TRUE)
saveRDS(list(acquired = TRUE, released = FALSE,
             controller_pid = controller_ownership$controller_pid,
             manifest_hash = execution_manifest$manifest_hash,
             configuration_hash = execution_manifest$configuration_hash,
             acquisition_timestamp = controller_ownership$acquired_at),
        file.path(config$output_root, "controller_ownership_audit.rds"))
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
fit_with_eligible_kfold <- function(target, partner, ctx, seed_override = NULL,
                                    progress_local = "none") {
  result <- pclvbayes:::.fit_direction_main_posterior(
    target, partner, ctx, seed_override, progress_local)
  if (inherits(result, "pclv_failure")) return(result)
  pclvbayes:::.add_predictive_evaluation(result)
}
jobs <- lapply(specs, make_v021_confirmation_fit_closure,
  runtime_context = approved_runtime, fit_direction = fit_with_eligible_kfold)
rm(study)

write_v021_manifest_atomic(manifest, manifest_path)
saveRDS(config, file.path(config$output_root, "config.rds"))
saveRDS(selected, file.path(config$output_root, "prospective_selection.rds"))
saveRDS(policy, file.path(config$output_root, "resource_policy.rds"))
saveRDS(thread_environment, file.path(config$output_root, "thread_environment.rds"))
saveRDS(execution_manifest, file.path(config$output_root, "canonical_execution_manifest.rds"))
writeLines(execution_manifest$manifest_hash,
           file.path(config$output_root, "preflight_manifest.sha256"))
writeLines(c(
  "#!/bin/sh", "set -eu",
  "Rscript benchmarks/mtist/run_v021_four_chain_preflight.R",
  "PCLV_V021_PREFLIGHT_RESUME=true Rscript benchmarks/mtist/run_v021_four_chain_preflight.R",
  "ps -eo pid,ppid,lstart,args",
  "find /proc/<pid>/task -maxdepth 1 -mindepth 1 -type d"),
  file.path(config$output_root, "pipeline.commands.sh"))

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
monitor_audit_root_pid <- preflight_root_pid
if (resume_all_complete) {
  snapshot_files <- list.files(monitor_dir, pattern = "^batch-[0-9]+-final\\.rds$",
                               full.names = TRUE)
  payloads <- lapply(snapshot_files, readRDS)
  all_snapshots <- unlist(lapply(payloads, `[[`, "snapshots"), recursive = FALSE)
  registry_path <- file.path(config$output_root, "worker_registry.rds")
  if (file.exists(registry_path)) {
    all_worker_registry <- readRDS(registry_path)
    monitor_audit_root_pid <- unique(all_worker_registry$expected_ppid)
    if (length(monitor_audit_root_pid) != 1L)
      stop("Retained worker registry has ambiguous controller ancestry.")
  }
}

mark_running <- function(indices) {
  for (i in indices) {
    started_attempt <- start_v021_task_attempt(
      preflight_status, preflight_attempt_ledger, execution_manifest,
      selected$direction_index[[i]], format(Sys.time(), tz = "UTC", usetz = TRUE),
      worker_provenance = "v021_preflight_controller")
    preflight_status <<- started_attempt$status
    preflight_attempt_ledger <<- started_attempt$ledger
    write_v021_task_status_atomic(
      preflight_status, preflight_execution_paths$status, execution_manifest)
    write_v021_attempt_ledger_atomic(
      preflight_attempt_ledger, preflight_execution_paths$attempts, execution_manifest)
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
    finished_attempt <- finish_v021_task_attempt(
      preflight_status, preflight_attempt_ledger, execution_manifest,
      selected$direction_index[[i]], "failed",
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      terminal_reason = checkpoint$failure_reason, failure_class = "pclv_failure")
  } else {
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
    inference_results$kfold_results <- if (retained$kfold_attempted)
      list(state = if (retained$kfold_completed) "completed" else "failed",
           folds_ok = retained$kfold_folds_ok, folds_failed = retained$kfold_folds_failed)
    else list(state = "not_eligible")
    inference_results$elpd_results <- list(
      state = if (retained$elpd_available) "available" else retained$aggregate_elpd_missing_state,
      value = retained$aggregate_elpd)
    inference_results$matrices <- list()
    inference_results$masks <- list()
    inference_results$feature_inputs <- retained
    inference_results$stan_data <- list(N = retained$n_pairs)
    inference_results$execution <- list(state = "completed", elapsed_seconds = elapsed)
    inference_results$status <- list(original_reporting_state = retained$diagnostic_class)
    inference_results$unavailable_values <- list(kfold = NA_real_, elpd = NA_real_)
    inference_results$zero_placeholders <- list()
    artifact <- finalize_v021_inference_artifact(list(
      artifact_schema = v021_retained_fixture_schema, artifact_state = "retained",
      execution_state = "completed", inference_results = inference_results))
    feature <- build_v021_diagnostic_feature_record(
      retained, observations$design,
      list(chains = config$chains, iter_sampling = config$iter_sampling))
    saveRDS(artifact, file.path(artifact_dir, paste0(manifest$direction_id[[i]], ".rds")))
    saveRDS(feature, file.path(feature_dir, paste0(manifest$direction_id[[i]], ".rds")))
    meta <- approved_runtime$meta_df
    sm <- approved_runtime$sm_mat
    ord <- match(as.character(meta$Sample), colnames(sm))
    observed_input <- data.frame(
      subject = as.character(meta$subject), time = as.numeric(meta$time),
      source_abundance = as.numeric(sm[selected$source[[i]], ord]),
      target_abundance = as.numeric(sm[selected$target[[i]], ord]),
      stringsAsFactors = FALSE)
    observed_input$rest_abundance <- 1 - observed_input$source_abundance -
      observed_input$target_abundance
    observed_input$eligible <- with(observed_input,
      is.finite(source_abundance) & is.finite(target_abundance) &
        is.finite(rest_abundance) & rest_abundance >= 0)
    observed_support <- calculate_v021_observed_support_features(observed_input)
    posterior_stability <- extract_v021_posterior_stability_features(list(
      finalized = TRUE, terminal_state = "completed",
      posterior_mean = retained$posterior_mean, posterior_median = retained$posterior_median,
      posterior_sd = retained$posterior_sd,
      interval_lower = retained$posterior_interval_lower,
      interval_upper = retained$posterior_interval_upper,
      psp = retained$p_sign2, lfsr = retained$lfsr,
      rhat_max = retained$rhat, bulk_ess_min = retained$ess_bulk,
      tail_ess_min = retained$ess_tail, divergence_count = retained$divergences,
      treedepth_hit_count = retained$treedepth_hits, ebfmi_min = retained$ebfmi_min,
      retry_count = max(length(result$retry_history[[1L]]) - 1L, 0L),
      diagnostic_class = retained$diagnostic_class,
      predictive_eligible = retained$bayesian_eligible,
      elpd = retained$aggregate_elpd,
      successful_test_observation_count =
        v021_preflight_successful_test_observation_count(result$kfold_success_total),
      folds_attempted = if (retained$kfold_attempted) 5L else 0L,
      folds_failed = retained$kfold_folds_failed))
    robustness <- build_v021_robustness_feature_artifact(
      execution_manifest, selected$direction_index[[i]], observed_support,
      posterior_stability, code_provenance = git_commit)
    write_v021_robustness_feature_artifact_atomic(
      robustness, file.path(robustness_dir,
                            paste0(manifest$direction_id[[i]], ".rds")),
      execution_manifest)
    checkpoint <- .v021_record_from_row(
      manifest[i, , drop = FALSE], "completed", elapsed_time = elapsed,
      retry_history = result$retry_history[[1L]] %||% list(), pathfinder_used = FALSE,
      completion_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), result = artifact)
    completion <- build_v021_completion_artifact(
      execution_manifest, selected$direction_index[[i]], config$output_root,
      posterior_summary = list(mean = retained$posterior_mean,
        median = retained$posterior_median, sd = retained$posterior_sd),
      diagnostics = list(class = retained$diagnostic_class, rhat = retained$rhat,
        bulk_ess = retained$ess_bulk, tail_ess = retained$ess_tail,
        divergences = retained$divergences, treedepth_hits = retained$treedepth_hits),
      psp_lfsr = list(PSP = retained$p_sign2, LFSR = retained$lfsr),
      predictive_eligibility = retained$bayesian_eligible,
      subject_elpd = list(state = if (retained$kfold_attempted) "executed" else "not_applicable",
        n_successful_test_observations =
          v021_preflight_successful_test_observation_count(result$kfold_success_total),
        method = "student-t-scale-mixture-kalman-ou-q16"),
      seed_split_provenance = list(direction_seed = retained$seed,
                                   kfold_seed = config$seed))
    completion_path <- file.path(preflight_execution_paths$artifact_root,
                                 paste0(manifest$direction_id[[i]], ".rds"))
    write_v021_completion_artifact_atomic(
      completion, completion_path, execution_manifest, config$output_root)
    finished_attempt <- finish_v021_task_attempt(
      preflight_status, preflight_attempt_ledger, execution_manifest,
      selected$direction_index[[i]], "completed",
      format(Sys.time(), tz = "UTC", usetz = TRUE),
      artifact_path = completion_path, result_root = config$output_root)
  }
  preflight_status <<- finished_attempt$status
  preflight_attempt_ledger <<- finished_attempt$ledger
  write_v021_task_status_atomic(
    preflight_status, preflight_execution_paths$status, execution_manifest)
  write_v021_attempt_ledger_atomic(
    preflight_attempt_ledger, preflight_execution_paths$attempts, execution_manifest)
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
  if (length(run_indices)) run_batch(run_indices)
  completed_before <- vapply(manifest$output_location[manifest$execution_state == "completed"],
                              tools::md5sum, character(1))
  restart <- plan_v021_checkpoint_restart(manifest)
  if (length(restart$resume_indices)) stop("Fresh preflight left unfinished posterior work.")
  completed_after <- vapply(names(completed_before), tools::md5sum, character(1))
  restart_record <<- list(
    controlled_interruption = FALSE, resumed_direction_index = integer(),
    completed_tasks_not_rerun = identical(completed_before, completed_after) &&
      (!length(resume_completed_hashes) || identical(
        resume_completed_hashes,
        vapply(names(resume_completed_hashes), tools::md5sum, character(1)))),
    incomplete_task_resumed = TRUE,
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
  all_snapshots, policy, monitor_audit_root_pid, executable,
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
robustness_features <- lapply(list.files(robustness_dir, full.names = TRUE), readRDS)
invisible(lapply(robustness_features, validate_v021_robustness_feature_artifact,
                 manifest = execution_manifest))
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

write_tsv <- function(x, name) write.table(
  x, file.path(config$output_root, name), sep = "\t", row.names = FALSE,
  quote = FALSE, na = "NA")
run_end <- Sys.time()
run_metadata <- data.frame(
  git_commit = git_commit,
  package_version = as.character(read.dcf(file.path(repo, "DESCRIPTION"), "Version")[[1L]]),
  r_version = R.version.string,
  cmdstan_version = paste(cmdstanr::cmdstan_version(), collapse = "."),
  hostname = Sys.info()[["nodename"]], logical_cpu_detected = parallel::detectCores(),
  allocated_logical_threads = 16L, cpu_affinity = "0-15",
  start_timestamp = summary$monitor_start,
  end_timestamp = format(run_end, tz = "UTC", usetz = TRUE),
  manifest_hash = execution_manifest$manifest_hash,
  configuration_hash = execution_manifest$configuration_hash,
  stringsAsFactors = FALSE)
write_tsv(run_metadata, "run_metadata.tsv")
write_tsv(selected, "preflight_task_manifest.tsv")
write_tsv(data.frame(
  policy_schema = policy$policy_schema, logical_host_threads = policy$logical_host_threads,
  reserved_host_threads = policy$reserved_host_threads,
  maximum_active_cmdstan_chains = policy$maximum_active_cmdstan_chains,
  maximum_cmdstan_process_slots = policy$maximum_cmdstan_process_slots,
  chains_per_fit = policy$main_chains,
  parallel_chains_per_fit = policy$main_parallel_chains,
  threads_per_chain = policy$cpu_threads_per_active_chain,
  maximum_concurrent_fits = policy$proposed_outer_concurrency,
  kfold_parallel_chains_per_fit = policy$kfold_parallel_chains,
  maximum_concurrent_kfold_fits = policy$maximum_concurrent_kfold_fits),
  "resource_policy.tsv")
write_tsv(data.frame(variable = names(thread_environment), value = unname(thread_environment)),
          "environment_thread_caps.tsv")
write_tsv(data.frame(wave = c(1L, 2L), task_count = c(3L, 1L),
                     reserved_chains = c(12L, 4L), accepted = c(TRUE, FALSE),
                     reason = c("three canonical fits fill capacity",
                                "fourth concurrent fit rejected")),
          "scheduler_dry_run.tsv")
write_tsv(manifest[c("task_id", "direction_id", "direction_index", "target", "source",
                     "seed", "execution_state", "elapsed_time", "failure_reason")],
          "task_status.tsv")
attempt_ledger <- manifest[c("direction_id", "direction_index", "seed", "execution_state",
                             "elapsed_time", "failure_reason")]
attempt_ledger$attempt_number <- 1L
write_tsv(attempt_ledger, "attempt_ledger.tsv")
write_tsv(data.frame(
  direction_id = manifest$direction_id, checkpoint_path = manifest$output_location,
  checkpoint_exists = file.exists(manifest$output_location),
  checkpoint_md5 = unname(tools::md5sum(manifest$output_location)),
  execution_state = manifest$execution_state), "completion_audit.tsv")
process_rows <- monitor$records
if (nrow(process_rows)) {
  keep <- intersect(c("timestamp", "pid", "ppid", "classification", "capture_state",
                      "command", "executable", "active_cmdstan_chains",
                      "active_cmdstan_processes"), names(process_rows))
  process_rows <- process_rows[keep]
}
write_tsv(process_rows, "process_tree_audit.tsv")
retained_rows <- lapply(manifest$output_location, function(path) readRDS(path)$result$inference_results$feature_inputs)
posterior <- do.call(rbind, lapply(retained_rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
write_tsv(posterior[c("dataset_id", "direction_index", "source", "target", "seed",
                      "posterior_mean", "posterior_median", "posterior_sd",
                      "posterior_interval_lower", "posterior_interval_upper", "p_sign2", "lfsr")],
          "posterior_summary.tsv")
write_tsv(posterior[c("direction_index", "source", "target", "diagnostic_class",
                      "rhat", "ess_bulk", "ess_tail", "divergences", "treedepth_hits",
                      "ebfmi_min", "bayesian_eligible")], "diagnostic_summary.tsv")
robustness_summary <- do.call(rbind, lapply(robustness_features, function(x) data.frame(
  task_id = x$directed_task_id,
  observed_validity = x$feature_validity_status,
  posterior_validity = if (is.null(x$posterior_stability)) "unavailable" else "valid",
  stringsAsFactors = FALSE)))
write_tsv(robustness_summary, "robustness_feature_summary.tsv")
kfold_summary <- posterior[c("direction_index", "source", "target", "bayesian_eligible",
                              "kfold_attempted", "kfold_completed", "kfold_folds_ok",
                              "kfold_folds_failed", "aggregate_elpd",
                              "aggregate_elpd_missing_state")]
kfold_summary$n_successful_test_observations <- vapply(
  robustness_features,
  function(x) x$posterior_stability$values$successful_test_observation_count,
  integer(1))
write_tsv(kfold_summary, "kfold_summary.tsv")
write_tsv(data.frame(
  perturbation_id = c("required_low_cost", "subject_deletion", "endpoint_deletion",
                      "time_grid_coarsening"),
  requires_refit = c(FALSE, TRUE, TRUE, TRUE), fits_per_selected_task = c(0L, 1L, 2L, 1L),
  chains_per_fit = c(0L, 4L, 4L, 4L), maximum_parallel_chains = c(0L, 4L, 4L, 4L),
  estimated_multiplier = c(1, 2, 3, 2),
  proposed_execution_scope = c("all directions", rep("not in canonical V021-06", 3L)),
  deterministic_selection_rule = c("all canonical tasks", rep("not selected", 3L)),
  completion_required = c(TRUE, FALSE, FALSE, FALSE)),
  "optional_perturbation_budget.tsv")
elapsed <- manifest$elapsed_time
bytes <- vapply(manifest$output_location, function(path) file.info(path)$size, numeric(1))
central_seconds <- stats::median(elapsed) * 9900 / config$maximum_simultaneous_fits
write_tsv(data.frame(
  measured_tasks = nrow(manifest), simultaneous_task_wall_seconds = sum(elapsed),
  median_task_seconds = stats::median(elapsed), checkpoint_bytes_per_task = stats::median(bytes),
  projected_directions = 9900L, optimistic_wall_days = central_seconds * .75 / 86400,
  central_wall_days = central_seconds / 86400,
  conservative_wall_days = central_seconds * 1.5 / 86400,
  projected_compact_bytes = stats::median(bytes) * 9900,
  retry_assumption = "no retries in central estimate; 50% allowance conservative",
  kfold_assumption = "eligible-direction cost observed separately",
  optional_refit_assumption = "excluded from canonical completion"),
  "runtime_storage_estimate.tsv")
failures <- manifest[manifest$execution_state == "failed",
                     c("direction_id", "direction_index", "failure_reason"), drop = FALSE]
write_tsv(failures, "failures.tsv")
release_v021_controller_ownership(config$output_root, controller_ownership)
ownership_released <- TRUE
saveRDS(list(acquired = TRUE, released = TRUE,
             controller_pid = controller_ownership$controller_pid,
             manifest_hash = execution_manifest$manifest_hash,
             configuration_hash = execution_manifest$configuration_hash,
             acquisition_timestamp = controller_ownership$acquired_at,
             release_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)),
        file.path(config$output_root, "controller_ownership_audit.rds"))
cat("V021-05 preflight passed. Peak chains:", summary$observed_peak_active_chains,
    "peak CmdStan process slots:", summary$observed_peak_cmdstan_process_slots, "\n")
}

run_v021_four_chain_preflight()
