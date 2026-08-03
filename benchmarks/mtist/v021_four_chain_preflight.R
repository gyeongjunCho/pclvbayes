# Benchmark-only bounded final-configuration preflight for ROADMAP V021-05.

v021_preflight_schema <- "v021_four_chain_preflight_v1"
v021_preflight_selection_rule <- paste(
  "MTIST dataset 37; canonical direction index from seed 20260802;",
  "first four direction_index rows; no truth or outcome fields inspected"
)

build_v021_four_chain_preflight_config <- function(output_root) {
  if (!is.character(output_root) || length(output_root) != 1L || is.na(output_root) ||
      !nzchar(output_root)) stop("output_root must be one path.")
  list(
    preflight_schema = v021_preflight_schema,
    dataset_id = "37", seed = 20260802L, selected_direction_count = 4L,
    selection_rule = v021_preflight_selection_rule,
    chains = 4L, iter_warmup = 2000L, iter_sampling = 2000L,
    nominal_retained_draws = 8000L, maximum_simultaneous_fits = 3L,
    run_kfold = FALSE, use_pathfinder = FALSE, controlled_interrupt_direction = 4L,
    output_root = normalizePath(output_root, mustWork = FALSE)
  )
}

validate_v021_four_chain_preflight_config <- function(config) {
  expected <- c(
    "preflight_schema", "dataset_id", "seed", "selected_direction_count",
    "selection_rule", "chains", "iter_warmup", "iter_sampling",
    "nominal_retained_draws", "maximum_simultaneous_fits", "run_kfold",
    "use_pathfinder", "controlled_interrupt_direction", "output_root"
  )
  if (!is.list(config) || !identical(names(config), expected) ||
      !identical(config$preflight_schema, v021_preflight_schema))
    stop("Invalid V021-05 preflight configuration.")
  if (!identical(config$dataset_id, "37") || !identical(config$seed, 20260802L) ||
      !identical(config$selected_direction_count, 4L) ||
      !identical(config$selection_rule, v021_preflight_selection_rule) ||
      !identical(config$chains, 4L) || !identical(config$iter_warmup, 2000L) ||
      !identical(config$iter_sampling, 2000L) ||
      !identical(config$nominal_retained_draws, 8000L) ||
      !identical(config$maximum_simultaneous_fits, 3L) ||
      !identical(config$run_kfold, FALSE) || !identical(config$use_pathfinder, FALSE) ||
      !identical(config$controlled_interrupt_direction, 4L))
    stop("Preflight configuration does not match the frozen V021-05 contract.")
  invisible(TRUE)
}

validate_v021_preflight_launch_capacity <- function(policy, fit_count) {
  fit_count <- .v021_resource_int(fit_count, "fit_count")
  derivation <- derive_safe_outer_concurrency(
    policy, build_v021_operation_spec("main_fit", fit_count))
  if (derivation$projected_active_cmdstan_chains >
      policy$maximum_active_cmdstan_chains ||
      derivation$projected_active_cmdstan_processes >
      policy$maximum_cmdstan_process_slots)
    stop("Preflight batch exceeds the configured chain or process-slot ceiling.")
  derivation
}

register_v021_preflight_worker <- function(pid, expected_ppid, worker_role,
                                            batch_id, task_id,
                                            proc_root = "/proc",
                                            stat_reader = function(path)
                                              readLines(path, warn = FALSE, n = 1L),
                                            clock = function()
                                              format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  pid <- .v021_resource_int(pid, "pid")
  expected_ppid <- .v021_resource_int(expected_ppid, "expected_ppid")
  stat <- tryCatch(stat_reader(file.path(proc_root, as.character(pid), "stat")),
                   error = identity)
  if (inherits(stat, "error") || !length(stat))
    stop("Outer-worker procfs identity is unavailable at registration.")
  identity <- .v021_parse_linux_stat(stat[[1L]], pid)
  if (!identical(identity$ppid, expected_ppid) ||
      !identity$process_state %in% .v021_linux_non_zombie_states)
    stop("Outer-worker PPID, ancestry, or process state is invalid at registration.")
  build_v021_worker_registry(
    pid, identity$start_time, expected_ppid, worker_role, batch_id, task_id,
    as.character(clock()))
}

select_v021_preflight_tasks <- function(task_table, count = 4L) {
  required <- c("task_id", "direction_index", "target", "source", "seed")
  if (!is.data.frame(task_table) || !all(required %in% names(task_table)))
    stop("task_table lacks canonical direction identity.")
  forbidden <- names(task_table)[.v021_truth_name(names(task_table)) |
    grepl("calibrat|withhold|same_nonzero|opposite_nonzero|absolute_zero",
          names(task_table), ignore.case = TRUE)]
  if (length(forbidden)) stop("Truth or future-policy fields are forbidden in preflight selection.")
  if (!identical(order(task_table$direction_index), seq_len(nrow(task_table))) ||
      anyDuplicated(task_table$direction_index))
    stop("task_table is not in canonical direction order.")
  count <- as.integer(count)
  if (length(count) != 1L || is.na(count) || count < 1L || count > nrow(task_table))
    stop("Invalid preflight task count.")
  out <- task_table[seq_len(count), required, drop = FALSE]
  out$dataset_id <- "37"
  out <- out[c("dataset_id", required)]
  rownames(out) <- NULL
  out
}

load_v021_preflight_observations <- function(root, dataset_id = 37L) {
  root <- normalizePath(root, mustWork = TRUE)
  path <- file.path(root, "mtist1.0", "mtist_datasets",
                    paste0("dataset_", as.integer(dataset_id), ".csv"))
  raw <- .read_indexed_csv(path)
  species <- names(raw)[grepl("^species_[0-9]+$", names(raw))]
  if (!length(species) || !all(c("time", "timeseries_id") %in% names(raw)))
    stop("Observational dataset lacks required fields.")
  abundance <- as.matrix(raw[species]); storage.mode(abundance) <- "double"
  if (any(!is.finite(abundance)) || any(abundance < 0) || any(rowSums(abundance) <= 0))
    stop("Invalid observational abundances.")
  relative <- abundance / rowSums(abundance)
  sample_id <- sprintf("did%s_ts%s_n%03d", dataset_id, raw$timeseries_id,
                       ave(seq_len(nrow(raw)), raw$timeseries_id, FUN = seq_along))
  rownames(relative) <- sample_id
  sample_meta <- data.frame(
    subject = as.character(raw$timeseries_id), time = as.numeric(raw$time),
    dataset_id = as.integer(dataset_id), row.names = sample_id,
    stringsAsFactors = FALSE
  )
  physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(t(relative), taxa_are_rows = TRUE),
    phyloseq::sample_data(sample_meta)
  )
  list(physeq = physeq, taxa = species,
       design = data.frame(subject = sample_meta$subject, time = sample_meta$time))
}

capture_v021_preflight_process_snapshot <- function(
    root_pid = Sys.getpid(), known_model_executables = character()) {
  capture_v021_linux_process_snapshot(
    root_pid = root_pid, known_model_executables = known_model_executables)
}

validate_v021_preflight_monitor <- function(monitor, policy) {
  if (!is.list(monitor) || !identical(monitor$monitoring_state, "verified"))
    stop("Complete process-tree monitoring was not verified.")
  if (!identical(monitor$compliance_status, "compliant") ||
      monitor$observed_peak_active_cmdstan_chains > policy$maximum_active_cmdstan_chains ||
      monitor$observed_peak_active_cmdstan_processes > policy$usable_chain_slots)
    stop("Observed CmdStan consumption exceeded the resource ceiling.")
  invisible(TRUE)
}

v021_preflight_retained_record <- function(task, fit_result) {
  if (!is.list(fit_result) || inherits(fit_result, "pclv_failure"))
    stop("A retained directional fit result is required.")
  diag <- fit_result$diag %||% list()
  list(
    dataset_id = as.character(task$dataset_id[[1L]]), pair_id = task$pair_id[[1L]],
    task_id = task$task_id[[1L]], direction_index = task$direction_index[[1L]],
    target = task$target[[1L]], source = task$source[[1L]], seed = task$seed[[1L]],
    posterior_mean = fit_result$a_mean, posterior_median = fit_result$a_median,
    posterior_sd = fit_result$a_sd, posterior_interval_lower = fit_result$a_q2.5,
    posterior_interval_upper = fit_result$a_q97.5,
    posterior_sign_probability = max(fit_result$positive_sign_probability,
                                     fit_result$negative_sign_probability),
    p_sign2 = fit_result$p_sign2, lfsr = fit_result$lfsr,
    rhat = diag$worst_rhat %||% NA_real_, ess_bulk = diag$min_ess_bulk %||% NA_real_,
    ess_tail = diag$min_ess_tail %||% NA_real_, divergences = diag$n_divergent %||% NA_integer_,
    treedepth_hits = diag$n_treedepth_hit %||% NA_integer_,
    ebfmi_min = diag$ebfmi_min %||% NA_real_,
    chain_sign_agreement = fit_result$chain_sign_agreement,
    diagnostic_class = fit_result$diagnostic_class,
    interaction_identifiable = fit_result$interaction_identifiable,
    residual_identifiable = fit_result$residual_identifiable,
    bayesian_eligible = identical(fit_result$diagnostic_class, "converged"),
    residual_regime_disagreement = fit_result$residual_regime_disagreement,
    kfold_attempted = FALSE, kfold_completed = FALSE, kfold_folds_ok = 0L,
    kfold_folds_failed = 0L, elpd_available = FALSE,
    aggregate_elpd = NA_real_, aggregate_elpd_missing_state = "not_executed",
    stacking_available = FALSE, stacking_weight = NA_real_,
    stacking_weight_missing_state = "not_executed", n_pairs = fit_result$n_pairs
  )
}

build_v021_preflight_summary <- function(config, selected, policy, monitor,
                                          manifest, executable_record,
                                          restart_record, feature_records,
                                          truth_absent = TRUE,
                                          orphan_process_check = NULL) {
  diagnostic_records <- if (nrow(monitor$records))
    monitor$records[monitor$records$classification == "cmdstan_diagnostic", , drop = FALSE]
  else monitor$records
  diagnostic_peak <- if (!nrow(diagnostic_records)) 0L else as.integer(max(
    vapply(split(diagnostic_records$pid, diagnostic_records$timestamp),
           function(pid) length(unique(pid)), integer(1))))
  vanished_records <- if (nrow(monitor$records))
    monitor$records[monitor$records$classification == "vanished_during_capture", , drop = FALSE]
  else monitor$records
  zombie_records <- if (nrow(monitor$records))
    monitor$records[monitor$records$classification == "zombie_process", , drop = FALSE]
  else monitor$records
  recapture_records <- if (nrow(monitor$records) &&
                             "capture_retry_count" %in% names(monitor$records))
    monitor$records[monitor$records$capture_retry_count > 0L, , drop = FALSE]
  else monitor$records[FALSE, , drop = FALSE]
  resource <- evaluate_v021_resource_preflight(
    policy, build_v021_operation_spec("main_fit", config$maximum_simultaneous_fits),
    v021_single_thread_environment(), monitor)
  terminal <- manifest$execution_state %in% c("completed", "failed", "skipped")
  passed <- identical(resource$state, "passed") && all(terminal) &&
    isTRUE(executable_record$reuse_verified) &&
    isFALSE(executable_record$worker_side_compilation) &&
    isTRUE(restart_record$completed_tasks_not_rerun) && isTRUE(truth_absent) &&
    length(feature_records) == sum(manifest$execution_state == "completed") &&
    isTRUE(orphan_process_check$passed)
  list(
    summary_schema = "v021_four_chain_preflight_summary_v2",
    resource_policy_schema = policy$policy_schema,
    logical_host_threads = policy$logical_host_threads,
    reserved_host_threads = policy$reserved_host_threads,
    usable_execution_capacity = policy$usable_chain_slots,
    maximum_active_cmdstan_chains = policy$maximum_active_cmdstan_chains,
    maximum_cmdstan_process_slots = policy$maximum_cmdstan_process_slots,
    attempt_id = basename(config$output_root),
    state = if (passed) "passed" else "failed",
    failure_reasons = if (passed) character() else c(resource$reasons,
      if (!all(terminal)) "nonterminal_checkpoint_state",
      if (!isTRUE(executable_record$reuse_verified)) "executable_reuse_unverified",
      if (isTRUE(executable_record$worker_side_compilation)) "worker_side_compilation",
      if (!isTRUE(restart_record$completed_tasks_not_rerun)) "completed_task_rerun",
      if (!isTRUE(truth_absent)) "truth_isolation_failed",
      if (!isTRUE(orphan_process_check$passed)) "orphan_preflight_process"),
    selected_task_ids = selected$task_id,
    selected_direction_indices = selected$direction_index,
    selected_directions = paste(selected$source, selected$target, sep = "->"),
    direction_seeds = selected$seed,
    chain_seeds = manifest$chain_seeds,
    selection_rule = config$selection_rule, chains = config$chains,
    iter_warmup = config$iter_warmup, iter_sampling = config$iter_sampling,
    nominal_retained_draws = config$nominal_retained_draws,
    maximum_simultaneous_fits = config$maximum_simultaneous_fits,
    projected_active_chains = as.integer(config$maximum_simultaneous_fits * config$chains),
    projected_cmdstan_process_slots = as.integer(
      config$maximum_simultaneous_fits * config$chains),
    output_root = config$output_root,
    dataset_id = config$dataset_id,
    executable_identity = executable_record$executable_path,
    compilation_count = executable_record$compilation_count,
    worker_compilation_count = executable_record$worker_compilation_count,
    thread_environment = v021_single_thread_environment(),
    monitor_start = if (nrow(monitor$records)) min(monitor$records$timestamp) else NA_character_,
    monitor_end = if (nrow(monitor$records)) max(monitor$records$timestamp) else NA_character_,
    monitoring_state = monitor$monitoring_state,
    observed_peak_active_chains = monitor$observed_peak_active_cmdstan_chains,
    observed_peak_cmdstan_process_slots = monitor$observed_peak_active_cmdstan_processes,
    observed_peak_cmdstan_diagnostic_processes = diagnostic_peak,
    cmdstan_diagnostic_process_ids = sort(unique(as.integer(diagnostic_records$pid))),
    cmdstan_diagnostic_executables = sort(unique(diagnostic_records$executable)),
    vanished_process_ids = sort(unique(as.integer(vanished_records$pid))),
    vanished_process_reasons = sort(unique(as.character(
      vanished_records$disappearance_reason))),
    zombie_process_ids = sort(unique(as.integer(zombie_records$pid))),
    zombie_process_reasons = sort(unique(as.character(zombie_records$zombie_reason))),
    transient_recapture_process_ids = sort(unique(as.integer(recapture_records$pid))),
    transient_recapture_event_count = as.integer(nrow(recapture_records)),
    transient_recapture_resolution_reasons = sort(unique(as.character(
      recapture_records$resolution_reason))),
    resource_state = resource$state,
    ceiling_respected = identical(resource$state, "passed"),
    all_thread_variables_one = TRUE,
    process_tree_monitoring_complete = identical(monitor$monitoring_state, "verified"),
    worker_side_compilation = executable_record$worker_side_compilation,
    restart_avoided_completed_tasks = restart_record$completed_tasks_not_rerun,
    executable_reuse_verified = executable_record$reuse_verified,
    checkpoint_states = manifest$execution_state,
    checkpoint_paths = manifest$output_location,
    manifest_path = file.path(config$output_root, "manifest.rds"),
    manifest_reconciled = isTRUE(restart_record$manifest_reconciled),
    interruption_record = restart_record[c("controlled_interruption", "resumed_direction_index")],
    restart_record = restart_record,
    incomplete_task_resumed = isTRUE(restart_record$incomplete_task_resumed),
    feature_records_validated = length(feature_records) ==
      sum(manifest$execution_state == "completed"),
    orphan_process_check = orphan_process_check,
    truth_absent = truth_absent,
    per_direction_states = data.frame(
      task_id = manifest$task_id, direction_index = manifest$direction_index,
      execution_state = manifest$execution_state,
      scientific_state = ifelse(manifest$execution_state == "failed",
                                "fit_failed", "retained_or_not_applicable"),
      stringsAsFactors = FALSE)
  )
}

validate_v021_preflight_summary <- function(summary) {
  required <- c(
    "summary_schema", "resource_policy_schema", "logical_host_threads",
    "reserved_host_threads", "usable_execution_capacity",
    "maximum_active_cmdstan_chains", "maximum_cmdstan_process_slots",
    "attempt_id", "state", "failure_reasons", "selected_task_ids",
    "selected_direction_indices", "selected_directions", "direction_seeds", "chain_seeds",
    "selection_rule", "chains",
    "iter_warmup", "iter_sampling", "nominal_retained_draws",
    "maximum_simultaneous_fits", "projected_active_chains",
    "projected_cmdstan_process_slots", "output_root", "dataset_id",
    "executable_identity", "compilation_count", "worker_compilation_count",
    "thread_environment", "monitor_start", "monitor_end", "monitoring_state",
    "observed_peak_active_chains", "observed_peak_cmdstan_process_slots",
    "observed_peak_cmdstan_diagnostic_processes", "cmdstan_diagnostic_process_ids",
    "cmdstan_diagnostic_executables", "vanished_process_ids",
    "vanished_process_reasons", "zombie_process_ids", "zombie_process_reasons",
    "transient_recapture_process_ids", "transient_recapture_event_count",
    "transient_recapture_resolution_reasons",
    "resource_state", "ceiling_respected", "all_thread_variables_one",
    "process_tree_monitoring_complete", "worker_side_compilation",
    "restart_avoided_completed_tasks", "executable_reuse_verified", "checkpoint_states",
    "checkpoint_paths", "manifest_path", "manifest_reconciled",
    "interruption_record", "restart_record", "incomplete_task_resumed",
    "feature_records_validated", "orphan_process_check", "truth_absent",
    "per_direction_states"
  )
  if (!is.list(summary) || !identical(names(summary), required) ||
      !identical(summary$summary_schema, "v021_four_chain_preflight_summary_v2") ||
      !identical(summary$resource_policy_schema, "v021_resource_policy_v2") ||
      !identical(summary$logical_host_threads, 16L) ||
      !identical(summary$reserved_host_threads, 4L) ||
      !identical(summary$usable_execution_capacity, 12L) ||
      !identical(summary$maximum_active_cmdstan_chains, 12L) ||
      !identical(summary$maximum_cmdstan_process_slots, 12L) ||
      !summary$state %in% c("passed", "failed")) stop("Invalid preflight summary.")
  if (identical(summary$state, "passed") && length(summary$failure_reasons))
    stop("Passed preflight cannot retain failure reasons.")
  if (identical(summary$state, "failed") && !length(summary$failure_reasons))
    stop("Failed preflight requires explicit reasons.")
  if (!is.integer(summary$compilation_count) || summary$compilation_count < 0L ||
      !is.integer(summary$worker_compilation_count) || summary$worker_compilation_count != 0L)
    stop("Worker-side compilation must be explicitly zero.")
  validate_v021_single_thread_environment(summary$thread_environment)
  invisible(TRUE)
}
