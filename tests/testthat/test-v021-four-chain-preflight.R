source(testthat::test_path("../../benchmarks/mtist/mtist_adapter.R"))
source(testthat::test_path("../../benchmarks/mtist/ten_species_helpers.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_truth_isolation.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_resource_policy.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_checkpoint_manifest.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_diagnostic_features.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_four_chain_preflight.R"))

test_that("preflight configuration and selection are deterministic and truth-free", {
  a <- build_v021_four_chain_preflight_config(tempdir())
  b <- build_v021_four_chain_preflight_config(tempdir())
  expect_identical(a, b)
  expect_silent(validate_v021_four_chain_preflight_config(a))
  index <- ten_species_direction_index(paste0("species_", 0:9), a$seed)
  selected <- select_v021_preflight_tasks(index)
  expect_identical(selected$direction_index, 1:4)
  expect_identical(selected$task_id, c(1L, 1L, 2L, 2L))
  bad <- index; bad$truth_sign <- 1
  expect_error(select_v021_preflight_tasks(bad), "Truth")
})

test_that("four-chain concurrency is bounded at three under policy v2", {
  policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
  derivation <- validate_v021_preflight_launch_capacity(policy, 3L)
  expect_identical(policy$policy_schema, "v021_resource_policy_v2")
  expect_identical(derivation$projected_active_cmdstan_chains, 12L)
  expect_identical(derivation$projected_active_cmdstan_processes, 12L)
  expect_error(validate_v021_preflight_launch_capacity(policy, 4L), "ceiling")
})

test_that("thread environment is scoped and restored", {
  old <- Sys.getenv(v021_thread_variables, unset = NA_character_)
  inside <- with_v021_single_thread_environment(function() Sys.getenv(v021_thread_variables))
  expect_true(all(inside == "1"))
  expect_identical(Sys.getenv(v021_thread_variables, unset = NA_character_), old)
})

.preflight_snapshot <- function(chains = 12L, unknown = FALSE, unreadable = FALSE) {
  n <- chains + 2L
  data.frame(
    timestamp = rep("2026-08-02 UTC", n), pid = seq_len(n),
    ppid = c(0L, 1L, rep(2L, chains)),
    command = c("R", "R worker", rep("model method=sample", chains)),
    executable = c("/R", "/R", rep("/model", chains)),
    potential_cmdstan = c(FALSE, FALSE, rep(TRUE, chains)),
    readable = c(TRUE, !unreadable, rep(TRUE, chains)), stringsAsFactors = FALSE
  ) |> within({ if (unknown) command[n] <- "unknown cmdstan" })
}

test_that("monitor traverses descendants and fails closed", {
  policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
  ok <- monitor_v021_process_snapshots(list(.preflight_snapshot(12L)), policy, 1L,
                                       known_model_executables = "/model")
  expect_identical(ok$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(ok$observed_peak_active_cmdstan_processes, 12L)
  expect_silent(validate_v021_preflight_monitor(ok, policy))
  unreadable <- monitor_v021_process_snapshots(list(.preflight_snapshot(12L, unreadable = TRUE)),
                                                policy, 1L, "/model")
  expect_error(validate_v021_preflight_monitor(unreadable, policy), "not verified")
  unknown <- monitor_v021_process_snapshots(list(.preflight_snapshot(12L, unknown = TRUE)),
                                             policy, 1L, "/model")
  expect_error(validate_v021_preflight_monitor(unknown, policy), "not verified")
  exceeded <- monitor_v021_process_snapshots(list(.preflight_snapshot(13L)), policy, 1L, "/model")
  expect_error(validate_v021_preflight_monitor(exceeded, policy), "exceeded")

  diagnose <- data.frame(
    timestamp = "2026-08-02 UTC", pid = 999L, ppid = 2L,
    command = "bin/diagnose /output/chain.csv",
    executable = "/opt/cmdstan-2.36.0/bin/diagnose",
    potential_cmdstan = TRUE, readable = TRUE, stringsAsFactors = FALSE)
  overlap <- monitor_v021_process_snapshots(
    list(rbind(.preflight_snapshot(12L), diagnose)), policy, 1L, "/model")
  expect_identical(overlap$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(overlap$observed_peak_active_cmdstan_processes, 13L)
  expect_identical(overlap$compliance_status, "exceeded")
})

test_that("outer-worker registration captures stable procfs identity", {
  stat <- paste(200L, "(R)", "S", paste(c(100L, rep("0", 17L), "2000"),
                                         collapse = " "))
  registry <- register_v021_preflight_worker(
    200L, 100L, "direction_fit_worker", "batch-01", "direction-000001",
    proc_root = "/fixture", stat_reader = function(path) stat,
    clock = function() "2026-08-02T00:00:00Z")
  expect_identical(registry$pid, 200L)
  expect_identical(registry$start_time, "2000")
  expect_identical(registry$expected_ppid, 100L)
  expect_identical(registry$worker_role, "direction_fit_worker")
  expect_identical(registry$batch_id, "batch-01")
  expect_identical(registry$task_id, "direction-000001")
  expect_error(register_v021_preflight_worker(
    200L, 101L, "direction_fit_worker", "batch-01", "direction-000001",
    proc_root = "/fixture", stat_reader = function(path) stat), "PPID")
})

test_that("checkpoint restart preserves terminal states and immutable completion", {
  root <- tempfile("v021-preflight-"); dir.create(root)
  tasks <- select_v021_preflight_tasks(ten_species_direction_index(paste0("species_", 0:2), 20260802L))
  manifest <- build_v021_checkpoint_manifest(
    tasks[c("task_id", "direction_index", "target", "source", "seed")],
    4L, file.path(root, "checkpoints"))
  calls <- integer()
  run <- function(row) {
    calls <<- c(calls, row$direction_index[[1L]])
    if (row$direction_index[[1L]] == 2L) return(list(execution_state = "incomplete",
      interrupted = TRUE, interruption_reason = "controlled_preflight_interrupt"))
    if (row$direction_index[[1L]] == 3L) return(list(execution_state = "skipped"))
    if (row$direction_index[[1L]] == 4L) return(list(execution_state = "failed",
      failure_reason = "fixed_fixture_failure"))
    list(execution_state = "completed", result = list(value = 1))
  }
  first <- execute_v021_checkpoint_restart(manifest, run, file.path(root, "manifest.rds"))
  expect_identical(first$execution_state, c("completed", "incomplete", "skipped", "failed"))
  second <- execute_v021_checkpoint_restart(first, function(row) {
    calls <<- c(calls, row$direction_index[[1L]])
    list(execution_state = "completed", result = list(value = 2))
  }, file.path(root, "manifest.rds"))
  expect_identical(calls, c(1L, 2L, 3L, 4L, 2L))
  expect_identical(second$execution_state, c("completed", "completed", "skipped", "failed"))
  expect_error(write_v021_checkpoint_atomic(readRDS(second$output_location[[1L]])),
               "overwrite")
  broken <- second; broken$execution_state[[1L]] <- "completed"
  unlink(broken$output_location[[1L]])
  expect_error(plan_v021_checkpoint_restart(broken), "no checkpoint")
})

test_that("feature and summary contracts remain truth-free and explicit", {
  task <- data.frame(dataset_id="37", pair_id="pair-000001", task_id=1L,
    direction_index=1L, target="a", source="b", seed=7L, stringsAsFactors=FALSE)
  fit <- list(a_mean=-.2, a_median=-.19, a_sd=.05, a_q2.5=-.3, a_q97.5=-.1,
    positive_sign_probability=.01, negative_sign_probability=.99, p_sign2=.02, lfsr=.01,
    diag=list(worst_rhat=1.001, min_ess_bulk=700, min_ess_tail=500,
              n_divergent=0L, n_treedepth_hit=0L, ebfmi_min=.8),
    chain_sign_agreement=TRUE, diagnostic_class="converged",
    interaction_identifiable=TRUE, residual_identifiable=TRUE,
    residual_regime_disagreement=FALSE, n_pairs=20L)
  retained <- v021_preflight_retained_record(task, fit)
  feature <- build_v021_diagnostic_feature_record(retained,
    run_metadata=list(chains=4L, iter_sampling=2000L))
  expect_identical(feature$values$nominal_retained_draws, 8000L)
  expect_false(identical(feature$values$nominal_retained_draws, feature$values$bulk_ess))
  expect_identical(feature$missingness$aggregate_elpd, "not_executed")
  expect_false(feature$values$alr_cap_exposure_available)
  expect_false(any(grepl("truth|withhold|calibrat", names(feature$values), ignore.case=TRUE)))

  policy <- build_v021_resource_policy(proposed_outer_concurrency=3L)
  monitor <- monitor_v021_process_snapshots(list(.preflight_snapshot(12L)), policy, 1L, "/model")
  manifest_tasks <- rbind(
    task[c("task_id", "pair_id", "direction_index", "target", "source", "seed")],
    data.frame(task_id=1L, pair_id="pair-000001", direction_index=2L,
      target="b", source="a", seed=8L, stringsAsFactors=FALSE))
  manifest <- build_v021_checkpoint_manifest(manifest_tasks, 4L, tempfile("cp-"))
  manifest$execution_state <- "completed"; manifest$completion_timestamp <- "2026-08-02 UTC"
  selected_summary <- cbind(dataset_id = "37", manifest_tasks, stringsAsFactors = FALSE)
  summary <- build_v021_preflight_summary(
    build_v021_four_chain_preflight_config(tempdir()), selected_summary, policy, monitor, manifest,
    list(executable_path="/model", reuse_verified=TRUE,
         worker_side_compilation=FALSE, compilation_count=0L,
         worker_compilation_count=0L),
    list(controlled_interruption=TRUE, resumed_direction_index=2L,
         completed_tasks_not_rerun=TRUE, incomplete_task_resumed=TRUE,
         manifest_reconciled=TRUE), list(feature, feature), TRUE,
    list(passed=TRUE, active_preflight_cmdstan_processes=0L))
  expect_identical(summary$state, "passed")
  expect_identical(summary$summary_schema, "v021_four_chain_preflight_summary_v2")
  expect_identical(summary$resource_policy_schema, "v021_resource_policy_v2")
  expect_identical(summary$logical_host_threads, 16L)
  expect_identical(summary$reserved_host_threads, 4L)
  expect_identical(summary$usable_execution_capacity, 12L)
  expect_identical(summary$attempt_id, basename(tempdir()))
  expect_identical(summary$worker_compilation_count, 0L)
  expect_identical(summary$projected_active_chains, 12L)
  expect_identical(summary$projected_cmdstan_process_slots, 12L)
  expect_identical(summary$observed_peak_cmdstan_diagnostic_processes, 0L)
  expect_identical(summary$cmdstan_diagnostic_process_ids, integer())
  expect_identical(summary$cmdstan_diagnostic_executables, character())
  expect_identical(summary$vanished_process_ids, integer())
  expect_identical(summary$vanished_process_reasons, character())
  expect_identical(summary$zombie_process_ids, integer())
  expect_identical(summary$zombie_process_reasons, character())
  expect_identical(summary$transient_recapture_process_ids, integer())
  expect_identical(summary$transient_recapture_event_count, 0L)
  expect_identical(summary$transient_recapture_resolution_reasons, character())
  expect_true(summary$manifest_reconciled)
  expect_true(summary$incomplete_task_resumed)
  expect_true(summary$orphan_process_check$passed)
  expect_silent(validate_v021_preflight_summary(summary))
  failed <- summary; failed$state <- "failed"; failed$failure_reasons <- "ceiling"
  expect_silent(validate_v021_preflight_summary(failed))
})

test_that("runner wires three-fit batches and registered monitoring exactly once", {
  source_lines <- readLines(testthat::test_path(
    "../../benchmarks/mtist/run_v021_four_chain_preflight.R"), warn = FALSE)
  expect_true(any(grepl("run_batch\\(c\\(1L, 2L, 3L\\)\\)", source_lines)))
  expect_true(any(grepl("worker_registry = batch_worker_registry", source_lines,
                        fixed = TRUE)))
  expect_true(any(grepl("worker_registry.rds", source_lines, fixed = TRUE)))
  expect_equal(sum(grepl("with_v021_single_thread_environment\\(jobs\\[\\[i\\]\\]\\)",
                         source_lines)), 1L)
})
