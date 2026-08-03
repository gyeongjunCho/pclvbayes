source(testthat::test_path("../../benchmarks/mtist/mtist_adapter.R"))
source(testthat::test_path("../../benchmarks/mtist/ten_species_helpers.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_truth_isolation.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_resource_policy.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_checkpoint_manifest.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_diagnostic_features.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_four_chain_preflight.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_full_100_species.R"))

test_that("full benchmark configuration is exact and deterministic", {
  a <- build_v021_full_config(tempdir())
  b <- build_v021_full_config(tempdir())
  expect_identical(a, b)
  expect_silent(validate_v021_full_config(a))
  expect_identical(a$benchmark_schema, "v021_full_100_species_v1")
  expect_identical(a$chains, 4L)
  expect_identical(a$iter_warmup, 2000L)
  expect_identical(a$iter_sampling, 2000L)
  expect_identical(a$nominal_retained_draws, 8000L)
  expect_identical(a$maximum_simultaneous_fits, 3L)
  expect_false(a$run_kfold)
  expect_false(a$use_pathfinder)
})

test_that("truth-free metadata selection chooses the canonical 100-species study", {
  metadata <- data.frame(
    did = c(7L, 1L, 3L, 37L), n_species = c(100L, 100L, 100L, 10L),
    noise = c(.01, .01, .01, .01), n_timeseries = c(25L, 10L, 10L, 10L),
    n_timepoints = c(15L, 15L, 15L, 15L),
    sampling_scheme = c("even", "even", "random", "even"),
    ground_truth = letters[1:4], stringsAsFactors = FALSE)
  selected <- select_v021_100_species_dataset(metadata, function(did) TRUE)
  expect_identical(selected$did, 1L)
  expect_false(any(grepl("truth", names(selected), ignore.case = TRUE)))
  expect_error(select_v021_100_species_dataset(metadata, function(did) FALSE),
               "No eligible")
})

test_that("canonical task generation yields 4950 pairs and 9900 directions", {
  taxa <- paste0("species_", 0:99)
  first <- build_v021_full_task_table(taxa, 20260802L)
  second <- build_v021_full_task_table(taxa, 20260802L)
  expect_identical(first, second)
  expect_identical(nrow(first), 9900L)
  expect_identical(length(unique(first$task_id)), 4950L)
  expect_identical(first$direction_index, seq_len(9900L))
  expect_true(all(table(first$task_id) == 2L))
  expect_identical(first$target[1:2], c("species_0", "species_1"))
  expect_identical(first$source[1:2], c("species_1", "species_0"))
  expect_identical(first$seed[1:2], c(20362803L, 20362804L))
  expect_error(build_v021_full_task_table(taxa[-1L], 20260802L), "100")
})

test_that("policy v2 bounds the full benchmark at three fits", {
  policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
  derivation <- validate_v021_preflight_launch_capacity(policy, 3L)
  expect_identical(policy$policy_schema, "v021_resource_policy_v2")
  expect_identical(policy$logical_host_threads, 16L)
  expect_identical(policy$reserved_host_threads, 4L)
  expect_identical(policy$usable_chain_slots, 12L)
  expect_identical(derivation$projected_active_cmdstan_chains, 12L)
  expect_identical(derivation$projected_active_cmdstan_processes, 12L)
  expect_error(validate_v021_preflight_launch_capacity(policy, 4L), "ceiling")
})

test_that("restart resumes incomplete work without rerunning terminal tasks", {
  root <- tempfile("v021-full-"); dir.create(root)
  tasks <- build_v021_full_task_table(paste0("species_", 0:99), 20260802L)[1:4, ]
  manifest <- build_v021_checkpoint_manifest(
    tasks[c("task_id", "direction_index", "target", "source", "seed")],
    4L, file.path(root, "checkpoints"))
  completed <- .v021_record_from_row(manifest[1, , drop = FALSE], "completed",
    completion_timestamp = "2026-08-03 UTC", result = list(value = 1))
  incomplete <- .v021_record_from_row(manifest[2, , drop = FALSE], "incomplete",
    interrupted = TRUE, interruption_reason = "fixture", failure_reason = "fixture")
  skipped <- .v021_record_from_row(manifest[3, , drop = FALSE], "skipped")
  failed <- .v021_record_from_row(manifest[4, , drop = FALSE], "failed",
    failure_reason = "fixture_failure")
  records <- list(completed, incomplete, skipped, failed)
  for (i in seq_along(records)) {
    write_v021_checkpoint_atomic(records[[i]])
    manifest <- .v021_apply_checkpoint(manifest, i, records[[i]])
  }
  plan <- validate_v021_full_restart_states(manifest)
  expect_identical(plan$resume_indices, 2L)
  expect_identical(plan$manifest$execution_state,
                   c("completed", "incomplete", "skipped", "failed"))
  expect_error(write_v021_checkpoint_atomic(completed), "overwrite")
})

test_that("paths and nohup launch command are deterministic", {
  config <- build_v021_full_config(file.path(tempdir(), "full"))
  paths <- v021_full_paths(config, "20260803-010203")
  expect_identical(paths$checkpoint_root, file.path(config$output_root, "checkpoints"))
  expect_identical(paths$manifest, file.path(config$output_root, "manifest.rds"))
  expect_match(paths$log_file, "20260803-010203\\.log$")
  command <- build_v021_full_launch_command(
    "benchmarks/mtist/run_v021_full_100_species.R",
    "benchmarks/mtist/configs/v021_full_100_species.R",
    paths$log_file, paths$pid_file)
  expect_match(command, "nohup Rscript")
  expect_match(command, "PCLV_V021_FULL_CONFIG")
  expect_match(command, "< /dev/null & echo \\$!")
})

test_that("runner is truth-free, checkpointed, monitored, and compilation-free", {
  runner <- readLines(testthat::test_path(
    "../../benchmarks/mtist/run_v021_full_100_species.R"), warn = FALSE)
  expect_false(any(grepl("load_mtist_study|truth_matrix|same_nonzero_sign|opposite_nonzero_sign|absolute_zero",
                         runner)))
  expect_true(any(grepl("build_v021_inference_spec", runner, fixed = TRUE)))
  expect_true(any(grepl("write_v021_checkpoint_atomic", runner, fixed = TRUE)))
  expect_true(any(grepl("write_v021_manifest_atomic", runner, fixed = TRUE)))
  expect_true(any(grepl("monitor_v021_process_snapshots", runner, fixed = TRUE)))
  expect_true(any(grepl("worker_compilation_count = 0L", runner, fixed = TRUE)))
  expect_false(any(grepl("cmdstan_model|compile\\(", runner)))
})
