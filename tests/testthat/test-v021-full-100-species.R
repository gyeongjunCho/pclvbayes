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
  expect_identical(a$benchmark_schema, "v021_full_100_species_v2")
  expect_identical(a$chains, 4L)
  expect_identical(a$iter_warmup, 2000L)
  expect_identical(a$iter_sampling, 2000L)
  expect_identical(a$nominal_retained_draws, 8000L)
  expect_identical(a$maximum_simultaneous_fits, 3L)
  expect_true(a$run_kfold)
  expect_false(a$use_pathfinder)
  expect_identical(a$cpu_affinity, "0-15")
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

test_that("three-direction debug preflight permits one incomplete pair only explicitly", {
  tasks <- build_v021_full_task_table(paste0("species_", 0:99), 20260802L)[1:3, ]
  input <- tasks[c("task_id", "direction_index", "target", "source", "seed")]
  expect_error(build_v021_checkpoint_manifest(input, 4L, tempfile("checkpoints-")),
               "both directions")
  manifest <- build_v021_checkpoint_manifest(
    input, 4L, tempfile("checkpoints-"), allow_incomplete_pairs = TRUE)
  expect_identical(manifest$direction_index, 1:3)
  expect_identical(manifest$seed, c(20362803L, 20362804L, 20363803L))
  runner <- readLines(testthat::test_path(
    "../../benchmarks/mtist/run_v021_full_100_species.R"), warn = FALSE)
  expect_true(any(grepl("runnable_ordinals %in% tasks$direction_index",
                         runner, fixed = TRUE)))
})

test_that("policy v3 bounds the full benchmark at three fits", {
  policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
  derivation <- validate_v021_preflight_launch_capacity(policy, 3L)
  expect_identical(policy$policy_schema, "v021_resource_policy_v3")
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
  expect_true(any(grepl("acquire_v021_controller_ownership", runner, fixed = TRUE)))
  expect_true(any(grepl("reserve_v021_capacity", runner, fixed = TRUE)))
  expect_true(any(grepl("release_v021_capacity", runner, fixed = TRUE)))
  expect_true(any(grepl("write_v021_reservations_atomic", runner, fixed = TRUE)))
  expect_true(any(grepl("write_v021_manifest_atomic", runner, fixed = TRUE)))
  expect_true(any(grepl("initialize_v021_execution_root", runner, fixed = TRUE)))
  expect_true(any(grepl("write_v021_task_status_atomic", runner, fixed = TRUE)))
  expect_true(any(grepl("write_v021_attempt_ledger_atomic", runner, fixed = TRUE)))
  expect_true(any(grepl(".add_predictive_evaluation", runner, fixed = TRUE)))
  expect_true(any(grepl("monitor_v021_process_snapshots", runner, fixed = TRUE)))
  expect_true(any(grepl("worker_compilation_count = 0L", runner, fixed = TRUE)))
  expect_false(any(grepl("cmdstan_model|compile\\(", runner)))
  expect_false(any(grepl("stacking_results|pseudo_BMA|elpd_pointwise", runner)))
})

test_that("full storage projection is conservative and complete", {
  projection <- build_v021_full_storage_projection(800 * 1024^3)
  expected <- c("cmdstan_csv_four_chain", "profile_diagnostic_sidecars",
    "task_logs", "process_monitor_records", "manifest_status_attempt_ledger",
    "completion_posterior_diagnostic", "robustness_features",
    "direction_checkpoints", "subject_level_kfold", "retry_overhead")
  expect_identical(projection$categories$category, expected)
  expect_gt(projection$permanent_conservative_bytes, projection$permanent_central_bytes)
  expect_gt(projection$safety_adjusted_requirement_bytes,
            projection$permanent_conservative_bytes)
  expect_true(projection$sufficient)
})

test_that("full dry run uses authoritative state retries and ownership without sampling", {
  root <- tempfile("v021-full-dry-")
  config <- build_v021_full_config(root)
  prepared <- prepare_v021_full_execution(
    config, paste0("species_", 0:99),
    list(code_commit = "fixture", benchmark_schema = config$benchmark_schema),
    initialize = FALSE)
  audit <- run_v021_full_dry_run_audit(
    prepared, build_v021_resource_policy(), root, 2L)
  expect_identical(audit$task_count, 9900L)
  expect_identical(audit$pair_count, 4950L)
  expect_true(audit$first_attempt_failed)
  expect_true(audit$second_attempt_completed)
  expect_true(audit$completed_task_skipped)
  expect_identical(audit$retry_count, 2L)
  expect_identical(audit$maximum_wave_slots, 12L)
  expect_true(audit$ownership_released)
  expect_false(audit$sampling_launched)
  expect_true(file.exists(file.path(root, "canonical_manifest.rds")))
  expect_true(file.exists(file.path(root, "task_status.rds")))
  expect_true(file.exists(file.path(root, "attempt_ledger.rds")))
})

test_that("parent failures retain exact calls, phases, identities, and safe metadata", {
  root <- tempfile("v021-traces-")
  context <- new_v021_failure_trace_context(
    root, "parent", "batch-00042", 42L, "direction-000084")
  set_v021_failure_trace_phase(context, "artifact_finalization")
  set_v021_failure_trace_phase(context, "checkpoint_write", completed = TRUE)
  set_v021_failure_trace_phase(context, "feature_generation")

  trigger <- function() {
    offending_object <- 7L
    harmless_label <- "kept"
    truth_matrix <- matrix(1, 2, 2)
    posterior_draws <- matrix(2, 2, 2)
    complete_study <- list(secret = 1)
    hidden_environment <- new.env()
    hidden_closure <- function() NULL
    (offending_object)()
  }
  expect_error(with_v021_failure_tracing(trigger, context),
               "attempt to apply non-function")
  captured <- context$last_trace
  expect_true(file.exists(captured$trace_path))
  trace <- readRDS(captured$trace_path)
  expect_identical(trace$trace_schema, "v021_full_failure_trace_v1")
  expect_identical(trace$origin, "parent")
  expect_identical(trace$execution_phase, "feature_generation")
  expect_identical(trace$last_completed_phase, "checkpoint_write")
  expect_identical(trace$batch_id, "batch-00042")
  expect_identical(trace$task_ids, 42L)
  expect_identical(trace$direction_ids, "direction-000084")
  expect_true(any(grepl("offending_object", trace$calls, fixed = TRUE)))
  metadata <- unlist(lapply(trace$frame_metadata, function(x)
    vapply(x$objects, `[[`, character(1), "object_name")), use.names = FALSE)
  expect_true("offending_object" %in% metadata)
  expect_false(any(grepl("truth|study|draw|environment|closure", metadata,
                         ignore.case = TRUE)))
  expect_false(any(grepl("\\.tmp-", list.files(root, all.files = TRUE))))
})

test_that("outer-worker and monitor-child failures retain their own traces", {
  root <- tempfile("v021-child-traces-")
  run_failure <- function(origin) {
    context <- new_v021_failure_trace_context(
      root, origin, "batch-00007", 7L, "direction-000014")
    set_v021_failure_trace_phase(context,
      if (origin == "outer_worker") "worker_execution" else "monitor_polling")
    run_v021_traced_child(function() {
      offending_callback <- 1L
      (offending_callback)()
    }, context)
  }
  if (.Platform$OS.type == "windows") skip("forked worker tracing requires Unix")
  worker_job <- parallel::mcparallel(run_failure("outer_worker"), silent = TRUE)
  worker <- unname(parallel::mccollect(worker_job))[[1L]]
  monitor <- run_failure("monitor_child")
  for (result in list(worker, monitor)) {
    expect_s3_class(result, "v021_traced_child_error")
    expect_match(result$condition_message, "attempt to apply non-function")
    expect_true(file.exists(result$trace$trace_path))
    persisted <- readRDS(result$trace$trace_path)
    expect_true(any(grepl("offending_callback", persisted$calls, fixed = TRUE)))
    expect_identical(persisted$task_ids, 7L)
    expect_identical(persisted$direction_ids, "direction-000014")
  }
  expect_identical(readRDS(worker$trace$trace_path)$origin, "outer_worker")
  expect_identical(readRDS(monitor$trace$trace_path)$origin, "monitor_child")
})

test_that("trace-write failure preserves the original condition", {
  context <- new_v021_failure_trace_context(
    tempfile("v021-trace-failure-"), "parent", "batch-1", 1L, "direction-1")
  bad_writer <- function(object, path) stop("injected trace writer failure")
  result <- tryCatch(with_v021_failure_tracing(function() {
    offending_value <- 2L
    (offending_value)()
  }, context, writer = bad_writer), error = identity)
  expect_s3_class(result, "error")
  expect_match(conditionMessage(result), "attempt to apply non-function")
  expect_match(context$last_trace$trace_persistence_error,
               "injected trace writer failure")
  expect_true(is.na(context$last_trace$trace_path))
  expect_true(any(grepl("offending_value", context$last_trace$trace$calls,
                         fixed = TRUE)))
})

test_that("runner marks and traces every V021-06 execution boundary", {
  runner <- readLines(testthat::test_path(
    "../../benchmarks/mtist/run_v021_full_100_species.R"), warn = FALSE)
  phases <- c("batch_launch", "worker_execution", "fit_return_handling",
    "scientific_state_classification", "artifact_finalization",
    "feature_generation", "checkpoint_write", "manifest_update",
    "summary_failure_persistence", "monitor_polling", "monitor_shutdown")
  for (phase in phases) expect_true(any(grepl(phase, runner, fixed = TRUE)))
  expect_true(any(grepl("run_v021_traced_child", runner, fixed = TRUE)))
  expect_true(any(grepl("persist_v021_failure_payload", runner, fixed = TRUE)))
})

.v021_cleanup_fixture <- function(root) {
  registry <- rbind(
    build_v021_worker_registry(200L, "20", 100L, "direction_fit_worker",
                               "batch-1", "direction-1", "2026-08-03 UTC"),
    build_v021_worker_registry(201L, "21", 100L, "resource_monitor",
                               "batch-1", "batch-1-monitor", "2026-08-03 UTC")
  )
  ownership <- new_v021_batch_ownership(
    "batch-1", 100L, registry, "/known/model", root,
    file.path(root, "cleanup-audit.rds"))
  records <- data.frame(
    pid = c(100L, 200L, 201L, 300L, 301L, 900L),
    ppid = c(1L, 100L, 100L, 200L, 200L, 1L),
    start_time = c("10", "20", "21", "30", "31", "90"),
    classification = c("parent_r", "registered_outer_worker",
      "registered_outer_worker", "cmdstan_chain", "cmdstan_diagnostic", "unrelated"),
    executable = c("/usr/bin/R", "/usr/bin/R", "/usr/bin/R", "/known/model",
                   "/cmdstan/bin/diagnose", "/unrelated"),
    stringsAsFactors = FALSE)
  list(ownership = ownership, records = records)
}

test_that("batch cleanup signals only proven owned descendants and escalates boundedly", {
  root <- tempfile("v021-cleanup-"); dir.create(root)
  fixture <- .v021_cleanup_fixture(root)
  live <- c(200L, 201L, 300L, 301L, 900L)
  signals <- list()
  signaler <- function(pids, signal) {
    signals[[length(signals) + 1L]] <<- list(pids = pids, signal = signal)
    if (signal == "SIGINT") live <<- setdiff(live, c(200L, 201L, 301L))
    if (signal == "SIGTERM") live <<- setdiff(live, pids)
  }
  alive <- function(pids, start_times) intersect(pids, live)
  audit <- cleanup_v021_owned_batch(
    fixture$ownership, simpleError("parent failure"), "parent", "fit_return_handling",
    function() fixture$records, signaler, alive, reaper = function() list(reaped = TRUE),
    sleeper = function(x) NULL, wait_attempts = 2L)
  expect_identical(sort(audit$owned_processes$pid), c(200L, 201L, 300L, 301L))
  expect_true(900L %in% audit$excluded_processes$pid)
  expect_identical(signals[[1L]]$signal, "SIGINT")
  expect_false(900L %in% signals[[1L]]$pids)
  expect_identical(signals[[2L]], list(pids = 300L, signal = "SIGTERM"))
  expect_length(audit$survivors, 0L)
  expect_false(audit$sigkill_used)
  expect_true(file.exists(fixture$ownership$cleanup_audit_path))
  expect_identical(readRDS(fixture$ownership$cleanup_audit_path)$cleanup_schema,
                   "v021_full_cleanup_audit_v1")
})

test_that("cleanup is idempotent and preserves original errors when audit writing fails", {
  root <- tempfile("v021-cleanup-idempotent-"); dir.create(root)
  fixture <- .v021_cleanup_fixture(root)
  writes <- 0L
  writer <- function(object, path) { writes <<- writes + 1L; stop("audit write failed") }
  live <- function(pids, start_times) integer()
  original <- simpleError("authoritative parent error")
  first <- cleanup_v021_owned_batch(
    fixture$ownership, original, "parent", "checkpoint_write",
    function() fixture$records, function(...) NULL, live,
    writer = writer, sleeper = function(x) NULL)
  second <- cleanup_v021_owned_batch(
    fixture$ownership, simpleError("replacement"), "parent", "manifest_update",
    function() stop("must not rerun"), function(...) stop("must not signal"), live,
    writer = writer)
  expect_identical(first, second)
  expect_identical(writes, 1L)
  expect_identical(first$trigger_condition, "authoritative parent error")
  expect_match(first$audit_persistence_error, "audit write failed")
})

test_that("normal completed cleanup is a no-op when owned children have exited", {
  root <- tempfile("v021-cleanup-success-"); dir.create(root)
  fixture <- .v021_cleanup_fixture(root)
  signal_count <- 0L
  audit <- cleanup_v021_owned_batch(
    fixture$ownership, NULL, "parent", "manifest_update",
    function() fixture$records,
    function(...) signal_count <<- signal_count + 1L,
    function(pids, start_times) integer(),
    reaper = function() list(reaped = TRUE), sleeper = function(x) NULL)
  expect_identical(signal_count, 0L)
  expect_length(audit$signals, 0L)
  expect_length(audit$survivors, 0L)
  expect_identical(audit$trigger_condition, "normal_return")
})

test_that("cleanup guard covers collection, monitor, storage, and persistence phases", {
  phases <- c("batch_launch", "fit_return_handling", "monitor_shutdown",
              "scientific_state_classification", "artifact_finalization",
              "feature_generation", "checkpoint_write", "manifest_update")
  for (phase in phases) {
    root <- tempfile("v021-phase-cleanup-"); dir.create(root)
    fixture <- .v021_cleanup_fixture(root)
    audit <- cleanup_v021_owned_batch(
      fixture$ownership, simpleError(paste("failure", phase)), "parent", phase,
      function() fixture$records, function(...) NULL,
      function(pids, start_times) integer(), sleeper = function(x) NULL)
    expect_identical(audit$phase, phase)
    expect_match(audit$trigger_condition, phase, fixed = TRUE)
  }
})

test_that("ordinary errors and interrupts retain traces without changing the condition", {
  root <- tempfile("v021-condition-traces-")
  error_context <- new_v021_failure_trace_context(root, "parent")
  error <- tryCatch(with_v021_failure_tracing(function() stop("ordinary error"),
                                               error_context), error = identity)
  expect_identical(conditionMessage(error), "ordinary error")
  expect_true(file.exists(error_context$last_trace$trace_path))

  interrupt_context <- new_v021_failure_trace_context(root, "parent")
  interrupt <- structure(list(message = "controlled interrupt", call = quote(worker_wait())),
                         class = c("interrupt", "condition"))
  observed <- tryCatch(with_v021_failure_tracing(function() stop(interrupt),
                                                  interrupt_context),
                       interrupt = identity, error = identity)
  expect_identical(conditionMessage(observed), "controlled interrupt")
  expect_true(file.exists(interrupt_context$last_trace$trace_path))
  expect_true("interrupt" %in%
    readRDS(interrupt_context$last_trace$trace_path)$condition_class)
})

test_that("a finished CSV cannot promote an unfinished direction", {
  root <- tempfile("v021-csv-no-promotion-"); dir.create(root)
  tasks <- data.frame(task_id = c(1L, 1L), direction_index = 1:2,
    target = c("a", "b"), source = c("b", "a"), seed = c(11L, 12L))
  manifest <- build_v021_checkpoint_manifest(tasks, 4L, file.path(root, "checkpoints"))
  running <- .v021_record_from_row(manifest[1, , drop = FALSE], "running")
  write_v021_checkpoint_atomic(running)
  manifest <- .v021_apply_checkpoint(manifest, 1L, running)
  csv <- file.path(root, "completed-looking.csv")
  writeLines(c("# num_samples = 2000", "lp__,a", "-1,0.2"), csv)
  plan <- plan_v021_checkpoint_restart(manifest)
  expect_true(file.exists(csv))
  expect_identical(plan$manifest$execution_state[[1L]], "incomplete")
  expect_false(plan$manifest$execution_state[[1L]] == "completed")
  expect_null(readRDS(plan$manifest$output_location[[1L]])$result)
})

test_that("cleanup never changes completed checkpoint content", {
  root <- tempfile("v021-cleanup-completed-"); dir.create(root)
  tasks <- data.frame(task_id = c(1L, 1L), direction_index = 1:2,
    target = c("a", "b"), source = c("b", "a"), seed = c(11L, 12L))
  manifest <- build_v021_checkpoint_manifest(tasks, 4L, file.path(root, "checkpoints"))
  completed <- .v021_record_from_row(manifest[1, , drop = FALSE], "completed",
    completion_timestamp = "2026-08-03 UTC", result = list(retained = TRUE))
  write_v021_checkpoint_atomic(completed)
  before <- unname(tools::md5sum(completed$output_location))
  fixture <- .v021_cleanup_fixture(root)
  cleanup_v021_owned_batch(fixture$ownership, simpleError("fixture"), "parent",
    "checkpoint_write", function() fixture$records, function(...) NULL,
    function(pids, start_times) integer(), sleeper = function(x) NULL)
  expect_identical(unname(tools::md5sum(completed$output_location)), before)
})

test_that("runner installs ownership cleanup before releasing workers", {
  runner <- readLines(testthat::test_path(
    "../../benchmarks/mtist/run_v021_full_100_species.R"), warn = FALSE)
  ownership_line <- grep("new_v021_batch_ownership", runner, fixed = TRUE)[[1L]]
  guard_line <- grep("on.exit({", runner, fixed = TRUE)[[1L]]
  release_line <- grep("file.create(start_file)", runner, fixed = TRUE)[[1L]]
  collect_line <- grep("collect_v021_tracked_jobs(fit_job_tracker)",
                       runner, fixed = TRUE)[[1L]]
  expect_lt(ownership_line, guard_line)
  expect_lt(guard_line, release_line)
  expect_lt(release_line, collect_line)
  expect_true(any(grepl("SIGKILL and uncatchable parent termination", runner,
                         fixed = TRUE)) ||
              any(grepl("SIGKILL and uncatchable parent termination",
                        readLines(testthat::test_path(
                          "../../benchmarks/mtist/v021_full_100_species.R")),
                        fixed = TRUE)))
})

test_that("intermediate sampler health never claims NA convergence passed", {
  helper <- paste(readLines(testthat::test_path("../../R/pclv_helpers.R"),
                            warn = FALSE), collapse = "\n")
  expect_match(helper, "sampler health ok; convergence summary pending", fixed = TRUE)
  expect_false(grepl("diag ok", helper, fixed = TRUE))
  expect_false(grepl("convergence ok", helper, fixed = TRUE))
})

test_that("outer worker preserves truth-free predictive context until controller evaluation", {
  skip_if_not_installed("phyloseq")
  physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(matrix(c(1, 2, 3, 4), 2, 2,
      dimnames = list(c("a", "b"), c("s1", "s2"))), taxa_are_rows = TRUE),
    phyloseq::sample_data(data.frame(subject = c("x", "x"), time = c(0, 1),
      row.names = c("s1", "s2"))))
  task <- data.frame(dataset_id = "dataset-1", task_id = 1L,
    direction_index = 1L, target = "a", source = "b", seed = 19L)
  fixture <- list(task = task, spec = build_v021_inference_spec(
    list(physeq = physeq, taxa = c("a", "b")), list(chains = 1L), task))
  ctx <- setNames(rep(list(1), length(v021_runtime_context_fields)),
                  v021_runtime_context_fields)
  ctx$meta_df <- data.frame(Sample = c("s1", "s2"), subject = c("x", "x"),
                            time = c(0, 1))
  ctx$sm_mat <- matrix(1:4, 2, 2,
    dimnames = list(c("a", "b"), c("s1", "s2")))
  ctx$mod_exe_file <- "/approved/canonical-model"
  runtime <- build_v021_runtime_context(ctx)
  context <- list(
    pair_in = data.frame(subject = "s1", time = 1), split_seed = 20260802L,
    sampling_seed = fixture$task$seed[[1L]], pair_tag = "fixture")
  fit_stub <- function(target, partner, ctx, seed_override, progress_local) {
    predictive <- list(
      pair_in = data.frame(subject = "s1", time = 1), split_seed = 20260802L,
      sampling_seed = seed_override, pair_tag = "fixture")
    list(target = target, partner = partner, seed = seed_override,
         diagnostic_failure = list(NULL), .predictive_context = predictive)
  }
  environment(fit_stub) <- baseenv()
  job <- make_v021_confirmation_fit_closure(fixture$spec, runtime, fit_stub)
  worker_result <- job()
  expect_identical(worker_result$.predictive_context, context)
  expect_silent(assert_v021_truth_free_schema(worker_result$.predictive_context,
                                              "worker predictive context"))
  evaluator <- function(x) {
    expect_identical(x$.predictive_context, context)
    x$.predictive_context <- NULL
    x$kfold_success_total <- 0L
    x
  }
  policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
  writes <- list()
  evaluated <- with_v021_predictive_reservation(
    worker_result, new_v021_reservations(), policy, "direction-000001", 1L,
    evaluator, persist = function(x) writes[[length(writes) + 1L]] <<- x)
  expect_false(".predictive_context" %in% names(evaluated$result))
  expect_identical(evaluated$result$kfold_success_total, 0L)
  expect_identical(vapply(writes, function(x) tail(x$state, 1L), character(1)),
                   c("reserved", "worker_terminated"))
})

test_that("predictive reservations release once on success and failure", {
  policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
  run <- function(evaluator) {
    writes <- list()
    value <- tryCatch(with_v021_predictive_reservation(
      list(.predictive_context = list()), new_v021_reservations(), policy,
      "direction-000001", 1L, evaluator,
      persist = function(x) writes[[length(writes) + 1L]] <<- x), error = identity)
    list(value = value, writes = writes)
  }
  ok <- run(function(x) x)
  expect_false(inherits(ok$value, "error"))
  expect_identical(length(ok$writes), 2L)
  expect_identical(tail(ok$writes[[2L]]$state, 1L), "worker_terminated")
  bad <- run(function(x) stop("predictive fixture failure"))
  expect_match(conditionMessage(bad$value), "predictive fixture failure")
  expect_identical(length(bad$writes), 2L)
  expect_identical(tail(bad$writes[[2L]]$state, 1L), "worker_terminated")
})

test_that("collected worker handles are not reaped twice after controller failure", {
  jobs <- lapply(1:3, function(pid) structure(list(pid = pid), class = "process"))
  fit_tracker <- new_v021_child_job_tracker(jobs)
  monitor_tracker <- new_v021_child_job_tracker()
  calls <- list()
  collector <- function(children, wait) {
    calls[[length(calls) + 1L]] <<- vapply(children, `[[`, integer(1), "pid")
    setNames(as.list(rep(TRUE, length(children))), calls[[length(calls)]])
  }
  expect_length(collect_v021_tracked_jobs(fit_tracker, collector), 3L)
  expect_identical(reap_v021_uncollected_jobs(
    fit_tracker, monitor_tracker, collector = collector)$class, "none")
  expect_identical(length(calls), 1L)
})

test_that("terminal controller children receive a final blocking reap", {
  calls <- list()
  collector <- function(wait = TRUE) {
    calls[[length(calls) + 1L]] <<- list(wait = wait)
    list()
  }
  expect_identical(reap_v021_terminal_children(collector), list())
  expect_identical(calls, list(list(wait = TRUE)))
})

test_that("empty durable monitor history starts with zero recovered peaks", {
  root <- tempfile("v021-monitor-peaks-")
  dir.create(file.path(root, "monitor"), recursive = TRUE)
  peaks <- recover_v021_monitor_peaks(
    root, build_v021_resource_policy(), "/models/pclv")
  expect_identical(peaks, list(chains = 0L, processes = 0L))
})

test_that("collected monitor exit waits for transient procfs zombie removal", {
  checks <- 0L
  path_exists <- function(path) {
    checks <<- checks + 1L
    checks < 3L
  }
  sleeps <- numeric()
  removed <- integer()
  expect_true(wait_v021_collected_child_exit(
    123L, remover = function(pid) removed <<- c(removed, pid),
    path_exists = path_exists,
    sleep = function(seconds) sleeps <<- c(sleeps, seconds)))
  expect_identical(removed, c(123L, 123L))
  expect_identical(sleeps, c(0.01, 0.01))
  expect_error(wait_v021_collected_child_exit(NA_integer_), "positive integer")
})

test_that("dead-controller restart reconciles stale reservations without completion", {
  policy <- build_v021_resource_policy(proposed_outer_concurrency = 3L)
  reservations <- reserve_v021_capacity(
    new_v021_reservations(), policy, "direction-000001", 1L, "kfold_fit")
  expect_error(reconcile_v021_dead_controller_reservations(
    reservations, policy, TRUE), "Live-controller")
  reconciled <- reconcile_v021_dead_controller_reservations(
    reservations, policy, FALSE)
  expect_identical(reconciled$state, "worker_terminated")
  expect_false(any(reconciled$state == "reserved"))
})

test_that("reservation attempt identity follows authoritative V021-03 status", {
  manifest <- build_v021_execution_manifest(
    dataset_id = "d", taxa_order = c("a", "b"),
    task_table = data.frame(task_id = c(1L, 1L), direction_index = 1:2,
      target = c("a", "b"), source = c("b", "a"), seed = c(11L, 12L)),
    public_seed = 1L,
    kfold_seed = 2L, preprocessing_config = list(id = "p"),
    posterior_config = list(id = "m"), kfold_config = list(id = "k"),
    predictive_config = list(id = "q"), provenance = list(commit = "x"))
  status <- new_v021_task_status(manifest)
  status$tasks$state <- "failed"
  status$tasks$attempt_count <- 1L
  expect_identical(v021_task_attempt_number(status, 1L, TRUE), 2L)
  status$tasks$state <- "running"
  status$tasks$attempt_count <- 2L
  expect_identical(v021_task_attempt_number(status, 1L), 2L)
  runner <- paste(readLines(testthat::test_path(
    "../../benchmarks/mtist/run_v021_full_100_species.R"), warn = FALSE),
    collapse = "\n")
  expect_false(grepl("manifest$retry_count[[i]] + 1L", runner, fixed = TRUE))
})

test_that("recording handlers rethrow once while trace, cleanup, and failure payload persist", {
  root <- tempfile("v021-resignal-"); dir.create(root)
  trace_root <- file.path(root, "traces")
  context <- new_v021_failure_trace_context(trace_root, "parent",
    "batch-9", 9L, "direction-18")
  trace_writes <- 0L
  trace_writer <- function(object, path) {
    trace_writes <<- trace_writes + 1L
    .v021_atomic_save_rds(object, path)
  }
  batch_records <- 0L
  execution_records <- 0L
  cleanup_runs <- 0L
  batch_recorder <- function(condition) batch_records <<- batch_records + 1L
  execution_recorder <- function(condition) execution_records <<- execution_records + 1L
  parent_handler <- function(condition) {
    if (is.null(context$last_trace))
      capture_v021_failure_trace(condition, context, writer = trace_writer)
  }
  failing_batch <- function() {
    on.exit(cleanup_runs <<- cleanup_runs + 1L, add = TRUE)
    offending_value <- 1L
    (offending_value)()
  }

  original <- tryCatch(
    withCallingHandlers(
      tryCatch(
        tryCatch(failing_batch(), error = function(condition)
          record_and_resignal_v021_condition(condition, batch_recorder)),
        error = function(condition)
          record_and_resignal_v021_condition(condition, execution_recorder)),
      error = parent_handler),
    error = identity)
  failure_path <- file.path(root, "failure_payload.rds")
  persisted <- persist_v021_failure_payload(list(
    failure_schema = "v021_full_100_species_failure_v1",
    reason = conditionMessage(original), parent_trace = context$last_trace),
    failure_path, trace_root)

  expect_match(conditionMessage(original), "attempt to apply non-function")
  expect_identical(batch_records, 1L)
  expect_identical(execution_records, 1L)
  expect_identical(cleanup_runs, 1L)
  expect_identical(trace_writes, 1L)
  expect_true(file.exists(context$last_trace$trace_path))
  expect_true(isTRUE(persisted$persisted))
  expect_true(file.exists(failure_path))
  expect_identical(readRDS(failure_path)$reason, conditionMessage(original))
  expect_identical(length(list.files(trace_root,
    pattern = "^v021_full_failure_trace_v1-.*[.]rds$")), 1L)
})
