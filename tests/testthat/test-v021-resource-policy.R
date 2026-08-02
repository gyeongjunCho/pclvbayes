source(testthat::test_path("../../benchmarks/mtist/v021_resource_policy.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_checkpoint_manifest.R"))

make_process_snapshot <- function(chain_count = 0L, timestamp = "2026-08-02T00:00:00Z",
                                  complete = TRUE, indirect = TRUE,
                                  unknown_potential = FALSE, pathfinder = FALSE) {
  root <- data.frame(
    timestamp = timestamp, pid = 100L, ppid = 1L,
    command = "Rscript preflight.R", executable = "/usr/bin/R",
    potential_cmdstan = FALSE, readable = TRUE, stringsAsFactors = FALSE
  )
  worker <- data.frame(
    timestamp = timestamp, pid = 200L, ppid = 100L,
    command = "R --worker", executable = "/usr/lib/R/bin/exec/R",
    potential_cmdstan = FALSE, readable = TRUE, stringsAsFactors = FALSE
  )
  rows <- list(root, worker)
  if (chain_count > 0L) {
    chains <- data.frame(
      timestamp = rep(timestamp, chain_count), pid = 300L + seq_len(chain_count),
      ppid = if (indirect) rep(200L, chain_count) else rep(100L, chain_count),
      command = paste0("/models/pclv id=", seq_len(chain_count), " method=sample"),
      executable = rep("/models/pclv", chain_count), potential_cmdstan = TRUE,
      readable = rep(TRUE, chain_count), stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- chains
  }
  if (unknown_potential) {
    rows[[length(rows) + 1L]] <- data.frame(
      timestamp = timestamp, pid = 999L, ppid = 200L,
      command = "/models/pclv method=mystery", executable = "/models/pclv",
      potential_cmdstan = TRUE, readable = TRUE, stringsAsFactors = FALSE
    )
  }
  if (pathfinder) {
    rows[[length(rows) + 1L]] <- data.frame(
      timestamp = timestamp, pid = 998L, ppid = 200L,
      command = "/models/pclv method=pathfinder", executable = "/models/pclv",
      potential_cmdstan = TRUE, readable = TRUE, stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  if (!complete) out$readable[out$pid == 200L] <- FALSE
  rownames(out) <- NULL
  out
}

primary_operation_spec <- function(workers = 2L) build_v021_operation_spec(
  c("pathfinder", "main_fit", "retry_fit", "kfold_fit"), workers)

raw_cmdline <- function(argv, terminal_nul = TRUE) {
  bytes <- lapply(argv, function(x) c(charToRaw(x), as.raw(0L)))
  out <- do.call(c, bytes)
  if (!terminal_nul) out <- head(out, -1L)
  out
}

snapshot_with_argv <- function(argv, executable = "/models/pclv", pid = 301L,
                               ppid = 100L, potential_cmdstan = TRUE,
                               readable = TRUE) {
  root <- data.frame(
    timestamp = "2026-08-02T00:00:00Z", pid = 100L, ppid = 1L,
    command = "Rscript preflight.R", argv = I(list(c("Rscript", "preflight.R"))),
    executable = "/usr/bin/R", potential_cmdstan = FALSE, readable = TRUE)
  child <- data.frame(
    timestamp = "2026-08-02T00:00:00Z", pid = pid, ppid = ppid,
    command = paste(argv, collapse = " "), argv = I(list(argv)),
    executable = executable, potential_cmdstan = potential_cmdstan,
    readable = readable)
  rbind(root, child)
}

test_that("Linux NUL-delimited cmdlines preserve exact argv boundaries", {
  argv <- c("/models/pclv", "method=sample", "num_samples=2000",
            "argument containing spaces", "", "id=1")
  decoded <- .v021_decode_linux_cmdline(raw_cmdline(argv))
  expect_identical(decoded$argv, argv)
  expect_identical(decoded$argv[[1L]], "/models/pclv")
  expect_identical(decoded$argv[[2L]], "method=sample")
  expect_match(decoded$command, "/models/pclv method=sample num_samples=2000")
  expect_match(decoded$command, "argument containing spaces  id=1", fixed = TRUE)
  expect_length(decoded$argv, length(argv))

  old <- paste(rawToChar(raw_cmdline(argv)[as.integer(raw_cmdline(argv)) != 0L],
                         multiple = TRUE), collapse = "")
  expect_match(old, "/models/pclvmethod=samplenum_samples=2000", fixed = TRUE)
  expect_false(grepl("(^|[[:space:]])method=sample([[:space:]]|$)", old))
})

test_that("Linux cmdline decoding fails closed for malformed or unreadable bytes", {
  expect_error(.v021_decode_linux_cmdline(character()), "raw vector")
  expect_error(.v021_decode_linux_cmdline(raw_cmdline(c("/models/pclv", "method=sample"),
                                                       terminal_nul = FALSE)), "NUL")
  expect_error(.v021_decode_linux_cmdline(raw_cmdline(c("", "method=sample"))),
               "executable")
})

test_that("bounded cmdline reads ignore reported size and reject unsafe input", {
  path <- tempfile("v021-cmdline-")
  on.exit(unlink(path), add = TRUE)
  bytes <- raw_cmdline(c("/models/pclv", "method=sample", "argument with spaces"))
  con <- file(path, open = "wb")
  writeBin(bytes, con)
  close(con)
  observed <- .v021_read_linux_cmdline(path)
  expect_identical(observed, bytes)
  expect_identical(.v021_decode_linux_cmdline(observed)$argv,
                   c("/models/pclv", "method=sample", "argument with spaces"))
  expect_identical(.v021_linux_cmdline_max_bytes, 1048576L)

  empty <- tempfile("v021-empty-cmdline-")
  on.exit(unlink(empty), add = TRUE)
  file.create(empty)
  expect_error(.v021_read_linux_cmdline(empty), "empty")
  expect_error(.v021_read_linux_cmdline(paste0(path, "-missing")), "unreadable")
  expect_error(.v021_read_linux_cmdline(path, max_bytes = 8L), "truncated")

  unterminated <- tempfile("v021-unterminated-cmdline-")
  on.exit(unlink(unterminated), add = TRUE)
  con <- file(unterminated, open = "wb")
  writeBin(raw_cmdline(c("/models/pclv", "method=sample"), terminal_nul = FALSE), con)
  close(con)
  expect_error(.v021_decode_linux_cmdline(.v021_read_linux_cmdline(unterminated)), "NUL")
})

test_that("live Linux procfs cmdline is read independently of pseudo-file size", {
  skip_if_not(Sys.info()[["sysname"]] == "Linux")
  path <- "/proc/self/cmdline"
  observed <- .v021_read_linux_cmdline(path)
  decoded <- .v021_decode_linux_cmdline(observed)
  expect_gt(length(observed), 0L)
  expect_gt(length(decoded$argv), 0L)
  expect_true(nzchar(decoded$argv[[1L]]))
  expect_identical(tail(observed, 1L), as.raw(0L))
})

test_that("argv classification recognizes known operations and rejects impostors", {
  sample <- snapshot_with_argv(c("/models/pclv", "method=sample", "id=1"))
  classified <- classify_v021_process_snapshot(sample, 100L, "/models/pclv")
  expect_identical(classified$classification[[2L]], "cmdstan_chain")

  pathfinder <- snapshot_with_argv(c("/models/pclv", "method=pathfinder", "num_paths=8"))
  classified <- classify_v021_process_snapshot(pathfinder, 100L, "/models/pclv")
  expect_identical(classified$classification[[2L]], "pathfinder_process")

  for (operation in c("retry", "kfold")) {
    command <- snapshot_with_argv(c("/models/pclv", "method=sample", paste0("operation=", operation)))
    classified <- classify_v021_process_snapshot(command, 100L, "/models/pclv")
    expect_identical(classified$classification[[2L]], "cmdstan_chain")
  }

  embedded <- snapshot_with_argv(c("/models/pclv", "note=method=sample in text"))
  classified <- classify_v021_process_snapshot(embedded, 100L, "/models/pclv")
  expect_identical(classified$classification[[2L]], "unknown_potential_cmdstan")

  unknown <- snapshot_with_argv(c("/models/pclv", "method=mystery"))
  classified <- classify_v021_process_snapshot(unknown, 100L, "/models/pclv")
  expect_identical(classified$classification[[2L]], "unknown_potential_cmdstan")

  unrelated <- snapshot_with_argv(c("/tmp/not-stan", "method=sample"),
                                  executable = "/tmp/not-stan",
                                  potential_cmdstan = FALSE)
  classified <- classify_v021_process_snapshot(unrelated, 100L, "/models/pclv")
  expect_identical(classified$classification[[2L]], "other_descendant")
})

test_that("an attempt3-shaped eight-chain snapshot is verified without sampling", {
  root <- snapshot_with_argv(c("R", "--worker"), executable = "/usr/lib/R/bin/exec/R",
                             pid = 200L, potential_cmdstan = FALSE)
  chains <- lapply(seq_len(8L), function(i) data.frame(
    timestamp = "2026-08-02T00:00:00Z", pid = 300L + i,
    ppid = if (i <= 4L) 201L else 202L,
    command = paste("/models/pclv", "method=sample", paste0("id=", i)),
    argv = I(list(c("/models/pclv", "method=sample", paste0("id=", i)))),
    executable = "/models/pclv", potential_cmdstan = TRUE, readable = TRUE))
  workers <- data.frame(
    timestamp = rep("2026-08-02T00:00:00Z", 2L), pid = 201:202,
    ppid = rep(100L, 2L), command = rep("R --worker", 2L),
    argv = I(rep(list(c("R", "--worker")), 2L)),
    executable = rep("/usr/lib/R/bin/exec/R", 2L),
    potential_cmdstan = FALSE, readable = TRUE)
  snapshot <- rbind(root[1L, ], workers, do.call(rbind, chains))
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 8L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 8L)
  expect_identical(monitor$compliance_status, "compliant")
})

test_that("the policy records the exact 12-thread, 2-reserved, 10-chain contract", {
  policy <- build_v021_resource_policy()
  expect_identical(policy$policy_schema, "v021_resource_policy_v1")
  expect_identical(policy$logical_host_threads, 12L)
  expect_identical(policy$reserved_host_threads, 2L)
  expect_identical(policy$usable_chain_slots, 10L)
  expect_identical(policy$maximum_active_cmdstan_chains, 10L)
  expect_identical(policy$cpu_threads_per_active_chain, 1L)
  expect_identical(policy$numerical_library_threads, 1L)
  expect_identical(policy$main_chains, 4L)
  expect_identical(policy$proposed_outer_concurrency, 2L)
  expect_false(policy$global_slot_scheduler)
  expect_silent(validate_v021_resource_policy(policy))
  expect_identical(policy, build_v021_resource_policy())
})

test_that("malformed and contradictory resource policies are rejected", {
  expect_error(build_v021_resource_policy(logical_host_threads = 11L), "Usable")
  expect_error(build_v021_resource_policy(reserved_host_threads = 3L), "Usable")
  expect_error(build_v021_resource_policy(cpu_threads_per_active_chain = 2L), "Exactly one")
  expect_error(build_v021_resource_policy(kfold_chains = 2L,
                                          kfold_parallel_chains = 3L), "cannot exceed")
  expect_error(build_v021_resource_policy(global_slot_scheduler = TRUE), "No verified")
  expect_error(build_v021_resource_policy(proposed_outer_concurrency = 3L), "exceeds")
  policy <- build_v021_resource_policy()
  policy$usable_chain_slots <- 9L
  expect_error(validate_v021_resource_policy(policy), "12/2/10")
  policy <- build_v021_resource_policy()
  policy$operation_slots$simultaneous_chain_slots[[1L]] <- 3L
  expect_error(validate_v021_resource_policy(policy), "contradicts")
  policy <- build_v021_resource_policy()
  policy$pathfinder_processes <- 0L
  policy$operation_slots$cmdstan_process_slots[[3L]] <- 0L
  expect_error(validate_v021_resource_policy(policy), "Pathfinder")
  expect_error(build_v021_resource_policy(retry_chains = 2L), "canonical")
})

test_that("four-chain directions derive at most two concurrent fits", {
  policy <- build_v021_resource_policy()
  derivation <- derive_safe_outer_concurrency(policy, primary_operation_spec(2L))
  expect_identical(derivation$per_fit_simultaneous_chain_slots, 4L)
  expect_identical(derivation$safe_outer_concurrency, 2L)
  expect_identical(derivation$projected_active_cmdstan_chains, 8L)
  expect_true(derivation$compliant)
  expect_error(derive_safe_outer_concurrency(policy, primary_operation_spec(3L)),
               "exceeds")
})

test_that("valid chain configurations use floor division without rounding upward", {
  cases <- data.frame(chains = c(1L, 2L, 3L, 5L, 6L, 10L),
                      expected = c(10L, 5L, 3L, 2L, 1L, 1L))
  for (i in seq_len(nrow(cases))) {
    policy <- build_v021_resource_policy(
      main_chains = cases$chains[[i]], retry_chains = cases$chains[[i]],
      kfold_chains = cases$chains[[i]], confirmation_chains = cases$chains[[i]],
      proposed_outer_concurrency = cases$expected[[i]])
    derivation <- derive_safe_outer_concurrency(
      policy, build_v021_operation_spec(c("main_fit", "retry_fit"), cases$expected[[i]]))
    expect_identical(derivation$safe_outer_concurrency, cases$expected[[i]])
    expect_lte(derivation$projected_active_cmdstan_chains, 10L)
  }
})

test_that("n_workers_outer cannot override the binding chain ceiling", {
  policy <- build_v021_resource_policy()
  expect_error(derive_safe_outer_concurrency(
    policy, build_v021_operation_spec("main_fit", requested_outer_concurrency = 10L)),
    "binding global ceiling")
})

test_that("every operation uses the same slot accounting model", {
  policy <- build_v021_resource_policy()
  slots <- policy$operation_slots
  expect_identical(slots$simultaneous_chain_slots,
                   c(4L, 4L, 0L, 1L, 4L))
  expect_identical(slots$cmdstan_process_slots,
                   c(4L, 4L, 1L, 1L, 4L))
  expect_identical(slots$total_chains, c(4L, 4L, 0L, 4L, 4L))
  expect_identical(policy$pathfinder_chain_slots, 0L)
  expect_identical(policy$pathfinder_processes, 1L)
  expect_identical(policy$pathfinder_num_paths, 8L)

  retry <- derive_safe_outer_concurrency(
    policy, build_v021_operation_spec("retry_fit", 2L))
  kfold <- derive_safe_outer_concurrency(
    policy, build_v021_operation_spec("kfold_fit", 10L))
  confirmation <- derive_safe_outer_concurrency(
    policy, build_v021_operation_spec("confirmation_fit", 2L))
  expect_identical(retry$projected_active_cmdstan_chains, 8L)
  expect_identical(kfold$projected_active_cmdstan_chains, 10L)
  expect_identical(confirmation$projected_active_cmdstan_chains, 8L)
})

test_that("overlapping operation plans fail when chains or processes exceed ten", {
  policy <- build_v021_resource_policy()
  safe <- data.frame(operation = c("main_fit", "pathfinder", "kfold_fit"),
                     concurrent_instances = c(2L, 1L, 1L))
  expect_identical(validate_v021_operation_plan(policy, safe)$active_cmdstan_chains, 9L)
  unsafe_chains <- data.frame(operation = c("main_fit", "retry_fit"),
                              concurrent_instances = c(2L, 1L))
  expect_error(validate_v021_operation_plan(policy, unsafe_chains), "exceeds")
  unsafe_processes <- data.frame(operation = c("main_fit", "pathfinder"),
                                 concurrent_instances = c(2L, 3L))
  expect_error(validate_v021_operation_plan(policy, unsafe_processes), "exceeds")
})

test_that("the exact numerical-thread contract is fixed to one and scoped", {
  expected_names <- c(
    "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "STAN_NUM_THREADS",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS"
  )
  contract <- v021_single_thread_environment()
  expect_identical(names(contract), expected_names)
  expect_true(all(contract == "1"))
  expect_silent(validate_v021_single_thread_environment(contract))
  bad <- contract
  bad[["OPENBLAS_NUM_THREADS"]] <- "2"
  expect_error(validate_v021_single_thread_environment(bad), "exactly one")
  expect_error(validate_v021_single_thread_environment(contract[-1L]), "missing")

  old <- Sys.getenv(expected_names, unset = NA_character_)
  observed <- with_v021_single_thread_environment(function()
    Sys.getenv(expected_names, unset = NA_character_))
  expect_identical(unname(observed), rep("1", length(expected_names)))
  expect_identical(Sys.getenv(expected_names, unset = NA_character_), old)
})

test_that("process accounting includes indirect descendants and operation identity", {
  policy <- build_v021_resource_policy()
  snapshot <- make_process_snapshot(4L, indirect = TRUE)
  map <- data.frame(pid = 301:304, operation = rep("main_fit", 4),
                    task_id = rep("task-1", 4), direction_id = rep("direction-1", 4),
                    stringsAsFactors = FALSE)
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), policy, root_pid = 100L,
    known_model_executables = "/models/pclv", operation_map = map)
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 4L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 4L)
  chains <- monitor$records[monitor$records$classification == "cmdstan_chain", ]
  expect_equal(nrow(chains), 4L)
  expect_true(all(chains$ppid == 200L))
  expect_true(all(chains$operation == "main_fit"))
  expect_true(all(chains$task_id == "task-1"))
  expect_true(all(chains$direction_id == "direction-1"))
})

test_that("unknown potential CmdStan descendants fail closed", {
  policy <- build_v021_resource_policy()
  monitor <- monitor_v021_process_snapshots(
    list(make_process_snapshot(2L, unknown_potential = TRUE)),
    policy, 100L, known_model_executables = "/models/pclv")
  expect_identical(monitor$monitoring_state, "unverified_process_tree")
  expect_identical(monitor$reason, "unknown_potential_cmdstan_descendant")
  expect_identical(monitor$compliance_status, "unverified")
})

test_that("process-tree peak ten passes and greater than ten fails", {
  policy <- build_v021_resource_policy()
  at_ceiling <- monitor_v021_process_snapshots(
    list(make_process_snapshot(4L, "t1"), make_process_snapshot(10L, "t2")),
    policy, 100L, known_model_executables = "/models/pclv")
  expect_identical(at_ceiling$observed_peak_active_cmdstan_chains, 10L)
  expect_identical(at_ceiling$compliance_status, "compliant")
  passed <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), at_ceiling)
  expect_identical(passed$state, "passed")
  expect_true(passed$passed)

  exceeded <- monitor_v021_process_snapshots(
    list(make_process_snapshot(11L)), policy, 100L,
    known_model_executables = "/models/pclv")
  failed <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), exceeded)
  expect_identical(exceeded$compliance_status, "exceeded")
  expect_identical(failed$state, "failed_ceiling_exceeded")
  expect_false(failed$passed)
  expect_match(failed$reasons, "exceeded")

  pathfinder_overlap <- monitor_v021_process_snapshots(
    list(make_process_snapshot(10L, pathfinder = TRUE)), policy, 100L,
    known_model_executables = "/models/pclv")
  expect_identical(pathfinder_overlap$observed_peak_active_cmdstan_chains, 10L)
  expect_identical(pathfinder_overlap$observed_peak_active_cmdstan_processes, 11L)
  overlap_failure <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(),
    pathfinder_overlap)
  expect_identical(overlap_failure$state, "failed_ceiling_exceeded")
})

test_that("unverified and monitoring-error process trees cannot pass preflight", {
  policy <- build_v021_resource_policy()
  incomplete <- monitor_v021_process_snapshots(
    list(make_process_snapshot(2L, complete = FALSE)), policy, 100L,
    known_model_executables = "/models/pclv")
  result <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), incomplete)
  expect_identical(result$state, "failed_unverified_process_tree")
  expect_match(result$reasons, "unreadable")

  monitor_error <- monitor_v021_process_snapshots(list(), policy, 100L)
  result <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), monitor_error)
  expect_identical(result$state, "failed_monitoring_error")
  expect_match(result$reasons, "no_process_snapshots")
})

test_that("incorrect thread environment and unsafe requests have explicit failures", {
  policy <- build_v021_resource_policy()
  monitor <- monitor_v021_process_snapshots(
    list(make_process_snapshot(8L)), policy, 100L,
    known_model_executables = "/models/pclv")
  bad_environment <- v021_single_thread_environment()
  bad_environment[["MKL_NUM_THREADS"]] <- "4"
  thread_failure <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), bad_environment, monitor)
  expect_identical(thread_failure$state, "failed_thread_environment")
  expect_match(thread_failure$reasons, "exactly one")

  unsafe <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(3L), v021_single_thread_environment(), monitor)
  expect_identical(unsafe$state, "failed_invalid_policy")
  expect_match(unsafe$reasons, "exceeds")
  expect_false(unsafe$passed)
})

test_that("preflight states and results are deterministic", {
  expect_setequal(v021_preflight_states, c(
    "passed", "failed_ceiling_exceeded", "failed_unverified_process_tree",
    "failed_thread_environment", "failed_invalid_policy", "failed_monitoring_error"
  ))
  policy <- build_v021_resource_policy()
  monitor <- monitor_v021_process_snapshots(
    list(make_process_snapshot(8L)), policy, 100L,
    known_model_executables = "/models/pclv")
  first <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), monitor)
  second <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), monitor)
  expect_identical(first, second)
  expect_identical(first$preflight_schema, "v021_resource_preflight_v1")
})

test_that("resource policy does not alter V021-03 identities or checkpoint semantics", {
  path <- tempfile("v021-resource-manifest-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  tasks <- data.frame(
    task_id = c(1L, 1L), direction_index = c(1L, 2L),
    target = c("a", "b"), source = c("b", "a"), seed = c(102018L, 102019L)
  )
  manifest <- build_v021_checkpoint_manifest(tasks, 4L, path)
  frozen <- serialize(manifest, NULL, version = 3)
  policy <- build_v021_resource_policy()
  derive_safe_outer_concurrency(policy, primary_operation_spec(2L))
  monitor_v021_process_snapshots(
    list(make_process_snapshot(8L)), policy, 100L,
    known_model_executables = "/models/pclv")
  expect_identical(serialize(manifest, NULL, version = 3), frozen)
  expect_identical(manifest$task_id, c(1L, 1L))
  expect_identical(manifest$direction_index, c(1L, 2L))
  expect_identical(manifest$seed, c(102018L, 102019L))
  expect_identical(manifest$chain_seeds,
                   build_v021_checkpoint_manifest(tasks, 4L, path)$chain_seeds)
  expect_true(all(manifest$execution_state == "pending"))
})

test_that("resource helpers neither modify inference values nor invoke sampling", {
  inference <- list(
    posterior = c(a_ij = -0.2), diagnostic_class = "converged",
    significant = TRUE, bayesian_eligible = TRUE,
    kfold = list(elpd = -3), stacking = c(weight = 1)
  )
  frozen <- serialize(inference, NULL, version = 3)
  policy <- build_v021_resource_policy()
  derive_safe_outer_concurrency(policy, primary_operation_spec(2L))
  expect_identical(serialize(inference, NULL, version = 3), frozen)

  helper_source <- readLines(testthat::test_path(
    "../../benchmarks/mtist/v021_resource_policy.R"), warn = FALSE)
  expect_false(any(grepl("cmdstan_model\\(|\\$sample\\(|\\$pathfinder\\(", helper_source)))
})
