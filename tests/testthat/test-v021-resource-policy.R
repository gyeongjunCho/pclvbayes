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

proc_stat <- function(pid, ppid, start_time, state = "S") paste(
  pid, "(fixture)", state, paste(c(ppid, rep("0", 17L), start_time), collapse = " "))

capture_race_fixture <- function(recheck = c("gone", "same", "reused", "ambiguous",
                                               "zombie"),
                                 malformed = FALSE) {
  recheck <- match.arg(recheck)
  root <- tempfile("v021-proc-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "100")); dir.create(file.path(root, "200"))
  stat_calls <- new.env(parent = emptyenv()); stat_calls$child <- 0L
  stat_reader <- function(path) {
    pid <- basename(dirname(path))
    if (pid == "100") return(proc_stat(100L, 1L, "1000"))
    stat_calls$child <- stat_calls$child + 1L
    if (stat_calls$child == 1L) return(proc_stat(200L, 100L, "2000"))
    if (recheck == "ambiguous") stop("fixture stat unreadable")
    proc_stat(200L, 100L, if (recheck == "reused") "3000" else "2000",
              if (recheck == "zombie") "Z" else "S")
  }
  cmdline_reader <- function(path) {
    pid <- basename(dirname(path))
    if (pid == "100") return(raw_cmdline(c("Rscript", "preflight.R")))
    if (malformed) return(raw_cmdline(c("", "method=sample")))
    stop("fixture cmdline disappeared")
  }
  path_exists <- function(path) {
    pid <- basename(path)
    if (pid == "200" && recheck == "gone") return(FALSE)
    TRUE
  }
  capture_v021_linux_process_snapshot(
    100L, proc_root = root, pid_lister = function(root) c("100", "200"),
    path_exists = path_exists, stat_reader = stat_reader,
    cmdline_reader = cmdline_reader,
    executable_reader = function(path) {
      if (basename(dirname(path)) == "100") "/usr/bin/R" else "/usr/bin/R"
    })
}

capture_recapture_fixture <- function(outcome = c("readable", "zombie", "gone",
                                                   "reused", "persistent")) {
  outcome <- match.arg(outcome)
  root <- tempfile("v021-proc-recapture-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "100")); dir.create(file.path(root, "200"))
  state <- new.env(parent = emptyenv())
  state$stat_calls <- 0L; state$capture_calls <- 0L
  state$rechecks <- 0L; state$delays <- numeric(); state$clock <- 0L
  stat_reader <- function(path) {
    pid <- basename(dirname(path))
    if (pid == "100") return(proc_stat(100L, 1L, "1000"))
    state$stat_calls <- state$stat_calls + 1L
    if (state$stat_calls <= 2L) return(proc_stat(200L, 100L, "2000"))
    if (outcome == "zombie") return(proc_stat(200L, 100L, "2000", "Z"))
    if (outcome == "reused") return(proc_stat(200L, 100L, "3000"))
    proc_stat(200L, 100L, "2000")
  }
  path_exists <- function(path) {
    if (basename(path) != "200") return(TRUE)
    state$rechecks <- state$rechecks + 1L
    !(outcome == "gone" && state$rechecks >= 2L)
  }
  cmdline_reader <- function(path) {
    if (basename(dirname(path)) == "100") return(raw_cmdline(c("Rscript", "preflight.R")))
    state$capture_calls <- state$capture_calls + 1L
    if (outcome == "readable" && state$capture_calls >= 2L)
      return(raw_cmdline(c("R", "--worker")))
    stop("transient cmdline")
  }
  snapshot <- capture_v021_linux_process_snapshot(
    100L, proc_root = root, pid_lister = function(root) c("100", "200"),
    path_exists = path_exists, stat_reader = stat_reader,
    cmdline_reader = cmdline_reader,
    executable_reader = function(path) {
      if (basename(dirname(path)) == "100") return("/usr/bin/R")
      if (outcome == "readable" && state$capture_calls >= 2L) return("/usr/bin/R")
      stop("transient executable")
    }, sleep = function(seconds) state$delays <- c(state$delays, seconds),
    clock = function() {
      state$clock <- state$clock + 1L
      sprintf("2026-08-02T00:00:0%dZ", state$clock)
    })
  list(snapshot = snapshot, state = state)
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

add_cmdstan_diagnostic <- function(snapshot, pid = 900L, ppid = 200L,
                                   executable = "/opt/cmdstan-2.36.0/bin/diagnose",
                                   argv = c("bin/diagnose", "/output/chain-1.csv",
                                            "/output/chain-2.csv"),
                                   potential_cmdstan = TRUE) {
  row <- data.frame(
    timestamp = snapshot$timestamp[[1L]], pid = pid, ppid = ppid,
    command = paste(argv, collapse = " "), executable = executable,
    potential_cmdstan = potential_cmdstan, readable = TRUE,
    stringsAsFactors = FALSE)
  if ("argv" %in% names(snapshot)) row$argv <- I(list(argv))
  row <- row[names(snapshot)]
  rbind(snapshot, row)
}

attempt9_worker_registry <- function(snapshot) {
  worker <- snapshot[snapshot$pid == 200L, , drop = FALSE]
  build_v021_worker_registry(
    pid = worker$pid, start_time = worker$start_time,
    expected_ppid = worker$ppid, worker_role = "direction_fit_worker",
    batch_id = "batch-01", task_id = "task-1",
    registration_timestamp = "2026-08-02T00:00:00Z"
  )
}

add_extended_sampling_chains <- function(snapshot, count, ppid = 200L) {
  template <- snapshot[snapshot$pid == 100L, , drop = FALSE]
  chains <- do.call(rbind, lapply(seq_len(count), function(i) {
    row <- template
    row$pid <- 300L + i; row$ppid <- ppid; row$start_time <- as.character(3000L + i)
    row$process_state <- "S"; row$initial_process_state <- "S"
    row$final_process_state <- "S"
    row$command <- paste("/models/pclv", "method=sample", paste0("id=", i))
    row$argv <- I(list(c("/models/pclv", "method=sample", paste0("id=", i))))
    row$executable <- "/models/pclv"; row$potential_cmdstan <- TRUE
    row$readable <- TRUE; row$capture_state <- "captured"
    row$disappearance_reason <- NA_character_; row$zombie_reason <- NA_character_
    row$capture_retry_count <- 0L; row$capture_retry_timestamps <- I(list(character()))
    row$resolution_reason <- "readable_initial_capture"
    row
  }))
  rbind(snapshot, chains)
}

add_extended_diagnostic <- function(
    snapshot, pid = 900L, ppid = 200L,
    executable = "/opt/cmdstan-2.36.0/bin/diagnose",
    argv = c("bin/diagnose", "/output/chain-1.csv"),
    potential_cmdstan = TRUE) {
  row <- snapshot[snapshot$pid == 100L, , drop = FALSE]
  row$pid <- pid; row$ppid <- ppid; row$start_time <- as.character(pid * 10L)
  row$process_state <- "S"; row$initial_process_state <- "S"
  row$final_process_state <- "S"; row$command <- paste(argv, collapse = " ")
  row$argv <- I(list(argv)); row$executable <- executable
  row$potential_cmdstan <- potential_cmdstan; row$readable <- TRUE
  row$capture_state <- "captured"; row$disappearance_reason <- NA_character_
  row$zombie_reason <- NA_character_; row$capture_retry_count <- 0L
  row$capture_retry_timestamps <- I(list(character()))
  row$resolution_reason <- "readable_initial_capture"
  rbind(snapshot, row)
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

test_that("verified procfs disappearance is retained and consumes no slots", {
  for (scenario in c("gone", "reused")) {
    snapshot <- capture_race_fixture(scenario)
    vanished <- snapshot[snapshot$pid == 200L, , drop = FALSE]
    expect_identical(vanished$capture_state, "vanished_during_capture")
    expect_identical(vanished$classification, NULL)
    expect_identical(vanished$disappearance_reason,
                     if (scenario == "gone") "pid_disappeared" else "pid_reused")
    expect_identical(vanished$start_time, "2000")
    expect_identical(vanished$ppid, 100L)
    expect_true(nzchar(vanished$discovery_time))
    expect_false(vanished$readable)
    expect_identical(vanished$argv[[1L]], character())
    monitor <- monitor_v021_process_snapshots(
      list(snapshot), build_v021_resource_policy(), 100L)
    record <- monitor$records[monitor$records$pid == 200L, , drop = FALSE]
    expect_identical(record$classification, "vanished_during_capture")
    expect_identical(monitor$monitoring_state, "verified")
    expect_identical(monitor$observed_peak_active_cmdstan_chains, 0L)
    expect_identical(monitor$observed_peak_active_cmdstan_processes, 0L)
  }
})

test_that("live unreadable and ambiguous procfs identities fail closed", {
  for (scenario in c("same", "ambiguous")) {
    snapshot <- capture_race_fixture(scenario)
    expect_identical(snapshot$capture_state[snapshot$pid == 200L],
                     if (scenario == "same") "unreadable_live" else "ambiguous_identity")
    monitor <- monitor_v021_process_snapshots(
      list(snapshot), build_v021_resource_policy(), 100L)
    expect_identical(monitor$monitoring_state, "monitoring_error")
    expect_match(monitor$reason, "Malformed decoded process argv")
  }

  malformed <- capture_race_fixture("same", malformed = TRUE)
  monitor <- monitor_v021_process_snapshots(
    list(malformed), build_v021_resource_policy(), 100L)
  expect_identical(monitor$monitoring_state, "monitoring_error")

  root <- tempfile("v021-proc-live-stat-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "100")); dir.create(file.path(root, "200"))
  expect_error(capture_v021_linux_process_snapshot(
    100L, proc_root = root, pid_lister = function(root) c("100", "200"),
    path_exists = function(path) TRUE,
    stat_reader = function(path) {
      if (basename(dirname(path)) == "100") proc_stat(100L, 1L, "1000")
      else stop("unreadable")
    }), "stat remained unreadable")
})

test_that("verified zombies are audited and consume no process capacity", {
  snapshot <- capture_race_fixture("zombie")
  zombie <- snapshot[snapshot$pid == 200L, , drop = FALSE]
  expect_identical(zombie$capture_state, "zombie_process")
  expect_identical(zombie$process_state, "Z")
  expect_identical(zombie$start_time, "2000")
  expect_identical(zombie$ppid, 100L)
  expect_identical(zombie$zombie_reason, "exited_unreaped")
  expect_true(nzchar(zombie$discovery_time))
  expect_false(zombie$readable)
  expect_identical(zombie$argv[[1L]], character())
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L)
  record <- monitor$records[monitor$records$pid == 200L, , drop = FALSE]
  expect_identical(record$classification, "zombie_process")
  expect_false(record$is_active_cmdstan_chain)
  expect_false(record$classification %in% c("completed", "cmdstan_chain"))
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 0L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 0L)

  non_zombie <- capture_race_fixture("same")
  expect_identical(non_zombie$process_state[non_zombie$pid == 200L], "S")
  expect_identical(monitor_v021_process_snapshots(
    list(non_zombie), build_v021_resource_policy(), 100L)$monitoring_state,
    "monitoring_error")

  reused <- capture_race_fixture("reused")
  expect_identical(reused$capture_state[reused$pid == 200L],
                   "vanished_during_capture")
  expect_identical(reused$disappearance_reason[reused$pid == 200L], "pid_reused")
})

test_that("attempt7-shaped zombie snapshot is verified with eight chains", {
  base <- make_process_snapshot(8L, indirect = TRUE)
  base$argv <- I(lapply(base$command,
                        function(x) strsplit(x, "[[:space:]]+")[[1L]]))
  snapshot <- data.frame(
    timestamp = base$timestamp, discovery_time = base$timestamp,
    pid = base$pid, ppid = base$ppid, start_time = as.character(base$pid * 10L),
    process_state = "S", command = base$command, argv = I(base$argv),
    executable = base$executable, potential_cmdstan = base$potential_cmdstan,
    readable = base$readable, capture_state = "captured",
    disappearance_reason = NA_character_, zombie_reason = NA_character_)
  zombie <- snapshot[snapshot$pid == 200L, , drop = FALSE]
  zombie$pid <- 201L; zombie$start_time <- "2010"; zombie$process_state <- "Z"
  zombie$command <- ""; zombie$argv <- I(list(character())); zombie$executable <- ""
  zombie$readable <- FALSE; zombie$capture_state <- "zombie_process"
  zombie$zombie_reason <- "exited_unreaped"
  snapshot <- rbind(snapshot, zombie)
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$compliance_status, "compliant")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 8L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 8L)
  expect_identical(monitor$records$classification[monitor$records$pid == 201L],
                   "zombie_process")
})

test_that("bounded identity-preserving recapture resolves transient procfs states", {
  readable <- capture_recapture_fixture("readable")
  row <- readable$snapshot[readable$snapshot$pid == 200L, , drop = FALSE]
  expect_identical(row$capture_state, "captured")
  expect_identical(row$capture_retry_count, 1L)
  expect_identical(row$resolution_reason, "readable_on_recapture")
  expect_identical(row$initial_process_state, "S")
  expect_identical(row$final_process_state, "S")
  expect_identical(row$capture_retry_timestamps[[1L]], "2026-08-02T00:00:01Z")
  expect_identical(readable$state$delays, .v021_procfs_recapture_delay_seconds)
  monitor <- monitor_v021_process_snapshots(
    list(readable$snapshot), build_v021_resource_policy(), 100L)
  expect_identical(monitor$monitoring_state, "verified")

  zombie <- capture_recapture_fixture("zombie")$snapshot
  expect_identical(zombie$capture_state[zombie$pid == 200L], "zombie_process")
  expect_identical(zombie$capture_retry_count[zombie$pid == 200L], 1L)
  expect_identical(zombie$final_process_state[zombie$pid == 200L], "Z")

  gone <- capture_recapture_fixture("gone")$snapshot
  expect_identical(gone$capture_state[gone$pid == 200L], "vanished_during_capture")
  expect_identical(gone$disappearance_reason[gone$pid == 200L], "pid_disappeared")

  reused <- capture_recapture_fixture("reused")$snapshot
  expect_identical(reused$capture_state[reused$pid == 200L],
                   "vanished_during_capture")
  expect_identical(reused$disappearance_reason[reused$pid == 200L], "pid_reused")
})

test_that("persistent live recapture failure remains unverified", {
  persistent <- capture_recapture_fixture("persistent")
  row <- persistent$snapshot[persistent$snapshot$pid == 200L, , drop = FALSE]
  expect_identical(row$capture_state, "unreadable_live")
  expect_identical(row$capture_retry_count, .v021_procfs_recapture_attempts)
  expect_length(row$capture_retry_timestamps[[1L]], .v021_procfs_recapture_attempts)
  expect_identical(row$resolution_reason, "recapture_limit_exhausted")
  expect_identical(persistent$state$delays,
                   rep(.v021_procfs_recapture_delay_seconds,
                       .v021_procfs_recapture_attempts))
  monitor <- monitor_v021_process_snapshots(
    list(persistent$snapshot), build_v021_resource_policy(), 100L)
  expect_identical(monitor$monitoring_state, "monitoring_error")
  expect_true(is.na(monitor$observed_peak_active_cmdstan_chains))
  expect_true(is.na(monitor$observed_peak_active_cmdstan_processes))
  expect_identical(monitor$compliance_status, "unverified")
})

test_that("registered unreadable outer workers retain identity without consuming slots", {
  persistent <- capture_recapture_fixture("persistent")$snapshot
  persistent$process_state[persistent$pid == 200L] <- "R"
  persistent$initial_process_state[persistent$pid == 200L] <- "R"
  persistent$final_process_state[persistent$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(persistent)
  expect_identical(names(registry), .v021_worker_registry_fields)
  expect_silent(validate_v021_worker_registry(registry))

  monitor <- monitor_v021_process_snapshots(
    list(persistent), build_v021_resource_policy(), 100L,
    worker_registry = registry)
  worker <- monitor$records[monitor$records$pid == 200L, , drop = FALSE]
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$compliance_status, "compliant")
  expect_identical(worker$classification, "registered_outer_worker")
  expect_identical(worker$worker_role, "direction_fit_worker")
  expect_identical(worker$worker_batch_id, "batch-01")
  expect_identical(worker$task_id, "task-1")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 0L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 0L)
})

test_that("registered worker descendants remain independently classified and counted", {
  snapshot <- capture_recapture_fixture("persistent")$snapshot
  snapshot$process_state[snapshot$pid == 200L] <- "R"
  snapshot$initial_process_state[snapshot$pid == 200L] <- "R"
  snapshot$final_process_state[snapshot$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(snapshot)
  snapshot <- add_extended_sampling_chains(snapshot, 8L)
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = registry)
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 8L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 8L)
  expect_equal(sum(monitor$records$classification == "cmdstan_chain"), 8L)
  expect_identical(monitor$records$classification[monitor$records$pid == 200L],
                   "registered_outer_worker")
})

test_that("registered worker identity and ancestry mismatches fail closed", {
  snapshot <- capture_recapture_fixture("persistent")$snapshot
  snapshot$process_state[snapshot$pid == 200L] <- "R"
  snapshot$initial_process_state[snapshot$pid == 200L] <- "R"
  snapshot$final_process_state[snapshot$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(snapshot)

  unregistered <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L)
  expect_identical(unregistered$monitoring_state, "monitoring_error")

  duplicate <- rbind(registry, registry)
  duplicate_result <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L,
    worker_registry = duplicate)
  expect_identical(duplicate_result$monitoring_state, "monitoring_error")
  expect_match(duplicate_result$reason, "duplicate PID")

  reused <- registry; reused$start_time <- "999999"
  reused_result <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L,
    worker_registry = reused)
  expect_identical(reused_result$monitoring_state, "monitoring_error")
  expect_match(reused_result$reason, "mismatch")

  wrong_parent <- registry; wrong_parent$expected_ppid <- 999L
  ancestry_result <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L,
    worker_registry = wrong_parent)
  expect_identical(ancestry_result$monitoring_state, "monitoring_error")
  expect_match(ancestry_result$reason, "mismatch")

  missing_identity <- monitor_v021_process_snapshots(
    list(make_process_snapshot()), build_v021_resource_policy(), 100L,
    worker_registry = registry)
  expect_identical(missing_identity$monitoring_state, "monitoring_error")
  expect_match(missing_identity$reason, "requires procfs identity")
})

test_that("registered worker exception never applies to unreadable CmdStan descendants", {
  snapshot <- capture_recapture_fixture("persistent")$snapshot
  snapshot$process_state[snapshot$pid == 200L] <- "R"
  snapshot$initial_process_state[snapshot$pid == 200L] <- "R"
  snapshot$final_process_state[snapshot$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(snapshot)
  child <- snapshot[snapshot$pid == 200L, , drop = FALSE]
  child$pid <- 301L; child$ppid <- 200L; child$start_time <- "3010"
  child$potential_cmdstan <- TRUE
  snapshot <- rbind(snapshot, child)
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = registry)
  expect_identical(monitor$monitoring_state, "monitoring_error")
  expect_match(monitor$reason, "Malformed decoded process argv")
})

test_that("attempt9-shaped registered worker snapshot is verified and compliant", {
  snapshot <- capture_recapture_fixture("persistent")$snapshot
  snapshot$process_state[snapshot$pid == 200L] <- "R"
  snapshot$initial_process_state[snapshot$pid == 200L] <- "R"
  snapshot$final_process_state[snapshot$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(snapshot)
  snapshot <- add_extended_sampling_chains(snapshot, 8L)
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = registry)
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$compliance_status, "compliant")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 8L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 8L)
})

test_that("canonical diagnose accepts verified readable or unreadable registered ancestry", {
  for (mode in c("readable", "persistent")) {
    snapshot <- capture_recapture_fixture(mode)$snapshot
    if (mode == "persistent") {
      snapshot$process_state[snapshot$pid == 200L] <- "R"
      snapshot$initial_process_state[snapshot$pid == 200L] <- "R"
      snapshot$final_process_state[snapshot$pid == 200L] <- "R"
    }
    registry <- attempt9_worker_registry(snapshot)
    snapshot <- add_extended_diagnostic(snapshot)
    monitor <- monitor_v021_process_snapshots(
      list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv",
      worker_registry = registry)
    parent <- monitor$records[monitor$records$pid == 200L, , drop = FALSE]
    diagnose <- monitor$records[monitor$records$pid == 900L, , drop = FALSE]
    expect_identical(monitor$monitoring_state, "verified")
    expect_identical(diagnose$classification, "cmdstan_diagnostic")
    expect_false(diagnose$is_active_cmdstan_chain)
    expect_identical(parent$classification, "registered_outer_worker")
    expect_false(parent$is_active_cmdstan_chain)
    expect_identical(monitor$observed_peak_active_cmdstan_chains, 0L)
    expect_identical(monitor$observed_peak_active_cmdstan_processes, 1L)
  }
})

test_that("registered diagnose ancestry mismatches remain fail closed", {
  snapshot <- capture_recapture_fixture("persistent")$snapshot
  snapshot$process_state[snapshot$pid == 200L] <- "R"
  snapshot$initial_process_state[snapshot$pid == 200L] <- "R"
  snapshot$final_process_state[snapshot$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(snapshot)
  diagnostic <- add_extended_diagnostic(snapshot)

  unregistered <- monitor_v021_process_snapshots(
    list(diagnostic), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(unregistered$monitoring_state, "monitoring_error")

  reused <- registry; reused$start_time <- "999999"
  mismatch <- monitor_v021_process_snapshots(
    list(diagnostic), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = reused)
  expect_identical(mismatch$monitoring_state, "monitoring_error")

  wrong_ppid <- add_extended_diagnostic(snapshot, ppid = 100L)
  wrong_parent <- monitor_v021_process_snapshots(
    list(wrong_ppid), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = registry)
  expect_identical(wrong_parent$monitoring_state, "verified")
  expect_identical(
    wrong_parent$records$classification[wrong_parent$records$pid == 900L],
    "cmdstan_diagnostic")

  duplicate <- monitor_v021_process_snapshots(
    list(diagnostic), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = rbind(registry, registry))
  expect_identical(duplicate$monitoring_state, "monitoring_error")
  expect_match(duplicate$reason, "duplicate PID")

  arbitrary <- add_extended_diagnostic(snapshot, executable = "/tmp/diagnose")
  arbitrary_result <- monitor_v021_process_snapshots(
    list(arbitrary), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = registry)
  expect_identical(arbitrary_result$monitoring_state, "unverified_process_tree")
})

test_that("registered monitoring is verified with strict twelve-slot utility accounting", {
  base <- capture_recapture_fixture("persistent")$snapshot
  base$process_state[base$pid == 200L] <- "R"
  base$initial_process_state[base$pid == 200L] <- "R"
  base$final_process_state[base$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(base)
  full_sampling <- add_extended_sampling_chains(base, 12L)
  utility_overlap <- add_extended_diagnostic(add_extended_sampling_chains(base, 11L))
  monitor <- monitor_v021_process_snapshots(
    list(full_sampling, utility_overlap), build_v021_resource_policy(), 100L,
    "/models/pclv", worker_registry = registry)
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$compliance_status, "compliant")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 12L)

  exceeded <- monitor_v021_process_snapshots(
    list(add_extended_diagnostic(full_sampling)), build_v021_resource_policy(),
    100L, "/models/pclv", worker_registry = registry)
  expect_identical(exceeded$monitoring_state, "verified")
  expect_identical(exceeded$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(exceeded$observed_peak_active_cmdstan_processes, 13L)
  expect_identical(exceeded$compliance_status, "exceeded")
})

test_that("verified terminal states precede live registered-worker validation", {
  live <- capture_recapture_fixture("persistent")$snapshot
  live$process_state[live$pid == 200L] <- "R"
  live$initial_process_state[live$pid == 200L] <- "R"
  live$final_process_state[live$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(live)

  vanished <- live
  row <- vanished$pid == 200L
  vanished$capture_state[row] <- "vanished_during_capture"
  vanished$disappearance_reason[row] <- "pid_disappeared"
  vanished$resolution_reason[row] <- "pid_disappeared"
  monitor <- monitor_v021_process_snapshots(
    list(vanished), build_v021_resource_policy(), 100L,
    worker_registry = registry)
  record <- monitor$records[monitor$records$pid == 200L, , drop = FALSE]
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(record$classification, "vanished_during_capture")
  expect_identical(record$disappearance_reason, "pid_disappeared")
  expect_false(record$classification %in% c("completed", "registered_outer_worker"))
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 0L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 0L)

  zombie <- capture_race_fixture("zombie")
  zombie_registry <- attempt9_worker_registry(zombie)
  zombie_monitor <- monitor_v021_process_snapshots(
    list(zombie), build_v021_resource_policy(), 100L,
    worker_registry = zombie_registry)
  expect_identical(zombie_monitor$monitoring_state, "verified")
  expect_identical(zombie_monitor$records$classification[
    zombie_monitor$records$pid == 200L], "zombie_process")

  live_monitor <- monitor_v021_process_snapshots(
    list(live), build_v021_resource_policy(), 100L,
    worker_registry = registry)
  expect_identical(live_monitor$monitoring_state, "verified")
  expect_identical(live_monitor$records$classification[
    live_monitor$records$pid == 200L], "registered_outer_worker")
})

test_that("ambiguous registered state remains fail closed after terminal precedence", {
  ambiguous <- capture_race_fixture("ambiguous")
  registry <- build_v021_worker_registry(
    200L, "2000", 100L, "direction_fit_worker", "batch-01", "task-1",
    "2026-08-02T00:00:00Z")
  monitor <- monitor_v021_process_snapshots(
    list(ambiguous), build_v021_resource_policy(), 100L,
    worker_registry = registry)
  expect_identical(monitor$monitoring_state, "monitoring_error")
  expect_match(monitor$reason, "mismatch")
})

test_that("registered exit preserves verified twelve-slot history", {
  live <- capture_recapture_fixture("persistent")$snapshot
  live$process_state[live$pid == 200L] <- "R"
  live$initial_process_state[live$pid == 200L] <- "R"
  live$final_process_state[live$pid == 200L] <- "R"
  registry <- attempt9_worker_registry(live)
  full <- add_extended_sampling_chains(live, 12L)
  vanished <- live
  row <- vanished$pid == 200L
  vanished$capture_state[row] <- "vanished_during_capture"
  vanished$disappearance_reason[row] <- "pid_disappeared"
  vanished$resolution_reason[row] <- "pid_disappeared"
  monitor <- monitor_v021_process_snapshots(
    list(full, vanished), build_v021_resource_policy(), 100L, "/models/pclv",
    worker_registry = registry)
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$compliance_status, "compliant")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 12L)
  expect_identical(tail(monitor$records$classification[
    monitor$records$pid == 200L], 1L), "vanished_during_capture")
})

test_that("attempt8-shaped transient worker resolves without undercounting", {
  base <- make_process_snapshot(8L, indirect = TRUE)
  base$argv <- I(lapply(base$command,
                        function(x) strsplit(x, "[[:space:]]+")[[1L]]))
  resolved <- data.frame(
    timestamp = base$timestamp, discovery_time = base$timestamp,
    pid = base$pid, ppid = base$ppid, start_time = as.character(base$pid * 10L),
    process_state = "S", initial_process_state = "S", final_process_state = "S",
    command = base$command, argv = I(base$argv), executable = base$executable,
    potential_cmdstan = base$potential_cmdstan, readable = TRUE,
    capture_state = "captured", disappearance_reason = NA_character_,
    zombie_reason = NA_character_, capture_retry_count = 0L,
    capture_retry_timestamps = I(rep(list(character()), nrow(base))),
    resolution_reason = "readable_initial_capture")
  retried <- capture_recapture_fixture("readable")$snapshot
  retried <- retried[retried$pid == 200L, , drop = FALSE]
  retried$pid <- 201L; retried$start_time <- "2010"
  resolved <- rbind(resolved, retried)
  monitor <- monitor_v021_process_snapshots(
    list(resolved), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$compliance_status, "compliant")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 8L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 8L)
  expect_identical(monitor$records$capture_retry_count[monitor$records$pid == 201L], 1L)
  expect_false(any(monitor$records$classification %in%
                     c("zombie_process", "vanished_during_capture")))
})

test_that("attempt6-shaped worker exit is verified and compliant", {
  base <- make_process_snapshot(8L, indirect = TRUE)
  base$argv <- I(lapply(base$command,
                        function(x) strsplit(x, "[[:space:]]+")[[1L]]))
  snapshot <- data.frame(
    timestamp = base$timestamp, discovery_time = base$timestamp,
    pid = base$pid, ppid = base$ppid, start_time = as.character(base$pid * 10L),
    command = base$command, argv = I(base$argv), executable = base$executable,
    potential_cmdstan = base$potential_cmdstan, readable = base$readable,
    capture_state = "captured", disappearance_reason = NA_character_)
  vanished <- snapshot[snapshot$pid == 200L, , drop = FALSE]
  vanished$pid <- 201L; vanished$start_time <- "2010"
  vanished$command <- ""; vanished$argv <- I(list(character()))
  vanished$executable <- ""; vanished$readable <- FALSE
  vanished$capture_state <- "vanished_during_capture"
  vanished$disappearance_reason <- "pid_disappeared"
  snapshot <- rbind(snapshot, vanished)
  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$compliance_status, "compliant")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 8L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 8L)
  expect_identical(monitor$records$classification[monitor$records$pid == 201L],
                   "vanished_during_capture")
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

test_that("canonical CmdStan diagnose is a non-chain process slot", {
  snapshot <- make_process_snapshot(8L, indirect = TRUE)
  snapshot$argv <- I(lapply(snapshot$command,
                            function(x) strsplit(x, "[[:space:]]+")[[1L]]))
  snapshot <- snapshot[c("timestamp", "pid", "ppid", "command", "argv",
                         "executable", "potential_cmdstan", "readable")]
  snapshot <- add_cmdstan_diagnostic(snapshot)
  classified <- classify_v021_process_snapshot(
    snapshot, 100L, known_model_executables = "/models/pclv")
  diagnostic <- classified[classified$pid == 900L, , drop = FALSE]
  expect_identical(diagnostic$classification, "cmdstan_diagnostic")
  expect_false(diagnostic$is_active_cmdstan_chain)

  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 8L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 9L)
  expect_identical(monitor$compliance_status, "compliant")

  exceeded <- add_cmdstan_diagnostic(make_process_snapshot(12L, indirect = TRUE))
  monitor <- monitor_v021_process_snapshots(
    list(exceeded), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(monitor$observed_peak_active_cmdstan_processes, 13L)
  expect_identical(monitor$compliance_status, "exceeded")
})

test_that("diagnose recognition accepts canonical controller-owned utility", {
  base <- make_process_snapshot(0L)
  unrelated <- add_cmdstan_diagnostic(
    base, executable = "/tmp/diagnose", potential_cmdstan = FALSE)
  classified <- classify_v021_process_snapshot(unrelated, 100L, "/models/pclv")
  expect_identical(classified$classification[classified$pid == 900L],
                   "other_descendant")

  wrong_parent <- add_cmdstan_diagnostic(base, ppid = 100L)
  monitor <- monitor_v021_process_snapshots(
    list(wrong_parent), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "verified")
  expect_identical(
    monitor$records$classification[monitor$records$pid == 900L],
    "cmdstan_diagnostic")

  unknown <- add_cmdstan_diagnostic(
    base, executable = "/opt/cmdstan-2.36.0/bin/stansummary",
    argv = c("bin/stansummary", "/output/chain-1.csv"))
  monitor <- monitor_v021_process_snapshots(
    list(unknown), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "unverified_process_tree")
  expect_identical(monitor$reason, "unknown_potential_cmdstan_descendant")
})

test_that("the policy records the exact 16-thread, 4-reserved, 12-slot contract", {
  policy <- build_v021_resource_policy()
  expect_identical(policy$policy_schema, "v021_resource_policy_v4")
  expect_identical(policy$logical_host_threads, 16L)
  expect_identical(policy$reserved_host_threads, 4L)
  expect_identical(policy$usable_chain_slots, 12L)
  expect_identical(policy$maximum_active_cmdstan_chains, 12L)
  expect_identical(policy$maximum_cmdstan_process_slots, 12L)
  expect_identical(policy$cpu_threads_per_active_chain, 1L)
  expect_identical(policy$numerical_library_threads, 1L)
  expect_identical(policy$main_chains, 4L)
  expect_identical(policy$proposed_outer_concurrency, 3L)
  expect_true(policy$global_slot_scheduler)
  expect_silent(validate_v021_resource_policy(policy))
  expect_identical(policy, build_v021_resource_policy())
})

test_that("malformed and contradictory resource policies are rejected", {
  expect_error(build_v021_resource_policy(logical_host_threads = 15L), "Usable")
  expect_error(build_v021_resource_policy(reserved_host_threads = 5L), "Usable")
  expect_error(build_v021_resource_policy(maximum_cmdstan_process_slots = 9L),
               "ceilings")
  expect_error(build_v021_resource_policy(cpu_threads_per_active_chain = 2L), "Exactly one")
  expect_error(build_v021_resource_policy(kfold_chains = 2L,
                                          kfold_parallel_chains = 3L), "cannot exceed")
  expect_error(build_v021_resource_policy(global_slot_scheduler = FALSE), "requires")
  expect_error(build_v021_resource_policy(proposed_outer_concurrency = 4L), "exceeds")
  policy <- build_v021_resource_policy()
  policy$usable_chain_slots <- 11L
  expect_error(validate_v021_resource_policy(policy), "16/4/12")
  policy <- build_v021_resource_policy()
  policy$operation_slots$simultaneous_chain_slots[[1L]] <- 3L
  expect_error(validate_v021_resource_policy(policy), "contradicts")
  policy <- build_v021_resource_policy()
  policy$pathfinder_processes <- 0L
  policy$operation_slots$cmdstan_process_slots[[3L]] <- 0L
  expect_error(validate_v021_resource_policy(policy), "Pathfinder")
  expect_error(build_v021_resource_policy(retry_chains = 2L),
               "canonical|Parallel chains")
})

test_that("four-chain directions derive at most three concurrent fits", {
  policy <- build_v021_resource_policy()
  derivation <- derive_safe_outer_concurrency(policy, primary_operation_spec(3L))
  expect_identical(derivation$per_fit_simultaneous_chain_slots, 4L)
  expect_identical(derivation$safe_outer_concurrency, 3L)
  expect_identical(derivation$projected_active_cmdstan_chains, 12L)
  expect_identical(derivation$projected_active_cmdstan_processes, 12L)
  expect_true(derivation$compliant)
  expect_identical(4L * policy$main_chains, 16L)
  expect_error(derive_safe_outer_concurrency(policy, primary_operation_spec(4L)),
               "exceeds")
})

test_that("valid chain configurations use floor division without rounding upward", {
  cases <- data.frame(chains = c(1L, 2L, 3L, 4L, 5L, 12L),
                      expected = c(12L, 6L, 4L, 3L, 2L, 1L))
  for (i in seq_len(nrow(cases))) {
    policy <- build_v021_resource_policy(
      main_chains = cases$chains[[i]], retry_chains = cases$chains[[i]],
      kfold_chains = cases$chains[[i]], confirmation_chains = cases$chains[[i]],
      proposed_outer_concurrency = cases$expected[[i]])
    derivation <- derive_safe_outer_concurrency(
      policy, build_v021_operation_spec(c("main_fit", "retry_fit"), cases$expected[[i]]))
    expect_identical(derivation$safe_outer_concurrency, cases$expected[[i]])
    expect_lte(derivation$projected_active_cmdstan_chains, 12L)
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
    policy, build_v021_operation_spec("kfold_fit", 12L))
  confirmation <- derive_safe_outer_concurrency(
    policy, build_v021_operation_spec("confirmation_fit", 2L))
  expect_identical(retry$projected_active_cmdstan_chains, 8L)
  expect_identical(kfold$projected_active_cmdstan_chains, 12L)
  expect_identical(confirmation$projected_active_cmdstan_chains, 8L)
})

test_that("overlapping operation plans fail when chains or processes exceed twelve", {
  policy <- build_v021_resource_policy()
  safe <- data.frame(operation = c("main_fit", "pathfinder", "kfold_fit"),
                     concurrent_instances = c(2L, 1L, 1L))
  expect_identical(validate_v021_operation_plan(policy, safe)$active_cmdstan_chains, 9L)
  unsafe_chains <- data.frame(operation = c("main_fit", "retry_fit"),
                              concurrent_instances = c(3L, 1L))
  expect_error(validate_v021_operation_plan(policy, unsafe_chains), "exceeds")
  unsafe_processes <- data.frame(operation = c("main_fit", "pathfinder"),
                                 concurrent_instances = c(3L, 1L))
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

test_that("unknown descendants retain a complete atomic post-mortem record", {
  snapshot <- capture_recapture_fixture("readable")$snapshot
  unknown <- snapshot[snapshot$pid == 100L, , drop = FALSE]
  unknown$pid <- 999L
  unknown$ppid <- 200L
  unknown$start_time <- "9990"
  unknown$process_state <- "S"
  unknown$initial_process_state <- "S"
  unknown$final_process_state <- "S"
  unknown$command <- "/models/pclv method=mystery id=9"
  unknown$argv <- I(list(c("/models/pclv", "method=mystery", "id=9")))
  unknown$executable <- "/models/pclv"
  unknown$potential_cmdstan <- TRUE
  unknown$readable <- TRUE
  unknown$capture_state <- "captured"
  unknown$disappearance_reason <- NA_character_
  unknown$zombie_reason <- NA_character_
  unknown$capture_retry_count <- 0L
  unknown$capture_retry_timestamps <- I(list(character()))
  unknown$resolution_reason <- "readable_initial_capture"
  snapshot <- rbind(snapshot, unknown)

  monitor <- monitor_v021_process_snapshots(
    list(snapshot), build_v021_resource_policy(), 100L, "/models/pclv")
  expect_identical(monitor$monitoring_state, "unverified_process_tree")
  expect_identical(monitor$reason, "unknown_potential_cmdstan_descendant")
  offender <- monitor$offending_processes
  expect_identical(names(offender), c(
    "pid", "ppid", "start_time", "argv", "command", "executable",
    "classification", "capture_state", "ancestry", "classifier_stage", "reason"))
  expect_identical(nrow(offender), 1L)
  expect_identical(offender$pid, 999L)
  expect_identical(offender$ppid, 200L)
  expect_identical(offender$start_time, "9990")
  expect_identical(offender$argv[[1L]], c("/models/pclv", "method=mystery", "id=9"))
  expect_identical(offender$command, "/models/pclv method=mystery id=9")
  expect_identical(offender$executable, "/models/pclv")
  expect_identical(offender$classification, "unknown_potential_cmdstan")
  expect_identical(offender$capture_state, "captured")
  expect_identical(offender$ancestry[[1L]], c(100L, 200L, 999L))
  expect_identical(offender$classifier_stage, "operation_classification")
  expect_identical(offender$reason, "unknown_potential_cmdstan_descendant")

  path <- tempfile("v021-monitor-failure-", fileext = ".rds")
  .v021_atomic_save_rds(list(latest_monitor = monitor), path)
  persisted <- readRDS(path)$latest_monitor$offending_processes
  expect_identical(persisted, offender)
  expect_false(any(grepl("tmp-", list.files(dirname(path), basename(path)))))
})

test_that("process-tree peak twelve passes and greater than twelve fails", {
  policy <- build_v021_resource_policy()
  at_ceiling <- monitor_v021_process_snapshots(
    list(make_process_snapshot(4L, "t1"), make_process_snapshot(12L, "t2")),
    policy, 100L, known_model_executables = "/models/pclv")
  expect_identical(at_ceiling$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(at_ceiling$compliance_status, "compliant")
  passed <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), at_ceiling)
  expect_identical(passed$state, "passed")
  expect_true(passed$passed)

  exceeded <- monitor_v021_process_snapshots(
    list(make_process_snapshot(13L)), policy, 100L,
    known_model_executables = "/models/pclv")
  failed <- evaluate_v021_resource_preflight(
    policy, primary_operation_spec(2L), v021_single_thread_environment(), exceeded)
  expect_identical(exceeded$compliance_status, "exceeded")
  expect_identical(failed$state, "failed_ceiling_exceeded")
  expect_false(failed$passed)
  expect_match(failed$reasons, "exceeded")

  pathfinder_overlap <- monitor_v021_process_snapshots(
    list(make_process_snapshot(12L, pathfinder = TRUE)), policy, 100L,
    known_model_executables = "/models/pclv")
  expect_identical(pathfinder_overlap$observed_peak_active_cmdstan_chains, 12L)
  expect_identical(pathfinder_overlap$observed_peak_active_cmdstan_processes, 13L)
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
    policy, primary_operation_spec(4L), v021_single_thread_environment(), monitor)
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
