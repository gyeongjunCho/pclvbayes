source(testthat::test_path("../../benchmarks/mtist/v021_checkpoint_manifest.R"))

make_v021_checkpoint_tasks <- function() {
  data.frame(
    task_id = rep(1:3, each = 2),
    direction_index = 1:6,
    target = c("a", "b", "a", "c", "b", "c"),
    source = c("b", "a", "c", "a", "c", "b"),
    seed = c(102018L, 102019L, 103018L, 103019L, 203018L, 203019L),
    stringsAsFactors = FALSE
  )
}

with_v021_checkpoint_dir <- function(code) {
  path <- tempfile("v021-checkpoints-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  force(code)(path)
}

retry_fixture <- function(seed, statuses = c("failed", "success")) {
  lapply(seq_along(statuses), function(i)
    list(attempt = as.integer(i), seed = as.integer(seed + i - 1L),
         status = statuses[[i]]))
}

test_that("checkpoint manifests and every identifier are deterministic", {
  with_v021_checkpoint_dir(function(path) {
    tasks <- make_v021_checkpoint_tasks()[c(5, 2, 6, 1, 4, 3), ]
    first <- build_v021_checkpoint_manifest(tasks, chains = 4L, checkpoint_dir = path)
    second <- build_v021_checkpoint_manifest(tasks, chains = 4L, checkpoint_dir = path)
    expect_identical(first, second)
    expect_identical(first$task_id, rep(1:3, each = 2))
    expect_identical(first$direction_index, 1:6)
    expect_identical(first$pair_id, rep(sprintf("pair-%06d", 1:3), each = 2))
    expect_identical(first$direction_id, sprintf("direction-%06d", 1:6))
    expect_identical(first$seed, make_v021_checkpoint_tasks()$seed)
    expect_identical(first$chain_seeds[[1L]],
                     first$seed[[1L]] + as.integer(10000000 * 1:4))
    expect_identical(anyDuplicated(unlist(first$chain_seeds, use.names = FALSE)), 0L)
    expect_true(all(first$execution_state == "pending"))
    expect_true(all(first$elapsed_time == 0))
    expect_silent(validate_v021_checkpoint_manifest(first))
  })
})

test_that("caller pair identifiers are retained without changing canonical directions", {
  with_v021_checkpoint_dir(function(path) {
    tasks <- make_v021_checkpoint_tasks()
    tasks$pair_id <- rep(c("dataset-A:a|b", "dataset-A:a|c", "dataset-A:b|c"), each = 2)
    manifest <- build_v021_checkpoint_manifest(tasks, 2L, path)
    expect_identical(manifest$pair_id, tasks$pair_id)
    expect_identical(manifest$target, tasks$target)
    expect_identical(manifest$source, tasks$source)
    expect_identical(manifest$direction_index, tasks$direction_index)
  })
})

test_that("manifest construction rejects duplicate and incomplete identities", {
  with_v021_checkpoint_dir(function(path) {
    tasks <- make_v021_checkpoint_tasks()
    duplicate <- tasks
    duplicate$direction_index[[2L]] <- duplicate$direction_index[[1L]]
    expect_error(build_v021_checkpoint_manifest(duplicate, 2L, path), "unique")
    expect_error(build_v021_checkpoint_manifest(tasks[-1, ], 2L, path), "both directions")
    conflict <- tasks
    conflict$pair_id <- c("same", "same", "same", "same", "pair-3", "pair-3")
    expect_error(build_v021_checkpoint_manifest(conflict, 2L, path), "conflicting")
  })
})

test_that("manifest and checkpoint writes are atomic and leave no temporary files", {
  with_v021_checkpoint_dir(function(path) {
    manifest <- build_v021_checkpoint_manifest(make_v021_checkpoint_tasks(), 2L, path)
    manifest_path <- file.path(path, "manifest.rds")
    write_v021_manifest_atomic(manifest, manifest_path)
    expect_identical(read_v021_checkpoint_manifest(manifest_path), manifest)
    expect_length(list.files(path, pattern = "\\.tmp-", all.files = TRUE), 0L)

    checkpoint <- .v021_record_from_row(
      manifest[1, , drop = FALSE], "completed", elapsed_time = 3.5,
      retry_history = retry_fixture(manifest$seed[[1L]]), pathfinder_used = TRUE,
      completion_timestamp = "2026-08-02 10:00:00 UTC",
      result = list(coefficient = -0.25)
    )
    write_v021_checkpoint_atomic(checkpoint)
    before <- unname(tools::md5sum(checkpoint$output_location))
    expect_error(write_v021_checkpoint_atomic(checkpoint), "completed")
    expect_identical(unname(tools::md5sum(checkpoint$output_location)), before)
    expect_identical(readRDS(checkpoint$output_location)$result,
                     list(coefficient = -0.25))
    expect_length(list.files(path, pattern = "\\.tmp-", all.files = TRUE), 0L)
  })
})

test_that("invalid checkpoint replacement cannot alter an existing output", {
  with_v021_checkpoint_dir(function(path) {
    manifest <- build_v021_checkpoint_manifest(make_v021_checkpoint_tasks(), 2L, path)
    running <- .v021_record_from_row(manifest[1, , drop = FALSE], "running")
    write_v021_checkpoint_atomic(running)
    before <- readBin(running$output_location, "raw", n = file.info(running$output_location)$size)
    invalid <- running
    invalid$direction_index <- 999L
    expect_error(write_v021_checkpoint_atomic(invalid), "conflicting identity")
    after <- readBin(running$output_location, "raw", n = file.info(running$output_location)$size)
    expect_identical(after, before)
  })
})

test_that("restart reconciles explicit states and resumes only unfinished tasks", {
  with_v021_checkpoint_dir(function(path) {
    manifest <- build_v021_checkpoint_manifest(make_v021_checkpoint_tasks(), 2L, path)
    completed <- .v021_record_from_row(
      manifest[1, , drop = FALSE], "completed", elapsed_time = 9,
      completion_timestamp = "2026-08-02 09:00:00 UTC", result = list(draws = c(1, 2))
    )
    failed <- .v021_record_from_row(
      manifest[2, , drop = FALSE], "failed", elapsed_time = 4,
      retry_history = retry_fixture(manifest$seed[[2L]], c("failed", "failed")),
      pathfinder_used = TRUE, failure_reason = "sampler_process_failed"
    )
    incomplete <- .v021_record_from_row(
      manifest[3, , drop = FALSE], "incomplete", elapsed_time = 2,
      interrupted = TRUE, interruption_reason = "host_shutdown",
      failure_reason = "host_shutdown"
    )
    running <- .v021_record_from_row(manifest[4, , drop = FALSE], "running", elapsed_time = 1)
    skipped <- .v021_record_from_row(
      manifest[5, , drop = FALSE], "skipped", failure_reason = "prospective_skip"
    )
    invisible(lapply(list(completed, failed, incomplete, running, skipped),
                     write_v021_checkpoint_atomic))

    completed_hash <- unname(tools::md5sum(completed$output_location))
    failed_hash <- unname(tools::md5sum(failed$output_location))
    plan <- plan_v021_checkpoint_restart(manifest)
    expect_identical(plan$manifest$execution_state,
                     c("completed", "failed", "incomplete", "incomplete", "skipped", "pending"))
    expect_identical(plan$resume_indices, c(3L, 4L, 6L))
    expect_identical(plan$completed_indices, 1L)
    expect_true(plan$manifest$interrupted[[4L]])
    expect_identical(plan$manifest$interruption_reason[[4L]], "interrupted_while_running")

    calls <- integer()
    worker <- function(row) {
      calls <<- c(calls, row$direction_index[[1L]])
      list(execution_state = "completed", elapsed_time = row$direction_index[[1L]] / 10,
           retry_history = list(list(attempt = 1L, seed = row$seed[[1L]], status = "success")),
           pathfinder_used = FALSE,
           result = list(direction_index = row$direction_index[[1L]], retained = TRUE))
    }
    manifest_path <- file.path(path, "manifest.rds")
    resumed <- execute_v021_checkpoint_restart(
      manifest, worker, manifest_path,
      timestamp = function() "2026-08-02 12:00:00 UTC"
    )
    expect_identical(calls, c(3L, 4L, 6L))
    expect_identical(resumed$direction_index, 1:6)
    expect_identical(resumed$execution_state,
                     c("completed", "failed", "completed", "completed", "skipped", "completed"))
    expect_identical(resumed$failure_reason[[2L]], "sampler_process_failed")
    expect_identical(resumed$retry_count[[2L]], 1L)
    expect_true(resumed$pathfinder_used[[2L]])
    expect_equal(resumed$elapsed_time[[3L]], 2.3)
    expect_equal(resumed$elapsed_time[[4L]], 1.4)
    expect_true(resumed$interrupted[[3L]])
    expect_true(resumed$interrupted[[4L]])
    expect_identical(resumed$interruption_reason[[3L]], "host_shutdown")
    expect_identical(resumed$interruption_reason[[4L]], "interrupted_while_running")
    expect_identical(unname(tools::md5sum(completed$output_location)), completed_hash)
    expect_identical(unname(tools::md5sum(failed$output_location)), failed_hash)
    expect_true(file.exists(failed$output_location))
    expect_identical(read_v021_checkpoint_manifest(manifest_path), resumed)
  })
})

test_that("a second restart never reruns or rewrites completed tasks", {
  with_v021_checkpoint_dir(function(path) {
    manifest <- build_v021_checkpoint_manifest(make_v021_checkpoint_tasks(), 2L, path)
    calls <- integer()
    worker <- function(row) {
      calls <<- c(calls, row$direction_index[[1L]])
      list(execution_state = "completed", elapsed_time = 0.1,
           result = list(value = row$direction_index[[1L]]))
    }
    first <- execute_v021_checkpoint_restart(
      manifest, worker, timestamp = function() "2026-08-02 12:00:00 UTC"
    )
    hashes <- unname(tools::md5sum(first$output_location))
    calls <- integer()
    second <- execute_v021_checkpoint_restart(
      first, worker, timestamp = function() "2026-08-03 12:00:00 UTC"
    )
    expect_length(calls, 0L)
    expect_identical(second, first)
    expect_identical(unname(tools::md5sum(second$output_location)), hashes)
  })
})

test_that("worker failures remain structured checkpoints and are not retried", {
  with_v021_checkpoint_dir(function(path) {
    manifest <- build_v021_checkpoint_manifest(make_v021_checkpoint_tasks()[1:2, ], 2L, path)
    calls <- integer()
    worker <- function(row) {
      calls <<- c(calls, row$direction_index[[1L]])
      if (row$direction_index[[1L]] == 1L) stop("fixture failure")
      list(execution_state = "incomplete", elapsed_time = 1.25,
           interrupted = TRUE, interruption_reason = "fixture interruption",
           failure_reason = "fixture interruption")
    }
    first <- execute_v021_checkpoint_restart(manifest, worker)
    expect_identical(first$execution_state, c("failed", "incomplete"))
    expect_identical(first$failure_reason, c("fixture failure", "fixture interruption"))
    expect_true(first$interrupted[[2L]])
    expect_equal(first$elapsed_time[[2L]], 1.25)
    failed_hash <- unname(tools::md5sum(first$output_location[[1L]]))

    calls <- integer()
    resumed <- execute_v021_checkpoint_restart(
      first,
      function(row) {
        calls <<- c(calls, row$direction_index[[1L]])
        list(execution_state = "completed", result = list(ok = TRUE))
      },
      timestamp = function() "2026-08-02 13:00:00 UTC"
    )
    expect_identical(calls, 2L)
    expect_identical(resumed$execution_state, c("failed", "completed"))
    expect_identical(unname(tools::md5sum(first$output_location[[1L]])), failed_hash)
  })
})

test_that("completed manifest rows without retained outputs fail closed", {
  with_v021_checkpoint_dir(function(path) {
    manifest <- build_v021_checkpoint_manifest(make_v021_checkpoint_tasks()[1:2, ], 2L, path)
    manifest$execution_state[[1L]] <- "completed"
    manifest$completion_timestamp[[1L]] <- "2026-08-02 12:00:00 UTC"
    expect_error(plan_v021_checkpoint_restart(manifest), "no checkpoint")
  })
})
