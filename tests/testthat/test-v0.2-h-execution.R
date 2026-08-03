make_task_success <- function(target, partner, ctx, seed_override, progress_local) {
  failure <- pclvbayes:::.pclv_failure("test", "unused", list())
  out <- pclvbayes:::.failed_direction_result(failure)
  out$ok <- TRUE; out$failure <- NULL
  out$n_pairs <- 5L
  out$a_mean <- seed_override / 1e6
  out$retry_history <- list(list(attempt = 1L, seed = seed_override, status = "success"))
  out$initialization_provenance <- list(list(actual = "scalar"))
  out$kfold_subject <- list(stats::setNames(c(-2, -1), c("S2", "S1")))
  out$kfold_subject_ppd <- list(stats::setNames(c(-1, -.5), c("S2", "S1")))
  out$kfold_subject_ids <- list(c("S2", "S1"))
  out$kfold_subject_counts <- list(c(S2 = 2L, S1 = 2L))
  out
}

run_task_fixture <- function(tasks, order = seq_along(tasks), progress = "none", scheduling = list(n_workers_kfold_eff = 1L), run_one = make_task_success) {
  outcomes <- lapply(order, function(k) pclvbayes:::.execute_pair_task(
    tasks[[k]], c("a", "b", "c", "d"), run_one,
    core_ctx = list(kfold_K = 2L, kfold_R = 1L, mod_exe_file = "/canonical/model"),
    scheduling = scheduling, progress = progress, mute_logs = TRUE
  ))
  pclvbayes:::.assemble_pair_outcomes(outcomes, tasks)
}

test_that("canonical pair tasks have stable order and direction seeds", {
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b", "c", "d"), 17L)
  expect_equal(vapply(tasks, `[[`, integer(1), "task_index"), 1:6)
  expect_identical(tasks[[1]]$direction_seeds,
                   c(ij = 17L + 100000L + 2000L + 1L,
                     ji = 17L + 100000L + 2000L + 2L))
  expect_identical(tasks[[1]]$taxon_i, "a")
  expect_identical(tasks[[1]]$taxon_j, "b")
})

test_that("sequential and reordered parallel completion assemble identically", {
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b", "c", "d"), 31L)
  sequential <- run_task_fixture(tasks)
  completed_reverse <- run_task_fixture(tasks, rev(seq_along(tasks)))
  expect_identical(sequential, completed_reverse)
  expect_identical(names(sequential), names(completed_reverse))
  expect_identical(vapply(sequential, typeof, character(1)),
                   vapply(completed_reverse, typeof, character(1)))
  expect_identical(sequential$i, c("a", "a", "a", "b", "b", "c"))
  expect_true(all(sequential$direction_ok_ij & sequential$direction_ok_ji))
  expect_identical(sequential$retry_history_ij, completed_reverse$retry_history_ij)
  expect_identical(sequential$initialization_provenance_ji,
                   completed_reverse$initialization_provenance_ji)
})

test_that("worker count and progress state cannot affect task results or seeds", {
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b", "c", "d"), 73L)
  one <- run_task_fixture(tasks, progress = "none", scheduling = list(n_workers_kfold_eff = 1L))
  many <- run_task_fixture(tasks, progress = "bar", scheduling = list(n_workers_kfold_eff = 3L))
  expect_identical(one, many)
  expect_identical(one$retry_history_ij, many$retry_history_ij)
})

test_that("a failed direction remains structured without deleting other pairs", {
  fail_one <- function(target, partner, ctx, seed_override, progress_local) {
    if (target == "b" && partner == "a")
      return(pclvbayes:::.pclv_failure("sampling", "worker_pair_failed", list()))
    make_task_success(target, partner, ctx, seed_override, progress_local)
  }
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b", "c", "d"), 91L)
  result <- run_task_fixture(tasks, rev(seq_along(tasks)), run_one = fail_one)
  expect_false(result$direction_ok_ji[[1]])
  expect_s3_class(result$failure_ji[[1]], "pclv_failure")
  expect_true(all(result$direction_ok_ij))
  expect_true(all(result$direction_ok_ji[-1]))
  expect_equal(nrow(result), length(tasks))
})

test_that("task executor passes the exact immutable model executable", {
  seen <- list()
  recorder <- function(target, partner, ctx, seed_override, progress_local) {
    seen[[length(seen) + 1L]] <<- list(exe = ctx$mod_exe_file, seed = seed_override,
                                       kfold_workers = ctx$n_workers_kfold_eff)
    make_task_success(target, partner, ctx, seed_override, progress_local)
  }
  task <- pclvbayes:::.make_pair_tasks(c("a", "b"), 5L)[[1L]]
  pclvbayes:::.execute_pair_task(task, c("a", "b"), recorder,
    core_ctx = list(kfold_K = 2L, kfold_R = 1L, mod_exe_file = "/parent/exact-model"),
    scheduling = list(n_workers_kfold_eff = 1L))
  expect_equal(vapply(seen, `[[`, character(1), "exe"), rep("/parent/exact-model", 2))
  expect_identical(vapply(seen, `[[`, integer(1), "seed"), unname(task$direction_seeds))
})

test_that("subject ELPD ordering follows canonical task order", {
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b", "c", "d"), 121L)
  seq_result <- run_task_fixture(tasks)
  rev_result <- run_task_fixture(tasks, rev(seq_along(tasks)))
  expect_identical(pclvbayes:::.expand_cross_subject_elpd(seq_result),
                   pclvbayes:::.expand_cross_subject_elpd(rev_result))
})
