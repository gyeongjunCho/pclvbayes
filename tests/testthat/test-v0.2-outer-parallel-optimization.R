test_that("outer dispatcher uses explicit globals and restores its plan", {
  state <- new.env(parent = emptyenv())
  state$current <- "caller"
  state$events <- list()
  plan <- function(strategy, workers) {
    if (missing(strategy)) return(state$current)
    state$events[[length(state$events) + 1L]] <- list(
      strategy = strategy, workers = if (missing(workers)) NULL else workers)
    state$current <- strategy
    invisible(strategy)
  }
  captured <- NULL
  map <- function(.x, .f, ..., .options) {
    captured <<- list(x = .x, f = .f, dots = list(...), options = .options)
    rev(.x)
  }
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b", "c"), 19L)
  out <- pclvbayes:::.outer_pair_map(
    tasks, c("a", "b", "c"), list(mod_exe_file = "/canonical/pclv"),
    list(n_workers_kfold_eff = 1L), "none", 2L, FALSE,
    .plan = plan, .map = map)
  expect_identical(out, rev(tasks))
  expect_false(captured$options$globals)
  expect_identical(captured$options$packages, "pclvbayes")
  expect_identical(captured$options$scheduling, Inf)
  expect_true(captured$options$seed)
  expect_identical(captured$dots$taxa_vec, c("a", "b", "c"))
  expect_identical(captured$dots$core_ctx$mod_exe_file, "/canonical/pclv")
  expect_identical(captured$dots$scheduling$n_workers_kfold_eff, 1L)
  expect_identical(state$current, "caller")
  expect_length(state$events, 2L)
  expect_identical(state$events[[1L]]$workers, 2L)
  expect_null(state$events[[2L]]$workers)
})

test_that("outer dispatcher restores caller plan after errors", {
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b"), 29L)
  run_case <- function(condition) {
    state <- new.env(parent = emptyenv()); state$current <- "caller"; state$n <- 0L
    plan <- function(strategy, workers) {
      if (missing(strategy)) return(state$current)
      state$n <- state$n + 1L; state$current <- strategy; invisible(strategy)
    }
    map <- function(...) stop(condition)
    expect_condition(pclvbayes:::.outer_pair_map(
      tasks, c("a", "b"), list(), list(n_workers_kfold_eff = 1L),
      "none", 2L, FALSE, .plan = plan, .map = map),
      class = class(condition)[[1L]])
    expect_identical(state$current, "caller")
    expect_identical(state$n, 2L)
  }
  run_case(simpleError("worker error"))
  run_case(structure(list(message = "simulated interrupt"),
                     class = c("simulated_interrupt", "error", "condition")))
})

test_that("progress dispatch retains task and seed payload", {
  state <- new.env(parent = emptyenv()); state$current <- "caller"
  plan <- function(strategy, workers) {
    if (missing(strategy)) return(state$current)
    state$current <- strategy; invisible(strategy)
  }
  captured <- NULL
  map <- function(.x, .f, ..., .options) {
    captured <<- list(x = .x, f = .f, dots = list(...), options = .options); .x
  }
  tasks <- pclvbayes:::.make_pair_tasks(c("a", "b", "c"), 37L)
  out <- pclvbayes:::.outer_pair_map(
    tasks, c("a", "b", "c"), list(), list(n_workers_kfold_eff = 1L),
    "bar", 2L, TRUE, .plan = plan, .map = map,
    .with_progress = function(expr) force(expr))
  expect_identical(out, tasks)
  expect_identical(captured$x, tasks)
  expect_identical(captured$f, pclvbayes:::.execute_pair_task_progress)
  expect_identical(captured$options$scheduling, Inf)
  expect_false(captured$options$globals)
  expect_true(is.function(captured$dots$progressor))
})

test_that("outer optimization does not add public arguments", {
  expect_false(any(c("outer_scheduling", "chunk_size", "globals") %in%
                     names(formals(fit_pclv_bayes))))
})
