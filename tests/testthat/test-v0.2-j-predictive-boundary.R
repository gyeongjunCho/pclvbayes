predictive_fixture <- function(diagnostic_failure = NULL) {
  list(
    a_mean = -.2, a_sd = .03, a_q2.5 = -.26, a_q97.5 = -.14,
    p_sign2 = .01, diagnostic_class = if (is.null(diagnostic_failure)) "converged" else "interaction_indeterminate",
    interaction_identifiable = is.null(diagnostic_failure),
    residual_identifiable = is.null(diagnostic_failure),
    diagnostic_failure = list(diagnostic_failure),
    retry_history = list(list(list(attempt = 0L, seed = 101L))),
    kfold = NULL, kfold_mean = NA_real_, kfold_method = NA_character_,
    kfold_outer_rounds = 0L, kfold_failed = TRUE,
    kfold_n_folds_ok = NA_integer_, kfold_n_folds_fail = NA_integer_,
    kfold_subject = list(NULL), kfold_subject_ppd = list(NULL),
    kfold_subject_ids = list(NULL), kfold_subject_counts = list(NULL),
    kfold_subject_success = list(NULL), kfold_subject_fail = list(NULL),
    kfold_success_total = NA_integer_, kfold_failures = list(NULL),
    kfold_splits = list(NULL), kfold_seed_used = NA_integer_,
    kfold_K = NA_integer_, kfold_R = NA_integer_, kfold_agg = "subject-uniform",
    kfold_sd = NA_real_, kfold_se = NA_real_, kfold_n_subjects = NA_integer_,
    kfold_retry_total = NA_integer_, kfold_retry_mean = NA_real_,
    kfold_nu_fold_means = list(NULL),
    .predictive_context = list(
      mod = structure(list(), class = "mock_model"), stan_list = list(N = 2L),
      sample_args = list(seed = 101L, step_size = 1, inv_metric = 2, metric_file = "x"),
      pair_in = data.frame(subject = c("a", "b")), K = 2L, R = 1L,
      seed = 101L, silent_sampler = TRUE, n_workers_kfold = 1L,
      max_retries = 0L, min_pairs = 1L, pair_tag = "b->a", progress = "none"
    )
  )
}

kfold_fixture <- function() list(
  K = 2L, R = 1L, elpd_mean = -1.5, elpd_method = "kalman-ou",
  n_folds_ok = 2L, n_folds_fail = 0L,
  elpd_subject = c(a = -1, b = -2), elpd_subject_ppd = c(a = -.5, b = -1),
  subject_test_counts = c(a = 1L, b = 1L),
  subject_success_counts = c(a = 1L, b = 1L),
  subject_failure_counts = c(a = 0L, b = 0L),
  total_successful_evaluations = 2L, failures = list(),
  splits_df = data.frame(r = 1L, k = 1L), retry_total = 0L,
  retry_mean = 0, nu_fold_means = c(4, 4)
)

test_that("predictive phase preserves Bayesian evidence and canonical K-fold payload", {
  calls <- 0L
  testthat::local_mocked_bindings(.repkfold_eval = function(...) {
    calls <<- calls + 1L
    args <- list(...)
    expect_identical(args$seed, 101L)
    expect_null(args$sample_args_base$step_size)
    expect_null(args$sample_args_base$inv_metric)
    expect_null(args$sample_args_base$metric_file)
    kfold_fixture()
  }, .package = "pclvbayes")
  main <- predictive_fixture()
  evidence <- main[c("a_mean", "a_sd", "a_q2.5", "a_q97.5", "p_sign2",
                     "diagnostic_class", "interaction_identifiable",
                     "residual_identifiable", "retry_history")]
  full <- pclvbayes:::.add_predictive_evaluation(main)
  expect_identical(calls, 1L)
  expect_identical(full[names(evidence)], evidence)
  expect_equal(full$kfold_mean, -1.5)
  expect_identical(full$kfold_seed_used, 101L)
  expect_equal(full$kfold_subject[[1L]], c(a = -1, b = -2))
})

test_that("main-posterior boundary contains no predictive invocation", {
  body_text <- paste(deparse(body(pclvbayes:::.fit_direction_main_posterior)), collapse = "\n")
  expect_false(grepl("repkfold_eval|kalman_ou_loglik|compute_stacking", body_text, ignore.case = TRUE))
})

test_that("ineligible directions bypass predictive evaluation unchanged", {
  calls <- 0L
  testthat::local_mocked_bindings(.repkfold_eval = function(...) {
    calls <<- calls + 1L; stop("must not run")
  }, .package = "pclvbayes")
  failure <- pclvbayes:::.pclv_failure("posterior_diagnostics", "interaction_indeterminate")
  main <- predictive_fixture(failure)
  full <- pclvbayes:::.add_predictive_evaluation(main)
  expect_identical(calls, 0L)
  expect_identical(full$diagnostic_class, main$diagnostic_class)
  expect_false(".predictive_context" %in% names(full))
})

test_that("full directed wrapper composes main and predictive phases once", {
  calls <- character()
  main <- predictive_fixture()
  testthat::local_mocked_bindings(
    .fit_direction_main_posterior = function(...) { calls <<- c(calls, "main"); main },
    .add_predictive_evaluation = function(x) { calls <<- c(calls, "predictive"); x },
    .package = "pclvbayes")
  out <- pclvbayes:::.run_one("a", "b", list(), 101L, "none")
  expect_identical(calls, c("main", "predictive"))
  expect_identical(out, main)
})

test_that("execution-boundary refactor adds no public arguments", {
  forbidden <- c("disable_kfold", "run_kfold", "posterior_only", "skip_elpd")
  expect_false(any(forbidden %in% names(formals(fit_pclv_bayes))))
})
