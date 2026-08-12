test_that("expected scientific failures use one structured schema", {
  failure <- pclvbayes:::.select_smoothed_pair_rows(
    matrix(c(.2, .3), nrow = 1, dimnames = list("x", c("s1", "s2"))),
    data.frame(Sample = c("s1", "s2"), subject = "A", time = c(0, 1)),
    j = "x", i = "x", min_pairs = 2
  )
  expect_s3_class(failure, "pclv_failure")
  expect_false(failure$ok)
  expect_identical(failure$stage, "preprocessing")
  expect_identical(failure$reason, "insufficient_rows")
  expect_type(failure$details, "list")

  split_failure <- pclvbayes:::.build_train_test(
    data.frame(subject = "A", time = 1, y = 1, xi = 1, xj = 1),
    train_subjects = "A", test_subjects = "B", min_pairs = 1
  )
  expect_s3_class(split_failure, "pclv_failure")
  expect_false(split_failure$ok)
  expect_identical(split_failure$reason, "empty_train_or_test")
})

test_that("retry exhaustion returns a structured computational failure", {
  local_mocked_bindings(
    .call_sample_silently = function(...) stop("synthetic CmdStan failure"),
    .package = "pclvbayes"
  )
  out <- pclvbayes:::.sample_with_retry(
    mod = NULL,
    base_args = list(seed = 1, chains = 1, parallel_chains = 1,
                     iter_warmup = 1, adapt_delta = .8, max_treedepth = 10,
                     metric = "diag_e", init = .2),
    stan_list = list(), max_retries = 1, silent_sampler = TRUE
  )
  expect_true(out$fit_failed)
  expect_equal(out$n_retries, 1L)
  expect_s3_class(out$failure, "pclv_failure")
  expect_identical(out$failure$stage, "sampling")
  expect_identical(out$failure$reason, "cmdstan_execution_failed")
  expect_match(out$failure$details$message, "synthetic CmdStan failure")
})

test_that("programmer invariants still stop immediately", {
  draws <- data.frame(r0 = 0, a_ii = 0, a_ij = 0, sigma = 1, sd_r0 = 0, nu = 5)
  held <- data.frame(subject = "A", time = 0, y = 0, xi = 0, xj = 0)
  missing_ou <- pclvbayes:::.proj_loglik_subject(draws, held)
  expect_s3_class(missing_ou, "pclv_failure")
  expect_identical(missing_ou$reason, "missing_ou_scale_draws")
  expect_error(
    pclvbayes:::.make_repkfold_splits(c("A", "B"), K = 3, R = 1),
    "K > #subjects"
  )
})

test_that("successful directions retain the existing pair summaries", {
  make_success <- function(value) {
    out <- pclvbayes:::.failed_direction_result(
      pclvbayes:::.pclv_failure("test", "unused", list())
    )
    out$ok <- TRUE; out$failure <- NULL
    out$n_pairs <- 7L; out$a_mean <- value
    out
  }
  run_one <- function(target, partner, ctx, seed_override, progress_local) {
    make_success(if (target == "a") .25 else -.4)
  }
  result <- pclvbayes:::.run_pair(
    1, 2, kfold_K = 2, kfold_R = 1, taxa_vec = c("a", "b"),
    .run_one = run_one, ctx = list(), progress = "none",
    mute_logs = TRUE, seed_base = 1
  )
  expect_true(result$direction_ok_ij)
  expect_true(result$direction_ok_ji)
  expect_null(result$failure_ij[[1]])
  expect_null(result$failure_ji[[1]])
  expect_equal(result$a_ij_mean, .25)
  expect_equal(result$a_ji_mean, -.4)
})

test_that("canonical expected-failure helpers never return bare NULL", {
  preprocessing <- pclvbayes:::.make_pair_inputs_glv(
    data.frame(subject = "A", time = 0, xi_raw = .2, xj_raw = .3),
    zero_mode_alr = "fixed", smooth_scale = "logra", nz_partner_min_frac = 0
  )
  expect_s3_class(preprocessing, "pclv_failure")
  expect_false(is.null(preprocessing))
})
