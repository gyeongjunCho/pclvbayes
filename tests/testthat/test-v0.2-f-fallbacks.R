test_that("canonical CV spline failures are explicit and never substituted", {
  too_short <- pclvbayes:::.smooth_spline_robust(1:2, c(0, 1))
  expect_s3_class(too_short, "pclv_failure")
  expect_identical(too_short$stage, "spline_smoothing")
  expect_identical(too_short$reason, "insufficient_unique_times")

  expect_error(
    pclvbayes:::.smooth_spline_robust(1:4, 1:4, spline_df = 3),
    "CV smoothing is required"
  )
})

test_that("invalid time intervals are rejected rather than repaired", {
  duplicate <- pclvbayes:::.delta_alr_over_dt(c(0, 1, 2), c(0, 1, 1), rep("A", 3))
  nonfinite <- pclvbayes:::.delta_alr_over_dt(c(0, 1), c(0, Inf), rep("A", 2))
  expect_s3_class(duplicate, "pclv_failure")
  expect_identical(duplicate$reason, "non_positive_dt")
  expect_s3_class(nonfinite, "pclv_failure")
  expect_identical(nonfinite$reason, "non_finite_time")
})

test_that("missing OU persistence and invalid nu cannot select another scorer", {
  held <- data.frame(subject = "A", time = c(0, 1), y = 0, xi = 0, xj = 0)
  base <- data.frame(r0 = 0, a_ii = 0, a_ij = 0, sigma = 1, sd_ou = 1, nu = 5)
  missing <- pclvbayes:::.proj_loglik_subject(base, held)
  expect_s3_class(missing, "pclv_failure")
  expect_identical(missing$reason, "missing_ou_persistence_draws")
  base$phi <- .8
  base$nu <- NA_real_
  invalid <- pclvbayes:::.proj_loglik_subject(base, held)
  expect_s3_class(invalid, "pclv_failure")
  expect_identical(invalid$reason, "invalid_nu_draws")
})

test_that("malformed initialization stops and retry exhaustion exposes provenance", {
  expect_error(
    pclvbayes:::.sample_with_retry(NULL, list(init = "repair-me"), list(), max_retries = 0),
    "malformed init"
  )

  local_mocked_bindings(
    .call_sample_silently = function(...) stop("synthetic failure"),
    .package = "pclvbayes"
  )
  out <- pclvbayes:::.sample_with_retry(
    NULL,
    list(seed = 9, chains = 1, parallel_chains = 1, iter_warmup = 10,
         adapt_delta = .8, max_treedepth = 10, metric = "diag_e", init = .2),
    list(), max_retries = 1, silent_sampler = TRUE
  )
  expect_s3_class(out$failure, "pclv_failure")
  expect_length(out$retry_history, 2)
  expect_equal(vapply(out$retry_history, `[[`, integer(1), "attempt"), 1:2)
  expect_identical(out$failure$details$attempt_history, out$retry_history)
  expect_equal(vapply(out$retry_history, `[[`, numeric(1), "seed"), c(9, 10))
  expect_equal(vapply(out$retry_history, `[[`, character(1), "initialization_method"),
               c("scalar", "scalar"))
  expect_true(all(vapply(out$retry_history, `[[`, character(1), "status") == "failed"))
})

test_that("missing compiled worker model is an explicit loading failure", {
  out <- pclvbayes:::.pclv_failure("model_loading", "compiled_model_unavailable",
                                   list(exe_file = NULL))
  expect_false(out$ok)
  expect_identical(out$stage, "model_loading")
  expect_identical(out$reason, "compiled_model_unavailable")
})
