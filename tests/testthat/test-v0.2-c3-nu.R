test_that("Stan estimates shared Student-t nu with the specified prior", {
  stan <- paste(readLines(test_path("..", "..", "inst", "stan", "pclv.stan"), warn = FALSE), collapse = "\n")
  expect_match(stan, "real log_nu_minus_two;", fixed = TRUE)
  expect_match(stan, "real<lower=2> nu = 2 + exp(log_nu_minus_two);", fixed = TRUE)
  expect_match(stan, "log_nu_minus_two ~ normal(log(3), 0.75);", fixed = TRUE)
  expect_match(stan, "y ~ student_t(nu, mu + e, sigma);", fixed = TRUE)
  expect_match(stan, "student_t_lpdf(y[n] | nu", fixed = TRUE)
  expect_false(grepl("nu_fixed|use_student_t", stan))
  expect_equal(2 + exp(log(3)), 5)
  expect_true(all(2 + exp(c(-20, -1, 0, 10)) > 2))
})

test_that("posterior nu extraction returns finite summaries", {
  draws <- data.frame(nu = c(2.5, 4, 6, 20))
  out <- pclvbayes:::.summarise_nu_draws(draws)
  expect_named(out, c("nu_mean", "nu_median", "nu_q05", "nu_q95"))
  expect_true(all(is.finite(unlist(out))))
  expect_equal(out$nu_mean, mean(draws$nu))
  expect_equal(out$nu_median, median(draws$nu))
})

test_that("posterior extraction rejects missing and invalid nu draws", {
  missing <- pclvbayes:::.summarise_nu_draws(data.frame(a_ij = 1:2))
  invalid <- pclvbayes:::.summarise_nu_draws(data.frame(nu = c(2, Inf)))
  expect_s3_class(missing, "pclv_failure")
  expect_identical(missing$stage, "posterior_extraction")
  expect_identical(missing$reason, "missing_nu_draws")
  expect_s3_class(invalid, "pclv_failure")
  expect_identical(invalid$reason, "invalid_nu_draws")
})

test_that("Kalman ELPD requires and uses draw-specific nu", {
  base <- data.frame(r0 = c(0, 0), a_ii = c(0, 0), a_ij = c(0, 0),
    sigma = c(0.4, 0.4), sd_ou = c(0.2, 0.2), lambda = c(0.7, 0.7), tau_r = c(0, 0))
  held <- data.frame(subject = c("A", "A"), time = c(0, 3),
    y = c(0.2, 1.5), xi = c(0, 0), xj = c(0, 0))
  missing <- pclvbayes:::.proj_loglik_subject(base, held)
  expect_s3_class(missing, "pclv_failure")
  expect_identical(missing$stage, "elpd_scoring")
  expect_identical(missing$reason, "missing_nu_draws")
  bad <- base; bad$nu <- c(2, NA_real_)
  invalid <- pclvbayes:::.proj_loglik_subject(bad, held)
  expect_s3_class(invalid, "pclv_failure")
  expect_identical(invalid$reason, "invalid_nu_draws")
  varied <- base; varied$nu <- c(2.5, 50)
  scored <- pclvbayes:::.proj_loglik_subject(varied, held)
  expect_true(all(is.finite(scored$full)))
  expect_false(isTRUE(all.equal(scored$full[1, ], scored$full[2, ])))
})

test_that("retry and K-fold APIs cannot override the nu model", {
  expect_false(any(c("nu_fix", "max_nu_cap") %in% names(formals(pclvbayes:::.sample_with_retry))))
  expect_false(any(c("nu_fixed_override", "nu_fixed_kfold") %in% names(formals(pclvbayes:::.fold_fit_and_score))))
  expect_false("nu_fixed_kfold" %in% names(formals(pclvbayes:::.repkfold_eval)))
})
