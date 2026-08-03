diag_ok_fixture <- function(...) {
  modifyList(list(worst_rhat = 1.001, min_ess_bulk = 900, min_ess_tail = 850,
                  n_divergent = 0L, n_treedepth_hit = 0L,
                  ebfmi_min = 0.7, ebfmi_med = 0.8, n_draws = 800L), list(...))
}

chain_draws <- function(aij, sigma_shift = c(0, 0), sdou_shift = c(0, 0)) {
  stopifnot(length(aij) == 2L)
  n <- 400L
  make <- function(k) {
    z <- seq(-1, 1, length.out = n)
    data.frame(
      a_ij = aij[[k]] + 0.015 * z,
      a_ii = -0.4 + 0.02 * z,
      sigma = 0.08 + sigma_shift[[k]] + 0.002 * z,
      sd_ou = 0.12 + sdou_shift[[k]] + 0.003 * z,
      phi = 0.8 + 0.01 * z,
      lambda = 0.1 - 0.002 * z,
      nu = 5 + 0.1 * z,
      .chain = k
    )
  }
  rbind(make(1L), make(2L))
}

classify_fixture <- function(draws, diag = diag_ok_fixture())
  pclvbayes:::.classify_chain_diagnostics(draws, diag)

test_that("positive and negative multi-chain interactions converge when diagnostics agree", {
  positive <- classify_fixture(chain_draws(c(0.3, 0.305)))
  negative <- classify_fixture(chain_draws(c(-0.3, -0.305)))
  expect_identical(positive$diagnostic_class, "converged")
  expect_identical(negative$diagnostic_class, "converged")
  expect_true(positive$chain_sign_agreement)
  expect_true(negative$chain_sign_agreement)
  expect_gt(positive$pooled_sign_probability, 0.99)
  expect_equal(names(positive$chain_aij_means), c("1", "2"))
})

test_that("stable interaction and separated residual scales are reported separately", {
  result <- classify_fixture(chain_draws(c(-0.3, -0.305), sigma_shift = c(0, 0.20)))
  expect_identical(result$diagnostic_class, "interaction_stable_residual_unstable")
  expect_true(result$interaction_identifiable)
  expect_false(result$residual_identifiable)
  expect_true(result$residual_regime_disagreement)
  expect_true(all(c("chain_specific_residual_regime", "residual_scale_nonidentifiability") %in%
                    result$indeterminate_reason))
})

test_that("near-zero chain regime is interaction magnitude disagreement", {
  draws <- chain_draws(c(-0.30, -0.005))
  result <- classify_fixture(draws)
  expect_identical(result$diagnostic_class, "interaction_indeterminate")
  expect_true("interaction_magnitude_disagreement" %in% result$indeterminate_reason)
  expect_false(result$interaction_identifiable)
})

test_that("opposite chain signs are interaction indeterminate", {
  result <- classify_fixture(chain_draws(c(-0.3, 0.3)))
  expect_identical(result$diagnostic_class, "interaction_indeterminate")
  expect_true("chain_sign_disagreement" %in% result$indeterminate_reason)
  expect_false(result$chain_sign_agreement)
})

test_that("stable pooled sign cannot override high R-hat or low E-BFMI", {
  high_rhat <- classify_fixture(chain_draws(c(0.3, 0.305)), diag_ok_fixture(worst_rhat = 1.2))
  low_energy <- classify_fixture(chain_draws(c(0.3, 0.305)), diag_ok_fixture(ebfmi_min = 0.1))
  expect_identical(high_rhat$diagnostic_class, "interaction_stable_residual_unstable")
  expect_true("high_rhat" %in% high_rhat$indeterminate_reason)
  expect_identical(low_energy$diagnostic_class, "interaction_stable_residual_unstable")
  expect_true("low_ebfmi" %in% low_energy$indeterminate_reason)
})

test_that("one-chain fits do not fabricate multi-chain agreement", {
  draws <- subset(chain_draws(c(-0.3, -0.3)), .chain == 1L)
  result <- classify_fixture(draws)
  expect_identical(result$diagnostic_class, "converged")
  expect_true(is.na(result$chain_sign_agreement))
  expect_true(is.na(result$residual_regime_disagreement))
})

test_that("a diagnostic failure preserves the successful opposite direction", {
  failed <- pclvbayes:::.failed_direction_result(
    pclvbayes:::.pclv_failure("posterior_diagnostics", "interaction_indeterminate")
  )
  failed$diagnostic_failure <- list(failed$failure)
  success <- failed
  success$ok <- TRUE; success$failure <- NULL; success$diagnostic_class <- "converged"
  success$diagnostic_failure <- list(NULL); success$interaction_identifiable <- TRUE
  runner <- function(target, partner, ctx, seed_override, progress_local) {
    if (target == "a") failed else success
  }
  result <- pclvbayes:::.run_pair(1L, 2L, taxa_vec = c("a", "b"), .run_one = runner,
                                  ctx = list(), progress = "none", seed_base = 1L)
  expect_false(result$direction_ok_ij)
  expect_true(result$direction_ok_ji)
  expect_s3_class(result$failure_ij[[1]], "pclv_failure")
})

test_that("chain sign certainty threshold is explicit and boundary-sensitive", {
  expect_identical(pclvbayes:::.PCLV_CHAIN_SIGN_PROB_MIN, 0.95)
  draws <- chain_draws(c(0.3, 0.3))
  # At the documented boundary, certainty is inclusive.
  draws$a_ij[draws$.chain == 2L][1:20] <- -0.01
  result <- classify_fixture(draws)
  expect_true(result$chain_aij_positive_probabilities[["2"]] >= 0.95)
})

test_that("conservative summarization excludes indeterminate directions", {
  failure <- pclvbayes:::.pclv_failure("posterior_diagnostics", "interaction_indeterminate")
  template <- pclvbayes:::.failed_direction_result(failure)
  indeterminate <- template
  indeterminate$diagnostic_class <- "interaction_indeterminate"
  indeterminate$diagnostic_failure <- list(failure)
  indeterminate$a_mean <- -0.4; indeterminate$a_q2.5 <- -0.5; indeterminate$a_q97.5 <- -0.3
  indeterminate$p_sign2 <- 0.01
  converged <- template
  converged$ok <- TRUE; converged$failure <- NULL; converged$diagnostic_class <- "converged"
  converged$diagnostic_failure <- list(NULL); converged$interaction_identifiable <- TRUE
  converged$residual_identifiable <- TRUE
  converged$a_mean <- 0.4; converged$a_q2.5 <- 0.3; converged$a_q97.5 <- 0.5
  converged$p_sign2 <- 0.01
  good_diag <- diag_ok_fixture()
  indeterminate$diag <- good_diag; converged$diag <- good_diag
  runner <- function(target, partner, ctx, seed_override, progress_local) {
    if (target == "a") indeterminate else converged
  }
  raw <- pclvbayes:::.run_pair(1L, 2L, taxa_vec = c("a", "b"), .run_one = runner,
                               ctx = list(), progress = "none", seed_base = 1L)
  fit_result <- list(cross = pclvbayes:::.mk_cross(raw), self = pclvbayes:::.mk_self(raw), raw = raw)
  summary <- summarize_bayes_pclv(fit_result)$cross
  bad <- summary[summary$diagnostic_class == "interaction_indeterminate", ]
  good <- summary[summary$diagnostic_class == "converged", ]
  expect_false(bad$diag_ok)
  expect_true(is.na(bad$bayes_FDR))
  expect_true(good$diag_ok)
  expect_true(is.finite(good$bayes_FDR))
})
