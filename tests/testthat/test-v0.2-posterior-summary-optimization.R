summary_draw_fixture <- function(chains = 2L, draws = 200L) {
  set.seed(812)
  out <- do.call(rbind, lapply(seq_len(chains), function(chain) data.frame(
    a_ij = rnorm(draws, -0.3, 0.04), a_ii = rnorm(draws, -0.5, 0.05),
    r0 = rnorm(draws), sigma = abs(rnorm(draws, .08, .005)),
    sd_ou = abs(rnorm(draws, .12, .008)), phi = pmin(pmax(rnorm(draws, .8, .01), .01), .99),
    lambda = abs(rnorm(draws, .1, .005)), nu = 2 + exp(rnorm(draws, log(3), .1)),
    log_nu_minus_two = rnorm(draws, log(3), .1), .chain = chain, .iteration = seq_len(draws),
    .draw = (chain - 1L) * draws + seq_len(draws)
  )))
  posterior::as_draws_df(out)
}

summary_diag_fixture <- function(n_draws = 400L) list(
  worst_rhat = 1.001, min_ess_bulk = 500, min_ess_tail = 480,
  n_divergent = 0L, n_treedepth_hit = 0L, ebfmi_min = .7,
  ebfmi_med = .75, ebfmi_chain = c(.7, .8), n_draws = n_draws
)

legacy_scientific_summary <- function(d) list(
  a_mean = mean(d$a_ij), a_median = median(d$a_ij),
  a_q = unname(quantile(d$a_ij, c(.025, .975))),
  self_mean = mean(d$a_ii), self_median = median(d$a_ii),
  p_sign2 = pclvbayes:::.clip01(2 * pmin(mean(d$a_ij > 0), mean(d$a_ij < 0))),
  lfsr = pclvbayes:::.lfsr_from_two_sided(
    pclvbayes:::.clip01(2 * pmin(mean(d$a_ij > 0), mean(d$a_ij < 0))))
)

test_that("one retained bundle preserves coefficient, sign, and LFSR conventions", {
  d <- summary_draw_fixture()
  old <- legacy_scientific_summary(d)
  bundle <- pclvbayes:::.build_posterior_summary_bundle(d, summary_diag_fixture())
  expect_equal(bundle$coefficients$interaction[["mean"]], old$a_mean)
  expect_equal(median(d$a_ij), old$a_median)
  expect_equal(unname(bundle$coefficients$interaction[c("q025", "q975")]), old$a_q)
  expect_equal(bundle$coefficients$self[["mean"]], old$self_mean)
  expect_equal(bundle$sign$interaction_p_two, old$p_sign2)
  expect_equal(bundle$sign$interaction_lfsr, old$lfsr)
})

test_that("bundle reuses unchanged nu, residual, and chain classification summaries", {
  d <- summary_draw_fixture()
  diag <- summary_diag_fixture()
  old_chain <- pclvbayes:::.classify_chain_diagnostics(d, diag)
  old_nu <- pclvbayes:::.summarise_nu_draws(d)
  bundle <- pclvbayes:::.build_posterior_summary_bundle(d, diag)
  expect_identical(bundle$chain, old_chain)
  expect_identical(bundle$nu, old_nu)
  expect_identical(bundle$chain$diagnostic_class, "converged")
  expect_identical(bundle$chain$chain_sign_agreement, old_chain$chain_sign_agreement)
  expect_identical(bundle$chain$chain_residual_summary, old_chain$chain_residual_summary)
  expect_identical(bundle$chain$indeterminate_reason, old_chain$indeterminate_reason)
})

test_that("one-chain bundle does not fabricate multi-chain agreement", {
  d <- summary_draw_fixture(chains = 1L)
  diag <- summary_diag_fixture(200L); diag$ebfmi_chain <- .7
  bundle <- pclvbayes:::.build_posterior_summary_bundle(d, diag)
  expect_true(is.na(bundle$chain$chain_sign_agreement))
  expect_true(is.na(bundle$chain$residual_regime_disagreement))
})

test_that("attempt sampler diagnostics avoid posterior draw materialization", {
  sampler <- posterior::as_draws_array(array(
    c(rep(0, 400), rep(5, 400), rnorm(400)), dim = c(200, 2, 3),
    dimnames = list(NULL, NULL, c("divergent__", "treedepth__", "energy__"))))
  posterior_calls <- 0L; sampler_calls <- 0L
  fit <- list(
    draws = function(...) { posterior_calls <<- posterior_calls + 1L; stop("posterior draws must not load") },
    sampler_diagnostics = function() { sampler_calls <<- sampler_calls + 1L; sampler }
  )
  diag <- pclvbayes:::.summarise_sampler_diag(fit, max_treedepth = 14L)
  expect_identical(posterior_calls, 0L)
  expect_identical(sampler_calls, 1L)
  expect_equal(diag$n_divergent, 0L)
  expect_equal(diag$n_treedepth_hit, 0L)
})

test_that("adding retained convergence diagnostics leaves sampler decisions unchanged", {
  d <- summary_draw_fixture()
  attempt <- summary_diag_fixture()
  retained <- pclvbayes:::.add_convergence_diag(attempt, d)
  fields <- c("n_divergent", "n_treedepth_hit", "ebfmi_min", "ebfmi_med", "n_draws")
  expect_identical(retained[fields], attempt[fields])
  expected <- posterior::summarise_draws(posterior::subset_draws(
    d, variable = intersect(c("a_ij", "a_ii", "r0", "sigma", "sd_ou", "phi", "nu", "log_nu_minus_two"),
                            posterior::variables(d))))
  expect_equal(retained$worst_rhat, max(expected$rhat, na.rm = TRUE))
  expect_equal(retained$min_ess_bulk, min(expected$ess_bulk, na.rm = TRUE))
  expect_equal(retained$min_ess_tail, min(expected$ess_tail, na.rm = TRUE))
})

test_that("posterior evidence remains immutable across downstream weight calculation", {
  d <- summary_draw_fixture()
  bundle <- pclvbayes:::.build_posterior_summary_bundle(d, summary_diag_fixture())
  evidence_before <- bundle[c("coefficients", "sign", "chain")]
  pointwise <- data.frame(from = rep(c("a", "b"), each = 3), to = "c",
                          subject = rep(1:3, 2), elpd = c(-1, -2, -1, -2, -1, -2),
                          n_test = 1L)
  invisible(pclvbayes:::.compute_pseudobma_weights_cross(
    list(elpd_pointwise_cross = pointwise), data.frame(from = c("a", "b"), to = "c"),
    min_models = 2L, min_subjects = 1L))
  expect_identical(bundle[c("coefficients", "sign", "chain")], evidence_before)
})

test_that("public API remains unchanged by private posterior summary reuse", {
  expect_false(any(c("posterior_bundle", "summary_cache", "profiling") %in%
                     names(formals(fit_pclv_bayes))))
})
