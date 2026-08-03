test_that("Student-t Kalman scoring reduces exactly to Student-t without latent-state uncertainty", {
  draws <- data.frame(
    r0 = 0,
    a_ii = 0,
    a_ij = 0,
    sigma = 1.25,
    nu = 5,
    sd_ou = 0,
    lambda = 1
  )
  held_out <- data.frame(
    y = 3,
    xi = 0,
    xj = 0,
    subject = "s1",
    time = 0
  )

  scored <- pclvbayes:::.proj_loglik_subject(draws, held_out)
  expected <- stats::dt(3 / 1.25, df = 5, log = TRUE) - log(1.25)

  expect_false(inherits(scored, "pclv_failure"))
  expect_equal(unname(scored$full[1, 1]), expected, tolerance = 1e-10)
  expect_identical(scored$method, "student-t-scale-mixture-kalman-ou-q16")
  expect_identical(scored$quadrature_nodes, 16L)
  expect_identical(scored$state_approximation, "gaussian-moment-collapse")
})

test_that("Student-t predictive scoring retains heavier tails than the old Gaussian moment match", {
  nu <- 3
  sigma <- 1
  y <- 8
  draws <- data.frame(
    r0 = 0,
    a_ii = 0,
    a_ij = 0,
    sigma = sigma,
    nu = nu,
    sd_ou = 0,
    lambda = 1
  )
  held_out <- data.frame(
    y = y,
    xi = 0,
    xj = 0,
    subject = "s1",
    time = 0
  )

  scored <- pclvbayes:::.proj_loglik_subject(draws, held_out)
  gaussian_moment_match <- stats::dnorm(
    y,
    mean = 0,
    sd = sqrt(sigma^2 * nu / (nu - 2)),
    log = TRUE
  )

  expect_gt(scored$full[1, 1], gaussian_moment_match)
})

test_that("Student-t scale-mixture Kalman scoring is deterministic", {
  draws <- data.frame(
    r0 = c(0.1, -0.2),
    a_ii = c(-0.5, -0.3),
    a_ij = c(0.2, 0.4),
    sigma = c(0.4, 0.7),
    nu = c(4, 12),
    sd_ou = c(0.3, 0.5),
    lambda = c(0.8, 1.1),
    tau_r = c(0.2, 0.1)
  )
  held_out <- data.frame(
    y = c(0.2, -0.1, 0.5),
    xi = c(-0.5, 0.1, 0.7),
    xj = c(0.3, -0.2, 0.4),
    subject = "s1",
    time = c(0, 1, 3)
  )

  first <- pclvbayes:::.proj_loglik_subject(draws, held_out)
  second <- pclvbayes:::.proj_loglik_subject(draws, held_out)

  expect_equal(first$full, second$full, tolerance = 0)
  expect_true(all(is.finite(first$full)))
})

test_that("Student-t Kalman scoring approaches Gaussian OU scoring as nu grows", {
  draws <- data.frame(
    r0 = 0,
    a_ii = 0,
    a_ij = 0,
    sigma = 1,
    nu = 1e6,
    sd_ou = 0.5,
    lambda = 1
  )
  held_out <- data.frame(
    y = 1.5,
    xi = 0,
    xj = 0,
    subject = "s1",
    time = 0
  )

  scored <- pclvbayes:::.proj_loglik_subject(draws, held_out)
  gaussian_limit <- stats::dnorm(
    1.5,
    mean = 0,
    sd = sqrt(1^2 + 0.5^2),
    log = TRUE
  )

  expect_equal(unname(scored$full[1, 1]), gaussian_limit, tolerance = 1e-5)
})

test_that("invalid random-intercept draws fail closed", {
  draws <- data.frame(
    r0 = 0,
    a_ii = 0,
    a_ij = 0,
    sigma = 1,
    nu = 5,
    sd_ou = 0.3,
    lambda = 1,
    tau_r = -0.1
  )
  held_out <- data.frame(
    y = 0,
    xi = 0,
    xj = 0,
    subject = "s1",
    time = 0
  )

  scored <- pclvbayes:::.proj_loglik_subject(draws, held_out)

  expect_s3_class(scored, "pclv_failure")
  expect_identical(scored$stage, "elpd_scoring")
  expect_identical(scored$reason, "invalid_ou_draws")
})
