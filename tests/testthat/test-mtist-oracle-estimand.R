source(testthat::test_path("../../benchmarks/mtist/oracle_estimand_audit.R"))

test_that("MTIST rows are targets and columns are sources", {
  A <- matrix(c(0, 2, -3, 0), 2, 2, byrow = TRUE,
              dimnames = list(c("i", "j"), c("i", "j")))
  expect_equal(A["i", "j"], 2)
  expect_equal(A["j", "i"], -3)
  expect_identical(oracle_rest_set(c("i", "j", "k"), "i", "j"), "k")
})

test_that("pair-to-rest derivative and contrast algebra are exact", {
  taxa <- c("i", "j", "k")
  x <- c(i = .2, j = .3, k = .5)
  growth <- c(i = .4, j = -.2, k = .1)
  A <- matrix(c(.1, .7, -.4, -.3, .2, .5, .6, -.8, .1), 3, 3,
              dimnames = list(taxa, taxa))
  g <- stats::setNames(oracle_glv_growth(x, growth, A), taxa)
  expected <- unname(g["i"] - g["k"])
  expect_equal(oracle_alr_derivative(x, "i", "j", growth, A, taxa), expected)
  expect_equal(oracle_contrast(x, "i", "j", A, taxa), A["i", "j"] - A["k", "j"])
  eps <- 1e-7
  x_next <- x + eps * (growth + as.numeric(A %*% x)) * x
  numeric_derivative <- unname((log(x_next["i"] / x_next["k"]) -
                           log(x["i"] / x["k"])) / eps)
  expect_equal(numeric_derivative, expected, tolerance = 1e-5)
})

test_that("projection signs are invariant to positive standardization", {
  dat <- data.frame(subject = rep(c("a", "b"), each = 5),
                    xi = seq_len(10), xj = c(2, 4, 1, 7, 3, 9, 5, 8, 6, 10))
  y <- .4 + .7 * dat$xi - .3 * dat$xj
  one <- oracle_projection(y, dat$xi, dat$xj, dat$subject)
  two <- oracle_projection(y, 10 * dat$xi + 2, 3 * dat$xj - 4, dat$subject)
  expect_true(one$identifiable && two$identifiable)
  expect_identical(one$sign, two$sign)
  expect_equal(one$rank, two$rank)
})

test_that("rank-deficient oracle projections return structured NA", {
  out <- oracle_projection(1:4, c(1, 1, 1, 1), c(2, 2, 2, 2),
                           rep("a", 4))
  expect_false(out$identifiable)
  expect_true(is.na(out$coefficient))
  expect_identical(out$reason, "rank_deficient")
})

test_that("subject boundaries produce no cross-subject predecessor", {
  y <- c(0, 1, 0, 1)
  subject <- c("a", "a", "b", "b")
  time <- c(0, 1, 0, 1)
  out <- pclvbayes:::.delta_alr_over_dt(y, time, subject)
  expect_true(is.na(out[[1]]) && is.na(out[[3]]))
  expect_equal(out[[2]], 1)
  expect_equal(out[[4]], 1)
})

test_that("oracle values are not converted to zero", {
  A <- matrix(c(0, .2, -.4, 0), 2, 2,
              dimnames = list(c("i", "j"), c("i", "j")))
  value <- oracle_contrast(c(i = .5, j = .5), "i", "j", A, c("i", "j"))
  expect_equal(value, A["i", "j"])
  expect_false(isTRUE(value == 0))
})
