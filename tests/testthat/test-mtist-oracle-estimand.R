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

oracle_transition_fixture <- function(subject = rep("a", 5), time = 0:4) {
  taxa <- c("species_4", "species_7", "rest")
  abundance <- cbind(
    species_4 = c(2, 3, 5, 8, 13),
    species_7 = c(3, 4, 6, 7, 9),
    rest = c(5, 6, 7, 9, 11)
  )
  list(
    abundance = abundance,
    metadata = data.frame(subject = subject, time = time),
    taxa = taxa
  )
}

test_that("raw oracle retains T minus one adjacent transitions including the final interval", {
  fixture <- oracle_transition_fixture()
  out <- oracle_raw_design(
    fixture$abundance, fixture$metadata,
    target = "species_4", source = "species_7", taxa = fixture$taxa
  )
  relative <- fixture$abundance / rowSums(fixture$abundance)
  zi <- log(relative[, "species_4"] / relative[, "rest"])
  zj <- log(relative[, "species_7"] / relative[, "rest"])

  expect_equal(nrow(out), 4L)
  expect_equal(out$predecessor_row, 1:4)
  expect_equal(out$outcome_row, 2:5)
  expect_equal(out$predecessor_time, 0:3)
  expect_equal(out$time, 1:4)
  expect_equal(out$dt, rep(1, 4))
  expect_equal(out$xi_unscaled, zi[1:4])
  expect_equal(out$xj_unscaled, zj[1:4])
  expect_equal(out$y, diff(zi))
  expect_identical(out$outcome_row[[4L]], 5L)
})

test_that("outcome-aligned derivative uses outcome rows", {
  fixture <- oracle_transition_fixture()
  derivative <- c(10, 20, 30, 40, 50)
  out <- oracle_raw_design(
    fixture$abundance, fixture$metadata,
    target = "species_4", source = "species_7", taxa = fixture$taxa,
    derivative = derivative
  )
  expect_equal(out$y, derivative[2:5])
  expect_equal(out$predecessor_row, 1:4)
  expect_equal(out$outcome_row, 2:5)
})

test_that("invalid observations remove adjacent intervals without bridging", {
  fixture <- oracle_transition_fixture()
  fixture$abundance[3, "species_7"] <- NA_real_
  out <- oracle_raw_design(
    fixture$abundance, fixture$metadata,
    target = "species_4", source = "species_7", taxa = fixture$taxa
  )
  expect_equal(out$predecessor_row, c(1L, 4L))
  expect_equal(out$outcome_row, c(2L, 5L))
  expect_false(any(out$predecessor_row == 2L & out$outcome_row == 4L))
})

test_that("raw oracle never bridges subjects", {
  one <- oracle_transition_fixture()
  abundance <- rbind(one$abundance[1:3, ], one$abundance[1:3, ])
  metadata <- data.frame(
    subject = c(rep("a", 3), rep("b", 3)),
    time = rep(0:2, 2)
  )
  out <- oracle_raw_design(
    abundance, metadata,
    target = "species_4", source = "species_7", taxa = one$taxa
  )
  expect_equal(as.integer(table(out$subject)), c(2L, 2L))
  expect_equal(names(table(out$subject)), c("a", "b"))
  expect_false(any(out$predecessor_row == 3L & out$outcome_row == 4L))
})

test_that("raw oracle rejects invalid adjacent dt", {
  fixture <- oracle_transition_fixture(time = c(0, 1, 1, 3, 4))
  expect_error(
    oracle_raw_design(
      fixture$abundance, fixture$metadata,
      target = "species_4", source = "species_7", taxa = fixture$taxa
    ),
    "strictly positive dt"
  )
  fixture$metadata$time[[3L]] <- NA_real_
  expect_error(
    oracle_raw_design(
      fixture$abundance, fixture$metadata,
      target = "species_4", source = "species_7", taxa = fixture$taxa
    ),
    "complete and finite"
  )
})

test_that("projection activates subject adjustment for character and factor subjects", {
  subject <- rep(c("a", "b"), each = 5)
  xi <- rep(c(-2, -1, 0, 1, 2), 2)
  xj <- c(2, -1, 1, -2, 0, -1, 2, -2, 1, 0)
  subject_shift <- ifelse(subject == "b", 20, 0)
  y <- 3 + subject_shift + .4 * xi - .7 * xj

  character_result <- oracle_projection(y, xi, xj, subject)
  factor_result <- oracle_projection(y, xi, xj, factor(subject))
  expect_true(character_result$identifiable)
  expect_equal(character_result$coefficient, -.7, tolerance = 1e-12)
  expect_equal(character_result$columns, 4L)
  expect_equal(factor_result, character_result)
})

test_that("single-subject projection retains the original branch", {
  xi <- c(-2, -1, 0, 1, 2)
  xj <- c(2, -1, 1, -2, 0)
  y <- 3 + .4 * xi - .7 * xj
  character_result <- oracle_projection(y, xi, xj, rep("a", 5))
  factor_result <- oracle_projection(y, xi, xj, factor(rep("a", 5)))
  expect_equal(character_result$coefficient, -.7, tolerance = 1e-12)
  expect_equal(character_result$columns, 3L)
  expect_equal(factor_result, character_result)
})

test_that("species_7 to species_4 focused path uses corrected transitions and projection", {
  fixture <- oracle_transition_fixture(
    subject = rep(c("a", "b"), each = 5), time = rep(0:4, 2)
  )
  fixture$abundance <- rbind(fixture$abundance, fixture$abundance * 1.1)
  design <- oracle_raw_design(
    fixture$abundance, fixture$metadata,
    target = "species_4", source = "species_7", taxa = fixture$taxa
  )
  projection <- oracle_projection(
    design$y, design$xi, design$xj, design$subject
  )
  expect_equal(nrow(design), 8L)
  expect_equal(as.integer(table(design$subject)), c(4L, 4L))
  expect_equal(names(table(design$subject)), c("a", "b"))
  expect_equal(max(design$outcome_row[design$subject == "a"]), 5L)
  expect_equal(max(design$outcome_row[design$subject == "b"]), 10L)
  expect_true(is.list(projection))
  expect_true(all(c("coefficient", "sign", "rank", "condition") %in%
                    names(projection)))
})

test_that("oracle values are not converted to zero", {
  A <- matrix(c(0, .2, -.4, 0), 2, 2,
              dimnames = list(c("i", "j"), c("i", "j")))
  value <- oracle_contrast(c(i = .5, j = .5), "i", "j", A, c("i", "j"))
  expect_equal(value, A["i", "j"])
  expect_false(isTRUE(value == 0))
})


test_that("empirical susceptibility uses the specified axis formula", {
  out <- oracle_empirical_susceptibility(
    canonical_sign = 1L,
    axis_classifications = list(a = c("same", "opposite", "ambiguous", "unavailable"),
                                b = c("same", "same")))
  expect_equal(out$axis$r_g[out$axis$axis == "a"], 0.5)
  expect_equal(out$axis$n_valid_comparisons[out$axis$axis == "a"], 3L)
  expect_equal(out$direction$empirical_sign_reversal_susceptibility, (0.5 + 0) / 2)
  expect_equal(out$direction$number_of_valid_comparisons, 5L)
})

test_that("missing and ambiguous comparisons are handled structurally", {
  out <- oracle_empirical_susceptibility(
    canonical_sign = 1L,
    axis_classifications = list(empty = "unavailable", amb = c("ambiguous", "unavailable")))
  expect_true(is.na(out$axis$r_g[out$axis$axis == "empty"]))
  expect_equal(out$axis$r_g[out$axis$axis == "amb"], 0.5)
  expect_equal(out$direction$number_of_valid_axes, 1L)
  missing <- oracle_empirical_susceptibility(NA_integer_, axis_classifications = list(a = "same"))
  expect_identical(missing$direction$score_status, "missing")
  expect_identical(missing$direction$score_reason, "missing_canonical_posterior_sign")
})

test_that("worst-axis ties are deterministic and absolute A is external", {
  out <- oracle_empirical_susceptibility(
    canonical_sign = -1L,
    axis_classifications = list(preprocessing = c("opposite"), subject = c("opposite"),
                                denominator = c("unavailable")))
  expect_equal(out$direction$worst_axis, "preprocessing;subject")
  expect_equal(out$direction$worst_axis_sign_reversal_susceptibility, 1)
  expect_equal(out$direction$number_of_valid_axes, 2L)
})
