# Regression/architecture tests for the shared longitudinal preprocessing layer.
#
# These tests intentionally separate two contracts:
#   (1) observed-support eligibility: sparse trajectories are excluded BEFORE
#       zero replacement; and
#   (2) triplet ALR numerics: for eligible trajectories, the current pcLV
#       replacement/rest-floor/closure/cap math remains equivalent.
#
# Future OU and fit_pclv_bayes() migrations should add tests here before deleting
# their local preprocessing implementations. In particular, pcLV support must be
# assessed from observed pre-spline RA, while its weak spline remains a separate
# model-specific temporal stage.

testthat::test_that("shared canonical triplet transform matches current pcLV row logic", {
  xi <- c(0.00, 0.20, 0.05, 0.00, 0.70, 0.01)
  xj <- c(0.10, 0.00, 0.20, 0.00, 0.20, 0.89)
  xr <- pmax(0, 1 - xi - xj)
  n <- length(xi)

  df <- data.frame(
    xi_raw = xi,
    xj_raw = xj,
    rest_raw = xr
  )

  old_trip <- t(vapply(
    seq_len(n),
    function(k) {
      .make_triplet_row(
        i = k,
        df = df,
        subj_minpos = rep(NA_real_, n),
        lib = rep(NA_real_, n),
        zero_mode_alr = "minpos_time",
        minpos_alpha = 0.5,
        eps_fixed = 1e-6,
        lib_eps_c = 0.65,
        rest_floor_frac = 1.0,
        minpos_base = "ij"
      )
    },
    numeric(4)
  ))
  colnames(old_trip) <- c("xi", "xj", "xr", "eps_star")
  old_alr <- .pair_alr(old_trip, alr_cap = 12)

  shared <- .triplet_alr_transform(
    xi_raw = xi,
    xj_raw = xj,
    rest_raw = xr,
    zero_mode_alr = "minpos_time",
    minpos_alpha = 0.5,
    minpos_base = "ij",
    eps_fixed = 1e-6,
    lib_eps_c = 0.65,
    rest_floor_frac = 1.0,
    alr_cap_mode = "fixed",
    alr_cap = 12
  )

  testthat::expect_equal(shared$xi, old_trip[, "xi"], tolerance = 1e-14)
  testthat::expect_equal(shared$xj, old_trip[, "xj"], tolerance = 1e-14)
  testthat::expect_equal(shared$xr, old_trip[, "xr"], tolerance = 1e-14)
  testthat::expect_equal(shared$eps_star, old_trip[, "eps_star"], tolerance = 1e-14)
  testthat::expect_equal(shared$alr_i, old_alr$i, tolerance = 1e-14)
  testthat::expect_equal(shared$alr_j, old_alr$j, tolerance = 1e-14)
  testthat::expect_equal(shared$xi + shared$xj + shared$xr, rep(1, n), tolerance = 1e-14)
  testthat::expect_lte(max(abs(shared$alr_i)), 12)
  testthat::expect_lte(max(abs(shared$alr_j)), 12)
})


testthat::test_that("matrix backend preserves the canonical element-wise transform", {
  xi <- matrix(c(0, 0.2, 0.1, 0), nrow = 2)
  xj <- matrix(c(0.1, 0, 0.2, 0), nrow = 2)
  xr <- 1 - xi - xj

  first_positive <- xi
  first_positive[first_positive <= 0] <- Inf
  second_positive <- xj
  second_positive[second_positive <= 0] <- Inf
  eps_t <- 0.5 * pmin(first_positive, second_positive)
  eps_t[!is.finite(eps_t)] <- 1e-6

  out <- .triplet_alr_apply(
    xi_raw = xi,
    xj_raw = xj,
    rest_raw = xr,
    eps_t = eps_t,
    rest_floor_frac = 1,
    alr_cap_mode = "fixed",
    alr_cap = 12
  )

  testthat::expect_identical(dim(out$alr_i), dim(xi))
  testthat::expect_identical(dim(out$alr_j), dim(xj))
  testthat::expect_equal(out$xi + out$xj + out$xr, matrix(1, nrow = 2, ncol = 2), tolerance = 1e-14)
})


testthat::test_that("effective n handles four observations and caps excessive lag requests", {
  x <- c(1, 3, 2, 4)
  y <- c(4, 2, 3, 1)

  ne_nw <- .eff_n(
    x, y,
    min_n = 4L,
    effn_method = "nw",
    nw_bw = 999L
  )
  ne_bartlett <- .eff_n(
    x, y,
    min_n = 4L,
    effn_method = "bartlett",
    L = 999L
  )

  testthat::expect_true(is.finite(ne_nw))
  testthat::expect_true(is.finite(ne_bartlett))
  testthat::expect_gte(ne_nw, 4)
  testthat::expect_lte(ne_nw, 4)
  testthat::expect_gte(ne_bartlett, 4)
  testthat::expect_lte(ne_bartlett, 4)
})


testthat::test_that("time aggregation reports the actual post-aggregation length", {
  mat <- cbind(
    x = c(1, 3, 5, 7, 9),
    y = c(2, 4, 6, 8, 10)
  )
  tt <- c(0, 0, 1, 2, 2)

  out <- .aggregate_by_time_series(mat, tt, mode = "mean")

  testthat::expect_equal(out$time, c(0, 1, 2))
  testthat::expect_equal(nrow(out$mat), 3L)
  testthat::expect_equal(out$mat[, "x"], c(2, 5, 8))
  testthat::expect_equal(out$mat[, "y"], c(3, 6, 9))
})


testthat::test_that("cor_meta_resid keeps exact subjectwise correlations without infinite meta effects", {
  testthat::skip_if_not_installed("phyloseq")
  testthat::skip_if_not_installed("metafor")

  subject <- rep(c("s1", "s2"), each = 5)
  time <- rep(seq_len(5), 2)
  sample_id <- paste0("s", seq_len(10))

  a <- rep(c(10, 20, 30, 40, 50), 2)
  b <- a
  c_rest <- rep(100, 10)
  otu <- rbind(A = a, B = b, C = c_rest)
  colnames(otu) <- sample_id

  md <- data.frame(subject = subject, time = time, row.names = sample_id)
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(md)
  )

  out <- cor_meta_resid(
    ps,
    subject_col = "subject",
    time_col = "time",
    taxa_vec = c("A", "B"),
    method = "pearson",
    transform = "raw",
    min_n = 5L,
    min_k = 2L,
    use_knha = FALSE,
    acf_correction = FALSE,
    effn_aggregate_by_time = "none"
  )

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$k, 2L)
  testthat::expect_true(is.finite(out$r_pooled))
  testthat::expect_gt(out$r_pooled, 0.999999)
})


testthat::test_that("cor_meta_resid applies min_n after time aggregation", {
  testthat::skip_if_not_installed("phyloseq")
  testthat::skip_if_not_installed("metafor")

  subject <- rep(c("s1", "s2"), each = 6)
  time <- rep(c(0, 0, 1, 1, 2, 2), 2)
  sample_id <- paste0("d", seq_len(12))

  otu <- rbind(
    A = rep(c(10, 11, 20, 21, 30, 31), 2),
    B = rep(c(30, 31, 20, 21, 10, 11), 2),
    C = rep(100, 12)
  )
  colnames(otu) <- sample_id

  md <- data.frame(subject = subject, time = time, row.names = sample_id)
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(md)
  )

  out <- cor_meta_resid(
    ps,
    subject_col = "subject",
    time_col = "time",
    taxa_vec = c("A", "B"),
    method = "pearson",
    transform = "raw",
    min_n = 4L,
    min_k = 2L,
    use_knha = FALSE,
    acf_correction = FALSE,
    effn_aggregate_by_time = "mean"
  )

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$k, 0L)
  testthat::expect_true(is.na(out$r_pooled))
})


testthat::test_that("observed-support gate excludes mostly-zero trajectories before replacement", {
  half <- .pair_observed_support(
    xi = c(1, 0, 1, 0, 1, 0, 1, 0),
    xj = c(0, 1, 0, 1, 0, 1, 0, 1),
    min_positive_frac = 0.5,
    min_positive_n = 4L
  )
  testthat::expect_true(half$keep)
  testthat::expect_equal(half$positive_frac_i, 0.5)
  testthat::expect_equal(half$zero_frac_i, 0.5)

  mostly_zero <- .pair_observed_support(
    xi = c(1, 0, 0, 0, 1, 0, 0, 0),
    xj = rep(1, 8),
    min_positive_frac = 0.5,
    min_positive_n = 4L
  )
  testthat::expect_false(mostly_zero$keep)
  testthat::expect_false(mostly_zero$keep_i)

  too_few_absolute <- .pair_observed_support(
    xi = c(1, 1, 1, 0, 0),
    xj = rep(1, 5),
    min_positive_frac = 0.5,
    min_positive_n = 4L
  )
  testthat::expect_false(too_few_absolute$keep)
  testthat::expect_gt(too_few_absolute$positive_frac_i, 0.5)
})


testthat::test_that("cor_meta_resid reports support-based subject exclusion", {
  testthat::skip_if_not_installed("phyloseq")
  testthat::skip_if_not_installed("metafor")

  subject <- rep(c("s1", "s2"), each = 8)
  time <- rep(seq_len(8), 2)
  sample_id <- paste0("z", seq_len(16))

  # s1 is well observed for both taxa. s2 taxon A is positive at only 2/8
  # time points and must be excluded BEFORE ALR zero replacement.
  A <- c(10, 12, 14, 16, 18, 20, 22, 24, 10, 0, 0, 0, 12, 0, 0, 0)
  B <- c(30, 28, 26, 24, 22, 20, 18, 16, 20, 18, 16, 14, 12, 10, 8, 6)
  C <- rep(100, 16)
  otu <- rbind(A = A, B = B, C = C)
  colnames(otu) <- sample_id

  md <- data.frame(subject = subject, time = time, row.names = sample_id)
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(md)
  )

  out <- cor_meta_resid(
    ps,
    subject_col = "subject",
    time_col = "time",
    taxa_vec = c("A", "B"),
    method = "pearson",
    transform = "alr",
    min_n = 5L,
    min_k = 2L,
    use_knha = FALSE,
    acf_correction = FALSE,
    effn_aggregate_by_time = "none",
    nz_partner_min_frac = 0.5,
    nz_partner_min_n = 4L
  )

  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$k, 1L)
  testthat::expect_equal(out$support_screened_subjects, 2L)
  testthat::expect_equal(out$support_excluded_subjects, 1L)
  testthat::expect_equal(out$nz_partner_min_frac, 0.5)
  testthat::expect_equal(out$nz_partner_min_n, 4L)
  testthat::expect_true(is.na(out$r_pooled))
})

