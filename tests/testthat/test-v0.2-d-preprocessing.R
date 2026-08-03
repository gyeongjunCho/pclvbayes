test_that("canonical preprocessing invokes each scientific transform once", {
  pair_df <- data.frame(
    subject = rep(c("A", "B"), each = 4),
    time = c(0, 2, 5, 9, 1, 4, 8, 13),
    xi_raw = c(.20, .30, .25, .35, .10, .20, .28, .32),
    xj_raw = c(.30, .20, .35, .25, .20, .30, .22, .18)
  )
  calls <- new.env(parent = emptyenv())
  calls$rest <- calls$alr <- calls$rate <- 0L
  rest_impl <- pclvbayes:::.pair_to_rest_abundance
  alr_impl <- pclvbayes:::.pair_alr
  rate_impl <- pclvbayes:::.delta_alr_over_dt
  local_mocked_bindings(
    .pair_to_rest_abundance = function(...) { calls$rest <- calls$rest + 1L; rest_impl(...) },
    .pair_alr = function(...) { calls$alr <- calls$alr + 1L; alr_impl(...) },
    .delta_alr_over_dt = function(...) { calls$rate <- calls$rate + 1L; rate_impl(...) },
    .package = "pclvbayes"
  )
  canonical <- pclvbayes:::.make_pair_inputs_glv(
    pair_df, zero_mode_alr = "fixed", eps_fixed = 1e-12,
    alr_cap = 12, smooth_scale = "logra", nz_partner_min_frac = 0
  )
  expect_equal(c(rest = calls$rest, alr = calls$alr, rate = calls$rate),
               c(rest = 1L, alr = 1L, rate = 1L))
  before_split <- c(calls$rest, calls$alr, calls$rate)
  split <- pclvbayes:::.build_train_test(canonical, "A", "B", min_pairs = 1)
  expect_false(is.null(split))
  expect_equal(c(calls$rest, calls$alr, calls$rate), before_split)
  expect_equal(split$train$y, canonical$y[canonical$subject == "A"])
  expect_equal(split$test$y, canonical$y[canonical$subject == "B"])
})

test_that("matrix and retained-pair entry paths produce identical canonical outputs", {
  sm_mat <- rbind(
    target = c(.20, .30, .25, .35, .10, .20, .28, .32),
    partner = c(.30, .20, .35, .25, .20, .30, .22, .18),
    other = c(.50, .50, .40, .40, .70, .50, .50, .50)
  )
  colnames(sm_mat) <- paste0("s", seq_len(ncol(sm_mat)))
  metadata <- data.frame(
    Sample = colnames(sm_mat), subject = rep(c("A", "B"), each = 4),
    time = c(0, 2, 5, 9, 1, 4, 8, 13)
  )
  retained <- data.frame(
    subject = rep(c("A", "B"), each = 4),
    time = c(0, 2, 5, 9, 1, 4, 8, 13),
    xi_raw = c(.20, .30, .25, .35, .10, .20, .28, .32),
    xj_raw = c(.30, .20, .35, .25, .20, .30, .22, .18)
  )
  args <- list(zero_mode_alr = "fixed", eps_fixed = 1e-12,
               alr_cap = 12, smooth_scale = "logra", nz_partner_min_frac = 0)
  from_matrix <- do.call(pclvbayes:::.make_pair_inputs_glv,
    c(list(sm_mat = sm_mat, meta_df = metadata, j = "partner", i = "target",
           min_pairs = 1, min_sd = 0), args))
  from_retained <- do.call(pclvbayes:::.make_pair_inputs_glv,
    c(list(pair_df = retained), args))
  expect_equal(from_matrix, from_retained, tolerance = 1e-12)
  expect_equal(from_matrix$y, from_retained$y, tolerance = 1e-12)
  expect_equal(from_matrix$xi_unscaled, from_retained$xi_unscaled, tolerance = 1e-12)
  expect_equal(from_matrix$xj_unscaled, from_retained$xj_unscaled, tolerance = 1e-12)
})

test_that("full-community smoothing recloses before pair extraction", {
  samples <- paste0("s", 0:5)
  mat <- rbind(
    a = c(.10, .13, .17, .22, .28, .35),
    b = c(.20, .24, .21, .18, .16, .14),
    c = c(.70, .63, .62, .60, .56, .51)
  )
  colnames(mat) <- samples
  metadata <- data.frame(
    Sample = samples, subject = "subject-1", time = 0:5
  )

  full <- pclvbayes:::.precompute_spline_smoothed(mat, metadata)
  pair <- pclvbayes:::.precompute_spline_smoothed(
    mat, metadata, taxa_list = c("a", "b")
  )

  expect_false(inherits(full, "pclv_failure"))
  expect_equal(unname(colSums(full)), rep(1, ncol(full)), tolerance = 1e-12)
  expect_equal(pair, full[c("a", "b"), , drop = FALSE], tolerance = 0)
  expect_equal(1 - pair["a", ] - pair["b", ], full["c", ], tolerance = 1e-12)

  forward <- pclvbayes:::.select_smoothed_pair_rows(
    full, metadata, j = "b", i = "a", min_pairs = 1
  )
  reverse <- pclvbayes:::.select_smoothed_pair_rows(
    full, metadata, j = "a", i = "b", min_pairs = 1
  )
  expect_equal(
    1 - forward$xi_raw - forward$xj_raw,
    1 - reverse$xi_raw - reverse$xj_raw,
    tolerance = 1e-15
  )
})

test_that("all observations produce T minus one transitions without bridging gaps", {
  pair <- data.frame(
    subject = "A", time = 0:4,
    xi_raw = c(.10, .13, .18, .24, .31),
    xj_raw = c(.20, .24, .21, .17, .14)
  )
  out <- pclvbayes:::.make_pair_inputs_glv(
    pair, zero_mode_alr = "fixed", eps_fixed = 1e-12,
    alr_cap = 12, smooth_scale = "logra", nz_partner_min_frac = 0
  )
  expect_equal(nrow(out), 4L)
  expect_equal(out$time, 1:4)

  rates <- pclvbayes:::.delta_alr_over_dt(
    c(0, NA_real_, 2, 3), 0:3, rep("A", 4)
  )
  expect_true(all(is.na(rates[1:3])))
  expect_equal(rates[[4]], 1)
})

test_that("ij and triplet minimum-positive bases are distinct", {
  row <- data.frame(xi_raw = 0, xj_raw = .2, rest_raw = .01)
  args <- list(
    i = 1L, df = row, subj_minpos = NA_real_, lib = NA_real_,
    zero_mode_alr = "minpos_time", minpos_alpha = .5,
    eps_fixed = 1e-8, lib_eps_c = .65, rest_floor_frac = 1
  )
  ij <- do.call(pclvbayes:::.make_triplet_row, c(args, list(minpos_base = "ij")))
  triplet <- do.call(
    pclvbayes:::.make_triplet_row, c(args, list(minpos_base = "triplet"))
  )
  expect_equal(unname(ij[[4]]), .1 / .4, tolerance = 1e-12)
  expect_equal(unname(triplet[[4]]), .005 / .215, tolerance = 1e-12)
  expect_gt(ij[[4]], triplet[[4]])
})

test_that("partner sparsity uses original abundance", {
  pair <- data.frame(
    subject = "A", time = 0:4,
    xi_raw = c(.10, .14, .19, .25, .32),
    xj_raw = c(0, 0, 0, .10, .20)
  )
  out <- pclvbayes:::.make_pair_inputs_glv(
    pair, zero_mode_alr = "fixed", eps_fixed = 1e-4,
    alr_cap = 12, smooth_scale = "logra", nz_partner_min_frac = .5
  )
  expect_s3_class(out, "pclv_failure")
  expect_identical(out$reason, "no_valid_lagged_rows")
})
