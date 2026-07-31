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
    subject = rep(c("A", "B"), each = 3),
    time = c(0, 2, 5, 1, 4, 8),
    xi_raw = c(.20, .30, .25, .10, .20, .28),
    xj_raw = c(.30, .20, .35, .20, .30, .22)
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
