test_that("current package identity and canonical exports are available", {
  expect_identical(
    unname(getNamespaceName(asNamespace("pclvbayes"))),
    "pclvbayes"
  )
  expect_true(all(c(
    "fit_pclv_bayes",
    "cor_meta_resid",
    "summarize_bayes_pclv"
  ) %in% getNamespaceExports("pclvbayes")))
})

test_that("cor_meta_resid keeps its screening output contract", {
  skip_if_not_installed("phyloseq")
  skip_if_not_installed("metafor")

  samples <- paste0("s", seq_len(18))
  subject <- rep(paste0("subject", 1:3), each = 6)
  time <- rep(c(0, 1, 3, 6, 10, 15), 3)
  phase <- rep(seq_len(6), 3) + rep(c(0, 2, 4), each = 6)

  counts <- rbind(
    taxon_a = 20 + 3 * phase,
    taxon_b = 75 - 2 * phase,
    taxon_c = 30 + (phase %% 4) * 5
  )
  colnames(counts) <- samples

  metadata <- data.frame(
    subject = subject,
    time = time,
    row.names = samples
  )
  physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(metadata)
  )

  screened <- cor_meta_resid(
    physeq,
    subject_col = "subject",
    time_col = "time",
    min_n = 5,
    min_k = 2,
    acf_correction = FALSE,
    return_subjectwise = TRUE
  )

  expect_named(screened, c("meta", "subjectwise"))
  expect_s3_class(screened$meta, "tbl_df")
  expect_s3_class(screened$subjectwise, "tbl_df")
  expect_true(all(c(
    "i", "j", "method", "transform", "k",
    "r_pooled", "ciL", "ciU", "pval", "qval",
    "tau2", "I2", "n_min", "n_median", "n_max"
  ) %in% names(screened$meta)))
  expect_true(all(c(
    "i", "j", "subject", "r", "n", "n_raw", "n_eff_in"
  ) %in% names(screened$subjectwise)))
})

test_that("canonical smoothed preprocessing preserves pair-to-rest ALR and delta-time rate", {
  sm_mat <- rbind(
    target = c(0.20, 0.30, 0.25, 0.35),
    partner = c(0.30, 0.30, 0.35, 0.25),
    other = c(0.50, 0.40, 0.40, 0.40)
  )
  colnames(sm_mat) <- paste0("s", 1:4)
  metadata <- data.frame(Sample = colnames(sm_mat), subject = "A", time = c(0, 2, 5, 9))

  out <- pclvbayes:::.make_pair_inputs_glv(
    sm_mat = sm_mat, meta_df = metadata, j = "partner", i = "target",
    min_pairs = 1, min_sd = 0, zero_mode_alr = "fixed", eps_fixed = 1e-12,
    alr_cap = 100, smooth_scale = "logra", nz_partner_min_frac = 0
  )
  retained_alr_i <- log(c(0.20 / 0.50, 0.30 / 0.40, 0.25 / 0.40))
  retained_alr_j <- log(c(0.30 / 0.50, 0.30 / 0.40, 0.35 / 0.40))
  expected_xi <- retained_alr_i[1:2]
  expected_xj <- retained_alr_j[1:2]
  expect_equal(out$y, diff(retained_alr_i) / c(2, 3), tolerance = 1e-9)
  expect_equal(out$xi_unscaled, expected_xi, tolerance = 1e-9)
  expect_equal(out$xj_unscaled, expected_xj, tolerance = 1e-9)
})

test_that("lagged predictors and delta ALR rates stay within subjects", {
  pair_df <- data.frame(
    subject = rep(c("A", "B"), each = 3),
    time = c(0, 2, 5, 1, 4, 10),
    xi_raw = c(0.20, 0.30, 0.40, 0.10, 0.20, 0.25),
    xj_raw = c(0.30, 0.20, 0.10, 0.20, 0.30, 0.35)
  )
  alr_i <- log(
    pair_df$xi_raw /
      (1 - pair_df$xi_raw - pair_df$xj_raw)
  )
  alr_j <- log(
    pair_df$xj_raw /
      (1 - pair_df$xi_raw - pair_df$xj_raw)
  )

  out <- pclvbayes:::.make_pair_inputs_glv(
    pair_df,
    zero_mode_alr = "fixed",
    eps_fixed = 1e-12,
    alr_cap = 100,
    smooth_scale = "logra",
    nz_partner_min_frac = 0
  )

  expect_equal(out$subject, c("A", "A", "B", "B"))
  expect_equal(out$time, c(2, 5, 4, 10))
  expected_xi <- alr_i[c(1, 2, 4, 5)]
  expect_equal(
    out$xi,
    (expected_xi - mean(expected_xi)) / stats::sd(expected_xi),
    tolerance = 1e-9
  )
  expected_xj <- alr_j[c(1, 2, 4, 5)]
  expect_equal(
    out$xj,
    (expected_xj - mean(expected_xj)) / stats::sd(expected_xj),
    tolerance = 1e-9
  )
  expect_equal(out$xi_unscaled, expected_xi, tolerance = 1e-9)
  expect_equal(out$xj_unscaled, expected_xj, tolerance = 1e-9)
  expect_equal(
    out$y,
    c(
      diff(alr_i[1:3]) / c(2, 3),
      diff(alr_i[4:6]) / c(3, 6)
    ),
    tolerance = 1e-9
  )
})

test_that("irregular-time OU indexes and dt reset by subject", {
  pair_in <- data.frame(
    subject = c("B", "A", "D", "C", "A", "B", "C", "D"),
    time = c(5, 2, 9, 4, 7, 11, 10, 12),
    y = seq_len(8),
    xi = seq_len(8) + 10,
    xj = seq_len(8) + 20
  )

  split <- pclvbayes:::.build_train_test(
    pair_in,
    train_subjects = c("A", "B"),
    test_subjects = c("C", "D"),
    min_pairs = 1
  )

  expect_equal(split$train$subject, c("A", "A", "B", "B"))
  expect_equal(split$prev_train, c(0L, 1L, 0L, 3L))
  expect_equal(split$dt_train, c(0, 5, 0, 6))
  expect_equal(split$train$y, c(2, 5, 1, 6))

  expect_equal(split$test$subject, c("C", "C", "D", "D"))
  expect_equal(split$prev_test, c(0L, 1L, 0L, 3L))
  expect_equal(split$dt_test, c(0, 6, 0, 3))
  expect_equal(split$test$y, c(4, 7, 3, 8))
})

test_that("Core scaling is global for full data and training-only for held-out subjects", {
  pair_in <- data.frame(
    subject = rep(c("train_a", "train_b", "test"), each = 2),
    time = rep(c(0, 3), 3),
    y = c(1, 2, 3, 4, 100, 200),
    xi = c(1, 3, 5, 7, 101, 103),
    xj = c(2, 4, 8, 10, 202, 204)
  )

  train_rows <- pair_in$subject != "test"
  canonical <- pair_in
  canonical$xi_unscaled <- canonical$xi
  canonical$xj_unscaled <- canonical$xj
  canonical$xi <- (canonical$xi - mean(canonical$xi)) / stats::sd(canonical$xi)
  canonical$xj <- (canonical$xj - mean(canonical$xj)) / stats::sd(canonical$xj)

  expect_equal(mean(canonical$xi), 0, tolerance = 1e-12)
  expect_equal(stats::sd(canonical$xi), 1, tolerance = 1e-12)
  expect_equal(mean(canonical$xj), 0, tolerance = 1e-12)
  expect_equal(stats::sd(canonical$xj), 1, tolerance = 1e-12)
  expect_equal(canonical$y, pair_in$y)

  xi_mean <- mean(pair_in$xi[train_rows])
  xi_sd <- stats::sd(pair_in$xi[train_rows])
  xj_mean <- mean(pair_in$xj[train_rows])
  xj_sd <- stats::sd(pair_in$xj[train_rows])

  split <- pclvbayes:::.build_train_test(
    canonical,
    train_subjects = c("train_a", "train_b"),
    test_subjects = "test",
    min_pairs = 1
  )

  expect_equal(split$train$xi, (c(1, 3, 5, 7) - xi_mean) / xi_sd)
  expect_equal(split$train$xj, (c(2, 4, 8, 10) - xj_mean) / xj_sd)
  expect_equal(split$test$xi, (c(101, 103) - xi_mean) / xi_sd)
  expect_equal(split$test$xj, (c(202, 204) - xj_mean) / xj_sd)
  expect_equal(split$train$y, c(1, 2, 3, 4))
  expect_equal(split$test$y, c(100, 200))

  failure <- pclvbayes:::.build_train_test(
    data.frame(subject = c("train", "train", "test"), time = c(0, 2, 5),
               y = c(2, 4, 8), xi = c(3, 3, 9), xj = c(1, 2, 3)),
    "train", "test", min_pairs = 1
  )
  expect_s3_class(failure, "pclv_failure")
  expect_identical(failure$stage, "kfold_training")
})


test_that("Core predictor variation validation is deterministic", {
  xi_failure <- pclvbayes:::.validate_predictor_variation(
    c(2, 2, 2), c(1, 2, 3), "full_data"
  )
  expect_s3_class(xi_failure, "pclv_failure")
  expect_identical(xi_failure$reason, "insufficient_predictor_variation")
  expect_identical(xi_failure$predictor, "xi")
  expect_equal(xi_failure$observed_sd, 0)
  expect_equal(xi_failure$required_sd, sqrt(.Machine$double.eps) * 2)

  xj_failure <- pclvbayes:::.validate_predictor_variation(
    c(1, 2, 3), c(4, 4, 4), "full_data"
  )
  expect_s3_class(xj_failure, "pclv_failure")
  expect_identical(xj_failure$predictor, "xj")

  near_failure <- pclvbayes:::.validate_predictor_variation(
    c(1, 1 + 1e-10, 1 + 2e-10), c(1, 2, 3), "full_data"
  )
  expect_s3_class(near_failure, "pclv_failure")
  expect_lte(near_failure$observed_sd, near_failure$required_sd)

  expect_null(pclvbayes:::.validate_predictor_variation(
    c(1, 2, 3), c(2, 4, 8), "full_data"
  ))
})

test_that("full-data validation occurs on post-lag ALR predictors", {
  make_pair <- function(constant) {
    varying <- seq(0.10, 0.20, length.out = 6)
    if (constant == "xi") {
      xj <- varying
      xi <- (1 - xj) / 3
    } else {
      xi <- varying
      xj <- (1 - xi) / 3
    }
    data.frame(subject = "A", time = 0:5, xi_raw = xi, xj_raw = xj)
  }
  build <- function(x) pclvbayes:::.make_pair_inputs_glv(
    x, zero_mode_alr = "fixed",
    eps_fixed = 1e-8, alr_cap = 12, smooth_scale = "logra",
    nz_partner_min_frac = 0
  )

  xi_failure <- build(make_pair("xi"))
  expect_s3_class(xi_failure, "pclv_failure")
  expect_identical(xi_failure$stage, "full_data")
  expect_identical(xi_failure$predictor, "xi")

  xj_failure <- build(make_pair("xj"))
  expect_s3_class(xj_failure, "pclv_failure")
  expect_identical(xj_failure$stage, "full_data")
  expect_identical(xj_failure$predictor, "xj")
})

test_that("fold validation uses training predictors without changing y", {
  pair_in <- data.frame(
    subject = c(rep("train", 3), rep("test", 2)),
    time = c(0, 2, 5, 0, 7), y = c(1, 2, 3, 20, 30),
    xi = c(1, 1 + 1e-10, 1 + 2e-10, 1e12, -1e12),
    xj = c(1, 2, 4, 1e12, -1e12)
  )
  failure <- pclvbayes:::.build_train_test(
    pair_in, "train", "test", min_pairs = 1
  )
  expect_s3_class(failure, "pclv_failure")
  expect_identical(failure$stage, "kfold_training")
  expect_identical(failure$predictor, "xi")

  valid <- pair_in
  valid$xi[valid$subject == "train"] <- c(1, 2, 4)
  split <- pclvbayes:::.build_train_test(valid, "train", "test", min_pairs = 1)
  expect_false(inherits(split, "pclv_failure"))
  expect_equal(split$train$y, c(1, 2, 3))
  expect_equal(split$test$y, c(20, 30))

  changed <- valid
  changed$xi[changed$subject == "test"] <- c(8e15, -8e15)
  changed$xj[changed$subject == "test"] <- c(-9e15, 9e15)
  changed <- pclvbayes:::.build_train_test(changed, "train", "test", min_pairs = 1)
  expect_equal(changed$train$xi, split$train$xi)
  expect_equal(changed$train$xj, split$train$xj)
})
test_that("directed fitting uses the fixed canonical Core choices", {
  expect_false("standardize_by_subject" %in% names(formals(pclvbayes::fit_pclv_bayes)))
  expect_false("alr_cap" %in% names(formals(pclvbayes::fit_pclv_bayes)))
  expect_equal(pclvbayes:::.PCLV_CORE_ALR_CAP, 12)
  deleted <- c(
    "transform", "lag", "resid_mode", "use_student_t", "nu_fix",
    "compute_elpd", "elpd_mode", "spline_df", "spline_spar", "spline_cv"
  )
  expect_length(intersect(deleted, names(formals(pclvbayes::fit_pclv_bayes))), 0)
  expect_identical(pclvbayes:::.PCLV_CORE_TRANSFORM, "alr")
  expect_identical(pclvbayes:::.PCLV_CORE_LAG, 1L)
  expect_identical(pclvbayes:::.PCLV_CORE_RESID_MODE, "ou")
  expect_true(pclvbayes:::.PCLV_CORE_USE_STUDENT_T)
  expect_false(exists(".PCLV_CORE_NU", envir = asNamespace("pclvbayes"), inherits = FALSE))
  expect_true(pclvbayes:::.PCLV_CORE_COMPUTE_ELPD)
  expect_identical(pclvbayes:::.PCLV_CORE_ELPD_MODE, "kalman")
  expect_null(pclvbayes:::.PCLV_CORE_SPLINE$df)
  expect_null(pclvbayes:::.PCLV_CORE_SPLINE$spar)
  expect_true(pclvbayes:::.PCLV_CORE_SPLINE$cv)

  expect_false(any(c("transform", "lag") %in%
                     names(formals(pclvbayes:::.make_pair_inputs_glv))))
  expect_false(any(c("resid_mode", "use_t", "elpd_mode") %in%
                     names(formals(pclvbayes:::.fold_fit_and_score))))
  expect_false(any(c("resid_mode", "use_t", "elpd_mode") %in%
                     names(formals(pclvbayes:::.repkfold_eval))))

  observed <- new.env(parent = emptyenv())

  local_mocked_bindings(
    get_pclv_model = function(...) {
      structure(list(), class = "mock_pclv_model")
    },
    .make_pair_inputs_glv = function(...) {
      args <- list(...)
      observed$alr_cap <- args$alr_cap
      pclvbayes:::.pclv_failure("test_seam", "captured", list())
    },
    .package = "pclvbayes"
  )

  pair_builder <- function(target, partner, ctx, eps, min_pairs) {
    data.frame(
      subject = rep("A", 4),
      time = 0:3,
      y = rep(0, 4),
      xi_raw = c(0.1, 0.2, 0.3, 0.4),
      xj_raw = c(0.2, 0.3, 0.2, 0.1)
    )
  }

  ctx <- list(
    mod_exe_file = NULL,
    meta_df = data.frame(),
    sm_mat = matrix(numeric(), 0, 0),
    eps = 1e-6,
    min_pairs = 1,
    zero_mode_alr = "fixed",
    minpos_alpha = 0.5,
    minpos_base = "ij",
    eps_fixed = 1e-6,
    lib_eps_c = 0.65,
    rest_floor_frac = 1,
    smooth_scale = "logra",
    alr_spline_df = NULL,
    alr_spline_spar = NULL,
    alr_spline_cv = TRUE,
    nz_partner_min_frac = 0,
    max_retries = 0,
    chains = 1,
    iter_warmup = 1,
    iter_sampling = 1,
    adapt_delta = 0.8,
    max_treedepth = 10,
    metric = "diag_e",
    init = 0.2,
    seed = 1,
    quiet = TRUE,
    silent_sampler = TRUE,
    n_workers_kfold_eff = 1,
    kfold_K = 2,
    kfold_R = 1,
    use_pathfinder_init = FALSE,
    pf_num_paths = 1,
    pf_draws = 1,
    pf_history_size = 1,
    pf_max_lbfgs_iters = 1,
    pf_psis_resample = FALSE,
    pair_builder = pair_builder
  )

  result <- pclvbayes:::.run_one(
    target = "target",
    partner = "partner",
    ctx = ctx,
    seed_override = 1,
    progress_local = "none"
  )

  expect_s3_class(result, "pclv_failure")
  expect_identical(result$reason, "captured")
  expect_equal(observed$alr_cap, 12)
})

test_that("a failed direction does not discard its opposite direction", {
  calls <- list()
  run_one_double <- function(
    target,
    partner,
    ctx,
    seed_override,
    progress_local
  ) {
    calls[[length(calls) + 1L]] <<- list(
      target = target,
      partner = partner,
      ctx = ctx,
      seed_override = seed_override,
      progress_local = progress_local
    )
    if (length(calls) == 1L) return(NULL)
    success <- pclvbayes:::.failed_direction_result(
      pclvbayes:::.pclv_failure("test", "unused", list())
    )
    success$ok <- TRUE
    success$failure <- NULL
    success
  }

  ctx <- list(marker = "ctx")
  result <- pclvbayes:::.run_pair(
    1,
    2,
    taxa_vec = c("a", "b"),
    .run_one = run_one_double,
    ctx = ctx,
    progress = "verbose",
    mute_logs = TRUE,
    seed_base = 7
  )

  expect_s3_class(result, "tbl_df")
  expect_false(result$direction_ok_ij)
  expect_true(result$direction_ok_ji)
  expect_s3_class(result$failure_ij[[1]], "pclv_failure")
  expect_identical(result$failure_ij[[1]]$reason, "missing_result")
  expect_length(calls, 2)
  expect_identical(calls[[1]]$target, "a")
  expect_identical(calls[[1]]$partner, "b")
  expect_identical(calls[[2]]$target, "b")
  expect_identical(calls[[2]]$partner, "a")
  expect_identical(calls[[1]]$ctx, ctx)
  expect_identical(calls[[2]]$ctx, ctx)
  expect_identical(calls[[1]]$progress_local, "none")
  expect_identical(calls[[2]]$progress_local, "none")
  expect_identical(calls[[1]]$seed_override, 102008L)
  expect_identical(calls[[2]]$seed_override, 102009L)
})

test_that("NULL pointwise payloads compact to an empty result", {
  payload <- tibble::tibble(
    i = "a",
    j = "b",
    kfold_subject_ids_ij = list(NULL),
    kfold_subject_ij = list(NULL),
    kfold_subject_ppd_ij = list(NULL),
    kfold_elpd_method_ij = NA_character_,
    kfold_subject_counts_ij = list(NULL),
    kfold_subject_ids_ji = list(NULL),
    kfold_subject_ji = list(NULL),
    kfold_subject_ppd_ji = list(NULL),
    kfold_elpd_method_ji = NA_character_,
    kfold_subject_counts_ji = list(NULL)
  )

  compacted <- pclvbayes:::.expand_cross_pw(payload)

  expect_s3_class(compacted, "tbl_df")
  expect_equal(nrow(compacted), 0)
  expect_equal(ncol(compacted), 0)
})

test_that("LFSR conversion and model weights retain normalization", {
  expect_equal(
    pclvbayes:::.lfsr_safe(c(-1, 0, 0.2, 1, 2, 3, NA_real_)),
    c(0, 0, 0.1, 0.5, 0.5, 0.5, NA_real_)
  )

  pointwise <- tibble::tibble(
    from = rep(c("a", "b", "c", "d"), each = 2),
    to = rep(
      c("target_1", "target_1", "target_2", "target_2"),
      each = 2
    ),
    subject = rep(c("s1", "s2"), 4),
    elpd = c(-1, -2, -3, -4, -2, -2, -1, -1),
    n_test = rep(1L, 8)
  )

  weights <- pclvbayes:::.compute_pseudobma_weights_cross(
    list(elpd_pointwise_cross = pointwise),
    plus = FALSE,
    min_models = 2,
    min_subjects = 2
  )

  sums <- tapply(weights$weight, weights$to, sum)
  expect_equal(as.numeric(sums), rep(1, length(sums)), tolerance = 1e-12)
  expect_true(all(is.finite(weights$weight)))
  expect_true(all(weights$weight >= 0 & weights$weight <= 1))
})

test_that("K-fold aggregation preserves unavailable and partial evidence", {
  calls <- new.env(parent = emptyenv())
  calls$n <- c(A = 0L, B = 0L)
  fold_double <- function(...) {
    args <- list(...)
    held_out <- args[[6]]
    calls$n[[held_out]] <- calls$n[[held_out]] + 1L
    if (held_out == "A" || calls$n[[held_out]] == 2L) {
      return(structure(
        list(
          stage = "kfold_training",
          reason = "insufficient_predictor_variation",
          predictor = "xj",
          observed_sd = 0,
          required_sd = sqrt(.Machine$double.eps)
        ),
        class = c("pclv_failure", "list")
      ))
    }
    list(
      elpd = stats::setNames(-5, held_out),
      elpd_ppd = stats::setNames(-2.5, held_out),
      n_obs = stats::setNames(2L, held_out),
      fold_diag = data.frame(
        n_retries = 0L, nu_mean = 5, ebfmi_min = 0.9,
        worst_rhat = 1, min_ess_bulk = 500,
        treedepth_hits = 0L, n_divergent = 0L
      )
    )
  }

  local_mocked_bindings(
    .fold_fit_and_score = fold_double,
    .package = "pclvbayes"
  )
  result <- pclvbayes:::.repkfold_eval(
    mod = NULL, stan_list_base = list(), sample_args_base = list(seed = 1),
    pair_in = data.frame(
      subject = rep(c("A", "B"), each = 2), time = rep(0:1, 2),
      y = 1:4, xi = c(1, 2, 3, 4), xj = c(2, 4, 6, 8)
    ),
    K = 2, R = 2, n_workers_kfold = 1,
    min_pairs = 1
  )

  expect_true(is.na(result$elpd_subject[["A"]]))
  expect_equal(result$elpd_subject[["B"]], -5)
  expect_equal(result$subject_success_counts, c(A = 0L, B = 1L))
  expect_equal(result$subject_failure_counts, c(A = 2L, B = 1L))
  expect_equal(result$subject_test_counts, c(A = 0L, B = 2L))
  expect_equal(result$total_successful_evaluations, 1L)
  expect_equal(result$nu_fold_means, 5)
  expect_equal(sum(is.finite(result$elpd_subject)), 1)
  expect_length(result$failures, 3)
  expect_true(all(c(
    "ok", "stage", "reason", "details", "predictor",
    "observed_sd", "required_sd", "repetition", "fold", "test_subjects"
  ) %in% names(result$failures[[1]])))
})

test_that("pointwise ELPD and weights require successful common evidence", {
  pair <- tibble::tibble(
    i = "target", j = "predator",
    kfold_subject_ids_ij = list(c("A", "B")),
    kfold_subject_ij = list(c(A = NA_real_, B = -5)),
    kfold_subject_ppd_ij = list(c(A = NA_real_, B = -2.5)),
    kfold_elpd_method_ij = "kalman-ou",
    kfold_subject_counts_ij = list(c(A = 0L, B = 2L)),
    kfold_subject_ids_ji = list(NULL),
    kfold_subject_ji = list(NULL),
    kfold_subject_ppd_ji = list(NULL),
    kfold_elpd_method_ji = NA_character_,
    kfold_subject_counts_ji = list(NULL)
  )
  pointwise <- pclvbayes:::.expand_cross_pw(pair)
  expect_true(is.na(pointwise$elpd[pointwise$subject == "A"]))
  expect_equal(pointwise$n_test[pointwise$subject == "A"], 0L)
  expect_equal(pointwise$n_test[pointwise$subject == "B"], 2L)

  evidence <- tibble::tibble(
    from = rep(c("p1", "p2"), each = 3), to = "target",
    subject = rep(c("A", "B", "C"), 2),
    elpd = c(NA, -2, -2.5, -1, -3, -3.5),
    elpd_ppd = c(NA, -2, -2.5, -1, -3, -3.5),
    n_test = c(0L, 2L, 2L, 1L, 2L, 2L)
  )
  unavailable <- pclvbayes:::.compute_pseudobma_weights_cross(
    list(elpd_pointwise_cross = evidence), min_models = 2, min_subjects = 3
  )
  expect_equal(nrow(unavailable), 0)

  weights <- pclvbayes:::.compute_pseudobma_weights_cross(
    list(elpd_pointwise_cross = evidence), min_models = 2, min_subjects = 2
  )
  expect_equal(sum(weights$weight), 1)
  expect_equal(sort(weights$from), c("p1", "p2"))

  stacking_unavailable <- pclvbayes:::.compute_stacking_weights_cross(
    list(elpd_pointwise_cross = evidence), min_models = 2, min_subjects = 3
  )
  expect_equal(nrow(stacking_unavailable), 0)
  stacking <- pclvbayes:::.compute_stacking_weights_cross(
    list(elpd_pointwise_cross = evidence), min_models = 2, min_subjects = 2
  )
  if (requireNamespace("loo", quietly = TRUE)) {
    expect_equal(sum(stacking$stacking), 1, tolerance = 1e-8)
    expect_equal(sort(stacking$from), c("p1", "p2"))
  } else {
    expect_equal(nrow(stacking), 0)
  }
})

test_that("canonical projection scoring is Kalman OU with draw-specific Student-t nu", {
  draws <- data.frame(
    r0 = c(0, 0.1), a_ii = c(-0.2, -0.1), a_ij = c(0.3, 0.2),
    sigma = c(0.4, 0.5), sd_ou = c(0.2, 0.25), lambda = c(0.7, 0.8),
    tau_r = c(0, 0), nu = c(4, 8)
  )
  held_out <- data.frame(
    subject = c("A", "A", "B", "B"), time = c(0, 2, 1, 5),
    y = c(0.1, 0.2, -0.1, 0.3), xi = c(-1, 0, 0.5, 1),
    xj = c(0.2, 0.4, -0.2, 0.1)
  )
  scored <- pclvbayes:::.proj_loglik_subject(draws, held_out)
  expect_equal(dim(scored$full), c(2, 2))
  expect_true(all(is.finite(scored$full)))
  expect_equal(scored$subjects, c("A", "B"))
  expect_equal(scored$n_obs, c(2L, 2L))
  expect_identical(names(formals(pclvbayes:::.proj_loglik_subject)),
                   c("draws_df", "pair_in"))
})

test_that("v0.2 public fit API is frozen and prototype is absent", {
  retained <- c(
    "physeq", "subject_col", "time_col", "taxa_vec", "nz_partner_min_frac",
    "min_unique_times", "min_pairs", "chains", "iter_warmup", "iter_sampling",
    "seed", "init", "adapt_delta", "max_treedepth", "progress",
    "n_workers_outer", "n_workers_kfold", "kfold_K", "kfold_R", "kfold_seed"
  )
  expect_identical(names(formals(pclvbayes::fit_pclv_bayes)), retained)
  removed <- c("zero_mode_alr", "minpos_alpha", "minpos_base", "smooth_scale",
               "alr_spline_df", "alr_spline_spar", "alr_spline_cv", "eps",
               "eps_fixed", "lib_eps_c", "rest_floor_frac", "metric", "quiet",
               "progress_every", "silent_sampler", "max_retries", "use_pathfinder_init",
               "pf_num_paths", "pf_draws", "pf_history_size", "pf_max_lbfgs_iters",
               "pf_psis_resample")
  expect_length(intersect(removed, names(formals(pclvbayes::fit_pclv_bayes))), 0L)
  expect_false("fit_pclv_bayes2" %in% getNamespaceExports("pclvbayes"))
  expect_false(file.exists(testthat::test_path("../../man/fit_pclv_bayes2.Rd")))
  expect_error(do.call(pclvbayes::fit_pclv_bayes,
                       c(list(physeq = list(), subject_col = "subject", time_col = "time"),
                         list(eps = 1e-6))), "unused argument")
})
