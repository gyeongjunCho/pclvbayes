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

test_that("smoothed pair builder uses pair-to-rest ALR and delta-time rate", {
  sm_mat <- rbind(
    target = c(0.20, 0.30, 0.25),
    partner = c(0.30, 0.30, 0.35),
    other = c(0.50, 0.40, 0.40)
  )
  colnames(sm_mat) <- c("s1", "s2", "s3")

  metadata <- data.frame(
    Sample = c("s1", "s2", "s3"),
    subject = "A",
    time = c(0, 2, 5)
  )

  out <- pclvbayes:::.build_pair_df_smoothed(
    sm_mat,
    metadata,
    j = "partner",
    i = "target",
    eps = 1e-12,
    min_pairs = 1,
    min_sd = 0
  )

  expected_alr <- log(c(
    0.20 / 0.50,
    0.30 / 0.40,
    0.25 / 0.40
  ))

  expect_equal(
    out$y,
    diff(expected_alr) / c(2, 3),
    tolerance = 1e-9
  )
  expect_equal(out$xi_raw, c(0.20, 0.30))
  expect_equal(out$xj_raw, c(0.30, 0.30))
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
    transform = "alr",
    lag = 1,
    zero_mode_alr = "fixed",
    eps_fixed = 1e-12,
    alr_cap = 100,
    smooth_scale = "logra",
    nz_partner_min_frac = 0,
    standardize_by_subject = FALSE,
    z_mode = "none"
  )

  expect_equal(out$subject, c("A", "A", "B", "B"))
  expect_equal(out$time, c(2, 5, 4, 10))
  expect_equal(
    out$xi,
    alr_i[c(1, 2, 4, 5)],
    tolerance = 1e-9
  )
  expect_equal(
    out$xj,
    alr_j[c(1, 2, 4, 5)],
    tolerance = 1e-9
  )
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
    use_global_scaling = FALSE,
    min_pairs = 1
  )

  expect_equal(split$train$subject, c("A", "A", "B", "B"))
  expect_equal(split$prev_train, c(0L, 1L, 0L, 3L))
  expect_equal(split$dt_train, c(0, 5, 0, 6))

  expect_equal(split$test$subject, c("C", "C", "D", "D"))
  expect_equal(split$prev_test, c(0L, 1L, 0L, 3L))
  expect_equal(split$dt_test, c(0, 6, 0, 3))
})

test_that("directed fitting currently applies a fixed effective ALR cap", {
  observed <- new.env(parent = emptyenv())

  local_mocked_bindings(
    get_pclv_model = function(...) {
      structure(list(), class = "mock_pclv_model")
    },
    .make_pair_inputs_glv = function(...) {
      args <- list(...)
      observed$alr_cap <- args$alr_cap
      NULL
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
    compute_elpd = FALSE,
    standardize_by_subject = FALSE,
    transform = "alr",
    lag = 1,
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
    nu_fix = 5,
    use_student_t = TRUE,
    elpd_mode = "kalman",
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
    resid_mode = "ou",
    use_pathfinder_init = FALSE,
    pf_num_paths = 1,
    pf_draws = 1,
    pf_history_size = 1,
    pf_max_lbfgs_iters = 1,
    pf_psis_resample = FALSE,
    pair_builder = pair_builder,
    alr_cap = 2
  )

  result <- pclvbayes:::.run_one(
    target = "target",
    partner = "partner",
    ctx = ctx,
    seed_override = 1,
    progress_local = "none"
  )

  expect_null(result)
  expect_equal(observed$alr_cap, 12)
  expect_false(identical(observed$alr_cap, ctx$alr_cap))
})

test_that("NULL directed results discard the whole pair", {
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
    if (length(calls) == 1L) NULL else list(n_pairs = 1)
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

  expect_null(result)
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
    elpd = c(-1, -2, -3, -4, -2, -2, -1, -1)
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
