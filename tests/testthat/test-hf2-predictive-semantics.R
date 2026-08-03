hf2_pair_payload <- function() {
  tibble::tibble(
    i = "target_i", j = "source_j",
    kfold_subject_ids_ij = list(c("s1", "s2")),
    kfold_subject_ij = list(c(s1 = -2, s2 = -4)),
    kfold_subject_ppd_ij = list(c(s1 = -1, s2 = -2)),
    kfold_elpd_method_ij = "student-t-scale-mixture-kalman-ou-q16",
    kfold_subject_counts_ij = list(c(s1 = 2L, s2 = 3L)),
    kfold_subject_ids_ji = list(c("s1", "s2")),
    kfold_subject_ji = list(c(s1 = -3, s2 = -5)),
    kfold_subject_ppd_ji = list(c(s1 = -1.5, s2 = -2.5)),
    kfold_elpd_method_ji = "student-t-scale-mixture-kalman-ou-q16",
    kfold_subject_counts_ji = list(c(s1 = 4L, s2 = 5L))
  )
}

hf2_run_kfold <- function(subjects, split_seed, sampling_seed) {
  observed_seeds <- integer()
  observed_train <- list()
  observed_test <- list()
  testthat::local_mocked_bindings(
    .fold_fit_and_score = function(mod, stan_list_base, sample_args_base,
                                   pair_in, train_subjects, test_subjects,
                                   max_retries = 3, silent_sampler = TRUE,
                                   sample_args_override = NULL,
                                   freeze_retry_hypers = FALSE,
                                   seed_override = NULL, min_pairs = 4L) {
      held_out <- as.character(test_subjects)
      observed_train[[length(observed_train) + 1L]] <<- as.character(train_subjects)
      observed_test[[length(observed_test) + 1L]] <<- held_out
      observed_seeds <<- c(observed_seeds, as.integer(seed_override))
      list(
        elpd = stats::setNames(rep(-1, length(held_out)), held_out),
        elpd_ppd = stats::setNames(rep(-0.5, length(held_out)), held_out),
        n_obs = stats::setNames(rep(1L, length(held_out)), held_out),
        fold_diag = data.frame(
          n_retries = 0L, nu_mean = 5, ebfmi_min = 1,
          worst_rhat = 1, min_ess_bulk = 1000,
          treedepth_hits = 0L, n_divergent = 0L)
      )
    },
    .package = "pclvbayes"
  )
  result <- pclvbayes:::.repkfold_eval(
    mod = NULL, stan_list_base = list(),
    sample_args_base = list(seed = as.integer(sampling_seed)),
    pair_in = data.frame(
      subject = subjects, time = seq_along(subjects), y = seq_along(subjects),
      xi = seq_along(subjects), xj = rev(seq_along(subjects))),
    K = 3L, R = 2L, seed = as.integer(split_seed),
    n_workers_kfold = 1L, min_pairs = 1L, has_progressr = FALSE)
  list(result = result, sampler_seeds = observed_seeds,
       train_subjects = observed_train, test_subjects = observed_test)
}

test_that("public kfold_seed is propagated into the runtime context", {
  controls <- list(
    eps = 1e-6, min_unique_times = 3L, kfold_seed = 9876L,
    min_pairs = 4L, zero_mode_alr = "minpos_time", minpos_alpha = 0.5,
    minpos_base = "ij", eps_fixed = 1e-6, lib_eps_c = 0.65,
    rest_floor_frac = 1, smooth_scale = "logra", alr_spline_df = NULL,
    alr_spline_spar = NULL, alr_spline_cv = TRUE, nz_partner_min_frac = .15,
    max_retries = 3L, chains = 4L, iter_warmup = 2000L,
    iter_sampling = 2000L, adapt_delta = .98, max_treedepth = 14L,
    metric = "diag_e", init = .2, seed = 123L, quiet = TRUE,
    silent_sampler = TRUE, kfold_K = 3L, kfold_R = 2L,
    use_pathfinder_init = TRUE, pf_num_paths = 8L, pf_draws = 1000L,
    pf_history_size = 50L, pf_max_lbfgs_iters = 200L, pf_psis_resample = TRUE)
  validated <- list(
    controls = controls, mat_rel = matrix(.5, 2, 2),
    meta_df = data.frame(subject = c("a", "b"), time = 1:2),
    taxa_vec = c("a", "b"))
  model <- list(exe_file = function() "/tmp/model")
  testthat::local_mocked_bindings(
    get_pclv_model = function(...) model,
    .precompute_spline_smoothed = function(...) matrix(.5, 2, 2),
    .package = "pclvbayes")
  runtime <- pclvbayes:::.prepare_fit_runtime(validated)
  expect_identical(runtime$ctx$kfold_seed, 9876L)
})

test_that("shared split seed is canonical and separate from sampler seeds", {
  subjects <- rep(sprintf("subject-%02d", 1:12), each = 2L)
  a <- hf2_run_kfold(subjects, split_seed = 700L, sampling_seed = 100L)
  b <- hf2_run_kfold(factor(rev(subjects), levels = rev(unique(subjects))),
                     split_seed = 700L, sampling_seed = 900L)
  c <- hf2_run_kfold(subjects, split_seed = 701L, sampling_seed = 100L)

  expect_identical(a$result$splits_df, b$result$splits_df)
  expect_false(identical(a$result$splits_df, c$result$splits_df))
  expect_identical(a$sampler_seeds, c(1101L, 1102L, 1103L, 2101L, 2102L, 2103L))
  expect_identical(b$sampler_seeds, c(1901L, 1902L, 1903L, 2901L, 2902L, 2903L))
  expect_identical(a$sampler_seeds, c$sampler_seeds)
  eligible <- sort(unique(as.character(subjects)))
  expect_true(all(vapply(seq_along(a$train_subjects), function(idx) {
    train <- a$train_subjects[[idx]]
    test <- a$test_subjects[[idx]]
    !length(intersect(train, test)) && identical(sort(c(train, test)), eligible)
  }, logical(1))))
})

test_that("predictive context records split and sampling seeds separately", {
  fixture <- list(
    diagnostic_failure = list(NULL), .predictive_context = list(
      mod = NULL, stan_list = list(), pair_in = data.frame(subject = "s"),
      sample_args = list(seed = 456L), split_seed = 123L, sampling_seed = 456L,
      silent_sampler = TRUE, n_workers_kfold = 1L, max_retries = 0L,
      min_pairs = 1L, K = 1L, R = 1L, pair_tag = "b->a", progress = "none"))
  testthat::local_mocked_bindings(
    .repkfold_eval = function(mod, stan_list_base, sample_args_base, pair_in,
                              K, R, seed, silent_sampler, max_retries,
                              n_workers_kfold, min_pairs,
                              freeze_retry_hypers = FALSE, progress = "none",
                              has_progressr = FALSE) {
      expect_identical(seed, 123L)
      expect_identical(sample_args_base$seed, 456L)
      list(K = 1L, R = 1L, elpd_mean = -1, elpd_method = "q16",
           n_folds_ok = 1L, n_folds_fail = 0L, elpd_subject = c(s = -1),
           elpd_subject_ppd = c(s = -1), subject_test_counts = c(s = 1L),
           subject_success_counts = c(s = 1L), subject_failure_counts = c(s = 0L),
           total_successful_evaluations = 1L, failures = list(),
           splits_df = data.frame(r = 1L, k = 1L), retry_total = 0L,
           retry_mean = 0, nu_fold_means = 5)
    }, .package = "pclvbayes")
  out <- pclvbayes:::.add_predictive_evaluation(fixture)
  expect_identical(out$kfold_seed_used, 123L)
})

test_that("public kfold_seed defaults to seed and is recorded", {
  captured <- NULL
  testthat::local_mocked_bindings(
    .validate_fit_pclv_inputs = function(physeq, subject_col, time_col,
                                         taxa_vec, controls) {
      captured <<- controls
      list(controls = controls, subject_col = subject_col, time_col = time_col,
           taxa_vec = c("a", "b"), meta_df = data.frame(), mat_rel = matrix(0, 0, 0))
    },
    .prepare_fit_runtime = function(validated) {
      pclvbayes:::.pclv_failure("test_boundary", "stop_before_runtime")
    },
    .package = "pclvbayes")
  result <- fit_pclv_bayes(NULL, "subject", "time", seed = 4321L,
                           progress = "none")
  expect_s3_class(result, "pclv_failure")
  expect_identical(captured$seed, 4321L)
  expect_identical(captured$kfold_seed, 4321L)

  fixture <- list(
    diagnostic_failure = list(NULL), .predictive_context = list(
      mod = NULL, stan_list = list(), pair_in = data.frame(subject = "s"),
      sample_args = list(seed = 4321L), split_seed = captured$kfold_seed,
      sampling_seed = captured$seed, silent_sampler = TRUE,
      n_workers_kfold = 1L, max_retries = 0L, min_pairs = 1L,
      K = 1L, R = 1L, pair_tag = "b->a", progress = "none"))
  testthat::local_mocked_bindings(
    .repkfold_eval = function(..., seed) {
      expect_identical(seed, 4321L)
      list(K = 1L, R = 1L, elpd_mean = -1, elpd_method = "q16",
           n_folds_ok = 1L, n_folds_fail = 0L, elpd_subject = c(s = -1),
           elpd_subject_ppd = c(s = -1), subject_test_counts = c(s = 1L),
           subject_success_counts = c(s = 1L), subject_failure_counts = c(s = 0L),
           total_successful_evaluations = 1L, failures = list(),
           splits_df = data.frame(r = 1L, k = 1L), retry_total = 0L,
           retry_mean = 0, nu_fold_means = 5)
    }, .package = "pclvbayes")
  evaluated <- pclvbayes:::.add_predictive_evaluation(fixture)
  expect_identical(evaluated$kfold_seed_used, 4321L)
})

test_that("split manifests represent their exact eligible subject universe", {
  universe_a <- sprintf("subject-%02d", 1:12)
  universe_b <- sprintf("subject-%02d", 1:9)
  a <- hf2_run_kfold(rep(universe_a, each = 2L), 811L, 100L)$result$splits_df
  a_reordered <- hf2_run_kfold(rep(rev(universe_a), each = 2L), 811L, 900L)$result$splits_df
  b <- hf2_run_kfold(rep(universe_b, each = 2L), 811L, 100L)$result$splits_df
  expect_identical(a, a_reordered)
  held_a <- unlist(a$test_subjects, use.names = FALSE)
  held_b <- unlist(b$test_subjects, use.names = FALSE)
  expect_setequal(held_a, universe_a)
  expect_setequal(held_b, universe_b)
  expect_false(setequal(held_a, held_b))
  expect_true(all(table(held_a) == 2L))
  expect_true(all(table(held_b) == 2L))
})

test_that("HF2 fitted-result schema exposes only directed subject ELPD", {
  payload <- hf2_pair_payload()
  testthat::local_mocked_bindings(
    .mk_cross = function(x) tibble::tibble(kind = "cross"),
    .mk_self = function(x) tibble::tibble(kind = "self"),
    .package = "pclvbayes")
  result <- pclvbayes:::.assemble_public_fit_result(payload)
  expect_named(result, c("cross", "self", "elpd_subject_cross", "raw"))
  expect_false(any(c("elpd_pointwise_cross", "elpd_pointwise_self") %in% names(result)))
  expect_named(result$elpd_subject_cross, c(
    "from", "to", "subject", "elpd", "elpd_per_observation",
    "n_successful_test_observations", "elpd_method"))
  expect_equal(unique(result$elpd_subject_cross[result$elpd_subject_cross$to == "target_i", c("from", "to")]),
               data.frame(from = "source_j", to = "target_i"), ignore_attr = TRUE)
  expect_equal(unique(result$elpd_subject_cross[result$elpd_subject_cross$to == "source_j", c("from", "to")]),
               data.frame(from = "target_i", to = "source_j"), ignore_attr = TRUE)
})

test_that("empty subject ELPD retains the exact typed public schema", {
  payload <- hf2_pair_payload()
  for (nm in grep("^kfold_subject", names(payload), value = TRUE)) payload[[nm]] <- list(NULL)
  payload$kfold_elpd_method_ij <- NA_character_
  payload$kfold_elpd_method_ji <- NA_character_
  empty <- pclvbayes:::.expand_cross_subject_elpd(payload)
  expect_s3_class(empty, "tbl_df")
  expect_equal(nrow(empty), 0L)
  expect_named(empty, c("from", "to", "subject", "elpd",
                        "elpd_per_observation", "n_successful_test_observations",
                        "elpd_method"))
  expect_identical(vapply(empty, typeof, character(1)), c(
    from = "character", to = "character", subject = "character",
    elpd = "double", elpd_per_observation = "double",
    n_successful_test_observations = "integer", elpd_method = "character"))

  testthat::local_mocked_bindings(
    .mk_cross = function(x) tibble::tibble(),
    .mk_self = function(x) tibble::tibble(), .package = "pclvbayes")
  result <- pclvbayes:::.assemble_public_fit_result(payload)
  expect_named(result, c("cross", "self", "elpd_subject_cross", "raw"))
  expect_identical(result$elpd_subject_cross, empty)
})

test_that("public raw table preserves existing input attributes", {
  payload <- hf2_pair_payload()
  attr(payload, "hf2_test_attribute") <- "preserve-me"
  testthat::local_mocked_bindings(
    .mk_cross = function(x) tibble::tibble(),
    .mk_self = function(x) tibble::tibble(), .package = "pclvbayes")
  result <- pclvbayes:::.assemble_public_fit_result(payload)
  expect_identical(attr(result$raw, "hf2_test_attribute"), "preserve-me")
  expect_identical(unclass(result$raw), unclass(payload))
})

test_that("invalid cross-pair weights and self predictive payload are removed", {
  forbidden_formals <- c("use_true_stacking", "stacking_use_ppd",
                         "stacking_min_models", "stacking_min_subjects")
  expect_false(any(forbidden_formals %in% names(formals(summarize_bayes_pclv))))
  ns <- asNamespace("pclvbayes")
  expect_false(exists(".compute_stacking_weights_cross", envir = ns, inherits = FALSE))
  expect_false(exists(".compute_pseudobma_weights_cross", envir = ns, inherits = FALSE))
  expect_false(exists(".softmax", envir = ns, inherits = FALSE))
  expect_false(exists(".expand_self_pw", envir = ns, inherits = FALSE))
  expect_false(exists(".expand_subject_cross", envir = ns, inherits = FALSE))
})

test_that("cross and self summaries omit invalid model-weight columns", {
  fit <- list(
    cross = data.frame(from = "b", to = "a", n_subjects = 5L,
      a_mean = .2, a_q2.5 = .1, a_q97.5 = .3, p_sign2 = .02,
      diagnostic_class = "converged"),
    self = data.frame(taxon = c("a", "b"),
      a_self_mean = c(-.2, -.3), a_self_q2.5 = c(-.3, -.4),
      a_self_q97.5 = c(-.1, -.2), p_sign2_self = c(.02, .02)),
    raw = data.frame(i = "a", j = "b", n_pairs_ij = 20L, n_pairs_ji = 20L,
      rhat_ij = 1, essb_ij = 1000, esst_ij = 1000, div_ij = 0L,
      tdhit_ij = 0L, ebfmi_min_ij = .9, diagnostic_class_ij = "converged",
      rhat_ji = 1, essb_ji = 1000, esst_ji = 1000, div_ji = 0L,
      tdhit_ji = 0L, ebfmi_min_ji = .9, diagnostic_class_ji = "converged"))
  cross <- summarize_bayes_pclv(fit, interaction = "cross")$cross
  self <- summarize_bayes_pclv(fit, interaction = "self")$self
  forbidden <- c("stacking", "pseudo_BMA", "pseudo_BMA_plus")
  expect_false(any(forbidden %in% names(cross)))
  expect_false(any(forbidden %in% names(self)))
})
