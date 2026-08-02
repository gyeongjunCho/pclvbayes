source(testthat::test_path("../../benchmarks/mtist/mtist_adapter.R"))
source(testthat::test_path("../../benchmarks/mtist/score_mtist.R"))
source(testthat::test_path("../../benchmarks/mtist/ten_species_helpers.R"))

test_that("ten-species selection is deterministic and metadata-driven", {
  catalog <- data.frame(
    did = c(90L, 37L, 36L, 55L), n_species = 10L,
    noise = c(.01, .01, .01, .1), n_timeseries = c(10L, 10L, 10L, 10L),
    n_timepoints = c(15L, 15L, 5L, 15L),
    sampling_scheme = c("even", "even", "even", "even"),
    compatible = TRUE)
  selected <- .select_mtist_10species_catalog(catalog)
  expect_identical(selected$did, 37L)
  expect_identical(selected$n_species, 10L)
})

test_that("ten taxa produce 45 canonical tasks and 90 directions", {
  taxa <- paste0("species_", 0:9)
  index <- ten_species_direction_index(taxa, 101L)
  expect_equal(length(unique(index$task_id)), 45L)
  expect_equal(nrow(index), 90L)
  expect_identical(index$direction_index, seq_len(90L))
  expect_identical(index$target[1:2], c("species_0", "species_1"))
  expect_identical(index$source[1:2], c("species_1", "species_0"))
  expect_identical(index$seed[1:2], unname(pclvbayes:::.pair_direction_seeds(101L, 1L, 2L)))
})

test_that("truth orientation is row target and column source", {
  truth <- matrix(0, 3L, 3L,
    dimnames = list(paste0("species_", 0:2), paste0("species_", 0:2)))
  truth["species_0", "species_1"] <- -0.7
  truth["species_1", "species_0"] <- 0.2
  expect_equal(mtist_truth_coefficient(truth, "species_0", "species_1"), -0.7)
  expect_equal(mtist_truth_coefficient(truth, "species_1", "species_0"), 0.2)
})

make_benchmark_directions <- function(taxa = paste0("species_", 0:9)) {
  x <- ten_species_direction_index(taxa, 101L)
  x$diagnostic_class <- rep(c("converged", "interaction_stable_residual_unstable",
                              "interaction_indeterminate", "sampler_diagnostics_failed"),
                            length.out = nrow(x))
  x$main_fit_retained <- x$diagnostic_class != "sampler_diagnostics_failed"
  x$posterior_mean <- ifelse(x$main_fit_retained, rep(c(-.2, .3), length.out = nrow(x)), NA_real_)
  x$final_significant <- x$diagnostic_class == "converged"
  x$elpd_available <- x$diagnostic_class == "converged"
  x$a_self_mean <- ifelse(x$diagnostic_class == "converged", -.5, NA_real_)
  x$bayesian_eligible <- x$diagnostic_class == "converged"
  x$predicted_sign <- ifelse(x$bayesian_eligible, sign(x$posterior_mean), NA_real_)
  x$truth_sign <- rep(c(-1, 1, 0), length.out = nrow(x))
  x
}

test_that("all diagnostic classes and opposite directions remain explicit", {
  x <- make_benchmark_directions()
  expect_equal(nrow(x), 90L)
  expect_setequal(unique(x$diagnostic_class), c(
    "converged", "interaction_stable_residual_unstable",
    "interaction_indeterminate", "sampler_diagnostics_failed"))
  pair <- x[x$task_id == 1L, ]
  expect_equal(nrow(pair), 2L)
  expect_false(identical(pair$diagnostic_class[[1L]], pair$diagnostic_class[[2L]]))
})

test_that("external zeros retain separate masks and do not mutate indeterminate coefficients", {
  taxa <- paste0("species_", 0:9)
  x <- make_benchmark_directions(taxa)
  indeterminate <- which(x$diagnostic_class == "interaction_indeterminate")[[1L]]
  value <- x$posterior_mean[[indeterminate]]
  matrices <- build_mtist_benchmark_matrices(x, taxa)
  expect_true(is.finite(value) && value != 0)
  expect_equal(x$posterior_mean[[indeterminate]], value)
  expect_equal(matrices$posterior[x$target[[indeterminate]], x$source[[indeterminate]]], 0)
  expect_true(matrices$masks$interaction_indeterminate[[indeterminate]])
  expect_equal(nrow(matrices$masks), 90L)
})

test_that("coverage metrics expose conditional denominators", {
  x <- make_benchmark_directions()
  metrics <- coverage_sign_metrics(x)
  eligible <- metrics[metrics$metric == "conditional_sign_accuracy", ]
  expect_equal(eligible$denominator, sum(x$bayesian_eligible))
  expect_equal(metrics$denominator[metrics$metric == "eligible_sign_coverage"], 90L)
})

test_that("absolute-A aliases preserve legacy benchmark metrics", {
  x <- make_benchmark_directions()
  metrics <- coverage_sign_metrics(x)
  old <- metrics[metrics$metric == "conditional_sign_accuracy", ]
  explicit <- metrics[metrics$metric == "absolute_A_cross_estimand_accuracy", ]
  expect_equal(explicit$value, old$value)
  expect_equal(explicit$numerator, old$numerator)
  expect_equal(explicit$denominator, old$denominator)
  zero_old <- metrics[metrics$metric == "truth_zero_false_positive", ]
  zero_new <- metrics[metrics$metric == "absolute_A_zero_to_nonzero", ]
  expect_equal(zero_new$value, zero_old$value)
  expect_equal(zero_new$denominator, zero_old$denominator)
})

test_that("focused transformed-oracle fixture preserves six-direction denominator", {
  focused <- data.frame(absolute_A_sign = c(-1, 1, -1, 0, -1, -1),
                        posterior_sign = c(1, 1, -1, -1, -1, 1),
                        transformed_oracle_sign = c(1, 1, -1, -1, -1, 1),
                        state_dependent = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE))
  expect_equal(nrow(focused), 6L)
  expect_equal(sum(focused$posterior_sign == focused$transformed_oracle_sign), 6L)
  expect_equal(sum(focused$posterior_sign == focused$absolute_A_sign), 3L)
  expect_equal(sum(focused$state_dependent), 5L)
})

test_that("diagonal uses median across pair-specific posterior self-effect means", {
  x <- make_benchmark_directions()
  n_self <- sum(x$target == "species_0" & x$diagnostic_class == "converged")
  x$a_self_mean[x$target == "species_0" & x$diagnostic_class == "converged"] <-
    seq(-.9, -.1, length.out = n_self)
  diagonal <- aggregate_mtist_diagonal(x, paste0("species_", 0:9))
  row <- diagonal[diagonal$taxon == "species_0", ]
  values <- x$a_self_mean[x$target == "species_0" & x$diagnostic_class == "converged" & is.finite(x$a_self_mean)]
  expect_equal(row$aggregated_value, median(values))
  expect_equal(row$n_pair_specific_means, length(values))
  expect_identical(row$aggregation, "median across pair-specific posterior self-effect means")
})

test_that("official ES identity and diagonal exclusion remain available", {
  skip_if_not(dir.exists(path.expand("~/mtist")))
  root <- mtist_root()
  study <- load_mtist_study(37L, root)
  expect_equal(official_mtist_es(study$truth, study$truth, root, FALSE), 1)
  expect_equal(official_mtist_es(study$truth, study$truth, root, TRUE), 1)
  changed_diagonal <- study$truth
  diag(changed_diagonal) <- 0
  expect_equal(official_mtist_es(study$truth, changed_diagonal, root, TRUE), 1)
  expect_lt(official_mtist_es(study$truth, changed_diagonal, root, FALSE), 1)
})

test_that("benchmark adds no public fit arguments", {
  expect_false(any(c("dataset_id", "truth", "benchmark_stage") %in%
                     names(formals(fit_pclv_bayes))))
})

test_that("confirmation selection is deterministic and truth-independent", {
  x <- make_benchmark_directions()
  x$direction_index <- seq_len(nrow(x))
  x$seed <- seq_len(nrow(x)) + 100L
  selected <- select_confirmation_directions(x)
  changed <- x
  changed$truth_sign <- rev(changed$truth_sign)
  changed$predicted_sign <- -changed$predicted_sign
  changed$posterior_mean <- changed$posterior_mean * 100
  expect_identical(
    selected[c("task_id", "direction_index", "target", "source", "seed")],
    select_confirmation_directions(changed)[c("task_id", "direction_index", "target", "source", "seed")]
  )
})

test_that("confirmation configuration fixes main-posterior-only effort", {
  config <- source(testthat::test_path("../../benchmarks/mtist/configs/ten_species_confirmation.R"))$value
  expect_identical(config$chains, 4L)
  expect_identical(config$iter_warmup, 2000L)
  expect_identical(config$iter_sampling, 2000L)
  expect_false(config$use_pathfinder_init)
  expect_identical(config$max_retries, 0L)
  expect_false(config$run_kfold)
  expect_false(config$calculate_elpd)
  expect_false(config$calculate_stacking)
  expect_lte(config$n_workers_outer, 3L)
})

test_that("confirmation outcomes preserve the prespecified denominator", {
  classes <- c("converged", "converged", "interaction_stable_residual_unstable",
               "interaction_indeterminate", "sampler_diagnostics_failed", "converged")
  signs <- c(-1, 1, -1, 1, -1, 1)
  stage_signs <- c(-1, -1, -1, 1, -1, 1)
  outcomes <- mapply(confirmation_outcome, "converged", stage_signs,
                     classes, signs, USE.NAMES = FALSE)
  expect_length(outcomes, 6L)
  expect_setequal(unique(outcomes), c(
    "confirmed_converged_same_sign", "confirmed_converged_sign_changed",
    "downgraded_residual_unstable", "downgraded_indeterminate",
    "long_run_sampler_failed"))
})

test_that("downgraded values are retained rather than converted to zero", {
  value <- -.27
  outcome <- confirmation_outcome("converged", -1,
                                  "interaction_indeterminate", sign(value))
  expect_identical(outcome, "downgraded_indeterminate")
  expect_equal(value, -.27)
})
