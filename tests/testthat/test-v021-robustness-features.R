source(testthat::test_path("../../benchmarks/mtist/v021_truth_isolation.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_checkpoint_manifest.R"))
source(testthat::test_path("../../benchmarks/mtist/ten_species_helpers.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_robustness_features.R"))

.robust_manifest <- function() {
  taxa <- c("a", "b", "c")
  tasks <- ten_species_direction_index(taxa, 101L)
  build_v021_execution_manifest(
    "fixture", taxa, tasks[c("task_id", "direction_index", "target", "source", "seed")],
    101L, 202L, list(id = "prep"), list(id = "posterior"),
    list(K = 2L, R = 1L), list(id = "q16"), list(code_commit = "fixture"))
}

.support_data <- function() data.frame(
  subject = c("s2", "s1", "s1", "s1", "s2"),
  time = c(0, 0, 1, 3, 2),
  source_abundance = c(.2, 0, .1, .4, .3),
  target_abundance = c(.1, .2, 0, .2, .3),
  rest_abundance = c(.7, .8, .9, .4, .4),
  eligible = c(TRUE, TRUE, TRUE, TRUE, TRUE), stringsAsFactors = FALSE)

.finalized <- function() list(
  finalized = TRUE, terminal_state = "completed", posterior_mean = .5,
  posterior_median = .45, posterior_sd = .2, interval_lower = .1,
  interval_upper = .9, psp = .98, lfsr = .02, rhat_max = 1.01,
  bulk_ess_min = 1200, tail_ess_min = 900, divergence_count = 0L,
  treedepth_hit_count = 1L, ebfmi_min = .4, retry_count = 1L,
  diagnostic_class = "identified", predictive_eligible = TRUE,
  elpd = -20, successful_test_observation_count = 10L,
  folds_attempted = 5L, folds_failed = 1L)

.comparison <- function(mean = .4, available = TRUE) list(
  dataset_id = "fixture", task_identity = "direction-000001", source = "b", target = "a",
  taxa_order = c("a", "b", "c"), estimand = "directed_pair_to_rest_a_ij",
  available = available, posterior_mean = mean, posterior_sd = .2,
  interval_lower = .1, interval_upper = .7, psp = .9, lfsr = .1,
  diagnostic_class = "identified", predictive_eligible = TRUE, elpd = -21,
  successful_test_observation_count = 10L, transformed_outcome_hash = "outcome",
  heldout_split_hash = "split")

test_that("feature dictionary is complete typed unique and explicit", {
  dictionary <- v021_robustness_feature_dictionary()
  expect_identical(attr(dictionary, "dictionary_version"),
                   v021_robustness_feature_dictionary_version)
  expect_identical(anyDuplicated(dictionary$name), 0L)
  expect_true(all(nzchar(dictionary$formula)))
  expect_true(all(nzchar(dictionary$denominator)))
  expect_true(all(nzchar(dictionary$missingness_rule)))
  expect_true(all(dictionary$type %in% c("integer", "double", "logical", "character")))
  expect_true(all(c("observed_support", "posterior_stability",
                    "perturbation_stability") %in% dictionary$group))
  expect_false(any(.v021_truth_name(dictionary$name)))
  expect_false(any(grepl("noise|risk|unsafe|withhold|sign_reversal|absolute_zero",
                         dictionary$name, ignore.case = TRUE)))
})

test_that("observed support uses explicit observation subject and transition denominators", {
  result <- calculate_v021_observed_support_features(.support_data())
  x <- result$values
  expect_identical(result$schema_version, v021_observed_support_features_schema)
  expect_identical(x$eligible_subject_count, 2L)
  expect_identical(x$eligible_observation_count, 5L)
  expect_identical(x$canonical_transition_count, 3L)
  expect_identical(x$min_transitions_per_subject, 1L)
  expect_identical(x$median_transitions_per_subject, 1.5)
  expect_identical(x$max_transitions_per_subject, 2L)
  expect_equal(c(x$min_time_gap, x$median_time_gap, x$max_time_gap), c(1, 2, 2))
  expect_gt(x$time_gap_cv, 0)
  expect_identical(x$effective_observation_count, 3L)
})

test_that("feature validators enforce frozen types ranges and missingness", {
  observed <- calculate_v021_observed_support_features(.support_data())
  bad <- observed; bad$values$source_prevalence <- 2
  expect_error(validate_v021_observed_support_features(bad), "outside")
  bad <- observed; bad$values$eligible_subject_count <- 2
  expect_error(validate_v021_observed_support_features(bad), "typed")
  posterior <- extract_v021_posterior_stability_features(.finalized())
  bad <- posterior; bad$values$psp <- 1.2
  expect_error(validate_v021_posterior_stability_features(bad), "outside")
  bad <- posterior; bad$values$lfsr <- NA_real_; bad$missingness$lfsr <- NA_character_
  expect_error(validate_v021_posterior_stability_features(bad), "missingness")
})

test_that("prevalence zero patterns abundance and rest support are exact", {
  x <- calculate_v021_observed_support_features(.support_data())$values
  expect_equal(x$source_prevalence, .8)
  expect_equal(x$target_prevalence, .8)
  expect_equal(x$pair_coprevalence, .6)
  expect_equal(x$source_near_zero_fraction, .2)
  expect_equal(x$target_near_zero_fraction, .2)
  expect_equal(x$both_zero_fraction, 0)
  expect_equal(x$source_only_zero_fraction, .2)
  expect_equal(x$target_only_zero_fraction, .2)
  expect_equal(x$rest_support_min, .4)
  expect_equal(x$rest_support_median, .7)
  expect_equal(x$weak_rest_support_fraction, 0)
  expect_true(x$source_abundance_q10 <= x$source_abundance_median)
  expect_true(x$target_abundance_median <= x$target_abundance_q90)
})

test_that("invalid observations are not bridged and zero denominators remain missing", {
  data <- .support_data()
  data$eligible[data$subject == "s1" & data$time == 1] <- FALSE
  x <- calculate_v021_observed_support_features(data)
  expect_identical(x$values$canonical_transition_count, 1L)
  expect_identical(x$values$min_transitions_per_subject, 0L)
  single <- data[data$subject == "s1" & data$time == 0, ]
  y <- calculate_v021_observed_support_features(single)
  expect_identical(y$values$canonical_transition_count, 0L)
  expect_true(is.na(y$values$min_time_gap))
  expect_identical(y$missingness$min_time_gap, "zero_denominator")
  none <- single; none$eligible <- FALSE
  expect_error(calculate_v021_observed_support_features(none), "zero eligible")
})

test_that("posterior stability preserves maintained PSP LFSR diagnostics and ELPD denominator", {
  result <- extract_v021_posterior_stability_features(.finalized())
  x <- result$values
  expect_identical(result$schema_version, v021_posterior_stability_features_schema)
  expect_equal(x$psp, .98)
  expect_equal(x$lfsr, .02)
  expect_equal(x$posterior_interval_width, .8)
  expect_equal(x$posterior_distance_zero, .5)
  expect_equal(x$posterior_distance_zero_sd, 2.5)
  expect_equal(x$elpd_per_successful_observation, -2)
  expect_equal(x$fold_failure_fraction, .2)
  expect_identical(x$diagnostic_class, "identified")
  expect_true(x$predictive_eligible)
})

test_that("posterior extraction requires finalization and preserves zero denominators", {
  bad <- .finalized(); bad$finalized <- FALSE
  expect_error(extract_v021_posterior_stability_features(bad), "finalized")
  zero <- .finalized(); zero$posterior_sd <- 0
  zero$successful_test_observation_count <- 0L; zero$folds_attempted <- 0L
  x <- extract_v021_posterior_stability_features(zero)
  expect_true(is.na(x$values$posterior_distance_zero_sd))
  expect_true(is.na(x$values$elpd_per_successful_observation))
  expect_true(is.na(x$values$fold_failure_fraction))
  expect_identical(x$missingness$elpd_per_successful_observation, "zero_denominator")
})

test_that("approved perturbations are deterministic bounded and estimand preserving", {
  a <- v021_approved_perturbations(c("s3", "s1", "s2", "s4"), 3L)
  b <- v021_approved_perturbations(factor(c("s4", "s2", "s1", "s3")), 3L)
  expect_identical(a, b)
  expect_identical(grep("^delete_subject__", names(a), value = TRUE),
                   paste0("delete_subject__", c("s1", "s2", "s3")))
  expect_true(all(vapply(a, function(x) isTRUE(x$response) && isTRUE(x$predictor), logical(1))))
  expect_lte(sum(grepl("^delete_subject__", names(a))), 3L)
})

test_that("perturbation identity is stable truth-free and does not alter manifest", {
  manifest <- .robust_manifest(); before <- manifest$manifest_hash
  approved <- v021_approved_perturbations(c("s2", "s1"))
  a <- build_v021_perturbation_identity(
    "remove_first_observation_per_subject", approved, "fixture", "direction-000001",
    "b", "a", c("a", "b", "c"), c("s2", "s1"), manifest$configuration_hash)
  b <- build_v021_perturbation_identity(
    "remove_first_observation_per_subject", approved, "fixture", "direction-000001",
    "b", "a", c("a", "b", "c"), c("s1", "s2"), manifest$configuration_hash)
  expect_identical(a, b)
  expect_identical(manifest$manifest_hash, before)
  expect_error(build_v021_perturbation_identity(
    "invented", approved, "fixture", "direction-000001", "b", "a",
    c("a", "b", "c"), c("s1"), manifest$configuration_hash), "Unapproved")
})

test_that("perturbation comparisons calculate shifts signs and interval overlap", {
  canonical <- .comparison(.5); perturbed <- .comparison(-.25)
  result <- compare_v021_perturbation_stability(canonical, perturbed)
  x <- result$values
  expect_equal(x$coefficient_shift_abs, .75)
  expect_equal(x$coefficient_shift_relative, 1.5)
  expect_equal(x$coefficient_shift_sd, 3.75)
  expect_false(x$posterior_sign_agreement)
  expect_true(x$interval_overlap_fraction >= 0 && x$interval_overlap_fraction <= 1)
  expect_equal(x$elpd_change_same_outcome, 0)
})

test_that("ELPD comparison requires the same outcome and held-out split", {
  canonical <- .comparison(.5); perturbed <- .comparison(.4)
  perturbed$heldout_split_hash <- "different"
  result <- compare_v021_perturbation_stability(canonical, perturbed)
  expect_true(is.na(result$values$elpd_change_same_outcome))
  expect_identical(result$missingness_reason, "identity_mismatch")
  perturbed <- .comparison(.4); perturbed$source <- "c"
  expect_error(compare_v021_perturbation_stability(canonical, perturbed), "identity mismatch")
})

test_that("unavailable perturbations remain optional and explicit", {
  unavailable <- .comparison(.4, available = FALSE)
  result <- compare_v021_perturbation_stability(.comparison(.5), unavailable)
  expect_identical(result$validity_status, "unavailable")
  expect_identical(result$missingness_reason, "perturbation_failed")
  expect_null(result$values)
})

test_that("feature artifacts are atomic provenance-bound and separate from task status", {
  manifest <- .robust_manifest()
  status <- new_v021_task_status(manifest)
  support <- calculate_v021_observed_support_features(.support_data())
  posterior <- extract_v021_posterior_stability_features(.finalized())
  artifact <- build_v021_robustness_feature_artifact(
    manifest, 1L, support, posterior, code_provenance = "fixture")
  expect_silent(validate_v021_robustness_feature_artifact(artifact, manifest))
  path <- tempfile("v021-feature-", fileext = ".rds")
  expect_silent(write_v021_robustness_feature_artifact_atomic(artifact, path, manifest))
  expect_identical(readRDS(path), artifact)
  expect_identical(status, new_v021_task_status(manifest))
  corrupt <- artifact; corrupt$source <- "c"
  expect_error(validate_v021_robustness_feature_artifact(corrupt, manifest), "orientation")
})

test_that("truth perturbation cannot alter features and recursive leakage fails closed", {
  support_a <- calculate_v021_observed_support_features(.support_data())
  absolute_matrix_a <- diag(3); absolute_matrix_b <- matrix(1, 3, 3)
  support_b <- calculate_v021_observed_support_features(.support_data())
  expect_identical(support_a, support_b)
  expect_false(identical(absolute_matrix_a, absolute_matrix_b))
  leaked <- .finalized(); leaked$metadata <- list(corrected_oracle = .2)
  expect_error(extract_v021_posterior_stability_features(leaked), "corrected_oracle")
  dictionary <- v021_robustness_feature_dictionary()
  bad <- dictionary; bad$name[[1L]] <- "sign_reversal"
  attr(bad, "dictionary_version") <- v021_robustness_feature_dictionary_version
  expect_error(validate_v021_robustness_feature_dictionary(bad), "truth-derived")
})

test_that("posterior sign stability is not absolute-gLV agreement", {
  features <- extract_v021_posterior_stability_features(.finalized())
  expect_equal(features$values$psp, .98)
  expect_equal(features$values$lfsr, .02)
  expect_false(any(c("absolute_zero_to_projected_nonzero", "sign_reversal",
                     "posterior_vs_truth", "absolute_A") %in%
                   v021_robustness_feature_dictionary()$name))
  expect_error(assert_v021_truth_free_schema(
    list(features = features, sign_reversal = TRUE), "V021-04"), "sign_reversal")
})
