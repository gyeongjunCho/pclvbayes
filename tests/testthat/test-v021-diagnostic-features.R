source(testthat::test_path("../../benchmarks/mtist/v021_diagnostic_features.R"))

v021_feature_fixture <- function() list(
  dataset_id = "dataset-37", pair_id = "species_0|species_1", task_id = 1L,
  direction_index = 1L, target = "species_0", source = "species_1", seed = 20362803L,
  posterior_mean = -0.25, posterior_median = -0.24, posterior_sd = 0.08,
  posterior_interval_lower = -0.41, posterior_interval_upper = -0.09,
  posterior_sign_probability = 0.995, p_sign2 = 0.01, lfsr = 0.005,
  rhat = 1.001, ess_bulk = 900, ess_tail = 740, divergences = 0L,
  treedepth_hits = 0L, ebfmi_min = 0.78, chain_sign_agreement = TRUE,
  diagnostic_class = "converged", interaction_identifiable = TRUE,
  residual_identifiable = TRUE, interaction_identifiability_class = "stable",
  residual_identifiability_class = "stable", bayesian_eligible = TRUE,
  residual_regime_disagreement = FALSE, kfold_attempted = TRUE,
  kfold_completed = TRUE, kfold_folds_ok = 5L, kfold_folds_failed = 0L,
  elpd_available = TRUE, aggregate_elpd = -12.5,
  n_pairs = 24L
)

test_that("V021-04 schema is versioned, ordered, unique, and reproducible", {
  first <- v021_diagnostic_feature_schema()
  second <- v021_diagnostic_feature_schema()
  expect_identical(first, second)
  expect_identical(attr(first, "schema_version"), "v021_diagnostic_features_v1")
  expect_identical(anyDuplicated(first$name), 0L)
  expect_identical(first$name[1:8], c(
    "schema_version", "dataset_id", "pair_id", "task_id", "direction_index",
    "target", "source", "seed"))
  expect_setequal(unique(first$group), c(
    "identity", "posterior_evidence", "sampler_diagnostics", "identifiability",
    "predictive_completion", "design_temporal", "abundance_sparsity", "alr_cap_exposure"))
})

test_that("schema validation rejects duplicates, unsupported types, and incomplete metadata", {
  schema <- v021_diagnostic_feature_schema()
  duplicate <- rbind(schema, schema[1, , drop = FALSE])
  attr(duplicate, "schema_version") <- v021_diagnostic_feature_schema_version
  expect_error(validate_v021_diagnostic_feature_schema(duplicate), "unique")
  bad_type <- schema; bad_type$type[[2L]] <- "object"
  expect_error(validate_v021_diagnostic_feature_schema(bad_type), "Unsupported")
  bad_provenance <- schema; bad_provenance$deterministic_provenance[[2L]] <- ""
  expect_error(validate_v021_diagnostic_feature_schema(bad_provenance), "complete")
  bad_missing <- schema; bad_missing$allowed_missingness[[2L]] <- "invented"
  expect_error(validate_v021_diagnostic_feature_schema(bad_missing), "missingness")
})

test_that("extraction is deterministic, non-mutating, and preserves canonical identity", {
  source <- v021_feature_fixture()
  before <- serialize(source, NULL, version = 3)
  design <- data.frame(subject = c("a", "a", "a", "b", "b"), time = c(0, 1, 3, 0, 2))
  a <- build_v021_diagnostic_feature_record(source, design, list(chains = 4L, iter_sampling = 2000L))
  b <- build_v021_diagnostic_feature_record(source, design, list(chains = 4L, iter_sampling = 2000L))
  expect_identical(a, b)
  expect_identical(serialize(source, NULL, version = 3), before)
  expect_identical(a$values$task_id, 1L)
  expect_identical(a$values$direction_index, 1L)
  expect_identical(a$values$target, "species_0")
  expect_identical(a$values$source, "species_1")
  expect_identical(a$values$seed, 20362803L)
  expect_equal(a$values$posterior_interval_width, 0.32)
  expect_identical(a$values$nominal_retained_draws, 8000L)
  expect_false(identical(a$values$nominal_retained_draws, a$values$bulk_ess))
})

test_that("diagnostic and predictive fields are copied without computation", {
  source <- v021_feature_fixture()
  record <- build_v021_diagnostic_feature_record(source)
  fields <- c("rhat", "bulk_ess", "tail_ess", "divergence_count",
              "treedepth_saturation_count", "ebfmi_min", "chain_sign_agreement",
              "diagnostic_class", "bayesian_eligible", "kfold_completed",
              "aggregate_elpd")
  expected <- c("rhat", "ess_bulk", "ess_tail", "divergences", "treedepth_hits",
                "ebfmi_min", "chain_sign_agreement", "diagnostic_class",
                "bayesian_eligible", "kfold_completed", "aggregate_elpd")
  for (i in seq_along(fields)) expect_identical(record$values[[fields[[i]]]], source[[expected[[i]]]])
})

test_that("truth, calibrated, and withholding inputs fail before extraction", {
  forbidden <- list(
    truth_sign = -1, truth_matrix = matrix(0, 2, 2), absolute_zero = TRUE,
    same_nonzero_sign = TRUE, calibrated_probability = 0.8,
    sign_withheld = TRUE, sign_reportable = FALSE
  )
  for (name in names(forbidden)) {
    source <- v021_feature_fixture(); source[[name]] <- forbidden[[name]]
    expect_error(build_v021_diagnostic_feature_record(source),
                 "Truth, calibration, or withholding", info = name)
  }
  nested <- v021_feature_fixture()
  nested$ignored_payload <- list(truth_sign = -1)
  expect_error(build_v021_diagnostic_feature_record(nested),
               "Truth, calibration, or withholding")
})

test_that("structured missingness never becomes zero or false", {
  source <- v021_feature_fixture()
  source$posterior_sd <- NA_real_
  source$posterior_sd_missing_state <- "indeterminate"
  source$aggregate_elpd <- NA_real_
  source$aggregate_elpd_missing_state <- "not_executed"
  record <- build_v021_diagnostic_feature_record(source)
  expect_true(is.na(record$values$posterior_sd))
  expect_identical(record$missingness$posterior_sd, "indeterminate")
  expect_true(is.na(record$values$aggregate_elpd))
  expect_identical(record$missingness$aggregate_elpd, "not_executed")
  expect_false(identical(record$values$posterior_sd, 0))
  expect_false(identical(record$values$aggregate_elpd, 0))
})

test_that("legacy unsupported fields and ALR-cap exposure remain explicit", {
  record <- build_v021_diagnostic_feature_record(v021_feature_fixture())
  expect_identical(record$missingness$zero_fraction, "absent_from_legacy_fixture")
  expect_true(is.na(record$values$zero_fraction))
  expect_false(record$values$alr_cap_exposure_available)
  expect_identical(record$missingness$alr_cap_exposure_available, "observed")
  expect_true(all(vapply(record$values[c("alr_cap_hit_count", "alr_cap_candidate_count",
                                         "alr_cap_hit_fraction")], is.na, logical(1))))
  expect_true(all(record$missingness[c("alr_cap_hit_count", "alr_cap_candidate_count",
                                       "alr_cap_hit_fraction")] == "absent_from_legacy_fixture"))
  expect_identical(record$values$alr_cap_unavailable_reason,
                   "canonical_pre_cap_accounting_not_retained")
})

test_that("design fixture formulas are deterministic and truth-free", {
  design <- data.frame(subject = c("a", "a", "a", "b", "b"), time = c(0, 1, 3, 0, 2))
  record <- build_v021_diagnostic_feature_record(v021_feature_fixture(), design)
  expect_identical(record$values$n_subjects, 2L)
  expect_identical(record$values$n_timepoints, 5L)
  expect_true(record$values$irregular_time)
  expect_equal(record$values$median_dt, 2)
  expect_equal(record$values$min_dt, 1)
  expect_equal(record$values$max_dt, 2)
})

test_that("record validation enforces types, states, and observed-value semantics", {
  record <- build_v021_diagnostic_feature_record(v021_feature_fixture())
  wrong_type <- record; wrong_type$values$task_id <- 1
  expect_error(validate_v021_diagnostic_feature_record(wrong_type), "type")
  wrong_state <- record; wrong_state$missingness$posterior_sd <- "scientific_zero"
  expect_error(validate_v021_diagnostic_feature_record(wrong_state), "missingness")
  false_observed <- record; false_observed$missingness$zero_fraction <- "observed"
  expect_error(validate_v021_diagnostic_feature_record(false_observed), "cannot be missing")
})
