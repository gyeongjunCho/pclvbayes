# Benchmark-only truth-free robustness contract for ROADMAP V021-04.

v021_observed_support_features_schema <- "v021_observed_support_features_v1"
v021_posterior_stability_features_schema <- "v021_posterior_stability_features_v1"
v021_perturbation_stability_features_schema <- "v021_perturbation_stability_features_v1"
v021_truth_free_robustness_features_schema <- "v021_truth_free_robustness_features_v1"
v021_robustness_feature_dictionary_version <- "v021_robustness_feature_dictionary_v1"
v021_robustness_feature_artifact_schema <- "v021_robustness_feature_artifact_v1"
v021_perturbation_schema <- "v021_truth_free_perturbation_v1"

v021_robustness_missing_reasons <- c(
  "not_applicable", "inference_not_eligible", "perturbation_failed",
  "zero_denominator", "identity_mismatch", "not_executed", "unavailable"
)

.v021_feature_row <- function(name, type, group, stage, formula, units,
                              valid_range, missingness_rule, denominator,
                              stability_direction = "unknown_before_calibration",
                              required = TRUE) {
  data.frame(
    name, type, group, availability_stage = stage, formula, units, valid_range,
    missingness_rule, denominator, stability_direction,
    required_low_cost = required, stringsAsFactors = FALSE)
}

v021_robustness_feature_dictionary <- function() {
  rows <- list(
    .v021_feature_row("eligible_subject_count", "integer", "observed_support", "pre_fit", "count unique eligible subject IDs", "subjects", "[1,Inf)", "zero_denominator when none", "eligible subjects"),
    .v021_feature_row("eligible_observation_count", "integer", "observed_support", "pre_fit", "count eligible observation rows", "observations", "[1,Inf)", "zero_denominator when none", "eligible observations"),
    .v021_feature_row("canonical_transition_count", "integer", "observed_support", "pre_fit", "count valid adjacent subject-local transitions", "transitions", "[0,Inf)", "observed even when zero", "candidate adjacent intervals"),
    .v021_feature_row("min_transitions_per_subject", "integer", "observed_support", "pre_fit", "minimum valid transitions by eligible subject", "transitions", "[0,Inf)", "zero_denominator when no subjects", "eligible subjects"),
    .v021_feature_row("median_transitions_per_subject", "double", "observed_support", "pre_fit", "median valid transitions by eligible subject", "transitions", "[0,Inf)", "zero_denominator when no subjects", "eligible subjects"),
    .v021_feature_row("max_transitions_per_subject", "integer", "observed_support", "pre_fit", "maximum valid transitions by eligible subject", "transitions", "[0,Inf)", "zero_denominator when no subjects", "eligible subjects"),
    .v021_feature_row("min_time_gap", "double", "observed_support", "pre_fit", "minimum valid adjacent time difference", "declared time units", "(0,Inf)", "zero_denominator when no valid transitions", "valid adjacent transitions"),
    .v021_feature_row("median_time_gap", "double", "observed_support", "pre_fit", "median valid adjacent time difference", "declared time units", "(0,Inf)", "zero_denominator when no valid transitions", "valid adjacent transitions"),
    .v021_feature_row("max_time_gap", "double", "observed_support", "pre_fit", "maximum valid adjacent time difference", "declared time units", "(0,Inf)", "zero_denominator when no valid transitions", "valid adjacent transitions"),
    .v021_feature_row("time_gap_cv", "double", "observed_support", "pre_fit", "sd(dt)/mean(dt)", "ratio", "[0,Inf)", "zero_denominator when fewer than two gaps or mean is zero", "valid adjacent transitions"),
    .v021_feature_row("source_prevalence", "double", "observed_support", "pre_fit", "mean(source abundance > 0)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("target_prevalence", "double", "observed_support", "pre_fit", "mean(target abundance > 0)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("pair_coprevalence", "double", "observed_support", "pre_fit", "mean(source > 0 and target > 0)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("source_abundance_q10", "double", "observed_support", "pre_fit", "10th percentile source abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("source_abundance_median", "double", "observed_support", "pre_fit", "median source abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("source_abundance_q90", "double", "observed_support", "pre_fit", "90th percentile source abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("target_abundance_q10", "double", "observed_support", "pre_fit", "10th percentile target abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("target_abundance_median", "double", "observed_support", "pre_fit", "median target abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("target_abundance_q90", "double", "observed_support", "pre_fit", "90th percentile target abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("source_near_zero_fraction", "double", "observed_support", "pre_fit", "mean(source <= declared near-zero threshold)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("target_near_zero_fraction", "double", "observed_support", "pre_fit", "mean(target <= declared near-zero threshold)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("both_zero_fraction", "double", "observed_support", "pre_fit", "mean(source == 0 and target == 0)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("source_only_zero_fraction", "double", "observed_support", "pre_fit", "mean(source == 0 and target > 0)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("target_only_zero_fraction", "double", "observed_support", "pre_fit", "mean(target == 0 and source > 0)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("rest_support_min", "double", "observed_support", "pre_fit", "minimum rest abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("rest_support_median", "double", "observed_support", "pre_fit", "median rest abundance", "relative abundance", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("weak_rest_support_fraction", "double", "observed_support", "pre_fit", "mean(rest <= declared weak-rest threshold)", "fraction", "[0,1]", "zero_denominator when no eligible observations", "eligible observations"),
    .v021_feature_row("subject_observation_cv", "double", "observed_support", "pre_fit", "sd(observation count by subject)/mean(count)", "ratio", "[0,Inf)", "zero_denominator when fewer than two subjects or mean is zero", "eligible subjects"),
    .v021_feature_row("effective_observation_count", "integer", "observed_support", "pre_fit", "observations represented by valid transition outcomes", "observations", "[0,Inf)", "observed even when zero", "eligible observations"),

    .v021_feature_row("posterior_mean", "double", "posterior_stability", "post_fit", "canonical posterior mean of a_ij", "coefficient", "(-Inf,Inf)", "unavailable when no finalized posterior", "posterior draws"),
    .v021_feature_row("posterior_median", "double", "posterior_stability", "post_fit", "canonical posterior median of a_ij", "coefficient", "(-Inf,Inf)", "unavailable when no finalized posterior", "posterior draws"),
    .v021_feature_row("posterior_sd", "double", "posterior_stability", "post_fit", "canonical posterior SD of a_ij", "coefficient", "[0,Inf)", "unavailable when no finalized posterior", "posterior draws"),
    .v021_feature_row("posterior_interval_width", "double", "posterior_stability", "post_fit", "upper 95% bound minus lower 95% bound", "coefficient", "[0,Inf)", "unavailable when interval absent", "posterior interval"),
    .v021_feature_row("posterior_distance_zero", "double", "posterior_stability", "post_fit", "abs(posterior mean)", "coefficient", "[0,Inf)", "unavailable when mean absent", "posterior draws"),
    .v021_feature_row("posterior_distance_zero_sd", "double", "posterior_stability", "post_fit", "abs(posterior mean)/posterior SD", "ratio", "[0,Inf)", "zero_denominator when posterior SD is zero", "posterior SD"),
    .v021_feature_row("psp", "double", "posterior_stability", "post_fit", "maintained posterior sign probability", "probability", "[0,1]", "unavailable when inference absent", "posterior draws"),
    .v021_feature_row("lfsr", "double", "posterior_stability", "post_fit", "maintained local false sign rate", "probability", "[0,1]", "unavailable when inference absent", "posterior draws"),
    .v021_feature_row("rhat_max", "double", "posterior_stability", "post_fit", "maximum maintained relevant R-hat", "ratio", "[0,Inf)", "unavailable when diagnostics absent", "diagnostic parameters"),
    .v021_feature_row("bulk_ess_min", "double", "posterior_stability", "post_fit", "minimum maintained relevant bulk ESS", "effective draws", "[0,Inf)", "unavailable when diagnostics absent", "diagnostic parameters"),
    .v021_feature_row("tail_ess_min", "double", "posterior_stability", "post_fit", "minimum maintained relevant tail ESS", "effective draws", "[0,Inf)", "unavailable when diagnostics absent", "diagnostic parameters"),
    .v021_feature_row("divergence_count", "integer", "posterior_stability", "post_fit", "sum divergent transitions", "events", "[0,Inf)", "unavailable when diagnostics absent", "sampled transitions"),
    .v021_feature_row("treedepth_hit_count", "integer", "posterior_stability", "post_fit", "sum maximum-treedepth hits", "events", "[0,Inf)", "unavailable when diagnostics absent", "sampled transitions"),
    .v021_feature_row("ebfmi_min", "double", "posterior_stability", "post_fit", "minimum chain E-BFMI", "ratio", "[0,Inf)", "unavailable when diagnostics absent", "chains"),
    .v021_feature_row("retry_count", "integer", "posterior_stability", "post_fit", "number of sampler retries", "attempts", "[0,Inf)", "unavailable when execution absent", "attempt ledger"),
    .v021_feature_row("diagnostic_class", "character", "posterior_stability", "post_fit", "maintained final diagnostic class", "class", "declared classes", "unavailable when inference absent", "finalized inference"),
    .v021_feature_row("predictive_eligible", "logical", "posterior_stability", "post_fit", "maintained predictive eligibility", "boolean", "TRUE|FALSE", "unavailable when inference absent", "finalized inference"),
    .v021_feature_row("elpd", "double", "posterior_stability", "post_fit", "pair-specific subject-level ELPD sum", "log predictive density", "(-Inf,Inf)", "inference_not_eligible or unavailable", "successful held-out observations"),
    .v021_feature_row("elpd_per_successful_observation", "double", "posterior_stability", "post_fit", "ELPD/n_successful_test_observations", "log predictive density", "(-Inf,Inf)", "zero_denominator when successful count is zero", "successful held-out observations"),
    .v021_feature_row("successful_test_observation_count", "integer", "posterior_stability", "post_fit", "count successful held-out observations", "observations", "[0,Inf)", "not_executed when K-fold absent", "attempted held-out observations"),
    .v021_feature_row("fold_failure_fraction", "double", "posterior_stability", "post_fit", "failed folds/attempted folds", "fraction", "[0,1]", "zero_denominator when no folds attempted", "attempted folds"),

    .v021_feature_row("coefficient_shift_abs", "double", "perturbation_stability", "optional_perturbation", "abs(perturbed mean-canonical mean)", "coefficient", "[0,Inf)", "perturbation_failed or unavailable", "one valid canonical-perturbed comparison", required = FALSE),
    .v021_feature_row("coefficient_shift_relative", "double", "perturbation_stability", "optional_perturbation", "absolute shift/max(abs(canonical mean), declared epsilon)", "ratio", "[0,Inf)", "zero_denominator only when epsilon invalid", "max(abs(canonical mean), epsilon)", required = FALSE),
    .v021_feature_row("coefficient_shift_sd", "double", "perturbation_stability", "optional_perturbation", "absolute shift/canonical posterior SD", "SD units", "[0,Inf)", "zero_denominator when canonical SD is zero", "canonical posterior SD", required = FALSE),
    .v021_feature_row("posterior_sign_agreement", "logical", "perturbation_stability", "optional_perturbation", "sign(canonical mean)==sign(perturbed mean), excluding zero", "boolean", "TRUE|FALSE", "not_applicable when either sign is zero", "one valid comparison", required = FALSE),
    .v021_feature_row("interval_overlap_fraction", "double", "perturbation_stability", "optional_perturbation", "intersection length/min(interval widths)", "fraction", "[0,1]", "zero_denominator for zero-width interval", "smaller interval width", required = FALSE),
    .v021_feature_row("psp_change", "double", "perturbation_stability", "optional_perturbation", "perturbed PSP-canonical PSP", "probability difference", "[-1,1]", "perturbation_failed or unavailable", "one valid comparison", required = FALSE),
    .v021_feature_row("lfsr_change", "double", "perturbation_stability", "optional_perturbation", "perturbed LFSR-canonical LFSR", "probability difference", "[-1,1]", "perturbation_failed or unavailable", "one valid comparison", required = FALSE),
    .v021_feature_row("diagnostic_class_changed", "logical", "perturbation_stability", "optional_perturbation", "canonical diagnostic class differs from perturbed class", "boolean", "TRUE|FALSE", "perturbation_failed or unavailable", "one valid comparison", required = FALSE),
    .v021_feature_row("predictive_eligibility_changed", "logical", "perturbation_stability", "optional_perturbation", "canonical predictive eligibility differs", "boolean", "TRUE|FALSE", "perturbation_failed or unavailable", "one valid comparison", required = FALSE),
    .v021_feature_row("elpd_change_same_outcome", "double", "perturbation_stability", "optional_perturbation", "perturbed ELPD-canonical ELPD for identical outcome and held-out split", "log predictive density", "(-Inf,Inf)", "identity_mismatch unless outcomes and folds match", "identical successful held-out observations", required = FALSE),
    .v021_feature_row("successful_observation_count_change", "integer", "perturbation_stability", "optional_perturbation", "perturbed successful count-canonical count", "observations", "(-Inf,Inf)", "perturbation_failed or unavailable", "attempted held-out observations", required = FALSE)
  )
  dictionary <- do.call(rbind, rows)
  attr(dictionary, "dictionary_version") <- v021_robustness_feature_dictionary_version
  validate_v021_robustness_feature_dictionary(dictionary)
  dictionary
}

validate_v021_robustness_feature_dictionary <- function(dictionary) {
  required <- c("name", "type", "group", "availability_stage", "formula", "units",
                "valid_range", "missingness_rule", "denominator",
                "stability_direction", "required_low_cost")
  if (!is.data.frame(dictionary) || !identical(names(dictionary), required) ||
      !nrow(dictionary) || anyNA(dictionary) || anyDuplicated(dictionary$name) ||
      !identical(attr(dictionary, "dictionary_version"),
                 v021_robustness_feature_dictionary_version))
    stop("Invalid V021 robustness feature dictionary.")
  if (any(!dictionary$type %in% c("integer", "double", "logical", "character")) ||
      any(!dictionary$group %in% c("observed_support", "posterior_stability",
                                   "perturbation_stability")) ||
      any(!nzchar(dictionary$denominator)) ||
      any(.v021_truth_name(dictionary$name)) ||
      any(grepl("false_sign_risk|reversal_probability|unsafe_edge|withholding_score|truth_susceptibility|sign_reversal|absolute_zero_to_projected_nonzero",
                dictionary$name, ignore.case = TRUE)))
    stop("Feature dictionary contains invalid or truth-derived semantics.")
  invisible(TRUE)
}

.v021_typed_na <- function(type) switch(type, integer = NA_integer_, double = NA_real_,
                                        logical = NA, character = NA_character_)

.v021_validate_feature_values <- function(record, group) {
  dictionary <- v021_robustness_feature_dictionary()
  rows <- dictionary[dictionary$group == group, , drop = FALSE]
  if (!is.list(record) || !identical(names(record$values), rows$name) ||
      !identical(names(record$missingness), rows$name))
    stop("Feature values do not match the frozen dictionary ordering.")
  for (i in seq_len(nrow(rows))) {
    value <- record$values[[i]]
    type_ok <- switch(rows$type[[i]], integer = is.integer(value),
      double = is.double(value), logical = is.logical(value), character = is.character(value))
    reason <- record$missingness[[i]]
    if (!type_ok || length(value) != 1L ||
        !is.character(reason) || length(reason) != 1L ||
        (!is.na(reason) && !reason %in% v021_robustness_missing_reasons) ||
        (is.na(value) != !is.na(reason)))
      stop("Invalid typed feature value or missingness for `", rows$name[[i]], "`.")
  }
  invisible(TRUE)
}

validate_v021_observed_support_features <- function(record) {
  if (!is.list(record) ||
      !identical(names(record), c("schema_version", "values", "missingness", "provenance")) ||
      !identical(record$schema_version, v021_observed_support_features_schema))
    stop("Invalid observed-support feature schema.")
  .v021_validate_feature_values(record, "observed_support")
  fractions <- c("source_prevalence", "target_prevalence", "pair_coprevalence",
                 "source_near_zero_fraction", "target_near_zero_fraction",
                 "both_zero_fraction", "source_only_zero_fraction",
                 "target_only_zero_fraction", "weak_rest_support_fraction")
  observed <- fractions[vapply(record$values[fractions], function(x) !is.na(x), logical(1))]
  if (length(observed) && any(unlist(record$values[observed]) < 0 |
                              unlist(record$values[observed]) > 1))
    stop("Observed-support fraction is outside [0,1].")
  assert_v021_truth_free_schema(record, "observed-support features")
  invisible(TRUE)
}

validate_v021_posterior_stability_features <- function(record) {
  if (!is.list(record) ||
      !identical(names(record), c("schema_version", "values", "missingness")) ||
      !identical(record$schema_version, v021_posterior_stability_features_schema))
    stop("Invalid posterior-stability feature schema.")
  .v021_validate_feature_values(record, "posterior_stability")
  for (name in c("psp", "lfsr", "fold_failure_fraction")) {
    value <- record$values[[name]]
    if (!is.na(value) && (value < 0 || value > 1))
      stop(name, " is outside [0,1].")
  }
  assert_v021_truth_free_schema(record, "posterior-stability features")
  invisible(TRUE)
}

.v021_adjacent_transitions <- function(data) {
  rows <- list()
  for (subject in sort(unique(as.character(data$subject)), method = "radix")) {
    idx <- which(as.character(data$subject) == subject)
    idx <- idx[order(data$time[idx], idx)]
    if (length(idx) < 2L) next
    for (j in seq_len(length(idx) - 1L)) {
      a <- idx[[j]]; b <- idx[[j + 1L]]
      valid <- isTRUE(data$eligible[[a]]) && isTRUE(data$eligible[[b]]) &&
        is.finite(data$time[[a]]) && is.finite(data$time[[b]]) &&
        data$time[[b]] > data$time[[a]]
      if (valid) rows[[length(rows) + 1L]] <- data.frame(
        subject = subject, predecessor = a, outcome = b,
        dt = data$time[[b]] - data$time[[a]], stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) data.frame(subject = character(), predecessor = integer(),
                                outcome = integer(), dt = double()) else do.call(rbind, rows)
}

calculate_v021_observed_support_features <- function(data, near_zero_threshold = 1e-6,
                                                       weak_rest_threshold = 0.05) {
  required <- c("subject", "time", "source_abundance", "target_abundance",
                "rest_abundance", "eligible")
  if (!is.numeric(near_zero_threshold) || length(near_zero_threshold) != 1L ||
      !is.finite(near_zero_threshold) || near_zero_threshold < 0 ||
      !is.numeric(weak_rest_threshold) || length(weak_rest_threshold) != 1L ||
      !is.finite(weak_rest_threshold) || weak_rest_threshold < 0 ||
      weak_rest_threshold > 1 ||
      !is.data.frame(data) || !identical(names(data), required) || !nrow(data) ||
      anyNA(data[c("subject", "eligible")]) ||
      any(!is.finite(as.matrix(data[c("source_abundance", "target_abundance",
                                      "rest_abundance")]))) ||
      any(as.matrix(data[c("source_abundance", "target_abundance", "rest_abundance")]) < 0) ||
      any(as.matrix(data[c("source_abundance", "target_abundance", "rest_abundance")]) > 1))
    stop("Invalid observed-support input.")
  assert_v021_truth_free_schema(data, "observed support input")
  eligible <- which(data$eligible & is.finite(data$time))
  x <- data[eligible, , drop = FALSE]
  transitions <- .v021_adjacent_transitions(data)
  subjects <- sort(unique(as.character(x$subject)), method = "radix")
  transition_counts <- table(factor(transitions$subject, levels = subjects))
  obs_counts <- table(factor(as.character(x$subject), levels = subjects))
  q <- function(z, p) as.numeric(stats::quantile(z, p, names = FALSE, type = 8))
  cv <- function(z) if (length(z) < 2L || mean(z) == 0) NA_real_ else stats::sd(z) / mean(z)
  if (!length(eligible)) stop("Observed-support denominator has zero eligible observations.")
  values <- list(
    eligible_subject_count = as.integer(length(subjects)),
    eligible_observation_count = as.integer(length(eligible)),
    canonical_transition_count = as.integer(nrow(transitions)),
    min_transitions_per_subject = as.integer(min(transition_counts)),
    median_transitions_per_subject = as.double(stats::median(transition_counts)),
    max_transitions_per_subject = as.integer(max(transition_counts)),
    min_time_gap = if (nrow(transitions)) min(transitions$dt) else NA_real_,
    median_time_gap = if (nrow(transitions)) stats::median(transitions$dt) else NA_real_,
    max_time_gap = if (nrow(transitions)) max(transitions$dt) else NA_real_,
    time_gap_cv = cv(transitions$dt),
    source_prevalence = mean(x$source_abundance > 0),
    target_prevalence = mean(x$target_abundance > 0),
    pair_coprevalence = mean(x$source_abundance > 0 & x$target_abundance > 0),
    source_abundance_q10 = q(x$source_abundance, .1),
    source_abundance_median = q(x$source_abundance, .5),
    source_abundance_q90 = q(x$source_abundance, .9),
    target_abundance_q10 = q(x$target_abundance, .1),
    target_abundance_median = q(x$target_abundance, .5),
    target_abundance_q90 = q(x$target_abundance, .9),
    source_near_zero_fraction = mean(x$source_abundance <= near_zero_threshold),
    target_near_zero_fraction = mean(x$target_abundance <= near_zero_threshold),
    both_zero_fraction = mean(x$source_abundance == 0 & x$target_abundance == 0),
    source_only_zero_fraction = mean(x$source_abundance == 0 & x$target_abundance > 0),
    target_only_zero_fraction = mean(x$target_abundance == 0 & x$source_abundance > 0),
    rest_support_min = min(x$rest_abundance),
    rest_support_median = stats::median(x$rest_abundance),
    weak_rest_support_fraction = mean(x$rest_abundance <= weak_rest_threshold),
    subject_observation_cv = cv(as.numeric(obs_counts)),
    effective_observation_count = as.integer(nrow(transitions))
  )
  missingness <- lapply(values, function(value)
    if (length(value) == 1L && !is.na(value)) NA_character_ else "zero_denominator")
  out <- list(schema_version = v021_observed_support_features_schema,
              values = values, missingness = missingness,
              provenance = list(near_zero_threshold = near_zero_threshold,
                                weak_rest_threshold = weak_rest_threshold,
                                subject_universe = subjects,
                                transition_rule = "subject_local_adjacent_predecessor_to_outcome_v1"))
  validate_v021_observed_support_features(out)
  out
}

extract_v021_posterior_stability_features <- function(finalized) {
  if (!is.list(finalized) || !isTRUE(finalized$finalized) ||
      !identical(finalized$terminal_state, "completed"))
    stop("Posterior stability requires a finalized completed inference record.")
  assert_v021_truth_free_schema(finalized, "finalized inference")
  required <- c("posterior_mean", "posterior_median", "posterior_sd", "interval_lower",
                "interval_upper", "psp", "lfsr", "rhat_max", "bulk_ess_min",
                "tail_ess_min", "divergence_count", "treedepth_hit_count",
                "ebfmi_min", "retry_count", "diagnostic_class", "predictive_eligible",
                "elpd", "successful_test_observation_count", "folds_attempted", "folds_failed")
  if (length(setdiff(required, names(finalized))))
    stop("Finalized inference lacks posterior-stability fields.")
  values <- finalized[required]
  values$posterior_interval_width <- values$interval_upper - values$interval_lower
  values$posterior_distance_zero <- abs(values$posterior_mean)
  values$posterior_distance_zero_sd <- if (values$posterior_sd > 0)
    abs(values$posterior_mean) / values$posterior_sd else NA_real_
  values$elpd_per_successful_observation <- if (values$successful_test_observation_count > 0)
    values$elpd / values$successful_test_observation_count else NA_real_
  values$fold_failure_fraction <- if (values$folds_attempted > 0)
    values$folds_failed / values$folds_attempted else NA_real_
  values$interval_lower <- values$interval_upper <- values$folds_attempted <- values$folds_failed <- NULL
  order <- v021_robustness_feature_dictionary()$name
  order <- order[v021_robustness_feature_dictionary()$group == "posterior_stability"]
  values <- values[order]
  missingness <- lapply(values, function(value)
    if (!is.na(value)) NA_character_ else "zero_denominator")
  out <- list(schema_version = v021_posterior_stability_features_schema,
              values = values, missingness = missingness)
  validate_v021_posterior_stability_features(out)
  out
}

v021_approved_perturbations <- function(subject_universe = character(), maximum_subject_deletions = 3L) {
  subjects <- sort(unique(as.character(subject_universe)), method = "radix")
  base <- list(
    no_smoothing = list(changes = "preprocessing_only", rule = "use reclosed observed community without spline smoothing where supported", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    smoothing_strength_lower = list(changes = "preprocessing_only", rule = "fixed spline spar minus 0.10 clipped to [0,1]; unavailable for incompatible smoothing mode", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    smoothing_strength_upper = list(changes = "preprocessing_only", rule = "fixed spline spar plus 0.10 clipped to [0,1]; unavailable for incompatible smoothing mode", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    zero_replacement_half = list(changes = "preprocessing_only", rule = "canonical replacement magnitude multiplied by 0.5", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    zero_replacement_double = list(changes = "preprocessing_only", rule = "canonical replacement magnitude multiplied by 2", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    minpos_triplet = list(changes = "preprocessing_only", rule = "replace canonical minpos_base ij with triplet", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    remove_first_observation_per_subject = list(changes = "input_observations", rule = "remove radix/time-ordered first observation of every subject", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    remove_last_observation_per_subject = list(changes = "input_observations", rule = "remove radix/time-ordered last observation of every subject", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    retain_early_half = list(changes = "input_observations", rule = "retain first ceiling(T/2) time-ordered observations per subject", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    retain_late_half = list(changes = "input_observations", rule = "retain last ceiling(T/2) time-ordered observations per subject", response = TRUE, predictor = TRUE, subject_universe = FALSE),
    coarsen_every_second_observation = list(changes = "input_observations", rule = "retain time-ordered subject-local positions 1,3,5,...", response = TRUE, predictor = TRUE, subject_universe = FALSE)
  )
  for (subject in head(subjects, maximum_subject_deletions))
    base[[paste0("delete_subject__", subject)]] <- list(
      changes = "input_observations", rule = paste0("remove subject ID exactly `", subject, "`"),
      response = TRUE, predictor = TRUE,
      subject_universe = TRUE, deleted_subject = subject)
  base
}

build_v021_perturbation_identity <- function(perturbation_name, approved,
                                              dataset_id, task_identity, source, target,
                                              taxa_order, subject_universe,
                                              canonical_configuration_hash) {
  if (!perturbation_name %in% names(approved)) stop("Unapproved perturbation.")
  body <- list(
    perturbation_schema = v021_perturbation_schema,
    perturbation_name = perturbation_name, changes = approved[[perturbation_name]],
    dataset_id = dataset_id, task_identity = task_identity, source = source, target = target,
    taxa_order = as.character(taxa_order),
    subject_universe = sort(unique(as.character(subject_universe)), method = "radix"),
    canonical_configuration_hash = canonical_configuration_hash,
    estimand = "directed_pair_to_rest_a_ij")
  assert_v021_truth_free_schema(body, "perturbation identity")
  c(body, list(perturbation_hash = v021_sha256(body)))
}

compare_v021_perturbation_stability <- function(canonical, perturbed,
                                                 relative_epsilon = 1e-8) {
  identity_fields <- c("dataset_id", "task_identity", "source", "target", "taxa_order",
                       "estimand")
  if (!all(vapply(identity_fields, function(x)
    identical(canonical[[x]], perturbed[[x]]), logical(1))))
    stop("Canonical and perturbation identity mismatch.")
  assert_v021_truth_free_schema(list(canonical = canonical, perturbed = perturbed),
                                "perturbation comparison")
  if (!isTRUE(perturbed$available)) return(list(
    schema_version = v021_perturbation_stability_features_schema,
    values = NULL, validity_status = "unavailable",
    missingness_reason = "perturbation_failed"))
  shift <- abs(perturbed$posterior_mean - canonical$posterior_mean)
  width_c <- canonical$interval_upper - canonical$interval_lower
  width_p <- perturbed$interval_upper - perturbed$interval_lower
  overlap <- max(0, min(canonical$interval_upper, perturbed$interval_upper) -
                   max(canonical$interval_lower, perturbed$interval_lower))
  same_elpd <- identical(canonical$transformed_outcome_hash,
                         perturbed$transformed_outcome_hash) &&
    identical(canonical$heldout_split_hash, perturbed$heldout_split_hash)
  values <- list(
    coefficient_shift_abs = shift,
    coefficient_shift_relative = shift / max(abs(canonical$posterior_mean), relative_epsilon),
    coefficient_shift_sd = if (canonical$posterior_sd > 0) shift / canonical$posterior_sd else NA_real_,
    posterior_sign_agreement = if (canonical$posterior_mean == 0 || perturbed$posterior_mean == 0)
      NA else sign(canonical$posterior_mean) == sign(perturbed$posterior_mean),
    interval_overlap_fraction = if (min(width_c, width_p) > 0)
      min(1, overlap / min(width_c, width_p)) else NA_real_,
    psp_change = perturbed$psp - canonical$psp,
    lfsr_change = perturbed$lfsr - canonical$lfsr,
    diagnostic_class_changed = !identical(perturbed$diagnostic_class, canonical$diagnostic_class),
    predictive_eligibility_changed = !identical(perturbed$predictive_eligible,
                                                 canonical$predictive_eligible),
    elpd_change_same_outcome = if (same_elpd) perturbed$elpd - canonical$elpd else NA_real_,
    successful_observation_count_change = as.integer(
      perturbed$successful_test_observation_count - canonical$successful_test_observation_count))
  list(schema_version = v021_perturbation_stability_features_schema,
       values = values, validity_status = "valid",
       missingness_reason = if (same_elpd) NA_character_ else "identity_mismatch")
}

build_v021_robustness_feature_artifact <- function(manifest, task_ordinal,
                                                    observed_support,
                                                    posterior_stability = NULL,
                                                    perturbation_stability = list(),
                                                    code_provenance) {
  validate_v021_execution_manifest(manifest)
  i <- match(as.integer(task_ordinal), manifest$tasks$task_ordinal)
  if (is.na(i)) stop("Unknown robustness-feature task identity.")
  task <- manifest$tasks[i, , drop = FALSE]
  body <- list(
    artifact_schema = v021_robustness_feature_artifact_schema,
    combined_schema = v021_truth_free_robustness_features_schema,
    dictionary_version = v021_robustness_feature_dictionary_version,
    manifest_hash = manifest$manifest_hash,
    configuration_hash = manifest$configuration_hash,
    dataset_id = manifest$dataset_id,
    directed_task_id = task$directed_task_id[[1L]],
    source = task$source[[1L]], target = task$target[[1L]],
    taxa_order_hash = v021_sha256(manifest$taxa_order),
    subject_universe_hash = v021_sha256(sort(unique(as.character(
      observed_support$provenance$subject_universe)), method = "radix")),
    feature_validity_status = "valid",
    missingness_reason = NA_character_,
    code_provenance = code_provenance,
    calculation_provenance = list(dictionary = v021_robustness_feature_dictionary_version),
    observed_support = observed_support,
    posterior_stability = posterior_stability,
    perturbation_stability = perturbation_stability)
  assert_v021_truth_free_schema(body, "robustness feature artifact")
  c(body, list(artifact_hash = v021_sha256(body)))
}

validate_v021_robustness_feature_artifact <- function(artifact, manifest) {
  validate_v021_execution_manifest(manifest)
  required <- c("artifact_schema", "combined_schema", "dictionary_version", "manifest_hash",
                "configuration_hash", "dataset_id", "directed_task_id", "source", "target",
                "taxa_order_hash", "subject_universe_hash", "feature_validity_status",
                "missingness_reason", "code_provenance", "calculation_provenance",
                "observed_support", "posterior_stability", "perturbation_stability",
                "artifact_hash")
  if (!is.list(artifact) || !identical(names(artifact), required) ||
      !identical(artifact$artifact_schema, v021_robustness_feature_artifact_schema) ||
      !identical(artifact$combined_schema, v021_truth_free_robustness_features_schema) ||
      !identical(artifact$dictionary_version, v021_robustness_feature_dictionary_version) ||
      !identical(artifact$manifest_hash, manifest$manifest_hash) ||
      !identical(artifact$configuration_hash, manifest$configuration_hash) ||
      !identical(artifact$taxa_order_hash, v021_sha256(manifest$taxa_order)))
    stop("Robustness feature artifact identity or schema mismatch.")
  task <- manifest$tasks[match(artifact$directed_task_id,
                               manifest$tasks$directed_task_id), , drop = FALSE]
  if (nrow(task) != 1L || !identical(artifact$source, task$source[[1L]]) ||
      !identical(artifact$target, task$target[[1L]]))
    stop("Robustness feature task orientation mismatch.")
  assert_v021_truth_free_schema(artifact, "robustness feature artifact")
  validate_v021_observed_support_features(artifact$observed_support)
  if (!is.null(artifact$posterior_stability))
    validate_v021_posterior_stability_features(artifact$posterior_stability)
  expected <- v021_sha256(artifact[setdiff(required, "artifact_hash")])
  if (!identical(artifact$artifact_hash, expected))
    stop("Robustness feature artifact hash mismatch.")
  invisible(TRUE)
}

write_v021_robustness_feature_artifact_atomic <- function(artifact, path, manifest) {
  validator <- function(x) validate_v021_robustness_feature_artifact(x, manifest)
  .v021_atomic_write_validated_rds(artifact, path, validator)
}
