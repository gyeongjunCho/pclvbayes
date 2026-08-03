# Benchmark-only diagnostic-feature contract for ROADMAP V021-04.

v021_diagnostic_feature_schema_version <- "v021_diagnostic_features_v1"

v021_feature_missing_states <- c(
  "observed", "unavailable", "not_applicable", "not_executed", "failed",
  "indeterminate", "absent_from_legacy_fixture"
)

.v021_feature_definitions <- list(
  # Identity is copied from the canonical task/direction records.
  c("schema_version", "identity", "character", "required", "identifier", "schema", "constant schema identifier"),
  c("dataset_id", "identity", "character", "required", "identifier", "retained_record", "copied without transformation"),
  c("pair_id", "identity", "character", "nullable", "identifier", "V021-03 manifest", "copied without transformation"),
  c("task_id", "identity", "integer", "required", "identifier", "canonical task record", "copied without transformation"),
  c("direction_index", "identity", "integer", "required", "identifier", "canonical direction record", "copied without transformation"),
  c("target", "identity", "character", "required", "taxon identifier", "canonical direction record", "copied without transformation"),
  c("source", "identity", "character", "required", "taxon identifier", "canonical direction record", "copied without transformation"),
  c("seed", "identity", "integer", "required", "seed", "canonical direction record", "copied without transformation"),

  c("posterior_mean", "posterior_evidence", "double", "nullable", "coefficient", "retained inference", "copied without transformation"),
  c("posterior_median", "posterior_evidence", "double", "nullable", "coefficient", "retained inference", "copied without transformation"),
  c("posterior_sd", "posterior_evidence", "double", "nullable", "coefficient", "retained inference", "copied without transformation"),
  c("posterior_interval_lower", "posterior_evidence", "double", "nullable", "coefficient", "retained inference", "copied without transformation"),
  c("posterior_interval_upper", "posterior_evidence", "double", "nullable", "coefficient", "retained inference", "copied without transformation"),
  c("posterior_interval_width", "posterior_evidence", "double", "nullable", "coefficient", "retained inference", "posterior_interval_upper - posterior_interval_lower"),
  c("posterior_sign_probability", "posterior_evidence", "double", "nullable", "probability", "retained inference", "copied without transformation"),
  c("psp_two_sided", "posterior_evidence", "double", "nullable", "probability", "retained inference p_sign2", "copied without transformation"),
  c("lfsr", "posterior_evidence", "double", "nullable", "probability", "retained inference", "copied without transformation"),
  c("nominal_retained_draws", "posterior_evidence", "integer", "nullable", "draws", "approved run metadata", "chains multiplied by iter_sampling; not ESS"),

  c("rhat", "sampler_diagnostics", "double", "nullable", "ratio", "retained inference", "copied without transformation"),
  c("bulk_ess", "sampler_diagnostics", "double", "nullable", "effective draws", "retained inference", "copied without transformation"),
  c("tail_ess", "sampler_diagnostics", "double", "nullable", "effective draws", "retained inference", "copied without transformation"),
  c("divergence_count", "sampler_diagnostics", "integer", "nullable", "events", "retained inference", "copied without transformation"),
  c("treedepth_saturation_count", "sampler_diagnostics", "integer", "nullable", "events", "retained inference", "copied without transformation"),
  c("ebfmi_min", "sampler_diagnostics", "double", "nullable", "ratio", "retained inference", "copied without transformation"),
  c("chain_sign_agreement", "sampler_diagnostics", "logical", "nullable", "boolean", "retained inference", "copied without transformation"),

  c("diagnostic_class", "identifiability", "character", "nullable", "class", "retained inference", "copied without reinterpretation"),
  c("interaction_identifiable", "identifiability", "logical", "nullable", "boolean", "retained inference", "copied without reinterpretation"),
  c("residual_identifiable", "identifiability", "logical", "nullable", "boolean", "retained inference", "copied without reinterpretation"),
  c("interaction_identifiability_class", "identifiability", "character", "nullable", "class", "retained inference", "copied without reinterpretation"),
  c("residual_identifiability_class", "identifiability", "character", "nullable", "class", "retained inference", "copied without reinterpretation"),
  c("bayesian_eligible", "identifiability", "logical", "nullable", "boolean", "retained inference", "copied without reinterpretation"),
  c("residual_allocation_disagreement", "identifiability", "logical", "nullable", "boolean", "retained inference", "copied without reinterpretation"),

  c("kfold_attempted", "predictive_completion", "logical", "nullable", "boolean", "retained inference", "copied; never executed here"),
  c("kfold_completed", "predictive_completion", "logical", "nullable", "boolean", "retained inference", "copied; never executed here"),
  c("kfold_folds_ok", "predictive_completion", "integer", "nullable", "folds", "retained inference", "copied; never recomputed"),
  c("kfold_folds_failed", "predictive_completion", "integer", "nullable", "folds", "retained inference", "copied; never recomputed"),
  c("elpd_available", "predictive_completion", "logical", "nullable", "boolean", "retained inference", "copied; never recomputed"),
  c("aggregate_elpd", "predictive_completion", "double", "nullable", "log predictive density", "retained inference", "copied; never recomputed"),

  c("n_subjects", "design_temporal", "integer", "nullable", "subjects", "approved truth-free design fixture", "number of unique non-missing subject identifiers"),
  c("n_timepoints", "design_temporal", "integer", "nullable", "observations", "approved truth-free design fixture", "number of finite time observations"),
  c("irregular_time", "design_temporal", "logical", "nullable", "boolean", "approved truth-free design fixture", "TRUE when positive within-subject time gaps are not all equal"),
  c("median_dt", "design_temporal", "double", "nullable", "input time units", "approved truth-free design fixture", "median positive within-subject time gap"),
  c("min_dt", "design_temporal", "double", "nullable", "input time units", "approved truth-free design fixture", "minimum positive within-subject time gap"),
  c("max_dt", "design_temporal", "double", "nullable", "input time units", "approved truth-free design fixture", "maximum positive within-subject time gap"),
  c("pair_sample_count", "design_temporal", "integer", "nullable", "analysis rows", "retained inference n_pairs", "copied without transformation"),

  c("zero_fraction", "abundance_sparsity", "double", "nullable", "fraction", "approved truth-free retained metadata", "copied without transformation"),
  c("dominance", "abundance_sparsity", "double", "nullable", "fraction", "approved truth-free retained metadata", "copied without transformation"),
  c("effective_sample_fraction", "abundance_sparsity", "double", "nullable", "fraction", "approved truth-free retained metadata", "copied without transformation"),
  c("usable_pair_fraction", "abundance_sparsity", "double", "nullable", "fraction", "approved truth-free retained metadata", "copied without transformation"),

  c("alr_cap_exposure_available", "alr_cap_exposure", "logical", "required", "boolean", "retained preprocessing metadata", "FALSE unless exact canonical pre-cap accounting is retained"),
  c("alr_cap_hit_count", "alr_cap_exposure", "integer", "nullable", "capped values", "retained preprocessing metadata", "copied only; equality at the canonical cap counts as a hit"),
  c("alr_cap_candidate_count", "alr_cap_exposure", "integer", "nullable", "candidate values", "retained preprocessing metadata", "number of ALR values assessed before canonical capping"),
  c("alr_cap_hit_fraction", "alr_cap_exposure", "double", "nullable", "fraction", "retained preprocessing metadata", "alr_cap_hit_count / alr_cap_candidate_count when denominator is positive"),
  c("alr_cap_unavailable_reason", "alr_cap_exposure", "character", "nullable", "reason", "feature extraction", "structured reason; legacy retained records lack pre-cap accounting")
)

v021_diagnostic_feature_schema <- function() {
  rows <- lapply(.v021_feature_definitions, function(x) data.frame(
    name = x[[1L]], group = x[[2L]], type = x[[3L]],
    nullable = identical(x[[4L]], "nullable"), units = x[[5L]],
    source = x[[6L]], deterministic_provenance = x[[7L]],
    stringsAsFactors = FALSE
  ))
  out <- do.call(rbind, rows)
  out$allowed_missingness <- I(lapply(out$nullable, function(nullable)
    if (nullable) v021_feature_missing_states else "observed"))
  rownames(out) <- NULL
  attr(out, "schema_version") <- v021_diagnostic_feature_schema_version
  validate_v021_diagnostic_feature_schema(out)
  out
}

validate_v021_diagnostic_feature_schema <- function(schema) {
  required <- c("name", "group", "type", "nullable", "units", "source",
                "deterministic_provenance", "allowed_missingness")
  if (!is.data.frame(schema) || !identical(names(schema), required) || !nrow(schema))
    stop("Invalid V021-04 feature schema structure.")
  if (!identical(attr(schema, "schema_version"), v021_diagnostic_feature_schema_version))
    stop("Invalid V021-04 feature schema version.")
  if (anyNA(schema$name) || any(!nzchar(schema$name)) || anyDuplicated(schema$name))
    stop("Feature names must be unique and non-empty.")
  if (any(!schema$type %in% c("character", "integer", "double", "logical")))
    stop("Unsupported feature type.")
  if (anyNA(schema[c("group", "units", "source", "deterministic_provenance")]) ||
      any(!nzchar(unlist(schema[c("group", "units", "source", "deterministic_provenance")]))))
    stop("Feature metadata and provenance must be complete.")
  valid_missing <- vapply(seq_len(nrow(schema)), function(i) {
    states <- schema$allowed_missingness[[i]]
    is.character(states) && length(states) > 0L && !anyNA(states) &&
      !anyDuplicated(states) && all(states %in% v021_feature_missing_states) &&
      (schema$nullable[[i]] || identical(states, "observed"))
  }, logical(1))
  if (!all(valid_missing)) stop("Invalid feature missingness declaration.")
  if (any(.v021_feature_forbidden_name(schema$name)))
    stop("Truth, calibration, or withholding fields are forbidden.")
  invisible(TRUE)
}

.v021_feature_forbidden_name <- function(x) grepl(
  "truth|absolute_[Aa]|same_nonzero_sign|opposite_nonzero_sign|absolute_zero|calibrat|withhold|withheld|sign_reportable",
  x, ignore.case = TRUE, perl = TRUE
)

.v021_assert_truth_free_feature_source <- function(x, path = "retained_record") {
  if (is.function(x) || is.environment(x) || typeof(x) == "externalptr" ||
      typeof(x) == "weakref" || typeof(x) == "language")
    stop("Executable or reference-bearing input is forbidden at ", path, ".")
  if (is.list(x)) {
    nms <- names(x)
    if (length(x) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms))))
      stop("Nested feature-source lists must be uniquely named at ", path, ".")
    if (length(nms) && (anyDuplicated(nms) || any(.v021_feature_forbidden_name(nms))))
      stop("Truth, calibration, or withholding fields are forbidden before feature construction.")
    for (nm in nms) .v021_assert_truth_free_feature_source(x[[nm]], paste0(path, "$", nm))
  }
  invisible(TRUE)
}

.v021_feature_na <- function(type) switch(type,
  character = NA_character_, integer = NA_integer_, double = NA_real_,
  logical = NA, stop("Unsupported feature type."))

.v021_feature_cast <- function(x, type, name) {
  if (length(x) != 1L) stop("Feature `", name, "` must be scalar.")
  value <- switch(type, character = as.character(x), integer = as.integer(x),
                  double = as.numeric(x), logical = as.logical(x))
  if (!is.na(x) && is.na(value)) stop("Feature `", name, "` cannot be cast to ", type, ".")
  value
}

.v021_design_summary <- function(design_fixture) {
  if (is.null(design_fixture)) return(NULL)
  if (!is.data.frame(design_fixture) ||
      !identical(sort(names(design_fixture)), c("subject", "time")))
    stop("design_fixture must contain exactly subject and time.")
  if (anyNA(design_fixture$subject) || any(!is.finite(design_fixture$time)))
    stop("design_fixture subject and time values must be observed and finite.")
  gaps <- unlist(lapply(split(design_fixture$time, design_fixture$subject), function(x) {
    x <- sort(unique(as.numeric(x)))
    if (length(x) < 2L) numeric() else diff(x)
  }), use.names = FALSE)
  gaps <- gaps[is.finite(gaps) & gaps > 0]
  list(
    n_subjects = length(unique(design_fixture$subject)),
    n_timepoints = nrow(design_fixture),
    irregular_time = if (length(gaps)) length(unique(gaps)) > 1L else NA,
    median_dt = if (length(gaps)) stats::median(gaps) else NA_real_,
    min_dt = if (length(gaps)) min(gaps) else NA_real_,
    max_dt = if (length(gaps)) max(gaps) else NA_real_
  )
}

build_v021_diagnostic_feature_record <- function(retained_record,
                                                  design_fixture = NULL,
                                                  run_metadata = NULL) {
  if (!is.list(retained_record) || is.object(retained_record) || is.null(names(retained_record)) ||
      anyDuplicated(names(retained_record))) stop("retained_record must be a uniquely named plain list.")
  .v021_assert_truth_free_feature_source(retained_record)
  if (!is.null(run_metadata) && (!is.list(run_metadata) ||
      !identical(sort(names(run_metadata)), c("chains", "iter_sampling"))))
    stop("run_metadata must contain exactly chains and iter_sampling.")

  schema <- v021_diagnostic_feature_schema()
  values <- setNames(lapply(schema$type, .v021_feature_na), schema$name)
  missingness <- setNames(as.list(rep("absent_from_legacy_fixture", nrow(schema))), schema$name)
  source_map <- c(
    dataset_id = "dataset_id", pair_id = "pair_id", task_id = "task_id",
    direction_index = "direction_index", target = "target", source = "source", seed = "seed",
    posterior_mean = "posterior_mean", posterior_median = "posterior_median",
    posterior_sd = "posterior_sd", posterior_interval_lower = "posterior_interval_lower",
    posterior_interval_upper = "posterior_interval_upper",
    posterior_sign_probability = "posterior_sign_probability", psp_two_sided = "p_sign2",
    lfsr = "lfsr", rhat = "rhat", bulk_ess = "ess_bulk", tail_ess = "ess_tail",
    divergence_count = "divergences", treedepth_saturation_count = "treedepth_hits",
    ebfmi_min = "ebfmi_min", chain_sign_agreement = "chain_sign_agreement",
    diagnostic_class = "diagnostic_class", interaction_identifiable = "interaction_identifiable",
    residual_identifiable = "residual_identifiable",
    interaction_identifiability_class = "interaction_identifiability_class",
    residual_identifiability_class = "residual_identifiability_class",
    bayesian_eligible = "bayesian_eligible",
    residual_allocation_disagreement = "residual_regime_disagreement",
    kfold_attempted = "kfold_attempted", kfold_completed = "kfold_completed",
    kfold_folds_ok = "kfold_folds_ok", kfold_folds_failed = "kfold_folds_failed",
    elpd_available = "elpd_available", aggregate_elpd = "aggregate_elpd",
    pair_sample_count = "n_pairs", zero_fraction = "zero_fraction", dominance = "dominance",
    effective_sample_fraction = "effective_sample_fraction", usable_pair_fraction = "usable_pair_fraction"
  )
  values$schema_version <- v021_diagnostic_feature_schema_version
  missingness$schema_version <- "observed"
  for (feature in names(source_map)) {
    source_name <- source_map[[feature]]
    if (source_name %in% names(retained_record)) {
      value <- retained_record[[source_name]]
      state_name <- paste0(source_name, "_missing_state")
      state <- if (state_name %in% names(retained_record)) retained_record[[state_name]] else
        if (length(value) == 1L && !is.na(value)) "observed" else "unavailable"
      values[[feature]] <- .v021_feature_cast(value, schema$type[match(feature, schema$name)], feature)
      missingness[[feature]] <- as.character(state)
    }
  }
  if (all(missingness[c("posterior_interval_lower", "posterior_interval_upper")] == "observed")) {
    values$posterior_interval_width <- values$posterior_interval_upper - values$posterior_interval_lower
    missingness$posterior_interval_width <- "observed"
  }
  if (!is.null(run_metadata)) {
    if (any(!vapply(run_metadata, function(x) is.numeric(x) && length(x) == 1L &&
                    is.finite(x) && x > 0 && x == as.integer(x), logical(1))))
      stop("run_metadata values must be positive integers.")
    values$nominal_retained_draws <- as.integer(run_metadata$chains * run_metadata$iter_sampling)
    missingness$nominal_retained_draws <- "observed"
  }
  design <- .v021_design_summary(design_fixture)
  if (!is.null(design)) for (feature in names(design)) {
    values[[feature]] <- .v021_feature_cast(design[[feature]], schema$type[match(feature, schema$name)], feature)
    missingness[[feature]] <- if (is.na(values[[feature]])) "not_applicable" else "observed"
  }

  # Canonical preprocessing does not retain pre-cap accounting in legacy records.
  values$alr_cap_exposure_available <- FALSE
  missingness$alr_cap_exposure_available <- "observed"
  values$alr_cap_unavailable_reason <- "canonical_pre_cap_accounting_not_retained"
  missingness$alr_cap_unavailable_reason <- "observed"
  missingness[c("alr_cap_hit_count", "alr_cap_candidate_count", "alr_cap_hit_fraction")] <-
    "absent_from_legacy_fixture"

  out <- list(schema_version = v021_diagnostic_feature_schema_version,
              values = values, missingness = missingness)
  validate_v021_diagnostic_feature_record(out, schema)
  unserialize(serialize(out, NULL, version = 3))
}

validate_v021_diagnostic_feature_record <- function(record,
                                                     schema = v021_diagnostic_feature_schema()) {
  validate_v021_diagnostic_feature_schema(schema)
  if (!is.list(record) || !identical(names(record), c("schema_version", "values", "missingness")) ||
      !identical(record$schema_version, v021_diagnostic_feature_schema_version) ||
      !identical(names(record$values), schema$name) ||
      !identical(names(record$missingness), schema$name))
    stop("Invalid V021-04 feature record structure or ordering.")
  for (i in seq_len(nrow(schema))) {
    name <- schema$name[[i]]
    state <- record$missingness[[name]]
    if (!is.character(state) || length(state) != 1L || is.na(state) ||
        !state %in% schema$allowed_missingness[[i]])
      stop("Invalid missingness for feature `", name, "`.")
    value <- record$values[[name]]
    expected <- schema$type[[i]]
    type_ok <- switch(expected, character = is.character(value), integer = is.integer(value),
                      double = is.double(value), logical = is.logical(value))
    if (!type_ok || length(value) != 1L) stop("Invalid type for feature `", name, "`.")
    if (identical(state, "observed") && is.na(value))
      stop("Observed feature `", name, "` cannot be missing.")
    if (!identical(state, "observed") && !is.na(value))
      stop("Non-observed feature `", name, "` must retain a typed missing value.")
  }
  if (any(.v021_feature_forbidden_name(names(record$values))))
    stop("Truth, calibration, or withholding fields are forbidden.")
  invisible(TRUE)
}
