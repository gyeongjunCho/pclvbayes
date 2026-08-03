# Benchmark-only truth-isolation helpers for ROADMAP V021-01.

v021_inference_spec_fields <- c(
  "physeq", "taxa_vec", "subject_col", "time_col", "config",
  "dataset_id", "task_id", "direction_index", "target", "source", "seed"
)

v021_required_inference_spec_fields <- c(
  "physeq", "taxa_vec", "subject_col", "time_col", "config",
  "task_id", "direction_index", "target", "source", "seed"
)

v021_inference_config_fields <- c(
  "nz_partner_min_frac", "min_unique_times", "min_pairs", "chains",
  "iter_warmup", "iter_sampling", "seed", "init", "adapt_delta",
  "max_treedepth", "progress", "n_workers_outer", "n_workers_kfold",
  "kfold_K", "kfold_R", "kfold_seed"
)

v021_runtime_context_fields <- c(
  "meta_df", "sm_mat", "eps", "min_pairs", "zero_mode_alr", "minpos_alpha",
  "minpos_base", "eps_fixed", "lib_eps_c", "rest_floor_frac", "smooth_scale",
  "alr_spline_df", "alr_spline_spar", "alr_spline_cv", "nz_partner_min_frac",
  "max_retries", "chains", "iter_warmup", "iter_sampling", "adapt_delta",
  "max_treedepth", "metric", "init", "seed", "quiet", "silent_sampler",
  "n_workers_kfold_eff", "kfold_K", "kfold_R", "kfold_seed", "use_pathfinder_init",
  "pf_num_paths", "pf_draws", "pf_history_size", "pf_max_lbfgs_iters",
  "pf_psis_resample", "mod_exe_file"
)

v021_inference_result_fields <- c(
  "posterior_coefficients", "posterior_summaries", "posterior_intervals",
  "posterior_sign_probabilities", "significance_decisions",
  "diagnostic_classes", "bayesian_eligibility", "kfold_results",
  "elpd_results", "matrices", "masks",
  "feature_inputs", "stan_data", "execution", "status",
  "unavailable_values", "zero_placeholders"
)

v021_retained_fixture_schema <- "v021_retained_inference_fixture_v1"

v021_truth_free_analysis_schema <- "v021_truth_free_analysis_specification_v1"
v021_finalized_inference_schema <- "v021_finalized_inference_record_v1"
v021_post_inference_truth_schema <- "v021_post_inference_truth_table_v1"

v021_prohibited_preinference_fields <- c(
  "absolute_a", "glv_a", "glv_a_target_source", "truth", "truth_sign",
  "expected_sign", "corrected_oracle", "oracle_sign", "structural_zero",
  "sign_reversal", "glv_to_oracle_class", "posterior_vs_truth",
  "posterior_vs_oracle", "posterior_vs_glv", "empirical_false_sign",
  "coverage", "withholding_truth_label", "absolute_zero",
  "absolute_zero_to_projected_nonzero", "projected_nonzero",
  "truth_derived_susceptibility", "truth_derived_eligibility",
  "truth_derived_priority", "truth_derived_seed"
)

.v021_copy <- function(x) unserialize(serialize(x, NULL, version = 3))

.v021_truth_name <- function(x) {
  normalized <- tolower(x)
  normalized %in% v021_prohibited_preinference_fields | grepl(
    "(^|[._])(ground[._]?truth|absolute[._]?a|truth[._]?(matrix|path|loader|coefficient|sign|label)|benchmark[._]?(outcome|label))($|[._])",
    normalized, perl = TRUE
  )
}

.v021_prohibited_field_paths <- function(x, path = "object") {
  found <- character()
  if (is.data.frame(x) || is.list(x)) {
    nms <- names(x)
    if (!is.null(nms)) {
      bad <- which(.v021_truth_name(nms))
      if (length(bad)) found <- c(found, paste0(path, "$", nms[bad]))
    }
    for (i in seq_along(x)) {
      child <- if (!is.null(nms) && nzchar(nms[[i]]))
        paste0(path, "$", nms[[i]]) else paste0(path, "[[", i, "]]")
      found <- c(found, .v021_prohibited_field_paths(x[[i]], child))
    }
  }
  unique(found)
}

assert_v021_truth_free_schema <- function(x, path = "object") {
  offending <- .v021_prohibited_field_paths(x, path)
  if (length(offending))
    stop("Prohibited pre-inference truth field(s): ", paste(offending, collapse = ", "))
  invisible(TRUE)
}

.v021_assert_names <- function(x, allowed, required = character(), path) {
  nms <- names(x)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms))
    stop(path, " must have unique, non-empty names.")
  unexpected <- setdiff(nms, allowed)
  if (length(unexpected))
    stop(path, " has unexpected field(s): ", paste(unexpected, collapse = ", "))
  missing <- setdiff(required, nms)
  if (length(missing))
    stop(path, " lacks required field(s): ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

.v021_assert_plain_value <- function(x, path, scan_truth_names = TRUE) {
  if (is.function(x) || is.environment(x) || typeof(x) == "externalptr" ||
      typeof(x) == "weakref" || typeof(x) == "symbol" || typeof(x) == "language")
    stop("Forbidden reference or executable object at ", path, ".")
  if (is.data.frame(x)) {
    if (scan_truth_names) assert_v021_truth_free_schema(x, path)
    if (length(setdiff(names(attributes(x)), c("names", "row.names", "class"))))
      stop("Arbitrary attributes are forbidden at ", path, ".")
    for (nm in names(x)) .v021_assert_plain_value(x[[nm]], paste0(path, "$", nm), scan_truth_names)
    return(invisible(TRUE))
  }
  if (is.object(x))
    stop("Unsupported S3/S4 object at ", path, ".")
  if (!is.null(attributes(x))) {
    allowed <- if (is.matrix(x) || is.array(x)) c("dim", "dimnames") else
      if (is.list(x) || is.atomic(x)) "names" else character()
    if (length(setdiff(names(attributes(x)), allowed)))
      stop("Arbitrary attributes are forbidden at ", path, ".")
  }
  if (is.list(x)) {
    if (is.null(names(x)) && length(x))
      stop("Unnamed nested lists are forbidden at ", path, ".")
    if (scan_truth_names && length(x)) assert_v021_truth_free_schema(x, path)
    for (nm in names(x))
      .v021_assert_plain_value(x[[nm]], paste0(path, "$", nm), scan_truth_names)
  }
  invisible(TRUE)
}

.v021_assert_physeq <- function(x) {
  if (!methods::is(x, "phyloseq") || !identical(as.character(class(x)), "phyloseq"))
    stop("physeq must be an approved phyloseq object.")
  attrs <- names(attributes(x))
  allowed_attrs <- c(methods::slotNames(x), "class")
  if (length(setdiff(attrs, allowed_attrs)))
    stop("physeq has arbitrary attributes.")
  sample_names <- colnames(as.data.frame(phyloseq::sample_data(x), stringsAsFactors = FALSE))
  taxonomy <- phyloseq::tax_table(x, errorIfNULL = FALSE)
  taxonomy_names <- if (is.null(taxonomy)) character() else colnames(as(taxonomy, "matrix"))
  if (any(.v021_truth_name(c(sample_names, taxonomy_names))))
    stop("physeq contains truth-bearing annotations.")
  invisible(TRUE)
}

.v021_validate_config <- function(x) {
  if (!is.list(x)) stop("config must be a named list.")
  if (!length(x)) return(invisible(TRUE))
  .v021_assert_names(x, v021_inference_config_fields, path = "config")
  for (nm in names(x)) {
    value <- x[[nm]]
    .v021_assert_plain_value(value, paste0("config$", nm))
    if (length(value) != 1L || !is.atomic(value) || is.matrix(value) || is.array(value))
      stop("config$", nm, " must be a plain scalar.")
  }
  invisible(TRUE)
}

build_v021_inference_spec <- function(study, config = list(), task) {
  if (!is.list(study) || is.object(study)) stop("study must be a plain source list.")
  if (!all(c("physeq", "taxa") %in% names(study)))
    stop("study must provide physeq and taxa.")
  if (!is.data.frame(task) || nrow(task) != 1L)
    stop("task must be one directed task row.")
  task_fields <- c("dataset_id", "task_id", "direction_index", "target", "source", "seed")
  required_task <- setdiff(task_fields, "dataset_id")
  if (!all(required_task %in% names(task))) stop("task lacks direction identity fields.")
  .v021_assert_physeq(study$physeq)
  .v021_validate_config(config)
  selected <- intersect(task_fields, names(task))
  spec <- c(list(
    physeq = .v021_copy(study$physeq),
    taxa_vec = as.character(.v021_copy(study$taxa)),
    subject_col = "subject", time_col = "time",
    config = .v021_copy(config)
  ), lapply(task[selected], function(value) .v021_copy(value[[1L]])))
  validate_v021_inference_spec(spec)
  spec
}

validate_v021_inference_spec <- function(x) {
  if (!is.list(x) || is.object(x)) stop("inference specification must be a plain list.")
  .v021_assert_names(
    x, v021_inference_spec_fields, v021_required_inference_spec_fields,
    "inference specification"
  )
  if (any(.v021_truth_name(names(x)))) stop("Truth-bearing inference field.")
  .v021_assert_physeq(x$physeq)
  if (!is.character(x$taxa_vec) || !length(x$taxa_vec) || anyNA(x$taxa_vec) ||
      any(!nzchar(x$taxa_vec)) || anyDuplicated(x$taxa_vec)) stop("Invalid taxa_vec.")
  if (!all(x$taxa_vec %in% phyloseq::taxa_names(x$physeq)))
    stop("taxa_vec is not contained in physeq.")
  for (nm in c("subject_col", "time_col", "target", "source"))
    if (!is.character(x[[nm]]) || length(x[[nm]]) != 1L || is.na(x[[nm]]) || !nzchar(x[[nm]]))
      stop(nm, " must be one non-missing string.")
  if (identical(x$target, x$source) || !all(c(x$target, x$source) %in% x$taxa_vec))
    stop("target and source must be distinct selected taxa.")
  for (nm in c("task_id", "direction_index", "seed"))
    if (length(x[[nm]]) != 1L || !is.numeric(x[[nm]]) || !is.finite(x[[nm]]))
      stop(nm, " must be one finite numeric identity.")
  if ("dataset_id" %in% names(x))
    .v021_assert_plain_value(x$dataset_id, "dataset_id")
  .v021_validate_config(x$config)
  allowed_attrs <- c("names")
  if (length(setdiff(names(attributes(x)), allowed_attrs)))
    stop("Inference specification has arbitrary attributes.")
  invisible(TRUE)
}

build_v021_confirmation_target <- function(stage_b_row) {
  fields <- c("task_id", "direction_index", "target", "source", "seed")
  if (!is.data.frame(stage_b_row) || nrow(stage_b_row) != 1L ||
      !all(fields %in% names(stage_b_row))) stop("Invalid Stage B direction row.")
  out <- stage_b_row[fields]
  rownames(out) <- NULL
  .v021_copy(out)
}

build_v021_runtime_context <- function(ctx) {
  if (!is.list(ctx) || is.object(ctx)) stop("Runtime context must be a plain list.")
  missing <- setdiff(v021_runtime_context_fields, names(ctx))
  if (length(missing))
    stop("Runtime context lacks approved field(s): ", paste(missing, collapse = ", "))
  if (!is.null(ctx$pair_builder))
    stop("Custom pair builders are outside the approved Core runtime context.")
  out <- lapply(ctx[v021_runtime_context_fields], .v021_copy)
  for (nm in names(out))
    .v021_assert_plain_value(out[[nm]], paste0("runtime_context$", nm))
  validate_v021_runtime_context(out)
  out
}

validate_v021_runtime_context <- function(ctx) {
  if (!is.list(ctx) || is.object(ctx)) stop("Runtime context must be a plain list.")
  .v021_assert_names(ctx, v021_runtime_context_fields, v021_runtime_context_fields,
                     "runtime context")
  for (nm in names(ctx))
    .v021_assert_plain_value(ctx[[nm]], paste0("runtime_context$", nm))
  if (!is.data.frame(ctx$meta_df) || !is.matrix(ctx$sm_mat) ||
      !is.numeric(ctx$sm_mat) || nrow(ctx$meta_df) != ncol(ctx$sm_mat))
    stop("Runtime inference data are structurally invalid.")
  taxa <- rownames(ctx$sm_mat)
  samples <- colnames(ctx$sm_mat)
  if (is.null(taxa) || anyNA(taxa) || any(!nzchar(taxa)) || anyDuplicated(taxa))
    stop("Runtime inference matrix lacks unique taxon rows.")
  if (is.null(samples) || anyNA(samples) || any(!nzchar(samples)) ||
      anyDuplicated(samples))
    stop("Runtime inference matrix lacks unique sample columns.")
  if (!"Sample" %in% names(ctx$meta_df))
    stop("Runtime metadata lacks canonical sample identifiers.")
  metadata_samples <- as.character(ctx$meta_df$Sample)
  if (anyNA(metadata_samples) || any(!nzchar(metadata_samples)) ||
      anyDuplicated(metadata_samples))
    stop("Runtime metadata lacks unique sample identifiers.")
  if (!setequal(metadata_samples, samples))
    stop("Runtime metadata and matrix sample identifier sets disagree.")
  sample_index <- match(metadata_samples, samples)
  if (length(sample_index) != nrow(ctx$meta_df) || anyNA(sample_index) ||
      anyDuplicated(sample_index))
    stop("Runtime metadata-to-matrix sample alignment is incomplete or not one-to-one.")
  if (!is.character(ctx$mod_exe_file) || length(ctx$mod_exe_file) != 1L ||
      is.na(ctx$mod_exe_file) || !nzchar(ctx$mod_exe_file))
    stop("Runtime executable identity is invalid.")
  invisible(TRUE)
}

make_v021_confirmation_fit_closure <- function(spec, runtime_context, fit_direction) {
  validate_v021_inference_spec(spec)
  validate_v021_runtime_context(runtime_context)
  if (!identical(as.character(spec$taxa_vec), rownames(runtime_context$sm_mat)))
    stop("Inference taxon declarations disagree with runtime matrix rows.")
  if (!is.function(fit_direction)) stop("fit_direction must be a function.")
  closure_env <- new.env(parent = baseenv())
  closure_env$spec <- .v021_copy(spec)
  closure_env$runtime_context <- runtime_context
  closure_env$fit_direction <- fit_direction
  fn <- function() {
    result <- fit_direction(
      target = spec$target, partner = spec$source, ctx = runtime_context,
      seed_override = as.integer(spec$seed), progress_local = "none"
    )
    result
  }
  environment(fn) <- closure_env
  fn
}

.v021_group_columns <- function(grouping_keys) {
  if (!length(grouping_keys)) return(c(dataset = NA_character_, interaction_matrix = NA_character_))
  if (is.null(names(grouping_keys)) || any(!names(grouping_keys) %in% c("dataset", "interaction_matrix")) ||
      anyDuplicated(names(grouping_keys)))
    stop("grouping_keys must be a named character vector using dataset and/or interaction_matrix.")
  out <- c(dataset = NA_character_, interaction_matrix = NA_character_)
  out[names(grouping_keys)] <- grouping_keys
  out
}

build_pair_group_manifest <- function(task_table, grouping_keys = character(), split_spec, seed = 1L) {
  required <- c("task_id", "direction_index", "target", "source", "seed")
  if (!is.data.frame(task_table) || !all(required %in% names(task_table)))
    stop("task_table lacks required direction identifiers.")
  if (!is.list(split_spec) || !identical(sort(names(split_spec)),
      sort(c("grouping_level", "allowed_partitions", "assignments"))))
    stop("split_spec must contain only grouping_level, allowed_partitions, and assignments.")
  level <- split_spec$grouping_level
  if (!is.character(level) || length(level) != 1L ||
      !level %in% c("unordered_pair", "dataset", "interaction_matrix"))
    stop("Invalid grouping_level.")
  allowed <- split_spec$allowed_partitions
  if (!is.character(allowed) || !length(allowed) || anyNA(allowed) ||
      any(!nzchar(allowed)) || anyDuplicated(allowed)) stop("Invalid allowed_partitions.")
  assignments <- split_spec$assignments
  if (!is.data.frame(assignments) || !identical(names(assignments), c("group_id", "partition")))
    stop("assignments must have exactly group_id and partition columns.")
  if (anyNA(assignments) || anyDuplicated(assignments$group_id))
    stop("Duplicate, incomplete, or conflicting assignments.")
  if (any(!assignments$partition %in% allowed)) stop("Assignment uses a disallowed partition.")
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed)) stop("Invalid manifest seed.")
  if (anyDuplicated(paste(task_table$task_id, task_table$direction_index, sep = "\r")))
    stop("Duplicate task/direction identifiers.")
  columns <- .v021_group_columns(grouping_keys)
  needed <- unname(columns[!is.na(columns)])
  if (!all(needed %in% names(task_table))) stop("Required grouping metadata is missing.")
  allowed_task_fields <- unique(c(required, needed, "pair_key"))
  if (length(setdiff(names(task_table), allowed_task_fields)))
    stop("task_table contains fields outside the manifest schema.")

  out <- task_table[intersect(allowed_task_fields, names(task_table))]
  out$dataset_group <- if (is.na(columns[["dataset"]])) NA_character_ else
    as.character(out[[columns[["dataset"]]]])
  out$matrix_group <- if (is.na(columns[["interaction_matrix"]])) NA_character_ else
    as.character(out[[columns[["interaction_matrix"]]]])
  pair_taxa <- vapply(seq_len(nrow(out)), function(i)
    paste(sort(c(as.character(out$target[[i]]), as.character(out$source[[i]]))), collapse = "|"),
    character(1))
  derived_pair <- vapply(seq_len(nrow(out)), function(i) {
    components <- c(out$matrix_group[[i]], out$dataset_group[[i]], pair_taxa[[i]])
    paste(components[!is.na(components) & nzchar(components)], collapse = "::")
  }, character(1))
  if ("pair_key" %in% names(out)) {
    if (anyNA(out$pair_key) || any(!nzchar(as.character(out$pair_key)))) stop("Invalid pair_key.")
    pair_map <- split(as.character(out$pair_key), derived_pair)
    if (any(vapply(pair_map, function(z) length(unique(z)) != 1L, logical(1))))
      stop("Conflicting pair_key values.")
    reverse_pair_map <- split(derived_pair, as.character(out$pair_key))
    if (any(vapply(reverse_pair_map, function(z) length(unique(z)) != 1L, logical(1))))
      stop("A pair_key identifies multiple unordered pairs.")
    out$pair_group <- as.character(out$pair_key)
  } else out$pair_group <- derived_pair
  if ((level == "dataset" && all(is.na(out$dataset_group))) ||
      (level == "interaction_matrix" && all(is.na(out$matrix_group))))
    stop("Grouping level metadata is missing.")
  group_id <- switch(level, unordered_pair = out$pair_group,
                     dataset = out$dataset_group, interaction_matrix = out$matrix_group)
  if (anyNA(group_id) || any(!nzchar(group_id))) stop("Incomplete grouping metadata.")
  missing_groups <- setdiff(unique(group_id), assignments$group_id)
  extra_groups <- setdiff(assignments$group_id, unique(group_id))
  if (length(missing_groups) || length(extra_groups))
    stop("Assignments are incomplete or contain unknown groups.")
  out$grouping_level <- level
  out$assignment_group <- group_id
  out$partition <- assignments$partition[match(group_id, assignments$group_id)]
  out$manifest_seed <- as.integer(seed)
  order_fields <- c("matrix_group", "dataset_group", "pair_group", "direction_index")
  out <- out[do.call(order, c(out[order_fields], list(na.last = TRUE))), , drop = FALSE]
  rownames(out) <- NULL
  assert_no_pair_split_leakage(out)
  out
}

assert_no_pair_split_leakage <- function(manifest) {
  required <- c("task_id", "direction_index", "target", "source", "pair_group", "partition")
  if (!is.data.frame(manifest) || !all(required %in% names(manifest)) || anyNA(manifest[required]))
    stop("Invalid or incomplete manifest.")
  canonical <- vapply(seq_len(nrow(manifest)), function(i)
    paste(sort(c(as.character(manifest$target[[i]]), as.character(manifest$source[[i]]))), collapse = "|"),
    character(1))
  pair_partition <- split(manifest$partition, manifest$pair_group)
  if (any(vapply(pair_partition, function(x) length(unique(x)) != 1L, logical(1))))
    stop("Unordered pair crosses partitions.")
  pair_taxa <- split(canonical, manifest$pair_group)
  if (any(vapply(pair_taxa, function(x) length(unique(x)) != 1L, logical(1))))
    stop("Pair group identifies conflicting unordered pairs.")
  rows <- split(seq_len(nrow(manifest)), manifest$pair_group)
  complete <- vapply(rows, function(ix) {
    if (length(ix) != 2L || length(unique(manifest$task_id[ix])) != 1L) return(FALSE)
    identical(as.character(manifest$target[ix]), rev(as.character(manifest$source[ix])))
  }, logical(1))
  if (any(!complete)) stop("Manifest contains an incomplete unordered pair.")
  invisible(TRUE)
}

.v021_validate_retained_artifact <- function(x, require_finalized = FALSE) {
  allowed <- c("artifact_schema", "artifact_state", "execution_state", "inference_results", ".v021_finalized")
  required <- c("artifact_schema", "artifact_state", "execution_state", "inference_results")
  if (!is.list(x) || is.object(x)) stop("Inference artifact must be a plain structured list.")
  if (length(setdiff(names(attributes(x)), "names")))
    stop("Inference artifact has arbitrary attributes.")
  .v021_assert_names(x, allowed, required, "inference artifact")
  if (!identical(x$artifact_schema, v021_retained_fixture_schema))
    stop("Unsupported inference artifact schema.")
  if (!identical(x$artifact_state, "retained") || !x$execution_state %in% c("completed", "retained"))
    stop("Inference artifact is unfinished, failed, incomplete, or unavailable.")
  if (!is.list(x$inference_results) || is.object(x$inference_results))
    stop("inference_results must be a plain list.")
  if (length(setdiff(names(attributes(x$inference_results)), "names")))
    stop("inference_results has arbitrary attributes.")
  .v021_assert_names(x$inference_results, v021_inference_result_fields,
                     v021_inference_result_fields, "inference_results")
  if (any(.v021_truth_name(names(x$inference_results)))) stop("Truth-bearing inference result field.")
  for (nm in names(x$inference_results))
    .v021_assert_plain_value(x$inference_results[[nm]], paste0("inference_results$", nm))
  if (require_finalized && !identical(x$.v021_finalized, TRUE))
    stop("Truth join requires an explicitly finalized retained inference artifact.")
  if (!require_finalized && ".v021_finalized" %in% names(x))
    stop("Artifact is already finalized.")
  invisible(TRUE)
}

finalize_v021_inference_artifact <- function(x) {
  .v021_validate_retained_artifact(x, require_finalized = FALSE)
  out <- .v021_copy(x)
  out$.v021_finalized <- TRUE
  out
}

classify_absolute_truth_outcome <- function(posterior_sign, absolute_A) {
  unavailable <- function(reason) list(
    absolute_truth_outcome = NA_character_,
    absolute_truth_outcome_status = "unavailable",
    absolute_truth_outcome_reason = list(code = reason)
  )
  if (length(posterior_sign) != 1L || !is.numeric(posterior_sign) ||
      is.na(posterior_sign) || !is.finite(posterior_sign))
    return(unavailable("posterior_sign_unavailable"))
  if (!posterior_sign %in% c(-1, 1)) return(unavailable("posterior_sign_indeterminate"))
  if (length(absolute_A) != 1L || !is.numeric(absolute_A) ||
      is.na(absolute_A) || !is.finite(absolute_A))
    return(unavailable("absolute_A_unavailable"))
  outcome <- if (absolute_A == 0) "absolute_zero" else if (posterior_sign == sign(absolute_A))
    "same_nonzero_sign" else "opposite_nonzero_sign"
  list(absolute_truth_outcome = outcome,
       absolute_truth_outcome_status = "classified",
       absolute_truth_outcome_reason = NULL)
}

join_v021_truth_post_inference <- function(inference_artifact, truth_artifact) {
  .v021_validate_retained_artifact(inference_artifact, require_finalized = TRUE)
  if (!is.list(truth_artifact) || is.object(truth_artifact) || is.null(names(truth_artifact)))
    stop("truth_artifact must be a named plain list.")
  .v021_assert_plain_value(truth_artifact, "truth_artifact", scan_truth_names = FALSE)
  out <- .v021_copy(inference_artifact)
  out$evaluation <- .v021_copy(truth_artifact)
  out
}

v021_analysis_specification_fields <- c(
  "analysis_schema", "dataset_id", "observed_data_id", "taxa_order",
  "task_manifest", "subject_universe", "time_metadata",
  "preprocessing_config", "posterior_config", "kfold_config", "public_seed",
  "kfold_seed", "direction_seeds", "scorer_name", "provenance"
)

build_v021_truth_free_analysis_specification <- function(
    dataset_id, observed_data_id, taxa_order, task_manifest, subject_universe,
    time_metadata, preprocessing_config, posterior_config, kfold_config,
    public_seed, kfold_seed, scorer_name, provenance) {
  spec <- list(
    analysis_schema = v021_truth_free_analysis_schema,
    dataset_id = as.character(dataset_id),
    observed_data_id = as.character(observed_data_id),
    taxa_order = as.character(taxa_order),
    task_manifest = .v021_copy(task_manifest),
    subject_universe = sort(unique(as.character(subject_universe)), method = "radix"),
    time_metadata = .v021_copy(time_metadata),
    preprocessing_config = .v021_copy(preprocessing_config),
    posterior_config = .v021_copy(posterior_config),
    kfold_config = .v021_copy(kfold_config),
    public_seed = as.integer(public_seed),
    kfold_seed = as.integer(kfold_seed),
    direction_seeds = as.integer(task_manifest$seed),
    scorer_name = as.character(scorer_name),
    provenance = .v021_copy(provenance)
  )
  validate_v021_truth_free_analysis_specification(spec)
  spec
}

validate_v021_truth_free_analysis_specification <- function(x) {
  if (!is.list(x) || is.object(x))
    stop("Truth-free analysis specification must be a plain list.")
  .v021_assert_names(x, v021_analysis_specification_fields,
                     v021_analysis_specification_fields, "analysis specification")
  assert_v021_truth_free_schema(x, "analysis specification")
  for (nm in names(x))
    .v021_assert_plain_value(x[[nm]], paste0("analysis specification$", nm))
  if (!identical(x$analysis_schema, v021_truth_free_analysis_schema))
    stop("Unsupported truth-free analysis schema.")
  scalar_text <- c("dataset_id", "observed_data_id", "scorer_name")
  if (any(vapply(x[scalar_text], function(z)
    !is.character(z) || length(z) != 1L || is.na(z) || !nzchar(z), logical(1))))
    stop("Analysis string identities must be non-empty scalars.")
  if (!is.character(x$taxa_order) || length(x$taxa_order) < 2L ||
      anyNA(x$taxa_order) || any(!nzchar(x$taxa_order)) || anyDuplicated(x$taxa_order))
    stop("Analysis taxa order must contain unique taxa.")
  required_task <- c("task_id", "direction_index", "target", "source", "seed")
  if (!is.data.frame(x$task_manifest) || !all(required_task %in% names(x$task_manifest)) ||
      !nrow(x$task_manifest) || anyDuplicated(x$task_manifest$direction_index) ||
      anyDuplicated(paste(x$task_manifest$task_id, x$task_manifest$direction_index)))
    stop("Analysis task manifest has invalid or ambiguous identities.")
  if (any(!x$task_manifest$target %in% x$taxa_order) ||
      any(!x$task_manifest$source %in% x$taxa_order) ||
      any(x$task_manifest$target == x$task_manifest$source))
    stop("Analysis task directions disagree with taxa order.")
  if (!identical(as.integer(x$direction_seeds), as.integer(x$task_manifest$seed)))
    stop("Direction seeds disagree with the task manifest.")
  if (!is.character(x$subject_universe) || !length(x$subject_universe) ||
      anyNA(x$subject_universe) || any(!nzchar(x$subject_universe)) ||
      anyDuplicated(x$subject_universe) ||
      !identical(x$subject_universe, sort(x$subject_universe, method = "radix")))
    stop("Subject universe must be canonical and non-empty.")
  if (!is.data.frame(x$time_metadata) ||
      !all(c("subject", "time") %in% names(x$time_metadata)))
    stop("Time metadata must contain subject and time.")
  for (nm in c("public_seed", "kfold_seed"))
    if (length(x[[nm]]) != 1L || is.na(x[[nm]]) || x[[nm]] < 0L)
      stop(nm, " must be one non-negative integer.")
  invisible(TRUE)
}

v021_finalized_inference_record_fields <- c(
  "record_schema", "finalized", "finalized_timestamp", "dataset_id",
  "taxa_order", "task_id", "direction_index", "target", "source",
  "direction_seed", "terminal_status", "posterior_summaries", "psp_lfsr",
  "diagnostics", "predictive_eligibility", "pair_specific_elpd",
  "failure_information", "seed_split_provenance"
)

build_v021_finalized_inference_record <- function(
    analysis_specification, task, terminal_status, posterior_summaries = list(),
    psp_lfsr = list(), diagnostics = list(), predictive_eligibility = NA,
    pair_specific_elpd = list(), failure_information = list(),
    seed_split_provenance = list(),
    finalized_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  validate_v021_truth_free_analysis_specification(analysis_specification)
  if (!is.data.frame(task) || nrow(task) != 1L)
    stop("Finalized record task must be one manifest row.")
  index <- match(task$direction_index, analysis_specification$task_manifest$direction_index)
  required <- c("task_id", "direction_index", "target", "source", "seed")
  if (is.na(index) || !all(required %in% names(task)) ||
      !identical(as.list(task[1L, required, drop = FALSE]),
                 as.list(analysis_specification$task_manifest[index, required, drop = FALSE])))
    stop("Finalized record task identity disagrees with the analysis manifest.")
  record <- list(
    record_schema = v021_finalized_inference_schema,
    finalized = TRUE,
    finalized_timestamp = as.character(finalized_timestamp),
    dataset_id = analysis_specification$dataset_id,
    taxa_order = analysis_specification$taxa_order,
    task_id = as.integer(task$task_id),
    direction_index = as.integer(task$direction_index),
    target = as.character(task$target), source = as.character(task$source),
    direction_seed = as.integer(task$seed), terminal_status = terminal_status,
    posterior_summaries = .v021_copy(posterior_summaries),
    psp_lfsr = .v021_copy(psp_lfsr), diagnostics = .v021_copy(diagnostics),
    predictive_eligibility = predictive_eligibility,
    pair_specific_elpd = .v021_copy(pair_specific_elpd),
    failure_information = .v021_copy(failure_information),
    seed_split_provenance = .v021_copy(seed_split_provenance)
  )
  validate_v021_finalized_inference_record(record)
  record
}

validate_v021_finalized_inference_record <- function(x) {
  if (!is.list(x) || is.object(x)) stop("Finalized inference record must be a plain list.")
  .v021_assert_names(x, v021_finalized_inference_record_fields,
                     v021_finalized_inference_record_fields, "finalized inference record")
  assert_v021_truth_free_schema(x, "finalized inference record")
  for (nm in names(x))
    .v021_assert_plain_value(x[[nm]], paste0("finalized inference record$", nm))
  if (!identical(x$record_schema, v021_finalized_inference_schema) ||
      !identical(x$finalized, TRUE) || !is.character(x$finalized_timestamp) ||
      length(x$finalized_timestamp) != 1L || is.na(x$finalized_timestamp) ||
      !nzchar(x$finalized_timestamp))
    stop("Inference record lacks an explicit finalized state.")
  if (!x$terminal_status %in% c("completed", "failed", "skipped", "unavailable"))
    stop("Inference record lacks a valid terminal status.")
  if (!is.character(x$taxa_order) || anyDuplicated(x$taxa_order) ||
      !all(c(x$target, x$source) %in% x$taxa_order) || identical(x$target, x$source))
    stop("Inference record direction or taxa order is invalid.")
  if (length(x$task_id) != 1L || length(x$direction_index) != 1L ||
      length(x$direction_seed) != 1L || anyNA(c(x$task_id, x$direction_index,
                                                x$direction_seed)))
    stop("Inference record identity is incomplete.")
  invisible(TRUE)
}

.v021_sign_with_tolerance <- function(x, tolerance) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) return(NA_integer_)
  if (abs(x) <= tolerance) 0L else as.integer(sign(x))
}

classify_v021_glv_to_oracle <- function(absolute_A, corrected_oracle,
                                         tolerance = 1e-10) {
  a_sign <- .v021_sign_with_tolerance(absolute_A, tolerance)
  oracle_sign <- .v021_sign_with_tolerance(corrected_oracle, tolerance)
  if (is.na(a_sign) || is.na(oracle_sign)) return(NA_character_)
  if (a_sign == 0L && oracle_sign == 0L) return("both_zero")
  if (a_sign == 0L) return("absolute_zero_to_projected_nonzero")
  if (oracle_sign == 0L) return("absolute_nonzero_to_projected_zero")
  if (a_sign == oracle_sign) "same_sign" else "sign_reversal"
}

build_v021_post_inference_truth_table <- function(
    inference_records, absolute_A, oracle_table, zero_tolerance = 1e-10) {
  if (!is.list(inference_records) || !length(inference_records))
    stop("At least one finalized inference record is required.")
  invisible(lapply(inference_records, validate_v021_finalized_inference_record))
  identities <- vapply(inference_records, function(x)
    paste(x$dataset_id, x$task_id, x$direction_index, x$target, x$source, sep = "\r"),
    character(1))
  if (anyDuplicated(identities)) stop("Finalized inference task identity is duplicated or ambiguous.")
  datasets <- unique(vapply(inference_records, `[[`, character(1), "dataset_id"))
  taxa <- lapply(inference_records, `[[`, "taxa_order")
  if (length(datasets) != 1L || !all(vapply(taxa, identical, logical(1), taxa[[1L]])))
    stop("Inference records have inconsistent dataset identity or taxa order.")
  taxa <- taxa[[1L]]
  if (!is.matrix(absolute_A) || !is.numeric(absolute_A) || any(!is.finite(absolute_A)) ||
      !identical(rownames(absolute_A), taxa) || !identical(colnames(absolute_A), taxa))
    stop("Absolute interaction matrix does not match the finalized taxa order.")
  oracle_required <- c("dataset_id", "task_id", "direction_index", "target", "source",
                       "corrected_oracle")
  if (!is.data.frame(oracle_table) || !all(oracle_required %in% names(oracle_table)))
    stop("Post-inference oracle table lacks required identity fields.")
  oracle_keys <- do.call(paste, c(oracle_table[oracle_required[1:5]], sep = "\r"))
  record_keys <- vapply(inference_records, function(x)
    paste(x$dataset_id, x$task_id, x$direction_index, x$target, x$source, sep = "\r"),
    character(1))
  if (anyDuplicated(oracle_keys) || anyNA(match(record_keys, oracle_keys)))
    stop("Post-inference oracle identity is duplicated, ambiguous, or incomplete.")
  matched <- oracle_table[match(record_keys, oracle_keys), , drop = FALSE]
  optional_oracle <- c(
    state_contrast_min = NA_real_, state_contrast_max = NA_real_,
    state_contrast_mean = NA_real_, state_contrast_median = NA_real_,
    state_sign_crossing = NA
  )
  for (nm in names(optional_oracle))
    if (!nm %in% names(matched)) matched[[nm]] <- optional_oracle[[nm]]
  rows <- lapply(seq_along(inference_records), function(i) {
    record <- inference_records[[i]]
    absolute <- as.numeric(absolute_A[record$target, record$source])
    oracle <- as.numeric(matched$corrected_oracle[[i]])
    posterior_sign <- record$psp_lfsr$posterior_sign
    if (is.null(posterior_sign)) posterior_sign <- NA_real_
    oracle_sign <- .v021_sign_with_tolerance(oracle, zero_tolerance)
    posterior_sign <- .v021_sign_with_tolerance(posterior_sign, 0)
    absolute_outcome <- classify_absolute_truth_outcome(posterior_sign, absolute)
    interval <- record$posterior_summaries$interval
    coverage <- if (is.numeric(interval) && length(interval) == 2L &&
                    all(is.finite(interval)))
      absolute >= min(interval) && absolute <= max(interval) else NA
    empirical_false_sign <- if (is.na(posterior_sign) || absolute == 0) NA else
      posterior_sign != sign(absolute)
    data.frame(
      truth_schema = v021_post_inference_truth_schema,
      dataset_id = record$dataset_id, task_id = record$task_id,
      direction_index = record$direction_index, target = record$target,
      source = record$source, absolute_A = absolute,
      structural_zero = identical(absolute, 0), corrected_oracle = oracle,
      oracle_sign = oracle_sign,
      state_contrast_min = as.numeric(matched$state_contrast_min[[i]]),
      state_contrast_max = as.numeric(matched$state_contrast_max[[i]]),
      state_contrast_mean = as.numeric(matched$state_contrast_mean[[i]]),
      state_contrast_median = as.numeric(matched$state_contrast_median[[i]]),
      state_sign_crossing = as.logical(matched$state_sign_crossing[[i]]),
      glv_to_oracle_class = classify_v021_glv_to_oracle(
        absolute, oracle, zero_tolerance),
      posterior_vs_oracle = if (is.na(posterior_sign) || is.na(oracle_sign))
        "not_reportable" else if (oracle_sign == 0L) "oracle_zero" else
        if (posterior_sign == oracle_sign) "agree" else "disagree",
      posterior_vs_absolute_A = absolute_outcome$absolute_truth_outcome,
      coverage = coverage, empirical_false_sign = empirical_false_sign,
      calibration_outcome = NA_character_,
      terminal_status = record$terminal_status,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
