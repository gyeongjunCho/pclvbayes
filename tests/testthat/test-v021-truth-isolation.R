source(testthat::test_path("../../benchmarks/mtist/v021_truth_isolation.R"))

make_v021_physeq <- function() {
  phyloseq::phyloseq(
    phyloseq::otu_table(
      matrix(c(1, 2, 3, 4), 2, 2,
             dimnames = list(c("a", "b"), c("s1", "s2"))),
      taxa_are_rows = TRUE
    ),
    phyloseq::sample_data(data.frame(
      subject = c("x", "x"), time = c(0, 1), row.names = c("s1", "s2")
    ))
  )
}

make_v021_spec <- function() {
  study <- list(
    physeq = make_v021_physeq(), taxa = c("a", "b"),
    truth = matrix(c(0, 1, -1, 0), 2), truth_path = "not-for-inference",
    truth_loader = function() stop("must not run")
  )
  task <- data.frame(
    dataset_id = "dataset-1", task_id = 1L, direction_index = 1L,
    target = "a", source = "b", seed = 19L
  )
  list(study = study, task = task,
       spec = build_v021_inference_spec(study, list(chains = 1L), task))
}

make_v021_results <- function() {
  list(
    posterior_coefficients = c(a_to_b = -0.25),
    posterior_summaries = data.frame(parameter = "a_ij", mean = -0.25),
    posterior_intervals = matrix(c(-0.4, -0.1), nrow = 1,
                                 dimnames = list("a_ij", c("lower", "upper"))),
    posterior_sign_probabilities = c(negative = 0.99, positive = 0.01),
    significance_decisions = c(a_to_b = TRUE),
    diagnostic_classes = c(a_to_b = "converged"),
    bayesian_eligibility = c(a_to_b = TRUE),
    kfold_results = list(completed = TRUE, folds_ok = 5L),
    elpd_results = list(available = TRUE, value = -12.5),
    matrices = list(posterior = matrix(c(0, -0.25, 0, 0), 2)),
    masks = list(unavailable = matrix(c(FALSE, FALSE, TRUE, FALSE), 2)),
    feature_inputs = data.frame(rhat = 1.001, ess_bulk = 850),
    stan_data = list(N = 2L, y = c(0.2, -0.1)),
    execution = list(state = "completed", attempts = 1L),
    status = list(original_reporting_state = "reportable"),
    unavailable_values = list(stacking = NA_real_),
    zero_placeholders = list(external_matrix = 0)
  )
}

make_v021_artifact <- function(artifact_state = "retained", execution_state = "completed") {
  list(
    artifact_schema = v021_retained_fixture_schema,
    artifact_state = artifact_state,
    execution_state = execution_state,
    inference_results = make_v021_results()
  )
}

make_v021_tasks <- function() {
  data.frame(
    task_id = rep(1:4, each = 2), direction_index = 1:8,
    target = c("a", "b", "a", "c", "a", "b", "a", "c"),
    source = c("b", "a", "c", "a", "b", "a", "c", "a"),
    seed = 101:108,
    dataset_id = rep(c("d1", "d1", "d2", "d2"), each = 2),
    matrix_id = rep(c("m1", "m1", "m2", "m2"), each = 2),
    stringsAsFactors = FALSE
  )
}

make_v021_runtime_source <- function(n_taxa = 2L, n_samples = 2L) {
  ctx <- setNames(rep(list(1), length(v021_runtime_context_fields)),
                  v021_runtime_context_fields)
  taxa <- if (n_taxa == 2L) c("a", "b") else paste0("taxon-", seq_len(n_taxa))
  samples <- paste0("sample-", seq_len(n_samples))
  ctx$meta_df <- data.frame(
    Sample = samples, subject = rep("x", n_samples), time = seq_len(n_samples) - 1L,
    stringsAsFactors = FALSE
  )
  ctx$sm_mat <- matrix(seq_len(n_taxa * n_samples), n_taxa, n_samples,
                       dimnames = list(taxa, samples))
  ctx$mod_exe_file <- "/approved/canonical-model"
  ctx
}

make_v021_runtime <- function() {
  ctx <- make_v021_runtime_source()
  build_v021_runtime_context(ctx)
}

test_that("runtime validation enforces canonical taxa-by-samples orientation", {
  canonical <- make_v021_runtime_source(10L, 150L)
  expect_identical(dim(canonical$sm_mat), c(10L, 150L))
  expect_identical(nrow(canonical$meta_df), 150L)
  expect_silent(validate_v021_runtime_context(canonical))
  expect_silent(build_v021_runtime_context(canonical))

  transposed <- canonical
  transposed$sm_mat <- t(transposed$sm_mat)
  expect_error(validate_v021_runtime_context(transposed), "structurally invalid")

  wrong_dimension <- canonical
  wrong_dimension$meta_df <- wrong_dimension$meta_df[-1L, , drop = FALSE]
  expect_error(validate_v021_runtime_context(wrong_dimension), "structurally invalid")

  wrong_samples <- canonical
  wrong_samples$meta_df$Sample[[1L]] <- "different-sample"
  expect_error(validate_v021_runtime_context(wrong_samples), "disagree")

  reordered_samples <- canonical
  reordered_samples$meta_df <- reordered_samples$meta_df[rev(seq_len(nrow(reordered_samples$meta_df))), ]
  expect_silent(validate_v021_runtime_context(reordered_samples))
  sample_index <- match(reordered_samples$meta_df$Sample,
                        colnames(reordered_samples$sm_mat))
  expect_length(sample_index, nrow(reordered_samples$meta_df))
  expect_false(anyNA(sample_index))
  expect_identical(anyDuplicated(sample_index), 0L)
  expect_identical(colnames(reordered_samples$sm_mat)[sample_index],
                   as.character(reordered_samples$meta_df$Sample))

  genuine_shape_permuted <- make_v021_runtime_source(10L, 150L)
  genuine_shape_permuted$meta_df <- genuine_shape_permuted$meta_df[
    c(121:150, 1:120), , drop = FALSE]
  expect_silent(validate_v021_runtime_context(genuine_shape_permuted))

  metadata_only <- canonical
  metadata_only$meta_df$Sample[[1L]] <- "metadata-only-sample"
  expect_error(validate_v021_runtime_context(metadata_only), "sets disagree")

  matrix_only <- canonical
  colnames(matrix_only$sm_mat)[[1L]] <- "matrix-only-sample"
  expect_error(validate_v021_runtime_context(matrix_only), "sets disagree")

  duplicate_metadata_sample <- canonical
  duplicate_metadata_sample$meta_df$Sample[[2L]] <- duplicate_metadata_sample$meta_df$Sample[[1L]]
  expect_error(validate_v021_runtime_context(duplicate_metadata_sample), "unique sample")

  missing_metadata_sample <- canonical
  missing_metadata_sample$meta_df$Sample[[1L]] <- NA_character_
  expect_error(validate_v021_runtime_context(missing_metadata_sample), "unique sample")

  empty_metadata_sample <- canonical
  empty_metadata_sample$meta_df$Sample[[1L]] <- ""
  expect_error(validate_v021_runtime_context(empty_metadata_sample), "unique sample")

  duplicate_matrix_sample <- canonical
  colnames(duplicate_matrix_sample$sm_mat)[[2L]] <- colnames(duplicate_matrix_sample$sm_mat)[[1L]]
  expect_error(validate_v021_runtime_context(duplicate_matrix_sample), "unique sample")

  missing_matrix_sample <- canonical
  colnames(missing_matrix_sample$sm_mat)[[1L]] <- NA_character_
  expect_error(validate_v021_runtime_context(missing_matrix_sample), "unique sample")

  empty_matrix_sample <- canonical
  colnames(empty_matrix_sample$sm_mat)[[1L]] <- ""
  expect_error(validate_v021_runtime_context(empty_matrix_sample), "unique sample")

  duplicate_taxon <- canonical
  rownames(duplicate_taxon$sm_mat)[[2L]] <- rownames(duplicate_taxon$sm_mat)[[1L]]
  expect_error(validate_v021_runtime_context(duplicate_taxon), "unique taxon")

  missing_taxon <- canonical
  rownames(missing_taxon$sm_mat)[[1L]] <- NA_character_
  expect_error(validate_v021_runtime_context(missing_taxon), "unique taxon")
})

test_that("approved runtime preserves the shared truth-free K-fold seed", {
  runtime <- make_v021_runtime_source()
  runtime$kfold_seed <- 1234L
  approved <- build_v021_runtime_context(runtime)
  expect_identical(approved$kfold_seed, 1234L)
  changed <- runtime
  changed$truth_matrix <- matrix(1, 2, 2)
  expect_identical(build_v021_runtime_context(changed)$kfold_seed,
                   approved$kfold_seed)
})

test_that("inference taxa declarations agree with runtime taxon rows", {
  skip_if_not_installed("phyloseq")
  fixture <- make_v021_spec()
  runtime <- make_v021_runtime()
  fit_stub <- function(...) list()
  environment(fit_stub) <- baseenv()
  expect_silent(make_v021_confirmation_fit_closure(fixture$spec, runtime, fit_stub))
  wrong <- runtime
  rownames(wrong$sm_mat) <- c("a", "different-taxon")
  expect_silent(validate_v021_runtime_context(wrong))
  expect_error(make_v021_confirmation_fit_closure(fixture$spec, wrong, fit_stub),
               "taxon declarations disagree")
})

make_split_spec <- function(level, ids, partitions) {
  list(
    grouping_level = level,
    allowed_partitions = c("development", "threshold_calibration", "locked_evaluation"),
    assignments = data.frame(group_id = ids, partition = partitions,
                             stringsAsFactors = FALSE)
  )
}

test_that("inference specifications use the exact approved top-level and config schemas", {
  skip_if_not_installed("phyloseq")
  fixture <- make_v021_spec()
  expect_identical(names(fixture$spec), v021_inference_spec_fields)
  expect_identical(names(fixture$spec$config), "chains")
  expect_silent(validate_v021_inference_spec(fixture$spec))

  expect_error(validate_v021_inference_spec(c(fixture$spec, list(extra = 1))), "unexpected")
  nested <- fixture$spec
  nested$config$unknown_option <- 1
  expect_error(validate_v021_inference_spec(nested), "unexpected")
})

test_that("truth-bearing and executable values cannot enter inference specifications", {
  skip_if_not_installed("phyloseq")
  spec <- make_v021_spec()$spec
  attacks <- list(
    truth = matrix(0, 2, 2), truth_coefficient = -0.2, truth_sign = -1,
    truth_path = "/secret", benchmark_outcome = "same", truth_loader = function() 1,
    callback = function() 1, environment = new.env(), closure = local({ z <- 1; function() z }),
    complete_study = list(physeq = spec$physeq, truth = matrix(0, 2, 2))
  )
  for (nm in names(attacks)) {
    attacked <- spec
    attacked[[nm]] <- attacks[[nm]]
    expect_error(validate_v021_inference_spec(attacked), info = nm)
  }

  nested_truth <- spec
  nested_truth$config <- list(chains = list(truth_matrix = matrix(0, 2, 2)))
  expect_error(validate_v021_inference_spec(nested_truth))
  disguised_matrix <- spec
  disguised_matrix$config <- list(chains = matrix(0, 1, 1))
  expect_error(validate_v021_inference_spec(disguised_matrix), "plain scalar")
  nested_object <- spec
  nested_object$config <- list(chains = structure(1, class = "unapproved"))
  expect_error(validate_v021_inference_spec(nested_object), "Unsupported")
  nested_s4 <- spec
  nested_s4$config <- list(chains = spec$physeq)
  expect_error(validate_v021_inference_spec(nested_s4), "Unsupported")
  nonempty_environment <- spec
  env <- new.env()
  env$truth <- matrix(0, 2, 2)
  nonempty_environment$config <- list(chains = env)
  expect_error(validate_v021_inference_spec(nonempty_environment), "Forbidden")
  hidden <- spec
  attr(hidden, "truth") <- matrix(0, 2, 2)
  expect_error(validate_v021_inference_spec(hidden), "attributes")
  pointer <- methods::new("externalptr")
  nested_pointer <- spec
  nested_pointer$config <- list(chains = pointer)
  expect_error(validate_v021_inference_spec(nested_pointer))

  annotated <- make_v021_physeq()
  phyloseq::sample_data(annotated)$truth_sign <- c(1, -1)
  bad_study <- list(physeq = annotated, taxa = c("a", "b"))
  expect_error(build_v021_inference_spec(bad_study, list(), make_v021_spec()$task),
               "truth-bearing")
})

test_that("specifications are independent copies and retain no source study", {
  skip_if_not_installed("phyloseq")
  fixture <- make_v021_spec()
  spec_before <- serialize(fixture$spec, NULL)
  fixture$study$truth[,] <- 99
  fixture$study$truth_path <- "changed"
  fixture$study$taxa[[1L]] <- "changed"
  phyloseq::otu_table(fixture$study$physeq)[1, 1] <- 999
  fixture$study$new_truth <- list(sign = 1)
  expect_identical(serialize(fixture$spec, NULL), spec_before)
  expect_false(any(c("study", "truth", "truth_path", "truth_loader") %in% names(fixture$spec)))
  expect_identical(fixture$spec$taxa_vec, c("a", "b"))
})

test_that("separate truth artifacts cannot alter feature inputs or Stan data", {
  artifact <- finalize_v021_inference_artifact(make_v021_artifact())
  feature_before <- artifact$inference_results$feature_inputs
  stan_before <- artifact$inference_results$stan_data
  truth <- list(absolute_A = -0.7, truth_sign = -1)
  joined <- join_v021_truth_post_inference(artifact, truth)
  expect_identical(joined$inference_results$feature_inputs, feature_before)
  expect_identical(joined$inference_results$stan_data, stan_before)
  expect_identical(artifact$inference_results$feature_inputs, feature_before)
  expect_identical(artifact$inference_results$stan_data, stan_before)
})

test_that("confirmation targets and fitting closures are physically truth-free", {
  skip_if_not_installed("phyloseq")
  fixture <- make_v021_spec()
  row <- cbind(fixture$task, diagnostic_class = "converged",
               truth_sign = -1, posterior_truth_agreement = TRUE)
  target <- build_v021_confirmation_target(row)
  expect_identical(names(target), c("task_id", "direction_index", "target", "source", "seed"))
  spec <- build_v021_inference_spec(fixture$study, list(chains = 1L), target)
  fit_stub <- function(target, partner, ctx, seed_override, progress_local)
    list(target = target, partner = partner, seed = seed_override, ctx = ctx,
         .predictive_context = list(pair_in = data.frame(subject = "s1"),
                                    split_seed = 11L, sampling_seed = seed_override))
  environment(fit_stub) <- baseenv()
  runtime <- make_v021_runtime()
  job <- make_v021_confirmation_fit_closure(spec, runtime, fit_stub)
  expect_setequal(ls(environment(job), all.names = TRUE),
                  c("spec", "runtime_context", "fit_direction"))
  expect_identical(parent.env(environment(job)), baseenv())
  expect_false(any(c("study", "stage_b", "row", "truth") %in%
                     ls(environment(job), all.names = TRUE)))
  expect_identical(job()$target, "a")
  expect_identical(job()$partner, "b")
  expect_true(is.list(job()$.predictive_context))
  expect_false(any(c("truth", "absolute_A", "oracle_sign") %in%
                     names(job()$.predictive_context)))
  runtime$truth_matrix <- matrix(1, 2, 2)
  expect_false("truth_matrix" %in% names(environment(job)$runtime_context))
})

test_that("the confirmation runner constructs isolated jobs before fitting", {
  runner <- readLines(testthat::test_path(
    "../../benchmarks/mtist/run_ten_species_confirmation.R"), warn = FALSE)
  target_line <- grep("fit_targets <-", runner, fixed = TRUE)
  spec_line <- grep("inference_specs <-", runner, fixed = TRUE)
  job_line <- grep("fit_jobs <-", runner, fixed = TRUE)
  map_line <- grep("results <- furrr::future_map", runner, fixed = TRUE)
  expect_true(target_line < spec_line && spec_line < job_line && job_line < map_line)
  fitting_block <- runner[map_line:grep("future::plan(future::sequential)", runner, fixed = TRUE)[2L]]
  expect_false(any(grepl("study|stage_b|evaluation_targets|targets\\[|row <-", fitting_block)))
  expect_true(any(grepl("b <- evaluation_targets", runner, fixed = TRUE)))
})

test_that("finalization requires the explicit retained fixture schema and state", {
  expect_error(finalize_v021_inference_artifact(list(value = 1)))
  expect_error(finalize_v021_inference_artifact(make_v021_artifact("unfinished")))
  for (state in c("failed", "incomplete", "unavailable", "unfinished", "pre_inference"))
    expect_error(finalize_v021_inference_artifact(make_v021_artifact(execution_state = state)),
                 info = state)
  incomplete <- make_v021_artifact()
  incomplete$inference_results$posterior_intervals <- NULL
  expect_error(finalize_v021_inference_artifact(incomplete), "lacks required")
  hidden <- make_v021_artifact()
  attr(hidden$inference_results, "truth") <- matrix(0, 2, 2)
  expect_error(finalize_v021_inference_artifact(hidden), "attributes")
  finalized <- finalize_v021_inference_artifact(make_v021_artifact())
  expect_identical(finalized$.v021_finalized, TRUE)
  expect_false(".v021_finalized" %in% names(make_v021_artifact()))
  expect_error(finalize_v021_inference_artifact(finalized), "already finalized")
})

test_that("truth joins are pure and preserve every frozen inference field", {
  original <- make_v021_artifact()
  finalized <- finalize_v021_inference_artifact(original)
  frozen <- .v021_copy(finalized$inference_results)
  joined <- join_v021_truth_post_inference(
    finalized,
    list(absolute_A = 0, absolute_truth_outcome = "absolute_zero",
         execution_state = "failed", original_reporting_state = "not-reportable")
  )
  expect_false("evaluation" %in% names(finalized))
  expect_false(".v021_finalized" %in% names(original))
  expect_identical(finalized$inference_results, frozen)
  expect_identical(joined$inference_results, frozen)
  for (nm in v021_inference_result_fields)
    expect_identical(joined$inference_results[[nm]], frozen[[nm]], info = nm)
  expect_identical(joined$evaluation$absolute_truth_outcome, "absolute_zero")
  expect_identical(joined$evaluation$execution_state, "failed")
  expect_identical(joined$evaluation$original_reporting_state, "not-reportable")
})

test_that("truth joins reject non-final and invalid inference states", {
  expect_error(join_v021_truth_post_inference(make_v021_artifact(), list(absolute_A = 1)),
               "finalized")
  for (state in c("failed", "incomplete", "unavailable", "unfinished", "pre_inference")) {
    x <- make_v021_artifact(execution_state = state)
    x$.v021_finalized <- TRUE
    expect_error(join_v021_truth_post_inference(x, list(absolute_A = 1)), info = state)
  }
})

test_that("truth-derived fields cannot become diagnostic features", {
  artifact <- make_v021_artifact()
  artifact$inference_results$feature_inputs$truth_sign <- -1
  expect_error(finalize_v021_inference_artifact(artifact), "truth_sign")
  artifact <- make_v021_artifact()
  artifact$inference_results$feature_inputs$benchmark_outcome <- "opposite"
  expect_error(finalize_v021_inference_artifact(artifact), "benchmark_outcome")
})

test_that("pair manifests preserve both directions and caller assignments", {
  tasks <- make_v021_tasks()[1:4, ]
  ids <- c("m1::d1::a|b", "m1::d1::a|c")
  split <- make_split_spec("unordered_pair", ids,
                           c("development", "locked_evaluation"))
  manifest <- build_pair_group_manifest(
    tasks, c(dataset = "dataset_id", interaction_matrix = "matrix_id"), split, 17L
  )
  expect_silent(assert_no_pair_split_leakage(manifest))
  expect_equal(unname(vapply(split(manifest$partition, manifest$pair_group),
                             function(x) length(unique(x)), integer(1))), c(1L, 1L))
  expect_setequal(names(manifest), c(
    names(tasks), "pair_group", "dataset_group", "matrix_group", "grouping_level",
    "assignment_group", "partition", "manifest_seed"
  ))
  expect_true(all(manifest$manifest_seed == 17L))
  expect_identical(manifest, build_pair_group_manifest(
    tasks, c(dataset = "dataset_id", interaction_matrix = "matrix_id"), split, 17L
  ))
  expect_identical(unique(manifest$partition[manifest$pair_group == ids[[1L]]]), "development")

  keyed <- tasks
  keyed$pair_key <- rep(c("existing-pair-1", "existing-pair-2"), each = 2)
  keyed_manifest <- build_pair_group_manifest(
    keyed, c(dataset = "dataset_id", interaction_matrix = "matrix_id"),
    make_split_spec("unordered_pair", c("existing-pair-1", "existing-pair-2"),
                    c("development", "locked_evaluation")), 17L
  )
  expect_setequal(unique(keyed_manifest$pair_group), c("existing-pair-1", "existing-pair-2"))
})

test_that("dataset and interaction-matrix grouping are preserved", {
  tasks <- make_v021_tasks()
  keys <- c(dataset = "dataset_id", interaction_matrix = "matrix_id")
  dataset_manifest <- build_pair_group_manifest(
    tasks, keys,
    make_split_spec("dataset", c("d1", "d2"), c("development", "locked_evaluation")), 4L
  )
  matrix_manifest <- build_pair_group_manifest(
    tasks, keys,
    make_split_spec("interaction_matrix", c("m1", "m2"),
                    c("threshold_calibration", "locked_evaluation")), 4L
  )
  expect_equal(unname(vapply(split(dataset_manifest$partition, dataset_manifest$dataset_group),
                             function(x) length(unique(x)), integer(1))), c(1L, 1L))
  expect_equal(unname(vapply(split(matrix_manifest$partition, matrix_manifest$matrix_group),
                             function(x) length(unique(x)), integer(1))), c(1L, 1L))
  expect_silent(assert_no_pair_split_leakage(dataset_manifest))
  expect_silent(assert_no_pair_split_leakage(matrix_manifest))
})

test_that("manifests reject leakage, duplicates, conflicts, and incomplete assignments", {
  tasks <- make_v021_tasks()[1:4, ]
  keys <- c(dataset = "dataset_id", interaction_matrix = "matrix_id")
  ids <- c("m1::d1::a|b", "m1::d1::a|c")
  split <- make_split_spec("unordered_pair", ids, c("development", "locked_evaluation"))
  manifest <- build_pair_group_manifest(tasks, keys, split, 1L)
  leaking <- manifest
  leaking$partition[[1L]] <- "locked_evaluation"
  expect_error(assert_no_pair_split_leakage(leaking), "crosses")
  duplicate_tasks <- rbind(tasks, tasks[1, ])
  expect_error(build_pair_group_manifest(duplicate_tasks, keys, split, 1L), "Duplicate")
  duplicate_assignment <- split
  duplicate_assignment$assignments <- rbind(split$assignments, split$assignments[1, ])
  expect_error(build_pair_group_manifest(tasks, keys, duplicate_assignment, 1L), "Duplicate")
  conflicting <- split
  conflicting$assignments <- rbind(
    split$assignments,
    data.frame(group_id = ids[[1L]], partition = "locked_evaluation")
  )
  expect_error(build_pair_group_manifest(tasks, keys, conflicting, 1L), "Duplicate")
  incomplete <- split
  incomplete$assignments <- incomplete$assignments[1, , drop = FALSE]
  expect_error(build_pair_group_manifest(tasks, keys, incomplete, 1L), "incomplete")
  expect_error(build_pair_group_manifest(tasks[-1, ], keys, split, 1L), "incomplete")
  truth_annotated <- tasks
  truth_annotated$truth_sign <- -1
  expect_error(build_pair_group_manifest(truth_annotated, keys, split, 1L), "outside")
})

test_that("V021-01 accepts assignments without selecting or freezing V021-08 splits", {
  tasks <- make_v021_tasks()[1:4, ]
  keys <- c(dataset = "dataset_id", interaction_matrix = "matrix_id")
  ids <- c("m1::d1::a|b", "m1::d1::a|c")
  first <- build_pair_group_manifest(
    tasks, keys, make_split_spec("unordered_pair", ids, c("development", "locked_evaluation")), 8L
  )
  second <- build_pair_group_manifest(
    tasks, keys, make_split_spec("unordered_pair", ids,
                                 c("threshold_calibration", "development")), 8L
  )
  expect_false(identical(first$partition, second$partition))
  expect_error(build_pair_group_manifest(tasks, keys, split_spec = NULL, seed = 8L))
})

test_that("absolute truth classification has exactly three substantive outcomes", {
  classified <- list(
    classify_absolute_truth_outcome(1, 2),
    classify_absolute_truth_outcome(-1, -2),
    classify_absolute_truth_outcome(1, -2),
    classify_absolute_truth_outcome(-1, 2),
    classify_absolute_truth_outcome(1, 0),
    classify_absolute_truth_outcome(-1, 0)
  )
  outcomes <- vapply(classified, `[[`, character(1), "absolute_truth_outcome")
  expect_setequal(unique(outcomes),
                  c("same_nonzero_sign", "opposite_nonzero_sign", "absolute_zero"))
  expect_false(any(grepl("failed|incomplete|unavailable|indeterminate|insignificant|reportable|omitted",
                         outcomes)))
  expect_identical(classify_absolute_truth_outcome(1, 0)$absolute_truth_outcome,
                   "absolute_zero")
  expect_identical(classify_absolute_truth_outcome(-1, 0)$absolute_truth_outcome,
                   "absolute_zero")
})

test_that("missing or invalid posterior signs are unavailable with structured reasons", {
  cases <- list(
    missing = NA_real_, unavailable = Inf, indeterminate = 0,
    invalid = 2, malformed = "positive", absent = numeric()
  )
  for (nm in names(cases)) {
    outcome <- classify_absolute_truth_outcome(cases[[nm]], 0)
    expect_true(is.na(outcome$absolute_truth_outcome), info = nm)
    expect_identical(outcome$absolute_truth_outcome_status, "unavailable", info = nm)
    expect_type(outcome$absolute_truth_outcome_reason, "list")
    expect_true(nzchar(outcome$absolute_truth_outcome_reason$code), info = nm)
  }
  missing_truth <- classify_absolute_truth_outcome(1, NA_real_)
  expect_true(is.na(missing_truth$absolute_truth_outcome))
  expect_identical(missing_truth$absolute_truth_outcome_reason$code, "absolute_A_unavailable")
})

test_that("execution and reporting states remain orthogonal to scientific zero", {
  artifact <- finalize_v021_inference_artifact(make_v021_artifact())
  states <- c("missing", "indeterminate", "failed", "unavailable", "not-reportable")
  for (state in states) {
    joined <- join_v021_truth_post_inference(
      artifact,
      list(absolute_truth_outcome = NA_character_,
           absolute_truth_outcome_status = "unavailable",
           execution_state = state, original_reporting_state = state)
    )
    expect_true(is.na(joined$evaluation$absolute_truth_outcome), info = state)
    expect_identical(joined$evaluation$execution_state, state, info = state)
    expect_identical(joined$evaluation$original_reporting_state, state, info = state)
    expect_false(identical(joined$evaluation$absolute_truth_outcome, "absolute_zero"), info = state)
  }
})

test_that("truth joining cannot flip coefficients or posterior signs", {
  artifact <- finalize_v021_inference_artifact(make_v021_artifact())
  coefficient <- artifact$inference_results$posterior_coefficients
  signs <- artifact$inference_results$posterior_sign_probabilities
  joined <- join_v021_truth_post_inference(
    artifact, list(absolute_A = 10, requested_sign = 1,
                   absolute_truth_outcome = "opposite_nonzero_sign")
  )
  expect_identical(joined$inference_results$posterior_coefficients, coefficient)
  expect_identical(joined$inference_results$posterior_sign_probabilities, signs)
})

make_v021_analysis_specification <- function(truth_payload = NULL) {
  tasks <- make_v021_tasks()[1:4, c(
    "task_id", "direction_index", "target", "source", "seed")]
  # truth_payload is deliberately external to every constructor argument.
  force(truth_payload)
  build_v021_truth_free_analysis_specification(
    dataset_id = "d1", observed_data_id = "observed-dataset-1",
    taxa_order = c("a", "b", "c"), task_manifest = tasks,
    subject_universe = c("subject-b", "subject-a"),
    time_metadata = data.frame(
      subject = c("subject-a", "subject-a", "subject-b", "subject-b"),
      time = c(0, 1, 0, 1), stringsAsFactors = FALSE),
    preprocessing_config = list(schema = "closure-v1", min_pairs = 4L),
    posterior_config = list(chains = 4L, retry_limit = 0L),
    kfold_config = list(K = 2L, R = 1L), public_seed = 101L,
    kfold_seed = 808L, scorer_name = "student-t-scale-mixture-kalman-ou-q16",
    provenance = list(git_commit = "fixture", config_id = "fixture-v1")
  )
}

make_v021_final_record <- function(spec, row = 1L, status = "completed") {
  build_v021_finalized_inference_record(
    spec, spec$task_manifest[row, , drop = FALSE], status,
    posterior_summaries = list(mean = -0.2, interval = c(-0.4, -0.1)),
    psp_lfsr = list(posterior_sign = -1L, PSP = 0.99, LFSR = 0.01),
    diagnostics = list(class = "converged", rhat = 1.001),
    predictive_eligibility = TRUE,
    pair_specific_elpd = list(value = -12, method = "student-t-scale-mixture-kalman-ou-q16"),
    failure_information = list(reason = NA_character_),
    seed_split_provenance = list(direction_seed = spec$direction_seeds[[row]],
                                 kfold_seed = spec$kfold_seed),
    finalized_timestamp = "2026-08-03 UTC"
  )
}

test_that("formal analysis specification is truth-free and truth perturbation invariant", {
  first <- make_v021_analysis_specification(matrix(c(0, 1, -1, 0), 2))
  second <- make_v021_analysis_specification(matrix(c(99, -8, 4, 2), 2))
  expect_identical(first, second)
  expect_silent(validate_v021_truth_free_analysis_specification(first))
  expect_identical(first$task_manifest, second$task_manifest)
  expect_identical(first$task_manifest$direction_index,
                   second$task_manifest$direction_index)
  expect_identical(first$direction_seeds, second$direction_seeds)
  expect_identical(first$kfold_seed, second$kfold_seed)
  expect_identical(first$preprocessing_config, second$preprocessing_config)
  expect_identical(first$posterior_config, second$posterior_config)
  expect_identical(first$kfold_config, second$kfold_config)
  expect_identical(first$subject_universe, c("subject-a", "subject-b"))
})

test_that("recursive prohibited fields identify exact offending paths", {
  allowed <- list(
    truth_free_note = "prose is not scanned by substring",
    filename = "oracle-report-is-not-a-field-name.tsv",
    metadata = data.frame(subject = "s1", coverage_note = "descriptive"),
    nested = list(predictive_eligibility = TRUE, PSP = .99, LFSR = .01)
  )
  expect_silent(assert_v021_truth_free_schema(allowed, "allowed"))

  attacks <- list(
    top = list(absolute_A = 1),
    nested = list(metadata = list(corrected_oracle = -.2)),
    list_column = data.frame(id = 1, payload = I(list(list(truth_sign = -1))))
  )
  expect_error(assert_v021_truth_free_schema(attacks$top, "spec"),
               "spec\\$absolute_A")
  expect_error(assert_v021_truth_free_schema(attacks$nested, "spec"),
               "spec\\$metadata\\$corrected_oracle")
  expect_error(assert_v021_truth_free_schema(attacks$list_column, "manifest"),
               "manifest\\$payload\\[\\[1\\]\\]\\$truth_sign")

  spec <- make_v021_analysis_specification()
  spec$posterior_config$withholding_truth_label <- "bad"
  expect_error(validate_v021_truth_free_analysis_specification(spec),
               "posterior_config\\$withholding_truth_label")
})

test_that("finalized inference records enforce terminal identity and exclude truth", {
  spec <- make_v021_analysis_specification()
  record <- make_v021_final_record(spec)
  expect_silent(validate_v021_finalized_inference_record(record))
  expect_identical(record$record_schema, v021_finalized_inference_schema)
  expect_true(record$finalized)

  unfinished <- record
  unfinished$finalized <- FALSE
  expect_error(validate_v021_finalized_inference_record(unfinished), "finalized state")
  running <- record
  running$terminal_status <- "running"
  expect_error(validate_v021_finalized_inference_record(running), "terminal status")
  missing <- record
  missing$posterior_summaries <- NULL
  expect_error(validate_v021_finalized_inference_record(missing), "lacks required")
  leaked <- record
  leaked$diagnostics$empirical_false_sign <- TRUE
  expect_error(validate_v021_finalized_inference_record(leaked),
               "diagnostics\\$empirical_false_sign")
  wrong_task <- spec$task_manifest[1, , drop = FALSE]
  wrong_task$source <- "c"
  expect_error(build_v021_finalized_inference_record(spec, wrong_task, "completed"),
               "disagrees")
})

test_that("post-inference truth table is separate, immutable, and uses A target source", {
  spec <- make_v021_analysis_specification()
  records <- list(make_v021_final_record(spec, 1L), make_v021_final_record(spec, 2L))
  records[[2L]]$psp_lfsr$posterior_sign <- 1L
  frozen <- serialize(records, NULL, version = 3)
  A <- matrix(c(0, 2, -4, -3, 0, 5, 7, -8, 0), 3, 3, byrow = TRUE,
              dimnames = list(spec$taxa_order, spec$taxa_order))
  oracle <- data.frame(
    dataset_id = "d1", task_id = c(1L, 1L), direction_index = c(1L, 2L),
    target = c("a", "b"), source = c("b", "a"),
    corrected_oracle = c(-.5, .4), stringsAsFactors = FALSE)
  truth <- build_v021_post_inference_truth_table(records, A, oracle)
  expect_identical(serialize(records, NULL, version = 3), frozen)
  expect_identical(truth$absolute_A, c(A["a", "b"], A["b", "a"]))
  expect_identical(truth$source, c("b", "a"))
  expect_identical(truth$target, c("a", "b"))
  expect_false(identical(truth$absolute_A[[1L]], truth$absolute_A[[2L]]))
  expect_identical(truth$posterior_vs_oracle, c("agree", "agree"))
  expect_false(any(c("absolute_A", "structural_zero", "corrected_oracle",
                     "glv_to_oracle_class", "posterior_vs_oracle") %in%
                   names(records[[1L]])))

  duplicated <- c(records, records[1L])
  expect_error(build_v021_post_inference_truth_table(duplicated, A, oracle),
               "duplicated or ambiguous")
  wrong_taxa <- A[c("b", "a", "c"), , drop = FALSE]
  expect_error(build_v021_post_inference_truth_table(records, wrong_taxa, oracle),
               "taxa order")
  incomplete <- records
  incomplete[[1L]]$terminal_status <- "incomplete"
  expect_error(build_v021_post_inference_truth_table(incomplete, A, oracle),
               "terminal status")
})

test_that("post-inference projection classes keep structural zero distinct", {
  expect_identical(classify_v021_glv_to_oracle(1, 2), "same_sign")
  expect_identical(classify_v021_glv_to_oracle(1, -2), "sign_reversal")
  expect_identical(classify_v021_glv_to_oracle(0, -2),
                   "absolute_zero_to_projected_nonzero")
  expect_identical(classify_v021_glv_to_oracle(2, 0),
                   "absolute_nonzero_to_projected_zero")
  expect_identical(classify_v021_glv_to_oracle(0, 0), "both_zero")
})

test_that("truth-free scientific decisions depend only on their declared inputs", {
  external_truth <- list(A = matrix(c(0, 1, -1, 0), 2), oracle_sign = -1L)
  draws <- c(-.4, -.2, -.1, .1)
  diagnostics <- list(worst_rhat = 1.001, min_ess_bulk = 800,
                      min_ess_tail = 700, ebfmi_min = .9,
                      n_divergent = 0L, n_treedepth_hit = 0L)
  split_one <- pclvbayes:::.make_repkfold_splits(
    factor(c("s3", "s1", "s2", "s4"), levels = c("s4", "s3", "s2", "s1")),
    K = 2L, R = 2L, seed = 44L)
  lfsr_one <- pclvbayes:::.lfsr_from_two_sided(2 * min(mean(draws > 0), mean(draws < 0)))
  external_truth$A[,] <- 999
  external_truth$oracle_sign <- 1L
  split_two <- pclvbayes:::.make_repkfold_splits(
    c("s1", "s2", "s3", "s4"), K = 2L, R = 2L, seed = 44L)
  lfsr_two <- pclvbayes:::.lfsr_from_two_sided(2 * min(mean(draws > 0), mean(draws < 0)))
  expect_identical(split_one, split_two)
  expect_identical(lfsr_one, lfsr_two)
  expect_identical(diagnostics, diagnostics)
})

