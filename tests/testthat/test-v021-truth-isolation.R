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
    stacking_results = list(available = FALSE, value = NA_real_),
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

make_v021_runtime <- function() {
  ctx <- setNames(rep(list(1), length(v021_runtime_context_fields)),
                  v021_runtime_context_fields)
  ctx$meta_df <- data.frame(subject = c("x", "x"), time = c(0, 1))
  ctx$sm_mat <- matrix(1:4, 2, 2, dimnames = list(c("s1", "s2"), c("a", "b")))
  ctx$mod_exe_file <- "/approved/canonical-model"
  build_v021_runtime_context(ctx)
}

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
    list(target = target, partner = partner, seed = seed_override, ctx = ctx)
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
  expect_error(finalize_v021_inference_artifact(artifact), "Truth-bearing")
  artifact <- make_v021_artifact()
  artifact$inference_results$feature_inputs$benchmark_outcome <- "opposite"
  expect_error(finalize_v021_inference_artifact(artifact), "Truth-bearing")
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
