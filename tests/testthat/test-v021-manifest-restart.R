source(testthat::test_path("../../benchmarks/mtist/mtist_adapter.R"))
source(testthat::test_path("../../benchmarks/mtist/ten_species_helpers.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_truth_isolation.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_checkpoint_manifest.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_full_100_species.R"))

v021_small_tasks <- function() {
  data.frame(
    task_id = rep(1:3, each = 2), direction_index = 1:6,
    target = c("a", "b", "a", "c", "b", "c"),
    source = c("b", "a", "c", "a", "c", "b"),
    seed = 101:106, stringsAsFactors = FALSE)
}

v021_small_manifest <- function(tasks = v021_small_tasks(),
                                 posterior = list(chains = 4L),
                                 provenance = list(code_commit = "fixture")) {
  build_v021_execution_manifest(
    dataset_id = "dataset-1", taxa_order = c("a", "b", "c"),
    task_table = tasks, public_seed = 11L, kfold_seed = 22L,
    preprocessing_config = list(schema = "closure-v1", alpha = .5),
    posterior_config = posterior,
    kfold_config = list(K = 2L, R = 1L),
    predictive_config = list(scorer = "student-t-scale-mixture-kalman-ou-q16"),
    provenance = provenance)
}

with_v021_execution_root <- function(code) {
  root <- tempfile("v021-execution-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  force(code)(root)
}

v021_fixture_artifact <- function(manifest, root, task = 1L) {
  build_v021_completion_artifact(
    manifest, task, root,
    posterior_summary = list(mean = -.2, interval = c(-.4, -.1)),
    diagnostics = list(class = "converged", rhat = 1.001),
    psp_lfsr = list(PSP = .99, LFSR = .01),
    predictive_eligibility = TRUE,
    subject_elpd = list(state = "not_executed"),
    seed_split_provenance = list(
      direction_seed = manifest$tasks$direction_seed[[task]],
      kfold_seed = manifest$tasks$kfold_seed[[task]]),
    failure_information = list(reason = NA_character_))
}

test_that("small canonical manifest has manually inspectable identities", {
  manifest <- v021_small_manifest()
  expect_silent(validate_v021_execution_manifest(manifest))
  expect_identical(manifest$unordered_pair_count, 3L)
  expect_identical(manifest$directed_task_count, 6L)
  expect_identical(manifest$tasks$pair_id, rep(sprintf("pair-%06d", 1:3), each = 2))
  expect_identical(manifest$tasks$directed_task_id, sprintf("direction-%06d", 1:6))
  expect_identical(manifest$tasks$source, c("b", "a", "c", "a", "c", "b"))
  expect_identical(manifest$tasks$target, c("a", "b", "a", "c", "b", "c"))
  expect_identical(manifest$tasks$source_index, c(2L, 1L, 3L, 1L, 3L, 2L))
  expect_identical(manifest$tasks$target_index, c(1L, 2L, 1L, 3L, 2L, 3L))
  expect_false(any(manifest$tasks$source == manifest$tasks$target))
  expect_identical(anyDuplicated(manifest$tasks$directed_task_id), 0L)
})

test_that("manifest ordering and hash ignore task row and factor-level order", {
  tasks <- v021_small_tasks()
  shuffled <- tasks[c(6, 2, 4, 1, 5, 3), ]
  shuffled$source <- factor(shuffled$source, levels = c("c", "b", "a"))
  shuffled$target <- factor(shuffled$target, levels = c("b", "a", "c"))
  first <- v021_small_manifest(tasks)
  second <- v021_small_manifest(shuffled)
  expect_identical(first, second)
  expect_identical(first$manifest_hash, second$manifest_hash)
})

test_that("full manifest contains exactly 4950 pairs and 9900 directions", {
  config <- build_v021_full_config(tempfile("v021-full-root-"))
  manifest <- build_v021_full_execution_manifest(
    config, paste0("species_", 0:99),
    provenance = list(code_commit = "fixture", benchmark_schema = config$benchmark_schema))
  expect_identical(nrow(manifest$tasks), 9900L)
  expect_identical(length(unique(manifest$tasks$pair_id)), 4950L)
  expect_identical(sum(manifest$tasks$source == manifest$tasks$target), 0L)
  expect_identical(anyDuplicated(manifest$tasks$directed_task_id), 0L)
  expect_identical(manifest$tasks$task_ordinal, 1:9900)
})

test_that("truth perturbation cannot alter manifest fields or hashes", {
  truth_one <- matrix(1, 3, 3)
  truth_two <- matrix(-99, 3, 3)
  build_without_truth <- function(truth) {
    force(truth)
    v021_small_manifest()
  }
  first <- build_without_truth(truth_one)
  second <- build_without_truth(truth_two)
  expect_identical(first, second)
  leaked <- v021_small_tasks()
  leaked$truth_sign <- -1L
  expect_error(v021_small_manifest(leaked), "Prohibited|schema")
  nested <- list(code_commit = "fixture", metadata = list(corrected_oracle = .2))
  expect_error(v021_small_manifest(provenance = nested), "corrected_oracle")
})

test_that("relevant configuration changes hashes but list ordering does not", {
  first <- v021_small_manifest(posterior = list(chains = 4L, warmup = 2000L))
  reordered <- v021_small_manifest(posterior = list(warmup = 2000L, chains = 4L))
  changed <- v021_small_manifest(posterior = list(chains = 4L, warmup = 1999L))
  expect_identical(first$configuration_hash, reordered$configuration_hash)
  expect_identical(first$manifest_hash, reordered$manifest_hash)
  expect_false(identical(first$configuration_hash, changed$configuration_hash))
  expect_false(identical(first$manifest_hash, changed$manifest_hash))
  expect_error(v021_small_manifest(
    provenance = list(code_commit = "fixture", start_timestamp = "now")),
    "Runtime")
})

test_that("immutable manifest rejects every material mismatch", {
  with_v021_execution_root(function(root) {
    first <- v021_small_manifest()
    path <- file.path(root, "canonical_manifest.rds")
    write_v021_execution_manifest_once(first, path)
    expect_identical(read_v021_execution_manifest(path), first)
    before <- unname(tools::md5sum(path))
    write_v021_execution_manifest_once(first, path)
    expect_identical(unname(tools::md5sum(path)), before)
    expect_error(write_v021_execution_manifest_once(
      v021_small_manifest(posterior = list(chains = 3L)), path),
      "configuration_hash|manifest_hash")
    changed <- first
    changed$dataset_id <- "different"
    changed$tasks$dataset_id <- "different"
    changed$manifest_hash <- v021_sha256(changed[setdiff(names(changed), "manifest_hash")])
    expect_error(compare_v021_execution_manifests(first, changed), "dataset_id")
    changed_taxa <- first
    changed_taxa$taxa_order <- rev(changed_taxa$taxa_order)
    changed_taxa$tasks$source_index <- match(changed_taxa$tasks$source,
                                             changed_taxa$taxa_order)
    changed_taxa$tasks$target_index <- match(changed_taxa$tasks$target,
                                             changed_taxa$taxa_order)
    changed_taxa$manifest_hash <- v021_sha256(
      changed_taxa[setdiff(names(changed_taxa), "manifest_hash")])
    expect_error(compare_v021_execution_manifests(first, changed_taxa), "taxa_order")
  })
})

test_that("legal task transitions are explicit and illegal rewrites fail", {
  expect_silent(validate_v021_state_transition("pending", "running"))
  expect_silent(validate_v021_state_transition("running", "completed"))
  expect_silent(validate_v021_state_transition("running", "failed"))
  expect_silent(validate_v021_state_transition("failed", "running"))
  expect_error(validate_v021_state_transition("completed", "pending"), "Illegal")
  expect_error(validate_v021_state_transition("pending", "completed"), "Illegal")
})

test_that("attempt ledger preserves identity, seeds, history, and artifact gate", {
  manifest <- v021_small_manifest()
  status <- new_v021_task_status(manifest)
  ledger <- new_v021_attempt_ledger(manifest)
  started <- start_v021_task_attempt(
    status, ledger, manifest, 1L, "2026-08-03 01:00 UTC", "worker-1")
  expect_identical(started$status$tasks$state[[1L]], "running")
  expect_identical(started$ledger$attempts$direction_seed[[1L]],
                   manifest$tasks$direction_seed[[1L]])
  expect_error(finish_v021_task_attempt(
    started$status, started$ledger, manifest, 1L, "completed",
    "2026-08-03 01:01 UTC"), "durable artifact")
  failed <- finish_v021_task_attempt(
    started$status, started$ledger, manifest, 1L, "failed",
    "2026-08-03 01:01 UTC", "interrupted", "interrupt")
  retry <- start_v021_task_attempt(
    failed$status, failed$ledger, manifest, 1L, "2026-08-03 01:02 UTC")
  expect_identical(retry$status$tasks$attempt_count[[1L]], 2L)
  expect_identical(retry$ledger$attempts$starting_state, c("pending", "failed"))
  expect_identical(unique(retry$ledger$attempts$direction_seed),
                   manifest$tasks$direction_seed[[1L]])
})

test_that("validated atomic writes reject malformed and partial state", {
  with_v021_execution_root(function(root) {
    manifest <- v021_small_manifest()
    initialize_v021_execution_root(manifest, root)
    paths <- v021_execution_paths(root)
    expect_silent(read_v021_task_status(paths$status, manifest))
    expect_silent(read_v021_attempt_ledger(paths$attempts, manifest))
    expect_length(list.files(root, pattern = "\\.tmp-", all.files = TRUE), 0L)
    stale <- file.path(root, ".task_status.rds.tmp-stale")
    writeBin(as.raw(1:3), stale)
    expect_silent(read_v021_task_status(paths$status, manifest))
    unlink(paths$status)
    expect_error(read_v021_task_status(paths$status, manifest), "missing")
    writeBin(as.raw(1:4), paths$status)
    expect_error(read_v021_task_status(paths$status, manifest), "malformed")
    saveRDS(list(status_schema = v021_task_status_schema), paths$status)
    expect_error(read_v021_task_status(paths$status, manifest), "Malformed")
  })
})

test_that("completion requires a valid immutable artifact with matching root", {
  with_v021_execution_root(function(root) {
    manifest <- v021_small_manifest()
    initialize_v021_execution_root(manifest, root)
    paths <- v021_execution_paths(root)
    artifact <- v021_fixture_artifact(manifest, root)
    path <- file.path(paths$artifact_root, "direction-000001.rds")
    write_v021_completion_artifact_atomic(artifact, path, manifest, root)
    expect_silent(validate_v021_completion_artifact(readRDS(path), manifest, root))
    expect_error(validate_v021_completion_artifact(artifact, manifest,
                                                    tempfile("other-root-")),
                 "result root")
    corrupt <- artifact
    corrupt$diagnostics <- list()
    expect_error(validate_v021_completion_artifact(corrupt, manifest, root),
                 "durable inference outputs|hash")
    missing <- artifact
    missing$posterior_summary <- NULL
    expect_error(validate_v021_completion_artifact(missing, manifest, root), "schema")
    fake_csv <- file.path(paths$artifact_root, "direction-000002.csv")
    writeLines("finished", fake_csv)
    expect_error(readRDS(fake_csv))
  })
})

test_that("resume skips completed, retries failed, and is idempotent", {
  with_v021_execution_root(function(root) {
    manifest <- v021_small_manifest()
    initialize_v021_execution_root(manifest, root)
    paths <- v021_execution_paths(root)
    status <- read_v021_task_status(paths$status, manifest)
    ledger <- read_v021_attempt_ledger(paths$attempts, manifest)

    artifact <- v021_fixture_artifact(manifest, root, 1L)
    artifact_path <- file.path(paths$artifact_root, "direction-000001.rds")
    write_v021_completion_artifact_atomic(artifact, artifact_path, manifest, root)
    one <- start_v021_task_attempt(status, ledger, manifest, 1L, "t1")
    one <- finish_v021_task_attempt(
      one$status, one$ledger, manifest, 1L, "completed", "t2",
      artifact_path = artifact_path, result_root = root)
    two <- start_v021_task_attempt(one$status, one$ledger, manifest, 2L, "t1")
    two <- finish_v021_task_attempt(
      two$status, two$ledger, manifest, 2L, "failed", "t2",
      terminal_reason = "fixture", failure_class = "fixture_error")
    write_v021_task_status_atomic(two$status, paths$status, manifest)
    write_v021_attempt_ledger_atomic(two$ledger, paths$attempts, manifest)

    first <- plan_v021_execution_resume(manifest, root, maximum_attempts = 2L)
    second <- plan_v021_execution_resume(manifest, root, maximum_attempts = 2L)
    expect_identical(first$runnable_indices, c(2:6))
    expect_identical(first$completed_indices, 1L)
    expect_identical(first$status, second$status)
    expect_identical(first$attempt_ledger, second$attempt_ledger)
    expect_identical(first$manifest$tasks$direction_seed,
                     manifest$tasks$direction_seed)
    expect_identical(first$manifest$tasks$kfold_seed,
                     manifest$tasks$kfold_seed)
    exhausted <- plan_v021_execution_resume(manifest, root, maximum_attempts = 1L)
    expect_false(2L %in% exhausted$runnable_indices)
  })
})

test_that("interrupted running is reconciled only after artifact validation", {
  with_v021_execution_root(function(root) {
    manifest <- v021_small_manifest()
    initialize_v021_execution_root(manifest, root)
    paths <- v021_execution_paths(root)
    status <- read_v021_task_status(paths$status, manifest)
    ledger <- read_v021_attempt_ledger(paths$attempts, manifest)
    running <- start_v021_task_attempt(status, ledger, manifest, 1L, "t1")
    write_v021_task_status_atomic(running$status, paths$status, manifest)
    write_v021_attempt_ledger_atomic(running$ledger, paths$attempts, manifest)
    plan <- plan_v021_execution_resume(
      manifest, root, maximum_attempts = 2L, persist_reconciliation = TRUE,
      timestamp = function() "reconciled")
    expect_identical(plan$status$tasks$state[[1L]], "failed")
    expect_identical(plan$status$tasks$terminal_reason[[1L]],
                     "interrupted_running_without_valid_artifact")
    expect_true(1L %in% plan$runnable_indices)
    expect_false(plan$status$tasks$artifact_validation[[1L]] == "validated")
    expect_identical(read_v021_task_status(paths$status, manifest), plan$status)
    expect_identical(plan$attempt_ledger$attempts$ending_state[[1L]], "failed")
    expect_identical(read_v021_attempt_ledger(paths$attempts, manifest),
                     plan$attempt_ledger)
  })
})

test_that("interrupted running with a validated artifact becomes completed", {
  with_v021_execution_root(function(root) {
    manifest <- v021_small_manifest()
    initialize_v021_execution_root(manifest, root)
    paths <- v021_execution_paths(root)
    status <- read_v021_task_status(paths$status, manifest)
    ledger <- read_v021_attempt_ledger(paths$attempts, manifest)
    running <- start_v021_task_attempt(status, ledger, manifest, 1L, "t1")
    artifact <- v021_fixture_artifact(manifest, root, 1L)
    artifact_path <- file.path(paths$artifact_root, "direction-000001.rds")
    write_v021_completion_artifact_atomic(artifact, artifact_path, manifest, root)
    running$status$tasks$artifact_path[[1L]] <- artifact_path
    write_v021_task_status_atomic(running$status, paths$status, manifest)
    write_v021_attempt_ledger_atomic(running$ledger, paths$attempts, manifest)
    plan <- plan_v021_execution_resume(manifest, root, maximum_attempts = 2L)
    expect_identical(plan$status$tasks$state[[1L]], "completed")
    expect_identical(plan$status$tasks$artifact_validation[[1L]],
                     "validated_after_interruption")
    expect_identical(plan$attempt_ledger$attempts$ending_state[[1L]], "completed")
    expect_false(1L %in% plan$runnable_indices)
  })
})

test_that("manifest status ledger and artifacts remain physically separate", {
  with_v021_execution_root(function(root) {
    manifest <- v021_small_manifest()
    paths <- initialize_v021_execution_root(manifest, root)
    expect_length(unique(c(paths$manifest, paths$status, paths$attempts,
                           paths$artifact_root)), 4L)
    expect_false("state" %in% names(manifest$tasks))
    expect_false(any(grepl("posterior|diagnostic|artifact", names(manifest$tasks))))
    expect_identical(nrow(read_v021_attempt_ledger(paths$attempts, manifest)$attempts), 0L)
  })
})

test_that("dry-run preparation creates no outputs and launches no sampler", {
  root <- tempfile("v021-dry-run-")
  config <- build_v021_full_config(root)
  prepared <- prepare_v021_full_execution(
    config, paste0("species_", 0:99),
    provenance = list(code_commit = "fixture", benchmark_schema = config$benchmark_schema),
    initialize = FALSE)
  expect_false(dir.exists(root))
  expect_false(prepared$sampling_launched)
  expect_true(prepared$dry_run)
  expect_identical(length(prepared$plan$runnable_indices), 9900L)
  expect_identical(prepared$status_summary$completed, 0L)
})
