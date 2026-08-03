source(testthat::test_path("../../benchmarks/mtist/v021_truth_isolation.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_resource_policy.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_checkpoint_manifest.R"))
source(testthat::test_path("../../benchmarks/mtist/ten_species_helpers.R"))
source(testthat::test_path("../../benchmarks/mtist/v021_full_100_species.R"))

.resource_manifest <- function(root, taxa = c("a", "b", "c")) {
  tasks <- ten_species_direction_index(taxa, 101L)
  build_v021_execution_manifest(
    "fixture", taxa, tasks[c("task_id", "direction_index", "target", "source", "seed")],
    101L, 202L, list(id = "prep"),
    list(id = "posterior", chains = 4L, parallel_chains = 4L,
         threads_per_chain = 1L),
    list(K = 2L, R = 1L, parallel_chains = 1L), list(id = "q16"),
    list(code_commit = "fixture"))
}

.fake_identity <- function(pid = Sys.getpid(), start = "123", host = "host-a",
                           boot = "boot-a") {
  list(pid = as.integer(pid), start_time = start, hostname = host, boot_id = boot)
}

test_that("canonical policy derives the twelve-chain capacity contract", {
  policy <- build_v021_resource_policy()
  expect_identical(policy$policy_schema, "v021_resource_policy_v3")
  expect_identical(policy$main_chains, 4L)
  expect_identical(policy$kfold_chains, 4L)
  expect_identical(policy$kfold_parallel_chains, 1L)
  expect_identical(policy$logical_host_threads, 16L)
  expect_identical(policy$reserved_host_threads, 4L)
  expect_identical(policy$maximum_active_cmdstan_chains, 12L)
  expect_identical(policy$proposed_outer_concurrency, 3L)
  expect_identical(policy$maximum_concurrent_kfold_fits, 1L)
  expect_identical(policy$controller_worker_limit, 3L)
  expect_true(policy$retries_share_chain_budget)
  expect_true(all(policy$environment_thread_caps == "1"))
  expect_identical(v021_job_demand(policy, "main_fit"), 4L)
  expect_identical(v021_job_demand(policy, "retry_fit"), 4L)
  expect_identical(v021_job_demand(policy, "kfold_fit"), 1L)
  expect_error(build_v021_resource_policy(main_chains = 13L,
                                           retry_chains = 13L,
                                           kfold_chains = 13L,
                                           confirmation_chains = 13L),
               "positive integer|fit within|unsafe|exceeds")
  expect_error(build_v021_resource_policy(main_chains = 4.5), "integer")
  expect_silent(assert_v021_truth_free_schema(policy, "resource policy"))
  expect_identical(v021_resource_policy_hash(policy), v021_resource_policy_hash(policy))
})

test_that("reservations share capacity across posterior retry and K-fold work", {
  policy <- build_v021_resource_policy()
  reservations <- new_v021_reservations()
  reservations <- reserve_v021_capacity(reservations, policy, "d1", 1L, "main_fit")
  reservations <- reserve_v021_capacity(reservations, policy, "d2", 1L, "retry_fit")
  expect_identical(sum(reservations$reserved_chain_slots), 8L)
  reservations <- reserve_v021_capacity(reservations, policy, "d3", 1L, "main_fit")
  expect_identical(sum(reservations$reserved_chain_slots), 12L)
  expect_error(reserve_v021_capacity(reservations, policy, "d4", 1L, "kfold_fit"),
               "Insufficient")
  expect_error(release_v021_capacity(reservations, policy, "d1", 1L, "main_fit", FALSE),
               "worker termination")
  reservations <- release_v021_capacity(reservations, policy, "d1", 1L,
                                         "main_fit", TRUE)
  expect_identical(sum(reservations$reserved_chain_slots[reservations$state == "reserved"]), 8L)
})

test_that("five four-chain jobs schedule deterministically as 12 8", {
  policy <- build_v021_resource_policy()
  ids <- sprintf("direction-%06d", c(5L, 1L, 4L, 2L, 3L))
  waves <- plan_v021_scheduler_waves(ids, policy)
  expect_identical(vapply(waves, `[[`, integer(1), "reserved_chain_slots"), c(12L, 8L))
  expect_identical(unlist(lapply(waves, `[[`, "task_identities")), sort(ids, method = "radix"))
  expect_identical(waves, plan_v021_scheduler_waves(rev(ids), policy))
  expect_identical(waves, plan_v021_scheduler_waves(factor(rev(ids), levels = ids), policy))
  expect_true(all(vapply(waves, `[[`, integer(1), "reserved_chain_slots") <= 12L))
})

test_that("thread limits are explicit scoped and reject incompatible values", {
  expected <- c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "STAN_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
                "RCPP_PARALLEL_NUM_THREADS")
  limits <- v021_single_thread_environment()
  expect_identical(names(limits), expected)
  expect_true(all(limits == "1"))
  bad <- limits; bad[["OMP_NUM_THREADS"]] <- "2"
  expect_error(validate_v021_single_thread_environment(bad), "exactly one")
})

test_that("controller ownership is exclusive and provenance-bound", {
  root <- tempfile("v021-owner-"); dir.create(root)
  identity <- .fake_identity()
  owner <- acquire_v021_controller_ownership(root, "manifest", "config", identity)
  expect_silent(validate_v021_controller_ownership(owner, root, "manifest", "config"))
  expect_error(acquire_v021_controller_ownership(root, "manifest", "config", identity),
               "already owned")
  expect_error(validate_v021_controller_ownership(owner, root, "other", "config"),
               "mismatch")
  expect_silent(release_v021_controller_ownership(root, owner))
  expect_false(dir.exists(v021_ownership_path(root)))
  expect_error(validate_v021_controller_ownership(owner, root, "manifest", "config"),
               "does not hold")
})

test_that("foreign and stale owners require explicit reconciliation", {
  root <- tempfile("v021-stale-"); dir.create(root)
  owner <- acquire_v021_controller_ownership(
    root, "manifest", "config", .fake_identity(pid = 2147480000L))
  expect_error(reconcile_v021_stale_ownership(
    root, .fake_identity(pid = 1L, host = "host-b")), "Foreign-host")
  expect_error(reconcile_v021_stale_ownership(root, .fake_identity(pid = 1L)),
               "administrative_override")
  expect_silent(reconcile_v021_stale_ownership(
    root, .fake_identity(pid = 1L), administrative_override = TRUE))
  expect_false(dir.exists(v021_ownership_path(root)))
  expect_true(is.list(owner))
})

test_that("same-host live ownership cannot be reconciled", {
  root <- tempfile("v021-live-"); dir.create(root)
  identity <- .v021_local_process_identity()
  owner <- acquire_v021_controller_ownership(root, "manifest", "config", identity)
  expect_error(reconcile_v021_stale_ownership(root, identity,
                                               administrative_override = TRUE),
               "still live")
  release_v021_controller_ownership(root, owner)
})

test_that("task claim is serialized by ownership and capacity", {
  root <- tempfile("v021-claim-"); dir.create(root)
  manifest <- .resource_manifest(root)
  status <- new_v021_task_status(manifest)
  ledger <- new_v021_attempt_ledger(manifest)
  policy <- build_v021_resource_policy()
  owner <- acquire_v021_controller_ownership(
    root, manifest$manifest_hash, manifest$configuration_hash, .fake_identity())
  claimed <- claim_v021_task_with_capacity(
    status, ledger, manifest, 1L, owner, root, policy, new_v021_reservations(),
    "2026-08-03T00:00:00Z")
  expect_identical(claimed$status$tasks$state[[1L]], "running")
  expect_identical(claimed$reservations$reserved_chain_slots, 4L)
  bad <- owner; bad$manifest_hash <- "wrong"
  expect_error(claim_v021_task_with_capacity(
    status, ledger, manifest, 1L, bad, root, policy, new_v021_reservations(),
    "2026-08-03T00:00:00Z"), "mismatch")
  release_v021_controller_ownership(root, owner)
})

test_that("resource audit distinguishes process exit from scientific completion", {
  policy <- build_v021_resource_policy()
  r <- reserve_v021_capacity(new_v021_reservations(), policy, "d1", 1L, "main_fit")
  audit <- audit_v021_resource_state(policy, r, TRUE, "worker-1", 4L)
  expect_identical(audit$total_reserved_chain_slots, 4L)
  expect_true(audit$ownership_valid)
  expect_error(audit_v021_resource_state(policy, r, TRUE, "worker-1", 5L),
               "exceed registered")
  r <- release_v021_capacity(r, policy, "d1", 1L, "main_fit", TRUE)
  expect_identical(sum(r$reserved_chain_slots[r$state == "reserved"]), 0L)
  # V021-03 completion remains artifact-gated; release does not change task status.
  expect_identical(new_v021_task_status(.resource_manifest(tempdir()))$tasks$state[[1L]],
                   "pending")
})

test_that("reservation persistence is atomic and malformed state fails closed", {
  policy <- build_v021_resource_policy()
  reservations <- reserve_v021_capacity(
    new_v021_reservations(), policy, "d1", 1L, "main_fit")
  path <- tempfile("v021-reservations-", fileext = ".rds")
  expect_silent(write_v021_reservations_atomic(reservations, path, policy))
  expect_identical(read_v021_reservations(path, policy), reservations)
  expect_false(any(grepl("\\.tmp-", list.files(dirname(path), all.files = TRUE))))
  saveRDS(list(partial = TRUE), path)
  expect_error(read_v021_reservations(path, policy), "schema")
})

test_that("resource planning is truth invariant and excludes scheduler metadata from manifest", {
  policy <- build_v021_resource_policy()
  expect_error(assert_v021_truth_free_schema(
    c(policy, list(truth_sign = 1L)), "resource policy"), "truth_sign")
  manifest <- .resource_manifest(tempdir())
  hash <- manifest$manifest_hash
  waves <- plan_v021_scheduler_waves(manifest$tasks$directed_task_id[1:5], policy)
  expect_identical(hash, manifest$manifest_hash)
  expect_false("worker_count" %in% names(manifest$configuration))
  expect_identical(vapply(waves, `[[`, integer(1), "reserved_chain_slots"), c(12L, 8L))
})

test_that("dry-run exposes policy without sampling", {
  config <- build_v021_full_config(tempfile("v021-dry-"))
  prepared <- prepare_v021_full_execution(
    config, paste0("species_", 0:99),
    list(code_commit = "fixture", benchmark_schema = config$benchmark_schema),
    initialize = FALSE)
  dry <- prepare_v021_resource_dry_run(prepared, build_v021_resource_policy())
  expect_false(dry$sampling_launched)
  expect_identical(dry$maximum_concurrent_tasks, 3L)
  expect_identical(dry$maximum_active_cmdstan_chains, 12L)
  expect_identical(dry$waves[[1L]]$reserved_chain_slots, 12L)
})
