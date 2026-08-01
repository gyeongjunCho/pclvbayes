profile_helpers_path <- testthat::test_path("..", "..", "benchmarks", "profiling", "profile_helpers.R")
source(profile_helpers_path, local = TRUE)

test_that("profiling timing records have deterministic schema", {
  row <- new_timing_row("x", "phase", 1, .2, .1)
  expect_identical(names(row), profiling_timing_schema())
  expect_identical(row$phase, "phase")
})

test_that("disabled profiling is observationally inert", {
  value <- profile_disabled({ list(a = 1, b = letters[1:2]) }, enabled = FALSE)
  expect_identical(value, list(a = 1, b = letters[1:2]))
})

test_that("unavailable platform metrics are represented by NA", {
  expect_true(is.na(read_proc_status_kb("DefinitelyNotAProcStatusField")))
  expect_true(is.na(system_time_available()) || nzchar(system_time_available()))
})

test_that("profiling errors retain their original failure semantics", {
  failure <- pclvbayes:::.pclv_failure("sampling", "cmdstan_execution_failed")
  observed <- profile_disabled(failure, enabled = FALSE)
  expect_identical(observed, failure)
  expect_s3_class(observed, "pclv_failure")
})
