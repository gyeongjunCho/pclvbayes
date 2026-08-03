test_that("control script accepts reviewed descendants and rejects unsafe repositories", {
  skip_if_not(nzchar(Sys.which("bash")))
  skip_if_not(nzchar(Sys.which("git")))
  script <- normalizePath(testthat::test_path(
    "../../benchmarks/mtist/v021_full_100_species_ctl.sh"))
  repo <- tempfile("v021-ctl-git-"); dir.create(repo)
  system2("git", c("-C", repo, "init", "-b", "v0.2.1-dev"), stdout = FALSE)
  system2("git", c("-C", repo, "config", "user.email", "fixture@example.test"))
  system2("git", c("-C", repo, "config", "user.name", "Fixture"))
  writeLines("one", file.path(repo, "fixture"))
  system2("git", c("-C", repo, "add", "fixture"))
  system2("git", c("-C", repo, "commit", "-m", "reviewed"), stdout = FALSE)
  reviewed <- system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE)
  writeLines("two", file.path(repo, "fixture"))
  system2("git", c("-C", repo, "commit", "-am", "descendant"), stdout = FALSE)
  descendant <- system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE)
  invoke <- function(command) suppressWarnings(system2("bash", c("-c", shQuote(sprintf(
    "source %s; %s", shQuote(script), command))), stdout = TRUE, stderr = TRUE))
  exit_status <- function(x) {
    status <- attr(x, "status")
    if (is.null(status)) 0L else as.integer(status)
  }
  expect_equal(exit_status(invoke(sprintf(
    "validate_reviewed_history %s %s %s %s", shQuote(repo), reviewed,
    reviewed, reviewed))), 0L)
  expect_equal(exit_status(invoke(sprintf(
    "validate_reviewed_history %s %s %s %s", shQuote(repo), descendant,
    reviewed, reviewed))), 0L)
  unrelated <- tempfile("v021-ctl-unrelated-"); dir.create(unrelated)
  system2("git", c("-C", unrelated, "init", "-b", "other"), stdout = FALSE)
  system2("git", c("-C", unrelated, "config", "user.email", "fixture@example.test"))
  system2("git", c("-C", unrelated, "config", "user.name", "Fixture"))
  writeLines("x", file.path(unrelated, "x")); system2("git", c("-C", unrelated, "add", "x"))
  system2("git", c("-C", unrelated, "commit", "-m", "unrelated"), stdout = FALSE)
  unrelated_head <- system2("git", c("-C", unrelated, "rev-parse", "HEAD"), stdout = TRUE)
  expect_true(exit_status(invoke(sprintf(
    "validate_reviewed_history %s %s %s %s", shQuote(unrelated), unrelated_head,
    reviewed, reviewed))) != 0L)
  expect_true(exit_status(invoke(sprintf(
    "validate_repository_state %s wrong-branch", shQuote(repo)))) != 0L)
  writeLines("dirty", file.path(repo, "fixture"))
  expect_true(exit_status(invoke(sprintf(
    "validate_repository_state %s v0.2.1-dev", shQuote(repo)))) != 0L)
})
