test_that("owned CmdStan output directories have an explicit lifecycle", {
  root <- withr::local_tempdir()
  withr::local_options(glvpair.output_root = root)

  fake_model <- list(sample = function(...) {
    args <- list(...)
    writeLines("sample", file.path(args$output_dir, "draws.csv"))
    new.env(parent = emptyenv())
  })

  fit <- pclvbayes:::.call_sample_silently(
    fake_model,
    list(chains = 1L, iter_warmup = 10L, iter_sampling = 10L),
    silent = TRUE
  )
  owned <- attr(fit, "pclv_owned_output_dir", exact = TRUE)
  expect_true(dir.exists(owned))
  expect_true(file.exists(file.path(owned, "draws.csv")))

  pclvbayes:::.cleanup_cmdstan_fit_output(fit)
  expect_false(dir.exists(owned))
})

test_that("failed sampling cleans package-owned output immediately", {
  root <- withr::local_tempdir()
  withr::local_options(glvpair.output_root = root)

  fake_model <- list(sample = function(...) {
    args <- list(...)
    writeLines("partial", file.path(args$output_dir, "partial.csv"))
    stop("synthetic sampler failure")
  })

  expect_error(
    pclvbayes:::.call_sample_silently(
      fake_model,
      list(chains = 1L, iter_warmup = 10L, iter_sampling = 10L),
      silent = FALSE
    ),
    "synthetic sampler failure"
  )
  expect_length(list.files(root, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("caller-owned CmdStan output directories are never deleted", {
  root <- withr::local_tempdir()
  supplied <- file.path(root, "caller-output")
  dir.create(supplied)

  fake_model <- list(sample = function(...) {
    args <- list(...)
    writeLines("sample", file.path(args$output_dir, "draws.csv"))
    new.env(parent = emptyenv())
  })

  fit <- pclvbayes:::.call_sample_silently(
    fake_model,
    list(output_dir = supplied, chains = 1L,
         iter_warmup = 10L, iter_sampling = 10L),
    silent = TRUE
  )
  expect_null(attr(fit, "pclv_owned_output_dir", exact = TRUE))
  pclvbayes:::.cleanup_cmdstan_fit_output(fit)
  expect_true(file.exists(file.path(supplied, "draws.csv")))
})

test_that("retry wrapper cleans a fit when diagnostics are interrupted", {
  root <- withr::local_tempdir()
  owned <- file.path(root, "owned-run")
  dir.create(owned)
  writeLines("sample", file.path(owned, "draws.csv"))

  testthat::local_mocked_bindings(
    .call_sample_silently = function(...) {
      fit <- new.env(parent = emptyenv())
      attr(fit, "pclv_owned_output_dir") <- owned
      fit
    },
    .summarise_diag = function(...) stop("diagnostic interruption"),
    .package = "pclvbayes"
  )

  expect_error(
    pclvbayes:::.sample_with_retry(
      mod = NULL,
      base_args = list(seed = 1L, chains = 1L, parallel_chains = 1L,
                       iter_warmup = 10L, iter_sampling = 10L,
                       adapt_delta = 0.8, max_treedepth = 10L,
                       metric = "diag_e", init = 0.2),
      stan_list = list(), max_retries = 0L, silent_sampler = TRUE
    ),
    "diagnostic interruption"
  )
  expect_false(dir.exists(owned))
})
