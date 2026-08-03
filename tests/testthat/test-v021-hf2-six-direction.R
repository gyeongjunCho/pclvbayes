test_that("HF2 six-direction manifest and runners are frozen and current", {
  repo <- normalizePath(file.path(testthat::test_path(), "../.."))
  root <- file.path(repo, "benchmarks", "mtist", "results",
                    "v021_hf2_six_direction_revalidation_v1")
  manifest_path <- file.path(root, "direction_manifest.csv")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  expected <- data.frame(
    source=c("species_1","species_2","species_3","species_7","species_8","species_5"),
    target=c("species_7","species_9","species_7","species_4","species_4","species_7"),
    task_id=c(15L,24L,28L,33L,34L,37L), direction_index=c(30L,48L,56L,65L,67L,74L),
    seed=c(20468804L,20570804L,20668804L,20768803L,20769803L,20868804L))
  expect_identical(manifest, expected)
  hash <- strsplit(system2("sha256sum", manifest_path, stdout=TRUE)[[1L]], " ")[[1L]][[1L]]
  expect_identical(hash, "9e6907aa47953146f56e08b7e78f3c6f57619e8bc4ee2f31f0190f3d6364c005")
  runners <- file.path(repo, "benchmarks", "mtist", c(
    "run_v021_hf_oracle_revalidation.R", "run_v021_hf_six_direction_inference.R",
    "run_v021_hf_six_direction_report.R"))
  expect_true(all(vapply(runners, function(path) !inherits(try(parse(path), silent=TRUE), "try-error"), logical(1))))
  text <- paste(unlist(lapply(runners, readLines, warn=FALSE)), collapse="\n")
  expect_match(text, "v021_hf2_six_direction_revalidation_v1", fixed=TRUE)
  expect_false(grepl("elpd_pointwise_cross|elpd_pointwise_self|pseudo_BMA|pseudo_BMA_plus", text))
  expect_match(text, "truth[d$target, d$source]", fixed=TRUE)
  expect_match(text, "approved_runtime$kfold_seed <- as.integer(runtime$ctx$kfold_seed)",
               fixed=TRUE)
})
