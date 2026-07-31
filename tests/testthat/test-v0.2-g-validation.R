make_core_physeq <- function(subjects = paste0("S", 1:6), times = 0:2) {
  ids <- as.vector(outer(subjects, times, function(s, t) paste(s, t, sep = "_")))
  md <- expand.grid(subject = subjects, time = times, KEEP.OUT.ATTRS = FALSE)
  rownames(md) <- ids
  mat <- rbind(a = seq_along(ids) + 1, b = rev(seq_along(ids)) + 2)
  colnames(mat) <- ids
  phyloseq::phyloseq(phyloseq::otu_table(mat, taxa_are_rows = TRUE),
                     phyloseq::sample_data(md))
}

core_controls <- function() {
  f <- formals(fit_pclv_bayes)
  env <- new.env(parent = environment())
  out <- list()
  for (nm in setdiff(names(f), c("physeq", "subject_col", "time_col", "taxa_vec"))) {
    out[[nm]] <- eval(f[[nm]], envir = env)
    assign(nm, out[[nm]], envir = env)
  }
  out
}

test_that("invalid public inputs fail before model loading", {
  called <- FALSE
  local_mocked_bindings(get_pclv_model = function(...) { called <<- TRUE; stop("loaded") },
                        .package = "pclvbayes")
  expect_error(fit_pclv_bayes(list(), "subject", "time"), "phyloseq")
  expect_false(called)
})

test_that("metadata columns, times, and taxa are validated centrally", {
  ps <- make_core_physeq()
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps, "missing", "time", NULL, core_controls()), "Missing sample metadata")
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps, "subject", "missing", NULL, core_controls()), "Missing sample metadata")
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps, "subject", "time", c("a", "unknown"), core_controls()), "Unknown taxa requested: unknown")

  md <- data.frame(phyloseq::sample_data(ps)); ix <- which(md$subject == md$subject[1]); md$time[ix[2]] <- md$time[ix[1]]
  ps_dup <- phyloseq::phyloseq(phyloseq::otu_table(ps), phyloseq::sample_data(md))
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps_dup, "subject", "time", NULL, core_controls()), "Duplicate time values")
  md$time[ix[2]] <- Inf
  ps_inf <- phyloseq::phyloseq(phyloseq::otu_table(ps), phyloseq::sample_data(md))
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps_inf, "subject", "time", NULL, core_controls()), "finite")
})

test_that("abundance and scalar controls are rejected at the boundary", {
  ps <- make_core_physeq(); ctl <- core_controls()
  bad <- as(phyloseq::otu_table(ps), "matrix"); bad[1, 1] <- -1
  ps_bad <- phyloseq::phyloseq(phyloseq::otu_table(bad, taxa_are_rows = TRUE), phyloseq::sample_data(ps))
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps_bad, "subject", "time", NULL, ctl), "nonnegative")
  bad[1, 1] <- Inf
  ps_bad <- phyloseq::phyloseq(phyloseq::otu_table(bad, taxa_are_rows = TRUE), phyloseq::sample_data(ps))
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps_bad, "subject", "time", NULL, ctl), "finite")

  ctl$quiet <- NA
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps, "subject", "time", NULL, ctl), "`quiet`")
  ctl <- core_controls(); ctl$chains <- 1.5
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps, "subject", "time", NULL, ctl), "`chains`")
  ctl <- core_controls(); ctl$adapt_delta <- 1
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps, "subject", "time", NULL, ctl), "`adapt_delta`")
  ctl <- core_controls(); ctl$kfold_K <- 7
  expect_error(pclvbayes:::.validate_fit_pclv_inputs(ps, "subject", "time", NULL, ctl), "cannot exceed")
})

test_that("valid public inputs normalize deterministically", {
  ps <- make_core_physeq(subjects = c("B", "A", "C", "D", "E", "F"))
  ctl <- core_controls(); ctl$chains <- 2; ctl$kfold_K <- 5
  out <- pclvbayes:::.validate_fit_pclv_inputs(ps, "subject", "time", c("b", "a"), ctl)
  expect_identical(out$taxa_vec, c("b", "a"))
  expect_type(out$controls$chains, "integer")
  expect_true(all(diff(out$meta_df$time[out$meta_df$subject == "A"]) > 0))
  expect_identical(unique(as.character(out$meta_df$subject)), sort(unique(as.character(out$meta_df$subject))))
})

test_that("scientific eligibility remains a structured runtime failure", {
  out <- pclvbayes:::.validate_predictor_variation(rep(1, 4), 1:4, "full_data")
  expect_s3_class(out, "pclv_failure")
  expect_identical(out$reason, "insufficient_predictor_variation")
})
