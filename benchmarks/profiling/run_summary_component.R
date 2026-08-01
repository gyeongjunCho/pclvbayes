source("benchmarks/profiling/profile_helpers.R")
devtools::load_all(quiet = TRUE)
set.seed(1)
n_iter <- 2000L; n_chain <- 4L
variables <- c("a_ij", "a_ii", "r0", "sigma", "sd_ou", "phi", "nu", "log_nu_minus_two")
draws <- posterior::as_draws_array(array(
  rnorm(n_iter * n_chain * length(variables)),
  dim = c(n_iter, n_chain, length(variables)), dimnames = list(NULL, NULL, variables)))
sampler <- posterior::as_draws_array(array(
  c(rbinom(n_iter * n_chain, 1, .001), sample(1:10, n_iter * n_chain, TRUE),
    rnorm(n_iter * n_chain)), dim = c(n_iter, n_chain, 3L),
  dimnames = list(NULL, NULL, c("divergent__", "treedepth__", "energy__"))))

draw_calls <- 0L
fit <- list(
  draws = function(variables = NULL) {
    draw_calls <<- draw_calls + 1L
    if (is.null(variables)) draws else posterior::subset_draws(draws, variable = variables)
  },
  sampler_diagnostics = function() sampler
)
legacy_attempt <- function() {
  all <- fit$draws(); use <- intersect(variables, posterior::variables(all))
  posterior::summarise_draws(fit$draws(use))
  posterior::as_draws_df(fit$sampler_diagnostics())
  posterior::ndraws(fit$draws())
}
optimized_attempt <- function() pclvbayes:::.summarise_sampler_diag(fit, 14L)
retained <- function() {
  d <- pclvbayes:::.safe_draws_df(fit)
  diag <- pclvbayes:::.add_convergence_diag(optimized_attempt(), d)
  pclvbayes:::.build_posterior_summary_bundle(d, diag)
}

invisible(legacy_attempt()); invisible(optimized_attempt()); invisible(retained())
measure <- function(fun, repetitions = 15L) {
  calls_before <- draw_calls
  times <- replicate(repetitions, system.time(fun())[["elapsed"]])
  data.frame(median_seconds = median(times), repetitions = repetitions,
             posterior_draw_calls = draw_calls - calls_before)
}
result <- rbind(
  cbind(path = "legacy_attempt", summarise_draws_calls_per_invocation = 1, measure(legacy_attempt)),
  cbind(path = "optimized_attempt", summarise_draws_calls_per_invocation = 0, measure(optimized_attempt)),
  cbind(path = "optimized_retained_final", summarise_draws_calls_per_invocation = 1, measure(retained))
)
out_dir <- file.path("benchmarks/profiling/results", paste0("summary-component-", format(Sys.time(), "%Y%m%d-%H%M%S")))
dir.create(out_dir, recursive = TRUE)
utils::write.csv(result, file.path(out_dir, "summary_component.csv"), row.names = FALSE)
cat("SUMMARY_COMPONENT_RESULT_DIR=", out_dir, "\n", sep = "")
print(result)
