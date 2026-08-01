args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
source(file.path(here, "..", "mtist_adapter.R"))
dir.create(file.path(here, "results"), recursive = TRUE, showWarnings = FALSE)

devtools::load_all(quiet = TRUE)
study <- load_mtist_study(361L)
taxa <- study$taxa
mat_rel <- t(study$relative)
meta <- data.frame(
  Sample = rownames(study$sample_metadata),
  subject = study$sample_metadata$subject,
  time = study$sample_metadata$time,
  stringsAsFactors = FALSE
)
sm <- pclvbayes:::.precompute_spline_smoothed(mat_rel, meta, taxa, 1e-6, 3L)
if (inherits(sm, "pclv_failure")) stop("Canonical pre-smoothing failed")
pair <- pclvbayes:::.make_pair_inputs_glv(
  sm_mat = sm, meta_df = meta, j = "species_2", i = "species_0",
  min_pairs = 4L, zero_mode_alr = "minpos_time", minpos_alpha = 0.5,
  minpos_base = "ij", eps_fixed = 1e-6, lib_eps_c = 0.65,
  rest_floor_frac = 1, alr_cap = 12, smooth_scale = "logra",
  alr_spline_cv = TRUE, nz_partner_min_frac = 0.15
)
if (inherits(pair, "pclv_failure")) stop("Canonical pair preprocessing failed")
pair <- pair[order(pair$subject, pair$time), , drop = FALSE]
n <- nrow(pair); prev <- integer(n); dt <- numeric(n)
for (ix in split(seq_len(n), pair$subject)) {
  if (length(ix) > 1L) {
    prev[ix[-1L]] <- ix[-length(ix)]
    dt[ix[-1L]] <- diff(pair$time[ix])
  }
}
stan_data <- list(N = n, y = pair$y, xi = pair$xi, xj = pair$xj,
                  S = length(unique(pair$subject)),
                  sid = as.integer(factor(pair$subject)), prev = prev, dt = dt)
utils::write.table(cbind(pair, sid = stan_data$sid, prev, dt),
                   file.path(here, "results", "fixture.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)

seed <- 20363802L
settings <- list(chains = 2L, parallel_chains = 2L, iter_warmup = 1000L,
                 iter_sampling = 2000L, seed = seed, init = 0.2,
                 adapt_delta = 0.98, max_treedepth = 14L, refresh = 0L,
                 show_messages = FALSE)
models <- list(
  A_current = cmdstanr::cmdstan_model(file.path(here, "..", "..", "..", "inst", "stan", "pclv.stan"), quiet = TRUE),
  B_decay_equivalent = cmdstanr::cmdstan_model(file.path(here, "decay_equivalent.stan"), quiet = TRUE),
  C_decay_regularized = cmdstanr::cmdstan_model(file.path(here, "decay_regularized.stan"), quiet = TRUE)
)

ebfmi <- function(e) mean(diff(e)^2) / stats::var(e)
summarize_fit <- function(fit, elapsed, variant) {
  vars <- c("a_ij", "a_ii", "sigma", "sd_ou", "phi", "lambda", "nu")
  smry <- fit$summary(variables = vars)
  draws <- posterior::as_draws_df(fit$draws(vars))
  diag <- posterior::as_draws_df(fit$sampler_diagnostics())
  chain_ids <- unique(diag$.chain)
  bfmi <- vapply(chain_ids, function(ch) ebfmi(diag$energy__[diag$.chain == ch]), numeric(1))
  cors <- data.frame(
    variant = variant,
    pair = c("sd_ou:phi", "sd_ou:lambda", "sigma:sd_ou", "sigma:nu", "a_ii:a_ij"),
    correlation = c(stats::cor(draws$sd_ou, draws$phi),
                    stats::cor(draws$sd_ou, draws$lambda),
                    stats::cor(draws$sigma, draws$sd_ou),
                    stats::cor(draws$sigma, draws$nu),
                    stats::cor(draws$a_ii, draws$a_ij))
  )
  detail <- transform(smry, variant = variant, elapsed_seconds = elapsed,
                      ess_per_second = ess_bulk / elapsed)
  overall <- data.frame(
    variant = variant, elapsed_seconds = elapsed,
    divergences = sum(diag$divergent__),
    max_treedepth_hits = sum(diag$treedepth__ >= 14),
    average_treedepth = mean(diag$treedepth__),
    step_size_chain1 = unique(diag$stepsize__[diag$.chain == chain_ids[1]])[1],
    step_size_chain2 = unique(diag$stepsize__[diag$.chain == chain_ids[2]])[1],
    ebfmi_chain1 = bfmi[1], ebfmi_chain2 = bfmi[2],
    max_rhat = max(smry$rhat, na.rm = TRUE), min_bulk_ess = min(smry$ess_bulk),
    min_tail_ess = min(smry$ess_tail),
    interaction_sign_probability = max(mean(draws$a_ij > 0), mean(draws$a_ij < 0)),
    mean_total_log_lik = mean(rowSums(posterior::as_draws_matrix(fit$draws("log_lik"))))
  )
  list(overall = overall, parameters = detail, correlations = cors)
}

all <- vector("list", length(models)); names(all) <- names(models)
fits <- vector("list", length(models)); names(fits) <- names(models)
for (variant in names(models)) {
  message("Sampling ", variant)
  started <- proc.time()[["elapsed"]]
  fits[[variant]] <- do.call(models[[variant]]$sample, c(list(data = stan_data), settings))
  elapsed <- proc.time()[["elapsed"]] - started
  all[[variant]] <- summarize_fit(fits[[variant]], elapsed, variant)
  variant_dir <- file.path(here, "results", variant)
  dir.create(variant_dir, recursive = TRUE, showWarnings = FALSE)
  fits[[variant]]$save_output_files(dir = variant_dir)
}
utils::write.table(do.call(rbind, lapply(all, `[[`, "overall")),
  file.path(here, "results", "diagnostic_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(do.call(rbind, lapply(all, `[[`, "parameters")),
  file.path(here, "results", "posterior_summaries.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(do.call(rbind, lapply(all, `[[`, "correlations")),
  file.path(here, "results", "posterior_correlations.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Quantitative prior implications for C, with no hard lower scale bound.
set.seed(seed)
m <- 200000L; dt_pos <- dt[dt > 0]; dt_unit <- mean(dt_pos)
prior <- data.frame(
  sigma = exp(stats::rnorm(m, log(0.02), 0.75)),
  sd_ou = exp(stats::rnorm(m, log(0.10), 0.75)),
  lambda = exp(stats::rnorm(m, log(-log(0.80) / dt_unit), 0.60))
)
prior$phi_dt_unit <- exp(-prior$lambda * dt_unit)
prior$rho_min_dt <- exp(-prior$lambda * min(dt_pos))
prior$rho_max_dt <- exp(-prior$lambda * max(dt_pos))
prior$residual_draw <- stats::rnorm(m, 0, prior$sd_ou) +
  prior$sigma * stats::rt(m, df = 2 + exp(stats::rnorm(m, log(3), 0.75)))
qfun <- function(x) unname(stats::quantile(x, c(.025, .5, .975)))
prior_summary <- do.call(rbind, lapply(names(prior), function(v)
  data.frame(quantity = v, q025 = qfun(prior[[v]])[1], median = qfun(prior[[v]])[2], q975 = qfun(prior[[v]])[3])))
empirical <- data.frame(quantity = c("response_sd", "response_mad", "response_min", "response_max", "dt_unit"),
                        value = c(sd(pair$y), mad(pair$y), min(pair$y), max(pair$y), dt_unit))
utils::write.table(prior_summary, file.path(here, "results", "prior_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(empirical, file.path(here, "results", "fixture_scale.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Simple identifiable synthetic recovery using the same irregular-time design.
set.seed(seed + 1L)
truth <- c(r0 = 0, sd_r0 = 0.03, a_ii = -0.15, a_ij = 0.10,
           sigma = 0.02, sd_ou = 0.08, phi = 0.80, nu = 5)
rsub <- stats::rnorm(stan_data$S, 0, truth[["sd_r0"]]); e <- numeric(n)
for (ii in seq_len(n)) {
  if (prev[ii] == 0L) e[ii] <- stats::rnorm(1, 0, truth[["sd_ou"]])
  else {
    rho <- truth[["phi"]]^(dt[ii] / dt_unit)
    e[ii] <- rho * e[prev[ii]] + stats::rnorm(1, 0, truth[["sd_ou"]] * sqrt(1-rho^2))
  }
}
syn <- stan_data
syn$y <- truth[["r0"]] + rsub[syn$sid] + truth[["a_ii"]] * syn$xi +
  truth[["a_ij"]] * syn$xj + e + truth[["sigma"]] * stats::rt(n, df = truth[["nu"]])
syn_settings <- settings; syn_settings$iter_warmup <- 500L; syn_settings$iter_sampling <- 1000L
started <- proc.time()[["elapsed"]]
syn_fit <- do.call(models$C_decay_regularized$sample, c(list(data = syn), syn_settings))
syn_elapsed <- proc.time()[["elapsed"]] - started
syn_sum <- syn_fit$summary(variables = c("a_ij", "a_ii", "sigma", "sd_ou", "phi", "nu"))
syn_draw <- posterior::as_draws_df(syn_fit$draws("a_ij"))
syn_diag <- posterior::as_draws_df(syn_fit$sampler_diagnostics())
synthetic_result <- data.frame(
  true_a_ij = truth[["a_ij"]], posterior_a_ij_median = syn_sum$median[syn_sum$variable == "a_ij"],
  posterior_a_ij_q05 = syn_sum$q5[syn_sum$variable == "a_ij"], posterior_a_ij_q95 = syn_sum$q95[syn_sum$variable == "a_ij"],
  probability_positive = mean(syn_draw$a_ij > 0), sign_recovered = mean(syn_draw$a_ij > 0) > 0.95,
  divergences = sum(syn_diag$divergent__), max_treedepth_hits = sum(syn_diag$treedepth__ >= 14),
  elapsed_seconds = syn_elapsed
)
utils::write.table(synthetic_result, file.path(here, "results", "synthetic_recovery.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(transform(syn_sum, truth = truth[match(variable, names(truth))]),
                   file.path(here, "results", "synthetic_posterior.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
synthetic_dir <- file.path(here, "results", "synthetic_C")
dir.create(synthetic_dir, recursive = TRUE, showWarnings = FALSE)
syn_fit$save_output_files(dir = synthetic_dir)
writeLines(capture.output(sessionInfo()), file.path(here, "results", "sessionInfo.txt"))
jsonlite::write_json(list(dataset_id = 361L, target = "species_0", source = "species_2",
  seed = seed, settings = settings, synthetic_settings = syn_settings,
  package_git_sha = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[1],
  mtist_git_sha = mtist_git_sha()), file.path(here, "results", "config.json"), auto_unbox = TRUE, pretty = TRUE)
print(do.call(rbind, lapply(all, `[[`, "overall")))
print(synthetic_result)
