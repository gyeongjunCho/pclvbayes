args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
source(file.path(here, "..", "mtist_adapter.R"))
out <- file.path(here, "results", "variance_allocation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
devtools::load_all(quiet = TRUE)

study <- load_mtist_study(361L)
mat_rel <- t(study$relative)
meta <- data.frame(Sample = rownames(study$sample_metadata),
                   subject = study$sample_metadata$subject,
                   time = study$sample_metadata$time, stringsAsFactors = FALSE)
sm <- pclvbayes:::.precompute_spline_smoothed(mat_rel, meta, study$taxa, 1e-6, 3L)
if (inherits(sm, "pclv_failure")) stop("Canonical pre-smoothing failed")

make_data <- function(target, source) {
  pair <- pclvbayes:::.make_pair_inputs_glv(
    sm_mat = sm, meta_df = meta, j = source, i = target, min_pairs = 4L,
    zero_mode_alr = "minpos_time", minpos_alpha = 0.5, minpos_base = "ij",
    eps_fixed = 1e-6, lib_eps_c = 0.65, rest_floor_frac = 1,
    alr_cap = 12, smooth_scale = "logra", alr_spline_cv = TRUE,
    nz_partner_min_frac = 0.15)
  if (inherits(pair, "pclv_failure")) stop("Canonical pair preprocessing failed")
  pair <- pair[order(pair$subject, pair$time), , drop = FALSE]
  n <- nrow(pair); prev <- integer(n); dt <- numeric(n)
  for (ix in split(seq_len(n), pair$subject)) if (length(ix) > 1L) {
    prev[ix[-1L]] <- ix[-length(ix)]; dt[ix[-1L]] <- diff(pair$time[ix])
  }
  list(pair = pair, stan = list(N = n, y = pair$y, xi = pair$xi, xj = pair$xj,
    S = length(unique(pair$subject)), sid = as.integer(factor(pair$subject)),
    prev = prev, dt = dt))
}

fixtures <- list(
  failed_species2_to_species0 = make_data("species_0", "species_2"),
  successful_species2_to_species1 = make_data("species_1", "species_2"))
seeds <- c(failed_species2_to_species0 = 20363802L,
           successful_species2_to_species1 = 20463802L)
models <- list(
  B_decay_equivalent = cmdstanr::cmdstan_model(file.path(here, "decay_equivalent.stan"), quiet = TRUE),
  D_variance_allocation = cmdstanr::cmdstan_model(file.path(here, "variance_allocation.stan"), quiet = TRUE))
base_settings <- list(chains = 2L, parallel_chains = 2L, iter_warmup = 1000L,
  iter_sampling = 2000L, init = 0.2, adapt_delta = 0.98,
  max_treedepth = 14L, refresh = 0L, show_messages = FALSE)

ebfmi <- function(x) mean(diff(x)^2) / stats::var(x)
summarize_one <- function(fit, elapsed, fixture, variant) {
  vars <- c("a_ij", "a_ii", "sigma", "sd_ou", "phi", "lambda", "nu")
  draws <- posterior::as_draws_df(fit$draws(vars))
  draws$s_total <- sqrt(draws$sigma^2 + draws$sd_ou^2)
  draws$omega <- draws$sd_ou^2 / draws$s_total^2
  smry <- posterior::summarise_draws(draws[, c(vars, "s_total", "omega")],
    mean, median, sd, ~posterior::quantile2(.x, probs = c(.05, .95)),
    posterior::rhat, posterior::ess_bulk, posterior::ess_tail)
  names(smry)[names(smry) == "5%"] <- "q05"
  names(smry)[names(smry) == "95%"] <- "q95"
  names(smry)[names(smry) == "posterior::rhat"] <- "rhat"
  names(smry)[names(smry) == "posterior::ess_bulk"] <- "ess_bulk"
  names(smry)[names(smry) == "posterior::ess_tail"] <- "ess_tail"
  diag <- posterior::as_draws_df(fit$sampler_diagnostics())
  chains <- unique(diag$.chain)
  bf <- vapply(chains, function(ch) ebfmi(diag$energy__[diag$.chain == ch]), numeric(1))
  corr <- function(a, b) stats::cor(draws[[a]], draws[[b]])
  corrs <- data.frame(fixture, variant,
    pair = c("sigma:sd_ou", "s_total:omega", "omega:lambda", "sd_ou:lambda",
             "sigma:nu", "a_ii:a_ij"),
    correlation = c(corr("sigma","sd_ou"), corr("s_total","omega"),
      corr("omega","lambda"), corr("sd_ou","lambda"), corr("sigma","nu"),
      corr("a_ii","a_ij")))
  overall <- data.frame(fixture, variant, elapsed_seconds = elapsed,
    divergences = sum(diag$divergent__),
    max_treedepth_hits = sum(diag$treedepth__ >= 14),
    average_treedepth = mean(diag$treedepth__),
    step_size_chain1 = unique(diag$stepsize__[diag$.chain == chains[1]])[1],
    step_size_chain2 = unique(diag$stepsize__[diag$.chain == chains[2]])[1],
    ebfmi_chain1 = bf[1], ebfmi_chain2 = bf[2],
    max_rhat = max(smry$rhat, na.rm = TRUE),
    min_bulk_ess = min(smry$ess_bulk), min_tail_ess = min(smry$ess_tail),
    min_ess_per_second = min(smry$ess_bulk) / elapsed,
    interaction_sign_probability = max(mean(draws$a_ij > 0), mean(draws$a_ij < 0)))
  list(overall = overall, summaries = transform(as.data.frame(smry), fixture = fixture, variant = variant),
       correlations = corrs)
}

results <- list(); fits <- list()
for (fixture in names(fixtures)) for (variant in names(models)) {
  key <- paste(fixture, variant, sep = "__")
  message("Sampling ", key)
  settings <- c(base_settings, list(seed = unname(seeds[[fixture]])))
  started <- proc.time()[["elapsed"]]
  fits[[key]] <- do.call(models[[variant]]$sample,
                         c(list(data = fixtures[[fixture]]$stan), settings))
  elapsed <- proc.time()[["elapsed"]] - started
  results[[key]] <- summarize_one(fits[[key]], elapsed, fixture, variant)
  fit_dir <- file.path(out, key); dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  fits[[key]]$save_output_files(dir = fit_dir)
}

bind_part <- function(part) do.call(rbind, lapply(results, `[[`, part))
utils::write.table(bind_part("overall"), file.path(out, "real_diagnostics.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(bind_part("summaries"), file.path(out, "real_posteriors.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(bind_part("correlations"), file.path(out, "real_correlations.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Exact implied canonical joint prior, evaluated in both coordinate systems.
set.seed(20363802L); m <- 300000L
prior <- data.frame(sigma = abs(stats::rnorm(m, 0, 0.5)),
                    sd_ou = abs(stats::rnorm(m, 0, 1)))
prior$s_total <- sqrt(prior$sigma^2 + prior$sd_ou^2)
prior$omega <- prior$sd_ou^2 / prior$s_total^2
prior$nu <- 2 + exp(stats::rnorm(m, log(3), 0.75))
prior$residual <- stats::rnorm(m, 0, prior$sd_ou) + prior$sigma * stats::rt(m, prior$nu)
q <- function(x) unname(stats::quantile(x, c(.025, .5, .975)))
prior_summary <- do.call(rbind, lapply(c("sigma","sd_ou","s_total","omega","residual"),
  function(v) data.frame(quantity = v, q025 = q(prior[[v]])[1],
    median = q(prior[[v]])[2], q975 = q(prior[[v]])[3])))
prior_extra <- data.frame(metric = c("cor_sigma_sd_ou", "residual_sd"),
                          value = c(cor(prior$sigma, prior$sd_ou), sd(prior$residual)))
utils::write.table(prior_summary, file.path(out, "prior_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(prior_extra, file.path(out, "prior_extra.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Three allocation regimes, with signs alternating across the two fixture designs.
cases <- data.frame(case = c("observation_dominant", "balanced", "ou_dominant"),
                    fixture = c(names(fixtures)[1], names(fixtures)[2], names(fixtures)[1]),
                    omega = c(.2, .5, .8), a_ij = c(.10, -.10, .10))
synthetic <- list()
for (cc in seq_len(nrow(cases))) {
  case <- cases[cc, ]; dat <- fixtures[[case$fixture]]$stan
  set.seed(30000000L + cc); dt_unit <- mean(dat$dt[dat$dt > 0]); total <- .10
  truth_sigma <- total * sqrt(1 - case$omega)
  truth_sd_ou <- total * sqrt(case$omega); truth_phi <- .8
  rsub <- stats::rnorm(dat$S, 0, .03); e <- numeric(dat$N)
  for (ii in seq_len(dat$N)) {
    if (dat$prev[ii] == 0L) e[ii] <- stats::rnorm(1, 0, truth_sd_ou)
    else {
      rho <- truth_phi^(dat$dt[ii] / dt_unit)
      e[ii] <- rho * e[dat$prev[ii]] + stats::rnorm(1, 0, truth_sd_ou * sqrt(1-rho^2))
    }
  }
  dat$y <- rsub[dat$sid] - .15 * dat$xi + case$a_ij * dat$xj + e +
    truth_sigma * stats::rt(dat$N, 5)
  syn_settings <- c(base_settings, list(seed = 30000000L + cc))
  started <- proc.time()[["elapsed"]]
  fit <- do.call(models$D_variance_allocation$sample, c(list(data = dat), syn_settings))
  elapsed <- proc.time()[["elapsed"]] - started
  ss <- summarize_one(fit, elapsed, case$case, "D_variance_allocation")
  ss$summaries$truth <- c(a_ij = case$a_ij, a_ii = -.15, sigma = truth_sigma,
    sd_ou = truth_sd_ou, phi = truth_phi, lambda = -log(truth_phi)/dt_unit,
    nu = 5, s_total = total, omega = case$omega)[ss$summaries$variable]
  ss$overall$true_interaction <- case$a_ij
  ss$overall$sign_recovered <- with(ss$summaries[ss$summaries$variable == "a_ij", ],
    sign(median) == sign(case$a_ij))
  synthetic[[case$case]] <- ss
  syn_dir <- file.path(out, paste0("synthetic_", case$case)); dir.create(syn_dir, recursive = TRUE, showWarnings = FALSE)
  fit$save_output_files(dir = syn_dir)
}
utils::write.table(do.call(rbind, lapply(synthetic, `[[`, "overall")), file.path(out, "synthetic_diagnostics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(do.call(rbind, lapply(synthetic, `[[`, "summaries")), file.path(out, "synthetic_recovery.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
jsonlite::write_json(list(dataset = 361L, directions = names(fixtures), seeds = seeds,
  settings = base_settings, prior = "exact_transformation_of_canonical_joint_half_normals"),
  file.path(out, "config.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(bind_part("overall")); print(do.call(rbind, lapply(synthetic, `[[`, "overall")))
