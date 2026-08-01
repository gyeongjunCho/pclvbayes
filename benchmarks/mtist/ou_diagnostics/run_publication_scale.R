args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
source(file.path(here, "..", "mtist_adapter.R"))
out <- file.path(here, "results", "publication_scale")
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
  list(N = n, y = pair$y, xi = pair$xi, xj = pair$xj,
       S = length(unique(pair$subject)), sid = as.integer(factor(pair$subject)),
       prev = prev, dt = dt)
}

fixtures <- list(
  failed_species2_to_species0 = make_data("species_0", "species_2"),
  successful_species2_to_species1 = make_data("species_1", "species_2"))
seeds <- c(failed_species2_to_species0 = 20363802L,
           successful_species2_to_species1 = 20463802L)

# Compile both models before starting any per-fit elapsed timer.
models <- list(
  A_current = cmdstanr::cmdstan_model(file.path(here, "..", "..", "..",
                                                "inst", "stan", "pclv.stan"), quiet = TRUE),
  B_decay_equivalent = cmdstanr::cmdstan_model(file.path(here, "decay_equivalent.stan"), quiet = TRUE))
settings <- list(chains = 1L, parallel_chains = 1L, iter_warmup = 3000L,
                 iter_sampling = 8000L, init = 0.2, adapt_delta = 0.98,
                 max_treedepth = 14L, refresh = 0L, show_messages = FALSE)

ebfmi <- function(x) mean(diff(x)^2) / stats::var(x)
summarize_fit <- function(fit, elapsed, fixture, variant) {
  vars <- c("a_ij", "a_ii", "sigma", "sd_ou", "phi", "lambda", "nu")
  draws <- posterior::as_draws_df(fit$draws(vars))
  smry <- fit$summary(variables = vars)
  diag <- posterior::as_draws_df(fit$sampler_diagnostics())
  overall <- data.frame(
    fixture, variant, elapsed_seconds = elapsed,
    divergences = sum(diag$divergent__),
    max_treedepth_hits = sum(diag$treedepth__ >= 14),
    average_treedepth = mean(diag$treedepth__),
    step_size = unique(diag$stepsize__)[1], ebfmi = ebfmi(diag$energy__),
    min_bulk_ess = min(smry$ess_bulk), min_tail_ess = min(smry$ess_tail),
    min_ess_per_second = min(smry$ess_bulk) / elapsed,
    interaction_sign_probability = max(mean(draws$a_ij > 0), mean(draws$a_ij < 0)))
  list(overall = overall,
       posterior = transform(smry, fixture = fixture, variant = variant,
                             ess_per_second = ess_bulk / elapsed))
}

results <- list()
for (fixture in names(fixtures)) for (variant in names(models)) {
  key <- paste(fixture, variant, sep = "__")
  message("Sampling ", key)
  started <- proc.time()[["elapsed"]]
  fit <- do.call(models[[variant]]$sample,
                 c(list(data = fixtures[[fixture]], seed = unname(seeds[[fixture]])), settings))
  elapsed <- proc.time()[["elapsed"]] - started
  results[[key]] <- summarize_fit(fit, elapsed, fixture, variant)
  fit_dir <- file.path(out, key); dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  fit$save_output_files(dir = fit_dir)
}
bind_part <- function(part) do.call(rbind, lapply(results, `[[`, part))
utils::write.table(bind_part("overall"), file.path(out, "diagnostics.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(bind_part("posterior"), file.path(out, "posteriors.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
jsonlite::write_json(list(dataset = 361L, directions = names(fixtures), seeds = seeds,
  variants = names(models), settings = settings, compilation_excluded = TRUE),
  file.path(out, "config.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(bind_part("overall"))
