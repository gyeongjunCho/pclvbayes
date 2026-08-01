args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
source(file.path(here, "..", "mtist_adapter.R"))
out <- file.path(here, "results", "multichain_canonical")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_con <- file(file.path(out, "console.log"), "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({sink(type = "message"); sink(type = "output"); close(log_con)}, add = TRUE)
devtools::load_all(quiet = TRUE)

study <- load_mtist_study(361L)
mat_rel <- t(study$relative)
meta <- data.frame(Sample = rownames(study$sample_metadata),
  subject = study$sample_metadata$subject, time = study$sample_metadata$time,
  stringsAsFactors = FALSE)
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
vars <- c("a_ij", "a_ii", "sigma", "sd_ou", "phi", "lambda", "nu")
settings <- list(chains = 4L, parallel_chains = 4L, iter_warmup = 3000L,
  iter_sampling = 2000L, init = 0.2, adapt_delta = 0.98,
  max_treedepth = 14L, refresh = 0L, show_messages = FALSE)

# Compile before either elapsed timer starts.
model <- cmdstanr::cmdstan_model(file.path(here, "..", "..", "..",
                                          "inst", "stan", "pclv.stan"), quiet = TRUE)
ebfmi <- function(x) mean(diff(x)^2) / stats::var(x)
summarise_draw_set <- function(x) {
  z <- posterior::summarise_draws(x, mean, median, sd,
    ~posterior::quantile2(.x, probs = c(.05, .95)), posterior::rhat,
    posterior::ess_bulk, posterior::ess_tail)
  names(z) <- sub("^posterior::", "", names(z))
  names(z)[names(z) == "5%"] <- "q05"; names(z)[names(z) == "95%"] <- "q95"
  as.data.frame(z)
}

all_results <- list()
for (fixture in names(fixtures)) {
  cat("Sampling ", fixture, "\n", sep = "")
  started <- proc.time()[["elapsed"]]
  fit <- do.call(model$sample, c(list(data = fixtures[[fixture]],
    seed = unname(seeds[[fixture]])), settings))
  elapsed <- proc.time()[["elapsed"]] - started
  draws <- posterior::as_draws_df(fit$draws(vars))
  diag <- posterior::as_draws_df(fit$sampler_diagnostics())
  pooled <- transform(summarise_draw_set(fit$draws(vars)), fixture = fixture)
  pooled$ess_per_second <- pooled$ess_bulk / elapsed

  chain_summaries <- list(); chain_diag <- list(); correlations <- list()
  for (ch in 1:4) {
    dc <- draws[draws$.chain == ch, , drop = FALSE]
    sc <- transform(summarise_draw_set(posterior::as_draws_array(
      posterior::subset_draws(fit$draws(vars), chain = ch))), fixture = fixture,
      chain = ch)
    sc$sign_probability <- vapply(sc$variable, function(v) {
      value <- dc[[v]]; max(mean(value > 0), mean(value < 0))
    }, numeric(1))
    chain_summaries[[ch]] <- sc
    dg <- diag[diag$.chain == ch, , drop = FALSE]
    n_draw <- nrow(dg); div <- sum(dg$divergent__); td <- sum(dg$treedepth__ >= 14)
    pass <- max(sc$rhat, na.rm = TRUE) < 1.05 && min(sc$ess_bulk) > 400 &&
      min(sc$ess_tail) > 400 && div <= ceiling(.001 * n_draw) &&
      td <= ceiling(.01 * n_draw) && ebfmi(dg$energy__) >= .30
    chain_diag[[ch]] <- data.frame(fixture, chain = ch, divergences = div,
      max_treedepth_hits = td, average_treedepth = mean(dg$treedepth__),
      step_size = unique(dg$stepsize__)[1], ebfmi = ebfmi(dg$energy__),
      max_split_rhat = max(sc$rhat, na.rm = TRUE), min_bulk_ess = min(sc$ess_bulk),
      min_tail_ess = min(sc$ess_tail), passes_all_gates = pass)
    pairs <- list(c("sigma","sd_ou"), c("sd_ou","phi"), c("sd_ou","lambda"),
                  c("sigma","nu"), c("a_ii","a_ij"))
    correlations[[ch]] <- do.call(rbind, lapply(pairs, function(p)
      data.frame(fixture, scope = paste0("chain_", ch),
        pair = paste(p, collapse = ":"), correlation = cor(dc[[p[1]]], dc[[p[2]]]))))
  }
  chain_summaries <- do.call(rbind, chain_summaries)
  chain_diag <- do.call(rbind, chain_diag)
  correlations <- do.call(rbind, correlations)
  pairs <- list(c("sigma","sd_ou"), c("sd_ou","phi"), c("sd_ou","lambda"),
                c("sigma","nu"), c("a_ii","a_ij"))
  correlations <- rbind(correlations, do.call(rbind, lapply(pairs, function(p)
    data.frame(fixture, scope = "pooled", pair = paste(p, collapse = ":"),
      correlation = cor(draws[[p[1]]], draws[[p[2]]])))))
  a <- draws$a_ij
  by_chain <- split(seq_along(a), draws$.chain)
  interaction <- data.frame(fixture, scope = c("pooled", paste0("chain_", 1:4)),
    mean = c(mean(a), vapply(by_chain, function(ix) mean(a[ix]), numeric(1))),
    median = c(median(a), vapply(by_chain, function(ix) median(a[ix]), numeric(1))),
    sign_probability = c(max(mean(a > 0), mean(a < 0)),
      vapply(by_chain, function(ix) max(mean(a[ix] > 0), mean(a[ix] < 0)), numeric(1))),
    positive_probability = c(mean(a > 0), vapply(by_chain, function(ix) mean(a[ix] > 0), numeric(1))))
  timing <- tryCatch(as.data.frame(fit$time()), error = function(e) data.frame())
  fit_dir <- file.path(out, fixture); dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  fit$save_output_files(dir = fit_dir)
  all_results[[fixture]] <- list(elapsed = data.frame(fixture, elapsed_seconds = elapsed,
    total_divergences = sum(chain_diag$divergences), total_treedepth_hits = sum(chain_diag$max_treedepth_hits),
    max_rhat = max(pooled$rhat), min_bulk_ess = min(pooled$ess_bulk),
    min_tail_ess = min(pooled$ess_tail), min_ess_per_second = min(pooled$ess_per_second),
    chains_passing_all_gates = sum(chain_diag$passes_all_gates)), pooled = pooled,
    chains = chain_summaries, diagnostics = chain_diag, correlations = correlations,
    interaction = interaction, timing = timing)
}

bind <- function(part) do.call(rbind, lapply(all_results, `[[`, part))
for (part in c("elapsed","pooled","chains","diagnostics","correlations","interaction"))
  utils::write.table(bind(part), file.path(out, paste0(part, ".tsv")), sep = "\t",
                     row.names = FALSE, quote = FALSE)
timing <- lapply(names(all_results), function(nm) {
  x <- all_results[[nm]]$timing; if (!nrow(x)) return(NULL); x$fixture <- nm; x
})
timing <- Filter(Negate(is.null), timing)
if (length(timing)) utils::write.table(do.call(rbind, timing), file.path(out, "timing.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
jsonlite::write_json(list(dataset = 361L, model = "canonical_A", seeds = seeds,
  chain_seed_policy = "CmdStan deterministic chain streams from recorded base seed",
  settings = settings, compilation_excluded = TRUE,
  gates = list(rhat_lt = 1.05, ess_gt = 400, ebfmi_gte = .30,
    divergence_max_fraction = .001, treedepth_max_fraction = .01)),
  file.path(out, "config.json"), auto_unbox = TRUE, pretty = TRUE)
jsonlite::write_json(list(timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  R_version = R.version.string, CmdStan_version = as.character(cmdstanr::cmdstan_version()),
  OS = unname(Sys.info()["sysname"]), release = unname(Sys.info()["release"]),
  logical_cores = parallel::detectCores(logical = TRUE)), file.path(out, "environment.json"),
  auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(bind("elapsed")); print(bind("diagnostics")); print(bind("interaction"))
