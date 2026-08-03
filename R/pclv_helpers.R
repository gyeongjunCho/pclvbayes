.PCLV_CORE_ALR_CAP <- 12
.PCLV_CORE_TRANSFORM <- "alr"
.PCLV_CORE_LAG <- 1L
.PCLV_CORE_RESID_MODE <- "ou"
.PCLV_CORE_USE_STUDENT_T <- TRUE
.PCLV_CORE_COMPUTE_ELPD <- TRUE
.PCLV_CORE_ELPD_MODE <- "student_t_scale_mixture_kalman"
.PCLV_CORE_T_KALMAN_QUADRATURE <- list(
  probability = c(
    0.005299532504175031, 0.027712488463383700,
    0.067184398806084122, 0.122297795822498500,
    0.191061877798678110, 0.270991611171386320,
    0.359198224610370540, 0.452493745081181290,
    0.547506254918818770, 0.640801775389629460,
    0.729008388828613630, 0.808938122201321890,
    0.877702204177501550, 0.932815601193915930,
    0.972287511536616300, 0.994700467495824970
  ),
  weight = c(
    0.013576229705877088, 0.031126761969323728,
    0.047579255841246303, 0.062314485627767036,
    0.074797994408288354, 0.084578259697501323,
    0.091301707522461820, 0.094725305227534320,
    0.094725305227534320, 0.091301707522461820,
    0.084578259697501323, 0.074797994408288354,
    0.062314485627767036, 0.047579255841246303,
    0.031126761969323728, 0.013576229705877088
  )
)
.PCLV_CORE_SPLINE <- list(df = NULL, spar = NULL, cv = TRUE)

.pclv_failure <- function(stage, reason, details = list()) {
  out <- c(list(ok = FALSE, stage = stage, reason = reason, details = details), details)
  out <- out[!duplicated(names(out))]
  structure(out, class = c("pclv_failure", "list"))
}

.is_pclv_failure <- function(x) inherits(x, "pclv_failure")

.predictor_variation_failure <- function(x, predictor, stage) {
  finite_x <- x[is.finite(x)]
  scale <- max(c(1, abs(finite_x)))
  required_sd <- sqrt(.Machine$double.eps) * scale
  observed_sd <- stats::sd(x)

  if (is.finite(observed_sd) && observed_sd > required_sd) {
    return(NULL)
  }

  .pclv_failure(stage, "insufficient_predictor_variation", list(
    predictor = predictor, observed_sd = observed_sd, required_sd = required_sd
  ))
}

.validate_predictor_variation <- function(xi, xj, stage) {
  failure <- .predictor_variation_failure(xi, "xi", stage)
  if (!is.null(failure)) {
    return(failure)
  }

  .predictor_variation_failure(xj, "xj", stage)
}

#' Per-subject spline smoothing of relative abundances (log-scale)
#'
#' Every taxon in the complete input community is smoothed independently within
#' subject on \code{log(pmax(x, 0) + eps)}. The reconstructed abundances are
#' then reclosed sample-wise over the complete community before the requested
#' taxa are returned. Reclosure is required because independent smoothing does
#' not preserve compositional column sums.
#'
#' @param mat_rel Taxa-by-samples relative abundance matrix (rows = taxa).
#' @param meta_df Data frame with columns \code{Sample}, \code{subject}, \code{time}.
#' @param taxa_list Character vector of taxa to return after full-community smoothing.
#' @param eps Small constant for log transform stability.
#' @param min_unique_times Minimum unique time points required to fit a spline.
#' @return A requested-taxa-by-samples numeric matrix whose values were
#'   normalized using the complete smoothed community.
#' @noRd
#' @keywords internal
.precompute_spline_smoothed <- function(mat_rel, meta_df,
                                        taxa_list = rownames(mat_rel),
                                        eps = 1e-6,
                                        min_unique_times = 3) {
  time_failure <- .validate_subject_times(
    meta_df$time, meta_df$subject, "spline_smoothing"
  )
  if (!is.null(time_failure)) return(time_failure)

  matrix_taxa <- rownames(mat_rel)
  matrix_samples <- colnames(mat_rel)
  metadata_samples <- as.character(meta_df$Sample)
  taxa_list <- as.character(taxa_list)

  if (is.null(matrix_taxa) || is.null(matrix_samples) ||
      anyNA(matrix_taxa) || anyNA(matrix_samples) ||
      anyDuplicated(matrix_taxa) || anyDuplicated(matrix_samples)) {
    stop("Validated abundance dimname invariant violated.")
  }
  if (!length(taxa_list) || anyNA(taxa_list) || any(!nzchar(taxa_list)) ||
      anyDuplicated(taxa_list) || !all(taxa_list %in% matrix_taxa)) {
    stop("Validated taxa invariant violated.")
  }
  if (length(metadata_samples) != length(matrix_samples) ||
      anyNA(metadata_samples) || any(!nzchar(metadata_samples)) ||
      anyDuplicated(metadata_samples) ||
      !setequal(metadata_samples, matrix_samples)) {
    return(.pclv_failure(
      "spline_smoothing", "sample_alignment_failed",
      list(
        metadata_n = length(metadata_samples),
        matrix_n = length(matrix_samples)
      )
    ))
  }

  col_idx_all <- match(metadata_samples, matrix_samples)
  if (anyNA(col_idx_all) || anyDuplicated(col_idx_all)) {
    return(.pclv_failure(
      "spline_smoothing", "sample_alignment_failed", list()
    ))
  }

  sm_mat <- matrix(
    NA_real_, nrow = nrow(mat_rel), ncol = ncol(mat_rel),
    dimnames = dimnames(mat_rel)
  )

  # The denominator represents the complete remaining community. Therefore all
  # taxa must participate in smoothing and reclosure, even when only a subset is
  # requested for directed fitting.
  taxa_use <- matrix_taxa

  for (tx in taxa_use) {
    vec_pred <- rep(NA_real_, nrow(meta_df))

    for (sb in unique(meta_df$subject)) {
      idx   <- which(meta_df$subject == sb)
      times <- meta_df$time[idx]
      cols  <- col_idx_all[idx]
      vals  <- as.numeric(mat_rel[tx, cols])

      ok <- is.finite(times) & is.finite(vals)
      if (!any(ok)) {
        return(.pclv_failure(
          "spline_smoothing", "no_finite_abundance",
          list(taxon = tx, subject = sb)
        ))
      }

      df   <- data.frame(time = times[ok], val = vals[ok])
      df2  <- stats::aggregate(val ~ time, df, mean)
      df2  <- df2[order(df2$time), , drop = FALSE]

      pred_log <- NULL
      if (nrow(df2) >= min_unique_times) {
        ylog <- log(pmax(df2$val, 0) + eps)
        rr <- .smooth_spline_robust(
          x = df2$time, y = ylog,
          spline_df = .PCLV_CORE_SPLINE$df,
          spline_spar = .PCLV_CORE_SPLINE$spar,
          use_cv = .PCLV_CORE_SPLINE$cv,
          min_unique = min_unique_times,
          min_df = 3.0,
          max_df = NULL,
          default_df = 4.0
        )
        if (.is_pclv_failure(rr)) {
          rr$details$taxon <- tx
          rr$details$subject <- sb
          return(rr)
        }
        pred_log <- rr$yhat[match(times, df2$time)]
      }

      if (is.null(pred_log)) {
        return(.pclv_failure(
          "spline_smoothing", "insufficient_unique_times",
          list(taxon = tx, subject = sb)
        ))
      }

      vec_pred[idx] <- pmax(exp(pred_log) - eps, 0)
    }

    sm_mat[tx, col_idx_all] <- vec_pred
  }

  invalid_value <- !is.finite(sm_mat) | sm_mat < 0
  if (any(invalid_value)) {
    bad <- which(invalid_value, arr.ind = TRUE)[1L, , drop = FALSE]
    return(.pclv_failure(
      "spline_smoothing", "invalid_smoothed_abundance",
      list(
        taxon = rownames(sm_mat)[bad[1L, "row"]],
        sample = colnames(sm_mat)[bad[1L, "col"]]
      )
    ))
  }

  column_totals <- colSums(sm_mat)
  bad_total <- !is.finite(column_totals) | column_totals <= 0
  if (any(bad_total)) {
    return(.pclv_failure(
      "spline_smoothing", "invalid_smoothed_composition_total",
      list(
        samples = names(column_totals)[bad_total],
        totals = unname(column_totals[bad_total])
      )
    ))
  }

  # Restore closure destroyed by independent taxon-wise smoothing.
  sm_mat <- sweep(sm_mat, 2L, column_totals, "/")
  closure_error <- abs(colSums(sm_mat) - 1)
  closure_tolerance <- max(1e-12, 64 * .Machine$double.eps * nrow(sm_mat))
  if (any(!is.finite(sm_mat)) || any(sm_mat < 0) ||
      any(!is.finite(closure_error)) ||
      max(closure_error) > closure_tolerance) {
    return(.pclv_failure(
      "spline_smoothing", "smoothed_composition_closure_failed",
      list(
        max_closure_error = suppressWarnings(max(closure_error)),
        tolerance = closure_tolerance
      )
    ))
  }

  # Subsetting happens only after normalization over the complete community.
  sm_mat[taxa_list, , drop = FALSE]
}

#'
#' Selects a subset of typical parameters if available (e.g., \code{a_ij},
#' \code{a_ii}, \code{r0}, noise and OU parameters) from a CmdStanR fit.
#'
#' @param fit A \pkg{cmdstanr} \code{CmdStanMCMC} fit.
#' @return A \pkg{posterior} draws data frame with available variables (possibly zero columns).
#' @noRd
#' @keywords internal
.safe_draws_df <- function(fit){
  drw <- fit$draws()
  keep <- intersect(
    c("a_ij","a_ii","r0",
      "sigma","sd_ou","phi",
      "sigma_ou","lambda","sigma_pred","tau_r","nu","log_nu_minus_two"),
    posterior::variables(drw)
  )
  if (!length(keep)) {
    return(posterior::as_draws_df(drw)[, 0, drop = FALSE])
  }
  posterior::as_draws_df(
    posterior::subset_draws(drw, variable = keep)
  )
}

#' Heuristic E-BFMI warning check via \code{cmdstan_diagnose()}
#'
#' Parses the output of \code{fit$cmdstan_diagnose()} and returns
#' \code{TRUE} if any chain’s E-BFMI appears below \code{thr}.
#'
#' @param fit A \pkg{cmdstanr} fit.
#' @param thr E-BFMI threshold (default 0.3).
#' @return Logical flag.
#' @noRd
#' @keywords internal
.ebfmi_warn_from_fit <- function(fit, thr = 0.3) {
  txt <- try(capture.output(fit$cmdstan_diagnose()), silent = TRUE)
  if (inherits(txt, "try-error") || is.null(txt)) return(FALSE)
  any(grepl(sprintf("E-BFMI .* less than %.1f", thr), txt))
}

#' Compute per-chain E-BFMI from sampler \code{energy__}
#'
#' Uses \code{mean(diff(E)^2)/var(E)} for each chain’s energy series.
#'
#' @param fit A \pkg{cmdstanr} \code{CmdStanMCMC} fit.
#' @return Numeric vector of E-BFMI per chain (may be empty).
#' @noRd
#' @keywords internal
.ebfmi_chainwise_from_energy <- function(fit) {
  sdiag <- try(fit$sampler_diagnostics(), silent = TRUE)
  if (inherits(sdiag, "try-error") || is.null(sdiag)) return(rep(NA_real_, 0))
  vars <- dimnames(sdiag)[[3]]
  pos  <- match("energy__", vars, nomatch = 0L)
  if (pos == 0L) return(rep(NA_real_, 0))
  E <- sdiag[, , pos, drop = FALSE]  # draws x chains x 1
  C <- dim(E)[2]
  eb <- rep(NA_real_, C)
  for (c in seq_len(C)) {
    ec <- as.numeric(E[, c, 1]); ec <- ec[is.finite(ec)]
    if (length(ec) < 3) next
    v <- stats::var(ec); if (!is.finite(v) || v <= 0) next
    d <- diff(ec)
    eb[c] <- mean(d * d) / v
  }
  eb
}

#' Summarise key MCMC diagnostics from a CmdStanR fit
#'
#' Reports worst \code{R-hat}, minimum \code{ESS} (bulk/tail), counts of
#' divergences and treedepth hits, per-chain E-BFMI summary, and total draws.
#'
#' @param fit A \pkg{cmdstanr} fit.
#' @param pars Parameter names to summarise if present.
#' @param max_treedepth Treedepth cap used by the sampler.
#' @return A named list of diagnostic summaries.
#' @noRd
#' @keywords internal
.add_convergence_diag <- function(diag, draws,
                                  pars = c("a_ij","a_ii","r0","sigma","sd_ou","phi","nu","log_nu_minus_two")) {
  avail <- posterior::variables(draws)
  use_pars <- intersect(pars, avail)
  if (!length(use_pars)) use_pars <- avail
  sdtab <- posterior::summarise_draws(posterior::subset_draws(draws, variable = use_pars))
  worst_rhat   <- suppressWarnings(max(sdtab$rhat,      na.rm = TRUE))
  min_ess_bulk <- suppressWarnings(min(sdtab$ess_bulk,  na.rm = TRUE))
  min_ess_tail <- suppressWarnings(min(sdtab$ess_tail,  na.rm = TRUE))
  if (!is.finite(worst_rhat))   worst_rhat   <- NA_real_
  if (!is.finite(min_ess_bulk)) min_ess_bulk <- NA_real_
  if (!is.finite(min_ess_tail)) min_ess_tail <- NA_real_

  diag$worst_rhat <- worst_rhat
  diag$min_ess_bulk <- min_ess_bulk
  diag$min_ess_tail <- min_ess_tail
  diag
}

.summarise_sampler_diag <- function(fit, max_treedepth = 12) {
  sdiag <- fit$sampler_diagnostics()
  sdiag_df <- posterior::as_draws_df(sdiag)

  n_div <- if ("divergent__" %in% names(sdiag_df)) {
    sum(as.integer(sdiag_df[["divergent__"]]), na.rm = TRUE)
  } else NA_integer_

  n_treedepth_hit <- if ("treedepth__" %in% names(sdiag_df)) {
    sum(as.integer(sdiag_df[["treedepth__"]] >= max_treedepth), na.rm = TRUE)
  } else NA_integer_

  # E-BFMI: mean(diff(E)^2) / var(E), once per chain and attempt.
  ebfmi_min <- ebfmi_med <- NA_real_
  eb <- numeric()
  if ("energy__" %in% names(sdiag_df)) {
    if (".chain" %in% names(sdiag_df)) {
      eb <- vapply(split(sdiag_df[["energy__"]], sdiag_df[[".chain"]]), function(ev) {
        e <- as.numeric(ev); e <- e[is.finite(e)]
        if (length(e) < 3) return(NA_real_)
        v <- stats::var(e); if (!is.finite(v) || v <= 0) return(NA_real_)
        mean(diff(e)^2) / v
      }, numeric(1))
      ebfmi_min <- suppressWarnings(min(eb, na.rm = TRUE)); if (!is.finite(ebfmi_min)) ebfmi_min <- NA_real_
      ebfmi_med <- suppressWarnings(stats::median(eb, na.rm = TRUE)); if (!is.finite(ebfmi_med)) ebfmi_med <- NA_real_
    } else {
      e <- as.numeric(sdiag_df[["energy__"]]); e <- e[is.finite(e)]
      if (length(e) >= 3) {
        v <- stats::var(e)
        if (is.finite(v) && v > 0) {
          val <- mean(diff(e)^2) / v
          eb <- val; ebfmi_min <- val; ebfmi_med <- val
        }
      }
    }
  }

  n_draws <- if (length(dim(sdiag)) >= 2L) prod(dim(sdiag)[1:2]) else nrow(sdiag_df)

  list(
    worst_rhat      = NA_real_,
    min_ess_bulk    = NA_real_,
    min_ess_tail    = NA_real_,
    n_divergent     = n_div,
    n_treedepth_hit = n_treedepth_hit,
    ebfmi_min       = ebfmi_min,
    ebfmi_med       = ebfmi_med,
    ebfmi_chain     = eb,
    n_draws         = as.integer(n_draws)
  )
}

# Backward-compatible complete diagnostic helper. Canonical retry execution
# uses .summarise_sampler_diag(); only a retained fit pays for convergence
# summaries after its selected draws have been materialized once.
.summarise_diag <- function(fit,
                            pars = c("a_ij","a_ii","r0","sigma","sd_ou","phi","nu","log_nu_minus_two"),
                            max_treedepth = 12) {
  draws <- .safe_draws_df(fit)
  .add_convergence_diag(.summarise_sampler_diag(fit, max_treedepth), draws, pars)
}

# A chain must put at least this much posterior mass on one sign before that
# sign is treated as identified. This matches the package's existing 95%
# posterior-sign interpretation while keeping LFSR filtering separate.
.PCLV_CHAIN_SIGN_PROB_MIN <- 0.95
.PCLV_EBFMI_MIN <- 0.30
`%||%` <- function(x, fallback) if (is.null(x)) fallback else x

.chain_parameter_summary <- function(draws, parameter) {
  if (!all(c(parameter, ".chain") %in% names(draws))) return(data.frame())
  pieces <- split(as.numeric(draws[[parameter]]), draws$.chain)
  rows <- lapply(names(pieces), function(chain) {
    x <- pieces[[chain]]
    x <- x[is.finite(x)]
    if (!length(x)) return(NULL)
    data.frame(
      chain = as.integer(chain), parameter = parameter,
      mean = mean(x), median = stats::median(x),
      q05 = unname(stats::quantile(x, 0.05)), q95 = unname(stats::quantile(x, 0.95)),
      positive_probability = mean(x > 0), negative_probability = mean(x < 0),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

.chain_intervals_overlap <- function(x) {
  nrow(x) > 0L && all(is.finite(x$q05)) && all(is.finite(x$q95)) &&
    max(x$q05) <= min(x$q95)
}

.classify_chain_diagnostics <- function(draws, diag) {
  empty <- list(
    diagnostic_class = "sampler_diagnostics_failed",
    interaction_identifiable = FALSE, residual_identifiable = FALSE,
    chain_sign_agreement = NA, pooled_sign_probability = NA_real_,
    chain_aij_means = numeric(), chain_aij_medians = numeric(),
    chain_aij_positive_probabilities = numeric(), chain_aij_negative_probabilities = numeric(),
    chain_aij_mean_range = NA_real_, chain_aij_median_range = NA_real_,
    chain_sign_probability_max_difference = NA_real_,
    chain_aii_means = numeric(), chain_aii_medians = numeric(),
    chain_aii_sign_probabilities = numeric(), chain_aii_sign_agreement = NA,
    chain_residual_summary = data.frame(), residual_median_ranges = numeric(),
    residual_regime_disagreement = NA,
    indeterminate_reason = "missing_or_invalid_posterior_draws"
  )
  if (!is.data.frame(draws) || !all(c("a_ij", "a_ii", ".chain") %in% names(draws)) ||
      any(!is.finite(draws$a_ij)) || any(!is.finite(draws$a_ii))) return(empty)

  aij <- .chain_parameter_summary(draws, "a_ij")
  aii <- .chain_parameter_summary(draws, "a_ii")
  if (!nrow(aij) || !nrow(aii)) return(empty)
  n_chains <- nrow(aij)
  dominant <- ifelse(aij$median > 0, 1L, ifelse(aij$median < 0, -1L, 0L))
  certain <- pmax(aij$positive_probability, aij$negative_probability) >= .PCLV_CHAIN_SIGN_PROB_MIN
  sign_agree <- if (n_chains > 1L) all(certain) && length(unique(dominant)) == 1L && dominant[[1L]] != 0L else NA
  magnitude_agree <- if (n_chains > 1L) .chain_intervals_overlap(aij) else NA
  pooled_pos <- mean(draws$a_ij > 0)
  pooled_neg <- mean(draws$a_ij < 0)

  self_dominant <- ifelse(aii$median > 0, 1L, ifelse(aii$median < 0, -1L, 0L))
  self_certain <- pmax(aii$positive_probability, aii$negative_probability) >= .PCLV_CHAIN_SIGN_PROB_MIN
  self_agree <- if (n_chains > 1L) all(self_certain) && length(unique(self_dominant)) == 1L && self_dominant[[1L]] != 0L else NA

  residual_parameters <- intersect(c("sigma", "sd_ou", "phi", "lambda", "nu"), names(draws))
  residual <- dplyr::bind_rows(lapply(residual_parameters, function(p) .chain_parameter_summary(draws, p)))
  residual_ranges <- numeric()
  separated <- character()
  if (nrow(residual)) {
    by_parameter <- split(residual, residual$parameter)
    residual_ranges <- vapply(by_parameter, function(x) diff(range(x$median)), numeric(1))
    if (n_chains > 1L) separated <- names(Filter(function(x) !.chain_intervals_overlap(x), by_parameter))
  }
  residual_disagree <- if (n_chains > 1L) length(separated) > 0L else NA

  diag_ok <- is.list(diag) && isTRUE(.ok_diag(
    diag$worst_rhat, diag$min_ess_bulk, diag$min_ess_tail,
    diag$n_divergent, diag$n_treedepth_hit, diag$n_draws
  )) && is.finite(diag$ebfmi_min) && diag$ebfmi_min >= .PCLV_EBFMI_MIN
  diagnostic_reasons <- character()
  if (!is.list(diag) || !is.finite(diag$worst_rhat) || diag$worst_rhat >= 1.05) diagnostic_reasons <- c(diagnostic_reasons, "high_rhat")
  if (!is.list(diag) || !is.finite(diag$min_ess_bulk) || !is.finite(diag$min_ess_tail) || diag$min_ess_bulk <= 400 || diag$min_ess_tail <= 400) diagnostic_reasons <- c(diagnostic_reasons, "low_ess")
  if (!is.list(diag) || !is.finite(diag$ebfmi_min) || diag$ebfmi_min < .PCLV_EBFMI_MIN) diagnostic_reasons <- c(diagnostic_reasons, "low_ebfmi")
  if (is.list(diag) && is.finite(diag$n_divergent) && diag$n_divergent > 0) diagnostic_reasons <- c(diagnostic_reasons, "divergences")
  if (is.list(diag) && is.finite(diag$n_treedepth_hit) && diag$n_treedepth_hit > 0) diagnostic_reasons <- c(diagnostic_reasons, "treedepth_saturation")

  interaction_reasons <- character()
  if (n_chains > 1L && (length(unique(dominant)) != 1L || dominant[[1L]] == 0L)) {
    interaction_reasons <- "chain_sign_disagreement"
  } else if (n_chains > 1L && (!all(certain) || !isTRUE(magnitude_agree))) {
    interaction_reasons <- "interaction_magnitude_disagreement"
  }
  residual_reasons <- character()
  if (length(separated)) {
    residual_reasons <- "chain_specific_residual_regime"
    if (any(c("sigma", "sd_ou") %in% separated)) residual_reasons <- c(residual_reasons, "residual_scale_nonidentifiability")
    if ("sd_ou" %in% separated && any(c("phi", "lambda") %in% separated)) residual_reasons <- c(residual_reasons, "ou_scale_decay_ridge")
  }

  if (n_chains == 1L) {
    diagnostic_class <- if (diag_ok) "converged" else "sampler_diagnostics_failed"
  } else if (length(interaction_reasons)) {
    diagnostic_class <- "interaction_indeterminate"
  } else if (!diag_ok || isTRUE(residual_disagree)) {
    diagnostic_class <- "interaction_stable_residual_unstable"
  } else diagnostic_class <- "converged"
  reasons <- unique(c(interaction_reasons, residual_reasons, diagnostic_reasons))

  list(
    diagnostic_class = diagnostic_class,
    interaction_identifiable = diagnostic_class %in% c("converged", "interaction_stable_residual_unstable"),
    residual_identifiable = identical(diagnostic_class, "converged"),
    chain_sign_agreement = sign_agree, pooled_sign_probability = max(pooled_pos, pooled_neg),
    chain_aij_means = stats::setNames(aij$mean, aij$chain), chain_aij_medians = stats::setNames(aij$median, aij$chain),
    chain_aij_positive_probabilities = stats::setNames(aij$positive_probability, aij$chain),
    chain_aij_negative_probabilities = stats::setNames(aij$negative_probability, aij$chain),
    chain_aij_mean_range = diff(range(aij$mean)), chain_aij_median_range = diff(range(aij$median)),
    chain_sign_probability_max_difference = diff(range(aij$positive_probability)),
    chain_aii_means = stats::setNames(aii$mean, aii$chain), chain_aii_medians = stats::setNames(aii$median, aii$chain),
    chain_aii_sign_probabilities = stats::setNames(pmax(aii$positive_probability, aii$negative_probability), aii$chain),
    chain_aii_sign_agreement = self_agree, chain_residual_summary = residual,
    residual_median_ranges = residual_ranges, residual_regime_disagreement = residual_disagree,
    indeterminate_reason = reasons
  )
}

.build_posterior_summary_bundle <- function(draws, diag) {
  if (!all(c("a_ij", "a_ii") %in% names(draws)))
    stop("Expected parameters a_ij and a_ii not found in draws.")
  nu <- .summarise_nu_draws(draws)
  if (.is_pclv_failure(nu)) return(nu)
  aij <- draws$a_ij
  aii <- draws$a_ii
  list(
    coefficients = list(
      interaction = c(mean = mean(aij), median = stats::median(aij), sd = stats::sd(aij),
                      q025 = unname(stats::quantile(aij, 0.025)),
                      q975 = unname(stats::quantile(aij, 0.975))),
      self = c(mean = mean(aii), sd = stats::sd(aii),
               q025 = unname(stats::quantile(aii, 0.025)),
               q975 = unname(stats::quantile(aii, 0.975)))
    ),
    sign = {
      interaction_positive_probability <- mean(aij > 0)
      interaction_negative_probability <- mean(aij < 0)
      interaction_p_two <- .clip01(2 * pmin(interaction_positive_probability, interaction_negative_probability))
      self_p_two <- .clip01(2 * pmin(mean(aii > 0), mean(aii < 0)))
      list(
        interaction_positive_probability = interaction_positive_probability,
        interaction_negative_probability = interaction_negative_probability,
        interaction_p_two = interaction_p_two,
        interaction_lfsr = .lfsr_from_two_sided(interaction_p_two),
        self_p_two = self_p_two,
        self_lfsr = .lfsr_from_two_sided(self_p_two)
      )
    },
    nu = nu,
    chain = .classify_chain_diagnostics(draws, diag)
  )
}

#' Local false sign rate from two-sided tail probability
#'
#' @param p_two Two-sided tail probability in \code{[0,1]}.
#' @return LFSR in \code{[0,0.5]}.
#' @noRd
#' @keywords internal
.lfsr_from_two_sided <- function(p_two) pmax(pmin(p_two / 2, 0.5), 0)

#' Clamp numeric vector to \code{[0,1]}
#' @param x Numeric vector.
#' @return Numeric vector clamped to \code{[0,1]}.
#' @noRd
#' @keywords internal
.clip01 <- function(x) pmin(pmax(x, 0), 1)

#' Alias of \code{.lfsr_from_two_sided()}
#' @inheritParams .lfsr_from_two_sided
#' @return LFSR in \code{[0,0.5]}.
#' @noRd
#' @keywords internal
.lfsr_from_two <- function(p_two) pmin(pmax(p_two/2, 0), 0.5)

#' Monotone cumulative q-value from LFSR
#'
#' Computes a cumulative average of sorted LFSR values, producing a
#' conservative, monotone \eqn{q}-like measure for sign error control.
#'
#' @param v Numeric vector of LFSR values.
#' @return Numeric vector of same length with cumulative \eqn{q}.
#' @noRd
#' @keywords internal
.q_from_lfsr <- function(v) {
  q <- rep(NA_real_, length(v))
  nn <- which(!is.na(v))
  if (length(nn)) {
    l  <- pmin(pmax(v[nn], 0), 0.5)
    oo <- nn[order(l)]
    q[oo] <- cumsum(l[order(l)]) / seq_along(oo)
  }
  q
}

#' Element-wise diagnostic pass/fail predicate
#'
#' Applies thresholds to vectors of diagnostics: \code{R-hat < 1.05},
#' \code{ESS > 400}, and small divergence/treedepth rates scaled by draws.
#'
#' @param rhat Worst R-hat.
#' @param essb Bulk ESS.
#' @param esst Tail ESS.
#' @param div Divergence counts.
#' @param tdhit Treedepth-hit counts.
#' @param n_draws Total draws (to scale tolerances).
#' @return Logical vector indicating OK diagnostics.
#' @noRd
#' @keywords internal
.ok_diag <- function(rhat, essb, esst, div, tdhit, n_draws = NA_real_) {
  # 길이 통일
  L <- max(length(rhat), length(essb), length(esst), length(div), length(tdhit), length(n_draws))
  rhat    <- rep_len(rhat,    L)
  essb    <- rep_len(essb,    L)
  esst    <- rep_len(esst,    L)
  div     <- rep_len(div,     L)
  tdhit   <- rep_len(tdhit,   L)
  n_draws <- rep_len(n_draws, L)

  # 임계값
  rhat_thr <- 1.05
  ess_thr  <- 400

  # n_draws가 스칼라가 아니어도 element-wise로 처리
  nd_ok      <- is.finite(n_draws) & n_draws > 0
  div_max    <- ifelse(nd_ok, ceiling(0.001 * n_draws), 8L)    # ~0.1%
  tdepth_max <- ifelse(nd_ok, ceiling(0.010 * n_draws), 80L)   # ~1%

  (is.finite(rhat)  & rhat < rhat_thr) &
    (is.finite(essb)  & essb > ess_thr)  &
    (is.finite(esst)  & esst > ess_thr)  &
    (is.finite(div)   & div  <= div_max) &
    (is.finite(tdhit) & tdhit<= tdepth_max)
}

# ----- Zero-aware ALR building & transforms -----
#' Construct zero-aware ALR triplet for a single row (explicit args)
#'
#' Same logic as before, but all required state is passed as parameters to avoid
#' free-variable lookups that can fail in parallel workers.
#'
#' @param i Row index.
#' @param df Data frame containing xi_raw, xj_raw, rest_raw.
#' @param subj_minpos Named/parallel vector of subject-level min positives.
#' @param lib Library-size vector (may be NA).
#' @param zero_mode_alr, minpos_alpha, eps_fixed, lib_eps_c, rest_floor_frac, minpos_base Settings.
#' @return Numeric vector c(xi, xj, xr, eps_star).
#' @noRd
#' @keywords internal
.make_triplet_row <- function(i, df, subj_minpos, lib,
                              zero_mode_alr, minpos_alpha, eps_fixed,
                              lib_eps_c, rest_floor_frac, minpos_base) {

  # 항상 numeric(1)로 강제 (NULL/list도 안전하게 NA로)
  xi <- as.numeric(df$xi_raw[i])[1]
  xj <- as.numeric(df$xj_raw[i])[1]
  xr <- as.numeric(df$rest_raw[i])[1]

  eps_t <- switch(zero_mode_alr,
                  "minpos_time" = {
                    base_pos <- c(
                      if (is.finite(xi) && xi > 0) xi,
                      if (is.finite(xj) && xj > 0) xj
                    )
                    if (identical(minpos_base, "triplet") &&
                        is.finite(xr) && xr > 0) {
                      base_pos <- c(base_pos, xr)
                    }
                    if (length(base_pos)) minpos_alpha * min(base_pos) else eps_fixed
                  },
                  "minpos_subject" = {
                    mp <- as.numeric(subj_minpos[i])[1]
                    if (is.finite(mp) && mp > 0) minpos_alpha * mp else eps_fixed
                  },
                  "lib" = {
                    L <- as.numeric(lib[i])[1]
                    if (!is.finite(L) || L <= 0) L <- 1
                    max(eps_fixed, lib_eps_c / L)
                  },
                  "fixed" = eps_fixed)

  xi <- if (is.finite(xi) && xi > 0) xi else eps_t
  xj <- if (is.finite(xj) && xj > 0) xj else eps_t
  xr <- if (is.finite(xr) && xr > 0) xr else eps_t
  xr <- max(xr, rest_floor_frac * eps_t)

  s <- xi + xj + xr
  if (!is.finite(s) || s <= 0) s <- 1

  xi_ <- xi / s; xj_ <- xj / s; xr_ <- xr / s
  eps_star <- eps_t / s

  c(as.numeric(xi_)[1], as.numeric(xj_)[1], as.numeric(xr_)[1], as.numeric(eps_star)[1])
}

.validate_subject_times <- function(time, subject, stage = "preprocessing") {
  if (any(!is.finite(time))) return(.pclv_failure(stage, "non_finite_time", list()))
  by_subject <- split(seq_along(time), subject)
  bad <- any(vapply(by_subject, function(ix) {
    length(ix) > 1L && any(diff(sort(time[ix])) <= 0)
  }, logical(1)))
  if (bad) return(.pclv_failure(stage, "non_positive_dt", list()))
  NULL
}

.select_smoothed_pair_rows <- function(sm_mat, meta_df, j, i,
                                       min_pairs = 4, min_dt = 1e-8, min_sd = 1e-12) {
  stopifnot(all(c("Sample", "subject", "time") %in% names(meta_df)))
  if (!all(c(i, j) %in% rownames(sm_mat))) {
    stop("Validated pair taxon invariant violated.")
  }

  metadata_samples <- as.character(meta_df$Sample)
  matrix_samples <- colnames(sm_mat)
  sample_idx <- match(metadata_samples, matrix_samples)
  if (length(sample_idx) != nrow(meta_df) || anyNA(sample_idx) ||
      anyDuplicated(sample_idx) || anyDuplicated(metadata_samples)) {
    return(.pclv_failure(
      "preprocessing", "sample_alignment_failed", list()
    ))
  }

  pair_df <- data.frame(
    subject = meta_df$subject,
    time = meta_df$time,
    xi_raw = pmax(as.numeric(sm_mat[i, sample_idx]), 0),
    xj_raw = pmax(as.numeric(sm_mat[j, sample_idx]), 0)
  )
  pair_df <- pair_df[order(pair_df$subject, pair_df$time), , drop = FALSE]

  if (any(!is.finite(pair_df$xi_raw)) || any(!is.finite(pair_df$xj_raw))) {
    return(.pclv_failure(
      "preprocessing", "non_finite_pair_abundance", list()
    ))
  }

  time_failure <- .validate_subject_times(pair_df$time, pair_df$subject)
  if (!is.null(time_failure)) return(time_failure)

  by_subject <- split(seq_len(nrow(pair_df)), pair_df$subject)
  transition_counts <- vapply(by_subject, function(ix) {
    if (length(ix) < 2L) return(0L)
    dt <- diff(as.numeric(pair_df$time[ix]))
    if (any(!is.finite(dt)) || any(dt <= min_dt)) return(NA_integer_)
    length(dt)
  }, integer(1))
  if (anyNA(transition_counts)) {
    return(.pclv_failure(
      "preprocessing", "dt_below_minimum", list(min_dt = min_dt)
    ))
  }

  # Preserve every raw observation. Lagging below removes only the first row per
  # subject. Keeping only predecessor rows here would systematically discard the
  # final valid transition of every subject.
  n_transitions <- sum(transition_counts)
  if (n_transitions < min_pairs) {
    return(.pclv_failure(
      "preprocessing", "insufficient_rows",
      list(observed_rows = n_transitions, required_rows = min_pairs)
    ))
  }

  sd_i <- stats::sd(pair_df$xi_raw)
  sd_j <- stats::sd(pair_df$xj_raw)
  if (!is.finite(sd_i) || !is.finite(sd_j) || sd_i < min_sd || sd_j < min_sd) {
    return(.pclv_failure(
      "preprocessing", "insufficient_raw_variation",
      list(observed_sd_i = sd_i, observed_sd_j = sd_j, required_sd = min_sd)
    ))
  }

  pair_df
}

.build_pair_df_smoothed <- function(sm_mat, meta_df, j, i, eps = 1e-8,
                                    min_pairs = 4, min_dt = 1e-8, min_sd = 1e-12) {
  .select_smoothed_pair_rows(sm_mat, meta_df, j, i, min_pairs, min_dt, min_sd)
}

.pair_to_rest_abundance <- function(xi, xj) pmax(0, 1 - xi - xj)

.pair_alr <- function(trip, alr_cap) {
  alr_i <- log(trip[, "xi"]) - log(trip[, "xr"])
  alr_j <- log(trip[, "xj"]) - log(trip[, "xr"])
  if (is.finite(alr_cap)) {
    alr_i <- pmax(pmin(alr_i, alr_cap), -alr_cap)
    alr_j <- pmax(pmin(alr_j, alr_cap), -alr_cap)
  } else {
    es <- pmax(trip[, "eps_star"], .Machine$double.eps)
    cap_vec <- pmin(8, pmax(0, log(pmax(1 - 2 * es, .Machine$double.eps) / es)))
    alr_i <- pmax(pmin(alr_i, cap_vec), -cap_vec)
    alr_j <- pmax(pmin(alr_j, cap_vec), -cap_vec)
  }
  list(i = alr_i, j = alr_j)
}

.delta_alr_over_dt <- function(v, time, subject) {
  time_failure <- .validate_subject_times(time, subject)
  if (!is.null(time_failure)) return(time_failure)
  by_subject <- split(seq_along(v), subject)
  unsplit(lapply(by_subject, function(ix) {
    vi <- v[ix]; ti <- time[ix]
    vi_lag <- dplyr::lag(vi, .PCLV_CORE_LAG)
    dt <- as.numeric(ti - dplyr::lag(ti, .PCLV_CORE_LAG))
    out <- (vi - vi_lag) / dt
    out[!is.finite(dt) | dt <= 0] <- NA_real_
    out
  }), subject)
}

#' Build model inputs for pairwise gLV regressions
#'
#' Creates lagged predictors and \code{ΔALR_i/Δt} response under either
#' ALR transformation with zero-aware ALR safeguards, optional
#' ALR smoothing, partner non-zero filters, and global predictor scaling.
#'
#' @param pair_df Optional preselected raw pair-abundance table; when NULL, rows are selected from sm_mat and meta_df.
#' @param lag Positive integer lag for predictors within subject.
#' @param zero_mode_alr Zero-handling mode for ALR (see code for options).
#' @param minpos_alpha Multiplier for data-driven epsilon.
#' @param minpos_base Whether min-positive search uses \code{"ij"} or \code{"triplet"}.
#' @param eps_fixed Fixed epsilon used when needed.
#' @param lib_eps_c Library-size epsilon coefficient (when \code{zero_mode_alr="lib"}).
#' @param rest_floor_frac Lower bound for \code{rest} after replacement.
#' @param alr_cap Finite cap for ALR magnitudes; \code{Inf} enables theory-based cap.
#' @param smooth_scale One of \code{"logra"} (pre-smoothed) or \code{"alr"} (inline).
#' @param alr_spline_df,alr_spline_spar,alr_spline_cv Spline controls for ALR smoothing.
#' @param nz_partner_min_frac Minimum fraction of non-zero partner entries per subject.
#' @return Data frame with columns \code{subject,time,y,xi,xj} and attributes
#'   \code{smooth_edf_mean}, \code{smooth_scale}, and \code{smoothed}.
#' @noRd
#' @keywords internal
.make_pair_inputs_glv <- function(pair_df = NULL,
                                  sm_mat = NULL, meta_df = NULL,
                                  j = NULL, i = NULL, min_pairs = 4,
                                  min_dt = 1e-8, min_sd = 1e-12,
                                  zero_mode_alr = c("minpos_time","minpos_subject","lib","fixed"),
                                  minpos_alpha = 0.5,
                                  minpos_base = c("ij","triplet"),
                                  eps_fixed = 1e-8,
                                  lib_eps_c = 0.65,
                                  rest_floor_frac = 1.0,
                                  alr_cap = Inf,
                                  smooth_scale = c("logra","alr"),
                                  alr_spline_df = NULL,
                                  alr_spline_spar = NULL,
                                  alr_spline_cv = TRUE,
                                  nz_partner_min_frac = 0.15) {
  zero_mode_alr <- zero_mode_alr[[1L]]
  minpos_base <- minpos_base[[1L]]
  smooth_scale <- smooth_scale[[1L]]
  if (!(zero_mode_alr %in% c("minpos_time", "minpos_subject", "lib", "fixed")) ||
      !(minpos_base %in% c("ij", "triplet")) ||
      !(smooth_scale %in% c("logra", "alr")))
    stop("Validated canonical preprocessing option invariant violated.")

  if (is.null(pair_df)) {
    pair_df <- .select_smoothed_pair_rows(
      sm_mat, meta_df, j, i, min_pairs, min_dt, min_sd
    )
    if (.is_pclv_failure(pair_df)) return(pair_df)
  }

  stopifnot(all(c("subject","time","xi_raw","xj_raw") %in% names(pair_df)))
  df <- pair_df[order(pair_df$subject, pair_df$time), , drop = FALSE]
  time_failure <- .validate_subject_times(df$time, df$subject)
  if (!is.null(time_failure)) return(time_failure)

  df$xi_raw <- pmax(as.numeric(df$xi_raw), 0)
  df$xj_raw <- pmax(as.numeric(df$xj_raw), 0)
  pair_total <- df$xi_raw + df$xj_raw
  composition_tolerance <- max(1e-12, 64 * .Machine$double.eps)
  if (any(!is.finite(pair_total)) ||
      any(pair_total > 1 + composition_tolerance)) {
    return(.pclv_failure(
      "preprocessing", "invalid_pair_composition",
      list(max_pair_total = suppressWarnings(max(pair_total)))
    ))
  }
  df$rest_raw <- .pair_to_rest_abundance(df$xi_raw, df$xj_raw)

  by_s <- split(seq_len(nrow(df)), df$subject)
  lagv <- function(v) unsplit(lapply(by_s, function(ix) dplyr::lag(v[ix], .PCLV_CORE_LAG)), df$subject)

  lib <- if ("libsize" %in% names(df)) df$libsize else NA_real_

  subj_minpos <- unsplit(lapply(by_s, function(ix) {
    base <- if (minpos_base == "ij") {
      c(df$xi_raw[ix][df$xi_raw[ix] > 0], df$xj_raw[ix][df$xj_raw[ix] > 0])
    } else {
      c(df$xi_raw[ix][df$xi_raw[ix] > 0],
        df$xj_raw[ix][df$xj_raw[ix] > 0],
        df$rest_raw[ix][df$rest_raw[ix] > 0])
    }
    if (length(base)) min(base) else NA_real_
  }), df$subject)

  trip <- t(vapply(
    seq_len(nrow(df)),
    function(ii)
      .make_triplet_row(ii, df, subj_minpos, lib,
                        zero_mode_alr, minpos_alpha, eps_fixed,
                        lib_eps_c, rest_floor_frac, minpos_base),
    numeric(4)
  ))
  colnames(trip) <- c("xi","xj","xr","eps_star")
  alr <- .pair_alr(trip, alr_cap)
  alr_i <- alr$i
  alr_j <- alr$j

  smooth_edf_mean <- NA_real_
  if (smooth_scale == "alr") {
    smooth_one <- function(v, t) {
      ok <- is.finite(v) & is.finite(t)
      if (sum(ok) < 3L || length(unique(t[ok])) < 3L) {
        return(.pclv_failure("spline_smoothing", "insufficient_unique_times", list()))
      }
      rr <- .smooth_spline_robust(
        x = t[ok], y = v[ok],
        spline_df = alr_spline_df,
        spline_spar = alr_spline_spar,
        use_cv = isTRUE(alr_spline_cv),
        min_unique = 3,
        min_df = 3.0,
        max_df = NULL,
        default_df = 4.0
      )
      if (.is_pclv_failure(rr)) return(rr)
      yhat <- rep(NA_real_, length(t))
      yhat[ok] <- rr$yhat
      list(y = yhat, df = rr$df)
    }
    edf_i <- edf_j <- rep(NA_real_, length(by_s))
    names(edf_i) <- names(edf_j) <- names(by_s)
    for (sname in names(by_s)) {
      ix <- by_s[[sname]]
      res_i <- smooth_one(alr_i[ix], df$time[ix])
      res_j <- smooth_one(alr_j[ix], df$time[ix])
      if (.is_pclv_failure(res_i)) return(res_i)
      if (.is_pclv_failure(res_j)) return(res_j)
      alr_i[ix] <- res_i$y
      alr_j[ix] <- res_j$y
      edf_i[sname] <- res_i$df
      edf_j[sname] <- res_j$df
    }
    edf_subj <- rowMeans(cbind(edf_i, edf_j), na.rm = TRUE)
    smooth_edf_mean <- mean(edf_subj, na.rm = TRUE)
    if (!is.finite(smooth_edf_mean)) smooth_edf_mean <- NA_real_
  }

  xi <- lagv(alr_i)
  xj <- lagv(alr_j)
  y <- .delta_alr_over_dt(alr_i, df$time, df$subject)
  if (.is_pclv_failure(y)) return(y)

  # 파트너 희소성 필터(옵션)
  if (is.finite(nz_partner_min_frac) && nz_partner_min_frac > 0) {
    by_s2 <- split(seq_len(nrow(df)), df$subject)
    keep_mask <- rep(TRUE, nrow(df))
    for (sname in names(by_s2)) {
      ix <- by_s2[[sname]]
      # Sparsity is a property of the raw partner abundance, not of its ALR.
      # A zero-replaced ALR is almost never exactly zero and therefore cannot be
      # used to detect an absent partner.
      frac_nz <- mean(is.finite(df$xj_raw[ix]) & df$xj_raw[ix] > 0)
      if (is.finite(frac_nz) && frac_nz < nz_partner_min_frac) {
        keep_mask[ix] <- FALSE
      }
    }
  } else {
    keep_mask <- rep(TRUE, nrow(df))
  }

  dat <- data.frame(
    subject = df$subject,
    time    = df$time,
    y = as.numeric(y),
    xi = as.numeric(xi),
    xj = as.numeric(xj)
  )

  ok <- is.finite(dat$y) & is.finite(dat$xi) & is.finite(dat$xj) & is.finite(dat$time) & keep_mask
  dat <- dat[ok, , drop = FALSE]
  if (!nrow(dat)) return(.pclv_failure("preprocessing", "no_valid_lagged_rows", list()))

  variation_failure <- .validate_predictor_variation(dat$xi, dat$xj, "full_data")
  if (!is.null(variation_failure)) {
    return(variation_failure)
  }


  mu_xi <- mean(dat$xi)
  sd_xi <- stats::sd(dat$xi)
  if (!is.finite(sd_xi) || sd_xi <= 0)
    return(.pclv_failure("standardization", "invalid_internal_scaling_state",
                         list(predictor = "xi", observed_sd = sd_xi)))
  mu_xj <- mean(dat$xj)
  sd_xj <- stats::sd(dat$xj)
  if (!is.finite(sd_xj) || sd_xj <= 0)
    return(.pclv_failure("standardization", "invalid_internal_scaling_state",
                         list(predictor = "xj", observed_sd = sd_xj)))

  dat$xi_unscaled <- dat$xi
  dat$xj_unscaled <- dat$xj
  dat$xi <- (dat$xi - mu_xi) / sd_xi
  dat$xj <- (dat$xj - mu_xj) / sd_xj

  attr(dat, "smooth_edf_mean") <- smooth_edf_mean
  attr(dat, "smooth_scale")    <- smooth_scale
  attr(dat, "smoothed")        <- smooth_scale %in% c("alr","logra")

  dat
}

# -------- cmdstanr call silencer & loglik extractor --------
#' Call \code{cmdstanr::sample()} with suppressed console output
#'
#' Silences stdout/messages (and sets \code{refresh=0}) to keep logs clean.
#'
#' @param mod A \pkg{cmdstanr} model.
#' @param args List of arguments forwarded to \code{mod$sample()}.
#' @param silent Logical; if \code{FALSE}, calls \code{sample()} verbatim.
#' @return A \pkg{cmdstanr} fit.
#' @noRd
#' @keywords internal
#' cmdstanr::sample()를 무음으로 호출 + 출력 경로/파일 충돌 방지
#'
#' - silent=TRUE: stdout/message 억제, refresh=0
#' - output_dir / output_basename 보장 및 basename에 UID 접미사 강제 부여
#' - 초미니 반복(iter_warmup/sampling <= 5)이면 짧은 랜덤 지터로 파일 I/O 경합 완화
#' - (옵션) args$.smoke_mode=TRUE면 parallel_chains <- 1 강제
#' - init 관련 오류가 나면 1회 폴백(init=NULL) + 새 디렉터리/베이스네임으로 재시도
.call_sample_silently <- function(mod, args, silent = TRUE) {
  # 1) 콘솔 무음 + 진행/메시지 억제
  if (isTRUE(silent)) {
    args$refresh <- 0L
    args$show_messages <- FALSE
  }

  # 유틸
  .make_uid <- function(n = 8L) paste(sample(c(letters, 0:9), n, TRUE), collapse = "")
  .now_tag  <- function() format(Sys.time(), "%Y%m%d%H%M%OS3")

  # 영속 루트: 옵션 없으면 사용자 캐시 디렉터리
  output_root <- getOption("glvpair.output_root",
                           tools::R_user_dir("glvpair", which = "cache"))
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  .ensure_outputs <- function(a, root = output_root,
                              dir_prefix = "run",
                              base_prefix = "glv_pairwise") {
    # output_dir 준비
    if (is.null(a$output_dir) || !nzchar(a$output_dir)) {
      subdir <- file.path(root, sprintf("%s_%s_%s", dir_prefix, .now_tag(), .make_uid(6)))
      dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      a$output_dir <- subdir
    } else {
      dir.create(a$output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    # output_basename: 항상 최종적으로 UID 접미사 부여(사용자 지정이어도 충돌 방지)
    if (is.null(a$output_basename) || !nzchar(a$output_basename)) {
      a$output_basename <- sprintf("%s", base_prefix)
    }
    a$output_basename <- sprintf("%s_%s", a$output_basename, .make_uid(6))

    # chain_ids 기본값
    if (is.null(a$chain_ids) && !is.null(a$chains)) {
      a$chain_ids <- seq_len(as.integer(a$chains))
    }
    a
  }

  owns_output_dir <- is.null(args$output_dir) || !nzchar(args$output_dir)
  args <- .ensure_outputs(args)
  owned_output_dir <- if (owns_output_dir) args$output_dir else NULL
  sample_completed <- FALSE
  on.exit({
    if (!sample_completed && !is.null(owned_output_dir))
      unlink(owned_output_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  # --- 스모크 모드(옵션): 초경량 테스트 시 체인 순차 실행 강제 ---
  if (isTRUE(args$.smoke_mode)) {
    args$parallel_chains <- 1L
  }

  # --- 초미니 반복에서 타임스탬프/파일 I/O 경합 완화용 짧은 지터 ---
  if (is.numeric(args$iter_warmup) && is.numeric(args$iter_sampling)) {
    if (args$iter_warmup <= 5L && args$iter_sampling <= 5L) {
      Sys.sleep(runif(1, 0, 0.25))  # 최대 0.25초 랜덤 대기
    }
  }

  # 3) stdout / message 무음 처리
  if (isTRUE(silent)) {
    out_con <- file(nullfile(), open = "wt")
    msg_con <- file(nullfile(), open = "wt")
    on.exit({
      try(sink(type = "message"), silent = TRUE)
      try(sink(), silent = TRUE)
      try(close(msg_con), silent = TRUE)
      try(close(out_con), silent = TRUE)
    }, add = TRUE)
    sink(out_con); sink(msg_con, type = "message")
  }

  # Sampling errors are handled by .sample_with_retry(); do not retry or alter init here.
  fit <- do.call(mod$sample, args)
  if (!is.null(owned_output_dir))
    attr(fit, "pclv_owned_output_dir") <- owned_output_dir
  sample_completed <- TRUE
  fit
}

# Remove only a sampler directory created and owned by .call_sample_silently().
# User-supplied output directories are never marked as owned and are preserved.
.cleanup_cmdstan_fit_output <- function(fit) {
  owned <- attr(fit, "pclv_owned_output_dir", exact = TRUE)
  if (is.character(owned) && length(owned) == 1L && nzchar(owned)) {
    unlink(owned, recursive = TRUE, force = TRUE)
    attr(fit, "pclv_owned_output_dir") <- NULL
  }
  invisible(NULL)
}


#' Convert Pathfinder draws into per-chain init lists
#'
#' @description
#' Takes posterior approximation draws returned by a
#' \code{CmdStanPathfinder} run and constructs a list of named
#' parameter initial values suitable for passing to
#' \code{cmdstanr::sample(init=...)}. Each chain receives a separate
#' named list of parameter values. If the draws are unusable (e.g.,
#' no matching parameters, NA/Inf values, or zero-length lists), the
#' function returns \code{NULL}.
#'
#' @param pf_fit A \code{CmdStanPathfinder} object (result of
#'   \code{mod$pathfinder()}).
#' @param mod The compiled \code{cmdstanr} model (not used currently,
#'   reserved for future extensions).
#' @param chains Integer number of chains to generate initial values
#'   for.
#' @param prefer Character vector of parameter names to prioritize
#'   when extracting from the Pathfinder draws. Defaults to common
#'   gLV parameters (\code{r0}, \code{a_ii}, \code{a_ij}, \code{sigma},
#'   \code{sd_ou}, \code{phi}, \code{lambda}, \code{sigma_ou},
#'   \code{tau_r}, \code{nu}).
#'
#' @return A list of length \code{chains}, each element being a named
#'   list of numeric initial values. Returns \code{NULL} if conversion
#'   fails.
#'
#' @examples
#' \dontrun{
#' mod <- cmdstanr::cmdstan_model("glv_pairwise.stan")
#' pf_fit <- mod$pathfinder(data = data_list)
#' inits <- .pf_inits_from_draws(pf_fit, mod, chains = 4)
#' fit <- mod$sample(data = data_list, chains = 4, init = inits)
#' }
#'
#' @keywords internal
#' @noRd
#' Convert Pathfinder draws into per-chain init lists
#'
#' Takes approximation draws from a CmdStanPathfinder fit and returns a
#' per-chain list of named numeric scalars suitable for `sample(init=...)`.
#' If no usable values are found, returns NULL; the caller records Pathfinder unavailability before using an approved ordinary initialization.
#' @noRd
#' Convert Pathfinder draws into per-chain init lists (strict)
#'
#' Takes approximation draws from a CmdStanPathfinder fit and returns a
#' per-chain list of named numeric scalars for `sample(init=...)`.
#' Returns `NULL` if no usable values can be constructed.
#' @noRd
.pf_inits_from_draws <- function(pf_fit, mod, chains = 4L,
                                 prefer = c("r0","a_ii","a_ij","sigma","sd_ou","phi",
                                            "lambda","sigma_ou","tau_r","nu")) {
  # 0) 모델 파라미터 집합
  model_params <- try(mod$variables()$parameters, silent = TRUE)
  if (inherits(model_params, "try-error") || is.null(model_params)) {
    model_params <- character(0)
  }

  # 1) draws 추출
  draws_df <- NULL
  try({
    dd <- pf_fit$draws()
    if (!is.null(dd)) draws_df <- posterior::as_draws_df(dd)
  }, silent = TRUE)
  if (is.null(draws_df)) return(NULL)

  # 2) PF draws 컬럼 ∩ 모델 파라미터 ∩ prefer
  vars_pf <- names(draws_df)
  pick <- intersect(prefer, intersect(vars_pf, model_params))
  if (!length(pick)) return(NULL)

  # 3) 체인 수만큼 행 선택 (부족하면 복제)
  nrow_use <- min(chains, nrow(draws_df))
  row_ids  <- seq_len(nrow_use)
  out <- vector("list", chains)
  for (k in seq_along(row_ids)) {
    r <- row_ids[k]
    vals <- as.list(draws_df[r, pick, drop = FALSE])
    keep <- vapply(vals, function(v) is.numeric(v) && length(v) == 1L && is.finite(v), logical(1))
    vals <- vals[keep]
    out[[k]] <- vals
  }
  # 4) 부족분은 마지막 유효 리스트로 채움 (모두 빈 경우는 NULL 반환)
  last_good <- NULL
  for (k in seq_len(nrow_use)) if (length(out[[k]]) > 0L) last_good <- out[[k]]
  if (is.null(last_good)) return(NULL)
  if (nrow_use < chains) for (k in (nrow_use + 1L):chains) out[[k]] <- last_good

  # 5) 모든 체인이 non-empty named list인지 확인
  ok <- all(vapply(out, function(li) is.list(li) && length(li) > 0L && length(names(li)) > 0L, logical(1)))
  if (!ok) return(NULL)
  out
}

#' Strict wrapper for \code{smooth.spline()} without scientific fallbacks
#'
#' Fits only the canonical cross-validated spline and returns a structured
#' failure instead of substituting another smoothing method.
#'
#' @param x,y Numeric vectors.
#' @param spline_df Optional effective degrees of freedom.
#' @param spline_spar Optional smoothing parameter.
#' @param use_cv Logical; try LOOCV before GCV.
#' @param min_unique Minimum distinct \code{x} needed for spline fit.
#' @param min_df,max_df Bounds for degrees of freedom; \code{NULL} auto-sets \code{max_df}.
#' @param default_df Default df if all attempts fail.
#' @return A list with \code{yhat} (fitted at \code{x}) and \code{df}.
#' @noRd
#' @keywords internal
.smooth_spline_robust <- function(x, y,
                                  spline_df = NULL, spline_spar = NULL,
                                  use_cv = TRUE, min_unique = 3,
                                  min_df = 3.0, max_df = NULL, default_df = 4.0) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (!all(ok)) return(.pclv_failure("spline_smoothing", "non_finite_spline_input", list()))
  if (!isTRUE(use_cv) || !is.null(spline_df) || !is.null(spline_spar)) {
    stop("Canonical Core spline invariant violated: CV smoothing is required.")
  }
  if (length(unique(x)) < min_unique) {
    return(.pclv_failure("spline_smoothing", "insufficient_unique_times",
                         list(observed = length(unique(x)), required = min_unique)))
  }
  ord <- order(x); x_ord <- x[ord]; y_ord <- y[ord]
  fit <- tryCatch(
    withCallingHandlers(stats::smooth.spline(x_ord, y_ord, cv = TRUE, all.knots = TRUE), warning = function(w) stop(conditionMessage(w), call. = FALSE)),
    error = function(e) .pclv_failure("spline_smoothing", "cv_spline_failed",
                                      list(message = conditionMessage(e)))
  )
  if (.is_pclv_failure(fit)) return(fit)
  yhat <- tryCatch(
    as.numeric(stats::predict(fit, x = x_ord)$y),
    error = function(e) .pclv_failure("spline_smoothing", "cv_spline_prediction_failed",
                                      list(message = conditionMessage(e)))
  )
  if (.is_pclv_failure(yhat)) return(yhat)
  if (length(yhat) != length(x_ord) || any(!is.finite(yhat))) {
    return(.pclv_failure("spline_smoothing", "invalid_cv_spline_prediction", list()))
  }
  inverse <- order(ord)
  list(yhat = yhat[inverse], df = as.numeric(fit$df), method = "cv")
}

#' Stable log-mean-exp
#'
#' Computes \eqn{\log(\mathrm{mean}(\exp(x)))} in a numerically stable manner.
#'
#' @param x Numeric vector.
#' @return A scalar on the log scale.
#' @noRd
#' @keywords internal
.log_mean_exp <- function(x) {
  xm <- max(x)
  xm + log(mean(exp(x - xm)))
}

#' Sampling wrapper with diagnostics-aware retries
#'
#' Runs \code{mod$sample()} and, if needed, retries with safer hyperparameters
#' (e.g., using higher \code{adapt_delta}, \code{dense_e}) until
#' diagnostics pass or \code{max_retries} is reached.
#'
#' @param mod A \pkg{cmdstanr} model.
#' @param base_args Baseline argument list for \code{mod$sample()}.
#' @param stan_list Data list; inserted into \code{base_args$data}.
#' @param max_retries Maximum retry attempts.
#' @param silent_sampler Logical; silence sampler output.
#' @param ebfmi_thresh E-BFMI threshold triggering retries.
#' @param tag Optional label for progress messages.
#' @param freeze_retry_hypers If \code{TRUE}, do not alter hyperparameters on retries.
#' @return A list with \code{fit}, \code{diag}, \code{n_retries},
#'   \code{fit_failed}, and \code{final_args}.
#' @noRd
#' @keywords internal
.sample_with_retry <- function(mod, base_args, stan_list,
                               max_retries = 3,
                               silent_sampler = FALSE,
                               ebfmi_thresh = 0.30,
                               tag = "",
                               freeze_retry_hypers = FALSE) {

  # 안전 기본값 헬퍼
  .or <- function(x, y) if (is.null(x)) y else x

  # 안전 init
  .init_safe <- function(stan_list, mod) {
    function(chain_id) {
      lst <- list()
      params <- try(names(mod$variables()$parameters), silent = TRUE)
      if (inherits(params, "try-error") || is.null(params)) params <- character(0)
      add <- function(n, v) if (n %in% params) lst[[n]] <<- v

      add("r0", 0); add("a_ii", 0); add("a_ij", 0)
      add("sigma", 0.20); add("sd_ou", 0.40); add("phi", 0.80)
      add("log_nu_minus_two", log(3))

      lst
    }
  }
  # 재시도 필요 여부 판단 (EBFMI + div + treedepth)
  .needs_retry <- function(diag, fit = NULL, thr = 0.30) {
    low_eb <- is.finite(diag$ebfmi_min) && (diag$ebfmi_min < thr)
    if (!low_eb && !is.null(fit)) {
      eb_vec <- diag$ebfmi_chain
      if (length(eb_vec)) low_eb <- any(is.finite(eb_vec) & (eb_vec < thr))
      if (!low_eb) low_eb <- .ebfmi_warn_from_fit(fit, thr = thr)
    }
    low_eb ||
      (is.finite(diag$n_divergent)     && diag$n_divergent > 0) ||
      (is.finite(diag$n_treedepth_hit) && diag$n_treedepth_hit > 0)
  }

  # 포맷터(로그용)
  .fmt_diag <- function(d)
    sprintf("div=%s, treedepth=%s, rhat=%s, ess_bulk=%s",
            ifelse(is.finite(d$n_divergent), d$n_divergent, NA),
            ifelse(is.finite(d$n_treedepth_hit), d$n_treedepth_hit, NA),
            ifelse(is.finite(d$worst_rhat), sprintf("%.3f", d$worst_rhat), NA),
            ifelse(is.finite(d$min_ess_bulk), d$min_ess_bulk, NA))

  # Malformed caller initialization is a programmer error, never repaired.
  init0 <- base_args$init
  valid_scalar <- is.numeric(init0) && length(init0) == 1L && is.finite(init0)
  valid_init <- is.null(init0) || valid_scalar || is.function(init0) || is.list(init0) || inherits(init0, "CmdStanPathfinder") || is.environment(init0)
  if (!valid_init) stop("Unsupported or malformed init specification.")

  # 준비
  base_args <- base_args
  base_args$data <- stan_list

  attempt <- 0L
  fit <- NULL; diag <- NULL; fit_failed <- FALSE
  final_args <- NULL
  attempt_history <- list()
  fit_owner_transferred <- FALSE
  on.exit({
    if (!fit_owner_transferred) .cleanup_cmdstan_fit_output(fit)
  }, add = TRUE)

  repeat {

    sample_args <- base_args
    if (attempt > 0L) {
      if (!freeze_retry_hypers) {
        sample_args$adapt_delta <- max(0.995, .or(base_args$adapt_delta, 0))
        sample_args$iter_warmup <- max(2000,  .or(base_args$iter_warmup, 1000))
        sample_args$metric      <- "dense_e"
        sample_args$step_size   <- 0.03
        if (is.numeric(base_args$init) && length(base_args$init) == 1L && is.finite(base_args$init)) {
          sample_args$init <- base_args$init
        } else {
          sample_args$init <- NULL  # 전 파라미터 자동 초기화
        }
      }
      # 리트라이에서도 seed 고정
      sample_args$seed <- as.integer(base_args$seed + attempt)
    }

    # --- init 가공: 타입 안전 처리 ---
    if (inherits(sample_args$init, "CmdStanPathfinder") || is.environment(sample_args$init)) {
      # Pathfinder 객체/환경은 여기서 손대지 않는다 (사전에 per-chain 변환하는게 원칙)
      # 단, 일부 cmdstanr 버전에서는 지원하지 않으므로 .call_sample_silently()에서 폴백 처리
    } else if (is.list(sample_args$init) && !is.function(sample_args$init)) {
      # Resolve number of chains (prefer 'chains', then fallback to 'parallel_chains')
      as_pos_int1 <- function(x) {
        if (is.null(x)) return(NA_integer_)
        x <- suppressWarnings(as.numeric(x))
        if (length(x) == 1L && is.finite(x) && !is.na(x) && x >= 1 && floor(x) == x) {
          return(as.integer(x))
        }
        NA_integer_
      }
      n <- as_pos_int1(sample_args$chains)
      if (is.na(n)) {
        n <- as_pos_int1(sample_args$parallel_chains)
        # keep internal consistency if chains was missing
        if (!is.na(n) && is.null(sample_args$chains)) sample_args$chains <- n
      }
      if (is.na(n)) n <- 1L

      init_obj <- sample_args$init

      # If init is already list-of-lists (per-chain), validate length & non-empty.
      # Else, replicate a single named list across chains.
      if (length(init_obj) > 0L && is.list(init_obj[[1L]]) && !is.null(names(init_obj[[1L]]))) {
        if (length(init_obj) != n) {
          stop("'init' length (", length(init_obj), ") must equal number of chains (", n, ").")
        }
        if (any(vapply(init_obj, function(x) length(x) == 0L, logical(1)))) {
          stop("'init' contains empty lists.")
        }
        # leave as-is
      } else {
        # 단일 파라미터 사전(named list)만 허용
        if (is.null(names(init_obj)) || !length(init_obj)) {
          stop("'init' must be a named list of parameter values or a per-chain list.")
        }
        sample_args$init <- replicate(n, init_obj, simplify = FALSE)
      }
    }
    message(sprintf("⏩ %sattempt %d/%d | metric=%s, adapt_delta=%s, warmup=%s, step_size=%s",
                    if (nzchar(tag)) paste0("[", tag, "] ") else "",
                    attempt, max_retries,
                    .or(sample_args$metric, "NA"),
                    ifelse(is.null(sample_args$adapt_delta), "NA", sprintf("%.3f", sample_args$adapt_delta)),
                    .or(sample_args$iter_warmup, "NA"),
                    .or(sample_args$step_size, "NA")))

    init_method <- if (is.null(sample_args$init)) "default" else if (is.numeric(sample_args$init)) "scalar" else if (is.function(sample_args$init)) "function" else if (is.list(sample_args$init)) "list" else "pathfinder"

    fit <- tryCatch(
      .call_sample_silently(mod, sample_args, silent = silent_sampler),
      error = function(e) .pclv_failure(
        "sampling", "cmdstan_execution_failed",
        list(message = conditionMessage(e), attempt = attempt)
      )
    )
    final_args <- sample_args
    if (.is_pclv_failure(fit)) {
      attempt_history[[length(attempt_history) + 1L]] <- list(
        attempt = attempt + 1L, seed = sample_args$seed, initialization_method = init_method,
        adapt_delta = sample_args$adapt_delta, max_treedepth = sample_args$max_treedepth,
        status = "failed", failure_reason = fit$reason
      )
      if (attempt >= max_retries) {
        fit$details$attempt_history <- attempt_history
        return(list(fit = NULL, diag = NULL, n_retries = attempt,
                    fit_failed = TRUE, final_args = final_args, failure = fit, retry_history = attempt_history))
      }
      attempt <- attempt + 1L
      next
    }

    diag <- .summarise_sampler_diag(
      fit,
      max_treedepth = .or(sample_args$max_treedepth, base_args$max_treedepth)
    )

    eb_vec <- diag$ebfmi_chain
    eb_str <- if (length(eb_vec)) paste(sprintf("%.3f", eb_vec), collapse=",") else "NA"
    ok <- !.needs_retry(diag, fit, thr = ebfmi_thresh)
    attempt_history[[length(attempt_history) + 1L]] <- list(
      attempt = attempt + 1L, seed = sample_args$seed, initialization_method = init_method,
      adapt_delta = sample_args$adapt_delta, max_treedepth = sample_args$max_treedepth,
      status = if (ok) "success" else "diagnostic_failure",
      failure_reason = if (ok) NA_character_ else "sampler_diagnostics_failed"
    )

    message(sprintf("%s%s | %s | E-BFMI chains=[%s]",
                    if (nzchar(tag)) paste0("[",tag,"] ") else "",
                    if (ok) "sampler health ok; convergence summary pending" else
                      "sampler health warning; convergence summary pending",
                    .fmt_diag(diag), eb_str))

    if (ok || attempt >= max_retries) {
      if (!ok) {
        fit_failed <- TRUE
      }
      break
    }
    .cleanup_cmdstan_fit_output(fit)
    attempt <- attempt + 1L
  }

  out <- list(
    fit = fit, diag = diag, n_retries = attempt,
    fit_failed = fit_failed, final_args = final_args,
    retry_history = attempt_history,
    failure = if (fit_failed) .pclv_failure(
      "sampling", "sampler_diagnostics_failed",
      list(n_retries = attempt, diagnostics = diag, attempt_history = attempt_history)
    ) else NULL
  )
  # A completed fit remains useful for chain-specific interpretation even when
  # retry diagnostic gates are exhausted. Its caller owns output cleanup.
  fit_owner_transferred <- !.is_pclv_failure(fit)
  out
}

.failed_direction_result <- function(failure) {
  na_diag <- list(worst_rhat = NA_real_, min_ess_bulk = NA_real_,
                  min_ess_tail = NA_real_, n_divergent = NA_integer_,
                  n_treedepth_hit = NA_integer_, ebfmi_min = NA_real_,
                  ebfmi_med = NA_real_)
  list(
    ok = FALSE, failure = failure, n_pairs = NA_integer_,
    a_mean = NA_real_, a_sd = NA_real_, a_q2.5 = NA_real_, a_q97.5 = NA_real_,
    p_sign2 = NA_real_, aii_mean = NA_real_, aii_sd = NA_real_,
    aii_q2.5 = NA_real_, aii_q97.5 = NA_real_, p_sign2_self = NA_real_,
    nu_mean = NA_real_, nu_median = NA_real_, nu_q05 = NA_real_, nu_q95 = NA_real_,
    diagnostic_class = "sampler_diagnostics_failed",
    interaction_identifiable = FALSE, residual_identifiable = FALSE,
    chain_sign_agreement = NA, pooled_sign_probability = NA_real_,
    chain_aij_means = list(numeric()), chain_aij_medians = list(numeric()),
    chain_aij_positive_probabilities = list(numeric()), chain_aij_negative_probabilities = list(numeric()),
    chain_aij_mean_range = NA_real_, chain_aij_median_range = NA_real_,
    chain_sign_probability_max_difference = NA_real_,
    chain_aii_means = list(numeric()), chain_aii_medians = list(numeric()),
    chain_aii_sign_probabilities = list(numeric()), chain_aii_sign_agreement = NA,
    chain_residual_summary = list(data.frame()), residual_median_ranges = list(numeric()),
    residual_regime_disagreement = NA, indeterminate_reason = list(failure$reason),
    diag = na_diag, retry_history = list(NULL), initialization_provenance = list(NULL),
    kfold_mean = NA_real_, kfold_method = NA_character_,
    kfold_subject = list(NULL), kfold_subject_ppd = list(NULL),
    kfold_subject_ids = list(NULL), kfold_subject_counts = list(NULL),
    kfold_subject_success = list(NULL), kfold_subject_fail = list(NULL),
    kfold_success_total = 0L, kfold_failures = list(list(failure)),
    kfold_splits = list(NULL), kfold_seed_used = NA_integer_,
    kfold_sd = NA_real_, kfold_se = NA_real_, kfold_n_subjects = 0L,
    kfold_retry_total = NA_integer_, kfold_retry_mean = NA_real_,
    kfold_nu_fold_means = list(NULL), kfold_outer_rounds = 0L,
    kfold_failed = TRUE, kfold_n_folds_ok = 0L, kfold_n_folds_fail = NA_integer_
  )
}

.pair_direction_seeds <- function(seed_base, idx_i, idx_j) {
  c(ij = as.integer(seed_base + 100000L * idx_i + 1000L * idx_j + 1L),
    ji = as.integer(seed_base + 100000L * idx_i + 1000L * idx_j + 2L))
}

.make_pair_tasks <- function(taxa_vec, seed_base) {
  idx <- utils::combn(seq_along(taxa_vec), 2L, simplify = FALSE)
  lapply(seq_along(idx), function(k) {
    ij <- idx[[k]]
    list(task_index = k, idx_i = ij[[1L]], idx_j = ij[[2L]],
         taxon_i = taxa_vec[[ij[[1L]]]], taxon_j = taxa_vec[[ij[[2L]]]],
         seed_base = as.integer(seed_base),
         direction_seeds = .pair_direction_seeds(seed_base, ij[[1L]], ij[[2L]]))
  })
}

.execute_pair_task <- function(task, taxa_vec, run_one, core_ctx, scheduling,
                               progress = "none", mute_logs = TRUE) {
  task_ctx <- core_ctx
  task_ctx$n_workers_kfold_eff <- scheduling$n_workers_kfold_eff
  result <- .run_pair(
    task$idx_i, task$idx_j, kfold_K = task_ctx$kfold_K, kfold_R = task_ctx$kfold_R,
    taxa_vec = taxa_vec, .run_one = run_one, ctx = task_ctx,
    progress = progress, mute_logs = mute_logs, seed_base = task$seed_base
  )
  list(task_index = task$task_index, result = result,
       direction_seeds = task$direction_seeds)
}

.execute_pair_task_progress <- function(task, taxa_vec, run_one, core_ctx,
                                        scheduling, progress, progressor) {
  result <- .execute_pair_task(
    task, taxa_vec, run_one, core_ctx, scheduling,
    progress = progress, mute_logs = TRUE
  )
  progressor(message = sprintf("pair %s-%s", task$taxon_i, task$taxon_j))
  result
}

.outer_pair_map <- function(tasks, taxa_vec, core_ctx, scheduling, progress,
                            workers, has_progressr,
                            .plan = future::plan,
                            .map = furrr::future_map,
                            .with_progress = progressr::with_progress) {
  previous_plan <- .plan()
  on.exit(.plan(previous_plan), add = TRUE)
  .plan(future::multisession, workers = workers)
  options <- furrr::furrr_options(
    seed = TRUE,
    globals = FALSE,
    packages = "pclvbayes",
    scheduling = Inf
  )
  common <- list(
    .x = tasks, taxa_vec = taxa_vec, run_one = .run_one,
    core_ctx = core_ctx, scheduling = scheduling, progress = progress,
    .options = options
  )
  if (isTRUE(has_progressr) && identical(progress, "bar")) {
    return(.with_progress({
      progressor <- progressr::progressor(steps = length(tasks))
      do.call(.map, c(common, list(
        .f = .execute_pair_task_progress,
        progressor = progressor
      )))
    }))
  }
  do.call(.map, c(common, list(
    .f = .execute_pair_task,
    mute_logs = TRUE
  )))
}

.assemble_pair_outcomes <- function(outcomes, tasks) {
  if (length(outcomes) != length(tasks)) stop("Pair outcome count invariant violated.")
  ord <- order(vapply(outcomes, `[[`, integer(1), "task_index"))
  observed <- vapply(outcomes[ord], `[[`, integer(1), "task_index")
  expected <- vapply(tasks, `[[`, integer(1), "task_index")
  if (!identical(observed, expected)) stop("Pair task ordering invariant violated.")
  dplyr::bind_rows(lapply(outcomes[ord], `[[`, "result"))
}

#' Run both directions (j→i and i→j) for a taxon pair
#'
#' Executes \code{.run_one()} twice with deterministic seeds and aggregates
#' posterior summaries, diagnostics, and repeated K-fold payloads into one row.
#'
#' @param idx_i,idx_j Integer indices into \code{taxa_vec}.
#' @param kfold_K,kfold_R Repeated K-fold settings (metadata pass-through).
#' @param taxa_vec Character vector of taxon names.
#' @param .run_one Callable for a single directed fit.
#' @param progress Progress mode string.
#' @param mute_logs Logical; suppress local progress.
#' @param seed_base Integer seed base for reproducibility.
#' @return A one-row tibble aggregating both directions, or \code{NULL}.
#' @noRd
#' @keywords internal
.run_pair <- function(idx_i, idx_j,kfold_K = NULL, kfold_R = NULL,
                      taxa_vec, .run_one, ctx, progress,
                      mute_logs = FALSE, seed_base) {
  ti <- taxa_vec[[idx_i]]
  pj <- taxa_vec[[idx_j]]

  direction_seeds <- .pair_direction_seeds(seed_base, idx_i, idx_j)
  seed_ij <- direction_seeds[["ij"]]
  seed_ji <- direction_seeds[["ji"]]

  res_ij <- .run_one(
    target = ti, partner = pj,
    ctx = ctx,
    seed_override = seed_ij,
    progress_local = if (mute_logs) "none" else progress
  )
  res_ji <- .run_one(
    target = pj, partner = ti,
    ctx = ctx,
    seed_override = seed_ji,
    progress_local = if (mute_logs) "none" else progress
  )

  if (is.null(res_ij)) res_ij <- .pclv_failure("directed_fit", "missing_result", list(direction = "ij"))
  if (is.null(res_ji)) res_ji <- .pclv_failure("directed_fit", "missing_result", list(direction = "ji"))

  failure_ij <- if (.is_pclv_failure(res_ij)) res_ij else NULL
  failure_ji <- if (.is_pclv_failure(res_ji)) res_ji else NULL
  if (!is.null(failure_ij)) res_ij <- .failed_direction_result(failure_ij)
  if (!is.null(failure_ji)) res_ji <- .failed_direction_result(failure_ji)

  # Older internal test doubles predate diagnostic_class; a completed result
  # without that private field retains the historical converged contract.
  class_ij <- if (isTRUE(res_ij$ok) && is.null(res_ij$failure)) "converged" else if (is.null(res_ij$diagnostic_class) || is.na(res_ij$diagnostic_class)) "converged" else res_ij$diagnostic_class
  class_ji <- if (isTRUE(res_ji$ok) && is.null(res_ji$failure)) "converged" else if (is.null(res_ji$diagnostic_class) || is.na(res_ji$diagnostic_class)) "converged" else res_ji$diagnostic_class
  diagnostic_failure_ij <- if (!is.null(failure_ij)) failure_ij else if (!identical(class_ij, "converged")) res_ij$diagnostic_failure[[1L]] else NULL
  diagnostic_failure_ji <- if (!is.null(failure_ji)) failure_ji else if (!identical(class_ji, "converged")) res_ji$diagnostic_failure[[1L]] else NULL

  tibble::as_tibble_row(list(
    i = ti, j = pj,
    direction_ok_ij = identical(class_ij, "converged") && is.null(failure_ij),
    direction_ok_ji = identical(class_ji, "converged") && is.null(failure_ji),
    failure_ij = list(diagnostic_failure_ij), failure_ji = list(diagnostic_failure_ji),
    diagnostic_class_ij = class_ij, diagnostic_class_ji = class_ji,
    interaction_identifiable_ij = res_ij$interaction_identifiable %||% TRUE,
    interaction_identifiable_ji = res_ji$interaction_identifiable %||% TRUE,
    residual_identifiable_ij = res_ij$residual_identifiable %||% TRUE,
    residual_identifiable_ji = res_ji$residual_identifiable %||% TRUE,
    chain_sign_agreement_ij = res_ij$chain_sign_agreement %||% NA,
    chain_sign_agreement_ji = res_ji$chain_sign_agreement %||% NA,
    pooled_sign_probability_ij = res_ij$pooled_sign_probability %||% NA_real_,
    pooled_sign_probability_ji = res_ji$pooled_sign_probability %||% NA_real_,
    chain_aij_means_ij = res_ij$chain_aij_means %||% list(numeric()),
    chain_aij_means_ji = res_ji$chain_aij_means %||% list(numeric()),
    chain_aij_medians_ij = res_ij$chain_aij_medians %||% list(numeric()),
    chain_aij_medians_ji = res_ji$chain_aij_medians %||% list(numeric()),
    chain_sign_probabilities_ij = res_ij$chain_aij_positive_probabilities %||% list(numeric()),
    chain_sign_probabilities_ji = res_ji$chain_aij_positive_probabilities %||% list(numeric()),
    chain_negative_sign_probabilities_ij = res_ij$chain_aij_negative_probabilities %||% list(numeric()),
    chain_negative_sign_probabilities_ji = res_ji$chain_aij_negative_probabilities %||% list(numeric()),
    chain_aij_mean_range_ij = res_ij$chain_aij_mean_range %||% NA_real_,
    chain_aij_mean_range_ji = res_ji$chain_aij_mean_range %||% NA_real_,
    chain_aij_median_range_ij = res_ij$chain_aij_median_range %||% NA_real_,
    chain_aij_median_range_ji = res_ji$chain_aij_median_range %||% NA_real_,
    chain_sign_probability_max_difference_ij = res_ij$chain_sign_probability_max_difference %||% NA_real_,
    chain_sign_probability_max_difference_ji = res_ji$chain_sign_probability_max_difference %||% NA_real_,
    chain_aii_means_ij = res_ij$chain_aii_means %||% list(numeric()),
    chain_aii_means_ji = res_ji$chain_aii_means %||% list(numeric()),
    chain_aii_medians_ij = res_ij$chain_aii_medians %||% list(numeric()),
    chain_aii_medians_ji = res_ji$chain_aii_medians %||% list(numeric()),
    chain_aii_sign_probabilities_ij = res_ij$chain_aii_sign_probabilities %||% list(numeric()),
    chain_aii_sign_probabilities_ji = res_ji$chain_aii_sign_probabilities %||% list(numeric()),
    chain_aii_sign_agreement_ij = res_ij$chain_aii_sign_agreement %||% NA,
    chain_aii_sign_agreement_ji = res_ji$chain_aii_sign_agreement %||% NA,
    chain_residual_summary_ij = res_ij$chain_residual_summary %||% list(data.frame()),
    chain_residual_summary_ji = res_ji$chain_residual_summary %||% list(data.frame()),
    residual_median_ranges_ij = res_ij$residual_median_ranges %||% list(numeric()),
    residual_median_ranges_ji = res_ji$residual_median_ranges %||% list(numeric()),
    residual_regime_disagreement_ij = res_ij$residual_regime_disagreement %||% NA,
    residual_regime_disagreement_ji = res_ji$residual_regime_disagreement %||% NA,
    indeterminate_reason_ij = res_ij$indeterminate_reason %||% list(character()),
    indeterminate_reason_ji = res_ji$indeterminate_reason %||% list(character()),
    retry_history_ij = res_ij$retry_history, retry_history_ji = res_ji$retry_history,
    initialization_provenance_ij = res_ij$initialization_provenance,
    initialization_provenance_ji = res_ji$initialization_provenance,
    n_pairs_ij = res_ij$n_pairs, n_pairs_ji = res_ji$n_pairs,
    a_ij_mean = res_ij$a_mean, a_ij_sd = res_ij$a_sd,
    a_ij_q2.5 = res_ij$a_q2.5, a_ij_q97.5 = res_ij$a_q97.5,
    p_sign2_ij = res_ij$p_sign2,
    a_ji_mean = res_ji$a_mean, a_ji_sd = res_ji$a_sd,
    a_ji_q2.5 = res_ji$a_q2.5, a_ji_q97.5 = res_ji$a_q97.5,
    p_sign2_ji = res_ji$p_sign2,
    nu_mean_ij = res_ij$nu_mean, nu_median_ij = res_ij$nu_median,
    nu_q05_ij = res_ij$nu_q05, nu_q95_ij = res_ij$nu_q95,
    nu_mean_ji = res_ji$nu_mean, nu_median_ji = res_ji$nu_median,
    nu_q05_ji = res_ji$nu_q05, nu_q95_ji = res_ji$nu_q95,
    a_ii_mean = res_ij$aii_mean, a_ii_sd = res_ij$aii_sd,
    a_ii_q2.5 = res_ij$aii_q2.5, a_ii_q97.5 = res_ij$aii_q97.5,
    p_sign2_ii = res_ij$p_sign2_self,
    a_jj_mean = res_ji$aii_mean, a_jj_sd = res_ji$aii_sd,
    a_jj_q2.5 = res_ji$aii_q2.5, a_jj_q97.5 = res_ji$aii_q97.5,
    p_sign2_jj = res_ji$p_sign2_self,
    rhat_ij  = res_ij$diag$worst_rhat, essb_ij = res_ij$diag$min_ess_bulk,
    esst_ij  = res_ij$diag$min_ess_tail, div_ij = res_ij$diag$n_divergent,
    tdhit_ij = res_ij$diag$n_treedepth_hit,
    ebfmi_min_ij = res_ij$diag$ebfmi_min, ebfmi_med_ij = res_ij$diag$ebfmi_med,
    rhat_ji  = res_ji$diag$worst_rhat, essb_ji = res_ji$diag$min_ess_bulk,
    esst_ji  = res_ji$diag$min_ess_tail, div_ji = res_ji$diag$n_divergent,
    tdhit_ji = res_ji$diag$n_treedepth_hit,
    ebfmi_min_ji = res_ji$diag$ebfmi_min, ebfmi_med_ji = res_ji$diag$ebfmi_med,
    kfold_elpd_mean_ij = res_ij$kfold_mean,
    kfold_elpd_mean_ji = res_ji$kfold_mean,
    kfold_elpd_method_ij = res_ij$kfold_method,
    kfold_elpd_method_ji = res_ji$kfold_method,
    kfold_K = if (is.null(kfold_K)) NA_integer_ else kfold_K, kfold_R = if (is.null(kfold_R)) NA_integer_ else kfold_R,
    kfold_subject_ij         = res_ij$kfold_subject,
    kfold_subject_ppd_ij     = res_ij$kfold_subject_ppd,
    kfold_subject_ids_ij     = res_ij$kfold_subject_ids,
    kfold_subject_counts_ij  = res_ij$kfold_subject_counts,
    kfold_subject_success_ij = res_ij$kfold_subject_success,
    kfold_subject_fail_ij    = res_ij$kfold_subject_fail,
    kfold_success_total_ij   = res_ij$kfold_success_total,
    kfold_failures_ij        = res_ij$kfold_failures,
    kfold_splits_ij          = res_ij$kfold_splits,
    kfold_seed_used_ij       = res_ij$kfold_seed_used,
    kfold_sd_ij              = res_ij$kfold_sd,
    kfold_se_ij              = res_ij$kfold_se,
    kfold_n_subjects_ij      = res_ij$kfold_n_subjects,
    kfold_retry_total_ij     = res_ij$kfold_retry_total,
    kfold_retry_mean_ij      = res_ij$kfold_retry_mean,
    kfold_nu_fold_means_ij  = res_ij$kfold_nu_fold_means,
    kfold_outer_rounds_ij    = res_ij$kfold_outer_rounds,
    kfold_failed_ij          = res_ij$kfold_failed,
    kfold_folds_ok_ij        = res_ij$kfold_n_folds_ok,
    kfold_folds_fail_ij      = res_ij$kfold_n_folds_fail,
    kfold_subject_ji         = res_ji$kfold_subject,
    kfold_subject_ppd_ji     = res_ji$kfold_subject_ppd,
    kfold_subject_ids_ji     = res_ji$kfold_subject_ids,
    kfold_subject_counts_ji  = res_ji$kfold_subject_counts,
    kfold_subject_success_ji = res_ji$kfold_subject_success,
    kfold_subject_fail_ji    = res_ji$kfold_subject_fail,
    kfold_success_total_ji   = res_ji$kfold_success_total,
    kfold_failures_ji        = res_ji$kfold_failures,
    kfold_splits_ji          = res_ji$kfold_splits,
    kfold_seed_used_ji       = res_ji$kfold_seed_used,
    kfold_sd_ji              = res_ji$kfold_sd,
    kfold_se_ji              = res_ji$kfold_se,
    kfold_n_subjects_ji      = res_ji$kfold_n_subjects,
    kfold_retry_total_ji     = res_ji$kfold_retry_total,
    kfold_retry_mean_ji      = res_ji$kfold_retry_mean,
    kfold_nu_fold_means_ji  = res_ji$kfold_nu_fold_means,
    kfold_outer_rounds_ji    = res_ji$kfold_outer_rounds,
    kfold_failed_ji          = res_ji$kfold_failed,
    kfold_folds_ok_ji        = res_ji$kfold_n_folds_ok,
    kfold_folds_fail_ji      = res_ji$kfold_n_folds_fail
  ))
}


.validate_nu_draws <- function(draws_df, stage) {
  if (!("nu" %in% names(draws_df))) {
    return(.pclv_failure(stage, "missing_nu_draws", list()))
  }
  nu <- as.numeric(draws_df$nu)
  if (length(nu) != nrow(draws_df) || any(!is.finite(nu)) || any(nu <= 2)) {
    return(.pclv_failure(stage, "invalid_nu_draws", list()))
  }
  NULL
}

.summarise_nu_draws <- function(draws_df) {
  failure <- .validate_nu_draws(draws_df, "posterior_extraction")
  if (!is.null(failure)) return(failure)
  nu <- as.numeric(draws_df$nu)
  list(nu_mean = mean(nu), nu_median = stats::median(nu),
       nu_q05 = unname(stats::quantile(nu, 0.05)),
       nu_q95 = unname(stats::quantile(nu, 0.95)))
}

#' Subject-level Student-t scale-mixture Kalman-OU log-likelihood
#'
#' Computes held-out subject log predictive densities under the canonical
#' irregular-time OU state process and Student-t observation likelihood. The
#' Student-t residual is represented as a Gamma-normal scale mixture. At each
#' observation, deterministic importance quadrature integrates the local
#' precision and the resulting Gaussian state mixture is collapsed by matching
#' its first two moments before the next OU prediction.
#'
#' This preserves the fitted Student-t tails and reduces exactly to the
#' Student-t density when latent-state uncertainty is zero. The moment collapse
#' is the only approximation; the previous variance-matched Gaussian scoring is
#' not used.
#'
#' @param draws_df Posterior draws data frame containing regression, Student-t,
#'   and OU parameters.
#' @param pair_in Test data with columns \code{y,xi,xj,subject,time}.
#' @return A list with matrix \code{full} (draws x subjects), subject IDs,
#'   observation counts, and scoring metadata.
#' @noRd
#' @keywords internal
.proj_loglik_subject <- function(draws_df, pair_in) {
  required_pair <- c("y", "xi", "xj", "subject", "time")
  if (!all(required_pair %in% names(pair_in))) {
    stop("pair_in lacks canonical predictive-scoring fields.")
  }
  if (!is.data.frame(draws_df) || !nrow(draws_df)) {
    return(.pclv_failure("elpd_scoring", "missing_posterior_draws", list()))
  }

  pair_in <- pair_in[order(pair_in$subject, pair_in$time), , drop = FALSE]
  time_failure <- .validate_subject_times(
    pair_in$time, pair_in$subject, "elpd_scoring"
  )
  if (!is.null(time_failure)) return(time_failure)

  required_draws <- c("r0", "a_ii", "a_ij", "sigma")
  if (!all(required_draws %in% names(draws_df))) {
    return(.pclv_failure(
      "elpd_scoring", "missing_predictive_draws",
      list(missing = setdiff(required_draws, names(draws_df)))
    ))
  }
  nu_failure <- .validate_nu_draws(draws_df, "elpd_scoring")
  if (!is.null(nu_failure)) return(nu_failure)
  core_draw_values <- unlist(draws_df[required_draws], use.names = FALSE)
  if (any(!is.finite(core_draw_values)) || any(draws_df$sigma <= 0)) {
    return(.pclv_failure("elpd_scoring", "invalid_predictive_draws", list()))
  }

  have_sdou <- "sd_ou" %in% names(draws_df)
  have_sigma_ou <- all(c("sigma_ou", "lambda") %in% names(draws_df))
  have_lambda <- "lambda" %in% names(draws_df)
  have_phi <- "phi" %in% names(draws_df)
  have_tau_r <- "tau_r" %in% names(draws_df)

  if (!have_sdou && !have_sigma_ou) {
    return(.pclv_failure(
      "elpd_scoring", "missing_ou_scale_draws", list()
    ))
  }
  if (!have_lambda && !have_phi) {
    return(.pclv_failure(
      "elpd_scoring", "missing_ou_persistence_draws", list()
    ))
  }

  invalid_ou <-
    (have_sdou && any(!is.finite(draws_df$sd_ou) | draws_df$sd_ou < 0)) ||
    (have_sigma_ou && any(
      !is.finite(draws_df$sigma_ou) | draws_df$sigma_ou < 0 |
        !is.finite(draws_df$lambda) | draws_df$lambda <= 0
    )) ||
    (have_lambda && any(!is.finite(draws_df$lambda) | draws_df$lambda <= 0)) ||
    (!have_lambda && have_phi && any(
      !is.finite(draws_df$phi) | draws_df$phi <= 0 | draws_df$phi >= 1
    )) ||
    (have_tau_r && any(!is.finite(draws_df$tau_r) | draws_df$tau_r < 0))
  if (invalid_ou) {
    return(.pclv_failure("elpd_scoring", "invalid_ou_draws", list()))
  }

  quadrature_probability <- .PCLV_CORE_T_KALMAN_QUADRATURE$probability
  quadrature_weight <- .PCLV_CORE_T_KALMAN_QUADRATURE$weight
  if (length(quadrature_probability) != length(quadrature_weight) ||
      !length(quadrature_probability) ||
      any(!is.finite(quadrature_probability)) ||
      any(quadrature_probability <= 0 | quadrature_probability >= 1) ||
      any(!is.finite(quadrature_weight)) || any(quadrature_weight <= 0) ||
      abs(sum(quadrature_weight) - 1) > 1e-12) {
    stop("Canonical Student-t quadrature invariant violated.")
  }

  D <- nrow(draws_df)
  Q <- length(quadrature_probability)
  nu_vec <- as.numeric(draws_df$nu)
  sigma2 <- as.numeric(draws_df$sigma)^2
  prior_shape <- nu_vec / 2
  prior_rate <- nu_vec / 2
  proposal_shape <- (nu_vec + 1) / 2

  row_log_sum_exp <- function(x) {
    row_max <- apply(x, 1L, max)
    shifted <- exp(x - row_max)
    row_max + log(rowSums(shifted))
  }

  # Deterministic importance quadrature for one observation. The proposal is
  # the conjugate local-precision posterior when state variance is zero and a
  # close approximation otherwise. This gives exact Student-t scoring in the
  # zero-state-uncertainty limit while retaining OU-state uncertainty.
  scale_mixture_weights <- function(innovation, state_variance) {
    state_variance <- pmax(as.numeric(state_variance), 0)
    proposal_rate <- (
      nu_vec + innovation^2 / pmax(sigma2 + state_variance, 1e-12)
    ) / 2

    probability_matrix <- matrix(
      rep(quadrature_probability, each = D), nrow = D, ncol = Q
    )
    shape_matrix <- matrix(
      rep(proposal_shape, times = Q), nrow = D, ncol = Q
    )
    rate_matrix <- matrix(
      rep(proposal_rate, times = Q), nrow = D, ncol = Q
    )
    precision <- stats::qgamma(
      probability_matrix,
      shape = shape_matrix,
      rate = rate_matrix
    )
    if (any(!is.finite(precision)) || any(precision <= 0)) {
      return(.pclv_failure(
        "elpd_scoring", "invalid_student_t_quadrature_precision", list()
      ))
    }

    observation_variance <- matrix(
      rep(state_variance, times = Q), nrow = D, ncol = Q
    ) + matrix(rep(sigma2, times = Q), nrow = D, ncol = Q) / precision
    observation_variance <- pmax(observation_variance, 1e-12)

    innovation_matrix <- matrix(
      rep(innovation, times = Q), nrow = D, ncol = Q
    )
    log_normal <- -0.5 * (
      log(2 * pi) + log(observation_variance) +
        innovation_matrix^2 / observation_variance
    )
    log_prior <- stats::dgamma(
      precision,
      shape = matrix(rep(prior_shape, times = Q), nrow = D, ncol = Q),
      rate = matrix(rep(prior_rate, times = Q), nrow = D, ncol = Q),
      log = TRUE
    )
    log_proposal <- stats::dgamma(
      precision,
      shape = shape_matrix,
      rate = rate_matrix,
      log = TRUE
    )
    log_quadrature_weight <- matrix(
      rep(log(quadrature_weight), each = D), nrow = D, ncol = Q
    )
    log_importance <-
      log_quadrature_weight + log_normal + log_prior - log_proposal

    log_predictive <- row_log_sum_exp(log_importance)
    normalized <- exp(log_importance - log_predictive)
    normalized <- normalized / rowSums(normalized)

    if (any(!is.finite(log_predictive)) || any(!is.finite(normalized))) {
      return(.pclv_failure(
        "elpd_scoring", "student_t_quadrature_failed", list()
      ))
    }

    list(
      log_predictive = log_predictive,
      weight = normalized,
      variance = observation_variance
    )
  }

  if (have_sdou) {
    stationary_variance <- as.numeric(draws_df$sd_ou)^2
  } else {
    stationary_variance <-
      as.numeric(draws_df$sigma_ou)^2 / (2 * as.numeric(draws_df$lambda))
  }
  stationary_variance <- pmax(stationary_variance, 0)

  use_random_intercept <- have_tau_r && any(draws_df$tau_r > 0)
  subjects <- unique(as.character(pair_in$subject))
  loglik_full <- matrix(
    NA_real_, nrow = D, ncol = length(subjects),
    dimnames = list(NULL, subjects)
  )
  n_obs <- integer(length(subjects))
  names(n_obs) <- subjects

  score_one_subject <- function(y, xi, xj, time) {
    J <- length(y)
    if (!J) return(rep(NA_real_, D))

    mu <- matrix(as.numeric(draws_df$r0), nrow = D, ncol = J) +
      tcrossprod(as.numeric(draws_df$a_ii), xi) +
      tcrossprod(as.numeric(draws_df$a_ij), xj)
    y_matrix <- matrix(rep(y, each = D), nrow = D, ncol = J)
    loglik <- rep(0, D)

    if (!use_random_intercept) {
      state_mean <- rep(0, D)
      state_variance <- stationary_variance
    } else {
      mean_ou <- rep(0, D)
      mean_intercept <- rep(0, D)
      var_ou <- stationary_variance
      var_intercept <- as.numeric(draws_df$tau_r)^2
      cov_ou_intercept <- rep(0, D)
    }

    delta_time <- c(0, diff(time))
    for (tt in seq_len(J)) {
      dt <- delta_time[[tt]]
      if (tt > 1L && (!is.finite(dt) || dt <= 0)) {
        stop("Validated OU time invariant violated.")
      }

      if (have_lambda) {
        persistence <- exp(-as.numeric(draws_df$lambda) * dt)
      } else {
        decay <- -log(as.numeric(draws_df$phi))
        persistence <- exp(-decay * dt)
      }
      process_variance <- pmax(
        stationary_variance * (1 - persistence^2), 0
      )

      if (!use_random_intercept) {
        predicted_mean <- persistence * state_mean
        predicted_variance <-
          persistence^2 * state_variance + process_variance
        predicted_variance <- pmax(predicted_variance, 0)
        innovation <- y_matrix[, tt] - mu[, tt] - predicted_mean

        mixture <- scale_mixture_weights(innovation, predicted_variance)
        if (.is_pclv_failure(mixture)) return(mixture)
        loglik <- loglik + mixture$log_predictive

        predicted_variance_matrix <- matrix(
          rep(predicted_variance, times = Q), nrow = D, ncol = Q
        )
        innovation_matrix <- matrix(
          rep(innovation, times = Q), nrow = D, ncol = Q
        )
        gain <- predicted_variance_matrix / mixture$variance
        component_mean <- matrix(
          rep(predicted_mean, times = Q), nrow = D, ncol = Q
        ) + gain * innovation_matrix
        component_variance <- pmax(
          predicted_variance_matrix -
            predicted_variance_matrix^2 / mixture$variance,
          0
        )

        state_mean <- rowSums(mixture$weight * component_mean)
        centered <- component_mean - state_mean
        state_variance <- rowSums(
          mixture$weight * (component_variance + centered^2)
        )
        state_variance <- pmax(state_variance, 1e-12)
      } else {
        predicted_mean_ou <- persistence * mean_ou
        predicted_mean_intercept <- mean_intercept
        predicted_var_ou <- persistence^2 * var_ou + process_variance
        predicted_cov <- persistence * cov_ou_intercept
        predicted_var_intercept <- var_intercept

        observed_state_variance <- pmax(
          predicted_var_ou + predicted_var_intercept + 2 * predicted_cov,
          0
        )
        innovation <-
          y_matrix[, tt] - mu[, tt] -
          predicted_mean_ou - predicted_mean_intercept

        mixture <- scale_mixture_weights(
          innovation, observed_state_variance
        )
        if (.is_pclv_failure(mixture)) return(mixture)
        loglik <- loglik + mixture$log_predictive

        cov_ou_observation <- predicted_var_ou + predicted_cov
        cov_intercept_observation <-
          predicted_var_intercept + predicted_cov
        cov_ou_matrix <- matrix(
          rep(cov_ou_observation, times = Q), nrow = D, ncol = Q
        )
        cov_intercept_matrix <- matrix(
          rep(cov_intercept_observation, times = Q),
          nrow = D, ncol = Q
        )
        innovation_matrix <- matrix(
          rep(innovation, times = Q), nrow = D, ncol = Q
        )

        gain_ou <- cov_ou_matrix / mixture$variance
        gain_intercept <- cov_intercept_matrix / mixture$variance
        component_mean_ou <- matrix(
          rep(predicted_mean_ou, times = Q), nrow = D, ncol = Q
        ) + gain_ou * innovation_matrix
        component_mean_intercept <- matrix(
          rep(predicted_mean_intercept, times = Q),
          nrow = D, ncol = Q
        ) + gain_intercept * innovation_matrix

        component_var_ou <- pmax(
          matrix(rep(predicted_var_ou, times = Q), nrow = D, ncol = Q) -
            cov_ou_matrix^2 / mixture$variance,
          0
        )
        component_var_intercept <- pmax(
          matrix(
            rep(predicted_var_intercept, times = Q),
            nrow = D, ncol = Q
          ) - cov_intercept_matrix^2 / mixture$variance,
          0
        )
        component_cov <-
          matrix(rep(predicted_cov, times = Q), nrow = D, ncol = Q) -
          cov_ou_matrix * cov_intercept_matrix / mixture$variance

        mean_ou <- rowSums(mixture$weight * component_mean_ou)
        mean_intercept <- rowSums(
          mixture$weight * component_mean_intercept
        )
        centered_ou <- component_mean_ou - mean_ou
        centered_intercept <- component_mean_intercept - mean_intercept

        var_ou <- rowSums(
          mixture$weight * (component_var_ou + centered_ou^2)
        )
        var_intercept <- rowSums(
          mixture$weight * (
            component_var_intercept + centered_intercept^2
          )
        )
        cov_ou_intercept <- rowSums(
          mixture$weight * (
            component_cov + centered_ou * centered_intercept
          )
        )

        var_ou <- pmax(var_ou, 1e-12)
        var_intercept <- pmax(var_intercept, 1e-12)
        covariance_bound <- sqrt(var_ou * var_intercept)
        cov_ou_intercept <- pmin(
          pmax(cov_ou_intercept, -covariance_bound),
          covariance_bound
        )
      }
    }

    loglik
  }

  for (s_idx in seq_along(subjects)) {
    subject_id <- subjects[[s_idx]]
    ix <- which(as.character(pair_in$subject) == subject_id)
    ix <- ix[order(pair_in$time[ix])]
    n_obs[[s_idx]] <- length(ix)
    subject_score <- score_one_subject(
      y = as.numeric(pair_in$y[ix]),
      xi = as.numeric(pair_in$xi[ix]),
      xj = as.numeric(pair_in$xj[ix]),
      time = as.numeric(pair_in$time[ix])
    )
    if (.is_pclv_failure(subject_score)) return(subject_score)
    loglik_full[, s_idx] <- subject_score
  }

  list(
    full = loglik_full,
    subjects = subjects,
    n_obs = n_obs,
    method = "student-t-scale-mixture-kalman-ou-q16",
    quadrature_nodes = Q,
    state_approximation = "gaussian-moment-collapse"
  )
}


#' Make repeated K-fold splits at the subject level
#'
#' Randomly partitions unique subjects into \code{K} folds, repeated \code{R} times.
#'
#' @param subject_vec Subject IDs aligned to rows of the data.
#' @param K Number of folds.
#' @param R Number of repetitions.
#' @param seed RNG seed used only for subject partitioning. With the same
#'   canonical subject universe, the same seed produces the same split manifest
#'   for every pair and direction.
#' @return A nested list \code{[[r]][[k]]} with \code{train_subjects}/\code{test_subjects}.
#' @noRd
#' @keywords internal
.make_repkfold_splits <- function(subject_vec,
                                  K = 5,
                                  R = 3,
                                  seed = 123) {
  set.seed(seed)
  # Canonicalize the subject universe before shuffling so row order or factor
  # level order cannot change the split manifest across pair directions.
  subs <- unique(as.character(subject_vec))
  if (anyNA(subs) || any(!nzchar(subs))) {
    stop("Subject IDs must be finite, non-missing, and non-empty.")
  }
  subs <- sort(subs, method = "radix")
  S <- length(subs)
  if (K > S)
    stop("K > #subjects")
  reps <- vector("list", R)
  for (r in seq_len(R)) {
    ord <- sample.int(S)
    folds <- split(subs[ord], rep(1:K, length.out = S))
    reps[[r]] <- lapply(seq_len(K), function(k) {
      list(test_subjects = folds[[k]],
           train_subjects = setdiff(subs, folds[[k]]))
    })
  }
  reps
}

#' Build train/test splits and prev/dt indices
#'
#' Creates train/test data frames with previous-row indices and inter-time \code{dt}
#' per subject and applies train-only global predictor scaling.
#'
#' @param pair_in Full input from \code{.make_pair_inputs_glv()}.
#' @param train_subjects,test_subjects Character vectors of subject IDs.
#' @return A list with \code{train}, \code{test}, and index vectors; or \code{NULL}.
#' @noRd
#' @keywords internal
.build_train_test <- function(pair_in,
                              train_subjects,
                              test_subjects,
                              min_pairs = 4) {
  mk_prev_dt <- function(df) {
    df <- df[order(df$subject, df$time), , drop = FALSE]
    N <- nrow(df)
    prev <- integer(N)
    dtv <- numeric(N)
    by_s <- split(seq_len(N), df$subject)
    for (sb in names(by_s)) {
      ix <- by_s[[sb]]
      prev[ix[1]] <- 0L
      dtv[ix[1]] <- 0
      if (length(ix) >= 2L) {
        for (k in 2:length(ix)) {
          prev[ix[k]] <- ix[k - 1]
          dt_k <- as.numeric(df$time[ix[k]] - df$time[ix[k - 1]])
          if (!is.finite(dt_k) || dt_k <= 0) stop("Canonical time invariant violated.")
          dtv[ix[k]] <- dt_k
        }
      }
    }
    list(df = df,
         prev = prev,
         dt = dtv)
  }

  if (all(c("xi_unscaled", "xj_unscaled") %in% names(pair_in))) {
    pair_in$xi <- pair_in$xi_unscaled
    pair_in$xj <- pair_in$xj_unscaled
    pair_in$xi_unscaled <- NULL
    pair_in$xj_unscaled <- NULL
  }
  tr <- subset(pair_in, subject %in% train_subjects)
  te <- subset(pair_in, subject %in% test_subjects)
  if (!nrow(tr) || !nrow(te))
    return(.pclv_failure("kfold_split", "empty_train_or_test", list()))
  trd <- mk_prev_dt(tr)
  ted <- mk_prev_dt(te)
  variation_failure <- .validate_predictor_variation(trd$df$xi, trd$df$xj, "kfold_training")
  if (!is.null(variation_failure)) {
    return(variation_failure)
  }

  sc <- list(
    xi_m = mean(trd$df$xi), xi_s = stats::sd(trd$df$xi),
    xj_m = mean(trd$df$xj), xj_s = stats::sd(trd$df$xj)
  )
  if (!is.finite(sc$xi_s) || sc$xi_s <= 0)
    return(.pclv_failure("kfold_training", "invalid_internal_scaling_state",
                         list(predictor = "xi", observed_sd = sc$xi_s)))
  if (!is.finite(sc$xj_s) || sc$xj_s <= 0)
    return(.pclv_failure("kfold_training", "invalid_internal_scaling_state",
                         list(predictor = "xj", observed_sd = sc$xj_s)))
  for (nm in c("xi", "xj")) {
    trd$df[[nm]] <- (trd$df[[nm]] - sc[[paste0(nm, "_m")]]) / sc[[paste0(nm, "_s")]]
    ted$df[[nm]] <- (ted$df[[nm]] - sc[[paste0(nm, "_m")]]) / sc[[paste0(nm, "_s")]]
  }

  # --- Guard: insufficient train samples per fold ---
  if (nrow(trd$df) < min_pairs)
    return(.pclv_failure("kfold_training", "insufficient_rows", list(observed_rows = nrow(trd$df), required_rows = min_pairs)))

  list(
    train = trd$df,
    test = ted$df,
    prev_train = trd$prev,
    dt_train = trd$dt,
    prev_test  = ted$prev,
    dt_test  = ted$dt
  )
}

#' Fit one (train/test) fold and compute subject-level ELPD
#'
#' Trains on \code{tr} subjects, scores on \code{te} subjects using
#' \code{.proj_loglik_subject()}, and returns per-subject ELPD vectors and fold diagnostics.
#'
#' @inheritParams .sample_with_retry
#' @param stan_list_base Base Stan data list (modified per fold).
#' @param sample_args_base Base sampler args (re-adapts per fold).
#' @param pair_in Full input data frame.
#' @param tr,te Character vectors of train/test subjects.
#' @param sample_args_override Optional override of sampler args.
#' @param seed_override Deterministic seed per fold.
#' @return A list with \code{elpd}, \code{elpd_ppd}, and \code{fold_diag}; or \code{NULL}.
#' @noRd
#' @keywords internal
.fold_fit_and_score <- function(mod,
                                stan_list_base,
                                sample_args_base,
                                pair_in,
                                train_subjects,
                                test_subjects,
                                max_retries = 3,
                                silent_sampler = TRUE,
                                sample_args_override = NULL,
                                freeze_retry_hypers = FALSE,
                                seed_override = NULL,
                                min_pairs = 4L) {
  pts <- .build_train_test(pair_in, train_subjects, test_subjects)
  if (inherits(pts, "pclv_failure")) {
    return(pts)
  }


  sl <- stan_list_base
  sl$N   <- nrow(pts$train)
  sl$y   <- pts$train$y
  sl$xi  <- pts$train$xi
  sl$xj  <- pts$train$xj
  sl$S   <- length(unique(pts$train$subject))
  sl$sid <- as.integer(factor(pts$train$subject))
  sl$prev <- pts$prev_train
  sl$dt  <- pts$dt_train


  # ★ 샘플러 제어는 공유하되, 적응 산출물(step_size/metric)은 폴드마다 재적응
  sa <- if (is.null(sample_args_override))
    sample_args_base
  else
    sample_args_override
  sa$data <- sl
  if (!is.null(seed_override))
    sa$seed <- seed_override   # 폴드별 결정적 seed
  # 폴드마다 재적응 유도: full-data 적응 결과 제거
  sa$step_size   <- NULL
  sa$inv_metric  <- NULL
  sa$metric_file <- NULL
  # ★ 풀런에서 넘어온 init(closure) 오염 차단: 폴드 데이터에 맞춰 재생성 또는 제거
  if (!(is.numeric(sa$init) &&
        length(sa$init) == 1L && is.finite(sa$init))) {
    sa$init <- NULL
  }

  # ★ 폴드 학습은 항상 체인 병렬 끔 (outer/메인 run과 중첩 병렬 방지)
  sa$parallel_chains <- 1L
  sa$chains <- if (is.null(sa$chains))
    4L
  else {
    x <- sa$chains
    if (is.character(x))
      x <- suppressWarnings(as.numeric(x))
    ok <- (length(x) == 1L) &&
      is.finite(x) && (x >= 1) && (floor(x) == x)
    if (ok)
      as.integer(x)
    else
      4L
  }

  tryfit <- .sample_with_retry(
    mod        = mod,
    base_args  = sa,
    stan_list  = sl,
    max_retries = max_retries,
    ebfmi_thresh = 0.30,
    silent_sampler = silent_sampler,
    tag = "kfold",
    freeze_retry_hypers = freeze_retry_hypers   # ★ 하이퍼를 바꾸지 않는 리트라이
  )
  fit <- tryfit$fit
  if (!is.null(tryfit$failure)) return(tryfit$failure)
  on.exit(.cleanup_cmdstan_fit_output(fit), add = TRUE)

  # 방어적 동일성 확인(디버그용; 필요시 주석 처리)
  fa <- tryfit$final_args

  d <- .safe_draws_df(fit)
  dg <- .add_convergence_diag(tryfit$diag, d)
  te_df <- pts$test
  llm <- .proj_loglik_subject(
    draws_df = d,
    pair_in = te_df
  )
  if (inherits(llm, "pclv_failure")) {
    return(llm)
  }
  # --- NEW: subject별 ELPD 벡터와 fold diagnostics 구성 ---
  # llm$full: (draws x n_test_subjects) 행렬, llm$subjects: 테스트 subject 벡터
  if (is.null(llm) ||
      is.null(llm$full) ||
      !is.matrix(llm$full) || ncol(llm$full) == 0) {
    # 테스트 데이터가 비정상인 폴드는 실패로 처리
    return(.pclv_failure("kfold_scoring", "invalid_pointwise_loglik", list()))
  }
  elpd_vec <- stats::setNames(apply(llm$full, 2, .log_mean_exp), llm$subjects)
  elpd_ppd_vec <- elpd_vec / pmax(llm$n_obs, 1L)

  fd <- data.frame(
    n_retries      = tryfit$n_retries,
    nu_mean        = mean(d$nu),
    ebfmi_min      = if (!is.null(dg$ebfmi_min))
      dg$ebfmi_min
    else
      NA_real_,
    worst_rhat     = if (!is.null(dg$worst_rhat))
      dg$worst_rhat
    else
      NA_real_,
    min_ess_bulk   = if (!is.null(dg$min_ess_bulk))
      dg$min_ess_bulk
    else
      NA_real_,
    treedepth_hits = if (!is.null(dg$n_treedepth_hit))
      dg$n_treedepth_hit
    else
      NA_integer_,
    n_divergent    = if (!is.null(dg$n_divergent))
      dg$n_divergent
    else
      NA_integer_
  )

  list(
    elpd = elpd_vec,
    elpd_ppd = elpd_ppd_vec,
    n_obs = stats::setNames(as.integer(llm$n_obs), llm$subjects),
    fold_diag = fd
  )
}

#' Repeated K-fold evaluation (subject-level) with diagnostics payload
#'
#' Runs \code{K × R} folds, aggregates subject-wise ELPD and PPD-normalized ELPD,
#' and returns fold diagnostics and split manifests for reproducibility.
#'
#' @inheritParams .fold_fit_and_score
#' @param seed Seed used only to build the repeated subject-level split manifest.
#'   Fold sampler seeds are derived separately from
#'   \code{sample_args_base$seed}.
#' @param n_workers_kfold Number of workers (multisession) for fold-level parallelism.
#' @return A list with aggregate ELPD summaries, diagnostics, splits, and counts.
#' @noRd
#' @keywords internal
.repkfold_eval <- function(mod,
                           stan_list_base,
                           sample_args_base,
                           pair_in,
                           K = 5,
                           R = 3,
                           seed = 123,
                           silent_sampler = TRUE,
                           max_retries = 3,
                           n_workers_kfold = 1L,
                           min_pairs = 4,
                           freeze_retry_hypers = FALSE,
                           progress = "none",
                           has_progressr = requireNamespace("progressr", quietly = TRUE)) { # retry computational settings remain fixed within folds
  # `seed` controls only the shared subject partition. Each directed model
  # retains its own sampler seed through sample_args_base$seed below.
  splits <- .make_repkfold_splits(pair_in$subject, K, R, seed)
  # --- compact split manifest for reproducibility (r,k,test_subjects)
  splits_df <- ({
    rows <- vector("list", K * R)
    t <- 0L
    for (r in seq_len(R)) {
      folds <- splits[[r]]
      for (k in seq_len(K)) {
        t <- t + 1L
        rows[[t]] <- data.frame(r = r,
                                k = k,
                                test_subjects = I(list(folds[[k]]$test_subjects)))
      }
    }
    do.call(rbind, rows)
  })
  # subject universe & how many times each was held out as test
  subs_all <- sort(unique(as.character(pair_in$subject)), method = "radix")
  test_counts_tab <- table(unlist(splits_df$test_subjects))
  subject_test_counts <- as.integer(test_counts_tab[subs_all])
  names(subject_test_counts) <- subs_all

  # task 리스트: (r, k, train_subjects, test_subjects, seed_offset)
  tasks <- list()
  for (r in seq_len(R)) {
    folds <- splits[[r]]
    for (k in seq_len(K)) {
      tr <- folds[[k]]$train_subjects
      te <- folds[[k]]$test_subjects
      base_seed <- if (!is.null(sample_args_base$seed))
        as.integer(sample_args_base$seed)
      else
        as.integer(seed)
      if (is.na(base_seed))
        base_seed <- 1L
      tasks[[length(tasks) + 1L]] <- list(
        r = r,
        k = k,
        tr = tr,
        te = te,
        seed = as.integer(base_seed + 1000L * r + k)
      )
    }
  }

  # 폴드별 실행 함수
  run_task <- function(task) {
    result <- .fold_fit_and_score(
      mod,
      stan_list_base,
      sample_args_base,
      pair_in,
      train_subjects = task$tr,
      test_subjects = task$te,
      silent_sampler = silent_sampler,
      max_retries = max_retries,
      freeze_retry_hypers = freeze_retry_hypers,
      seed_override = task$seed,
      min_pairs = min_pairs
    )
    if (is.null(result)) {
      result <- .pclv_failure("kfold_scoring", "fold_evaluation_failed",
                              list(predictor = NA_character_))
    }
    if (.is_pclv_failure(result)) {
      if (is.null(result$ok) || is.null(result$details)) {
        legacy_details <- result[setdiff(names(result), c("ok", "stage", "reason", "details"))]
        result <- .pclv_failure(result$stage, result$reason, legacy_details)
      }
      if (is.null(result$predictor)) result$predictor <- NA_character_
      result$repetition <- task$r
      result$fold <- task$k
      result$test_subjects <- task$te
      result$details$predictor <- result$predictor
      result$details$repetition <- task$r
      result$details$fold <- task$k
      result$details$test_subjects <- task$te
    }
    result
  }

  # 실행: 순차 또는 병렬
  if (n_workers_kfold > 1L) {
    # 안전하게 PSOCK(멀티세션) 사용 권장 (fork 이슈 회피)
    op <- future::plan()
    on.exit(future::plan(op), add = TRUE)
    future::plan(future::multisession, workers = n_workers_kfold)
    if (isTRUE(has_progressr) && identical(progress, "bar")) {
      progressr::with_progress({
        p <- progressr::progressor(steps = length(tasks))
        res_list <- furrr::future_map(tasks, function(task) {
          res <- run_task(task)
          p(message = sprintf("kfold r=%d k=%d", task$r, task$k))
          res
        }, .options = furrr::furrr_options(seed = TRUE))
      })
    } else {
      res_list <- furrr::future_map(tasks, run_task, .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    res_list <- purrr::map(tasks, run_task)
  }

  subs <- sort(unique(pair_in$subject))
  agg <- setNames(numeric(length(subs)), subs)
  cnt <- setNames(integer(length(subs)), subs)
  agg_ppd <- setNames(numeric(length(subs)), subs)
  cnt_ppd <- setNames(integer(length(subs)), subs)
  fail_cnt <- setNames(integer(length(subs)), subs)
  successful_obs <- setNames(integer(length(subs)), subs)

  fold_diag_df <- list()
  fold_failures <- list()

  n_ok <- 0L
  n_fail <- 0L

  for (el in res_list) {
    if (inherits(el, "pclv_failure")) {
      n_fail <- n_fail + 1L
      fold_failures[[length(fold_failures) + 1L]] <- el
      failed_subjects <- intersect(el$test_subjects, subs)
      fail_cnt[failed_subjects] <- fail_cnt[failed_subjects] + 1L
    } else {
      n_ok <- n_ok + 1L
      idx <- names(el$elpd)
      agg[idx] <- agg[idx] + el$elpd
      cnt[idx] <- cnt[idx] + 1L
      if (!is.null(el$elpd_ppd)) {
        agg_ppd[idx] <- agg_ppd[idx] + el$elpd_ppd
        cnt_ppd[idx] <- cnt_ppd[idx] + 1L
      }
      if (!is.null(el$n_obs)) {
        obs_idx <- intersect(names(el$n_obs), subs)
        successful_obs[obs_idx] <- successful_obs[obs_idx] + as.integer(el$n_obs[obs_idx])
      }
      fold_diag_df[[length(fold_diag_df) + 1L]] <- data.frame(
        n_retries = el$fold_diag$n_retries,
        nu_mean = el$fold_diag$nu_mean,
        ebfmi_min = el$fold_diag$ebfmi_min,
        worst_rhat = el$fold_diag$worst_rhat,
        min_ess_bulk = el$fold_diag$min_ess_bulk,
        treedepth_hits = el$fold_diag$treedepth_hits,
        n_divergent = el$fold_diag$n_divergent
      )
    }
  }
  elpd_subject <- rep(NA_real_, length(subs))
  names(elpd_subject) <- subs
  has_elpd <- cnt > 0L
  elpd_subject[has_elpd] <- agg[has_elpd] / cnt[has_elpd]
  elpd_subject_ppd <- rep(NA_real_, length(subs))
  names(elpd_subject_ppd) <- subs
  has_ppd <- cnt_ppd > 0L
  elpd_subject_ppd[has_ppd] <- agg_ppd[has_ppd] / cnt_ppd[has_ppd]
  elpd_mean_ppd <- if (length(elpd_subject_ppd)) {
    m <- mean(elpd_subject_ppd[is.finite(elpd_subject_ppd)], na.rm = TRUE)
    if (is.finite(m)) m else NA_real_
  } else NA_real_

  elpd_mean <- if (length(elpd_subject)) {
    m <- mean(elpd_subject[is.finite(elpd_subject)], na.rm = TRUE)
    if (is.finite(m)) m else NA_real_
  } else NA_real_


  fd <- if (length(fold_diag_df))
    do.call(rbind, fold_diag_df)
  else
    data.frame(
      n_retries = integer(),
      nu_mean = double(),
      ebfmi_min = double(),
      worst_rhat = double(),
      min_ess_bulk = double(),
      treedepth_hits = integer(),
      n_divergent = integer()
    )
  # 리트라이 요약
  retry_total <- if (nrow(fd))
    sum(fd$n_retries, na.rm = TRUE)
  else
    0L
  retry_mean  <- if (nrow(fd))
    mean(fd$n_retries, na.rm = TRUE)
  else
    NA_real_
  nu_fold_means <- if (nrow(fd)) fd$nu_mean else numeric()


  list(
    type = "repeated-kfold",
    K = K,
    R = R,
    elpd_subject = elpd_subject,
    elpd_subject_ppd = elpd_subject_ppd,
    elpd_mean    = elpd_mean,
    elpd_mean_ppd= elpd_mean_ppd,
    elpd_method = "student-t-scale-mixture-kalman-ou-q16",
    n_folds_ok = n_ok,
    n_folds_fail = n_fail,
    # reproducibility payload for subject-level predictive evidence
    splits_df = splits_df,
    subject_ids = subs_all,
    # Number of successful held-out observations contributing to each subject's
    # ELPD summary (legacy field retained for compatibility).
    subject_test_counts = successful_obs,
    # Number of folds in which each subject was scheduled for holdout.
    subject_holdout_counts = subject_test_counts,
    subject_success_counts = cnt,
    subject_failure_counts = fail_cnt,
    total_successful_evaluations = sum(cnt),
    fold_diag = fd,
    retry_total = retry_total,
    retry_mean  = retry_mean,
    nu_fold_means = nu_fold_means,
    failures = fold_failures
  )
}

#' Run a single directed regression (j → i)
#'
#' End-to-end pipeline for one direction: builds inputs, fits the model with
#' retry logic, computes sign probabilities, and (optionally) performs
#' repeated K-fold with train-only scaling and ν-consistency policy.
#'
#' @param target Target taxon name (i).
#' @param partner Partner taxon name (j).
#' @param seed_override Optional integer seed for reproducibility.
#' @param progress_local Progress mode string; inherits external \code{progress}.
#' @return A list with posterior summaries, diagnostics, and K-fold payload; or \code{NULL}.
#' @noRd
#' @keywords internal
.accepted_scale_audit <- function(draws, variable = "sigma") {
  if (!is.data.frame(draws) || !variable %in% names(draws))
    return(list(variable = variable, available = FALSE, reason = "variable_unavailable"))
  values <- as.numeric(draws[[variable]])
  chains <- if (".chain" %in% names(draws)) as.integer(draws$.chain) else rep(1L, length(values))
  chain_ids <- sort(unique(chains))
  chain_minimum <- setNames(vapply(chain_ids, function(chain) {
    x <- values[chains == chain & is.finite(values)]
    if (length(x)) min(x) else NA_real_
  }, numeric(1)), as.character(chain_ids))
  finite <- values[is.finite(values)]
  list(
    variable = variable,
    available = TRUE,
    minimum = if (length(finite)) min(finite) else NA_real_,
    maximum = if (length(finite)) max(finite) else NA_real_,
    zero_count = as.integer(sum(values == 0, na.rm = TRUE)),
    nonfinite_count = as.integer(sum(!is.finite(values))),
    chain_minimum = chain_minimum
  )
}

.fit_direction_main_posterior <- function(target,
                                          partner,
                                          ctx,
                                          seed_override = NULL,
                                          progress_local = progress) {

  # ---- 컨텍스트 바인딩(워커 환경 내에서 사용) ----
  meta_df              <- ctx$meta_df
  sm_mat               <- ctx$sm_mat
  eps                  <- ctx$eps
  min_pairs            <- ctx$min_pairs
  zero_mode_alr        <- ctx$zero_mode_alr
  minpos_alpha         <- ctx$minpos_alpha
  minpos_base          <- ctx$minpos_base
  eps_fixed            <- ctx$eps_fixed
  lib_eps_c            <- ctx$lib_eps_c
  rest_floor_frac      <- ctx$rest_floor_frac
  smooth_scale         <- ctx$smooth_scale
  alr_spline_df        <- ctx$alr_spline_df
  alr_spline_spar      <- ctx$alr_spline_spar
  alr_spline_cv        <- ctx$alr_spline_cv
  nz_partner_min_frac  <- ctx$nz_partner_min_frac
  max_retries          <- ctx$max_retries
  chains               <- ctx$chains
  iter_warmup          <- ctx$iter_warmup
  iter_sampling        <- ctx$iter_sampling
  adapt_delta          <- ctx$adapt_delta
  max_treedepth        <- ctx$max_treedepth
  metric               <- ctx$metric
  init                 <- ctx$init
  seed                 <- ctx$seed
  quiet                <- ctx$quiet
  silent_sampler       <- ctx$silent_sampler
  n_workers_kfold_eff  <- ctx$n_workers_kfold_eff
  kfold_K              <- ctx$kfold_K
  kfold_R              <- ctx$kfold_R
  kfold_seed           <- ctx$kfold_seed

  # PF 옵션
  use_pathfinder_init <- isTRUE(ctx$use_pathfinder_init)
  pf_num_paths        <- ctx$pf_num_paths
  pf_draws            <- ctx$pf_draws
  pf_history_size     <- ctx$pf_history_size
  pf_max_lbfgs_iters  <- ctx$pf_max_lbfgs_iters
  pf_psis_resample    <- ctx$pf_psis_resample

  # The exact compiled canonical model is loaded after preprocessing succeeds.

  # 안전 초기화(ELPD 비활성 시도 포함)
  kfold_outer_rounds_local <- 0L

  # 0) 페어 데이터 구성 ---------------------------------------------------------
  # allow custom pair builder (e.g., cross-kingdom mixing)
if (!is.null(ctx$pair_builder) && is.function(ctx$pair_builder)) {
  pair_df <- ctx$pair_builder(
    target = target,
    partner = partner,
    ctx = ctx,
    eps = eps,
    min_pairs = min_pairs
  )
  if (is.null(pair_df)) return(.pclv_failure("preprocessing", "pair_builder_no_data", list()))
} else {
  pair_df <- NULL
}

  pair_in <- .make_pair_inputs_glv(
    pair_df = pair_df,
    sm_mat = sm_mat, meta_df = meta_df,
    j = partner, i = target, min_pairs = min_pairs,
    zero_mode_alr = zero_mode_alr,
    minpos_alpha = minpos_alpha,
    minpos_base  = minpos_base,
    eps_fixed = eps_fixed,
    lib_eps_c = lib_eps_c,
    rest_floor_frac = rest_floor_frac,
    alr_cap = .PCLV_CORE_ALR_CAP,
    smooth_scale = smooth_scale,
    alr_spline_df = alr_spline_df,
    alr_spline_spar = alr_spline_spar,
    alr_spline_cv = alr_spline_cv,
    nz_partner_min_frac = nz_partner_min_frac
  )
  if (inherits(pair_in, "pclv_failure")) {
    return(pair_in)
  }
  if (nrow(pair_in) < min_pairs)
    return(.pclv_failure("preprocessing", "insufficient_analysis_rows", list(observed_rows = nrow(pair_in), required_rows = min_pairs)))

  # 결측/비유한 제거
  pair_in <- pair_in[is.finite(pair_in$y) &
                       is.finite(pair_in$xi) &
                       is.finite(pair_in$xj) &
                       is.finite(pair_in$time), , drop = FALSE]
  if (!nrow(pair_in))
    return(.pclv_failure("preprocessing", "no_finite_analysis_rows", list()))

  # 1) prev/dt 구성(AR(1)용) ----------------------------------------------------
  pair_in <- pair_in[order(pair_in$subject, pair_in$time), , drop = FALSE]
  N <- nrow(pair_in)
  prev <- integer(N)
  dtv <- numeric(N)
  by_s <- split(seq_len(N), pair_in$subject)
  for (sb in names(by_s)) {
    ix <- by_s[[sb]]
    prev[ix[1]] <- 0L
    dtv[ix[1]] <- 0
    if (length(ix) >= 2L) {
      for (k in 2:length(ix)) {
        prev[ix[k]] <- ix[k - 1]
        dt_k <- as.numeric(pair_in$time[ix[k]] - pair_in$time[ix[k - 1]])
        if (!is.finite(dt_k) || dt_k <= 0) stop("Canonical time invariant violated.")
        dtv[ix[k]] <- dt_k
      }
    }
  }

  n_subj <- length(unique(pair_in$subject))
  if (n_subj < 2L) {
    warning(
      "Fewer than 2 subjects for pair {",
      partner,
      "->",
      target,
      "} — AR(1) persistence may be weakly identified."
    )
  }

  sedf  <- attr(pair_in, "smooth_edf_mean")
  ssc   <- attr(pair_in, "smooth_scale")
  sflag <- isTRUE(attr(pair_in, "smoothed"))

  # 2) Stan data (기본) ---------------------------------------------------------
  stan_list <- c(
    list(
      N   = N,
      y   = pair_in$y,
      xi  = pair_in$xi,
      xj  = pair_in$xj,
      S   = n_subj,
      sid = as.integer(factor(pair_in$subject)),
      prev = prev,
      dt   = dtv
    )
  )

  # 앞쪽 값 우선으로 중복 키 제거
  nm <- names(stan_list)
  if (anyDuplicated(nm))
    stan_list <- stan_list[match(unique(nm), nm)]

  # Load only the exact compiled canonical model; workers never compile variants.
  if (is.null(ctx$mod_exe_file) || !file.exists(ctx$mod_exe_file)) {
    return(.pclv_failure("model_loading", "compiled_model_unavailable",
                         list(exe_file = ctx$mod_exe_file)))
  }
  mod <- tryCatch(
    cmdstanr::cmdstan_model(stan_file = NULL, exe_file = ctx$mod_exe_file),
    error = function(e) .pclv_failure("model_loading", "compiled_model_load_failed",
                                      list(message = conditionMessage(e), exe_file = ctx$mod_exe_file))
  )
  if (.is_pclv_failure(mod)) return(mod)

  # --- seed 결정(쌍별/방향별 재현성) -------------------------------------------
  seed_main <- if (!is.null(seed_override))
    as.integer(seed_override)
  else
    as.integer(seed)

  # 3) 샘플링 베이스 인자 --------------------------------------------------------
  base_args <- list(
    data = stan_list,
    seed = seed_main,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    metric = metric,
    init = init,
    refresh = if (quiet)
      0
    else
      100
  )

  # ----(NEW) Pathfinder로 좋은 초기값 만들기----
  # CmdStanR는 pathfinder fit을 init에 바로 받을 수 있음.
  # 체인 수보다 PF draw가 적으면 자동으로(with replacement) 뽑아 씀.
  # 참고: cmdstanr reference (model-method-pathfinder, init에 CmdStanPathfinder 허용). :contentReference[oaicite:1]{index=1}
  initialization_provenance <- list(
    requested = if (use_pathfinder_init) "pathfinder" else "user_or_default",
    pathfinder_status = if (use_pathfinder_init) "not_run" else "not_requested",
    pathfinder_failure = NULL,
    actual = if (is.null(init)) "default" else if (is.numeric(init)) "scalar" else if (is.function(init)) "function" else "list"
  )
  if (use_pathfinder_init) {
    pf_dir <- tempfile(pattern = "glvpf_", tmpdir = tempdir())
    dir.create(pf_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(if (dir.exists(pf_dir)) unlink(pf_dir, recursive = TRUE, force = TRUE), add = TRUE)
    pf_fit <- tryCatch(mod$pathfinder(
      data = stan_list, seed = seed_main, init = 0.05,
      num_paths = pf_num_paths, draws = pf_draws,
      history_size = pf_history_size, max_lbfgs_iters = pf_max_lbfgs_iters,
      psis_resample = pf_psis_resample, show_messages = !silent_sampler,
      output_dir = pf_dir, output_basename = paste0("glv_pf_", target, "_", partner)
    ), error = function(e) .pclv_failure("initialization", "pathfinder_failed",
                                         list(message = conditionMessage(e))))
    if (.is_pclv_failure(pf_fit)) {
      initialization_provenance$pathfinder_status <- "failed"
      initialization_provenance$pathfinder_failure <- pf_fit
    } else {
      pf_inits <- .pf_inits_from_draws(pf_fit, mod, chains = chains)
      if (is.null(pf_inits)) {
        initialization_provenance$pathfinder_status <- "failed"
        initialization_provenance$pathfinder_failure <- .pclv_failure(
          "initialization", "pathfinder_draws_unusable", list()
        )
      } else {
        base_args$init <- pf_inits
        initialization_provenance$pathfinder_status <- "success"
        initialization_provenance$actual <- "pathfinder"
        if (!is.null(iter_warmup) && iter_warmup >= 1500)
          base_args$iter_warmup <- max(500L, floor(iter_warmup / 2))
      }
    }
    if (dir.exists(pf_dir)) unlink(pf_dir, recursive = TRUE, force = TRUE)
  }
  # 4) 공통 래퍼로 샘플 + 리트라이 ----------------------------------------------
  tag_lbl  <- sprintf("%s\u2192%s main", partner, target)
  pair_tag <- sprintf("%s\u2192%s", partner, target)

  # --- main run start ---
  if (progress_local != "none")
    cat(sprintf("%s pair: main run started\n", pair_tag))

  res_try <- .sample_with_retry(
    mod = mod,
    base_args = base_args,
    stan_list = stan_list,
    max_retries = max_retries,
    silent_sampler = silent_sampler,
    ebfmi_thresh = 0.30,
    tag = tag_lbl
  )

  # --- main run done ---
  if (progress_local != "none")
    cat(sprintf("%s pair: main run completed\n", pair_tag))

  fit_full   <- res_try$fit
  if (!is.null(res_try$failure) && is.null(fit_full)) {
    res_try$failure$details$retry_history <- res_try$retry_history
    res_try$failure$details$initialization_provenance <- initialization_provenance
    return(res_try$failure)
  }
  on.exit(.cleanup_cmdstan_fit_output(fit_full), add = TRUE)
  diag       <- res_try$diag
  n_retries  <- res_try$n_retries
  fit_failed <- res_try$fit_failed
  final_args_main <- res_try$final_args

  # 5) 드로우 요약 ---------------------------------------------------------------
  if (quiet) {
    d <- suppressWarnings(.safe_draws_df(fit_full))
  } else {
    d <- .safe_draws_df(fit_full)  # a_ij, a_ii, r0, tau_r, sigma, sigma_ou, lambda, ...
  }
  diag <- .add_convergence_diag(diag, d)
  posterior_summary <- .build_posterior_summary_bundle(d, diag)
  if (.is_pclv_failure(posterior_summary)) return(posterior_summary)
  chain_diagnostics <- posterior_summary$chain
  if (!is.null(res_try$failure) && identical(chain_diagnostics$diagnostic_class, "converged")) {
    multi_chain <- length(chain_diagnostics$chain_aij_means) > 1L
    chain_diagnostics$diagnostic_class <- if (multi_chain)
      "interaction_stable_residual_unstable" else "sampler_diagnostics_failed"
    chain_diagnostics$residual_identifiable <- FALSE
    if (!multi_chain) chain_diagnostics$interaction_identifiable <- FALSE
    chain_diagnostics$indeterminate_reason <- unique(c(
      chain_diagnostics$indeterminate_reason, "sampler_diagnostics_failed"
    ))
  }
  diagnostic_failure <- NULL
  if (!identical(chain_diagnostics$diagnostic_class, "converged")) {
    diagnostic_failure <- .pclv_failure(
      "posterior_diagnostics", chain_diagnostics$diagnostic_class,
      list(
        reasons = chain_diagnostics$indeterminate_reason,
        interaction_identifiable = chain_diagnostics$interaction_identifiable,
        residual_identifiable = chain_diagnostics$residual_identifiable,
        retry_history = res_try$retry_history
      )
    )
  }

  interaction_summary <- posterior_summary$coefficients$interaction
  self_summary <- posterior_summary$coefficients$self
  p_sign2_dir <- posterior_summary$sign$interaction_p_two
  p_sign2_self <- posterior_summary$sign$self_p_two
  nu_summary <- posterior_summary$nu
  accepted_scale_audit <- .accepted_scale_audit(d, "sigma")
  rm(d)

  # Predictive evaluation is a separate downstream phase.
  predictive_context <- list(
    mod = mod, stan_list = stan_list, pair_in = pair_in,
    sample_args = final_args_main,
    # Split construction is shared across pair directions. The sampler seed
    # remains direction-specific inside sample_args.
    split_seed = as.integer(kfold_seed),
    sampling_seed = if (is.null(final_args_main$seed)) seed_main else final_args_main$seed,
    silent_sampler = silent_sampler, n_workers_kfold = n_workers_kfold_eff,
    max_retries = max_retries, min_pairs = min_pairs, K = kfold_K, R = kfold_R,
    pair_tag = pair_tag, progress = progress_local
  )
  kfold <- NULL

  # 7) 리턴 ----------------------------------------------------------------------
  list(
    n_pairs = nrow(pair_in),

    a_mean  = interaction_summary[["mean"]],
    a_median = interaction_summary[["median"]],
    a_sd = interaction_summary[["sd"]],
    a_q2.5  = interaction_summary[["q025"]],
    a_q97.5 = interaction_summary[["q975"]],
    positive_sign_probability = posterior_summary$sign$interaction_positive_probability,
    negative_sign_probability = posterior_summary$sign$interaction_negative_probability,
    lfsr = posterior_summary$sign$interaction_lfsr,
    p_sign2 = p_sign2_dir,

    aii_mean = self_summary[["mean"]],
    aii_sd = self_summary[["sd"]],
    aii_q2.5 = self_summary[["q025"]],
    aii_q97.5 = self_summary[["q975"]],
    p_sign2_self = p_sign2_self,

    nu_mean = nu_summary$nu_mean,
    nu_median = nu_summary$nu_median,
    nu_q05 = nu_summary$nu_q05,
    nu_q95 = nu_summary$nu_q95,

    diagnostic_class = chain_diagnostics$diagnostic_class,
    interaction_identifiable = chain_diagnostics$interaction_identifiable,
    residual_identifiable = chain_diagnostics$residual_identifiable,
    chain_sign_agreement = chain_diagnostics$chain_sign_agreement,
    pooled_sign_probability = chain_diagnostics$pooled_sign_probability,
    chain_aij_means = list(chain_diagnostics$chain_aij_means),
    chain_aij_medians = list(chain_diagnostics$chain_aij_medians),
    chain_aij_positive_probabilities = list(chain_diagnostics$chain_aij_positive_probabilities),
    chain_aij_negative_probabilities = list(chain_diagnostics$chain_aij_negative_probabilities),
    chain_aij_mean_range = chain_diagnostics$chain_aij_mean_range,
    chain_aij_median_range = chain_diagnostics$chain_aij_median_range,
    chain_sign_probability_max_difference = chain_diagnostics$chain_sign_probability_max_difference,
    chain_aii_means = list(chain_diagnostics$chain_aii_means),
    chain_aii_medians = list(chain_diagnostics$chain_aii_medians),
    chain_aii_sign_probabilities = list(chain_diagnostics$chain_aii_sign_probabilities),
    chain_aii_sign_agreement = chain_diagnostics$chain_aii_sign_agreement,
    chain_residual_summary = list(chain_diagnostics$chain_residual_summary),
    residual_median_ranges = list(chain_diagnostics$residual_median_ranges),
    residual_regime_disagreement = chain_diagnostics$residual_regime_disagreement,
    indeterminate_reason = list(chain_diagnostics$indeterminate_reason),
    diagnostic_failure = list(diagnostic_failure),

    # 원시 진단치/리트라이 메타(후처리 summariser에서 필터 예정)
    diag = diag,
    n_retries = n_retries,
    fit_failed = fit_failed,
    retry_history = list(res_try$retry_history),
    initialization_provenance = list(initialization_provenance),

    # 스무딩 메타
    smoothed = sflag,
    smooth_scale = ssc,
    smooth_edf_mean = sedf,
    subjects = unique(pair_in$subject),

    # Repeated K-fold result and reproducibility payload
    kfold = kfold,
    kfold_mean = if (is.null(kfold)) NA_real_ else kfold$elpd_mean,
    kfold_method = if (is.null(kfold)) NA_character_ else kfold$elpd_method,

    # 용어 명확화: 바깥 재평가 라운드(outer)만 카운트
    kfold_outer_rounds = kfold_outer_rounds_local,
    kfold_failed = is.null(kfold) ||
      (!is.null(kfold$n_folds_fail) && kfold$n_folds_fail > 0L),
    kfold_n_folds_ok   = if (is.null(kfold))
      NA_integer_
    else
      kfold$n_folds_ok,
    kfold_n_folds_fail = if (is.null(kfold))
      NA_integer_
    else
      kfold$n_folds_fail,
    # subject-level out-of-fold ELPD vector (named)
    kfold_subject       = list(if (is.null(kfold))
      NULL
      else
        kfold$elpd_subject),
    kfold_subject_ppd   = list(if (is.null(kfold))
      NULL
      else
        kfold$elpd_subject_ppd),
    kfold_subject_ids   = list(if (is.null(kfold))
      NULL
      else
        names(kfold$elpd_subject)),
    kfold_subject_counts = list(if (is.null(kfold))
      NULL
      else
        kfold$subject_test_counts),
    kfold_subject_success = list(if (is.null(kfold))
      NULL
      else
        kfold$subject_success_counts),
    kfold_subject_fail = list(if (is.null(kfold))
      NULL
      else
        kfold$subject_failure_counts),
    kfold_success_total = if (is.null(kfold))
      NA_integer_
    else
      kfold$total_successful_evaluations,
    kfold_failures = list(if (is.null(kfold))
      NULL
      else
        kfold$failures),
    # split manifest (r,k,test_subjects) & seed to ensure same partitions across models
    kfold_splits        = list(if (is.null(kfold))
      NULL
      else
        kfold$splits_df),
    kfold_seed_used     = NA_integer_,
    kfold_K             = if (is.null(kfold))
      NA_integer_
    else
      kfold$K,
    kfold_R             = if (is.null(kfold))
      NA_integer_
    else
      kfold$R,
    kfold_agg           = "subject-uniform",
    # convenient dispersion summaries
    kfold_sd            = if (is.null(kfold))
      NA_real_
    else
      stats::sd(kfold$elpd_subject, na.rm = TRUE),
    kfold_se            = if (is.null(kfold))
      NA_real_
    else {
      nn <- sum(is.finite(kfold$elpd_subject))
      stats::sd(kfold$elpd_subject, na.rm = TRUE) / sqrt(pmax(nn, 1L))
    },
    kfold_n_subjects    = if (is.null(kfold))
      NA_integer_
    else
      sum(is.finite(kfold$elpd_subject)),

    # 폴드-내 리트라이/nu 사용 요약
    kfold_retry_total  = if (is.null(kfold) ||
                             is.null(kfold$retry_total))
      NA_integer_
    else
      kfold$retry_total,
    kfold_retry_mean   = if (is.null(kfold) ||
                             is.null(kfold$retry_mean))
      NA_real_
    else
      kfold$retry_mean,
    kfold_nu_fold_means = list(if (is.null(kfold) ||
                                    is.null(kfold$nu_fold_means))
      NULL
      else
        kfold$nu_fold_means),
    .predictive_context = predictive_context,
    .accepted_scale_audit = accepted_scale_audit
  )
}

.add_predictive_evaluation <- function(main_result) {
  if (.is_pclv_failure(main_result)) return(main_result)
  predictive <- main_result$.predictive_context
  if (is.null(predictive)) stop("Main-posterior result lacks predictive context.")
  main_result$.predictive_context <- NULL
  if (!is.null(main_result$diagnostic_failure[[1L]])) return(main_result)

  if (predictive$progress != "none")
    cat(sprintf("%s pair: k-fold evaluation started (K=%d, R=%d)\n",
                predictive$pair_tag, predictive$K, predictive$R))
  sample_args <- predictive$sample_args
  sample_args$seed <- as.integer(predictive$sampling_seed)
  sample_args$step_size <- NULL
  sample_args$inv_metric <- NULL
  sample_args$metric_file <- NULL
  kfold <- .repkfold_eval(
    mod = predictive$mod,
    stan_list_base = predictive$stan_list,
    sample_args_base = sample_args,
    pair_in = predictive$pair_in,
    K = predictive$K,
    R = predictive$R,
    seed = predictive$split_seed,
    silent_sampler = predictive$silent_sampler,
    n_workers_kfold = predictive$n_workers_kfold,
    max_retries = predictive$max_retries,
    freeze_retry_hypers = TRUE,
    min_pairs = predictive$min_pairs,
    progress = predictive$progress
  )
  if (predictive$progress != "none")
    cat(sprintf("%s pair: k-fold evaluation completed\n", predictive$pair_tag))

  main_result$kfold <- kfold
  main_result$kfold_mean <- if (is.null(kfold)) NA_real_ else kfold$elpd_mean
  main_result$kfold_method <- if (is.null(kfold)) NA_character_ else kfold$elpd_method
  main_result$kfold_outer_rounds <- 1L
  main_result$kfold_failed <- is.null(kfold) ||
    (!is.null(kfold$n_folds_fail) && kfold$n_folds_fail > 0L)
  main_result$kfold_n_folds_ok <- if (is.null(kfold)) NA_integer_ else kfold$n_folds_ok
  main_result$kfold_n_folds_fail <- if (is.null(kfold)) NA_integer_ else kfold$n_folds_fail
  main_result$kfold_subject <- list(if (is.null(kfold)) NULL else kfold$elpd_subject)
  main_result$kfold_subject_ppd <- list(if (is.null(kfold)) NULL else kfold$elpd_subject_ppd)
  main_result$kfold_subject_ids <- list(if (is.null(kfold)) NULL else names(kfold$elpd_subject))
  main_result$kfold_subject_counts <- list(if (is.null(kfold)) NULL else kfold$subject_test_counts)
  main_result$kfold_subject_success <- list(if (is.null(kfold)) NULL else kfold$subject_success_counts)
  main_result$kfold_subject_fail <- list(if (is.null(kfold)) NULL else kfold$subject_failure_counts)
  main_result$kfold_success_total <- if (is.null(kfold)) NA_integer_ else kfold$total_successful_evaluations
  main_result$kfold_failures <- list(if (is.null(kfold)) NULL else kfold$failures)
  main_result$kfold_splits <- list(if (is.null(kfold)) NULL else kfold$splits_df)
  main_result$kfold_seed_used <- if (is.null(kfold)) NA_integer_ else predictive$split_seed
  main_result$kfold_K <- if (is.null(kfold)) NA_integer_ else kfold$K
  main_result$kfold_R <- if (is.null(kfold)) NA_integer_ else kfold$R
  main_result$kfold_sd <- if (is.null(kfold)) NA_real_ else stats::sd(kfold$elpd_subject, na.rm = TRUE)
  main_result$kfold_se <- if (is.null(kfold)) NA_real_ else {
    n <- sum(is.finite(kfold$elpd_subject))
    stats::sd(kfold$elpd_subject, na.rm = TRUE) / sqrt(pmax(n, 1L))
  }
  main_result$kfold_n_subjects <- if (is.null(kfold)) NA_integer_ else sum(is.finite(kfold$elpd_subject))
  main_result$kfold_retry_total <- if (is.null(kfold) || is.null(kfold$retry_total)) NA_integer_ else kfold$retry_total
  main_result$kfold_retry_mean <- if (is.null(kfold) || is.null(kfold$retry_mean)) NA_real_ else kfold$retry_mean
  main_result$kfold_nu_fold_means <- list(if (is.null(kfold) || is.null(kfold$nu_fold_means)) NULL else kfold$nu_fold_means)
  main_result
}

.run_one <- function(target, partner, ctx, seed_override = NULL,
                     progress_local = progress) {
  main_result <- .fit_direction_main_posterior(
    target = target, partner = partner, ctx = ctx,
    seed_override = seed_override, progress_local = progress_local
  )
  .add_predictive_evaluation(main_result)
}


#' Expand pairwise results into directed edge table
#'
#' Produces a long table of directed edges with means, intervals, sign
#' probabilities, and ELPD summaries for \code{j→i} and \code{i→j}.
#'
#' @param df Tibble produced by \code{.run_pair()} over many pairs.
#' @return A tibble with columns \code{from,to,a_mean,a_q2.5,a_q97.5,p_sign2,elpd_*}.
#' @noRd
#' @keywords internal
.mk_cross <- function(df) {
  ij <- df |>
    dplyr::transmute(
      from = .data$j, to = .data$i,
      a_mean = .data$a_ij_mean,
      a_q2.5 = .data$a_ij_q2.5, a_q97.5 = .data$a_ij_q97.5,
      p_sign2 = .data$p_sign2_ij,
      diagnostic_class = .data$diagnostic_class_ij,
      interaction_identifiable = .data$interaction_identifiable_ij,
      residual_identifiable = .data$residual_identifiable_ij,
      elpd_mean = .data$kfold_elpd_mean_ij,
      elpd_sd   = .data$kfold_sd_ij,
      elpd_se   = .data$kfold_se_ij,
      n_subjects = .data$kfold_n_subjects_ij
    )
  ji <- df |>
    dplyr::transmute(
      from = .data$i, to = .data$j,
      a_mean = .data$a_ji_mean,
      a_q2.5 = .data$a_ji_q2.5, a_q97.5 = .data$a_ji_q97.5,
      p_sign2 = .data$p_sign2_ji,
      diagnostic_class = .data$diagnostic_class_ji,
      interaction_identifiable = .data$interaction_identifiable_ji,
      residual_identifiable = .data$residual_identifiable_ji,
      elpd_mean = .data$kfold_elpd_mean_ji,
      elpd_sd   = .data$kfold_sd_ji,
      elpd_se   = .data$kfold_se_ji,
      n_subjects = .data$kfold_n_subjects_ji
    )
  dplyr::bind_rows(ij, ji)
}

#' Extract self-effect summaries
#'
#' Builds a table of \code{a_self_*} summaries and \code{p_sign2_self} by taxon.
#'
#' @param df Tibble produced by \code{.run_pair()} over many pairs.
#' @return A tibble with columns \code{taxon,a_self_mean,a_self_q2.5,a_self_q97.5,p_sign2_self}.
#' @noRd
#' @keywords internal
.mk_self <- function(df) {
  ii <- df |>
    dplyr::transmute(
      taxon = .data$i,
      a_self_mean = .data$a_ii_mean,
      a_self_q2.5 = .data$a_ii_q2.5, a_self_q97.5 = .data$a_ii_q97.5,
      p_sign2_self = .data$p_sign2_ii
    )
  jj <- df |>
    dplyr::transmute(
      taxon = .data$j,
      a_self_mean = .data$a_jj_mean,
      a_self_q2.5 = .data$a_jj_q2.5, a_self_q97.5 = .data$a_jj_q97.5,
      p_sign2_self = .data$p_sign2_jj
    )
  dplyr::bind_rows(ii, jj)
}


#' Build a two-column (subject, value) tibble from named vectors
#'
#' @param keys Character vector of subject IDs.
#' @param vals Numeric vector of values aligned to \code{keys}.
#' @return A tibble with \code{subject} and \code{elpd} (or \code{NULL}).
#' @noRd
#' @keywords internal
.as_named_frame <- function(keys, vals) {
  if (is.null(vals) || length(vals) == 0) return(NULL)
  tibble::tibble(subject = as.character(keys), elpd = as.numeric(vals))
}


#' Build a three-column (subject, elpd, elpd_ppd) tibble
#'
#' @param keys Character vector of subject IDs.
#' @param v1 Numeric vector of ELPD values.
#' @param v2 Optional numeric vector of PPD-normalized ELPD values.
#' @return A tibble with \code{subject}, \code{elpd}, \code{elpd_ppd}.
#' @noRd
#' @keywords internal
.as_named_frame2 <- function(keys, v1, v2) {
  if (is.null(v1) || length(v1) == 0) {
    tibble::tibble(subject = character(0), elpd = numeric(0), elpd_ppd = numeric(0))
  } else {
    tibble::tibble(
      subject = as.character(keys),
      elpd = as.numeric(v1),
      elpd_ppd = if (is.null(v2)) NA_real_ else as.numeric(v2)
    )
  }
}

#' Construct an empty directed subject-level ELPD table
#' @return A zero-row tibble with the stable public column schema and types.
#' @noRd
.empty_cross_subject_elpd <- function() {
  tibble::tibble(
    from = character(),
    to = character(),
    subject = character(),
    elpd = double(),
    elpd_per_observation = double(),
    n_successful_test_observations = integer(),
    elpd_method = character()
  )
}

#' Expand per-pair subject ELPD payloads to long directed format
#'
#' Unnests subject-level full pair-model ELPD for both directions into a long
#' tibble of \code{from → to} rows. It is not self-only predictive evidence.
#'
#' @param df Tibble with nested K-fold payload columns.
#' @return A tibble with the stable subject-level predictive schema.
#' @noRd
#' @keywords internal
.expand_cross_subject_elpd <- function(df) {
  # ij
  ij_rows <- purrr::pmap_dfr(
    list(df$i, df$j, df$kfold_subject_ids_ij, df$kfold_subject_ij, df$kfold_subject_ppd_ij,
         df$kfold_elpd_method_ij, df$kfold_subject_counts_ij),
    function(i, j, ids, elv, elv_ppd, method, cnts) {
      if (is.null(ids) || is.null(elv)) return(NULL)
      tab <- .as_named_frame2(ids, elv, elv_ppd)
      tab$n_successful_test_observations <- as.integer(if (is.null(cnts)) NA else cnts[tab$subject])
      tab$from <- j; tab$to <- i
      tab$elpd_method <- if (is.null(method)) NA_character_ else method
      tab$elpd_per_observation <- tab$elpd_ppd
      tab[, c("from", "to", "subject", "elpd", "elpd_per_observation",
              "n_successful_test_observations", "elpd_method")]
    }
  )
  # ji
  ji_rows <- purrr::pmap_dfr(
    list(df$j, df$i, df$kfold_subject_ids_ji, df$kfold_subject_ji, df$kfold_subject_ppd_ji,
         df$kfold_elpd_method_ji, df$kfold_subject_counts_ji),
    function(j, i, ids, elv, elv_ppd, method, cnts) {
      if (is.null(ids) || is.null(elv)) return(NULL)
      tab <- .as_named_frame2(ids, elv, elv_ppd)
      tab$n_successful_test_observations <- as.integer(if (is.null(cnts)) NA else cnts[tab$subject])
      tab$from <- i; tab$to <- j
      tab$elpd_method <- if (is.null(method)) NA_character_ else method
      tab$elpd_per_observation <- tab$elpd_ppd
      tab[, c("from", "to", "subject", "elpd", "elpd_per_observation",
              "n_successful_test_observations", "elpd_method")]
    }
  )
  out <- dplyr::bind_rows(ij_rows, ji_rows)
  if (!nrow(out)) .empty_cross_subject_elpd() else out
}

#' Assemble the breaking HF2 public predictive-result schema
#' @param res Canonical bidirectional pair-result table.
#' @return A four-element fitted result list.
#' @noRd
.assemble_public_fit_result <- function(res) {
  list(
    cross = .mk_cross(res),
    self = .mk_self(res),
    elpd_subject_cross = .expand_cross_subject_elpd(res),
    raw = res
  )
}
