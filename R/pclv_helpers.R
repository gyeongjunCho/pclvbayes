##### R/glv_helpers.R

#' Get abundance matrix with taxa as rows
#'
#' Coerces the \code{otu_table} from a \pkg{phyloseq} object to a numeric
#' matrix and ensures taxa are on rows (transposes when samples are rows).
#'
#' @param physeq A \pkg{phyloseq} object.
#' @return A numeric matrix with taxa as rows and samples as columns.
#' @noRd
#' @keywords internal
.get_abund_matrix_precomputed <- function(physeq) {
  mat <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) mat <- t(mat)
  mat
}

#' Extract minimal sample metadata (Sample/subject/time)
#'
#' Builds a compact data frame containing sample IDs and user-specified
#' subject and time columns, coercing \code{time} to numeric when possible.
#'
#' @param physeq A \pkg{phyloseq} object.
#' @param subject_col Column name in \code{sample_data} indicating subject ID.
#' @param time_col Column name in \code{sample_data} indicating time.
#' @return A data frame with columns \code{Sample}, \code{subject}, \code{time}.
#' @noRd
#' @keywords internal
.get_sample_meta <- function(physeq, subject_col, time_col) {
  md <- data.frame(
    sample = sample_names(physeq),
    as(sample_data(physeq), "data.frame"),
    check.names = FALSE
  )
  md <- md |>
    dplyr::transmute(
      Sample  = sample,
      subject = .data[[subject_col]],
      time    = .data[[time_col]]
    )
  if (!is.numeric(md$time)) suppressWarnings(md$time <- as.numeric(md$time))
  md
}



#' Per-subject spline smoothing of relative abundances (log-scale)
#'
#' For each taxon and subject, fits a robust smoothing spline on
#' \code{log(pmax(x,0)+eps)} against time, predicts back to the original
#' time points, and returns smoothed abundances on the original scale
#' (with \code{eps} subtraction and non-negativity guard).
#'
#' @param mat_rel Taxa-by-samples relative abundance matrix (rows = taxa).
#' @param meta_df Data frame with columns \code{Sample}, \code{subject}, \code{time}.
#' @param taxa_list Optional character vector of taxa to smooth.
#' @param eps Small constant for log transform stability.
#' @param spline_df Optional effective degrees of freedom for \code{smooth.spline}.
#' @param spline_spar Optional \code{spar} parameter for \code{smooth.spline}.
#' @param spline_cv Logical; use leave-one-out CV (fallbacks to GCV).
#' @param min_unique_times Minimum unique time points required to fit a spline.
#' @return A taxa-by-samples numeric matrix of smoothed abundances.
#' @noRd
#' @keywords internal
.precompute_spline_smoothed <- function(mat_rel, meta_df,
                                        taxa_list = rownames(mat_rel),
                                        eps = 1e-6,
                                        spline_df = NULL,
                                        spline_spar = NULL,
                                        spline_cv = TRUE,
                                        min_unique_times = 3) {
  sm_mat <- matrix(NA_real_, nrow = nrow(mat_rel), ncol = ncol(mat_rel),
                   dimnames = dimnames(mat_rel))
  col_idx_all <- match(meta_df$Sample, colnames(mat_rel))

  taxa_use <- intersect(taxa_list, rownames(mat_rel))
  if (!length(taxa_use)) return(sm_mat)

  for (tx in taxa_use) {
    vec_pred <- rep(NA_real_, nrow(meta_df))

    for (sb in unique(meta_df$subject)) {
      idx   <- which(meta_df$subject == sb)
      times <- meta_df$time[idx]
      cols  <- col_idx_all[idx]
      vals  <- as.numeric(mat_rel[tx, cols])

      ok <- is.finite(times) & is.finite(vals)
      if (!any(ok)) next

      df   <- data.frame(time = times[ok], val = vals[ok])
      df2  <- stats::aggregate(val ~ time, df, mean)
      df2  <- df2[order(df2$time), , drop = FALSE]

      pred_log <- NULL
      if (nrow(df2) >= min_unique_times) {
        ylog <- log(pmax(df2$val, 0) + eps)
        rr <- .smooth_spline_robust(
          x = df2$time, y = ylog,
          spline_df = spline_df,
          spline_spar = spline_spar,
          use_cv = isTRUE(spline_cv),
          min_unique = min_unique_times,
          min_df = 3.0,
          max_df = NULL,
          default_df = 4.0
        )
        pred_log <- try(as.numeric(stats::approx(
          x = df2$time, y = rr$yhat, xout = times, rule = 2
        )$y), silent = TRUE)
        if (inherits(pred_log, "try-error") || !length(pred_log)) pred_log <- rr$yhat
      }

      if (is.null(pred_log)) {
        const_log <- mean(log(pmax(vals, 0) + eps))
        pred_log  <- rep(const_log, length(times))
      }

      vec_pred[idx] <- pmax(exp(pred_log) - eps, 0)
    }

    sm_mat[tx, col_idx_all] <- vec_pred
  }

  sm_mat
}

#' Build smoothed pair table and ΔALR/Δt response
#'
#' Constructs a per-interval table for a directed pair \code{j -> i}:
#' computes ALR for taxon \code{i} versus \code{rest} at consecutive times
#' and returns the per-interval rate \code{(ALR_next - ALR_now)/dt},
#' alongside raw smoothed abundances \code{xi_raw}, \code{xj_raw}.
#'
#' @param sm_mat Smoothed abundance matrix from \code{.precompute_spline_smoothed()}.
#' @param meta_df Metadata with \code{Sample}, \code{subject}, \code{time}.
#' @param j Partner taxon name (from-index).
#' @param i Target taxon name (to-index).
#' @param eps Small constant for ALR stability.
#' @param min_pairs Minimum number of valid intervals to return a table.
#' @param min_dt Minimum positive time difference to accept an interval.
#' @param min_sd Minimum SD threshold for \code{xi_raw}/\code{xj_raw}.
#' @return A tibble with columns \code{subject}, \code{time}, \code{y}, \code{xi_raw}, \code{xj_raw}, or \code{NULL}.
#' @noRd
#' @keywords internal
.build_pair_df_smoothed <- function(sm_mat, meta_df, j, i,
                                    eps = 1e-8,
                                    min_pairs = 4,
                                    min_dt = 1e-8,
                                    min_sd = 1e-12) {
  stopifnot(all(c("Sample","subject","time") %in% names(meta_df)))
  idx_i <- match(meta_df$Sample, colnames(sm_mat))
  xi_all <- as.numeric(sm_mat[i, idx_i])
  xj_all <- as.numeric(sm_mat[j, idx_i])

  df <- meta_df |>
    dplyr::mutate(
      xi_raw = pmax(xi_all, 0),
      xj_raw = pmax(xj_all, 0)
    ) |>
    dplyr::arrange(subject, time)

  pair_df <- df |>
    dplyr::group_by(subject) |>
    dplyr::arrange(time, .by_group = TRUE) |>
    dplyr::mutate(
      time_next = dplyr::lead(time),
      xi_next   = dplyr::lead(xi_raw),
      xj_next   = dplyr::lead(xj_raw),
      dt        = time_next - time
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      is.finite(time), is.finite(time_next),
      is.finite(xi_raw), is.finite(xj_raw),
      is.finite(xi_next), is.finite(xj_next),
      is.finite(dt), dt > min_dt
    ) |>
    dplyr::mutate(
      rest_now  = pmax(0, 1 - xi_raw - xj_raw),
      rest_next = pmax(0, 1 - xi_next - xj_next),
      alr_i_now  = log(pmax(xi_raw,  0) + eps) - log(pmax(rest_now,  0) + eps),
      alr_i_next = log(pmax(xi_next, 0) + eps) - log(pmax(rest_next, 0) + eps),
      y = (alr_i_next - alr_i_now) / dt
    ) |>
    dplyr::transmute(subject, time, y, xi_raw, xj_raw) |>
    dplyr::filter(is.finite(y), is.finite(xi_raw), is.finite(xj_raw))

  if (nrow(pair_df) < min_pairs) return(NULL)
  if (stats::sd(pair_df$xi_raw) < min_sd || stats::sd(pair_df$xj_raw) < min_sd) return(NULL)
  pair_df
}

#' Robust z-score with caps and NA/Inf guards
#'
#' Centers and scales a numeric vector using finite entries only, replaces
#' non-finite results with zero, and caps by \code{[-cap, cap]}.
#'
#' @param v Numeric vector.
#' @param cap Positive cap for absolute z-scores.
#' @return A numeric vector of capped z-scores.
#' @noRd
#' @keywords internal
.z <- function(v, cap = 5){
  v_fin <- v[is.finite(v)]
  m <- if (length(v_fin)) mean(v_fin) else 0
  s <- if (length(v_fin) > 1) sd(v_fin) else 1
  if (!is.finite(s) || s == 0) s <- 1
  z <- (v - m) / s
  z[!is.finite(z)] <- 0
  z[z >  cap] <-  cap
  z[z < -cap] <- -cap
  z
}

#' Safely extract common parameter draws as a data frame
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
      "sigma_ou","lambda","sigma_pred","tau_r","nu"),
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
.summarise_diag <- function(fit,
                            pars = c("a_ij","a_ii","r0","sigma","sd_ou","phi"),
                            max_treedepth = 12) {
  # (1) 파라미터 요약치: R-hat / ESS (가용한 것만 선택)
  all_dd <- fit$draws()
  avail  <- dimnames(all_dd)$variables
  use_pars <- intersect(pars, avail)
  if (!length(use_pars)) use_pars <- avail
  sdtab <- posterior::summarise_draws(fit$draws(use_pars))
  worst_rhat   <- suppressWarnings(max(sdtab$rhat,      na.rm = TRUE))
  min_ess_bulk <- suppressWarnings(min(sdtab$ess_bulk,  na.rm = TRUE))
  min_ess_tail <- suppressWarnings(min(sdtab$ess_tail,  na.rm = TRUE))
  if (!is.finite(worst_rhat))   worst_rhat   <- NA_real_
  if (!is.finite(min_ess_bulk)) min_ess_bulk <- NA_real_
  if (!is.finite(min_ess_tail)) min_ess_tail <- NA_real_

  # (2) 샘플러 진단: divergent__, treedepth__, energy__
  sdiag_df <- posterior::as_draws_df(fit$sampler_diagnostics())

  n_div <- if ("divergent__" %in% names(sdiag_df)) {
    sum(as.integer(sdiag_df[["divergent__"]]), na.rm = TRUE)
  } else NA_integer_

  n_treedepth_hit <- if ("treedepth__" %in% names(sdiag_df)) {
    sum(as.integer(sdiag_df[["treedepth__"]] >= max_treedepth), na.rm = TRUE)
  } else NA_integer_

  # (3) E-BFMI: mean(diff(E)^2) / var(E) (체인별 계산 후 요약)
  ebfmi_min <- ebfmi_med <- NA_real_
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
          ebfmi_min <- val; ebfmi_med <- val
        }
      }
    }
  }

  # (4) 총 저장 draw 수(스칼라)
  n_draws <- posterior::ndraws(fit$draws())

  list(
    worst_rhat      = worst_rhat,
    min_ess_bulk    = min_ess_bulk,
    min_ess_tail    = min_ess_tail,
    n_divergent     = n_div,
    n_treedepth_hit = n_treedepth_hit,
    ebfmi_min       = ebfmi_min,
    ebfmi_med       = ebfmi_med,
    n_draws         = as.integer(n_draws)
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
                    base_pos <- c(if (is.finite(xi) && xi > 0) xi,
                                  if (is.finite(xj) && xj > 0) xj,
                                  if (is.finite(xr) && xr > 0) xr)
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

#' Build model inputs for pairwise gLV regressions
#'
#' Creates lagged predictors and \code{ΔALR_i/Δt} response under either
#' ALR or raw-RA transforms, with zero-aware ALR safeguards, optional
#' ALR smoothing, partner non-zero filters, and subject/global z-scaling.
#'
#' @param pair_df Tibble from \code{.build_pair_df_smoothed()}.
#' @param transform Either \code{"alr"} or \code{"raw"} for predictors.
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
#' @param standardize_by_subject Logical; subject-wise z-standardization.
#' @param z_mode One of \code{"subject"}, \code{"global"}, or \code{"none"}.
#' @param z_external Optional list of global means/SDs when \code{z_mode="global"}.
#' @return Data frame with columns \code{subject,time,y,xi,xj} and attributes
#'   \code{smooth_edf_mean}, \code{smooth_scale}, \code{smoothed}, optionally \code{z_stats}.
#' @noRd
#' @keywords internal
.make_pair_inputs_glv <- function(pair_df,
                                  transform = c("alr","raw"),
                                  lag = 1,
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
                                  nz_partner_min_frac = 0.15,
                                  standardize_by_subject = TRUE,
                                  z_mode = c("subject","global","none"),
                                  z_external = NULL) {
  transform     <- match.arg(transform)
  zero_mode_alr <- match.arg(zero_mode_alr)
  minpos_base   <- match.arg(minpos_base)
  smooth_scale  <- match.arg(smooth_scale)
  z_mode        <- match.arg(z_mode)

  stopifnot(all(c("subject","time","xi_raw","xj_raw") %in% names(pair_df)))
  df <- pair_df[order(pair_df$subject, pair_df$time), , drop = FALSE]

  df$xi_raw   <- pmax(as.numeric(df$xi_raw), 0)
  df$xj_raw   <- pmax(as.numeric(df$xj_raw), 0)
  df$rest_raw <- pmax(0, 1 - df$xi_raw - df$xj_raw)

  by_s <- split(seq_len(nrow(df)), df$subject)
  lagv <- function(v) unsplit(lapply(by_s, function(ix) dplyr::lag(v[ix], lag)), df$subject)
  diff_over_dt <- function(v, t) {
    unsplit(lapply(by_s, function(ix) {
      vi <- v[ix]; ti <- t[ix]
      vi_lag <- dplyr::lag(vi, lag)
      dt     <- as.numeric(ti - dplyr::lag(ti, lag))
      # 하한: subject 내 양의 dt들의 중앙값의 25% (예)
      dt_pos <- dt[is.finite(dt) & dt > 0]
      dt_min <- if (length(dt_pos)) 0.25 * stats::median(dt_pos) else 1
      dt_adj <- pmax(dt, dt_min)
      out <- (vi - vi_lag) / dt_adj
      out[!is.finite(dt) | dt <= 0] <- NA_real_
      out
    }), df$subject)
  }

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
  alr_i <- log(trip[, "xi"]) - log(trip[, "xr"])
  alr_j <- log(trip[, "xj"]) - log(trip[, "xr"])

  if (is.finite(alr_cap)) {
    alr_i <- pmax(pmin(alr_i,  alr_cap), -alr_cap)
    alr_j <- pmax(pmin(alr_j,  alr_cap), -alr_cap)
  } else {
    es <- pmax(trip[, "eps_star"], .Machine$double.eps)
    num <- pmax(1 - 2 * es, .Machine$double.eps)
    B_theory <- log(num / es)
    cap_vec  <- pmin(8, pmax(0, B_theory))
    alr_i <- pmax(pmin(alr_i,  cap_vec), -cap_vec)
    alr_j <- pmax(pmin(alr_j,  cap_vec), -cap_vec)
  }

  smooth_edf_mean <- NA_real_
  if (smooth_scale == "alr") {
    smooth_one <- function(v, t) {
      ok <- is.finite(v) & is.finite(t)
      if (sum(ok) < 3L || length(unique(t[ok])) < 3L) {
        return(list(y = v, df = NA_real_))
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
      yhat <- rep(NA_real_, length(t))
      yhat[ok] <- as.numeric(stats::approx(
        x = t[ok], y = rr$yhat, xout = t[ok], rule = 2
      )$y)
      list(y = yhat, df = rr$df)
    }
    edf_i <- edf_j <- rep(NA_real_, length(by_s))
    names(edf_i) <- names(edf_j) <- names(by_s)
    for (sname in names(by_s)) {
      ix <- by_s[[sname]]
      res_i <- smooth_one(alr_i[ix], df$time[ix])
      res_j <- smooth_one(alr_j[ix], df$time[ix])
      alr_i[ix] <- res_i$y
      alr_j[ix] <- res_j$y
      edf_i[sname] <- res_i$df
      edf_j[sname] <- res_j$df
    }
    edf_subj <- rowMeans(cbind(edf_i, edf_j), na.rm = TRUE)
    smooth_edf_mean <- mean(edf_subj, na.rm = TRUE)
    if (!is.finite(smooth_edf_mean)) smooth_edf_mean <- NA_real_
  }

  if (transform == "alr") {
    xi <- lagv(alr_i)
    xj <- lagv(alr_j)
    y  <- diff_over_dt(alr_i, df$time)
  } else {
    xi <- lagv(df$xi_raw)
    xj <- lagv(df$xj_raw)
    y  <- diff_over_dt(alr_i, df$time)  # response stays ΔALR_i/Δt
  }

  # 파트너 희소성 필터(옵션)
  if (is.finite(nz_partner_min_frac) && nz_partner_min_frac > 0) {
    by_s2 <- split(seq_len(nrow(df)), df$subject)
    keep_mask <- rep(TRUE, nrow(df))
    for (sname in names(by_s2)) {
      ix <- by_s2[[sname]]
      frac_nz <- mean(is.finite(xj[ix]) & xj[ix] != 0, na.rm = TRUE)
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

  z_stats <- NULL
  if (z_mode == "none") {
    # no standardization
  } else if (!is.null(z_external) && z_mode == "global") {
    mu_y  <- z_external$mu_y;  sd_y  <- z_external$sd_y
    mu_xi <- z_external$mu_xi; sd_xi <- z_external$sd_xi
    mu_xj <- z_external$mu_xj; sd_xj <- z_external$sd_xj
    if (!is.finite(sd_y)  || sd_y  <= 0) sd_y  <- 1
    if (!is.finite(sd_xi) || sd_xi <= 0) sd_xi <- 1
    if (!is.finite(sd_xj) || sd_xj <= 0) sd_xj <- 1
    dat$y  <- (dat$y  - mu_y)  / sd_y
    dat$xi <- (dat$xi - mu_xi) / sd_xi
    dat$xj <- (dat$xj - mu_xj) / sd_xj
  } else {
    if (isTRUE(standardize_by_subject) && z_mode == "subject") {
      by_s2 <- split(seq_len(nrow(dat)), dat$subject)
      dat$y  <- unsplit(lapply(by_s2, function(ix) .z(dat$y[ix])),  dat$subject)
      dat$xi <- unsplit(lapply(by_s2, function(ix) .z(dat$xi[ix])), dat$subject)
      dat$xj <- unsplit(lapply(by_s2, function(ix) .z(dat$xj[ix])), dat$subject)
    } else {
      mu_y  <- mean(dat$y[is.finite(dat$y)],  na.rm = TRUE); sd_y  <- stats::sd(dat$y,  na.rm = TRUE); if (!is.finite(sd_y)  || sd_y  == 0) sd_y  <- 1
      mu_xi <- mean(dat$xi[is.finite(dat$xi)], na.rm = TRUE); sd_xi <- stats::sd(dat$xi, na.rm = TRUE); if (!is.finite(sd_xi) || sd_xi == 0) sd_xi <- 1
      mu_xj <- mean(dat$xj[is.finite(dat$xj)], na.rm = TRUE); sd_xj <- stats::sd(dat$xj, na.rm = TRUE); if (!is.finite(sd_xj) || sd_xj == 0) sd_xj <- 1
      dat$y  <- (dat$y  - mu_y)  / sd_y
      dat$xi <- (dat$xi - mu_xi) / sd_xi
      dat$xj <- (dat$xj - mu_xj) / sd_xj
      z_stats <- list(mu_y = mu_y, sd_y = sd_y,
                      mu_xi = mu_xi, sd_xi = sd_xi,
                      mu_xj = mu_xj, sd_xj = sd_xj)
    }
  }

  ok <- is.finite(dat$y) & is.finite(dat$xi) & is.finite(dat$xj) & is.finite(dat$time) & keep_mask
  dat <- dat[ok, , drop = FALSE]
  if (!nrow(dat)) return(NULL)

  attr(dat, "smooth_edf_mean") <- smooth_edf_mean
  attr(dat, "smooth_scale")    <- smooth_scale
  attr(dat, "smoothed")        <- smooth_scale %in% c("alr","logra")
  if (!is.null(z_stats)) attr(dat, "z_stats") <- z_stats

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
  # 1) silent=FALSE면 그대로 호출
  if (!isTRUE(silent)) return(do.call(mod$sample, args))

  # 2) 콘솔 무음 + 진행/메시지 억제
  args$refresh <- 0L
  args$show_messages <- FALSE

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

  args <- .ensure_outputs(args)

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
  out_con <- file(nullfile(), open = "wt")
  msg_con <- file(nullfile(), open = "wt")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(sink(), silent = TRUE)
    try(close(msg_con), silent = TRUE)
    try(close(out_con), silent = TRUE)
  }, add = TRUE)
  sink(out_con); sink(msg_con, type = "message")

  # 4) 1차 시도
  fit <- try(suppressWarnings(suppressMessages(do.call(mod$sample, args))), silent = TRUE)

  # 5) 실패 시 재시도: (a) init 관련 오류면 폴백, (b) 새 디렉터리/베이스네임 + 병렬체인 1
  if (inherits(fit, "try-error")) {
    msg <- tryCatch(as.character(attr(fit, "condition")), error = function(e) "")
    # init 관련 에러면 init=NULL로 폴백
    if (length(msg) && any(grepl("'init' contains empty lists|invalid init|bad init", msg, ignore.case = TRUE))) {
      message("Pathfinder/init not accepted; falling back to default init.")
      args$init <- NULL
    }
    # 재시도 인자: 새 디렉터리/베이스네임 + parallel_chains=1
    args_retry <- args
    args_retry$output_dir      <- file.path(output_root, sprintf("run_retry_%s_%s", .now_tag(), .make_uid(6)))
    dir.create(args_retry$output_dir, recursive = TRUE, showWarnings = FALSE)
    args_retry$output_basename <- sprintf("glv_pairwise_retry_%s", .make_uid(6))
    if (!is.null(args_retry$chains)) args_retry$parallel_chains <- 1L

    # 초미니면 재시도 전에도 짧은 지터 한 번 더
    if (is.numeric(args_retry$iter_warmup) && is.numeric(args_retry$iter_sampling)) {
      if (args_retry$iter_warmup <= 5L && args_retry$iter_sampling <= 5L) {
        Sys.sleep(runif(1, 0, 0.25))
      }
    }

    fit <- suppressWarnings(suppressMessages(do.call(mod$sample, args_retry)))
  }

  fit
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
#' If no usable values are found, returns NULL (caller should fallback).
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

#' Make init from PF or fallback to numeric scalar (never returns empty lists)
#' @noRd
.init_from_pf_or_scalar <- function(pf_fit, mod, chains, fallback_numeric = 0.2, verbose = FALSE, tag = NULL) {
  say <- function(...) if (verbose) cat(sprintf("%s %s\n", if (!is.null(tag)) tag else "", sprintf(...)))
  pf_inits <- NULL
  try({ pf_inits <- .pf_inits_from_draws(pf_fit, mod, chains = chains) }, silent = TRUE)
  if (is.null(pf_inits)) {
    say("Pathfinder draws unusable → fallback numeric init=%.3f", fallback_numeric)
    return(fallback_numeric)
  }
  # 최종 검증: 체인 수/빈리스트/이름 보유
  if (!is.list(pf_inits) || length(pf_inits) != chains) {
    say("PF init shape invalid → fallback numeric init=%.3f", fallback_numeric)
    return(fallback_numeric)
  }
  good <- TRUE
  for (k in seq_len(chains)) {
    li <- pf_inits[[k]]
    if (!is.list(li) || length(li) == 0L || !length(names(li))) { good <- FALSE; break }
    # 모든 값 numeric scalar & finite
    if (!all(vapply(li, function(v) is.numeric(v) && length(v) == 1L && is.finite(v), logical(1)))) {
      good <- FALSE; break
    }
  }
  if (!good) {
    say("PF init contains empty/invalid lists → fallback numeric init=%.3f", fallback_numeric)
    return(fallback_numeric)
  }
  pf_inits
}

#' Robust wrapper for \code{smooth.spline()} with safe fallbacks
#'
#' Sorts and aggregates duplicate \code{x}, tries \code{df}/\code{spar}/CV/GCV
#' strategies in order, and falls back to constant fits when needed.
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
                                  spline_df = NULL,
                                  spline_spar = NULL,
                                  use_cv = TRUE,
                                  min_unique = 3,
                                  min_df = 3.0,
                                  max_df = NULL,
                                  default_df = 4.0) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (!length(x)) return(list(yhat = y, df = NA_real_))

  # 정렬 + 중복 x 평균 집계(LOOCV 안정화)
  ord <- order(x); x <- x[ord]; y <- y[ord]
  df2 <- stats::aggregate(y ~ x, data.frame(x = x, y = y), mean)
  x_u <- as.numeric(df2$x); y_u <- as.numeric(df2$y)

  # 유니크 시점 부족 → 상수
  if (length(x_u) < min_unique) {
    yhat <- rep(stats::mean(y_u), length(x))
    return(list(yhat = yhat[order(ord)], df = NA_real_))
  }

  # df 상한 자동
  if (is.null(max_df)) {
    max_df <- max(min(length(x_u) - 1L, 10L), min_df)
  }

  try_fit <- function(mode = c("df","spar","cv","gcv")) {
    mode <- match.arg(mode)
    if (mode == "df" && !is.null(spline_df)) {
      df_use <- min(max(spline_df, min_df), max_df)
      return(suppressWarnings(stats::smooth.spline(x_u, y_u, df = df_use)))
    }
    if (mode == "spar" && !is.null(spline_spar)) {
      return(suppressWarnings(stats::smooth.spline(x_u, y_u, spar = spline_spar)))
    }
    if (mode == "cv") {
      return(suppressWarnings(stats::smooth.spline(x_u, y_u, cv = TRUE, all.knots = TRUE)))
    }
    return(suppressWarnings(stats::smooth.spline(x_u, y_u, cv = FALSE, all.knots = TRUE)))
  }

  fit <- NULL
  for (m in c(if (!is.null(spline_df)) "df" else NULL,
              if (!is.null(spline_spar)) "spar" else NULL,
              if (isTRUE(use_cv)) "cv" else NULL,
              "gcv")) {
    fit <- try(try_fit(m), silent = TRUE)
    if (!inherits(fit, "try-error")) break
  }

  # 실패 → 기본 df
  if (inherits(fit, "try-error") || is.null(fit)) {
    fit <- try(suppressWarnings(stats::smooth.spline(
      x_u, y_u, df = min(max(default_df, min_df), max_df)
    )), silent = TRUE)
  }

  # 그래도 실패 → 상수
  if (inherits(fit, "try-error") || is.null(fit)) {
    yhat <- rep(stats::mean(y_u), length(x))
    return(list(yhat = yhat[order(ord)], df = NA_real_))
  }

  # 예측 (외삽 가드)
  yhat <- try(suppressWarnings(as.numeric(stats::predict(fit, x = x)$y)), silent = TRUE)
  if (inherits(yhat, "try-error") || !length(yhat) || any(!is.finite(yhat))) {
    yhat <- suppressWarnings(as.numeric(stats::approx(
      x = x_u,
      y = as.numeric(stats::predict(fit, x = x_u)$y),
      xout = x, rule = 2
    )$y))
  }

  df_out <- suppressWarnings(as.numeric(fit$df))
  list(yhat = yhat[order(ord)], df = ifelse(is.finite(df_out), df_out, NA_real_))
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
#' (e.g., bumping \code{nu}, higher \code{adapt_delta}, \code{dense_e}) until
#' diagnostics pass or \code{max_retries} is reached.
#'
#' @param mod A \pkg{cmdstanr} model.
#' @param base_args Baseline argument list for \code{mod$sample()}.
#' @param stan_list Data list; inserted into \code{base_args$data}.
#' @param max_retries Maximum retry attempts.
#' @param nu_fix Base fixed-\eqn{\nu} for t-noise; may be bumped up to \code{max_nu_cap}.
#' @param silent_sampler Logical; silence sampler output.
#' @param ebfmi_thresh E-BFMI threshold triggering retries.
#' @param max_nu_cap Maximum \eqn{\nu} allowed.
#' @param tag Optional label for progress messages.
#' @param freeze_retry_hypers If \code{TRUE}, do not alter hyperparameters on retries.
#' @return A list with \code{fit}, \code{diag}, \code{n_retries}, \code{used_nu},
#'   \code{fit_failed}, and \code{final_args}.
#' @noRd
#' @keywords internal
.sample_with_retry <- function(mod, base_args, stan_list,
                               max_retries = 3,
                               nu_fix,
                               silent_sampler = FALSE,
                               ebfmi_thresh = 0.30,
                               max_nu_cap = 7,
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

      lst
    }
  }
  # 재시도 필요 여부 판단 (EBFMI + div + treedepth)
  .needs_retry <- function(diag, fit = NULL, thr = 0.30) {
    low_eb <- is.finite(diag$ebfmi_min) && (diag$ebfmi_min < thr)
    if (!low_eb && !is.null(fit)) {
      eb_vec <- .ebfmi_chainwise_from_energy(fit)
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

  # 준비
  base_args <- base_args
  base_args$data <- stan_list
  nu_fix <- .or(nu_fix, .or(stan_list$nu_fixed, 4L))

  attempt <- 0L
  used_nu <- nu_fix
  fit <- NULL; diag <- NULL; fit_failed <- FALSE
  final_args <- NULL

  repeat {
    # nu bump (freeze면 bump 안 함)
    if (attempt == 0L) {
      used_nu <- min(max_nu_cap, nu_fix)
    } else if (!freeze_retry_hypers) {
      used_nu <- min(max_nu_cap, nu_fix + min(attempt, 2L))  # +1, +2까지
    } else {
      used_nu <- nu_fix
    }
    base_args$data$nu_fixed <- used_nu

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
      sample_args$seed <- base_args$seed
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
    message(sprintf("⏩ %sattempt %d/%d | nu_fixed=%d | metric=%s, adapt_delta=%s, warmup=%s, step_size=%s",
                    if (nzchar(tag)) paste0("[",tag,"] ") else "",
                    attempt, max_retries, used_nu,
                    .or(sample_args$metric, "NA"),
                    ifelse(is.null(sample_args$adapt_delta),"NA",sprintf("%.3f", sample_args$adapt_delta)),
                    .or(sample_args$iter_warmup, "NA"),
                    .or(sample_args$step_size, "NA")))

    fit <- .call_sample_silently(mod, sample_args, silent = silent_sampler)

    diag <- .summarise_diag(
      fit,
      max_treedepth = .or(sample_args$max_treedepth, base_args$max_treedepth)
    )

    final_args <- sample_args

    eb_vec <- .ebfmi_chainwise_from_energy(fit)
    eb_str <- if (length(eb_vec)) paste(sprintf("%.3f", eb_vec), collapse=",") else "NA"
    ok <- !.needs_retry(diag, fit, thr = ebfmi_thresh)

    message(sprintf("%s%s | %s | E-BFMI chains=[%s]",
                    if (nzchar(tag)) paste0("[",tag,"] ") else "",
                    if (ok) "✅ diag ok" else "⚠️ diag warn",
                    .fmt_diag(diag), eb_str))

    if (ok || attempt >= max_retries) { if (!ok) fit_failed <- TRUE; break }
    attempt <- attempt + 1L
  }

  list(
    fit = fit,
    diag = diag,
    n_retries = attempt,
    used_nu = used_nu,
    fit_failed = fit_failed,
    final_args = final_args
  )
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

  seed_ij <- as.integer(seed_base + 100000L * idx_i + 1000L * idx_j + 1L)
  seed_ji <- as.integer(seed_base + 100000L * idx_i + 1000L * idx_j + 2L)

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

  if (is.null(res_ij) || is.null(res_ji)) return(NULL)

  tibble::as_tibble_row(list(
    i = ti, j = pj,
    n_pairs_ij = res_ij$n_pairs, n_pairs_ji = res_ji$n_pairs,
    a_ij_mean = res_ij$a_mean, a_ij_sd = res_ij$a_sd,
    a_ij_q2.5 = res_ij$a_q2.5, a_ij_q97.5 = res_ij$a_q97.5,
    p_sign2_ij = res_ij$p_sign2,
    a_ji_mean = res_ji$a_mean, a_ji_sd = res_ji$a_sd,
    a_ji_q2.5 = res_ji$a_q2.5, a_ji_q97.5 = res_ji$a_q97.5,
    p_sign2_ji = res_ji$p_sign2,
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
    kfold_K = kfold_K, kfold_R = kfold_R,
    kfold_subject_ij         = res_ij$kfold_subject,
    kfold_subject_ppd_ij     = res_ij$kfold_subject_ppd,
    kfold_subject_ids_ij     = res_ij$kfold_subject_ids,
    kfold_subject_counts_ij  = res_ij$kfold_subject_counts,
    kfold_splits_ij          = res_ij$kfold_splits,
    kfold_seed_used_ij       = res_ij$kfold_seed_used,
    kfold_sd_ij              = res_ij$kfold_sd,
    kfold_se_ij              = res_ij$kfold_se,
    kfold_n_subjects_ij      = res_ij$kfold_n_subjects,
    kfold_retry_total_ij     = res_ij$kfold_retry_total,
    kfold_retry_mean_ij      = res_ij$kfold_retry_mean,
    kfold_nu_used_counts_ij  = res_ij$kfold_nu_used_counts,
    kfold_outer_rounds_ij    = res_ij$kfold_outer_rounds,
    kfold_failed_ij          = res_ij$kfold_failed,
    kfold_folds_ok_ij        = res_ij$kfold_n_folds_ok,
    kfold_folds_fail_ij      = res_ij$kfold_n_folds_fail,
    kfold_subject_ji         = res_ji$kfold_subject,
    kfold_subject_ppd_ji     = res_ji$kfold_subject_ppd,
    kfold_subject_ids_ji     = res_ji$kfold_subject_ids,
    kfold_subject_counts_ji  = res_ji$kfold_subject_counts,
    kfold_splits_ji          = res_ji$kfold_splits,
    kfold_seed_used_ji       = res_ji$kfold_seed_used,
    kfold_sd_ji              = res_ji$kfold_sd,
    kfold_se_ji              = res_ji$kfold_se,
    kfold_n_subjects_ji      = res_ji$kfold_n_subjects,
    kfold_retry_total_ji     = res_ji$kfold_retry_total,
    kfold_retry_mean_ji      = res_ji$kfold_retry_mean,
    kfold_nu_used_counts_ji  = res_ji$kfold_nu_used_counts,
    kfold_outer_rounds_ji    = res_ji$kfold_outer_rounds,
    kfold_failed_ji          = res_ji$kfold_failed,
    kfold_folds_ok_ji        = res_ji$kfold_n_folds_ok,
    kfold_folds_fail_ji      = res_ji$kfold_n_folds_fail
  ))
}

#' Subject-level projection log-likelihood (indep or Kalman-OU)
#'
#' Computes per-subject log-likelihood of held-out data given posterior draws,
#' either under an independent Gaussian/t-model or via a Kalman filter for OU residuals.
#'
#' @param draws_df Posterior draws data frame (at least \code{r0,a_ii,a_ij} and noise params).
#' @param pair_in Test data with columns \code{y,xi,xj,subject,time}.
#' @param resid_mode Residual mode, e.g. \code{"ou"} or \code{"wn"}.
#' @param use_t Logical; enable Student-\eqn{t} scoring if \eqn{\nu} available.
#' @param fit_full Optional full fit to extract \eqn{\nu} when needed.
#' @param nu_scalar Optional scalar \eqn{\nu} to use if draws lack \eqn{\nu}.
#' @param elpd_mode One of \code{"indep"} or \code{"kalman"} (default depends on \code{resid_mode}).
#' @return A list with matrix \code{full} (draws × subjects), \code{subjects}, and \code{n_obs}.
#' @noRd
#' @keywords internal
.proj_loglik_subject <- function(draws_df,
                                 pair_in,
                                 resid_mode,
                                 use_t = FALSE,
                                 fit_full = NULL,
                                 nu_scalar = NULL,
                                 elpd_mode = NULL) {
  # 자동 디폴트: OU면 kalman, 아니면 indep
  if (is.null(elpd_mode)) {
    elpd_mode <- if (identical(resid_mode, "ou")) "kalman" else "indep"
  } else {
    elpd_mode <- match.arg(elpd_mode, c("indep","kalman"))
  }
  stopifnot(all(c("y", "xi", "xj", "subject") %in% names(pair_in)))
  subs <- unique(pair_in$subject)
  D    <- nrow(draws_df)
  Ssub <- length(subs)

  # 불규칙 간격 dt/subject-wise 인덱스 구성
  mk_prev_dt <- function(df) {
    df <- df[order(df$subject, df$time), , drop = FALSE]
    by_s <- split(seq_len(nrow(df)), df$subject)
    dtv <- numeric(nrow(df)); dtv[] <- NA_real_
    for (sb in names(by_s)) {
      ix <- by_s[[sb]]; dtv[ix[1]] <- NA_real_
      if (length(ix) >= 2L) for (k in 2:length(ix)) {
        dtv[ix[k]] <- as.numeric(df$time[ix[k]] - df$time[ix[k-1]])
      }
    }
    dtv
  }
  dt_all <- mk_prev_dt(pair_in)

  # 우선순위: (1) sigma_pred 생성량 → (2) sigma & sd_ou → (3) sigma, sigma_ou, lambda → (4) sigma
  if (resid_mode == "ou") {
    if ("sigma_pred" %in% names(draws_df)) {
      sigma_pred <- draws_df$sigma_pred
    } else if (all(c("sigma", "sd_ou") %in% names(draws_df))) {
      sigma_pred <- sqrt(draws_df$sigma^2 + draws_df$sd_ou^2)
    } else if (all(c("sigma", "sigma_ou", "lambda") %in% names(draws_df))) {
      sigma_pred <- sqrt(draws_df$sigma^2 + draws_df$sigma_ou^2 / (2 * pmax(draws_df$lambda, 1e-8)))
    } else {
      sigma_pred <- draws_df$sigma
    }
  } else {
    # "wn"
    sigma_pred <- draws_df$sigma
  }

  nu_vec <- NULL
  if (isTRUE(use_t)) {
    if ("nu" %in% names(draws_df)) {
      nu_vec <- draws_df$nu                     # 1) prefer estimated df (if modeled)
    } else if (!is.null(nu_scalar) && is.finite(nu_scalar)) {
      nu_vec <- rep(nu_scalar, nrow(draws_df))  # 2) fallback to fixed df from wrapper
    } else if (!is.null(fit_full)) {
      nu_draws <- try(fit_full$draws("nu"), silent = TRUE)
      if (!inherits(nu_draws, "try-error") &&
          !is.null(nu_draws)) {
        nu_vec <- as.numeric(posterior::as_draws_df(nu_draws)[["nu"]])
      }
    }
  }
  # If df unavailable, fall back to Normal scoring downstream
  use_t <- isTRUE(use_t) && !is.null(nu_vec)

  loglik_full <- matrix(NA_real_, nrow = D, ncol = Ssub)
  n_obs_vec   <- integer(Ssub)

  # --- 칼만필터 기반 OU 주변우도 (resid_mode="ou" & elpd_mode="kalman") ---
  .ou_kf_ll_one_subject <- function(draws_df, y, xi, xj, tvec) {
    J <- length(y)
    if (J == 0) return(rep(NA_real_, nrow(draws_df)))
    # mean: mu = r0 + a_ii * xi + a_ij * xj
    r0_mat <- matrix(draws_df$r0, nrow = nrow(draws_df), ncol = J)
    aii_xi <- tcrossprod(draws_df$a_ii, xi)
    aij_xj <- tcrossprod(draws_df$a_ij, xj)
    mu     <- r0_mat + aii_xi + aij_xj     # (D x J)

    y_rep  <- matrix(rep(y, each = nrow(draws_df)), nrow = nrow(draws_df))

    # 측정잡음 분산 R (t면 ν/(ν-2) 팩터로 가우시안 근사)
    if (isTRUE(use_t)) {
      t_fac <- ifelse(nu_vec > 2, nu_vec/(nu_vec - 2), 5.0)
      R <- (draws_df$sigma^2) * t_fac
    } else {
      R <- (draws_df$sigma^2)
    }
    R <- pmax(R, 1e-10)

    # OU 이산화: a_t, q_t
    have_lambda <- "lambda" %in% names(draws_df)
    have_phi    <- "phi"    %in% names(draws_df)
    have_sdou   <- "sd_ou"  %in% names(draws_df)

    # stationary var of OU
    if (have_sdou) {
      sd2 <- (draws_df$sd_ou)^2
    } else if (all(c("sigma_ou","lambda") %in% names(draws_df))) {
      sd2 <- (draws_df$sigma_ou^2) / (2 * pmax(draws_df$lambda, 1e-8))
    } else {
      # 정보 부족 → 독립 근사 fallback
      return(matrixStats::rowSums2(
        dnorm(y_rep, mean = mu, sd = matrix(sigma_pred, nrow = nrow(draws_df), ncol = J), log = TRUE)
      ))
    }

    # 필터 초기화: e0 ~ N(0, sd2)
    D <- nrow(draws_df)
    ll <- rep(0.0, D)

    ll <- rep(0.0, D)

    use_ri <- ("tau_r" %in% names(draws_df)) && any(is.finite(draws_df$tau_r))
    if (!use_ri) {
      m  <- rep(0.0, D)     # e-state mean per draw
      P  <- sd2             # e-state var per draw
    } else {
      # 2D state: [e_t, r_s]^T, with r_s time-invariant
      m_e <- rep(0.0, D); m_r <- rep(0.0, D)
      P_ee <- sd2
      P_rr <- pmax(draws_df$tau_r^2, 0)
      P_er <- rep(0.0, D)
    }

    # 시간 루프
    # 주의: 첫 관측의 dt가 NA일 수 있으므로 dt<=0 → 작은 값으로 보정
    dtv <- as.numeric(c(NA, diff(tvec)))
    l2pi <- log(2*pi)
    for (t in seq_len(J)) {
      dt <- dtv[t]
      if (!is.finite(dt) || dt < 0) dt <- 0        # 첫 관측은 dt=0로 처리
      if (have_lambda) {
        lam <- pmax(draws_df$lambda, 1e-8)
        a_t <- exp(-lam * dt)
      } else if (have_phi) {
        phi <- pmin(pmax(draws_df$phi, 1e-8), 0.999999)
        lam <- -log(phi)
        a_t <- exp(-lam * dt)
      } else {
        a_t <- rep(0.0, D)
      }
      q_t <- pmax(sd2 * (1 - a_t^2), 1e-16)
      if (!use_ri) {
        # ---- 1D (OU only) ----
        # 예측
        m_pred <- a_t * m
        P_pred <- (a_t^2) * P + q_t
        # 혁신
        v_t <- y_rep[, t] - mu[, t] - m_pred
        F_t <- pmax(P_pred + R, 1e-10)
        # loglik
        ll <- ll + (-0.5 * (l2pi + log(F_t) + (v_t * v_t) / F_t))
        # 업데이트
        K_t <- P_pred / F_t
        m   <- m_pred + K_t * v_t
        P   <- pmax((1 - K_t) * P_pred, 1e-12)
      } else {
        # ---- 2D (OU + random intercept) ----
        # 예측
        m_e_pred <- a_t * m_e
        m_r_pred <- m_r                   # r_s is time-invariant
        P_ee_p <- (a_t^2) * P_ee + q_t
        P_er_p <- a_t * P_er
        P_rr_p <- P_rr                    # no process noise on r_s
        # 혁신
        v_t <- y_rep[, t] - mu[, t] - (m_e_pred + m_r_pred)
        F_t <- pmax(P_ee_p + P_rr_p + 2*P_er_p + R, 1e-10)  # H=[1,1]
        # loglik
        ll <- ll + (-0.5 * (l2pi + log(F_t) + (v_t * v_t) / F_t))
        # 칼만 이득 (2x1)
        K_e <- (P_ee_p + P_er_p) / F_t
        K_r <- (P_er_p + P_rr_p) / F_t
        # 업데이트
        m_e <- m_e_pred + K_e * v_t
        m_r <- m_r_pred + K_r * v_t
        # 공분산 업데이트: P = (I-KH)P_pred
        # 요소별 전개 (H=[1,1])
        P_ee <- pmax(P_ee_p - K_e*(P_ee_p + P_er_p), 1e-12)
        P_rr <- pmax(P_rr_p - K_r*(P_er_p + P_rr_p), 1e-12)
        P_er <-       P_er_p - K_e*(P_er_p + P_rr_p)
      }
    }
    ll
  }

  for (s_idx in seq_along(subs)) {
    sb <- subs[s_idx]
    ix <- which(pair_in$subject == sb)
    y  <- pair_in$y[ix]
    xi <- pair_in$xi[ix]
    xj <- pair_in$xj[ix]
    tvec <- pair_in$time[ix]
    J  <- length(ix)
    n_obs_vec[s_idx] <- J

    if (identical(resid_mode, "ou") && identical(elpd_mode, "kalman")) {
      ll_f <- .ou_kf_ll_one_subject(draws_df, y, xi, xj, tvec)
    } else {
      # 독립 근사(레거시)
      r0_mat  <- matrix(draws_df$r0, nrow = D, ncol = J)
      aii_xi  <- tcrossprod(draws_df$a_ii, xi)
      aij_xj  <- tcrossprod(draws_df$a_ij, xj)
      mu_full <- r0_mat + aii_xi + aij_xj
      sig_mat <- matrix(sigma_pred, nrow = D, ncol = J)
      y_mat   <- matrix(rep(y, each = D), nrow = D, ncol = J)
      if (!is.null(nu_vec)) {
        nu_mat <- matrix(nu_vec, nrow = D, ncol = J)
        ll_f <- rowSums(stats::dt((y_mat - mu_full) / sig_mat,
                                  df = nu_mat,
                                  log = TRUE
        ) - log(sig_mat))
      } else {
        ll_f <- rowSums(stats::dnorm(
          y_mat,
          mean = mu_full,
          sd = sig_mat,
          log = TRUE
        ))
      }
    }
    loglik_full[, s_idx] <- ll_f
  }
  list(full = loglik_full, subjects = subs, n_obs = n_obs_vec)
}

#' Make repeated K-fold splits at the subject level
#'
#' Randomly partitions unique subjects into \code{K} folds, repeated \code{R} times.
#'
#' @param subject_vec Subject IDs aligned to rows of the data.
#' @param K Number of folds.
#' @param R Number of repetitions.
#' @param seed RNG seed.
#' @return A nested list \code{[[r]][[k]]} with \code{train_subjects}/\code{test_subjects}.
#' @noRd
#' @keywords internal
.make_repkfold_splits <- function(subject_vec,
                                  K = 5,
                                  R = 3,
                                  seed = 123) {
  set.seed(seed)
  subs <- unique(subject_vec)
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
#' per subject; optionally applies train-only global scaling to avoid leakage.
#'
#' @param pair_in Full input from \code{.make_pair_inputs_glv()}.
#' @param train_subjects,test_subjects Character vectors of subject IDs.
#' @param use_global_scaling Logical; apply scaling based on train only.
#' @return A list with \code{train}, \code{test}, and index vectors; or \code{NULL}.
#' @noRd
#' @keywords internal
.build_train_test <- function(pair_in,
                              train_subjects,
                              test_subjects,
                              use_global_scaling = TRUE,
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
          dtv[ix[k]] <- if (is.finite(dt_k) &&
                            dt_k > 0)
            dt_k
          else
            1e-6
        }
      }
    }
    list(df = df,
         prev = prev,
         dt = dtv)
  }
  tr <- subset(pair_in, subject %in% train_subjects)
  te <- subset(pair_in, subject %in% test_subjects)
  if (!nrow(tr) || !nrow(te))
    return(NULL)
  trd <- mk_prev_dt(tr)
  ted <- mk_prev_dt(te)
  # --- ★ Train-only global scaling (OOF leakage guard for K-fold) ---
  if (isTRUE(use_global_scaling)) {
    sc <- list(
      y_m  = mean(trd$df$y),
      y_s  = sd(trd$df$y),
      xi_m = mean(trd$df$xi),
      xi_s = sd(trd$df$xi),
      xj_m = mean(trd$df$xj),
      xj_s = sd(trd$df$xj)
    )
    zs <- function(x, m, s)
      (x - m) / ifelse(is.finite(s) && s > 0, s, 1)
    for (nm in c("y", "xi", "xj")) {
      trd$df[[nm]] <- zs(trd$df[[nm]], sc[[paste0(nm, "_m")]], sc[[paste0(nm, "_s")]])
      ted$df[[nm]] <- zs(ted$df[[nm]], sc[[paste0(nm, "_m")]], sc[[paste0(nm, "_s")]])
    }
  }

  # --- Guard: insufficient train samples per fold ---
  if (nrow(trd$df) < min_pairs)
    return(NULL)

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
#' @param resid_mode Residual mode string.
#' @param use_t Logical; enable t-scoring if df is available.
#' @param elpd_mode \code{"indep"} or \code{"kalman"}.
#' @param nu_fixed_override Optional fixed \eqn{\nu} for all folds.
#' @param sample_args_override Optional override of sampler args.
#' @param seed_override Deterministic seed per fold.
#' @return A list with \code{elpd}, \code{elpd_ppd}, and \code{fold_diag}; or \code{NULL}.
#' @noRd
#' @keywords internal
.fold_fit_and_score <- function(mod,
                                stan_list_base,
                                sample_args_base,
                                pair_in,
                                tr,
                                te,
                                resid_mode,
                                max_retries = 3,
                                use_t = FALSE,
                                elpd_mode = NULL,
                                silent_sampler = TRUE,
                                nu_fixed_override = NULL,
                                sample_args_override = NULL,
                                freeze_retry_hypers = FALSE,
                                seed_override = NULL,
                                min_pairs = min_pairs) {
  # 자동 디폴트: OU면 kalman, 아니면 indep
  if (is.null(elpd_mode)) {
    elpd_mode <- if (identical(resid_mode, "ou")) "kalman" else "indep"
  } else {
    elpd_mode <- match.arg(elpd_mode, c("indep","kalman"))
  }
  pts <- .build_train_test(pair_in, tr, te, use_global_scaling = TRUE)
  if (is.null(pts))
    return(NULL)

  sl <- stan_list_base
  sl$N   <- nrow(pts$train)
  sl$y   <- pts$train$y
  sl$xi  <- pts$train$xi
  sl$xj  <- pts$train$xj
  sl$S   <- length(unique(pts$train$subject))
  sl$sid <- as.integer(factor(pts$train$subject))
  sl$prev <- pts$prev_train
  sl$dt  <- pts$dt_train

  # K-fold에서도 ν 고정 주입을 허용 (훈련 진단 기반 글로벌 규칙; 누출 아님)
  # 이유: override 값은 오직 "훈련 단계의 샘플러 진단"으로만 결정되며,
  # 테스트 ELPD/관측에는 의존하지 않는다. 폴드 간 동일 ν로 비교 일관성 확보.

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
    # ν 고정 정책: override가 주어지면 그 값으로 **고정**(bump 금지)
    # (훈련 진단만으로 결정된 사전 규칙 → 누출 아님)
    nu_fix      = if (is.null(nu_fixed_override)) nu_fix else nu_fixed_override,
    max_nu_cap  = if (is.null(nu_fixed_override)) 7      else nu_fixed_override,
    ebfmi_thresh = 0.30,
    silent_sampler = silent_sampler,
    tag = "kfold",
    freeze_retry_hypers = freeze_retry_hypers   # ★ 하이퍼를 바꾸지 않는 리트라이
  )
  fit <- tryfit$fit

  # 방어적 동일성 확인(디버그용; 필요시 주석 처리)
  fa <- tryfit$final_args

  d <- .safe_draws_df(fit)
  te_df <- pts$test
  # Resolve df for Student-t scoring without relying on rlang::`%||%`
  use_t_flag <- isTRUE(use_t)
  nu_for_score <- NULL
  if (use_t_flag) {
    if (!is.null(tryfit$used_nu) && is.finite(tryfit$used_nu)) {
      nu_for_score <- tryfit$used_nu
    } else if (!is.null(stan_list_base$nu_fixed) &&
               is.finite(stan_list_base$nu_fixed)) {
      nu_for_score <- stan_list_base$nu_fixed
    } else if (exists("nu_fix") && is.finite(nu_fix)) {
      nu_for_score <- nu_fix
    }
  }
  llm <- .proj_loglik_subject(
    draws_df   = d,
    pair_in    = te_df,
    resid_mode = resid_mode,
    use_t      = use_t_flag,
    fit_full   = fit,
    nu_scalar  = nu_for_score,
    elpd_mode  = elpd_mode
  )
  # --- NEW: subject별 ELPD 벡터와 fold diagnostics 구성 ---
  # llm$full: (draws x n_test_subjects) 행렬, llm$subjects: 테스트 subject 벡터
  if (is.null(llm) ||
      is.null(llm$full) ||
      !is.matrix(llm$full) || ncol(llm$full) == 0) {
    # 테스트 데이터가 비정상인 폴드는 실패로 처리
    return(NULL)
  }
  elpd_vec <- stats::setNames(apply(llm$full, 2, .log_mean_exp), llm$subjects)
  elpd_ppd_vec <- elpd_vec / pmax(llm$n_obs, 1L)

  dg <- tryfit$diag
  fd <- data.frame(
    n_retries      = tryfit$n_retries,
    nu_used        = tryfit$used_nu,
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

  list(elpd = elpd_vec, elpd_ppd = elpd_ppd_vec, fold_diag = fd)
}

#' Repeated K-fold evaluation (subject-level) with diagnostics payload
#'
#' Runs \code{K × R} folds, aggregates subject-wise ELPD and PPD-normalized ELPD,
#' and returns fold diagnostics and split manifests for reproducibility.
#'
#' @inheritParams .fold_fit_and_score
#' @param n_workers_kfold Number of workers (multisession) for fold-level parallelism.
#' @param nu_fixed_kfold Fixed \eqn{\nu} to use in all folds (no bump).
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
                           resid_mode,
                           use_t = FALSE,
                           elpd_mode = NULL,
                           silent_sampler = TRUE,
                           max_retries = 3,
                           n_workers_kfold = 1L,
                           min_pairs = 4,
                           nu_fixed_kfold = NULL,         # ← 1차 K-fold에서 고정할 ν (예: 본계산 used_nu)
                           freeze_retry_hypers = FALSE) { # ← 폴드 내 bump 금지 여부
  # 자동 디폴트: OU → "kalman", WN → "indep"
  if (is.null(elpd_mode)) {
    elpd_mode <- if (identical(resid_mode, "ou")) "kalman" else "indep"
  } else {
    elpd_mode <- match.arg(elpd_mode, c("indep","kalman"))
  }
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
  subs_all <- sort(unique(pair_in$subject))
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
    .fold_fit_and_score(
      mod,
      stan_list_base,
      sample_args_base,
      pair_in,
      task$tr,
      task$te,
      resid_mode,
      use_t,
      elpd_mode = elpd_mode,
      silent_sampler,
      max_retries = max_retries,
      nu_fixed_override   = nu_fixed_kfold,     # ← ν 고정 주입
      freeze_retry_hypers = freeze_retry_hypers,
      seed_override = task$seed,
      min_pairs = min_pairs
    )
  }

  # 실행: 순차 또는 병렬
  if (n_workers_kfold > 1L) {
    # 안전하게 PSOCK(멀티세션) 사용 권장 (fork 이슈 회피)
    op <- future::plan()
    on.exit(future::plan(op), add = TRUE)
    future::plan(future::multisession, workers = n_workers_kfold)
    if (has_progressr && identical(progress, "bar")) {
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
  agg      <- setNames(numeric(length(subs)), subs)
  cnt      <- setNames(integer(length(subs)), subs)
  agg_ppd  <- setNames(numeric(length(subs)), subs)
  cnt_ppd  <- setNames(integer(length(subs)), subs)

  fold_diag_df <- list()

  n_ok <- 0L
  n_fail <- 0L

  for (el in res_list) {
    if (is.null(el)) {
      n_fail <- n_fail + 1L
    } else {
      n_ok <- n_ok + 1L
      # elpd
      idx <- names(el$elpd)
      agg[idx]     <- agg[idx]     + el$elpd
      cnt[idx]     <- cnt[idx]     + 1L
      if (!is.null(el$elpd_ppd)) {
        agg_ppd[idx] <- agg_ppd[idx] + el$elpd_ppd
        cnt_ppd[idx] <- cnt_ppd[idx] + 1L
      }
      # fold diag
      fold_diag_df[[length(fold_diag_df) + 1L]] <- data.frame(
        n_retries      = el$fold_diag$n_retries,
        nu_used        = el$fold_diag$nu_used,
        ebfmi_min      = el$fold_diag$ebfmi_min,
        worst_rhat     = el$fold_diag$worst_rhat,
        min_ess_bulk   = el$fold_diag$min_ess_bulk,
        treedepth_hits = el$fold_diag$treedepth_hits,
        n_divergent    = el$fold_diag$n_divergent
      )
    }
  }

  elpd_subject     <- agg / pmax(cnt, 1L)
  elpd_subject_ppd <- agg_ppd / pmax(cnt_ppd, 1L)
  elpd_mean    <- if (length(elpd_subject)) mean(elpd_subject[is.finite(elpd_subject)], na.rm = TRUE) else NA_real_
  elpd_mean_ppd<- if (length(elpd_subject_ppd)) mean(elpd_subject_ppd[is.finite(elpd_subject_ppd)], na.rm = TRUE) else NA_real_

  elpd_mean    <- if (length(elpd_subject)) {
    m <- mean(elpd_subject[is.finite(elpd_subject)], na.rm = TRUE)
    if (is.finite(m)) m else NA_real_
  } else NA_real_


  fd <- if (length(fold_diag_df))
    do.call(rbind, fold_diag_df)
  else
    data.frame(
      n_retries = integer(),
      nu_used = integer(),
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
  nu_counts   <- if (nrow(fd))
    as.integer(table(fd$nu_used))
  else
    integer()
  names(nu_counts) <- if (nrow(fd))
    names(table(fd$nu_used))
  else
    character()


  list(
    type = "repeated-kfold",
    K = K,
    R = R,
    elpd_subject = elpd_subject,
    elpd_subject_ppd = elpd_subject_ppd,
    elpd_mean    = elpd_mean,
    elpd_mean_ppd= elpd_mean_ppd,
    elpd_method  = if (identical(elpd_mode,"kalman") && identical(resid_mode,"ou")) "kalman-ou" else "indep",
    n_folds_ok = n_ok,
    n_folds_fail = n_fail,
    # reproducibility payload for pseudo-BMA/stacking
    splits_df = splits_df,
    subject_ids = subs_all,
    subject_test_counts = subject_test_counts,
    fold_diag = fd,
    retry_total = retry_total,
    retry_mean  = retry_mean,
    nu_used_counts = nu_counts
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
.run_one <- function(target,
                     partner,
                     ctx,
                     seed_override = NULL,
                     progress_local = progress) {

  # ---- 컨텍스트 바인딩(워커 환경 내에서 사용) ----
  meta_df              <- ctx$meta_df
  sm_mat               <- ctx$sm_mat
  eps                  <- ctx$eps
  min_pairs            <- ctx$min_pairs
  compute_elpd         <- ctx$compute_elpd
  standardize_by_subject <- ctx$standardize_by_subject
  transform            <- ctx$transform
  lag                  <- ctx$lag
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
  nu_fix               <- ctx$nu_fix
  use_student_t        <- ctx$use_student_t
  elpd_mode            <- ctx$elpd_mode
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
  resid_mode           <- ctx$resid_mode

  # PF 옵션
  use_pathfinder_init <- isTRUE(ctx$use_pathfinder_init)
  pf_num_paths        <- ctx$pf_num_paths
  pf_draws            <- ctx$pf_draws
  pf_history_size     <- ctx$pf_history_size
  pf_max_lbfgs_iters  <- ctx$pf_max_lbfgs_iters
  pf_psis_resample    <- ctx$pf_psis_resample

  # ---- 병렬 워커에서 Stan 모델 핸들 확보: 이미 빌드된 exe 재사용 ----
  # cmdstanr는 R6 클래스를 직접 new()로 만들지 말고, cmdstan_model()을 사용해야 함.
  if (!is.null(ctx$mod_exe_file) && file.exists(ctx$mod_exe_file)) {
    mod <- cmdstanr::cmdstan_model(stan_file = NULL, exe_file = ctx$mod_exe_file)
  } else {
    # 안전한 폴백: exe 경로가 누락/손상되면 캐시에서 로드(재컴파일은 캐시 히트 시 생략됨)
    mod <- get_pclv_model(quiet = quiet)
  }

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
} else {
  pair_df <- .build_pair_df_smoothed(sm_mat, meta_df, partner, target, eps, min_pairs)
}
  if (is.null(pair_df))
    return(NULL)
  std_by_subj_eff <- if (isTRUE(compute_elpd))
    FALSE
  else
    standardize_by_subject
  if (isTRUE(compute_elpd) && isTRUE(standardize_by_subject)) {
    if (!quiet) {
      message(
        "compute_elpd=TRUE with standardize_by_subject=TRUE: ",
        "disabling subject-wise standardization to avoid OOF leakage; ",
        "K-fold will apply train-only scaling."
      )
    }
  }

  pair_in <- .make_pair_inputs_glv(
    pair_df = pair_df,
    transform = transform,
    lag = lag,
    zero_mode_alr = zero_mode_alr,
    minpos_alpha = minpos_alpha,
    minpos_base  = minpos_base,
    eps_fixed = eps_fixed,
    lib_eps_c = lib_eps_c,
    rest_floor_frac = rest_floor_frac,
    alr_cap = 12,
    smooth_scale = smooth_scale,
    alr_spline_df = alr_spline_df,
    alr_spline_spar = alr_spline_spar,
    alr_spline_cv = alr_spline_cv,
    nz_partner_min_frac = nz_partner_min_frac,
    standardize_by_subject = std_by_subj_eff
  )
  if (is.null(pair_in) || nrow(pair_in) < min_pairs)
    return(NULL)

  # 결측/비유한 제거
  pair_in <- pair_in[is.finite(pair_in$y) &
                       is.finite(pair_in$xi) &
                       is.finite(pair_in$xj) &
                       is.finite(pair_in$time), , drop = FALSE]
  if (!nrow(pair_in))
    return(NULL)

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
        dtv[ix[k]] <- if (is.finite(dt_k) &&
                          dt_k > 0)
          dt_k
        else
          1e-6
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
      dt   = dtv,
      use_student_t = as.integer(isTRUE(use_student_t)),
      # 0/1
      nu_fixed      = nu_fix                              # 시도별 래퍼에서 업데이트됨
    )
  )

  # 앞쪽 값 우선으로 중복 키 제거
  nm <- names(stan_list)
  if (anyDuplicated(nm))
    stan_list <- stan_list[match(unique(nm), nm)]

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
  if (use_pathfinder_init) {
    if (progress_local != "none")
      cat(sprintf("%s pair: running Pathfinder for inits\n", pair_tag))
    pf_dir <- file.path(tempdir(), sprintf("glvpf_%s", format(Sys.time(), "%Y%m%d%H%M%S")))
    dir.create(pf_dir, recursive = TRUE, showWarnings = FALSE)
    pf_init_scalar <- 0.05
    pf_fit <- mod$pathfinder(
      data    = stan_list,
      seed    = seed_main,
      init    = pf_init_scalar,
      num_paths = pf_num_paths,
      draws     = pf_draws,
      history_size    = pf_history_size,
      max_lbfgs_iters = pf_max_lbfgs_iters,
      psis_resample   = pf_psis_resample,
      show_messages   = !silent_sampler,
      output_dir      = pf_dir,
      output_basename = paste0("glv_pf_", target, "_", partner)
    )
    # ▶ PF 결과를 per-chain init으로 만들되, 한 체인이라도 비면 숫자 스칼라로 강제 폴백
    base_args$init <- .init_from_pf_or_scalar(
      pf_fit  = pf_fit,
      mod     = mod,
      chains  = chains,
      # 메인 인자로 숫자 스칼라 init이 왔으면 그 값을 우선 사용, 아니면 0.2
      fallback_numeric = if (is.numeric(init) && length(init) == 1L && is.finite(init)) init else 0.2,
      verbose = (progress_local != "none"),
      tag     = pair_tag
    )
    # Pathfinder 사용 시 warmup 단축
    if (!is.null(iter_warmup) && iter_warmup >= 1500) {
      base_args$iter_warmup <- max(500L, floor(iter_warmup/2))
    }
  }

  # ---[ NEW ] 추가 유틸: init 유효성 강검증 & 안전 폴백 -------------------------
  # list-of-lists 형식 점검, 빈 리스트 제거, 체인 수 보정, 전부 실패 시 숫자 init 반환
  .as_chain_inits_safe <- function(init_obj, chains,
                                   fallback_numeric = 0.2,
                                   verbose_tag = NULL) {
    say <- function(msg) {
      if (!is.null(verbose_tag)) cat(sprintf("%s %s\n", verbose_tag, msg))
    }
    # 1) 숫자 스칼라/함수면 그대로 허용(샘플러가 해석)
    if (is.numeric(init_obj) && length(init_obj) == 1L && is.finite(init_obj)) {
      return(init_obj)
    }
    if (is.function(init_obj)) {
      return(init_obj)
    }
    # 2) named list면 체인수만큼 복제
    if (is.list(init_obj) && length(init_obj) > 0L && !is.list(init_obj[[1L]])) {
      if (!length(names(init_obj))) {
        say("init named list missing names → fallback numeric init")
        return(fallback_numeric)
      }
      return(rep(list(init_obj), chains))
    }
    # 3) list-of-lists면 각 체인 점검
    if (is.list(init_obj) && length(init_obj) > 0L && is.list(init_obj[[1L]])) {
      L <- length(init_obj)
      # 길이 보정
      if (L != chains) {
        if (L == 1L) {
          init_obj <- rep(init_obj, chains)
        } else {
          say(sprintf("init length (%d) != chains (%d) → fallback numeric init", L, chains))
          return(fallback_numeric)
        }
      }
      # 체인별 빈 리스트/NA 제거
      ok_all <- TRUE
      for (k in seq_len(chains)) {
        li <- init_obj[[k]]
        # 비었거나 이름 없음 → 실패
        if (!is.list(li) || length(li) == 0L || !length(names(li))) { ok_all <- FALSE; break }
        # 숫자 아닌 값/NA/Inf 제거
        keep <- vapply(li, function(v) is.numeric(v) && length(v) == 1L && is.finite(v), logical(1))
        li <- li[keep]
        init_obj[[k]] <- li
        if (length(li) == 0L) { ok_all <- FALSE; break }
      }
      if (!ok_all) {
        say("init contains empty lists after validation → fallback numeric init")
        return(fallback_numeric)
      }
      return(init_obj)
    }
    # 4) 나머지 전부 실패 → 숫자 스칼라 폴백
    say("init type not recognized → fallback numeric init")
    fallback_numeric
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
    nu_fix = nu_fix,
    silent_sampler = silent_sampler,
    ebfmi_thresh = 0.30,
    max_nu_cap = 7,
    tag = tag_lbl
  )

  # --- main run done ---
  if (progress_local != "none")
    cat(sprintf("%s pair: main run completed\n", pair_tag))

  fit_full   <- res_try$fit
  diag       <- res_try$diag
  n_retries  <- res_try$n_retries
  fit_failed <- res_try$fit_failed
  used_nu    <- res_try$used_nu
  final_args_main <- res_try$final_args
  nu_final        <- used_nu

  # 5) 드로우 요약 ---------------------------------------------------------------
  if (quiet) {
    d <- suppressWarnings(.safe_draws_df(fit_full))
  } else {
    d <- .safe_draws_df(fit_full)  # a_ij, a_ii, r0, tau_r, sigma, sigma_ou, lambda, ...
  }
  if (isTRUE(use_student_t) && !("nu" %in% names(d))) {
    d$nu <- rep(used_nu, nrow(d))
  }

  # 방어적 체크
  if (!all(c("a_ij", "a_ii") %in% names(d))) {
    stop("Expected parameters a_ij and a_ii not found in draws.")
  }

  aij <- d$a_ij
  aii <- d$a_ii
  p_sign2_dir  <- .clip01(2 * pmin(mean(aij > 0), mean(aij < 0)))
  p_sign2_self <- .clip01(2 * pmin(mean(aii > 0), mean(aii < 0)))

  # ----- 6) Repeated K-fold with outer retries ----------------------------------
  kfold <- NULL

  # ▼ k-fold 메타 기본값(ELPD 계산 안 할 때도 반환값 보장)
  kfold_rounds_attempted <- 0L
  kfold_failed_overall   <- NA
  kfold_n_ok   <- NA_integer_
  kfold_n_fail <- NA_integer_

  if (isTRUE(compute_elpd)) {
    seed_for_kfold <- if (is.null(final_args_main$seed))
      seed_main
    else
      final_args_main$seed
    if (progress_local != "none")
      cat(
        sprintf(
          "%s pair: k-fold evaluation started (K=%d, R=%d)\n",
          pair_tag,
          kfold_K,
          kfold_R
        )
      )
    sargs_kfold <- final_args_main
    sargs_kfold$step_size   <- NULL   # 초기 step size 강제 금지 (폴드별 adapt)
    sargs_kfold$inv_metric  <- NULL   # metric 파일/행렬 강제 금지
    sargs_kfold$metric_file <- NULL


    # 본계산에서 사용된 ν로 K-fold 전 폴드 **고정** 실행 (bump 금지) ---
    nu_kfold_main <- nu_final

    kfold <- .repkfold_eval(
      mod = mod,
      stan_list_base   = stan_list,
      sample_args_base = sargs_kfold,
      pair_in = pair_in,
      K = kfold_K,
      R = kfold_R,
      seed = seed_for_kfold,
      resid_mode = resid_mode,
      use_t = isTRUE(use_student_t),
      elpd_mode = elpd_mode,
      silent_sampler = silent_sampler,
      n_workers_kfold = n_workers_kfold_eff,
      max_retries = max_retries,
      nu_fixed_kfold = nu_kfold_main,          # ← ν 고정
      freeze_retry_hypers = TRUE,
      min_pairs = min_pairs
    )
    if (progress_local != "none")
      cat(sprintf("%s pair: k-fold evaluation completed\n", pair_tag))
    kfold_n_ok   <- if (is.null(kfold))
      NA_integer_
    else
      kfold$n_folds_ok
    kfold_n_fail <- if (is.null(kfold))
      NA_integer_
    else
      kfold$n_folds_fail
    # 폴드 훈련-진단 실패가 하나라도 있으면 ν=7로 K-fold **전부 재실행** ---
    # (정보 누출 아님: 테스트 ELPD/관측은 보지 않고, 훈련 진단만 사용)
    kfold_outer_rounds_local <- 1L
    bad_count <- 0L
    if (!is.null(kfold) && is.data.frame(kfold$fold_diag) && nrow(kfold$fold_diag)) {
      fd <- kfold$fold_diag
      bad_flag <- with(fd,
                       (is.finite(n_divergent)     & n_divergent     > 0) |
                         (is.finite(ebfmi_min)       & ebfmi_min       < 0.30) |
                         (is.finite(treedepth_hits)  & treedepth_hits  > 0) |
                         (is.finite(worst_rhat)      & worst_rhat      >= 1.05) |
                         (is.finite(min_ess_bulk)    & min_ess_bulk    < 400)
      )
      bad_count <- sum(bad_flag, na.rm = TRUE)
    }
    if (isTRUE(bad_count > 0) && is.finite(nu_kfold_main) && nu_kfold_main < 7) {
      if (progress_local != "none")
        cat(sprintf("%s pair: k-fold retrial triggered by TRAIN diagnostics → nu=7 (global)\n", pair_tag))
      kfold <- .repkfold_eval(
        mod = mod,
        stan_list_base   = stan_list,
        sample_args_base = sargs_kfold,
        pair_in = pair_in,
        K = kfold_K,
        R = kfold_R,
        seed = seed_for_kfold,
        resid_mode = resid_mode,
        use_t = isTRUE(use_student_t),
        elpd_mode = elpd_mode,
        silent_sampler = silent_sampler,
        n_workers_kfold = n_workers_kfold_eff,
        max_retries = max_retries,
        nu_fixed_kfold = 7,                   # ← 글로벌 한 단계 상승(최대 1회)
        freeze_retry_hypers = TRUE
      )
      kfold_outer_rounds_local <- 2L
      if (progress_local != "none")
        cat(sprintf("%s pair: k-fold evaluation (nu=7) completed\n", pair_tag))
      kfold_n_ok   <- if (is.null(kfold)) NA_integer_ else kfold$n_folds_ok
      kfold_n_fail <- if (is.null(kfold)) NA_integer_ else kfold$n_folds_fail


      # --- 본계산도 ν=7로 일관화: 테스트 점수는 쓰지 않고, K-fold '훈련 진단' 트리거만으로 재적합 ---
      # 누출이 아닌 이유:
      #  (1) ν=7 결정은 오직 K-fold '훈련' 진단(다이버전스/E-BFMI/treedepth/R-hat/ESS)으로만 내림
      #  (2) 테스트 ELPD/관측치는 보지 않으며, 선택/재적합 판단에 사용되지 않음
      #  (3) 최종 보고는 사양 확정 후 단일 CV 및 본계산 결과만 사용
      if (is.finite(nu_final) && nu_final < 7) {
        if (progress_local != "none")
          cat(sprintf("%s pair: main refit with nu=7 for spec consistency\n", pair_tag))
        base_args_refit <- final_args_main
        # 리핏 시 bump/하이퍼 변경 금지, seed/체인 동일, metric 그대로
        base_args_refit$data$nu_fixed <- 7
        base_args_refit$adapt_delta   <- final_args_main$adapt_delta
        base_args_refit$max_treedepth <- final_args_main$max_treedepth
        base_args_refit$metric        <- final_args_main$metric
        # 리핏은 한 번만, bump 금지
        refit_try <- .sample_with_retry(
          mod = mod,
          base_args = base_args_refit,
          stan_list = base_args_refit$data,
          max_retries = 0,
          nu_fix = 7,
          max_nu_cap = 7,
          silent_sampler = silent_sampler,
          freeze_retry_hypers = TRUE,
          tag = sprintf("%s\u2192%s main-refit-nu7", partner, target)
        )
        fit_full <- refit_try$fit
        diag     <- refit_try$diag
        nu_final <- 7
        d <- .safe_draws_df(fit_full)
        aij <- d$a_ij; aii <- d$a_ii
        p_sign2_dir  <- .clip01(2 * pmin(mean(aij > 0), mean(aij < 0)))
        p_sign2_self <- .clip01(2 * pmin(mean(aii > 0), mean(aii < 0)))
      }
    }
  }




  # 7) 리턴 ----------------------------------------------------------------------
  list(
    n_pairs = nrow(pair_in),

    a_mean  = mean(aij),
    a_sd = stats::sd(aij),
    a_q2.5  = stats::quantile(aij, 0.025),
    a_q97.5 = stats::quantile(aij, 0.975),
    p_sign2 = p_sign2_dir,

    aii_mean = mean(aii),
    aii_sd = stats::sd(aii),
    aii_q2.5 = stats::quantile(aii, 0.025),
    aii_q97.5 = stats::quantile(aii, 0.975),
    p_sign2_self = p_sign2_self,

    # 원시 진단치/리트라이 메타(후처리 summariser에서 필터 예정)
    diag = diag,
    n_retries = n_retries,
    fit_failed = fit_failed,
    nu_fixed_used = nu_final,

    # 스무딩 메타
    smoothed = sflag,
    smooth_scale = ssc,
    smooth_edf_mean = sedf,
    subjects = unique(pair_in$subject),

    # Repeated K-fold 결과 요약 (+ pseudo-BMA/stacking용 payload)
    kfold = kfold,
    kfold_mean = if (is.null(kfold)) NA_real_ else kfold$elpd_mean,
    kfold_method = if (is.null(kfold)) NA_character_ else kfold$elpd_method,

    # 용어 명확화: 바깥 재평가 라운드(outer)만 카운트
    kfold_outer_rounds = if (isTRUE(compute_elpd)) kfold_outer_rounds_local else 0L,
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
    # subject-level OOF elpd vector (named) for model averaging/stacking
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
    # split manifest (r,k,test_subjects) & seed to ensure same partitions across models
    kfold_splits        = list(if (is.null(kfold))
      NULL
      else
        kfold$splits_df),
    kfold_seed_used     = if (is.null(kfold))
      NA_integer_
    else
      seed_for_kfold,
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
    kfold_nu_used_counts = list(if (is.null(kfold) ||
                                    is.null(kfold$nu_used_counts))
      NULL
      else
        kfold$nu_used_counts)
  )
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

#' Expand per-pair self ELPD payloads to long format
#'
#' Unnests subject-level self-effect ELPD (and PPD-normalized ELPD) from
#' pairwise results into a long tibble.
#'
#' @param df Tibble with nested K-fold payload columns for \code{i} and \code{j}.
#' @return A tibble with \code{taxon,subject,elpd,elpd_ppd,n_test,elpd_method}.
#' @noRd
#' @keywords internal
.expand_self_pw <- function(df) {
  ii_rows <- purrr::pmap_dfr(
    list(df$i, df$kfold_subject_ids_ij, df$kfold_subject_ij, df$kfold_subject_ppd_ij,
         df$kfold_elpd_method_ij, df$kfold_subject_counts_ij),
    function(i, ids, elv, elv_ppd, method, cnts) {
      # 방어적 언래핑(리스트 한 겹만 더 들어온 경우)
      if (is.list(ids) && length(ids) == 1)       ids <- ids[[1]]
      if (is.list(elv) && length(elv) == 1)       elv <- elv[[1]]
      if (is.list(elv_ppd) && length(elv_ppd) == 1) elv_ppd <- elv_ppd[[1]]
      if (is.null(ids) || is.null(elv)) return(NULL)
      tab <- .as_named_frame2(ids, elv, elv_ppd)
      tab$n_test <- as.integer(if (is.null(cnts)) NA else cnts[tab$subject])
      tab$taxon <- i
      tab$elpd_method <- if (is.null(method)) NA_character_ else method
      tab[, c("taxon","subject","elpd","elpd_ppd","n_test","elpd_method")]
    }
  )
  jj_rows <- purrr::pmap_dfr(
    list(df$j, df$kfold_subject_ids_ji, df$kfold_subject_ji, df$kfold_subject_ppd_ji,
         df$kfold_elpd_method_ji, df$kfold_subject_counts_ji),
    function(j, ids, elv, elv_ppd, method, cnts) {
      # 방어적 언래핑
      if (is.list(ids) && length(ids) == 1)       ids <- ids[[1]]
      if (is.list(elv) && length(elv) == 1)       elv <- elv[[1]]
      if (is.list(elv_ppd) && length(elv_ppd) == 1) elv_ppd <- elv_ppd[[1]]
      if (is.null(ids) || is.null(elv)) return(NULL)
      tab <- .as_named_frame2(ids, elv, elv_ppd)
      tab$n_test <- as.integer(if (is.null(cnts)) NA else cnts[tab$subject])
      tab$taxon <- j
      tab$elpd_method <- if (is.null(method)) NA_character_ else method
      tab[, c("taxon","subject","elpd","elpd_ppd","n_test","elpd_method")]
    }
  )
  dplyr::bind_rows(ii_rows, jj_rows)
}

#' Expand per-pair cross ELPD payloads to long directed format
#'
#' Unnests subject-level cross-edge ELPD (and PPD-normalized ELPD) for both
#' directions into a long tibble of \code{from → to} rows.
#'
#' @param df Tibble with nested K-fold payload columns.
#' @return A tibble with \code{from,to,subject,elpd,elpd_ppd,n_test,elpd_method}.
#' @noRd
#' @keywords internal
.expand_cross_pw <- function(df) {
  # ij
  ij_rows <- purrr::pmap_dfr(
    list(df$i, df$j, df$kfold_subject_ids_ij, df$kfold_subject_ij, df$kfold_subject_ppd_ij,
         df$kfold_elpd_method_ij, df$kfold_subject_counts_ij),
    function(i, j, ids, elv, elv_ppd, method, cnts) {
      if (is.null(ids) || is.null(elv)) return(NULL)
      tab <- .as_named_frame2(ids, elv, elv_ppd)
      tab$n_test <- as.integer(if (is.null(cnts)) NA else cnts[tab$subject])
      tab$from <- j; tab$to <- i
      tab$elpd_method <- if (is.null(method)) NA_character_ else method
      tab[, c("from","to","subject","elpd","elpd_ppd","n_test","elpd_method")]
    }
  )
  # ji
  ji_rows <- purrr::pmap_dfr(
    list(df$j, df$i, df$kfold_subject_ids_ji, df$kfold_subject_ji, df$kfold_subject_ppd_ji,
         df$kfold_elpd_method_ji, df$kfold_subject_counts_ji),
    function(j, i, ids, elv, elv_ppd, method, cnts) {
      if (is.null(ids) || is.null(elv)) return(NULL)
      tab <- .as_named_frame2(ids, elv, elv_ppd)
      tab$n_test <- as.integer(if (is.null(cnts)) NA else cnts[tab$subject])
      tab$from <- i; tab$to <- j
      tab$elpd_method <- if (is.null(method)) NA_character_ else method
      tab[, c("from","to","subject","elpd","elpd_ppd","n_test","elpd_method")]
    }
  )
  dplyr::bind_rows(ij_rows, ji_rows)
}
