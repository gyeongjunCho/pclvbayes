# R/summarize_bayes_pglv.R

#' Convert wide (i–j) gLV table into a long directed edge list (Bayes-only)
#'
#' @description
#' One row per **directed** edge (j -> i and i -> j).
#' Bayesian-only: uses p_sign2, lfsr, q_bayes, psp. No frequentist p/q.
#' PLOTTING TIP: Do NOT map |a_mean| to edge width. Use `delta_elpd_per_obs`
#' as edge thickness (importance) and `a_sign` for color. The absolute magnitude
#' of a_mean from pairwise fits is not comparable across pairs.
#'
#' @param df Output of [fit_glv_pairwise()].
#' @param alpha Optional numeric in (0,1]; if given, keep edges with q_bayes ≤ alpha.
#' @param p_sign2_max Optional numeric in (0,1]; if given, keep edges with p_sign2 ≤ p_sign2_max.
#' @return data.frame with:
#'   from,to,n_pairs,a_mean,a_sign("pos/neg/zero"),a_q2.5,a_q97.5,
#'   p_sign2,lfsr,q_bayes,psp,rhat,essb,esst,div,tdhit,ok_diag,keep_final,
#'   delta_elpd,delta_elpd_se,delta_elpd_per_obs,n_heldout,
#'   smoothed/smooth_* (if present), resid_mode/prior_mode (if present).
#' @export
summarize_bayes_pglv <- function(df, alpha = NULL, p_sign2_max = NULL) {
  stopifnot(is.data.frame(df))

  req_cols <- c("i","j",
                "n_pairs_ij","a_ij_mean","a_ij_sd","a_ij_q2.5","a_ij_q97.5",
                "n_pairs_ji","a_ji_mean","a_ji_sd","a_ji_q2.5","a_ji_q97.5")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "), call. = FALSE)

  has_bayes_ij <- any(c("lfsr_ij","p_sign2_ij") %in% names(df))
  has_bayes_ji <- any(c("lfsr_ji","p_sign2_ji") %in% names(df))
  if (!has_bayes_ij || !has_bayes_ji)
    stop("`df` must contain Bayesian significance columns (lfsr_* or p_sign2_*).", call. = FALSE)

  .as_lfsr <- function(lfsr_col, ps2_col) {
    if (!is.null(lfsr_col)) return(pmin(pmax(lfsr_col, 0), 0.5))
    if (!is.null(ps2_col))  return(pmin(pmax(ps2_col/2, 0), 0.5))
    rep(NA_real_, nrow(df))
  }
  .q_from_l <- function(v) {
    q <- rep(NA_real_, length(v)); nn <- which(!is.na(v))
    if (length(nn)) { l <- pmin(pmax(v[nn], 0), 0.5); oo <- nn[order(l)]; q[oo] <- cumsum(l[order(l)]) / seq_along(oo) }
    q
  }
  .sign_str <- function(x) ifelse(x > 0, "pos", ifelse(x < 0, "neg", "zero"))

  # --- IJ: j -> i ---
  lfsr_ij <- .as_lfsr(if ("lfsr_ij" %in% names(df)) df$lfsr_ij else NULL,
                      if ("p_sign2_ij" %in% names(df)) df$p_sign2_ij else NULL)
  qbij_ij <- if ("q_bayes_ij" %in% names(df)) df$q_bayes_ij else .q_from_l(lfsr_ij)
  psp_ij  <- if ("psp_ij" %in% names(df)) df$psp_ij else (1 - lfsr_ij)
  ps2_ij  <- if ("p_sign2_ij" %in% names(df)) df$p_sign2_ij else (2 * lfsr_ij)

  df_ij <- df %>%
    dplyr::transmute(
      from = j, to = i,
      n_pairs = n_pairs_ij,
      a_mean = a_ij_mean,
      a_sign = .sign_str(a_ij_mean),                # sign-only (no magnitude use)
      a_q2.5 = a_ij_q2.5, a_q97.5 = a_ij_q97.5,
      p_sign2 = ps2_ij, lfsr = lfsr_ij, q_bayes = qbij_ij, psp = psp_ij,
      rhat = if ("rhat_ij" %in% names(df)) rhat_ij else NA_real_,
      essb = if ("essb_ij" %in% names(df)) essb_ij else NA_real_,
      esst = if ("esst_ij" %in% names(df)) esst_ij else NA_real_,
      div  = if ("div_ij"  %in% names(df)) div_ij  else NA_real_,
      tdhit= if ("tdhit_ij"%in% names(df)) tdhit_ij else NA_real_,
      ok_diag   = if ("ok_diag_ij" %in% names(df)) ok_diag_ij else NA,
      keep_final= if ("keep_ij_final" %in% names(df)) keep_ij_final else NA,
      delta_elpd         = if ("delta_elpd_ij" %in% names(df)) df$delta_elpd_ij else NA_real_,
      delta_elpd_se      = if ("delta_elpd_se_ij" %in% names(df)) df$delta_elpd_se_ij else NA_real_,
      delta_elpd_per_obs = if ("delta_elpd_per_obs_ij" %in% names(df)) df$delta_elpd_per_obs_ij else NA_real_,
      n_heldout          = if ("n_heldout_ij" %in% names(df)) df$n_heldout_ij else NA_integer_,
      resid_mode = if ("resid_mode" %in% names(df)) df$resid_mode else NA_character_,
      prior_mode = if ("prior_mode" %in% names(df)) df$prior_mode else NA_character_,
      smoothed   = if ("smoothed" %in% names(df)) df$smoothed else NA,
      smooth_scale = if ("smooth_scale" %in% names(df)) df$smooth_scale else NA_character_,
      smooth_df     = if ("smooth_df" %in% names(df)) df$smooth_df else NA_real_,
      smooth_spar   = if ("smooth_spar" %in% names(df)) df$smooth_spar else NA_real_,
      smooth_cv     = if ("smooth_cv" %in% names(df)) df$smooth_cv else NA
    )

  # --- JI: i -> j ---
  lfsr_ji <- .as_lfsr(if ("lfsr_ji" %in% names(df)) df$lfsr_ji else NULL,
                      if ("p_sign2_ji" %in% names(df)) df$p_sign2_ji else NULL)
  qbij_ji <- if ("q_bayes_ji" %in% names(df)) df$q_bayes_ji else .q_from_l(lfsr_ji)
  psp_ji  <- if ("psp_ji" %in% names(df)) df$psp_ji else (1 - lfsr_ji)
  ps2_ji  <- if ("p_sign2_ji" %in% names(df)) df$p_sign2_ji else (2 * lfsr_ji)

  df_ji <- df %>%
    dplyr::transmute(
      from = i, to = j,
      n_pairs = n_pairs_ji,
      a_mean = a_ji_mean,
      a_sign = .sign_str(a_ji_mean),
      a_q2.5 = a_ji_q2.5, a_q97.5 = a_ji_q97.5,
      p_sign2 = ps2_ji, lfsr = lfsr_ji, q_bayes = qbij_ji, psp = psp_ji,
      rhat = if ("rhat_ji" %in% names(df)) rhat_ji else NA_real_,
      essb = if ("essb_ji" %in% names(df)) essb_ji else NA_real_,
      esst = if ("esst_ji" %in% names(df)) esst_ji else NA_real_,
      div  = if ("div_ji"  %in% names(df)) div_ji  else NA_real_,
      tdhit= if ("tdhit_ji"%in% names(df)) tdhit_ji else NA_real_,
      ok_diag   = if ("ok_diag_ji" %in% names(df)) ok_diag_ji else NA,
      keep_final= if ("keep_ji_final" %in% names(df)) keep_ji_final else NA,
      delta_elpd         = if ("delta_elpd_ji" %in% names(df)) df$delta_elpd_ji else NA_real_,
      delta_elpd_se      = if ("delta_elpd_se_ji" %in% names(df)) df$delta_elpd_se_ji else NA_real_,
      delta_elpd_per_obs = if ("delta_elpd_per_obs_ji" %in% names(df)) df$delta_elpd_per_obs_ji else NA_real_,
      n_heldout          = if ("n_heldout_ji" %in% names(df)) df$n_heldout_ji else NA_integer_,
      resid_mode = if ("resid_mode" %in% names(df)) df$resid_mode else NA_character_,
      prior_mode = if ("prior_mode" %in% names(df)) df$prior_mode else NA_character_,
      smoothed   = if ("smoothed" %in% names(df)) df$smoothed else NA,
      smooth_scale = if ("smooth_scale" %in% names(df)) df$smooth_scale else NA_character_,
      smooth_df     = if ("smooth_df" %in% names(df)) df$smooth_df else NA_real_,
      smooth_spar   = if ("smooth_spar" %in% names(df)) df$smooth_spar else NA_real_,
      smooth_cv     = if ("smooth_cv" %in% names(df)) df$smooth_cv else NA
    )

  out <- dplyr::bind_rows(df_ij, df_ji)

  # ---- Optional Bayes filters (intersection if both given) -------------------
  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || alpha <= 0 || alpha > 1) {
      stop("`alpha` must be in (0,1].", call. = FALSE)
    }
    out <- dplyr::filter(out, is.finite(q_bayes) & q_bayes <= alpha)
  }
  if (!is.null(p_sign2_max)) {
    if (!is.numeric(p_sign2_max) || p_sign2_max <= 0 || p_sign2_max > 1) {
      stop("`p_sign2_max` must be in (0,1].", call. = FALSE)
    }
    out <- dplyr::filter(out, is.finite(p_sign2) & p_sign2 <= p_sign2_max)
  }

  # Legacy heads-up
  if (any(c("p_ij","q_ij","p_ji","q_ji") %in% names(df))) {
    warning("Legacy frequentist columns (p_*, q_*) are ignored. This is Bayes-only.", call. = FALSE)
  }

  out
}
