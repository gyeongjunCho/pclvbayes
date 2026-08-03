#' Summarize Bayesian pcLV results
#'
#' @description
#' Given the result list from \code{fit_pclv_bayes()}, returns a summary table
#' for directed cross effects or self effects.
#'
#' Directed cross coefficients are pair-to-rest dynamic coefficients and are
#' not generally absolute direct-gLV effects. Posterior sign support and MCMC
#' diagnostics are the primary edge evidence. Pair-specific repeated K-fold
#' ELPD may be retained as supporting predictive evidence in the fit object.
#'
#' Cross-pair model weights are not returned because distinct pair-to-rest
#' models generally predict different transformed outcomes and therefore do
#' not form a common-outcome stacking or pseudo-BMA model set.
#'
#' @details
#' Diagnostic criteria depend on \code{diag_mode}:
#' \itemize{
#'   \item \strong{strict}: \code{rhat < 1.01}, \code{essb > 1000},
#'     \code{esst > 1000}, \code{div == 0}, \code{tdhit == 0}, and
#'     \code{ebfmi_min >= 0.40}.
#'   \item \strong{moderate}: \code{rhat < 1.05}, \code{essb > 400},
#'     \code{esst > 400}, \code{div <= 8}, \code{tdhit <= 80}, and
#'     \code{ebfmi_min >= 0.30}.
#' }
#'
#' Indeterminate directions remain explicit and are not interpreted as zero.
#' Tables are sorted by \code{bayes_FDR} in ascending order with missing values
#' last.
#'
#' @param df A list returned by \code{fit_pclv_bayes()} containing
#'   \code{$cross}, \code{$self}, and \code{$raw}.
#' @param alpha Optional threshold for Bayes-FDR/LFSR. When supplied, adds
#'   \code{pass_bayes_fdr}.
#' @param diag_mode Either \code{"moderate"} or \code{"strict"}.
#' @param interaction One of \code{"cross"} or \code{"self"}.
#'
#' @return A list containing the requested \code{$cross} or \code{$self}
#'   summary table.
#'
#' @export
summarize_bayes_pclv <- function(df,
                                 alpha = NULL,
                                 diag_mode = c("moderate","strict"),
                                 interaction = "cross") {
  # deps
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("Packages 'dplyr' and 'tibble' are required.")
  }

  diag_mode <- match.arg(diag_mode, c("moderate","strict"))
  interaction <- match.arg(interaction, c("cross","self"))
  do_cross <- (interaction == "cross")
  do_self  <- (interaction == "self")

  # thresholds
  thr <- switch(diag_mode,
                "strict"   = list(rhat = 1.01, ess = 1000, div = 0L,  tdhit = 0L,  ebfmi = 0.40),
                "moderate" = list(rhat = 1.05, ess =  400, div = 8L,  tdhit = 80L, ebfmi = 0.30)
  )

  # inputs
  if (!is.list(df) || !all(c("cross","self","raw") %in% names(df))) {
    stop("`df` must include $cross, $self, and $raw.")
  }
  cross <- tibble::as_tibble(df$cross)
  self  <- tibble::as_tibble(df$self)
  raw   <- tibble::as_tibble(df$raw)

  # Results created before chain-specific reporting are interpreted using the
  # historical direction_ok contract. New results always carry these fields.
  for (suffix in c("ij", "ji")) {
    class_col <- paste0("diagnostic_class_", suffix)
    ok_col <- paste0("direction_ok_", suffix)
    if (!class_col %in% names(raw)) {
      raw[[class_col]] <- if (ok_col %in% names(raw))
        ifelse(raw[[ok_col]], "converged", "sampler_diagnostics_failed") else "converged"
    }
  }

  # diagnostics per direction from 'raw'
  cross_diag <- dplyr::bind_rows(
    dplyr::transmute(
      raw,
      from = .data$j, to = .data$i,
      n_pairs = .data$n_pairs_ij,
      rhat = .data$rhat_ij, essb = .data$essb_ij, esst = .data$esst_ij,
      div = .data$div_ij, tdhit = .data$tdhit_ij,
      ebfmi_min = .data$ebfmi_min_ij,
      diagnostic_class_diag = .data$diagnostic_class_ij
    ),
    dplyr::transmute(
      raw,
      from = .data$i, to = .data$j,
      n_pairs = .data$n_pairs_ji,
      rhat = .data$rhat_ji, essb = .data$essb_ji, esst = .data$esst_ji,
      div = .data$div_ji, tdhit = .data$tdhit_ji,
      ebfmi_min = .data$ebfmi_min_ji,
      diagnostic_class_diag = .data$diagnostic_class_ji
    )
  )
  self_diag <- dplyr::bind_rows(
    dplyr::transmute(
      raw,
      taxon = .data$i,
      n_pairs = .data$n_pairs_ij,
      rhat = .data$rhat_ij, essb = .data$essb_ij, esst = .data$esst_ij,
      div = .data$div_ij, tdhit = .data$tdhit_ij,
      ebfmi_min = .data$ebfmi_min_ij,
      diagnostic_class = .data$diagnostic_class_ij
    ),
    dplyr::transmute(
      raw,
      taxon = .data$j,
      n_pairs = .data$n_pairs_ji,
      rhat = .data$rhat_ji, essb = .data$essb_ji, esst = .data$esst_ji,
      div = .data$div_ji, tdhit = .data$tdhit_ji,
      ebfmi_min = .data$ebfmi_min_ji,
      diagnostic_class = .data$diagnostic_class_ji
    )
  )

  # ---------------------------
  # CROSS summary (diag_ok gating)
  # ---------------------------
  if (do_cross) {
    # 1) 진단 붙이고 diag_ok 먼저 계산
    cross_merged <- cross |>
      dplyr::left_join(cross_diag, by = c("from","to")) |>
      dplyr::mutate(
        sampler_diag_ok = .diag_ok_fun(.data$rhat, .data$essb, .data$esst,
                                       .data$div, .data$tdhit, .data$ebfmi_min, thr = thr),
        diagnostic_class = dplyr::coalesce(.data$diagnostic_class, .data$diagnostic_class_diag),
        diag_ok = .data$sampler_diag_ok & .data$diagnostic_class == "converged",
        a_sign = dplyr::case_when(
          is.finite(.data$a_mean) & .data$a_mean >  0 ~ "+",
          is.finite(.data$a_mean) & .data$a_mean <  0 ~ "-",
          TRUE ~ "0"
        ),
        bayes_FDR = ifelse(.data$diag_ok, .lfsr_safe(.data$p_sign2), NA_real_)
      )

    # Cross-pair pseudo-BMA and stacking are intentionally disabled.
    # Distinct incoming edges generally use different pair-to-rest responses,
    # so they are not common-outcome candidate models.
    cross_sum <- cross_merged |>
      dplyr::transmute(
        from, to,
        n_subjects = .data$n_subjects,
        n_pairs = .data$n_pairs,
        a_sign, a_mean, a_q2.5, a_q97.5,
        p_sign2, bayes_FDR,
        rhat, essb, esst, div, tdhit, ebfmi_min,
        diagnostic_class, sampler_diag_ok, diag_ok
      )

    if (!is.null(alpha)) {
      cross_sum <- dplyr::mutate(
        cross_sum,
        pass_bayes_fdr = is.finite(.data$bayes_FDR) & .data$bayes_FDR <= alpha
      )
    }
    cross_sum <- dplyr::arrange(cross_sum, .data$bayes_FDR)
  }

  # ---------------------------
  # SELF summary (diag_ok gating)
  # ---------------------------
  if (do_self) {
    self_sum <- self |>
      dplyr::left_join(self_diag, by = c("taxon")) |>
      dplyr::mutate(
        from = .data$taxon,
        to   = .data$taxon,
        a_sign = dplyr::case_when(
          is.finite(.data$a_self_mean) & .data$a_self_mean >  0 ~ "+",
          is.finite(.data$a_self_mean) & .data$a_self_mean <  0 ~ "-",
          TRUE ~ "0"
        ),
        sampler_diag_ok = .diag_ok_fun(.data$rhat, .data$essb, .data$esst,
                                       .data$div, .data$tdhit, .data$ebfmi_min, thr = thr),
        diag_ok = .data$sampler_diag_ok & .data$diagnostic_class == "converged",
        bayes_FDR = ifelse(.data$diag_ok, .lfsr_safe(.data$p_sign2_self), NA_real_)
      ) |>
      dplyr::transmute(
        from, to,
        n_pairs = .data$n_pairs,
        a_sign,
        a_mean  = .data$a_self_mean,
        a_q2.5  = .data$a_self_q2.5,
        a_q97.5 = .data$a_self_q97.5,
        p_sign2 = .data$p_sign2_self,
        bayes_FDR,
        rhat, essb, esst, div, tdhit, ebfmi_min,
        diagnostic_class, sampler_diag_ok, diag_ok
      )

    if (!is.null(alpha)) {
      self_sum <- dplyr::mutate(
        self_sum,
        pass_bayes_fdr = is.finite(.data$bayes_FDR) & .data$bayes_FDR <= alpha
      )
    }
    self_sum <- dplyr::arrange(self_sum, .data$bayes_FDR)
  }

  out <- list()
  if (do_cross) out$cross <- cross_sum
  if (do_self)  out$self  <- self_sum
  out
}
