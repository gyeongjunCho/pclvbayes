#' Summarize Bayesian pcLV results (cross/self, diagnostics, weights)
#'
#' @description
#' Given the result list from \code{fit_glv_pairwise()}, returns separate
#' summary tables for \strong{cross} (from→to) edges and \strong{self} effects.
#'
#' Included fields:
#' - Key intuitive metrics: \code{from}, \code{to}, \code{n_pairs}, \code{a_sign},
#'   \code{a_mean}, \code{a_q2.5}, \code{a_q97.5}
#' - Reliability metrics: \code{p_sign2}, \code{bayes_FDR} (=\code{LFSR}),
#'   optional \code{pass_bayes_fdr}
#' - Weight-like metrics: \code{pseudo_BMA} (no bootstrap), \code{pseudo_BMA_plus}
#'   (a.k.a. pseudo-BMA+; Bayesian bootstrap), and \code{stacking}
#'   (true stacking if available; otherwise omitted)
#' - Diagnostics: \code{rhat}, \code{essb}, \code{esst}, \code{div},
#'   \code{tdhit}, \code{ebfmi_min}, \code{diag_ok}
#'
#' @details
#' \itemize{
#' \item Diagnostic criteria depend on \code{diag_mode}:
#'   \itemize{
#'     \item \strong{strict}: \code{rhat < 1.01}, \code{essb > 1000}, \code{esst > 1000},
#'           \code{div == 0}, \code{tdhit == 0}, \code{ebfmi_min >= 0.40}
#'     \item \strong{moderate} (default): \code{rhat < 1.05}, \code{essb > 400}, \code{esst > 400},
#'           \code{div <= 8}, \code{tdhit <= 80}, \code{ebfmi_min >= 0.30}
#'   }
#' \item \code{pseudo_BMA}: a simple weight computed by taking the sum of
#'   pointwise ELPD (\eqn{\widehat{\mathrm{elpd}}}) and normalizing its
#'   \code{exp}-transformed values.
#' \item \code{pseudo_BMA_plus}: the same idea with an added Bayesian bootstrap
#'   (Dirichlet(1,…,1)) to regularize the weights.
#' \item \code{stacking}: if \code{use_true_stacking = TRUE}, pointwise ELPD
#'   (\code{df$elpd_pointwise_cross}) is available, and the \pkg{loo} package
#'   is installed, then true stacking weights are computed per target
#'   (each \code{to}) using the internal helper
#'   \code{.compute_stacking_weights_cross()} (columns sum to 1). Omitted
#'   otherwise.
#'
#' Each table is returned sorted by \strong{bayes_FDR} in ascending order (NA last).
#'
#' @param df A list returned by \code{fit_glv_pairwise()}, which must include
#'   at least \code{$cross}, \code{$self}, and \code{$raw}. To compute true stacking,
#'   it should also include \code{$elpd_pointwise_cross}.
#' @param alpha (optional) Threshold for Bayes-FDR (LFSR). If provided,
#'   adds a \code{pass_bayes_fdr} column.
#' @param diag_mode Either \code{"moderate"} or \code{"strict"}.
#' @param use_true_stacking Logical; if \code{TRUE} (default) attempts to use
#'   \code{.compute_stacking_weights_cross()} when feasible (otherwise omitted).
#' @param stacking_use_ppd Logical; pass-through to \code{.compute_stacking_weights_cross()}.
#' @param stacking_min_models Integer; pass-through.
#' @param stacking_min_subjects Integer; pass-through.
#'
#' @param interaction Which interaction summary to return. One of
#'   \code{"cross"} or \code{"self"}; default \code{"cross"}.
#'
#' @return A list containing the requested table (\code{$cross} or \code{$self}).
#'   The table is sorted by \strong{bayes_FDR} (ascending; \code{NA} last).
#'
#' @examples
#'
#' @seealso \code{fit_glv_pairwise()}
#' @export
summarize_bayes_pclv <- function(df,
                                 alpha = NULL,
                                 diag_mode = c("moderate","strict"),
                                 use_true_stacking = TRUE,
                                 stacking_use_ppd = FALSE,
                                 stacking_min_models = 2L,
                                 stacking_min_subjects = 1L,
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

  # diagnostics per direction from 'raw'
  cross_diag <- dplyr::bind_rows(
    dplyr::transmute(
      raw,
      from = .data$j, to = .data$i,
      n_pairs = .data$n_pairs_ij,
      rhat = .data$rhat_ij, essb = .data$essb_ij, esst = .data$esst_ij,
      div = .data$div_ij, tdhit = .data$tdhit_ij,
      ebfmi_min = .data$ebfmi_min_ij
    ),
    dplyr::transmute(
      raw,
      from = .data$i, to = .data$j,
      n_pairs = .data$n_pairs_ji,
      rhat = .data$rhat_ji, essb = .data$essb_ji, esst = .data$esst_ji,
      div = .data$div_ji, tdhit = .data$tdhit_ji,
      ebfmi_min = .data$ebfmi_min_ji
    )
  )
  self_diag <- dplyr::bind_rows(
    dplyr::transmute(
      raw,
      taxon = .data$i,
      n_pairs = .data$n_pairs_ij,
      rhat = .data$rhat_ij, essb = .data$essb_ij, esst = .data$esst_ij,
      div = .data$div_ij, tdhit = .data$tdhit_ij,
      ebfmi_min = .data$ebfmi_min_ij
    ),
    dplyr::transmute(
      raw,
      taxon = .data$j,
      n_pairs = .data$n_pairs_ji,
      rhat = .data$rhat_ji, essb = .data$essb_ji, esst = .data$esst_ji,
      div = .data$div_ji, tdhit = .data$tdhit_ji,
      ebfmi_min = .data$ebfmi_min_ji
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
        diag_ok = .diag_ok_fun(.data$rhat, .data$essb, .data$esst,
                               .data$div, .data$tdhit, .data$ebfmi_min, thr = thr),
        a_sign = dplyr::case_when(
          is.finite(.data$a_mean) & .data$a_mean >  0 ~ "+",
          is.finite(.data$a_mean) & .data$a_mean <  0 ~ "-",
          TRUE ~ "0"
        ),
        bayes_FDR = ifelse(.data$diag_ok, .lfsr_safe(.data$p_sign2), NA_real_)
      )

    # 2) diag_ok==TRUE 엣지 목록
    ok_edges <- cross_merged |>
      dplyr::filter(.data$diag_ok) |>
      dplyr::distinct(.data$from, .data$to)

    # 3) pseudo-BMA (BB=FALSE) / pseudo-BMA+ (BB=TRUE) 계산 (pointwise elpd 기반)
    w_pbma <- .compute_pseudobma_weights_cross(
      df, edges = ok_edges, use_ppd = stacking_use_ppd,
      plus = FALSE, min_models = stacking_min_models,
      min_subjects = stacking_min_subjects
    )
    w_pbma_plus <- .compute_pseudobma_weights_cross(
      df, edges = ok_edges, use_ppd = stacking_use_ppd,
      plus = TRUE, min_models = stacking_min_models,
      min_subjects = stacking_min_subjects
    )

    cross_w <- cross_merged |>
      dplyr::left_join(dplyr::rename(w_pbma,  pseudo_BMA      = .data$weight), by = c("from","to")) |>
      dplyr::left_join(dplyr::rename(w_pbma_plus, pseudo_BMA_plus = .data$weight), by = c("from","to"))


    cross_sum <- cross_w |>
      dplyr::transmute(
        from, to,
        n_subjects = .data$n_subjects,
        n_pairs = .data$n_pairs,
        a_sign, a_mean, a_q2.5, a_q97.5,
        p_sign2, bayes_FDR,
        pseudo_BMA,
        pseudo_BMA_plus,
        stacking = NA_real_,
        rhat, essb, esst, div, tdhit, ebfmi_min,
        diag_ok
      )

    # 4) TRUE stacking도 diag_ok==TRUE 엣지만 사용
    if (isTRUE(use_true_stacking) &&
        ("elpd_pointwise_cross" %in% names(df)) &&
        !is.null(df$elpd_pointwise_cross) &&
        nrow(df$elpd_pointwise_cross)) {
      if (nrow(ok_edges)) {
        pw_ok <- tibble::as_tibble(df$elpd_pointwise_cross) |>
          dplyr::semi_join(ok_edges, by = c("from","to"))
        if (nrow(pw_ok)) {
          df2 <- df
          df2$elpd_pointwise_cross <- pw_ok
          sw <- try(
            .compute_stacking_weights_cross(
              df2,
              use_ppd = stacking_use_ppd,
              min_models = stacking_min_models,
              min_subjects = stacking_min_subjects
            ),
            silent = TRUE
          )
          if (!inherits(sw, "try-error") && is.data.frame(sw) && nrow(sw)) {
            cross_sum <- dplyr::left_join(
              cross_sum,
              dplyr::rename(sw, stacking_true = .data$stacking),
              by = c("from","to")
            )
            cross_sum$stacking <- ifelse(
              is.finite(cross_sum$stacking_true),
              cross_sum$stacking_true,
              cross_sum$stacking
            )
            cross_sum$stacking_true <- NULL
          }
        }
      }
    }

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
        diag_ok = .diag_ok_fun(.data$rhat, .data$essb, .data$esst,
                               .data$div, .data$tdhit, .data$ebfmi_min, thr = thr),
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
        pseudo_BMA = NA_real_,
        pseudo_BMA_plus = NA_real_,
        stacking = NA_real_,
        rhat, essb, esst, div, tdhit, ebfmi_min,
        diag_ok
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
