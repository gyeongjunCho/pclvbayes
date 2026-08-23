#' Summarize Bayesian pcLV results
#'
#' @description
#' Given the result list from \code{fit_pclv_bayes()}, returns a summary table
#' for directed cross effects or pair-context self effects.
#'
#' Directed cross coefficients are pair-to-rest dynamic coefficients and are
#' not generally absolute direct-gLV effects. Posterior sign support, MCMC
#' diagnostics, and pair-specific held-out predictive evidence are retained as
#' distinct local evidence fingerprints.
#'
#' When \code{df$elpd_pointwise_cross} is available, repeated-K-fold ELPD is
#' summarized per directed pair and joined to the cross table. These ELPD
#' summaries describe predictive adequacy for that pair-specific transformed
#' outcome. They are not cross-edge model probabilities or common-outcome
#' model weights.
#'
#' Cross-pair stacking, pseudo-BMA, and pseudo-BMA+ weights are deliberately
#' not returned because distinct pair-to-rest models generally predict
#' different transformed outcomes and therefore do not form a common-outcome
#' model set.
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
#' Missing diagnostic values fail closed. Results without a diagnostic class
#' are marked \code{"diagnostics_unavailable"} unless a historical
#' \code{direction_ok_*} field explicitly establishes their status.
#'
#' Predictive evidence is optional. If an edge has no valid predictive score,
#' predictive summary fields remain missing and \code{predictive_status}
#' records whether evidence was unavailable, partial, or invalid/failed.
#'
#' Self coefficients are pair-context estimands. The returned self table keeps
#' both \code{taxon} and \code{partner}, preventing unrelated pair-specific
#' self estimates and diagnostics from being joined by taxon alone.
#'
#' Indeterminate directions remain explicit and are not interpreted as zero.
#' Non-finite coefficient means receive a missing sign rather than \code{"0"}.
#' Tables are sorted by \code{bayes_FDR} in ascending order with missing values
#' last.
#'
#' @param df A list returned by \code{fit_pclv_bayes()} containing
#'   \code{$cross}, \code{$self}, and \code{$raw}. It may additionally
#'   contain \code{$elpd_pointwise_cross}.
#' @param alpha Optional threshold for Bayes-FDR/LFSR. When supplied, adds
#'   \code{pass_bayes_fdr}.
#' @param diag_mode Either \code{"moderate"} or \code{"strict"}.
#' @param interaction One of \code{"cross"} or \code{"self"}.
#'
#' @return A list containing the requested \code{$cross} or \code{$self}
#'   summary table. The cross table includes pair-specific predictive
#'   fingerprint fields when available. The self table includes \code{taxon}
#'   and \code{partner}.
#'
#' @export
summarize_bayes_pclv <- function(df,
                                 alpha = NULL,
                                 diag_mode = c("moderate", "strict"),
                                 interaction = "cross") {
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("Packages 'dplyr' and 'tibble' are required.")
  }

  diag_mode <- match.arg(diag_mode, c("moderate", "strict"))
  interaction <- match.arg(interaction, c("cross", "self"))
  do_cross <- identical(interaction, "cross")
  do_self <- identical(interaction, "self")

  if (!is.null(alpha)) {
    if (!is.numeric(alpha) ||
        length(alpha) != 1L ||
        !is.finite(alpha) ||
        alpha < 0 ||
        alpha > 1) {
      stop("`alpha` must be NULL or a finite scalar in [0, 1].")
    }
    alpha <- as.numeric(alpha)
  }

  thr <- switch(
    diag_mode,
    strict = list(
      rhat = 1.01,
      ess = 1000,
      div = 0L,
      tdhit = 0L,
      ebfmi = 0.40
    ),
    moderate = list(
      rhat = 1.05,
      ess = 400,
      div = 8L,
      tdhit = 80L,
      ebfmi = 0.30
    )
  )

  if (!is.list(df) ||
      !all(c("cross", "self", "raw") %in% names(df))) {
    stop("`df` must include $cross, $self, and $raw.")
  }

  cross <- tibble::as_tibble(df$cross)
  self <- tibble::as_tibble(df$self)
  raw <- tibble::as_tibble(df$raw)

  # Validate only the columns required by the requested summary. The raw table
  # supplies directional diagnostics and pair identity; self posterior
  # summaries are taken from `df$self`, preserving compatibility with
  # historical/test fixtures whose raw table does not repeat those columns.
  required_raw <- c(
    "i", "j",
    "n_pairs_ij", "n_pairs_ji",
    "rhat_ij", "essb_ij", "esst_ij", "div_ij", "tdhit_ij",
    "ebfmi_min_ij",
    "rhat_ji", "essb_ji", "esst_ji", "div_ji", "tdhit_ji",
    "ebfmi_min_ji"
  )

  missing_raw <- setdiff(required_raw, names(raw))
  if (length(missing_raw)) {
    stop(
      sprintf(
        "`df$raw` lacks required column(s): %s.",
        paste(missing_raw, collapse = ", ")
      )
    )
  }

  if (do_cross) {
    required_cross <- c(
      "from", "to", "n_subjects",
      "a_mean", "a_q2.5", "a_q97.5", "p_sign2"
    )
    missing_cross <- setdiff(required_cross, names(cross))
    if (length(missing_cross)) {
      stop(
        sprintf(
          "`df$cross` lacks required column(s): %s.",
          paste(missing_cross, collapse = ", ")
        )
      )
    }
  }

  if (do_self) {
    required_self <- c(
      "taxon",
      "a_self_mean", "a_self_q2.5", "a_self_q97.5",
      "p_sign2_self"
    )
    missing_self <- setdiff(required_self, names(self))
    if (length(missing_self)) {
      stop(
        sprintf(
          "`df$self` lacks required column(s): %s.",
          paste(missing_self, collapse = ", ")
        )
      )
    }
  }

  # Historical results may lack chain-specific diagnostic classes. They are
  # accepted only when the corresponding historical direction_ok field exists.
  for (suffix in c("ij", "ji")) {
    class_col <- paste0("diagnostic_class_", suffix)
    ok_col <- paste0("direction_ok_", suffix)

    if (!class_col %in% names(raw)) {
      if (ok_col %in% names(raw)) {
        ok_value <- raw[[ok_col]]
        raw[[class_col]] <- ifelse(
          is.na(ok_value),
          "diagnostics_unavailable",
          ifelse(
            ok_value,
            "converged",
            "sampler_diagnostics_failed"
          )
        )
      } else {
        raw[[class_col]] <- rep(
          "diagnostics_unavailable",
          nrow(raw)
        )
      }
    }

    raw[[class_col]] <- as.character(raw[[class_col]])
    raw[[class_col]][
      is.na(raw[[class_col]]) | !nzchar(raw[[class_col]])
    ] <- "diagnostics_unavailable"
  }

  if (!"diagnostic_class" %in% names(cross)) {
    cross$diagnostic_class <- NA_character_
  } else {
    cross$diagnostic_class <- as.character(cross$diagnostic_class)
  }

  cross_diag <- dplyr::bind_rows(
    dplyr::transmute(
      raw,
      from = as.character(.data$j),
      to = as.character(.data$i),
      n_pairs = as.integer(.data$n_pairs_ij),
      rhat = .data$rhat_ij,
      essb = .data$essb_ij,
      esst = .data$esst_ij,
      div = .data$div_ij,
      tdhit = .data$tdhit_ij,
      ebfmi_min = .data$ebfmi_min_ij,
      diagnostic_class_diag = .data$diagnostic_class_ij
    ),
    dplyr::transmute(
      raw,
      from = as.character(.data$i),
      to = as.character(.data$j),
      n_pairs = as.integer(.data$n_pairs_ji),
      rhat = .data$rhat_ji,
      essb = .data$essb_ji,
      esst = .data$esst_ji,
      div = .data$div_ji,
      tdhit = .data$tdhit_ji,
      ebfmi_min = .data$ebfmi_min_ji,
      diagnostic_class_diag = .data$diagnostic_class_ji
    )
  )

  cross_key <- paste(cross_diag$from, cross_diag$to, sep = "\r")
  if (anyDuplicated(cross_key)) {
    stop(
      "`df$raw` contains duplicate directed identities; cross diagnostics cannot be joined one-to-one."
    )
  }

  if (do_cross) {
    n_cross_before <- nrow(cross)

    cross_merged <- cross |>
      dplyr::left_join(
        cross_diag,
        by = c("from", "to")
      )

    if (nrow(cross_merged) != n_cross_before) {
      stop(
        "Cross diagnostic join changed the row count; directed identities are not one-to-one."
      )
    }

    # Pair-specific predictive evidence is summarized independently of the
    # diagnostic gate. This preserves the predictive fingerprint even when a
    # direction later fails MCMC reportability.
    pred_cross <- .summarize_elpd_cross(df)

    if (nrow(pred_cross)) {
      pred_key <- paste(pred_cross$from, pred_cross$to, sep = "\r")
      if (anyDuplicated(pred_key)) {
        stop(
          "Pair-specific predictive evidence contains duplicate directed identities."
        )
      }

      cross_merged <- cross_merged |>
        dplyr::left_join(
          pred_cross,
          by = c("from", "to")
        )

      if (nrow(cross_merged) != n_cross_before) {
        stop(
          "Predictive-evidence join changed the row count; directed identities are not one-to-one."
        )
      }
    }

    # Typed defaults keep the public schema stable when predictive evidence is
    # globally absent or absent for a particular edge.
    predictive_numeric <- c(
      "elpd_total",
      "elpd_per_test",
      "elpd_subject_mean",
      "elpd_subject_sd",
      "elpd_ppd_mean",
      "elpd_split_sd",
      "n_test_total",
      "predictive_success_fraction"
    )
    predictive_integer <- c(
      "n_predictive_subjects",
      "n_predictive_rows",
      "n_successful_splits"
    )

    for (nm in predictive_numeric) {
      if (!nm %in% names(cross_merged)) {
        cross_merged[[nm]] <- NA_real_
      }
    }
    for (nm in predictive_integer) {
      if (!nm %in% names(cross_merged)) {
        cross_merged[[nm]] <- NA_integer_
      }
    }
    if (!"predictive_available" %in% names(cross_merged)) {
      cross_merged$predictive_available <- FALSE
    }
    if (!"predictive_status" %in% names(cross_merged)) {
      cross_merged$predictive_status <- "predictive_unavailable"
    }

    cross_merged$predictive_available <- dplyr::coalesce(
      as.logical(cross_merged$predictive_available),
      FALSE
    )
    cross_merged$predictive_status <- dplyr::coalesce(
      as.character(cross_merged$predictive_status),
      "predictive_unavailable"
    )

    cross_sum <- cross_merged |>
      dplyr::mutate(
        diagnostic_class = dplyr::coalesce(
          .data$diagnostic_class,
          .data$diagnostic_class_diag,
          "diagnostics_unavailable"
        ),
        sampler_diag_ok = .diag_ok_fun(
          .data$rhat,
          .data$essb,
          .data$esst,
          .data$div,
          .data$tdhit,
          .data$ebfmi_min,
          thr = thr
        ),
        diag_ok =
          .data$sampler_diag_ok &
          .data$diagnostic_class == "converged",
        a_sign = dplyr::case_when(
          !is.finite(.data$a_mean) ~ NA_character_,
          .data$a_mean > 0 ~ "+",
          .data$a_mean < 0 ~ "-",
          TRUE ~ "0"
        ),
        bayes_FDR = ifelse(
          .data$diag_ok,
          .lfsr_safe(.data$p_sign2),
          NA_real_
        )
      ) |>
      dplyr::transmute(
        from = as.character(.data$from),
        to = as.character(.data$to),
        n_subjects = .data$n_subjects,
        n_pairs = .data$n_pairs,
        a_sign = .data$a_sign,
        a_mean = .data$a_mean,
        a_q2.5 = .data$a_q2.5,
        a_q97.5 = .data$a_q97.5,
        p_sign2 = .data$p_sign2,
        bayes_FDR = .data$bayes_FDR,
        rhat = .data$rhat,
        essb = .data$essb,
        esst = .data$esst,
        div = .data$div,
        tdhit = .data$tdhit,
        ebfmi_min = .data$ebfmi_min,
        diagnostic_class = .data$diagnostic_class,
        sampler_diag_ok = .data$sampler_diag_ok,
        diag_ok = .data$diag_ok,
        elpd_total = .data$elpd_total,
        elpd_per_test = .data$elpd_per_test,
        elpd_subject_mean = .data$elpd_subject_mean,
        elpd_subject_sd = .data$elpd_subject_sd,
        elpd_ppd_mean = .data$elpd_ppd_mean,
        elpd_split_sd = .data$elpd_split_sd,
        n_test_total = .data$n_test_total,
        n_predictive_subjects = .data$n_predictive_subjects,
        n_predictive_rows = .data$n_predictive_rows,
        n_successful_splits = .data$n_successful_splits,
        predictive_success_fraction = .data$predictive_success_fraction,
        predictive_available = .data$predictive_available,
        predictive_status = .data$predictive_status
      )

    if (!is.null(alpha)) {
      cross_sum <- dplyr::mutate(
        cross_sum,
        pass_bayes_fdr =
          is.finite(.data$bayes_FDR) &
          .data$bayes_FDR <= alpha
      )
    }

    cross_sum <- dplyr::arrange(
      cross_sum,
      .data$bayes_FDR
    )
  }

  if (do_self) {
    # New results carry `partner` directly in df$self. Historical package
    # results omitted it, but .mk_self() emitted rows deterministically as all
    # i-side rows followed by all j-side rows. Reconstruct partner only when
    # that exact invariant is verifiable; otherwise fail rather than performing
    # an ambiguous taxon-only join.
    if (!"partner" %in% names(self)) {
      expected_n <- 2L * nrow(raw)

      if (nrow(self) != expected_n) {
        stop(
          paste(
            "`df$self` lacks `partner`, and its row count does not match",
            "the historical two-rows-per-pair layout; self contexts are",
            "ambiguous."
          )
        )
      }

      expected_taxon <- c(
        as.character(raw$i),
        as.character(raw$j)
      )
      reconstructed_partner <- c(
        as.character(raw$j),
        as.character(raw$i)
      )

      if (!identical(as.character(self$taxon), expected_taxon)) {
        stop(
          paste(
            "`df$self` lacks `partner`, and its row order does not match",
            "the historical .mk_self() layout; self contexts cannot be",
            "reconstructed safely."
          )
        )
      }

      self$partner <- reconstructed_partner
    }

    self$taxon <- as.character(self$taxon)
    self$partner <- as.character(self$partner)

    if (anyNA(self$taxon) ||
        anyNA(self$partner) ||
        any(!nzchar(self$taxon)) ||
        any(!nzchar(self$partner))) {
      stop(
        "`df$self` contains missing or empty taxon-partner identities."
      )
    }

    self_diag <- dplyr::bind_rows(
      dplyr::transmute(
        raw,
        taxon = as.character(.data$i),
        partner = as.character(.data$j),
        n_pairs = as.integer(.data$n_pairs_ij),
        rhat = .data$rhat_ij,
        essb = .data$essb_ij,
        esst = .data$esst_ij,
        div = .data$div_ij,
        tdhit = .data$tdhit_ij,
        ebfmi_min = .data$ebfmi_min_ij,
        diagnostic_class = .data$diagnostic_class_ij
      ),
      dplyr::transmute(
        raw,
        taxon = as.character(.data$j),
        partner = as.character(.data$i),
        n_pairs = as.integer(.data$n_pairs_ji),
        rhat = .data$rhat_ji,
        essb = .data$essb_ji,
        esst = .data$esst_ji,
        div = .data$div_ji,
        tdhit = .data$tdhit_ji,
        ebfmi_min = .data$ebfmi_min_ji,
        diagnostic_class = .data$diagnostic_class_ji
      )
    )

    self_key <- paste(
      self$taxon,
      self$partner,
      sep = "\r"
    )
    self_diag_key <- paste(
      self_diag$taxon,
      self_diag$partner,
      sep = "\r"
    )

    if (anyDuplicated(self_key)) {
      stop(
        "`df$self` contains duplicate taxon-partner self identities."
      )
    }
    if (anyDuplicated(self_diag_key)) {
      stop(
        "`df$raw` contains duplicate taxon-partner diagnostic identities."
      )
    }

    n_self_before <- nrow(self)

    self_merged <- self |>
      dplyr::left_join(
        self_diag,
        by = c("taxon", "partner")
      )

    if (nrow(self_merged) != n_self_before) {
      stop(
        "Self diagnostic join changed the row count; identities are not one-to-one."
      )
    }

    self_sum <- self_merged |>
      dplyr::mutate(
        diagnostic_class = dplyr::coalesce(
          as.character(.data$diagnostic_class),
          "diagnostics_unavailable"
        ),
        sampler_diag_ok = .diag_ok_fun(
          .data$rhat,
          .data$essb,
          .data$esst,
          .data$div,
          .data$tdhit,
          .data$ebfmi_min,
          thr = thr
        ),
        diag_ok =
          .data$sampler_diag_ok &
          .data$diagnostic_class == "converged",
        from = .data$taxon,
        to = .data$taxon,
        a_sign = dplyr::case_when(
          !is.finite(.data$a_self_mean) ~ NA_character_,
          .data$a_self_mean > 0 ~ "+",
          .data$a_self_mean < 0 ~ "-",
          TRUE ~ "0"
        ),
        bayes_FDR = ifelse(
          .data$diag_ok,
          .lfsr_safe(.data$p_sign2_self),
          NA_real_
        )
      ) |>
      dplyr::transmute(
        taxon = .data$taxon,
        partner = .data$partner,
        from = .data$from,
        to = .data$to,
        n_pairs = .data$n_pairs,
        a_sign = .data$a_sign,
        a_mean = .data$a_self_mean,
        a_q2.5 = .data$a_self_q2.5,
        a_q97.5 = .data$a_self_q97.5,
        p_sign2 = .data$p_sign2_self,
        bayes_FDR = .data$bayes_FDR,
        rhat = .data$rhat,
        essb = .data$essb,
        esst = .data$esst,
        div = .data$div,
        tdhit = .data$tdhit,
        ebfmi_min = .data$ebfmi_min,
        diagnostic_class = .data$diagnostic_class,
        sampler_diag_ok = .data$sampler_diag_ok,
        diag_ok = .data$diag_ok
      )

    if (!is.null(alpha)) {
      self_sum <- dplyr::mutate(
        self_sum,
        pass_bayes_fdr =
          is.finite(.data$bayes_FDR) &
          .data$bayes_FDR <= alpha
      )
    }

    self_sum <- dplyr::arrange(
      self_sum,
      .data$bayes_FDR,
      .data$taxon,
      .data$partner
    )
  }

  out <- list()
  if (do_cross) {
    out$cross <- cross_sum
  }
  if (do_self) {
    out$self <- self_sum
  }
  out
}
