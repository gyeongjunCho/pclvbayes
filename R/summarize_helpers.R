#' @noRd
#' @keywords internal
.as_int <- function(x) suppressWarnings(as.integer(x))

#' @noRd
#' @keywords internal
.as_num <- function(x) suppressWarnings(as.numeric(x))

#' @noRd
#' @keywords internal
.lfsr_safe <- function(p_two) {
  if (exists(".lfsr_from_two", mode = "function", inherits = TRUE)) {
    return(.lfsr_from_two(p_two))
  }
  pmin(pmax(p_two / 2, 0), 0.5)
}

#' Reconstruct between-chain directional stability from stored medians
#'
#' Returns TRUE when all available chains have finite, non-zero posterior
#' medians with the same sign. Returns NA when chain information is unavailable.
#'
#' @noRd
#' @keywords internal
.chain_direction_stable_safe <- function(medians) {
  x <- suppressWarnings(
    as.numeric(unlist(medians, use.names = FALSE))
  )

  if (!length(x) || any(!is.finite(x))) {
    return(NA)
  }

  s <- sign(x)
  if (any(s == 0)) {
    return(FALSE)
  }

  length(unique(s)) == 1L
}

#' Descriptive all-chain posterior sign-confidence flag
#'
#' This is not an interaction-identifiability gate. It records whether every
#' chain puts at least `threshold` posterior mass on one sign. Cross effects
#' provide separate positive and negative probabilities; self effects may
#' provide the already-maximized dominant-sign probability.
#'
#' @noRd
#' @keywords internal
.all_chain_sign_confident_safe <- function(positive = NULL,
                                            negative = NULL,
                                            dominant = NULL,
                                            threshold = 0.95) {
  if (!is.null(dominant)) {
    p <- suppressWarnings(
      as.numeric(unlist(dominant, use.names = FALSE))
    )

    if (!length(p) || any(!is.finite(p))) {
      return(NA)
    }

    return(all(p >= threshold))
  }

  pp <- suppressWarnings(
    as.numeric(unlist(positive, use.names = FALSE))
  )
  pn <- suppressWarnings(
    as.numeric(unlist(negative, use.names = FALSE))
  )

  if (!length(pp) ||
      length(pp) != length(pn) ||
      any(!is.finite(pp)) ||
      any(!is.finite(pn))) {
    return(NA)
  }

  all(pmax(pp, pn) >= threshold)
}

#' Element-wise diagnostic pass/fail predicate used by the public summarizer
#'
#' Missing, non-finite, negative, non-integral, or malformed diagnostic counts
#' fail closed. A diagnostic is never treated as passing merely because its
#' divergence or treedepth count is unavailable.
#'
#' @noRd
#' @keywords internal
.diag_ok_fun <- function(rhat, essb, esst, div, tdhit, ebfmi_min, thr) {
  rhat_num <- .as_num(rhat)
  essb_num <- .as_num(essb)
  esst_num <- .as_num(esst)
  div_num <- .as_num(div)
  tdhit_num <- .as_num(tdhit)
  ebfmi_num <- .as_num(ebfmi_min)

  div_int <- .as_int(div_num)
  tdhit_int <- .as_int(tdhit_num)

  valid_div <-
    is.finite(div_num) &
    !is.na(div_int) &
    div_num >= 0 &
    div_num == div_int &
    div_int <= thr$div

  valid_tdhit <-
    is.finite(tdhit_num) &
    !is.na(tdhit_int) &
    tdhit_num >= 0 &
    tdhit_num == tdhit_int &
    tdhit_int <= thr$tdhit

  (is.finite(rhat_num) & rhat_num < thr$rhat) &
    (is.finite(essb_num) & essb_num > thr$ess) &
    (is.finite(esst_num) & esst_num > thr$ess) &
    valid_div &
    valid_tdhit &
    (is.finite(ebfmi_num) & ebfmi_num >= thr$ebfmi)
}

#' Empty pair-specific predictive-evidence summary
#'
#' @noRd
#' @keywords internal
.empty_elpd_cross_summary <- function() {
  tibble::tibble(
    from = character(),
    to = character(),
    elpd_total = numeric(),
    elpd_per_test = numeric(),
    elpd_subject_mean = numeric(),
    elpd_subject_sd = numeric(),
    elpd_ppd_mean = numeric(),
    elpd_split_sd = numeric(),
    n_test_total = numeric(),
    n_predictive_subjects = integer(),
    n_predictive_rows = integer(),
    n_successful_splits = integer(),
    predictive_success_fraction = numeric(),
    predictive_available = logical(),
    predictive_status = character()
  )
}

#' Summarize pair-specific repeated-K-fold predictive evidence
#'
#' @description
#' Summarizes \code{df$elpd_pointwise_cross} for each directed pair without
#' comparing distinct pair-to-rest models as if they predicted a common
#' outcome. The resulting quantities are predictive fingerprints for the same
#' pair-specific transformed outcome, not stacking or pseudo-BMA model weights.
#'
#' Repeated scores for the same subject are averaged first so that increasing
#' the number of repeated CV splits does not mechanically multiply
#' \code{elpd_total}. Thus \code{elpd_total} is the sum of subject-level mean
#' held-out ELPD values, and \code{elpd_per_test} divides that quantity by the
#' corresponding sum of subject-level mean test counts.
#'
#' If \code{fold} and/or \code{repeat} are present, \code{elpd_split_sd} is the
#' standard deviation of split-level ELPD per test observation. Otherwise it is
#' missing. Rows with non-finite ELPD or non-positive/non-finite \code{n_test}
#' are retained in the availability accounting but excluded from numerical
#' summaries.
#'
#' @param df A result list that may contain \code{$elpd_pointwise_cross}.
#'
#' @return A tibble with one row per directed pair represented in
#'   \code{$elpd_pointwise_cross}. Returns an empty typed tibble when predictive
#'   evidence is absent.
#'
#' @noRd
#' @keywords internal
.summarize_elpd_cross <- function(df) {
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("Packages 'dplyr' and 'tibble' are required.")
  }

  if (!is.list(df) ||
      !("elpd_pointwise_cross" %in% names(df)) ||
      is.null(df$elpd_pointwise_cross) ||
      !NROW(df$elpd_pointwise_cross)) {
    return(.empty_elpd_cross_summary())
  }

  pw <- tibble::as_tibble(df$elpd_pointwise_cross)

  required <- c("from", "to", "subject", "elpd", "n_test")
  missing_required <- setdiff(required, names(pw))
  if (length(missing_required)) {
    stop(
      sprintf(
        "`df$elpd_pointwise_cross` lacks required column(s): %s.",
        paste(missing_required, collapse = ", ")
      )
    )
  }

  pw <- pw |>
    dplyr::mutate(
      from = as.character(.data$from),
      to = as.character(.data$to),
      subject = as.character(.data$subject),
      elpd_num = .as_num(.data$elpd),
      n_test_num = .as_num(.data$n_test),
      valid_predictive =
        is.finite(.data$elpd_num) &
        is.finite(.data$n_test_num) &
        .data$n_test_num > 0
    )

  if (anyNA(pw$from) ||
      anyNA(pw$to) ||
      any(!nzchar(pw$from)) ||
      any(!nzchar(pw$to))) {
    stop(
      "`df$elpd_pointwise_cross` contains missing or empty directed identities."
    )
  }

  # Missing subject labels cannot define subject-level predictive units. Such
  # rows remain part of failure accounting but are not numerically summarized.
  valid_subject <-
    !is.na(pw$subject) &
    nzchar(pw$subject)

  pw$valid_predictive <-
    pw$valid_predictive &
    valid_subject

  if ("elpd_ppd" %in% names(pw)) {
    pw$elpd_ppd_num <- .as_num(pw$elpd_ppd)
  } else {
    pw$elpd_ppd_num <- NA_real_
  }

  # Preserve every pair represented in the predictive table, including pairs
  # for which all predictive rows failed.
  pair_status <- pw |>
    dplyr::group_by(.data$from, .data$to) |>
    dplyr::summarise(
      n_rows_all = dplyr::n(),
      n_rows_valid = sum(.data$valid_predictive),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      predictive_success_fraction =
        .data$n_rows_valid / .data$n_rows_all,
      predictive_available = .data$n_rows_valid > 0L,
      predictive_status = dplyr::case_when(
        .data$n_rows_valid <= 0L ~ "predictive_invalid_or_failed",
        .data$n_rows_valid < .data$n_rows_all ~ "predictive_partial",
        TRUE ~ "predictive_available"
      )
    )

  valid <- pw |>
    dplyr::filter(.data$valid_predictive)

  if (!nrow(valid)) {
    return(
      pair_status |>
        dplyr::transmute(
          from = .data$from,
          to = .data$to,
          elpd_total = NA_real_,
          elpd_per_test = NA_real_,
          elpd_subject_mean = NA_real_,
          elpd_subject_sd = NA_real_,
          elpd_ppd_mean = NA_real_,
          elpd_split_sd = NA_real_,
          n_test_total = NA_real_,
          n_predictive_subjects = 0L,
          n_predictive_rows = 0L,
          n_successful_splits = NA_integer_,
          predictive_success_fraction =
            .data$predictive_success_fraction,
          predictive_available =
            .data$predictive_available,
          predictive_status =
            .data$predictive_status
        )
    )
  }

  # First collapse repeated CV scores for the same subject. This prevents R
  # repeats from multiplying the apparent amount of predictive information.
  subject_summary <- valid |>
    dplyr::group_by(.data$from, .data$to, .data$subject) |>
    dplyr::summarise(
      elpd_subject = mean(.data$elpd_num),
      n_test_subject = mean(.data$n_test_num),
      elpd_ppd_subject = if (any(is.finite(.data$elpd_ppd_num))) {
        mean(.data$elpd_ppd_num[is.finite(.data$elpd_ppd_num)])
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  pair_numeric <- subject_summary |>
    dplyr::group_by(.data$from, .data$to) |>
    dplyr::summarise(
      elpd_total = sum(.data$elpd_subject),
      n_test_total = sum(.data$n_test_subject),
      elpd_subject_mean = mean(.data$elpd_subject),
      elpd_subject_sd = if (dplyr::n() > 1L) {
        stats::sd(.data$elpd_subject)
      } else {
        NA_real_
      },
      elpd_ppd_mean = if (any(is.finite(.data$elpd_ppd_subject))) {
        mean(.data$elpd_ppd_subject[is.finite(.data$elpd_ppd_subject)])
      } else {
        NA_real_
      },
      n_predictive_subjects = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      elpd_per_test = dplyr::if_else(
        is.finite(.data$n_test_total) & .data$n_test_total > 0,
        .data$elpd_total / .data$n_test_total,
        NA_real_
      )
    )

  row_counts <- valid |>
    dplyr::group_by(.data$from, .data$to) |>
    dplyr::summarise(
      n_predictive_rows = dplyr::n(),
      .groups = "drop"
    )

  split_cols <- intersect(c("repeat", "fold"), names(valid))
  if (length(split_cols)) {
    split_id <- do.call(
      paste,
      c(
        lapply(
          split_cols,
          function(nm) {
            x <- as.character(valid[[nm]])
            x[is.na(x) | !nzchar(x)] <- "<NA>"
            x
          }
        ),
        sep = "\r"
      )
    )
    valid$split_id <- split_id

    split_summary <- valid |>
      dplyr::group_by(.data$from, .data$to, .data$split_id) |>
      dplyr::summarise(
        split_elpd = sum(.data$elpd_num),
        split_n_test = sum(.data$n_test_num),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        split_elpd_per_test = dplyr::if_else(
          is.finite(.data$split_n_test) & .data$split_n_test > 0,
          .data$split_elpd / .data$split_n_test,
          NA_real_
        )
      )

    split_pair <- split_summary |>
      dplyr::group_by(.data$from, .data$to) |>
      dplyr::summarise(
        elpd_split_sd = {
          z <- .data$split_elpd_per_test[
            is.finite(.data$split_elpd_per_test)
          ]
          if (length(z) > 1L) stats::sd(z) else NA_real_
        },
        n_successful_splits = sum(
          is.finite(.data$split_elpd_per_test)
        ),
        .groups = "drop"
      )
  } else {
    split_pair <- pair_numeric |>
      dplyr::transmute(
        from = .data$from,
        to = .data$to,
        elpd_split_sd = NA_real_,
        n_successful_splits = NA_integer_
      )
  }

  out <- pair_status |>
    dplyr::left_join(
      pair_numeric,
      by = c("from", "to")
    ) |>
    dplyr::left_join(
      row_counts,
      by = c("from", "to")
    ) |>
    dplyr::left_join(
      split_pair,
      by = c("from", "to")
    ) |>
    dplyr::transmute(
      from = .data$from,
      to = .data$to,
      elpd_total = .data$elpd_total,
      elpd_per_test = .data$elpd_per_test,
      elpd_subject_mean = .data$elpd_subject_mean,
      elpd_subject_sd = .data$elpd_subject_sd,
      elpd_ppd_mean = .data$elpd_ppd_mean,
      elpd_split_sd = .data$elpd_split_sd,
      n_test_total = .data$n_test_total,
      n_predictive_subjects = dplyr::coalesce(
        as.integer(.data$n_predictive_subjects),
        0L
      ),
      n_predictive_rows = dplyr::coalesce(
        as.integer(.data$n_predictive_rows),
        0L
      ),
      n_successful_splits = as.integer(.data$n_successful_splits),
      predictive_success_fraction =
        .data$predictive_success_fraction,
      predictive_available =
        .data$predictive_available,
      predictive_status =
        .data$predictive_status
    )

  key <- paste(out$from, out$to, sep = "\r")
  if (anyDuplicated(key)) {
    stop(
      "Predictive-evidence summarization produced duplicate directed identities."
    )
  }

  out
}
