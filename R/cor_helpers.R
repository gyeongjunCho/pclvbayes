# Correlation/residual helper functions -----------------------------------------
#
# Architecture note (2026-08): observed-positive support and triplet ALR
# construction live in compositional_helpers.R so correlation, future OU, and
# eventually fit_pclv_bayes() can share one compositional policy. This file owns
# correlation-specific time aggregation and discrete-lag effective-n utilities.
# A future OU implementation should NOT extend .eff_n() into a pseudo-continuous
# correction; it should use actual dt in the OU model while reusing the shared
# observed-support and triplet helpers.

#' Effective sample size with autocorrelation correction
#'
#' Methods:
#' - \code{effn_method = "ar1"}: AR(1) approximation
#'   \code{n * (1 - rho_x(1) rho_y(1)) / (1 + rho_x(1) rho_y(1))}.
#' - \code{effn_method = "bartlett"}: Bartlett-window sum of lagged ACF products.
#' - \code{effn_method = "nw"}: Newey-West automatic bandwidth with a Bartlett
#'   kernel (default).
#'
#' If \code{t} and time aggregation are supplied, duplicate/near-duplicate time
#' points are aggregated before the ACF calculation.  The ACF itself is indexed
#' by observation lag, not physical time lag; therefore this correction assumes
#' observations are approximately equally spaced after any aggregation.  A
#' future continuous-time OU correction should use actual time differences
#' instead of this discrete-lag helper.
#'
#' For Spearman analyses, callers should supply the same scale used for the
#' correlation itself (ranks or rank residuals).
#'
#' @param x,y Numeric series of equal length.
#' @param t Optional numeric time vector.
#' @param aggregate_by_time One of \code{"none"}, \code{"mean"}, or
#'   \code{"median"}.
#' @param min_n Minimum number of observations and minimum clamped effective n.
#'   Must be at least 4 for downstream Fisher-z variance formulas.
#' @param clamp Whether to constrain effective n to \code{[min_n, n]}.
#' @param time_tol Optional non-negative time tolerance for aggregation.
#' @param effn_method Effective-n method.
#' @param L Optional Bartlett lag cutoff.
#' @param nw_bw Optional Newey-West bandwidth.
#' @return Numeric effective sample size, or \code{NA_real_} when too short.
#' @noRd
.eff_n <- function(
    x,
    y,
    t = NULL,
    aggregate_by_time = c("none", "mean", "median"),
    min_n = 4L,
    clamp = TRUE,
    time_tol = NULL,
    effn_method = c("nw", "bartlett", "ar1"),
    L = NULL,
    nw_bw = NULL) {

  aggregate_by_time <- match.arg(aggregate_by_time)
  effn_method <- match.arg(effn_method)

  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }

  if (!is.null(t)) {
    t <- suppressWarnings(as.numeric(t))
    if (length(t) != length(x)) {
      stop("t must be NULL or have the same length as x and y.", call. = FALSE)
    }
  }

  min_n <- suppressWarnings(as.integer(min_n))
  if (length(min_n) != 1L || is.na(min_n) || min_n < 4L) {
    stop("min_n must be a single integer >= 4.", call. = FALSE)
  }
  if (!is.null(time_tol)) {
    time_tol <- suppressWarnings(as.numeric(time_tol))
    if (length(time_tol) != 1L || !is.finite(time_tol) || time_tol < 0) {
      stop("time_tol must be NULL or a finite non-negative scalar.", call. = FALSE)
    }
  }

  ok <- is.finite(x) & is.finite(y)
  if (!is.null(t)) ok <- ok & is.finite(t)
  if (!any(ok)) return(NA_real_)

  x <- x[ok]
  y <- y[ok]
  if (!is.null(t)) t <- t[ok]

  if (!is.null(t) && aggregate_by_time != "none") {
    ag <- .aggregate_by_time_series(
      mat = cbind(x = x, y = y),
      tt = t,
      mode = aggregate_by_time,
      tol = time_tol
    )
    x <- as.numeric(ag$mat[, 1L])
    y <- as.numeric(ag$mat[, 2L])
    t <- ag$time
  } else if (!is.null(t) && anyDuplicated(t)) {
    warning(
      "Duplicate time stamps detected but aggregate_by_time='none'. ",
      "Consider mean/median aggregation or a positive time_tol.",
      call. = FALSE
    )
  }

  n <- min(length(x), length(y))
  if (n < min_n) return(NA_real_)

  # lag-k ACF vector, excluding lag 0 ------------------------------------------
  acf_vec <- function(v, lag.max) {
    lag.max <- min(as.integer(lag.max), length(v) - 1L)
    if (!is.finite(lag.max) || lag.max <= 0L) return(numeric(0))

    out <- try(
      stats::acf(
        v,
        plot = FALSE,
        lag.max = lag.max,
        na.action = stats::na.fail
      )$acf,
      silent = TRUE
    )
    if (inherits(out, "try-error") || !is.numeric(out)) {
      return(rep(0, lag.max))
    }

    a <- as.numeric(out[-1L])
    if (length(a) < lag.max) a <- c(a, rep(0, lag.max - length(a)))
    a <- a[seq_len(lag.max)]
    a[!is.finite(a)] <- 0
    pmax(pmin(a, 0.99), -0.99)
  }

  if (identical(effn_method, "ar1")) {
    lag1_acf <- function(v) {
      if (length(v) < 2L) return(0)
      out <- try(
        stats::acf(
          v,
          plot = FALSE,
          lag.max = 1L,
          na.action = stats::na.fail
        )$acf[2L],
        silent = TRUE
      )
      if (inherits(out, "try-error") || !is.finite(out)) 0 else as.numeric(out)
    }

    ax <- max(min(lag1_acf(x), 0.99), -0.99)
    ay <- max(min(lag1_acf(y), 0.99), -0.99)
    axy <- max(min(ax * ay, 0.95), -0.95)
    neff <- n * (1 - axy) / (1 + axy)
  } else {
    if (identical(effn_method, "nw")) {
      if (is.null(nw_bw)) {
        L_use <- max(1L, floor(4 * (n / 100)^(2 / 9)))
      } else {
        nw_bw_num <- suppressWarnings(as.numeric(nw_bw))
        if (length(nw_bw_num) != 1L || !is.finite(nw_bw_num) || nw_bw_num <= 0) {
          stop("nw_bw must be NULL or a finite positive scalar.", call. = FALSE)
        }
        L_use <- floor(nw_bw_num)
      }
    } else {
      if (is.null(L)) {
        L_use <- max(1L, floor(sqrt(n)))
      } else {
        L_num <- suppressWarnings(as.numeric(L))
        if (length(L_num) != 1L || !is.finite(L_num) || L_num <= 0) {
          stop("L must be NULL or a finite positive scalar.", call. = FALSE)
        }
        L_use <- floor(L_num)
      }
    }

    # ACF lags beyond n - 1 are undefined and must not be requested.
    L_use <- min(as.integer(L_use), n - 1L)
    if (L_use <= 0L) return(if (isTRUE(clamp)) as.numeric(n) else as.numeric(n))

    w <- 1 - seq_len(L_use) / (L_use + 1)
    ax <- acf_vec(x, L_use)
    ay <- acf_vec(y, L_use)
    gam <- sum(w * ax * ay, na.rm = TRUE)
    denom <- 1 + 2 * gam

    if (!is.finite(denom) || denom <= 1e-6) denom <- 1e-6
    neff <- n / denom
  }

  if (!is.finite(neff)) neff <- n
  if (isTRUE(clamp)) neff <- max(min_n, min(n, neff))
  as.numeric(neff)
}


#' Aggregate multivariate series at duplicate or near-duplicate times
#'
#' The same aggregation helper is used by correlation preprocessing and is
#' intentionally generic enough for future OU/residual code.  With a positive
#' tolerance, adjacent sorted observations separated by at most \code{tol} are
#' placed in the same contiguous time group.
#'
#' @param mat Matrix-like object with observations in rows.
#' @param tt Numeric time vector of length \code{nrow(mat)}.
#' @param mode Aggregation function, \code{"mean"} or \code{"median"}.
#' @param tol Optional non-negative grouping tolerance.
#' @return List with aggregated \code{mat} and representative mean \code{time}.
#' @noRd
.aggregate_by_time_series <- function(
    mat,
    tt,
    mode = c("mean", "median"),
    tol = NULL) {

  mode <- match.arg(mode)
  mat <- as.matrix(mat)

  if (is.null(tt)) return(list(mat = mat, time = NULL))

  tt <- suppressWarnings(as.numeric(tt))
  if (length(tt) != nrow(mat)) {
    stop("tt must have one value per matrix row.", call. = FALSE)
  }
  if (any(!is.finite(tt))) {
    stop("tt must contain only finite times.", call. = FALSE)
  }
  if (!is.null(tol)) {
    tol <- suppressWarnings(as.numeric(tol))
    if (length(tol) != 1L || !is.finite(tol) || tol < 0) {
      stop("tol must be NULL or a finite non-negative scalar.", call. = FALSE)
    }
  }

  if (nrow(mat) == 0L) return(list(mat = mat, time = tt))

  ord <- order(tt)
  tt2 <- tt[ord]
  m2 <- mat[ord, , drop = FALSE]

  if (!is.null(tol) && tol > 0) {
    grp <- c(1L, 1L + cumsum(diff(tt2) > tol))
  } else {
    grp <- as.integer(factor(tt2, levels = unique(tt2)))
  }

  idx_list <- split(seq_along(tt2), grp)
  agg_fun <- if (identical(mode, "median")) stats::median else base::mean

  mat_agg <- vapply(
    seq_len(ncol(m2)),
    function(j) {
      vapply(
        idx_list,
        function(ix) agg_fun(m2[ix, j], na.rm = TRUE),
        numeric(1)
      )
    },
    numeric(length(idx_list))
  )

  if (!is.matrix(mat_agg)) {
    mat_agg <- matrix(mat_agg, nrow = length(idx_list), ncol = ncol(m2))
  }
  colnames(mat_agg) <- colnames(m2)
  rownames(mat_agg) <- NULL

  tt_agg <- vapply(
    idx_list,
    function(ix) mean(tt2[ix], na.rm = TRUE),
    numeric(1)
  )
  names(tt_agg) <- NULL

  list(mat = mat_agg, time = tt_agg)
}

