# Shared longitudinal support + pair-to-rest triplet ALR preprocessing --------
#
# DESIGN PHILOSOPHY (2026-08)
# ---------------------------
# This file is intended to become the single deterministic preprocessing layer
# shared by cor_meta_resid(), future continuous-time OU analyses, and, after
# equivalence/benchmark checks, fit_pclv_bayes().
#
# The current philosophy is deliberately conservative:
#   1) Do not reconstruct a longitudinal trajectory that is mostly unobserved.
#      A taxon must have enough observed-positive time points before any zero
#      replacement is allowed.  Sparse subject x taxon trajectories are excluded
#      rather than imputed or rescued by smoothing.
#   2) Zero replacement is retained only as a finite-logratio safeguard for an
#      otherwise sufficiently observed trajectory.  It is NOT interpreted as a
#      biological estimate of the abundance at a zero observation.
#   3) The support filter is applied to taxa i and j, not to rest.  The rest is a
#      derived community remainder and is protected separately by the rest floor.
#   4) Temporal processing remains model-specific. cor_meta_resid() does not
#      spline the trajectory. A future OU model should model actual irregular dt
#      directly. pcLV may retain its weak full-community CV spline because the
#      derivative target benefits from mild smoothing, but the observed-support
#      decision should ultimately be made from pre-spline observed abundances.
#
# CURRENT SHARED SUPPORT POLICY CANDIDATE
# ---------------------------------------
# At the defaults below, exactly 50% positive observations passes; >50% zeros
# fails.  The absolute positive-count requirement prevents a very short series
# from passing on fraction alone.  These values are intentionally centralized so
# OU and pcLV can later adopt the same rule after validation.
#
# MIGRATION PLAN
# --------------
#   Phase 1 (this refactor): cor_meta_resid() uses .pair_observed_support() and
#     .triplet_alr_transform() from this file.
#   Phase 2: future OU code reuses the same support + triplet helpers, while the
#     OU likelihood itself uses actual time intervals. No spline-based zero
#     imputation is planned as the default.
#   Phase 3 (current): pcLV reuses the same low-level triplet ALR kernel while
#     retaining its own weak full-community spline and current 0.15 partner
#     support rule. The stricter observed-support candidate below is NOT applied
#     to pcLV by this refactor; changing support remains a separate scientific
#     policy decision.
#
# The minpos_time x 0.5 zero rule below is retained for compatibility for now.
# With the support filter in front of it, it should act only on intermittent zeros
# rather than determine a mostly-zero trajectory. Its exact epsilon policy can be
# revisited independently without changing the observed-support principle.
#
# Canonical triplet policy:
#   zero_mode_alr  = "minpos_time"
#   minpos_alpha   = 0.5
#   minpos_base    = "ij"
#   eps_fixed      = 1e-6
#   lib_eps_c      = 0.65
#   rest_floor_frac= 1.0
#   alr_cap_mode   = "fixed"
#   alr_cap        = 12
#
# Important: pcLV additionally smooths the full relative-abundance community
# within subject and recloses it before this triplet transform.  Smoothing is
# intentionally not part of this file; this file owns only {i,j,rest} handling.


.PCLV_OBSERVED_SUPPORT_POLICY <- list(
  nz_partner_min_frac = 0.5,
  nz_partner_min_n = 4L
)


#' Check whether a subject-level pair has enough observed-positive support
#'
#' This check must occur before zero replacement.  It deliberately treats zeros
#' as observed lack of positive support rather than values to be freely imputed.
#' The caller should supply the abundance vectors on the time grid used to define
#' an observed time point (for example, after duplicate-time aggregation).
#'
#' At the default fraction threshold, exactly half of the resolved observations
#' may be positive; a trajectory with more than half zeros is excluded.  Both
#' taxa must also satisfy the absolute positive-count threshold.
#'
#' @param xi,xj Non-negative observed abundance vectors of equal length.
#' @param min_positive_frac Minimum positive fraction required for each taxon.
#' @param min_positive_n Minimum number of positive observations per taxon.
#' @return A list with keep flag and support diagnostics for i and j.
#' @noRd
.pair_observed_support <- function(
    xi,
    xj,
    min_positive_frac = .PCLV_OBSERVED_SUPPORT_POLICY$nz_partner_min_frac,
    min_positive_n = .PCLV_OBSERVED_SUPPORT_POLICY$nz_partner_min_n) {

  xi <- suppressWarnings(as.numeric(xi))
  xj <- suppressWarnings(as.numeric(xj))
  if (!length(xi) || length(xi) != length(xj)) {
    stop("xi and xj must be non-empty vectors of equal length.", call. = FALSE)
  }
  if (any(!is.finite(xi)) || any(!is.finite(xj)) || any(xi < 0) || any(xj < 0)) {
    stop("Observed-support inputs must be finite and non-negative.", call. = FALSE)
  }

  min_positive_frac <- suppressWarnings(as.numeric(min_positive_frac))
  if (length(min_positive_frac) != 1L || !is.finite(min_positive_frac) ||
      min_positive_frac < 0 || min_positive_frac > 1) {
    stop("min_positive_frac must be a finite scalar in [0, 1].", call. = FALSE)
  }
  min_positive_n <- suppressWarnings(as.integer(min_positive_n))
  if (length(min_positive_n) != 1L || is.na(min_positive_n) || min_positive_n < 1L) {
    stop("min_positive_n must be a single integer >= 1.", call. = FALSE)
  }

  n <- length(xi)
  npos_i <- sum(xi > 0)
  npos_j <- sum(xj > 0)
  frac_i <- npos_i / n
  frac_j <- npos_j / n
  keep_i <- frac_i >= min_positive_frac && npos_i >= min_positive_n
  keep_j <- frac_j >= min_positive_frac && npos_j >= min_positive_n

  list(
    keep = isTRUE(keep_i && keep_j),
    keep_i = isTRUE(keep_i),
    keep_j = isTRUE(keep_j),
    n = as.integer(n),
    n_positive_i = as.integer(npos_i),
    n_positive_j = as.integer(npos_j),
    positive_frac_i = as.numeric(frac_i),
    positive_frac_j = as.numeric(frac_j),
    zero_frac_i = as.numeric(1 - frac_i),
    zero_frac_j = as.numeric(1 - frac_j),
    min_positive_frac = as.numeric(min_positive_frac),
    min_positive_n = as.integer(min_positive_n)
  )
}

.PCLV_TRIPLET_ALR_POLICY <- list(
  zero_mode_alr = "minpos_time",
  minpos_alpha = 0.5,
  minpos_base = "ij",
  eps_fixed = 1e-6,
  lib_eps_c = 0.65,
  rest_floor_frac = 1.0,
  alr_cap_mode = "fixed",
  alr_cap = 12
)


#' Build the canonical per-observation minpos-time epsilon
#'
#' Low-level vector/matrix helper shared by the high-level compositional
#' transform and pcLV's all-pair vectorized parent-process precompute.  It owns
#' only the deterministic epsilon rule; support filtering remains caller-specific.
#'
#' @param xi_raw,xj_raw Non-negative objects with identical shape.
#' @param rest_raw Optional non-negative object of the same shape; required when
#'   \code{minpos_base = "triplet"}.
#' @param minpos_alpha Positive multiplier for the local minimum positive value.
#' @param minpos_base Either \code{"ij"} or \code{"triplet"}.
#' @param eps_fixed Positive fallback when no positive base value exists.
#' @return Positive epsilon object with the same shape as \code{xi_raw}.
#' @noRd
.triplet_minpos_time_epsilon <- function(
    xi_raw,
    xj_raw,
    rest_raw = NULL,
    minpos_alpha = .PCLV_TRIPLET_ALR_POLICY$minpos_alpha,
    minpos_base = .PCLV_TRIPLET_ALR_POLICY$minpos_base,
    eps_fixed = .PCLV_TRIPLET_ALR_POLICY$eps_fixed) {

  minpos_base <- match.arg(minpos_base, c("ij", "triplet"))
  same_shape <- function(a, b) {
    length(a) == length(b) && identical(dim(a), dim(b))
  }
  if (!length(xi_raw) || !same_shape(xi_raw, xj_raw)) {
    stop("xi_raw and xj_raw must be non-empty objects with identical shape.",
         call. = FALSE)
  }
  if (any(!is.finite(xi_raw)) || any(!is.finite(xj_raw)) ||
      any(xi_raw < 0) || any(xj_raw < 0)) {
    stop("minpos-time inputs must be finite and non-negative.", call. = FALSE)
  }
  if (identical(minpos_base, "triplet")) {
    if (is.null(rest_raw) || !same_shape(xi_raw, rest_raw) ||
        any(!is.finite(rest_raw)) || any(rest_raw < 0)) {
      stop("rest_raw must be finite, non-negative, and shape-matched for triplet minpos.",
           call. = FALSE)
    }
  }

  minpos_alpha <- suppressWarnings(as.numeric(minpos_alpha))
  eps_fixed <- suppressWarnings(as.numeric(eps_fixed))
  if (length(minpos_alpha) != 1L || !is.finite(minpos_alpha) || minpos_alpha <= 0) {
    stop("minpos_alpha must be a finite positive scalar.", call. = FALSE)
  }
  if (length(eps_fixed) != 1L || !is.finite(eps_fixed) || eps_fixed <= 0) {
    stop("eps_fixed must be a finite positive scalar.", call. = FALSE)
  }

  xi_positive <- xi_raw
  xj_positive <- xj_raw
  xi_positive[xi_positive <= 0] <- Inf
  xj_positive[xj_positive <= 0] <- Inf
  min_positive <- pmin(xi_positive, xj_positive)

  if (identical(minpos_base, "triplet")) {
    rest_positive <- rest_raw
    rest_positive[rest_positive <= 0] <- Inf
    min_positive <- pmin(min_positive, rest_positive)
  }

  eps_t <- minpos_alpha * min_positive
  eps_t[!is.finite(eps_t) | eps_t <= 0] <- eps_fixed
  eps_t
}


#' Apply zero replacement, rest floor, closure, and ALR capping
#'
#' Low-level backend for the shared triplet transform.  Inputs may be vectors
#' or matrices as long as all objects have identical shape.  This separation is
#' deliberate: pairwise correlation/OU callers can use the high-level vector
#' wrapper below, while pcLV's all-pair vectorized precompute can later call this
#' backend without giving up matrix-level vectorization.
#'
#' @param xi_raw,xj_raw,rest_raw Non-negative objects with identical shape.
#' @param eps_t Positive replacement epsilon with identical shape.
#' @param rest_floor_frac Non-negative multiplier applied before closure.
#' @param alr_cap_mode One of \code{"fixed"}, \code{"dynamic"}, or \code{"none"}.
#' @param alr_cap Positive finite cap for fixed mode and ceiling for dynamic mode.
#' @return A list containing closed triplet components, \code{eps_star}, ALRs,
#'   and safeguard counts.
#' @noRd
.triplet_alr_apply <- function(
    xi_raw,
    xj_raw,
    rest_raw,
    eps_t,
    rest_floor_frac = .PCLV_TRIPLET_ALR_POLICY$rest_floor_frac,
    alr_cap_mode = .PCLV_TRIPLET_ALR_POLICY$alr_cap_mode,
    alr_cap = .PCLV_TRIPLET_ALR_POLICY$alr_cap) {

  alr_cap_mode <- match.arg(alr_cap_mode, c("fixed", "dynamic", "none"))

  same_shape <- function(a, b) {
    length(a) == length(b) && identical(dim(a), dim(b))
  }
  if (!same_shape(xi_raw, xj_raw) ||
      !same_shape(xi_raw, rest_raw) ||
      !same_shape(xi_raw, eps_t)) {
    stop("Triplet inputs and eps_t must have identical shape.", call. = FALSE)
  }
  if (!length(xi_raw)) {
    stop("Triplet inputs must contain at least one observation.", call. = FALSE)
  }
  if (any(!is.finite(xi_raw)) || any(!is.finite(xj_raw)) ||
      any(!is.finite(rest_raw)) || any(!is.finite(eps_t)) ||
      any(xi_raw < 0) || any(xj_raw < 0) || any(rest_raw < 0) ||
      any(eps_t <= 0)) {
    stop(
      "Triplet inputs must be finite/non-negative and eps_t must be positive.",
      call. = FALSE
    )
  }

  rest_floor_frac <- suppressWarnings(as.numeric(rest_floor_frac))
  if (length(rest_floor_frac) != 1L || !is.finite(rest_floor_frac) ||
      rest_floor_frac < 0) {
    stop("rest_floor_frac must be a finite non-negative scalar.", call. = FALSE)
  }
  if (!identical(alr_cap_mode, "none")) {
    alr_cap <- suppressWarnings(as.numeric(alr_cap))
    if (length(alr_cap) != 1L || !is.finite(alr_cap) || alr_cap <= 0) {
      stop("alr_cap must be a finite positive scalar.", call. = FALSE)
    }
  }

  # Exact Core order: replacement -> rest floor -> closure -> ALR -> cap.
  xi <- xi_raw
  xj <- xj_raw
  xr <- rest_raw
  xi[xi <= 0] <- eps_t[xi <= 0]
  xj[xj <= 0] <- eps_t[xj <= 0]
  xr[xr <= 0] <- eps_t[xr <= 0]

  rest_floor <- rest_floor_frac * eps_t
  floor_mask <- xr < rest_floor
  if (any(floor_mask)) xr[floor_mask] <- rest_floor[floor_mask]

  triplet_sum <- xi + xj + xr
  if (any(!is.finite(triplet_sum)) || any(triplet_sum <= 0)) {
    stop("Triplet closure failed because a row has a non-positive total.", call. = FALSE)
  }

  xi <- xi / triplet_sum
  xj <- xj / triplet_sum
  xr <- xr / triplet_sum
  eps_star <- eps_t / triplet_sum

  closure_error <- abs(xi + xj + xr - 1)
  closure_tolerance <- max(1e-12, 64 * .Machine$double.eps)
  if (any(!is.finite(closure_error)) ||
      max(closure_error) > closure_tolerance) {
    stop("Triplet closure failed numerical validation.", call. = FALSE)
  }

  alr_i_raw <- log(xi) - log(xr)
  alr_j_raw <- log(xj) - log(xr)

  if (identical(alr_cap_mode, "fixed")) {
    cap_vec <- eps_t
    cap_vec[] <- alr_cap
  } else if (identical(alr_cap_mode, "dynamic")) {
    es <- pmax(eps_star, .Machine$double.eps)
    theory_cap <- pmax(
      0,
      log(pmax(1 - 2 * es, .Machine$double.eps) / es)
    )
    cap_vec <- pmin(alr_cap, theory_cap)
  } else {
    cap_vec <- eps_t
    cap_vec[] <- Inf
  }

  cap_i <- is.finite(cap_vec) & abs(alr_i_raw) > cap_vec
  cap_j <- is.finite(cap_vec) & abs(alr_j_raw) > cap_vec

  alr_i <- alr_i_raw
  alr_j <- alr_j_raw
  finite_cap <- is.finite(cap_vec)
  if (any(finite_cap)) {
    alr_i[finite_cap] <- pmax(
      pmin(alr_i_raw[finite_cap], cap_vec[finite_cap]),
      -cap_vec[finite_cap]
    )
    alr_j[finite_cap] <- pmax(
      pmin(alr_j_raw[finite_cap], cap_vec[finite_cap]),
      -cap_vec[finite_cap]
    )
  }

  if (any(!is.finite(alr_i)) || any(!is.finite(alr_j))) {
    stop("Non-finite ALR values remained after triplet preprocessing.", call. = FALSE)
  }

  list(
    xi = xi,
    xj = xj,
    xr = xr,
    eps_star = eps_star,
    alr_i = alr_i,
    alr_j = alr_j,
    floor_count = as.integer(sum(floor_mask)),
    cap_count = as.integer(sum(cap_i) + sum(cap_j)),
    cap_count_i = as.integer(sum(cap_i)),
    cap_count_j = as.integer(sum(cap_j))
  )
}


#' Build zero-aware pair-to-rest triplet ALRs
#'
#' High-level vector wrapper around \code{.triplet_alr_apply()}.  It constructs
#' epsilon according to the selected policy and validates that the raw pair and
#' rest come from one closed relative-abundance composition.
#'
#' @param xi_raw,xj_raw Non-negative relative-abundance vectors for taxa i and j.
#' @param rest_raw Optional non-negative rest vector.  When \code{NULL}, it is
#'   computed as \code{pmax(0, 1 - xi_raw - xj_raw)}.
#' @param subject Optional subject vector. Required for
#'   \code{zero_mode_alr = "minpos_subject"}.
#' @param lib Optional positive library-size vector. Required for
#'   \code{zero_mode_alr = "lib"}.
#' @param zero_mode_alr Zero-replacement mode.
#' @param minpos_alpha Multiplier for minimum-positive replacement.
#' @param minpos_base Whether minimum-positive search uses \code{"ij"} or
#'   \code{"triplet"}.
#' @param eps_fixed Positive fallback/fixed epsilon.
#' @param lib_eps_c Positive coefficient for library-size replacement.
#' @param rest_floor_frac Non-negative multiplier applied before triplet closure.
#' @param alr_cap_mode One of \code{"fixed"}, \code{"dynamic"}, or \code{"none"}.
#' @param alr_cap Positive finite ALR cap for fixed mode, and ceiling for dynamic mode.
#' @return A list from \code{.triplet_alr_apply()} plus the epsilon vector and
#'   settings used.
#' @noRd
.triplet_alr_transform <- function(
    xi_raw,
    xj_raw,
    rest_raw = NULL,
    subject = NULL,
    lib = NULL,
    zero_mode_alr = c("minpos_time", "minpos_subject", "lib", "fixed"),
    minpos_alpha = .PCLV_TRIPLET_ALR_POLICY$minpos_alpha,
    minpos_base = c("ij", "triplet"),
    eps_fixed = .PCLV_TRIPLET_ALR_POLICY$eps_fixed,
    lib_eps_c = .PCLV_TRIPLET_ALR_POLICY$lib_eps_c,
    rest_floor_frac = .PCLV_TRIPLET_ALR_POLICY$rest_floor_frac,
    alr_cap_mode = c("fixed", "dynamic", "none"),
    alr_cap = .PCLV_TRIPLET_ALR_POLICY$alr_cap) {

  zero_mode_alr <- match.arg(zero_mode_alr)
  minpos_base <- match.arg(minpos_base)
  alr_cap_mode <- match.arg(alr_cap_mode)

  xi_raw <- suppressWarnings(as.numeric(xi_raw))
  xj_raw <- suppressWarnings(as.numeric(xj_raw))
  n <- length(xi_raw)

  if (length(xj_raw) != n) {
    stop("xi_raw and xj_raw must have the same length.", call. = FALSE)
  }
  if (!n) {
    stop("Triplet ALR input must contain at least one observation.", call. = FALSE)
  }
  if (any(!is.finite(xi_raw)) || any(!is.finite(xj_raw)) ||
      any(xi_raw < 0) || any(xj_raw < 0)) {
    stop("xi_raw and xj_raw must be finite and non-negative.", call. = FALSE)
  }

  closure_tol <- max(1e-10, 64 * .Machine$double.eps)
  pair_total <- xi_raw + xj_raw
  if (any(pair_total > 1 + closure_tol)) {
    stop(
      "xi_raw + xj_raw exceeds 1; pair inputs must come from a common closed composition.",
      call. = FALSE
    )
  }

  if (is.null(rest_raw)) {
    rest_raw <- pmax(0, 1 - pair_total)
  } else {
    rest_raw <- suppressWarnings(as.numeric(rest_raw))
    if (length(rest_raw) != n || any(!is.finite(rest_raw)) || any(rest_raw < 0)) {
      stop("rest_raw must be finite, non-negative, and match input length.", call. = FALSE)
    }
    raw_closure_error <- abs(pair_total + rest_raw - 1)
    if (any(!is.finite(raw_closure_error)) ||
        max(raw_closure_error) > closure_tol) {
      stop(
        "xi_raw, xj_raw, and rest_raw must form a closed relative-abundance triplet.",
        call. = FALSE
      )
    }
  }

  positive_scalar <- function(x, name, allow_zero = FALSE) {
    z <- suppressWarnings(as.numeric(x))
    ok <- length(z) == 1L && is.finite(z) &&
      (if (allow_zero) z >= 0 else z > 0)
    if (!ok) {
      stop(
        name,
        if (allow_zero) " must be a finite non-negative scalar." else
          " must be a finite positive scalar.",
        call. = FALSE
      )
    }
    z
  }

  minpos_alpha <- positive_scalar(minpos_alpha, "minpos_alpha")
  eps_fixed <- positive_scalar(eps_fixed, "eps_fixed")
  lib_eps_c <- positive_scalar(lib_eps_c, "lib_eps_c")
  rest_floor_frac <- positive_scalar(
    rest_floor_frac, "rest_floor_frac", allow_zero = TRUE
  )
  if (!identical(alr_cap_mode, "none")) {
    alr_cap <- positive_scalar(alr_cap, "alr_cap")
  }

  if (!is.null(subject) && length(subject) != n) {
    stop("subject must be NULL or match input length.", call. = FALSE)
  }
  if (identical(zero_mode_alr, "minpos_subject") &&
      (is.null(subject) || anyNA(subject))) {
    stop(
      "subject is required and must not contain NA for minpos_subject.",
      call. = FALSE
    )
  }

  if (identical(zero_mode_alr, "lib")) {
    lib <- suppressWarnings(as.numeric(lib))
    if (length(lib) != n || any(!is.finite(lib)) || any(lib <= 0)) {
      stop(
        "lib must contain positive finite library sizes for zero_mode_alr='lib'.",
        call. = FALSE
      )
    }
  }

  # Build epsilon_t using the policy selected by the caller.
  if (identical(zero_mode_alr, "minpos_time")) {
    eps_t <- .triplet_minpos_time_epsilon(
      xi_raw = xi_raw,
      xj_raw = xj_raw,
      rest_raw = rest_raw,
      minpos_alpha = minpos_alpha,
      minpos_base = minpos_base,
      eps_fixed = eps_fixed
    )
  } else if (identical(zero_mode_alr, "minpos_subject")) {
    subject_chr <- as.character(subject)
    by_subject <- split(seq_len(n), subject_chr)
    subject_min <- vapply(by_subject, function(ix) {
      positive <- c(
        xi_raw[ix][xi_raw[ix] > 0],
        xj_raw[ix][xj_raw[ix] > 0]
      )
      if (identical(minpos_base, "triplet")) {
        positive <- c(positive, rest_raw[ix][rest_raw[ix] > 0])
      }
      if (length(positive)) min(positive) else NA_real_
    }, numeric(1))
    min_pos <- unname(subject_min[subject_chr])
    eps_t <- ifelse(
      is.finite(min_pos) & min_pos > 0,
      minpos_alpha * min_pos,
      eps_fixed
    )
  } else if (identical(zero_mode_alr, "lib")) {
    eps_t <- pmax(eps_fixed, lib_eps_c / lib)
  } else {
    eps_t <- rep(eps_fixed, n)
  }

  eps_t[!is.finite(eps_t) | eps_t <= 0] <- eps_fixed

  out <- .triplet_alr_apply(
    xi_raw = xi_raw,
    xj_raw = xj_raw,
    rest_raw = rest_raw,
    eps_t = eps_t,
    rest_floor_frac = rest_floor_frac,
    alr_cap_mode = alr_cap_mode,
    alr_cap = alr_cap
  )

  out$eps_t <- eps_t
  out$settings <- list(
    zero_mode_alr = zero_mode_alr,
    minpos_alpha = minpos_alpha,
    minpos_base = minpos_base,
    eps_fixed = eps_fixed,
    lib_eps_c = lib_eps_c,
    rest_floor_frac = rest_floor_frac,
    alr_cap_mode = alr_cap_mode,
    alr_cap = if (identical(alr_cap_mode, "none")) NA_real_ else alr_cap
  )
  out
}
