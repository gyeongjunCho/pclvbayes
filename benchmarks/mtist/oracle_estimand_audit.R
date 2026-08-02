oracle_sign <- function(x, tolerance = 1e-10) {
  if (!is.finite(x)) return(NA_integer_)
  if (abs(x) <= tolerance) 0L else as.integer(sign(x))
}

oracle_rest_set <- function(taxa, target, source)
  setdiff(taxa, c(target, source))

oracle_glv_growth <- function(x, growth, A)
  as.numeric(growth + A %*% x)

oracle_alr_derivative <- function(x, target, source, growth, A, taxa) {
  rest <- oracle_rest_set(taxa, target, source)
  names(x) <- taxa
  g <- stats::setNames(oracle_glv_growth(x, growth, A), taxa)
  weights <- x[rest] / sum(x[rest])
  as.numeric(g[[target]] - sum(weights * g[rest]))
}

oracle_contrast <- function(x, target, source, A, taxa) {
  rest <- oracle_rest_set(taxa, target, source)
  names(x) <- taxa
  weights <- x[rest] / sum(x[rest])
  as.numeric(A[target, source] - sum(weights * A[rest, source]))
}

oracle_projection <- function(y, xi, xj, subject, tolerance = 1e-10) {
  dat <- data.frame(y = y, xi = xi, xj = xj,
                    subject = factor(subject), stringsAsFactors = FALSE)
  if (any(!is.finite(dat$y)) || any(!is.finite(dat$xi)) ||
      any(!is.finite(dat$xj)))
    return(list(coefficient = NA_real_, sign = NA_integer_, rank = NA_integer_,
                columns = NA_integer_, condition = NA_real_,
                identifiable = FALSE, reason = "non_finite_design"))
  X <- if (nlevels(dat) >= 2L) stats::model.matrix(~ subject + xi + xj, dat) else stats::model.matrix(~ xi + xj, dat)
  qr_x <- qr(X, tol = tolerance)
  if (qr_x$rank < ncol(X))
    return(list(coefficient = NA_real_, sign = NA_integer_, rank = qr_x$rank,
                columns = ncol(X), condition = kappa(X),
                identifiable = FALSE, reason = "rank_deficient"))
  coefficient <- unname(qr.coef(qr_x, dat$y)[["xj"]])
  list(coefficient = coefficient, sign = oracle_sign(coefficient),
       rank = qr_x$rank, columns = ncol(X), condition = kappa(X),
       identifiable = TRUE, reason = NA_character_)
}

oracle_raw_design <- function(abundance, metadata, target, source, taxa,
                              derivative = NULL) {
  rows <- vector("list", length(unique(metadata$subject)))
  subjects <- unique(metadata$subject)
  for (s in seq_along(subjects)) {
    subject <- subjects[[s]]
    ix <- which(metadata$subject == subject)
    ix <- ix[order(metadata$time[ix])]
    x <- abundance[ix, taxa, drop = FALSE]
    relative <- x / rowSums(x)
    rest <- setdiff(taxa, c(target, source))
    zi <- log(relative[, target] / rowSums(relative[, rest, drop = FALSE]))
    zj <- log(relative[, source] / rowSums(relative[, rest, drop = FALSE]))
    current <- seq_len(nrow(x) - 1L)
    retained <- current[-1L]
    if (is.null(derivative)) {
      response <- diff(zi[current]) / diff(metadata$time[ix][current])
    } else {
      response <- derivative[ix[current[-length(current)]]]
    }
    rows[[s]] <- data.frame(
      subject = subject, time = metadata$time[ix][retained],
      y = response, xi_unscaled = zi[current[-length(current)]],
      xj_unscaled = zj[current[-length(current)]])
  }
  dat <- do.call(rbind, rows)
  dat$xi_scale <- stats::sd(dat$xi_unscaled)
  dat$xj_scale <- stats::sd(dat$xj_unscaled)
  dat$xi <- (dat$xi_unscaled - mean(dat$xi_unscaled)) / dat$xi_scale
  dat$xj <- (dat$xj_unscaled - mean(dat$xj_unscaled)) / dat$xj_scale
  rownames(dat) <- NULL
  dat
}

oracle_canonical_design <- function(relative, metadata, target, source, taxa,
                                    eps = 1e-6) {
  mat <- t(relative[, taxa, drop = FALSE])
  colnames(mat) <- metadata$Sample
  smoothed <- pclvbayes:::.precompute_spline_smoothed(
    mat, metadata[c("Sample", "subject", "time")], taxa, eps, 3L)
  if (inherits(smoothed, "pclv_failure")) return(smoothed)
  pclvbayes:::.make_pair_inputs_glv(
    sm_mat = smoothed, meta_df = metadata[c("Sample", "subject", "time")],
    j = source, i = target, min_pairs = 4L,
    zero_mode_alr = "minpos_time", minpos_alpha = 0.5,
    minpos_base = "ij", eps_fixed = 1e-6, lib_eps_c = 0.65,
    rest_floor_frac = 1, alr_cap = 12, smooth_scale = "logra",
    alr_spline_df = NULL, alr_spline_spar = NULL, alr_spline_cv = TRUE,
    nz_partner_min_frac = 0.15)
}

oracle_contrast_summary <- function(values, tolerance = 1e-10) {
  data.frame(
    contrast_mean = mean(values), contrast_median = stats::median(values),
    contrast_min = min(values), contrast_max = max(values),
    contrast_q05 = unname(stats::quantile(values, .05)),
    contrast_q95 = unname(stats::quantile(values, .95)),
    contrast_fraction_positive = mean(values > tolerance),
    contrast_fraction_negative = mean(values < -tolerance),
    contrast_fraction_near_zero = mean(abs(values) <= tolerance),
    contrast_constant_sign = length(unique(vapply(values, oracle_sign, integer(1),
                                                   tolerance = tolerance))) == 1L,
    contrast_changes_sign = any(values > tolerance) && any(values < -tolerance))
}

oracle_primary_category <- function(A_value, contrast_changes,
                                    instant_sign, finite_sign,
                                    noiseless_sign, observed_sign,
                                    posterior_sign) {
  A_sign <- oracle_sign(A_value)
  if (A_sign == 0L && any(c(finite_sign, noiseless_sign, observed_sign,
                             posterior_sign) != 0L, na.rm = TRUE))
    return("truth_zero_induced_compositional_effect")
  if (isTRUE(contrast_changes))
    return("state_dependent_no_single_sign")
  if (is.finite(instant_sign) && instant_sign == A_sign &&
      is.finite(finite_sign) && finite_sign != instant_sign)
    return("finite_interval_sign_flip")
  if (is.finite(finite_sign) && is.finite(noiseless_sign) &&
      noiseless_sign != finite_sign)
    return("preprocessing_sign_flip")
  if (is.finite(noiseless_sign) && is.finite(observed_sign) &&
      observed_sign != noiseless_sign)
    return("observation_noise_sign_flip")
  if (is.finite(observed_sign) && is.finite(posterior_sign) &&
      posterior_sign != observed_sign)
    return("bayesian_fit_sign_flip")
  if (is.finite(posterior_sign) &&
      posterior_sign %in% c(noiseless_sign, observed_sign) &&
      posterior_sign != A_sign)
    return("posterior_matches_oracle_estimand")
  "unresolved"
}

oracle_design_from_smoothed <- function(smoothed, metadata, target, source) {
  pclvbayes:::.make_pair_inputs_glv(
    sm_mat = smoothed, meta_df = metadata[c("Sample", "subject", "time")],
    j = source, i = target, min_pairs = 4L,
    zero_mode_alr = "minpos_time", minpos_alpha = 0.5,
    minpos_base = "ij", eps_fixed = 1e-6, lib_eps_c = 0.65,
    rest_floor_frac = 1, alr_cap = 12, smooth_scale = "logra",
    alr_spline_df = NULL, alr_spline_spar = NULL, alr_spline_cv = TRUE,
    nz_partner_min_frac = 0.15)
}

oracle_sign_text <- function(x) {
  value <- oracle_sign(x)
  if (is.na(value)) "NA" else as.character(value)
}
