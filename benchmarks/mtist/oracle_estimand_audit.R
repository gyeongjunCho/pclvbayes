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
  subject_value <- as.character(subject)
  dat <- data.frame(y = y, xi = xi, xj = xj,
                    subject = factor(subject_value), stringsAsFactors = FALSE)
  if (any(!is.finite(dat$y)) || any(!is.finite(dat$xi)) ||
      any(!is.finite(dat$xj)) || anyNA(subject_value) ||
      any(!nzchar(subject_value)))
    return(list(coefficient = NA_real_, sign = NA_integer_, rank = NA_integer_,
                columns = NA_integer_, condition = NA_real_,
                identifiable = FALSE, reason = "non_finite_design"))
  n_subjects <- length(unique(subject_value))
  X <- if (n_subjects >= 2L)
    stats::model.matrix(~ subject + xi + xj, dat)
  else
    stats::model.matrix(~ xi + xj, dat)
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
  if (!is.matrix(abundance) || nrow(abundance) != nrow(metadata) ||
      !all(c("subject", "time") %in% names(metadata)) ||
      !all(taxa %in% colnames(abundance)) ||
      !all(c(target, source) %in% taxa) ||
      !length(setdiff(taxa, c(target, source))))
    stop("Malformed oracle raw-design inputs.")
  if (!is.null(derivative) && length(derivative) != nrow(metadata))
    stop("Oracle derivative must align with observation rows.")
  subject_values <- as.character(metadata$subject)
  if (anyNA(subject_values) || any(!nzchar(subject_values)) ||
      any(!is.finite(metadata$time)))
    stop("Oracle subject and time values must be complete and finite.")

  rows <- vector("list", length(unique(subject_values)))
  subjects <- unique(subject_values)
  for (s in seq_along(subjects)) {
    subject <- subjects[[s]]
    ix <- which(subject_values == subject)
    ix <- ix[order(metadata$time[ix])]
    if (length(ix) < 2L) next
    predecessor <- seq_len(length(ix) - 1L)
    outcome <- predecessor + 1L
    time <- as.numeric(metadata$time[ix])
    dt <- time[outcome] - time[predecessor]
    if (any(!is.finite(dt)) || any(dt <= 0))
      stop("Oracle adjacent intervals require finite, strictly positive dt.")

    x <- abundance[ix, taxa, drop = FALSE]
    total <- rowSums(x)
    valid_observation <- apply(is.finite(x) & x >= 0, 1L, all) &
      is.finite(total) & total > 0
    relative <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x),
                       dimnames = dimnames(x))
    relative[valid_observation, ] <-
      x[valid_observation, , drop = FALSE] / total[valid_observation]
    rest <- setdiff(taxa, c(target, source))
    rest_total <- rowSums(relative[, rest, drop = FALSE])
    zi <- log(relative[, target] / rest_total)
    zj <- log(relative[, source] / rest_total)
    if (is.null(derivative)) {
      response <- (zi[outcome] - zi[predecessor]) / dt
    } else {
      response <- as.numeric(derivative[ix[outcome]])
    }
    keep <- valid_observation[predecessor] & valid_observation[outcome] &
      is.finite(zi[predecessor]) & is.finite(zi[outcome]) &
      is.finite(zj[predecessor]) & is.finite(response)
    rows[[s]] <- data.frame(
      subject = subject,
      predecessor_row = ix[predecessor][keep],
      outcome_row = ix[outcome][keep],
      predecessor_time = time[predecessor][keep],
      time = time[outcome][keep],
      dt = dt[keep],
      y = response[keep],
      xi_unscaled = zi[predecessor][keep],
      xj_unscaled = zj[predecessor][keep])
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  dat <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(dat)) stop("No valid adjacent oracle transitions.")
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


oracle_classify_sign_comparison <- function(comparison_sign, canonical_sign,
                                            tolerance = 1e-10,
                                            identifiable = TRUE,
                                            attempted = TRUE) {
  if (!attempted) return("unavailable")
  if (!is.finite(canonical_sign) || canonical_sign == 0L) return("unavailable")
  if (!isTRUE(identifiable) || !is.finite(comparison_sign)) return("ambiguous")
  value <- oracle_sign(comparison_sign, tolerance)
  if (is.na(value) || value == 0L) return("ambiguous")
  if (value == canonical_sign) "same" else "opposite"
}

oracle_axis_summary <- function(axis, classifications, notes = NA_character_) {
  classifications <- as.character(classifications)
  n <- table(factor(classifications,
                    levels = c("same", "opposite", "ambiguous", "unavailable")))
  valid <- unname(n[["same"]] + n[["opposite"]] + n[["ambiguous"]])
  risk <- if (valid > 0L)
    unname((n[["opposite"]] + 0.5 * n[["ambiguous"]]) / valid) else NA_real_
  data.frame(axis = axis, n_same = unname(n[["same"]]),
             n_opposite = unname(n[["opposite"]]),
             n_ambiguous = unname(n[["ambiguous"]]),
             n_unavailable = unname(n[["unavailable"]]),
             n_valid_comparisons = valid, r_g = risk,
             notes = notes, stringsAsFactors = FALSE)
}

oracle_empirical_susceptibility <- function(canonical_sign,
                                             canonical_sign_probability = NA_real_,
                                             canonical_lfsr = NA_real_,
                                             diagnostic_class = NA_character_,
                                             axis_classifications = list(),
                                             axis_notes = list(),
                                             tolerance = 1e-10) {
  valid_sign <- is.finite(canonical_sign) && canonical_sign %in% c(-1, 1)
  if (!valid_sign) return(list(
    axis = data.frame(), direction = data.frame(
      canonical_posterior_sign = NA_integer_,
      canonical_posterior_sign_probability = canonical_sign_probability,
      canonical_lfsr = canonical_lfsr, diagnostic_class = diagnostic_class,
      empirical_sign_reversal_susceptibility = NA_real_,
      worst_axis_sign_reversal_susceptibility = NA_real_, worst_axis = NA_character_,
      number_of_valid_axes = 0L, number_of_valid_comparisons = 0L,
      score_status = "missing", score_reason = "missing_canonical_posterior_sign",
      stringsAsFactors = FALSE)))
  axes <- names(axis_classifications)
  if (!length(axes)) return(list(
    axis = data.frame(), direction = data.frame(
      canonical_posterior_sign = canonical_sign,
      canonical_posterior_sign_probability = canonical_sign_probability,
      canonical_lfsr = canonical_lfsr, diagnostic_class = diagnostic_class,
      empirical_sign_reversal_susceptibility = NA_real_,
      worst_axis_sign_reversal_susceptibility = NA_real_, worst_axis = NA_character_,
      number_of_valid_axes = 0L, number_of_valid_comparisons = 0L,
      score_status = "missing", score_reason = "no_robustness_axes_available",
      stringsAsFactors = FALSE)))
  axis <- do.call(rbind, lapply(axes, function(a) {
    note <- axis_notes[[a]]
    if (is.null(note)) note <- NA_character_
    oracle_axis_summary(a, axis_classifications[[a]], note)
  }))
  valid <- is.finite(axis$r_g)
  if (!any(valid)) status <- "missing" else status <- "available"
  if (!any(valid)) score <- worst <- NA_real_ else {
    score <- mean(axis$r_g[valid])
    worst <- max(axis$r_g[valid])
  }
  worst_axes <- if (any(valid)) axis$axis[valid & axis$r_g == worst] else character()
  direction <- data.frame(
    canonical_posterior_sign = canonical_sign,
    canonical_posterior_sign_probability = canonical_sign_probability,
    canonical_lfsr = canonical_lfsr, diagnostic_class = diagnostic_class,
    empirical_sign_reversal_susceptibility = score,
    worst_axis_sign_reversal_susceptibility = worst,
    worst_axis = if (length(worst_axes)) paste(worst_axes, collapse = ";") else NA_character_,
    number_of_valid_axes = sum(valid),
    number_of_valid_comparisons = sum(axis$n_valid_comparisons[valid]),
    score_status = status,
    score_reason = if (status == "available") NA_character_ else "no_valid_axis",
    stringsAsFactors = FALSE)
  list(axis = axis, direction = direction)
}
