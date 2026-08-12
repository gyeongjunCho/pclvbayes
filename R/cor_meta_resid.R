#' Pairwise correlations with subject-wise random-effects meta-analysis
#'
#' @description
#' Computes within-subject correlations for each unordered taxon pair and pools
#' them with random-effects meta-analysis using \pkg{metafor}.  CLR and ILR are
#' not used.
#'
#' With \code{transform = "alr"} (recommended), each pair is represented as the
#' three-part composition \code{{i, j, rest}} and transformed to pair-to-rest
#' ALRs.  The shared \code{.triplet_alr_transform()} helper is used so that
#' correlation/residual screening and future OU-residual analyses can use the
#' same zero replacement, rest floor, closure, and ALR cap logic.
#'
#' Before any zero replacement, each subject-level taxon pair must have enough
#' observed-positive longitudinal support.  By default each taxon must be
#' positive at >=50% of resolved time points and at >=4 time points. Thus a
#' trajectory with >50% zeros is excluded, while exactly 50% positive can pass
#' if the absolute-count criterion is also met. This is intentionally a
#' conservative "do not infer what was mostly unobserved" rule.
#'
#' Zero replacement is then used only as a finite-logratio safeguard for the
#' intermittent zeros that remain in an otherwise supported trajectory. It is
#' not interpreted as biological imputation. The default ALR contract currently
#' retains the pcLV-compatible rules: \code{zero_mode_alr = "minpos_time"},
#' \code{zero_minpos_base = "ij"}, \code{zero_minpos_alpha = 0.5},
#' \code{rest_floor_frac = 1}, and a fixed ALR cap of 12.
#'
#' This function does not perform the full-community spline smoothing used by
#' the pcLV dynamic model. Correlation keeps observed trajectories unsmoothed; a
#' future OU implementation should use actual time intervals. After pcLV
#' equivalence/benchmark checks, fit_pclv_bayes() is intended to reuse the same
#' observed-support and triplet helpers while retaining its weak spline as a
#' model-specific temporal stage.
#'
#' With \code{transform = "raw"}, correlations are computed directly on
#' relative abundance.  Optional partial correlation can control for
#' \code{rest = 1 - x_i - x_j}.  For Spearman partial correlation, the function
#' operates on ranks and computes effective sample size on rank residuals when
#' rest is controlled, including when \code{partial_method = "matrix"}.
#'
#' Pearson correlations use Fisher-z variance by default. Spearman correlations
#' use Bonett-Wright variance by default, with FHP and Fisher-z alternatives.
#' Optional Knapp-Hartung adjustment is supported.
#'
#' Effective sample-size correction uses discrete observation-lag ACFs.  Time
#' values are used for ordering and optional duplicate/near-duplicate time
#' aggregation, but the ACF is not a continuous-time model.  Therefore the
#' \code{nw}, \code{bartlett}, and \code{ar1} corrections assume approximately
#' equally spaced observations after aggregation.  A future OU implementation
#' should use actual time intervals while reusing the same triplet ALR helper.
#'
#' @param physeq A phyloseq object.
#' @param subject_col Sample metadata column identifying subjects.
#' @param time_col Sample metadata column identifying numeric sampling times.
#' @param taxa_vec Optional taxa identifiers to analyze.
#' @param interesting_taxa Optional taxa subset for pair screening.
#' @param method Correlation method, \code{"pearson"} or \code{"spearman"}.
#' @param transform Transformation, \code{"alr"} or \code{"raw"}.
#' @param min_n Minimum observations per subject after time aggregation.
#' @param min_k Minimum number of eligible subjects for meta-analysis.
#' @param meta_method Meta-analysis method, \code{"DL"} or \code{"REML"}.
#' @param use_knha Whether to use Knapp-Hartung adjustment.
#' @param prevalence_cut Optional prevalence filter in [0, 1].
#' @param mean_ra_cut Optional mean relative-abundance filter in [0, 1].
#' @param detrend Whether to remove a linear time trend within subject.
#' @param q_method Multiple-testing correction method.
#' @param q_weight_by Weighting scheme for weighted q-value correction.
#' @param k_filter Optional minimum subject count filter.
#' @param r_abs_min Optional minimum absolute pooled correlation in [0, 1].
#' @param var_method Variance estimator method.
#' @param return_subjectwise Whether to return subject-level results.
#' @param partial_rest Whether to adjust raw correlations for the rest component.
#' @param partial_method Spearman partial-correlation method.
#' @param zero_minpos_alpha Multiplier for minimum-positive zero replacement.
#' @param zero_minpos_base Base used for minimum-positive replacement.
#' @param zero_mode_alr Zero-replacement mode for ALR inputs.
#' @param zero_eps_fixed Positive fixed/fallback epsilon.
#' @param zero_lib_pseudo Positive library-size epsilon coefficient.
#' @param rest_floor_frac Non-negative rest-floor multiplier.
#' @param alr_cap_mode ALR capping mode. The canonical default is \code{"fixed"}.
#' @param alr_cap_value Positive ALR cap, or dynamic-cap ceiling.
#' @param effn_method Effective-sample-size correction method.
#' @param effn_L Optional Bartlett lag cutoff.
#' @param effn_bw Optional Newey-West bandwidth.
#' @param effn_aggregate_by_time Time aggregation for correlation and effective n.
#' @param effn_time_tol Optional non-negative time tolerance for aggregation.
#' @param effn_phi_cap Optional AR(1) product cap in [0, 1).
#' @param effn_blend Optional blend in [0, 1] toward ACF-corrected effective n.
#' @param effn_min_frac Optional effective-n floor as a fraction of input n.
#' @param cozero_check Whether to screen co-zero patterns.
#' @param cozero_tol Subject-level co-zero threshold in [0, 1].
#' @param cozero_action Action for high co-zero subjects.
#' @param cozero_min_subject_frac Pair-level flag threshold in [0, 1].
#' @param acf_correction Whether to apply autocorrelation correction.
#' @param nz_partner_min_frac Minimum observed-positive fraction required for
#'   each taxon within subject before zero replacement. Default 0.5.
#' @param nz_partner_min_n Minimum absolute number of observed-positive time
#'   points required for each taxon within subject. Default 4.
#' @param warn_alr Whether to emit one safeguard warning per affected pair.
#' @param all_columns Whether to retain setting columns irrelevant to the chosen
#'   effective-n and partial-correlation modes.
#' @return A tibble, or a list with \code{meta} and \code{subjectwise} when
#'   \code{return_subjectwise = TRUE}.
#' @export
cor_meta_resid <- function(
    physeq,
    subject_col,
    time_col,
    taxa_vec = NULL,
    interesting_taxa = NULL,
    method = c("pearson", "spearman"),
    transform = c("alr", "raw"),
    min_n = 5L,
    min_k = 2L,
    meta_method = c("DL", "REML"),
    use_knha = TRUE,
    prevalence_cut = NULL,
    mean_ra_cut = NULL,
    detrend = FALSE,
    q_method = "BH",
    q_weight_by = c("k", "n_median", "n_effin_median", "none"),
    k_filter = NULL,
    r_abs_min = NULL,
    var_method = c("auto", "fisherz", "bonett", "fhp"),
    return_subjectwise = FALSE,
    partial_rest = FALSE,
    partial_method = c("resid", "matrix"),
    zero_minpos_alpha = 0.5,
    zero_minpos_base = c("ij", "triplet"),
    zero_mode_alr = c("minpos_time", "minpos_subject", "lib", "fixed"),
    zero_eps_fixed = 1e-6,
    zero_lib_pseudo = 0.65,
    rest_floor_frac = 1.0,
    alr_cap_mode = c("fixed", "dynamic", "none"),
    alr_cap_value = 12,
    effn_method = c("nw", "bartlett", "ar1"),
    effn_L = NULL,
    effn_bw = NULL,
    effn_aggregate_by_time = c("median", "mean", "none"),
    effn_time_tol = NULL,
    effn_phi_cap = 0.6,
    effn_blend = 0.5,
    effn_min_frac = 0.6,
    cozero_check = FALSE,
    cozero_tol = 0.5,
    cozero_action = c("exclude", "flag", "warn"),
    cozero_min_subject_frac = 0.4,
    acf_correction = TRUE,
    nz_partner_min_frac = .PCLV_OBSERVED_SUPPORT_POLICY$nz_partner_min_frac,
    nz_partner_min_n = .PCLV_OBSERVED_SUPPORT_POLICY$nz_partner_min_n,
    warn_alr = FALSE,
    all_columns = FALSE) {

  stopifnot(requireNamespace("metafor", quietly = TRUE))
  stopifnot(requireNamespace("phyloseq", quietly = TRUE))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))

  method <- match.arg(method)
  transform <- match.arg(transform)
  meta_method <- match.arg(meta_method)
  var_method <- match.arg(var_method)
  zero_minpos_base <- match.arg(zero_minpos_base)
  zero_mode_alr <- match.arg(zero_mode_alr)
  alr_cap_mode <- match.arg(alr_cap_mode)
  effn_method <- match.arg(effn_method)
  effn_aggregate_by_time <- match.arg(effn_aggregate_by_time)
  partial_method <- match.arg(partial_method)
  q_method <- match.arg(q_method, c("BH", "BY", "qvalue", "wBH", "wqvalue"))
  q_weight_by <- match.arg(q_weight_by)
  cozero_action <- match.arg(cozero_action)

  if (identical(method, "pearson") &&
      !var_method %in% c("auto", "fisherz")) {
    warning(
      "Pearson correlations use Fisher-z variance; var_method='",
      var_method, "' will be treated as 'fisherz'.",
      call. = FALSE
    )
  }

  # --- Argument validation ----------------------------------------------------
  .scalar_num <- function(x, name, lower = -Inf, upper = Inf,
                          lower_open = FALSE, upper_open = FALSE,
                          allow_null = FALSE, integer = FALSE) {
    if (is.null(x)) {
      if (allow_null) return(NULL)
      stop(name, " must not be NULL.", call. = FALSE)
    }
    z <- suppressWarnings(as.numeric(x))
    if (length(z) != 1L || !is.finite(z)) {
      stop(name, " must be a finite scalar.", call. = FALSE)
    }
    lower_ok <- if (lower_open) z > lower else z >= lower
    upper_ok <- if (upper_open) z < upper else z <= upper
    if (!lower_ok || !upper_ok) {
      stop(name, " is outside its allowed range.", call. = FALSE)
    }
    if (integer && z != floor(z)) {
      stop(name, " must be an integer.", call. = FALSE)
    }
    if (integer) as.integer(z) else z
  }

  .scalar_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(name, " must be TRUE or FALSE.", call. = FALSE)
    }
    x
  }

  min_n <- .scalar_num(min_n, "min_n", lower = 4, integer = TRUE)
  min_k <- .scalar_num(min_k, "min_k", lower = 2, integer = TRUE)
  use_knha <- .scalar_flag(use_knha, "use_knha")
  detrend <- .scalar_flag(detrend, "detrend")
  return_subjectwise <- .scalar_flag(return_subjectwise, "return_subjectwise")
  partial_rest <- .scalar_flag(partial_rest, "partial_rest")
  cozero_check <- .scalar_flag(cozero_check, "cozero_check")
  acf_correction <- .scalar_flag(acf_correction, "acf_correction")
  warn_alr <- .scalar_flag(warn_alr, "warn_alr")
  all_columns <- .scalar_flag(all_columns, "all_columns")

  prevalence_cut <- .scalar_num(
    prevalence_cut, "prevalence_cut", 0, 1, allow_null = TRUE
  )
  mean_ra_cut <- .scalar_num(
    mean_ra_cut, "mean_ra_cut", 0, 1, allow_null = TRUE
  )
  k_filter <- .scalar_num(
    k_filter, "k_filter", lower = 1, allow_null = TRUE, integer = TRUE
  )
  r_abs_min <- .scalar_num(
    r_abs_min, "r_abs_min", 0, 1, allow_null = TRUE
  )

  zero_minpos_alpha <- .scalar_num(
    zero_minpos_alpha, "zero_minpos_alpha", lower = 0, lower_open = TRUE
  )
  zero_eps_fixed <- .scalar_num(
    zero_eps_fixed, "zero_eps_fixed", lower = 0, lower_open = TRUE
  )
  zero_lib_pseudo <- .scalar_num(
    zero_lib_pseudo, "zero_lib_pseudo", lower = 0, lower_open = TRUE
  )
  rest_floor_frac <- .scalar_num(
    rest_floor_frac, "rest_floor_frac", lower = 0
  )
  if (!identical(alr_cap_mode, "none")) {
    alr_cap_value <- .scalar_num(
      alr_cap_value, "alr_cap_value", lower = 0, lower_open = TRUE
    )
  }

  effn_L <- .scalar_num(
    effn_L, "effn_L", lower = 1, allow_null = TRUE, integer = TRUE
  )
  effn_bw <- .scalar_num(
    effn_bw, "effn_bw", lower = 1, allow_null = TRUE, integer = TRUE
  )
  effn_time_tol <- .scalar_num(
    effn_time_tol, "effn_time_tol", lower = 0, allow_null = TRUE
  )
  effn_phi_cap <- .scalar_num(
    effn_phi_cap, "effn_phi_cap", 0, 1,
    upper_open = TRUE, allow_null = TRUE
  )
  effn_blend <- .scalar_num(
    effn_blend, "effn_blend", 0, 1, allow_null = TRUE
  )
  effn_min_frac <- .scalar_num(
    effn_min_frac, "effn_min_frac", 0, 1, allow_null = TRUE
  )

  cozero_tol <- .scalar_num(cozero_tol, "cozero_tol", 0, 1)
  cozero_min_subject_frac <- .scalar_num(
    cozero_min_subject_frac, "cozero_min_subject_frac", 0, 1
  )
  nz_partner_min_frac <- .scalar_num(
    nz_partner_min_frac, "nz_partner_min_frac", 0, 1
  )
  nz_partner_min_n <- .scalar_num(
    nz_partner_min_n, "nz_partner_min_n", lower = 1, integer = TRUE
  )

  # --- 0) abundance matrix and metadata ---------------------------------------
  otu_mat <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  storage.mode(otu_mat) <- "double"

  if (is.null(rownames(otu_mat)) || is.null(colnames(otu_mat)) ||
      anyDuplicated(rownames(otu_mat)) || anyDuplicated(colnames(otu_mat))) {
    stop("otu_table must have unique taxon and sample names.", call. = FALSE)
  }
  if (any(!is.finite(otu_mat))) {
    stop("otu_table contains non-finite abundance values.", call. = FALSE)
  }
  if (any(otu_mat < 0)) {
    stop("otu_table contains negative abundance values.", call. = FALSE)
  }

  all_taxa <- rownames(otu_mat)
  if (is.null(taxa_vec)) taxa_vec <- all_taxa
  taxa_vec <- intersect(as.character(taxa_vec), all_taxa)
  if (length(taxa_vec) < 2L) {
    stop("taxa_vec must contain at least 2 taxa present in physeq.", call. = FALSE)
  }

  samp_df <- data.frame(
    sample = phyloseq::sample_names(physeq),
    as(phyloseq::sample_data(physeq), "data.frame"),
    check.names = FALSE
  )
  if (!all(c(subject_col, time_col) %in% colnames(samp_df))) {
    stop("subject_col/time_col not found in sample_data.", call. = FALSE)
  }

  meta_df <- samp_df |>
    dplyr::transmute(
      Sample = as.character(sample),
      subject = .data[[subject_col]],
      time = suppressWarnings(as.numeric(.data[[time_col]]))
    )

  if (anyNA(meta_df$Sample) || anyDuplicated(meta_df$Sample) ||
      anyNA(meta_df$subject) || any(!is.finite(meta_df$time))) {
    stop("Sample IDs, subjects, and numeric times must be complete and finite.", call. = FALSE)
  }
  if (!setequal(colnames(otu_mat), meta_df$Sample)) {
    stop("otu_table and sample_data sample IDs are not aligned.", call. = FALSE)
  }

  sample_idx <- match(meta_df$Sample, colnames(otu_mat))
  if (anyNA(sample_idx) || anyDuplicated(sample_idx)) {
    stop("Failed to align otu_table with sample_data.", call. = FALSE)
  }
  otu_mat <- otu_mat[, sample_idx, drop = FALSE]
  subjects <- unique(meta_df$subject)

  # Normalize over the complete community before any taxon filtering.
  lib <- colSums(otu_mat)
  if (any(!is.finite(lib)) || any(lib <= 0)) {
    stop("All samples must have positive finite library totals.", call. = FALSE)
  }
  ra_mat <- sweep(otu_mat, 2L, lib, "/")
  if (any(!is.finite(ra_mat)) || any(ra_mat < 0)) {
    stop("Relative-abundance normalization produced invalid values.", call. = FALSE)
  }

  # Optional RA filters ---------------------------------------------------------
  if (!is.null(mean_ra_cut) || !is.null(prevalence_cut)) {
    keep_mean <- if (is.null(mean_ra_cut)) {
      rep(TRUE, nrow(ra_mat))
    } else {
      rowMeans(ra_mat) >= mean_ra_cut
    }
    keep_prev <- if (is.null(prevalence_cut)) {
      rep(TRUE, nrow(ra_mat))
    } else {
      rowMeans(ra_mat > 0) >= prevalence_cut
    }
    taxa_vec <- intersect(taxa_vec, rownames(ra_mat)[keep_mean & keep_prev])
  }
  if (length(taxa_vec) < 2L) {
    stop("After filtering, fewer than 2 taxa remain.", call. = FALSE)
  }

  if (!is.null(interesting_taxa)) {
    interesting_taxa <- intersect(as.character(interesting_taxa), taxa_vec)
    if (!length(interesting_taxa)) {
      stop("interesting_taxa has no overlap with taxa_vec.", call. = FALSE)
    }
  }

  # --- 1) unordered pairs ------------------------------------------------------
  pair_mat <- utils::combn(taxa_vec, 2L)
  if (!is.null(interesting_taxa)) {
    keep_cols <- apply(
      pair_mat,
      2L,
      function(z) any(z %in% interesting_taxa)
    )
    pair_mat <- pair_mat[, keep_cols, drop = FALSE]
    if (!ncol(pair_mat)) {
      stop("No pairs satisfy interesting_taxa constraint.", call. = FALSE)
    }
  }

  npairs <- ncol(pair_mat)
  meta_rows <- vector("list", npairs)
  subjectwise_rows <- if (return_subjectwise) vector("list", npairs) else NULL

  use_partial <- isTRUE(partial_rest)
  if (identical(transform, "alr") && use_partial) {
    warning(
      "partial_rest is redundant for transform='alr' and will be disabled.",
      call. = FALSE
    )
    use_partial <- FALSE
  }

  .aggregate_series <- function(mat, tt) {
    if (identical(effn_aggregate_by_time, "none")) {
      return(list(mat = as.matrix(mat), time = tt))
    }
    .aggregate_by_time_series(
      mat = mat,
      tt = tt,
      mode = effn_aggregate_by_time,
      tol = effn_time_tol
    )
  }

  # Series supplied here are already on the exact time resolution used for the
  # correlation itself. Avoid a second aggregation inside .eff_n().
  .neff_and_nin <- function(xv, yv, tt) {
    list(
      neff = .eff_n(
        xv,
        yv,
        t = tt,
        aggregate_by_time = "none",
        min_n = 4L,
        time_tol = NULL,
        effn_method = effn_method,
        L = effn_L,
        nw_bw = effn_bw
      ),
      nin = length(xv)
    )
  }

  .linear_residual <- function(v, tt) {
    out <- try(stats::residuals(stats::lm(v ~ tt)), silent = TRUE)
    if (inherits(out, "try-error") || length(out) != length(v) ||
        any(!is.finite(out))) NULL else as.numeric(out)
  }

  .pair_cozero_summary <- function(coz_list) {
    fin <- coz_list[is.finite(coz_list)]
    list(
      median = if (length(fin)) stats::median(fin) else NA_real_,
      max = if (length(fin)) max(fin) else NA_real_,
      flag = if (isTRUE(cozero_check) && length(fin)) {
        mean(fin >= cozero_tol) >= cozero_min_subject_frac
      } else {
        FALSE
      }
    )
  }

  .cap_value_for_output <- function() {
    if (identical(transform, "alr") && !identical(alr_cap_mode, "none")) {
      as.numeric(alr_cap_value)
    } else {
      NA_real_
    }
  }

  # --- 2) pair and subject loops ----------------------------------------------
  #
  # SHARED-PREPROCESSING ROADMAP:
  # - cor_meta_resid(): observed support -> minimal zero safeguard -> triplet ALR.
  # - future OU: same two shared helpers, then a continuous-time model using dt.
  # - fit_pclv_bayes(): after current benchmarks/equivalence tests, evaluate
  #   observed support on pre-spline RA, keep its weak full-community spline as a
  #   pcLV-specific temporal step, then call the same triplet ALR backend.
  # Do not move spline smoothing into compositional_helpers.R.
  for (p in seq_len(npairs)) {
    i <- pair_mat[1L, p]
    j <- pair_mat[2L, p]

    r_list <- numeric(0)
    n_list <- integer(0)
    n_raw_list <- integer(0)
    n_effin_list <- integer(0)
    s_list <- character(0)
    coz_list <- numeric(0)
    floor_count_pair <- 0L
    cap_count_pair <- 0L
    support_screened_pair <- 0L
    support_excluded_pair <- 0L

    for (sb in subjects) {
      idx <- which(meta_df$subject == sb)
      if (length(idx) < min_n) next

      o <- order(meta_df$time[idx])
      idx <- idx[o]
      n_raw <- length(idx)

      Xi_raw <- as.numeric(otu_mat[i, idx, drop = TRUE])
      Xj_raw <- as.numeric(otu_mat[j, idx, drop = TRUE])
      lib_t <- as.numeric(lib[idx])

      Xi_ra <- as.numeric(ra_mat[i, idx, drop = TRUE])
      Xj_ra <- as.numeric(ra_mat[j, idx, drop = TRUE])
      rest_ra <- pmax(0, 1 - Xi_ra - Xj_ra)
      tt <- as.numeric(meta_df$time[idx])

      # Defensive complete-case guard. Upstream validation should make this a
      # no-op, but missing values are never recoded as biological zeros.
      cc <- stats::complete.cases(
        Xi_raw, Xj_raw, Xi_ra, Xj_ra, rest_ra, lib_t, tt
      )
      if (!any(cc)) next
      Xi_raw <- Xi_raw[cc]
      Xj_raw <- Xj_raw[cc]
      Xi_ra <- Xi_ra[cc]
      Xj_ra <- Xj_ra[cc]
      rest_ra <- rest_ra[cc]
      lib_t <- lib_t[cc]
      tt <- tt[cc]

      if (length(tt) < min_n) next

      # co-zero check -----------------------------------------------------------
      if (isTRUE(cozero_check)) {
        cozero_frac_sb <- mean((Xi_raw == 0) & (Xj_raw == 0))
        if (identical(cozero_action, "exclude") &&
            cozero_frac_sb >= cozero_tol) {
          next
        }
        if (identical(cozero_action, "warn") &&
            cozero_frac_sb >= cozero_tol) {
          warning(
            sprintf(
              "High co-zero in subject %s for pair (%s,%s): %.2f",
              as.character(sb), i, j, cozero_frac_sb
            ),
            call. = FALSE
          )
        }
      } else {
        cozero_frac_sb <- NA_real_
      }

      # Observed-support gate --------------------------------------------------
      #
      # IMPORTANT: support is evaluated BEFORE zero replacement.  When duplicate
      # or near-duplicate times are being aggregated for this analysis, use that
      # same resolved time grid here so technical replication cannot inflate the
      # apparent longitudinal support.  This gate is deliberately shared in
      # spirit with future OU/pcLV preprocessing: sparse trajectories are excluded
      # rather than reconstructed by spline or zero imputation.
      support_series <- .aggregate_series(
        cbind(xi = Xi_ra, xj = Xj_ra),
        tt
      )
      support_screened_pair <- support_screened_pair + 1L
      support <- .pair_observed_support(
        xi = support_series$mat[, 1L],
        xj = support_series$mat[, 2L],
        min_positive_frac = nz_partner_min_frac,
        min_positive_n = nz_partner_min_n
      )
      if (!isTRUE(support$keep)) {
        support_excluded_pair <- support_excluded_pair + 1L
        next
      }

      if (!is.finite(stats::var(Xi_ra)) || !is.finite(stats::var(Xj_ra)) ||
          stats::var(Xi_ra) < 1e-16 || stats::var(Xj_ra) < 1e-16) {
        next
      }

      if (identical(transform, "alr")) {
        trip <- .triplet_alr_transform(
          xi_raw = Xi_ra,
          xj_raw = Xj_ra,
          rest_raw = rest_ra,
          subject = rep(as.character(sb), length(Xi_ra)),
          lib = lib_t,
          zero_mode_alr = zero_mode_alr,
          minpos_alpha = zero_minpos_alpha,
          minpos_base = zero_minpos_base,
          eps_fixed = zero_eps_fixed,
          lib_eps_c = zero_lib_pseudo,
          rest_floor_frac = rest_floor_frac,
          alr_cap_mode = alr_cap_mode,
          alr_cap = alr_cap_value
        )

        floor_count_pair <- floor_count_pair + trip$floor_count
        cap_count_pair <- cap_count_pair + trip$cap_count

        out <- .aggregate_series(
          cbind(alr_i = trip$alr_i, alr_j = trip$alr_j),
          tt
        )
        alr_i <- as.numeric(out$mat[, 1L])
        alr_j <- as.numeric(out$mat[, 2L])
        tt <- as.numeric(out$time)

        if (isTRUE(detrend)) {
          if (length(alr_i) < 3L) next
          ri <- .linear_residual(alr_i, tt)
          rj <- .linear_residual(alr_j, tt)
          if (is.null(ri) || is.null(rj)) next
          alr_i <- ri
          alr_j <- rj
        }

        n_used <- length(alr_i)
        if (n_used < min_n || length(alr_j) != n_used) next

        if (identical(method, "spearman")) {
          Ri <- rank(alr_i, ties.method = "average")
          Rj <- rank(alr_j, ties.method = "average")
          r_s <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
          if (!is.finite(r_s)) next
          ne <- .neff_and_nin(Ri, Rj, tt)
        } else {
          r_s <- suppressWarnings(stats::cor(alr_i, alr_j, method = "pearson"))
          if (!is.finite(r_s)) next
          ne <- .neff_and_nin(alr_i, alr_j, tt)
        }

      } else {
        out <- .aggregate_series(
          cbind(xi = Xi_ra, xj = Xj_ra, xr = rest_ra),
          tt
        )
        xi <- as.numeric(out$mat[, 1L])
        xj <- as.numeric(out$mat[, 2L])
        xr <- as.numeric(out$mat[, 3L])
        tt <- as.numeric(out$time)

        if (isTRUE(detrend)) {
          if (length(xi) < 3L) next
          xi_r <- .linear_residual(xi, tt)
          xj_r <- .linear_residual(xj, tt)
          xr_r <- .linear_residual(xr, tt)
          if (is.null(xi_r) || is.null(xj_r) || is.null(xr_r)) next
          xi <- xi_r
          xj <- xj_r
          xr <- xr_r
        }

        n_used <- length(xi)
        if (n_used < min_n || length(xj) != n_used || length(xr) != n_used) next

        if (use_partial) {
          var_xr <- stats::var(xr)
          if (!is.finite(var_xr) || var_xr < 1e-12) {
            # Degenerate control: explicit fallback to simple correlation.
            if (identical(method, "spearman")) {
              Ri <- rank(xi, ties.method = "average")
              Rj <- rank(xj, ties.method = "average")
              r_s <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
              if (!is.finite(r_s)) next
              ne <- .neff_and_nin(Ri, Rj, tt)
            } else {
              r_s <- suppressWarnings(stats::cor(xi, xj, method = "pearson"))
              if (!is.finite(r_s)) next
              ne <- .neff_and_nin(xi, xj, tt)
            }
          } else if (identical(method, "spearman")) {
            Ri <- rank(xi, ties.method = "average")
            Rj <- rank(xj, ties.method = "average")
            Rr <- rank(xr, ties.method = "average")

            res_i <- .linear_residual(Ri, Rr)
            res_j <- .linear_residual(Rj, Rr)
            if (is.null(res_i) || is.null(res_j)) {
              r_s <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
              if (!is.finite(r_s)) next
              ne <- .neff_and_nin(Ri, Rj, tt)
            } else {
              if (identical(partial_method, "matrix")) {
                r_ij <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
                r_ir <- suppressWarnings(stats::cor(Ri, Rr, method = "pearson"))
                r_jr <- suppressWarnings(stats::cor(Rj, Rr, method = "pearson"))
                if (!all(is.finite(c(r_ij, r_ir, r_jr)))) next
                den <- sqrt(
                  max(1e-12, 1 - r_ir^2) *
                    max(1e-12, 1 - r_jr^2)
                )
                r_s <- (r_ij - r_ir * r_jr) / den
                r_s <- max(-1, min(1, r_s))
              } else {
                r_s <- suppressWarnings(stats::cor(res_i, res_j, method = "pearson"))
              }
              if (!is.finite(r_s)) next

              # Both partial methods use the same adjusted series for n_eff.
              ne <- .neff_and_nin(res_i, res_j, tt)
            }
          } else {
            res_i <- .linear_residual(xi, xr)
            res_j <- .linear_residual(xj, xr)
            if (is.null(res_i) || is.null(res_j)) {
              r_s <- suppressWarnings(stats::cor(xi, xj, method = "pearson"))
              if (!is.finite(r_s)) next
              ne <- .neff_and_nin(xi, xj, tt)
            } else {
              r_s <- suppressWarnings(stats::cor(res_i, res_j, method = "pearson"))
              if (!is.finite(r_s)) next
              ne <- .neff_and_nin(res_i, res_j, tt)
            }
          }
        } else {
          if (identical(method, "spearman")) {
            Ri <- rank(xi, ties.method = "average")
            Rj <- rank(xj, ties.method = "average")
            r_s <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
            if (!is.finite(r_s)) next
            ne <- .neff_and_nin(Ri, Rj, tt)
          } else {
            r_s <- suppressWarnings(stats::cor(xi, xj, method = "pearson"))
            if (!is.finite(r_s)) next
            ne <- .neff_and_nin(xi, xj, tt)
          }
        }
      }

      # Effective-n policy ------------------------------------------------------
      n_in <- if (is.finite(ne$nin)) as.numeric(ne$nin) else as.numeric(n_used)
      if (!is.finite(n_in) || n_in < 4) next

      n_eff_use <- if (isTRUE(acf_correction) && is.finite(ne$neff)) {
        as.numeric(ne$neff)
      } else {
        n_in
      }

      if (isTRUE(acf_correction)) {
        if (identical(effn_method, "ar1") && !is.null(effn_phi_cap)) {
          phi_hat <- (n_in - n_eff_use) / (n_in + n_eff_use)
          if (!is.finite(phi_hat)) phi_hat <- 0
          phi_hat <- max(0, min(0.999, phi_hat))
          phi_hat <- min(phi_hat, effn_phi_cap)
          n_eff_use <- n_in * (1 - phi_hat) / (1 + phi_hat)
        }

        if (!is.null(effn_blend)) {
          n_eff_use <- (1 - effn_blend) * n_in + effn_blend * n_eff_use
        }
        if (!is.null(effn_min_frac)) {
          n_eff_use <- max(n_eff_use, effn_min_frac * n_in)
        }
      }

      # Fisher-z variance requires n_eff > 3. Never let correction create more
      # information than the series actually contains.
      n_eff_use <- min(n_in, max(4, floor(n_eff_use)))

      r_list <- c(r_list, max(-1, min(1, as.numeric(r_s))))
      n_list <- c(n_list, as.integer(n_eff_use))
      n_raw_list <- c(n_raw_list, as.integer(n_raw))
      n_effin_list <- c(n_effin_list, as.integer(n_in))
      s_list <- c(s_list, as.character(sb))
      coz_list <- c(coz_list, cozero_frac_sb)
    }

    if (isTRUE(warn_alr) &&
        (floor_count_pair > 0L || cap_count_pair > 0L)) {
      warning(
        sprintf(
          "ALR safeguards applied for pair (%s, %s): rest-floor=%d, alr-cap=%d.",
          i, j, floor_count_pair, cap_count_pair
        ),
        call. = FALSE
      )
    }

    k_eff <- length(r_list)
    cozero_summary <- .pair_cozero_summary(coz_list)
    effective_var_method <- if (identical(method, "pearson")) {
      "fisherz"
    } else if (identical(var_method, "auto")) {
      "bonett"
    } else {
      var_method
    }

    .summary_value <- function(x, fun, default = NA_real_) {
      if (length(x)) fun(x) else default
    }

    build_na_row <- function() {
      tibble::tibble(
        i = i,
        j = j,
        method = method,
        transform = transform,
        partial_rest = use_partial,
        partial_method = if (use_partial && identical(method, "spearman")) {
          partial_method
        } else {
          NA_character_
        },
        k = k_eff,
        r_pooled = NA_real_,
        ciL = NA_real_,
        ciU = NA_real_,
        pval = NA_real_,
        tau2 = NA_real_,
        I2 = NA_real_,
        Q = NA_real_,
        Q_p = NA_real_,
        n_min = .summary_value(n_list, min, NA_integer_),
        n_median = .summary_value(n_list, stats::median),
        n_max = .summary_value(n_list, max, NA_integer_),
        n_raw_min = .summary_value(n_raw_list, min, NA_integer_),
        n_raw_median = .summary_value(n_raw_list, stats::median),
        n_raw_max = .summary_value(n_raw_list, max, NA_integer_),
        n_effin_min = .summary_value(n_effin_list, min, NA_integer_),
        n_effin_median = .summary_value(n_effin_list, stats::median),
        n_effin_max = .summary_value(n_effin_list, max, NA_integer_),
        cozero_frac_median = cozero_summary$median,
        cozero_frac_max = cozero_summary$max,
        cozero_flag_pair = cozero_summary$flag,
        support_screened_subjects = support_screened_pair,
        support_excluded_subjects = support_excluded_pair,
        nz_partner_min_frac = nz_partner_min_frac,
        nz_partner_min_n = nz_partner_min_n,
        meta_method = meta_method,
        knha = use_knha,
        var_method = effective_var_method,
        alr_cap_mode = if (identical(transform, "alr")) alr_cap_mode else NA_character_,
        alr_cap_value = .cap_value_for_output(),
        effn_method = effn_method,
        effn_L = if (identical(effn_method, "bartlett")) {
          if (is.null(effn_L)) NA_real_ else effn_L
        } else {
          NA_real_
        },
        effn_bw = if (identical(effn_method, "nw")) {
          if (is.null(effn_bw)) NA_real_ else effn_bw
        } else {
          NA_real_
        },
        effn_phi_cap = if (identical(effn_method, "ar1")) {
          if (is.null(effn_phi_cap)) NA_real_ else effn_phi_cap
        } else {
          NA_real_
        },
        effn_blend = if (is.null(effn_blend)) NA_real_ else effn_blend,
        effn_min_frac = if (is.null(effn_min_frac)) NA_real_ else effn_min_frac,
        effn_aggregate_by_time = effn_aggregate_by_time,
        effn_time_tol = if (is.null(effn_time_tol)) NA_real_ else effn_time_tol
      )
    }

    if (k_eff < min_k) {
      meta_rows[[p]] <- build_na_row()
      if (return_subjectwise) {
        subjectwise_rows[[p]] <- tibble::tibble(
          i = i,
          j = j,
          subject = s_list,
          r = r_list,
          n = n_list,
          n_raw = n_raw_list,
          n_eff_in = n_effin_list
        )
      }
      next
    }

    # --- 3) random-effects meta-analysis --------------------------------------
    # Exact +/-1 correlations are valid finite-sample outcomes but atanh(+/-1)
    # is infinite. Clamp only for the Fisher-z transformation; subjectwise r is
    # retained exactly in r_list.
    fisher_eps <- 1e-12
    r_for_meta <- pmax(pmin(r_list, 1 - fisher_eps), -1 + fisher_eps)
    yi <- atanh(r_for_meta)

    if (identical(method, "pearson")) {
      vi <- 1 / (n_list - 3)
    } else if (identical(effective_var_method, "bonett")) {
      vi <- (1 + r_for_meta^2 / 2) / (n_list - 3)
    } else if (identical(effective_var_method, "fhp")) {
      vi <- 1.06 / (n_list - 3)
    } else if (identical(effective_var_method, "fisherz")) {
      vi <- 1 / (n_list - 3)
    } else {
      stop("Unknown var_method for Spearman: ", effective_var_method, call. = FALSE)
    }

    if (any(!is.finite(yi)) || any(!is.finite(vi)) || any(vi <= 0)) {
      fit <- structure("invalid effect/variance", class = "try-error")
    } else {
      fit <- try(
        metafor::rma.uni(
          yi = yi,
          vi = vi,
          method = meta_method,
          test = if (use_knha) "knha" else "z"
        ),
        silent = TRUE
      )
    }

    if (inherits(fit, "try-error")) {
      meta_rows[[p]] <- build_na_row()
    } else {
      meta_rows[[p]] <- tibble::tibble(
        i = i,
        j = j,
        method = method,
        transform = transform,
        partial_rest = use_partial,
        partial_method = if (use_partial && identical(method, "spearman")) {
          partial_method
        } else {
          NA_character_
        },
        k = k_eff,
        r_pooled = tanh(as.numeric(fit$b[1L, 1L])),
        ciL = tanh(as.numeric(fit$ci.lb)),
        ciU = tanh(as.numeric(fit$ci.ub)),
        pval = as.numeric(fit$pval),
        tau2 = as.numeric(fit$tau2),
        I2 = as.numeric(fit$I2),
        Q = as.numeric(fit$QE),
        Q_p = as.numeric(fit$QEp),
        n_min = min(n_list),
        n_median = stats::median(n_list),
        n_max = max(n_list),
        n_raw_min = min(n_raw_list),
        n_raw_median = stats::median(n_raw_list),
        n_raw_max = max(n_raw_list),
        n_effin_min = min(n_effin_list),
        n_effin_median = stats::median(n_effin_list),
        n_effin_max = max(n_effin_list),
        cozero_frac_median = cozero_summary$median,
        cozero_frac_max = cozero_summary$max,
        cozero_flag_pair = cozero_summary$flag,
        support_screened_subjects = support_screened_pair,
        support_excluded_subjects = support_excluded_pair,
        nz_partner_min_frac = nz_partner_min_frac,
        nz_partner_min_n = nz_partner_min_n,
        meta_method = meta_method,
        knha = use_knha,
        var_method = effective_var_method,
        alr_cap_mode = if (identical(transform, "alr")) alr_cap_mode else NA_character_,
        alr_cap_value = .cap_value_for_output(),
        effn_method = effn_method,
        effn_L = if (identical(effn_method, "bartlett")) {
          if (is.null(effn_L)) NA_real_ else effn_L
        } else {
          NA_real_
        },
        effn_bw = if (identical(effn_method, "nw")) {
          if (is.null(effn_bw)) NA_real_ else effn_bw
        } else {
          NA_real_
        },
        effn_phi_cap = if (identical(effn_method, "ar1")) {
          if (is.null(effn_phi_cap)) NA_real_ else effn_phi_cap
        } else {
          NA_real_
        },
        effn_blend = if (is.null(effn_blend)) NA_real_ else effn_blend,
        effn_min_frac = if (is.null(effn_min_frac)) NA_real_ else effn_min_frac,
        effn_aggregate_by_time = effn_aggregate_by_time,
        effn_time_tol = if (is.null(effn_time_tol)) NA_real_ else effn_time_tol
      )
    }

    if (return_subjectwise) {
      subjectwise_rows[[p]] <- tibble::tibble(
        i = i,
        j = j,
        subject = s_list,
        r = r_list,
        n = n_list,
        n_raw = n_raw_list,
        n_eff_in = n_effin_list
      )
    }
  }

  meta_tbl <- dplyr::bind_rows(meta_rows)
  rownames(meta_tbl) <- NULL

  if (!all_columns && nrow(meta_tbl) > 0L) {
    uniq_meth <- unique(meta_tbl$effn_method)
    if (length(uniq_meth) == 1L) {
      if (identical(uniq_meth, "ar1")) {
        meta_tbl <- dplyr::select(meta_tbl, -effn_L, -effn_bw)
      } else if (identical(uniq_meth, "bartlett")) {
        meta_tbl <- dplyr::select(meta_tbl, -effn_bw, -effn_phi_cap)
      } else if (identical(uniq_meth, "nw")) {
        meta_tbl <- dplyr::select(meta_tbl, -effn_L, -effn_phi_cap)
      }
    }
    if (all(!meta_tbl$partial_rest)) {
      meta_tbl <- dplyr::select(meta_tbl, -partial_method)
    }
    if (all(meta_tbl$effn_aggregate_by_time == "none")) {
      meta_tbl <- dplyr::select(meta_tbl, -effn_time_tol)
    }
  }

  # --- 4) multiple-testing correction -----------------------------------------
  if (nrow(meta_tbl)) {
    .get_weights <- function(df, by) {
      if (identical(by, "none")) return(rep(1, nrow(df)))
      w <- switch(
        by,
        "n_median" = df$n_median,
        "n_effin_median" = df$n_effin_median,
        "k" = df$k,
        rep(1, nrow(df))
      )
      w[!is.finite(w) | w <= 0] <- 1
      w / mean(w)
    }

    .safe_qvalue <- function(p) {
      q <- rep(NA_real_, length(p))
      idx <- which(is.finite(p) & p >= 0 & p <= 1)
      if (!length(idx)) return(list(q = q, pi0 = NA_real_))

      fallback <- stats::p.adjust(p[idx], method = "BH")
      if (!requireNamespace("qvalue", quietly = TRUE)) {
        q[idx] <- fallback
        return(list(q = q, pi0 = NA_real_))
      }

      qq <- try(qvalue::qvalue(p[idx]), silent = TRUE)
      if (inherits(qq, "try-error")) {
        q[idx] <- fallback
        return(list(q = q, pi0 = NA_real_))
      }

      q[idx] <- qq$qvalues
      list(q = q, pi0 = as.numeric(qq$pi0))
    }

    meta_tbl$q_weight <- if (q_method %in% c("wBH", "wqvalue")) {
      q_weight_by
    } else {
      "none"
    }

    if (q_method %in% c("wBH", "wqvalue")) {
      w <- .get_weights(meta_tbl, q_weight_by)
      p_tilde <- meta_tbl$pval / w
      p_tilde[is.finite(p_tilde)] <- pmin(1, pmax(0, p_tilde[is.finite(p_tilde)]))

      if (identical(q_method, "wBH")) {
        meta_tbl$qval <- stats::p.adjust(p_tilde, method = "BH")
        meta_tbl$pi0 <- NA_real_
      } else {
        qq <- .safe_qvalue(p_tilde)
        meta_tbl$qval <- qq$q
        meta_tbl$pi0 <- qq$pi0
      }
    } else if (identical(q_method, "qvalue")) {
      qq <- .safe_qvalue(meta_tbl$pval)
      meta_tbl$qval <- qq$q
      meta_tbl$pi0 <- qq$pi0
    } else {
      meta_tbl$qval <- stats::p.adjust(meta_tbl$pval, method = q_method)
      meta_tbl$pi0 <- NA_real_
    }

    if (!is.null(k_filter)) {
      meta_tbl <- dplyr::filter(meta_tbl, k >= k_filter)
    }
    if (!is.null(r_abs_min)) {
      meta_tbl <- dplyr::filter(
        meta_tbl,
        is.finite(r_pooled),
        abs(r_pooled) >= r_abs_min
      )
    }
  }

  if (return_subjectwise) {
    list(
      meta = meta_tbl,
      subjectwise = dplyr::bind_rows(subjectwise_rows)
    )
  } else {
    meta_tbl
  }
}

