##### R/fit_pclv_bayes2.R
# Cross-kingdom extension: accept TWO phyloseq objects and fit ONLY cross pairs
# by mixing kingdom-level compositions with a user-tunable ratio (default 5:5).

#' @noRd
#' @keywords internal
.normalize_ratio <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 1L) x <- c(x, 1 - x)
  if (length(x) != 2L) stop("mix_ratio must be length 1 or 2.")
  if (any(!is.finite(x)) || any(x < 0)) stop("mix_ratio must be finite and non-negative.")
  s <- sum(x)
  if (!is.finite(s) || s <= 0) stop("mix_ratio sum must be > 0.")
  x / s
}

#' @noRd
#' @keywords internal
.align_two_phyloseq <- function(phy1, phy2, prefer_order = c("phy1","phy2")) {
  prefer_order <- match.arg(prefer_order)
  s1 <- phyloseq::sample_names(phy1)
  s2 <- phyloseq::sample_names(phy2)
  common <- intersect(s1, s2)
  if (length(common) < 3L) {
    stop("Need at least 3 common samples between the two phyloseq objects.")
  }
  # prune to common samples
  phy1 <- phyloseq::prune_samples(common, phy1)
  phy2 <- phyloseq::prune_samples(common, phy2)

  # keep the sample order of the preferred physeq
  ord <- if (prefer_order == "phy1") phyloseq::sample_names(phy1) else phyloseq::sample_names(phy2)

  phy1 <- phyloseq::prune_samples(ord, phy1)
  phy2 <- phyloseq::prune_samples(ord, phy2)
  list(phy1 = phy1, phy2 = phy2, samples = ord)
}

#' @noRd
#' @keywords internal
.pair_builder_mixed_triplet <- function(target, partner, ctx, eps, min_pairs) {
  # Build a 2-row "mini sm_mat" with already-scaled xi/xj, then reuse the
  # existing .build_pair_df_smoothed() logic (dt filtering, etc.).
  meta_df <- ctx$meta_df
  stopifnot(!is.null(meta_df), "Sample" %in% names(meta_df))

  km <- ctx$kingdom_map
  if (is.null(km)) stop("ctx$kingdom_map is required for mixed pair builder.")

  k_t <- unname(km[[target]])
  k_p <- unname(km[[partner]])
  if (!is.character(k_t) || !nzchar(k_t) || !is.character(k_p) || !nzchar(k_p)) {
    stop("Target/partner not found in ctx$kingdom_map.")
  }

  # pick the right smoothed matrix and weight for each taxon
  get_scaled <- function(taxon, k) {
    if (identical(k, ctx$kingdom_names[1])) {
      mat <- ctx$sm_mat_1
      w   <- ctx$mix_w[1]
    } else if (identical(k, ctx$kingdom_names[2])) {
      mat <- ctx$sm_mat_2
      w   <- ctx$mix_w[2]
    } else {
      stop("Unknown kingdom label in ctx$kingdom_map: ", k)
    }
    if (is.null(mat) || !is.matrix(mat)) stop("Smoothed matrix missing in ctx.")
    if (!taxon %in% rownames(mat)) stop("Taxon not in smoothed matrix: ", taxon)

    # meta_df sample order
    idx <- match(meta_df$Sample, colnames(mat))
    if (any(!is.finite(idx))) stop("Sample alignment failed when building pair_df.")
    as.numeric(w) * as.numeric(mat[taxon, idx])
  }

  xi <- get_scaled(target,  k_t)
  xj <- get_scaled(partner, k_p)

  mini <- rbind(xi, xj)
  rownames(mini) <- c(target, partner)
  colnames(mini) <- meta_df$Sample

  .build_pair_df_smoothed(mini, meta_df, j = partner, i = target, eps = eps, min_pairs = min_pairs)
}

#' Pairwise Bayesian pcLV for two phyloseq objects (cross-kingdom only)
#'
#' This is a cross-kingdom wrapper around \code{fit_pclv_bayes()} that:
#' \itemize{
#'   \item Pre-smooths each kingdom separately on \eqn{\log(\mathrm{RA}+\epsilon)}.
#'   \item Builds triplet-ALR using a mixed composition
#'         \eqn{p^{mix} = w_1 p^{(1)} \oplus w_2 p^{(2)}} where \eqn{w_1+w_2=1}.
#'   \item Fits ONLY cross pairs (one taxon from each phyloseq).
#' }
#'
#' Default mixing ratio is 5:5, but users can explore sensitivity (e.g., 3:7, 7:3).
#'
#' @param physeq_1,physeq_2 Two phyloseq objects (e.g., 16S and ITS) sharing sample IDs.
#' @param subject_col,time_col Columns in sample_data denoting subject and time.
#' @param taxa_vec_1,taxa_vec_2 Optional taxa subsets per kingdom.
#' @param mix_ratio Length-2 numeric ratio (default c(5,5)); normalized to weights.
#'        (e.g., c(3,7), c(7,3)). If length-1, treated as w1 and w2=1-w1.
#' @param kingdom_names Length-2 labels for the two kingdoms (default c("k1","k2")).
#' @param prefer_sample_order Which physeq decides sample ordering ("phy1" or "phy2").
#' @param ... Other arguments forwarded to the internal pcLV fitting machinery.
#' @return A tibble (same schema as \code{fit_pclv_bayes()}), with extra attributes:
#'   \code{mix_ratio}, \code{mix_weights}, \code{kingdom_names}, \code{common_samples}.
#' @return tibble/data.frame of fits
#' @export
fit_pclv_bayes2 <- function(physeq_1,
                            physeq_2,
                            subject_col,
                            time_col,
                            taxa_vec_1 = NULL,
                            taxa_vec_2 = NULL,
                            mix_ratio  = c(5, 5),
                            kingdom_names = c("k1","k2"),
                            prefer_sample_order = c("phy1","phy2"),
                            ...) {

  prefer_sample_order <- match.arg(prefer_sample_order)
  if (length(kingdom_names) != 2L) stop("kingdom_names must be length 2.")
  kingdom_names <- as.character(kingdom_names)

  # ---- align samples ----
  al <- .align_two_phyloseq(physeq_1, physeq_2, prefer_order = prefer_sample_order)
  physeq_1 <- al$phy1
  physeq_2 <- al$phy2

  # ---- weights (default 5:5) ----
  w <- .normalize_ratio(mix_ratio)
  # keep the user-facing ratio for printing (do not normalize for display)
  mix_ratio <- as.numeric(mix_ratio)
  if (length(mix_ratio) == 1L) mix_ratio <- c(mix_ratio, 1 - mix_ratio)

  # ---- meta + matrices ----
  meta_df  <- .get_sample_meta(physeq_1, subject_col, time_col)
  meta_df2 <- .get_sample_meta(physeq_2, subject_col, time_col)

  # reorder meta_df2 to meta_df order (defensive)
  if (!all(meta_df$Sample == meta_df2$Sample)) {
    meta_df2 <- meta_df2[match(meta_df$Sample, meta_df2$Sample), , drop = FALSE]
  }
  if (!isTRUE(all.equal(meta_df[, c("subject","time")],
                        meta_df2[, c("subject","time")],
                        check.attributes = FALSE))) {
    warning("subject/time columns differ between physeq_1 and physeq_2 after sample alignment; using physeq_1 metadata.")
  }

  mat_rel_1 <- .get_abund_matrix_precomputed(physeq_1)
  mat_rel_2 <- .get_abund_matrix_precomputed(physeq_2)

  if (is.null(taxa_vec_1)) taxa_vec_1 <- phyloseq::taxa_names(physeq_1)
  if (is.null(taxa_vec_2)) taxa_vec_2 <- phyloseq::taxa_names(physeq_2)
  taxa_vec_1 <- intersect(as.character(taxa_vec_1), rownames(mat_rel_1))
  taxa_vec_2 <- intersect(as.character(taxa_vec_2), rownames(mat_rel_2))
  if (length(taxa_vec_1) < 1L || length(taxa_vec_2) < 1L) {
    stop("Need at least 1 taxon in each phyloseq (taxa_vec_1 and taxa_vec_2).")
  }

  # ---- pull common '...' args used by smoothing (use defaults in fit_pclv_bayes.R if missing) ----
  dots <- list(...)
  eps <- if (!is.null(dots$eps)) dots$eps else 1e-8
  spline_df      <- if (!is.null(dots$spline_df)) dots$spline_df else NULL
  spline_spar    <- if (!is.null(dots$spline_spar)) dots$spline_spar else NULL
  spline_cv      <- if (!is.null(dots$spline_cv)) dots$spline_cv else TRUE
  min_unique_times <- if (!is.null(dots$min_unique_times)) dots$min_unique_times else 3L
  min_pairs      <- if (!is.null(dots$min_pairs)) dots$min_pairs else 12L

  # ---- legacy smoothing on each kingdom separately ----
  sm_mat_1 <- .precompute_spline_smoothed(
    mat_rel_1, meta_df, taxa_vec_1, eps,
    spline_df, spline_spar, spline_cv, min_unique_times
  )
  sm_mat_2 <- .precompute_spline_smoothed(
    mat_rel_2, meta_df, taxa_vec_2, eps,
    spline_df, spline_spar, spline_cv, min_unique_times
  )

  # ---- build ctx (borrow the same ctx keys as fit_pclv_bayes) ----
  # We reuse get_glv_model() to obtain the compiled exe, then point workers at it.
  # (The patched glv_helpers2 step-0 will call ctx$pair_builder if present.)
  mod <- get_glv_model(quiet = isTRUE(dots$quiet))

  # Taxa vector is the concatenation; we will fit ONLY cross pairs via pair_idx.
  taxa_vec <- c(taxa_vec_1, taxa_vec_2)
  kingdom_map <- setNames(
    c(rep(kingdom_names[1], length(taxa_vec_1)),
      rep(kingdom_names[2], length(taxa_vec_2))),
    taxa_vec
  )

  # Base ctx: copy most knobs from ...
  # (If a knob isn't provided, the downstream code in .run_one already has defaults.)
  ctx <- dots
  ctx$mod_exe_file <- mod$exe_file()
  ctx$meta_df      <- meta_df
  ctx$sm_mat       <- NULL  # unused when pair_builder is set
  ctx$eps          <- eps
  ctx$min_pairs    <- min_pairs
  ctx$min_unique_times <- min_unique_times

  # cross-kingdom specific
  ctx$sm_mat_1 <- sm_mat_1
  ctx$sm_mat_2 <- sm_mat_2
  ctx$kingdom_names <- kingdom_names
  ctx$kingdom_map   <- kingdom_map
  ctx$mix_ratio     <- mix_ratio
  ctx$mix_w         <- w

  # custom builder
  ctx$pair_builder  <- .pair_builder_mixed_triplet

  # ---- build cross-only pair list (i in 1..n1, j in (n1+1)..(n1+n2)) ----
  n1 <- length(taxa_vec_1); n2 <- length(taxa_vec_2)
  idx_i <- seq_len(n1)
  idx_j <- (n1 + 1L):(n1 + n2)
  pair_idx <- vector("list", length = n1 * n2)
  kk <- 0L
  for (ii in idx_i) {
    for (jj in idx_j) {
      kk <- kk + 1L
      pair_idx[[kk]] <- c(ii, jj)
    }
  }

  # ---- outer loop (same pattern as fit_pclv_bayes, but using pair_idx) ----
  n_workers_outer <- if (!is.null(dots$n_workers_outer)) as.integer(dots$n_workers_outer) else 1L
  progress <- if (!is.null(dots$progress)) dots$progress else "bar"
  progress_every <- if (!is.null(dots$progress_every)) as.integer(dots$progress_every) else 1L
  seed <- if (!is.null(dots$seed)) as.integer(dots$seed) else 1L
  n_tasks <- 2L * length(pair_idx)   # two directions per pair
  task_k  <- 0L
  out <- list(); k <- 0L

  has_progressr <- requireNamespace("progressr", quietly = TRUE)

  if (n_workers_outer <= 1L) {
    if (identical(progress, "bar")) {
      pb <- utils::txtProgressBar(min = 0, max = n_tasks, style = 3)
      on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }
    bump <- function() {
      task_k <<- task_k + 1L
      if (identical(progress, "bar") && task_k %% progress_every == 0L) {
        utils::setTxtProgressBar(pb, task_k)
      }
    }

    for (idx in pair_idx) {
      row <- .run_pair(
        idx[1], idx[2],
        taxa_vec = taxa_vec,
        .run_one = .run_one,
        ctx = ctx,
        progress = progress,
        mute_logs = FALSE,
        seed_base = seed
      )
      bump(); bump()
      if (!is.null(row)) {
        # attach kingdom labels + weights (for downstream filtering / reporting)
        row$kingdom_i <- kingdom_map[[row$i]]
        row$kingdom_j <- kingdom_map[[row$j]]
        row$mix_w_i   <- ifelse(row$kingdom_i == kingdom_names[1], w[1], w[2])
        row$mix_w_j   <- ifelse(row$kingdom_j == kingdom_names[1], w[1], w[2])
        k <- k + 1L
        out[[k]] <- row
      }
    }
  } else {
    # parallel mode: pair-level parallelization with future + furrr (same as fit_pclv_bayes)
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("furrr", quietly = TRUE)) {
      stop("Packages 'future' and 'furrr' are required for n_workers_outer > 1.")
    }
    op <- future::plan()
    on.exit(future::plan(op), add = TRUE)
    future::plan(future::multisession, workers = n_workers_outer)

    if (has_progressr && identical(progress, "bar")) {
      progressr::with_progress({
        p <- progressr::progressor(steps = length(pair_idx))
        out <- furrr::future_map(
          pair_idx,
          function(idx) {
            res <- .run_pair(
              idx[1], idx[2],
              taxa_vec = taxa_vec,
              .run_one = .run_one,
              ctx = ctx,
              progress = progress,
              mute_logs = TRUE,
              seed_base = seed
            )
            p(message = sprintf("cross pair %s-%s", taxa_vec[idx[1]], taxa_vec[idx[2]]))
            if (!is.null(res)) {
              res$kingdom_i <- kingdom_map[[res$i]]
              res$kingdom_j <- kingdom_map[[res$j]]
              res$mix_w_i   <- ifelse(res$kingdom_i == kingdom_names[1], w[1], w[2])
              res$mix_w_j   <- ifelse(res$kingdom_j == kingdom_names[1], w[1], w[2])
            }
            res
          },
          .options = furrr::furrr_options(seed = TRUE)
        )
      })
    } else {
      out <- furrr::future_map(
        pair_idx,
        function(idx) {
          res <- .run_pair(
            idx[1], idx[2],
            taxa_vec = taxa_vec,
            .run_one = .run_one,
            ctx = ctx,
            progress = progress,
            mute_logs = TRUE,
            seed_base = seed
          )
          if (!is.null(res)) {
            res$kingdom_i <- kingdom_map[[res$i]]
            res$kingdom_j <- kingdom_map[[res$j]]
            res$mix_w_i   <- ifelse(res$kingdom_i == kingdom_names[1], w[1], w[2])
            res$mix_w_j   <- ifelse(res$kingdom_j == kingdom_names[1], w[1], w[2])
          }
          res
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
    }
  }

  res <- dplyr::bind_rows(out)
  attr(res, "mix_ratio") <- mix_ratio
  attr(res, "mix_weights") <- setNames(as.numeric(w), kingdom_names)
  attr(res, "kingdom_names") <- kingdom_names
  attr(res, "common_samples") <- meta_df$Sample
  res
}
