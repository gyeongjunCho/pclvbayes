#' Pairwise correlations (pair-to-rest ALR or raw RA) with subject-wise random-effects meta-analysis
#'
#' @description
#' Computes within-subject correlations for each unordered taxon pair and pools them
#' via random-effects meta-analysis (REML or DL) using \pkg{metafor}. **CLR and ILR are not used.**
#'
#' Two transformation modes are supported:
#'
#' - `transform = "alr"` (recommended): To mitigate spurious near-1 correlations from
#'   mathematical coupling, use **pair-to-rest ALR** logratios after zero-aware
#'   multiplicative replacement on the 3-part ``composition `{i, j, rest}```:
#'     - `alr_i = log(x_i / x_rest)`, `alr_j = log(x_j / x_rest)`.
#'   We then compute within-subject correlations between `alr_i` and `alr_j`
#'   (Spearman or Pearson) and aggregate across subjects.
#'
#'   **ALR safeguards (important):**
#'   When `rest` is extremely small, ALR values can explode. This function provides
#'   (i) time-varying zero replacement (`zero_mode_alr`), (ii) a **rest floor** (`rest_floor_frac`),
#'   and (iii) optional **winsorization** controlled by `alr_cap_mode`.
#'
#' - `transform = "raw"`: Work directly on relative abundance (RA) time series.
#'   - `partial_rest = FALSE`: simple correlation (Spearman or Pearson).
#'   - `partial_rest = TRUE`: partial correlation controlling for `rest = 1 - x_i - x_j`.
#'       * Spearman: rank-transform `x_i, x_j, rest` to `R_i, R_j, R_r`,
#'         then either regress residuals (`partial_method = "resid"`, default) or
#'         use the rank-correlation matrix formula (`partial_method = "matrix"`) to obtain the partial r.
#'       * Pearson: regress `x_i ~ rest`, `x_j ~ rest`, then correlate residuals.
#'
#' Variance scaling for meta-analysis defaults to **Fisher's z** for Pearson and
#' **Bonett–Wright** for Spearman; **Fieller–Hartley–Pearson (FHP)** is available.
#' Optional Knapp–Hartung adjustment is supported.
#'
#' **Effective n (neff):**
#' - `effn_method = c("nw","bartlett","ar1")` (default `"nw"`).
#' - For short series, you can weaken ACF correction via:
#'   * `effn_phi_cap` (AR(1) only): cap `|phi_x * phi_y|`.
#'   * `effn_blend`: blend `n_eff <- (1 - b) * n + b * n_eff` (e.g., `b = 0.5`).
#'   * `effn_min_frac`: floor `n_eff >= f * n` (e.g., `f = 0.6`).
#' - Optional time aggregation for `n_eff`: `effn_aggregate_by_time` (`"none"|"mean"|"median"`)
#'   with tolerance `effn_time_tol`.
#' - For Spearman, `n_eff` is computed on the **same scale used for correlation**:
#'   ranks (simple Spearman) or residuals of ranks (partial-resid Spearman).
#'
#' @return
#' If `return_subjectwise = TRUE`, a list with `meta` and `subjectwise`.
#' Otherwise, a data.frame; the meta table includes `piL/piU` (prediction interval) and
#' the `n_eff` settings used.
#'
#' @export
cor_meta_resid <- function(physeq,
                           subject_col,
                           time_col,
                           taxa_vec = NULL,
                           interesting_taxa = NULL,
                           method = c("pearson", "spearman"),
                           transform = c("alr", "raw"),
                           min_n = 5,
                           min_k = 2,
                           meta_method = c("DL", "REML"),
                           use_knha = TRUE,
                           prevalence_cut = NULL,
                           mean_ra_cut = NULL,
                           detrend = FALSE,
                           q_method = "BH",
                           q_weight_by = c("k","n_median","n_effin_median","none"),
                           k_filter = NULL,
                           r_abs_min = NULL,
                           var_method = c("auto", "fisherz", "bonett", "fhp"),
                           return_subjectwise = FALSE,
                           partial_rest = FALSE,
                           partial_method = c("resid","matrix"),
                           zero_minpos_alpha = 0.5,
                           zero_minpos_base = c("ij", "triplet"),
                           zero_mode_alr = c("minpos_subject", "minpos_time", "lib"),
                           zero_lib_pseudo = 0.5,
                           rest_floor_frac = 1.0,
                           # ALR capping
                           alr_cap_mode = c("dynamic","fixed","none"),
                           alr_cap_value = 12,
                           # neff options
                           effn_method = c("nw","bartlett","ar1"),
                           effn_L = NULL,
                           effn_bw = NULL,
                           effn_aggregate_by_time = c("median","mean","none"),
                           effn_time_tol = NULL,
                           # weak ACF-correction controls (optional)
                           effn_phi_cap = 0.6,   # e.g., 0.6 (cap on |phi_x * phi_y|)
                           effn_blend   = 0.5,   # e.g., 0.5 (blend toward raw n)
                           effn_min_frac = 0.6,  # e.g., 0.6 (floor for n_eff relative to n)
                           # misc
                           acf_correction = TRUE,
                           nz_partner_min_frac = 0.15,
                           warn_alr = FALSE) {
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
  q_method <- match.arg(q_method, c("BH","BY","qvalue","wBH","wqvalue"))
  q_weight_by <- match.arg(q_weight_by)

  # --- 0) matrix & metadata ---------------------------------------------------
  otu_mat <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq))
    otu_mat <- t(otu_mat)
  all_taxa <- rownames(otu_mat)
  if (is.null(taxa_vec)) taxa_vec <- all_taxa
  taxa_vec <- intersect(taxa_vec, all_taxa)
  if (length(taxa_vec) < 2) stop("taxa_vec must contain at least 2 taxa.")

  samp_df <- data.frame(
    sample = phyloseq::sample_names(physeq),
    as(phyloseq::sample_data(physeq), "data.frame"),
    check.names = FALSE
  )
  if (!all(c(subject_col, time_col) %in% colnames(samp_df)))
    stop("subject_col/time_col not found in sample_data.")

  meta_df <- samp_df |>
    dplyr::transmute(Sample = sample,
                     subject = .data[[subject_col]],
                     time = suppressWarnings(as.numeric(.data[[time_col]])))

  common_samples <- intersect(colnames(otu_mat), meta_df$Sample)
  otu_mat <- otu_mat[, common_samples, drop = FALSE]
  meta_df <- meta_df[match(common_samples, meta_df$Sample), , drop = FALSE]
  subjects <- unique(meta_df$subject)

  # Normalize to relative abundance
  lib <- colSums(otu_mat, na.rm = TRUE)
  lib[!is.finite(lib) | lib <= 0] <- 1
  ra_mat <- sweep(otu_mat, 2, lib, "/")
  ra_mat[!is.finite(ra_mat)] <- 0

  # Optional RA filters
  if (!is.null(mean_ra_cut) || !is.null(prevalence_cut)) {
    keep_mean <- if (is.null(mean_ra_cut)) rep(TRUE, nrow(ra_mat)) else rowMeans(ra_mat, na.rm = TRUE) >= mean_ra_cut
    keep_prev <- if (is.null(prevalence_cut)) rep(TRUE, nrow(ra_mat)) else rowMeans(ra_mat > 0, na.rm = TRUE) >= prevalence_cut
    taxa_vec <- intersect(taxa_vec, rownames(ra_mat)[keep_mean & keep_prev])
  }
  if (length(taxa_vec) < 2) stop("After filtering, fewer than 2 taxa remain.")

  # Restrict to interesting_taxa if provided
  if (!is.null(interesting_taxa)) {
    interesting_taxa <- intersect(interesting_taxa, taxa_vec)
    if (length(interesting_taxa) == 0)
      stop("interesting_taxa has no overlap with taxa_vec.")
  }

  # --- 1) generate all unordered pairs ---------------------------------------
  pair_mat <- utils::combn(taxa_vec, 2)
  if (!is.null(interesting_taxa)) {
    keep_cols <- apply(pair_mat, 2, function(col)
      (col[1] %in% interesting_taxa) || (col[2] %in% interesting_taxa))
    pair_mat <- pair_mat[, keep_cols, drop = FALSE]
    if (ncol(pair_mat) == 0) stop("No pairs satisfy interesting_taxa constraint.")
  }

  npairs <- ncol(pair_mat)
  meta_rows <- vector("list", npairs)
  subjectwise_rows <- if (return_subjectwise) vector("list", npairs) else NULL

  # --- 1.5) partial_rest 플래그 (lag 제거됨) ----------------------------------
  use_partial <- isTRUE(partial_rest)
  if (transform == "alr" && isTRUE(partial_rest)) {
    warning("partial_rest is redundant in transform='alr' and will be disabled.", call. = FALSE)
    use_partial <- FALSE
  }

  # --- 내부 헬퍼: neff 계산 + 입력길이(집계후) 추정 ----------------------------
  .neff_and_nin <- function(xv, yv, tt) {
    nin <- length(xv)
    if (!is.null(tt) && effn_aggregate_by_time != "none" && length(tt) == length(xv)) {
      ord <- order(tt); tt2 <- tt[ord]
      if (!is.null(effn_time_tol) && is.finite(effn_time_tol) && effn_time_tol > 0) {
        grp <- c(1L, 1L + cumsum(diff(tt2) > effn_time_tol))
        nin <- length(unique(grp))
      } else {
        nin <- length(unique(tt2))
      }
    }
    neff <- .eff_n(
      xv, yv,
      t = if (!is.null(tt) && length(tt) == length(xv)) tt else NULL,,
      aggregate_by_time = effn_aggregate_by_time,
      time_tol = effn_time_tol,
      effn_method = effn_method,
      L = effn_L,
      nw_bw = effn_bw
    )
    list(neff = neff, nin = nin)
  }

  # --- 2) loop over pairs -----------------------------------------------------
  for (p in seq_len(npairs)) {
    i <- pair_mat[1, p]
    j <- pair_mat[2, p]

    r_list <- numeric(0)
    n_list <- integer(0)
    n_raw_list <- integer(0)     # 집계/전처리 전
    n_effin_list <- integer(0)   # neff 입력 길이(집계후 추정)
    s_list <- character(0)
    floor_count_pair <- 0L
    cap_count_pair <- 0L

    for (sb in subjects) {
      idx <- which(meta_df$subject == sb)
      if (length(idx) < min_n) next
      n_raw <- length(idx)  # 집계/정렬 이전 길이 기록

      o <- order(meta_df$time[idx])
      idx <- idx[o]

      # series (raw counts & RA)
      Xi_raw <- as.numeric(otu_mat[i, idx, drop = TRUE])  # counts
      Xj_raw <- as.numeric(otu_mat[j, idx, drop = TRUE])  # counts
      lib_t  <- as.numeric(lib[idx])                      # library size per time

      Xi_ra  <- as.numeric(ra_mat[i, idx, drop = TRUE])   # RA
      Xj_ra  <- as.numeric(ra_mat[j, idx, drop = TRUE])   # RA
      rest_ra <- pmax(0, 1 - Xi_ra - Xj_ra)               # RA
      tt <- meta_df$time[idx]

      cc <- stats::complete.cases(Xi_raw, Xj_raw, Xi_ra, Xj_ra, rest_ra, lib_t, tt)
      if (!any(cc)) next
      Xi_raw <- Xi_raw[cc]; Xj_raw <- Xj_raw[cc]
      Xi_ra  <- Xi_ra[cc];  Xj_ra  <- Xj_ra[cc]
      rest_ra <- rest_ra[cc]
      lib_t <- lib_t[cc]
      tt <- tt[cc]

      # i, j 모두 변이/비영점 체크
      if (var(Xi_ra, na.rm = TRUE) < 1e-16 || var(Xj_ra, na.rm = TRUE) < 1e-16) next
      if (mean(Xi_ra > 0) < nz_partner_min_frac || mean(Xj_ra > 0) < nz_partner_min_frac) next

      # ---- per-subject correlation ------------------------------------------
      if (transform == "alr") {
        # 1) ε_t 선택
        if (zero_mode_alr == "lib") {
          eps_t <- pmax(zero_lib_pseudo / pmax(lib_t, 1), 1e-12)
        } else if (zero_mode_alr == "minpos_time") {
          if (zero_minpos_base == "ij") {
            pos_min <- pmin(ifelse(Xi_ra > 0, Xi_ra, Inf),
                            ifelse(Xj_ra > 0, Xj_ra, Inf))
          } else {
            pos_min <- pmin(
              ifelse(Xi_ra > 0, Xi_ra, Inf),
              ifelse(Xj_ra > 0, Xj_ra, Inf),
              ifelse(rest_ra > 0, rest_ra, Inf)
            )
          }
          pos_min[!is.finite(pos_min)] <- NA_real_
          eps_t <- zero_minpos_alpha * ifelse(is.na(pos_min), 1e-8, pos_min)
          eps_t <- pmax(eps_t, 1e-12)
        } else { # "minpos_subject"
          base_pos <- if (zero_minpos_base == "ij")
            c(Xi_ra[Xi_ra > 0], Xj_ra[Xj_ra > 0])
          else
            c(Xi_ra[Xi_ra > 0], Xj_ra[Xj_ra > 0], rest_ra[rest_ra > 0])
          eps0 <- if (length(base_pos)) zero_minpos_alpha * min(base_pos) else 1e-8
          eps_t <- rep(pmax(eps0, 1e-12), length(Xi_ra))
        }

        # 2) zero-replacement + closure({i,j,rest})
        xi <- ifelse(Xi_ra > 0, Xi_ra, eps_t)
        xj <- ifelse(Xj_ra > 0, Xj_ra, eps_t)
        xr <- ifelse(rest_ra > 0, rest_ra, eps_t)

        s_pre <- xi + xj + xr
        s_pre[!is.finite(s_pre) | s_pre <= 0] <- 1
        xi <- xi / s_pre; xj <- xj / s_pre; xr <- xr / s_pre
        eps_star <- eps_t / s_pre

        # 3) rest floor
        if (is.finite(rest_floor_frac) && rest_floor_frac > 0) {
          rest_min <- rest_floor_frac * eps_star
          mask <- xr < rest_min
          if (any(mask)) {
            sc <- (1 - rest_min[mask]) / pmax(xi[mask] + xj[mask], .Machine$double.eps)
            xi[mask] <- xi[mask] * sc
            xj[mask] <- xj[mask] * sc
            xr[mask] <- rest_min[mask]
            floor_count_pair <- floor_count_pair + sum(mask)
          }
        }

        # 4) pair-to-rest ALR
        alr_i <- log(xi) - log(xr)
        alr_j <- log(xj) - log(xr)

        # 4.5) (옵션) 동일 시간 집계: 상관 & neff & 동적캡 모두 동일 시간해상도로 맞춤
        if (effn_aggregate_by_time != "none" && length(tt) > 1L && all(is.finite(tt))) {
          out <- .aggregate_by_time_series(
            mat = cbind(alr_i, alr_j, eps_star),
            tt  = tt,
            mode = effn_aggregate_by_time, tol = effn_time_tol
          )
          alr_i    <- out$mat[, 1]
          alr_j    <- out$mat[, 2]
          eps_star <- out$mat[, 3]
          tt       <- out$time
        }

        # 5) winsorization (alr_cap_mode)
        if (alr_cap_mode == "fixed" && is.finite(alr_cap_value)) {
          cap_before_i <- sum(abs(alr_i) > alr_cap_value, na.rm = TRUE)
          cap_before_j <- sum(abs(alr_j) > alr_cap_value, na.rm = TRUE)
          if (cap_before_i + cap_before_j > 0) {
            alr_i <- pmax(pmin(alr_i,  alr_cap_value), -alr_cap_value)
            alr_j <- pmax(pmin(alr_j,  alr_cap_value), -alr_cap_value)
            cap_count_pair <- cap_count_pair + cap_before_i + cap_before_j
          }
        } else if (alr_cap_mode == "dynamic") {
          # 안전장치: eps_star와 ALR 길이 불일치 시 동적 캡 스킵(경고만)
          if (length(eps_star) != length(alr_i) || length(eps_star) != length(alr_j)) {
            warning("length mismatch between eps_star and ALR; skipping dynamic cap for this subject.", call. = FALSE)
          } else {
            es  <- pmax(eps_star, .Machine$double.eps)
            num <- pmax(1 - 2 * es, .Machine$double.eps)
            B_theory <- log(num / es)                 # 이론 경계
            cap_vec  <- pmin(alr_cap_value, pmax(0, B_theory))  # 상한을 인자로
            cap_before_i <- sum(abs(alr_i) > cap_vec, na.rm = TRUE)
            cap_before_j <- sum(abs(alr_j) > cap_vec, na.rm = TRUE)
            if (cap_before_i + cap_before_j > 0) {
              alr_i <- pmax(pmin(alr_i, cap_vec), -cap_vec)
              alr_j <- pmax(pmin(alr_j, cap_vec), -cap_vec)
              cap_count_pair <- cap_count_pair + cap_before_i + cap_before_j
            }
          }
        } # "none"이면 아무 것도 안 함

        # 6) (옵션) 선형 detrend
        if (isTRUE(detrend) && length(alr_i) >= 3 && all(is.finite(tt))) {
          alr_i <- residuals(lm(alr_i ~ tt))
          alr_j <- residuals(lm(alr_j ~ tt))
        }

        # 7) 상관 & neff
        if (method == "spearman") {
          Ri <- rank(alr_i); Rj <- rank(alr_j)
          r_s <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
          if (!is.finite(r_s)) next
          ne <- .neff_and_nin(Ri, Rj, tt)
        } else {
          r_s <- suppressWarnings(stats::cor(alr_i, alr_j, method = "pearson"))
          if (!is.finite(r_s)) next
          ne <- .neff_and_nin(alr_i, alr_j, tt)
        }

      } else {
        # transform == "raw"
        xi <- Xi_ra; xj <- Xj_ra; xr <- rest_ra

        # (옵션) 동일 시간 집계: RAW 스케일에서도 일관화
        if (effn_aggregate_by_time != "none" && length(tt) > 1L && all(is.finite(tt))) {
          out <- .aggregate_by_time_series(
            mat = cbind(xi, xj, xr),
            tt  = tt,
            mode = effn_aggregate_by_time, tol = effn_time_tol
          )
          xi <- out$mat[, 1]; xj <- out$mat[, 2]; xr <- out$mat[, 3]; tt <- out$time
        }

        # (옵션) 선형 detrend (RA 스케일)
        if (isTRUE(detrend) && length(xi) >= 3 && all(is.finite(tt))) {
          xi <- residuals(lm(xi ~ tt))
          xj <- residuals(lm(xj ~ tt))
          xr <- residuals(lm(xr ~ tt))
        }

        if (use_partial) {
          var_xr <- stats::var(xr, na.rm = TRUE)
          if (!is.finite(var_xr) || var_xr < 1e-12 || length(xi) < 3L) {
            # fallback: simple correlation
            if (method == "spearman") {
              Ri <- rank(xi); Rj <- rank(xj)
              r_s <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
              if (!is.finite(r_s)) next
              ne <- .neff_and_nin(Ri, Rj, tt)
            } else {
              r_s <- suppressWarnings(stats::cor(xi, xj, method = "pearson"))
              if (!is.finite(r_s)) next
              ne <- .neff_and_nin(xi, xj, tt)
            }
          } else {
            if (method == "spearman") {
              Ri <- rank(xi, ties.method = "average")
              Rj <- rank(xj, ties.method = "average")
              Rr <- rank(xr, ties.method = "average")

              if (partial_method == "matrix") {
                r_ij <- suppressWarnings(stats::cor(Ri, Rj, method="pearson"))
                r_ir <- suppressWarnings(stats::cor(Ri, Rr, method="pearson"))
                r_jr <- suppressWarnings(stats::cor(Rj, Rr, method="pearson"))
                num <- r_ij - r_ir * r_jr
                den <- sqrt(max(1e-12, 1 - r_ir^2) * max(1e-12, 1 - r_jr^2))
                r_s <- if (den > 0) num / den else r_ij
                ne <- .neff_and_nin(Ri, Rj, tt)  # 랭크 스케일로 neff
              } else {
                res_i <- try(residuals(lm(Ri ~ Rr)), silent = TRUE)
                res_j <- try(residuals(lm(Rj ~ Rr)), silent = TRUE)
                if (inherits(res_i,"try-error") || inherits(res_j,"try-error")) {
                  r_s <- suppressWarnings(stats::cor(Ri, Rj, method = "pearson"))
                  ne <- .neff_and_nin(Ri, Rj, tt)
                } else {
                  r_s <- suppressWarnings(stats::cor(res_i, res_j, method = "pearson"))
                  ne <- .neff_and_nin(res_i, res_j, tt)  # 랭크-잔차 스케일로 neff
                }
              }
            } else {
              res_i <- try(residuals(lm(xi ~ xr)), silent = TRUE)
              res_j <- try(residuals(lm(xj ~ xr)), silent = TRUE)
              if (inherits(res_i,"try-error") || inherits(res_j,"try-error")) {
                r_s <- suppressWarnings(stats::cor(xi, xj, method = "pearson"))
                ne <- .neff_and_nin(xi, xj, tt)
              } else {
                r_s <- suppressWarnings(stats::cor(res_i, res_j, method = "pearson"))
                ne <- .neff_and_nin(res_i, res_j, tt)
              }
            }
          }
        } else {
          if (method == "spearman") {
            Ri <- rank(xi); Rj <- rank(xj)
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

      n_s <- length(xi) # 현재 사용 길이
      if (n_s < min_n) next

      # ---- AR(1)-lite: 약한 보정 적용 (φ-cap → blend → relative floor) ----
      # n_in: neff 계산의 입력 길이(집계 후), 없으면 현재 길이
      n_in <- ifelse(is.finite(ne$nin), as.numeric(ne$nin), as.numeric(n_s))
      # 기본 neff (보정 off 또는 실패 시 n_s 사용)
      n_eff_use <- if (isTRUE(acf_correction) && is.finite(ne$neff)) as.numeric(ne$neff) else as.numeric(n_s)

      if (isTRUE(acf_correction)) {
        # 1) φ-cap (AR(1)에서만): neff와 n_in으로 역산한 φ_hat을 상한으로 캡
        if (identical(effn_method, "ar1") &&
            !is.null(effn_phi_cap) && is.finite(effn_phi_cap) &&
            effn_phi_cap >= 0 && effn_phi_cap < 1 &&
            is.finite(n_eff_use) && is.finite(n_in) && n_in > 0) {
          # invert: n_eff = n * (1 - φ) / (1 + φ)  ⇒  φ_hat = (n - n_eff) / (n + n_eff)
          phi_hat <- (n_in - n_eff_use) / (n_in + n_eff_use)
          if (!is.finite(phi_hat)) phi_hat <- 0
          phi_hat <- max(0, min(0.999, phi_hat))              # 안정화
          phi_hat <- min(phi_hat, effn_phi_cap)               # cap
          n_eff_use <- n_in * (1 - phi_hat) / (1 + phi_hat)   # 재계산
        }
        # 2) blend toward raw n: n_eff <- (1-b) * n_in + b * n_eff
        if (!is.null(effn_blend)) {
          b <- as.numeric(effn_blend)
          if (is.finite(b) && b > 0 && b < 1) {
            n_eff_use <- (1 - b) * n_in + b * n_eff_use
          }
        }
        # 3) relative floor: n_eff >= f * n_in
        if (!is.null(effn_min_frac)) {
          f <- as.numeric(effn_min_frac)
          if (is.finite(f) && f > 0 && f <= 1) {
            n_eff_use <- max(n_eff_use, f * n_in)
          }
        }
      }
      # 절대 바닥 및 정수화
      n_eff_use <- max(3, floor(n_eff_use))

      r_list       <- c(r_list, r_s)
      n_list       <- c(n_list, n_eff_use)
      n_raw_list   <- c(n_raw_list, n_raw)
      n_effin_list <- c(n_effin_list, n_in)
      s_list       <- c(s_list, as.character(sb))
    } # end subject loop

    # Optional one-line warning per pair if safeguards triggered
    if (warn_alr && (floor_count_pair > 0 || cap_count_pair > 0)) {
      warning(sprintf("ALR safeguards applied for pair (%s, %s): rest-floor=%d, alr-cap=%d.",
                      i, j, floor_count_pair, cap_count_pair), call. = FALSE)
    }

    k_eff <- length(r_list)

    # --- 도우미: NA 행 빌더 ---------------------------------------------------
    build_na_row <- function() {
      tibble::tibble(
        i = i, j = j,
        method = method,
        transform = transform,
        partial_rest = use_partial,
        partial_method = if (use_partial && method == "spearman") partial_method else NA_character_,
        k = k_eff,
        r_pooled = NA_real_,
        ciL = NA_real_, ciU = NA_real_,
        piL = NA_real_, piU = NA_real_,
        pval = NA_real_,
        tau2 = NA_real_, I2 = NA_real_,
        Q = NA_real_, Q_p = NA_real_,
        n_min = ifelse(length(n_list) > 0, min(n_list), NA_integer_),
        n_median = ifelse(length(n_list) > 0, stats::median(n_list), NA_real_),
        n_max = ifelse(length(n_list) > 0, max(n_list), NA_integer_),
        n_raw_min = ifelse(length(n_raw_list) > 0, min(n_raw_list), NA_integer_),
        n_raw_median = ifelse(length(n_raw_list) > 0, stats::median(n_raw_list), NA_real_),
        n_raw_max = ifelse(length(n_raw_list) > 0, max(n_raw_list), NA_integer_),
        n_effin_min = ifelse(length(n_effin_list) > 0, min(n_effin_list), NA_integer_),
        n_effin_median = ifelse(length(n_effin_list) > 0, stats::median(n_effin_list), NA_real_),
        n_effin_max = ifelse(length(n_effin_list) > 0, max(n_effin_list), NA_integer_),
        meta_method = meta_method,
        knha = use_knha,
        var_method = if (var_method == "auto") (if (method == "spearman") "bonett" else "fisherz") else var_method,
        alr_cap_mode = if (transform == "alr") alr_cap_mode else NA_character_,
        alr_cap_value = if (transform == "alr" && alr_cap_mode == "fixed") alr_cap_value else NA_real_,
        effn_method = effn_method,
        effn_L = ifelse(effn_method == "bartlett", if (is.null(effn_L)) NA_real_ else effn_L, NA_real_),
        effn_bw = ifelse(effn_method == "nw", if (is.null(effn_bw)) NA_real_ else effn_bw, NA_real_),
        effn_phi_cap = ifelse(effn_method == "ar1", if (is.null(effn_phi_cap)) NA_real_ else effn_phi_cap, NA_real_),
        effn_blend = ifelse(is.null(effn_blend), NA_real_, effn_blend),
        effn_min_frac = ifelse(is.null(effn_min_frac), NA_real_, effn_min_frac),
        effn_aggregate_by_time = effn_aggregate_by_time,
        effn_time_tol = ifelse(is.null(effn_time_tol), NA_real_, effn_time_tol)
      )
    }

    if (k_eff < min_k) {
      meta_rows[[p]] <- build_na_row()
      if (return_subjectwise) {
        subjectwise_rows[[p]] <- tibble::tibble(
          i = i, j = j, subject = s_list,
          r = r_list, n = n_list,
          n_raw = n_raw_list, n_eff_in = n_effin_list
        )
      }
      next
    }

    # --- 3) random-effects meta-analysis -------------------------------------
    eff_var_method <- if (var_method == "auto") {
      if (method == "spearman") "bonett" else "fisherz"
    } else var_method

    # Build effect and variance on the correct scale
    if (method == "pearson") {
      yi <- atanh(r_list)
      vi <- 1 / (n_list - 3)
      scale <- "z"
    } else {
      yi <- atanh(r_list)
      if (eff_var_method == "bonett") {
        vi <- (1 + (r_list^2)/2) / (n_list - 3)
      } else if (eff_var_method == "fhp") {
        vi <- 1.06 / (n_list - 3)
      } else if (eff_var_method == "fisherz") {
        vi <- 1 / (n_list - 3)
      } else {
        stop("Unknown var_method for Spearman: ", eff_var_method)
      }
      scale <- "z"
    }

    fit <- try(metafor::rma.uni(
      yi = yi, vi = vi,
      method = meta_method,
      test = if (use_knha) "knha" else "z"
    ), silent = TRUE)

    if (inherits(fit, "try-error")) {
      meta_rows[[p]] <- build_na_row()
    } else {
      # 평균 효과 & CI
      if (scale == "z") {
        r_hat <- tanh(fit$b[1, 1])
        ciL <- tanh(fit$ci.lb); ciU <- tanh(fit$ci.ub)
        # 예측구간(PI)
        pr <- try(metafor::predict(fit, transf = tanh), silent = TRUE)
        if (inherits(pr, "try-error")) {
          piL <- NA_real_; piU <- NA_real_
        } else {
          piL <- pr$pi.lb; piU <- pr$pi.ub
        }
      } else {
        r_hat <- fit$b[1, 1]; ciL <- fit$ci.lb; ciU <- fit$ci.ub
        pr <- try(metafor::predict(fit), silent = TRUE)
        if (inherits(pr, "try-error")) {
          piL <- NA_real_; piU <- NA_real_
        } else {
          piL <- pr$pi.lb; piU <- pr$pi.ub
        }
      }

      meta_rows[[p]] <- tibble::tibble(
        i = i, j = j,
        method = method,
        transform = transform,
        partial_rest = use_partial,
        partial_method = if (use_partial && method == "spearman") partial_method else NA_character_,
        k = k_eff,
        r_pooled = r_hat,
        ciL = ciL, ciU = ciU,
        piL = piL, piU = piU,
        pval = fit$pval,
        tau2 = fit$tau2, I2 = fit$I2,
        Q = fit$QE, Q_p = fit$QEp,
        n_min = min(n_list),
        n_median = stats::median(n_list),
        n_max = max(n_list),
        n_raw_min = min(n_raw_list),
        n_raw_median = stats::median(n_raw_list),
        n_raw_max = max(n_raw_list),
        n_effin_min = min(n_effin_list),
        n_effin_median = stats::median(n_effin_list),
        n_effin_max = max(n_effin_list),
        meta_method = meta_method,
        knha = use_knha,
        var_method = eff_var_method,
        alr_cap_mode = if (transform == "alr") alr_cap_mode else NA_character_,
        alr_cap_value = if (transform == "alr" && alr_cap_mode == "fixed") alr_cap_value else NA_real_,
        effn_method = effn_method,
        effn_L = ifelse(effn_method == "bartlett", if (is.null(effn_L)) NA_real_ else effn_L, NA_real_),
        effn_bw = ifelse(effn_method == "nw", if (is.null(effn_bw)) NA_real_ else effn_bw, NA_real_),
        effn_phi_cap = ifelse(effn_method == "ar1", if (is.null(effn_phi_cap)) NA_real_ else effn_phi_cap, NA_real_),
        effn_blend = ifelse(is.null(effn_blend), NA_real_, effn_blend),
        effn_min_frac = ifelse(is.null(effn_min_frac), NA_real_, effn_min_frac),
        effn_aggregate_by_time = effn_aggregate_by_time,
        effn_time_tol = ifelse(is.null(effn_time_tol), NA_real_, effn_time_tol)
      )
    }

    if (return_subjectwise) {
      subjectwise_rows[[p]] <- tibble::tibble(
        i = i, j = j, subject = s_list,
        r = r_list, n = n_list,
        n_raw = n_raw_list, n_eff_in = n_effin_list
      )
    }
  } # end pair loop

  meta_tbl <- dplyr::bind_rows(meta_rows)
  rownames(meta_tbl) <- NULL

  if (nrow(meta_tbl)) {
    # -- helper: compute weights (mean-normalized to 1) --
    .get_weights <- function(df, by) {
      if (by == "none") return(rep(1, nrow(df)))
      w <- switch(by,
                  "n_median"        = df$n_median,
                  "n_effin_median"  = df$n_effin_median,
                  "k"               = df$k,
                  rep(1, nrow(df))
      )
      w[!is.finite(w) | w <= 0] <- 1
      w / mean(w, na.rm = TRUE)
    }

    if (q_method %in% c("wBH","wqvalue")) {
      w <- .get_weights(meta_tbl, q_weight_by)
      p_tilde <- pmin(1, meta_tbl$pval / w)  # weighted p-values
      meta_tbl$q_weight <- q_weight_by
      if (q_method == "wBH") {
        meta_tbl$qval <- stats::p.adjust(p_tilde, method = "BH")
        meta_tbl$pi0  <- NA_real_
      } else {  # "wqvalue"
        if (requireNamespace("qvalue", quietly = TRUE)) {
          qq <- try(qvalue::qvalue(p_tilde), silent = TRUE)
          if (!inherits(qq, "try-error")) {
            meta_tbl$qval <- qq$qvalues
            meta_tbl$pi0  <- qq$pi0
          } else {
            meta_tbl$qval <- stats::p.adjust(p_tilde, method = "BH")
            meta_tbl$pi0  <- NA_real_
          }
        } else {
          meta_tbl$qval <- stats::p.adjust(p_tilde, method = "BH")
          meta_tbl$pi0  <- NA_real_
        }
      }
    } else if (q_method == "qvalue") {
      if (requireNamespace("qvalue", quietly = TRUE)) {
        qq <- try(qvalue::qvalue(meta_tbl$pval), silent = TRUE)
        if (!inherits(qq, "try-error")) {
          meta_tbl$qval <- qq$qvalues
          meta_tbl$pi0  <- qq$pi0
        } else {
          meta_tbl$qval <- stats::p.adjust(meta_tbl$pval, method = "BH")
          meta_tbl$pi0  <- NA_real_
        }
      } else {
        meta_tbl$qval <- stats::p.adjust(meta_tbl$pval, method = "BH")
        meta_tbl$pi0  <- NA_real_
      }
    } else {
      meta_tbl$qval <- stats::p.adjust(meta_tbl$pval, method = q_method)
      meta_tbl$pi0  <- NA_real_
      meta_tbl$q_weight <- "none"
    }
    if (!is.null(k_filter)) meta_tbl <- dplyr::filter(meta_tbl, k >= k_filter)
    if (!is.null(r_abs_min)) meta_tbl <- dplyr::filter(meta_tbl, is.finite(r_pooled), abs(r_pooled) >= r_abs_min)
  }

  if (return_subjectwise) {
    subj_tbl <- dplyr::bind_rows(subjectwise_rows)
    return(list(meta = meta_tbl, subjectwise = subj_tbl))
  } else {
    return(meta_tbl)
  }
}
