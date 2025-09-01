#' Pairwise Bayesian gLV fitting (bidirectional) with OU/white-noise residuals
#'
#' @description
#' For each unordered taxon ``pair `{i, j}```, fits two ordered regressions (j -> i and i -> j)
#' with a Bayesian, gLV-inspired linear model. Residuals can be time-correlated via an
#' OU (approx AR(1)) process or set to white-noise. The response is the per-interval change
#' of the target's ALR (`Delta ALR_i / Delta t`, which can be standardized), and predictors
#' are **lagged states** of x_i, x_j (either ALR or raw relative abundance).
#'
#' @details
#' **Transforms**
#' - `transform = "alr"` (recommended): perform zero-aware replacement on the triplet
#'   (i, j, rest), then compute pair-to-rest ALR:
#'   `x_i = log(x_i / x_rest)`, `x_j = log(x_j / x_rest)`.
#'   Response is `y = Delta ALR_i / Delta t` (standardized if configured); predictors are **lagged** ALR.
#' - `transform = "raw"`: predictors use **lagged relative abundances**; the response remains
#'   `Delta ALR_i / Delta t`.
#'
#' **Validation**
#' - When `compute_elpd = TRUE`, runs **Repeated K-fold (K x R)** at the **subject level** and computes
#'   per-subject projection log-likelihood under the **fold-train posterior**, returning the **mean
#'   ELPD per subject**. **Train-only scaling** is used to avoid out-of-fold leakage.
#'
#' **ELPD computation modes (automatic defaults)**
#' - If `elpd_mode` is `NULL`, the scorer is chosen automatically:
#'   - `resid_mode = "ou"`  -> `kalman` (exact OU state-space marginal log-likelihood via a Kalman filter).
#'   - `resid_mode = "wn"`  -> `indep`  (i.i.d. Gaussian likelihood; no state-space).
#' - You can override with `elpd_mode = "kalman"` or `"indep"` explicitly.
#' - `kalman` (OU only): irregular-interval OU innovations are integrated as a **linear Gaussian state-space**
#'   model; measurement variance uses `sigma^2` (or a t->Gaussian variance expansion when `use_student_t = TRUE`).
#'   When required OU parameters are unavailable, it **falls back** to `indep`.
#' - `indep`: independent Gaussian scoring. For OU-like dispersion, predictive SD is
#'   `sd_pred = sqrt(sigma^2 + sd_ou^2)` (or `sigma_pred` if provided); for white-noise, `sd_pred = sqrt(sigma^2)`.
#'
#' **Noise model & scoring**
#' - Optional **Student-t** observation noise via `use_student_t = TRUE` (also used in K-fold scoring).
#' - For `elpd_mode = "kalman"` (OU), predictive uncertainty is handled **inside** the Kalman filter with
#'   OU transition `a_t = exp(-lambda * delta_t)` and process noise `q_t = sd_ou^2 * (1 - a_t^2)`;
#'   measurement variance uses `sigma^2` (or its t-variance expansion when `nu > 2`).
#'
#' **Sampling robustness**
#' - A diagnostics-aware **retry** policy (up to `max_retries`) monitors divergences, tree depth,
#'   and E-BFMI; if needed, it **escalates `nu_fix` within [4, 7]** while keeping other sampler controls
#'   stable for fold runs. Diagnostics are returned for downstream filtering.
#'
#' **Parallelism & progress**
#' - The **outer pair loop** can run in parallel with `n_workers_outer > 1` (PSOCK via **future/furrr**).
#' - **K-fold** can run in parallel with `n_workers_kfold > 1`. To avoid nested parallelism,
#'   when `n_workers_outer > 1` the effective K-fold workers are forced to **1**.
#' - If **progressr** is installed and `progress = "bar"`, progress bars are shown for both outer
#'   pairs and K-fold tasks; otherwise it falls back silently.
#'
#' **Interpretation**
#' - For edges, prioritize **signs** (PSP/LFSR derived from posterior draws) and **MCMC diagnostics**;
#'   optionally use **ELPD per subject** (from K x R) as supporting evidence. Aggregation is **subject-uniform**.
#'
#' @inheritParams get_sample_meta
#' @inheritParams precompute_spline_smoothed
#' @inheritParams .make_pair_inputs_glv
#'
#' @param physeq A \code{phyloseq} object with taxa in rows (will be transposed if needed).
#' @param taxa_vec Optional character vector of taxa IDs to subset; defaults to all taxa.
#' @param eps Small constant for legacy log-RA smoothing at the builder stage; default \code{1e-6}.
#' @param min_pairs Minimum number of valid time-adjacent pairs per direction; default \code{4}.
#'
#' @param chains,iter_warmup,iter_sampling,seed,quiet,adapt_delta,max_treedepth,metric,init
#' Sampling controls passed to \pkg{cmdstanr}.
#'
#' @param compute_elpd Logical; if \code{TRUE}, run **Repeated K-fold** (subject-level) and return mean
#' ELPD per subject; default \code{TRUE}.
#' @param transform One of \code{"alr"} (recommended) or \code{"raw"}.
#' @param lag Integer lag for predictors (within subject); default \code{1}.
#' @param resid_mode One of \code{"ou"}, \code{"wn"}; sets the residual process (OU or white-noise).
#' @param nz_partner_min_frac Drop subjects whose partner predictor is non-zero in fewer than this
#' fraction of rows; default \code{0.15}.
#' @param standardize_by_subject If \code{TRUE}, z-standardize within subject (train-only); default \code{TRUE}.
#'
#' @param zero_mode_alr Mode for zero-aware replacement in the ALR triplet; one of
#' \code{"minpos_time"}, \code{"minpos_subject"}, \code{"lib"}, \code{"fixed"}.
#' @param minpos_alpha Multiplier for minimal positive pseudo-count when applicable; default \code{0.5}.
#' @param minpos_base Basis for minimal positive detection, \code{"ij"} or \code{"triplet"}.
#' @param eps_fixed Fixed pseudo-count floor; default \code{1e-6}.
#' @param lib_eps_c Library-size scaled pseudo-count constant; default \code{0.65}.
#' @param rest_floor_frac Floor fraction for the rest component; default \code{1.0}.
#' @param alr_cap Hard cap for ALR magnitude; \code{Inf} to use adaptive theoretical caps.
#' @param smooth_scale Smoothing scale for optional per-subject spline on predictors;
#' one of \code{"logra"}, \code{"alr"}.
#' @param alr_spline_df,alr_spline_spar,alr_spline_cv Controls for spline smoothing if used.
#'
#' @param silent_sampler If \code{TRUE}, silence cmdstanr progress/messages in the wrapper; default \code{FALSE}.
#' @param progress Progress output: \code{"bar"}, \code{"verbose"}, or \code{"none"}.
#' When \code{"bar"} and \pkg{progressr} is available, outer/K-fold progress bars are displayed.
#' @param progress_every Update frequency for the legacy sequential progress bar; default \code{1}.
#'
#' @param use_student_t If \code{TRUE}, use Student-t observation noise (also used in K-fold scoring);
#' in `kalman` scoring, this is approximated via a variance expansion (`nu/(nu - 2)`) when `nu > 2`.
#'
#' @param elpd_mode Character or \code{NULL}. One of \code{"kalman"}, \code{"indep"}.
#' If \code{NULL} (recommended), defaults to \code{"kalman"} when \code{resid_mode = "ou"} and to
#' \code{"indep"} when \code{resid_mode = "wn"}. You may override explicitly if needed.
#'
#' @param kfold_K Number of folds K; default \code{5}.
#' @param kfold_R Number of repetitions R; default \code{3}.
#' @param kfold_seed Seed for K-fold split reproducibility (defaults to `seed`).
#' @param n_workers_kfold Number of parallel workers for K-fold (PSOCK via future/furrr); default \code{1}.
#' @param n_workers_outer Number of parallel workers for the outer pair loop; default \code{1}.
#' If \code{> 1}, K-fold parallelism is automatically disabled to avoid nested parallelism.
#'
#' @param max_retries Maximum sampler retry attempts (diagnostics-aware); default \code{3}.
#' @param nu_fix Initial fixed `nu` for Student-t; the wrapper may escalate up to 7 on retries; default \code{4}.
#'
#' @return
#' A \code{tibble} with one row per unordered ``pair `{i, j}```, summarizing **both directions**
#' (j -> i and i -> j) and self-effects: posterior summaries, sign probabilities (PSP/LFSR, reported as
#' two-sided \code{p_sign2}), key MCMC diagnostics, retry metadata, and (if enabled) **Repeated K-fold**
#' mean ELPD per subject plus basic fold counts. The ELPD scorer used is reported via \code{kfold$elpd_method}
#' (either \code{"kalman-ou"} or \code{"indep"}). Aggregation is **subject-uniform** by default.
#'
#' @section Progress & Parallel
#' - Set global handlers once for pretty bars: \preformatted{
#'   if (requireNamespace("progressr", quietly = TRUE)) {
#'     progressr::handlers(global = TRUE); progressr::handlers("cli")
#'   }}
#' - Reproducibility: seeds are deterministic per pair/fold (`seed`, `kfold_seed`) and
#'   \code{furrr::future_map(..., .options = furrr::furrr_options(seed = TRUE))} is used.
#'
#' @seealso \code{get_glv_model()}, \code{get_sample_meta()}, \code{precompute_spline_smoothed()},
#' \code{build_pair_df_smoothed()}, \code{.make_pair_inputs_glv()}.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("progressr", quietly = TRUE)) {
#'   progressr::handlers(global = TRUE); progressr::handlers("cli")
#' }
#' EP_res <- fit_glv_pairwise(
#'   physeq = EP_phy_obj, subject_col = "plot2", time_col = "week",
#'   transform = "alr", resid_mode = "ou",  # elpd_mode auto -> "kalman"
#'   use_student_t = TRUE,
#'   chains = 4, iter_warmup = 1500, iter_sampling = 1500,
#'   n_workers_outer = 4,   # outer pair loop in parallel
#'   n_workers_kfold = 1,   # k-fold parallel disabled when outer > 1
#'   progress = "bar"
#' )
#' }
#'
#' @importFrom progressr with_progress progressor
#' @export

fit_glv_pairwise <- function(# --- 필수 입력 ---
  physeq,
  subject_col,
  time_col,
  taxa_vec = NULL,

  # --- 모델/변환 설정 ---
  transform = c("alr", "raw"),
  lag = 1,
  resid_mode = c("ou", "wn"),
  use_student_t = TRUE,
  standardize_by_subject = TRUE,
  nz_partner_min_frac = 0.15,

  # --- 빌더 단계(로그-RA) 스무딩 & 최소 요구량 ---
  eps = 1e-6,
  spline_df = NULL,
  spline_spar = NULL,
  spline_cv = TRUE,
  min_unique_times = 3,
  min_pairs = 4,

  # --- ALR zero 처리 & (선택) ALR 스무딩 ---
  zero_mode_alr = c("minpos_time", "minpos_subject", "lib", "fixed"),
  minpos_alpha = 0.5,
  minpos_base = c("ij", "triplet"),
  eps_fixed = 1e-6,
  lib_eps_c = 0.65,
  rest_floor_frac = 1.0,
  alr_cap = Inf,
  smooth_scale = c("logra", "alr"),
  alr_spline_df = NULL,
  alr_spline_spar = NULL,
  alr_spline_cv = TRUE,

  # --- Stan 샘플링 제어 ---
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  seed = 1234,
  init = 0.2,
  adapt_delta = 0.98,
  max_treedepth = 14,
  metric = "diag_e",
  quiet = FALSE,

  # --- UX/진행 & 바깥 병렬 ---
  progress = c("bar", "verbose", "none"),
  progress_every = 1,
  silent_sampler = FALSE,
  n_workers_outer = 1L,

  # --- ELPD / K-fold ---
  compute_elpd = TRUE,
  elpd_mode = NULL,
  kfold_K = 5,
  kfold_R = 3,
  kfold_seed = seed,
  n_workers_kfold = 1L,

  # --- 리트라이 & t-꼬리 제어 ---
  max_retries = 3,
  nu_fix = 5)
{
  transform     <- match.arg(transform)
  resid_mode    <- match.arg(resid_mode)
  # elpd_mode 자동 디폴트: OU → "kalman", WN → "indep"
  if (is.null(elpd_mode)) {
    elpd_mode <- if (identical(resid_mode, "ou")) "kalman" else "indep"
  } else {
    elpd_mode <- match.arg(elpd_mode, c("kalman","indep"))
  }
  zero_mode_alr <- match.arg(zero_mode_alr)
  minpos_base   <- match.arg(minpos_base)
  smooth_scale  <- match.arg(smooth_scale)
  progress      <- match.arg(progress)

  # Stan model (nu_fixed + use_student_t 지원 가정)
  mod <- get_glv_model(quiet = quiet)

  # metadata & matrix
  meta_df <- get_sample_meta(physeq, subject_col, time_col)
  mat_rel <- get_abund_matrix_precomputed(physeq)
  if (is.null(taxa_vec))
    taxa_vec <- phyloseq::taxa_names(physeq)
  taxa_vec <- intersect(taxa_vec, rownames(mat_rel))
  if (length(taxa_vec) < 2)
    stop("taxa_vec must contain at least 2 taxa.")

  # legacy smoothing on log(RA+eps)
  sm_mat <- precompute_spline_smoothed(
    mat_rel,
    meta_df,
    taxa_vec,
    eps,
    spline_df,
    spline_spar,
    spline_cv,
    min_unique_times
  )

  # subject-level projection loglik (full) for any subset (K-fold scoring에 사용)
  .proj_loglik_subject <- function(draws_df,
                                   pair_in,
                                   resid_mode,
                                   use_t = FALSE,
                                   fit_full = NULL,
                                   nu_scalar = NULL,
                                   elpd_mode = NULL) {
    # 자동 디폴트: OU면 kalman, 아니면 indep
    if (is.null(elpd_mode)) {
      elpd_mode <- if (identical(resid_mode, "ou")) "kalman" else "indep"
    } else {
      elpd_mode <- match.arg(elpd_mode, c("indep","kalman"))
    }
    stopifnot(all(c("y", "xi", "xj", "subject") %in% names(pair_in)))
    subs <- unique(pair_in$subject)
    D    <- nrow(draws_df)
    Ssub <- length(subs)

    # 불규칙 간격 dt/subject-wise 인덱스 구성
    mk_prev_dt <- function(df) {
      df <- df[order(df$subject, df$time), , drop = FALSE]
      by_s <- split(seq_len(nrow(df)), df$subject)
      dtv <- numeric(nrow(df)); dtv[] <- NA_real_
      for (sb in names(by_s)) {
        ix <- by_s[[sb]]; dtv[ix[1]] <- NA_real_
        if (length(ix) >= 2L) for (k in 2:length(ix)) {
          dtv[ix[k]] <- as.numeric(df$time[ix[k]] - df$time[ix[k-1]])
        }
      }
      dtv
    }
    dt_all <- mk_prev_dt(pair_in)

    # 우선순위: (1) sigma_pred 생성량 → (2) sigma & sd_ou → (3) sigma, sigma_ou, lambda → (4) sigma
    if (resid_mode == "ou") {
      if ("sigma_pred" %in% names(draws_df)) {
        sigma_pred <- draws_df$sigma_pred
      } else if (all(c("sigma", "sd_ou") %in% names(draws_df))) {
        sigma_pred <- sqrt(draws_df$sigma^2 + draws_df$sd_ou^2)
      } else if (all(c("sigma", "sigma_ou", "lambda") %in% names(draws_df))) {
        sigma_pred <- sqrt(draws_df$sigma^2 + draws_df$sigma_ou^2 / (2 * pmax(draws_df$lambda, 1e-8)))
      } else {
        sigma_pred <- draws_df$sigma
      }
    } else {
      # "wn"
      sigma_pred <- draws_df$sigma
    }

    nu_vec <- NULL
    if (isTRUE(use_t)) {
      if ("nu" %in% names(draws_df)) {
        nu_vec <- draws_df$nu                     # 1) prefer estimated df (if modeled)
      } else if (!is.null(nu_scalar) && is.finite(nu_scalar)) {
        nu_vec <- rep(nu_scalar, nrow(draws_df))  # 2) fallback to fixed df from wrapper
      } else if (!is.null(fit_full)) {
        nu_draws <- try(fit_full$draws("nu"), silent = TRUE)
        if (!inherits(nu_draws, "try-error") &&
            !is.null(nu_draws)) {
          nu_vec <- as.numeric(posterior::as_draws_df(nu_draws)[["nu"]])
        }
      }
    }
    # If df unavailable, fall back to Normal scoring downstream
    use_t <- isTRUE(use_t) && !is.null(nu_vec)

    loglik_full <- matrix(NA_real_, nrow = D, ncol = Ssub)
    n_obs_vec   <- integer(Ssub)

    # --- 칼만필터 기반 OU 주변우도 (resid_mode="ou" & elpd_mode="kalman") ---
    .ou_kf_ll_one_subject <- function(draws_df, y, xi, xj, tvec) {
      J <- length(y)
      if (J == 0) return(rep(NA_real_, nrow(draws_df)))
      # mean: mu = r0 + a_ii * xi + a_ij * xj
      r0_mat <- matrix(draws_df$r0, nrow = nrow(draws_df), ncol = J)
      aii_xi <- tcrossprod(draws_df$a_ii, xi)
      aij_xj <- tcrossprod(draws_df$a_ij, xj)
      mu     <- r0_mat + aii_xi + aij_xj     # (D x J)

      y_rep  <- matrix(rep(y, each = nrow(draws_df)), nrow = nrow(draws_df))

      # 측정잡음 분산 R (t면 ν/(ν-2) 팩터로 가우시안 근사)
      if (isTRUE(use_t)) {
        t_fac <- ifelse(nu_vec > 2, nu_vec/(nu_vec - 2), 5.0)
        R <- (draws_df$sigma^2) * t_fac
      } else {
        R <- (draws_df$sigma^2)
      }
      R <- pmax(R, 1e-10)

      # OU 이산화: a_t, q_t
      have_lambda <- "lambda" %in% names(draws_df)
      have_phi    <- "phi"    %in% names(draws_df)
      have_sdou   <- "sd_ou"  %in% names(draws_df)

      # stationary var of OU
      if (have_sdou) {
        sd2 <- (draws_df$sd_ou)^2
      } else if (all(c("sigma_ou","lambda") %in% names(draws_df))) {
        sd2 <- (draws_df$sigma_ou^2) / (2 * pmax(draws_df$lambda, 1e-8))
      } else {
        # 정보 부족 → 독립 근사 fallback
        return(matrixStats::rowSums2(
          dnorm(y_rep, mean = mu, sd = matrix(sigma_pred, nrow = nrow(draws_df), ncol = J), log = TRUE)
        ))
      }

      # 필터 초기화: e0 ~ N(0, sd2)
      D <- nrow(draws_df)
      ll <- rep(0.0, D)

      ll <- rep(0.0, D)

      use_ri <- ("tau_r" %in% names(draws_df)) && any(is.finite(draws_df$tau_r))
      if (!use_ri) {
        m  <- rep(0.0, D)     # e-state mean per draw
        P  <- sd2             # e-state var per draw
      } else {
        # 2D state: [e_t, r_s]^T, with r_s time-invariant
        m_e <- rep(0.0, D); m_r <- rep(0.0, D)
        P_ee <- sd2
        P_rr <- pmax(draws_df$tau_r^2, 0)
        P_er <- rep(0.0, D)
      }

      # 시간 루프
      # 주의: 첫 관측의 dt가 NA일 수 있으므로 dt<=0 → 작은 값으로 보정
      dtv <- as.numeric(c(NA, diff(tvec)))
      l2pi <- log(2*pi)
      for (t in seq_len(J)) {
        dt <- dtv[t]
        if (!is.finite(dt) || dt < 0) dt <- 0        # 첫 관측은 dt=0로 처리
        if (have_lambda) {
          lam <- pmax(draws_df$lambda, 1e-8)
          a_t <- exp(-lam * dt)
        } else if (have_phi) {
          phi <- pmin(pmax(draws_df$phi, 1e-8), 0.999999)
          lam <- -log(phi)
          a_t <- exp(-lam * dt)
        } else {
          a_t <- rep(0.0, D)
        }
        q_t <- pmax(sd2 * (1 - a_t^2), 1e-16)
        if (!use_ri) {
          # ---- 1D (OU only) ----
          # 예측
          m_pred <- a_t * m
          P_pred <- (a_t^2) * P + q_t
          # 혁신
          v_t <- y_rep[, t] - mu[, t] - m_pred
          F_t <- pmax(P_pred + R, 1e-10)
          # loglik
          ll <- ll + (-0.5 * (l2pi + log(F_t) + (v_t * v_t) / F_t))
          # 업데이트
          K_t <- P_pred / F_t
          m   <- m_pred + K_t * v_t
          P   <- pmax((1 - K_t) * P_pred, 1e-12)
        } else {
          # ---- 2D (OU + random intercept) ----
          # 예측
          m_e_pred <- a_t * m_e
          m_r_pred <- m_r                   # r_s is time-invariant
          P_ee_p <- (a_t^2) * P_ee + q_t
          P_er_p <- a_t * P_er
          P_rr_p <- P_rr                    # no process noise on r_s
          # 혁신
          v_t <- y_rep[, t] - mu[, t] - (m_e_pred + m_r_pred)
          F_t <- pmax(P_ee_p + P_rr_p + 2*P_er_p + R, 1e-10)  # H=[1,1]
          # loglik
          ll <- ll + (-0.5 * (l2pi + log(F_t) + (v_t * v_t) / F_t))
          # 칼만 이득 (2x1)
          K_e <- (P_ee_p + P_er_p) / F_t
          K_r <- (P_er_p + P_rr_p) / F_t
          # 업데이트
          m_e <- m_e_pred + K_e * v_t
          m_r <- m_r_pred + K_r * v_t
          # 공분산 업데이트: P = (I-KH)P_pred
          # 요소별 전개 (H=[1,1])
          P_ee <- pmax(P_ee_p - K_e*(P_ee_p + P_er_p), 1e-12)
          P_rr <- pmax(P_rr_p - K_r*(P_er_p + P_rr_p), 1e-12)
          P_er <-       P_er_p - K_e*(P_er_p + P_rr_p)
        }
      }
      ll
    }

    for (s_idx in seq_along(subs)) {
      sb <- subs[s_idx]
      ix <- which(pair_in$subject == sb)
      y  <- pair_in$y[ix]
      xi <- pair_in$xi[ix]
      xj <- pair_in$xj[ix]
      tvec <- pair_in$time[ix]
      J  <- length(ix)
      n_obs_vec[s_idx] <- J

      if (identical(resid_mode, "ou") && identical(elpd_mode, "kalman")) {
        ll_f <- .ou_kf_ll_one_subject(draws_df, y, xi, xj, tvec)
      } else {
        # 독립 근사(레거시)
        r0_mat  <- matrix(draws_df$r0, nrow = D, ncol = J)
        aii_xi  <- tcrossprod(draws_df$a_ii, xi)
        aij_xj  <- tcrossprod(draws_df$a_ij, xj)
        mu_full <- r0_mat + aii_xi + aij_xj
        sig_mat <- matrix(sigma_pred, nrow = D, ncol = J)
        y_mat   <- matrix(rep(y, each = D), nrow = D, ncol = J)
        if (!is.null(nu_vec)) {
          nu_mat <- matrix(nu_vec, nrow = D, ncol = J)
          ll_f <- rowSums(stats::dt((y_mat - mu_full) / sig_mat,
                                    df = nu_mat,
                                    log = TRUE
          ) - log(sig_mat))
        } else {
          ll_f <- rowSums(stats::dnorm(
            y_mat,
            mean = mu_full,
            sd = sig_mat,
            log = TRUE
          ))
        }
      }
      loglik_full[, s_idx] <- ll_f
    }
    list(full = loglik_full, subjects = subs, n_obs = n_obs_vec)
  }

  # ---------- Repeated K-fold helpers ----------
  .make_repkfold_splits <- function(subject_vec,
                                    K = 5,
                                    R = 3,
                                    seed = 123) {
    set.seed(seed)
    subs <- unique(subject_vec)
    S <- length(subs)
    if (K > S)
      stop("K > #subjects")
    reps <- vector("list", R)
    for (r in seq_len(R)) {
      ord <- sample.int(S)
      folds <- split(subs[ord], rep(1:K, length.out = S))
      reps[[r]] <- lapply(seq_len(K), function(k) {
        list(test_subjects = folds[[k]],
             train_subjects = setdiff(subs, folds[[k]]))
      })
    }
    reps
  }

  .build_train_test <- function(pair_in,
                                train_subjects,
                                test_subjects,
                                use_global_scaling = TRUE) {
    mk_prev_dt <- function(df) {
      df <- df[order(df$subject, df$time), , drop = FALSE]
      N <- nrow(df)
      prev <- integer(N)
      dtv <- numeric(N)
      by_s <- split(seq_len(N), df$subject)
      for (sb in names(by_s)) {
        ix <- by_s[[sb]]
        prev[ix[1]] <- 0L
        dtv[ix[1]] <- 0
        if (length(ix) >= 2L) {
          for (k in 2:length(ix)) {
            prev[ix[k]] <- ix[k - 1]
            dt_k <- as.numeric(df$time[ix[k]] - df$time[ix[k - 1]])
            dtv[ix[k]] <- if (is.finite(dt_k) &&
                              dt_k > 0)
              dt_k
            else
              1e-6
          }
        }
      }
      list(df = df,
           prev = prev,
           dt = dtv)
    }
    tr <- subset(pair_in, subject %in% train_subjects)
    te <- subset(pair_in, subject %in% test_subjects)
    if (!nrow(tr) || !nrow(te))
      return(NULL)
    trd <- mk_prev_dt(tr)
    ted <- mk_prev_dt(te)
    # --- ★ Train-only global scaling (OOF leakage guard for K-fold) ---
    if (isTRUE(use_global_scaling)) {
      sc <- list(
        y_m  = mean(trd$df$y),
        y_s  = sd(trd$df$y),
        xi_m = mean(trd$df$xi),
        xi_s = sd(trd$df$xi),
        xj_m = mean(trd$df$xj),
        xj_s = sd(trd$df$xj)
      )
      zs <- function(x, m, s)
        (x - m) / ifelse(is.finite(s) && s > 0, s, 1)
      for (nm in c("y", "xi", "xj")) {
        trd$df[[nm]] <- zs(trd$df[[nm]], sc[[paste0(nm, "_m")]], sc[[paste0(nm, "_s")]])
        ted$df[[nm]] <- zs(ted$df[[nm]], sc[[paste0(nm, "_m")]], sc[[paste0(nm, "_s")]])
      }
    }

    # --- Guard: insufficient train samples per fold ---
    if (nrow(trd$df) < min_pairs)
      return(NULL)

    list(
      train = trd$df,
      test = ted$df,
      prev_train = trd$prev,
      dt_train = trd$dt,
      prev_test  = ted$prev,
      dt_test  = ted$dt
    )
  }

  .fold_fit_and_score <- function(mod,
                                  stan_list_base,
                                  sample_args_base,
                                  pair_in,
                                  tr,
                                  te,
                                  resid_mode,
                                  use_t = FALSE,
                                  elpd_mode = NULL,
                                  silent_sampler = TRUE,
                                  nu_fixed_override = NULL,
                                  sample_args_override = NULL,
                                  freeze_retry_hypers = FALSE,
                                  seed_override = NULL) {
    # 자동 디폴트: OU면 kalman, 아니면 indep
    if (is.null(elpd_mode)) {
      elpd_mode <- if (identical(resid_mode, "ou")) "kalman" else "indep"
    } else {
      elpd_mode <- match.arg(elpd_mode, c("indep","kalman"))
    }
    pts <- .build_train_test(pair_in, tr, te, use_global_scaling = TRUE)
    if (is.null(pts))
      return(NULL)

    sl <- stan_list_base
    sl$N   <- nrow(pts$train)
    sl$y   <- pts$train$y
    sl$xi  <- pts$train$xi
    sl$xj  <- pts$train$xj
    sl$S   <- length(unique(pts$train$subject))
    sl$sid <- as.integer(factor(pts$train$subject))
    sl$prev <- pts$prev_train
    sl$dt  <- pts$dt_train

    # K-fold에서도 ν 고정 주입을 허용 (훈련 진단 기반 글로벌 규칙; 누출 아님)
    # 이유: override 값은 오직 "훈련 단계의 샘플러 진단"으로만 결정되며,
    # 테스트 ELPD/관측에는 의존하지 않는다. 폴드 간 동일 ν로 비교 일관성 확보.

    # ★ 샘플러 제어는 공유하되, 적응 산출물(step_size/metric)은 폴드마다 재적응
    sa <- if (is.null(sample_args_override))
      sample_args_base
    else
      sample_args_override
    sa$data <- sl
    if (!is.null(seed_override))
      sa$seed <- seed_override   # 폴드별 결정적 seed
    # 폴드마다 재적응 유도: full-data 적응 결과 제거
    sa$step_size   <- NULL
    sa$inv_metric  <- NULL
    sa$metric_file <- NULL
    # ★ 풀런에서 넘어온 init(closure) 오염 차단: 폴드 데이터에 맞춰 재생성 또는 제거
    if (!(is.numeric(sa$init) &&
          length(sa$init) == 1L && is.finite(sa$init))) {
      sa$init <- NULL
    }

    # ★ 폴드 학습은 항상 체인 병렬 끔 (outer/메인 run과 중첩 병렬 방지)
    sa$parallel_chains <- 1L
    sa$chains <- if (is.null(sa$chains))
      4L
    else {
      x <- sa$chains
      if (is.character(x))
        x <- suppressWarnings(as.numeric(x))
      ok <- (length(x) == 1L) &&
        is.finite(x) && (x >= 1) && (floor(x) == x)
      if (ok)
        as.integer(x)
      else
        4L
    }

    tryfit <- .sample_with_retry(
      mod        = mod,
      base_args  = sa,
      stan_list  = sl,
      max_retries = max_retries,
      # ν 고정 정책: override가 주어지면 그 값으로 **고정**(bump 금지)
      # (훈련 진단만으로 결정된 사전 규칙 → 누출 아님)
      nu_fix      = if (is.null(nu_fixed_override)) nu_fix else nu_fixed_override,
      max_nu_cap  = if (is.null(nu_fixed_override)) 7      else nu_fixed_override,
      ebfmi_thresh = 0.30,
      silent_sampler = silent_sampler,
      tag = "kfold",
      freeze_retry_hypers = freeze_retry_hypers   # ★ 하이퍼를 바꾸지 않는 리트라이
    )
    fit <- tryfit$fit

    # 방어적 동일성 확인(디버그용; 필요시 주석 처리)
    fa <- tryfit$final_args

    d <- safe_draws_df(fit)
    te_df <- pts$test
    # Resolve df for Student-t scoring without relying on rlang::`%||%`
    use_t_flag <- isTRUE(use_t)
    nu_for_score <- NULL
    if (use_t_flag) {
      if (!is.null(tryfit$used_nu) && is.finite(tryfit$used_nu)) {
        nu_for_score <- tryfit$used_nu
      } else if (!is.null(stan_list_base$nu_fixed) &&
                 is.finite(stan_list_base$nu_fixed)) {
        nu_for_score <- stan_list_base$nu_fixed
      } else if (exists("nu_fix") && is.finite(nu_fix)) {
        nu_for_score <- nu_fix
      }
    }
    llm <- .proj_loglik_subject(
      draws_df   = d,
      pair_in    = te_df,
      resid_mode = resid_mode,
      use_t      = use_t_flag,
      fit_full   = fit,
      nu_scalar  = nu_for_score,
      elpd_mode  = elpd_mode
    )
    # --- NEW: subject별 ELPD 벡터와 fold diagnostics 구성 ---
    # llm$full: (draws x n_test_subjects) 행렬, llm$subjects: 테스트 subject 벡터
    if (is.null(llm) ||
        is.null(llm$full) ||
        !is.matrix(llm$full) || ncol(llm$full) == 0) {
      # 테스트 데이터가 비정상인 폴드는 실패로 처리
      return(NULL)
    }
    elpd_vec <- stats::setNames(apply(llm$full, 2, .log_mean_exp), llm$subjects)
    elpd_ppd_vec <- elpd_vec / pmax(llm$n_obs, 1L)

    dg <- tryfit$diag
    fd <- data.frame(
      n_retries      = tryfit$n_retries,
      nu_used        = tryfit$used_nu,
      ebfmi_min      = if (!is.null(dg$ebfmi_min))
        dg$ebfmi_min
      else
        NA_real_,
      worst_rhat     = if (!is.null(dg$worst_rhat))
        dg$worst_rhat
      else
        NA_real_,
      min_ess_bulk   = if (!is.null(dg$min_ess_bulk))
        dg$min_ess_bulk
      else
        NA_real_,
      treedepth_hits = if (!is.null(dg$n_treedepth_hit))
        dg$n_treedepth_hit
      else
        NA_integer_,
      n_divergent    = if (!is.null(dg$n_divergent))
        dg$n_divergent
      else
        NA_integer_
    )

    list(elpd = elpd_vec, fold_diag = fd)
    list(elpd = elpd_vec, elpd_ppd = elpd_ppd_vec, fold_diag = fd)
  }

  .repkfold_eval <- function(mod,
                             stan_list_base,
                             sample_args_base,
                             pair_in,
                             K = 5,
                             R = 3,
                             seed = 123,
                             resid_mode,
                             use_t = FALSE,
                             elpd_mode = NULL,
                             silent_sampler = TRUE,
                             n_workers_kfold = 1L,
                             nu_fixed_kfold = NULL,         # ← 1차 K-fold에서 고정할 ν (예: 본계산 used_nu)
                             freeze_retry_hypers = FALSE) { # ← 폴드 내 bump 금지 여부
    # 자동 디폴트: OU → "kalman", WN → "indep"
    if (is.null(elpd_mode)) {
      elpd_mode <- if (identical(resid_mode, "ou")) "kalman" else "indep"
    } else {
      elpd_mode <- match.arg(elpd_mode, c("indep","kalman"))
    }
    splits <- .make_repkfold_splits(pair_in$subject, K, R, seed)
    # --- compact split manifest for reproducibility (r,k,test_subjects)
    splits_df <- ({
      rows <- vector("list", K * R)
      t <- 0L
      for (r in seq_len(R)) {
        folds <- splits[[r]]
        for (k in seq_len(K)) {
          t <- t + 1L
          rows[[t]] <- data.frame(r = r,
                                  k = k,
                                  test_subjects = I(list(folds[[k]]$test_subjects)))
        }
      }
      do.call(rbind, rows)
    })
    # subject universe & how many times each was held out as test
    subs_all <- sort(unique(pair_in$subject))
    test_counts_tab <- table(unlist(splits_df$test_subjects))
    subject_test_counts <- as.integer(test_counts_tab[subs_all])
    names(subject_test_counts) <- subs_all

    # task 리스트: (r, k, train_subjects, test_subjects, seed_offset)
    tasks <- list()
    for (r in seq_len(R)) {
      folds <- splits[[r]]
      for (k in seq_len(K)) {
        tr <- folds[[k]]$train_subjects
        te <- folds[[k]]$test_subjects
        base_seed <- if (!is.null(sample_args_base$seed))
          as.integer(sample_args_base$seed)
        else
          as.integer(seed)
        if (is.na(base_seed))
          base_seed <- 1L
        tasks[[length(tasks) + 1L]] <- list(
          r = r,
          k = k,
          tr = tr,
          te = te,
          seed = as.integer(base_seed + 1000L * r + k)
        )
      }
    }

    # 폴드별 실행 함수
    run_task <- function(task) {
      .fold_fit_and_score(
        mod,
        stan_list_base,
        sample_args_base,
        pair_in,
        task$tr,
        task$te,
        resid_mode,
        use_t,
        elpd_mode = elpd_mode,
        silent_sampler,
        nu_fixed_override   = nu_fixed_kfold,     # ← ν 고정 주입
        freeze_retry_hypers = freeze_retry_hypers,
        seed_override = task$seed
      )
    }

    # 실행: 순차 또는 병렬
    if (n_workers_kfold > 1L) {
      # 안전하게 PSOCK(멀티세션) 사용 권장 (fork 이슈 회피)
      op <- future::plan()
      on.exit(future::plan(op), add = TRUE)
      future::plan(future::multisession, workers = n_workers_kfold)
      if (has_progressr && identical(progress, "bar")) {
        progressr::with_progress({
          p <- progressr::progressor(steps = length(tasks))
          res_list <- furrr::future_map(tasks, function(task) {
            res <- run_task(task)
            p(message = sprintf("kfold r=%d k=%d", task$r, task$k))
            res
          }, .options = furrr::furrr_options(seed = TRUE))
        })
      } else {
        res_list <- furrr::future_map(tasks, run_task, .options = furrr::furrr_options(seed = TRUE))
      }
    } else {
      res_list <- purrr::map(tasks, run_task)
    }

    subs <- sort(unique(pair_in$subject))
    agg      <- setNames(numeric(length(subs)), subs)
    cnt      <- setNames(integer(length(subs)), subs)
    agg_ppd  <- setNames(numeric(length(subs)), subs)
    cnt_ppd  <- setNames(integer(length(subs)), subs)

    fold_diag_df <- list()

    n_ok <- 0L
    n_fail <- 0L

    for (el in res_list) {
      if (is.null(el)) {
        n_fail <- n_fail + 1L
      } else {
        n_ok <- n_ok + 1L
        # elpd
        idx <- names(el$elpd)
        agg[idx]     <- agg[idx]     + el$elpd
        cnt[idx]     <- cnt[idx]     + 1L
        if (!is.null(el$elpd_ppd)) {
          agg_ppd[idx] <- agg_ppd[idx] + el$elpd_ppd
          cnt_ppd[idx] <- cnt_ppd[idx] + 1L
        }
        # fold diag
        fold_diag_df[[length(fold_diag_df) + 1L]] <- data.frame(
          n_retries      = el$fold_diag$n_retries,
          nu_used        = el$fold_diag$nu_used,
          ebfmi_min      = el$fold_diag$ebfmi_min,
          worst_rhat     = el$fold_diag$worst_rhat,
          min_ess_bulk   = el$fold_diag$min_ess_bulk,
          treedepth_hits = el$fold_diag$treedepth_hits,
          n_divergent    = el$fold_diag$n_divergent
        )
      }
    }

    elpd_subject     <- agg / pmax(cnt, 1L)
    elpd_subject_ppd <- agg_ppd / pmax(cnt_ppd, 1L)
    elpd_mean    <- if (length(elpd_subject)) mean(elpd_subject[is.finite(elpd_subject)], na.rm = TRUE) else NA_real_
    elpd_mean_ppd<- if (length(elpd_subject_ppd)) mean(elpd_subject_ppd[is.finite(elpd_subject_ppd)], na.rm = TRUE) else NA_real_

    elpd_mean    <- if (length(elpd_subject)) {
      m <- mean(elpd_subject[is.finite(elpd_subject)], na.rm = TRUE)
      if (is.finite(m)) m else NA_real_
    } else NA_real_


    fd <- if (length(fold_diag_df))
      do.call(rbind, fold_diag_df)
    else
      data.frame(
        n_retries = integer(),
        nu_used = integer(),
        ebfmi_min = double(),
        worst_rhat = double(),
        min_ess_bulk = double(),
        treedepth_hits = integer(),
        n_divergent = integer()
      )
    # 리트라이 요약
    retry_total <- if (nrow(fd))
      sum(fd$n_retries, na.rm = TRUE)
    else
      0L
    retry_mean  <- if (nrow(fd))
      mean(fd$n_retries, na.rm = TRUE)
    else
      NA_real_
    nu_counts   <- if (nrow(fd))
      as.integer(table(fd$nu_used))
    else
      integer()
    names(nu_counts) <- if (nrow(fd))
      names(table(fd$nu_used))
    else
      character()


    list(
      type = "repeated-kfold",
      K = K,
      R = R,
      elpd_subject = elpd_subject,
      elpd_mean    = elpd_mean,
      elpd_subject_ppd = elpd_subject_ppd,
      elpd_mean    = elpd_mean,
      elpd_mean_ppd= elpd_mean_ppd,
      elpd_method  = if (identical(elpd_mode,"kalman") && identical(resid_mode,"ou")) "kalman-ou" else "indep",
      n_folds_ok = n_ok,
      n_folds_fail = n_fail,
      # reproducibility payload for pseudo-BMA/stacking
      splits_df = splits_df,
      subject_ids = subs_all,
      subject_test_counts = subject_test_counts,
      fold_diag = fd,
      retry_total = retry_total,
      retry_mean  = retry_mean,
      nu_used_counts = nu_counts
    )
  }

  # --------------------------------------------------------
  # --- Avoid nested parallelism: if outer-parallel, disable kfold-parallel
  n_workers_kfold_eff <- if (n_workers_outer > 1L)
    1L
  else
    n_workers_kfold

  # progressr availability (optional)
  has_progressr <- requireNamespace("progressr", quietly = TRUE)

  # run one ordered direction: partner j -> target i
  run_one <- function(target,
                      partner,
                      seed_override = NULL,
                      progress_local = progress) {
    # 안전 초기화(ELPD 비활성 시도 포함)
    kfold_outer_rounds_local <- 0L

    # 0) 페어 데이터 구성 ---------------------------------------------------------
    pair_df <- build_pair_df_smoothed(sm_mat, meta_df, partner, target, eps, min_pairs)
    if (is.null(pair_df))
      return(NULL)
    std_by_subj_eff <- if (isTRUE(compute_elpd))
      FALSE
    else
      standardize_by_subject
    if (isTRUE(compute_elpd) && isTRUE(standardize_by_subject)) {
      warning(
        "compute_elpd=TRUE with standardize_by_subject=TRUE: ",
        "disabling subject-wise standardization to avoid OOF leakage; ",
        "K-fold will apply train-only scaling."
      )
    }

    pair_in <- .make_pair_inputs_glv(
      pair_df = pair_df,
      transform = transform,
      lag = lag,
      zero_mode_alr = zero_mode_alr,
      minpos_alpha = minpos_alpha,
      minpos_base  = minpos_base,
      eps_fixed = eps_fixed,
      lib_eps_c = lib_eps_c,
      rest_floor_frac = rest_floor_frac,
      alr_cap = 12,
      smooth_scale = smooth_scale,
      alr_spline_df = alr_spline_df,
      alr_spline_spar = alr_spline_spar,
      alr_spline_cv = alr_spline_cv,
      nz_partner_min_frac = nz_partner_min_frac,
      standardize_by_subject = std_by_subj_eff
    )
    if (is.null(pair_in) || nrow(pair_in) < min_pairs)
      return(NULL)

    # 결측/비유한 제거
    pair_in <- pair_in[is.finite(pair_in$y) &
                         is.finite(pair_in$xi) &
                         is.finite(pair_in$xj) &
                         is.finite(pair_in$time), , drop = FALSE]
    if (!nrow(pair_in))
      return(NULL)

    # 1) prev/dt 구성(AR(1)용) ----------------------------------------------------
    pair_in <- pair_in[order(pair_in$subject, pair_in$time), , drop = FALSE]
    N <- nrow(pair_in)
    prev <- integer(N)
    dtv <- numeric(N)
    by_s <- split(seq_len(N), pair_in$subject)
    for (sb in names(by_s)) {
      ix <- by_s[[sb]]
      prev[ix[1]] <- 0L
      dtv[ix[1]] <- 0
      if (length(ix) >= 2L) {
        for (k in 2:length(ix)) {
          prev[ix[k]] <- ix[k - 1]
          dt_k <- as.numeric(pair_in$time[ix[k]] - pair_in$time[ix[k - 1]])
          dtv[ix[k]] <- if (is.finite(dt_k) &&
                            dt_k > 0)
            dt_k
          else
            1e-6
        }
      }
    }

    n_subj <- length(unique(pair_in$subject))
    if (n_subj < 2L) {
      warning(
        "Fewer than 2 subjects for pair {",
        partner,
        "->",
        target,
        "} — AR(1) persistence may be weakly identified."
      )
    }

    sedf  <- attr(pair_in, "smooth_edf_mean")
    ssc   <- attr(pair_in, "smooth_scale")
    sflag <- isTRUE(attr(pair_in, "smoothed"))

    # 2) Stan data (기본) ---------------------------------------------------------
    stan_list <- c(
      list(
        N   = N,
        y   = pair_in$y,
        xi  = pair_in$xi,
        xj  = pair_in$xj,
        S   = n_subj,
        sid = as.integer(factor(pair_in$subject)),
        prev = prev,
        dt   = dtv,
        use_student_t = as.integer(isTRUE(use_student_t)),
        # 0/1
        nu_fixed      = nu_fix                              # 시도별 래퍼에서 업데이트됨
      )
    )

    # 앞쪽 값 우선으로 중복 키 제거
    nm <- names(stan_list)
    if (anyDuplicated(nm))
      stan_list <- stan_list[match(unique(nm), nm)]

    # --- seed 결정(쌍별/방향별 재현성) -------------------------------------------
    seed_main <- if (!is.null(seed_override))
      as.integer(seed_override)
    else
      as.integer(seed)

    # 3) 샘플링 베이스 인자 --------------------------------------------------------
    base_args <- list(
      data = stan_list,
      seed = seed_main,
      chains = chains,
      parallel_chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      metric = metric,
      init = init,
      refresh = if (quiet)
        0
      else
        100
    )

    # 4) 공통 래퍼로 샘플 + 리트라이 ----------------------------------------------
    tag_lbl  <- sprintf("%s\u2192%s main", partner, target)
    pair_tag <- sprintf("%s\u2192%s", partner, target)

    # --- main run start ---
    if (progress_local != "none")
      cat(sprintf("%s pair: main run started\n", pair_tag))

    res_try <- .sample_with_retry(
      mod = mod,
      base_args = base_args,
      stan_list = stan_list,
      max_retries = max_retries,
      nu_fix = nu_fix,
      silent_sampler = silent_sampler,
      ebfmi_thresh = 0.30,
      max_nu_cap = 7,
      tag = tag_lbl
    )

    # --- main run done ---
    if (progress_local != "none")
      cat(sprintf("%s pair: main run completed\n", pair_tag))

    fit_full   <- res_try$fit
    diag       <- res_try$diag
    n_retries  <- res_try$n_retries
    fit_failed <- res_try$fit_failed
    used_nu    <- res_try$used_nu
    final_args_main <- res_try$final_args
    nu_final        <- used_nu

    # 5) 드로우 요약 ---------------------------------------------------------------
    d <- safe_draws_df(fit_full)  # a_ij, a_ii, r0, tau_r, sigma, sigma_ou, lambda, ...
    if (isTRUE(use_student_t) && !("nu" %in% names(d))) {
      d$nu <- rep(used_nu, nrow(d))
    }

    # 방어적 체크
    if (!all(c("a_ij", "a_ii") %in% names(d))) {
      stop("Expected parameters a_ij and a_ii not found in draws.")
    }

    aij <- d$a_ij
    aii <- d$a_ii
    p_sign2_dir  <- .clip01(2 * pmin(mean(aij > 0), mean(aij < 0)))
    p_sign2_self <- .clip01(2 * pmin(mean(aii > 0), mean(aii < 0)))

    # ----- 6) Repeated K-fold with outer retries ----------------------------------
    kfold <- NULL

    # ▼ k-fold 메타 기본값(ELPD 계산 안 할 때도 반환값 보장)
    kfold_rounds_attempted <- 0L
    kfold_failed_overall   <- NA
    kfold_n_ok   <- NA_integer_
    kfold_n_fail <- NA_integer_

    if (isTRUE(compute_elpd)) {
      seed_for_kfold <- if (is.null(final_args_main$seed))
        seed_main
      else
        final_args_main$seed
      if (progress_local != "none")
        cat(
          sprintf(
            "%s pair: k-fold evaluation started (K=%d, R=%d)\n",
            pair_tag,
            kfold_K,
            kfold_R
          )
        )
      sargs_kfold <- final_args_main
      sargs_kfold$step_size   <- NULL   # 초기 step size 강제 금지 (폴드별 adapt)
      sargs_kfold$inv_metric  <- NULL   # metric 파일/행렬 강제 금지
      sargs_kfold$metric_file <- NULL


      # 본계산에서 사용된 ν로 K-fold 전 폴드 **고정** 실행 (bump 금지) ---
      nu_kfold_main <- nu_final

      kfold <- .repkfold_eval(
        mod = mod,
        stan_list_base   = stan_list,
        sample_args_base = sargs_kfold,
        pair_in = pair_in,
        K = kfold_K,
        R = kfold_R,
        seed = seed_for_kfold,
        resid_mode = resid_mode,
        use_t = isTRUE(use_student_t),
        elpd_mode = elpd_mode,
        silent_sampler = silent_sampler,
        n_workers_kfold = n_workers_kfold_eff,
        nu_fixed_kfold = nu_kfold_main,          # ← ν 고정
        freeze_retry_hypers = TRUE
      )
      if (progress_local != "none")
        cat(sprintf("%s pair: k-fold evaluation completed\n", pair_tag))
      kfold_n_ok   <- if (is.null(kfold))
        NA_integer_
      else
        kfold$n_folds_ok
      kfold_n_fail <- if (is.null(kfold))
        NA_integer_
      else
        kfold$n_folds_fail
      # 폴드 훈련-진단 실패가 하나라도 있으면 ν=7로 K-fold **전부 재실행** ---
      # (정보 누출 아님: 테스트 ELPD/관측은 보지 않고, 훈련 진단만 사용)
      kfold_outer_rounds_local <- 1L
      bad_count <- 0L
      if (!is.null(kfold) && is.data.frame(kfold$fold_diag) && nrow(kfold$fold_diag)) {
        fd <- kfold$fold_diag
        bad_flag <- with(fd,
                         (is.finite(n_divergent)     & n_divergent     > 0) |
                           (is.finite(ebfmi_min)       & ebfmi_min       < 0.30) |
                           (is.finite(treedepth_hits)  & treedepth_hits  > 0) |
                           (is.finite(worst_rhat)      & worst_rhat      >= 1.05) |
                           (is.finite(min_ess_bulk)    & min_ess_bulk    < 400)
        )
        bad_count <- sum(bad_flag, na.rm = TRUE)
      }
      if (isTRUE(bad_count > 0) && is.finite(nu_kfold_main) && nu_kfold_main < 7) {
        if (progress_local != "none")
          cat(sprintf("%s pair: k-fold retrial triggered by TRAIN diagnostics → nu=7 (global)\n", pair_tag))
        kfold <- .repkfold_eval(
          mod = mod,
          stan_list_base   = stan_list,
          sample_args_base = sargs_kfold,
          pair_in = pair_in,
          K = kfold_K,
          R = kfold_R,
          seed = seed_for_kfold,
          resid_mode = resid_mode,
          use_t = isTRUE(use_student_t),
          elpd_mode = elpd_mode,
          silent_sampler = silent_sampler,
          n_workers_kfold = n_workers_kfold_eff,
          nu_fixed_kfold = 7,                   # ← 글로벌 한 단계 상승(최대 1회)
          freeze_retry_hypers = TRUE
        )
        kfold_outer_rounds_local <- 2L
        if (progress_local != "none")
          cat(sprintf("%s pair: k-fold evaluation (nu=7) completed\n", pair_tag))
        kfold_n_ok   <- if (is.null(kfold)) NA_integer_ else kfold$n_folds_ok
        kfold_n_fail <- if (is.null(kfold)) NA_integer_ else kfold$n_folds_fail


        # --- 본계산도 ν=7로 일관화: 테스트 점수는 쓰지 않고, K-fold '훈련 진단' 트리거만으로 재적합 ---
        # 누출이 아닌 이유:
        #  (1) ν=7 결정은 오직 K-fold '훈련' 진단(다이버전스/E-BFMI/treedepth/R-hat/ESS)으로만 내림
        #  (2) 테스트 ELPD/관측치는 보지 않으며, 선택/재적합 판단에 사용되지 않음
        #  (3) 최종 보고는 사양 확정 후 단일 CV 및 본계산 결과만 사용
        if (is.finite(nu_final) && nu_final < 7) {
          if (progress_local != "none")
            cat(sprintf("%s pair: main refit with nu=7 for spec consistency\n", pair_tag))
          base_args_refit <- final_args_main
          # 리핏 시 bump/하이퍼 변경 금지, seed/체인 동일, metric 그대로
          base_args_refit$data$nu_fixed <- 7
          base_args_refit$adapt_delta   <- final_args_main$adapt_delta
          base_args_refit$max_treedepth <- final_args_main$max_treedepth
          base_args_refit$metric        <- final_args_main$metric
          # 리핏은 한 번만, bump 금지
          refit_try <- .sample_with_retry(
            mod = mod,
            base_args = base_args_refit,
            stan_list = base_args_refit$data,
            max_retries = 0,
            nu_fix = 7,
            max_nu_cap = 7,
            silent_sampler = silent_sampler,
            freeze_retry_hypers = TRUE,
            tag = sprintf("%s\u2192%s main-refit-nu7", partner, target)
          )
          fit_full <- refit_try$fit
          diag     <- refit_try$diag
          nu_final <- 7
          d <- safe_draws_df(fit_full)
          aij <- d$a_ij; aii <- d$a_ii
          p_sign2_dir  <- .clip01(2 * pmin(mean(aij > 0), mean(aij < 0)))
          p_sign2_self <- .clip01(2 * pmin(mean(aii > 0), mean(aii < 0)))
        }
      }
    }




    # 7) 리턴 ----------------------------------------------------------------------
    list(
      n_pairs = nrow(pair_in),

      a_mean  = mean(aij),
      a_sd = stats::sd(aij),
      a_q2.5  = stats::quantile(aij, 0.025),
      a_q97.5 = stats::quantile(aij, 0.975),
      p_sign2 = p_sign2_dir,

      aii_mean = mean(aii),
      aii_sd = stats::sd(aii),
      aii_q2.5 = stats::quantile(aii, 0.025),
      aii_q97.5 = stats::quantile(aii, 0.975),
      p_sign2_self = p_sign2_self,

      # 원시 진단치/리트라이 메타(후처리 summariser에서 필터 예정)
      diag = diag,
      n_retries = n_retries,
      fit_failed = fit_failed,
      nu_fixed_used = nu_final,

      # 스무딩 메타
      smoothed = sflag,
      smooth_scale = ssc,
      smooth_edf_mean = sedf,
      subjects = unique(pair_in$subject),

      # Repeated K-fold 결과 요약 (+ pseudo-BMA/stacking용 payload)
      kfold = kfold,
      kfold_mean = if (is.null(kfold)) NA_real_ else kfold$elpd_mean,

      # 용어 명확화: 바깥 재평가 라운드(outer)만 카운트
      kfold_outer_rounds = if (isTRUE(compute_elpd)) kfold_outer_rounds_local else 0L,
      kfold_failed = is.null(kfold) ||
        (!is.null(kfold$n_folds_fail) && kfold$n_folds_fail > 0L),
      kfold_n_folds_ok   = if (is.null(kfold))
        NA_integer_
      else
        kfold$n_folds_ok,
      kfold_n_folds_fail = if (is.null(kfold))
        NA_integer_
      else
        kfold$n_folds_fail,
      # subject-level OOF elpd vector (named) for model averaging/stacking
      kfold_subject       = list(if (is.null(kfold))
        NULL
        else
          kfold$elpd_subject),
      kfold_subject_ids   = list(if (is.null(kfold))
        NULL
        else
          names(kfold$elpd_subject)),
      kfold_subject_counts = list(if (is.null(kfold))
        NULL
        else
          kfold$subject_test_counts),
      # split manifest (r,k,test_subjects) & seed to ensure same partitions across models
      kfold_splits        = list(if (is.null(kfold))
        NULL
        else
          kfold$splits_df),
      kfold_seed_used     = if (is.null(kfold))
        NA_integer_
      else
        seed_for_kfold,
      kfold_K             = if (is.null(kfold))
        NA_integer_
      else
        kfold$K,
      kfold_R             = if (is.null(kfold))
        NA_integer_
      else
        kfold$R,
      kfold_agg           = "subject-uniform",
      # convenient dispersion summaries
      kfold_sd            = if (is.null(kfold))
        NA_real_
      else
        stats::sd(kfold$elpd_subject, na.rm = TRUE),
      kfold_se            = if (is.null(kfold))
        NA_real_
      else {
        nn <- sum(is.finite(kfold$elpd_subject))
        stats::sd(kfold$elpd_subject, na.rm = TRUE) / sqrt(pmax(nn, 1L))
      },
      kfold_n_subjects    = if (is.null(kfold))
        NA_integer_
      else
        sum(is.finite(kfold$elpd_subject)),

      # 폴드-내 리트라이/nu 사용 요약
      kfold_retry_total  = if (is.null(kfold) ||
                               is.null(kfold$retry_total))
        NA_integer_
      else
        kfold$retry_total,
      kfold_retry_mean   = if (is.null(kfold) ||
                               is.null(kfold$retry_mean))
        NA_real_
      else
        kfold$retry_mean,
      kfold_nu_used_counts = list(if (is.null(kfold) ||
                                      is.null(kfold$nu_used_counts))
        NULL
        else
          kfold$nu_used_counts)
    )
  }

  # outer loop (both directions per pair) --------------------------------
  out <- list()
  k <- 0L
  n_taxa  <- length(taxa_vec)
  n_tasks <- 2L * base::choose(n_taxa, 2L)
  task_k  <- 0L

  if (n_workers_outer <= 1L) {
    # ====== 순차 모드: 기존 진행바 유지 ======
    if (progress == "bar") {
      pb <- utils::txtProgressBar(min = 0,
                                  max = n_tasks,
                                  style = 3)
      on.exit(try(close(pb), silent = TRUE)
              , add = TRUE)
    }
    bump <- function() {
      task_k <<- task_k + 1L
      if (progress == "bar" && task_k %% progress_every == 0L) {
        utils::setTxtProgressBar(pb, task_k)
      }
    }
    log_once <- function(tag, res, i_nm, j_nm) {
      if (progress != "verbose")
        return(invisible())
      if (is.null(res)) {
        cat(sprintf(
          "[%d/%d] %s: %s → %s | skipped\n",
          task_k,
          n_tasks,
          tag,
          j_nm,
          i_nm
        ))
        return(invisible())
      }
      d <- res$diag
      cat(
        sprintf(
          "[%d/%d] %s: %s → %s | div=%s ebfmi_min=%s rhat=%s ess_bulk=%s\n",
          task_k,
          n_tasks,
          tag,
          j_nm,
          i_nm,
          ifelse(is.finite(d$n_divergent), d$n_divergent, NA),
          ifelse(
            is.finite(d$ebfmi_min),
            sprintf("%.3f", d$ebfmi_min),
            NA
          ),
          ifelse(
            is.finite(d$worst_rhat),
            sprintf("%.3f", d$worst_rhat),
            NA
          ),
          ifelse(is.finite(d$min_ess_bulk), d$min_ess_bulk, NA)
        )
      )
    }

    for (i in 1:(n_taxa - 1L)) {
      for (j in (i + 1L):n_taxa) {
        row <- .run_pair(i, j, mute_logs = FALSE, seed_base = seed)
        # run_one 내부에서 두 번 호출되므로 2 step 증가
        bump()
        bump()
        if (!is.null(row)) {
          k <- k + 1L
          out[[k]] <- row
          if (progress == "verbose") {
            cat(sprintf("pair %s-%s done\n", taxa_vec[i], taxa_vec[j]))
          }
        } else {
          if (progress == "verbose") {
            cat(sprintf("pair %s-%s skipped\n", taxa_vec[i], taxa_vec[j]))
          }
        }
      }
    }

  } else {
    # ====== 병렬 모드: furrr로 쌍 병렬 ======
    pair_idx <- utils::combn(n_taxa, 2, simplify = FALSE)

    # furrr plan 설정 (PSOCK 권장)
    op <- future::plan()
    on.exit(future::plan(op), add = TRUE)
    future::plan(future::multisession, workers = n_workers_outer)

    # 병렬에서는 내부 로그를 끄고, progressr가 있으면 진행바 표시
    if (has_progressr && identical(progress, "bar")) {
      progressr::with_progress({
        p <- progressr::progressor(steps = length(pair_idx))
        out <- furrr::future_map(pair_idx, function(idx) {
          res <- .run_pair(idx[1],
                           idx[2],
                           mute_logs = TRUE,
                           seed_base = seed)
          p(message = sprintf("pair %s-%s", taxa_vec[idx[1]], taxa_vec[idx[2]]))
          res
        }, .options = furrr::furrr_options(
          seed = TRUE,
          globals = list(
            .run_pair = .run_pair,
            # 워커에 강제 전달
            run_one   = run_one,
            # ← run_one 누락 방지 (핵심)
            taxa_vec  = taxa_vec,
            kfold_K   = kfold_K,
            kfold_R   = kfold_R,
            progress  = progress,
            seed      = seed
          )
        ))
      })
    } else {
      out <- furrr::future_map(pair_idx, function(idx)
        .run_pair(idx[1], idx[2], mute_logs = TRUE, seed_base = seed), .options = furrr::furrr_options(
          seed = TRUE,
          globals = list(
            .run_pair = .run_pair,
            run_one   = run_one,
            taxa_vec  = taxa_vec,
            kfold_K   = kfold_K,
            kfold_R   = kfold_R,
            progress  = progress,
            seed      = seed
          )
        ))
    }
  }

  res <- dplyr::bind_rows(purrr::compact(out))
  rownames(res) <- NULL
  return(res)
}
