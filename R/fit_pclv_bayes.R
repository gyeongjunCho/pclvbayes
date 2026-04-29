##### R/fit_pclv_bayes.R

#' Pairwise Bayesian **pcLV** fitting (bidirectional) with OU/white-noise residuals
#'
#' @description
#' This function belongs to the **pGLVbayes** framework, which extends the
#' classical generalized Lotka–Volterra (gLV) model to a **Bayesian,
#' pairwise compositional (pcLV)** formulation suitable for microbiome
#' time-series data.
#'
#' Conceptually, the gLV equations are reformulated on the additive
#' log-ratio (ALR) scale to handle compositional constraints, and are then
#' estimated in a **pairwise Bayesian regression** structure with OU or
#' white-noise residuals. Each unordered taxon pair `{i, j}` is fitted in
#' both directions (j → i and i → j), providing directional interaction
#' coefficients that can be aggregated across subjects.
#'
#' Residual processes may follow an OU (Ornstein–Uhlenbeck, ≈ AR(1))
#' correlation structure or be modeled as independent Gaussian noise.
#' Optional Student-t noise and repeated K×R cross-validation with Kalman
#' ELPD scoring provide robust inference and model comparison.
#'
#' In short, **pGLVbayes** generalizes the traditional gLV by introducing
#' (i) compositional transformation, (ii) pairwise modularization, and
#' (iii) fully Bayesian estimation with diagnostic-aware sampling.
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
#'   and E-BFMI; if needed, it **escalates `nu_fix` within 4 ~ 7** while keeping other sampler controls
#'   stable for fold runs. Diagnostics are returned for downstream filtering.
#' - Optional **Pathfinder** initialization (`use_pathfinder_init = TRUE`) uses
#'   `cmdstanr::pathfinder()` to obtain near-posterior inits (and mass-matrix info) before HMC.
#'   When enabled and `iter_warmup >= 1500`, the warmup is **auto-shortened** to
#'   `max(500, floor(iter_warmup/2))` to exploit the better starting point. You can fine-tune
#'   Pathfinder via `pf_num_paths`, `pf_draws`, `pf_history_size`, `pf_max_lbfgs_iters`,
#'   and `pf_psis_resample`.
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
#' @param physeq A \code{phyloseq} object with taxa in rows (will be transposed if needed).
#' @param taxa_vec Optional character vector of taxa IDs to subset; defaults to all taxa.
#' @param eps Small constant for legacy log-RA smoothing at the builder stage; default \code{1e-6}.
#' @param min_pairs Minimum number of valid time-adjacent pairs per direction; default \code{4}.
#'
#' @param subject_col Column name in sample metadata indicating subjects (e.g., plant/plot).
#' @param time_col    Column name in sample metadata indicating (numeric) time within subject.
#'
#' @param chains,iter_warmup,iter_sampling,seed,quiet,adapt_delta,max_treedepth,metric,init
#' Sampling controls passed to \pkg{cmdstanr}.
#' @param use_pathfinder_init Logical; whether to run \code{cmdstanr::pathfinder()} to
#' obtain Pathfinder-based initial draws before HMC/NUTS. **Default = TRUE (recommended)**.
#' When enabled and \code{iter_warmup >= 1500}, warmup is automatically reduced to
#' \code{max(500, floor(iter_warmup/2))}. This typically preserves or improves convergence
#' while reducing runtime.
#' @param pf_num_paths Integer; number of quasi-Newton paths for Pathfinder (default e.g. \code{8}).
#' @param pf_draws Integer; number of draws to sample from the Pathfinder approximation (default e.g. \code{2000}).
#' @param pf_history_size Integer; limited-memory BFGS history size (default e.g. \code{10}).
#' @param pf_max_lbfgs_iters Integer; maximum L-BFGS iterations per path (default e.g. \code{1000}).
#' @param pf_psis_resample Logical; if \code{TRUE}, use PSIS-importance resampling inside Pathfinder.
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
#' @section Progress & Parallel:
#' - Set global handlers once for pretty bars: \preformatted{
#'   if (requireNamespace("progressr", quietly = TRUE)) {
#'     progressr::handlers(global = TRUE); progressr::handlers("cli")
#'   }}
#' - Reproducibility: seeds are deterministic per pair/fold (`seed`, `kfold_seed`) and
#'   \code{furrr::future_map(..., .options = furrr::furrr_options(seed = TRUE))} is used.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("progressr", quietly = TRUE)) {
#'   progressr::handlers(global = TRUE); progressr::handlers("cli")
#' }
#' EP_res <- fit_pclv_bayes(
#'   physeq = EP_phy_obj, subject_col = "plot2", time_col = "week",
#'   transform = "alr", resid_mode = "ou",  # elpd_mode auto -> "kalman"
#'   use_student_t = TRUE,
#'   use_pathfinder_init = TRUE,
#'   chains = 4, iter_warmup = 1000, iter_sampling = 1500,
#'   n_workers_outer = 4,   # outer pair loop in parallel
#'   n_workers_kfold = 1,   # k-fold parallel disabled when outer > 1
#'   progress = "bar"
#' )
#' }
#'
#' @name fit_pclv_bayes
#' @rdname fit_pclv_bayes
#' @export

fit_pclv_bayes <- function(# --- 필수 입력 ---
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
  init = 0.2, # pathfinder가 있으면 무시됨
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
  nu_fix = 5,

  # --- Pathfinder (옵션) ---
  use_pathfinder_init = TRUE,
  pf_num_paths = 8,
  pf_draws = 1000,
  pf_history_size = 50,
  pf_max_lbfgs_iters = 200,
  pf_psis_resample = TRUE
)
{

  try({
    .cleanup_stale_csv_start <- function(older_than_hours = 24) {
      cutoff <- Sys.time() - older_than_hours * 3600
      roots <- unique(c(
        getOption("glvpair.output_root",
                  tools::R_user_dir("glvpair", which = "cache")),
        tempdir()
      ))
      for (rt in roots) {
        if (!dir.exists(rt)) next
        files <- list.files(rt, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
        if (!length(files)) next
        # glvpair/cmdstan 관련 경로만, csv/txt/json만 대상으로 제한
        sel <- grepl("(glvpair|cmdstan)", files) &
          grepl("\\.(csv|txt|json)$", files, ignore.case = TRUE)
        if (!any(sel)) next
        info <- suppressWarnings(file.info(files[sel]))
        stale <- rownames(info)[is.finite(info$mtime) & info$mtime < cutoff]
        if (length(stale)) unlink(stale, recursive = TRUE)
      }
    }
    .cleanup_stale_csv_start(older_than_hours = 24)
  }, silent = TRUE)


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
  mod <- get_pclv_model(quiet = quiet)

  # metadata & matrix
  meta_df <- .get_sample_meta(physeq, subject_col, time_col)
  mat_rel <- .get_abund_matrix_precomputed(physeq)
  if (is.null(taxa_vec))
    taxa_vec <- phyloseq::taxa_names(physeq)
  taxa_vec <- intersect(taxa_vec, rownames(mat_rel))
  if (length(taxa_vec) < 2)
    stop("taxa_vec must contain at least 2 taxa.")

  # legacy smoothing on log(RA+eps)
  sm_mat <- .precompute_spline_smoothed(
    mat_rel,
    meta_df,
    taxa_vec,
    eps,
    spline_df,
    spline_spar,
    spline_cv,
    min_unique_times
  )

  # --------------------------------------------------------
  # --- Avoid nested parallelism: if outer-parallel, disable kfold-parallel
  n_workers_kfold_eff <- if (n_workers_outer > 1L)
    1L
  else
    n_workers_kfold

  # progressr availability (optional)
  has_progressr <- requireNamespace("progressr", quietly = TRUE)

  # ---------- 런타임 컨텍스트를 한데 모아 병렬 워커로 전달 ----------
  ctx <- list(
    mod_exe_file = mod$exe_file(),
    # 데이터
    meta_df = meta_df,
    sm_mat = sm_mat,
    eps = eps,
    min_pairs = min_pairs,
    min_unique_times = min_unique_times,
    # 입력/전처리 설정
    compute_elpd = compute_elpd,
    standardize_by_subject = standardize_by_subject,
    transform = transform,
    lag = lag,
    zero_mode_alr = zero_mode_alr,
    minpos_alpha = minpos_alpha,
    minpos_base  = minpos_base,
    eps_fixed = eps_fixed,
    lib_eps_c = lib_eps_c,
    rest_floor_frac = rest_floor_frac,
    smooth_scale = smooth_scale,
    alr_spline_df = alr_spline_df,
    alr_spline_spar = alr_spline_spar,
    alr_spline_cv = alr_spline_cv,
    nz_partner_min_frac = nz_partner_min_frac,
    # 모델/샘플러 설정
    max_retries = max_retries,
    nu_fix = nu_fix,
    use_student_t = use_student_t,
    elpd_mode = elpd_mode,
    resid_mode = resid_mode,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    metric = metric,
    init = init,
    seed = seed,
    quiet = quiet,
    silent_sampler = silent_sampler,
    # K-fold/병렬
    n_workers_kfold_eff = n_workers_kfold_eff,
    kfold_K = kfold_K,
    kfold_R = kfold_R,

    # PF 옵션 전달
    use_pathfinder_init = use_pathfinder_init,
    pf_num_paths = pf_num_paths,
    pf_draws = pf_draws,
    pf_history_size = pf_history_size,
    pf_max_lbfgs_iters = pf_max_lbfgs_iters,
    pf_psis_resample = pf_psis_resample
  )


  if (has_progressr && identical(progress, "bar")) {
    old_opt <- options(
      progressr.enable   = TRUE,  # with_progress 안에서 진행바 활성
      progressr.clear    = FALSE, # 완료 후 바를 지우지 않음
      progressr.handlers = list(
        progressr::handler_progress(
          format = "|:bar| :percent :current of :total",
          clear  = FALSE,
          width  = 40
        )
      )
    )
    on.exit(options(old_opt), add = TRUE)  # 함수 종료 시 원복
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

    say_every_pairs <- 10L

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
        row <- .run_pair(
          i, j,
          taxa_vec = taxa_vec,
          .run_one  = .run_one,
          ctx = ctx,
          progress = progress,
          mute_logs = FALSE, seed_base = seed,
        )
        # .run_one 내부에서 두 번 호출되므로 2 step 증가
        bump()
        bump()

        # 진행바 한 줄 위에 누적 메시지(요약) 출력 후, 바를 다시 그려서 레이아웃 복구
        if (progress == "bar" && (task_k %% (2L * say_every_pairs) == 0L)) {
          pct <- 100 * task_k / n_tasks
          cat(sprintf(
            "\ncompleted %d / %d tasks (%.1f%%) — last pair %s-%s\n",
            task_k / 2L, n_tasks / 2L, pct, taxa_vec[i], taxa_vec[j]
          ))
          utils::setTxtProgressBar(pb, task_k)
        }

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
        out <- furrr::future_map(
          pair_idx,
          function(idx) {
            res <- .run_pair(
              idx[1], idx[2],
              taxa_vec = taxa_vec,
              .run_one  = .run_one,
              ctx = ctx,
              kfold_K  = kfold_K,
              kfold_R  = kfold_R,
              progress = progress,
              mute_logs = TRUE, seed_base = seed
            )
            # 개별 페어 진행 메시지
            p(message = sprintf("pair %s-%s", taxa_vec[idx[1]], taxa_vec[idx[2]]),
              amount = 1)
            res
          },
          .options = furrr::furrr_options(
            seed = TRUE,
            globals = TRUE
          )
        )
      })
    } else {
      out <- furrr::future_map(
        pair_idx,
        function(idx) .run_pair(
          idx[1], idx[2],
          ctx = ctx,
          taxa_vec = taxa_vec,
          .run_one  = .run_one,
          progress = progress,
          mute_logs = TRUE, seed_base = seed
        ),
        .options = furrr::furrr_options(
          seed = TRUE,
          globals = TRUE
        )
      )
    }
  }

  res <- dplyr::bind_rows(purrr::compact(out))
  rownames(res) <- NULL

  cross_tbl <- .mk_cross(res)
  self_tbl <- .mk_self(res)
  # ---------- ③ elpd_pointwise_cross ----------
  elpd_pointwise_cross_tbl <- .expand_cross_pw(res)
  elpd_pointwise_self_tbl <- .expand_self_pw(res)
  # ---------- ⑤ raw: 디버그용 와이드 테이블 ----------
  raw_tbl <- res

  if ("fit" %in% names(raw_tbl)) {
    try({
      ff <- try(raw_tbl$fit$output_files(), silent = TRUE)
      if (!inherits(ff, "try-error") && length(ff) > 0) {
        invisible(lapply(ff, function(f) if (file.exists(f)) unlink(f)))
      }
    }, silent = TRUE)
  }

  # 최종 리스트로 반환
  out_obj <- list(
    cross = cross_tbl,
    self  = self_tbl,
    elpd_pointwise_cross = elpd_pointwise_cross_tbl,
    elpd_pointwise_self  = elpd_pointwise_self_tbl,
    raw = raw_tbl
  )
  return(out_obj)
}
