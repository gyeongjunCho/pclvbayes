##### R/fit_pclv_bayes.R
#' Pairwise Bayesian **pcLV** Core fitting with irregular-time OU residuals
#'
#' @description
#' This function belongs to the **pGLVbayes** framework, which extends the
#' classical generalized Lotka–Volterra (gLV) model to a **Bayesian,
#' pairwise compositional (pcLV)** formulation suitable for microbiome
#' time-series data.
#'
#' Conceptually, the gLV equations are reformulated on the additive
#' log-ratio (ALR) scale to handle compositional constraints, and are then
#' estimated in a **pairwise Bayesian regression** structure with irregular-time
#' OU residuals. Each unordered taxon pair `{i, j}` is fitted in
#' both directions (j → i and i → j), providing directional interaction
#' coefficients that can be aggregated across subjects.
#'
#' Residual dependence follows an irregular-time OU process. Student-t observation
#' noise and repeated subject-level K-fold Kalman ELPD scoring are fixed Core choices.
#'
#' In short, **pGLVbayes** generalizes the traditional gLV by introducing
#' (i) compositional transformation, (ii) pairwise modularization, and
#' (iii) fully Bayesian estimation with diagnostic-aware sampling.
#'
#' @details
#' The canonical Core uses pair-to-rest ALR predictors at lag 1 and the response
#' `y = Delta ALR_i / Delta t`. Predictors are globally standardized using full
#' analysis data for the main fit and training-only statistics within each fold;
#' the response is not standardized.
#'
#' Repeated subject-level K-fold ELPD is always computed. Held-out observations
#' are scored with the irregular-time OU Kalman path using a Student-t observation
#' likelihood with one posterior-estimated `nu > 2` per directed fit. Spline smoothing is selected by cross-validation.
#' Failed folds remain explicit and contribute no zero-valued ELPD placeholders.
#'
#' - A diagnostics-aware **retry** policy (up to `max_retries`) monitors divergences,
#'   tree depth, and E-BFMI without changing the Student-t model or its `nu` prior.
#'   Diagnostics are returned for downstream filtering.
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
#' @param nz_partner_min_frac Drop subjects whose partner predictor is non-zero in fewer than this
#' fraction of rows; default \code{0.15}.
#'
#' @param zero_mode_alr Mode for zero-aware replacement in the ALR triplet; one of
#' \code{"minpos_time"}, \code{"minpos_subject"}, \code{"lib"}, \code{"fixed"}.
#' @param minpos_alpha Multiplier for minimal positive pseudo-count when applicable; default \code{0.5}.
#' @param minpos_base Basis for minimal positive detection, \code{"ij"} or \code{"triplet"}.
#' @param eps_fixed Fixed pseudo-count floor; default \code{1e-6}.
#' @param lib_eps_c Library-size scaled pseudo-count constant; default \code{0.65}.
#' @param rest_floor_frac Floor fraction for the rest component; default \code{1.0}.
#' @param smooth_scale Smoothing scale for optional per-subject spline on predictors;
#' one of \code{"logra"}, \code{"alr"}.
#' @param alr_spline_df,alr_spline_spar,alr_spline_cv Controls for spline smoothing if used.
#'
#' @param silent_sampler If \code{TRUE}, silence cmdstanr progress/messages in the wrapper; default \code{FALSE}.
#' @param progress Progress output: \code{"bar"}, \code{"verbose"}, or \code{"none"}.
#' When \code{"bar"} and \pkg{progressr} is available, outer/K-fold progress bars are displayed.
#' @param progress_every Update frequency for the legacy sequential progress bar; default \code{1}.
#'
#'
#' @param kfold_K Number of folds K; default \code{5}.
#' @param kfold_R Number of repetitions R; default \code{3}.
#' @param kfold_seed Seed for K-fold split reproducibility (defaults to `seed`).
#' @param n_workers_kfold Number of parallel workers for K-fold (PSOCK via future/furrr); default \code{1}.
#' @param n_workers_outer Number of parallel workers for the outer pair loop; default \code{1}.
#' If \code{> 1}, K-fold parallelism is automatically disabled to avoid nested parallelism.
#'
#' @param max_retries Maximum sampler retry attempts (diagnostics-aware); default \code{3}.
#'
#' @return
#' A list containing directional and self-effect summaries, pointwise repeated
#' K-fold ELPD, and the raw bidirectional pair table. ELPD uses the
#' `kalman-ou` method and includes fold-evidence counts.
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
  nz_partner_min_frac = 0.15,

  # --- 빌더 단계(로그-RA) 스무딩 & 최소 요구량 ---
  eps = 1e-6,
  min_unique_times = 3,
  min_pairs = 4,

  # --- ALR zero 처리 & (선택) ALR 스무딩 ---
  zero_mode_alr = c("minpos_time", "minpos_subject", "lib", "fixed"),
  minpos_alpha = 0.5,
  minpos_base = c("ij", "triplet"),
  eps_fixed = 1e-6,
  lib_eps_c = 0.65,
  rest_floor_frac = 1.0,
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
  kfold_K = 5,
  kfold_R = 3,
  kfold_seed = seed,
  n_workers_kfold = 1L,

  # --- 리트라이 & t-꼬리 제어 ---
  max_retries = 3,

  # --- Pathfinder (옵션) ---
  use_pathfinder_init = TRUE,
  pf_num_paths = 8,
  pf_draws = 1000,
  pf_history_size = 50,
  pf_max_lbfgs_iters = 200,
  pf_psis_resample = TRUE
)
{
  public_control_names <- setdiff(names(formals(sys.function())),
                                  c("physeq", "subject_col", "time_col", "taxa_vec"))
  validated <- .validate_fit_pclv_inputs(
    physeq, subject_col, time_col, taxa_vec,
    controls = mget(public_control_names, envir = environment(), inherits = FALSE)
  )
  for (nm in names(validated$controls))
    assign(nm, validated$controls[[nm]], envir = environment())
  subject_col <- validated$subject_col
  time_col <- validated$time_col
  taxa_vec <- validated$taxa_vec
  meta_df <- validated$meta_df
  mat_rel <- validated$mat_rel

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

  # Stan model (canonical Student-t likelihood)
  mod <- tryCatch(
    get_pclv_model(quiet = quiet),
    error = function(e) .pclv_failure("model_loading", "canonical_model_compilation_failed",
                                      list(message = conditionMessage(e)))
  )
  if (.is_pclv_failure(mod)) return(mod)

  # legacy smoothing on log(RA+eps)
  sm_mat <- .precompute_spline_smoothed(
    mat_rel,
    meta_df,
    taxa_vec,
    eps,
    min_unique_times
  )
  if (inherits(sm_mat, "pclv_failure")) return(sm_mat)

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

  # Canonical unordered-pair tasks; scheduling is the only mode-specific layer.
  tasks <- .make_pair_tasks(taxa_vec, seed)
  scheduling <- list(n_workers_kfold_eff = n_workers_kfold_eff)
  run_task <- function(task, mute_logs) {
    .execute_pair_task(
      task = task, taxa_vec = taxa_vec, run_one = .run_one,
      core_ctx = ctx, scheduling = scheduling,
      progress = progress, mute_logs = mute_logs
    )
  }

  if (n_workers_outer <= 1L) {
    pb <- NULL
    if (progress == "bar") {
      pb <- utils::txtProgressBar(min = 0, max = length(tasks), style = 3)
      on.exit(try(close(pb), silent = TRUE), add = TRUE)
    }
    out <- vector("list", length(tasks))
    for (task_index in seq_along(tasks)) {
      out[[task_index]] <- run_task(tasks[[task_index]], mute_logs = FALSE)
      if (!is.null(pb) && task_index %% progress_every == 0L)
        utils::setTxtProgressBar(pb, task_index)
      if (progress == "verbose")
        cat(sprintf("pair %s-%s done\n", tasks[[task_index]]$taxon_i, tasks[[task_index]]$taxon_j))
    }
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_workers_outer)
    if (has_progressr && identical(progress, "bar")) {
      progressr::with_progress({
        progressor <- progressr::progressor(steps = length(tasks))
        out <- furrr::future_map(
          tasks,
          function(task) {
            result <- run_task(task, mute_logs = TRUE)
            progressor(message = sprintf("pair %s-%s", task$taxon_i, task$taxon_j))
            result
          },
          .options = furrr::furrr_options(seed = TRUE, globals = TRUE)
        )
      })
    } else {
      out <- furrr::future_map(
        tasks, function(task) run_task(task, mute_logs = TRUE),
        .options = furrr::furrr_options(seed = TRUE, globals = TRUE)
      )
    }
  }

  res <- .assemble_pair_outcomes(out, tasks)
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
