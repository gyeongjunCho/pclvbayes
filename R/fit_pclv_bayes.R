.prepare_fit_runtime <- function(validated) {
  controls <- validated$controls
  mod <- tryCatch(
    get_pclv_model(quiet = controls$quiet),
    error = function(e) .pclv_failure(
      "model_loading", "canonical_model_compilation_failed",
      list(message = conditionMessage(e))
    )
  )
  if (.is_pclv_failure(mod)) return(mod)
  sm_mat <- .precompute_spline_smoothed(
    validated$mat_rel, validated$meta_df, validated$taxa_vec,
    controls$eps, controls$min_unique_times
  )
  if (.is_pclv_failure(sm_mat)) return(sm_mat)
  context_names <- c(
    "eps", "min_pairs", "min_unique_times", "zero_mode_alr",
    "minpos_alpha", "minpos_base", "eps_fixed", "lib_eps_c",
    "rest_floor_frac", "smooth_scale", "alr_spline_df",
    "alr_spline_spar", "alr_spline_cv", "nz_partner_min_frac",
    "max_retries", "chains", "iter_warmup", "iter_sampling",
    "adapt_delta", "max_treedepth", "metric", "init", "seed",
    "quiet", "silent_sampler", "kfold_K", "kfold_R", "kfold_seed",
    "use_pathfinder_init", "pf_num_paths", "pf_draws",
    "pf_history_size", "pf_max_lbfgs_iters", "pf_psis_resample"
  )
  list(
    ctx = c(list(
      mod_exe_file = mod$exe_file(),
      meta_df = validated$meta_df,
      sm_mat = sm_mat
    ), controls[context_names]),
    model = mod
  )
}

##### R/fit_pclv_bayes.R
#' Bayesian directed pair-to-rest fitting with irregular-time OU residuals
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
#' both directions (j → i and i → j), providing directed pair-to-rest dynamic coefficients that can be aggregated
#' across subjects; these are not generally absolute direct-gLV coefficients.
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
#' The validated \code{kfold_seed} is used only to construct subject-level
#' split manifests and is shared across pair directions with the same subject
#' universe. Fold sampler seeds remain direction-specific.
#'
#' - A diagnostics-aware retry policy and Pathfinder initialization use the fixed
#'   canonical v0.2 settings; their implementation controls are private.
#'
#' **Parallelism & progress**
#' - The **outer pair loop** can run in parallel with `n_workers_outer > 1` (PSOCK via **future/furrr**).
#' - Repeated K-fold fold tasks run sequentially within each directed fit to
#'   avoid nested process-level parallelism.
#' - Main and K-fold fits run their MCMC chains in parallel.
#' - If **progressr** is installed and `progress = "bar"`, outer-pair progress
#'   is displayed; otherwise execution falls back silently.
#'
#' **Interpretation**
#' - For edges, prioritize **posterior-supported pair-to-rest directions** (PSP/LFSR
#'   derived from posterior draws) and **MCMC diagnostics**. A sign is not by
#'   itself evidence of direct biological facilitation or inhibition.
#'   optionally use **ELPD per subject** (from K x R) as supporting evidence. Aggregation is **subject-uniform**.
#'
#' @param physeq A \code{phyloseq} object with taxa in rows (will be transposed if needed).
#' @param taxa_vec Optional character vector of taxa IDs to subset; defaults to all taxa.
#' @param nz_partner_min_frac Minimum partner nonzero fraction for eligibility; default \code{0.15}.
#' @param min_unique_times Minimum unique times per subject; default \code{3}.
#' @param min_pairs Minimum valid adjacent pairs per direction; default \code{4}.
#' @param chains Number of MCMC chains; default \code{4}.
#' @param iter_warmup Warmup iterations per chain; default \code{2000}.
#' @param iter_sampling Post-warmup sampling iterations per chain; default \code{2000}.
#' @param seed Reproducibility seed for fitting; default \code{1234}.
#' @param init Initial value or initializer specification passed to the sampler; default \code{0.2}.
#' @param adapt_delta Target HMC acceptance probability; default \code{0.98}.
#' @param max_treedepth Maximum HMC tree depth; default \code{14}.
#'
#' @param subject_col Column name in sample metadata indicating subjects (e.g., plant/plot).
#' @param time_col    Column name in sample metadata indicating (numeric) time within subject.
#'
#'
#' @param nz_partner_min_frac Drop subjects whose partner predictor is non-zero in fewer than this
#' fraction of rows; default \code{0.15}.
#'
#'
#' @param progress Progress output: \code{"bar"}, \code{"verbose"}, or \code{"none"}.
#' When \code{"bar"} and \pkg{progressr} is available, outer/K-fold progress bars are displayed.
#'
#'
#' @param kfold_K Number of folds K; default \code{5}.
#' @param kfold_R Number of repetitions R; default \code{3}.
#' @param kfold_seed Seed used only for repeated K-fold subject splits; defaults to \code{seed}. The same value is applied across all pair directions, while fold sampler seeds remain direction-specific.
#' @param n_workers_outer Number of parallel workers for the outer pair loop;
#' default \code{1}.
#'
#'
#' @return
#' A list containing directional and self-effect summaries, directed
#' subject-level repeated K-fold ELPD, and the raw bidirectional pair table.
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
#'   # Pathfinder initialization is fixed by the private v0.2 policy.
#'   chains = 4, iter_warmup = 1000, iter_sampling = 1500,
#'   n_workers_outer = 4,   # outer pair loop in parallel
#'   progress = "bar"
#' )
#' }
#'#' @rdname fit_pclv_bayes
#' @export

fit_pclv_bayes <- function(
  physeq, subject_col, time_col, taxa_vec = NULL,
  nz_partner_min_frac = 0.15,
  min_unique_times = 3,
  min_pairs = 4,
  chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  seed = 1234,
  init = 0.2,
  adapt_delta = 0.98,
  max_treedepth = 14,
  progress = c("bar", "verbose", "none"),
  n_workers_outer = 1L,
  kfold_K = 5,
  kfold_R = 3,
  kfold_seed = seed
){
  # Frozen v0.2 canonical policy. These values are intentionally private and
  # are not accepted through ... or a public control list.
  eps <- 1e-6
  zero_mode_alr <- "minpos_time"
  minpos_alpha <- 0.5
  minpos_base <- "ij"
  eps_fixed <- 1e-6
  lib_eps_c <- 0.65
  rest_floor_frac <- 1.0
  smooth_scale <- "logra"
  alr_spline_df <- NULL
  alr_spline_spar <- NULL
  alr_spline_cv <- TRUE
  metric <- "diag_e"
  quiet <- FALSE
  progress_every <- 1L
  silent_sampler <- FALSE
  max_retries <- 3L
  n_workers_kfold <- as.integer(chains)
  use_pathfinder_init <- TRUE
  pf_num_paths <- 8L
  pf_draws <- 1000L
  pf_history_size <- 50L
  pf_max_lbfgs_iters <- 200L
  pf_psis_resample <- TRUE
  public_control_names <- setdiff(names(formals(sys.function())),
                                  c("physeq", "subject_col", "time_col", "taxa_vec"))
  internal_control_names <- c(
    "eps", "zero_mode_alr", "minpos_alpha", "minpos_base", "eps_fixed",
    "lib_eps_c", "rest_floor_frac", "smooth_scale", "alr_spline_df",
    "alr_spline_spar", "alr_spline_cv", "metric", "quiet", "progress_every",
    "silent_sampler", "max_retries", "n_workers_kfold",
    "use_pathfinder_init", "pf_num_paths",
    "pf_draws", "pf_history_size", "pf_max_lbfgs_iters", "pf_psis_resample"
  )
  validated <- .validate_fit_pclv_inputs(
    physeq, subject_col, time_col, taxa_vec,
    controls = mget(c(public_control_names, internal_control_names),
                    envir = environment(), inherits = FALSE)
  )
  for (nm in names(validated$controls))
    assign(nm, validated$controls[[nm]], envir = environment())
  subject_col <- validated$subject_col
  time_col <- validated$time_col
  taxa_vec <- validated$taxa_vec
  meta_df <- validated$meta_df
  mat_rel <- validated$mat_rel


  runtime <- .prepare_fit_runtime(validated)
  if (.is_pclv_failure(runtime)) return(runtime)
  ctx <- runtime$ctx
  # K-fold task-level parallelism is a private fixed policy. Individual
  # fold fits still execute their requested MCMC chains in parallel.
  n_workers_kfold_eff <- n_workers_kfold
  ctx$n_workers_kfold_eff <- n_workers_kfold_eff
  has_progressr <- requireNamespace("progressr", quietly = TRUE)

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
    out <- .outer_pair_map(
      tasks = tasks, taxa_vec = taxa_vec, core_ctx = ctx,
      scheduling = scheduling, progress = progress,
      workers = n_workers_outer, has_progressr = has_progressr
    )
  }

  res <- .assemble_pair_outcomes(out, tasks)
  rownames(res) <- NULL

  # Assemble from the finalized raw result so that any existing attributes on
  # the canonical result table are preserved.
  raw_tbl <- res
  out_obj <- .assemble_public_fit_result(raw_tbl)

  if ("fit" %in% names(raw_tbl)) {
    try({
      ff <- try(raw_tbl$fit$output_files(), silent = TRUE)
      if (!inherits(ff, "try-error") && length(ff) > 0) {
        invisible(lapply(ff, function(f) if (file.exists(f)) unlink(f)))
      }
    }, silent = TRUE)
  }

  return(out_obj)
}
