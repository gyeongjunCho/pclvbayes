# R/fit_pclv_laplace.R

# -------------------------------------------------------------------------
# pcLV Laplace-only audit fit
#
# Purpose:
#   Run the canonical pcLV preprocessing and Stan model, but perform only
#   MAP optimization + Laplace approximation. No HMC, Pathfinder, retry,
#   diagnostics, or K-fold evaluation is run.
#
# Intended use:
#   Calibrate Laplace sign/LFSR against an existing HMC reference.
# -------------------------------------------------------------------------


.laplace_cleanup_fit <- function(fit) {
  if (is.null(fit)) return(invisible(NULL))

  files <- tryCatch(
    fit$output_files(),
    error = function(e) character()
  )

  if (length(files)) {
    invisible(lapply(
      files,
      function(f) {
        if (is.character(f) &&
            length(f) == 1L &&
            !is.na(f) &&
            file.exists(f)) {
          try(unlink(f), silent = TRUE)
        }
      }
    ))
  }

  invisible(NULL)
}


.laplace_failure_result <- function(
    target,
    partner,
    direction,
    stage,
    reason,
    message = NA_character_,
    n_pairs = NA_integer_,
    N = NA_integer_,
    S = NA_integer_,
    map_ok = FALSE,
    map_return_code = NA_integer_,
    laplace_return_code = NA_integer_,
    optimize_sec = NA_real_,
    laplace_sec = NA_real_,
    runtime_sec = NA_real_,
    laplace_draws_requested = NA_integer_) {

  list(
    from = partner,
    to = target,
    direction = direction,

    laplace_ok = FALSE,
    map_ok = isTRUE(map_ok),

    failure_stage = as.character(stage),
    failure_reason = as.character(reason),
    failure_message = as.character(message),

    n_pairs = as.integer(n_pairs),
    N = as.integer(N),
    S = as.integer(S),

    map_return_code = as.integer(map_return_code),
    laplace_return_code = as.integer(laplace_return_code),

    map_aij = NA_real_,

    a_mean = NA_real_,
    a_median = NA_real_,
    a_sd = NA_real_,
    a_q2.5 = NA_real_,
    a_q97.5 = NA_real_,

    positive_sign_probability = NA_real_,
    negative_sign_probability = NA_real_,
    dominant_sign_probability = NA_real_,
    lfsr = NA_real_,
    p_sign2 = NA_real_,

    laplace_draws_requested =
      as.integer(laplace_draws_requested),
    laplace_draws_finite = 0L,
    finite_draw_fraction = NA_real_,

    mode_mean_shift_sd = NA_real_,

    optimize_sec = as.numeric(optimize_sec),
    laplace_sec = as.numeric(laplace_sec),
    runtime_sec = as.numeric(runtime_sec)
  )
}


.fit_direction_laplace <- function(
    target,
    partner,
    ctx,
    pair_state,
    direction,
    seed_override,
    laplace_draws = 1000L) {

  total_start <- proc.time()[["elapsed"]]

  fail <- function(
      stage,
      reason,
      message = NA_character_,
      n_pairs = NA_integer_,
      N = NA_integer_,
      S = NA_integer_,
      map_ok = FALSE,
      map_return_code = NA_integer_,
      laplace_return_code = NA_integer_,
      optimize_sec = NA_real_,
      laplace_sec = NA_real_) {

    .laplace_failure_result(
      target = target,
      partner = partner,
      direction = direction,
      stage = stage,
      reason = reason,
      message = message,
      n_pairs = n_pairs,
      N = N,
      S = S,
      map_ok = map_ok,
      map_return_code = map_return_code,
      laplace_return_code = laplace_return_code,
      optimize_sec = optimize_sec,
      laplace_sec = laplace_sec,
      runtime_sec =
        proc.time()[["elapsed"]] - total_start,
      laplace_draws_requested = laplace_draws
    )
  }

  if (is.null(pair_state)) {
    return(fail(
      "preprocessing",
      "missing_precomputed_pair_state"
    ))
  }

  if (!(direction %in% c("ij", "ji"))) {
    return(fail(
      "preprocessing",
      "invalid_direction",
      direction
    ))
  }

  # -----------------------------------------------------------------------
  # 1. Canonical precomputed directed input
  # -----------------------------------------------------------------------

  pair_in <- tryCatch(
    .make_pair_inputs_precomputed(
      pair_state = pair_state,
      subject = ctx$pair_subject,
      time = ctx$pair_time,
      direction = direction
    ),
    error = function(e) {
      .pclv_failure(
        "preprocessing",
        "pair_input_construction_failed",
        list(message = conditionMessage(e))
      )
    }
  )

  if (.is_pclv_failure(pair_in)) {
    msg <- if (!is.null(pair_in$details$message)) {
      pair_in$details$message
    } else {
      NA_character_
    }

    return(fail(
      stage = pair_in$stage %||% "preprocessing",
      reason = pair_in$reason %||% "pair_input_failed",
      message = msg
    ))
  }

  if (nrow(pair_in) < ctx$min_pairs) {
    return(fail(
      "preprocessing",
      "insufficient_analysis_rows",
      n_pairs = nrow(pair_in)
    ))
  }

  pair_in <- pair_in[
    is.finite(pair_in$y) &
      is.finite(pair_in$xi) &
      is.finite(pair_in$xj) &
      is.finite(pair_in$time),
    ,
    drop = FALSE
  ]

  if (!nrow(pair_in)) {
    return(fail(
      "preprocessing",
      "no_finite_analysis_rows"
    ))
  }

  pair_in <- pair_in[
    order(pair_in$subject, pair_in$time),
    ,
    drop = FALSE
  ]

  N <- nrow(pair_in)

  # -----------------------------------------------------------------------
  # 2. Same prev/dt construction as canonical HMC
  # -----------------------------------------------------------------------

  prev <- integer(N)
  dtv <- numeric(N)

  by_s <- split(seq_len(N), pair_in$subject)

  for (sb in names(by_s)) {
    ix <- by_s[[sb]]

    prev[ix[[1L]]] <- 0L
    dtv[ix[[1L]]] <- 0

    if (length(ix) >= 2L) {
      for (k in 2:length(ix)) {
        prev[ix[[k]]] <- ix[[k - 1L]]

        dt_k <- as.numeric(
          pair_in$time[ix[[k]]] -
            pair_in$time[ix[[k - 1L]]]
        )

        if (!is.finite(dt_k) || dt_k <= 0) {
          return(fail(
            "preprocessing",
            "invalid_time_interval",
            n_pairs = N,
            N = N
          ))
        }

        dtv[ix[[k]]] <- dt_k
      }
    }
  }

  n_subj <- length(unique(pair_in$subject))

  stan_list <- list(
    N = N,
    y = pair_in$y,
    xi = pair_in$xi,
    xj = pair_in$xj,
    S = n_subj,
    sid = as.integer(factor(pair_in$subject)),
    prev = prev,
    dt = dtv
  )

  if (.is_pclv_failure(ctx$mod)) {
    return(fail(
      "model_loading",
      "compiled_model_unavailable",
      n_pairs = N,
      N = N,
      S = n_subj
    ))
  }

  if (is.null(ctx$mod)) {
    return(fail(
      "model_loading",
      "missing_worker_model",
      n_pairs = N,
      N = N,
      S = n_subj
    ))
  }

  mod <- ctx$mod
  seed_main <- as.integer(seed_override)

  fit_map <- NULL
  fit_laplace <- NULL

  on.exit({
    .laplace_cleanup_fit(fit_laplace)
    .laplace_cleanup_fit(fit_map)
  }, add = TRUE)

  # -----------------------------------------------------------------------
  # 3. MAP
  #
  # jacobian = TRUE is essential here: this is the posterior mode used by
  # the Laplace approximation, not the unconstrained penalized MLE.
  # -----------------------------------------------------------------------

  opt_start <- proc.time()[["elapsed"]]

  fit_map <- tryCatch(
    mod$optimize(
      data = stan_list,
      seed = seed_main,
      init = ctx$init,
      jacobian = TRUE,
      refresh = 0,
      show_messages = FALSE,
      show_exceptions = FALSE
    ),
    error = function(e) e
  )

  optimize_sec <-
    proc.time()[["elapsed"]] - opt_start

  if (inherits(fit_map, "error")) {
    return(fail(
      "optimization",
      "map_exception",
      message = conditionMessage(fit_map),
      n_pairs = N,
      N = N,
      S = n_subj,
      optimize_sec = optimize_sec
    ))
  }

  map_rc <- tryCatch(
    as.integer(fit_map$return_codes()[[1L]]),
    error = function(e) NA_integer_
  )

  if (!is.finite(map_rc) || map_rc != 0L) {
    output_message <- tryCatch(
      paste(fit_map$output(), collapse = "\n"),
      error = function(e) NA_character_
    )

    return(fail(
      "optimization",
      "map_nonzero_return_code",
      message = output_message,
      n_pairs = N,
      N = N,
      S = n_subj,
      map_ok = FALSE,
      map_return_code = map_rc,
      optimize_sec = optimize_sec
    ))
  }

  map_aij <- tryCatch(
    as.numeric(fit_map$mle("a_ij")[[1L]]),
    error = function(e) NA_real_
  )

  if (!is.finite(map_aij)) {
    return(fail(
      "optimization",
      "invalid_map_aij",
      n_pairs = N,
      N = N,
      S = n_subj,
      map_ok = TRUE,
      map_return_code = map_rc,
      optimize_sec = optimize_sec
    ))
  }

  # -----------------------------------------------------------------------
  # 4. Laplace approximation
  # -----------------------------------------------------------------------

  lap_start <- proc.time()[["elapsed"]]

  fit_laplace <- tryCatch(
    mod$laplace(
      data = stan_list,
      seed = seed_main,
      mode = fit_map,
      jacobian = TRUE,
      draws = as.integer(laplace_draws),
      refresh = 0,
      show_messages = FALSE,
      show_exceptions = FALSE
    ),
    error = function(e) e
  )

  laplace_sec <-
    proc.time()[["elapsed"]] - lap_start

  if (inherits(fit_laplace, "error")) {
    return(fail(
      "laplace",
      "laplace_exception",
      message = conditionMessage(fit_laplace),
      n_pairs = N,
      N = N,
      S = n_subj,
      map_ok = TRUE,
      map_return_code = map_rc,
      optimize_sec = optimize_sec,
      laplace_sec = laplace_sec
    ))
  }

  lap_rc <- tryCatch(
    as.integer(fit_laplace$return_codes()[[1L]]),
    error = function(e) NA_integer_
  )

  if (!is.finite(lap_rc) || lap_rc != 0L) {
    output_message <- tryCatch(
      paste(fit_laplace$output(), collapse = "\n"),
      error = function(e) NA_character_
    )

    return(fail(
      "laplace",
      "laplace_nonzero_return_code",
      message = output_message,
      n_pairs = N,
      N = N,
      S = n_subj,
      map_ok = TRUE,
      map_return_code = map_rc,
      laplace_return_code = lap_rc,
      optimize_sec = optimize_sec,
      laplace_sec = laplace_sec
    ))
  }

  # Read only a_ij into R. CmdStan itself generated the full Laplace draw file,
  # but we do not retain any other parameter draws.
  a_draws <- tryCatch(
    {
      x <- fit_laplace$draws(
        variables = "a_ij",
        format = "matrix"
      )
      as.numeric(x[, "a_ij"])
    },
    error = function(e) e
  )

  if (inherits(a_draws, "error")) {
    return(fail(
      "posterior_extraction",
      "aij_draw_extraction_failed",
      message = conditionMessage(a_draws),
      n_pairs = N,
      N = N,
      S = n_subj,
      map_ok = TRUE,
      map_return_code = map_rc,
      laplace_return_code = lap_rc,
      optimize_sec = optimize_sec,
      laplace_sec = laplace_sec
    ))
  }

  finite <- is.finite(a_draws)
  n_finite <- sum(finite)

  if (!n_finite) {
    return(fail(
      "posterior_extraction",
      "no_finite_aij_draws",
      n_pairs = N,
      N = N,
      S = n_subj,
      map_ok = TRUE,
      map_return_code = map_rc,
      laplace_return_code = lap_rc,
      optimize_sec = optimize_sec,
      laplace_sec = laplace_sec
    ))
  }

  a <- a_draws[finite]

  a_mean <- mean(a)
  a_median <- stats::median(a)
  a_sd <- stats::sd(a)

  q <- stats::quantile(
    a,
    probs = c(0.025, 0.975),
    names = FALSE,
    type = 8
  )

  p_pos <- mean(a > 0)
  p_neg <- mean(a < 0)

  # With a continuous Laplace approximation exact zeros should be negligible,
  # but normalize defensively if they occur.
  sign_mass <- p_pos + p_neg
  if (is.finite(sign_mass) && sign_mass > 0) {
    p_pos <- p_pos / sign_mass
    p_neg <- p_neg / sign_mass
  } else {
    p_pos <- NA_real_
    p_neg <- NA_real_
  }

  lfsr <- if (
    is.finite(p_pos) &&
      is.finite(p_neg)
  ) {
    min(p_pos, p_neg)
  } else {
    NA_real_
  }

  p_sign2 <- if (is.finite(lfsr)) {
    min(1, 2 * lfsr)
  } else {
    NA_real_
  }

  mode_mean_shift_sd <- if (
    is.finite(a_sd) &&
      a_sd > 0
  ) {
    (a_mean - map_aij) / a_sd
  } else {
    NA_real_
  }

  list(
    from = partner,
    to = target,
    direction = direction,

    laplace_ok = TRUE,
    map_ok = TRUE,

    failure_stage = NA_character_,
    failure_reason = NA_character_,
    failure_message = NA_character_,

    n_pairs = as.integer(N),
    N = as.integer(N),
    S = as.integer(n_subj),

    map_return_code = map_rc,
    laplace_return_code = lap_rc,

    map_aij = map_aij,

    a_mean = a_mean,
    a_median = a_median,
    a_sd = a_sd,
    a_q2.5 = q[[1L]],
    a_q97.5 = q[[2L]],

    positive_sign_probability = p_pos,
    negative_sign_probability = p_neg,
    dominant_sign_probability = max(p_pos, p_neg),
    lfsr = lfsr,
    p_sign2 = p_sign2,

    laplace_draws_requested =
      as.integer(laplace_draws),
    laplace_draws_finite =
      as.integer(n_finite),
    finite_draw_fraction =
      n_finite / length(a_draws),

    # Useful audit quantity:
    # how far the approximate posterior mean is from the MAP in SD units.
    mode_mean_shift_sd = mode_mean_shift_sd,

    optimize_sec = optimize_sec,
    laplace_sec = laplace_sec,
    runtime_sec =
      proc.time()[["elapsed"]] - total_start
  )
}


.run_pair_laplace <- function(
    task,
    taxa_vec,
    core_ctx,
    laplace_draws,
    progress = "none",
    mute_logs = TRUE) {

  ti <- taxa_vec[[task$idx_i]]
  pj <- taxa_vec[[task$idx_j]]

  task_ctx <- core_ctx

  pair_failure <-
    !is.null(task$pair_state) &&
    .is_pclv_failure(task$pair_state$failure)

  if (pair_failure) {
    f <- task$pair_state$failure

    ij <- .laplace_failure_result(
      target = ti,
      partner = pj,
      direction = "ij",
      stage = f$stage %||% "preprocessing",
      reason = f$reason %||% "pair_precompute_failed",
      laplace_draws_requested = laplace_draws
    )

    ji <- .laplace_failure_result(
      target = pj,
      partner = ti,
      direction = "ji",
      stage = f$stage %||% "preprocessing",
      reason = f$reason %||% "pair_precompute_failed",
      laplace_draws_requested = laplace_draws
    )

    return(dplyr::bind_rows(ij, ji))
  }

  have_executable <-
    is.character(core_ctx$mod_exe_file) &&
    length(core_ctx$mod_exe_file) == 1L &&
    nzchar(core_ctx$mod_exe_file) &&
    file.exists(core_ctx$mod_exe_file)

  if (!have_executable) {
    ij <- .laplace_failure_result(
      target = ti,
      partner = pj,
      direction = "ij",
      stage = "model_loading",
      reason = "compiled_model_unavailable",
      laplace_draws_requested = laplace_draws
    )

    ji <- .laplace_failure_result(
      target = pj,
      partner = ti,
      direction = "ji",
      stage = "model_loading",
      reason = "compiled_model_unavailable",
      laplace_draws_requested = laplace_draws
    )

    return(dplyr::bind_rows(ij, ji))
  }

  task_ctx$mod <- tryCatch(
    cmdstanr::cmdstan_model(
      stan_file = NULL,
      exe_file = core_ctx$mod_exe_file
    ),
    error = function(e) {
      .pclv_failure(
        "model_loading",
        "compiled_model_load_failed",
        list(message = conditionMessage(e))
      )
    }
  )

  if (.is_pclv_failure(task_ctx$mod)) {
    f <- task_ctx$mod

    ij <- .laplace_failure_result(
      target = ti,
      partner = pj,
      direction = "ij",
      stage = f$stage,
      reason = f$reason,
      laplace_draws_requested = laplace_draws
    )

    ji <- .laplace_failure_result(
      target = pj,
      partner = ti,
      direction = "ji",
      stage = f$stage,
      reason = f$reason,
      laplace_draws_requested = laplace_draws
    )

    return(dplyr::bind_rows(ij, ji))
  }

  seeds <- task$direction_seeds

  if (!mute_logs && progress == "verbose") {
    cat(sprintf(
      "Laplace pair %s-%s started\n",
      ti, pj
    ))
  }

  ij <- .fit_direction_laplace(
    target = ti,
    partner = pj,
    ctx = task_ctx,
    pair_state = task$pair_state,
    direction = "ij",
    seed_override = seeds[["ij"]],
    laplace_draws = laplace_draws
  )

  ji <- .fit_direction_laplace(
    target = pj,
    partner = ti,
    ctx = task_ctx,
    pair_state = task$pair_state,
    direction = "ji",
    seed_override = seeds[["ji"]],
    laplace_draws = laplace_draws
  )

  if (!mute_logs && progress == "verbose") {
    cat(sprintf(
      "Laplace pair %s-%s completed\n",
      ti, pj
    ))
  }

  dplyr::bind_rows(ij, ji)
}


.execute_laplace_pair_task <- function(
    task,
    taxa_vec,
    core_ctx,
    laplace_draws,
    progress = "none",
    mute_logs = TRUE) {

  result <- .run_pair_laplace(
    task = task,
    taxa_vec = taxa_vec,
    core_ctx = core_ctx,
    laplace_draws = laplace_draws,
    progress = progress,
    mute_logs = mute_logs
  )

  list(
    task_index = task$task_index,
    result = result
  )
}


.execute_laplace_pair_task_progress <- function(
    task,
    taxa_vec,
    core_ctx,
    laplace_draws,
    progress,
    progressor) {

  result <- .execute_laplace_pair_task(
    task = task,
    taxa_vec = taxa_vec,
    core_ctx = core_ctx,
    laplace_draws = laplace_draws,
    progress = progress,
    mute_logs = TRUE
  )

  progressor(
    message = sprintf(
      "pair %s-%s",
      task$taxon_i,
      task$taxon_j
    )
  )

  result
}


.outer_laplace_pair_map <- function(
    tasks,
    taxa_vec,
    core_ctx,
    laplace_draws,
    progress,
    workers,
    has_progressr) {

  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)

  future::plan(
    future::multisession,
    workers = workers
  )

  opts <- furrr::furrr_options(
    seed = TRUE,
    globals = FALSE,
    packages = "pclvbayes",
    scheduling = Inf
  )

  common <- list(
    .x = tasks,
    taxa_vec = taxa_vec,
    core_ctx = core_ctx,
    laplace_draws = laplace_draws,
    progress = progress,
    .options = opts
  )

  if (
    isTRUE(has_progressr) &&
    identical(progress, "bar")
  ) {
    return(
      progressr::with_progress({
        p <- progressr::progressor(
          steps = length(tasks)
        )

        do.call(
          furrr::future_map,
          c(
            common,
            list(
              .f = .execute_laplace_pair_task_progress,
              progressor = p
            )
          )
        )
      })
    )
  }

  do.call(
    furrr::future_map,
    c(
      common,
      list(
        .f = .execute_laplace_pair_task,
        mute_logs = TRUE
      )
    )
  )
}


#' Fast Laplace approximation of canonical pcLV directions
#'
#' @description
#' Runs the same deterministic preprocessing and canonical Stan model used by
#' fit_pclv_bayes(), but performs MAP optimization followed by a Laplace
#' approximation only. It never runs NUTS-HMC, Pathfinder, retries, or K-fold.
#'
#' This function is intended for calibration and future pre-HMC screening.
#'
#' @param physeq A phyloseq object.
#' @param subject_col Subject identifier column.
#' @param time_col Numeric time column.
#' @param taxa_vec Optional taxa subset.
#' @param nz_partner_min_frac Same canonical pcLV eligibility threshold.
#' @param min_unique_times Same canonical pcLV minimum unique times.
#' @param min_pairs Same canonical minimum adjacent-pair count.
#' @param seed Reproducibility seed.
#' @param init Initial value passed to MAP optimization.
#' @param laplace_draws Number of approximate posterior draws per direction.
#' @param progress "bar", "verbose", or "none".
#' @param n_workers_outer Number of concurrent unordered-pair workers.
#'
#' @return A list with directed Laplace results in `cross` and configuration.
#'
#' @export
fit_pclv_laplace <- function(
    physeq,
    subject_col,
    time_col,
    taxa_vec = NULL,
    nz_partner_min_frac = 0.15,
    min_unique_times = 3,
    min_pairs = 4,
    seed = 1234,
    init = 0.2,
    laplace_draws = 1000L,
    progress = c("bar", "verbose", "none"),
    n_workers_outer = 1L) {

  progress <- match.arg(progress)

  if (
    length(laplace_draws) != 1L ||
    !is.numeric(laplace_draws) ||
    !is.finite(laplace_draws) ||
    laplace_draws < 100 ||
    floor(laplace_draws) != laplace_draws
  ) {
    stop(
      "`laplace_draws` must be an integer >= 100.",
      call. = FALSE
    )
  }

  laplace_draws <- as.integer(laplace_draws)

  # -----------------------------------------------------------------------
  # Reuse the canonical pcLV validator.
  #
  # HMC/K-fold controls below are supplied only because the existing validator
  # validates the full canonical control schema. They are NEVER executed by
  # this function.
  # -----------------------------------------------------------------------

  controls <- list(
    eps = 1e-6,
    min_pairs = min_pairs,
    min_unique_times = min_unique_times,
    zero_mode_alr = "minpos_time",
    minpos_alpha = 0.5,
    minpos_base = "ij",
    eps_fixed = 1e-6,
    lib_eps_c = 0.65,
    rest_floor_frac = 1.0,
    smooth_scale = "logra",
    alr_spline_df = NULL,
    alr_spline_spar = NULL,
    alr_spline_cv = TRUE,
    nz_partner_min_frac = nz_partner_min_frac,

    # Unused HMC placeholders required by canonical validation:
    max_retries = 0L,
    chains = 1L,
    iter_warmup = 0L,
    iter_sampling = 1L,
    adapt_delta = 0.98,
    max_treedepth = 14L,
    metric = "diag_e",

    init = init,
    seed = seed,
    quiet = TRUE,
    progress = progress,
    progress_every = 1L,
    silent_sampler = TRUE,
    n_workers_outer = n_workers_outer,

    # Unused predictive placeholders:
    kfold_K = 2L,
    kfold_R = 1L,
    kfold_seed = seed,

    # Explicitly disabled:
    use_pathfinder_init = FALSE,
    pf_num_paths = 1L,
    pf_draws = 1L,
    pf_history_size = 1L,
    pf_max_lbfgs_iters = 1L,
    pf_psis_resample = TRUE
  )

  validated <- .validate_fit_pclv_inputs(
    physeq = physeq,
    subject_col = subject_col,
    time_col = time_col,
    taxa_vec = taxa_vec,
    controls = controls
  )

  controls <- validated$controls
  taxa_vec <- validated$taxa_vec

  old_thread_env <-
    .pclv_force_single_thread_libraries()

  on.exit(
    .pclv_restore_thread_libraries(
      old_thread_env
    ),
    add = TRUE
  )

  if (!identical(progress, "none")) {
    cat(
      "Preparing canonical pcLV trajectories for Laplace-only fitting.\n"
    )
  }

  runtime <- .prepare_fit_runtime(validated)

  if (.is_pclv_failure(runtime)) {
    return(runtime)
  }

  ctx <- runtime$ctx

  pair_precompute <-
    .precompute_pair_states_vectorized(
      sm_mat = ctx$sm_mat,
      meta_df = ctx$meta_df,
      taxa_vec = taxa_vec,
      zero_mode_alr = controls$zero_mode_alr,
      minpos_alpha = controls$minpos_alpha,
      minpos_base = controls$minpos_base,
      eps_fixed = controls$eps_fixed,
      rest_floor_frac = controls$rest_floor_frac,
      alr_cap = .PCLV_CORE_ALR_CAP,
      smooth_scale = controls$smooth_scale,
      nz_partner_min_frac =
        controls$nz_partner_min_frac,
      min_pairs = controls$min_pairs
    )

  if (.is_pclv_failure(pair_precompute)) {
    return(pair_precompute)
  }

  pair_states <- pair_precompute$states

  ctx$pair_subject <-
    pair_precompute$subject
  ctx$pair_time <-
    pair_precompute$time

  # Workers do not need the large originals.
  ctx$sm_mat <- NULL
  ctx$meta_df <- NULL

  tasks <- .make_pair_tasks(
    taxa_vec,
    controls$seed,
    pair_states = pair_states
  )

  rm(
    pair_states,
    pair_precompute,
    runtime,
    validated,
    physeq
  )

  if (!identical(progress, "none")) {
    cat(sprintf(
      paste0(
        "Starting Laplace-only fitting: ",
        "%d unordered pairs / %d directed models, ",
        "%d approximate draws each.\n"
      ),
      length(tasks),
      2L * length(tasks),
      laplace_draws
    ))
  }

  overall_start <-
    proc.time()[["elapsed"]]

  has_progressr <-
    requireNamespace(
      "progressr",
      quietly = TRUE
    )

  if (controls$n_workers_outer <= 1L) {

    pb <- NULL

    if (identical(progress, "bar")) {
      pb <- utils::txtProgressBar(
        min = 0,
        max = length(tasks),
        style = 3
      )

      on.exit(
        try(close(pb), silent = TRUE),
        add = TRUE
      )
    }

    out <- vector(
      "list",
      length(tasks)
    )

    for (k in seq_along(tasks)) {

      out[[k]] <-
        .execute_laplace_pair_task(
          task = tasks[[k]],
          taxa_vec = taxa_vec,
          core_ctx = ctx,
          laplace_draws = laplace_draws,
          progress = progress,
          mute_logs = FALSE
        )

      if (!is.null(pb)) {
        utils::setTxtProgressBar(
          pb,
          k
        )
      }
    }

  } else {

    out <- .outer_laplace_pair_map(
      tasks = tasks,
      taxa_vec = taxa_vec,
      core_ctx = ctx,
      laplace_draws = laplace_draws,
      progress = progress,
      workers =
        controls$n_workers_outer,
      has_progressr =
        has_progressr
    )
  }

  ord <- order(
    vapply(
      out,
      `[[`,
      integer(1),
      "task_index"
    )
  )

  cross <- dplyr::bind_rows(
    lapply(
      out[ord],
      `[[`,
      "result"
    )
  )

  overall_sec <-
    proc.time()[["elapsed"]] -
    overall_start

  configuration <- list(
    inference = "laplace_only",
    canonical_model = "pclv.stan",
    taxa_count = length(taxa_vec),
    directed_model_count = nrow(cross),
    unordered_pair_count = length(tasks),
    laplace_draws = laplace_draws,
    seed = controls$seed,
    init = controls$init,
    n_workers_outer =
      controls$n_workers_outer,
    nz_partner_min_frac =
      controls$nz_partner_min_frac,
    min_unique_times =
      controls$min_unique_times,
    min_pairs =
      controls$min_pairs,
    jacobian = TRUE,
    hmc_run = FALSE,
    pathfinder_run = FALSE,
    kfold_run = FALSE,
    runtime_sec = overall_sec
  )

  structure(
    list(
      cross = cross,
      configuration = configuration
    ),
    class = c(
      "pclv_laplace_fit",
      "list"
    )
  )
}
