# Benchmark-only full 100-species execution contract for ROADMAP V021-06.

v021_full_schema <- "v021_full_100_species_v5"
v021_full_status_schema <- "v021_full_100_species_status_v1"
v021_full_failure_trace_schema <- "v021_full_failure_trace_v1"
v021_full_cleanup_audit_schema <- "v021_full_cleanup_audit_v1"
v021_full_selection_rule <- paste(
  "MTIST 100-species metadata; noise=0.01; even sampling; at least 10 series;",
  "at least 15 timepoints; lowest dataset ID; no truth fields inspected"
)

new_v021_child_job_tracker <- function(jobs = list()) {
  tracker <- new.env(parent = emptyenv())
  tracker$jobs <- jobs
  tracker
}

collect_v021_tracked_jobs <- function(tracker, collector = parallel::mccollect,
                                      wait = TRUE) {
  if (!is.environment(tracker) || is.null(tracker$jobs))
    stop("Invalid V021 child-job tracker.")
  jobs <- tracker$jobs
  if (!length(jobs)) return(list())
  value <- collector(jobs, wait = wait)
  tracker$jobs <- list()
  value
}

reap_v021_uncollected_jobs <- function(..., collector = parallel::mccollect) {
  trackers <- list(...)
  jobs <- unlist(lapply(trackers, function(x) {
    if (!is.environment(x) || is.null(x$jobs)) stop("Invalid V021 child-job tracker.")
    x$jobs
  }), recursive = FALSE)
  if (!length(jobs)) return(list(class = "none", child_pids = integer()))
  pids <- vapply(jobs, function(job) as.integer(job$pid), integer(1))
  collected <- tryCatch(collector(jobs, wait = FALSE), error = identity)
  lapply(trackers, function(x) x$jobs <- list())
  list(class = class(collected), child_pids = pids)
}

v021_proc_identity <- function(
    pid,
    proc_root = "/proc",
    stat_reader = function(path) readLines(path, warn = FALSE, n = 1L)) {
  pid <- as.integer(pid)
  if (length(pid) != 1L || is.na(pid) || pid < 1L)
    stop("Process PID must be one positive integer.")

  stat_path <- file.path(proc_root, as.character(pid), "stat")
  stat <- tryCatch(stat_reader(stat_path), error = function(e) character())
  if (!length(stat)) return(NULL)

  parsed <- tryCatch(
    .v021_parse_linux_stat(stat, pid),
    error = function(e) NULL
  )
  if (is.null(parsed)) {
    return(list(
      pid = pid,
      start_time = NA_character_,
      process_state = NA_character_,
      readable = FALSE
    ))
  }

  list(
    pid = pid,
    start_time = as.character(parsed$start_time),
    process_state = as.character(parsed$process_state),
    readable = TRUE
  )
}

wait_v021_collected_child_exit <- function(
    pid,
    expected_start_time,
    timeout_seconds = 10,
    poll_seconds = 0.01,
    identity_reader = v021_proc_identity,
    sleep = Sys.sleep,
    clock = function() proc.time()[["elapsed"]]) {
  pid <- as.integer(pid)
  if (length(pid) != 1L || is.na(pid) || pid < 1L)
    stop("Collected child PID must be one positive integer.")
  if (!is.character(expected_start_time) ||
      length(expected_start_time) != 1L ||
      is.na(expected_start_time) ||
      !nzchar(expected_start_time))
    stop("Collected child start time must be one non-empty character value.")
  if (!is.numeric(timeout_seconds) ||
      length(timeout_seconds) != 1L ||
      !is.finite(timeout_seconds) ||
      timeout_seconds < 0)
    stop("timeout_seconds must be one finite nonnegative value.")
  if (!is.numeric(poll_seconds) ||
      length(poll_seconds) != 1L ||
      !is.finite(poll_seconds) ||
      poll_seconds <= 0)
    stop("poll_seconds must be one finite positive value.")
  if (!is.function(identity_reader) ||
      !is.function(sleep) ||
      !is.function(clock))
    stop("Child-exit readers and timing hooks must be functions.")

  deadline <- clock() + timeout_seconds

  repeat {
    identity <- identity_reader(pid)

    # Missing /proc identity means that the collected process has exited.
    if (is.null(identity))
      return(invisible(TRUE))

    if (!is.list(identity) ||
        is.null(identity$start_time) ||
        length(identity$start_time) != 1L ||
        is.na(identity$start_time) ||
        !nzchar(as.character(identity$start_time)) ||
        is.null(identity$process_state) ||
        length(identity$process_state) != 1L ||
        is.na(identity$process_state) ||
        !nzchar(as.character(identity$process_state)))
      stop("Invalid collected-child process identity.")

    # The PID may have been reused after the monitor exited. A different
    # start-time identity is not the collected child and must not fail the run.
    if (!identical(as.character(identity$start_time), expected_start_time))
      return(invisible(TRUE))

    # A zombie has already terminated. It may remain briefly visible in /proc
    # until its parent completes reaping, but it is not a live child.
    if (identical(as.character(identity$process_state), "Z"))
      return(invisible(TRUE))

    if (clock() >= deadline)
      break

    sleep(poll_seconds)
  }

  stop("Collected child identity remained in the process table after collection.")
}


recover_v021_monitor_peaks <- function(output_root, policy, executable) {
  validate_v021_resource_policy(policy)
  monitor_root <- file.path(output_root, "monitor")
  files <- list.files(monitor_root, pattern = "^batch-[0-9]+(-kfold)?\\.rds$",
                      full.names = TRUE)
  if (!length(files)) return(list(chains = 0L, processes = 0L))
  peaks <- lapply(files, function(path) {
    batch_id <- sub("\\.rds$", "", basename(path))
    payload <- readRDS(path)
    ownership <- readRDS(file.path(output_root, "batch_ownership",
                                   paste0(batch_id, ".rds")))
    if (!is.list(payload) || !is.list(payload$snapshots) ||
        !length(payload$snapshots) || !is.list(ownership) ||
        is.null(ownership$parent_pid) || is.null(ownership$worker_registry))
      stop("Invalid durable V021 monitor evidence.")
    monitor <- monitor_v021_process_snapshots(
      payload$snapshots, policy, as.integer(ownership$parent_pid), executable,
      worker_registry = ownership$worker_registry)
    validate_v021_preflight_monitor(monitor, policy)
    c(chains = monitor$observed_peak_active_cmdstan_chains,
      processes = monitor$observed_peak_active_cmdstan_processes)
  })
  peaks <- do.call(rbind, peaks)
  list(chains = as.integer(max(peaks[, "chains"])),
       processes = as.integer(max(peaks[, "processes"])))
}

reconcile_v021_dead_controller_reservations <- function(reservations, policy,
                                                         controller_alive) {
  validate_v021_reservations(reservations, policy)
  if (isTRUE(controller_alive))
    stop("Live-controller reservations cannot be administratively reconciled.")
  reservations$state[reservations$state == "reserved"] <- "worker_terminated"
  validate_v021_reservations(reservations, policy)
  reservations
}

v021_task_attempt_number <- function(status, task_ordinal, next_attempt = FALSE) {
  if (!is.list(status) || !is.data.frame(status$tasks))
    stop("Invalid V021 task status for reservation identity.")
  i <- match(as.integer(task_ordinal), status$tasks$task_ordinal)
  if (is.na(i)) stop("Unknown task ordinal for reservation identity.")
  value <- as.integer(status$tasks$attempt_count[[i]])
  if (isTRUE(next_attempt)) value <- value + 1L
  if (value < 1L) stop("Reservation attempt identity is not active.")
  value
}

with_v021_predictive_reservation <- function(
    result, reservations, policy, task_identity, attempt_number,
    evaluator, persist = function(x) invisible(NULL)) {
  if (!is.function(evaluator) || !is.function(persist))
    stop("Predictive evaluator and reservation persistence must be functions.")
  current <- reserve_v021_capacity(
    reservations, policy, task_identity, attempt_number, "kfold_fit")
  persist(current)
  released <- FALSE
  on.exit({
    if (!released) {
      current <- release_v021_capacity(
        current, policy, task_identity, attempt_number, "kfold_fit",
        worker_terminated = TRUE)
      persist(current)
      released <- TRUE
    }
  }, add = TRUE)
  value <- evaluator(result)
  current <- release_v021_capacity(
    current, policy, task_identity, attempt_number, "kfold_fit",
    worker_terminated = TRUE)
  persist(current)
  released <- TRUE
  list(result = value, reservations = current)
}


# Build exact repeated-K-fold tasks for one directed main-posterior result.
# The scientific split, fold seed, fit, and aggregation contracts mirror
# pclvbayes:::.repkfold_eval(); only the execution topology is moved to the
# benchmark controller so fold fits from different directions share one
# global 12-slot queue.
prepare_v021_predictive_fold_plan <- function(main_result, position,
                                               direction_id) {
  if (!is.numeric(position) || length(position) != 1L ||
      is.na(position) || position < 1L ||
      position != as.integer(position))
    stop("Predictive direction position must be one positive integer.")
  if (!is.character(direction_id) || length(direction_id) != 1L ||
      is.na(direction_id) || !nzchar(direction_id))
    stop("Predictive direction identity must be one non-empty string.")
  if (inherits(main_result, c("try-error", "pclv_failure",
                              "v021_traced_child_error")))
    stop("Only successful main-posterior results can build predictive folds.")

  predictive <- main_result$.predictive_context
  if (is.null(predictive))
    stop("Main-posterior result lacks predictive context.")
  result <- main_result
  result$.predictive_context <- NULL

  diagnostic_failure <- result$diagnostic_failure
  if (is.list(diagnostic_failure) && length(diagnostic_failure) >= 1L &&
      !is.null(diagnostic_failure[[1L]])) {
    return(list(
      position = as.integer(position), direction_id = direction_id,
      eligible = FALSE, result = result, predictive = predictive,
      splits_df = NULL, subject_ids = character(),
      subject_holdout_counts = integer(), tasks = list(),
      queue_indices = integer()
    ))
  }

  required <- c(
    "mod", "stan_list", "pair_in", "sample_args", "split_seed",
    "sampling_seed", "silent_sampler", "max_retries", "min_pairs",
    "K", "R", "pair_tag", "progress"
  )
  if (!is.list(predictive) || !all(required %in% names(predictive)))
    stop("Predictive context is incomplete for fold-level scheduling.")
  if (!is.data.frame(predictive$pair_in) ||
      !"subject" %in% names(predictive$pair_in))
    stop("Predictive pair input lacks subject identities.")

  K <- as.integer(predictive$K)
  R <- as.integer(predictive$R)
  if (length(K) != 1L || is.na(K) || K < 1L ||
      length(R) != 1L || is.na(R) || R < 1L)
    stop("Predictive K and R must be positive integers.")

  sample_args <- predictive$sample_args
  sample_args$seed <- as.integer(predictive$sampling_seed)
  sample_args$step_size <- NULL
  sample_args$inv_metric <- NULL
  sample_args$metric_file <- NULL

  splits <- pclvbayes:::.make_repkfold_splits(
    predictive$pair_in$subject, K, R, predictive$split_seed)
  rows <- vector("list", K * R)
  tasks <- vector("list", K * R)
  task_index <- 0L
  base_seed <- if (!is.null(sample_args$seed))
    as.integer(sample_args$seed) else as.integer(predictive$split_seed)
  if (is.na(base_seed)) base_seed <- 1L

  for (r in seq_len(R)) {
    folds <- splits[[r]]
    for (k in seq_len(K)) {
      task_index <- task_index + 1L
      train_subjects <- folds[[k]]$train_subjects
      test_subjects <- folds[[k]]$test_subjects
      fold_identity <- sprintf(
        "%s-r%02d-k%02d", direction_id, as.integer(r), as.integer(k))
      rows[[task_index]] <- data.frame(
        r = as.integer(r), k = as.integer(k),
        test_subjects = I(list(test_subjects))
      )
      tasks[[task_index]] <- list(
        position = as.integer(position),
        direction_id = direction_id,
        fold_identity = fold_identity,
        repetition = as.integer(r), fold = as.integer(k),
        train_subjects = train_subjects, test_subjects = test_subjects,
        seed = as.integer(base_seed + 1000L * r + k),
        mod = predictive$mod,
        stan_list = predictive$stan_list,
        sample_args = sample_args,
        pair_in = predictive$pair_in,
        silent_sampler = predictive$silent_sampler,
        max_retries = predictive$max_retries,
        min_pairs = predictive$min_pairs
      )
    }
  }

  splits_df <- do.call(rbind, rows)
  subjects <- sort(
    unique(as.character(predictive$pair_in$subject)), method = "radix")
  counts <- table(unlist(splits_df$test_subjects, use.names = FALSE))
  holdout_counts <- as.integer(counts[subjects])
  holdout_counts[is.na(holdout_counts)] <- 0L
  names(holdout_counts) <- subjects

  list(
    position = as.integer(position), direction_id = direction_id,
    eligible = TRUE, result = result, predictive = predictive,
    splits_df = splits_df, subject_ids = subjects,
    subject_holdout_counts = holdout_counts,
    tasks = tasks, queue_indices = integer()
  )
}

build_v021_predictive_fold_queue <- function(plans) {
  if (!is.list(plans)) stop("Predictive direction plans must be a list.")
  queue <- list()
  updated <- plans
  for (plan_index in seq_along(updated)) {
    plan <- updated[[plan_index]]
    if (!is.list(plan) || is.null(plan$tasks) || is.null(plan$eligible))
      stop("Invalid predictive direction plan.")
    if (!isTRUE(plan$eligible)) {
      updated[[plan_index]]$queue_indices <- integer()
      next
    }
    begin <- length(queue) + 1L
    queue <- c(queue, plan$tasks)
    finish <- length(queue)
    updated[[plan_index]]$queue_indices <- seq.int(begin, finish)
  }
  identities <- vapply(queue, `[[`, character(1), "fold_identity")
  if (length(identities) && (anyNA(identities) || anyDuplicated(identities)))
    stop("Predictive fold queue identities must be unique.")
  list(plans = updated, tasks = queue)
}

run_v021_predictive_fold_task <- function(task) {
  required <- c(
    "repetition", "fold", "train_subjects", "test_subjects", "seed",
    "mod", "stan_list", "sample_args", "pair_in", "silent_sampler",
    "max_retries", "min_pairs"
  )
  if (!is.list(task) || !all(required %in% names(task)))
    stop("Invalid predictive fold task.")

  result <- pclvbayes:::.fold_fit_and_score(
    mod = task$mod,
    stan_list_base = task$stan_list,
    sample_args_base = task$sample_args,
    pair_in = task$pair_in,
    train_subjects = task$train_subjects,
    test_subjects = task$test_subjects,
    silent_sampler = task$silent_sampler,
    max_retries = task$max_retries,
    freeze_retry_hypers = TRUE,
    seed_override = task$seed,
    min_pairs = task$min_pairs
  )
  if (is.null(result)) {
    result <- pclvbayes:::.pclv_failure(
      "kfold_scoring", "fold_evaluation_failed",
      list(predictor = NA_character_))
  }
  if (inherits(result, "pclv_failure")) {
    if (is.null(result$ok) || is.null(result$details)) {
      legacy_details <- result[setdiff(
        names(result), c("ok", "stage", "reason", "details"))]
      result <- pclvbayes:::.pclv_failure(
        result$stage, result$reason, legacy_details)
    }
    if (is.null(result$predictor)) result$predictor <- NA_character_
    result$repetition <- task$repetition
    result$fold <- task$fold
    result$test_subjects <- task$test_subjects
    result$details$predictor <- result$predictor
    result$details$repetition <- task$repetition
    result$details$fold <- task$fold
    result$details$test_subjects <- task$test_subjects
  }
  result
}

.v021_predictive_worker_failure <- function(plan, value) {
  message <- if (inherits(value, "v021_traced_child_error") &&
                 !is.null(value$condition_message)) {
    as.character(value$condition_message)
  } else if (inherits(value, "try-error")) {
    as.character(value)
  } else if (inherits(value, "error")) {
    conditionMessage(value)
  } else {
    "invalid predictive fold worker result"
  }
  pclvbayes:::.pclv_failure(
    "kfold_evaluation", "unexpected_kfold_error",
    list(
      message = message,
      pair_tag = plan$predictive$pair_tag,
      K = plan$predictive$K, R = plan$predictive$R,
      split_seed = plan$predictive$split_seed,
      sampling_seed = plan$predictive$sampling_seed
    )
  )
}

aggregate_v021_predictive_fold_results <- function(plan, fold_results) {
  if (!is.list(plan) || !isTRUE(plan$eligible) || !is.list(fold_results) ||
      length(fold_results) != length(plan$tasks))
    stop("Predictive fold results do not match their direction plan.")

  invalid <- vapply(fold_results, function(value) {
    inherits(value, c("try-error", "v021_traced_child_error", "error"))
  }, logical(1))
  if (any(invalid))
    return(.v021_predictive_worker_failure(
      plan, fold_results[[which(invalid)[[1L]]]]))

  subs <- sort(unique(plan$predictive$pair_in$subject))
  agg <- setNames(numeric(length(subs)), subs)
  cnt <- setNames(integer(length(subs)), subs)
  agg_ppd <- setNames(numeric(length(subs)), subs)
  cnt_ppd <- setNames(integer(length(subs)), subs)
  fail_cnt <- setNames(integer(length(subs)), subs)
  successful_obs <- setNames(integer(length(subs)), subs)
  fold_diag_df <- list()
  fold_failures <- list()
  n_ok <- 0L
  n_fail <- 0L

  for (element in fold_results) {
    if (inherits(element, "pclv_failure")) {
      n_fail <- n_fail + 1L
      fold_failures[[length(fold_failures) + 1L]] <- element
      failed_subjects <- intersect(element$test_subjects, subs)
      fail_cnt[failed_subjects] <- fail_cnt[failed_subjects] + 1L
    } else {
      if (!is.list(element) || is.null(element$elpd) ||
          is.null(element$fold_diag))
        return(.v021_predictive_worker_failure(plan, element))
      n_ok <- n_ok + 1L
      idx <- names(element$elpd)
      agg[idx] <- agg[idx] + element$elpd
      cnt[idx] <- cnt[idx] + 1L
      if (!is.null(element$elpd_ppd)) {
        agg_ppd[idx] <- agg_ppd[idx] + element$elpd_ppd
        cnt_ppd[idx] <- cnt_ppd[idx] + 1L
      }
      if (!is.null(element$n_obs)) {
        obs_idx <- intersect(names(element$n_obs), subs)
        successful_obs[obs_idx] <- successful_obs[obs_idx] +
          as.integer(element$n_obs[obs_idx])
      }
      fold_diag_df[[length(fold_diag_df) + 1L]] <- data.frame(
        n_retries = element$fold_diag$n_retries,
        nu_mean = element$fold_diag$nu_mean,
        ebfmi_min = element$fold_diag$ebfmi_min,
        worst_rhat = element$fold_diag$worst_rhat,
        min_ess_bulk = element$fold_diag$min_ess_bulk,
        treedepth_hits = element$fold_diag$treedepth_hits,
        n_divergent = element$fold_diag$n_divergent
      )
    }
  }

  elpd_subject <- rep(NA_real_, length(subs)); names(elpd_subject) <- subs
  has_elpd <- cnt > 0L
  elpd_subject[has_elpd] <- agg[has_elpd] / cnt[has_elpd]
  elpd_subject_ppd <- rep(NA_real_, length(subs)); names(elpd_subject_ppd) <- subs
  has_ppd <- cnt_ppd > 0L
  elpd_subject_ppd[has_ppd] <- agg_ppd[has_ppd] / cnt_ppd[has_ppd]
  elpd_mean <- if (length(elpd_subject)) {
    value <- mean(elpd_subject[is.finite(elpd_subject)], na.rm = TRUE)
    if (is.finite(value)) value else NA_real_
  } else NA_real_
  elpd_mean_ppd <- if (length(elpd_subject_ppd)) {
    value <- mean(
      elpd_subject_ppd[is.finite(elpd_subject_ppd)], na.rm = TRUE)
    if (is.finite(value)) value else NA_real_
  } else NA_real_

  fold_diag <- if (length(fold_diag_df)) do.call(rbind, fold_diag_df) else
    data.frame(
      n_retries = integer(), nu_mean = double(), ebfmi_min = double(),
      worst_rhat = double(), min_ess_bulk = double(),
      treedepth_hits = integer(), n_divergent = integer()
    )
  retry_total <- if (nrow(fold_diag))
    sum(fold_diag$n_retries, na.rm = TRUE) else 0L
  retry_mean <- if (nrow(fold_diag))
    mean(fold_diag$n_retries, na.rm = TRUE) else NA_real_
  nu_fold_means <- if (nrow(fold_diag)) fold_diag$nu_mean else numeric()

  list(
    type = "repeated-kfold",
    K = as.integer(plan$predictive$K),
    R = as.integer(plan$predictive$R),
    elpd_subject = elpd_subject,
    elpd_subject_ppd = elpd_subject_ppd,
    elpd_mean = elpd_mean,
    elpd_mean_ppd = elpd_mean_ppd,
    elpd_method = "student-t-scale-mixture-kalman-ou-q16",
    n_folds_ok = n_ok,
    n_folds_fail = n_fail,
    splits_df = plan$splits_df,
    subject_ids = plan$subject_ids,
    subject_test_counts = successful_obs,
    subject_holdout_counts = plan$subject_holdout_counts,
    subject_success_counts = cnt,
    subject_failure_counts = fail_cnt,
    total_successful_evaluations = sum(cnt),
    fold_diag = fold_diag,
    retry_total = retry_total,
    retry_mean = retry_mean,
    nu_fold_means = nu_fold_means,
    failures = fold_failures
  )
}

.v021_apply_predictive_failure <- function(main_result, predictive, failure) {
  main_result$kfold <- NULL
  main_result$kfold_mean <- NA_real_
  main_result$kfold_method <- NA_character_
  main_result$kfold_outer_rounds <- 1L
  main_result$kfold_failed <- TRUE
  main_result$kfold_n_folds_ok <- 0L
  main_result$kfold_n_folds_fail <- NA_integer_
  main_result$kfold_subject <- list(NULL)
  main_result$kfold_subject_ppd <- list(NULL)
  main_result$kfold_subject_ids <- list(NULL)
  main_result$kfold_subject_counts <- list(NULL)
  main_result$kfold_subject_success <- list(NULL)
  main_result$kfold_subject_fail <- list(NULL)
  main_result$kfold_success_total <- 0L
  main_result$kfold_failures <- list(failure)
  main_result$kfold_splits <- list(NULL)
  main_result$kfold_seed_used <- as.integer(predictive$split_seed)
  main_result$kfold_K <- as.integer(predictive$K)
  main_result$kfold_R <- as.integer(predictive$R)
  main_result$kfold_sd <- NA_real_
  main_result$kfold_se <- NA_real_
  main_result$kfold_n_subjects <- 0L
  main_result$kfold_retry_total <- NA_integer_
  main_result$kfold_retry_mean <- NA_real_
  main_result$kfold_nu_fold_means <- list(NULL)
  main_result
}

finalize_v021_predictive_direction <- function(plan, fold_results = list()) {
  if (!is.list(plan) || is.null(plan$result) || is.null(plan$eligible))
    stop("Invalid predictive direction plan.")
  if (!isTRUE(plan$eligible)) return(plan$result)

  kfold <- tryCatch(
    aggregate_v021_predictive_fold_results(plan, fold_results),
    error = function(error) pclvbayes:::.pclv_failure(
      "kfold_evaluation", "unexpected_kfold_error",
      list(
        message = conditionMessage(error),
        pair_tag = plan$predictive$pair_tag,
        K = plan$predictive$K, R = plan$predictive$R,
        split_seed = plan$predictive$split_seed,
        sampling_seed = plan$predictive$sampling_seed
      )
    )
  )
  if (inherits(kfold, "pclv_failure"))
    return(.v021_apply_predictive_failure(
      plan$result, plan$predictive, kfold))
  if (!is.list(kfold) || is.null(kfold$elpd_subject) ||
      is.null(kfold$elpd_method) || is.null(kfold$n_folds_ok) ||
      is.null(kfold$n_folds_fail)) {
    failure <- pclvbayes:::.pclv_failure(
      "kfold_evaluation", "invalid_kfold_result",
      list(pair_tag = plan$predictive$pair_tag))
    return(.v021_apply_predictive_failure(
      plan$result, plan$predictive, failure))
  }

  finite_elpd <- is.finite(kfold$elpd_subject)
  n_finite_subjects <- sum(finite_elpd)
  kfold_sd <- if (n_finite_subjects >= 2L)
    stats::sd(kfold$elpd_subject[finite_elpd]) else NA_real_
  kfold_se <- if (n_finite_subjects >= 2L)
    kfold_sd / sqrt(n_finite_subjects) else NA_real_
  result <- plan$result
  result$kfold <- kfold
  result$kfold_mean <- if (is.null(kfold$elpd_mean))
    NA_real_ else as.numeric(kfold$elpd_mean)
  result$kfold_method <- as.character(kfold$elpd_method)
  result$kfold_outer_rounds <- 1L
  result$kfold_failed <- isTRUE(kfold$n_folds_fail > 0L)
  result$kfold_n_folds_ok <- as.integer(kfold$n_folds_ok)
  result$kfold_n_folds_fail <- as.integer(kfold$n_folds_fail)
  result$kfold_subject <- list(kfold$elpd_subject)
  result$kfold_subject_ppd <- list(kfold$elpd_subject_ppd)
  result$kfold_subject_ids <- list(names(kfold$elpd_subject))
  result$kfold_subject_counts <- list(kfold$subject_test_counts)
  result$kfold_subject_success <- list(kfold$subject_success_counts)
  result$kfold_subject_fail <- list(kfold$subject_failure_counts)
  result$kfold_success_total <- as.integer(
    kfold$total_successful_evaluations)
  result$kfold_failures <- list(kfold$failures)
  result$kfold_splits <- list(kfold$splits_df)
  result$kfold_seed_used <- as.integer(plan$predictive$split_seed)
  result$kfold_K <- as.integer(kfold$K)
  result$kfold_R <- as.integer(kfold$R)
  result$kfold_sd <- kfold_sd
  result$kfold_se <- kfold_se
  result$kfold_n_subjects <- as.integer(n_finite_subjects)
  result$kfold_retry_total <- if (is.null(kfold$retry_total))
    NA_integer_ else as.integer(kfold$retry_total)
  result$kfold_retry_mean <- if (is.null(kfold$retry_mean))
    NA_real_ else as.numeric(kfold$retry_mean)
  result$kfold_nu_fold_means <- list(
    if (is.null(kfold$nu_fold_means)) NULL else kfold$nu_fold_means)
  result
}

build_v021_full_config <- function(output_root) {
  if (!is.character(output_root) || length(output_root) != 1L || is.na(output_root) ||
      !nzchar(output_root)) stop("output_root must be one path.")
  list(
    benchmark_schema = v021_full_schema,
    dataset_id = "1",
    seed = 20260802L,
    taxa_count = 100L,
    unordered_pair_count = 4950L,
    directed_task_count = 9900L,
    selection_rule = v021_full_selection_rule,
    chains = 4L,
    main_parallel_chains = 1L,
    kfold_parallel_chains = 1L,
    iter_warmup = 2000L,
    iter_sampling = 2000L,
    nominal_retained_draws = 8000L,
    maximum_simultaneous_fits = 12L,
    maximum_simultaneous_kfold_fits = 12L,
    rolling_window_size = 60L,
    maximum_task_attempts = 2L,
    kfold_seed = 20260802L,
    preprocessing_config_id = "pclv_smoothed_full_composition_closure_v1",
    posterior_config_id = "student_t_irregular_time_ou_4x2000_v1",
    predictive_config_id = "student-t-scale-mixture-kalman-ou-q16_k5_r1_v1",
    run_kfold = TRUE,
    use_pathfinder = FALSE,
    cpu_affinity = "0-15",
    output_root = normalizePath(output_root, mustWork = FALSE)
  )
}

validate_v021_full_config <- function(config) {
  fields <- c(
    "benchmark_schema", "dataset_id", "seed", "taxa_count",
    "unordered_pair_count", "directed_task_count", "selection_rule", "chains",
    "main_parallel_chains", "kfold_parallel_chains",
    "iter_warmup", "iter_sampling", "nominal_retained_draws",
    "maximum_simultaneous_fits", "maximum_simultaneous_kfold_fits",
    "rolling_window_size",
    "maximum_task_attempts", "kfold_seed",
    "preprocessing_config_id", "posterior_config_id", "predictive_config_id",
    "run_kfold", "use_pathfinder", "cpu_affinity", "output_root"
  )
  if (!is.list(config) || !identical(names(config), fields) ||
      !identical(config$benchmark_schema, v021_full_schema) ||
      !identical(config$dataset_id, "1") || !identical(config$seed, 20260802L) ||
      !identical(config$taxa_count, 100L) ||
      !identical(config$unordered_pair_count, 4950L) ||
      !identical(config$directed_task_count, 9900L) ||
      !identical(config$selection_rule, v021_full_selection_rule) ||
      !identical(config$chains, 4L) ||
      !identical(config$main_parallel_chains, 1L) ||
      !identical(config$kfold_parallel_chains, 1L) ||
      !identical(config$iter_warmup, 2000L) ||
      !identical(config$iter_sampling, 2000L) ||
      !identical(config$nominal_retained_draws, 8000L) ||
      !identical(config$maximum_simultaneous_fits, 12L) ||
      !identical(config$maximum_simultaneous_kfold_fits, 12L) ||
      !identical(config$rolling_window_size, 60L) ||
      !identical(config$maximum_task_attempts, 2L) ||
      !identical(config$kfold_seed, 20260802L) ||
      !identical(config$preprocessing_config_id,
                 "pclv_smoothed_full_composition_closure_v1") ||
      !identical(config$posterior_config_id,
                 "student_t_irregular_time_ou_4x2000_v1") ||
      !identical(config$predictive_config_id,
                 "student-t-scale-mixture-kalman-ou-q16_k5_r1_v1") ||
      !identical(config$run_kfold, TRUE) || !identical(config$use_pathfinder, FALSE) ||
      !identical(config$cpu_affinity, "0-15"))
    stop("Invalid V021-06 full benchmark configuration.")
  invisible(TRUE)
}

canonicalize_v021_full_restart_config <- function(config) {
  validate_v021_full_config(config)
  fields <- names(config)
  canonical <- setNames(lapply(fields, function(name) config[[name]]), fields)
  output_root <- canonical$output_root
  if (!is.character(output_root) || length(output_root) != 1L ||
      is.na(output_root) || !nzchar(output_root))
    stop("V021-06 restart configuration has an invalid output_root.")
  canonical$output_root <- normalizePath(output_root, mustWork = FALSE)
  canonical
}


compare_v021_checkpoint_manifest_identity <- function(stored, requested) {
  validate_v021_checkpoint_manifest(stored)
  validate_v021_checkpoint_manifest(requested)

  identity <- c(
    "task_id", "pair_id", "direction_id", "direction_index",
    "target", "source", "seed", "chain_seeds"
  )
  mismatch <- identity[!vapply(
    identity,
    function(name) identical(stored[[name]], requested[[name]]),
    logical(1)
  )]
  if (length(mismatch))
    stop(
      "V021-06 restart manifest identity changed: ",
      paste(mismatch, collapse = ", "),
      "."
    )

  canonical_path <- function(path) {
    if (!is.character(path) || anyNA(path) || any(!nzchar(path)))
      stop("V021-06 restart manifest contains an invalid output_location.")
    parent <- dirname(path)
    canonical_parent <- vapply(
      parent,
      function(directory) normalizePath(
        directory, winslash = "/", mustWork = TRUE
      ),
      character(1)
    )
    unname(file.path(canonical_parent, basename(path)))
  }
  stored_path <- canonical_path(stored$output_location)
  requested_path <- canonical_path(requested$output_location)
  if (!identical(stored_path, requested_path))
    stop("V021-06 restart manifest identity changed: output_location.")

  invisible(TRUE)
}

compare_v021_full_restart_configs <- function(stored, requested) {
  stored <- canonicalize_v021_full_restart_config(stored)
  requested <- canonicalize_v021_full_restart_config(requested)
  fields <- names(stored)
  if (!identical(fields, names(requested)))
    stop("V021-06 restart configuration schema changed.")
  mismatch <- fields[!vapply(
    fields,
    function(name) identical(stored[[name]], requested[[name]]),
    logical(1)
  )]
  if (length(mismatch))
    stop(
      "V021-06 restart configuration changed: ",
      paste(mismatch, collapse = ", "), "."
    )
  invisible(TRUE)
}

build_v021_full_resource_policy <- function(config) {
  validate_v021_full_config(config)
  build_v021_resource_policy(
    main_chains = config$chains,
    main_parallel_chains = config$main_parallel_chains,
    retry_chains = config$chains,
    retry_parallel_chains = config$main_parallel_chains,
    kfold_chains = config$chains,
    kfold_parallel_chains = config$kfold_parallel_chains,
    confirmation_chains = config$chains,
    confirmation_parallel_chains = config$main_parallel_chains,
    proposed_outer_concurrency = config$maximum_simultaneous_fits,
    maximum_concurrent_kfold_fits =
      config$maximum_simultaneous_kfold_fits
  )
}

v021_current_cpu_affinity <- function(path = "/proc/self/status") {
  lines <- readLines(path, warn = FALSE)
  value <- sub("^Cpus_allowed_list:[[:space:]]*", "",
               grep("^Cpus_allowed_list:", lines, value = TRUE))
  if (length(value) != 1L || !nzchar(value)) stop("CPU affinity is unavailable.")
  value
}

select_v021_100_species_dataset <- function(metadata, dataset_exists) {
  required <- c("did", "n_species", "noise", "n_timeseries", "n_timepoints",
                "sampling_scheme")
  if (!is.data.frame(metadata) || !all(required %in% names(metadata)))
    stop("MTIST metadata lacks required truth-free selection fields.")
  forbidden <- names(metadata)[.v021_truth_name(names(metadata))]
  if (length(forbidden)) metadata <- metadata[setdiff(names(metadata), forbidden)]
  eligible <- metadata$n_species == 100L & metadata$noise == 0.01 &
    metadata$sampling_scheme == "even" & metadata$n_timeseries >= 10L &
    metadata$n_timepoints >= 15L & vapply(metadata$did, dataset_exists, logical(1))
  if (!any(eligible)) stop("No eligible local 100-species observational dataset exists.")
  metadata[which(eligible)[order(metadata$did[eligible])][[1L]], required, drop = FALSE]
}

build_v021_full_task_table <- function(taxa, seed, dataset_id = "1") {
  if (!is.character(taxa) || length(taxa) != 100L || anyNA(taxa) ||
      any(!nzchar(taxa)) || anyDuplicated(taxa)) stop("Exactly 100 unique taxa are required.")
  index <- ten_species_direction_index(taxa, seed)
  if (nrow(index) != 9900L || length(unique(index$task_id)) != 4950L ||
      !identical(index$direction_index, seq_len(9900L)))
    stop("Canonical 100-species direction enumeration failed.")
  index$dataset_id <- as.character(dataset_id)
  index <- index[c("dataset_id", "task_id", "direction_index", "target", "source", "seed")]
  rownames(index) <- NULL
  index
}

v021_full_paths <- function(config, timestamp = format(Sys.time(), "%Y%m%d-%H%M%S")) {
  validate_v021_full_config(config)
  root <- config$output_root
  list(
    output_root = root,
    checkpoint_root = file.path(root, "checkpoints"),
    manifest = file.path(root, "manifest.rds"),
    monitor_root = file.path(root, "monitor"),
    status = file.path(root, "status.rds"),
    summary = file.path(root, "summary.rds"),
    pid_file = file.path(root, "v021_full_100_species.pid"),
    log_file = file.path(root, paste0("v021_full_100_species_", timestamp, ".log"))
  )
}

build_v021_full_launch_command <- function(runner, config_path, log_file, pid_file) {
  values <- c(runner, config_path, log_file, pid_file)
  if (anyNA(values) || any(!nzchar(values))) stop("Launch paths must be non-empty.")
  sprintf("PCLV_V021_FULL_CONFIG=%s nohup Rscript %s > %s 2>&1 < /dev/null & echo $! > %s",
          shQuote(config_path), shQuote(runner), shQuote(log_file), shQuote(pid_file))
}

write_v021_full_status <- function(path, state, config, completed = 0L,
                                    running = 0L, failed = 0L, skipped = 0L,
                                    peak_chains = 0L, peak_slots = 0L,
                                    reason = NA_character_) {
  record <- list(
    status_schema = v021_full_status_schema,
    state = state,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    parent_pid = as.integer(Sys.getpid()),
    output_root = config$output_root,
    completed = as.integer(completed), running = as.integer(running),
    failed = as.integer(failed), skipped = as.integer(skipped),
    peak_chains = as.integer(peak_chains), peak_slots = as.integer(peak_slots),
    reason = reason
  )
  .v021_atomic_save_rds(record, path)
  invisible(record)
}

validate_v021_full_restart_states <- function(manifest) {
  plan <- plan_v021_checkpoint_restart(manifest)
  terminal <- plan$manifest$execution_state %in% c("completed", "failed", "skipped")
  if (any(plan$resume_indices %in% which(terminal)))
    stop("Terminal tasks cannot enter the V021-06 resume plan.")
  plan
}

build_v021_full_execution_manifest <- function(config, taxa, provenance) {
  validate_v021_full_config(config)
  tasks <- build_v021_full_task_table(taxa, config$seed, config$dataset_id)
  build_v021_execution_manifest(
    dataset_id = config$dataset_id, taxa_order = taxa,
    task_table = tasks[c("task_id", "direction_index", "target", "source", "seed")],
    public_seed = config$seed, kfold_seed = config$kfold_seed,
    preprocessing_config = list(identity = config$preprocessing_config_id),
    posterior_config = list(
      identity = config$posterior_config_id, chains = config$chains,
      parallel_chains = config$main_parallel_chains,
      iter_warmup = config$iter_warmup, iter_sampling = config$iter_sampling),
    kfold_config = list(
      K = 5L, R = 1L, enabled = config$run_kfold,
      scheduler = "global_direction_repetition_fold_queue",
      scheduling_unit = "fold_fit",
      parallel_chains = config$kfold_parallel_chains,
      maximum_concurrent_fits = config$maximum_simultaneous_kfold_fits),
    predictive_config = list(identity = config$predictive_config_id),
    provenance = provenance)
}


v021_operational_resume_bridge_schema <-
  "v021_operational_resume_bridge_v1"

v021_operational_resume_allowed_paths <- function() {
  c(
    "benchmarks/mtist/run_v021_full_100_species.R",
    "benchmarks/mtist/v021_resource_policy.R",
    "benchmarks/mtist/v021_full_100_species.R",
    "benchmarks/mtist/v021_full_100_species_v5_ctl.sh",
    "tests/testthat/test-v021-full-100-species.R",
    "tests/testthat/test-v021-resource-policy.R",
    "tests/testthat/test-v021-resource-concurrency.R"
  )
}

.v021_git_capture <- function(repository_root, arguments) {
  repository_root <- normalizePath(repository_root, mustWork = TRUE)
  output <- suppressWarnings(system2(
    "git", c("-C", repository_root, arguments),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(output = as.character(output), status = as.integer(status))
}

.v021_one_commit <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !grepl("^[0-9a-f]{40}$", value))
    stop(name, " must be one full lowercase Git commit SHA.")
  value
}

validate_v021_operational_resume_bridge <- function(bridge) {
  fields <- c(
    "bridge_schema", "authorization", "canonical_manifest_hash",
    "canonical_configuration_hash", "canonical_commit", "runtime_commit",
    "requested_manifest_hash", "unchanged_manifest_fields", "changed_files",
    "changed_file_status", "allowed_paths", "created_timestamp", "bridge_hash"
  )
  if (!is.list(bridge) || !identical(names(bridge), fields) ||
      !identical(bridge$bridge_schema, v021_operational_resume_bridge_schema) ||
      !identical(bridge$authorization, "operational_hotfix_only"))
    stop("Invalid V021 operational-resume bridge schema.")
  for (name in c("canonical_manifest_hash", "canonical_configuration_hash",
                 "requested_manifest_hash", "created_timestamp")) {
    if (!is.character(bridge[[name]]) || length(bridge[[name]]) != 1L ||
        is.na(bridge[[name]]) || !nzchar(bridge[[name]]))
      stop("Operational-resume bridge has invalid ", name, ".")
  }
  .v021_one_commit(bridge$canonical_commit, "canonical_commit")
  .v021_one_commit(bridge$runtime_commit, "runtime_commit")
  if (!is.character(bridge$unchanged_manifest_fields) ||
      !length(bridge$unchanged_manifest_fields) ||
      !is.character(bridge$changed_files) || !length(bridge$changed_files) ||
      anyNA(bridge$changed_files) || any(!nzchar(bridge$changed_files)) ||
      anyDuplicated(bridge$changed_files) ||
      !is.character(bridge$changed_file_status) ||
      length(bridge$changed_file_status) != length(bridge$changed_files) ||
      any(!bridge$changed_file_status %in% c("A", "M")) ||
      !is.character(bridge$allowed_paths) || !length(bridge$allowed_paths) ||
      any(!bridge$changed_files %in% bridge$allowed_paths))
    stop("Operational-resume bridge contains invalid changed-file evidence.")
  expected_hash <- v021_sha256(bridge[setdiff(fields, "bridge_hash")])
  if (!identical(bridge$bridge_hash, expected_hash))
    stop("Operational-resume bridge hash mismatch.")
  invisible(TRUE)
}

build_v021_operational_resume_bridge <- function(
    stored, requested, repository_root, runtime_commit,
    created_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    allowed_paths = v021_operational_resume_allowed_paths()) {
  validate_v021_execution_manifest(stored)
  validate_v021_execution_manifest(requested)
  repository_root <- normalizePath(repository_root, mustWork = TRUE)
  runtime_commit <- .v021_one_commit(runtime_commit, "runtime_commit")

  immutable_fields <- c(
    "manifest_schema", "dataset_id", "taxa_order", "unordered_pair_count",
    "directed_task_count", "configuration", "configuration_hash", "tasks"
  )
  immutable_mismatch <- immutable_fields[!vapply(
    immutable_fields,
    function(name) identical(stored[[name]], requested[[name]]),
    logical(1)
  )]
  if (length(immutable_mismatch))
    stop(
      "Operational resume rejected scientific manifest changes: ",
      paste(immutable_mismatch, collapse = ", "), "."
    )

  stored_provenance <- stored$provenance
  requested_provenance <- requested$provenance
  if (!is.list(stored_provenance) || !is.list(requested_provenance) ||
      is.null(stored_provenance$code_commit) ||
      is.null(requested_provenance$code_commit))
    stop("Operational resume requires code_commit provenance.")
  canonical_commit <- .v021_one_commit(
    stored_provenance$code_commit, "canonical provenance code_commit")
  requested_commit <- .v021_one_commit(
    requested_provenance$code_commit, "requested provenance code_commit")
  if (!identical(requested_commit, runtime_commit))
    stop("Requested manifest commit does not match the runtime commit.")
  stored_provenance$code_commit <- NULL
  requested_provenance$code_commit <- NULL
  if (!identical(stored_provenance, requested_provenance))
    stop("Operational resume rejected non-commit provenance changes.")

  head <- .v021_git_capture(repository_root, c("rev-parse", "HEAD"))
  if (head$status != 0L || length(head$output) != 1L ||
      !identical(trimws(head$output[[1L]]), runtime_commit))
    stop("Runtime commit is not the repository HEAD.")
  tracked_status <- .v021_git_capture(
    repository_root, c("status", "--short", "--untracked-files=no"))
  if (tracked_status$status != 0L || length(tracked_status$output))
    stop("Operational resume requires a clean tracked working tree.")

  ancestor <- .v021_git_capture(
    repository_root,
    c("merge-base", "--is-ancestor", canonical_commit, runtime_commit)
  )
  if (ancestor$status != 0L)
    stop("Canonical manifest commit is not an ancestor of runtime HEAD.")

  diff <- .v021_git_capture(
    repository_root,
    c("diff", "--name-status", "--no-renames",
      paste0(canonical_commit, "..", runtime_commit), "--")
  )
  if (diff$status != 0L)
    stop("Could not inspect operational-hotfix changed files.")
  lines <- diff$output[nzchar(diff$output)]
  if (!length(lines))
    stop("Operational resume requested without an auditable code difference.")
  split <- strsplit(lines, "\t", fixed = TRUE)
  if (any(lengths(split) != 2L))
    stop("Operational-hotfix Git diff has an unsupported name-status shape.")
  statuses <- vapply(split, `[[`, character(1), 1L)
  changed_files <- vapply(split, `[[`, character(1), 2L)
  if (any(!statuses %in% c("A", "M")))
    stop("Operational resume rejects deleted, renamed, or type-changed files.")
  if (anyDuplicated(changed_files) || any(!changed_files %in% allowed_paths)) {
    rejected <- setdiff(changed_files, allowed_paths)
    stop(
      "Operational resume rejected non-allowlisted changes: ",
      paste(rejected, collapse = ", "), "."
    )
  }

  body <- list(
    bridge_schema = v021_operational_resume_bridge_schema,
    authorization = "operational_hotfix_only",
    canonical_manifest_hash = stored$manifest_hash,
    canonical_configuration_hash = stored$configuration_hash,
    canonical_commit = canonical_commit,
    runtime_commit = runtime_commit,
    requested_manifest_hash = requested$manifest_hash,
    unchanged_manifest_fields = immutable_fields,
    changed_files = changed_files,
    changed_file_status = statuses,
    allowed_paths = allowed_paths,
    created_timestamp = as.character(created_timestamp)
  )
  bridge <- c(body, list(bridge_hash = v021_sha256(body)))
  validate_v021_operational_resume_bridge(bridge)
  bridge
}

write_v021_operational_resume_bridge <- function(bridge, result_root) {
  validate_v021_operational_resume_bridge(bridge)
  result_root <- normalizePath(result_root, mustWork = TRUE)
  bridge_root <- file.path(result_root, "operational_resume_bridges")
  dir.create(bridge_root, recursive = TRUE, showWarnings = FALSE)
  stamp <- gsub("[^0-9]", "", bridge$created_timestamp)
  if (!nzchar(stamp)) stamp <- "timestamp-unavailable"
  path <- file.path(
    bridge_root,
    sprintf(
      "bridge-%s-%s-%s.rds", stamp,
      substr(bridge$runtime_commit, 1L, 12L),
      substr(bridge$bridge_hash, 1L, 12L)
    )
  )
  if (file.exists(path)) {
    existing <- tryCatch(readRDS(path), error = identity)
    if (inherits(existing, "error"))
      stop("Existing operational-resume bridge is unreadable.")
    validate_v021_operational_resume_bridge(existing)
    if (!identical(existing$bridge_hash, bridge$bridge_hash))
      stop("Operational-resume bridge path conflicts with existing evidence.")
    return(invisible(path))
  }
  .v021_atomic_save_rds(bridge, path)
  invisible(path)
}

prepare_v021_full_execution <- function(
    config, taxa, provenance, initialize = FALSE,
    operational_resume = FALSE, repository_root = NULL,
    runtime_commit = NULL,
    bridge_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  requested_manifest <- build_v021_full_execution_manifest(
    config, taxa, provenance)
  manifest <- requested_manifest
  paths <- v021_execution_paths(config$output_root)
  bridge <- NULL

  if (file.exists(paths$manifest)) {
    stored <- read_v021_execution_manifest(paths$manifest)
    if (isTRUE(operational_resume) &&
        !identical(stored$manifest_hash, requested_manifest$manifest_hash)) {
      if (is.null(repository_root) || is.null(runtime_commit))
        stop("Operational resume requires repository_root and runtime_commit.")
      bridge <- build_v021_operational_resume_bridge(
        stored = stored,
        requested = requested_manifest,
        repository_root = repository_root,
        runtime_commit = runtime_commit,
        created_timestamp = bridge_timestamp
      )
      manifest <- stored
    } else {
      compare_v021_execution_manifests(stored, requested_manifest)
      manifest <- stored
    }
  } else if (isTRUE(operational_resume)) {
    stop("Operational resume requires an existing canonical result root.")
  }

  if (isTRUE(initialize))
    initialize_v021_execution_root(manifest, config$output_root)
  if (file.exists(paths$manifest)) {
    plan <- plan_v021_execution_resume(
      manifest, config$output_root, config$maximum_task_attempts,
      persist_reconciliation = FALSE)
  } else {
    plan <- list(
      manifest = manifest, status = new_v021_task_status(manifest),
      attempt_ledger = new_v021_attempt_ledger(manifest),
      runnable_indices = seq_len(nrow(manifest$tasks)), completed_indices = integer(),
      reconciliation_changed = FALSE, paths = paths)
  }
  list(
    manifest = manifest,
    requested_manifest = requested_manifest,
    operational_resume_bridge = bridge,
    plan = plan,
    status_summary = compact_v021_execution_status(plan), paths = paths,
    dry_run = !isTRUE(initialize), sampling_launched = FALSE
  )
}

prepare_v021_resource_dry_run <- function(prepared, policy) {
  if (!is.list(prepared) || is.null(prepared$manifest) || is.null(prepared$plan))
    stop("prepared execution is malformed.")
  validate_v021_execution_manifest(prepared$manifest)
  validate_v021_resource_policy(policy)
  ids <- prepared$manifest$tasks$directed_task_id[prepared$plan$runnable_indices]
  list(
    policy_schema = policy$policy_schema,
    policy_hash = v021_resource_policy_hash(policy),
    waves = if (length(ids)) plan_v021_scheduler_waves(ids, policy) else list(),
    maximum_concurrent_tasks = policy$proposed_outer_concurrency,
    maximum_active_cmdstan_chains = policy$maximum_active_cmdstan_chains,
    sampling_launched = FALSE)
}

build_v021_full_storage_projection <- function(
    free_bytes, directions = 9900L, kfold_eligible_fraction = 0.5,
    retry_fraction = 0.1, safety_factor = 1.25) {
  if (!is.numeric(free_bytes) || length(free_bytes) != 1L ||
      !is.finite(free_bytes) || free_bytes <= 0) stop("free_bytes must be positive.")
  categories <- data.frame(
    category = c("cmdstan_csv_four_chain", "profile_diagnostic_sidecars",
      "task_logs", "process_monitor_records", "manifest_status_attempt_ledger",
      "completion_posterior_diagnostic", "robustness_features",
      "direction_checkpoints", "subject_level_kfold", "retry_overhead"),
    retained = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    central_bytes_per_direction = c(16, 1, .25, 1, .02, .05, .02, .02, .10, 0) * 1024^2,
    conservative_bytes_per_direction = c(32, 2, .5, 2, .04, .10, .04, .04, .25, 0) * 1024^2,
    stringsAsFactors = FALSE)
  base_central <- sum(categories$central_bytes_per_direction) * directions
  base_conservative <- sum(categories$conservative_bytes_per_direction) * directions
  kfold_central <- directions * kfold_eligible_fraction * 0.10 * 1024^2
  kfold_conservative <- directions * kfold_eligible_fraction * 0.25 * 1024^2
  retry_central <- base_central * retry_fraction
  retry_conservative <- base_conservative * retry_fraction
  permanent_central <- base_central + kfold_central + retry_central
  permanent_conservative <- base_conservative + kfold_conservative + retry_conservative
  temporary_peak <- 12 * 32 * 1024^2
  safety_requirement <- (permanent_conservative + temporary_peak) * safety_factor
  list(
    projection_schema = "v021_full_storage_projection_v1",
    categories = categories,
    directions = as.integer(directions),
    kfold_eligible_fraction = kfold_eligible_fraction,
    retry_fraction = retry_fraction,
    temporary_peak_bytes = temporary_peak,
    permanent_central_bytes = permanent_central,
    permanent_conservative_bytes = permanent_conservative,
    safety_factor = safety_factor,
    safety_adjusted_requirement_bytes = safety_requirement,
    free_bytes = free_bytes,
    sufficient = free_bytes >= safety_requirement)
}

run_v021_full_dry_run_audit <- function(prepared, policy, result_root,
                                         maximum_attempts = 2L) {
  validate_v021_execution_manifest(prepared$manifest)
  initialize_v021_execution_root(prepared$manifest, result_root)
  owner <- acquire_v021_controller_ownership(
    result_root, prepared$manifest$manifest_hash,
    prepared$manifest$configuration_hash)
  released <- FALSE
  on.exit(if (!released) release_v021_controller_ownership(result_root, owner), add = TRUE)
  plan <- plan_v021_execution_resume(
    prepared$manifest, result_root, maximum_attempts, persist_reconciliation = TRUE)
  first <- start_v021_task_attempt(
    plan$status, plan$attempt_ledger, prepared$manifest, 1L, "dry-t0", "dry-run")
  failed <- finish_v021_task_attempt(
    first$status, first$ledger, prepared$manifest, 1L, "failed", "dry-t1",
    terminal_reason = "injected_dry_run_failure", failure_class = "dry_run")
  write_v021_task_status_atomic(failed$status, plan$paths$status, prepared$manifest)
  write_v021_attempt_ledger_atomic(failed$ledger, plan$paths$attempts, prepared$manifest)
  retry_plan <- plan_v021_execution_resume(
    prepared$manifest, result_root, maximum_attempts, persist_reconciliation = TRUE)
  retry <- start_v021_task_attempt(
    retry_plan$status, retry_plan$attempt_ledger, prepared$manifest, 1L,
    "dry-t2", "dry-run-retry")
  artifact <- build_v021_completion_artifact(
    prepared$manifest, 1L, result_root,
    posterior_summary = list(mean = 0), diagnostics = list(class = "dry_run"),
    psp_lfsr = list(PSP = NA_real_, LFSR = NA_real_),
    predictive_eligibility = FALSE,
    subject_elpd = list(state = "not_applicable",
                        n_successful_test_observations = 0L),
    seed_split_provenance = list(
      direction_seed = prepared$manifest$tasks$direction_seed[[1L]],
      kfold_seed = prepared$manifest$configuration$kfold_seed))
  artifact_path <- file.path(retry_plan$paths$artifact_root, "dry-run-direction-000001.rds")
  write_v021_completion_artifact_atomic(
    artifact, artifact_path, prepared$manifest, result_root)
  completed <- finish_v021_task_attempt(
    retry$status, retry$ledger, prepared$manifest, 1L, "completed", "dry-t3",
    artifact_path = artifact_path, result_root = result_root)
  write_v021_task_status_atomic(completed$status, plan$paths$status, prepared$manifest)
  write_v021_attempt_ledger_atomic(completed$ledger, plan$paths$attempts, prepared$manifest)
  final_plan <- plan_v021_execution_resume(
    prepared$manifest, result_root, maximum_attempts, persist_reconciliation = TRUE)
  waves <- plan_v021_scheduler_waves(
    prepared$manifest$tasks$directed_task_id[final_plan$runnable_indices], policy)
  release_v021_controller_ownership(result_root, owner)
  released <- TRUE
  list(
    audit_schema = "v021_full_dry_run_audit_v1",
    task_count = nrow(prepared$manifest$tasks),
    pair_count = length(unique(prepared$manifest$tasks$pair_id)),
    first_attempt_failed = TRUE, second_attempt_completed = TRUE,
    completed_task_skipped = !1L %in% final_plan$runnable_indices,
    retry_count = completed$status$tasks$attempt_count[[1L]],
    maximum_wave_slots = max(vapply(waves, `[[`, integer(1), "reserved_chain_slots")),
    ownership_released = !dir.exists(v021_ownership_path(result_root)),
    kfold_planning = list(
      priority = "canonical completed-main direction order",
      scheduler = "kfold_global_fold_slot_rolling",
      parallel_chains = policy$kfold_parallel_chains,
      scheduling_unit = "direction_repetition_fold",
      maximum_concurrent_fits = policy$maximum_concurrent_kfold_fits,
      shared_budget = TRUE, main_kfold_overlap = FALSE),
    obsolete_weight_fields_absent = TRUE,
    sampling_launched = FALSE)
}

v021_full_observed_support_input <- function(runtime_context, task) {
  if (!is.list(runtime_context) ||
      !all(c("meta_df", "sm_mat") %in% names(runtime_context)) ||
      !is.data.frame(task) || nrow(task) != 1L)
    stop("Invalid V021 observed-support runtime input.")
  meta <- runtime_context$meta_df
  sm <- runtime_context$sm_mat
  sample_order <- match(as.character(meta$Sample), colnames(sm))
  if (anyNA(sample_order) || anyDuplicated(sample_order) ||
      !all(c(task$source[[1L]], task$target[[1L]]) %in% rownames(sm)))
    stop("Observed-support sample or taxon identity mismatch.")
  source <- as.numeric(sm[task$source[[1L]], sample_order])
  target <- as.numeric(sm[task$target[[1L]], sample_order])
  rest <- 1 - source - target
  data.frame(
    subject = as.character(meta$subject), time = as.numeric(meta$time),
    source_abundance = source, target_abundance = target,
    rest_abundance = rest,
    eligible = is.finite(source) & is.finite(target) & is.finite(rest) & rest >= 0,
    stringsAsFactors = FALSE)
}

.v021_trace_sensitive_name <- function(x) {
  grepl("truth|study|draw|posterior|sample_matrix|coefficient_matrix",
        x, ignore.case = TRUE)
}

.v021_safe_object_metadata <- function(frame) {
  object_names <- ls(frame, all.names = TRUE)
  object_names <- object_names[!.v021_trace_sensitive_name(object_names)]
  if (length(object_names)) {
    lazy <- tryCatch(rlang::env_binding_are_lazy(frame, object_names),
                     error = function(e) rep(TRUE, length(object_names)))
    active <- tryCatch(rlang::env_binding_are_active(frame, object_names),
                       error = function(e) rep(TRUE, length(object_names)))
    object_names <- object_names[!lazy & !active]
  }
  rows <- lapply(object_names, function(object_name) {
    value <- tryCatch(get(object_name, envir = frame, inherits = FALSE),
                      error = identity)
    if (inherits(value, "error") || is.environment(value) || is.function(value) ||
        typeof(value) == "externalptr") return(NULL)
    value_names <- attr(value, "names", exact = TRUE)
    if (!is.null(value_names))
      value_names <- value_names[!.v021_trace_sensitive_name(value_names)]
    list(
      object_name = object_name,
      typeof = typeof(value),
      class = as.character(attr(value, "class", exact = TRUE) %||% typeof(value)),
      length = as.integer(length(value)),
      dimensions = as.integer(attr(value, "dim", exact = TRUE) %||% integer()),
      names = head(as.character(value_names %||% character()), 100L)
    )
  })
  Filter(Negate(is.null), rows)
}

new_v021_failure_trace_context <- function(trace_root, origin,
                                            batch_id = NA_character_,
                                            task_ids = integer(),
                                            direction_ids = character()) {
  if (!is.character(trace_root) || length(trace_root) != 1L || is.na(trace_root) ||
      !nzchar(trace_root)) stop("trace_root must be one path.")
  if (!origin %in% c("parent", "outer_worker", "monitor_child"))
    stop("Invalid failure-trace origin.")
  context <- new.env(parent = emptyenv())
  context$trace_root <- trace_root
  context$origin <- origin
  context$batch_id <- as.character(batch_id)
  context$task_ids <- as.integer(task_ids)
  context$direction_ids <- as.character(direction_ids)
  context$phase <- "initialization"
  context$last_completed_phase <- NA_character_
  context$last_trace <- NULL
  context$cleanup_audit_path <- NA_character_
  context$monitor_state <- NULL
  context
}

set_v021_failure_trace_phase <- function(context, phase, completed = FALSE) {
  if (!is.environment(context) || !is.character(phase) || length(phase) != 1L ||
      is.na(phase) || !nzchar(phase)) stop("Invalid failure-trace phase.")
  context$phase <- phase
  if (isTRUE(completed)) context$last_completed_phase <- phase
  invisible(context)
}

.v021_current_ppid <- function() {
  stat <- tryCatch(readLines("/proc/self/stat", warn = FALSE, n = 1L),
                   error = function(e) character())
  if (!length(stat)) return(NA_integer_)
  parsed <- regexec("^[0-9]+ \\((.*)\\) [^ ] ([0-9]+) ", stat)
  fields <- regmatches(stat, parsed)[[1L]]
  if (length(fields) != 3L) NA_integer_ else as.integer(fields[[3L]])
}

capture_v021_failure_trace <- function(condition, context,
                                        writer = .v021_atomic_save_rds,
                                        calls = sys.calls(), frames = sys.frames(),
                                        timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  if (!inherits(condition, "condition") || !is.environment(context))
    stop("Invalid failure-trace capture input.")
  call_text <- vapply(calls, function(x) paste(deparse(x, width.cutoff = 500L),
                                                collapse = " "), character(1))
  frame_metadata <- lapply(seq_along(frames), function(i) list(
    frame_index = as.integer(i),
    call = if (i <= length(call_text)) call_text[[i]] else NA_character_,
    objects = .v021_safe_object_metadata(frames[[i]])
  ))
  condition_call <- conditionCall(condition)
  trace <- list(
    trace_schema = v021_full_failure_trace_schema,
    timestamp = as.character(timestamp),
    origin = context$origin,
    pid = as.integer(Sys.getpid()),
    ppid = .v021_current_ppid(),
    execution_phase = context$phase,
    last_completed_phase = context$last_completed_phase,
    batch_id = context$batch_id,
    task_ids = context$task_ids,
    direction_ids = context$direction_ids,
    cleanup_audit_path = context$cleanup_audit_path,
    monitor_state = context$monitor_state,
    condition_class = class(condition),
    condition_message = conditionMessage(condition),
    condition_call = if (is.null(condition_call)) NA_character_ else
      paste(deparse(condition_call, width.cutoff = 500L), collapse = " "),
    calls = call_text,
    frame_metadata = frame_metadata
  )
  stamp <- gsub("[^0-9]", "", as.character(timestamp))
  path <- file.path(context$trace_root, sprintf(
    "%s-%s-pid-%d-%s.rds", v021_full_failure_trace_schema, context$origin,
    Sys.getpid(), stamp))
  persistence_error <- tryCatch({
    writer(trace, path)
    NA_character_
  }, error = conditionMessage)
  result <- list(trace = trace,
                 trace_path = if (is.na(persistence_error)) path else NA_character_,
                 trace_persistence_error = persistence_error)
  context$last_trace <- result
  result
}

with_v021_failure_tracing <- function(code, context,
                                       writer = .v021_atomic_save_rds) {
  if (!is.function(code)) stop("code must be a function.")
  handler <- function(condition) {
    if (is.null(context$last_trace)) {
      capture_error <- tryCatch({
        capture_v021_failure_trace(condition, context, writer = writer)
        NA_character_
      }, error = conditionMessage)
      if (!is.na(capture_error)) context$last_trace <- list(
        trace = list(
          trace_schema = v021_full_failure_trace_schema,
          timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
          origin = context$origin,
          condition_class = class(condition),
          condition_message = conditionMessage(condition),
          condition_call = if (is.null(conditionCall(condition))) NA_character_ else
            paste(deparse(conditionCall(condition)), collapse = " ")
        ),
        trace_path = NA_character_,
        trace_persistence_error = paste("trace_capture_failed", capture_error, sep = ": ")
      )
    }
  }
  withCallingHandlers(code(), error = handler, interrupt = handler)
}

run_v021_traced_child <- function(code, context,
                                   writer = .v021_atomic_save_rds) {
  tryCatch(
    with_v021_failure_tracing(code, context, writer),
    error = function(condition) structure(list(
      condition_class = class(condition),
      condition_message = conditionMessage(condition),
      condition_call = if (is.null(conditionCall(condition))) NA_character_ else
        paste(deparse(conditionCall(condition), width.cutoff = 500L), collapse = " "),
      trace = context$last_trace
    ), class = c("v021_traced_child_error", "list"))
  )
}

record_and_resignal_v021_condition <- function(condition, recorder) {
  if (!inherits(condition, "condition") || !is.function(recorder))
    stop("Invalid condition re-signal input.")
  recorder(condition)
  stop(condition)
}

persist_v021_failure_payload <- function(payload, path, trace_root,
                                          writer = .v021_atomic_save_rds) {
  if (file.exists(path)) {
    archive_error <- tryCatch({
      prior <- readRDS(path)
      archive <- file.path(trace_root, sprintf(
        "prior-failure-payload-%s.rds",
        format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC")))
      writer(prior, archive)
      NA_character_
    }, error = conditionMessage)
    if (!is.na(archive_error)) return(list(
      persisted = FALSE,
      persistence_error = paste("prior_failure_archive_failed", archive_error, sep = ": ")
    ))
  }
  persistence_error <- tryCatch({
    writer(payload, path)
    NA_character_
  }, error = conditionMessage)
  list(persisted = is.na(persistence_error), persistence_error = persistence_error)
}

new_v021_batch_ownership <- function(batch_id, parent_pid, worker_registry,
                                      known_model_executable, output_root,
                                      cleanup_audit_path) {
  validate_v021_worker_registry(worker_registry)
  values <- c(batch_id, known_model_executable, output_root, cleanup_audit_path)
  if (anyNA(values) || any(!nzchar(values))) stop("Invalid V021-06 batch ownership paths.")
  ownership <- new.env(parent = emptyenv())
  ownership$batch_id <- batch_id
  ownership$parent_pid <- as.integer(parent_pid)
  ownership$worker_registry <- worker_registry
  ownership$known_model_executable <- normalizePath(known_model_executable, mustWork = FALSE)
  ownership$output_root <- normalizePath(output_root, mustWork = FALSE)
  ownership$cleanup_audit_path <- cleanup_audit_path
  ownership$cleanup_started <- FALSE
  ownership$cleanup_result <- NULL
  ownership
}

.v021_owned_batch_rows <- function(records, ownership) {
  required <- c("pid", "ppid", "start_time", "classification", "executable")
  if (!is.data.frame(records) || !all(required %in% names(records)) ||
      anyDuplicated(records$pid)) stop("Invalid cleanup process inventory.")
  registry <- ownership$worker_registry
  registry_match <- match(records$pid, registry$pid)
  registered <- !is.na(registry_match)
  verified_registered <- registered
  verified_registered[registered] <-
    records$start_time[registered] == registry$start_time[registry_match[registered]] &
    records$ppid[registered] == registry$expected_ppid[registry_match[registered]] &
    records$classification[registered] == "registered_outer_worker"
  verified_worker_pids <- records$pid[verified_registered]
  ancestry_owned <- vapply(records$pid, function(pid) {
    ancestry <- .v021_process_ancestry(records, pid, ownership$parent_pid)
    any(verified_worker_pids %in% ancestry)
  }, logical(1))
  known_descendant <- records$classification %in%
    c("cmdstan_chain", "pathfinder_process", "cmdstan_diagnostic")
  exact_model <- normalizePath(records$executable, mustWork = FALSE) ==
    ownership$known_model_executable
  owned <- verified_registered | (ancestry_owned & known_descendant &
    (exact_model | records$classification %in% c("pathfinder_process", "cmdstan_diagnostic")))
  list(
    owned = records[owned, , drop = FALSE],
    excluded = records[!owned & records$pid != ownership$parent_pid, , drop = FALSE],
    registry_mismatch = records[registered & !verified_registered, , drop = FALSE]
  )
}

cleanup_v021_owned_batch <- function(
    ownership, trigger_condition = NULL, origin = "parent", phase = NA_character_,
    inventory_reader, signaler, alive_reader, reaper = function() list(),
    sleeper = Sys.sleep, wait_attempts = 5L, wait_seconds = 0.1,
    writer = .v021_atomic_save_rds,
    clock = function() format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  if (!is.environment(ownership)) stop("Invalid V021-06 batch ownership record.")
  if (isTRUE(ownership$cleanup_started)) return(ownership$cleanup_result)
  ownership$cleanup_started <- TRUE
  original_message <- if (inherits(trigger_condition, "condition"))
    conditionMessage(trigger_condition) else as.character(trigger_condition %||% "normal_return")
  audit <- list(
    cleanup_schema = v021_full_cleanup_audit_schema,
    batch_id = ownership$batch_id,
    trigger_condition = original_message,
    origin = origin,
    phase = as.character(phase),
    parent_pid = ownership$parent_pid,
    start_timestamp = clock(),
    output_root = ownership$output_root,
    owned_processes = data.frame(),
    excluded_processes = data.frame(),
    registry_mismatches = data.frame(),
    signals = list(),
    reap_result = NULL,
    survivors = integer(),
    sigkill_used = FALSE,
    sigkill_limitation = "SIGKILL and uncatchable parent termination cannot guarantee R-level tracing.",
    cleanup_error = NA_character_,
    audit_persistence_error = NA_character_,
    end_timestamp = NA_character_
  )
  cleanup_error <- tryCatch({
    inventory <- inventory_reader()
    selected <- .v021_owned_batch_rows(inventory, ownership)
    audit$owned_processes <- selected$owned
    audit$excluded_processes <- selected$excluded
    audit$registry_mismatches <- selected$registry_mismatch
    if (nrow(selected$registry_mismatch))
      stop("Registered worker identity changed during cleanup.")
    owned_identity <- selected$owned[rev(seq_len(nrow(selected$owned))), , drop = FALSE]
    owned_pids <- owned_identity$pid
    owned_start_times <- owned_identity$start_time
    live <- alive_reader(owned_pids, owned_start_times)
    if (length(live)) {
      audit$signals[[length(audit$signals) + 1L]] <- list(
        signal = "SIGINT", pids = as.integer(live), timestamp = clock())
      signaler(live, "SIGINT")
    }
    for (attempt in seq_len(as.integer(wait_attempts))) {
      live <- alive_reader(owned_pids, owned_start_times)
      if (!length(live)) break
      sleeper(wait_seconds)
    }
    live <- alive_reader(owned_pids, owned_start_times)
    if (length(live)) {
      audit$signals[[length(audit$signals) + 1L]] <- list(
        signal = "SIGTERM", pids = as.integer(live), timestamp = clock())
      signaler(live, "SIGTERM")
      for (attempt in seq_len(as.integer(wait_attempts))) {
        live <- alive_reader(owned_pids, owned_start_times)
        if (!length(live)) break
        sleeper(wait_seconds)
      }
    }
    audit$survivors <- as.integer(alive_reader(owned_pids, owned_start_times))
    NA_character_
  }, error = conditionMessage)
  audit$cleanup_error <- cleanup_error
  audit$reap_result <- tryCatch(reaper(), error = function(e)
    list(error = conditionMessage(e)))
  audit$end_timestamp <- clock()
  persistence_error <- tryCatch({
    writer(audit, ownership$cleanup_audit_path)
    NA_character_
  }, error = conditionMessage)
  audit$audit_persistence_error <- persistence_error
  ownership$cleanup_result <- audit
  audit
}
