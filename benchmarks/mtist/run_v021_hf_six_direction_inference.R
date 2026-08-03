args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
repo <- normalizePath(file.path(script_dir, "../.."))
pkgload::load_all(repo, quiet = TRUE)
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "ten_species_helpers.R"))
source(file.path(script_dir, "v021_truth_isolation.R"))
source(file.path(script_dir, "v021_resource_policy.R"))
source(file.path(script_dir, "v021_checkpoint_manifest.R"))
source(file.path(script_dir, "v021_diagnostic_features.R"))
source(file.path(script_dir, "v021_four_chain_preflight.R"))

output_root <- file.path(script_dir, "results", "v021_hf2_six_direction_revalidation_v1")
manifest_csv <- file.path(output_root, "direction_manifest.csv")
expected_hash <- "9e6907aa47953146f56e08b7e78f3c6f57619e8bc4ee2f31f0190f3d6364c005"
sha256 <- system2("sha256sum", shQuote(manifest_csv), stdout = TRUE)
if (!startsWith(sha256[[1L]], expected_hash)) stop("Frozen manifest hash changed.")
inference_root <- file.path(output_root, "inference")
if (dir.exists(inference_root) && length(list.files(
    inference_root, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, include.dirs = FALSE, no.. = TRUE)))
  stop("Fresh inference root is not empty; refusing checkpoint reuse.")
dir.create(inference_root, showWarnings = FALSE)
checkpoint_root <- file.path(inference_root, "checkpoints"); dir.create(checkpoint_root, showWarnings = FALSE)
artifact_root <- file.path(inference_root, "artifacts"); dir.create(artifact_root, showWarnings = FALSE)

atomic_save <- function(x, path) {
  tmp <- tempfile(paste0(".", basename(path), "-"), dirname(path))
  saveRDS(x, tmp, version = 3)
  if (!file.rename(tmp, path)) stop("Atomic output promotion failed: ", path)
}
write_tsv <- function(x, path) utils::write.table(
  x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

started <- Sys.time()
manifest <- utils::read.csv(manifest_csv, stringsAsFactors = FALSE)
manifest$dataset_id <- "37"
manifest$pair_id <- sprintf("task-%02d", manifest$task_id)
manifest <- manifest[c("dataset_id", "pair_id", "task_id", "direction_index",
                       "target", "source", "seed")]

policy <- build_v021_resource_policy(proposed_outer_concurrency = 1L)
capacity <- validate_v021_preflight_launch_capacity(policy, 1L)
thread_environment <- v021_single_thread_environment()
validate_v021_single_thread_environment(thread_environment)

root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
observations <- load_v021_preflight_observations(root, 37L)
study <- list(physeq = observations$physeq, taxa = observations$taxa)

fit_formals <- formals(fit_pclv_bayes)
control_names <- setdiff(names(fit_formals), c("physeq", "subject_col", "time_col", "taxa_vec"))
evaluation <- new.env(parent = environment(fit_pclv_bayes))
controls <- list()
for (name in control_names) {
  controls[[name]] <- eval(fit_formals[[name]], evaluation)
  assign(name, controls[[name]], evaluation)
}
controls[c(
  "eps", "zero_mode_alr", "minpos_alpha", "minpos_base", "eps_fixed",
  "lib_eps_c", "rest_floor_frac", "smooth_scale", "alr_spline_df",
  "alr_spline_spar", "alr_spline_cv", "metric", "quiet", "progress_every",
  "silent_sampler", "max_retries", "use_pathfinder_init", "pf_num_paths",
  "pf_draws", "pf_history_size", "pf_max_lbfgs_iters", "pf_psis_resample"
)] <- list(
  1e-6, "minpos_time", 0.5, "ij", 1e-6, 0.65, 1.0, "logra", NULL,
  NULL, TRUE, "diag_e", FALSE, 1L, FALSE, 0L, FALSE,
  8L, 1000L, 50L, 200L, TRUE
)
controls[c("chains", "iter_warmup", "iter_sampling", "seed", "kfold_seed", "progress",
           "n_workers_outer", "n_workers_kfold", "kfold_K", "kfold_R")] <- list(
  4L, 2000L, 2000L, 20260802L, 20260802L, "none", 1L, 1L, 5L, 1L)
validated <- pclvbayes:::.validate_fit_pclv_inputs(
  observations$physeq, "subject", "time", observations$taxa, controls)
runtime <- pclvbayes:::.prepare_fit_runtime(validated)
if (inherits(runtime, "pclv_failure")) stop(runtime$reason)
runtime$ctx$n_workers_kfold_eff <- 1L
runtime$ctx$max_retries <- 0L
runtime$ctx$use_pathfinder_init <- FALSE
approved_runtime <- build_v021_runtime_context(runtime$ctx)
# HF2 adds a validated scalar split seed to the production runtime. The older
# benchmark projection predates that field, so restore only this truth-free
# scalar after projecting the otherwise approved context.
approved_runtime$kfold_seed <- as.integer(runtime$ctx$kfold_seed)
if (length(approved_runtime$kfold_seed) != 1L ||
    is.na(approved_runtime$kfold_seed) || approved_runtime$kfold_seed < 0L)
  stop("Validated HF2 K-fold split seed is unavailable.")
executable <- normalizePath(approved_runtime$mod_exe_file, mustWork = TRUE)
executable_before <- unclass(file.info(executable)[c("size", "mtime")])
inference_config <- controls[intersect(names(controls), v021_inference_config_fields)]
specs <- lapply(seq_len(nrow(manifest)), function(i)
  build_v021_inference_spec(study, inference_config, manifest[i, , drop = FALSE]))
rm(study)

checkpoint_manifest <- manifest
checkpoint_manifest$chain_seeds <- I(lapply(checkpoint_manifest$seed, v021_chain_seeds,
                                            chains = 4L))
checkpoint_manifest$execution_state <- "pending"
checkpoint_manifest$elapsed_time <- 0
checkpoint_manifest$failure_reason <- NA_character_
checkpoint_manifest$output_location <- file.path(
  checkpoint_root, sprintf("direction-%06d.rds", checkpoint_manifest$direction_index))
atomic_save(checkpoint_manifest, file.path(inference_root, "manifest.rds"))
atomic_save(list(
  schema = "v021_hf2_six_direction_inference_config_v1",
  manifest_sha256 = expected_hash, chains = 4L, iter_warmup = 2000L,
  iter_sampling = 2000L, nominal_draws = 8000L,
  resource_policy = policy, capacity = capacity,
  thread_environment = thread_environment,
  preprocessing_schema = "pclv_smoothed_full_composition_closure_v1",
  posterior_likelihood = "student_t", residual_model = "irregular_time_ou",
  predictive_scorer = "student-t-scale-mixture-kalman-ou-q16",
  executable = executable, started = format(started, tz = "UTC", usetz = TRUE)
), file.path(inference_root, "config.rds"))

results <- vector("list", nrow(manifest))
states <- vector("list", nrow(manifest))
run_all <- function() {
  for (i in seq_len(nrow(manifest))) {
    task_started <- Sys.time()
    checkpoint_manifest$execution_state[[i]] <<- "running"
    atomic_save(checkpoint_manifest, file.path(inference_root, "manifest.rds"))
    result <- pclvbayes:::.fit_direction_main_posterior(
      target = specs[[i]]$target, partner = specs[[i]]$source,
      ctx = approved_runtime, seed_override = as.integer(specs[[i]]$seed),
      progress_local = "none")
    if (!inherits(result, "pclv_failure") &&
        is.null(result$diagnostic_failure[[1L]])) {
      result <- pclvbayes:::.add_predictive_evaluation(result)
    } else if (!inherits(result, "pclv_failure")) {
      result$.predictive_context <- NULL
    }
    elapsed <- as.numeric(difftime(Sys.time(), task_started, units = "secs"))
    if (inherits(result, "pclv_failure")) {
      state <- "failed"
      reason <- paste(result$stage, result$reason, sep = ":")
      checkpoint <- list(
        schema = "v021_hf_direction_checkpoint_v1",
        identity = manifest[i, , drop = FALSE],
        chain_seeds = checkpoint_manifest$chain_seeds[[i]],
        execution_state = state, elapsed_time = elapsed,
        failure_reason = reason,
        retry_history = result$details$retry_history %||% list(),
        completion_timestamp = NA_character_, artifact_path = NA_character_)
    } else {
      state <- "completed"
      reason <- NA_character_
      subject_elpd <- unlist(result$kfold_subject[[1L]], use.names = TRUE)
      subject_counts <- unlist(result$kfold_subject_counts[[1L]], use.names = TRUE)
      ids <- union(names(subject_elpd), names(subject_counts))
      result$elpd_subject_cross <- tibble::tibble(
        from = rep(specs[[i]]$source, length(ids)),
        to = rep(specs[[i]]$target, length(ids)),
        subject = as.character(ids),
        elpd = as.double(subject_elpd[ids]),
        elpd_per_observation = as.double(subject_elpd[ids] / subject_counts[ids]),
        n_successful_test_observations = as.integer(subject_counts[ids]),
        elpd_method = rep(as.character(result$kfold_method), length(ids)))
      result$.predictive_context <- NULL
      atomic_save(result, file.path(artifact_root,
        sprintf("direction-%06d.rds", manifest$direction_index[[i]])))
      checkpoint <- list(
        schema = "v021_hf_direction_checkpoint_v1",
        identity = manifest[i, , drop = FALSE],
        chain_seeds = checkpoint_manifest$chain_seeds[[i]],
        execution_state = state, elapsed_time = elapsed,
        failure_reason = NA_character_,
        retry_history = result$retry_history[[1L]] %||% list(),
        pathfinder_used = FALSE,
        completion_timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
        artifact_path = file.path(artifact_root,
          sprintf("direction-%06d.rds", manifest$direction_index[[i]])))
    }
    atomic_save(checkpoint, checkpoint_manifest$output_location[[i]])
    checkpoint_manifest$execution_state[[i]] <<- state
    checkpoint_manifest$elapsed_time[[i]] <<- elapsed
    checkpoint_manifest$failure_reason[[i]] <<- reason
    atomic_save(checkpoint_manifest, file.path(inference_root, "manifest.rds"))
    results[[i]] <<- result
    states[[i]] <<- data.frame(
      task_id = manifest$task_id[[i]], direction_index = manifest$direction_index[[i]],
      source = manifest$source[[i]], target = manifest$target[[i]], seed = manifest$seed[[i]],
      execution_state = state, failure_reason = reason, elapsed_seconds = elapsed,
      stringsAsFactors = FALSE)
    cat(sprintf("DIRECTION %d %s->%s %s elapsed=%.2f\n",
      manifest$direction_index[[i]], manifest$source[[i]], manifest$target[[i]],
      state, elapsed))
  }
}
withr::with_envvar(thread_environment, run_all())

state_table <- do.call(rbind, states)
write_tsv(state_table, file.path(inference_root, "terminal_states.tsv"))

split_rows <- list(); seed_rows <- list()
for (i in which(state_table$execution_state == "completed")) {
  result <- results[[i]]
  splits <- result$kfold_splits[[1L]]
  if (is.null(splits) || !nrow(splits)) next
  for (j in seq_len(nrow(splits))) {
    eligible <- sort(unique(as.character(result$subjects)))
    held_out <- as.character(splits$test_subjects[[j]])
    split_rows[[length(split_rows) + 1L]] <- data.frame(
      source = manifest$source[[i]], target = manifest$target[[i]],
      eligible_subject_universe = paste(eligible, collapse = ";"),
      repetition = splits$r[[j]], fold = splits$k[[j]],
      train_subjects = paste(setdiff(eligible, held_out), collapse = ";"),
      test_subjects = paste(held_out, collapse = ";"),
      kfold_seed_used = result$kfold_seed_used,
      direction_sampling_seed = manifest$seed[[i]],
      fold_sampler_seed = as.integer(manifest$seed[[i]] + 1000L * splits$r[[j]] + splits$k[[j]]),
      stringsAsFactors = FALSE)
  }
  seed_rows[[length(seed_rows) + 1L]] <- data.frame(
    source = manifest$source[[i]], target = manifest$target[[i]],
    direction_index = manifest$direction_index[[i]],
    kfold_seed_used = result$kfold_seed_used,
    direction_sampling_seed = manifest$seed[[i]],
    chain_seeds = paste(checkpoint_manifest$chain_seeds[[i]], collapse = ";"),
    stringsAsFactors = FALSE)
}
empty_splits <- data.frame(source=character(), target=character(), eligible_subject_universe=character(),
  repetition=integer(), fold=integer(), train_subjects=character(), test_subjects=character(),
  kfold_seed_used=integer(), direction_sampling_seed=integer(), fold_sampler_seed=integer())
empty_seeds <- data.frame(source=character(), target=character(), direction_index=integer(),
  kfold_seed_used=integer(), direction_sampling_seed=integer(), chain_seeds=character())
write_tsv(if (length(split_rows)) do.call(rbind, split_rows) else empty_splits,
          file.path(output_root, "six_direction_kfold_split_audit.tsv"))
write_tsv(if (length(seed_rows)) do.call(rbind, seed_rows) else empty_seeds,
          file.path(output_root, "six_direction_seed_audit.tsv"))
executable_after <- unclass(file.info(executable)[c("size", "mtime")])
summary <- list(
  schema = "v021_hf2_six_direction_inference_summary_v1",
  states = table(state_table$execution_state),
  executable = executable,
  executable_reused = identical(executable_before, executable_after),
  worker_compilation_count = 0L,
  truth_loaded = FALSE,
  manifest_sha256 = expected_hash,
  started = format(started, tz = "UTC", usetz = TRUE),
  finished = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
atomic_save(summary, file.path(inference_root, "summary.rds"))
cat("INFERENCE_ROOT=", inference_root, "\n", sep = "")
print(state_table, row.names = FALSE)
