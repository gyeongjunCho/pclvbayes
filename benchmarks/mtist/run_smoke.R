args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "score_mtist.R"))
config_path <- Sys.getenv("PCLV_MTIST_CONFIG", unset = file.path(script_dir, "configs", "smoke.R"))
config <- source(config_path)$value
config$n_workers_outer <- config$n_workers_outer %||% config$n_workers %||% 1L
config$n_workers_kfold <- config$n_workers_kfold %||% config$n_workers %||% 1L

root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
study <- load_mtist_study(config$dataset_id, root)
n_subjects <- length(unique(study$sample_metadata$subject))
if (n_subjects < 2L) stop("Smoke benchmark requires at least two subjects.")
config$kfold_K <- min(5L, n_subjects)
config$mtist_root <- root
config$conservative_filter_scope <- "off_diagonal_only"
config$primary_interpretation <- "off_diagonal_diagonal_excluded"

result_dir <- file.path(script_dir, "results", config$result_label %||% as.character(config$dataset_id))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(result_dir, "console.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE); sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) }, add = TRUE)

cat("MTIST pcLVbayes smoke benchmark\n")
cat("dataset:", config$dataset_id, " subjects:", n_subjects, " taxa:", paste(study$taxa, collapse = ","), "\n")
devtools::load_all(quiet = TRUE)
invisible(get_pclv_model(quiet = TRUE))
started <- Sys.time()
fit <- fit_pclv_bayes(
  study$physeq, subject_col = "subject", time_col = "time", taxa_vec = study$taxa,
  chains = config$chains, iter_warmup = config$iter_warmup,
  iter_sampling = config$iter_sampling, seed = config$seed,
  n_workers_outer = config$n_workers_outer, n_workers_kfold = config$n_workers_kfold,
  kfold_K = config$kfold_K, kfold_R = config$kfold_R
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
if (inherits(fit, "pclv_failure") || !is.list(fit) || is.null(fit$raw))
  stop("Top-level fit failed: ", paste(capture.output(str(fit)), collapse = " "))

taxa <- study$taxa; n <- length(taxa)
posterior <- matrix(0, n, n, dimnames = list(taxa, taxa))
conservative <- posterior
direction_rows <- list(); failure_rows <- list(); self_rows <- list(); k <- 0L
failure_text <- function(x, field) if (is.null(x)) NA_character_ else as.character(x[[field]] %||% NA_character_)
for (r in seq_len(nrow(fit$raw))) {
  z <- fit$raw[r, , drop = FALSE]
  specs <- list(
    list(source = z$j[[1]], target = z$i[[1]], ok = z$direction_ok_ij[[1]],
         mean = z$a_ij_mean[[1]], p = z$p_sign2_ij[[1]], self = z$a_ii_mean[[1]], failure = z$failure_ij[[1]],
         retries = max(NROW(z$retry_history_ij[[1]]) - 1L, 0L), rhat = z$rhat_ij[[1]], ess_bulk = z$essb_ij[[1]],
         ess_tail = z$esst_ij[[1]], divergences = z$div_ij[[1]], treedepth_hits = z$tdhit_ij[[1]],
         ebfmi_min = z$ebfmi_min_ij[[1]], kfold_ok = z$kfold_folds_ok_ij[[1]], kfold_failed = z$kfold_folds_fail_ij[[1]]),
    list(source = z$i[[1]], target = z$j[[1]], ok = z$direction_ok_ji[[1]],
         mean = z$a_ji_mean[[1]], p = z$p_sign2_ji[[1]], self = z$a_jj_mean[[1]], failure = z$failure_ji[[1]],
         retries = max(NROW(z$retry_history_ji[[1]]) - 1L, 0L), rhat = z$rhat_ji[[1]], ess_bulk = z$essb_ji[[1]],
         ess_tail = z$esst_ji[[1]], divergences = z$div_ji[[1]], treedepth_hits = z$tdhit_ji[[1]],
         ebfmi_min = z$ebfmi_min_ji[[1]], kfold_ok = z$kfold_folds_ok_ji[[1]], kfold_failed = z$kfold_folds_fail_ji[[1]])
  )
  for (d in specs) {
    k <- k + 1L; ok <- isTRUE(d$ok) && is.finite(d$mean); lfsr <- if (ok) d$p / 2 else NA_real_
    posterior[d$target, d$source] <- if (ok) d$mean else 0
    conservative[d$target, d$source] <- if (ok && is.finite(lfsr) && lfsr <= config$conservative_threshold) d$mean else 0
    direction_rows[[k]] <- data.frame(source = d$source, target = d$target, direction_ok = ok,
      posterior_mean = if (ok) d$mean else NA_real_, p_sign2 = if (ok) d$p else NA_real_,
      lfsr = lfsr, conservative_pass = ok && is.finite(lfsr) && lfsr <= config$conservative_threshold,
      retries = d$retries, rhat = d$rhat, ess_bulk = d$ess_bulk, ess_tail = d$ess_tail,
      divergences = d$divergences, treedepth_hits = d$treedepth_hits, ebfmi_min = d$ebfmi_min,
      kfold_folds_ok = d$kfold_ok, kfold_folds_failed = d$kfold_failed,
      stage = failure_text(d$failure, "stage"), reason = failure_text(d$failure, "reason"))
    failure_rows[[k]] <- data.frame(source = d$source, target = d$target, failed = !ok,
      matrix_value_on_failure = if (!ok) 0 else NA_real_, stage = failure_text(d$failure, "stage"),
      reason = failure_text(d$failure, "reason"))
    self_rows[[k]] <- data.frame(target = d$target, direction_ok = ok,
                                 a_self_mean = if (ok && is.finite(d$self)) d$self else NA_real_)
  }
}
directions <- do.call(rbind, direction_rows); failures <- do.call(rbind, failure_rows)
self_evidence <- do.call(rbind, self_rows)
diagonal <- do.call(rbind, lapply(taxa, function(taxon) {
  values <- self_evidence$a_self_mean[self_evidence$target == taxon & is.finite(self_evidence$a_self_mean)]
  unavailable <- !length(values)
  value <- if (unavailable) 0 else stats::median(values)
  posterior[taxon, taxon] <<- value; conservative[taxon, taxon] <<- value
  data.frame(taxon = taxon, n_pair_specific_means = length(values), aggregated_value = value,
             diagonal_unavailable = unavailable,
             reason = if (unavailable) "no_successful_pair_specific_self_effect_mean" else NA_character_,
             aggregation = "median across pair-specific posterior self-effect means")
}))

validate_mtist_prediction(posterior, taxa); validate_mtist_prediction(conservative, taxa)
if (!identical(rownames(study$truth), taxa) || !identical(colnames(study$truth), taxa)) stop("Truth order mismatch.")
n_success <- sum(directions$direction_ok); n_failed <- nrow(directions) - n_success
scores <- c(
  posterior_mean_es_with_diagonal = official_mtist_es(study$truth, posterior, root, FALSE),
  posterior_mean_es_without_diagonal = official_mtist_es(study$truth, posterior, root, TRUE),
  conservative_es_with_diagonal = official_mtist_es(study$truth, conservative, root, FALSE),
  conservative_es_without_diagonal = official_mtist_es(study$truth, conservative, root, TRUE)
)
metrics <- rbind(sign_metrics(study$truth, posterior, "posterior_mean", n_success, n_failed),
                 sign_metrics(study$truth, conservative, "conservative", n_success, n_failed),
                 data.frame(matrix = "official_es", metric = names(scores), value = as.numeric(scores),
                            numerator = NA_real_, denominator = NA_real_, undefined_reason = NA_character_,
                            definition = "official MTIST calculate_es_score"))

write_matrix <- function(x, name) utils::write.table(x, file.path(result_dir, name), sep = ",",
  row.names = FALSE, col.names = FALSE, quote = FALSE, na = "NA")
write_matrix(posterior, "posterior_mean_matrix.csv"); write_matrix(conservative, "conservative_matrix.csv")
write_matrix(study$truth, "truth_matrix.csv")
utils::write.table(failures, file.path(result_dir, "failure_mask.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
utils::write.table(directions, file.path(result_dir, "direction_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
utils::write.table(diagonal, file.path(result_dir, "diagonal_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
utils::write.table(metrics, file.path(result_dir, "metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

pc_sha <- system2("git", c("-C", shQuote(normalizePath(file.path(script_dir, "../.."))), "rev-parse", "HEAD"), stdout = TRUE)[[1L]]
config$dataset_id <- as.integer(config$dataset_id); config$truth_matrix_id <- study$record$ground_truth[[1L]]
config$taxon_order <- taxa; config$n_subjects <- n_subjects; config$n_timepoints_per_subject <- as.integer(study$record$n_timepoints[[1L]])
config$elapsed_wall_seconds <- elapsed
jsonlite::write_json(config, file.path(result_dir, "config.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
cpu <- tryCatch(readLines("/proc/cpuinfo"), error = function(e) character())
cpu_model <- sub("^[^:]+: *", "", grep("^model name", cpu, value = TRUE)[1L] %||% NA_character_)
environment <- list(timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), pcLVbayes_git_sha = pc_sha,
  mtist_git_sha = mtist_git_sha(root), dataset_id = config$dataset_id, truth_matrix_id = config$truth_matrix_id,
  seed = config$seed, chains = config$chains, iter_warmup = config$iter_warmup,
  iter_sampling = config$iter_sampling, kfold_K = config$kfold_K, kfold_R = config$kfold_R,
  n_workers_outer = config$n_workers_outer, n_workers_kfold_requested = config$n_workers_kfold,
  n_workers_kfold_effective = if (config$n_workers_outer > 1L) 1L else config$n_workers_kfold,
  elapsed_wall_seconds = elapsed, R_version = R.version.string,
  CmdStan_version = as.character(cmdstanr::cmdstan_version()), operating_system = Sys.info()[["sysname"]],
  os_release = Sys.info()[["release"]], cpu_model = cpu_model, logical_core_count = parallel::detectCores(logical = TRUE))
jsonlite::write_json(environment, file.path(result_dir, "environment.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
cat("elapsed_seconds:", elapsed, " successes:", n_success, " failures:", n_failed, "\n")
print(scores); print(metrics)
