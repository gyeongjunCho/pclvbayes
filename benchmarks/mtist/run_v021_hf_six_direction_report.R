args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
repo <- normalizePath(file.path(script_dir, "../.."))
pkgload::load_all(repo, quiet = TRUE)
source(file.path(script_dir, "mtist_adapter.R"))

output_root <- file.path(script_dir, "results", "v021_hf2_six_direction_revalidation_v1")
manifest_path <- file.path(output_root, "direction_manifest.csv")
expected_hash <- "9e6907aa47953146f56e08b7e78f3c6f57619e8bc4ee2f31f0190f3d6364c005"
actual_hash <- strsplit(system2("sha256sum", shQuote(manifest_path), stdout = TRUE)[[1L]], " ")[[1L]][[1L]]
if (!identical(actual_hash, expected_hash)) stop("Frozen manifest hash changed.")
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
states <- utils::read.delim(file.path(output_root, "inference", "terminal_states.tsv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
if (any(states$execution_state %in% c("pending", "running", "incomplete")))
  stop("Every direction must have a terminal state before truth joining.")

write_tsv <- function(x, path) utils::write.table(
  x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
collapse_named <- function(x) {
  x <- unlist(x, use.names = TRUE)
  paste(names(x), format(x, digits = 10, scientific = FALSE), sep = "=", collapse = ";")
}
zero_tolerance <- 1e-10
sign_tol <- function(x) if (!is.finite(x)) NA_integer_ else if (abs(x) <= zero_tolerance) 0L else as.integer(sign(x))

# Absolute MTIST truth is first loaded here, after inference reaches terminal states.
root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
study <- load_mtist_study(37L, root)
truth <- study$truth
oracle_all <- utils::read.delim(file.path(
  output_root, "deterministic_audit", "corrected_oracle_all_90.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE)
oracle <- oracle_all[match(manifest$direction_index, oracle_all$direction_index), , drop = FALSE]
if (anyNA(oracle$direction_index)) stop("Corrected oracle rows are incomplete.")

rows <- vector("list", nrow(manifest))
failures <- list()
for (i in seq_len(nrow(manifest))) {
  d <- manifest[i, ]
  state <- states$execution_state[match(d$direction_index, states$direction_index)]
  A_value <- as.numeric(truth[d$target, d$source])
  A_zero <- identical(A_value, 0)
  A_sign <- sign_tol(A_value)
  oracle_value <- oracle$observed_coefficient[[i]]
  oracle_sign <- sign_tol(oracle_value)
  glv_class <- if (A_zero && oracle_sign != 0L) "absolute_zero_to_projected_nonzero" else
    if (!A_zero && oracle_sign == 0L) "absolute_nonzero_to_projected_zero" else
    if (A_zero && oracle_sign == 0L) "both_zero" else
    if (A_sign == oracle_sign) "same_sign" else "sign_reversal"
  state_class <- if (oracle$contrast_fraction_near_zero[[i]] == 1) "state_all_zero" else
    if (isTRUE(oracle$contrast_changes_sign[[i]])) "state_sign_crossing" else "state_sign_stable"

  artifact <- file.path(output_root, "inference", "artifacts",
                        sprintf("direction-%06d.rds", d$direction_index))
  if (!identical(state, "completed") || !file.exists(artifact)) {
    failures[[length(failures) + 1L]] <- data.frame(
      source=d$source, target=d$target, direction_index=d$direction_index,
      execution_state=state, reason=states$failure_reason[match(d$direction_index, states$direction_index)])
    fit <- NULL
  } else fit <- readRDS(artifact)

  reportable <- !is.null(fit) && identical(fit$diagnostic_class, "converged")
  posterior_sign <- if (is.null(fit)) NA_integer_ else sign_tol(fit$a_mean)
  posterior_vs_oracle <- if (is.null(fit)) "not_reportable" else
    if (!reportable && identical(fit$diagnostic_class, "interaction_indeterminate")) "posterior_indeterminate" else
    if (!reportable) "not_reportable" else if (oracle_sign == 0L) "oracle_zero" else
    if (posterior_sign == oracle_sign) "agree" else "disagree"
  posterior_vs_glv <- if (is.null(fit)) "not_reportable" else
    if (!reportable && identical(fit$diagnostic_class, "interaction_indeterminate")) "posterior_indeterminate" else
    if (!reportable) "not_reportable" else if (A_zero) "absolute_zero" else
    if (posterior_sign == A_sign) "same_sign" else "opposite_sign"
  subject_elpd <- if (is.null(fit)) NULL else fit$elpd_subject_cross
  successful_obs <- if (is.null(subject_elpd) || !nrow(subject_elpd)) 0L else
    sum(subject_elpd$n_successful_test_observations, na.rm = TRUE)
  total_subject_elpd <- if (is.null(subject_elpd) || !nrow(subject_elpd)) NA_real_ else
    sum(subject_elpd$elpd[is.finite(subject_elpd$elpd)])

  rows[[i]] <- data.frame(
    source=d$source, target=d$target,
    glv_matrix_row=d$target, glv_matrix_column=d$source,
    glv_A_target_source=A_value, glv_A_structural_zero=A_zero,
    glv_A_sign=A_sign, corrected_oracle=oracle_value,
    corrected_oracle_sign=oracle_sign, glv_to_oracle_class=glv_class,
    state_contrast_min=oracle$contrast_min[[i]], state_contrast_max=oracle$contrast_max[[i]],
    state_contrast_mean=oracle$contrast_mean[[i]], state_contrast_median=oracle$contrast_median[[i]],
    fraction_state_positive=oracle$contrast_fraction_positive[[i]],
    fraction_state_negative=oracle$contrast_fraction_negative[[i]],
    fraction_state_zero=oracle$contrast_fraction_near_zero[[i]],
    fraction_state_opposite_to_glv=oracle$contrast_fraction_opposite_to_A[[i]],
    state_contrast_class=state_class, state_crosses_zero=oracle$contrast_changes_sign[[i]],
    execution_state=state,
    posterior_mean=if (is.null(fit)) NA_real_ else fit$a_mean,
    posterior_median=if (is.null(fit)) NA_real_ else fit$a_median,
    posterior_q2.5=if (is.null(fit)) NA_real_ else fit$a_q2.5,
    posterior_q97.5=if (is.null(fit)) NA_real_ else fit$a_q97.5,
    PSP=if (is.null(fit)) NA_real_ else max(fit$positive_sign_probability, fit$negative_sign_probability),
    LFSR=if (is.null(fit)) NA_real_ else fit$lfsr,
    diagnostic_class=if (is.null(fit)) NA_character_ else fit$diagnostic_class,
    predictive_eligible=reportable,
    elpd=if (is.null(fit)) NA_real_ else fit$kfold_mean,
    elpd_per_observation=if (successful_obs > 0L) total_subject_elpd / successful_obs else NA_real_,
    n_successful_test_observations=successful_obs,
    elpd_method=if (is.null(fit)) NA_character_ else fit$kfold_method,
    posterior_vs_oracle=posterior_vs_oracle, posterior_vs_glv=posterior_vs_glv,
    rhat=if (is.null(fit)) NA_real_ else fit$diag$worst_rhat,
    bulk_ESS=if (is.null(fit)) NA_real_ else fit$diag$min_ess_bulk,
    tail_ESS=if (is.null(fit)) NA_real_ else fit$diag$min_ess_tail,
    divergences=if (is.null(fit)) NA_integer_ else fit$diag$n_divergent,
    treedepth_hits=if (is.null(fit)) NA_integer_ else fit$diag$n_treedepth_hit,
    chain_EBFMI=if (is.null(fit)) NA_character_ else collapse_named(fit$diag$ebfmi_chain),
    chain_means=if (is.null(fit)) NA_character_ else collapse_named(fit$chain_aij_means),
    chain_probability_positive=if (is.null(fit)) NA_character_ else collapse_named(fit$chain_aij_positive_probabilities),
    chain_probability_negative=if (is.null(fit)) NA_character_ else collapse_named(fit$chain_aij_negative_probabilities),
    retry_count=if (is.null(fit)) NA_integer_ else fit$n_retries,
    stringsAsFactors=FALSE)
}
summary <- do.call(rbind, rows)
write_tsv(summary, file.path(output_root, "six_direction_hf2_summary.tsv"))
write_tsv(summary[c("source","target","glv_matrix_row","glv_matrix_column",
  "glv_A_target_source","glv_A_structural_zero","glv_A_sign","corrected_oracle",
  "corrected_oracle_sign","glv_to_oracle_class")],
  file.path(output_root, "six_direction_glv_projection_audit.tsv"))
write_tsv(summary[c("source","target","state_contrast_min","state_contrast_max",
  "state_contrast_mean","state_contrast_median","fraction_state_positive",
  "fraction_state_negative","fraction_state_zero","fraction_state_opposite_to_glv",
  "state_contrast_class","state_crosses_zero")],
  file.path(output_root, "six_direction_state_contrast.tsv"))
write_tsv(summary[c("source","target","posterior_vs_oracle","posterior_vs_glv",
  "diagnostic_class","predictive_eligible")],
  file.path(output_root, "six_direction_posterior_oracle_agreement.tsv"))

# The production split manifest stores held-out subjects. Materialize training
# subjects for the audit as the deterministic complement within each recorded
# eligible universe; do not infer or alter any split assignment.
split_path <- file.path(output_root, "six_direction_kfold_split_audit.tsv")
split_audit <- utils::read.delim(split_path, stringsAsFactors=FALSE, check.names=FALSE)
if (nrow(split_audit)) {
  split_audit$train_subjects <- vapply(seq_len(nrow(split_audit)), function(i) {
    universe <- strsplit(split_audit$eligible_subject_universe[[i]], ";", fixed=TRUE)[[1L]]
    test <- strsplit(split_audit$test_subjects[[i]], ";", fixed=TRUE)[[1L]]
    paste(setdiff(universe, test), collapse=";")
  }, character(1))
  write_tsv(split_audit, split_path)
}
empty_failures <- data.frame(source=character(),target=character(),direction_index=integer(),
                             execution_state=character(),reason=character())
write_tsv(if (length(failures)) do.call(rbind, failures) else empty_failures,
          file.path(output_root, "six_direction_failures.tsv"))

focus <- summary[summary$source == "species_7" & summary$target == "species_4", , drop=FALSE]
write_tsv(focus, file.path(output_root, "species_7_to_species_4_focused_audit.tsv"))
historical <- data.frame(
  claim=c("converged","transformed_oracle_agreement","absolute_A_agreement",
          "state_dependent_contrast_changes","absolute_zero_induced_transformed_effects"),
  historical=c("6/6","6/6","3/6","44/90","16/90"),
  corrected=c(sprintf("%d/6",sum(summary$diagnostic_class=="converged",na.rm=TRUE)),
    sprintf("%d/%d",sum(summary$posterior_vs_oracle=="agree"),sum(summary$posterior_vs_oracle %in% c("agree","disagree"))),
    sprintf("%d/%d",sum(summary$posterior_vs_glv=="same_sign"),sum(summary$posterior_vs_glv %in% c("same_sign","opposite_sign"))),
    sprintf("%d/90",sum(oracle_all$truth_only_contrast_changes_sign)),
    sprintf("%d/90",sum(oracle_all$A_zero & oracle_all$observed_sign != 0L))),
  stringsAsFactors=FALSE)
historical$classification <- ifelse(historical$historical==historical$corrected,"reproduced","changed")
write_tsv(historical, file.path(output_root,"historical_comparison.tsv"))

metadata <- data.frame(
  key=c("git_commit","package_version","R_version","cmdstanr_version","cmdstan_version",
        "manifest_sha256","zero_tolerance","source_target_orientation","public_seed",
        "shared_kfold_seed","K","R","predictive_scorer","executable_path",
        "resource_policy_schema","run_started_utc","finished_utc","working_tree_diff_stat"),
  value=c(system2("git",c("-C",shQuote(repo),"rev-parse","HEAD"),stdout=TRUE)[[1L]],
    as.character(utils::packageVersion("pclvbayes")),R.version.string,
    as.character(utils::packageVersion("cmdstanr")),as.character(cmdstanr::cmdstan_version()),
    actual_hash,format(zero_tolerance,scientific=TRUE),"source -> target = A[target, source]",
    "20260802","20260802","5","1","student-t-scale-mixture-kalman-ou-q16",
    readRDS(file.path(output_root,"inference","config.rds"))$executable,
    readRDS(file.path(output_root,"inference","config.rds"))$resource_policy$policy_schema,
    readRDS(file.path(output_root,"inference","config.rds"))$started,
    format(Sys.time(),tz="UTC",usetz=TRUE),
    paste(system2("git",c("-C",shQuote(repo),"diff","--stat"),stdout=TRUE),collapse=" | ")),
  stringsAsFactors=FALSE)
write_tsv(metadata,file.path(output_root,"run_metadata.tsv"))
cat("REPORT_ROOT=",output_root,"\n",sep="")
print(summary,row.names=FALSE)
print(historical,row.names=FALSE)
