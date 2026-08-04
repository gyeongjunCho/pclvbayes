args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
repo <- normalizePath(file.path(script_dir, "../.."))
pkgload::load_all(repo, quiet = TRUE)
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "ten_species_helpers.R"))
source(file.path(script_dir, "oracle_estimand_audit.R"))

output_root <- file.path(
  script_dir, "results", "v021_hf2_six_direction_revalidation_v2"
)
manifest_path <- file.path(output_root, "direction_manifest.csv")
expected_hash <- "9e6907aa47953146f56e08b7e78f3c6f57619e8bc4ee2f31f0190f3d6364c005"
actual_hash <- unname(tools::md5sum(manifest_path))
sha256 <- system2("sha256sum", shQuote(manifest_path), stdout = TRUE)
if (!startsWith(sha256[[1L]], expected_hash)) stop("Frozen manifest hash changed.")
audit_root <- file.path(output_root, "deterministic_audit")
if (dir.exists(audit_root) && length(list.files(audit_root, all.files = TRUE, no.. = TRUE)))
  stop("Corrected deterministic audit output already exists.")
dir.create(audit_root, recursive = TRUE)

started <- Sys.time()
zero_tolerance <- 1e-10
root <- mtist_root(Sys.getenv("MTIST_ROOT", unset = "~/mtist"))
study <- load_mtist_study(37L, root)
taxa <- study$taxa
A <- study$truth
growth <- scan(file.path(root, "mtist1.0", "ground_truths", "growth_rates",
                         "10_sp_gr_1.csv"), quiet = TRUE)
observed_meta <- study$sample_metadata
observed_meta$Sample <- rownames(observed_meta)
observed_absolute <- study$absolute
colnames(observed_absolute) <- taxa

noiseless_root <- file.path(output_root, "deterministic_inputs")
noiseless_path <- file.path(noiseless_root, "official_noiseless_100.csv")
if (!file.exists(noiseless_path)) {
  status <- system2("python3", c(
    file.path(script_dir, "regenerate_noiseless_mtist.py"),
    "--mtist-root", root, "--output", noiseless_root))
  if (!identical(status, 0L) || !file.exists(noiseless_path))
    stop("Canonical noiseless trajectory regeneration failed.")
}
noiseless100 <- utils::read.csv(noiseless_path, check.names = FALSE)
noiseless_rows <- lapply(unique(observed_meta$subject), function(subject) {
  wanted <- observed_meta[observed_meta$subject == subject, , drop = FALSE]
  available <- noiseless100[as.character(noiseless100$subject) == subject, , drop = FALSE]
  nearest <- vapply(wanted$time, function(t) which.min(abs(available$time - t)), integer(1))
  if (any(abs(available$time[nearest] - wanted$time) > 1e-10))
    stop("Official noiseless times do not reproduce benchmark times.")
  out <- available[nearest, c("subject", "time", taxa)]
  out$Sample <- wanted$Sample
  out
})
noiseless <- do.call(rbind, noiseless_rows)
rownames(noiseless) <- NULL
noiseless_meta <- noiseless[c("Sample", "subject", "time")]
noiseless_absolute <- as.matrix(noiseless[, taxa])
relative_noiseless <- noiseless_absolute / rowSums(noiseless_absolute)
relative_observed <- observed_absolute / rowSums(observed_absolute)
mat_noiseless <- t(relative_noiseless); colnames(mat_noiseless) <- noiseless_meta$Sample
mat_observed <- t(relative_observed); colnames(mat_observed) <- observed_meta$Sample
smoothed_noiseless <- suppressWarnings(pclvbayes:::.precompute_spline_smoothed(
  mat_noiseless, noiseless_meta, taxa, 1e-6, 3L))
smoothed_observed <- suppressWarnings(pclvbayes:::.precompute_spline_smoothed(
  mat_observed, observed_meta, taxa, 1e-6, 3L))
if (inherits(smoothed_noiseless, "pclv_failure") ||
    inherits(smoothed_observed, "pclv_failure")) stop("Canonical smoothing failed.")
closure_noiseless <- max(abs(colSums(smoothed_noiseless) - 1))
closure_observed <- max(abs(colSums(smoothed_observed) - 1))

directions <- ten_species_direction_index(taxa, 20260802L)
rows <- vector("list", nrow(directions))
sensitivity_rows <- list()
for (r in seq_len(nrow(directions))) {
  target <- directions$target[[r]]
  source <- directions$source[[r]]
  truth_contrast <- apply(noiseless_absolute, 1L, oracle_contrast,
                          target = target, source = source, A = A, taxa = taxa)
  contrast <- apply(t(smoothed_observed), 1L, oracle_contrast,
                    target = target, source = source, A = A, taxa = taxa)
  derivative <- apply(noiseless_absolute, 1L, oracle_alr_derivative,
                      target = target, source = source, growth = growth,
                      A = A, taxa = taxa)
  instantaneous <- oracle_raw_design(
    noiseless_absolute, noiseless_meta, target, source, taxa, derivative)
  finite <- oracle_raw_design(
    noiseless_absolute, noiseless_meta, target, source, taxa)
  noiseless_design <- oracle_design_from_smoothed(
    smoothed_noiseless, noiseless_meta, target, source)
  observed_design <- oracle_design_from_smoothed(
    smoothed_observed, observed_meta, target, source)
  instant_projection <- oracle_projection(
    instantaneous$y, instantaneous$xi, instantaneous$xj, instantaneous$subject)
  finite_projection <- oracle_projection(finite$y, finite$xi, finite$xj, finite$subject)
  noiseless_projection <- oracle_projection(
    noiseless_design$y, noiseless_design$xi, noiseless_design$xj,
    noiseless_design$subject)
  observed_projection <- oracle_projection(
    observed_design$y, observed_design$xi, observed_design$xj,
    observed_design$subject)
  rest <- setdiff(taxa, c(target, source))
  complement_error <- max(abs(
    1 - smoothed_observed[target, ] - smoothed_observed[source, ] -
      colSums(smoothed_observed[rest, , drop = FALSE])
  ))
  cs <- oracle_contrast_summary(contrast, tolerance = zero_tolerance)
  truth_cs <- oracle_contrast_summary(truth_contrast, tolerance = zero_tolerance)
  A_value <- A[target, source]
  A_sign <- oracle_sign(A_value, zero_tolerance)
  rows[[r]] <- cbind(directions[r, ], data.frame(
    A_value = A_value, A_sign = A_sign,
    A_zero = A_value == 0,
    cs,
    truth_only_contrast_changes_sign = truth_cs$contrast_changes_sign,
    contrast_fraction_opposite_to_A = if (A_sign == 0L) NA_real_ else
      mean(vapply(contrast, oracle_sign, integer(1), tolerance = zero_tolerance) == -A_sign),
    raw_observations = nrow(observed_meta),
    instantaneous_transitions = nrow(instantaneous),
    finite_transitions = nrow(finite),
    observed_transitions = nrow(observed_design),
    complement_rest_max_error = complement_error,
    predictor_xi_mean = mean(observed_design$xi_unscaled),
    predictor_xi_scale = stats::sd(observed_design$xi_unscaled),
    predictor_xj_mean = mean(observed_design$xj_unscaled),
    predictor_xj_scale = stats::sd(observed_design$xj_unscaled),
    instantaneous_coefficient = instant_projection$coefficient,
    instantaneous_sign = instant_projection$sign,
    instantaneous_rank = instant_projection$rank,
    instantaneous_condition = instant_projection$condition,
    finite_coefficient = finite_projection$coefficient,
    finite_sign = finite_projection$sign,
    finite_rank = finite_projection$rank,
    finite_condition = finite_projection$condition,
    noiseless_coefficient = noiseless_projection$coefficient,
    noiseless_sign = noiseless_projection$sign,
    noiseless_rank = noiseless_projection$rank,
    noiseless_condition = noiseless_projection$condition,
    observed_coefficient = observed_projection$coefficient,
    observed_sign = observed_projection$sign,
    observed_rank = observed_projection$rank,
    observed_condition = observed_projection$condition,
    preprocessing_sign_changed = finite_projection$sign != observed_projection$sign,
    stringsAsFactors = FALSE))

  if (directions$direction_index[[r]] %in% c(30L, 48L, 56L, 65L, 67L, 74L)) {
    subjects <- unique(observed_design$subject)
    subject_signs <- vapply(split(seq_len(nrow(observed_design)), observed_design$subject),
      function(ix) oracle_projection(observed_design$y[ix], observed_design$xi[ix],
        observed_design$xj[ix], observed_design$subject[ix])$sign, integer(1))
    loso_signs <- vapply(subjects, function(sb) {
      keep <- observed_design$subject != sb
      oracle_projection(observed_design$y[keep], observed_design$xi[keep],
        observed_design$xj[keep], observed_design$subject[keep])$sign
    }, integer(1))
    sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.frame(
      direction_index = directions$direction_index[[r]],
      subject_signs = paste(names(subject_signs), subject_signs, sep = "=", collapse = ";"),
      loso_signs = paste(subjects, loso_signs, sep = "=", collapse = ";"),
      stringsAsFactors = FALSE)
  }
}
audit <- do.call(rbind, rows)
sensitivity <- do.call(rbind, sensitivity_rows)
summary <- data.frame(
  metric = c(
    "directions", "contrast_sign_changing",
    "A_zero_to_instantaneous_nonzero", "A_zero_to_finite_nonzero",
    "A_zero_to_noiseless_nonzero", "A_zero_to_observed_nonzero"
  ),
  value = c(
    nrow(audit), sum(audit$truth_only_contrast_changes_sign),
    sum(audit$A_zero & audit$instantaneous_sign != 0L),
    sum(audit$A_zero & audit$finite_sign != 0L),
    sum(audit$A_zero & audit$noiseless_sign != 0L),
    sum(audit$A_zero & audit$observed_sign != 0L)
  ),
  denominator = 90L
)
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
focused <- merge(manifest, audit, by = c("task_id", "direction_index", "target", "source", "seed"),
                 all.x = TRUE, sort = FALSE)
focused <- focused[match(manifest$direction_index, focused$direction_index), ]
focused <- merge(focused, sensitivity, by = "direction_index", all.x = TRUE, sort = FALSE)
focused <- focused[match(manifest$direction_index, focused$direction_index), ]

write_tsv <- function(x, path) utils::write.table(
  x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(audit, file.path(audit_root, "corrected_oracle_all_90.tsv"))
write_tsv(summary, file.path(audit_root, "corrected_oracle_summary.tsv"))
write_tsv(focused, file.path(audit_root, "corrected_six_direction_oracle.tsv"))
write_tsv(sensitivity, file.path(audit_root, "corrected_sensitivity.tsv"))
metadata <- list(
  schema = "v021_hf2_corrected_oracle_audit_v1",
  git_commit = system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
  dataset_id = 37L, direction_count = 90L,
  preprocessing_schema = "pclv_smoothed_full_composition_closure_v1",
  transition_schema = "subject_adjacent_predecessor_outcome_v1",
  projection_schema = "subject_adjusted_pair_to_rest_qr_v1",
  zero_tolerance = zero_tolerance,
  closure_error_noiseless = closure_noiseless,
  closure_error_observed = closure_observed,
  manifest_sha256 = expected_hash,
  started = format(started, tz = "UTC", usetz = TRUE),
  finished = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(metadata, file.path(audit_root, "audit_metadata.rds"))
cat("AUDIT_ROOT=", audit_root, "\n", sep = "")
print(summary, row.names = FALSE)
print(focused[, c("source", "target", "direction_index", "finite_coefficient",
                  "finite_sign", "observed_coefficient", "observed_sign",
                  "A_value", "A_sign")], row.names = FALSE)
