args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
repo <- normalizePath(file.path(script_dir, "../.."))
source(file.path(script_dir, "mtist_adapter.R"))
source(file.path(script_dir, "ten_species_helpers.R"))
source(file.path(script_dir, "oracle_estimand_audit.R"))
devtools::load_all(repo, quiet = TRUE)
root <- mtist_root()
result_dir <- file.path(script_dir, "results", "oracle_estimand_audit")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
noiseless_path <- file.path(result_dir, "official_noiseless_100.csv")
if (!file.exists(noiseless_path)) stop("Official noiseless regeneration is missing.")

study <- load_mtist_study(37L, root)
taxa <- study$taxa
A <- study$truth
growth <- scan(file.path(root, "mtist1.0", "ground_truths", "growth_rates",
                         "10_sp_gr_1.csv"), quiet = TRUE)
observed_meta <- study$sample_metadata
observed_meta$Sample <- rownames(observed_meta)
observed_absolute <- study$absolute
colnames(observed_absolute) <- taxa

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
initial_observed <- observed_absolute[ave(seq_len(nrow(observed_meta)), observed_meta$subject,
                                         FUN = seq_along) == 1L, , drop = FALSE]
initial_noiseless <- noiseless_absolute[ave(seq_len(nrow(noiseless_meta)), noiseless_meta$subject,
                                           FUN = seq_along) == 1L, , drop = FALSE]
if (!isTRUE(all.equal(unname(initial_observed), unname(initial_noiseless), tolerance = 0)))
  stop("Regenerated official initial states do not match dataset 37.")

relative_noiseless <- noiseless_absolute / rowSums(noiseless_absolute)
relative_observed <- observed_absolute / rowSums(observed_absolute)
mat_noiseless <- t(relative_noiseless); colnames(mat_noiseless) <- noiseless_meta$Sample
mat_observed <- t(relative_observed); colnames(mat_observed) <- observed_meta$Sample
smoothed_noiseless <- suppressWarnings(pclvbayes:::.precompute_spline_smoothed(
  mat_noiseless, noiseless_meta, taxa, 1e-6, 3L))
smoothed_observed <- suppressWarnings(pclvbayes:::.precompute_spline_smoothed(
  mat_observed, observed_meta, taxa, 1e-6, 3L))
if (inherits(smoothed_noiseless, "pclv_failure") || inherits(smoothed_observed, "pclv_failure"))
  stop("Canonical smoothing failed.")

directions <- ten_species_direction_index(taxa, 20260802L)
rows <- vector("list", nrow(directions))
contrast_lag <- vector("list", nrow(directions))
sensitivities <- vector("list", nrow(directions))
for (r in seq_len(nrow(directions))) {
  target <- directions$target[[r]]; source <- directions$source[[r]]
  contrast <- apply(noiseless_absolute, 1L, oracle_contrast,
                    target = target, source = source, A = A, taxa = taxa)
  derivative <- apply(noiseless_absolute, 1L, oracle_alr_derivative,
                      target = target, source = source, growth = growth,
                      A = A, taxa = taxa)
  instant_design <- oracle_raw_design(noiseless_absolute, noiseless_meta,
                                      target, source, taxa, derivative)
  finite_design <- oracle_raw_design(noiseless_absolute, noiseless_meta,
                                     target, source, taxa)
  instant_projection <- oracle_projection(instant_design$y, instant_design$xi,
                                          instant_design$xj, instant_design$subject)
  finite_projection <- oracle_projection(finite_design$y, finite_design$xi,
                                         finite_design$xj, finite_design$subject)
  noiseless_design <- oracle_design_from_smoothed(smoothed_noiseless, noiseless_meta,
                                                  target, source)
  observed_design <- oracle_design_from_smoothed(smoothed_observed, observed_meta,
                                                 target, source)
  noiseless_projection <- if (inherits(noiseless_design, "pclv_failure"))
    list(coefficient = NA_real_, sign = NA_integer_, rank = NA_integer_,
         columns = NA_integer_, condition = NA_real_, identifiable = FALSE,
         reason = noiseless_design$reason) else
    oracle_projection(noiseless_design$y, noiseless_design$xi,
                      noiseless_design$xj, noiseless_design$subject)
  observed_projection <- if (inherits(observed_design, "pclv_failure"))
    list(coefficient = NA_real_, sign = NA_integer_, rank = NA_integer_,
         columns = NA_integer_, condition = NA_real_, identifiable = FALSE,
         reason = observed_design$reason) else
    oracle_projection(observed_design$y, observed_design$xi,
                      observed_design$xj, observed_design$subject)
  cs <- oracle_contrast_summary(contrast)
  rows[[r]] <- cbind(directions[r, ], data.frame(
    A_value = A[target, source], A_sign = oracle_sign(A[target, source]),
    A_zero = A[target, source] == 0,
    cs,
    instantaneous_coefficient = instant_projection$coefficient,
    instantaneous_sign = instant_projection$sign,
    instantaneous_rank = instant_projection$rank,
    instantaneous_columns = instant_projection$columns,
    instantaneous_condition = instant_projection$condition,
    instantaneous_identifiable = instant_projection$identifiable,
    instantaneous_reason = instant_projection$reason,
    finite_coefficient = finite_projection$coefficient,
    finite_sign = finite_projection$sign,
    finite_rank = finite_projection$rank,
    finite_condition = finite_projection$condition,
    noiseless_coefficient = noiseless_projection$coefficient,
    noiseless_sign = noiseless_projection$sign,
    noiseless_rank = noiseless_projection$rank,
    noiseless_condition = noiseless_projection$condition,
    noiseless_reason = noiseless_projection$reason,
    observed_coefficient = observed_projection$coefficient,
    observed_sign = observed_projection$sign,
    observed_rank = observed_projection$rank,
    observed_condition = observed_projection$condition,
    observed_reason = observed_projection$reason,
    stringsAsFactors = FALSE))
  contrast_lag[[r]] <- data.frame(direction_index = directions$direction_index[[r]],
                                  subject = noiseless_meta$subject,
                                  time = noiseless_meta$time, contrast = contrast)
  if (directions$direction_index[[r]] %in% c(30L, 48L, 56L, 65L, 67L, 74L)) {
    subject_signs <- vapply(split(seq_len(nrow(observed_design)), observed_design$subject),
      function(ix) oracle_projection(observed_design$y[ix], observed_design$xi[ix],
                                     observed_design$xj[ix], observed_design$subject[ix])$sign,
      integer(1))
    subjects <- unique(observed_design$subject)
    loso_signs <- vapply(subjects, function(s)
      oracle_projection(observed_design$y[observed_design$subject != s],
                        observed_design$xi[observed_design$subject != s],
                        observed_design$xj[observed_design$subject != s],
                        observed_design$subject[observed_design$subject != s])$sign,
      integer(1))
    sensitivities[[r]] <- data.frame(
      direction_index = directions$direction_index[[r]],
      subject_signs = paste(names(subject_signs), subject_signs, sep = "=", collapse = ";"),
      loso_signs = paste(subjects, loso_signs, sep = "=", collapse = ";"),
      rest_weight_min = min(contrast), rest_weight_max = max(contrast),
      unsmoothed_sign = finite_projection$sign,
      smoothed_noiseless_sign = noiseless_projection$sign,
      stringsAsFactors = FALSE)
  }
}
audit <- do.call(rbind, rows)
contrasts <- do.call(rbind, contrast_lag)
sensitivity <- do.call(rbind, sensitivities[!vapply(sensitivities, is.null, logical(1))])
confirmation <- utils::read.delim(file.path(script_dir, "results",
  "ten_species_37_confirmation_4x2000_2000", "confirmation_comparison.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE)
focused <- merge(audit, confirmation, by = c("task_id", "direction_index", "target", "source", "seed"),
                 all = FALSE, sort = FALSE)
focused$primary_category <- mapply(oracle_primary_category,
  focused$A_value, focused$contrast_changes_sign, focused$instantaneous_sign,
  focused$finite_sign, focused$noiseless_sign, focused$observed_sign,
  focused$long_sign, USE.NAMES = FALSE)

sign_columns <- c("A_sign", "instantaneous_sign", "finite_sign",
                  "noiseless_sign", "observed_sign")
cross <- do.call(rbind, lapply(combn(sign_columns, 2L, simplify = FALSE), function(pair) {
  tab <- table(factor(audit[[pair[[1L]]]], levels = -1:1),
               factor(audit[[pair[[2L]]]], levels = -1:1))
  data.frame(from = pair[[1L]], to = pair[[2L]],
             from_sign = rep(-1:1, each = 3), to_sign = rep(-1:1, 3),
             count = as.integer(tab))
}))
agreement <- do.call(rbind, lapply(sign_columns[-1L], function(column) data.frame(
  comparison = paste("A_sign", column, sep = "_vs_"),
  agreement_n = sum(audit$A_sign == audit[[column]], na.rm = TRUE),
  denominator = sum(is.finite(audit$A_sign) & is.finite(audit[[column]])),
  truth_positive_agreement = sum(audit$A_sign == 1L & audit[[column]] == 1L),
  truth_positive_denominator = sum(audit$A_sign == 1L),
  truth_negative_agreement = sum(audit$A_sign == -1L & audit[[column]] == -1L),
  truth_negative_denominator = sum(audit$A_sign == -1L))))
summary <- data.frame(
  metric = c("directions", "contrast_sign_changing", "instantaneous_rank_deficient",
             "finite_rank_deficient", "noiseless_rank_deficient", "observed_rank_deficient",
             "A_zero_to_instantaneous_nonzero", "A_zero_to_finite_nonzero",
             "A_zero_to_noiseless_nonzero", "A_zero_to_observed_nonzero"),
  value = c(nrow(audit), sum(audit$contrast_changes_sign),
    sum(!audit$instantaneous_identifiable), sum(!is.finite(audit$finite_coefficient)),
    sum(!is.finite(audit$noiseless_coefficient)), sum(!is.finite(audit$observed_coefficient)),
    sum(audit$A_sign == 0 & audit$instantaneous_sign != 0),
    sum(audit$A_sign == 0 & audit$finite_sign != 0),
    sum(audit$A_sign == 0 & audit$noiseless_sign != 0),
    sum(audit$A_sign == 0 & audit$observed_sign != 0)),
  denominator = c(90L, 90L, rep(90L, 8L)))

write_tsv <- function(x, name) utils::write.table(x, file.path(result_dir, name),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write_tsv(audit, "all_direction_oracles.tsv")
write_tsv(focused, "focused_six.tsv")
write_tsv(contrasts, "instantaneous_contrasts.tsv")
write_tsv(sensitivity, "focused_sensitivity.tsv")
write_tsv(cross, "sign_cross_tabs.tsv")
write_tsv(agreement, "sign_agreement.tsv")
write_tsv(summary, "audit_summary.tsv")
cat("INITIAL_STATES_EXACT=TRUE\n")
print(focused[c("source", "target", "A_sign", "contrast_changes_sign",
  "instantaneous_sign", "finite_sign", "noiseless_sign", "observed_sign",
  "long_sign", "primary_category")], row.names = FALSE)
print(agreement, row.names = FALSE)
print(summary, row.names = FALSE)
