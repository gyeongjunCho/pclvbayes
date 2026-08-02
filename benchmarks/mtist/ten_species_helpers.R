ten_species_direction_index <- function(taxa, seed) {
  tasks <- pclvbayes:::.make_pair_tasks(taxa, seed)
  rows <- lapply(tasks, function(task) {
    data.frame(
      task_id = rep(task$task_index, 2L),
      direction_index = c(2L * task$task_index - 1L, 2L * task$task_index),
      target = c(task$taxon_i, task$taxon_j),
      source = c(task$taxon_j, task$taxon_i),
      seed = unname(task$direction_seeds), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

mtist_truth_coefficient <- function(truth, target, source) {
  if (!all(c(target, source) %in% rownames(truth))) stop("Unknown truth taxon.")
  as.numeric(truth[target, source])
}

.bench_scalar <- function(z, name, default = NA) {
  if (!name %in% names(z)) return(default)
  value <- z[[name]][[1L]]
  if (is.null(value) || !length(value)) default else value
}

.bench_failure <- function(x, field) {
  if (is.null(x) || !length(x)) return(NA_character_)
  as.character(x[[field]] %||% NA_character_)
}

.bench_list_value <- function(z, name) {
  if (!name %in% names(z)) return(NULL)
  z[[name]][[1L]]
}

.residual_median <- function(summary, parameter) {
  if (!is.data.frame(summary) || !nrow(summary) ||
      !all(c("parameter", "median") %in% names(summary))) return(NA_real_)
  x <- summary$median[summary$parameter == parameter]
  if (!length(x) || !any(is.finite(x))) NA_real_ else stats::median(x[is.finite(x)])
}

expand_mtist_directions <- function(fit, truth, dataset_id, seed,
                                    alpha = 0.05) {
  raw <- fit$raw
  taxa <- rownames(truth)
  expected <- ten_species_direction_index(taxa, seed)
  cross <- summarize_bayes_pclv(fit, alpha = alpha, interaction = "cross")$cross
  lookup_cross <- function(target, source, field, default = NA) {
    hit <- cross[cross$to == target & cross$from == source, , drop = FALSE]
    if (nrow(hit) != 1L || !field %in% names(hit)) default else hit[[field]][[1L]]
  }
  rows <- vector("list", 2L * nrow(raw))
  k <- 0L
  for (r in seq_len(nrow(raw))) {
    z <- raw[r, , drop = FALSE]
    specs <- list(
      list(suffix = "ij", target = z$i[[1L]], source = z$j[[1L]],
           mean = "a_ij_mean", self = "a_ii_mean"),
      list(suffix = "ji", target = z$j[[1L]], source = z$i[[1L]],
           mean = "a_ji_mean", self = "a_jj_mean"))
    for (d in specs) {
      k <- k + 1L
      suffix <- d$suffix
      failure <- .bench_list_value(z, paste0("failure_", suffix))
      history <- .bench_list_value(z, paste0("retry_history_", suffix))
      residual <- .bench_list_value(z, paste0("chain_residual_summary_", suffix))
      cls <- as.character(.bench_scalar(z, paste0("diagnostic_class_", suffix),
                                        "sampler_diagnostics_failed"))
      coefficient <- as.numeric(.bench_scalar(z, d$mean, NA_real_))
      p_sign2 <- as.numeric(.bench_scalar(z, paste0("p_sign2_", suffix), NA_real_))
      lfsr <- if (is.finite(p_sign2)) p_sign2 / 2 else NA_real_
      folds_ok <- as.integer(.bench_scalar(z, paste0("kfold_folds_ok_", suffix), 0L))
      folds_fail <- as.integer(.bench_scalar(z, paste0("kfold_folds_fail_", suffix), 0L))
      if (!is.finite(folds_ok)) folds_ok <- 0L
      if (!is.finite(folds_fail)) folds_fail <- 0L
      elpd <- as.numeric(.bench_scalar(z, paste0("kfold_elpd_mean_", suffix), NA_real_))
      diag_ok <- isTRUE(lookup_cross(d$target, d$source, "diag_ok", FALSE))
      bayes_fdr <- as.numeric(lookup_cross(d$target, d$source, "bayes_FDR", NA_real_))
      stacking <- as.numeric(lookup_cross(d$target, d$source, "stacking", NA_real_))
      rows[[k]] <- data.frame(
        dataset_id = as.integer(dataset_id), task_id = r,
        direction_index = 2L * r - as.integer(suffix == "ij"),
        target = d$target, source = d$source,
        seed = expected$seed[expected$target == d$target & expected$source == d$source],
        main_fit_attempted = length(history) > 0L,
        main_fit_retained = is.finite(coefficient), retry_count = max(length(history) - 1L, 0L),
        diagnostic_class = cls,
        interaction_identifiable = isTRUE(.bench_scalar(z, paste0("interaction_identifiable_", suffix), FALSE)),
        residual_identifiable = isTRUE(.bench_scalar(z, paste0("residual_identifiable_", suffix), FALSE)),
        diagnostic_reasons = paste(.bench_list_value(z, paste0("indeterminate_reason_", suffix)) %||% character(), collapse = ";"),
        posterior_mean = coefficient, posterior_median = NA_real_,
        posterior_median_note = "not_exposed_in_canonical_raw_schema",
        posterior_sign_probability = if (is.finite(p_sign2)) 1 - p_sign2 / 2 else NA_real_,
        p_sign2 = p_sign2, lfsr = lfsr,
        chain_sign_agreement = as.logical(.bench_scalar(z, paste0("chain_sign_agreement_", suffix), NA)),
        a_self_mean = as.numeric(.bench_scalar(z, d$self, NA_real_)),
        sigma_median = .residual_median(residual, "sigma"),
        sd_ou_median = .residual_median(residual, "sd_ou"),
        phi_median = .residual_median(residual, "phi"),
        lambda_median = .residual_median(residual, "lambda"),
        nu_median = as.numeric(.bench_scalar(z, paste0("nu_median_", suffix), NA_real_)),
        rhat = as.numeric(.bench_scalar(z, paste0("rhat_", suffix), NA_real_)),
        ess_bulk = as.numeric(.bench_scalar(z, paste0("essb_", suffix), NA_real_)),
        ess_tail = as.numeric(.bench_scalar(z, paste0("esst_", suffix), NA_real_)),
        divergences = as.integer(.bench_scalar(z, paste0("div_", suffix), NA_integer_)),
        treedepth_hits = as.integer(.bench_scalar(z, paste0("tdhit_", suffix), NA_integer_)),
        ebfmi_min = as.numeric(.bench_scalar(z, paste0("ebfmi_min_", suffix), NA_real_)),
        bayesian_eligible = diag_ok && is.finite(bayes_fdr),
        kfold_attempted = folds_ok + folds_fail > 0L,
        kfold_completed = folds_ok > 0L,
        kfold_folds_ok = folds_ok, kfold_folds_failed = folds_fail,
        elpd_available = is.finite(elpd), aggregate_elpd = elpd,
        stacking_available = is.finite(stacking), stacking_weight = stacking,
        final_significant = diag_ok && is.finite(bayes_fdr) && bayes_fdr <= alpha,
        truth_coefficient = mtist_truth_coefficient(truth, d$target, d$source),
        truth_sign = sign(mtist_truth_coefficient(truth, d$target, d$source)),
        predicted_sign = if (diag_ok && is.finite(coefficient)) sign(coefficient) else NA_real_,
        elapsed_seconds = NA_real_, failure_stage = .bench_failure(failure, "stage"),
        failure_reason = .bench_failure(failure, "reason"), stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out <- out[match(paste(expected$target, expected$source),
                   paste(out$target, out$source)), , drop = FALSE]
  rownames(out) <- NULL
  if (nrow(out) != nrow(expected) || anyNA(out$task_id))
    stop("Directed benchmark record invariant violated.")
  out
}

build_mtist_benchmark_matrices <- function(directions, taxa, alpha = 0.05) {
  posterior <- conservative <- matrix(0, length(taxa), length(taxa),
                                       dimnames = list(taxa, taxa))
  mask <- directions[c("task_id", "direction_index", "target", "source")]
  mask$execution_failure <- !directions$main_fit_retained
  mask$diagnostic_failure <- directions$diagnostic_class == "sampler_diagnostics_failed"
  mask$interaction_indeterminate <- directions$diagnostic_class == "interaction_indeterminate"
  mask$residual_unstable <- directions$diagnostic_class == "interaction_stable_residual_unstable"
  mask$not_significant <- !directions$final_significant
  mask$elpd_unavailable <- !directions$elpd_available
  for (r in seq_len(nrow(directions))) {
    d <- directions[r, ]
    usable <- identical(d$diagnostic_class, "converged") && is.finite(d$posterior_mean)
    posterior[d$target, d$source] <- if (usable) d$posterior_mean else 0
    conservative[d$target, d$source] <- if (usable && d$final_significant) d$posterior_mean else 0
  }
  list(posterior = posterior, conservative = conservative, masks = mask)
}

aggregate_mtist_diagonal <- function(directions, taxa) {
  do.call(rbind, lapply(taxa, function(taxon) {
    values <- directions$a_self_mean[
      directions$target == taxon & directions$diagnostic_class == "converged" &
        is.finite(directions$a_self_mean)]
    data.frame(taxon = taxon, n_pair_specific_means = length(values),
      aggregated_value = if (length(values)) stats::median(values) else 0,
      diagonal_unavailable = !length(values),
      reason = if (length(values)) NA_character_ else "no_converged_pair_specific_self_effect_mean",
      aggregation = "median across pair-specific posterior self-effect means")
  }))
}

coverage_sign_metrics <- function(directions) {
  eligible <- directions$bayesian_eligible & is.finite(directions$predicted_sign)
  significant <- directions$final_significant & is.finite(directions$predicted_sign)
  truth_nz <- directions$truth_sign != 0
  metric <- function(name, numerator, denominator) data.frame(
    metric = name, value = if (denominator > 0) numerator / denominator else NA_real_,
    numerator = numerator, denominator = denominator,
    undefined_reason = if (denominator > 0) NA_character_ else "denominator_is_zero")
  rbind(
    metric("converged_coverage", sum(directions$diagnostic_class == "converged"), nrow(directions)),
    metric("eligible_sign_coverage", sum(eligible), nrow(directions)),
    metric("interaction_indeterminate_fraction", sum(directions$diagnostic_class == "interaction_indeterminate"), nrow(directions)),
    metric("diagnostic_failure_fraction", sum(directions$diagnostic_class == "sampler_diagnostics_failed"), nrow(directions)),
    metric("conditional_sign_accuracy", sum(directions$predicted_sign[eligible] == directions$truth_sign[eligible]), sum(eligible)),
    metric("absolute_A_cross_estimand_accuracy", sum(directions$predicted_sign[eligible] == directions$truth_sign[eligible]), sum(eligible)),
    metric("conditional_false_sign", sum(directions$predicted_sign[eligible] != directions$truth_sign[eligible]), sum(eligible)),
    metric("positive_truth_coverage", sum(eligible & directions$truth_sign > 0), sum(directions$truth_sign > 0)),
    metric("negative_truth_coverage", sum(eligible & directions$truth_sign < 0), sum(directions$truth_sign < 0)),
    metric("significant_count", sum(significant), 1L),
    metric("significant_sign_accuracy", sum(directions$predicted_sign[significant] == directions$truth_sign[significant]), sum(significant)),
    metric("truth_nonzero_omitted", sum(truth_nz & !significant), sum(truth_nz)),
    metric("truth_zero_false_positive", sum(!truth_nz & significant), sum(!truth_nz)),
    metric("absolute_A_zero_to_nonzero", sum(!truth_nz & significant), sum(!truth_nz)))
}

direction_asymmetry <- function(directions) {
  keys <- unique(vapply(seq_len(nrow(directions)), function(i)
    paste(sort(c(directions$target[[i]], directions$source[[i]])), collapse = "|"), character(1)))
  rows <- lapply(keys, function(key) {
    taxa <- strsplit(key, "|", fixed = TRUE)[[1L]]
    x <- directions[directions$target %in% taxa & directions$source %in% taxa, ]
    ok <- x$diagnostic_class == "converged"
    data.frame(pair = key, success_pattern = c("neither", "one", "both")[[sum(ok) + 1L]],
               different_diagnostic_class = length(unique(x$diagnostic_class)) > 1L)
  })
  do.call(rbind, rows)
}

stage_counts <- function(directions) data.frame(
  stage = c("directed_records", "main_fit_attempted", "main_fit_retained",
            "diagnostically_converged", "bayesian_eligible", "kfold_attempted",
            "kfold_completed", "elpd_available", "stacking_available", "final_significant"),
  count = c(nrow(directions), sum(directions$main_fit_attempted, na.rm = TRUE),
            sum(directions$main_fit_retained, na.rm = TRUE), sum(directions$diagnostic_class == "converged"),
            sum(directions$bayesian_eligible, na.rm = TRUE), sum(directions$kfold_attempted, na.rm = TRUE),
            sum(directions$kfold_completed, na.rm = TRUE), sum(directions$elpd_available, na.rm = TRUE),
            sum(directions$stacking_available, na.rm = TRUE), sum(directions$final_significant, na.rm = TRUE)))

select_confirmation_directions <- function(stage_b) {
  required <- c("task_id", "direction_index", "target", "source", "seed",
                "diagnostic_class")
  if (!all(required %in% names(stage_b))) stop("Stage B records lack selection fields.")
  selected <- stage_b[stage_b$diagnostic_class == "converged", , drop = FALSE]
  selected <- selected[order(selected$direction_index), , drop = FALSE]
  rownames(selected) <- NULL
  selected
}

confirmation_outcome <- function(stage_b_class, stage_b_sign,
                                 long_class, long_sign) {
  if (identical(long_class, "converged"))
    return(if (identical(as.numeric(stage_b_sign), as.numeric(long_sign)))
      "confirmed_converged_same_sign" else "confirmed_converged_sign_changed")
  switch(long_class,
    interaction_stable_residual_unstable = "downgraded_residual_unstable",
    interaction_indeterminate = "downgraded_indeterminate",
    sampler_diagnostics_failed = "long_run_sampler_failed",
    stop("Unknown long-run diagnostic class: ", long_class)
  )
}
