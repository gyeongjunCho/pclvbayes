official_mtist_es <- function(truth, prediction, mtist_root, exclude_diagonal = FALSE) {
  truth_file <- tempfile(fileext = ".csv")
  pred_file <- tempfile(fileext = ".csv")
  on.exit(unlink(c(truth_file, pred_file)), add = TRUE)
  utils::write.table(truth, truth_file, sep = ",", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)
  utils::write.table(prediction, pred_file, sep = ",", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)
  code <- paste(
    "import sys, numpy as np",
    "sys.path.insert(0, sys.argv[1])",
    "from mtist.mtist_utils import calculate_es_score",
    "t=np.loadtxt(sys.argv[2], delimiter=',')",
    "p=np.loadtxt(sys.argv[3], delimiter=',')",
    "print(repr(float(calculate_es_score(t,p,exclude_self_interaction=(sys.argv[4]=='1')))))",
    sep = ";"
  )
  output <- system2("python3", c("-c", shQuote(code), shQuote(mtist_root),
                                 shQuote(truth_file), shQuote(pred_file),
                                 if (exclude_diagonal) "1" else "0"),
                    stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  score <- suppressWarnings(as.numeric(tail(output, 1L)))
  if (status != 0L || length(score) != 1L || !is.finite(score))
    stop("Official MTIST ES scoring failed: ", paste(output, collapse = "\n"))
  score
}

.metric <- function(name, numerator, denominator, definition) {
  if (!is.finite(denominator) || denominator <= 0) {
    data.frame(metric = name, value = NA_real_, numerator = numerator,
               denominator = denominator, undefined_reason = "denominator_is_zero",
               definition = definition)
  } else {
    data.frame(metric = name, value = numerator / denominator, numerator = numerator,
               denominator = denominator, undefined_reason = NA_character_,
               definition = definition)
  }
}

sign_metrics <- function(truth, prediction, matrix_name, n_success, n_failed) {
  off <- row(truth) != col(truth)
  ts <- sign(truth[off]); ps <- sign(prediction[off])
  truth_nz <- ts != 0; truth_pos <- ts > 0; truth_neg <- ts < 0; truth_zero <- ts == 0
  pred_nz <- ps != 0
  rows <- list(
    .metric("sign_accuracy_truth_nonzero", sum(ps[truth_nz] == ts[truth_nz]), sum(truth_nz),
            "legacy alias: correct predicted sign / nonzero absolute-A truth edges"),
    .metric("absolute_A_sign_agreement", sum(ps[truth_nz] == ts[truth_nz]), sum(truth_nz),
            "correct predicted sign / nonzero absolute-A truth edges"),
    .metric("positive_recall", sum(ps[truth_pos] > 0), sum(truth_pos),
            "predicted positive / positive truth edges"),
    .metric("negative_recall", sum(ps[truth_neg] < 0), sum(truth_neg),
            "predicted negative / negative truth edges"),
    .metric("signed_precision", sum(pred_nz & ps == ts), sum(pred_nz),
            "correct signed predictions / predicted nonzero edges"),
    .metric("false_positive_rate_truth_zero", sum(pred_nz & truth_zero), sum(truth_zero),
            "legacy alias: predicted nonzero / zero absolute-A truth edges"),
    .metric("absolute_A_zero_to_nonzero", sum(pred_nz & truth_zero), sum(truth_zero),
            "predicted nonzero / zero absolute-A truth edges"),
    .metric("predicted_zero_fraction", sum(!pred_nz), length(ps),
            "predicted zero / all off-diagonal entries")
  )
  pos_tp <- sum(ps > 0 & ts > 0); pos_fp <- sum(ps > 0 & ts <= 0); pos_fn <- sum(ps <= 0 & ts > 0)
  neg_tp <- sum(ps < 0 & ts < 0); neg_fp <- sum(ps < 0 & ts >= 0); neg_fn <- sum(ps >= 0 & ts < 0)
  f1 <- function(tp, fp, fn) if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_
  recalls <- c(if (sum(truth_pos)) sum(ps[truth_pos] > 0) / sum(truth_pos) else NA_real_,
               if (sum(truth_neg)) sum(ps[truth_neg] < 0) / sum(truth_neg) else NA_real_)
  f1s <- c(f1(pos_tp, pos_fp, pos_fn), f1(neg_tp, neg_fp, neg_fn))
  macro_row <- function(name, values, definition) data.frame(
    metric = name, value = if (all(is.finite(values))) mean(values) else NA_real_,
    numerator = NA_real_, denominator = sum(is.finite(values)),
    undefined_reason = if (all(is.finite(values))) NA_character_ else "one_or_more_sign_classes_undefined",
    definition = definition
  )
  rows <- c(rows, list(
    macro_row("macro_recall_positive_negative", recalls, "mean of positive and negative recall"),
    macro_row("macro_f1_positive_negative", f1s, "mean one-vs-rest F1 for positive and negative signs"),
    data.frame(metric = c("successful_directions", "failed_directions", "failure_rate"),
               value = c(n_success, n_failed, n_failed / (n_success + n_failed)),
               numerator = c(n_success, n_failed, n_failed),
               denominator = c(1, 1, n_success + n_failed), undefined_reason = NA_character_,
               definition = c("successful off-diagonal directed fits", "failed off-diagonal directed fits",
                              "failed / all off-diagonal directed fits"))
  ))
  out <- do.call(rbind, rows); out$matrix <- matrix_name
  out[, c("matrix", "metric", "value", "numerator", "denominator", "undefined_reason", "definition")]
}
