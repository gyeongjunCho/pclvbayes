
#' @noRd
#' @keywords internal
.as_int <- function(x) suppressWarnings(as.integer(x))

#' @noRd
#' @keywords internal
.lfsr_safe <- function(p_two) {
  if (exists(".lfsr_from_two", mode = "function", inherits = TRUE)) {
    return(.lfsr_from_two(p_two))
  }
  pmin(pmax(p_two/2, 0), 0.5)
}

#' @noRd
#' @keywords internal
.softmax <- function(z) {
  z <- as.numeric(z)
  if (!length(z) || !any(is.finite(z))) return(rep(NA_real_, length(z)))
  z0 <- z - max(z[is.finite(z)], na.rm = TRUE)
  w <- exp(z0)
  s <- sum(w[is.finite(w)], na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(z)))
  w / s
}

#' @noRd
#' @keywords internal
.diag_ok_fun <- function(rhat, essb, esst, div, tdhit, ebfmi_min, thr) {
  (is.finite(rhat) & rhat < thr$rhat) &
    (is.finite(essb) & essb > thr$ess) &
    (is.finite(esst) & esst > thr$ess) &
    (is.na(div)   | .as_int(div)   <= thr$div) &
    (is.na(tdhit) | .as_int(tdhit) <= thr$tdhit) &
    (is.finite(ebfmi_min) & ebfmi_min >= thr$ebfmi)
}


#' Compute true stacking weights (by target) from pointwise ELPD
#'
#' @description
#' Uses \code{df$elpd_pointwise_cross} from \code{fit_glv_pairwise()} to compute
#' proper stacking weights (per \strong{to}) across incoming edges (\strong{from}→\strong{to}).
#' Each subject acts as a "data point". Only subjects with finite ELPD across
#' \emph{all} candidate edges for the same target are kept (intersection).
#'
#' @param df A list returned by \code{fit_glv_pairwise()} containing \code{$elpd_pointwise_cross}.
#' @param use_ppd Logical; if \code{TRUE}, uses per-point ELPD (\code{elpd_ppd})
#'   multiplied by \code{n_test} to align with aggregate ELPD scale.
#' @param min_models Minimum number of competing edges per target. Default \code{2}.
#' @param min_subjects Minimum number of common subjects. Default \code{1}.
#'
#' @return A tibble with columns \code{from}, \code{to}, \code{stacking}
#'   (weights sum to 1 within each \code{to}). If stacking is not feasible
#'   (no \code{elpd_pointwise_cross}, no \pkg{loo}, insufficient overlap),
#'   returns an empty tibble.
#' @noRd
#' @keywords internal
.compute_stacking_weights_cross <- function(df,
                                            use_ppd = FALSE,
                                            min_models = 2L,
                                            min_subjects = 1L) {
  if (!("elpd_pointwise_cross" %in% names(df)) ||
      is.null(df$elpd_pointwise_cross) ||
      !nrow(df$elpd_pointwise_cross)) {
    return(tibble::tibble(from = character(), to = character(), stacking = numeric()))
  }
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("Packages 'dplyr', 'tidyr', and 'tibble' are required.")
  }

  pw <- tibble::as_tibble(df$elpd_pointwise_cross)
  val_col <- if (isTRUE(use_ppd)) "elpd_ppd" else "elpd"

  out_list <- list()

  for (tg in sort(unique(pw$to))) {
    sub_pw <- pw[pw$to == tg, , drop = FALSE]
    if (!nrow(sub_pw)) next

    if (isTRUE(use_ppd) && all(c("elpd_ppd", "n_test") %in% names(sub_pw))) {
      sub_pw$val <- as.numeric(sub_pw$elpd_ppd) * as.numeric(sub_pw$n_test)
    } else {
      sub_pw$val <- as.numeric(sub_pw[[val_col]])
    }

    wide <- tidyr::pivot_wider(
      sub_pw[, c("subject","from","val")],
      names_from = "from", values_from = "val"
    )

    m <- as.matrix(wide[ , setdiff(names(wide), "subject"), drop = FALSE])
    row_ok <- apply(m, 1L, function(r) all(is.finite(r)))
    m <- m[row_ok, , drop = FALSE]

    K <- ncol(m); N <- nrow(m)
    if (is.null(K) || K < min_models || N < min_subjects) next

    if (!requireNamespace("loo", quietly = TRUE)) next

    # 지터: RNG 오염 없이 결정적으로 분리
    col_all_same <- apply(m, 2L, function(v) length(unique(v[is.finite(v)])) <= 1L)
    if (any(col_all_same)) {
      jitter_vals <- matrix(seq_len(N), nrow = N, ncol = sum(col_all_same))
      m[, col_all_same] <- m[, col_all_same, drop = FALSE] + jitter_vals * 1e-12
    }

    w <- try(loo::stacking_weights(m), silent = TRUE)
    if (inherits(w, "try-error") || is.null(w) || !length(w)) next

    nm <- colnames(m)
    out_list[[length(out_list) + 1L]] <- tibble::tibble(
      from = nm,
      to   = tg,
      stacking = as.numeric(w)
    )
  }

  if (!length(out_list)) {
    return(tibble::tibble(from = character(), to = character(), stacking = numeric()))
  }
  dplyr::bind_rows(out_list)
}

#' @noRd
#' @keywords internal
.compute_pseudobma_weights_cross <- function(df,
                                             edges = NULL,
                                             use_ppd = FALSE,
                                             plus = FALSE,
                                             BB_n = 1000L,
                                             min_models = 2L,
                                             min_subjects = 1L) {
  if (!("elpd_pointwise_cross" %in% names(df)) ||
      is.null(df$elpd_pointwise_cross) || !nrow(df$elpd_pointwise_cross)) {
    return(tibble::tibble(from = character(), to = character(), weight = numeric()))
  }
  if (!requireNamespace("tibble", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Packages 'tibble','tidyr','dplyr' are required.")
  }
  # plus=TRUE를 쓰려면 loo 필요
  if (isTRUE(plus) && !requireNamespace("loo", quietly = TRUE)) {
    return(tibble::tibble(from = character(), to = character(), weight = numeric()))
  }

  pw <- tibble::as_tibble(df$elpd_pointwise_cross)
  if (!is.null(edges) && nrow(edges)) {
    pw <- dplyr::semi_join(pw, edges, by = c("from","to"))
  }
  if (!nrow(pw)) {
    return(tibble::tibble(from = character(), to = character(), weight = numeric()))
  }

  out <- list()
  for (tg in sort(unique(pw$to))) {
    sub <- pw[pw$to == tg, , drop = FALSE]
    if (!nrow(sub)) next
    # 값 선택: elpd 또는 (elpd_ppd * n_test)
    if (isTRUE(use_ppd) && all(c("elpd_ppd","n_test") %in% names(sub))) {
      sub$val <- as.numeric(sub$elpd_ppd) * as.numeric(sub$n_test)
    } else {
      sub$val <- as.numeric(sub$elpd)
    }
    wide <- tidyr::pivot_wider(sub[, c("subject","from","val")],
                               names_from = "from", values_from = "val")
    X <- as.matrix(wide[ , setdiff(names(wide), "subject"), drop = FALSE])
    row_ok <- apply(X, 1L, function(r) all(is.finite(r)))
    X <- X[row_ok, , drop = FALSE]
    K <- ncol(X); N <- nrow(X)
    if (is.null(K) || K < min_models || N < min_subjects) next

    # 열 전부 동일하면 수치불안 방지용 미세 지터
    col_all_same <- apply(X, 2L, function(v) length(unique(v[is.finite(v)])) <= 1L)
    if (any(col_all_same)) {
      jitter_vals <- matrix(seq_len(N), nrow = N, ncol = sum(col_all_same))
      X[, col_all_same] <- X[, col_all_same, drop = FALSE] + jitter_vals * 1e-12
    }

    if (isTRUE(plus)) {
      w <- try(loo::pseudobma_weights(X, BB = TRUE, BB_n = BB_n), silent = TRUE)
    } else {
      # pseudo-BMA (BB=FALSE): exp(sum elpd) 정규화
      s <- colSums(X)
      s0 <- s - max(s[is.finite(s)], na.rm = TRUE)
      v <- exp(s0); v[!is.finite(v)] <- NA_real_
      if (all(!is.finite(v))) v[] <- NA_real_
      w <- v / sum(v, na.rm = TRUE)
    }
    if (inherits(w, "try-error") || !length(w)) next
    nm <- colnames(X)
    out[[length(out)+1L]] <- tibble::tibble(from = nm, to = tg, weight = as.numeric(w))
  }

  if (!length(out)) {
    return(tibble::tibble(from = character(), to = character(), weight = numeric()))
  }
  dplyr::bind_rows(out)
}
