##### R/cor_helpers.R
#' @keywords internal
#' Effective sample size with autocorrelation correction.
#'
#' Methods:
#' - effn_method = "ar1": n * (1 - ρx(1)ρy(1)) / (1 + ρx(1)ρy(1))  (AR(1) 근사)
#' - effn_method = "bartlett": n / { 1 + 2 * sum_{k=1}^L w_k ρx(k)ρy(k) },
#'     L은 사용자가 지정(기본 sqrt(n)), w_k = 1 - k/(L+1) (Bartlett window)
#' - effn_method = "nw": Newey–West 자동 대역폭(bw = floor(4*(n/100)^(2/9)))을 사용,
#'     L = bw, w_k 동일(Bartlett 커널). (기본)
#'
#' Notes:
#' - Spearman 상관이면 상관에 사용한 것과 **동일 스케일**(랭크 또는 랭크-잔차)로 x,y를 넣으세요.
#' - time 정보와 집계 옵션을 같이 주면, neff 계산을 위한 ACF 산출 전에
#'   동일/근접 타임스탬프를 mean/median으로 집계하여 n 팽창을 방지합니다.
.eff_n <- function(x, y,
                   t = NULL,
                   aggregate_by_time = c("none", "mean", "median"),
                   min_n = 4,
                   clamp = TRUE,
                   time_tol = NULL,
                   effn_method = c("nw","bartlett","ar1"),
                   L = NULL,
                   nw_bw = NULL) {
  aggregate_by_time <- match.arg(aggregate_by_time)
  effn_method <- match.arg(effn_method)

  # --- 기본 전처리 ------------------------------------------------------------
  x <- as.numeric(x); y <- as.numeric(y)
  if (!is.null(t)) t <- suppressWarnings(as.numeric(t))

  ok <- is.finite(x) & is.finite(y)
  if (!is.null(t)) ok <- ok & is.finite(t)
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; y <- y[ok]; if (!is.null(t)) t <- t[ok]

  # --- 동일/근접 시간 스탬프 집계 (neff 계산용) ------------------------------
  if (!is.null(t) && aggregate_by_time != "none") {
    ord <- order(t); t <- t[ord]; x <- x[ord]; y <- y[ord]
    if (!is.null(time_tol) && is.finite(time_tol) && time_tol > 0) {
      grp <- c(1L, 1L + cumsum(diff(t) > time_tol))
      sp  <- split(seq_along(t), grp)
    } else {
      sp  <- split(seq_along(t), t) # 정확히 동일한 시간만 병합
    }
    agg_fun <- if (aggregate_by_time == "mean") base::mean else stats::median
    x <- as.numeric(vapply(sp, function(ix) agg_fun(x[ix], na.rm = TRUE), numeric(1)))
    y <- as.numeric(vapply(sp, function(ix) agg_fun(y[ix], na.rm = TRUE), numeric(1)))
  } else if (!is.null(t) && any(duplicated(t))) {
    warning("Duplicate time stamps detected but aggregate_by_time='none'. ",
            "Consider aggregating (mean/median) or set a small time_tol.",
            call. = FALSE)
  }

  n <- min(length(x), length(y))
  if (!is.finite(n) || n < 5) return(NA_real_)

  # --- 보조: lag-k ACF 벡터 ---------------------------------------------------
  acf_vec <- function(v, lag.max) {
    if (lag.max <= 0) return(numeric(0))
    out <- try(stats::acf(v, plot = FALSE, lag.max = lag.max,
                          na.action = stats::na.fail)$acf, silent = TRUE)
    if (inherits(out, "try-error") || !is.numeric(out)) return(rep(0, lag.max))
    a <- as.numeric(out[-1]) # drop lag 0
    a[!is.finite(a)] <- 0
    pmax(pmin(a, 0.99), -0.99)
  }

  if (effn_method == "ar1") {
    lag1_acf <- function(v) {
      out <- try(stats::acf(v, plot = FALSE, lag.max = 1,
                            na.action = stats::na.fail)$acf[2],
                 silent = TRUE)
      if (inherits(out, "try-error") || !is.finite(out)) 0 else as.numeric(out)
    }
    ax <- max(min(lag1_acf(x), 0.99), -0.99)
    ay <- max(min(lag1_acf(y), 0.99), -0.99)
    axy <- max(min(ax * ay, 0.95), -0.95)
    neff <- n * (1 - axy) / (1 + axy)
  } else {
    # Bartlett / Newey–West
    if (effn_method == "nw") {
      bw <- if (is.null(nw_bw) || !is.finite(nw_bw) || nw_bw <= 0)
        max(1L, floor(4 * (n / 100)^(2/9)))
      else floor(nw_bw)
      L_use <- bw
      w <- 1 - (seq_len(L_use)) / (L_use + 1)  # Bartlett kernel
    } else { # "bartlett"
      L_use <- if (is.null(L) || !is.finite(L) || L <= 0) max(1L, floor(sqrt(n))) else floor(L)
      w <- 1 - (seq_len(L_use)) / (L_use + 1)  # Bartlett window
    }
    ax <- acf_vec(x, L_use)
    ay <- acf_vec(y, L_use)
    gam <- sum(w * (ax * ay), na.rm = TRUE)
    denom <- 1 + 2 * gam
    if (!is.finite(denom) || denom <= 1e-6) denom <- 1e-6
    neff <- n / denom
  }

  if (!is.finite(neff)) neff <- n
  if (isTRUE(clamp)) neff <- max(min_n, min(n, neff))
  as.numeric(neff)
}


#' @keywords internal
.aggregate_by_time_series <- function(mat, tt, mode = c("mean","median"), tol = NULL) {
  mode <- match.arg(mode)
  if (is.null(tt) || length(tt) != nrow(mat)) return(list(mat = mat, time = tt))
  ord <- order(tt); tt2 <- tt[ord]; m2 <- as.matrix(mat)[ord, , drop = FALSE]
  if (!is.null(tol) && is.finite(tol) && tol > 0) {
    grp <- c(1L, 1L + cumsum(diff(tt2) > tol))
  } else {
    # 정확히 동일한 시간만 묶기
    grp <- as.integer(factor(tt2, levels = unique(tt2)))
  }
  idx_list <- split(seq_along(tt2), grp)
  agg_fun <- if (mode == "median") stats::median else base::mean
  mat_agg <- vapply(
    X = seq_len(ncol(m2)),
    FUN = function(j) vapply(idx_list, function(ix) agg_fun(m2[ix, j], na.rm = TRUE), numeric(1)),
    FUN.VALUE = numeric(length(idx_list))
  )
  if (!is.matrix(mat_agg)) mat_agg <- matrix(mat_agg, ncol = ncol(m2))
  # 대표 시간은 그룹 평균으로
  tt_agg <- vapply(idx_list, function(ix) mean(tt2[ix], na.rm = TRUE), numeric(1))
  list(mat = mat_agg, time = tt_agg)
}
