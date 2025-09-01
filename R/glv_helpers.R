#' @keywords internal
get_abund_matrix_precomputed <- function(physeq) {
  mat <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) mat <- t(mat)
  mat
}

#' @keywords internal
get_sample_meta <- function(physeq, subject_col, time_col) {
  md <- data.frame(
    sample = sample_names(physeq),
    as(sample_data(physeq), "data.frame"),
    check.names = FALSE
  )
  md <- md |>
    dplyr::transmute(
      Sample  = sample,
      subject = .data[[subject_col]],
      time    = .data[[time_col]]
    )
  if (!is.numeric(md$time)) suppressWarnings(md$time <- as.numeric(md$time))
  md
}


#' @keywords internal
precompute_spline_smoothed <- function(mat_rel, meta_df,
                                       taxa_list = rownames(mat_rel),
                                       eps = 1e-6,
                                       spline_df = NULL,
                                       spline_spar = NULL,
                                       spline_cv = TRUE,
                                       min_unique_times = 3) {
  sm_mat <- matrix(NA_real_, nrow = nrow(mat_rel), ncol = ncol(mat_rel),
                   dimnames = dimnames(mat_rel))
  col_idx_all <- match(meta_df$Sample, colnames(mat_rel))

  taxa_use <- intersect(taxa_list, rownames(mat_rel))
  if (!length(taxa_use)) return(sm_mat)

  for (tx in taxa_use) {
    vec_pred <- rep(NA_real_, nrow(meta_df))

    for (sb in unique(meta_df$subject)) {
      idx   <- which(meta_df$subject == sb)
      times <- meta_df$time[idx]
      cols  <- col_idx_all[idx]
      vals  <- as.numeric(mat_rel[tx, cols])

      ok <- is.finite(times) & is.finite(vals)
      if (!any(ok)) next

      df   <- data.frame(time = times[ok], val = vals[ok])
      df2  <- stats::aggregate(val ~ time, df, mean)
      df2  <- df2[order(df2$time), , drop = FALSE]

      pred_log <- NULL
      if (nrow(df2) >= min_unique_times) {
        ylog <- log(pmax(df2$val, 0) + eps)
        rr <- .smooth_spline_robust(
          x = df2$time, y = ylog,
          spline_df = spline_df,
          spline_spar = spline_spar,
          use_cv = isTRUE(spline_cv),
          min_unique = min_unique_times,
          min_df = 3.0,
          max_df = NULL,
          default_df = 4.0
        )
        pred_log <- try(as.numeric(stats::approx(
          x = df2$time, y = rr$yhat, xout = times, rule = 2
        )$y), silent = TRUE)
        if (inherits(pred_log, "try-error") || !length(pred_log)) pred_log <- rr$yhat
      }

      if (is.null(pred_log)) {
        const_log <- mean(log(pmax(vals, 0) + eps))
        pred_log  <- rep(const_log, length(times))
      }

      vec_pred[idx] <- pmax(exp(pred_log) - eps, 0)
    }

    sm_mat[tx, col_idx_all] <- vec_pred
  }

  sm_mat
}


#' @keywords internal
build_pair_df_smoothed <- function(sm_mat, meta_df, j, i,
                                   eps = 1e-8,
                                   min_pairs = 4,
                                   min_dt = 1e-8,
                                   min_sd = 1e-12) {
  stopifnot(all(c("Sample","subject","time") %in% names(meta_df)))
  idx_i <- match(meta_df$Sample, colnames(sm_mat))
  xi_all <- as.numeric(sm_mat[i, idx_i])
  xj_all <- as.numeric(sm_mat[j, idx_i])

  df <- meta_df |>
    dplyr::mutate(
      xi_raw = pmax(xi_all, 0),
      xj_raw = pmax(xj_all, 0)
    ) |>
    dplyr::arrange(subject, time)

  pair_df <- df |>
    dplyr::group_by(subject) |>
    dplyr::arrange(time, .by_group = TRUE) |>
    dplyr::mutate(
      time_next = dplyr::lead(time),
      xi_next   = dplyr::lead(xi_raw),
      xj_next   = dplyr::lead(xj_raw),
      dt        = time_next - time
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      is.finite(time), is.finite(time_next),
      is.finite(xi_raw), is.finite(xj_raw),
      is.finite(xi_next), is.finite(xj_next),
      is.finite(dt), dt > min_dt
    ) |>
    dplyr::mutate(
      rest_now  = pmax(0, 1 - xi_raw - xj_raw),
      rest_next = pmax(0, 1 - xi_next - xj_next),
      alr_i_now  = log(pmax(xi_raw,  0) + eps) - log(pmax(rest_now,  0) + eps),
      alr_i_next = log(pmax(xi_next, 0) + eps) - log(pmax(rest_next, 0) + eps),
      y = (alr_i_next - alr_i_now) / dt
    ) |>
    dplyr::transmute(subject, time, y, xi_raw, xj_raw) |>
    dplyr::filter(is.finite(y), is.finite(xi_raw), is.finite(xj_raw))

  if (nrow(pair_df) < min_pairs) return(NULL)
  if (stats::sd(pair_df$xi_raw) < min_sd || stats::sd(pair_df$xj_raw) < min_sd) return(NULL)
  pair_df
}

#' @keywords internal
.z <- function(v, cap = 5){
  v_fin <- v[is.finite(v)]
  m <- if (length(v_fin)) mean(v_fin) else 0
  s <- if (length(v_fin) > 1) sd(v_fin) else 1
  if (!is.finite(s) || s == 0) s <- 1
  z <- (v - m) / s
  z[!is.finite(z)] <- 0
  z[z >  cap] <-  cap
  z[z < -cap] <- -cap
  z
}

#' @keywords internal
safe_draws_df <- function(fit){
  dd <- posterior::as_draws_df(fit$draws())
  keep <- intersect(
    c("a_ij","a_ii","r0",
      "sigma","sd_ou","phi",
      "sigma_ou","lambda","sigma_pred","tau_r","nu"),
    names(dd)
  )
  # keep가 비어도 0-열 data.frame을 돌려 안정성 확보
  dd[, keep, drop = FALSE]
}

#' @keywords internal
.ebfmi_warn_from_fit <- function(fit, thr = 0.3) {
  txt <- try(capture.output(fit$cmdstan_diagnose()), silent = TRUE)
  if (inherits(txt, "try-error") || is.null(txt)) return(FALSE)
  any(grepl(sprintf("E-BFMI .* less than %.1f", thr), txt))
}

#' @keywords internal
.ebfmi_chainwise_from_energy <- function(fit) {
  sdiag <- try(fit$sampler_diagnostics(), silent = TRUE)
  if (inherits(sdiag, "try-error") || is.null(sdiag)) return(rep(NA_real_, 0))
  vars <- dimnames(sdiag)[[3]]
  pos  <- match("energy__", vars, nomatch = 0L)
  if (pos == 0L) return(rep(NA_real_, 0))
  E <- sdiag[, , pos, drop = FALSE]  # draws x chains x 1
  C <- dim(E)[2]
  eb <- rep(NA_real_, C)
  for (c in seq_len(C)) {
    ec <- as.numeric(E[, c, 1]); ec <- ec[is.finite(ec)]
    if (length(ec) < 3) next
    v <- stats::var(ec); if (!is.finite(v) || v <= 0) next
    d <- diff(ec)
    eb[c] <- mean(d * d) / v
  }
  eb
}

#' @keywords internal
summarise_diag <- function(fit,
                           pars = c("a_ij","a_ii","r0","sigma","sd_ou","phi"),
                           max_treedepth = 12) {
  # (1) 파라미터 요약치: R-hat / ESS (가용한 것만 선택)
  all_dd <- fit$draws()
  avail  <- dimnames(all_dd)$variables
  use_pars <- intersect(pars, avail)
  if (!length(use_pars)) use_pars <- avail
  sdtab <- posterior::summarise_draws(fit$draws(use_pars))
  worst_rhat   <- suppressWarnings(max(sdtab$rhat,      na.rm = TRUE))
  min_ess_bulk <- suppressWarnings(min(sdtab$ess_bulk,  na.rm = TRUE))
  min_ess_tail <- suppressWarnings(min(sdtab$ess_tail,  na.rm = TRUE))
  if (!is.finite(worst_rhat))   worst_rhat   <- NA_real_
  if (!is.finite(min_ess_bulk)) min_ess_bulk <- NA_real_
  if (!is.finite(min_ess_tail)) min_ess_tail <- NA_real_

  # (2) 샘플러 진단: divergent__, treedepth__, energy__
  sdiag_df <- posterior::as_draws_df(fit$sampler_diagnostics())

  n_div <- if ("divergent__" %in% names(sdiag_df)) {
    sum(as.integer(sdiag_df[["divergent__"]]), na.rm = TRUE)
  } else NA_integer_

  n_treedepth_hit <- if ("treedepth__" %in% names(sdiag_df)) {
    sum(as.integer(sdiag_df[["treedepth__"]] >= max_treedepth), na.rm = TRUE)
  } else NA_integer_

  # (3) E-BFMI: mean(diff(E)^2) / var(E) (체인별 계산 후 요약)
  ebfmi_min <- ebfmi_med <- NA_real_
  if ("energy__" %in% names(sdiag_df)) {
    if (".chain" %in% names(sdiag_df)) {
      eb <- vapply(split(sdiag_df[["energy__"]], sdiag_df[[".chain"]]), function(ev) {
        e <- as.numeric(ev); e <- e[is.finite(e)]
        if (length(e) < 3) return(NA_real_)
        v <- stats::var(e); if (!is.finite(v) || v <= 0) return(NA_real_)
        mean(diff(e)^2) / v
      }, numeric(1))
      ebfmi_min <- suppressWarnings(min(eb, na.rm = TRUE)); if (!is.finite(ebfmi_min)) ebfmi_min <- NA_real_
      ebfmi_med <- suppressWarnings(stats::median(eb, na.rm = TRUE)); if (!is.finite(ebfmi_med)) ebfmi_med <- NA_real_
    } else {
      e <- as.numeric(sdiag_df[["energy__"]]); e <- e[is.finite(e)]
      if (length(e) >= 3) {
        v <- stats::var(e)
        if (is.finite(v) && v > 0) {
          val <- mean(diff(e)^2) / v
          ebfmi_min <- val; ebfmi_med <- val
        }
      }
    }
  }

  # (4) 총 저장 draw 수(스칼라)
  n_draws <- posterior::ndraws(fit$draws())

  list(
    worst_rhat      = worst_rhat,
    min_ess_bulk    = min_ess_bulk,
    min_ess_tail    = min_ess_tail,
    n_divergent     = n_div,
    n_treedepth_hit = n_treedepth_hit,
    ebfmi_min       = ebfmi_min,
    ebfmi_med       = ebfmi_med,
    n_draws         = as.integer(n_draws)
  )
}

#' @keywords internal
.lfsr_from_two_sided <- function(p_two) pmax(pmin(p_two / 2, 0.5), 0)

#' @keywords internal
add_bayes_fdr <- function(df, p_cols, alpha = 0.10) {
  stopifnot(all(p_cols %in% names(df)))
  for (nm in p_cols) {
    lfsr <- .lfsr_from_two_sided(df[[nm]])
    q <- .q_from_lfsr(lfsr)
    base <- sub("^p_", "", nm)
    df[[paste0("lfsr_", base)]]    <- lfsr
    df[[paste0("q_bayes_", base)]] <- q
    df[[paste0("keep_", base)]]    <- !is.na(q) & q <= alpha
    df[[paste0("psp_", base)]]     <- 1 - lfsr
  }
  df
}

#' @keywords internal
.clip01 <- function(x) pmin(pmax(x, 0), 1)

#' @keywords internal
.lfsr_from_two <- function(p_two) pmin(pmax(p_two/2, 0), 0.5)

#' @keywords internal
.q_from_lfsr <- function(v) {
  q <- rep(NA_real_, length(v))
  nn <- which(!is.na(v))
  if (length(nn)) {
    l  <- pmin(pmax(v[nn], 0), 0.5)
    oo <- nn[order(l)]
    q[oo] <- cumsum(l[order(l)]) / seq_along(oo)
  }
  q
}

#' @keywords internal
.ok_diag <- function(rhat, essb, esst, div, tdhit, n_draws = NA_real_) {
  # 길이 통일
  L <- max(length(rhat), length(essb), length(esst), length(div), length(tdhit), length(n_draws))
  rhat    <- rep_len(rhat,    L)
  essb    <- rep_len(essb,    L)
  esst    <- rep_len(esst,    L)
  div     <- rep_len(div,     L)
  tdhit   <- rep_len(tdhit,   L)
  n_draws <- rep_len(n_draws, L)

  # 임계값
  rhat_thr <- 1.05
  ess_thr  <- 400

  # n_draws가 스칼라가 아니어도 element-wise로 처리
  nd_ok      <- is.finite(n_draws) & n_draws > 0
  div_max    <- ifelse(nd_ok, ceiling(0.001 * n_draws), 8L)    # ~0.1%
  tdepth_max <- ifelse(nd_ok, ceiling(0.010 * n_draws), 80L)   # ~1%

  (is.finite(rhat)  & rhat < rhat_thr) &
    (is.finite(essb)  & essb > ess_thr)  &
    (is.finite(esst)  & esst > ess_thr)  &
    (is.finite(div)   & div  <= div_max) &
    (is.finite(tdhit) & tdhit<= tdepth_max)
}

# ----- Zero-aware ALR building & transforms -----

#' @keywords internal
.make_pair_inputs_glv <- function(pair_df,
                                  transform = c("alr","raw"),
                                  lag = 1,
                                  zero_mode_alr = c("minpos_time","minpos_subject","lib","fixed"),
                                  minpos_alpha = 0.5,
                                  minpos_base = c("ij","triplet"),
                                  eps_fixed = 1e-8,
                                  lib_eps_c = 0.65,
                                  rest_floor_frac = 1.0,
                                  alr_cap = Inf,
                                  smooth_scale = c("logra","alr"),
                                  alr_spline_df = NULL,
                                  alr_spline_spar = NULL,
                                  alr_spline_cv = TRUE,
                                  nz_partner_min_frac = 0.15,
                                  standardize_by_subject = TRUE,
                                  z_mode = c("subject","global","none"),
                                  z_external = NULL) {
  transform     <- match.arg(transform)
  zero_mode_alr <- match.arg(zero_mode_alr)
  minpos_base   <- match.arg(minpos_base)
  smooth_scale  <- match.arg(smooth_scale)
  z_mode        <- match.arg(z_mode)

  stopifnot(all(c("subject","time","xi_raw","xj_raw") %in% names(pair_df)))
  df <- pair_df[order(pair_df$subject, pair_df$time), , drop = FALSE]

  df$xi_raw   <- pmax(as.numeric(df$xi_raw), 0)
  df$xj_raw   <- pmax(as.numeric(df$xj_raw), 0)
  df$rest_raw <- pmax(0, 1 - df$xi_raw - df$xj_raw)

  by_s <- split(seq_len(nrow(df)), df$subject)
  lagv <- function(v) unsplit(lapply(by_s, function(ix) dplyr::lag(v[ix], lag)), df$subject)
  diff_over_dt <- function(v, t) {
    unsplit(lapply(by_s, function(ix) {
      vi <- v[ix]; ti <- t[ix]
      vi_lag <- dplyr::lag(vi, lag)
      dt     <- as.numeric(ti - dplyr::lag(ti, lag))
      # 하한: subject 내 양의 dt들의 중앙값의 25% (예)
      dt_pos <- dt[is.finite(dt) & dt > 0]
      dt_min <- if (length(dt_pos)) 0.25 * stats::median(dt_pos) else 1
      dt_adj <- pmax(dt, dt_min)
      out <- (vi - vi_lag) / dt_adj
      out[!is.finite(dt) | dt <= 0] <- NA_real_
      out
    }), df$subject)
  }

  lib <- if ("libsize" %in% names(df)) df$libsize else NA_real_

  subj_minpos <- unsplit(lapply(by_s, function(ix) {
    base <- if (minpos_base == "ij") {
      c(df$xi_raw[ix][df$xi_raw[ix] > 0], df$xj_raw[ix][df$xj_raw[ix] > 0])
    } else {
      c(df$xi_raw[ix][df$xi_raw[ix] > 0],
        df$xj_raw[ix][df$xj_raw[ix] > 0],
        df$rest_raw[ix][df$rest_raw[ix] > 0])
    }
    if (length(base)) min(base) else NA_real_
  }), df$subject)

  make_triplet_row <- function(i) {
    xi <- df$xi_raw[i]; xj <- df$xj_raw[i]; xr <- df$rest_raw[i]
    eps_t <- switch(zero_mode_alr,
                    "minpos_time" = {
                      base_pos <- c(if (xi > 0) xi, if (xj > 0) xj, if (xr > 0) xr)
                      if (length(base_pos)) minpos_alpha * min(base_pos) else eps_fixed
                    },
                    "minpos_subject" = {
                      mp <- subj_minpos[i]
                      if (is.finite(mp) && mp > 0) minpos_alpha * mp else eps_fixed
                    },
                    "lib" = {
                      L <- lib[i]; if (!is.finite(L) || L <= 0) L <- 1
                      max(eps_fixed, lib_eps_c / L)
                    },
                    "fixed" = eps_fixed)
    xi <- if (xi > 0) xi else eps_t
    xj <- if (xj > 0) xj else eps_t
    xr <- if (xr > 0) xr else eps_t
    xr <- max(xr, rest_floor_frac * eps_t)

    s <- xi + xj + xr
    if (!is.finite(s) || s <= 0) s <- 1
    xi_ <- xi / s; xj_ <- xj / s; xr_ <- xr / s

    eps_star <- eps_t / s
    c(xi_, xj_, xr_, eps_star)
  }

  trip <- t(vapply(seq_len(nrow(df)), make_triplet_row, numeric(4)))
  colnames(trip) <- c("xi","xj","xr","eps_star")
  alr_i <- log(trip[, "xi"]) - log(trip[, "xr"])
  alr_j <- log(trip[, "xj"]) - log(trip[, "xr"])

  if (is.finite(alr_cap)) {
    alr_i <- pmax(pmin(alr_i,  alr_cap), -alr_cap)
    alr_j <- pmax(pmin(alr_j,  alr_cap), -alr_cap)
  } else {
    es <- pmax(trip[, "eps_star"], .Machine$double.eps)
    num <- pmax(1 - 2 * es, .Machine$double.eps)
    B_theory <- log(num / es)
    cap_vec  <- pmin(8, pmax(0, B_theory))
    alr_i <- pmax(pmin(alr_i,  cap_vec), -cap_vec)
    alr_j <- pmax(pmin(alr_j,  cap_vec), -cap_vec)
  }

  smooth_edf_mean <- NA_real_
  if (smooth_scale == "alr") {
    smooth_one <- function(v, t) {
      ok <- is.finite(v) & is.finite(t)
      if (sum(ok) < 3L || length(unique(t[ok])) < 3L) {
        return(list(y = v, df = NA_real_))
      }
      rr <- .smooth_spline_robust(
        x = t[ok], y = v[ok],
        spline_df = alr_spline_df,
        spline_spar = alr_spline_spar,
        use_cv = isTRUE(alr_spline_cv),
        min_unique = 3,
        min_df = 3.0,
        max_df = NULL,
        default_df = 4.0
      )
      yhat <- rep(NA_real_, length(t))
      yhat[ok] <- as.numeric(stats::approx(
        x = t[ok], y = rr$yhat, xout = t[ok], rule = 2
      )$y)
      list(y = yhat, df = rr$df)
    }
    edf_i <- edf_j <- rep(NA_real_, length(by_s))
    names(edf_i) <- names(edf_j) <- names(by_s)
    for (sname in names(by_s)) {
      ix <- by_s[[sname]]
      res_i <- smooth_one(alr_i[ix], df$time[ix])
      res_j <- smooth_one(alr_j[ix], df$time[ix])
      alr_i[ix] <- res_i$y
      alr_j[ix] <- res_j$y
      edf_i[sname] <- res_i$df
      edf_j[sname] <- res_j$df
    }
    edf_subj <- rowMeans(cbind(edf_i, edf_j), na.rm = TRUE)
    smooth_edf_mean <- mean(edf_subj, na.rm = TRUE)
    if (!is.finite(smooth_edf_mean)) smooth_edf_mean <- NA_real_
  }

  if (transform == "alr") {
    xi <- lagv(alr_i)
    xj <- lagv(alr_j)
    y  <- diff_over_dt(alr_i, df$time)
  } else {
    xi <- lagv(df$xi_raw)
    xj <- lagv(df$xj_raw)
    y  <- diff_over_dt(alr_i, df$time)  # response stays ΔALR_i/Δt
  }

  # 파트너 희소성 필터(옵션)
  if (is.finite(nz_partner_min_frac) && nz_partner_min_frac > 0) {
    by_s2 <- split(seq_len(nrow(df)), df$subject)
    keep_mask <- rep(TRUE, nrow(df))
    for (sname in names(by_s2)) {
      ix <- by_s2[[sname]]
      frac_nz <- mean(is.finite(xj[ix]) & xj[ix] != 0, na.rm = TRUE)
      if (is.finite(frac_nz) && frac_nz < nz_partner_min_frac) {
        keep_mask[ix] <- FALSE
      }
    }
  } else {
    keep_mask <- rep(TRUE, nrow(df))
  }

  dat <- data.frame(
    subject = df$subject,
    time    = df$time,
    y = as.numeric(y),
    xi = as.numeric(xi),
    xj = as.numeric(xj)
  )

  z_stats <- NULL
  if (z_mode == "none") {
    # no standardization
  } else if (!is.null(z_external) && z_mode == "global") {
    mu_y  <- z_external$mu_y;  sd_y  <- z_external$sd_y
    mu_xi <- z_external$mu_xi; sd_xi <- z_external$sd_xi
    mu_xj <- z_external$mu_xj; sd_xj <- z_external$sd_xj
    if (!is.finite(sd_y)  || sd_y  <= 0) sd_y  <- 1
    if (!is.finite(sd_xi) || sd_xi <= 0) sd_xi <- 1
    if (!is.finite(sd_xj) || sd_xj <= 0) sd_xj <- 1
    dat$y  <- (dat$y  - mu_y)  / sd_y
    dat$xi <- (dat$xi - mu_xi) / sd_xi
    dat$xj <- (dat$xj - mu_xj) / sd_xj
  } else {
    if (isTRUE(standardize_by_subject) && z_mode == "subject") {
      by_s2 <- split(seq_len(nrow(dat)), dat$subject)
      dat$y  <- unsplit(lapply(by_s2, function(ix) .z(dat$y[ix])),  dat$subject)
      dat$xi <- unsplit(lapply(by_s2, function(ix) .z(dat$xi[ix])), dat$subject)
      dat$xj <- unsplit(lapply(by_s2, function(ix) .z(dat$xj[ix])), dat$subject)
    } else {
      mu_y  <- mean(dat$y[is.finite(dat$y)],  na.rm = TRUE); sd_y  <- stats::sd(dat$y,  na.rm = TRUE); if (!is.finite(sd_y)  || sd_y  == 0) sd_y  <- 1
      mu_xi <- mean(dat$xi[is.finite(dat$xi)], na.rm = TRUE); sd_xi <- stats::sd(dat$xi, na.rm = TRUE); if (!is.finite(sd_xi) || sd_xi == 0) sd_xi <- 1
      mu_xj <- mean(dat$xj[is.finite(dat$xj)], na.rm = TRUE); sd_xj <- stats::sd(dat$xj, na.rm = TRUE); if (!is.finite(sd_xj) || sd_xj == 0) sd_xj <- 1
      dat$y  <- (dat$y  - mu_y)  / sd_y
      dat$xi <- (dat$xi - mu_xi) / sd_xi
      dat$xj <- (dat$xj - mu_xj) / sd_xj
      z_stats <- list(mu_y = mu_y, sd_y = sd_y,
                      mu_xi = mu_xi, sd_xi = sd_xi,
                      mu_xj = mu_xj, sd_xj = sd_xj)
    }
  }

  ok <- is.finite(dat$y) & is.finite(dat$xi) & is.finite(dat$xj) & is.finite(dat$time) & keep_mask
  dat <- dat[ok, , drop = FALSE]
  if (!nrow(dat)) return(NULL)

  attr(dat, "smooth_edf_mean") <- smooth_edf_mean
  attr(dat, "smooth_scale")    <- smooth_scale
  attr(dat, "smoothed")        <- smooth_scale %in% c("alr","logra")
  if (!is.null(z_stats)) attr(dat, "z_stats") <- z_stats

  dat
}

# -------- cmdstanr call silencer & loglik extractor --------

#' @keywords internal
.call_sample_silently <- function(mod, args, silent = TRUE) {
  if (!isTRUE(silent)) return(do.call(mod$sample, args))

  # 강제: Stan 쪽 진행표시/메시지 끄기
  args$refresh <- 0L
  args$show_messages <- FALSE

  # stdout / message를 null로 잠시 리다이렉트
  out_con <- file(nullfile(), open = "wt")
  msg_con <- file(nullfile(), open = "wt")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(sink(), silent = TRUE)
    try(close(msg_con), silent = TRUE)
    try(close(out_con), silent = TRUE)
  }, add = TRUE)

  sink(out_con)                  # stdout
  sink(msg_con, type = "message")# messages

  # 남은 warning/message까지 무음
  fit <- suppressWarnings(suppressMessages(
    do.call(mod$sample, args)
  ))
  fit
}

#' @keywords internal
.smooth_spline_robust <- function(x, y,
                                  spline_df = NULL,
                                  spline_spar = NULL,
                                  use_cv = TRUE,
                                  min_unique = 3,
                                  min_df = 3.0,
                                  max_df = NULL,
                                  default_df = 4.0) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (!length(x)) return(list(yhat = y, df = NA_real_))

  # 정렬 + 중복 x 평균 집계(LOOCV 안정화)
  ord <- order(x); x <- x[ord]; y <- y[ord]
  df2 <- stats::aggregate(y ~ x, data.frame(x = x, y = y), mean)
  x_u <- as.numeric(df2$x); y_u <- as.numeric(df2$y)

  # 유니크 시점 부족 → 상수
  if (length(x_u) < min_unique) {
    yhat <- rep(stats::mean(y_u), length(x))
    return(list(yhat = yhat[order(ord)], df = NA_real_))
  }

  # df 상한 자동
  if (is.null(max_df)) {
    max_df <- max(min(length(x_u) - 1L, 10L), min_df)
  }

  try_fit <- function(mode = c("df","spar","cv","gcv")) {
    mode <- match.arg(mode)
    if (mode == "df" && !is.null(spline_df)) {
      df_use <- min(max(spline_df, min_df), max_df)
      return(suppressWarnings(stats::smooth.spline(x_u, y_u, df = df_use)))
    }
    if (mode == "spar" && !is.null(spline_spar)) {
      return(suppressWarnings(stats::smooth.spline(x_u, y_u, spar = spline_spar)))
    }
    if (mode == "cv") {
      return(suppressWarnings(stats::smooth.spline(x_u, y_u, cv = TRUE, all.knots = TRUE)))
    }
    return(suppressWarnings(stats::smooth.spline(x_u, y_u, cv = FALSE, all.knots = TRUE)))
  }

  fit <- NULL
  for (m in c(if (!is.null(spline_df)) "df" else NULL,
              if (!is.null(spline_spar)) "spar" else NULL,
              if (isTRUE(use_cv)) "cv" else NULL,
              "gcv")) {
    fit <- try(try_fit(m), silent = TRUE)
    if (!inherits(fit, "try-error")) break
  }

  # 실패 → 기본 df
  if (inherits(fit, "try-error") || is.null(fit)) {
    fit <- try(suppressWarnings(stats::smooth.spline(
      x_u, y_u, df = min(max(default_df, min_df), max_df)
    )), silent = TRUE)
  }

  # 그래도 실패 → 상수
  if (inherits(fit, "try-error") || is.null(fit)) {
    yhat <- rep(stats::mean(y_u), length(x))
    return(list(yhat = yhat[order(ord)], df = NA_real_))
  }

  # 예측 (외삽 가드)
  yhat <- try(as.numeric(stats::predict(fit, x = x)$y), silent = TRUE)
  if (inherits(yhat, "try-error") || !length(yhat) || any(!is.finite(yhat))) {
    yhat <- as.numeric(stats::approx(
      x = x_u,
      y = as.numeric(stats::predict(fit, x = x_u)$y),
      xout = x, rule = 2
    )$y)
  }

  df_out <- suppressWarnings(as.numeric(fit$df))
  list(yhat = yhat[order(ord)], df = ifelse(is.finite(df_out), df_out, NA_real_))
}


#' @keywords internal
.log_mean_exp <- function(x) {
  xm <- max(x)
  xm + log(mean(exp(x - xm)))
}

#' @keywords internal
.sample_with_retry <- function(mod, base_args, stan_list,
                               max_retries = 3,
                               nu_fix,
                               silent_sampler = FALSE,
                               ebfmi_thresh = 0.30,
                               max_nu_cap = 7,
                               tag = "",
                               freeze_retry_hypers = FALSE) {

  # 안전 기본값 헬퍼
  .or <- function(x, y) if (is.null(x)) y else x

  # 안전 init
  .init_safe <- function(stan_list, mod) {
    function(chain_id) {
      lst <- list()
      params <- try(names(mod$variables()$parameters), silent = TRUE)
      if (inherits(params, "try-error") || is.null(params)) params <- character(0)
      add <- function(n, v) if (n %in% params) lst[[n]] <<- v

      add("r0", 0); add("a_ii", 0); add("a_ij", 0)
      add("sigma", 0.20); add("sd_ou", 0.40); add("phi", 0.80)

      lst
    }
  }
  # 재시도 필요 여부 판단 (EBFMI + div + treedepth)
  .needs_retry <- function(diag, fit = NULL, thr = 0.30) {
    low_eb <- is.finite(diag$ebfmi_min) && (diag$ebfmi_min < thr)
    if (!low_eb && !is.null(fit)) {
      eb_vec <- .ebfmi_chainwise_from_energy(fit)
      if (length(eb_vec)) low_eb <- any(is.finite(eb_vec) & (eb_vec < thr))
      if (!low_eb) low_eb <- .ebfmi_warn_from_fit(fit, thr = thr)
    }
    low_eb ||
      (is.finite(diag$n_divergent)     && diag$n_divergent > 0) ||
      (is.finite(diag$n_treedepth_hit) && diag$n_treedepth_hit > 0)
  }

  # 포맷터(로그용)
  .fmt_diag <- function(d)
    sprintf("div=%s, treedepth=%s, rhat=%s, ess_bulk=%s",
            ifelse(is.finite(d$n_divergent), d$n_divergent, NA),
            ifelse(is.finite(d$n_treedepth_hit), d$n_treedepth_hit, NA),
            ifelse(is.finite(d$worst_rhat), sprintf("%.3f", d$worst_rhat), NA),
            ifelse(is.finite(d$min_ess_bulk), d$min_ess_bulk, NA))

  # 준비
  base_args <- base_args
  base_args$data <- stan_list
  nu_fix <- .or(nu_fix, .or(stan_list$nu_fixed, 4L))

  attempt <- 0L
  used_nu <- nu_fix
  fit <- NULL; diag <- NULL; fit_failed <- FALSE
  final_args <- NULL

  repeat {
    # nu bump (freeze면 bump 안 함)
    if (attempt == 0L) {
      used_nu <- min(max_nu_cap, nu_fix)
    } else if (!freeze_retry_hypers) {
      used_nu <- min(max_nu_cap, nu_fix + min(attempt, 2L))  # +1, +2까지
    } else {
      used_nu <- nu_fix
    }
    base_args$data$nu_fixed <- used_nu

    sample_args <- base_args
    if (attempt > 0L) {
      if (!freeze_retry_hypers) {
        sample_args$adapt_delta <- max(0.995, .or(base_args$adapt_delta, 0))
        sample_args$iter_warmup <- max(2000,  .or(base_args$iter_warmup, 1000))
        sample_args$metric      <- "dense_e"
        sample_args$step_size   <- 0.03
        if (is.numeric(base_args$init) && length(base_args$init) == 1L && is.finite(base_args$init)) {
          sample_args$init <- base_args$init
        } else {
          sample_args$init <- NULL  # 전 파라미터 자동 초기화
        }
      }
      # 리트라이에서도 seed 고정
      sample_args$seed <- base_args$seed
    }

    if (is.list(sample_args$init) && !is.function(sample_args$init)) {
      # Resolve number of chains (prefer `chains`, then fallback to `parallel_chains`)
      as_pos_int1 <- function(x) {
        if (is.null(x)) return(NA_integer_)
        x <- suppressWarnings(as.numeric(x))
        if (length(x) == 1L && is.finite(x) && !is.na(x) && x >= 1 && floor(x) == x) {
          return(as.integer(x))
        }
        NA_integer_
      }
      n <- as_pos_int1(sample_args$chains)
      if (is.na(n)) {
        n <- as_pos_int1(sample_args$parallel_chains)
        # keep internal consistency if chains was missing
        if (!is.na(n) && is.null(sample_args$chains)) sample_args$chains <- n
      }
      if (is.na(n)) n <- 1L

      init_obj <- sample_args$init
      # If init is already list-of-lists (per-chain), validate length.
      # Else, replicate a single parameter list across chains.
      if (length(init_obj) > 0L && is.list(init_obj[[1L]])) {
        if (length(init_obj) != n) {
          stop("`init` length (", length(init_obj), ") must equal number of chains (", n, ").")
        }
        # leave as-is
      } else {
        sample_args$init <- rep(list(init_obj), n)
      }
    }

    message(sprintf("⏩ %sattempt %d/%d | nu_fixed=%d | metric=%s, adapt_delta=%s, warmup=%s, step_size=%s",
                    if (nzchar(tag)) paste0("[",tag,"] ") else "",
                    attempt, max_retries, used_nu,
                    .or(sample_args$metric, "NA"),
                    ifelse(is.null(sample_args$adapt_delta),"NA",sprintf("%.3f", sample_args$adapt_delta)),
                    .or(sample_args$iter_warmup, "NA"),
                    .or(sample_args$step_size, "NA")))

    fit <- .call_sample_silently(mod, sample_args, silent = silent_sampler)

    diag <- summarise_diag(
      fit,
      max_treedepth = .or(sample_args$max_treedepth, base_args$max_treedepth)
    )

    final_args <- sample_args

    eb_vec <- .ebfmi_chainwise_from_energy(fit)
    eb_str <- if (length(eb_vec)) paste(sprintf("%.3f", eb_vec), collapse=",") else "NA"
    ok <- !.needs_retry(diag, fit, thr = ebfmi_thresh)

    message(sprintf("%s%s | %s | E-BFMI chains=[%s]",
                    if (nzchar(tag)) paste0("[",tag,"] ") else "",
                    if (ok) "✅ diag ok" else "⚠️ diag warn",
                    .fmt_diag(diag), eb_str))

    if (ok || attempt >= max_retries) { if (!ok) fit_failed <- TRUE; break }
    attempt <- attempt + 1L
  }

  list(
    fit = fit,
    diag = diag,
    n_retries = attempt,
    used_nu = used_nu,
    fit_failed = fit_failed,
    final_args = final_args
  )
}

#' @keywords internal
.run_pair <- function(idx_i, idx_j, mute_logs = FALSE, seed_base) {
  # 인덱스 → taxon 이름
  ti <- taxa_vec[[idx_i]]
  pj <- taxa_vec[[idx_j]]
  # 쌍별 seed 파생 (항상 seed_base 사용; 전역 참조 없음)
  seed_ij <- as.integer(seed_base + 100000L * idx_i + 1000L * idx_j + 1L)
  seed_ji <- as.integer(seed_base + 100000L * idx_i + 1000L * idx_j + 2L)

  res_ij <- run_one(
    target = ti, partner = pj,
    seed_override = seed_ij,
    progress_local = if (mute_logs) "none" else progress
  )
  res_ji <- run_one(
    target = pj, partner = ti,
    seed_override = seed_ji,
    progress_local = if (mute_logs) "none" else progress
  )

  # 둘 중 하나라도 완전 NULL이면 스킵
  if (is.null(res_ij) || is.null(res_ji)) return(NULL)

  tibble::as_tibble_row(list(
    i = ti, j = pj,
    # n_pairs
    n_pairs_ij = res_ij$n_pairs, n_pairs_ji = res_ji$n_pairs,
    # directional effects
    a_ij_mean = res_ij$a_mean, a_ij_sd = res_ij$a_sd,
    a_ij_q2.5 = res_ij$a_q2.5, a_ij_q97.5 = res_ij$a_q97.5,
    p_sign2_ij = res_ij$p_sign2,
    a_ji_mean = res_ji$a_mean, a_ji_sd = res_ji$a_sd,
    a_ji_q2.5 = res_ji$a_q2.5, a_ji_q97.5 = res_ji$a_q97.5,
    p_sign2_ji = res_ji$p_sign2,
    # self-effects
    a_ii_mean = res_ij$aii_mean, a_ii_sd = res_ij$aii_sd,
    a_ii_q2.5 = res_ij$aii_q2.5, a_ii_q97.5 = res_ij$aii_q97.5,
    p_sign2_ii = res_ij$p_sign2_self,
    a_jj_mean = res_ji$aii_mean, a_jj_sd = res_ji$aii_sd,
    a_jj_q2.5 = res_ji$aii_q2.5, a_jj_q97.5 = res_ji$aii_q97.5,
    p_sign2_jj = res_ji$p_sign2_self,
    # diagnostics
    rhat_ij  = res_ij$diag$worst_rhat, essb_ij = res_ij$diag$min_ess_bulk,
    esst_ij  = res_ij$diag$min_ess_tail, div_ij = res_ij$diag$n_divergent,
    tdhit_ij = res_ij$diag$n_treedepth_hit,
    ebfmi_min_ij = res_ij$diag$ebfmi_min, ebfmi_med_ij = res_ij$diag$ebfmi_med,
    rhat_ji  = res_ji$diag$worst_rhat, essb_ji = res_ji$diag$min_ess_bulk,
    esst_ji  = res_ji$diag$min_ess_tail, div_ji = res_ji$diag$n_divergent,
    tdhit_ji = res_ji$diag$n_treedepth_hit,
    ebfmi_min_ji = res_ji$diag$ebfmi_min, ebfmi_med_ji = res_ji$diag$ebfmi_med,
    # Repeated K-fold 요약 (mean)
    kfold_elpd_mean_ij = res_ij$kfold_mean,
    kfold_elpd_mean_ji = res_ji$kfold_mean,
    kfold_K = kfold_K, kfold_R = kfold_R,
    # pseudo-BMA/stacking용 payload (list-columns)
    kfold_subject_ij         = list(res_ij$kfold_subject),
    kfold_subject_ids_ij     = list(res_ij$kfold_subject_ids),
    kfold_subject_counts_ij  = list(res_ij$kfold_subject_counts),
    kfold_splits_ij          = list(res_ij$kfold_splits),
    kfold_seed_used_ij       = res_ij$kfold_seed_used,
    kfold_sd_ij              = res_ij$kfold_sd,
    kfold_se_ij              = res_ij$kfold_se,
    kfold_n_subjects_ij      = res_ij$kfold_n_subjects,
    kfold_retry_total_ij     = res_ij$kfold_retry_total,
    kfold_retry_mean_ij      = res_ij$kfold_retry_mean,
    kfold_nu_used_counts_ij  = list(res_ij$kfold_nu_used_counts),
    kfold_outer_rounds_ij    = res_ij$kfold_outer_rounds,
    kfold_failed_ij          = res_ij$kfold_failed,
    kfold_folds_ok_ij        = res_ij$kfold_n_folds_ok,
    kfold_folds_fail_ij      = res_ij$kfold_n_folds_fail,
    #
    kfold_subject_ji         = list(res_ji$kfold_subject),
    kfold_subject_ids_ji     = list(res_ji$kfold_subject_ids),
    kfold_subject_counts_ji  = list(res_ji$kfold_subject_counts),
    kfold_splits_ji          = list(res_ji$kfold_splits),
    kfold_seed_used_ji       = res_ji$kfold_seed_used,
    kfold_sd_ji              = res_ji$kfold_sd,
    kfold_se_ji              = res_ji$kfold_se,
    kfold_n_subjects_ji      = res_ji$kfold_n_subjects,
    kfold_retry_total_ji     = res_ji$kfold_retry_total,
    kfold_retry_mean_ji      = res_ji$kfold_retry_mean,
    kfold_nu_used_counts_ji  = list(res_ji$kfold_nu_used_counts),
    kfold_outer_rounds_ji    = res_ji$kfold_outer_rounds,
    kfold_failed_ji          = res_ji$kfold_failed,
    kfold_folds_ok_ji        = res_ji$kfold_n_folds_ok,
    kfold_folds_fail_ji      = res_ji$kfold_n_folds_fail
  ))
}

