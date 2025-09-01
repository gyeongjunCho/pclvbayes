.onLoad <- function(libname, pkgname) {
  # Be quiet & never error on load (CRAN friendly)
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(invisible())

  ok <- FALSE
  ver <- try(cmdstanr::cmdstan_version(), silent = TRUE)  # no args
  if (!inherits(ver, "try-error") && !is.na(ver)) ok <- TRUE

  # Only nudge interactive users; never stop() on load
  if (!ok && interactive()) {
    packageStartupMessage(
      "CmdStan not found. Install with: cmdstanr::install_cmdstan()"
    )
  }
}


# R/zzz.R (하단에 추가)

utils::globalVariables(c(
  # dplyr 파이프라인 내 컬럼/심볼
  ".data",
  "xi_raw","xj_raw","xj_next","dt","rest_now","rest_next",
  "alr_i_next","alr_i_now","y",
  "i","j","k",
  "rhat_ij","essb_ij","esst_ij","div_ij","tdhit_ij","ok_diag_ij","keep_ij_final",
  "rhat_ji","essb_ji","esst_ji","div_ji","tdhit_ji","ok_diag_ji","keep_ji_final",
  "q_bayes","p_sign2","r_pooled",
  # 파이프라인 외부에서 NOTE로 잡힌 심볼들
  "taxa_vec","progress","kfold_K","kfold_R"
))
