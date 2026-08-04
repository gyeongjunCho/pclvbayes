#' @noRd
#' @keywords internal
.as_int <- function(x) suppressWarnings(as.integer(x))

#' @noRd
#' @keywords internal
.as_num <- function(x) suppressWarnings(as.numeric(x))

#' @noRd
#' @keywords internal
.lfsr_safe <- function(p_two) {
  if (exists(".lfsr_from_two", mode = "function", inherits = TRUE)) {
    return(.lfsr_from_two(p_two))
  }
  pmin(pmax(p_two / 2, 0), 0.5)
}

#' Element-wise diagnostic pass/fail predicate used by the public summarizer
#'
#' Missing, non-finite, negative, non-integral, or malformed diagnostic counts
#' fail closed. A diagnostic is never treated as passing merely because its
#' divergence or treedepth count is unavailable.
#'
#' @noRd
#' @keywords internal
.diag_ok_fun <- function(rhat, essb, esst, div, tdhit, ebfmi_min, thr) {
  rhat_num <- .as_num(rhat)
  essb_num <- .as_num(essb)
  esst_num <- .as_num(esst)
  div_num <- .as_num(div)
  tdhit_num <- .as_num(tdhit)
  ebfmi_num <- .as_num(ebfmi_min)

  div_int <- .as_int(div_num)
  tdhit_int <- .as_int(tdhit_num)

  valid_div <-
    is.finite(div_num) &
    !is.na(div_int) &
    div_num >= 0 &
    div_num == div_int &
    div_int <= thr$div

  valid_tdhit <-
    is.finite(tdhit_num) &
    !is.na(tdhit_int) &
    tdhit_num >= 0 &
    tdhit_num == tdhit_int &
    tdhit_int <= thr$tdhit

  (is.finite(rhat_num) & rhat_num < thr$rhat) &
    (is.finite(essb_num) & essb_num > thr$ess) &
    (is.finite(esst_num) & esst_num > thr$ess) &
    valid_div &
    valid_tdhit &
    (is.finite(ebfmi_num) & ebfmi_num >= thr$ebfmi)
}
