#' @noRd
#' @keywords internal
.as_int <- function(x) suppressWarnings(as.integer(x))

#' @noRd
#' @keywords internal
.lfsr_safe <- function(p_two) {
  if (exists(".lfsr_from_two", mode = "function", inherits = TRUE)) {
    return(.lfsr_from_two(p_two))
  }
  pmin(pmax(p_two / 2, 0), 0.5)
}

#' @noRd
#' @keywords internal
.diag_ok_fun <- function(rhat, essb, esst, div, tdhit, ebfmi_min, thr) {
  (is.finite(rhat) & rhat < thr$rhat) &
    (is.finite(essb) & essb > thr$ess) &
    (is.finite(esst) & esst > thr$ess) &
    (is.na(div) | .as_int(div) <= thr$div) &
    (is.na(tdhit) | .as_int(tdhit) <= thr$tdhit) &
    (is.finite(ebfmi_min) & ebfmi_min >= thr$ebfmi)
}
