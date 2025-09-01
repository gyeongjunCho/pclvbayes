# R/stan_model.R
#' @keywords internal
.pGLV_env <- new.env(parent = emptyenv())

#' @keywords internal
#' Load & cache the single AR(1) cmdstan model
get_glv_model <- function(quiet = TRUE, rebuild = FALSE) {
  if (!is.null(.pGLV_env$mod) && !isTRUE(rebuild)) return(.pGLV_env$mod)

  pkg <- tryCatch(utils::packageName(), error = function(e) "pglvbayes")
  stan_file <- system.file("stan", "glv_pairwise.stan", package = pkg, mustWork = TRUE)

  has_cmdstan <- !is.na(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NA))
  if (!has_cmdstan) {
    stop("CmdStan not found. Please run cmdstanr::install_cmdstan() and try again.", call. = FALSE)
  }

  mod <- cmdstanr::cmdstan_model(stan_file, quiet = quiet)
  if (isTRUE(rebuild)) mod$compile(force_recompile = TRUE, quiet = quiet)
  .pGLV_env$mod <- mod
  mod
}
