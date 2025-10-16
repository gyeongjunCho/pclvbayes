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


