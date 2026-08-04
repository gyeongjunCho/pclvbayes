# R/stan_model.R
#' @keywords internal
.pGLV_env <- new.env(parent = emptyenv())

#' @keywords internal
#' Load and cache the canonical pcLV CmdStan model.
get_pclv_model <- function(quiet = TRUE, rebuild = FALSE) {
  if (!is.logical(quiet) ||
      length(quiet) != 1L ||
      is.na(quiet)) {
    stop(
      "`quiet` must be a non-missing scalar logical value.",
      call. = FALSE
    )
  }

  if (!is.logical(rebuild) ||
      length(rebuild) != 1L ||
      is.na(rebuild)) {
    stop(
      "`rebuild` must be a non-missing scalar logical value.",
      call. = FALSE
    )
  }

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Package 'cmdstanr' is required.",
      call. = FALSE
    )
  }

  package_name <- tryCatch(
    utils::packageName(),
    error = function(e) NULL
  )

  if (is.null(package_name) ||
      length(package_name) != 1L ||
      is.na(package_name) ||
      !nzchar(package_name)) {
    package_name <- "pclvbayes"
  }

  stan_file <- system.file(
    "stan",
    "pclv.stan",
    package = package_name,
    mustWork = FALSE
  )

  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stop(
      sprintf(
        paste(
          "Canonical Stan source was not found:",
          "inst/stan/pclv.stan in package '%s'."
        ),
        package_name
      ),
      call. = FALSE
    )
  }

  stan_file <- normalizePath(
    stan_file,
    winslash = "/",
    mustWork = TRUE
  )

  cached_model_is_valid <- function(model) {
    if (is.null(model)) {
      return(FALSE)
    }

    exe_file <- tryCatch(
      model$exe_file(),
      error = function(e) NULL
    )

    if (!is.character(exe_file) ||
        length(exe_file) != 1L ||
        is.na(exe_file) ||
        !nzchar(exe_file) ||
        !file.exists(exe_file)) {
      return(FALSE)
    }

    cached_stan_file <- tryCatch(
      model$stan_file(),
      error = function(e) NULL
    )

    # Older CmdStanR objects may not expose stan_file(). In that case,
    # executable existence remains the available cache-validity criterion.
    if (is.null(cached_stan_file) ||
        !is.character(cached_stan_file) ||
        length(cached_stan_file) != 1L ||
        is.na(cached_stan_file) ||
        !nzchar(cached_stan_file)) {
      return(TRUE)
    }

    cached_stan_file <- tryCatch(
      normalizePath(
        cached_stan_file,
        winslash = "/",
        mustWork = TRUE
      ),
      error = function(e) NULL
    )

    !is.null(cached_stan_file) &&
      identical(cached_stan_file, stan_file)
  }

  cached_model <- if (exists(
    "mod",
    envir = .pGLV_env,
    inherits = FALSE
  )) {
    .pGLV_env$mod
  } else {
    NULL
  }

  if (!isTRUE(rebuild) &&
      cached_model_is_valid(cached_model)) {
    return(cached_model)
  }

  cmdstan_version <- tryCatch(
    cmdstanr::cmdstan_version(),
    error = function(e) NULL
  )

  if (is.null(cmdstan_version) ||
      !length(cmdstan_version) ||
      anyNA(cmdstan_version)) {
    stop(
      paste(
        "CmdStan was not found.",
        "Run cmdstanr::install_cmdstan() and try again."
      ),
      call. = FALSE
    )
  }

  # cmdstan_model() performs the compile itself. Passing force_recompile here
  # avoids constructing once and then compiling the same model a second time.
  model <- tryCatch(
    cmdstanr::cmdstan_model(
      stan_file = stan_file,
      compile = TRUE,
      force_recompile = isTRUE(rebuild),
      quiet = quiet
    ),
    error = function(e) {
      stop(
        sprintf(
          "Canonical pcLV Stan model compilation failed: %s",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  exe_file <- tryCatch(
    model$exe_file(),
    error = function(e) NULL
  )

  if (!is.character(exe_file) ||
      length(exe_file) != 1L ||
      is.na(exe_file) ||
      !nzchar(exe_file) ||
      !file.exists(exe_file)) {
    stop(
      paste(
        "Canonical pcLV model was constructed,",
        "but its compiled executable is unavailable."
      ),
      call. = FALSE
    )
  }

  # Replace the cache only after a fully valid model has been obtained.
  .pGLV_env$mod <- model
  model
}
