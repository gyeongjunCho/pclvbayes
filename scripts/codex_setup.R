required <- c(
  "devtools",
  "testthat",
  "roxygen2",
  "cmdstanr",
  "posterior",
  "loo"
)

installed <- rownames(installed.packages())
missing <- setdiff(required, installed)

if (length(missing) > 0L) {
  install.packages(
    missing,
    repos = c(
      CRAN = "https://cloud.r-project.org",
      STAN = "https://mc-stan.org/r-packages/"
    )
  )
}

devtools::install_deps(
  dependencies = TRUE,
  upgrade = "never"
)
