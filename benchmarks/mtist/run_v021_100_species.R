stopifnot(file.exists("DESCRIPTION"))

repo <- normalizePath(getwd())

# -------------------------------------------------------------------------
# Run configuration: 여기만 수정하면 benchmark와 실제 실행이 함께 바뀜
# -------------------------------------------------------------------------

dataset_id <- 1L

run_seed <- 20260802L

chains <- 4L
iter_warmup <- 2000L
iter_sampling <- 1000L

init <- 0.2
adapt_delta <- 0.98
max_treedepth <- 14L

n_workers_outer <- 15L

kfold_K <- 5L
kfold_R <- 1L
kfold_seed <- run_seed

nz_partner_min_frac <- 0.15
min_unique_times <- 3L
min_pairs <- 4L


# -------------------------------------------------------------------------
# Install current package source into isolated library
# -------------------------------------------------------------------------

lib <- file.path(
  repo,
  "benchmarks",
  "mtist",
  ".lib-v021-installed"
)

unlink(lib, recursive = TRUE, force = TRUE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

# worker / CmdStan / BLAS의 암묵적 다중 스레딩 방지
Sys.setenv(
  R_LIBS_USER = lib,
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

build_dir <- tempfile("pclv-build-")
dir.create(build_dir)

pkg_tar <- devtools::build(
  pkg = repo,
  path = build_dir,
  binary = FALSE,
  vignettes = FALSE,
  manual = FALSE
)

install.packages(
  pkg_tar,
  repos = NULL,
  type = "source",
  lib = lib,
  INSTALL_opts = c("--preclean", "--clean")
)

unlink(build_dir, recursive = TRUE, force = TRUE)

.libPaths(c(lib, .libPaths()))

library(pclvbayes, lib.loc = lib)
library(phyloseq)

cat(
  "package version: ",
  as.character(packageVersion("pclvbayes", lib.loc = lib)),
  "\n",
  sep = ""
)

cat(
  "package path:    ",
  find.package("pclvbayes", lib.loc = lib),
  "\n",
  sep = ""
)


# -------------------------------------------------------------------------
# Load MTIST dataset
# -------------------------------------------------------------------------

mtist_root <- normalizePath(
  path.expand(Sys.getenv("MTIST_ROOT", unset = "~/mtist")),
  mustWork = TRUE
)

dataset_csv <- file.path(
  mtist_root,
  "mtist1.0",
  "mtist_datasets",
  sprintf("dataset_%d.csv", dataset_id)
)

stopifnot(file.exists(dataset_csv))

raw <- read.csv(
  dataset_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_meta <- c("timeseries_id", "time")
stopifnot(all(required_meta %in% names(raw)))

species <- names(raw)[
  grepl("^species_[0-9]+$", names(raw))
]

stopifnot(
  length(species) == 100L,
  !anyDuplicated(species)
)

# subject별 시간순 정렬
raw <- raw[
  order(raw$timeseries_id, raw$time),
  ,
  drop = FALSE
]

abundance <- as.matrix(raw[species])
storage.mode(abundance) <- "double"

stopifnot(
  all(is.finite(abundance)),
  all(abundance >= 0),
  all(rowSums(abundance) > 0),
  all(is.finite(raw$time)),
  !anyDuplicated(
    paste(raw$timeseries_id, raw$time, sep = "\r")
  )
)

# MTIST absolute density -> relative abundance
relative <- abundance / rowSums(abundance)

sample_index <- ave(
  seq_len(nrow(raw)),
  raw$timeseries_id,
  FUN = seq_along
)

sample_id <- sprintf(
  "did%d_ts%s_n%03d",
  dataset_id,
  raw$timeseries_id,
  sample_index
)

rownames(relative) <- sample_id

sample_meta <- data.frame(
  subject = as.character(raw$timeseries_id),
  time = as.numeric(raw$time),
  dataset_id = dataset_id,
  row.names = sample_id,
  stringsAsFactors = FALSE
)

taxonomy <- matrix(
  species,
  ncol = 1L,
  dimnames = list(species, "Taxon")
)

ps100 <- phyloseq::phyloseq(
  phyloseq::otu_table(
    t(relative),
    taxa_are_rows = TRUE
  ),
  phyloseq::sample_data(sample_meta),
  phyloseq::tax_table(taxonomy)
)

stopifnot(
  phyloseq::ntaxa(ps100) == 100L,
  phyloseq::nsamples(ps100) == nrow(raw),
  identical(phyloseq::taxa_names(ps100), species)
)


# -------------------------------------------------------------------------
# Output / provenance
# -------------------------------------------------------------------------

out_root <- file.path(
  repo,
  "benchmarks",
  "mtist",
  "results",
  "v021_installed_package_100_species_pairparallel_v2"
)

if (
  dir.exists(out_root) &&
  length(list.files(out_root, all.files = TRUE, no.. = TRUE))
) {
  stop("Output root is already non-empty: ", out_root)
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

git_commit <- trimws(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
)

configuration <- list(
  execution = "installed-package-public-api",

  dataset_id = dataset_id,
  dataset_path = normalizePath(dataset_csv),
  dataset_md5 = unname(tools::md5sum(dataset_csv)),
  taxa_count = phyloseq::ntaxa(ps100),
  sample_count = phyloseq::nsamples(ps100),

  nz_partner_min_frac = nz_partner_min_frac,
  min_unique_times = min_unique_times,
  min_pairs = min_pairs,

  seed = run_seed,
  init = init,
  chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,

  # Parallel architecture
  parallel_unit = "unordered_pair",
  n_workers_outer = n_workers_outer,
  parallel_chains = 1L,

  kfold_K = kfold_K,
  kfold_R = kfold_R,
  kfold_seed = kfold_seed,

  package_commit = git_commit,
  package_version = as.character(
    packageVersion("pclvbayes", lib.loc = lib)
  ),
  package_path = find.package(
    "pclvbayes",
    lib.loc = lib
  ),
  cmdstanr_version = as.character(
    packageVersion("cmdstanr")
  ),
  cmdstan_version = as.character(
    cmdstanr::cmdstan_version()
  )
)

saveRDS(
  ps100,
  file.path(out_root, "input_phyloseq.rds")
)

saveRDS(
  configuration,
  file.path(out_root, "configuration.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_root, "sessionInfo.txt")
)

print(configuration)


# -------------------------------------------------------------------------
# Run fit
# -------------------------------------------------------------------------

status_path <- file.path(
  out_root,
  "status.rds"
)

result_path <- file.path(
  out_root,
  "fit_pclv_bayes_result.rds"
)

started_at <- Sys.time()

saveRDS(
  list(
    state = "running",
    started_at = started_at,
    configuration = configuration
  ),
  status_path
)

fit100 <- tryCatch(
  {
    pclvbayes::fit_pclv_bayes(
      physeq = ps100,
      subject_col = "subject",
      time_col = "time",
      taxa_vec = species,

      nz_partner_min_frac = nz_partner_min_frac,
      min_unique_times = min_unique_times,
      min_pairs = min_pairs,

      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,

      seed = run_seed,
      init = init,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,

      progress = "bar",

      n_workers_outer = n_workers_outer,

      kfold_K = kfold_K,
      kfold_R = kfold_R,
      kfold_seed = kfold_seed
    )
  },

  error = function(e) {
    failed_at <- Sys.time()

    saveRDS(
      list(
        state = "failed",
        started_at = started_at,
        failed_at = failed_at,
        elapsed_seconds = as.numeric(
          difftime(
            failed_at,
            started_at,
            units = "secs"
          )
        ),
        error_class = class(e),
        error_message = conditionMessage(e),
        calls = sys.calls(),
        configuration = configuration
      ),
      status_path
    )

    stop(e)
  }
)


# -------------------------------------------------------------------------
# Preserve result
# -------------------------------------------------------------------------

completed_at <- Sys.time()

saveRDS(
  fit100,
  result_path,
  compress = FALSE
)

saveRDS(
  list(
    state = "completed",
    started_at = started_at,
    completed_at = completed_at,
    elapsed_seconds = as.numeric(
      difftime(
        completed_at,
        started_at,
        units = "secs"
      )
    ),
    result_path = result_path,
    result_names = names(fit100),
    configuration = configuration
  ),
  status_path
)

cat(
  "Completed in ",
  round(
    as.numeric(
      difftime(
        completed_at,
        started_at,
        units = "hours"
      )
    ),
    3
  ),
  " hours\n",
  sep = ""
)

print(names(fit100))
