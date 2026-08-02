mtist_root <- function(path = NULL) {
  root <- normalizePath(if (is.null(path)) "~/mtist" else path,
                        mustWork = FALSE)
  if (!dir.exists(root)) stop("MTIST repository not found: ", root)
  root
}

mtist_git_sha <- function(root = mtist_root()) {
  sha <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"),
                 stdout = TRUE, stderr = TRUE)
  if (!length(sha) || attr(sha, "status") %||% 0L != 0L)
    stop("Could not determine the MTIST Git SHA.")
  trimws(sha[[1L]])
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.mtist_paths <- function(root) {
  base <- file.path(root, "mtist1.0")
  list(
    datasets = file.path(base, "mtist_datasets"),
    metadata = file.path(base, "mtist_datasets", "mtist_metadata.csv"),
    truths = file.path(base, "ground_truths", "interaction_coefficients")
  )
}

.drop_confirmed_csv_index <- function(x, path) {
  if (!length(names(x)) || !identical(names(x)[[1L]], ""))
    stop("Expected a leading unnamed CSV index column in ", path)
  index <- x[[1L]]
  if (!is.numeric(index) || anyNA(index) || any(!is.finite(index)))
    stop("The unnamed CSV index column is malformed in ", path)
  x[-1L]
}

.read_indexed_csv <- function(path) {
  if (!file.exists(path)) stop("Missing MTIST file: ", path)
  .drop_confirmed_csv_index(utils::read.csv(path, check.names = FALSE), path)
}

enumerate_mtist_studies <- function(n_species, root = mtist_root()) {
  paths <- .mtist_paths(root)
  meta <- .read_indexed_csv(paths$metadata)
  meta <- meta[meta$n_species == as.integer(n_species), , drop = FALSE]
  meta$dataset_path <- file.path(paths$datasets, paste0("dataset_", meta$did, ".csv"))
  meta$truth_path <- file.path(
    paths$truths,
    paste0(sub("_gt_", "_aij_", meta$ground_truth, fixed = TRUE), ".csv")
  )
  meta$compatible <- file.exists(meta$dataset_path) & file.exists(meta$truth_path) &
    meta$n_timeseries >= 2L
  meta
}


enumerate_mtist_3species <- function(root = mtist_root())
  enumerate_mtist_studies(3L, root)

enumerate_mtist_10species <- function(root = mtist_root())
  enumerate_mtist_studies(10L, root)

.select_mtist_10species_catalog <- function(candidates) {
  eligible <- candidates[candidates$compatible & candidates$noise == 0.01 &
    candidates$sampling_scheme == "even" & candidates$n_timeseries >= 10L &
    candidates$n_timepoints >= 15L, , drop = FALSE]
  if (!nrow(eligible)) stop("No compatible representative 10-species study found.")
  eligible[order(eligible$did), , drop = FALSE][1L, , drop = FALSE]
}

select_mtist_10species <- function(root = mtist_root())
  .select_mtist_10species_catalog(enumerate_mtist_10species(root))

load_mtist_study <- function(dataset_id, root = mtist_root()) {
  paths <- .mtist_paths(root)
  metadata <- .read_indexed_csv(paths$metadata)
  metadata_record <- metadata[metadata$did == dataset_id, , drop = FALSE]
  if (nrow(metadata_record) != 1L) stop("Dataset ID must identify one MTIST study: ", dataset_id)
  catalog <- enumerate_mtist_studies(metadata_record$n_species[[1L]], root)
  record <- catalog[catalog$did == dataset_id, , drop = FALSE]
  if (nrow(record) != 1L) stop("Dataset ID must identify one MTIST study: ", dataset_id)
  if (!isTRUE(record$compatible[[1L]])) stop("Dataset is not compatible: ", dataset_id)

  data <- .read_indexed_csv(record$dataset_path[[1L]])
  species <- names(data)[grepl("^species_[0-9]+$", names(data))]
  if (length(species) != record$n_species[[1L]]) stop("Species-column count disagrees with metadata.")

  fields <- c("did", "n_species", "ground_truth", "noise", "n_timeseries",
              "n_timepoints", "sampling_scheme")
  for (field in fields) {
    values <- unique(data[[field]])
    if (length(values) != 1L || is.na(values)) stop("Dataset field is not constant: ", field)
    expected <- record[[field]][[1L]]
    if (is.numeric(expected)) {
      if (!isTRUE(all.equal(as.numeric(values), as.numeric(expected))))
        stop("Dataset field disagrees with metadata: ", field)
    } else if (!identical(as.character(values), as.character(expected))) {
      stop("Dataset field disagrees with metadata: ", field)
    }
  }

  abundance <- as.matrix(data[, species, drop = FALSE])
  storage.mode(abundance) <- "double"
  if (any(!is.finite(abundance))) stop("Non-finite abundance values detected.")
  if (any(abundance < 0)) stop("Negative abundance values detected.")
  totals <- rowSums(abundance)
  if (any(!is.finite(totals) | totals <= 0)) stop("Non-positive or non-finite abundance row sum detected.")
  relative <- abundance / totals
  if (any(abs(rowSums(relative) - 1) > 1e-12)) stop("Relative-abundance normalization failed.")

  if (!is.numeric(data$time) || any(!is.finite(data$time))) stop("Invalid numeric time values.")
  if (anyNA(data$timeseries_id)) stop("Missing timeseries_id values.")
  duplicate_time <- duplicated(data[, c("timeseries_id", "time")])
  if (any(duplicate_time)) stop("Duplicate times within a time series.")
  observed_counts <- table(data$timeseries_id)
  if (length(observed_counts) != record$n_timeseries[[1L]] ||
      any(observed_counts != record$n_timepoints[[1L]]))
    stop("Observed time-series dimensions disagree with metadata.")

  sample_id <- sprintf("did%s_ts%s_n%03d", dataset_id, data$timeseries_id,
                       ave(seq_len(nrow(data)), data$timeseries_id, FUN = seq_along))
  if (anyDuplicated(sample_id)) stop("Generated sample IDs are not unique.")
  rownames(relative) <- sample_id
  sample_meta <- data.frame(
    subject = as.character(data$timeseries_id), time = as.numeric(data$time),
    dataset_id = as.integer(data$did), truth_name = as.character(data$ground_truth),
    noise = as.numeric(data$noise), sampling_scheme = as.character(data$sampling_scheme),
    row.names = sample_id, check.names = FALSE
  )
  otu <- t(relative)
  if (!identical(rownames(otu), species) || !identical(colnames(otu), rownames(sample_meta)))
    stop("OTU and metadata names are inconsistent.")
  physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::sample_data(sample_meta)
  )
  truth <- as.matrix(utils::read.csv(record$truth_path[[1L]], header = FALSE,
                                     check.names = FALSE))
  storage.mode(truth) <- "double"
  if (!identical(dim(truth), c(length(species), length(species))) || any(!is.finite(truth)))
    stop("Truth matrix has invalid dimensions or values.")
  dimnames(truth) <- list(species, species)

  list(physeq = physeq, absolute = abundance, relative = relative,
       sample_metadata = sample_meta, taxa = species, truth = truth,
       record = record[1L, , drop = FALSE])
}

validate_mtist_prediction <- function(x, taxa) {
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(length(taxa), length(taxa))))
    stop("Prediction matrix must be a square numeric matrix in MTIST taxon order.")
  if (any(!is.finite(x))) stop("Prediction matrix contains non-finite values.")
  if (!identical(rownames(x), taxa) || !identical(colnames(x), taxa))
    stop("Prediction matrix taxon order is incorrect.")
  invisible(TRUE)
}
