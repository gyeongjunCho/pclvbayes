# Central boundary validator for fit_pclv_bayes().
.validate_fit_pclv_inputs <- function(
  physeq,
  subject_col,
  time_col,
  taxa_vec = NULL,
  controls
) {
  chr1 <- function(x, name) {
    if (!is.character(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !nzchar(x)) {
      stop(
        sprintf(
          "`%s` must be a non-missing scalar character value.",
          name
        ),
        call. = FALSE
      )
    }
    x
  }

  log1 <- function(x, name) {
    if (!is.logical(x) ||
        length(x) != 1L ||
        is.na(x)) {
      stop(
        sprintf(
          "`%s` must be a non-missing scalar logical value.",
          name
        ),
        call. = FALSE
      )
    }
    x
  }

  int1 <- function(
    x,
    name,
    lower = 1,
    upper = .Machine$integer.max
  ) {
    valid <-
      is.numeric(x) &&
      length(x) == 1L &&
      is.finite(x) &&
      x == floor(x) &&
      x >= lower &&
      x <= upper

    if (!valid) {
      stop(
        sprintf(
          "`%s` must be a finite integer between %s and %s.",
          name,
          format(lower, scientific = FALSE),
          format(upper, scientific = FALSE)
        ),
        call. = FALSE
      )
    }

    as.integer(x)
  }

  num1 <- function(
    x,
    name,
    lower = -Inf,
    upper = Inf,
    open_lower = FALSE,
    open_upper = FALSE
  ) {
    valid <-
      is.numeric(x) &&
      length(x) == 1L &&
      is.finite(x)

    if (valid) {
      valid <- if (open_lower) x > lower else x >= lower
    }
    if (valid) {
      valid <- if (open_upper) x < upper else x <= upper
    }

    if (!valid) {
      stop(
        sprintf(
          "`%s` must be a finite scalar in the required range.",
          name
        ),
        call. = FALSE
      )
    }

    as.numeric(x)
  }

  choice1 <- function(x, name, choices) {
    # Preserve ordinary match.arg()-style defaults such as
    # progress = c("bar", "verbose", "none").
    if (length(x) > 1L) {
      if (!identical(x, choices)) {
        stop(
          sprintf("`%s` must be a scalar choice.", name),
          call. = FALSE
        )
      }
      x <- choices[[1L]]
    }

    x <- chr1(x, name)

    if (!(x %in% choices)) {
      stop(
        sprintf(
          "`%s` must be one of: %s.",
          name,
          paste(choices, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    x
  }

  validate_named_init <- function(x, label) {
    valid_names <-
      is.list(x) &&
      length(x) > 0L &&
      !is.null(names(x)) &&
      !anyNA(names(x)) &&
      all(nzchar(names(x))) &&
      !anyDuplicated(names(x))

    if (!valid_names) {
      stop(
        sprintf(
          "`%s` must be a non-empty named parameter list.",
          label
        ),
        call. = FALSE
      )
    }

    valid_values <- vapply(
      x,
      function(value) {
        is.numeric(value) &&
          length(value) > 0L &&
          all(is.finite(value))
      },
      logical(1)
    )

    if (!all(valid_values)) {
      bad <- names(x)[!valid_values]
      stop(
        sprintf(
          "`%s` contains invalid, empty, or non-finite parameter values: %s.",
          label,
          paste(bad, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  # ---------------------------------------------------------------------------
  # Control-list contract
  # ---------------------------------------------------------------------------

  if (!is.list(controls) ||
      is.null(names(controls)) ||
      anyNA(names(controls)) ||
      any(!nzchar(names(controls))) ||
      anyDuplicated(names(controls))) {
    stop(
      "`controls` must be a uniquely named list.",
      call. = FALSE
    )
  }

  required_controls <- c(
    "eps",
    "min_pairs",
    "min_unique_times",
    "zero_mode_alr",
    "minpos_alpha",
    "minpos_base",
    "eps_fixed",
    "lib_eps_c",
    "rest_floor_frac",
    "smooth_scale",
    "alr_spline_df",
    "alr_spline_spar",
    "alr_spline_cv",
    "nz_partner_min_frac",
    "max_retries",
    "chains",
    "iter_warmup",
    "iter_sampling",
    "adapt_delta",
    "max_treedepth",
    "metric",
    "init",
    "seed",
    "quiet",
    "progress",
    "progress_every",
    "silent_sampler",
    "n_workers_outer",
    "kfold_K",
    "kfold_R",
    "kfold_seed",
    "use_pathfinder_init",
    "pf_num_paths",
    "pf_draws",
    "pf_history_size",
    "pf_max_lbfgs_iters",
    "pf_psis_resample"
  )

  missing_controls <- setdiff(required_controls, names(controls))

  if (length(missing_controls)) {
    stop(
      sprintf(
        "Missing required control value(s): %s.",
        paste(missing_controls, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # phyloseq and abundance-table validation
  # ---------------------------------------------------------------------------

  if (!inherits(physeq, "phyloseq")) {
    stop(
      "`physeq` must inherit from class phyloseq.",
      call. = FALSE
    )
  }

  otu <- tryCatch(
    phyloseq::otu_table(physeq),
    error = function(e) NULL
  )
  md <- tryCatch(
    data.frame(
      phyloseq::sample_data(physeq),
      check.names = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(otu) ||
      is.null(md) ||
      phyloseq::ntaxa(physeq) < 1L ||
      phyloseq::nsamples(physeq) < 1L) {
    stop(
      "`physeq` must contain a non-empty OTU table and sample metadata.",
      call. = FALSE
    )
  }

  mat <- as(otu, "matrix")

  if (!phyloseq::taxa_are_rows(physeq)) {
    mat <- t(mat)
  }

  sample_names <- as.character(phyloseq::sample_names(physeq))
  taxa_names <- as.character(phyloseq::taxa_names(physeq))

  valid_sample_names <-
    length(sample_names) == ncol(mat) &&
    !anyNA(sample_names) &&
    all(nzchar(sample_names)) &&
    !anyDuplicated(sample_names)

  if (!valid_sample_names) {
    stop(
      paste(
        "`physeq` sample names must be present, unique,",
        "and match OTU columns."
      ),
      call. = FALSE
    )
  }

  valid_taxa_names <-
    length(taxa_names) == nrow(mat) &&
    !anyNA(taxa_names) &&
    all(nzchar(taxa_names)) &&
    !anyDuplicated(taxa_names)

  if (!valid_taxa_names) {
    stop(
      paste(
        "`physeq` taxa names must be present, unique,",
        "and match OTU rows."
      ),
      call. = FALSE
    )
  }

  if (!identical(colnames(mat), sample_names) ||
      !identical(rownames(mat), taxa_names)) {
    stop(
      "`physeq` OTU dimnames are inconsistent with phyloseq names.",
      call. = FALSE
    )
  }

  if (nrow(md) != length(sample_names) ||
      is.null(rownames(md)) ||
      anyNA(rownames(md)) ||
      any(!nzchar(rownames(md))) ||
      anyDuplicated(rownames(md)) ||
      !setequal(rownames(md), sample_names)) {
    stop(
      paste(
        "`physeq` OTU and sample-metadata dimensions",
        "or names are inconsistent."
      ),
      call. = FALSE
    )
  }

  # Canonical alignment before extracting metadata columns.
  md <- md[sample_names, , drop = FALSE]

  if (!is.numeric(mat) ||
      any(!is.finite(mat)) ||
      any(mat < 0)) {
    stop(
      paste(
        "`physeq` abundance values must be numeric,",
        "finite, and nonnegative."
      ),
      call. = FALSE
    )
  }

  storage.mode(mat) <- "double"

  sample_totals <- colSums(mat)

  if (any(!is.finite(sample_totals)) ||
      any(sample_totals <= 0)) {
    bad_samples <- sample_names[
      !is.finite(sample_totals) | sample_totals <= 0
    ]

    stop(
      sprintf(
        paste(
          "`physeq` contains sample(s) with a non-finite",
          "or zero total abundance: %s."
        ),
        paste(bad_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Canonical sample-wise closure. This supports both raw counts and
  # already-normalized inputs without allowing library size to enter smoothing.
  mat_rel <- sweep(mat, 2L, sample_totals, "/")

  if (any(!is.finite(mat_rel)) ||
      any(mat_rel < 0)) {
    stop(
      "Relative-abundance normalization produced invalid values.",
      call. = FALSE
    )
  }

  closure_error <- abs(colSums(mat_rel) - 1)
  closure_tolerance <- max(
    1e-12,
    64 * .Machine$double.eps * nrow(mat_rel)
  )

  if (any(!is.finite(closure_error)) ||
      max(closure_error) > closure_tolerance) {
    stop(
      "Relative-abundance normalization failed the closure invariant.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Metadata
  # ---------------------------------------------------------------------------

  subject_col <- chr1(subject_col, "subject_col")
  time_col <- chr1(time_col, "time_col")

  if (identical(subject_col, time_col)) {
    stop(
      "`subject_col` and `time_col` must identify different columns.",
      call. = FALSE
    )
  }

  missing_metadata <- setdiff(
    c(subject_col, time_col),
    names(md)
  )

  if (length(missing_metadata)) {
    stop(
      sprintf(
        "Missing sample metadata column(s): %s.",
        paste(missing_metadata, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  subject_raw <- md[[subject_col]]
  subject <- as.character(subject_raw)

  if (!length(subject) ||
      length(subject) != nrow(md) ||
      anyNA(subject) ||
      any(!nzchar(subject))) {
    stop(
      sprintf(
        "`%s` subject values must be present and non-missing.",
        subject_col
      ),
      call. = FALSE
    )
  }

  if (is.numeric(subject_raw) &&
      any(!is.finite(subject_raw))) {
    stop(
      sprintf(
        "`%s` numeric subject values must be finite.",
        subject_col
      ),
      call. = FALSE
    )
  }

  time_raw <- md[[time_col]]

  time <- if (is.numeric(time_raw)) {
    as.numeric(time_raw)
  } else {
    suppressWarnings(as.numeric(as.character(time_raw)))
  }

  if (length(time) != nrow(md) ||
      any(!is.finite(time))) {
    stop(
      sprintf(
        "`%s` time values must be numeric or safely coercible and finite.",
        time_col
      ),
      call. = FALSE
    )
  }

  time_by_subject <- split(time, subject)

  duplicate_time_subjects <- names(
    Filter(
      function(x) anyDuplicated(x) > 0L,
      time_by_subject
    )
  )

  if (length(duplicate_time_subjects)) {
    stop(
      sprintf(
        "Duplicate time values detected within subject(s): %s.",
        paste(duplicate_time_subjects, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Taxa selection
  # ---------------------------------------------------------------------------

  if (is.null(taxa_vec)) {
    taxa_vec <- taxa_names
  }

  if (!is.character(taxa_vec) ||
      !length(taxa_vec) ||
      anyNA(taxa_vec) ||
      any(!nzchar(taxa_vec)) ||
      anyDuplicated(taxa_vec)) {
    stop(
      "`taxa_vec` must be NULL or a non-empty unique character vector.",
      call. = FALSE
    )
  }

  unknown_taxa <- setdiff(taxa_vec, taxa_names)

  if (length(unknown_taxa)) {
    stop(
      sprintf(
        "Unknown taxa requested: %s.",
        paste(unknown_taxa, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(taxa_vec) < 2L) {
    stop(
      "`taxa_vec` must contain at least two taxa.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Scalar controls
  # ---------------------------------------------------------------------------

  logical_controls <- c(
    "alr_spline_cv",
    "quiet",
    "silent_sampler",
    "use_pathfinder_init",
    "pf_psis_resample"
  )

  for (name in logical_controls) {
    controls[[name]] <- log1(controls[[name]], name)
  }

  integer_minimums <- c(
    chains = 1,
    iter_warmup = 0,
    iter_sampling = 1,
    seed = 0,
    progress_every = 1,
    n_workers_outer = 1,
    kfold_K = 2,
    kfold_R = 1,
    kfold_seed = 0,
    max_retries = 0,
    pf_num_paths = 1,
    pf_draws = 1,
    pf_history_size = 1,
    pf_max_lbfgs_iters = 1,
    min_unique_times = 3,
    min_pairs = 1,
    max_treedepth = 1
  )

  for (name in names(integer_minimums)) {
    controls[[name]] <- int1(
      controls[[name]],
      name,
      lower = integer_minimums[[name]]
    )
  }

  controls$adapt_delta <- num1(
    controls$adapt_delta,
    "adapt_delta",
    lower = 0,
    upper = 1,
    open_lower = TRUE,
    open_upper = TRUE
  )

  for (name in c(
    "eps",
    "minpos_alpha",
    "eps_fixed",
    "lib_eps_c"
  )) {
    controls[[name]] <- num1(
      controls[[name]],
      name,
      lower = 0,
      open_lower = TRUE
    )
  }

  controls$rest_floor_frac <- num1(
    controls$rest_floor_frac,
    "rest_floor_frac",
    lower = 0,
    upper = 1,
    open_lower = TRUE
  )

  controls$nz_partner_min_frac <- num1(
    controls$nz_partner_min_frac,
    "nz_partner_min_frac",
    lower = 0,
    upper = 1
  )

  if (!is.null(controls$alr_spline_df) ||
      !is.null(controls$alr_spline_spar) ||
      !isTRUE(controls$alr_spline_cv)) {
    stop(
      paste(
        "Canonical CV smoothing requires",
        "`alr_spline_df = NULL`,",
        "`alr_spline_spar = NULL`, and",
        "`alr_spline_cv = TRUE`."
      ),
      call. = FALSE
    )
  }

  controls$zero_mode_alr <- choice1(
    controls$zero_mode_alr,
    "zero_mode_alr",
    c("minpos_time", "minpos_subject", "lib", "fixed")
  )

  controls$minpos_base <- choice1(
    controls$minpos_base,
    "minpos_base",
    c("ij", "triplet")
  )

  controls$smooth_scale <- choice1(
    controls$smooth_scale,
    "smooth_scale",
    c("logra", "alr")
  )

  controls$progress <- choice1(
    controls$progress,
    "progress",
    c("bar", "verbose", "none")
  )

  controls$metric <- choice1(
    controls$metric,
    "metric",
    c("diag_e", "dense_e", "unit_e")
  )

  subject_count <- length(unique(subject))

  if (controls$kfold_K > subject_count) {
    stop(
      "`kfold_K` cannot exceed the number of unique subjects.",
      call. = FALSE
    )
  }

  observations_per_subject <- lengths(time_by_subject)
  insufficient_time_subjects <- names(
    observations_per_subject[
      observations_per_subject < controls$min_unique_times
    ]
  )

  if (length(insufficient_time_subjects)) {
    stop(
      sprintf(
        paste(
          "Subject(s) have fewer than `min_unique_times = %d`",
          "distinct observations: %s."
        ),
        controls$min_unique_times,
        paste(insufficient_time_subjects, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Initialization contract
  # ---------------------------------------------------------------------------

  init <- controls$init

  valid_init_type <-
    is.null(init) ||
    is.function(init) ||
    (
      is.numeric(init) &&
      length(init) == 1L &&
      is.finite(init)
    ) ||
    inherits(init, "CmdStanPathfinder") ||
    is.environment(init) ||
    (
      is.list(init) &&
      length(init) > 0L
    )

  if (!valid_init_type) {
    stop(
      "`init` is malformed or unsupported.",
      call. = FALSE
    )
  }

  if (is.list(init) && !is.function(init)) {
    per_chain <- is.list(init[[1L]])

    if (per_chain) {
      if (length(init) != controls$chains) {
        stop(
          sprintf(
            paste(
              "`init` per-chain list must contain exactly",
              "%d entries."
            ),
            controls$chains
          ),
          call. = FALSE
        )
      }

      for (chain in seq_along(init)) {
        validate_named_init(
          init[[chain]],
          sprintf("init[[%d]]", chain)
        )
      }
    } else {
      validate_named_init(init, "init")
    }
  }

  # ---------------------------------------------------------------------------
  # Derived seed-range contract
  # ---------------------------------------------------------------------------

  n_taxa <- length(taxa_vec)

  # For i < j, the largest pair-direction offset occurs at
  # i = n_taxa - 1 and j = n_taxa.
  max_direction_offset <-
    100000 * (n_taxa - 1) +
    1000 * n_taxa +
    2

  # Main-fit retries may advance the direction seed. K-fold then derives
  # another offset and its own retries may advance it again.
  max_derived_sampler_seed <-
    as.double(controls$seed) +
    max_direction_offset +
    controls$max_retries +
    1000 * controls$kfold_R +
    controls$kfold_K +
    controls$max_retries

  if (!is.finite(max_derived_sampler_seed) ||
      max_derived_sampler_seed > .Machine$integer.max) {
    stop(
      paste(
        "`seed` is too large for the selected taxa count,",
        "K-fold settings, and retry policy;",
        "a derived CmdStan seed would exceed the R integer range."
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Canonical aligned outputs
  # ---------------------------------------------------------------------------

  meta_df <- data.frame(
    Sample = sample_names,
    subject = subject,
    time = time,
    stringsAsFactors = FALSE
  )

  meta_df <- meta_df[
    order(meta_df$subject, meta_df$time, method = "radix"),
    ,
    drop = FALSE
  ]
  rownames(meta_df) <- NULL

  list(
    meta_df = meta_df,
    mat_rel = mat_rel,
    taxa_vec = taxa_vec,
    subject_col = subject_col,
    time_col = time_col,
    controls = controls
  )
}
