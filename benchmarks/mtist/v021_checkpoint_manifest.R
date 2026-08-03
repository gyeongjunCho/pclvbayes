# Benchmark-only checkpoint and restart infrastructure for ROADMAP V021-03.

v021_checkpoint_states <- c(
  "pending", "running", "completed", "failed", "skipped", "incomplete"
)

v021_manifest_fields <- c(
  "task_id", "pair_id", "direction_id", "direction_index", "target", "source",
  "seed", "chain_seeds", "execution_state", "elapsed_time", "retry_count",
  "pathfinder_used", "output_location", "completion_timestamp", "failure_reason",
  "retry_history", "interrupted", "interruption_reason"
)

v021_checkpoint_fields <- c(
  "checkpoint_version", v021_manifest_fields, "result"
)

.v021_checkpoint_version <- "v021_checkpoint_v1"
.v021_manifest_version <- "v021_manifest_v1"

.v021_checkpoint_copy <- function(x) unserialize(serialize(x, NULL, version = 3))

.v021_one_string <- function(x, name, missing_ok = FALSE) {
  if (!is.character(x) || length(x) != 1L || (!missing_ok && is.na(x)) ||
      (!missing_ok && !nzchar(x))) stop(name, " must be one non-empty string.")
  invisible(TRUE)
}

.v021_one_number <- function(x, name, integer = FALSE, nonnegative = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      (integer && x != as.integer(x)) || (nonnegative && x < 0))
    stop(name, " is invalid.")
  invisible(TRUE)
}

.v021_validate_retry_history <- function(x) {
  if (!is.list(x)) stop("retry_history must be a list.")
  if (!length(x)) return(invisible(TRUE))
  valid <- vapply(x, function(attempt) {
    is.list(attempt) && !is.null(names(attempt)) &&
      all(c("attempt", "seed", "status") %in% names(attempt)) &&
      is.numeric(attempt$attempt) && length(attempt$attempt) == 1L &&
      is.numeric(attempt$seed) && length(attempt$seed) == 1L &&
      is.character(attempt$status) && length(attempt$status) == 1L
  }, logical(1))
  if (!all(valid)) stop("retry_history contains an invalid attempt record.")
  invisible(TRUE)
}

v021_chain_seeds <- function(direction_seed, chains) {
  .v021_one_number(direction_seed, "direction_seed", integer = TRUE, nonnegative = TRUE)
  .v021_one_number(chains, "chains", integer = TRUE, nonnegative = TRUE)
  if (chains < 1L) stop("chains must be positive.")
  seeds <- as.double(direction_seed) + 10000000 * seq_len(as.integer(chains))
  if (any(seeds > .Machine$integer.max)) stop("Deterministic chain seed exceeds integer range.")
  as.integer(seeds)
}

.v021_manifest_pair_ids <- function(task_table) {
  unordered <- vapply(seq_len(nrow(task_table)), function(i)
    paste(sort(c(as.character(task_table$target[[i]]),
                 as.character(task_table$source[[i]]))), collapse = "|"),
    character(1))
  if ("pair_id" %in% names(task_table)) {
    pair_id <- as.character(task_table$pair_id)
    if (anyNA(pair_id) || any(!nzchar(pair_id))) stop("pair_id is incomplete.")
  } else {
    pair_id <- sprintf("pair-%06d", as.integer(task_table$task_id))
  }
  by_pair <- split(unordered, pair_id)
  if (any(vapply(by_pair, function(x) length(unique(x)) != 1L, logical(1))))
    stop("A pair_id identifies conflicting unordered pairs.")
  by_unordered_pair <- split(pair_id, unordered)
  if (any(vapply(by_unordered_pair, function(x) length(unique(x)) != 1L, logical(1))))
    stop("An unordered pair has conflicting pair identifiers.")
  by_task <- split(pair_id, task_table$task_id)
  if (any(vapply(by_task, function(x) length(unique(x)) != 1L, logical(1))))
    stop("A task_id identifies conflicting pair identifiers.")
  pair_id
}

build_v021_checkpoint_manifest <- function(task_table, chains, checkpoint_dir) {
  required <- c("task_id", "direction_index", "target", "source", "seed")
  allowed <- c(required, "pair_id")
  if (!is.data.frame(task_table) || !all(required %in% names(task_table)) ||
      length(setdiff(names(task_table), allowed)))
    stop("task_table does not match the checkpoint manifest input schema.")
  if (!nrow(task_table)) stop("task_table must not be empty.")
  .v021_one_number(chains, "chains", integer = TRUE, nonnegative = TRUE)
  if (chains < 1L) stop("chains must be positive.")
  .v021_one_string(checkpoint_dir, "checkpoint_dir")
  for (nm in c("task_id", "direction_index", "seed")) {
    if (!is.numeric(task_table[[nm]]) || anyNA(task_table[[nm]]) ||
        any(!is.finite(task_table[[nm]])) || any(task_table[[nm]] != as.integer(task_table[[nm]])))
      stop(nm, " must contain finite integers.")
  }
  if (anyDuplicated(task_table$direction_index)) stop("direction_index must be unique.")
  if (anyNA(task_table[c("target", "source")]) ||
      any(!nzchar(as.character(task_table$target))) ||
      any(!nzchar(as.character(task_table$source))) ||
      any(as.character(task_table$target) == as.character(task_table$source)))
    stop("Invalid target/source identity.")

  pair_id <- .v021_manifest_pair_ids(task_table)
  pair_rows <- split(seq_len(nrow(task_table)), pair_id)
  complete <- vapply(pair_rows, function(ix) {
    length(ix) == 2L && length(unique(task_table$task_id[ix])) == 1L &&
      identical(as.character(task_table$target[ix]),
                rev(as.character(task_table$source[ix])))
  }, logical(1))
  if (any(!complete)) stop("Every checkpoint pair must contain both directions.")

  ord <- order(as.integer(task_table$task_id), as.integer(task_table$direction_index))
  task_table <- task_table[ord, , drop = FALSE]
  pair_id <- pair_id[ord]
  direction_index <- as.integer(task_table$direction_index)
  direction_id <- sprintf("direction-%06d", direction_index)
  chain_seeds <- lapply(as.integer(task_table$seed), v021_chain_seeds, chains = chains)
  all_chain_seeds <- unlist(chain_seeds, use.names = FALSE)
  if (anyDuplicated(all_chain_seeds)) stop("Deterministic chain seeds are not unique.")

  manifest <- data.frame(
    task_id = as.integer(task_table$task_id), pair_id = pair_id,
    direction_id = direction_id, direction_index = direction_index,
    target = as.character(task_table$target), source = as.character(task_table$source),
    seed = as.integer(task_table$seed), execution_state = "pending",
    elapsed_time = 0, retry_count = 0L, pathfinder_used = FALSE,
    output_location = file.path(checkpoint_dir, paste0(direction_id, ".rds")),
    completion_timestamp = NA_character_, failure_reason = NA_character_,
    interrupted = FALSE, interruption_reason = NA_character_,
    stringsAsFactors = FALSE
  )
  manifest$chain_seeds <- I(chain_seeds)
  manifest$retry_history <- I(rep(list(list()), nrow(manifest)))
  manifest <- manifest[v021_manifest_fields]
  attr(manifest, "manifest_version") <- .v021_manifest_version
  validate_v021_checkpoint_manifest(manifest)
  manifest
}

validate_v021_checkpoint_manifest <- function(manifest) {
  if (!is.data.frame(manifest) || !identical(names(manifest), v021_manifest_fields))
    stop("Invalid checkpoint manifest schema.")
  if (!identical(attr(manifest, "manifest_version"), .v021_manifest_version))
    stop("Invalid checkpoint manifest version.")
  if (!nrow(manifest) || anyDuplicated(manifest$direction_index) ||
      anyDuplicated(manifest$direction_id) || anyDuplicated(manifest$output_location))
    stop("Checkpoint manifest identities are empty or duplicated.")
  expected_order <- order(manifest$task_id, manifest$direction_index)
  if (!identical(expected_order, seq_len(nrow(manifest))))
    stop("Checkpoint manifest is not in deterministic task order.")
  if (any(!manifest$execution_state %in% v021_checkpoint_states))
    stop("Checkpoint manifest contains an invalid execution state.")
  if (anyNA(manifest[c("task_id", "pair_id", "direction_id", "direction_index",
                       "target", "source", "seed", "execution_state", "elapsed_time",
                       "retry_count", "pathfinder_used", "output_location", "interrupted")]))
    stop("Checkpoint manifest contains missing required metadata.")
  if (any(manifest$elapsed_time < 0) || any(manifest$retry_count < 0))
    stop("Checkpoint manifest contains negative execution metadata.")
  if (any(manifest$execution_state == "completed" & is.na(manifest$completion_timestamp)))
    stop("Completed manifest rows require completion_timestamp.")
  if (any(manifest$execution_state != "completed" & !is.na(manifest$completion_timestamp)))
    stop("Only completed manifest rows may have completion_timestamp.")
  if (any(manifest$execution_state == "failed" & is.na(manifest$failure_reason)))
    stop("Failed manifest rows require failure_reason.")
  invisible(lapply(manifest$retry_history, .v021_validate_retry_history))
  chain_values <- unlist(manifest$chain_seeds, use.names = FALSE)
  if (!length(chain_values) || anyNA(chain_values) || anyDuplicated(chain_values))
    stop("Checkpoint manifest chain seeds are incomplete or duplicated.")
  invisible(TRUE)
}

.v021_record_from_row <- function(row, execution_state, elapsed_time = 0,
                                  retry_history = list(), pathfinder_used = FALSE,
                                  completion_timestamp = NA_character_,
                                  failure_reason = NA_character_, interrupted = FALSE,
                                  interruption_reason = NA_character_, result = NULL) {
  if (!is.data.frame(row) || nrow(row) != 1L) stop("row must be one manifest row.")
  if (!execution_state %in% v021_checkpoint_states) stop("Invalid execution_state.")
  record <- list(
    checkpoint_version = .v021_checkpoint_version,
    task_id = row$task_id[[1L]], pair_id = row$pair_id[[1L]],
    direction_id = row$direction_id[[1L]], direction_index = row$direction_index[[1L]],
    target = row$target[[1L]], source = row$source[[1L]], seed = row$seed[[1L]],
    chain_seeds = .v021_checkpoint_copy(row$chain_seeds[[1L]]),
    execution_state = execution_state, elapsed_time = as.numeric(elapsed_time),
    retry_count = as.integer(max(length(retry_history) - 1L, 0L)),
    pathfinder_used = isTRUE(pathfinder_used), output_location = row$output_location[[1L]],
    completion_timestamp = completion_timestamp, failure_reason = failure_reason,
    retry_history = .v021_checkpoint_copy(retry_history), interrupted = isTRUE(interrupted),
    interruption_reason = interruption_reason, result = .v021_checkpoint_copy(result)
  )
  validate_v021_checkpoint(record)
  record
}

validate_v021_checkpoint <- function(checkpoint) {
  if (!is.list(checkpoint) || !identical(names(checkpoint), v021_checkpoint_fields))
    stop("Invalid checkpoint schema.")
  if (!identical(checkpoint$checkpoint_version, .v021_checkpoint_version))
    stop("Invalid checkpoint version.")
  .v021_one_number(checkpoint$task_id, "task_id", integer = TRUE, nonnegative = TRUE)
  .v021_one_number(checkpoint$direction_index, "direction_index", integer = TRUE,
                   nonnegative = TRUE)
  .v021_one_number(checkpoint$seed, "seed", integer = TRUE, nonnegative = TRUE)
  for (nm in c("pair_id", "direction_id", "target", "source", "output_location"))
    .v021_one_string(checkpoint[[nm]], nm)
  if (!checkpoint$execution_state %in% v021_checkpoint_states)
    stop("Invalid checkpoint execution_state.")
  .v021_one_number(checkpoint$elapsed_time, "elapsed_time", nonnegative = TRUE)
  .v021_one_number(checkpoint$retry_count, "retry_count", integer = TRUE,
                   nonnegative = TRUE)
  if (!is.logical(checkpoint$pathfinder_used) || length(checkpoint$pathfinder_used) != 1L ||
      is.na(checkpoint$pathfinder_used)) stop("Invalid pathfinder_used.")
  if (!is.logical(checkpoint$interrupted) || length(checkpoint$interrupted) != 1L ||
      is.na(checkpoint$interrupted)) stop("Invalid interrupted flag.")
  if (!is.integer(checkpoint$chain_seeds) || !length(checkpoint$chain_seeds) ||
      anyNA(checkpoint$chain_seeds) || anyDuplicated(checkpoint$chain_seeds))
    stop("Invalid checkpoint chain_seeds.")
  .v021_validate_retry_history(checkpoint$retry_history)
  if (!identical(checkpoint$retry_count,
                 as.integer(max(length(checkpoint$retry_history) - 1L, 0L))))
    stop("retry_count does not match retry_history.")
  completed <- identical(checkpoint$execution_state, "completed")
  if (completed && (is.na(checkpoint$completion_timestamp) || is.null(checkpoint$result)))
    stop("Completed checkpoints require a timestamp and retained result.")
  if (!completed && !is.na(checkpoint$completion_timestamp))
    stop("Only completed checkpoints may have completion_timestamp.")
  if (identical(checkpoint$execution_state, "failed") && is.na(checkpoint$failure_reason))
    stop("Failed checkpoints require failure_reason.")
  if (isTRUE(checkpoint$interrupted) && is.na(checkpoint$interruption_reason))
    stop("Interrupted checkpoints require interruption_reason.")
  invisible(TRUE)
}

.v021_atomic_save_rds <- function(object, path) {
  .v021_one_string(path, "path")
  directory <- dirname(path)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE))
    stop("Could not create checkpoint directory: ", directory)
  temporary <- tempfile(paste0(".", basename(path), ".tmp-"), tmpdir = directory)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  connection <- file(temporary, open = "wb")
  tryCatch({
    serialize(object, connection, version = 3)
    flush(connection)
  }, finally = close(connection))
  if (!file.rename(temporary, path)) stop("Atomic checkpoint rename failed for: ", path)
  invisible(path)
}

write_v021_checkpoint_atomic <- function(checkpoint, path = checkpoint$output_location) {
  validate_v021_checkpoint(checkpoint)
  if (!identical(path, checkpoint$output_location))
    stop("Checkpoint path does not match output_location.")
  if (file.exists(path)) {
    existing <- readRDS(path)
    validate_v021_checkpoint(existing)
    if (identical(existing$execution_state, "completed"))
      stop("Refusing to overwrite a completed checkpoint.")
    identity_fields <- c("task_id", "pair_id", "direction_id", "direction_index",
                         "target", "source", "seed", "chain_seeds", "output_location")
    if (!identical(existing[identity_fields], checkpoint[identity_fields]))
      stop("Refusing to replace a checkpoint with conflicting identity.")
  }
  .v021_atomic_save_rds(checkpoint, path)
}

write_v021_manifest_atomic <- function(manifest, path) {
  validate_v021_checkpoint_manifest(manifest)
  .v021_atomic_save_rds(manifest, path)
}

read_v021_checkpoint_manifest <- function(path) {
  manifest <- readRDS(path)
  validate_v021_checkpoint_manifest(manifest)
  manifest
}

.v021_identity_matches <- function(row, checkpoint) {
  fields <- c("task_id", "pair_id", "direction_id", "direction_index",
              "target", "source", "seed", "chain_seeds", "output_location")
  manifest_values <- lapply(fields, function(nm) {
    if (nm == "chain_seeds") row[[nm]][[1L]] else row[[nm]][[1L]]
  })
  names(manifest_values) <- fields
  identical(manifest_values, checkpoint[fields])
}

.v021_apply_checkpoint <- function(manifest, index, checkpoint) {
  fields <- c("execution_state", "elapsed_time", "retry_count", "pathfinder_used",
              "completion_timestamp", "failure_reason", "retry_history",
              "interrupted", "interruption_reason")
  for (nm in fields) {
    if (nm == "retry_history") manifest[[nm]][[index]] <- checkpoint[[nm]]
    else manifest[[nm]][[index]] <- checkpoint[[nm]]
  }
  manifest
}

plan_v021_checkpoint_restart <- function(manifest) {
  validate_v021_checkpoint_manifest(manifest)
  reconciled <- .v021_checkpoint_copy(manifest)
  for (i in seq_len(nrow(reconciled))) {
    path <- reconciled$output_location[[i]]
    if (!file.exists(path)) {
      if (identical(reconciled$execution_state[[i]], "completed"))
        stop("Completed manifest row has no checkpoint: ", path)
      if (identical(reconciled$execution_state[[i]], "running")) {
        reconciled$execution_state[[i]] <- "incomplete"
        reconciled$interrupted[[i]] <- TRUE
        reconciled$interruption_reason[[i]] <- "checkpoint_missing_after_running_state"
      }
      next
    }
    checkpoint <- readRDS(path)
    validate_v021_checkpoint(checkpoint)
    if (!.v021_identity_matches(reconciled[i, , drop = FALSE], checkpoint))
      stop("Checkpoint identity conflicts with manifest: ", path)
    if (identical(reconciled$execution_state[[i]], "completed") &&
        !identical(checkpoint$execution_state, "completed"))
      stop("Completed manifest row has a non-completed checkpoint: ", path)
    if (identical(checkpoint$execution_state, "running")) {
      checkpoint$execution_state <- "incomplete"
      checkpoint$interrupted <- TRUE
      checkpoint$interruption_reason <- "interrupted_while_running"
      write_v021_checkpoint_atomic(checkpoint)
    }
    reconciled <- .v021_apply_checkpoint(reconciled, i, checkpoint)
  }
  validate_v021_checkpoint_manifest(reconciled)
  resume_states <- c("pending", "incomplete")
  list(
    manifest = reconciled,
    resume_indices = which(reconciled$execution_state %in% resume_states),
    completed_indices = which(reconciled$execution_state == "completed")
  )
}

.v021_outcome_value <- function(outcome, name, default) {
  if (!name %in% names(outcome) || is.null(outcome[[name]])) default else outcome[[name]]
}

execute_v021_checkpoint_restart <- function(manifest, run_task, manifest_path = NULL,
                                            timestamp = function()
                                              format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  if (!is.function(run_task)) stop("run_task must be a function.")
  plan <- plan_v021_checkpoint_restart(manifest)
  current <- plan$manifest
  if (!is.null(manifest_path)) write_v021_manifest_atomic(current, manifest_path)
  for (i in plan$resume_indices) {
    row <- current[i, , drop = FALSE]
    running <- .v021_record_from_row(
      row, "running", elapsed_time = row$elapsed_time[[1L]],
      retry_history = row$retry_history[[1L]],
      pathfinder_used = row$pathfinder_used[[1L]],
      interrupted = row$interrupted[[1L]],
      interruption_reason = row$interruption_reason[[1L]]
    )
    write_v021_checkpoint_atomic(running)
    current <- .v021_apply_checkpoint(current, i, running)
    if (!is.null(manifest_path)) write_v021_manifest_atomic(current, manifest_path)

    started <- proc.time()[["elapsed"]]
    outcome <- tryCatch(
      run_task(.v021_checkpoint_copy(row)),
      interrupt = function(e) list(
        execution_state = "incomplete", interrupted = TRUE,
        interruption_reason = conditionMessage(e), failure_reason = conditionMessage(e)
      ),
      error = function(e) list(
        execution_state = "failed", failure_reason = conditionMessage(e)
      )
    )
    elapsed <- proc.time()[["elapsed"]] - started
    if (!is.list(outcome) || !outcome$execution_state %in%
        c("completed", "failed", "skipped", "incomplete"))
      stop("run_task returned an invalid checkpoint outcome.")
    state <- outcome$execution_state
    retries <- c(row$retry_history[[1L]],
                 .v021_outcome_value(outcome, "retry_history", list()))
    completion <- if (identical(state, "completed")) as.character(timestamp()) else NA_character_
    checkpoint <- .v021_record_from_row(
      row, state,
      elapsed_time = row$elapsed_time[[1L]] +
        .v021_outcome_value(outcome, "elapsed_time", elapsed),
      retry_history = retries,
      pathfinder_used = row$pathfinder_used[[1L]] ||
        isTRUE(.v021_outcome_value(outcome, "pathfinder_used", FALSE)),
      completion_timestamp = completion,
      failure_reason = .v021_outcome_value(outcome, "failure_reason", NA_character_),
      interrupted = row$interrupted[[1L]] ||
        isTRUE(.v021_outcome_value(outcome, "interrupted", FALSE)),
      interruption_reason = .v021_outcome_value(
        outcome, "interruption_reason", row$interruption_reason[[1L]]),
      result = .v021_outcome_value(outcome, "result", NULL)
    )
    write_v021_checkpoint_atomic(checkpoint)
    current <- .v021_apply_checkpoint(current, i, checkpoint)
    if (!is.null(manifest_path)) write_v021_manifest_atomic(current, manifest_path)
  }
  validate_v021_checkpoint_manifest(current)
  current
}

# V021-03 canonical execution contract. The v1 helpers above remain only for
# the already archived preflight format; full V021 execution uses these v2
# objects so immutable identity is never mixed with mutable attempt state.

v021_execution_manifest_schema <- "v021_execution_manifest_v2"
v021_task_status_schema <- "v021_task_status_v1"
v021_attempt_ledger_schema <- "v021_attempt_ledger_v1"
v021_completion_artifact_schema <- "v021_completion_artifact_v1"
v021_run_provenance_schema <- "v021_run_provenance_v1"

v021_execution_states <- c("pending", "running", "completed", "failed")
v021_legal_state_transitions <- list(
  pending = "running", running = c("completed", "failed"),
  completed = character(), failed = "running"
)

.v021_canonicalize_hash_value <- function(x) {
  if (is.data.frame(x)) {
    out <- x
    names(out) <- as.character(names(out))
    return(out)
  }
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms)) {
      if (anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms))
        stop("Canonical hash input has invalid list names.")
      x <- x[order(nms, method = "radix")]
    }
    return(lapply(x, .v021_canonicalize_hash_value))
  }
  if (is.factor(x)) return(as.character(x))
  x
}

v021_sha256 <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The digest package is required for benchmark provenance hashing.")
  digest::digest(serialize(.v021_canonicalize_hash_value(x), NULL, version = 3),
                 algo = "sha256", serialize = FALSE)
}

.v021_execution_configuration <- function(preprocessing_config,
                                           posterior_config,
                                           kfold_config,
                                           predictive_config,
                                           public_seed, kfold_seed) {
  config <- list(
    preprocessing = .v021_checkpoint_copy(preprocessing_config),
    posterior = .v021_checkpoint_copy(posterior_config),
    kfold = .v021_checkpoint_copy(kfold_config),
    predictive = .v021_checkpoint_copy(predictive_config),
    public_seed = as.integer(public_seed),
    kfold_seed = as.integer(kfold_seed)
  )
  assert_v021_truth_free_schema(config, "execution configuration")
  .v021_assert_plain_value(config, "execution configuration")
  list(value = config, hash = v021_sha256(config))
}

build_v021_execution_manifest <- function(
    dataset_id, taxa_order, task_table, public_seed, kfold_seed,
    preprocessing_config, posterior_config, kfold_config, predictive_config,
    provenance) {
  if (!is.character(dataset_id) || length(dataset_id) != 1L || is.na(dataset_id) ||
      !nzchar(dataset_id)) stop("dataset_id must be one non-empty string.")
  taxa_order <- as.character(taxa_order)
  if (length(taxa_order) < 2L || anyNA(taxa_order) || any(!nzchar(taxa_order)) ||
      anyDuplicated(taxa_order)) stop("taxa_order must contain unique taxa.")
  required <- c("task_id", "direction_index", "target", "source", "seed")
  if (!is.data.frame(task_table) || !all(required %in% names(task_table)))
    stop("task_table lacks canonical direction identity.")
  assert_v021_truth_free_schema(task_table, "task manifest")
  if (length(setdiff(names(task_table), required)))
    stop("task_table contains fields outside the canonical manifest input schema.")
  task_table <- task_table[required]
  task_table$target <- as.character(task_table$target)
  task_table$source <- as.character(task_table$source)
  task_table <- task_table[order(as.integer(task_table$direction_index)), , drop = FALSE]
  rownames(task_table) <- NULL
  if (anyNA(task_table) || anyDuplicated(task_table$direction_index) ||
      any(task_table$target == task_table$source) ||
      any(!task_table$target %in% taxa_order) || any(!task_table$source %in% taxa_order))
    stop("task_table has invalid, duplicate, self, or unknown directions.")
  expected_count <- length(taxa_order) * (length(taxa_order) - 1L)
  if (nrow(task_table) != expected_count ||
      !identical(as.integer(task_table$direction_index), seq_len(expected_count)))
    stop("task_table does not contain every canonical directed non-self task.")
  source_index <- match(task_table$source, taxa_order)
  target_index <- match(task_table$target, taxa_order)
  pair_key <- vapply(seq_len(nrow(task_table)), function(i) {
    indices <- sort(c(source_index[[i]], target_index[[i]]))
    paste(taxa_order[indices], collapse = "|")
  }, character(1))
  pair_levels <- unique(pair_key)
  pair_ordinal <- match(pair_key, pair_levels)
  if (length(pair_levels) != choose(length(taxa_order), 2L) ||
      any(table(pair_key) != 2L))
    stop("Canonical unordered-pair coverage is incomplete.")
  tasks <- data.frame(
    dataset_id = rep(dataset_id, nrow(task_table)),
    pair_id = sprintf("pair-%06d", pair_ordinal),
    unordered_pair = pair_key,
    source = task_table$source, target = task_table$target,
    source_index = as.integer(source_index), target_index = as.integer(target_index),
    directed_task_id = sprintf("direction-%06d", task_table$direction_index),
    task_ordinal = as.integer(task_table$direction_index),
    public_seed = rep(as.integer(public_seed), nrow(task_table)),
    direction_seed = as.integer(task_table$seed),
    kfold_seed = rep(as.integer(kfold_seed), nrow(task_table)),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(tasks$directed_task_id) ||
      anyDuplicated(paste(tasks$source, tasks$target, sep = "\r")))
    stop("Directed task identity is duplicated.")
  configuration <- .v021_execution_configuration(
    preprocessing_config, posterior_config, kfold_config, predictive_config,
    public_seed, kfold_seed)
  if (!is.list(provenance) || is.null(names(provenance)))
    stop("provenance must be a named truth-free list.")
  assert_v021_truth_free_schema(provenance, "manifest provenance")
  .v021_assert_plain_value(provenance, "manifest provenance")
  if (any(grepl("timestamp|wall_clock|pid|result_root", names(provenance),
                ignore.case = TRUE)))
    stop("Runtime or result-root metadata is forbidden in canonical provenance.")
  provenance <- .v021_canonicalize_hash_value(provenance)
  manifest_without_hash <- list(
    manifest_schema = v021_execution_manifest_schema,
    dataset_id = dataset_id,
    taxa_order = taxa_order,
    unordered_pair_count = as.integer(length(pair_levels)),
    directed_task_count = as.integer(nrow(tasks)),
    configuration = configuration$value,
    configuration_hash = configuration$hash,
    provenance = .v021_checkpoint_copy(provenance),
    tasks = tasks
  )
  manifest <- c(manifest_without_hash,
                list(manifest_hash = v021_sha256(manifest_without_hash)))
  validate_v021_execution_manifest(manifest)
  manifest
}

validate_v021_execution_manifest <- function(manifest) {
  fields <- c("manifest_schema", "dataset_id", "taxa_order",
              "unordered_pair_count", "directed_task_count", "configuration",
              "configuration_hash", "provenance", "tasks", "manifest_hash")
  if (!is.list(manifest) || !identical(names(manifest), fields) ||
      !identical(manifest$manifest_schema, v021_execution_manifest_schema))
    stop("Invalid V021 execution manifest schema.")
  assert_v021_truth_free_schema(manifest, "execution manifest")
  taxa <- manifest$taxa_order
  tasks <- manifest$tasks
  required_tasks <- c("dataset_id", "pair_id", "unordered_pair", "source", "target",
                      "source_index", "target_index", "directed_task_id",
                      "task_ordinal", "public_seed", "direction_seed", "kfold_seed")
  if (!is.character(taxa) || anyNA(taxa) || anyDuplicated(taxa) ||
      !is.data.frame(tasks) || !identical(names(tasks), required_tasks))
    stop("Execution manifest taxa or task schema is malformed.")
  expected <- length(taxa) * (length(taxa) - 1L)
  if (!identical(nrow(tasks), expected) ||
      !identical(tasks$task_ordinal, seq_len(expected)) ||
      any(tasks$source == tasks$target) || anyDuplicated(tasks$directed_task_id) ||
      anyDuplicated(paste(tasks$source, tasks$target, sep = "\r")) ||
      !identical(tasks$source_index, as.integer(match(tasks$source, taxa))) ||
      !identical(tasks$target_index, as.integer(match(tasks$target, taxa))) ||
      !all(tasks$dataset_id == manifest$dataset_id) ||
      !identical(manifest$unordered_pair_count, as.integer(choose(length(taxa), 2L))) ||
      !identical(manifest$directed_task_count, as.integer(expected)))
    stop("Execution manifest canonical task identities or ordering are invalid.")
  if (!identical(manifest$configuration_hash,
                 v021_sha256(manifest$configuration)))
    stop("Execution manifest configuration hash mismatch.")
  expected_hash <- v021_sha256(manifest[setdiff(fields, "manifest_hash")])
  if (!identical(manifest$manifest_hash, expected_hash))
    stop("Execution manifest hash mismatch.")
  invisible(TRUE)
}

compare_v021_execution_manifests <- function(stored, requested) {
  validate_v021_execution_manifest(stored)
  validate_v021_execution_manifest(requested)
  fields <- c("manifest_schema", "dataset_id", "taxa_order", "configuration_hash",
              "provenance", "tasks", "manifest_hash")
  mismatch <- fields[!vapply(fields, function(nm)
    identical(stored[[nm]], requested[[nm]]), logical(1))]
  if (length(mismatch))
    stop("Stored immutable manifest mismatch: ", paste(mismatch, collapse = ", "),
         ". Use a new result root.")
  invisible(TRUE)
}

v021_execution_paths <- function(result_root) {
  root <- normalizePath(result_root, mustWork = FALSE)
  list(
    result_root = root,
    manifest = file.path(root, "canonical_manifest.rds"),
    status = file.path(root, "task_status.rds"),
    attempts = file.path(root, "attempt_ledger.rds"),
    artifact_root = file.path(root, "artifacts")
  )
}

.v021_atomic_write_validated_rds <- function(object, path, validator) {
  .v021_one_string(path, "path")
  if (!is.function(validator)) stop("validator must be a function.")
  validator(object)
  directory <- dirname(path)
  if (!dir.exists(directory) &&
      !dir.create(directory, recursive = TRUE, showWarnings = FALSE))
    stop("Could not create state directory: ", directory)
  temporary <- tempfile(paste0(".", basename(path), ".tmp-"), tmpdir = directory)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  connection <- file(temporary, open = "wb")
  tryCatch({ serialize(object, connection, version = 3); flush(connection) },
           finally = close(connection))
  candidate <- tryCatch(readRDS(temporary), error = identity)
  if (inherits(candidate, "error")) stop("Atomic temporary state is unreadable.")
  validator(candidate)
  if (!file.rename(temporary, path)) stop("Atomic state rename failed for: ", path)
  invisible(path)
}

write_v021_execution_manifest_once <- function(manifest, path) {
  validate_v021_execution_manifest(manifest)
  if (file.exists(path)) {
    stored <- read_v021_execution_manifest(path)
    compare_v021_execution_manifests(stored, manifest)
    return(invisible(path))
  }
  .v021_atomic_write_validated_rds(manifest, path,
                                   validate_v021_execution_manifest)
}

read_v021_execution_manifest <- function(path) {
  if (!file.exists(path)) stop("Canonical manifest is missing: ", path)
  value <- tryCatch(readRDS(path), error = identity)
  if (inherits(value, "error")) stop("Canonical manifest is malformed: ", path)
  validate_v021_execution_manifest(value)
  value
}

new_v021_task_status <- function(manifest) {
  validate_v021_execution_manifest(manifest)
  table <- data.frame(
    directed_task_id = manifest$tasks$directed_task_id,
    task_ordinal = manifest$tasks$task_ordinal,
    state = "pending", attempt_count = 0L,
    last_attempt_id = NA_character_, updated_timestamp = NA_character_,
    terminal_reason = NA_character_, artifact_path = NA_character_,
    artifact_validation = "not_checked", stringsAsFactors = FALSE)
  structure(list(status_schema = v021_task_status_schema,
                 manifest_hash = manifest$manifest_hash, tasks = table),
            class = "v021_task_status")
}

validate_v021_task_status <- function(status, manifest = NULL) {
  fields <- c("status_schema", "manifest_hash", "tasks")
  task_fields <- c("directed_task_id", "task_ordinal", "state", "attempt_count",
                   "last_attempt_id", "updated_timestamp", "terminal_reason",
                   "artifact_path", "artifact_validation")
  if (!is.list(status) || !identical(names(status), fields) ||
      !identical(status$status_schema, v021_task_status_schema) ||
      !is.data.frame(status$tasks) || !identical(names(status$tasks), task_fields) ||
      any(!status$tasks$state %in% v021_execution_states) ||
      any(status$tasks$attempt_count < 0L) || anyDuplicated(status$tasks$directed_task_id))
    stop("Malformed V021 task status.")
  assert_v021_truth_free_schema(status, "task status")
  if (!is.null(manifest)) {
    validate_v021_execution_manifest(manifest)
    if (!identical(status$manifest_hash, manifest$manifest_hash) ||
        !identical(status$tasks$directed_task_id, manifest$tasks$directed_task_id) ||
        !identical(status$tasks$task_ordinal, manifest$tasks$task_ordinal))
      stop("Task status provenance conflicts with the immutable manifest.")
  }
  invisible(TRUE)
}

new_v021_attempt_ledger <- function(manifest) {
  validate_v021_execution_manifest(manifest)
  rows <- data.frame(
    attempt_id = character(), directed_task_id = character(),
    task_ordinal = integer(), attempt_number = integer(),
    starting_state = character(), ending_state = character(),
    start_timestamp = character(), end_timestamp = character(),
    direction_seed = integer(), worker_provenance = character(),
    terminal_reason = character(), failure_class = character(),
    artifact_path = character(), artifact_validation = character(),
    stringsAsFactors = FALSE)
  structure(list(ledger_schema = v021_attempt_ledger_schema,
                 manifest_hash = manifest$manifest_hash, attempts = rows),
            class = "v021_attempt_ledger")
}

validate_v021_attempt_ledger <- function(ledger, manifest = NULL) {
  attempt_fields <- c(
    "attempt_id", "directed_task_id", "task_ordinal", "attempt_number",
    "starting_state", "ending_state", "start_timestamp", "end_timestamp",
    "direction_seed", "worker_provenance", "terminal_reason", "failure_class",
    "artifact_path", "artifact_validation")
  if (!is.list(ledger) ||
      !identical(names(ledger), c("ledger_schema", "manifest_hash", "attempts")) ||
      !identical(ledger$ledger_schema, v021_attempt_ledger_schema) ||
      !is.data.frame(ledger$attempts) ||
      !identical(names(ledger$attempts), attempt_fields) ||
      anyDuplicated(ledger$attempts$attempt_id))
    stop("Malformed V021 attempt ledger.")
  assert_v021_truth_free_schema(ledger, "attempt ledger")
  if (!is.null(manifest) && !identical(ledger$manifest_hash, manifest$manifest_hash))
    stop("Attempt ledger provenance conflicts with the immutable manifest.")
  invisible(TRUE)
}

validate_v021_state_transition <- function(from, to) {
  if (!from %in% names(v021_legal_state_transitions) ||
      !to %in% v021_legal_state_transitions[[from]])
    stop("Illegal V021 task-state transition: ", from, " -> ", to)
  invisible(TRUE)
}

write_v021_task_status_atomic <- function(status, path, manifest) {
  validator <- function(x) validate_v021_task_status(x, manifest)
  .v021_atomic_write_validated_rds(status, path, validator)
}

write_v021_attempt_ledger_atomic <- function(ledger, path, manifest) {
  validator <- function(x) validate_v021_attempt_ledger(x, manifest)
  .v021_atomic_write_validated_rds(ledger, path, validator)
}

read_v021_task_status <- function(path, manifest) {
  if (!file.exists(path)) stop("Task status is missing: ", path)
  x <- tryCatch(readRDS(path), error = identity)
  if (inherits(x, "error")) stop("Task status is malformed: ", path)
  validate_v021_task_status(x, manifest); x
}

read_v021_attempt_ledger <- function(path, manifest) {
  if (!file.exists(path)) stop("Attempt ledger is missing: ", path)
  x <- tryCatch(readRDS(path), error = identity)
  if (inherits(x, "error")) stop("Attempt ledger is malformed: ", path)
  validate_v021_attempt_ledger(x, manifest); x
}

build_v021_completion_artifact <- function(
    manifest, task_ordinal, result_root, posterior_summary, diagnostics,
    psp_lfsr, predictive_eligibility, subject_elpd, seed_split_provenance,
    failure_information = list()) {
  validate_v021_execution_manifest(manifest)
  task <- manifest$tasks[match(as.integer(task_ordinal), manifest$tasks$task_ordinal),,
                         drop = FALSE]
  if (nrow(task) != 1L) stop("Unknown completion task ordinal.")
  body <- list(
    artifact_schema = v021_completion_artifact_schema,
    manifest_hash = manifest$manifest_hash,
    configuration_hash = manifest$configuration_hash,
    result_root = normalizePath(result_root, mustWork = FALSE),
    dataset_id = manifest$dataset_id,
    directed_task_id = task$directed_task_id[[1L]],
    task_ordinal = task$task_ordinal[[1L]], source = task$source[[1L]],
    target = task$target[[1L]], direction_seed = task$direction_seed[[1L]],
    finalized = TRUE, posterior_summary = .v021_checkpoint_copy(posterior_summary),
    diagnostics = .v021_checkpoint_copy(diagnostics),
    psp_lfsr = .v021_checkpoint_copy(psp_lfsr),
    predictive_eligibility = predictive_eligibility,
    subject_elpd = .v021_checkpoint_copy(subject_elpd),
    seed_split_provenance = .v021_checkpoint_copy(seed_split_provenance),
    failure_information = .v021_checkpoint_copy(failure_information)
  )
  artifact <- c(body, list(artifact_hash = v021_sha256(body)))
  validate_v021_completion_artifact(artifact, manifest, result_root)
  artifact
}

validate_v021_completion_artifact <- function(artifact, manifest, result_root) {
  validate_v021_execution_manifest(manifest)
  required <- c("artifact_schema", "manifest_hash", "configuration_hash",
                "result_root", "dataset_id", "directed_task_id", "task_ordinal",
                "source", "target", "direction_seed", "finalized",
                "posterior_summary", "diagnostics", "psp_lfsr",
                "predictive_eligibility", "subject_elpd", "seed_split_provenance",
                "failure_information", "artifact_hash")
  if (!is.list(artifact) || !identical(names(artifact), required) ||
      !identical(artifact$artifact_schema, v021_completion_artifact_schema) ||
      !identical(artifact$finalized, TRUE))
    stop("Invalid V021 completion artifact schema or finalization state.")
  assert_v021_truth_free_schema(artifact, "completion artifact")
  if (!identical(artifact$manifest_hash, manifest$manifest_hash) ||
      !identical(artifact$configuration_hash, manifest$configuration_hash) ||
      !identical(artifact$result_root, normalizePath(result_root, mustWork = FALSE)))
    stop("Completion artifact provenance conflicts with manifest or result root.")
  task <- manifest$tasks[match(artifact$directed_task_id,
                               manifest$tasks$directed_task_id), , drop = FALSE]
  if (nrow(task) != 1L ||
      !identical(as.integer(artifact$task_ordinal), task$task_ordinal[[1L]]) ||
      !identical(artifact$source, task$source[[1L]]) ||
      !identical(artifact$target, task$target[[1L]]) ||
      !identical(as.integer(artifact$direction_seed), task$direction_seed[[1L]]))
    stop("Completion artifact task identity conflicts with manifest.")
  if (!is.list(artifact$posterior_summary) || !length(artifact$posterior_summary) ||
      !is.list(artifact$diagnostics) || !length(artifact$diagnostics) ||
      !is.list(artifact$psp_lfsr) || !all(c("PSP", "LFSR") %in% names(artifact$psp_lfsr)) ||
      !is.logical(artifact$predictive_eligibility) ||
      length(artifact$predictive_eligibility) != 1L ||
      !is.list(artifact$seed_split_provenance) ||
      !all(c("direction_seed", "kfold_seed") %in% names(artifact$seed_split_provenance)))
    stop("Completion artifact lacks required durable inference outputs.")
  expected_hash <- v021_sha256(artifact[setdiff(required, "artifact_hash")])
  if (!identical(artifact$artifact_hash, expected_hash))
    stop("Completion artifact hash mismatch.")
  invisible(TRUE)
}

write_v021_completion_artifact_atomic <- function(artifact, path, manifest,
                                                   result_root) {
  validator <- function(x) validate_v021_completion_artifact(x, manifest, result_root)
  if (file.exists(path)) {
    existing <- tryCatch(readRDS(path), error = identity)
    if (inherits(existing, "error")) stop("Existing completion artifact is malformed.")
    validator(existing)
    if (!identical(existing$artifact_hash, artifact$artifact_hash))
      stop("Refusing to overwrite an immutable completed artifact.")
    return(invisible(path))
  }
  .v021_atomic_write_validated_rds(artifact, path, validator)
}

initialize_v021_execution_root <- function(manifest, result_root) {
  paths <- v021_execution_paths(result_root)
  dir.create(paths$result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$artifact_root, recursive = TRUE, showWarnings = FALSE)
  write_v021_execution_manifest_once(manifest, paths$manifest)
  stored <- read_v021_execution_manifest(paths$manifest)
  compare_v021_execution_manifests(stored, manifest)
  if (!file.exists(paths$status))
    write_v021_task_status_atomic(new_v021_task_status(manifest), paths$status, manifest)
  if (!file.exists(paths$attempts))
    write_v021_attempt_ledger_atomic(new_v021_attempt_ledger(manifest),
                                     paths$attempts, manifest)
  invisible(paths)
}

plan_v021_execution_resume <- function(manifest, result_root, maximum_attempts,
                                        persist_reconciliation = FALSE,
                                        timestamp = function()
                                          format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  paths <- v021_execution_paths(result_root)
  stored <- read_v021_execution_manifest(paths$manifest)
  compare_v021_execution_manifests(stored, manifest)
  status <- read_v021_task_status(paths$status, stored)
  ledger <- read_v021_attempt_ledger(paths$attempts, stored)
  maximum_attempts <- as.integer(maximum_attempts)
  if (length(maximum_attempts) != 1L || is.na(maximum_attempts) || maximum_attempts < 1L)
    stop("maximum_attempts must be a positive integer.")
  changed <- FALSE
  for (i in seq_len(nrow(status$tasks))) {
    state <- status$tasks$state[[i]]
    artifact_path <- status$tasks$artifact_path[[i]]
    artifact_valid <- FALSE
    if (!is.na(artifact_path) && nzchar(artifact_path) && file.exists(artifact_path)) {
      artifact <- tryCatch(readRDS(artifact_path), error = identity)
      artifact_valid <- !inherits(artifact, "error") && isTRUE(tryCatch({
        validate_v021_completion_artifact(artifact, stored, result_root); TRUE
      }, error = function(e) FALSE))
    }
    if (identical(state, "completed") && !artifact_valid)
      stop("Completed task lacks a valid durable artifact: ",
           status$tasks$directed_task_id[[i]])
    if (identical(state, "running")) {
      attempt_id <- status$tasks$last_attempt_id[[i]]
      ledger_index <- match(attempt_id, ledger$attempts$attempt_id)
      if (is.na(ledger_index) || !is.na(ledger$attempts$ending_state[[ledger_index]]))
        stop("Running task has ambiguous or missing attempt history: ",
             status$tasks$directed_task_id[[i]])
      if (artifact_valid) {
        status$tasks$state[[i]] <- "completed"
        status$tasks$artifact_validation[[i]] <- "validated_after_interruption"
        ledger$attempts$ending_state[[ledger_index]] <- "completed"
        ledger$attempts$artifact_path[[ledger_index]] <- artifact_path
        ledger$attempts$artifact_validation[[ledger_index]] <-
          "validated_after_interruption"
      } else {
        status$tasks$state[[i]] <- "failed"
        status$tasks$terminal_reason[[i]] <- "interrupted_running_without_valid_artifact"
        status$tasks$artifact_validation[[i]] <- "absent_or_invalid"
        ledger$attempts$ending_state[[ledger_index]] <- "failed"
        ledger$attempts$terminal_reason[[ledger_index]] <-
          "interrupted_running_without_valid_artifact"
        ledger$attempts$failure_class[[ledger_index]] <- "interrupted_attempt"
        ledger$attempts$artifact_validation[[ledger_index]] <- "absent_or_invalid"
      }
      status$tasks$updated_timestamp[[i]] <- timestamp()
      ledger$attempts$end_timestamp[[ledger_index]] <- status$tasks$updated_timestamp[[i]]
      changed <- TRUE
    }
  }
  runnable <- which(status$tasks$state == "pending" |
                    (status$tasks$state == "failed" &
                     status$tasks$attempt_count < maximum_attempts))
  completed <- which(status$tasks$state == "completed")
  if (persist_reconciliation && changed) {
    write_v021_task_status_atomic(status, paths$status, stored)
    write_v021_attempt_ledger_atomic(ledger, paths$attempts, stored)
  }
  list(manifest = stored, status = status, attempt_ledger = ledger,
       runnable_indices = runnable, completed_indices = completed,
       reconciliation_changed = changed, paths = paths)
}

compact_v021_execution_status <- function(plan) {
  states <- table(factor(plan$status$tasks$state, levels = v021_execution_states))
  data.frame(
    manifest_hash = plan$manifest$manifest_hash,
    pending = unname(states[["pending"]]), running = unname(states[["running"]]),
    completed = unname(states[["completed"]]), failed = unname(states[["failed"]]),
    runnable = length(plan$runnable_indices), stringsAsFactors = FALSE)
}

start_v021_task_attempt <- function(status, ledger, manifest, task_ordinal,
                                     timestamp, worker_provenance = NA_character_) {
  validate_v021_task_status(status, manifest)
  validate_v021_attempt_ledger(ledger, manifest)
  i <- match(as.integer(task_ordinal), status$tasks$task_ordinal)
  if (is.na(i)) stop("Unknown task ordinal.")
  from <- status$tasks$state[[i]]
  validate_v021_state_transition(from, "running")
  attempt_number <- status$tasks$attempt_count[[i]] + 1L
  attempt_id <- sprintf("%s-attempt-%03d",
                        status$tasks$directed_task_id[[i]], attempt_number)
  if (attempt_id %in% ledger$attempts$attempt_id)
    stop("Attempt ledger already contains this attempt identity.")
  task <- manifest$tasks[i, , drop = FALSE]
  row <- data.frame(
    attempt_id = attempt_id, directed_task_id = task$directed_task_id,
    task_ordinal = task$task_ordinal, attempt_number = attempt_number,
    starting_state = from, ending_state = NA_character_,
    start_timestamp = as.character(timestamp), end_timestamp = NA_character_,
    direction_seed = task$direction_seed,
    worker_provenance = as.character(worker_provenance),
    terminal_reason = NA_character_, failure_class = NA_character_,
    artifact_path = NA_character_, artifact_validation = "not_checked",
    stringsAsFactors = FALSE)
  ledger$attempts <- rbind(ledger$attempts, row)
  status$tasks$state[[i]] <- "running"
  status$tasks$attempt_count[[i]] <- attempt_number
  status$tasks$last_attempt_id[[i]] <- attempt_id
  status$tasks$updated_timestamp[[i]] <- as.character(timestamp)
  status$tasks$terminal_reason[[i]] <- NA_character_
  status$tasks$artifact_path[[i]] <- NA_character_
  status$tasks$artifact_validation[[i]] <- "not_checked"
  validate_v021_task_status(status, manifest)
  validate_v021_attempt_ledger(ledger, manifest)
  list(status = status, ledger = ledger, attempt_id = attempt_id)
}

finish_v021_task_attempt <- function(status, ledger, manifest, task_ordinal,
                                      ending_state, timestamp,
                                      terminal_reason = NA_character_,
                                      failure_class = NA_character_,
                                      artifact_path = NA_character_,
                                      artifact_validation = "not_checked",
                                      result_root = NULL) {
  validate_v021_task_status(status, manifest)
  validate_v021_attempt_ledger(ledger, manifest)
  i <- match(as.integer(task_ordinal), status$tasks$task_ordinal)
  if (is.na(i)) stop("Unknown task ordinal.")
  validate_v021_state_transition(status$tasks$state[[i]], ending_state)
  attempt_id <- status$tasks$last_attempt_id[[i]]
  j <- match(attempt_id, ledger$attempts$attempt_id)
  if (is.na(j) || !is.na(ledger$attempts$ending_state[[j]]))
    stop("Current attempt ledger row is missing or already terminal.")
  if (identical(ending_state, "completed")) {
    if (is.null(result_root) || is.na(artifact_path) || !file.exists(artifact_path))
      stop("Completion requires an existing durable artifact and result root.")
    artifact <- tryCatch(readRDS(artifact_path), error = identity)
    if (inherits(artifact, "error")) stop("Completion artifact is unreadable.")
    validate_v021_completion_artifact(artifact, manifest, result_root)
    artifact_validation <- "validated"
  }
  if (identical(ending_state, "failed") &&
      (is.na(terminal_reason) || !nzchar(terminal_reason)))
    stop("Failed attempts require a terminal reason.")
  ledger$attempts$ending_state[[j]] <- ending_state
  ledger$attempts$end_timestamp[[j]] <- as.character(timestamp)
  ledger$attempts$terminal_reason[[j]] <- terminal_reason
  ledger$attempts$failure_class[[j]] <- failure_class
  ledger$attempts$artifact_path[[j]] <- artifact_path
  ledger$attempts$artifact_validation[[j]] <- artifact_validation
  status$tasks$state[[i]] <- ending_state
  status$tasks$updated_timestamp[[i]] <- as.character(timestamp)
  status$tasks$terminal_reason[[i]] <- terminal_reason
  status$tasks$artifact_path[[i]] <- artifact_path
  status$tasks$artifact_validation[[i]] <- artifact_validation
  validate_v021_task_status(status, manifest)
  validate_v021_attempt_ledger(ledger, manifest)
  list(status = status, ledger = ledger)
}
