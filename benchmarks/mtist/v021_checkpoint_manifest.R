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
