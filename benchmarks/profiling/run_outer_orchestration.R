args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))) else getwd()
repo <- normalizePath(file.path(script_dir, "../.."))
profile_library <- tempfile("pclv-outer-library-")
dir.create(profile_library)
on.exit(unlink(profile_library, recursive = TRUE), add = TRUE)
withr::with_libpaths(profile_library, action = "prefix",
  devtools::install(repo, quiet = TRUE, upgrade = FALSE,
                    dependencies = FALSE, keep_source = TRUE))
.libPaths(c(profile_library, .libPaths()))
library(pclvbayes)

workers <- 2L
repetitions <- 5L
tasks <- pclvbayes:::.make_pair_tasks(c("species_0", "species_1", "species_2"), 20260801L)
task_payload_bytes <- length(serialize(tasks, NULL))

# Deterministic CPU work accompanies the wait representing an external CmdStan
# process. Production benefit is assessed separately by the real MTIST profile.
work <- function(task, durations) {
  started <- proc.time()[["elapsed"]]
  x <- matrix(seq_len(1600L) + task$task_index, 40L)
  checksum <- sum(crossprod(x))
  Sys.sleep(durations[[task$task_index]])
  list(task_index = task$task_index, seed = task$direction_seeds,
       duration = proc.time()[["elapsed"]] - started, checksum = checksum)
}

legacy_map <- function(tasks, durations) {
  run_task <- function(task) work(task, durations)
  furrr::future_map(
    tasks, function(task) run_task(task),
    .options = furrr::furrr_options(seed = TRUE, globals = TRUE, scheduling = 1))
}

optimized_map <- function(tasks, durations) {
  furrr::future_map(
    tasks, work, durations = durations,
    .options = furrr::furrr_options(
      seed = TRUE, globals = FALSE, packages = "pclvbayes", scheduling = Inf))
}

measure <- function(label, fun, durations) {
  elapsed <- numeric(repetitions)
  signatures <- vector("list", repetitions)
  for (r in seq_len(repetitions)) {
    timing <- system.time(value <- fun(tasks, durations))
    elapsed[[r]] <- timing[["elapsed"]]
    signatures[[r]] <- lapply(
      value[order(vapply(value, `[[`, integer(1), "task_index"))],
      function(x) x[c("task_index", "seed", "checksum")])
  }
  data.frame(
    case = label, median_elapsed_seconds = median(elapsed),
    min_elapsed_seconds = min(elapsed), max_elapsed_seconds = max(elapsed),
    repetitions = repetitions,
    deterministic = all(vapply(signatures[-1L], identical, logical(1), signatures[[1L]])))
}

probe_env <- list2env(list(tasks = tasks, durations = c(0.10, 0.10, 0.65),
                                legacy_map = legacy_map), parent = environment())
legacy_probe <- future::getGlobalsAndPackages(
  quote(legacy_map(tasks, durations)), envir = probe_env, globals = TRUE)
legacy_global_bytes <- attr(legacy_probe$globals, "total_size")
if (is.null(legacy_global_bytes) || !is.finite(legacy_global_bytes))
  legacy_global_bytes <- length(serialize(as.list(legacy_probe$globals), NULL))
legacy_global_names <- paste(names(legacy_probe$globals), collapse = ",")
explicit_payload_bytes <- length(serialize(list(
  tasks = tasks, durations = probe_env$durations, work = work), NULL))

old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::multisession, workers = workers)
legacy_map(tasks, c(0.01, 0.01, 0.01))
optimized_map(tasks, c(0.01, 0.01, 0.01))
results <- rbind(
  measure("legacy_equal", legacy_map, c(0.15, 0.15, 0.15)),
  measure("optimized_equal", optimized_map, c(0.15, 0.15, 0.15)),
  measure("legacy_uneven", legacy_map, c(0.10, 0.65, 0.65)),
  measure("optimized_uneven", optimized_map, c(0.10, 0.65, 0.65)))
results$task_payload_bytes <- task_payload_bytes
results$legacy_discovered_global_bytes <- legacy_global_bytes
results$optimized_explicit_payload_bytes <- explicit_payload_bytes
results$legacy_discovered_globals <- legacy_global_names; results$workers <- workers
results$submitted_tasks <- length(tasks)

out <- file.path(script_dir, "results", format(Sys.time(), "%Y%m%d-%H%M%S-outer"))
dir.create(out, recursive = TRUE)
utils::write.csv(results, file.path(out, "orchestration.csv"), row.names = FALSE)
jsonlite::write_json(list(
  repository_sha = system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
  timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  workers = workers, repetitions = repetitions, task_count = length(tasks),
  legacy = list(globals = TRUE, scheduling = 1),
  optimized = list(globals = FALSE, packages = "pclvbayes", scheduling = "Inf"),
  R = R.version.string, session = capture.output(sessionInfo())),
  file.path(out, "metadata.json"), pretty = TRUE, auto_unbox = TRUE)
cat("ORCHESTRATION_RESULT_DIR=", out, "\n", sep = "")
print(results)
