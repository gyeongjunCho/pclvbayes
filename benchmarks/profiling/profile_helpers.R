profiling_timing_schema <- function() {
  c("run", "phase", "calls", "elapsed_seconds", "user_seconds", "system_seconds",
    "measurement_scope", "available", "note")
}

new_timing_row <- function(run, phase, elapsed = NA_real_, user = NA_real_,
                           system = NA_real_, calls = 1L, scope = "exclusive",
                           available = is.finite(elapsed), note = NA_character_) {
  out <- data.frame(
    run = as.character(run), phase = as.character(phase), calls = as.integer(calls),
    elapsed_seconds = as.numeric(elapsed), user_seconds = as.numeric(user),
    system_seconds = as.numeric(system), measurement_scope = as.character(scope),
    available = as.logical(available), note = as.character(note),
    stringsAsFactors = FALSE
  )
  out[profiling_timing_schema()]
}

time_expression <- function(run, phase, expr, scope = "exclusive") {
  started <- proc.time()
  value <- force(expr)
  elapsed <- proc.time() - started
  list(value = value, timing = new_timing_row(
    run, phase, elapsed[["elapsed"]], elapsed[["user.self"]], elapsed[["sys.self"]],
    scope = scope
  ))
}

profile_disabled <- function(expr, enabled = FALSE) {
  if (!isTRUE(enabled)) return(force(expr))
  time_expression("internal", "expression", expr)$value
}

read_proc_status_kb <- function(field) {
  path <- "/proc/self/status"
  if (!file.exists(path)) return(NA_real_)
  line <- grep(paste0("^", field, ":"), readLines(path, warn = FALSE), value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(sub(paste0("^", field, ":\\s*([0-9]+).*"), "\\1", line[[1L]]))
}

platform_memory <- function() {
  data.frame(rss_kb = read_proc_status_kb("VmRSS"), peak_rss_kb = read_proc_status_kb("VmHWM"))
}

directory_snapshot <- function(path) {
  if (!dir.exists(path)) return(data.frame(path = path, files = 0L, directories = 0L, bytes = 0))
  entries <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                        include.dirs = TRUE, no.. = TRUE)
  info <- if (length(entries)) file.info(entries) else data.frame(size = numeric(), isdir = logical())
  data.frame(path = normalizePath(path), files = sum(!info$isdir), directories = sum(info$isdir),
             bytes = sum(info$size[!info$isdir], na.rm = TRUE))
}

system_time_available <- function() {
  path <- Sys.which("time")
  if (!nzchar(path) || !identical(normalizePath(path), "/usr/bin/time")) NA_character_ else path
}

object_sizes <- function(objects) {
  data.frame(object = names(objects), bytes = vapply(objects, object.size, numeric(1)))
}

summarise_rprof <- function(path, top_n = 20L) {
  if (!file.exists(path) || !file.info(path)$size) return(data.frame())
  tab <- summaryRprof(path)$by.total
  if (!nrow(tab)) return(data.frame())
  tab$function_name <- rownames(tab)
  rownames(tab) <- NULL
  utils::head(tab[order(tab$total.time, decreasing = TRUE), ], top_n)
}

scientific_result_signature <- function(fit) {
  raw <- fit$raw
  cols <- intersect(c("i", "j", "direction_ok_ij", "direction_ok_ji",
                      "a_ij_mean", "a_ji_mean", "diagnostic_class_ij", "diagnostic_class_ji"), names(raw))
  raw[, cols, drop = FALSE]
}

