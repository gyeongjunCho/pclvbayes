# Benchmark-only resource and process-monitoring contract for ROADMAP V021-02.

v021_resource_policy_schema <- "v021_resource_policy_v1"
v021_preflight_states <- c(
  "passed", "failed_ceiling_exceeded", "failed_unverified_process_tree",
  "failed_thread_environment", "failed_invalid_policy", "failed_monitoring_error"
)

v021_thread_variables <- c(
  "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "STAN_NUM_THREADS",
  "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
  "RCPP_PARALLEL_NUM_THREADS"
)

v021_operation_types <- c(
  "main_fit", "retry_fit", "pathfinder", "kfold_fit", "confirmation_fit"
)

.v021_linux_cmdline_max_bytes <- 1024L * 1024L

.v021_resource_int <- function(x, name, positive = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || (positive && x < 1L) || (!positive && x < 0L))
    stop(name, " must be ", if (positive) "a positive" else "a non-negative",
         " integer.")
  as.integer(x)
}

.v021_operation_table <- function(main_chains, retry_chains, pathfinder_processes,
                                  kfold_chains, kfold_parallel_chains,
                                  confirmation_chains) {
  data.frame(
    operation = v021_operation_types,
    total_chains = as.integer(c(
      main_chains, retry_chains, 0L, kfold_chains, confirmation_chains
    )),
    simultaneous_chain_slots = as.integer(c(
      main_chains, retry_chains, 0L, kfold_parallel_chains, confirmation_chains
    )),
    cmdstan_process_slots = as.integer(c(
      main_chains, retry_chains, pathfinder_processes,
      kfold_parallel_chains, confirmation_chains
    )),
    execution_relationship = c(
      "direction_main", "sequential_retry_of_direction",
      "sequential_before_direction_main", "predictive_fold",
      "confirmation_main_only"
    ),
    stringsAsFactors = FALSE
  )
}

build_v021_resource_policy <- function(
    logical_host_threads = 12L, reserved_host_threads = 2L,
    maximum_active_cmdstan_chains = 10L, cpu_threads_per_active_chain = 1L,
    main_chains = 4L, retry_chains = main_chains,
    pathfinder_processes = 1L, pathfinder_num_paths = 8L,
    kfold_chains = main_chains, kfold_parallel_chains = 1L,
    confirmation_chains = main_chains, proposed_outer_concurrency = NULL,
    global_slot_scheduler = FALSE) {
  logical_host_threads <- .v021_resource_int(logical_host_threads, "logical_host_threads")
  reserved_host_threads <- .v021_resource_int(
    reserved_host_threads, "reserved_host_threads", positive = FALSE)
  maximum_active_cmdstan_chains <- .v021_resource_int(
    maximum_active_cmdstan_chains, "maximum_active_cmdstan_chains")
  cpu_threads_per_active_chain <- .v021_resource_int(
    cpu_threads_per_active_chain, "cpu_threads_per_active_chain")
  main_chains <- .v021_resource_int(main_chains, "main_chains")
  retry_chains <- .v021_resource_int(retry_chains, "retry_chains")
  pathfinder_processes <- .v021_resource_int(
    pathfinder_processes, "pathfinder_processes", positive = FALSE)
  pathfinder_num_paths <- .v021_resource_int(pathfinder_num_paths, "pathfinder_num_paths")
  kfold_chains <- .v021_resource_int(kfold_chains, "kfold_chains")
  kfold_parallel_chains <- .v021_resource_int(
    kfold_parallel_chains, "kfold_parallel_chains")
  confirmation_chains <- .v021_resource_int(confirmation_chains, "confirmation_chains")
  usable <- logical_host_threads - reserved_host_threads
  if (usable < 1L || usable != maximum_active_cmdstan_chains)
    stop("Usable host threads must equal the global CmdStan-chain ceiling.")
  if (cpu_threads_per_active_chain != 1L)
    stop("Exactly one CPU thread per active chain is required.")
  if (kfold_parallel_chains > kfold_chains)
    stop("kfold_parallel_chains cannot exceed kfold_chains.")
  if (!is.logical(global_slot_scheduler) || length(global_slot_scheduler) != 1L ||
      is.na(global_slot_scheduler)) stop("global_slot_scheduler must be one logical value.")
  if (isTRUE(global_slot_scheduler))
    stop("No verified global slot scheduler exists for V021-02.")

  operation_slots <- .v021_operation_table(
    main_chains, retry_chains, pathfinder_processes, kfold_chains,
    kfold_parallel_chains, confirmation_chains
  )
  binding_slots <- max(operation_slots$simultaneous_chain_slots,
                       operation_slots$cmdstan_process_slots)
  safe_outer <- floor(usable / binding_slots)
  if (is.null(proposed_outer_concurrency)) proposed_outer_concurrency <- safe_outer
  proposed_outer_concurrency <- .v021_resource_int(
    proposed_outer_concurrency, "proposed_outer_concurrency")
  if (proposed_outer_concurrency > safe_outer)
    stop("Proposed outer concurrency exceeds the global slot ceiling.")

  policy <- list(
    policy_schema = v021_resource_policy_schema,
    logical_host_threads = logical_host_threads,
    reserved_host_threads = reserved_host_threads,
    usable_chain_slots = as.integer(usable),
    maximum_active_cmdstan_chains = maximum_active_cmdstan_chains,
    cpu_threads_per_active_chain = cpu_threads_per_active_chain,
    numerical_library_threads = 1L,
    main_chains = main_chains,
    retry_chains = retry_chains,
    pathfinder_processes = pathfinder_processes,
    pathfinder_chain_slots = 0L,
    pathfinder_num_paths = pathfinder_num_paths,
    kfold_chains = kfold_chains,
    kfold_parallel_chains = kfold_parallel_chains,
    confirmation_chains = confirmation_chains,
    proposed_outer_concurrency = proposed_outer_concurrency,
    global_slot_scheduler = FALSE,
    operation_slots = operation_slots
  )
  validate_v021_resource_policy(policy)
  policy
}

validate_v021_resource_policy <- function(policy) {
  fields <- c(
    "policy_schema", "logical_host_threads", "reserved_host_threads",
    "usable_chain_slots", "maximum_active_cmdstan_chains",
    "cpu_threads_per_active_chain", "numerical_library_threads", "main_chains",
    "retry_chains", "pathfinder_processes", "pathfinder_chain_slots",
    "pathfinder_num_paths", "kfold_chains", "kfold_parallel_chains",
    "confirmation_chains", "proposed_outer_concurrency", "global_slot_scheduler",
    "operation_slots"
  )
  if (!is.list(policy) || !identical(names(policy), fields))
    stop("Invalid V021 resource-policy schema.")
  if (!identical(policy$policy_schema, v021_resource_policy_schema))
    stop("Invalid V021 resource-policy version.")
  for (nm in setdiff(fields, c("policy_schema", "global_slot_scheduler", "operation_slots")))
    .v021_resource_int(policy[[nm]], nm, positive = nm != "pathfinder_chain_slots" &&
                         nm != "pathfinder_processes" && nm != "reserved_host_threads")
  if (policy$logical_host_threads != 12L || policy$reserved_host_threads != 2L ||
      policy$usable_chain_slots != 10L || policy$maximum_active_cmdstan_chains != 10L)
    stop("Policy violates the exact 12/2/10 host contract.")
  if (policy$logical_host_threads - policy$reserved_host_threads !=
      policy$usable_chain_slots)
    stop("Policy host-thread arithmetic is contradictory.")
  if (policy$cpu_threads_per_active_chain != 1L ||
      policy$numerical_library_threads != 1L)
    stop("Policy requires exactly one thread per chain and numerical library.")
  if (!identical(policy$global_slot_scheduler, FALSE))
    stop("V021-02 has no verified global slot scheduler.")
  expected_operations <- .v021_operation_table(
    policy$main_chains, policy$retry_chains, policy$pathfinder_processes,
    policy$kfold_chains, policy$kfold_parallel_chains, policy$confirmation_chains
  )
  if (!identical(policy$operation_slots, expected_operations))
    stop("Operation-slot table contradicts the policy.")
  if (policy$pathfinder_chain_slots != 0L)
    stop("Inspected Pathfinder behavior consumes no MCMC chain slots.")
  if (policy$pathfinder_processes != 1L || policy$pathfinder_num_paths != 8L)
    stop("Policy contradicts inspected canonical Pathfinder behavior.")
  if (policy$retry_chains != policy$main_chains ||
      policy$kfold_chains != policy$main_chains ||
      policy$kfold_parallel_chains != 1L ||
      policy$confirmation_chains != policy$main_chains)
    stop("Policy contradicts canonical main, retry, K-fold, or confirmation chain behavior.")
  if (policy$proposed_outer_concurrency *
      max(policy$operation_slots$simultaneous_chain_slots,
          policy$operation_slots$cmdstan_process_slots) > policy$usable_chain_slots)
    stop("Proposed outer concurrency is mathematically unsafe.")
  invisible(TRUE)
}

build_v021_operation_spec <- function(operations, requested_outer_concurrency) {
  if (!is.character(operations) || !length(operations) || anyNA(operations) ||
      any(!operations %in% v021_operation_types)) stop("Invalid operation types.")
  requested_outer_concurrency <- .v021_resource_int(
    requested_outer_concurrency, "requested_outer_concurrency")
  list(operations = unique(operations),
       requested_outer_concurrency = requested_outer_concurrency)
}

derive_safe_outer_concurrency <- function(policy, operation_spec) {
  validate_v021_resource_policy(policy)
  if (!is.list(operation_spec) ||
      !identical(names(operation_spec), c("operations", "requested_outer_concurrency")))
    stop("Invalid operation specification.")
  spec <- build_v021_operation_spec(
    operation_spec$operations, operation_spec$requested_outer_concurrency)
  slots <- policy$operation_slots[
    match(spec$operations, policy$operation_slots$operation), , drop = FALSE]
  per_fit_chain_slots <- max(slots$simultaneous_chain_slots)
  per_fit_process_slots <- max(slots$cmdstan_process_slots)
  chain_capacity <- if (per_fit_chain_slots > 0L)
    floor(policy$maximum_active_cmdstan_chains / per_fit_chain_slots) else Inf
  process_capacity <- if (per_fit_process_slots > 0L)
    floor(policy$usable_chain_slots / per_fit_process_slots) else Inf
  safe <- as.integer(min(chain_capacity, process_capacity))
  if (safe < 1L) stop("Operation cannot fit within the global resource ceiling.")
  requested <- spec$requested_outer_concurrency
  projected_chains <- requested * per_fit_chain_slots
  projected_processes <- requested * per_fit_process_slots
  if (requested > safe || projected_chains > policy$maximum_active_cmdstan_chains ||
      projected_processes > policy$usable_chain_slots)
    stop("Requested outer concurrency exceeds the binding global ceiling.")
  list(
    policy_schema = policy$policy_schema,
    operations = spec$operations,
    per_fit_simultaneous_chain_slots = as.integer(per_fit_chain_slots),
    per_fit_cmdstan_process_slots = as.integer(per_fit_process_slots),
    maximum_active_cmdstan_chains = policy$maximum_active_cmdstan_chains,
    usable_process_slots = policy$usable_chain_slots,
    safe_outer_concurrency = safe,
    requested_outer_concurrency = requested,
    projected_active_cmdstan_chains = as.integer(projected_chains),
    projected_active_cmdstan_processes = as.integer(projected_processes),
    global_slot_scheduler = FALSE,
    compliant = TRUE
  )
}

validate_v021_operation_plan <- function(policy, operation_plan) {
  validate_v021_resource_policy(policy)
  if (!is.data.frame(operation_plan) ||
      !identical(names(operation_plan), c("operation", "concurrent_instances")) ||
      !nrow(operation_plan) || anyNA(operation_plan))
    stop("Invalid overlapping operation plan.")
  if (any(!operation_plan$operation %in% v021_operation_types))
    stop("Unknown operation in overlapping plan.")
  if (!is.numeric(operation_plan$concurrent_instances) ||
      any(operation_plan$concurrent_instances < 0L) ||
      any(operation_plan$concurrent_instances != as.integer(operation_plan$concurrent_instances)))
    stop("Invalid concurrent_instances.")
  slots <- policy$operation_slots[
    match(operation_plan$operation, policy$operation_slots$operation), , drop = FALSE]
  active_chains <- sum(operation_plan$concurrent_instances * slots$simultaneous_chain_slots)
  active_processes <- sum(operation_plan$concurrent_instances * slots$cmdstan_process_slots)
  if (active_chains > policy$maximum_active_cmdstan_chains ||
      active_processes > policy$usable_chain_slots)
    stop("Overlapping operation plan exceeds the global ceiling.")
  list(active_cmdstan_chains = as.integer(active_chains),
       active_cmdstan_processes = as.integer(active_processes), compliant = TRUE)
}

v021_single_thread_environment <- function() {
  stats::setNames(rep("1", length(v021_thread_variables)), v021_thread_variables)
}

validate_v021_single_thread_environment <- function(environment) {
  required <- v021_single_thread_environment()
  if (is.null(names(environment)) || anyDuplicated(names(environment)) ||
      !all(names(required) %in% names(environment)))
    stop("Required numerical-thread environment variables are missing.")
  observed <- as.character(environment[names(required)])
  if (anyNA(observed) || any(observed != "1"))
    stop("Every governed numerical-thread variable must equal exactly one.")
  invisible(TRUE)
}

with_v021_single_thread_environment <- function(code) {
  if (!is.function(code)) stop("code must be a function.")
  environment <- v021_single_thread_environment()
  validate_v021_single_thread_environment(environment)
  withr::with_envvar(environment, code())
}

.v021_descendant_pids <- function(snapshot, root_pid) {
  descendants <- as.integer(root_pid)
  repeat {
    children <- snapshot$pid[snapshot$ppid %in% descendants]
    expanded <- unique(c(descendants, children))
    if (identical(sort(expanded), sort(descendants))) break
    descendants <- expanded
  }
  descendants
}

.v021_decode_linux_cmdline <- function(command_raw) {
  if (!is.raw(command_raw) || !length(command_raw))
    stop("Linux cmdline must be a non-empty raw vector.")
  nul <- which(command_raw == as.raw(0L))
  if (!length(nul) || tail(nul, 1L) != length(command_raw))
    stop("Linux cmdline must end at a NUL argument boundary.")
  starts <- c(1L, head(nul, -1L) + 1L)
  ends <- nul - 1L
  argv <- Map(function(start, end) {
    if (end < start) return("")
    rawToChar(command_raw[start:end])
  }, starts, ends)
  argv <- unlist(argv, use.names = FALSE)
  if (!is.character(argv) || !length(argv) || anyNA(argv) || !nzchar(argv[[1L]]))
    stop("Linux cmdline does not contain a valid executable argument.")
  list(argv = argv, command = paste(argv, collapse = " "))
}

.v021_read_linux_cmdline <- function(path,
                                      max_bytes = .v021_linux_cmdline_max_bytes) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("Linux cmdline path must be one non-empty string.")
  if (!is.numeric(max_bytes) || length(max_bytes) != 1L || is.na(max_bytes) ||
      !is.finite(max_bytes) || max_bytes < 1L || max_bytes != as.integer(max_bytes))
    stop("Linux cmdline maximum must be a positive integer.")
  connection <- tryCatch(suppressWarnings(file(path, open = "rb")), error = identity)
  if (inherits(connection, "error"))
    stop("Linux cmdline is unreadable: ", conditionMessage(connection))
  on.exit(close(connection), add = TRUE)
  chunks <- list()
  total <- 0L
  repeat {
    remaining <- as.integer(max_bytes) - total
    request <- min(4096L, remaining + 1L)
    chunk <- tryCatch(readBin(connection, what = "raw", n = request), error = identity)
    if (inherits(chunk, "error"))
      stop("Linux cmdline is unreadable: ", conditionMessage(chunk))
    if (!length(chunk)) break
    total <- total + length(chunk)
    if (total > max_bytes)
      stop("Linux cmdline exceeded the bounded read limit and may be truncated.")
    chunks[[length(chunks) + 1L]] <- chunk
  }
  command_raw <- if (length(chunks)) do.call(c, chunks) else raw()
  if (!length(command_raw)) stop("Linux cmdline is empty.")
  command_raw
}

.v021_snapshot_argv <- function(snapshot) {
  if ("argv" %in% names(snapshot)) {
    valid <- vapply(snapshot$argv, function(x)
      is.character(x) && length(x) > 0L && !anyNA(x) && nzchar(x[[1L]]), logical(1))
    if (!all(valid)) stop("Malformed decoded process argv.")
    return(snapshot$argv)
  }
  lapply(snapshot$command, function(x) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) character()
    else strsplit(x, "[[:space:]]+", perl = TRUE)[[1L]]
  })
}

classify_v021_process_snapshot <- function(snapshot, root_pid,
                                            known_model_executables = character(),
                                            operation_map = NULL) {
  required <- c("timestamp", "pid", "ppid", "command", "executable",
                "potential_cmdstan", "readable")
  allowed <- append(required, "argv", after = 4L)
  if (!is.data.frame(snapshot) ||
      !(identical(names(snapshot), required) || identical(names(snapshot), allowed)) ||
      anyDuplicated(snapshot$pid) || !root_pid %in% snapshot$pid)
    stop("Invalid process snapshot.")
  descendants <- .v021_descendant_pids(snapshot, root_pid)
  out <- snapshot[snapshot$pid %in% descendants, , drop = FALSE]
  out$is_descendant <- out$pid != root_pid
  argv <- .v021_snapshot_argv(out)
  sample_argument <- vapply(argv, function(x) "method=sample" %in% x, logical(1))
  pathfinder_argument <- vapply(argv, function(x) "method=pathfinder" %in% x, logical(1))
  known_executable <- out$executable %in% known_model_executables
  sample_process <- known_executable & sample_argument
  pathfinder_process <- known_executable & pathfinder_argument
  r_process <- grepl("(^|/)(R|Rscript)([[:space:]]|$)", out$command)
  out$classification <- "other_descendant"
  out$classification[out$pid == root_pid] <- "parent_r"
  out$classification[out$is_descendant & r_process] <- "outer_worker"
  out$classification[out$is_descendant & pathfinder_process] <- "pathfinder_process"
  out$classification[out$is_descendant & sample_process] <- "cmdstan_chain"
  unknown <- out$is_descendant & (out$potential_cmdstan | known_executable) &
    !sample_process & !pathfinder_process
  out$classification[unknown] <- "unknown_potential_cmdstan"
  out$is_active_cmdstan_chain <- out$classification == "cmdstan_chain"
  out$operation <- NA_character_
  out$task_id <- NA_character_
  out$direction_id <- NA_character_
  if (!is.null(operation_map)) {
    map_fields <- c("pid", "operation", "task_id", "direction_id")
    if (!is.data.frame(operation_map) || !identical(names(operation_map), map_fields) ||
        anyDuplicated(operation_map$pid)) stop("Invalid process operation map.")
    matched <- match(out$pid, operation_map$pid)
    present <- !is.na(matched)
    out$operation[present] <- operation_map$operation[matched[present]]
    out$task_id[present] <- as.character(operation_map$task_id[matched[present]])
    out$direction_id[present] <- as.character(operation_map$direction_id[matched[present]])
  }
  out
}

monitor_v021_process_snapshots <- function(snapshots, policy, root_pid,
                                            known_model_executables = character(),
                                            operation_map = NULL) {
  validate_v021_resource_policy(policy)
  if (!is.list(snapshots) || !length(snapshots))
    return(list(monitoring_state = "monitoring_error", reason = "no_process_snapshots",
                records = data.frame(), observed_peak_active_cmdstan_chains = NA_integer_,
                observed_peak_active_cmdstan_processes = NA_integer_,
                configured_ceiling = policy$maximum_active_cmdstan_chains,
                compliance_status = "unverified"))
  classified <- tryCatch(lapply(snapshots, classify_v021_process_snapshot,
    root_pid = root_pid, known_model_executables = known_model_executables,
    operation_map = operation_map), error = identity)
  if (inherits(classified, "error"))
    return(list(monitoring_state = "monitoring_error",
                reason = conditionMessage(classified), records = data.frame(),
                observed_peak_active_cmdstan_chains = NA_integer_,
                observed_peak_active_cmdstan_processes = NA_integer_,
                configured_ceiling = policy$maximum_active_cmdstan_chains,
                compliance_status = "unverified"))
  records <- do.call(rbind, classified)
  counts <- vapply(classified, function(x) sum(x$is_active_cmdstan_chain), integer(1))
  process_counts <- vapply(classified, function(x)
    sum(x$classification %in% c("cmdstan_chain", "pathfinder_process")), integer(1))
  peak <- max(counts)
  process_peak <- max(process_counts)
  complete <- all(vapply(snapshots, function(x) all(x$readable), logical(1)))
  unknown <- any(records$classification == "unknown_potential_cmdstan")
  state <- if (!complete || unknown) "unverified_process_tree" else "verified"
  reason <- if (!complete) "unreadable_descendant_process" else if (unknown)
    "unknown_potential_cmdstan_descendant" else NA_character_
  records$active_cmdstan_chain_count <- rep(counts, vapply(classified, nrow, integer(1)))
  records$active_cmdstan_process_count <- rep(
    process_counts, vapply(classified, nrow, integer(1)))
  records$observed_peak_active_cmdstan_chain_count <- peak
  records$observed_peak_active_cmdstan_process_count <- process_peak
  records$configured_ceiling <- policy$maximum_active_cmdstan_chains
  records$compliance_status <- if (state != "verified") "unverified" else if (
    peak <= policy$maximum_active_cmdstan_chains &&
      process_peak <= policy$usable_chain_slots) "compliant" else "exceeded"
  list(
    monitoring_state = state, reason = reason, records = records,
    observed_peak_active_cmdstan_chains = as.integer(peak),
    observed_peak_active_cmdstan_processes = as.integer(process_peak),
    configured_ceiling = policy$maximum_active_cmdstan_chains,
    compliance_status = unique(records$compliance_status)[[1L]]
  )
}

capture_v021_linux_process_snapshot <- function(
    root_pid = Sys.getpid(), known_model_executables = character(),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), proc_root = "/proc") {
  if (!dir.exists(proc_root)) stop("Linux /proc process monitoring is unavailable.")
  entries <- list.files(proc_root, pattern = "^[0-9]+$", full.names = FALSE)
  rows <- lapply(entries, function(pid_text) {
    stat_path <- file.path(proc_root, pid_text, "stat")
    cmd_path <- file.path(proc_root, pid_text, "cmdline")
    stat <- tryCatch(readLines(stat_path, warn = FALSE, n = 1L), error = function(e) character())
    if (!length(stat)) return(NULL)
    ppid <- suppressWarnings(as.integer(sub("^[0-9]+ \\(.*\\) [A-Z] ([0-9]+).*$", "\\1", stat)))
    command_raw <- tryCatch(.v021_read_linux_cmdline(cmd_path),
                            error = function(e) raw())
    decoded <- tryCatch(.v021_decode_linux_cmdline(command_raw), error = function(e) NULL)
    readable <- !is.null(decoded)
    argv <- if (readable) decoded$argv else character()
    command <- if (readable) decoded$command else ""
    executable <- tryCatch(Sys.readlink(file.path(proc_root, pid_text, "exe")),
                           error = function(e) "")
    potential <- any(grepl("^method=", argv)) || executable %in% known_model_executables ||
      grepl("cmdstan", command, ignore.case = TRUE)
    data.frame(timestamp = as.character(timestamp), pid = as.integer(pid_text),
               ppid = ppid, command = command, argv = I(list(argv)),
               executable = executable, potential_cmdstan = potential,
               readable = readable, stringsAsFactors = FALSE)
  })
  snapshot <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(snapshot) || !root_pid %in% snapshot$pid)
    stop("Root process could not be established in /proc.")
  descendants <- .v021_descendant_pids(snapshot, root_pid)
  snapshot <- snapshot[snapshot$pid %in% descendants, , drop = FALSE]
  rownames(snapshot) <- NULL
  snapshot
}

monitor_v021_process_tree <- function(policy, root_pid = Sys.getpid(), samples = 1L,
                                       interval_seconds = 0.1,
                                       known_model_executables = character(),
                                       operation_map = NULL,
                                       snapshot_reader = capture_v021_linux_process_snapshot) {
  samples <- .v021_resource_int(samples, "samples")
  if (!is.numeric(interval_seconds) || length(interval_seconds) != 1L ||
      is.na(interval_seconds) || !is.finite(interval_seconds) || interval_seconds < 0)
    stop("interval_seconds must be non-negative.")
  snapshots <- vector("list", samples)
  for (i in seq_len(samples)) {
    snapshots[[i]] <- snapshot_reader(
      root_pid = root_pid, known_model_executables = known_model_executables)
    if (i < samples && interval_seconds > 0) Sys.sleep(interval_seconds)
  }
  monitor_v021_process_snapshots(
    snapshots, policy, root_pid, known_model_executables, operation_map)
}

.v021_preflight_result <- function(state, reasons, policy = NULL, derivation = NULL,
                                    monitor = NULL, thread_environment = NULL) {
  if (!state %in% v021_preflight_states) stop("Invalid preflight state.")
  list(
    preflight_schema = "v021_resource_preflight_v1", state = state,
    passed = identical(state, "passed"), reasons = as.character(reasons),
    policy = policy, derivation = derivation, monitor = monitor,
    thread_environment = thread_environment
  )
}

evaluate_v021_resource_preflight <- function(policy, operation_spec,
                                               thread_environment, monitor) {
  policy_error <- tryCatch({
    validate_v021_resource_policy(policy)
    NULL
  }, error = identity)
  if (!is.null(policy_error)) return(.v021_preflight_result(
    "failed_invalid_policy", conditionMessage(policy_error), policy = policy))
  derivation <- tryCatch(derive_safe_outer_concurrency(policy, operation_spec),
                         error = identity)
  if (inherits(derivation, "error")) return(.v021_preflight_result(
    "failed_invalid_policy", conditionMessage(derivation), policy = policy))
  thread_error <- tryCatch({
    validate_v021_single_thread_environment(thread_environment)
    NULL
  }, error = identity)
  if (!is.null(thread_error)) return(.v021_preflight_result(
    "failed_thread_environment", conditionMessage(thread_error), policy, derivation,
    thread_environment = thread_environment))
  if (!is.list(monitor) || is.null(monitor$monitoring_state))
    return(.v021_preflight_result(
      "failed_monitoring_error", "Invalid monitoring result.", policy, derivation,
      monitor, thread_environment))
  if (identical(monitor$monitoring_state, "monitoring_error"))
    return(.v021_preflight_result(
      "failed_monitoring_error", monitor$reason, policy, derivation,
      monitor, thread_environment))
  if (!identical(monitor$monitoring_state, "verified"))
    return(.v021_preflight_result(
      "failed_unverified_process_tree", monitor$reason, policy, derivation,
      monitor, thread_environment))
  if (!is.numeric(monitor$observed_peak_active_cmdstan_chains) ||
      length(monitor$observed_peak_active_cmdstan_chains) != 1L ||
      is.na(monitor$observed_peak_active_cmdstan_chains))
    return(.v021_preflight_result(
      "failed_monitoring_error", "Observed peak is unavailable.", policy, derivation,
      monitor, thread_environment))
  if (!identical(as.integer(monitor$configured_ceiling),
                 policy$maximum_active_cmdstan_chains))
    return(.v021_preflight_result(
      "failed_monitoring_error", "Monitor ceiling does not match policy.", policy,
      derivation, monitor, thread_environment))
  if (!is.numeric(monitor$observed_peak_active_cmdstan_processes) ||
      length(monitor$observed_peak_active_cmdstan_processes) != 1L ||
      is.na(monitor$observed_peak_active_cmdstan_processes))
    return(.v021_preflight_result(
      "failed_monitoring_error", "Observed CmdStan process peak is unavailable.",
      policy, derivation, monitor, thread_environment))
  if (monitor$observed_peak_active_cmdstan_chains >
      policy$maximum_active_cmdstan_chains ||
      monitor$observed_peak_active_cmdstan_processes > policy$usable_chain_slots)
    return(.v021_preflight_result(
      "failed_ceiling_exceeded",
      "Observed active CmdStan chain or process slots exceeded the ceiling.",
      policy, derivation, monitor, thread_environment))
  .v021_preflight_result("passed", character(), policy, derivation, monitor,
                         thread_environment)
}
