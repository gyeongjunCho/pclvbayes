# Benchmark-only resource and process-monitoring contract for ROADMAP V021-02.

v021_resource_policy_schema <- "v021_resource_policy_v5"
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
.v021_procfs_recapture_attempts <- 2L
.v021_procfs_recapture_delay_seconds <- 0.01
.v021_linux_non_zombie_states <- c("R", "S", "D", "T", "t", "X", "x", "K", "W", "P", "I")
.v021_worker_registry_fields <- c(
  "pid", "start_time", "expected_ppid", "worker_role", "batch_id", "task_id",
  "registration_timestamp"
)

.v021_resource_int <- function(x, name, positive = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || (positive && x < 1L) || (!positive && x < 0L))
    stop(name, " must be ", if (positive) "a positive" else "a non-negative",
         " integer.")
  as.integer(x)
}

build_v021_worker_registry <- function(pid, start_time, expected_ppid, worker_role,
                                       batch_id, task_id, registration_timestamp) {
  registry <- data.frame(
    pid = as.integer(pid), start_time = as.character(start_time),
    expected_ppid = as.integer(expected_ppid), worker_role = as.character(worker_role),
    batch_id = as.character(batch_id), task_id = as.character(task_id),
    registration_timestamp = as.character(registration_timestamp),
    stringsAsFactors = FALSE
  )
  validate_v021_worker_registry(registry)
  registry
}

validate_v021_worker_registry <- function(registry) {
  if (!is.data.frame(registry) ||
      !identical(names(registry), .v021_worker_registry_fields) || !nrow(registry))
    stop("Invalid registered outer-worker registry schema.")
  integer_fields <- c("pid", "expected_ppid")
  if (any(vapply(registry[integer_fields], function(x)
    !is.integer(x) || anyNA(x), logical(1))) || any(registry$pid < 1L) ||
      any(registry$expected_ppid < 1L))
    stop("Registered outer-worker PID and PPID identities must be positive integers.")
  text_fields <- setdiff(.v021_worker_registry_fields, integer_fields)
  if (any(vapply(registry[text_fields], function(x)
    !is.character(x) || anyNA(x) || any(!nzchar(x)), logical(1))))
    stop("Registered outer-worker identity metadata must be non-empty strings.")
  if (any(!grepl("^[0-9]+$", registry$start_time)))
    stop("Registered outer-worker start times must be procfs clock-tick identities.")
  if (anyDuplicated(registry$pid))
    stop("Registered outer-worker registry contains duplicate PID entries.")
  invisible(TRUE)
}

.v021_operation_table <- function(
    main_chains, main_parallel_chains,
    retry_chains, retry_parallel_chains,
    pathfinder_processes,
    kfold_chains, kfold_parallel_chains,
    confirmation_chains, confirmation_parallel_chains) {
  data.frame(
    operation = v021_operation_types,
    total_chains = as.integer(c(
      main_chains, retry_chains, 0L, kfold_chains, confirmation_chains
    )),
    simultaneous_chain_slots = as.integer(c(
      main_parallel_chains, retry_parallel_chains, 0L,
      kfold_parallel_chains, confirmation_parallel_chains
    )),
    cmdstan_process_slots = as.integer(c(
      main_parallel_chains, retry_parallel_chains, pathfinder_processes,
      kfold_parallel_chains, confirmation_parallel_chains
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
    logical_host_threads = 16L, reserved_host_threads = 4L,
    maximum_active_cmdstan_chains = 12L, maximum_cmdstan_process_slots = 12L,
    cpu_threads_per_active_chain = 1L,
    main_chains = 4L, main_parallel_chains = main_chains,
    retry_chains = main_chains, retry_parallel_chains = main_parallel_chains,
    pathfinder_processes = 1L, pathfinder_num_paths = 8L,
    kfold_chains = main_chains, kfold_parallel_chains = 1L,
    confirmation_chains = main_chains,
    confirmation_parallel_chains = main_parallel_chains,
    proposed_outer_concurrency = NULL,
    maximum_concurrent_kfold_fits = 1L,
    global_slot_scheduler = TRUE) {
  logical_host_threads <- .v021_resource_int(logical_host_threads, "logical_host_threads")
  reserved_host_threads <- .v021_resource_int(
    reserved_host_threads, "reserved_host_threads", positive = FALSE)
  maximum_active_cmdstan_chains <- .v021_resource_int(
    maximum_active_cmdstan_chains, "maximum_active_cmdstan_chains")
  maximum_cmdstan_process_slots <- .v021_resource_int(
    maximum_cmdstan_process_slots, "maximum_cmdstan_process_slots")
  cpu_threads_per_active_chain <- .v021_resource_int(
    cpu_threads_per_active_chain, "cpu_threads_per_active_chain")
  main_chains <- .v021_resource_int(main_chains, "main_chains")
  main_parallel_chains <- .v021_resource_int(
    main_parallel_chains, "main_parallel_chains")
  retry_chains <- .v021_resource_int(retry_chains, "retry_chains")
  retry_parallel_chains <- .v021_resource_int(
    retry_parallel_chains, "retry_parallel_chains")
  pathfinder_processes <- .v021_resource_int(
    pathfinder_processes, "pathfinder_processes", positive = FALSE)
  pathfinder_num_paths <- .v021_resource_int(pathfinder_num_paths, "pathfinder_num_paths")
  kfold_chains <- .v021_resource_int(kfold_chains, "kfold_chains")
  kfold_parallel_chains <- .v021_resource_int(
    kfold_parallel_chains, "kfold_parallel_chains")
  confirmation_chains <- .v021_resource_int(confirmation_chains, "confirmation_chains")
  confirmation_parallel_chains <- .v021_resource_int(
    confirmation_parallel_chains, "confirmation_parallel_chains")
  usable <- logical_host_threads - reserved_host_threads
  if (usable < 1L || usable != maximum_active_cmdstan_chains ||
      usable != maximum_cmdstan_process_slots)
    stop("Usable host threads must equal both global CmdStan ceilings.")
  if (cpu_threads_per_active_chain != 1L)
    stop("Exactly one CPU thread per active chain is required.")
  if (main_parallel_chains > main_chains ||
      retry_parallel_chains > retry_chains ||
      kfold_parallel_chains > kfold_chains ||
      confirmation_parallel_chains > confirmation_chains)
    stop("Parallel chains cannot exceed total chains for an operation.")
  if (!is.logical(global_slot_scheduler) || length(global_slot_scheduler) != 1L ||
      is.na(global_slot_scheduler)) stop("global_slot_scheduler must be one logical value.")
  if (!isTRUE(global_slot_scheduler))
    stop("V021-02 requires the verified global slot scheduler.")

  operation_slots <- .v021_operation_table(
    main_chains, main_parallel_chains,
    retry_chains, retry_parallel_chains,
    pathfinder_processes,
    kfold_chains, kfold_parallel_chains,
    confirmation_chains, confirmation_parallel_chains
  )
  binding_slots <- max(operation_slots$simultaneous_chain_slots,
                       operation_slots$cmdstan_process_slots)
  safe_outer <- floor(usable / binding_slots)
  if (is.null(proposed_outer_concurrency)) proposed_outer_concurrency <- safe_outer
  proposed_outer_concurrency <- .v021_resource_int(
    proposed_outer_concurrency, "proposed_outer_concurrency")
  maximum_concurrent_kfold_fits <- .v021_resource_int(
    maximum_concurrent_kfold_fits, "maximum_concurrent_kfold_fits")
  if (proposed_outer_concurrency > safe_outer)
    stop("Proposed outer concurrency exceeds the global slot ceiling.")
  kfold_row <- operation_slots[
    operation_slots$operation == "kfold_fit", , drop = FALSE]
  kfold_safe <- min(
    floor(maximum_active_cmdstan_chains /
            kfold_row$simultaneous_chain_slots[[1L]]),
    floor(maximum_cmdstan_process_slots /
            kfold_row$cmdstan_process_slots[[1L]])
  )
  if (maximum_concurrent_kfold_fits > kfold_safe)
    stop("Maximum concurrent K-fold fits exceed the global slot ceiling.")

  policy <- list(
    policy_schema = v021_resource_policy_schema,
    logical_host_threads = logical_host_threads,
    reserved_host_threads = reserved_host_threads,
    usable_chain_slots = as.integer(usable),
    maximum_active_cmdstan_chains = maximum_active_cmdstan_chains,
    maximum_cmdstan_process_slots = maximum_cmdstan_process_slots,
    cpu_threads_per_active_chain = cpu_threads_per_active_chain,
    numerical_library_threads = 1L,
    main_chains = main_chains,
    main_parallel_chains = main_parallel_chains,
    retry_chains = retry_chains,
    retry_parallel_chains = retry_parallel_chains,
    pathfinder_processes = pathfinder_processes,
    pathfinder_chain_slots = 0L,
    pathfinder_num_paths = pathfinder_num_paths,
    kfold_chains = kfold_chains,
    kfold_parallel_chains = kfold_parallel_chains,
    confirmation_chains = confirmation_chains,
    confirmation_parallel_chains = confirmation_parallel_chains,
    proposed_outer_concurrency = proposed_outer_concurrency,
    maximum_concurrent_kfold_fits = maximum_concurrent_kfold_fits,
    controller_worker_limit = max(
      proposed_outer_concurrency, maximum_concurrent_kfold_fits),
    retries_share_chain_budget = TRUE,
    environment_thread_caps = v021_single_thread_environment(),
    result_root_ownership_policy = v021_controller_ownership_schema,
    scheduler_poll_interval_seconds = 0.25,
    global_slot_scheduler = TRUE,
    operation_slots = operation_slots
  )
  validate_v021_resource_policy(policy)
  policy
}

validate_v021_resource_policy <- function(policy) {
  fields <- c(
    "policy_schema", "logical_host_threads", "reserved_host_threads",
    "usable_chain_slots", "maximum_active_cmdstan_chains",
    "maximum_cmdstan_process_slots",
    "cpu_threads_per_active_chain", "numerical_library_threads", "main_chains",
    "main_parallel_chains", "retry_chains", "retry_parallel_chains",
    "pathfinder_processes", "pathfinder_chain_slots",
    "pathfinder_num_paths", "kfold_chains", "kfold_parallel_chains",
    "confirmation_chains", "confirmation_parallel_chains",
    "proposed_outer_concurrency",
    "maximum_concurrent_kfold_fits", "controller_worker_limit",
    "retries_share_chain_budget", "environment_thread_caps",
    "result_root_ownership_policy", "scheduler_poll_interval_seconds",
    "global_slot_scheduler",
    "operation_slots"
  )
  if (!is.list(policy) || !identical(names(policy), fields))
    stop("Invalid V021 resource-policy schema.")
  if (!identical(policy$policy_schema, v021_resource_policy_schema))
    stop("Invalid V021 resource-policy version.")
  for (nm in setdiff(fields, c("policy_schema", "global_slot_scheduler", "operation_slots",
                               "retries_share_chain_budget", "environment_thread_caps",
                               "result_root_ownership_policy",
                               "scheduler_poll_interval_seconds")))
    .v021_resource_int(policy[[nm]], nm, positive = nm != "pathfinder_chain_slots" &&
                         nm != "pathfinder_processes" && nm != "reserved_host_threads")
  if (policy$logical_host_threads != 16L || policy$reserved_host_threads != 4L ||
      policy$usable_chain_slots != 12L || policy$maximum_active_cmdstan_chains != 12L ||
      policy$maximum_cmdstan_process_slots != 12L)
    stop("Policy violates the exact 16/4/12 host contract.")
  if (policy$logical_host_threads - policy$reserved_host_threads !=
      policy$usable_chain_slots)
    stop("Policy host-thread arithmetic is contradictory.")
  if (policy$cpu_threads_per_active_chain != 1L ||
      policy$numerical_library_threads != 1L)
    stop("Policy requires exactly one thread per chain and numerical library.")
  if (!identical(policy$global_slot_scheduler, TRUE))
    stop("V021-02 requires the verified global slot scheduler.")
  if (!identical(
        policy$controller_worker_limit,
        max(policy$proposed_outer_concurrency,
            policy$maximum_concurrent_kfold_fits)) ||
      !identical(policy$retries_share_chain_budget, TRUE) ||
      !identical(policy$result_root_ownership_policy, v021_controller_ownership_schema) ||
      !is.numeric(policy$scheduler_poll_interval_seconds) ||
      length(policy$scheduler_poll_interval_seconds) != 1L ||
      !is.finite(policy$scheduler_poll_interval_seconds) ||
      policy$scheduler_poll_interval_seconds <= 0)
    stop("Policy scheduler, retry, or ownership settings are invalid.")
  validate_v021_single_thread_environment(policy$environment_thread_caps)
  expected_operations <- .v021_operation_table(
    policy$main_chains, policy$main_parallel_chains,
    policy$retry_chains, policy$retry_parallel_chains,
    policy$pathfinder_processes,
    policy$kfold_chains, policy$kfold_parallel_chains,
    policy$confirmation_chains, policy$confirmation_parallel_chains
  )
  if (!identical(policy$operation_slots, expected_operations))
    stop("Operation-slot table contradicts the policy.")
  if (policy$pathfinder_chain_slots != 0L)
    stop("Inspected Pathfinder behavior consumes no MCMC chain slots.")
  if (policy$pathfinder_processes != 1L || policy$pathfinder_num_paths != 8L)
    stop("Policy contradicts inspected canonical Pathfinder behavior.")
  if (policy$retry_chains != policy$main_chains ||
      policy$retry_parallel_chains != policy$main_parallel_chains ||
      policy$kfold_chains != policy$main_chains ||
      policy$kfold_parallel_chains != 1L ||
      policy$confirmation_chains != policy$main_chains ||
      policy$confirmation_parallel_chains != policy$main_parallel_chains)
    stop("Policy contradicts canonical main, retry, K-fold, or confirmation chain behavior.")
  if (policy$proposed_outer_concurrency *
      max(policy$operation_slots$simultaneous_chain_slots,
          policy$operation_slots$cmdstan_process_slots) > policy$usable_chain_slots)
    stop("Proposed outer concurrency is mathematically unsafe.")
  kfold_row <- policy$operation_slots[
    policy$operation_slots$operation == "kfold_fit", , drop = FALSE]
  if (policy$maximum_concurrent_kfold_fits *
        kfold_row$simultaneous_chain_slots[[1L]] >
        policy$maximum_active_cmdstan_chains ||
      policy$maximum_concurrent_kfold_fits *
        kfold_row$cmdstan_process_slots[[1L]] >
        policy$maximum_cmdstan_process_slots)
    stop("Maximum concurrent K-fold fits are mathematically unsafe.")
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
    floor(policy$maximum_cmdstan_process_slots / per_fit_process_slots) else Inf
  safe <- as.integer(min(chain_capacity, process_capacity))
  # An operation specification describes operations that can occur sequentially
  # within one direction, not stages that overlap.  Therefore a mixed
  # main/retry/K-fold specification is bound by the main outer-concurrency
  # ceiling; the K-fold-specific ceiling binds only a K-fold-only stage.
  has_primary_stage <- any(
    spec$operations %in% c("main_fit", "retry_fit", "confirmation_fit")
  )
  if (has_primary_stage) {
    safe <- min(safe, policy$proposed_outer_concurrency)
  } else if ("kfold_fit" %in% spec$operations) {
    safe <- min(safe, policy$maximum_concurrent_kfold_fits)
  }
  if (safe < 1L) stop("Operation cannot fit within the global resource ceiling.")
  requested <- spec$requested_outer_concurrency
  projected_chains <- requested * per_fit_chain_slots
  projected_processes <- requested * per_fit_process_slots
  if (requested > safe || projected_chains > policy$maximum_active_cmdstan_chains ||
      projected_processes > policy$maximum_cmdstan_process_slots)
    stop("Requested outer concurrency exceeds the binding global ceiling.")
  list(
    policy_schema = policy$policy_schema,
    operations = spec$operations,
    per_fit_simultaneous_chain_slots = as.integer(per_fit_chain_slots),
    per_fit_cmdstan_process_slots = as.integer(per_fit_process_slots),
    maximum_active_cmdstan_chains = policy$maximum_active_cmdstan_chains,
    usable_process_slots = policy$maximum_cmdstan_process_slots,
    safe_outer_concurrency = safe,
    requested_outer_concurrency = requested,
    projected_active_cmdstan_chains = as.integer(projected_chains),
    projected_active_cmdstan_processes = as.integer(projected_processes),
    global_slot_scheduler = TRUE,
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
      active_processes > policy$maximum_cmdstan_process_slots)
    stop("Overlapping operation plan exceeds the global ceiling.")
  list(active_cmdstan_chains = as.integer(active_chains),
       active_cmdstan_processes = as.integer(active_processes), compliant = TRUE)
}

v021_controller_ownership_schema <- "v021_controller_ownership_v1"
v021_reservation_schema <- "v021_chain_reservations_v1"

v021_resource_policy_hash <- function(policy) {
  validate_v021_resource_policy(policy)
  if (!exists("v021_sha256", mode = "function"))
    stop("V021-03 hashing must be loaded before policy hashing.")
  v021_sha256(policy)
}

v021_job_demand <- function(policy, job_type) {
  validate_v021_resource_policy(policy)
  if (!is.character(job_type) || length(job_type) != 1L ||
      !job_type %in% policy$operation_slots$operation)
    stop("Unknown V021 job type.")
  row <- policy$operation_slots[match(job_type, policy$operation_slots$operation), ]
  as.integer(max(row$simultaneous_chain_slots, row$cmdstan_process_slots))
}

new_v021_reservations <- function() {
  data.frame(
    reservation_schema = character(), task_identity = character(),
    attempt_number = integer(), job_type = character(), reserved_chain_slots = integer(),
    worker_identity = character(), acquired_at = character(), state = character(),
    stringsAsFactors = FALSE)
}

validate_v021_reservations <- function(reservations, policy) {
  expected <- names(new_v021_reservations())
  if (!is.data.frame(reservations) || !identical(names(reservations), expected))
    stop("Invalid V021 reservation schema.")
  if (!nrow(reservations)) return(invisible(TRUE))
  if (anyNA(reservations) || any(reservations$reservation_schema != v021_reservation_schema) ||
      any(!reservations$job_type %in% policy$operation_slots$operation) ||
      any(!reservations$state %in% c("reserved", "worker_terminated")) ||
      any(reservations$attempt_number < 1L) || any(reservations$reserved_chain_slots < 1L) ||
      anyDuplicated(paste(reservations$task_identity, reservations$attempt_number,
                          reservations$job_type, sep = "\r")))
    stop("Invalid V021 reservation record.")
  active <- reservations$state == "reserved"
  if (sum(reservations$reserved_chain_slots[active]) > policy$maximum_active_cmdstan_chains)
    stop("Active reservations exceed the CmdStan chain ceiling.")
  invisible(TRUE)
}

write_v021_reservations_atomic <- function(reservations, path, policy) {
  validate_v021_reservations(reservations, policy)
  validator <- function(x) validate_v021_reservations(x, policy)
  if (!exists(".v021_atomic_write_validated_rds", mode = "function"))
    stop("V021-03 atomic persistence must be loaded before reservations are written.")
  .v021_atomic_write_validated_rds(reservations, path, validator)
}

read_v021_reservations <- function(path, policy) {
  value <- tryCatch(readRDS(path), error = identity)
  if (inherits(value, "error")) stop("Reservation state is malformed or unreadable.")
  validate_v021_reservations(value, policy)
  value
}

reserve_v021_capacity <- function(reservations, policy, task_identity, attempt_number,
                                  job_type, worker_identity = "not_started",
                                  acquired_at = format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  validate_v021_reservations(reservations, policy)
  demand <- v021_job_demand(policy, job_type)
  used <- sum(reservations$reserved_chain_slots[reservations$state == "reserved"])
  if (used + demand > policy$maximum_active_cmdstan_chains)
    stop("Insufficient V021 CmdStan chain capacity.")
  row <- data.frame(
    reservation_schema = v021_reservation_schema,
    task_identity = as.character(task_identity), attempt_number = as.integer(attempt_number),
    job_type = job_type, reserved_chain_slots = demand,
    worker_identity = as.character(worker_identity), acquired_at = as.character(acquired_at),
    state = "reserved", stringsAsFactors = FALSE)
  out <- rbind(reservations, row)
  validate_v021_reservations(out, policy)
  out
}

release_v021_capacity <- function(reservations, policy, task_identity, attempt_number,
                                  job_type, worker_terminated) {
  validate_v021_reservations(reservations, policy)
  if (!isTRUE(worker_terminated))
    stop("Capacity cannot be released before worker termination is established.")
  hit <- which(reservations$task_identity == task_identity &
                 reservations$attempt_number == as.integer(attempt_number) &
                 reservations$job_type == job_type & reservations$state == "reserved")
  if (length(hit) != 1L) stop("Reservation identity is missing or ambiguous.")
  reservations$state[[hit]] <- "worker_terminated"
  validate_v021_reservations(reservations, policy)
  reservations
}

plan_v021_scheduler_waves <- function(task_identities, policy,
                                       job_type = "main_fit") {
  validate_v021_resource_policy(policy)
  task_identities <- sort(unique(as.character(task_identities)), method = "radix")
  if (!length(task_identities) || anyNA(task_identities) || any(!nzchar(task_identities)))
    stop("Runnable task identities must be non-empty and unique.")
  demand <- v021_job_demand(policy, job_type)
  per_wave <- floor(policy$maximum_active_cmdstan_chains / demand)
  if (per_wave < 1L) stop("One job exceeds the total CmdStan chain budget.")
  groups <- split(task_identities, ceiling(seq_along(task_identities) / per_wave))
  lapply(seq_along(groups), function(i) list(
    wave = as.integer(i), task_identities = unname(groups[[i]]),
    reserved_chain_slots = as.integer(length(groups[[i]]) * demand)))
}

.v021_local_process_identity <- function(pid = Sys.getpid()) {
  stat_path <- file.path("/proc", as.character(pid), "stat")
  stat <- if (file.exists(stat_path)) tryCatch(readLines(stat_path, warn = FALSE), error = identity)
  else character()
  parsed <- if (length(stat) == 1L && !inherits(stat, "error"))
    tryCatch(.v021_parse_linux_stat(stat), error = function(e) NULL) else NULL
  boot <- tryCatch(readLines("/proc/sys/kernel/random/boot_id", warn = FALSE),
                   error = function(e) NA_character_)
  list(pid = as.integer(pid),
       start_time = if (is.null(parsed)) NA_character_ else parsed$start_time,
       hostname = unname(Sys.info()[["nodename"]]),
       boot_id = if (length(boot)) boot[[1L]] else NA_character_)
}

v021_ownership_path <- function(result_root) file.path(result_root, ".v021-controller-lock")

acquire_v021_controller_ownership <- function(result_root, manifest_hash,
                                               configuration_hash,
                                               identity = .v021_local_process_identity()) {
  lock <- v021_ownership_path(result_root)
  if (!dir.exists(result_root) && !dir.create(result_root, recursive = TRUE))
    stop("Could not create result root for ownership.")
  if (!dir.create(lock, showWarnings = FALSE))
    stop("V021 result root is already owned; reconcile explicitly.")
  record <- list(
    ownership_schema = v021_controller_ownership_schema,
    result_root = normalizePath(result_root, mustWork = TRUE),
    manifest_hash = manifest_hash, configuration_hash = configuration_hash,
    hostname = identity$hostname, controller_pid = as.integer(identity$pid),
    process_start = identity$start_time, boot_id = identity$boot_id,
    acquired_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  tryCatch({
    temporary <- file.path(lock, ".owner.rds.tmp")
    saveRDS(record, temporary)
    if (!identical(readRDS(temporary), record) ||
        !file.rename(temporary, file.path(lock, "owner.rds")))
      stop("Atomic ownership metadata promotion failed.")
  }, error = function(e) {
    unlink(lock, recursive = TRUE); stop(e)
  })
  record
}

validate_v021_controller_ownership <- function(record, result_root,
                                                manifest_hash, configuration_hash) {
  required <- c("ownership_schema", "result_root", "manifest_hash", "configuration_hash",
                "hostname", "controller_pid", "process_start", "boot_id", "acquired_at")
  if (!is.list(record) || !identical(names(record), required) ||
      !identical(record$ownership_schema, v021_controller_ownership_schema) ||
      !identical(record$result_root, normalizePath(result_root, mustWork = TRUE)) ||
      !identical(record$manifest_hash, manifest_hash) ||
      !identical(record$configuration_hash, configuration_hash))
    stop("Controller ownership identity or provenance mismatch.")
  stored <- tryCatch(suppressWarnings(readRDS(
    file.path(v021_ownership_path(result_root), "owner.rds"))),
                     error = identity)
  if (inherits(stored, "error") || !identical(stored, record))
    stop("Controller does not hold the active result-root ownership record.")
  invisible(TRUE)
}

release_v021_controller_ownership <- function(result_root, record) {
  lock <- v021_ownership_path(result_root)
  stored <- tryCatch(readRDS(file.path(lock, "owner.rds")), error = identity)
  if (inherits(stored, "error") || !identical(stored, record))
    stop("Refusing to release ownership not held by this controller.")
  unlink(lock, recursive = TRUE)
  invisible(TRUE)
}

reconcile_v021_stale_ownership <- function(result_root, current_identity,
                                            administrative_override = FALSE) {
  lock <- v021_ownership_path(result_root)
  stored <- tryCatch(readRDS(file.path(lock, "owner.rds")), error = identity)
  if (inherits(stored, "error")) stop("Ownership record is malformed.")
  if (!identical(stored$hostname, current_identity$hostname) ||
      !identical(stored$boot_id, current_identity$boot_id))
    stop("Foreign-host or prior-boot ownership requires administrative reconciliation.")
  live <- dir.exists(file.path("/proc", as.character(stored$controller_pid))) &&
    identical(.v021_local_process_identity(stored$controller_pid)$start_time,
              stored$process_start)
  if (live) stop("Controller ownership is still live.")
  if (!isTRUE(administrative_override))
    stop("Stale ownership requires explicit administrative_override.")
  unlink(lock, recursive = TRUE)
  invisible(stored)
}

audit_v021_resource_state <- function(policy, reservations, ownership_valid,
                                       worker_identities = character(),
                                       actual_cmdstan_chains = NA_integer_) {
  validate_v021_reservations(reservations, policy)
  reserved <- sum(reservations$reserved_chain_slots[reservations$state == "reserved"])
  if (!is.na(actual_cmdstan_chains) && actual_cmdstan_chains > reserved)
    stop("Observed CmdStan descendants exceed registered reservations.")
  list(controller_pid = as.integer(Sys.getpid()), worker_identities = worker_identities,
       total_reserved_chain_slots = as.integer(reserved),
       maximum_allowed_chain_slots = policy$maximum_active_cmdstan_chains,
       configured_thread_limits = v021_single_thread_environment(),
       oversubscribed = reserved > policy$maximum_active_cmdstan_chains,
       ownership_valid = isTRUE(ownership_valid))
}

claim_v021_task_with_capacity <- function(status, ledger, manifest, task_ordinal,
                                           ownership_record, result_root, policy,
                                           reservations, timestamp,
                                           job_type = "main_fit") {
  validate_v021_controller_ownership(
    ownership_record, result_root, manifest$manifest_hash,
    manifest$configuration_hash)
  i <- match(as.integer(task_ordinal), manifest$tasks$task_ordinal)
  if (is.na(i)) stop("Unknown task ordinal.")
  attempt_number <- status$tasks$attempt_count[[i]] + 1L
  reserved <- reserve_v021_capacity(
    reservations, policy, manifest$tasks$directed_task_id[[i]],
    attempt_number, job_type)
  started <- tryCatch(start_v021_task_attempt(
    status, ledger, manifest, task_ordinal, timestamp,
    worker_provenance = "controller_reserved"), error = identity)
  if (inherits(started, "error")) stop(conditionMessage(started))
  list(status = started$status, ledger = started$ledger,
       reservations = reserved, attempt_id = started$attempt_id)
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

.v021_snapshot_argv <- function(snapshot, allow_empty = rep(FALSE, nrow(snapshot))) {
  if (!is.logical(allow_empty) || length(allow_empty) != nrow(snapshot) || anyNA(allow_empty))
    stop("Invalid process argv exception mask.")
  if ("argv" %in% names(snapshot)) {
    terminal <- "capture_state" %in% names(snapshot) &
      snapshot$capture_state %in% c("vanished_during_capture", "zombie_process")
    valid <- vapply(snapshot$argv, function(x)
      is.character(x) && length(x) > 0L && !anyNA(x) && nzchar(x[[1L]]), logical(1))
    if (!all(valid | terminal | allow_empty)) stop("Malformed decoded process argv.")
    return(snapshot$argv)
  }
  lapply(snapshot$command, function(x) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) character()
    else strsplit(x, "[[:space:]]+", perl = TRUE)[[1L]]
  })
}

.v021_process_ancestry <- function(records, pid, root_pid) {
  path <- as.integer(pid)
  seen <- integer()
  current <- as.integer(pid)
  repeat {
    if (current %in% seen) break
    seen <- c(seen, current)
    if (identical(current, as.integer(root_pid))) break
    row <- match(current, records$pid)
    if (is.na(row)) break
    parent <- as.integer(records$ppid[[row]])
    if (is.na(parent) || !parent %in% records$pid) break
    path <- c(parent, path)
    current <- parent
  }
  path
}

.v021_empty_offending_processes <- function() {
  data.frame(
    pid = integer(), ppid = integer(), start_time = character(),
    argv = I(list()), command = character(), executable = character(),
    classification = character(), capture_state = character(),
    ancestry = I(list()), classifier_stage = character(), reason = character(),
    stringsAsFactors = FALSE)
}

.v021_offending_processes <- function(classified, reason) {
  if (!is.list(classified) || !length(classified))
    return(.v021_empty_offending_processes())
  rows <- list()
  for (records in classified) {
    hit <- which(records$classification == "unknown_potential_cmdstan")
    if (!length(hit)) next
    for (i in hit) {
      rows[[length(rows) + 1L]] <- data.frame(
        pid = as.integer(records$pid[[i]]),
        ppid = as.integer(records$ppid[[i]]),
        start_time = if ("start_time" %in% names(records))
          as.character(records$start_time[[i]]) else NA_character_,
        argv = I(list(records$argv[[i]])),
        command = as.character(records$command[[i]]),
        executable = as.character(records$executable[[i]]),
        classification = as.character(records$classification[[i]]),
        capture_state = if ("capture_state" %in% names(records))
          as.character(records$capture_state[[i]]) else "captured",
        ancestry = I(list(.v021_process_ancestry(
          records, records$pid[[i]], records$pid[records$classification == "parent_r"][[1L]]))),
        classifier_stage = "operation_classification",
        reason = as.character(reason), stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) .v021_empty_offending_processes() else do.call(rbind, rows)
}

classify_v021_process_snapshot <- function(snapshot, root_pid,
                                            known_model_executables = character(),
                                            operation_map = NULL,
                                            worker_registry = NULL) {
  required <- c("timestamp", "pid", "ppid", "command", "executable",
                "potential_cmdstan", "readable")
  allowed <- append(required, "argv", after = 4L)
  extended_v1 <- c("timestamp", "discovery_time", "pid", "ppid", "start_time",
                "command", "argv", "executable", "potential_cmdstan", "readable",
                "capture_state", "disappearance_reason")
  extended <- c("timestamp", "discovery_time", "pid", "ppid", "start_time",
                "process_state", "command", "argv", "executable",
                "potential_cmdstan", "readable", "capture_state",
                "disappearance_reason", "zombie_reason")
  extended_retry <- c("timestamp", "discovery_time", "pid", "ppid", "start_time",
                      "process_state", "initial_process_state", "final_process_state",
                      "command", "argv", "executable", "potential_cmdstan", "readable",
                      "capture_state", "disappearance_reason", "zombie_reason",
                      "capture_retry_count", "capture_retry_timestamps",
                      "resolution_reason")
  if (!is.data.frame(snapshot) ||
      !(identical(names(snapshot), required) || identical(names(snapshot), allowed) ||
        identical(names(snapshot), extended_v1) || identical(names(snapshot), extended) ||
        identical(names(snapshot), extended_retry)) ||
      anyDuplicated(snapshot$pid) || !root_pid %in% snapshot$pid)
    stop("Invalid process snapshot.")
  if (identical(names(snapshot), extended_v1) || identical(names(snapshot), extended) ||
      identical(names(snapshot), extended_retry)) {
    states <- c("captured", "vanished_during_capture", "unreadable_live",
                "ambiguous_identity", "zombie_process")
    if (anyNA(snapshot$capture_state) || !all(snapshot$capture_state %in% states) ||
        any(snapshot$readable != (snapshot$capture_state == "captured")))
      stop("Invalid process capture state.")
    vanished_rows <- snapshot$capture_state == "vanished_during_capture"
    valid_vanished <- !snapshot$readable & !nzchar(snapshot$command) &
      !nzchar(snapshot$executable) &
      vapply(snapshot$argv, function(x) is.character(x) && !length(x), logical(1)) &
      snapshot$disappearance_reason %in% c("pid_disappeared", "pid_reused")
    if (any(vanished_rows & !valid_vanished))
      stop("Invalid vanished process record.")
    if (identical(names(snapshot), extended) || identical(names(snapshot), extended_retry)) {
      zombie_rows <- snapshot$capture_state == "zombie_process"
      valid_zombie <- !snapshot$readable & snapshot$process_state == "Z" &
        !nzchar(snapshot$command) & !nzchar(snapshot$executable) &
        vapply(snapshot$argv, function(x) is.character(x) && !length(x), logical(1)) &
        snapshot$zombie_reason == "exited_unreaped"
      if (any(zombie_rows & !valid_zombie)) stop("Invalid zombie process record.")
    }
    if (identical(names(snapshot), extended_retry)) {
      retry_valid <- is.integer(snapshot$capture_retry_count) &
        snapshot$capture_retry_count >= 0L &
        vapply(seq_len(nrow(snapshot)), function(i) {
          timestamps <- snapshot$capture_retry_timestamps[[i]]
          is.character(timestamps) && length(timestamps) == snapshot$capture_retry_count[[i]] &&
            !anyNA(timestamps) && all(nzchar(timestamps))
        }, logical(1)) & !is.na(snapshot$initial_process_state) &
        !is.na(snapshot$final_process_state) & !is.na(snapshot$resolution_reason)
      if (!all(retry_valid)) stop("Invalid procfs recapture audit metadata.")
    }
  }
  descendants <- .v021_descendant_pids(snapshot, root_pid)
  vanished <- if ("capture_state" %in% names(snapshot))
    snapshot$capture_state == "vanished_during_capture" else rep(FALSE, nrow(snapshot))
  zombie <- if ("capture_state" %in% names(snapshot))
    snapshot$capture_state == "zombie_process" else rep(FALSE, nrow(snapshot))
  out <- snapshot[snapshot$pid %in% descendants | vanished | zombie, , drop = FALSE]
  out$is_descendant <- out$pid != root_pid
  registered <- rep(FALSE, nrow(out))
  registry_match <- rep(NA_integer_, nrow(out))
  if (!is.null(worker_registry)) {
    validate_v021_worker_registry(worker_registry)
    if (!all(c("start_time", "process_state", "capture_state") %in% names(out)))
      stop("Registered outer-worker classification requires procfs identity data.")
    registry_match <- match(out$pid, worker_registry$pid)
    terminal_capture <- out$capture_state %in%
      c("vanished_during_capture", "zombie_process")
    candidate <- !is.na(registry_match) & !terminal_capture
    if (any(candidate)) {
      registry_rows <- worker_registry[registry_match[candidate], , drop = FALSE]
      observed <- out[candidate, , drop = FALSE]
      valid_state <- !is.na(observed$process_state) &
        observed$process_state %in% .v021_linux_non_zombie_states
      valid_identity <- observed$start_time == registry_rows$start_time &
        observed$ppid == registry_rows$expected_ppid & valid_state &
        observed$ppid %in% descendants & observed$is_descendant &
        observed$capture_state %in% c("captured", "unreadable_live") &
        !observed$potential_cmdstan
      if (anyNA(valid_identity) || !all(valid_identity))
        stop("Registered outer-worker PID, start-time, PPID, or ancestry mismatch.")
      registered[candidate] <- TRUE
    }
  }
  registered_unreadable <- registered
  if ("capture_state" %in% names(out))
    registered_unreadable <- registered & out$capture_state == "unreadable_live"
  argv <- .v021_snapshot_argv(out, allow_empty = registered_unreadable)
  if (!"argv" %in% names(out)) out$argv <- I(argv)
  sample_argument <- vapply(argv, function(x) "method=sample" %in% x, logical(1))
  pathfinder_argument <- vapply(argv, function(x) "method=pathfinder" %in% x, logical(1))
  known_executable <- out$executable %in% known_model_executables
  sample_process <- known_executable & sample_argument
  pathfinder_process <- known_executable & pathfinder_argument
  r_process <- grepl("(^|/)(R|Rscript)([[:space:]]|$)", out$command)
  parent_index <- match(out$ppid, out$pid)
  valid_parent <- !is.na(parent_index)
  parent_is_r_worker <- parent_is_registered_worker <- rep(FALSE, nrow(out))
  parent_is_r_worker[valid_parent] <- r_process[parent_index[valid_parent]] &
    out$pid[parent_index[valid_parent]] != root_pid
  parent_is_registered_worker[valid_parent] <- registered[parent_index[valid_parent]]
  canonical_diagnose_executable <- grepl(
    "/cmdstan-[^/]+/bin/diagnose$", out$executable)
  canonical_diagnose_argv <- vapply(argv, function(x) {
    length(x) >= 2L && identical(x[[1L]], "bin/diagnose") &&
      all(grepl("\\.csv$", x[-1L]))
  }, logical(1))
  registered_diagnostic_identity <- rep(FALSE, nrow(out))
  if (all(c("start_time", "process_state", "capture_state") %in% names(out)))
    registered_diagnostic_identity <- out$readable & out$capture_state == "captured" &
      !is.na(out$start_time) & grepl("^[0-9]+$", out$start_time) &
      out$process_state %in% .v021_linux_non_zombie_states
  parent_is_controller <- valid_parent &
    out$pid[parent_index] == root_pid & r_process[parent_index]
  diagnostic_process <- out$is_descendant & canonical_diagnose_executable &
    canonical_diagnose_argv & (parent_is_controller | parent_is_r_worker |
      (parent_is_registered_worker & registered_diagnostic_identity))
  out_vanished <- if ("capture_state" %in% names(out))
    out$capture_state == "vanished_during_capture" else rep(FALSE, nrow(out))
  out$classification <- "other_descendant"
  out$classification[out_vanished] <- "vanished_during_capture"
  if ("capture_state" %in% names(out))
    out$classification[out$capture_state == "zombie_process"] <- "zombie_process"
  out$classification[out$pid == root_pid] <- "parent_r"
  out$classification[out$is_descendant & r_process] <- "outer_worker"
  out$classification[registered] <- "registered_outer_worker"
  out$classification[out$is_descendant & pathfinder_process] <- "pathfinder_process"
  out$classification[out$is_descendant & sample_process] <- "cmdstan_chain"
  out$classification[diagnostic_process] <- "cmdstan_diagnostic"
  unknown <- out$is_descendant & (out$potential_cmdstan | known_executable) &
    !sample_process & !pathfinder_process & !diagnostic_process
  out$classification[unknown] <- "unknown_potential_cmdstan"
  out$is_active_cmdstan_chain <- out$classification == "cmdstan_chain"
  out$operation <- NA_character_
  out$task_id <- NA_character_
  out$direction_id <- NA_character_
  out$worker_role <- NA_character_
  out$worker_batch_id <- NA_character_
  out$worker_registration_timestamp <- NA_character_
  if (any(registered)) {
    rows <- worker_registry[registry_match[registered], , drop = FALSE]
    out$worker_role[registered] <- rows$worker_role
    out$worker_batch_id[registered] <- rows$batch_id
    out$task_id[registered] <- rows$task_id
    out$worker_registration_timestamp[registered] <- rows$registration_timestamp
  }
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
                                            operation_map = NULL,
                                            worker_registry = NULL) {
  validate_v021_resource_policy(policy)
  if (!is.list(snapshots) || !length(snapshots))
    return(list(monitoring_state = "monitoring_error", reason = "no_process_snapshots",
                records = data.frame(),
                offending_processes = .v021_empty_offending_processes(),
                observed_peak_active_cmdstan_chains = NA_integer_,
                observed_peak_active_cmdstan_processes = NA_integer_,
                configured_ceiling = policy$maximum_active_cmdstan_chains,
                compliance_status = "unverified"))
  classified <- tryCatch(lapply(snapshots, classify_v021_process_snapshot,
    root_pid = root_pid, known_model_executables = known_model_executables,
    operation_map = operation_map, worker_registry = worker_registry), error = identity)
  if (inherits(classified, "error"))
    return(list(monitoring_state = "monitoring_error",
                reason = conditionMessage(classified), records = data.frame(),
                offending_processes = .v021_empty_offending_processes(),
                observed_peak_active_cmdstan_chains = NA_integer_,
                observed_peak_active_cmdstan_processes = NA_integer_,
                configured_ceiling = policy$maximum_active_cmdstan_chains,
                compliance_status = "unverified"))
  records <- do.call(rbind, classified)
  counts <- vapply(classified, function(x) sum(x$is_active_cmdstan_chain), integer(1))
  process_counts <- vapply(classified, function(x)
    sum(x$classification %in%
          c("cmdstan_chain", "pathfinder_process", "cmdstan_diagnostic")), integer(1))
  peak <- max(counts)
  process_peak <- max(process_counts)
  main_slot_row <- policy$operation_slots[
    policy$operation_slots$operation == "main_fit", , drop = FALSE]
  kfold_slot_row <- policy$operation_slots[
    policy$operation_slots$operation == "kfold_fit", , drop = FALSE]
  configured_launch_chain_ceiling <- min(
    policy$maximum_active_cmdstan_chains,
    max(
      policy$proposed_outer_concurrency *
        main_slot_row$simultaneous_chain_slots[[1L]],
      policy$maximum_concurrent_kfold_fits *
        kfold_slot_row$simultaneous_chain_slots[[1L]]
    ))
  complete <- all(vapply(classified, function(x) {
    terminal <- if ("capture_state" %in% names(x))
      x$capture_state %in% c("vanished_during_capture", "zombie_process")
    else rep(FALSE, nrow(x))
    all(x$readable | terminal | x$classification == "registered_outer_worker")
  }, logical(1)))
  unknown <- any(records$classification == "unknown_potential_cmdstan")
  state <- if (!complete || unknown) "unverified_process_tree" else "verified"
  reason <- if (!complete) "unreadable_descendant_process" else if (unknown)
    "unknown_potential_cmdstan_descendant" else NA_character_
  offending_processes <- .v021_offending_processes(classified, reason)
  records$active_cmdstan_chain_count <- rep(counts, vapply(classified, nrow, integer(1)))
  records$active_cmdstan_process_count <- rep(
    process_counts, vapply(classified, nrow, integer(1)))
  records$observed_peak_active_cmdstan_chain_count <- peak
  records$observed_peak_active_cmdstan_process_count <- process_peak
  records$configured_ceiling <- policy$maximum_active_cmdstan_chains
  records$configured_launch_chain_ceiling <- configured_launch_chain_ceiling
  records$compliance_status <- if (state != "verified") "unverified" else if (
    peak <= policy$maximum_active_cmdstan_chains &&
      process_peak <= policy$maximum_cmdstan_process_slots) "compliant" else "exceeded"
  list(
    monitoring_state = state, reason = reason, records = records,
    offending_processes = offending_processes,
    observed_peak_active_cmdstan_chains = as.integer(peak),
    observed_peak_active_cmdstan_processes = as.integer(process_peak),
    configured_ceiling = policy$maximum_active_cmdstan_chains,
    configured_launch_chain_ceiling = as.integer(configured_launch_chain_ceiling),
    compliance_status = unique(records$compliance_status)[[1L]]
  )
}

.v021_parse_linux_stat <- function(stat, expected_pid) {
  if (!is.character(stat) || length(stat) != 1L || is.na(stat) || !nzchar(stat))
    stop("Linux process stat is unreadable.")
  matched <- regexec("^([0-9]+) \\((.*)\\) ([^ ]) (.*)$", stat)
  fields <- regmatches(stat, matched)[[1L]]
  if (length(fields) != 5L || !identical(as.integer(fields[[2L]]), as.integer(expected_pid)))
    stop("Linux process stat identity is malformed.")
  remainder <- strsplit(fields[[5L]], " ", fixed = TRUE)[[1L]]
  if (length(remainder) < 19L) stop("Linux process stat is incomplete.")
  ppid <- suppressWarnings(as.integer(remainder[[1L]]))
  start_time <- remainder[[19L]]
  if (is.na(ppid) || !nzchar(start_time) || !grepl("^[0-9]+$", start_time))
    stop("Linux process stat identity is malformed.")
  list(ppid = ppid, start_time = start_time, process_state = fields[[4L]])
}

.v021_recheck_linux_process <- function(pid_dir, pid, start_time,
                                         path_exists, stat_reader) {
  if (!isTRUE(path_exists(pid_dir)))
    return(list(state = "vanished_during_capture", reason = "pid_disappeared"))
  observed <- tryCatch(
    .v021_parse_linux_stat(stat_reader(file.path(pid_dir, "stat")), pid),
    error = identity)
  if (inherits(observed, "error"))
    return(list(state = "ambiguous_identity", reason = "identity_recheck_unreadable"))
  if (!identical(observed$start_time, start_time))
    return(list(state = "vanished_during_capture", reason = "pid_reused"))
  if (identical(observed$process_state, "Z"))
    return(list(state = "zombie_process", reason = "exited_unreaped",
                process_state = "Z"))
  list(state = "identity_stable", reason = "identity_stable",
       process_state = observed$process_state)
}

capture_v021_linux_process_snapshot <- function(
    root_pid = Sys.getpid(), known_model_executables = character(),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), proc_root = "/proc",
    pid_lister = function(root) list.files(root, pattern = "^[0-9]+$", full.names = FALSE),
    path_exists = dir.exists,
    stat_reader = function(path) readLines(path, warn = FALSE, n = 1L),
    cmdline_reader = .v021_read_linux_cmdline,
    executable_reader = Sys.readlink,
    sleep = Sys.sleep,
    clock = function() format(Sys.time(), tz = "UTC", usetz = TRUE)) {
  if (!dir.exists(proc_root)) stop("Linux /proc process monitoring is unavailable.")
  discovery_time <- as.character(timestamp)
  entries <- pid_lister(proc_root)
  identities <- lapply(entries, function(pid_text) {
    stat_path <- file.path(proc_root, pid_text, "stat")
    stat <- tryCatch(stat_reader(stat_path), error = function(e) character())
    if (!length(stat)) {
      if (isTRUE(path_exists(file.path(proc_root, pid_text))))
        stop("Linux process stat remained unreadable for live PID ", pid_text, ".")
      return(NULL)
    }
    parsed <- tryCatch(.v021_parse_linux_stat(stat, as.integer(pid_text)), error = identity)
    if (inherits(parsed, "error")) return(NULL)
    data.frame(pid = as.integer(pid_text), ppid = parsed$ppid,
               start_time = parsed$start_time, process_state = parsed$process_state,
               stringsAsFactors = FALSE)
  })
  inventory <- do.call(rbind, identities[!vapply(identities, is.null, logical(1))])
  if (is.null(inventory) || !root_pid %in% inventory$pid)
    stop("Root process could not be established in /proc.")
  descendants <- .v021_descendant_pids(inventory, root_pid)
  inventory <- inventory[inventory$pid %in% descendants, , drop = FALSE]
  rows <- lapply(seq_len(nrow(inventory)), function(i) {
    pid <- inventory$pid[[i]]; pid_text <- as.character(pid)
    pid_dir <- file.path(proc_root, pid_text)
    retry_timestamps <- character()
    cmd_result <- executable <- NULL
    recheck <- NULL
    for (attempt in 0:.v021_procfs_recapture_attempts) {
      if (attempt > 0L) {
        sleep(.v021_procfs_recapture_delay_seconds)
        retry_timestamps <- c(retry_timestamps, as.character(clock()))
      }
      cmd_result <- tryCatch({
        raw <- cmdline_reader(file.path(pid_dir, "cmdline"))
        .v021_decode_linux_cmdline(raw)
      }, error = identity)
      executable <- tryCatch(executable_reader(file.path(pid_dir, "exe")), error = identity)
      capture_ok <- !inherits(cmd_result, "error") && !inherits(executable, "error") &&
        is.character(executable) && length(executable) == 1L && nzchar(executable)
      recheck <- .v021_recheck_linux_process(pid_dir, pid, inventory$start_time[[i]],
                                             path_exists, stat_reader)
      if (!identical(recheck$state, "identity_stable")) break
      if (capture_ok) {
        recheck$state <- "captured"
        recheck$reason <- if (attempt) "readable_on_recapture" else "readable_initial_capture"
        break
      }
      if (attempt == .v021_procfs_recapture_attempts) {
        recheck$state <- "unreadable_live"
        recheck$reason <- "recapture_limit_exhausted"
      }
    }
    readable <- identical(recheck$state, "captured")
    argv <- if (readable) cmd_result$argv else character()
    command <- if (readable) cmd_result$command else ""
    executable <- if (readable) executable else ""
    potential <- any(grepl("^method=", argv)) || executable %in% known_model_executables ||
      grepl("cmdstan", command, ignore.case = TRUE)
    data.frame(timestamp = as.character(timestamp), discovery_time = discovery_time,
               pid = pid, ppid = inventory$ppid[[i]],
               start_time = inventory$start_time[[i]],
               process_state = if (!is.null(recheck$process_state))
                 recheck$process_state else inventory$process_state[[i]],
               initial_process_state = inventory$process_state[[i]],
               final_process_state = if (!is.null(recheck$process_state))
                 recheck$process_state else inventory$process_state[[i]],
               command = command, argv = I(list(argv)),
               executable = executable, potential_cmdstan = potential,
               readable = readable, capture_state = recheck$state,
               disappearance_reason = if (recheck$state == "vanished_during_capture")
                 recheck$reason else NA_character_,
               zombie_reason = if (recheck$state == "zombie_process") recheck$reason
                 else NA_character_, capture_retry_count = as.integer(length(retry_timestamps)),
               capture_retry_timestamps = I(list(retry_timestamps)),
               resolution_reason = recheck$reason, stringsAsFactors = FALSE)
  })
  snapshot <- do.call(rbind, rows)
  rownames(snapshot) <- NULL
  snapshot
}

monitor_v021_process_tree <- function(policy, root_pid = Sys.getpid(), samples = 1L,
                                       interval_seconds = 0.1,
                                       known_model_executables = character(),
                                       operation_map = NULL,
                                       worker_registry = NULL,
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
    snapshots, policy, root_pid, known_model_executables, operation_map,
    worker_registry)
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
      monitor$observed_peak_active_cmdstan_processes >
        policy$maximum_cmdstan_process_slots ||
      identical(monitor$compliance_status, "exceeded"))
    return(.v021_preflight_result(
      "failed_ceiling_exceeded",
      "Observed active CmdStan chain or process slots exceeded the ceiling.",
      policy, derivation, monitor, thread_environment))
  .v021_preflight_result("passed", character(), policy, derivation, monitor,
                         thread_environment)
}
