#!/usr/bin/env bash
set -euo pipefail

REPO="/home/heuklang/pclvbayes"
RESULT_ROOT="$REPO/benchmarks/mtist/results/v021_full_100_species_v5"
CONFIG="$REPO/benchmarks/mtist/configs/v021_full_100_species_v5.R"
RUNNER="$REPO/benchmarks/mtist/run_v021_full_100_species.R"
EXECUTABLE="$REPO/inst/stan/pclv"
LOG_DIR="$REPO/benchmarks/mtist/logs"
LOG="$LOG_DIR/v021_full_100_species_v5.nohup.log"
PIDFILE="$LOG_DIR/v021_full_100_species_v5.controller.pid"
LAUNCH_META="$LOG_DIR/v021_full_100_species_v5.launch.tsv"
EXPECTED_BRANCH="v0.2.1-dev"
REQUIRED_CONTROL_SCRIPT_COMMIT="32d983bd7282db1d5f4297576ffd612213b69338"
REQUIRED_FULL_RUN_COMMIT="8f1571a7cfedabfaef0389ed61944ff1cc19f97e"
EXPECTED_AFFINITY="0-15"
PROCESS_PATTERN='run_v021_full_100_species|(^|[[:space:]])(\./|/[^[:space:]]*/)?pclv([[:space:]]|$)'

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

script_commits() {
  SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
  SCRIPT_REL="${SCRIPT_PATH#"$REPO/"}"
  SCRIPT_COMMIT="$(git -C "$REPO" log -1 --format=%H -- "$SCRIPT_REL" 2>/dev/null || true)"
  CURRENT_HEAD="$(git -C "$REPO" rev-parse HEAD)"
}

validate_reviewed_history() {
  local repo="$1" current_head="$2" required_control="$3" required_integration="$4"
  git -C "$repo" merge-base --is-ancestor "$required_control" "$current_head" ||
    die "reviewed control-script commit is not an ancestor of HEAD"
  git -C "$repo" merge-base --is-ancestor "$required_integration" "$current_head" ||
    die "reviewed full-run integration commit is not an ancestor of HEAD"
}

validate_repository_state() {
  local repo="$1" expected_branch="$2"
  test "$(git -C "$repo" branch --show-current)" = "$expected_branch" ||
    die "current branch is not $expected_branch"
  test -z "$(git -C "$repo" status --short --untracked-files=no)" ||
    die "tracked working tree is not clean"
}

controller_pid() {
  test -f "$PIDFILE" || return 1
  local pid
  pid="$(tr -d '[:space:]' < "$PIDFILE" 2>/dev/null || true)"
  [[ "$pid" =~ ^[0-9]+$ ]] || return 1
  printf '%s\n' "$pid"
}

controller_alive() {
  local pid
  pid="$(controller_pid)" || return 1
  kill -0 "$pid" 2>/dev/null
}

common_start_validation() {
  test -d "$REPO/.git" || die "repository is missing: $REPO"
  cd "$REPO"
  validate_repository_state "$REPO" "$EXPECTED_BRANCH"
  test -f "$RUNNER" || die "runner is missing: $RUNNER"
  test -f "$CONFIG" || die "configuration is missing: $CONFIG"
  test -x "$EXECUTABLE" || die "compiled executable is missing: $EXECUTABLE"
  command -v taskset >/dev/null || die "taskset is unavailable"
  test "$(nproc)" -ge 16 || die "fewer than 16 logical CPUs are available"
  taskset -c "$EXPECTED_AFFINITY" true 2>/dev/null ||
    die "CPU affinity $EXPECTED_AFFINITY is unavailable"
  script_commits
  test -n "$SCRIPT_COMMIT" || die "control script is uncommitted"
  validate_reviewed_history "$REPO" "$CURRENT_HEAD" \
    "$REQUIRED_CONTROL_SCRIPT_COMMIT" "$REQUIRED_FULL_RUN_COMMIT"
  if controller_alive; then
    die "controller PID $(controller_pid) is already alive"
  fi
}

archive_file() {
  local path="$1" stamp="$2"
  test ! -e "$path" || mv "$path" "${path}.preserved.${stamp}"
}

start_controller() {
  local output_log="$1" metadata="$2" mode="$3"
  printf '%s_commit\t%s\n%s_time\t%s\ncpu_affinity\t0-15\nactive_chain_cap\t12\nconcurrent_fits\t3\n' \
    "$mode" "$CURRENT_HEAD" "$mode" "$(date --iso-8601=seconds)" > "$metadata"

  nohup env \
    MTIST_ROOT=/home/heuklang/mtist \
    PCLV_V021_FULL_CONFIG="$CONFIG" \
    STAN_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    OMP_THREAD_LIMIT=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    BLIS_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    RCPP_PARALLEL_NUM_THREADS=1 \
    taskset -c "$EXPECTED_AFFINITY" \
    Rscript "$RUNNER" \
    > "$output_log" 2>&1 < /dev/null &

  local pid=$!
  printf '%s\n' "$pid" > "$PIDFILE"
  disown "$pid" 2>/dev/null || true
  printf 'V021-06 controller PID: %s\nCommit: %s\nResult root: %s\nLog: %s\nPID file: %s\n' \
    "$pid" "$CURRENT_HEAD" "$RESULT_ROOT" "$output_log" "$PIDFILE"
}

action_launch() {
  common_start_validation
  test ! -e "$RESULT_ROOT" || die "first-launch result root already exists: $RESULT_ROOT"
  mkdir -p "$LOG_DIR"
  local stamp
  stamp="$(date +%Y%m%d-%H%M%S)"
  archive_file "$PIDFILE" "$stamp"
  archive_file "$LOG" "$stamp"
  archive_file "$LAUNCH_META" "$stamp"
  start_controller "$LOG" "$LAUNCH_META" launch
}

action_status() {
  test -d "$REPO/.git" || die "repository is missing: $REPO"
  script_commits
  printf 'Current HEAD: %s\nScript-owning commit: %s\n' \
    "$CURRENT_HEAD" "${SCRIPT_COMMIT:-uncommitted}"
  if test -z "$(git -C "$REPO" status --short)"; then
    echo "Worktree: clean"
  else
    echo "Worktree: dirty"
    git -C "$REPO" status --short
  fi
  test -e "$RESULT_ROOT" && echo "Result root: exists" || echo "Result root: absent"
  local pid=""
  pid="$(controller_pid 2>/dev/null || true)"
  if test -n "$pid"; then
    echo "Controller PID: $pid"
    if kill -0 "$pid" 2>/dev/null; then
      echo "Controller: alive"
      ps -fp "$pid" || true
    else
      echo "Controller: not alive"
    fi
  else
    echo "Controller PID: unavailable"
  fi
  if test -f "$LOG"; then
    echo "Recent log:"
    tail -n 30 "$LOG"
  else
    echo "Log: absent"
  fi
  if test -e "$RESULT_ROOT"; then
    du -sh "$RESULT_ROOT" || true
  fi
}

action_tail() {
  test -f "$LOG" || die "log does not exist: $LOG"
  tail -n 100 -f "$LOG"
}

action_processes() {
  pgrep -af "$PROCESS_PATTERN" || true
  ps -eo pid,ppid,psr,%cpu,%mem,rss,etime,stat,cmd --forest |
    grep -E "$PROCESS_PATTERN" | grep -v grep || true
}

relevant_processes() {
  local current="$$" parent pid line
  declare -A excluded=()
  while [[ "$current" =~ ^[0-9]+$ ]] && test "$current" -gt 1; do
    excluded["$current"]=1
    parent="$(ps -o ppid= -p "$current" 2>/dev/null | tr -d '[:space:]' || true)"
    test -n "$parent" || break
    current="$parent"
  done
  while read -r pid line; do
    test -n "$pid" || continue
    test -n "${excluded[$pid]+x}" || printf '%s %s\n' "$pid" "$line"
  done < <(pgrep -af "$PROCESS_PATTERN" 2>/dev/null || true)
}

descendant_pids() {
  local root_pid="$1" current child
  local -a queue=("$root_pid")
  declare -A seen=(["$root_pid"]=1)
  while ((${#queue[@]})); do
    current="${queue[0]}"
    queue=("${queue[@]:1}")
    while read -r child; do
      test -n "$child" || continue
      if test -z "${seen[$child]+x}"; then
        seen["$child"]=1
        printf '%s\n' "$child"
        queue+=("$child")
      fi
    done < <(pgrep -P "$current" 2>/dev/null || true)
  done
}

action_affinity() {
  local pid
  pid="$(controller_pid)" || die "controller PID is unavailable"
  kill -0 "$pid" 2>/dev/null || die "controller PID $pid is not alive"
  echo "Expected affinity: $EXPECTED_AFFINITY"
  taskset -pc "$pid"
  while read -r child; do
    kill -0 "$child" 2>/dev/null && taskset -pc "$child" || true
  done < <(descendant_pids "$pid")
}

action_progress() {
  V021_RESULT_ROOT="$RESULT_ROOT" Rscript - <<'RS'
root <- Sys.getenv("V021_RESULT_ROOT")
null_default <- function(x, default) {
  if (is.null(x) || length(x) == 0L) default else x
}
safe_read <- function(path, label) {
  if (!file.exists(path)) {
    cat(label, ": not created\n", sep = "")
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) {
    warning(sprintf("Could not read %s (%s): %s", label, path, conditionMessage(e)))
    NULL
  })
}

manifest <- safe_read(file.path(root, "canonical_manifest.rds"), "manifest")
status <- safe_read(file.path(root, "task_status.rds"), "task status")
ledger <- safe_read(file.path(root, "attempt_ledger.rds"), "attempt ledger")

if (!is.null(manifest)) {
  cat("manifest_hash:", null_default(manifest$manifest_hash, "unavailable"), "\n")
  cat("configuration_hash:", null_default(manifest$configuration_hash, "unavailable"), "\n")
  cat("total_tasks:", if (is.data.frame(manifest$tasks)) nrow(manifest$tasks) else NA_integer_, "\n")
}

if (!is.null(status) && is.data.frame(status$tasks)) {
  counts <- table(factor(status$tasks$state,
    levels = c("pending", "running", "completed", "failed")))
  for (name in names(counts)) cat(name, ": ", unname(counts[[name]]), "\n", sep = "")
  cat("exhausted_attempts:", sum(status$tasks$state == "failed" &
    status$tasks$attempt_count >= 2L, na.rm = TRUE), "\n")
  running <- status$tasks[status$tasks$state == "running", , drop = FALSE]
  if (nrow(running) && !is.null(manifest) && is.data.frame(manifest$tasks)) {
    running <- merge(
      running[c("directed_task_id", "task_ordinal", "attempt_count")],
      manifest$tasks[c("directed_task_id", "source", "target", "direction_seed")],
      by = "directed_task_id", sort = FALSE)
    cat("currently_running_tasks:\n")
    print(running, row.names = FALSE)
  } else cat("currently_running_tasks: none\n")
}

if (!is.null(ledger) && is.data.frame(ledger$attempts)) {
  cat("attempt_ledger_rows:", nrow(ledger$attempts), "\n")
  failures <- ledger$attempts[ledger$attempts$ending_state == "failed", , drop = FALSE]
  if (nrow(failures)) {
    cat("recent_failures:\n")
    print(tail(failures, 20L), row.names = FALSE)
  } else cat("recent_failures: none\n")
}

artifact_dir <- file.path(root, "artifacts")
artifact_files <- if (dir.exists(artifact_dir))
  list.files(artifact_dir, pattern = "\\.rds$", full.names = TRUE) else character()
artifacts <- list()
for (path in artifact_files) {
  value <- tryCatch(readRDS(path), error = function(e) {
    warning(sprintf("Skipping unreadable artifact %s: %s", path, conditionMessage(e)))
    NULL
  })
  if (!is.null(value) && is.list(value) &&
      identical(value$artifact_schema, "v021_completion_artifact_v1"))
    artifacts[[length(artifacts) + 1L]] <- value
  else if (!is.null(value)) warning(sprintf("Skipping non-completion artifact %s", path))
}
cat("predictive_eligible_artifacts:", sum(vapply(artifacts,
  function(x) isTRUE(x$predictive_eligibility), logical(1))), "\n")
kfold_states <- vapply(artifacts, function(x)
  as.character(null_default(x$subject_elpd$state, "unavailable")), character(1))
kfold_counts <- table(factor(kfold_states,
  levels = c("pending", "running", "completed", "failed", "not_applicable", "unavailable")))
cat("kfold_states:\n")
print(kfold_counts)

owner <- safe_read(file.path(root, ".v021-controller-lock", "owner.rds"), "ownership")
if (!is.null(owner)) print(owner)
reservations <- safe_read(file.path(root, "chain_reservations.rds"), "chain reservations")
if (!is.null(reservations)) print(reservations)
RS
}

action_storage() {
  local target="$REPO"
  if test -e "$RESULT_ROOT"; then
    target="$RESULT_ROOT"
    du -sh "$RESULT_ROOT"
  else
    echo "Result root: absent"
  fi
  df -h "$target"
  df -i "$target"
  echo "Central retained estimate: 211.3 GB"
  echo "Conservative retained estimate: 423.5 GB"
  echo "Safety-adjusted requirement: 529.4 GB"
  local free_bytes
  free_bytes="$(df -B1 --output=avail "$target" | tail -n 1 | tr -d '[:space:]')"
  if [[ "$free_bytes" =~ ^[0-9]+$ ]]; then
    awk -v free="$free_bytes" -v required=529400000000 \
      'BEGIN { printf "Current free space: %.1f GB\nHeadroom: %.1f GB\n", free/1e9, (free-required)/1e9 }'
  fi
}

action_orphans() {
  if controller_alive; then
    echo "Controller $(controller_pid) is active."
    return 0
  fi
  local matches
  matches="$(relevant_processes)"
  if test -n "$matches"; then
    echo "WARNING: possible V021 orphan processes:" >&2
    printf '%s\n' "$matches" >&2
  else
    echo "No V021 orphan processes found."
  fi
}

action_restart() {
  local reconcile="false"
  case "${1:-}" in
    "") ;;
    --reconcile-stale-lock) reconcile="true" ;;
    *) die "unknown restart option: $1" ;;
  esac
  common_start_validation
  test -d "$RESULT_ROOT" || die "restart result root is missing: $RESULT_ROOT"
  local required
  for required in canonical_manifest.rds task_status.rds attempt_ledger.rds; do
    test -f "$RESULT_ROOT/$required" || die "missing restart state: $required"
  done
  mkdir -p "$LOG_DIR"
  local stamp restart_log restart_meta old_pid matches
  stamp="$(date +%Y%m%d-%H%M%S)"
  restart_log="$LOG_DIR/v021_full_100_species_v5.restart.$stamp.log"
  restart_meta="$LOG_DIR/v021_full_100_species_v5.restart.$stamp.tsv"
  old_pid="$(controller_pid 2>/dev/null || true)"
  if test -n "$old_pid" && kill -0 "$old_pid" 2>/dev/null; then
    die "existing controller is still live: $old_pid"
  fi
  archive_file "$PIDFILE" "$stamp"

  if test -d "$RESULT_ROOT/.v021-controller-lock"; then
    matches="$(relevant_processes)"
    test -z "$matches" || die "V021 processes remain active; refusing stale-lock reconciliation: $matches"
    if test "$reconcile" != "true"; then
      die "stale ownership exists; run: bash $SCRIPT_PATH restart --reconcile-stale-lock"
    fi
    V021_RESULT_ROOT="$RESULT_ROOT" Rscript - <<'RS'
source("/home/heuklang/pclvbayes/benchmarks/mtist/v021_resource_policy.R")
root <- Sys.getenv("V021_RESULT_ROOT")
reconcile_v021_stale_ownership(
  root,
  .v021_local_process_identity(),
  administrative_override = TRUE
)
RS
  fi
  start_controller "$restart_log" "$restart_meta" restart
}

action_help() {
  cat <<'EOF'
Usage:
  v021_full_100_species_v5_ctl.sh launch
  v021_full_100_species_v5_ctl.sh status
  v021_full_100_species_v5_ctl.sh tail
  v021_full_100_species_v5_ctl.sh processes
  v021_full_100_species_v5_ctl.sh affinity
  v021_full_100_species_v5_ctl.sh progress
  v021_full_100_species_v5_ctl.sh storage
  v021_full_100_species_v5_ctl.sh orphans
  v021_full_100_species_v5_ctl.sh restart
  v021_full_100_species_v5_ctl.sh restart --reconcile-stale-lock
  v021_full_100_species_v5_ctl.sh help
EOF
}

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  return 0
fi

case "${1:-help}" in
  launch) shift; test "$#" -eq 0 || die "launch accepts no options"; action_launch ;;
  restart) shift; test "$#" -le 1 || die "too many restart options"; action_restart "${1:-}" ;;
  status) shift; test "$#" -eq 0 || die "status accepts no options"; action_status ;;
  tail) shift; test "$#" -eq 0 || die "tail accepts no options"; action_tail ;;
  processes) shift; test "$#" -eq 0 || die "processes accepts no options"; action_processes ;;
  affinity) shift; test "$#" -eq 0 || die "affinity accepts no options"; action_affinity ;;
  progress) shift; test "$#" -eq 0 || die "progress accepts no options"; action_progress ;;
  storage) shift; test "$#" -eq 0 || die "storage accepts no options"; action_storage ;;
  orphans) shift; test "$#" -eq 0 || die "orphans accepts no options"; action_orphans ;;
  help|-h|--help)
    test "$#" -eq 0 || shift
    test "$#" -eq 0 || die "help accepts no options"
    action_help
    ;;
  *) die "unknown action: $1" ;;
esac
