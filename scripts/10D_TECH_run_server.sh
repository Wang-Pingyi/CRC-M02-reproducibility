#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

ROOT="${1:-${CRC_PROJECT_ROOT}}"
RUN_ID="stage10d_tech_$(date +%Y%m%d_%H%M%S)"
STATUS_DIR="$ROOT/status"
LOG_DIR="$ROOT/logs/stage10d_tech"
RESULT_DIR="$ROOT/results/stage10d_tech"
REPORT_DIR="$ROOT/reports"
LOG_FILE="$LOG_DIR/${RUN_ID}.log"
RUN_RECORD="$LOG_DIR/${RUN_ID}.tsv"

mkdir -p "$STATUS_DIR/history" "$LOG_DIR" "$REPORT_DIR" "$RESULT_DIR"

for protected in \
  STAGE10D_TECH_PRIMARY_SCORE_RESULTS.tsv \
  STAGE10D_TECH_COVERAGE_RESULTS.tsv \
  STAGE10D_TECH_SENSITIVITY_RESULTS.tsv \
  STAGE10D_TECH_NULL_MODULE_RESULTS.tsv \
  STAGE10D_TECH_DECISION.md; do
  if [[ -e "$RESULT_DIR/$protected" ]]; then
    printf 'Refusing to overwrite existing Stage 10D-TECH result: %s\n' "$RESULT_DIR/$protected" >&2
    exit 73
  fi
done

for state in RUNNING SUCCESS FAILED; do
  marker="$STATUS_DIR/STAGE10D_TECH.$state"
  if [[ -e "$marker" ]]; then
    mv "$marker" "$STATUS_DIR/history/STAGE10D_TECH.${RUN_ID}.${state}"
  fi
done

{
  printf 'field\tvalue\n'
  printf 'run_id\t%s\n' "$RUN_ID"
  printf 'status\tRUNNING\n'
  printf 'started_at\t%s\n' "$(date --iso-8601=seconds)"
  printf 'pid\t%s\n' "$$"
  printf 'root\t%s\n' "$ROOT"
  printf 'log\t%s\n' "$LOG_FILE"
  printf 'cpu_thread_cap\t24\n'
  printf 'virtual_memory_limit_kb\t180000000\n'
  printf 'seed\t20260731\n'
  printf 'expected_primary_output\t%s\n' "$RESULT_DIR/STAGE10D_TECH_PRIMARY_SCORE_RESULTS.tsv"
  printf 'completion_check\ttest -f %s\n' "$STATUS_DIR/STAGE10D_TECH.SUCCESS"
  printf 'failure_check\ttest -f %s\n' "$STATUS_DIR/STAGE10D_TECH.FAILED"
  printf 'recovery\treview log and submit a new run ID; never overwrite completed outputs\n'
} > "$RUN_RECORD"
cp "$RUN_RECORD" "$STATUS_DIR/STAGE10D_TECH.RUNNING"

finish() {
  code=$?
  if [[ $code -eq 0 ]]; then
    mv "$STATUS_DIR/STAGE10D_TECH.RUNNING" "$STATUS_DIR/STAGE10D_TECH.SUCCESS"
    {
      printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
      printf 'validation\tPASS\n'
      printf 'next_stage_started\tNO\n'
    } >> "$STATUS_DIR/STAGE10D_TECH.SUCCESS"
  else
    mv "$STATUS_DIR/STAGE10D_TECH.RUNNING" "$STATUS_DIR/STAGE10D_TECH.FAILED" 2>/dev/null || true
    {
      printf 'failed_at\t%s\n' "$(date --iso-8601=seconds)"
      printf 'exit_code\t%s\n' "$code"
      printf 'log\t%s\n' "$LOG_FILE"
      printf 'next_action\tCodex should inspect only the final log, repair the failed step and use a new run ID\n'
    } >> "$STATUS_DIR/STAGE10D_TECH.FAILED"
    {
      printf '# Stage 10D-TECH server summary\n\n'
      printf -- '- Status: FAILED\n'
      printf -- '- Run ID: %s\n' "$RUN_ID"
      printf -- '- Log: %s\n' "$LOG_FILE"
      printf -- '- No later stage was started.\n'
    } > "$REPORT_DIR/STAGE10D_TECH_SUMMARY.md"
  fi
  exit "$code"
}
trap finish EXIT

exec > >(tee -a "$LOG_FILE") 2>&1

unset R_LIBS_USER
export OMP_NUM_THREADS=24
export OPENBLAS_NUM_THREADS=24
export MKL_NUM_THREADS=24
export VECLIB_MAXIMUM_THREADS=24
ulimit -v 180000000

cd "$ROOT"
printf 'Stage 10D-TECH started: %s\n' "$(date --iso-8601=seconds)"
python3 scripts/10C2_SP_validate_outputs.py
test -f cache/stage10d_tech_smoke/SMOKE.SUCCESS
Rscript --vanilla scripts/10D_TECH_pseudospot_validation.R --root "$ROOT"
python3 scripts/10D_TECH_validate_outputs.py "$ROOT"
printf 'Stage 10D-TECH completed: %s\n' "$(date --iso-8601=seconds)"
