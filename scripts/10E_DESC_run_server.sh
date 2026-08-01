#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:?project root required}"
RUN_ID="${2:?run ID required}"
LOG_DIR="$ROOT/logs/stage10e_desc/$RUN_ID"
STATUS_DIR="$ROOT/status"
mkdir -p "$LOG_DIR" "$STATUS_DIR" "$ROOT/results/stage10e_desc" "$ROOT/figures/stage10e_desc/source_data"

RUNNING="$STATUS_DIR/STAGE10E_DESC.RUNNING"
SUCCESS="$STATUS_DIR/STAGE10E_DESC.SUCCESS"
FAILED="$STATUS_DIR/STAGE10E_DESC.FAILED"
rm -f "$SUCCESS" "$FAILED"
printf 'stage=STAGE10E_DESC\nrun_id=%s\nstarted_utc=%s\n' "$RUN_ID" "$(date -u +%FT%TZ)" > "$RUNNING"

on_fail() {
  code=$?
  printf 'stage=STAGE10E_DESC\nrun_id=%s\nfailed_utc=%s\nexit_code=%s\n' "$RUN_ID" "$(date -u +%FT%TZ)" "$code" > "$FAILED"
  rm -f "$RUNNING"
  exit "$code"
}
trap on_fail ERR

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

Rscript "$ROOT/scripts/10E_DESC_case4_descriptive.R" "$ROOT" "$RUN_ID" formal \
  > "$LOG_DIR/formal.stdout.log" 2> "$LOG_DIR/formal.stderr.log"
python3 "$ROOT/scripts/10E_DESC_validate_outputs.py" --root "$ROOT" --run-id "$RUN_ID" \
  > "$LOG_DIR/validator.stdout.log" 2> "$LOG_DIR/validator.stderr.log"

printf 'stage=STAGE10E_DESC\nrun_id=%s\ncompleted_utc=%s\nvalidator=PASS\n' "$RUN_ID" "$(date -u +%FT%TZ)" > "$SUCCESS"
rm -f "$RUNNING"
trap - ERR
