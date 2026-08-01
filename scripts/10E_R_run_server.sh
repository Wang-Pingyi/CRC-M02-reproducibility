#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:?project root required}"
RUN_ID="${2:?run ID required}"
LOG_DIR="$ROOT/logs/stage10e_roi_remediation/$RUN_ID"
OUT_DIR="$ROOT/results/stage10e_roi_remediation"
STATUS_DIR="$ROOT/status"
mkdir -p "$LOG_DIR" "$OUT_DIR" "$STATUS_DIR"

RUNNING="$STATUS_DIR/STAGE10E_R.RUNNING"
SUCCESS="$STATUS_DIR/STAGE10E_R.SUCCESS"
FAILED="$STATUS_DIR/STAGE10E_R.FAILED"
rm -f "$SUCCESS" "$FAILED"
printf 'run_id=%s\nstarted_utc=%s\n' "$RUN_ID" "$(date -u +%FT%TZ)" > "$RUNNING"

on_fail() {
  code=$?
  printf 'run_id=%s\nfailed_utc=%s\nexit_code=%s\n' "$RUN_ID" "$(date -u +%FT%TZ)" "$code" > "$FAILED"
  rm -f "$RUNNING"
  exit "$code"
}
trap on_fail ERR

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

python3 "$ROOT/scripts/10E_R_roi_registration.py" --root "$ROOT" --mode formal --run-id "$RUN_ID" \
  > "$LOG_DIR/formal.stdout.log" 2> "$LOG_DIR/formal.stderr.log"
python3 "$ROOT/scripts/10E_R_validate_outputs.py" --root "$ROOT" --run-id "$RUN_ID" \
  > "$LOG_DIR/validator.stdout.log" 2> "$LOG_DIR/validator.stderr.log"

printf 'run_id=%s\ncompleted_utc=%s\nvalidator=PASS\n' "$RUN_ID" "$(date -u +%FT%TZ)" > "$SUCCESS"
rm -f "$RUNNING"
trap - ERR
