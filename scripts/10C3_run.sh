#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-${CRC_PROJECT_ROOT}}"
MODE="${2:-full}"
RUN_ID="stage10c3_$(date +%Y%m%d_%H%M%S)"
STATUS_DIR="$ROOT/status"
LOG_DIR="$ROOT/logs/stage10c3"
R_LIB="$ROOT/environment/R-library"

mkdir -p "$STATUS_DIR/history" "$LOG_DIR"
for old_state in RUNNING SUCCESS FAILED; do
  if [[ -f "$STATUS_DIR/STAGE10C3.$old_state" ]]; then
    mv "$STATUS_DIR/STAGE10C3.$old_state" "$STATUS_DIR/history/STAGE10C3.${RUN_ID}.${old_state}"
  fi
done
{
  printf 'run_id=%s\n' "$RUN_ID"
  printf 'started_at=%s\n' "$(date --iso-8601=seconds)"
  printf 'pid=%s\n' "$$"
  printf 'root=%s\n' "$ROOT"
  printf 'mode=%s\n' "$MODE"
} > "$STATUS_DIR/STAGE10C3.RUNNING"

finish() {
  code=$?
  if [[ $code -eq 0 ]]; then
    mv "$STATUS_DIR/STAGE10C3.RUNNING" "$STATUS_DIR/STAGE10C3.SUCCESS"
    printf 'completed_at=%s\n' "$(date --iso-8601=seconds)" >> "$STATUS_DIR/STAGE10C3.SUCCESS"
  else
    mv "$STATUS_DIR/STAGE10C3.RUNNING" "$STATUS_DIR/STAGE10C3.FAILED" 2>/dev/null || true
    printf 'failed_at=%s\nexit_code=%s\n' "$(date --iso-8601=seconds)" "$code" >> "$STATUS_DIR/STAGE10C3.FAILED"
  fi
  exit "$code"
}
trap finish EXIT

export R_LIBS_USER="$R_LIB"
export OMP_NUM_THREADS=24
export OPENBLAS_NUM_THREADS=24
export MKL_NUM_THREADS=24

score_extra=()
if [[ "$MODE" == "full" ]]; then
  Rscript "$ROOT/scripts/10C3_qc_annotation.R" --root "$ROOT"
elif [[ "$MODE" == "score-only" ]]; then
  test -f "$ROOT/results/stage10c3/QC_PREM02.SUCCESS"
elif [[ "$MODE" == "score-resume" ]]; then
  test -f "$ROOT/results/stage10c3/QC_PREM02.SUCCESS"
  score_extra+=(--resume-after-report-failure)
elif [[ "$MODE" == "finalize-only" ]]; then
  test -f "$ROOT/results/stage10c3/STAGE10C3_PATIENT_SAMPLE_SCORES.tsv"
  Rscript "$ROOT/scripts/10C3_finalize_outputs.R" "$ROOT"
  exit 0
else
  printf 'Unknown mode: %s\n' "$MODE" >&2
  exit 2
fi
Rscript "$ROOT/scripts/10C3_score_locked_m02.R" --root "$ROOT" "${score_extra[@]}"
