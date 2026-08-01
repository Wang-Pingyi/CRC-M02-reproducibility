#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

ROOT="${1:-${CRC_PROJECT_ROOT}}"
RUN_ID="stage10e_$(date +%Y%m%d_%H%M%S)"
STAGE="STAGE10E"
GIT_COMMIT="${STAGE10E_GIT_COMMIT:-UNRECORDED}"
STATUS_DIR="$ROOT/status"
LOG_DIR="$ROOT/logs/stage10e"
RESULT_DIR="$ROOT/results/stage10e"
DERIVED_DIR="$ROOT/data_processed/stage10e/$RUN_ID/derived"
LOG_FILE="$LOG_DIR/${RUN_ID}.log"
RUN_RECORD="$STATUS_DIR/${STAGE}_SUBMISSION.tsv"
RAW_ZIP="$ROOT/data_raw/stage10B_20260730_205808/E-GEAD-622/E-GEAD-622.processed.zip"

mkdir -p "$STATUS_DIR/history" "$LOG_DIR" "$RESULT_DIR" "$DERIVED_DIR" "$ROOT/reports"
for protected in STAGE10E_SPATIAL_QC.tsv STAGE10E_EXCLUSION_LOG.tsv STAGE10E_DECONVOLUTION_QC.tsv STAGE10E_GENESET_INTEGRITY_CHECK.tsv STAGE10E_ANALYSIS_READY_MANIFEST.tsv STAGE10E_DECISION.md; do
  [[ ! -e "$RESULT_DIR/$protected" ]] || { echo "Refusing to overwrite $RESULT_DIR/$protected" >&2; exit 73; }
done
for state in RUNNING SUCCESS FAILED; do
  [[ -e "$STATUS_DIR/$STAGE.$state" ]] && mv "$STATUS_DIR/$STAGE.$state" "$STATUS_DIR/history/$STAGE.$RUN_ID.$state"
done

{
  printf 'field\tvalue\n'
  printf 'stage\t%s\nrun_id\t%s\nscheduler\ttmux\njob_id_or_pid\t%s\nlaunch_time\t%s\n' "$STAGE" "$RUN_ID" "${TMUX_PANE:-tmux}" "$(date --iso-8601=seconds)"
  printf 'git_commit\t%s\n' "$GIT_COMMIT"
  printf 'launch_command\tbash scripts/10E_run_server.sh %s\n' "$ROOT"
  printf 'wrapper_script\tscripts/10E_run_server.sh\n'
  printf 'expected_outputs\tresults/stage10e/STAGE10E_DECISION.md;data/metadata/spatial_patient_slide_roi_manifest.tsv\n'
  printf 'full_log\t%s\n' "$LOG_FILE"
  printf 'status_check_command\ttest -f %s/%s.SUCCESS && cat %s/%s.SUCCESS\n' "$STATUS_DIR" "$STAGE" "$STATUS_DIR" "$STAGE"
  printf 'failure_recovery_command\treview %s then submit a new tmux run; preserve this run-specific derived directory\n' "$LOG_FILE"
  printf 'resource_request\t24 CPU threads; virtual memory 180000000 KB; GPU none\n'
  printf 'smoke_test\tcache/stage10e_smoke/STAGE10E_SMOKE.SUCCESS\n'
} > "$RUN_RECORD"
cp "$RUN_RECORD" "$STATUS_DIR/$STAGE.RUNNING"
finish() {
  code=$?
  if [[ $code -eq 0 ]]; then
    mv "$STATUS_DIR/$STAGE.RUNNING" "$STATUS_DIR/$STAGE.SUCCESS"
    printf 'completed_at\t%s\nvalidation\tPASS\nnext_stage_started\tNO\n' "$(date --iso-8601=seconds)" >> "$STATUS_DIR/$STAGE.SUCCESS"
  else
    mv "$STATUS_DIR/$STAGE.RUNNING" "$STATUS_DIR/$STAGE.FAILED" 2>/dev/null || true
    printf 'failed_at\t%s\nexit_code\t%s\nlog\t%s\nrecovery\trepair only the failing step in a new run ID; do not overwrite outputs\n' "$(date --iso-8601=seconds)" "$code" "$LOG_FILE" >> "$STATUS_DIR/$STAGE.FAILED"
  fi
  exit "$code"
}
trap finish EXIT
exec > >(tee -a "$LOG_FILE") 2>&1
unset R_LIBS_USER
export OMP_NUM_THREADS=24 OPENBLAS_NUM_THREADS=24 MKL_NUM_THREADS=24 VECLIB_MAXIMUM_THREADS=24
ulimit -v 180000000
cd "$ROOT"
[[ -f "$RAW_ZIP" ]] || { echo "Missing raw archive" >&2; exit 2; }
for slide in Rectum_kyudai_Beppu_20200303 Ascending_kyudai_Beppu_20200430 Sigmoid_kyudai_Beppu_20210602 Transverse_kyudai_Beppu_20211111; do
  mkdir -p "$DERIVED_DIR/$slide"
  unzip -p "$RAW_ZIP" "$slide.tar" | tar -xf - -C "$DERIVED_DIR"
done
mkdir -p "$DERIVED_DIR/source_render"
pdftoppm -f 1 -l 1 -png -r 120 "$ROOT/metadata/stage10C2_SP/sources/mmc1.pdf" "$DERIVED_DIR/source_render/mmc1"
Rscript --vanilla scripts/10E_spatial_blind_qc.R --root "$ROOT" --run-id "$RUN_ID" --input-dir "$DERIVED_DIR"
python3 scripts/10E_validate_outputs.py "$ROOT"
