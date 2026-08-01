#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT="${1:?project path required}"
RUN_ID="${2:?run id required}"
SETUP_COMMIT="${3:-not_recorded}"
LOG_DIR="$PROJECT/logs/08B_bulk_validation"
RUN_LOG="$LOG_DIR/stage_8B_${RUN_ID}.log"
READY="$LOG_DIR/READY_FOR_CODEX_QC"
ATTENTION="$LOG_DIR/NEEDS_CODEX_ATTENTION"
mkdir -p "$LOG_DIR"
rm -f "$READY" "$ATTENTION"

on_error() {
  local status=$?
  printf 'Stage 8B failed (exit %s) at %s. See %s\n' "$status" "$(date --iso-8601=seconds)" "$RUN_LOG" > "$ATTENTION"
  exit "$status"
}
trap on_error ERR

exec > >(tee -a "$RUN_LOG") 2>&1
echo "Stage 8B start: $(date --iso-8601=seconds)"
echo "Run ID: $RUN_ID"
echo "Local setup Git commit: $SETUP_COMMIT"

R_LIB_USER="${R_LIBS_USER:-${SERVER_STORAGE_ROOT}/R/x86_64-pc-linux-gnu-library/4.3}"
export R_LIBS_USER
if ! Rscript -e 'quit(status=ifelse(requireNamespace("CMScaller", quietly=TRUE),0,1))'; then
  R CMD INSTALL --library="$R_LIB_USER" "$PROJECT/environment/sources/08B/CMScaller-main.tar.gz"
fi

Rscript "$PROJECT/scripts/08B_bulk_cohorts.R" "$PROJECT" "$RUN_ID"
Rscript "$PROJECT/scripts/08B_tcga_auxiliary.R" "$PROJECT" "$RUN_ID"
Rscript "$PROJECT/scripts/08B_meta_gate_report.R" "$PROJECT" "$RUN_ID"
Rscript "$PROJECT/scripts/08B_validate_outputs.R" "$PROJECT" "$RUN_ID"
echo "Stage 8B complete: $(date --iso-8601=seconds)"
trap - ERR
