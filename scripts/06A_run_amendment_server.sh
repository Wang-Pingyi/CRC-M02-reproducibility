#!/usr/bin/env bash

set -Eeuo pipefail
umask 002

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_DIR="${PROJECT_DIR}/logs/06A_amendment"
RESULT_DIR="${PROJECT_DIR}/results/06A_amendment"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RUN_DIR}/stage_6A_amendment_${RUN_ID}.log"
STATE_FILE="${RUN_DIR}/stage_6A_amendment_state.tsv"
COMPLETE_MARKER="${RUN_DIR}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_DIR}/NEEDS_CODEX_ATTENTION"

mkdir -p "${RUN_DIR}" "${RESULT_DIR}"

archive_marker() {
  local marker="$1"
  if [[ -e "${marker}" ]]; then
    mv "${marker}" "${marker}.${RUN_ID}.previous"
  fi
}
archive_marker "${COMPLETE_MARKER}"
archive_marker "${FAILED_MARKER}"

on_error() {
  local exit_code=$?
  {
    printf 'status\tfailed\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'exit_code\t%s\n' "${exit_code}"
    printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
    printf 'log\t%s\n' "${LOG_FILE}"
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_6A_AMENDMENT "Stage 6A amendment failed; review ${FAILED_MARKER}"
  exit "${exit_code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_6A_amendment\n'
  printf 'log\t%s\n' "${LOG_FILE}"
} > "${STATE_FILE}"

cd "${PROJECT_DIR}"

PRIMARY_FILES=(
  "results/06A_pseudobulk/pseudobulk_results.tsv"
  "results/06A_pseudobulk/candidate_programs.tsv"
  "reports/stage_6A_pseudobulk_discovery.md"
)
sha256sum "${PRIMARY_FILES[@]}" > "${RESULT_DIR}/primary_freeze_sha256.before.tsv"

echo "Stage 6A amendment started: $(date --iso-8601=seconds)"
Rscript scripts/06A_amendment_exploratory.R "${PROJECT_DIR}"

sha256sum "${PRIMARY_FILES[@]}" > "${RESULT_DIR}/primary_freeze_sha256.after.tsv"
diff -u \
  "${RESULT_DIR}/primary_freeze_sha256.before.tsv" \
  "${RESULT_DIR}/primary_freeze_sha256.after.tsv"

Rscript scripts/06A_validate_amendment.R "${PROJECT_DIR}"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'validation\t%s\n' "${RESULT_DIR}/validation_checks.tsv"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_6A_exploratory_amendment.md"
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_6A_AMENDMENT "Stage 6A amendment completed and is ready for Codex QC"
echo "Stage 6A amendment completed: $(date --iso-8601=seconds)"
