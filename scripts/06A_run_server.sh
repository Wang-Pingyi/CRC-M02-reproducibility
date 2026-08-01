#!/usr/bin/env bash

set -Eeuo pipefail
umask 002

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_DIR="${PROJECT_DIR}/logs/06A_pseudobulk"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RUN_DIR}/stage_6A_${RUN_ID}.log"
STATE_FILE="${RUN_DIR}/stage_6A_state.tsv"
COMPLETE_MARKER="${RUN_DIR}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_DIR}/NEEDS_CODEX_ATTENTION"

mkdir -p "${RUN_DIR}"

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
  {
    printf 'status\tfailed\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
    printf 'log\t%s\n' "${LOG_FILE}"
  } > "${STATE_FILE}"
  logger -t CRC_6A "Stage 6A failed; review ${FAILED_MARKER}"
  exit "${exit_code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1

{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_6A_pseudobulk\n'
  printf 'log\t%s\n' "${LOG_FILE}"
} > "${STATE_FILE}"

cd "${PROJECT_DIR}"
echo "Stage 6A server run started: $(date --iso-8601=seconds)"
Rscript scripts/06A_pseudobulk_discovery.R "${PROJECT_DIR}"
Rscript scripts/06A_validate_outputs.R "${PROJECT_DIR}"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'validation\t%s\n' "${PROJECT_DIR}/results/06A_pseudobulk/validation_checks.tsv"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_6A_pseudobulk_discovery.md"
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_6A "Stage 6A completed and is ready for Codex QC"
echo "Stage 6A server run completed: $(date --iso-8601=seconds)"
