#!/usr/bin/env bash

set -Eeuo pipefail
umask 002

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_DIR="${PROJECT_DIR}/logs/06B_regulatory_inference"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RUN_DIR}/stage_6B_qc_fix_${RUN_ID}.log"
STATE_FILE="${RUN_DIR}/stage_6B_qc_fix_state.tsv"
READY_MARKER="${RUN_DIR}/READY_FOR_CODEX_FINAL_ACCEPTANCE"
FAILED_MARKER="${RUN_DIR}/QC_FIX_NEEDS_CODEX_ATTENTION"

mkdir -p "${RUN_DIR}"
for marker in "${READY_MARKER}" "${FAILED_MARKER}"; do
  if [[ -e "${marker}" ]]; then
    mv "${marker}" "${marker}.${RUN_ID}.previous"
  fi
done

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
  logger -t CRC_STAGE_6B "Stage 6B communication QC correction failed"
  exit "${exit_code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_6B_qcfix\n'
  printf 'log\t%s\n' "${LOG_FILE}"
} > "${STATE_FILE}"

cd "${PROJECT_DIR}"
cp \
  results/06B_regulatory_inference/prioritized_ligand_receptor.tsv \
  results/06B_regulatory_inference/qc_rejected_axes_pre_category_filter.tsv

Rscript scripts/06B_syntax_check.R
Rscript scripts/06B_liana_donor_communication.R "${PROJECT_DIR}"
Rscript scripts/06B_finalize_report.R "${PROJECT_DIR}"
Rscript scripts/06B_validate_outputs.R "${PROJECT_DIR}"
Rscript scripts/06B_acceptance_qc.R "${PROJECT_DIR}"

{
  printf 'status\tready_for_final_acceptance\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'acceptance\t%s\n' \
    "${PROJECT_DIR}/results/06B_regulatory_inference/independent_acceptance_checks.tsv"
  printf 'instruction\tStop; do not enter external validation\n'
} > "${READY_MARKER}"
cp "${READY_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_6B "Stage 6B communication QC correction ready for final acceptance"
echo "Stage 6B communication QC correction completed"
