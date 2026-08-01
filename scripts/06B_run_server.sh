#!/usr/bin/env bash

set -Eeuo pipefail
umask 002

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_DIR="${PROJECT_DIR}/logs/06B_regulatory_inference"
RESULT_DIR="${PROJECT_DIR}/results/06B_regulatory_inference"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RUN_DIR}/stage_6B_${RUN_ID}.log"
STATE_FILE="${RUN_DIR}/stage_6B_state.tsv"
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
    printf 'next_action\tCodex should inspect the final log and modify the failing step before retrying\n'
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_STAGE_6B "Stage 6B failed; review ${FAILED_MARKER}"
  exit "${exit_code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_6B\n'
  printf 'log\t%s\n' "${LOG_FILE}"
} > "${STATE_FILE}"

cd "${PROJECT_DIR}"

LOCKED_FILES=(
  "results/06A_amendment/exploratory_candidate_modules.tsv"
  "results/06A_amendment/stage_blind_module_membership.tsv"
  "results/06A_amendment/source_data/stage_blind_module_scores.tsv"
  "reports/stage_6A_exploratory_amendment.md"
  "reports/stage_6A_pseudobulk_discovery.md"
)
sha256sum "${LOCKED_FILES[@]}" > "${RESULT_DIR}/stage_6A_locked_inputs.before.sha256"

echo "Stage 6B started: $(date --iso-8601=seconds)"
Rscript scripts/06B_install_dependencies.R "${PROJECT_DIR}"
Rscript scripts/06B_syntax_check.R
Rscript scripts/06B_smoke_test.R "${PROJECT_DIR}"
Rscript scripts/06B_preflight.R "${PROJECT_DIR}"
Rscript scripts/06B_real_data_pilot.R "${PROJECT_DIR}"
Rscript scripts/06B_trajectory_tradeSeq.R "${PROJECT_DIR}"
Rscript scripts/06B_regulator_pathway.R "${PROJECT_DIR}"
Rscript scripts/06B_liana_donor_communication.R "${PROJECT_DIR}"
Rscript scripts/06B_finalize_report.R "${PROJECT_DIR}"

sha256sum "${LOCKED_FILES[@]}" > "${RESULT_DIR}/stage_6A_locked_inputs.after.sha256"
diff -u \
  "${RESULT_DIR}/stage_6A_locked_inputs.before.sha256" \
  "${RESULT_DIR}/stage_6A_locked_inputs.after.sha256"

Rscript scripts/06B_validate_outputs.R "${PROJECT_DIR}"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'validation\t%s\n' "${RESULT_DIR}/validation_checks.tsv"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_6B_regulatory_inference.md"
  printf 'instruction\tStop; do not enter external validation\n'
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_6B "Stage 6B completed and is ready for Codex QC"
echo "Stage 6B completed: $(date --iso-8601=seconds)"
