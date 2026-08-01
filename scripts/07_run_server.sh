#!/usr/bin/env bash

set -Eeuo pipefail
umask 002

PROJECT_DIR="${CRC_PROJECT_ROOT}"
SETUP_COMMIT="0300a21"
RUN_DIR="${PROJECT_DIR}/logs/07_singlecell_replication"
RESULT_DIR="${PROJECT_DIR}/results/07_singlecell_replication"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RUN_DIR}/stage_7_${RUN_ID}.log"
STATE_FILE="${RUN_DIR}/stage_7_state.tsv"
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
    printf 'next_action\tCodex should inspect the final log and repair only the failed step\n'
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_STAGE_7 "Stage 7 failed; review ${FAILED_MARKER}"
  exit "${exit_code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_stage7\n'
  printf 'setup_commit\t%s\n' "${SETUP_COMMIT}"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'monitoring\tserver_background_only\n'
} > "${STATE_FILE}"

cd "${PROJECT_DIR}"

LOCKED_FILES=(
  "results/06A_pseudobulk/candidate_programs.tsv"
  "results/06A_amendment/exploratory_candidate_modules.tsv"
  "results/06A_amendment/stage_blind_module_membership.tsv"
  "reports/stage_6A_pseudobulk_discovery.md"
  "reports/stage_6A_exploratory_amendment.md"
)
mapfile -t RAW_FILES < <(
  find data_raw/GSE161277 data_raw/GSE132465 -maxdepth 1 -type f -print | sort
)

sha256sum "${LOCKED_FILES[@]}" > "${RESULT_DIR}/locked_inputs.before.sha256"
sha256sum "${RAW_FILES[@]}" > "${RESULT_DIR}/raw_inputs.before.sha256"

echo "Stage 7 started: $(date --iso-8601=seconds)"
Rscript scripts/07_install_dependencies.R "${PROJECT_DIR}"
Rscript scripts/07_syntax_check.R
Rscript scripts/07_smoke_test.R "${PROJECT_DIR}"
Rscript scripts/07_preflight.R "${PROJECT_DIR}"
Rscript scripts/07_real_data_pilot.R "${PROJECT_DIR}"
Rscript scripts/07_prepare_GSE161277.R "${PROJECT_DIR}"
Rscript scripts/07_prepare_GSE132465.R "${PROJECT_DIR}"
Rscript scripts/07_run_replication_models.R "${PROJECT_DIR}"
Rscript scripts/07_finalize_report.R "${PROJECT_DIR}"

sha256sum "${LOCKED_FILES[@]}" > "${RESULT_DIR}/locked_inputs.after.sha256"
sha256sum "${RAW_FILES[@]}" > "${RESULT_DIR}/raw_inputs.after.sha256"
diff -u \
  "${RESULT_DIR}/locked_inputs.before.sha256" \
  "${RESULT_DIR}/locked_inputs.after.sha256"
diff -u \
  "${RESULT_DIR}/raw_inputs.before.sha256" \
  "${RESULT_DIR}/raw_inputs.after.sha256"

Rscript scripts/07_validate_outputs.R "${PROJECT_DIR}"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'validation\t%s\n' "${RESULT_DIR}/validation_checks.tsv"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_7_singlecell_replication.md"
  printf 'instruction\tStop; do not enter Stage 8\n'
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_7 "Stage 7 completed and is ready for Codex QC"
echo "Stage 7 completed: $(date --iso-8601=seconds)"
