#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

if [[ $# -ne 1 ]]; then
  printf 'Usage: %s <run_id>\n' "$0" >&2
  exit 2
fi

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_ID="$1"
RUN_ROOT="${PROJECT_DIR}/logs/09B_model_training"
RESULT_DIR="${PROJECT_DIR}/results/09B_model_training/${RUN_ID}"
LOG_FILE="${RUN_ROOT}/stage_9B_${RUN_ID}_validation_resume.log"
STATE_FILE="${RUN_ROOT}/stage_9B_state.tsv"
COMPLETE_MARKER="${RUN_ROOT}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_ROOT}/NEEDS_CODEX_ATTENTION"
MODEL_PATH="${PROJECT_DIR}/objects/locked_stool_model.rds"

on_error() {
  code=$?
  {
    printf 'status\tfailed\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'exit_code\t%s\n' "${code}"
    printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
    printf 'log\t%s\n' "${LOG_FILE}"
    printf 'test_expression_accessed\tFALSE\n'
    printf 'next_action\tCodex should inspect validation/report failure; do not run the test set\n'
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  exit "${code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
cd "${PROJECT_DIR}"
test -d "${RESULT_DIR}"
test -r "${MODEL_PATH}"
test -r "${RESULT_DIR}/nested_cv_predictions_training_only.tsv"
test -r "${RESULT_DIR}/nested_cv_performance.tsv"
test -r "${RESULT_DIR}/locked_stool_model.sha256"

rm -f "${COMPLETE_MARKER}"
if [[ -e "${FAILED_MARKER}" ]]; then
  mv "${FAILED_MARKER}" "${FAILED_MARKER}.${RUN_ID}.validation_logic.previous"
fi
{
  printf 'status\tvalidating_existing_locked_model\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'model_retrained\tFALSE\n'
  printf 'test_expression_accessed\tFALSE\n'
} > "${STATE_FILE}"

Rscript scripts/09B_validate_locked_model.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/09B_finalize_report.R "${PROJECT_DIR}" "${RUN_ID}"
sha256sum -c "${RESULT_DIR}/locked_stool_model.sha256"
test ! -w "${MODEL_PATH}"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_9B_model_training.md"
  printf 'model\t%s\n' "${MODEL_PATH}"
  printf 'model_sha256\t%s\n' "$(awk '{print $1}' "${RESULT_DIR}/locked_stool_model.sha256")"
  printf 'model_retrained\tFALSE\n'
  printf 'test_expression_accessed\tFALSE\n'
  printf 'git_tag\tpending_desktop_qc_stool-model-locked\n'
  printf 'next_action\tCodex independent QC only; do not run the test set\n'
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
