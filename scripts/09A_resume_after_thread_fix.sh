#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_ID="20260729_213907"
RUN_ROOT="${PROJECT_DIR}/logs/09A_stool_feasibility"
RESULT_DIR="${PROJECT_DIR}/results/09A_stool_feasibility/${RUN_ID}"
CEL_DIR="${PROJECT_DIR}/data_processed/09A_stool_feasibility/${RUN_ID}/training_cel"
LOG_FILE="${RUN_ROOT}/stage_9A_${RUN_ID}_thread_fix_resume.log"
STATE_FILE="${RUN_ROOT}/stage_9A_state.tsv"
COMPLETE_MARKER="${RUN_ROOT}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_ROOT}/NEEDS_CODEX_ATTENTION"
RAW_TAR="${PROJECT_DIR}/data_raw/GSE99573/GSE99573_RAW.tar"
RAW_SOFT="${PROJECT_DIR}/data_raw/GSE99573/GSE99573_series_metadata.soft.txt"

on_error() {
  code=$?
  {
    printf 'status\tfailed\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'exit_code\t%s\n' "${code}"
    printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
    printf 'log\t%s\n' "${LOG_FILE}"
    printf 'resume\tthread_fix\n'
    printf 'next_action\tCodex should inspect this failure; do not access test CEL files or start model training\n'
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_STAGE_9A "Stage 9A thread-fix resume failed"
  exit "${code}"
}
trap on_error ERR

mkdir -p "${RUN_ROOT}"
if [[ -e "${FAILED_MARKER}" ]]; then
  mv "${FAILED_MARKER}" "${FAILED_MARKER}.${RUN_ID}.pthread.previous"
fi
exec > >(tee -a "${LOG_FILE}") 2>&1
cd "${PROJECT_DIR}"

test "$(find "${CEL_DIR}" -type f -iname '*.CEL' | wc -l)" -eq 265
test "$(find "${CEL_DIR}" -type f -iname '*.CEL.gz' | wc -l)" -eq 265
{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'resumed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_stage9A\n'
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'resume_reason\tpreprocessCore pthread_create EINVAL\n'
  printf 'test_expression_access\tforbidden\n'
} > "${STATE_FILE}"

bash scripts/09A_install_r_packages.sh "${PROJECT_DIR}"
export R_LIBS_USER="${PROJECT_DIR}/environment/Rlib_stage9A"
Rscript scripts/09A_stool_feasibility.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/09A_finalize_report.R "${PROJECT_DIR}" "${RUN_ID}"

if find "${CEL_DIR}" -type f | grep -F -f <(
  awk -F '\t' 'NR > 1 && $3 == "testing" {print $1 "_"}' \
    "${RESULT_DIR}/GSE99573_sample_inventory.tsv"
) | grep -q .; then
  echo "Test-set file found in training work directory"
  exit 14
fi
sha256sum "${RAW_TAR}" "${RAW_SOFT}" > "${RESULT_DIR}/raw_inputs.after.sha256"
diff -u \
  "${RESULT_DIR}/raw_inputs.before.sha256" \
  "${RESULT_DIR}/raw_inputs.after.sha256"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_9A_stool_feasibility.md"
  printf 'test_expression_accessed\tFALSE\n'
  printf 'resume\tthread_fix\n'
  printf 'next_action\tCodex independent QC only; do not start model training\n'
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_9A "Stage 9A completed after single-thread preprocessCore fix"
