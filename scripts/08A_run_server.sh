#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

# Keep RMA/preprocessCore and linked numerical libraries within the server's
# per-process thread allowance. Stage 8A is sample-light and gains little from
# nested BLAS/OpenMP threads.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_DIR="${PROJECT_DIR}/logs/08A_bulk_preprocessing"
RESULT_BASE="${PROJECT_DIR}/results/08A_bulk_preprocessing"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RUN_DIR}/stage_8A_${RUN_ID}.log"
STATE_FILE="${RUN_DIR}/stage_8A_state.tsv"
COMPLETE_MARKER="${RUN_DIR}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_DIR}/NEEDS_CODEX_ATTENTION"

mkdir -p "${RUN_DIR}" "${RESULT_BASE}"
for marker in "${COMPLETE_MARKER}" "${FAILED_MARKER}"; do
  [[ -e "${marker}" ]] && mv "${marker}" "${marker}.${RUN_ID}.previous"
done

on_error() {
  code=$?
  {
    printf 'status\tfailed\nrun_id\t%s\nexit_code\t%s\ntimestamp\t%s\nlog\t%s\nnext_action\tCodex should inspect only the failed preprocessing step; do not modify data_raw\n' "${RUN_ID}" "${code}" "$(date --iso-8601=seconds)" "${LOG_FILE}"
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_STAGE_8A "Stage 8A failed; see ${FAILED_MARKER}"
  exit "${code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
cd "${PROJECT_DIR}"
{
  printf 'status\trunning\nrun_id\t%s\nstarted\t%s\ntmux_session\tcrc_stage8A\nlog\t%s\nmonitoring\tserver_background_only\n' "${RUN_ID}" "$(date --iso-8601=seconds)" "${LOG_FILE}"
} > "${STATE_FILE}"

mapfile -t RAW_FILES < <(find data_raw/GSE41657 data_raw/GSE100179 data_raw/GSE8671 -maxdepth 1 -type f -print | sort)
sha256sum "${RAW_FILES[@]}" > "${RESULT_BASE}/raw_inputs_${RUN_ID}.before.sha256"
Rscript scripts/08A_GSE41657_preprocess.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/08A_GSE100179_preprocess.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/08A_GSE8671_preprocess.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/08A_finalize_report.R "${PROJECT_DIR}" "${RUN_ID}"
sha256sum "${RAW_FILES[@]}" > "${RESULT_BASE}/raw_inputs_${RUN_ID}.after.sha256"
diff -u "${RESULT_BASE}/raw_inputs_${RUN_ID}.before.sha256" "${RESULT_BASE}/raw_inputs_${RUN_ID}.after.sha256"
{
  printf 'status\tready_for_qc\nrun_id\t%s\ncompleted\t%s\nlog\t%s\nreport\t%s\ninstruction\tStop after Stage 8A; do not run trend testing or Stage 8B\n' "${RUN_ID}" "$(date --iso-8601=seconds)" "${LOG_FILE}" "${PROJECT_DIR}/reports/stage_8A_bulk_preprocessing.md"
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_8A "Stage 8A completed and is ready for Codex QC"
