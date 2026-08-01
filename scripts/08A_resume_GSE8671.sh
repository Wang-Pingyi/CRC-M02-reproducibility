#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_ID="20260729_060546"
RUN_DIR="${PROJECT_DIR}/logs/08A_bulk_preprocessing"
RESULT_BASE="${PROJECT_DIR}/results/08A_bulk_preprocessing"
WORK_DIR="${PROJECT_DIR}/data_processed/stage_8A_bulk_preprocessing/${RUN_ID}/work/GSE8671"
GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8671/matrix/GSE8671_series_matrix.txt.gz"
GEO_FILE="${WORK_DIR}/GSE8671_series_matrix.txt.gz"
PART_FILE="${GEO_FILE}.part"
LOG_FILE="${RUN_DIR}/stage_8A_${RUN_ID}_resume_GSE8671.log"
STATE_FILE="${RUN_DIR}/stage_8A_state.tsv"
COMPLETE_MARKER="${RUN_DIR}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_DIR}/NEEDS_CODEX_ATTENTION"

mkdir -p "${WORK_DIR}"
rm -f "${COMPLETE_MARKER}" "${FAILED_MARKER}"

on_error() {
  code=$?
  {
    printf 'status\tfailed\nrun_id\t%s\nresume_step\tGSE8671\nexit_code\t%s\ntimestamp\t%s\nlog\t%s\nnext_action\tCodex should inspect the GSE8671 resume only; do not modify data_raw\n' "${RUN_ID}" "${code}" "$(date --iso-8601=seconds)" "${LOG_FILE}"
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  exit "${code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
{
  printf 'status\trunning\nrun_id\t%s\nresume_step\tGSE8671\nstarted\t%s\ntmux_session\tcrc_stage8A\nlog\t%s\n' "${RUN_ID}" "$(date --iso-8601=seconds)" "${LOG_FILE}"
} > "${STATE_FILE}"

cd "${PROJECT_DIR}"
if [[ -s "${GEO_FILE}" ]] && ! gzip -t "${GEO_FILE}" 2>/dev/null; then
  if [[ ! -e "${PART_FILE}" ]]; then
    mv "${GEO_FILE}" "${PART_FILE}"
  fi
fi
if [[ ! -s "${GEO_FILE}" ]] || ! gzip -t "${GEO_FILE}" 2>/dev/null; then
  curl --fail --location --retry 8 --retry-all-errors --retry-delay 10 \
    --connect-timeout 30 --max-time 900 --continue-at - \
    --output "${PART_FILE}" "${GEO_URL}"
  gzip -t "${PART_FILE}"
  mv "${PART_FILE}" "${GEO_FILE}"
fi

Rscript scripts/08A_GSE8671_preprocess.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/08A_finalize_report.R "${PROJECT_DIR}" "${RUN_ID}"

mapfile -t RAW_FILES < <(find data_raw/GSE41657 data_raw/GSE100179 data_raw/GSE8671 -maxdepth 1 -type f -print | sort)
sha256sum "${RAW_FILES[@]}" > "${RESULT_BASE}/raw_inputs_${RUN_ID}.after.sha256"
diff -u "${RESULT_BASE}/raw_inputs_${RUN_ID}.before.sha256" "${RESULT_BASE}/raw_inputs_${RUN_ID}.after.sha256"

{
  printf 'status\tready_for_qc\nrun_id\t%s\ncompleted\t%s\nlog\t%s\nreport\t%s\ninstruction\tStop after Stage 8A; do not run trend testing or Stage 8B\n' "${RUN_ID}" "$(date --iso-8601=seconds)" "${LOG_FILE}" "${PROJECT_DIR}/reports/stage_8A_bulk_preprocessing.md"
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_8A "Stage 8A completed after GSE8671 resume and is ready for Codex QC"
