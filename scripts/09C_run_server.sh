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
RUN_ROOT="${PROJECT_DIR}/logs/09C_external_test"
RESULT_ROOT="${PROJECT_DIR}/results/09C_external_test"
WORK_ROOT="${PROJECT_DIR}/data_processed/09C_external_test"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RESULT_DIR="${RESULT_ROOT}/${RUN_ID}"
WORK_DIR="${WORK_ROOT}/${RUN_ID}"
COMPRESSED_DIR="${WORK_DIR}/test_cel_gz"
CEL_DIR="${WORK_DIR}/test_cel"
LOG_FILE="${RUN_ROOT}/stage_9C_${RUN_ID}.log"
STATE_FILE="${RUN_ROOT}/stage_9C_state.tsv"
STARTED_MARKER="${RUN_ROOT}/ONE_TIME_TEST_STARTED"
COMPLETE_MARKER="${RUN_ROOT}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_ROOT}/NEEDS_CODEX_ATTENTION"
SETUP_GIT_COMMIT="${STAGE9C_SETUP_COMMIT:-UNKNOWN}"
MODEL_PATH="${PROJECT_DIR}/objects/locked_stool_model.rds"
RAW_ARCHIVE="${PROJECT_DIR}/data_raw/GSE99573/GSE99573_RAW.tar"
ACCESS_AUDIT="${PROJECT_DIR}/results_final/stage_9A_locked_test_access_audit.tsv"
EXPECTED_MODEL_SHA="3d4cd825b81fada4d0a3a92907d766dba4cdfea8e3ef7ce926b6790339c861fb"
EXPECTED_RAW_SHA="80cc17d2317d26da0dd739fe76f292b3c2212f5b6ced4c85611f94f3fe2dcdff"

mkdir -p "${RUN_ROOT}"

on_error() {
  code=$?
  {
    printf 'status\tfailed\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'exit_code\t%s\n' "${code}"
    printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
    printf 'log\t%s\n' "${LOG_FILE}"
    printf 'one_time_test_started\t%s\n' "$([[ -e "${STARTED_MARKER}" ]] && printf TRUE || printf FALSE)"
    printf 'model_retrained\tFALSE\n'
    printf 'next_action\tCodex should inspect this same run; do not create a new test evaluation\n'
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_STAGE_9C "Stage 9C failed; see ${FAILED_MARKER}"
  exit "${code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
cd "${PROJECT_DIR}"

# Complete every preflight before the irreversible test-access marker.
test -r AGENTS.md
test -r protocol/stool_model_test_once_plan.md
test -r config/stage_9C_lock_verification.tsv
test -r config/stool_test_operating_points.tsv
test -r "${MODEL_PATH}"
test -r "${RAW_ARCHIVE}"
test -r "${ACCESS_AUDIT}"
test ! -e "${STARTED_MARKER}"
test ! -e "${COMPLETE_MARKER}"
test ! -e "${RESULT_DIR}"
test ! -e "${WORK_DIR}"
test "$(stat -c '%a' "${MODEL_PATH}")" = "444"
test "$(sha256sum "${MODEL_PATH}" | awk '{print $1}')" = "${EXPECTED_MODEL_SHA}"
test "$(sha256sum "${RAW_ARCHIVE}" | awk '{print $1}')" = "${EXPECTED_RAW_SHA}"

MEMBER_LIST="${RUN_ROOT}/test_members_${RUN_ID}.txt"
awk -F $'\t' '
  NR > 1 && $2 == "testing" {
    gsub(/\r$/, "", $3)
    print $3
  }
' "${ACCESS_AUDIT}" > "${MEMBER_LIST}"
test "$(wc -l < "${MEMBER_LIST}")" -eq 65
test "$(sort -u "${MEMBER_LIST}" | wc -l)" -eq 65
tar -tf "${RAW_ARCHIVE}" > "${RUN_ROOT}/raw_members_${RUN_ID}.txt"
while IFS= read -r member; do
  grep -Fxq "${member}" "${RUN_ROOT}/raw_members_${RUN_ID}.txt"
done < "${MEMBER_LIST}"

{
  printf 'status\tpreflight_passed\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_stage9C\n'
  printf 'setup_git_commit\t%s\n' "${SETUP_GIT_COMMIT}"
  printf 'model_sha256\t%s\n' "${EXPECTED_MODEL_SHA}"
  printf 'stool_model_git_tag\tstool-model-locked\n'
  printf 'stool_model_tag_commit\tba5e02da27f964076ac272bf4c23904d0ab892a2\n'
  printf 'test_samples\t65\n'
  printf 'model_retraining\tforbidden\n'
  printf 'monitoring\tserver_background_only\n'
} > "${STATE_FILE}"

set -o noclobber
{
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'irreversible_test_access_authorized\tTRUE\n'
  printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
  printf 'setup_git_commit\t%s\n' "${SETUP_GIT_COMMIT}"
  printf 'model_sha256\t%s\n' "${EXPECTED_MODEL_SHA}"
  printf 'test_samples\t65\n'
} > "${STARTED_MARKER}"
set +o noclobber

mkdir -p "${RESULT_DIR}" "${COMPRESSED_DIR}" "${CEL_DIR}"
tar -xf "${RAW_ARCHIVE}" -C "${COMPRESSED_DIR}" -T "${MEMBER_LIST}"
test "$(find "${COMPRESSED_DIR}" -maxdepth 1 -type f -name '*.CEL.gz' | wc -l)" -eq 65

while IFS= read -r member; do
  source_path="${COMPRESSED_DIR}/${member}"
  output_path="${CEL_DIR}/${member%.gz}"
  test -s "${source_path}"
  gzip -t "${source_path}"
  gzip -dc "${source_path}" > "${output_path}"
  test -s "${output_path}"
done < "${MEMBER_LIST}"
test "$(find "${CEL_DIR}" -maxdepth 1 -type f -name '*.CEL' | wc -l)" -eq 65

Rscript scripts/09C_test_once.R \
  "${PROJECT_DIR}" "${RUN_ID}" "${SETUP_GIT_COMMIT}"
Rscript scripts/09C_validate_outputs.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/09C_finalize_report.R "${PROJECT_DIR}" "${RUN_ID}"

test "$(sha256sum "${MODEL_PATH}" | awk '{print $1}')" = "${EXPECTED_MODEL_SHA}"
test "$(sha256sum "${RAW_ARCHIVE}" | awk '{print $1}')" = "${EXPECTED_RAW_SHA}"
test "$(stat -c '%a' "${MODEL_PATH}")" = "444"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_9C_external_test.md"
  printf 'results\t%s\n' "${RESULT_DIR}/stool_test_results.tsv"
  printf 'model_sha256\t%s\n' "${EXPECTED_MODEL_SHA}"
  printf 'test_samples_evaluated\t65\n'
  printf 'one_time_test_completed\tTRUE\n'
  printf 'model_retrained\tFALSE\n'
  printf 'next_action\tCodex independent QC only; do not reuse the test set\n'
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_9C "Stage 9C completed and is ready for Codex QC"
