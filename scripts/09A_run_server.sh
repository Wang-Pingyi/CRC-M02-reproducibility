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
RUN_ROOT="${PROJECT_DIR}/logs/09A_stool_feasibility"
RESULT_ROOT="${PROJECT_DIR}/results/09A_stool_feasibility"
WORK_ROOT="${PROJECT_DIR}/data_processed/09A_stool_feasibility"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RUN_ROOT}/stage_9A_${RUN_ID}.log"
STATE_FILE="${RUN_ROOT}/stage_9A_state.tsv"
COMPLETE_MARKER="${RUN_ROOT}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_ROOT}/NEEDS_CODEX_ATTENTION"
RAW_TAR="${PROJECT_DIR}/data_raw/GSE99573/GSE99573_RAW.tar"
RAW_SOFT="${PROJECT_DIR}/data_raw/GSE99573/GSE99573_series_metadata.soft.txt"
RESULT_DIR="${RESULT_ROOT}/${RUN_ID}"
WORK_DIR="${WORK_ROOT}/${RUN_ID}"
CEL_DIR="${WORK_DIR}/training_cel"

mkdir -p "${RUN_ROOT}" "${RESULT_DIR}" "${CEL_DIR}"
for marker in "${COMPLETE_MARKER}" "${FAILED_MARKER}"; do
  if [[ -e "${marker}" ]]; then
    mv "${marker}" "${marker}.${RUN_ID}.previous"
  fi
done

on_error() {
  code=$?
  {
    printf 'status\tfailed\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'exit_code\t%s\n' "${code}"
    printf 'timestamp\t%s\n' "$(date --iso-8601=seconds)"
    printf 'log\t%s\n' "${LOG_FILE}"
    printf 'next_action\tCodex should inspect this failure; do not access test CEL files or start model training\n'
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_STAGE_9A "Stage 9A failed; see ${FAILED_MARKER}"
  exit "${code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
cd "${PROJECT_DIR}"
{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_stage9A\n'
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'test_expression_access\tforbidden\n'
  printf 'monitoring\tserver_background_only\n'
} > "${STATE_FILE}"

test -r AGENTS.md
test -r protocol/analysis_protocol_v1.md
test -r protocol/stool_model_protocol.md
test -r reports/stage_4A_metadata_audit.md
test -r reports/stage_8B_bulk_validation.md
test -r "${RAW_TAR}"
test -r "${RAW_SOFT}"
if find data_raw/GSE99573 -maxdepth 1 -type f -perm /222 | grep -q .; then
  echo "GSE99573 raw files are writable; aborting"
  exit 11
fi

sha256sum "${RAW_TAR}" "${RAW_SOFT}" > "${RESULT_DIR}/raw_inputs.before.sha256"
python3 scripts/09A_prepare_split_audit.py "${PROJECT_DIR}" "${RUN_ID}"

# Install from pre-transferred official source tarballs. This avoids dependence
# on an unstable server-side Bioconductor connection during the long job.
bash scripts/09A_install_r_packages.sh "${PROJECT_DIR}"
export R_LIBS_USER="${PROJECT_DIR}/environment/Rlib_stage9A"

# Selective archive extraction is the test firewall: the file list contains
# exactly the 265 frozen training members and no test or Not Used CEL.
tar -xf "${RAW_TAR}" \
  -C "${CEL_DIR}" \
  --files-from="${RESULT_DIR}/training_tar_members.txt"
extracted_gz="$(find "${CEL_DIR}" -type f -iname '*.CEL.gz' | wc -l)"
if [[ "${extracted_gz}" -ne 265 ]]; then
  echo "Expected 265 extracted training CEL.gz files, found ${extracted_gz}"
  exit 12
fi
find "${CEL_DIR}" -type f -iname '*.CEL.gz' -print0 \
  | xargs -0 -n 1 -P 4 gzip -dk
extracted_cel="$(find "${CEL_DIR}" -type f -iname '*.CEL' | wc -l)"
if [[ "${extracted_cel}" -ne 265 ]]; then
  echo "Expected 265 uncompressed training CEL files, found ${extracted_cel}"
  exit 13
fi

Rscript scripts/09A_stool_feasibility.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/09A_finalize_report.R "${PROJECT_DIR}" "${RUN_ID}"

# Firewall and raw-data invariants are checked after analysis.
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
  printf 'next_action\tCodex independent QC only; do not start model training\n'
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_9A "Stage 9A completed and is ready for Codex QC"
