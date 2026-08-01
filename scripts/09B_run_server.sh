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
RUN_ROOT="${PROJECT_DIR}/logs/09B_model_training"
RESULT_ROOT="${PROJECT_DIR}/results/09B_model_training"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RESULT_DIR="${RESULT_ROOT}/${RUN_ID}"
LOG_FILE="${RUN_ROOT}/stage_9B_${RUN_ID}.log"
STATE_FILE="${RUN_ROOT}/stage_9B_state.tsv"
COMPLETE_MARKER="${RUN_ROOT}/READY_FOR_CODEX_QC"
FAILED_MARKER="${RUN_ROOT}/NEEDS_CODEX_ATTENTION"
SETUP_GIT_COMMIT="${STAGE9B_SETUP_COMMIT:-UNKNOWN}"
MODEL_PATH="${PROJECT_DIR}/objects/locked_stool_model.rds"

mkdir -p "${RUN_ROOT}" "${RESULT_DIR}"
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
    printf 'test_expression_accessed\tFALSE\n'
    printf 'next_action\tCodex should inspect this failure; do not run the test set\n'
  } > "${FAILED_MARKER}"
  cp "${FAILED_MARKER}" "${STATE_FILE}"
  logger -t CRC_STAGE_9B "Stage 9B failed; see ${FAILED_MARKER}"
  exit "${code}"
}
trap on_error ERR

exec > >(tee -a "${LOG_FILE}") 2>&1
cd "${PROJECT_DIR}"
{
  printf 'status\trunning\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'started\t%s\n' "$(date --iso-8601=seconds)"
  printf 'tmux_session\tcrc_stage9B\n'
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'setup_git_commit\t%s\n' "${SETUP_GIT_COMMIT}"
  printf 'test_expression_access\tforbidden\n'
  printf 'monitoring\tserver_background_only\n'
} > "${STATE_FILE}"

test -r AGENTS.md
test -r protocol/stool_model_protocol.md
test -r protocol/stool_model_training_lock_plan.md
test -r reports/stage_9A_stool_feasibility.md
test -r reports/stage_9A_acceptance_audit.md
test -r objects/GSE99573_9A_training_RMA_20260729_213907.rds
test ! -e "${MODEL_PATH}"

# Preflight proves that the only expression object contains the 265 training
# IDs and no test ID. No raw TAR or test CEL is opened in Stage 9B.
Rscript - <<'RS'
project <- "${CRC_PROJECT_ROOT}"
object <- readRDS(file.path(
  project, "objects", "GSE99573_9A_training_RMA_20260729_213907.rds"
))
manifest <- read.delim(
  file.path(project, "metadata", "dataset_manifest.tsv"),
  check.names = FALSE
)
manifest <- manifest[manifest$accession == "GSE99573", ]
training <- manifest$sample_id[manifest$validation_split == "training"]
testing <- manifest$sample_id[manifest$validation_split == "testing"]
stopifnot(
  identical(object$test_expression_accessed, FALSE),
  ncol(object$expression) == 265L,
  setequal(colnames(object$expression), training),
  !any(colnames(object$expression) %in% testing)
)
cat("STAGE9B_TEST_FIREWALL_PREFLIGHT_OK\n")
RS

Rscript scripts/09B_train_lock_model.R \
  "${PROJECT_DIR}" "${RUN_ID}" "${SETUP_GIT_COMMIT}"
Rscript scripts/09B_validate_locked_model.R "${PROJECT_DIR}" "${RUN_ID}"
Rscript scripts/09B_finalize_report.R "${PROJECT_DIR}" "${RUN_ID}"

sha256sum -c <(
  awk '{print $1 "  " $2}' \
    "${RESULT_DIR}/locked_stool_model.sha256"
)
test ! -w "${MODEL_PATH}"

{
  printf 'status\tready_for_qc\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'completed\t%s\n' "$(date --iso-8601=seconds)"
  printf 'log\t%s\n' "${LOG_FILE}"
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_9B_model_training.md"
  printf 'model\t%s\n' "${MODEL_PATH}"
  printf 'model_sha256\t%s\n' "$(awk '{print $1}' "${RESULT_DIR}/locked_stool_model.sha256")"
  printf 'test_expression_accessed\tFALSE\n'
  printf 'git_tag\tpending_desktop_qc_stool-model-locked\n'
  printf 'next_action\tCodex independent QC only; do not run the test set\n'
} > "${COMPLETE_MARKER}"
cp "${COMPLETE_MARKER}" "${STATE_FILE}"
logger -t CRC_STAGE_9B "Stage 9B completed and is ready for Codex QC"
