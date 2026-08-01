#!/usr/bin/env bash
set -Eeuo pipefail
umask 002

PROJECT_DIR="${CRC_PROJECT_ROOT}"
RUN_ID="${1:-20260729_105540}"
RESULT_DIR="${PROJECT_DIR}/results/08A_bulk_preprocessing/${RUN_ID}"
LOG_DIR="${PROJECT_DIR}/logs/08A_bulk_preprocessing"
CHECKS="${RESULT_DIR}/validation_checks.tsv"
MARKER="${LOG_DIR}/STAGE_8A_ACCEPTED"
STATE="${LOG_DIR}/stage_8A_state.tsv"

[[ -s "${CHECKS}" ]]
awk -F '\t' 'NR > 1 && $2 != "TRUE" { exit 1 }' "${CHECKS}"
grep -Fq 'Stage 8A is technically and biologically accepted.' \
  "${PROJECT_DIR}/reports/stage_8A_bulk_preprocessing.md"
grep -Fq 'Stage 8B remains unauthorized' \
  "${PROJECT_DIR}/reports/stage_8A_acceptance_audit.md"

{
  printf 'status\taccepted\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'accepted\t%s\n' "$(date --iso-8601=seconds)"
  printf 'validation\t33/33 passed\n'
  printf 'report\t%s\n' "${PROJECT_DIR}/reports/stage_8A_bulk_preprocessing.md"
  printf 'acceptance_audit\t%s\n' "${PROJECT_DIR}/reports/stage_8A_acceptance_audit.md"
  printf 'instruction\tStop after Stage 8A; Stage 8B is not authorized\n'
} > "${MARKER}"
cp "${MARKER}" "${STATE}"
logger -t CRC_STAGE_8A "Stage 8A accepted; stopped before Stage 8B"
