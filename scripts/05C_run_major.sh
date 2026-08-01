#!/usr/bin/env bash
set -uo pipefail

PROJECT_ROOT="${1:-$(pwd)}"
LOG_DIR="${PROJECT_ROOT}/logs/05C_annotation"
RESULT_DIR="${PROJECT_ROOT}/results/05C_annotation"
mkdir -p "${LOG_DIR}" "${RESULT_DIR}"
rm -f \
  "${LOG_DIR}/major_discovery.exit_code" \
  "${LOG_DIR}/major_discovery.finished_at"

STARTED_AT="$(date -Is)"
printf '%s\n' "${STARTED_AT}" >"${LOG_DIR}/major_discovery.started_at"

/usr/bin/time -v \
  -o "${LOG_DIR}/major_discovery.resources.txt" \
  Rscript "${PROJECT_ROOT}/scripts/05C_major_annotation_discovery.R" \
    "${PROJECT_ROOT}" \
  >"${LOG_DIR}/major_discovery.log" 2>&1
EXIT_CODE=$?

printf '%s\n' "${EXIT_CODE}" >"${LOG_DIR}/major_discovery.exit_code"
date -Is >"${LOG_DIR}/major_discovery.finished_at"
exit "${EXIT_CODE}"
