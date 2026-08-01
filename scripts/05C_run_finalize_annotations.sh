#!/usr/bin/env bash
set -uo pipefail

PROJECT_ROOT="${1:-$(pwd)}"
LOG_DIR="${PROJECT_ROOT}/logs/05C_annotation"
mkdir -p "${LOG_DIR}"
rm -f \
  "${LOG_DIR}/finalize_annotations.exit_code" \
  "${LOG_DIR}/finalize_annotations.finished_at"

date -Is >"${LOG_DIR}/finalize_annotations.started_at"
/usr/bin/time -v \
  -o "${LOG_DIR}/finalize_annotations.resources.txt" \
  Rscript "${PROJECT_ROOT}/scripts/05C_finalize_annotations.R" \
    "${PROJECT_ROOT}" \
  >"${LOG_DIR}/finalize_annotations.log" 2>&1
EXIT_CODE=$?

printf '%s\n' "${EXIT_CODE}" \
  >"${LOG_DIR}/finalize_annotations.exit_code"
date -Is >"${LOG_DIR}/finalize_annotations.finished_at"
exit "${EXIT_CODE}"
