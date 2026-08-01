#!/usr/bin/env bash
set -uo pipefail

PROJECT_ROOT="${1:-$(pwd)}"
cd "${PROJECT_ROOT}"

LOG_DIR="logs/05C_annotation"
mkdir -p "${LOG_DIR}"
rm -f "${LOG_DIR}/validation.exit" "${LOG_DIR}/validation.finished"

/usr/bin/time -v \
  -o "${LOG_DIR}/validation.resources.txt" \
  Rscript scripts/05C_validate_outputs.R . \
  > "${LOG_DIR}/validation.log" 2>&1
rc=$?

printf '%s\n' "${rc}" > "${LOG_DIR}/validation.exit"
date -Is > "${LOG_DIR}/validation.finished"
exit "${rc}"
