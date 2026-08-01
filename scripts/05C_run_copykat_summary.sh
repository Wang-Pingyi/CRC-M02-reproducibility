#!/usr/bin/env bash
set -uo pipefail

PROJECT_ROOT="${1:-$(pwd)}"
cd "${PROJECT_ROOT}"

LOG_DIR="logs/05C_copykat"
mkdir -p "${LOG_DIR}"
rm -f "${LOG_DIR}/summary.exit" "${LOG_DIR}/summary.finished"

/usr/bin/time -v \
  -o "${LOG_DIR}/summary.resources.txt" \
  Rscript scripts/05C_summarize_copykat.R . \
  > "${LOG_DIR}/summarize_copykat.log" 2>&1
rc=$?

printf '%s\n' "${rc}" > "${LOG_DIR}/summary.exit"
date -Is > "${LOG_DIR}/summary.finished"
exit "${rc}"
