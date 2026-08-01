#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_DIR="${CRC_PROJECT_ROOT}"
SESSION="crc_stage9C"

cd "${PROJECT_DIR}"
if tmux has-session -t "${SESSION}" 2>/dev/null; then
  printf 'Refusing to replace active tmux session: %s\n' "${SESSION}" >&2
  exit 1
fi
if [[ -e logs/09C_external_test/ONE_TIME_TEST_STARTED ]]; then
  printf 'Refusing a second one-time test: start marker already exists\n' >&2
  exit 1
fi
tmux new-session -d -s "${SESSION}" \
  "cd '${PROJECT_DIR}' && STAGE9C_SETUP_COMMIT='${STAGE9C_SETUP_COMMIT:-UNKNOWN}' bash scripts/09C_run_server.sh"
printf 'Submitted tmux session: %s\n' "${SESSION}"
