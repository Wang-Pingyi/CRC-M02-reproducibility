#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_DIR="${CRC_PROJECT_ROOT}"
SESSION="crc_stage9A"

if tmux has-session -t "${SESSION}" 2>/dev/null; then
  echo "tmux session ${SESSION} already exists"
  exit 2
fi

tmux new-session -d -s "${SESSION}" \
  "cd \"${PROJECT_DIR}\" && exec bash scripts/09A_run_server.sh"
echo "STAGE9A_SUBMITTED tmux_session=${SESSION}"
