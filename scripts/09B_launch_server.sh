#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_DIR="${CRC_PROJECT_ROOT}"
SESSION="crc_stage9B"
SETUP_COMMIT="${1:?Usage: 09B_launch_server.sh <setup_git_commit>}"

if tmux has-session -t "${SESSION}" 2>/dev/null; then
  echo "tmux session ${SESSION} already exists"
  exit 2
fi
tmux new-session -d -s "${SESSION}" \
  "cd \"${PROJECT_DIR}\" && STAGE9B_SETUP_COMMIT=\"${SETUP_COMMIT}\" exec bash scripts/09B_run_server.sh"
echo "STAGE9B_SUBMITTED tmux_session=${SESSION} setup_commit=${SETUP_COMMIT}"
