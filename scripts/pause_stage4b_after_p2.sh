#!/usr/bin/env bash
set -uo pipefail

DOWNLOAD_SESSION="${1:-crc_4B_download_20260724T210501Z}"
PROJECT_ROOT="${2:-${CRC_PROJECT_ROOT}}"
MAX_WAIT_SECONDS="${3:-2700}"
TARGET="${PROJECT_ROOT}/data_raw/GSE226997/GSM7089856_Ajou_Visium_P2.tar.gz"
LOG="${PROJECT_ROOT}/logs/stage_4B_pause_guard_20260725.log"

start="$(date +%s)"
deadline="$((start + MAX_WAIT_SECONDS))"
reason="max_wait_reached"

while tmux has-session -t "${DOWNLOAD_SESSION}" 2>/dev/null; do
  if [[ -s "${TARGET}" ]]; then
    reason="P2_completed"
    break
  fi
  if [[ "$(date +%s)" -ge "${deadline}" ]]; then
    break
  fi
  sleep 10
done

if tmux has-session -t "${DOWNLOAD_SESSION}" 2>/dev/null; then
  tmux kill-session -t "${DOWNLOAD_SESSION}"
fi

{
  printf "paused_at_utc=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf "reason=%s\n" "${reason}"
  printf "completed_files=%s\n" \
    "$(find "${PROJECT_ROOT}/data_raw" -type f ! -name "*.part" | wc -l)"
  printf "partial_files=%s\n" \
    "$(find "${PROJECT_ROOT}/data_raw" -type f -name "*.part" | wc -l)"
} >"${LOG}"
