#!/usr/bin/env bash
set -Eeuo pipefail
cd "${CRC_PROJECT_ROOT}"
exec bash scripts/08A_run_server.sh
