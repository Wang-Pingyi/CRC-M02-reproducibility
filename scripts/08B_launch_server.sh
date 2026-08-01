#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT="${CRC_PROJECT_ROOT}"
SESSION="crc_stage8B"
SETUP_COMMIT="${1:-not_recorded}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
SUBMISSION="$PROJECT/logs/08B_bulk_validation/SUBMISSION.tsv"
mkdir -p "$(dirname "$SUBMISSION")"

if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "Refusing to replace active tmux session: $SESSION" >&2
  exit 1
fi

printf -v command '%q ' bash "$PROJECT/scripts/08B_run_server.sh" "$PROJECT" "$RUN_ID" "$SETUP_COMMIT"
tmux new-session -d -s "$SESSION" "$command"
printf 'run_id\ttmux_session\tsetup_commit\tsubmitted_at\n%s\t%s\t%s\t%s\n' \
  "$RUN_ID" "$SESSION" "$SETUP_COMMIT" "$(date --iso-8601=seconds)" > "$SUBMISSION"
printf 'RUN_ID=%s\nTMUX_SESSION=%s\n' "$RUN_ID" "$SESSION"
