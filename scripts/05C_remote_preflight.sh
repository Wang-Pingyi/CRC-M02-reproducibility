#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:?Project root is required}"

echo "ROOT"
printf '%s\n' "${PROJECT_ROOT}"
if [[ -d "${PROJECT_ROOT}" ]]; then
  echo "ROOT_EXISTS=yes"
else
  echo "ROOT_EXISTS=no"
  exit 1
fi

find "${PROJECT_ROOT}" -maxdepth 2 -type f \
  \( -name AGENTS.md \
     -o -name stage_5B_full_qc_integration.md \
     -o -name GSE201348_harmony_integrated.rds \) \
  -printf '%p\t%s bytes\n'

echo "SOFTWARE"
Rscript --version 2>&1
Rscript -e '
  cat(
    "Seurat=", as.character(packageVersion("Seurat")),
    "\nHarmony=", as.character(packageVersion("harmony")),
    "\n",
    sep = ""
  )
'

echo "RESOURCE"
free -h
df -h "${PROJECT_ROOT}"

echo "PROCESSES"
pgrep -af '05C_|Rscript' || true
tmux ls 2>/dev/null || true
