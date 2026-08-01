#!/usr/bin/env bash
set -uo pipefail

PROJECT_ROOT="${1:-$(pwd)}"
cd "${PROJECT_ROOT}"

mkdir -p logs/05C_annotation
rm -f logs/05C_prepare_copykat.exit logs/05C_prepare_copykat.finished

/usr/bin/time -v \
  -o logs/05C_annotation/prepare_copykat.resources.txt \
  Rscript scripts/05C_prepare_copykat_inputs.R . \
  > logs/05C_prepare_copykat.log 2>&1
rc=$?

printf '%s\n' "${rc}" > logs/05C_prepare_copykat.exit
date -Is > logs/05C_prepare_copykat.finished
exit "${rc}"
