#!/usr/bin/env bash

set -uo pipefail

project_dir="${1:-$(pwd)}"
run_name="05B_biological_qc_revision"
log_dir="${project_dir}/logs/${run_name}"
mkdir -p "${log_dir}"
export R_LIBS_USER="${project_dir}/environment/R-library${R_LIBS_USER:+:${R_LIBS_USER}}"

started_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
/usr/bin/time -v Rscript \
  "${project_dir}/scripts/05B_full_qc_integration.R" \
  "${project_dir}" \
  >"${log_dir}/run.out" \
  2>"${log_dir}/run.err"
status=$?
finished_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

{
  printf 'exit_code=%s\n' "${status}"
  printf 'started_utc=%s\n' "${started_utc}"
  printf 'finished_utc=%s\n' "${finished_utc}"
} >"${log_dir}/run.status"

exit "${status}"
