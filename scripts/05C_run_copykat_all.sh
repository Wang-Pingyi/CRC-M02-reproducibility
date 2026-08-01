#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${1:-$(pwd)}"
PARAM_FILE="${PROJECT_ROOT}/config/annotation_parameters.tsv"
PLAN_FILE="${PROJECT_ROOT}/results/05C_annotation/copykat_sample_plan.tsv"
LOG_DIR="${PROJECT_ROOT}/logs/05C_copykat"
R_LIBRARY="${PROJECT_ROOT}/environment/R_library"
mkdir -p "${LOG_DIR}" "${R_LIBRARY}"
export R_LIBS_USER="${R_LIBRARY}"

EXIT_FILE="${LOG_DIR}/pipeline.exit"
FINISHED_FILE="${LOG_DIR}/pipeline.finished"
DRIVER_LOG="${LOG_DIR}/pipeline_driver.log"
rm -f "${EXIT_FILE}" "${FINISHED_FILE}"
: > "${DRIVER_LOG}"
exec > >(tee -a "${DRIVER_LOG}") 2>&1
record_exit() {
  local rc=$?
  trap - EXIT
  printf '%s\n' "${rc}" > "${EXIT_FILE}"
  date -Is > "${FINISHED_FILE}"
  exit "${rc}"
}
trap record_exit EXIT

Rscript "${PROJECT_ROOT}/scripts/05C_install_copykat.R" "${PROJECT_ROOT}" \
  >"${LOG_DIR}/install_copykat.log" 2>&1

PARALLEL_SAMPLES="$(
  awk -F '\t' '$1=="cnv" && $2=="parallel_samples" {print $3}' \
    "${PARAM_FILE}"
)"
if [[ -z "${PARALLEL_SAMPLES}" ]]; then
  echo "Missing cnv/parallel_samples parameter" >&2
  exit 1
fi

mapfile -t SAMPLE_KEYS < <(
  awk -F '\t' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        header_value=$i
        gsub(/"/, "", header_value)
        if (header_value=="sample_key") key_col=i
        if (header_value=="status") status_col=i
      }
      next
    }
    {
      status_value=$status_col
      gsub(/"/, "", status_value)
      if (status_value=="ready") {
        gsub(/"/, "", $key_col)
        print $key_col
      }
    }
  ' "${PLAN_FILE}"
)

if [[ "${#SAMPLE_KEYS[@]}" -eq 0 ]]; then
  echo "No ready CopyKAT inputs found" >&2
  exit 1
fi

run_one() {
  local sample_key="$1"
  local status_file="${PROJECT_ROOT}/cache/05C_copykat_outputs/${sample_key}/run_status.tsv"
  if [[ -s "${status_file}" ]] && awk -F '\t' '
    NR>1 {
      status_value=$2
      gsub(/"/, "", status_value)
      if (status_value=="completed") found=1
    }
    END {exit !found}
  ' "${status_file}"; then
    echo "${sample_key}: already completed"
    return 0
  fi
  Rscript "${PROJECT_ROOT}/scripts/05C_run_copykat_sample.R" \
    "${PROJECT_ROOT}" "${sample_key}" \
    >"${LOG_DIR}/${sample_key}.log" 2>&1
}
export -f run_one
export PROJECT_ROOT LOG_DIR

printf '%s\n' "${SAMPLE_KEYS[@]}" |
  xargs -r -P "${PARALLEL_SAMPLES}" -I '{}' \
    bash -c 'run_one "$@"' _ '{}'

Rscript "${PROJECT_ROOT}/scripts/05C_summarize_copykat.R" "${PROJECT_ROOT}" \
  >"${LOG_DIR}/summarize_copykat.log" 2>&1
