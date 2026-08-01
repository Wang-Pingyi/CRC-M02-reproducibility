#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${CRC_PROJECT_ROOT:-}" ]]; then
  PROJECT_ROOT="${CRC_PROJECT_ROOT}"
elif [[ -d "${CRC_PROJECT_ROOT}" ]]; then
  PROJECT_ROOT="${CRC_PROJECT_ROOT}"
else
  PROJECT_ROOT="${CRC_PROJECT_ROOT}"
fi

PLAN="${PROJECT_ROOT}/metadata/download_plan.tsv"
DATA_RAW="${PROJECT_ROOT}/data_raw"
LOG_DIR="${PROJECT_ROOT}/logs"
REPORT_DIR="${PROJECT_ROOT}/reports"
META_DIR="${PROJECT_ROOT}/metadata"
HEADER_DIR="${LOG_DIR}/stage_4B_headers"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_LOG="${LOG_DIR}/stage_4B_download_${RUN_ID}.log"
FAIL_LOG="${LOG_DIR}/stage_4B_failures_${RUN_ID}.tsv"
FILE_MANIFEST="${META_DIR}/file_manifest.tsv"
REPORT="${REPORT_DIR}/stage_4B_download_integrity.md"
STATUS_FILE="${PROJECT_ROOT}/STATUS.md"
INCLUDE_OPTIONAL_SPATIAL="${STAGE4B_INCLUDE_OPTIONAL_SPATIAL:-1}"
STAGE4B_PROXY="${STAGE4B_PROXY:-}"

if [[ -n "${STAGE4B_PROXY}" ]]; then
  export ALL_PROXY="${STAGE4B_PROXY}"
fi

mkdir -p "${DATA_RAW}" "${LOG_DIR}" "${REPORT_DIR}" "${META_DIR}" "${HEADER_DIR}"
exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Stage 4B run ID: ${RUN_ID}"
echo "Project root: ${PROJECT_ROOT}"
echo "Download plan: ${PLAN}"
echo "Optional spatial enabled: ${INCLUDE_OPTIONAL_SPATIAL}"
echo "Proxy configured: $([[ -n "${STAGE4B_PROXY}" ]] && echo yes || echo no)"

if [[ -n "${STAGE4B_PROXY}" ]]; then
  if ! curl --proxy "${STAGE4B_PROXY}" --fail --location --silent --show-error \
      --connect-timeout 20 --max-time 45 -o /dev/null \
      -w 'Preflight NCBI HTTP %{http_code}, total %{time_total}s\n' \
      'https://www.ncbi.nlm.nih.gov/'; then
    echo "ERROR: proxy preflight to NCBI failed; refusing to start downloads" >&2
    exit 4
  fi
else
  echo "WARNING: no proxy configured; external downloads may fail" >&2
fi

if [[ ! -s "${PLAN}" ]]; then
  echo "ERROR: download plan missing or empty: ${PLAN}" >&2
  exit 2
fi
for command_name in curl sha256sum stat awk sed find flock; do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    echo "ERROR: required command missing: ${command_name}" >&2
    exit 2
  fi
done

exec 9>"${LOG_DIR}/stage_4B_download.lock"
if ! flock -n 9; then
  echo "ERROR: another Stage 4B downloader holds the lock" >&2
  exit 3
fi

# data_raw is immutable between runs. Temporarily restore owner write permission
# only for this controlled downloader.
chmod u+rwx "${DATA_RAW}"
find "${DATA_RAW}" -type d -exec chmod u+rwx {} +
find "${DATA_RAW}" -type f -exec chmod u+rw {} +

printf "accession\trelative_path\tofficial_url\treason\n" >"${FAIL_LOG}"

safe_header_name() {
  printf "%s" "$1" | sed 's#[^A-Za-z0-9._-]#_#g'
}

download_one() {
  local accession="$1"
  local relative_path="$2"
  local official_url="$3"
  local data_class="$4"
  local priority="$5"
  local source_version="$6"
  local destination="${DATA_RAW}/${relative_path}"
  local partial="${destination}.part"
  local header_file="${HEADER_DIR}/$(safe_header_name "${relative_path}").headers.txt"

  if [[ "${priority}" == "optional_spatial" && "${INCLUDE_OPTIONAL_SPATIAL}" != "1" ]]; then
    echo "SKIP optional spatial: ${relative_path}"
    return 0
  fi

  mkdir -p "$(dirname -- "${destination}")"
  if [[ -s "${destination}" ]]; then
    echo "EXISTS: ${relative_path}"
    return 0
  fi

  echo "DOWNLOAD: ${relative_path}"
  if curl \
    --fail \
    --location \
    --http1.1 \
    --continue-at - \
    --retry 8 \
    --retry-delay 10 \
    --retry-max-time 0 \
    --retry-all-errors \
    --connect-timeout 60 \
    --speed-time 300 \
    --speed-limit 1024 \
    --remote-time \
    --user-agent "CRC-carcinogenesis-stage4B/1.0" \
    --dump-header "${header_file}.tmp" \
    --output "${partial}" \
    "${official_url}"; then
    mv -f -- "${header_file}.tmp" "${header_file}"
    mv -f -- "${partial}" "${destination}"
    chmod a-w "${destination}"
    echo "DONE: ${relative_path} ($(stat -c %s "${destination}") bytes)"
  else
    local rc=$?
    rm -f -- "${header_file}.tmp"
    printf "%s\t%s\t%s\tcurl_exit_%s\n" \
      "${accession}" "${relative_path}" "${official_url}" "${rc}" >>"${FAIL_LOG}"
    echo "FAILED: ${relative_path}; curl exit ${rc}" >&2
  fi
}

tail -n +2 "${PLAN}" |
while IFS=$'\t' read -r accession relative_path official_url data_class priority source_version; do
  [[ -n "${accession}" ]] || continue
  download_one \
    "${accession}" \
    "${relative_path}" \
    "${official_url}" \
    "${data_class}" \
    "${priority}" \
    "${source_version}"
done

printf "accession\trelative_path\tfile_name\tdata_class\tpriority\tsource_version\tofficial_url\tdownload_date_utc\tsource_last_modified\tsize_bytes\tsha256\tstatus\n" \
  >"${FILE_MANIFEST}.tmp"

planned=0
verified=0
missing=0
optional_skipped=0
total_bytes=0

tail -n +2 "${PLAN}" |
while IFS=$'\t' read -r accession relative_path official_url data_class priority source_version; do
  [[ -n "${accession}" ]] || continue
  planned=$((planned + 1))
  destination="${DATA_RAW}/${relative_path}"
  header_file="${HEADER_DIR}/$(safe_header_name "${relative_path}").headers.txt"
  file_name="$(basename -- "${relative_path}")"
  if [[ "${priority}" == "optional_spatial" && "${INCLUDE_OPTIONAL_SPATIAL}" != "1" ]]; then
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\toptional_not_requested\n" \
      "${accession}" "${relative_path}" "${file_name}" "${data_class}" "${priority}" \
      "${source_version}" "${official_url}" >>"${FILE_MANIFEST}.tmp"
    optional_skipped=$((optional_skipped + 1))
  elif [[ -s "${destination}" ]]; then
    size_bytes="$(stat -c %s "${destination}")"
    sha256="$(sha256sum "${destination}" | awk '{print $1}')"
    download_date_utc="$(date -u -r "${destination}" +%Y-%m-%dT%H:%M:%SZ)"
    source_last_modified="NA"
    if [[ -s "${header_file}" ]]; then
      source_last_modified="$(
        awk 'BEGIN{IGNORECASE=1} /^last-modified:/ {$1=""; sub(/^ /,""); gsub(/\r/,""); value=$0} END{if(value!="") print value}' \
          "${header_file}"
      )"
      [[ -n "${source_last_modified}" ]] || source_last_modified="NA"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tverified\n" \
      "${accession}" "${relative_path}" "${file_name}" "${data_class}" "${priority}" \
      "${source_version}" "${official_url}" "${download_date_utc}" \
      "${source_last_modified}" "${size_bytes}" "${sha256}" >>"${FILE_MANIFEST}.tmp"
    verified=$((verified + 1))
    total_bytes=$((total_bytes + size_bytes))
  else
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tmissing_or_failed\n" \
      "${accession}" "${relative_path}" "${file_name}" "${data_class}" "${priority}" \
      "${source_version}" "${official_url}" >>"${FILE_MANIFEST}.tmp"
    missing=$((missing + 1))
  fi
done

# The counters above were updated in a pipeline subshell; derive final values
# from the completed manifest instead.
mv -f -- "${FILE_MANIFEST}.tmp" "${FILE_MANIFEST}"
planned="$(tail -n +2 "${FILE_MANIFEST}" | wc -l | awk '{print $1}')"
verified="$(awk -F '\t' 'NR>1 && $12=="verified"{n++} END{print n+0}' "${FILE_MANIFEST}")"
missing="$(awk -F '\t' 'NR>1 && $12=="missing_or_failed"{n++} END{print n+0}' "${FILE_MANIFEST}")"
optional_skipped="$(awk -F '\t' 'NR>1 && $12=="optional_not_requested"{n++} END{print n+0}' "${FILE_MANIFEST}")"
total_bytes="$(awk -F '\t' 'NR>1 && $12=="verified"{sum+=$10} END{printf "%.0f",sum+0}' "${FILE_MANIFEST}")"
failure_rows="$(tail -n +2 "${FAIL_LOG}" | wc -l | awk '{print $1}')"

cat >"${REPORT}.tmp" <<EOF
# Stage 4B download integrity report

## Run

- Run ID: ${RUN_ID}
- Project root: ${PROJECT_ROOT}
- Started/submitted under tmux; run log: ${RUN_LOG}
- Download plan: ${PLAN}
- FASTQ policy: not downloaded
- Optional GSE226997 spatial bundles enabled: ${INCLUDE_OPTIONAL_SPATIAL}
- SSH relay proxy used: $([[ -n "${STAGE4B_PROXY}" ]] && echo yes || echo no)

## Integrity summary

- Planned files: ${planned}
- SHA256-verified files: ${verified}
- Missing or failed files: ${missing}
- Optional files skipped by policy: ${optional_skipped}
- Bytes verified: ${total_bytes}
- Curl failure records in this run: ${failure_rows}
- File manifest: ${FILE_MANIFEST}
- Failure log: ${FAIL_LOG}

## Immutability

Original archives were not unpacked or overwritten. Completed files have write
permission removed. Processed copies, when later authorized, must be written to
\`data_processed/\`.

## Completion decision

EOF

if [[ "${missing}" -eq 0 && "${failure_rows}" -eq 0 ]]; then
  cat >>"${REPORT}.tmp" <<EOF
Stage 4B download and SHA256 verification completed successfully.
EOF
  find "${DATA_RAW}" -type f -exec chmod a-w {} +
  find "${DATA_RAW}" -depth -type d -exec chmod a-w {} +
  completion_status="complete"
  exit_code=0
else
  cat >>"${REPORT}.tmp" <<EOF
Stage 4B is incomplete. Re-run the same script to resume partial downloads and
retry failed files. The raw-data tree remains owner-writable only so the
controlled downloader can resume; completed individual files are read-only.
EOF
  completion_status="incomplete_retry_required"
  exit_code=1
fi
mv -f -- "${REPORT}.tmp" "${REPORT}"

cat >"${STATUS_FILE}.tmp" <<EOF
# CRC Carcinogenesis Remote Status

- Stage: 4B — public-data download and integrity validation
- Status: ${completion_status}
- Run ID: ${RUN_ID}
- Project root: ${PROJECT_ROOT}
- Planned files: ${planned}
- SHA256-verified files: ${verified}
- Missing or failed files: ${missing}
- Optional files skipped: ${optional_skipped}
- File manifest: ${FILE_MANIFEST}
- Integrity report: ${REPORT}
- Run log: ${RUN_LOG}
- Next stage authorization: not granted
EOF
mv -f -- "${STATUS_FILE}.tmp" "${STATUS_FILE}"

echo "Stage 4B status: ${completion_status}"
echo "Verified: ${verified}/${planned}; missing: ${missing}"
exit "${exit_code}"
