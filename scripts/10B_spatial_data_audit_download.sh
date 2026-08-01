#!/usr/bin/env bash
# Stage 10B: metadata audit and immutable processed-data download only.
# This script deliberately does not unpack matrices, run Space Ranger, score
# modules, or perform any spatial/single-cell biological analysis.
set -Eeuo pipefail

STAGE="stage10B"
MODE=""
RUN_ID=""
PROJECT_ROOT=""
SOURCES_FILE=""
JOB_ID_OR_PID="${JOB_ID_OR_PID:-not_submitted}"
LAUNCH_COMMAND="${LAUNCH_COMMAND:-not_recorded}"
GIT_COMMIT="${GIT_COMMIT:-not_recorded}"
RESUME_FAILED=0
RECOVERY_ID=""
RECOVERY_EVENT=""

usage() {
  cat <<'USAGE'
Usage:
  10B_spatial_data_audit_download.sh --smoke --root PATH --run-id ID --sources FILE
  10B_spatial_data_audit_download.sh --full  --root PATH --run-id ID --sources FILE
  10B_spatial_data_audit_download.sh --full --resume-failed --root PATH --run-id ID --sources FILE

The smoke test fetches only metadata plus a one-mebibyte HTTP range request.
The full mode downloads only listed processed archives/matrices and metadata.
It never downloads raw sequencing .tar files or FASTQ files.
USAGE
}

while (($#)); do
  case "$1" in
    --smoke) MODE="smoke" ;;
    --full) MODE="full" ;;
    --resume-failed) RESUME_FAILED=1 ;;
    --root) PROJECT_ROOT="$2"; shift ;;
    --run-id) RUN_ID="$2"; shift ;;
    --sources) SOURCES_FILE="$2"; shift ;;
    --help|-h) usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

[[ "$MODE" == "smoke" || "$MODE" == "full" ]] || { usage >&2; exit 2; }
[[ -n "$PROJECT_ROOT" && -n "$RUN_ID" && -n "$SOURCES_FILE" ]] || { usage >&2; exit 2; }
[[ -d "$PROJECT_ROOT" && -r "$SOURCES_FILE" ]] || { printf 'Missing root or source plan.\n' >&2; exit 2; }

META_DIR="$PROJECT_ROOT/metadata/stage10B/$RUN_ID"
RAW_DIR="$PROJECT_ROOT/data_raw/stage10B_$RUN_ID"
LOG_DIR="$PROJECT_ROOT/logs/10B_data_audit/$RUN_ID"
STATUS_DIR="$PROJECT_ROOT/status"
REPORTS_DIR="$PROJECT_ROOT/reports"
SMOKE_DIR="$META_DIR/smoke_test"
DOWNLOAD_MANIFEST="$META_DIR/DOWNLOAD_MANIFEST.tsv"

mkdir -p "$META_DIR" "$LOG_DIR" "$STATUS_DIR" "$REPORTS_DIR"

timestamp() { date --iso-8601=seconds; }
trim() { local x="$1"; x="${x#${x%%[![:space:]]*}}"; x="${x%${x##*[![:space:]]}}"; printf '%s' "$x"; }
atomic_text() {
  local target="$1"; shift
  local tmp="${target}.tmp.$$"
  printf '%s\n' "$@" > "$tmp"
  mv -f "$tmp" "$target"
}

raw_unlocked=0
lock_raw_dir() {
  if [[ -d "$RAW_DIR" ]]; then
    chmod -R a-w "$RAW_DIR" || true
  fi
  chmod a-w "$PROJECT_ROOT/data_raw" || true
  raw_unlocked=0
}
unlock_stage_raw_dir() {
  # Only the new Stage 10B directory is made writable. Existing raw data are
  # never modified. The parent is immediately returned to read-only mode.
  chmod u+w "$PROJECT_ROOT/data_raw"
  mkdir -p "$RAW_DIR"
  chmod -R u+w "$RAW_DIR"
  chmod a-w "$PROJECT_ROOT/data_raw"
  raw_unlocked=1
}
write_failure_marker() {
  local code="$1" line="$2"
  [[ "$MODE" == "full" ]] || return 0
  lock_raw_dir
  atomic_text "$STATUS_DIR/$STAGE.FAILED" \
    "stage=$STAGE" "run_id=$RUN_ID" "exit_code=$code" "line=$line" \
    "time=$(timestamp)" "log=$LOG_DIR/full.log" \
    "preserved_outputs=$META_DIR;$RAW_DIR" \
    "recovery=JOB_ID_OR_PID=RETRY GIT_COMMIT=$GIT_COMMIT bash scripts/10B_spatial_data_audit_download.sh --full --resume-failed --root '$PROJECT_ROOT' --run-id '$RUN_ID' --sources '$SOURCES_FILE'"
  rm -f "$STATUS_DIR/$STAGE.RUNNING"
}
on_error() {
  local code="$?" line="$1"
  set +e
  printf '[%s] ERROR line=%s exit=%s\n' "$(timestamp)" "$line" "$code" >> "$LOG_DIR/${MODE}.log"
  write_failure_marker "$code" "$line"
  exit "$code"
}
trap 'on_error $LINENO' ERR

for tool in curl sha256sum md5sum unzip zipinfo python3 gzip tar awk; do
  command -v "$tool" >/dev/null 2>&1 || { printf 'Missing required tool: %s\n' "$tool" >&2; exit 127; }
done

target_path() {
  local dataset="$1" filename="$2"
  if [[ "$MODE" == "smoke" ]]; then
    # A smoke test must not unlock or populate the immutable raw-data tree.
    printf '%s/downloads/%s/%s' "$SMOKE_DIR" "$dataset" "$filename"
    return
  fi
  case "$dataset" in
    PMC11121171|Cell_Reports_111929) printf '%s/articles/%s' "$META_DIR" "$filename" ;;
    *) printf '%s/%s/%s' "$RAW_DIR" "$dataset" "$filename" ;;
  esac
}

download_one() {
  local dataset="$1" asset_type="$2" filename="$3" url="$4" expected_bytes="$5" expected_md5="$6"
  local target part observed_md5 observed_bytes
  local -a timeout_args=()
  [[ "$expected_bytes" == "NA" ]] && expected_bytes=""
  [[ "$expected_md5" == "NA" ]] && expected_md5=""
  [[ "$MODE" == "smoke" ]] && timeout_args=(--max-time 60)
  target="$(target_path "$dataset" "$filename")"
  mkdir -p "$(dirname "$target")"
  if [[ -s "$target" ]]; then
    observed_bytes="$(stat -c '%s' "$target")"
    if [[ -n "$expected_md5" ]]; then
      observed_md5="$(md5sum "$target" | awk '{print tolower($1)}')"
      [[ "$observed_md5" == "${expected_md5,,}" ]] && return 0
    elif [[ -z "$expected_bytes" || "$observed_bytes" == "$expected_bytes" ]]; then
      return 0
    fi
  fi
  part="${target}.part"
  if [[ -s "$part" ]]; then
    curl --fail --location --retry 5 --retry-all-errors --connect-timeout 30 "${timeout_args[@]}" --continue-at - --output "$part" "$url"
  else
    curl --fail --location --retry 5 --retry-all-errors --connect-timeout 30 "${timeout_args[@]}" --output "$part" "$url"
  fi
  observed_bytes="$(stat -c '%s' "$part")"
  [[ -z "$expected_bytes" || "$observed_bytes" == "$expected_bytes" ]] || { printf 'Unexpected size for %s: %s != %s\n' "$filename" "$observed_bytes" "$expected_bytes" >&2; return 1; }
  if [[ -n "$expected_md5" ]]; then
    observed_md5="$(md5sum "$part" | awk '{print tolower($1)}')"
    [[ "$observed_md5" == "${expected_md5,,}" ]] || { printf 'MD5 mismatch for %s\n' "$filename" >&2; return 1; }
  fi
  mv -f "$part" "$target"
}

fetch_metadata_only() {
  local dataset priority asset_type filename url expected_bytes expected_md5 required notes
  while IFS=$'\t' read -r dataset priority asset_type filename url expected_bytes expected_md5 required notes; do
    [[ "$dataset" == "dataset_id" || -z "$dataset" ]] && continue
    [[ "$MODE" == "smoke" && !( "$dataset" == "E-GEAD-622" && "$asset_type" =~ ^(filelist|idf|sdrf)$ ) ]] && continue
    [[ "$expected_bytes" == "NA" ]] && expected_bytes=""
    [[ "$expected_md5" == "NA" ]] && expected_md5=""
    required="$(trim "$required")"
    case "$asset_type" in
      filelist|idf|sdrf|series_metadata|article_metadata)
        download_one "$dataset" "$asset_type" "$filename" "$url" "$expected_bytes" "$expected_md5"
        ;;
    esac
  done < "$SOURCES_FILE"
}

download_all_planned_assets() {
  local dataset priority asset_type filename url expected_bytes expected_md5 required notes
  while IFS=$'\t' read -r dataset priority asset_type filename url expected_bytes expected_md5 required notes; do
    [[ "$dataset" == "dataset_id" || -z "$dataset" ]] && continue
    download_one "$dataset" "$asset_type" "$filename" "$url" "$expected_bytes" "$expected_md5"
  done < "$SOURCES_FILE"
}

write_manifest() {
  local tmp="$DOWNLOAD_MANIFEST.tmp.$$"
  printf 'dataset_id\tpriority\tasset_type\tpath\tofficial_url\texpected_bytes\texpected_md5\tobserved_bytes\tobserved_md5\tsha256\tstatus\n' > "$tmp"
  local dataset priority asset_type filename url expected_bytes expected_md5 required notes target obs_bytes obs_md5 sha status
  while IFS=$'\t' read -r dataset priority asset_type filename url expected_bytes expected_md5 required notes; do
    [[ "$dataset" == "dataset_id" || -z "$dataset" ]] && continue
    [[ "$expected_bytes" == "NA" ]] && expected_bytes=""
    [[ "$expected_md5" == "NA" ]] && expected_md5=""
    target="$(target_path "$dataset" "$filename")"
    status="missing"; obs_bytes="NA"; obs_md5="NA"; sha="NA"
    if [[ -s "$target" ]]; then
      obs_bytes="$(stat -c '%s' "$target")"
      obs_md5="$(md5sum "$target" | awk '{print tolower($1)}')"
      sha="$(sha256sum "$target" | awk '{print $1}')"
      status="downloaded"
      if [[ -n "$expected_bytes" && "$obs_bytes" != "$expected_bytes" ]]; then status="size_mismatch"; fi
      if [[ -n "$expected_md5" && "$obs_md5" != "${expected_md5,,}" ]]; then status="md5_mismatch"; fi
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$dataset" "$priority" "$asset_type" "${target#$PROJECT_ROOT/}" "$url" \
      "${expected_bytes:-NA}" "${expected_md5:-NA}" "$obs_bytes" "$obs_md5" "$sha" "$status" >> "$tmp"
  done < "$SOURCES_FILE"

  # Existing GSE226997 bundles are reused. The Stage 4B SHA256 values are
  # copied as audited manifest evidence; these large archives are not copied
  # or re-downloaded in Stage 10B.
  python3 - "$PROJECT_ROOT" >> "$tmp" <<'PY'
import csv, os, sys
root = sys.argv[1]
manifest = os.path.join(root, 'metadata', 'file_manifest.tsv')
for row in csv.DictReader(open(manifest, encoding='utf-8'), delimiter='\t'):
    data_class = row.get('data_class') or row.get('file_type')
    official_url = row.get('official_url') or row.get('source_url') or 'NA'
    if row.get('accession') == 'GSE226997' and data_class in {'processed_spatial_bundle', 'series_metadata'}:
        path = os.path.join(root, 'data_raw', row['relative_path'])
        exists = os.path.isfile(path) and os.path.getsize(path) == int(row['size_bytes'])
        status = 'reused_existing_size_matches_stage4b_manifest' if exists else 'missing_or_size_changed'
        print('\t'.join([
            'GSE226997', 'secondary_existing', data_class, row['relative_path'], official_url,
            row['size_bytes'], 'NA', str(os.path.getsize(path)) if os.path.isfile(path) else 'NA', 'NA',
            row['sha256'], status
        ]))
PY
  mv -f "$tmp" "$DOWNLOAD_MANIFEST"
  sha256sum "$DOWNLOAD_MANIFEST" | awk '{print $1 "  DOWNLOAD_MANIFEST.tsv"}' > "$META_DIR/DOWNLOAD_MANIFEST.sha256"
}

discover_cellreports_supplements() {
  # Publisher-generated supplement URLs are discovered at run time rather than
  # hard-coded. An unavailable publisher asset is recorded but never invalidates
  # the independently archived public matrices.
  local articles="$META_DIR/articles"
  local landing="$articles/Cell_Reports_111929.landing.html"
  local discovered="$articles/Cell_Reports_111929.discovered_urls.tsv"
  local manifest="$articles/ARTICLE_SUPPLEMENT_MANIFEST.tsv"
  mkdir -p "$articles"
  printf 'asset\tofficial_url\tobserved_bytes\tsha256\tstatus\n' > "$manifest"
  if ! curl --fail --location --retry 3 --retry-all-errors --connect-timeout 30 --max-time 120 \
      --output "$landing.part" 'https://www.cell.com/cell-reports/fulltext/S2211-1247(22)01830-7'; then
    rm -f "$landing.part"
    printf 'publisher_landing_page\thttps://www.cell.com/cell-reports/fulltext/S2211-1247(22)01830-7\tNA\tNA\tnot_retrieved\n' >> "$manifest"
    return 0
  fi
  mv -f "$landing.part" "$landing"
  python3 - "$landing" "$discovered" <<'PY'
import html, re, sys
from pathlib import Path
from urllib.parse import urljoin
src, out = map(Path, sys.argv[1:])
text = src.read_text(encoding='utf-8', errors='replace')
hrefs = re.findall(r'''href=["']([^"']+)["']''', text, flags=re.I)
urls = []
for href in hrefs:
    url = urljoin('https://www.cell.com', html.unescape(href))
    low = url.lower()
    if ('attachment' in low or re.search(r'/mmc\\d+', low)) and re.search(r'\\.(pdf|xlsx?|csv|tsv|docx?|zip|gz)(?:[?#]|$)', low):
        urls.append(url)
with out.open('w', encoding='utf-8') as h:
    h.write('asset\tofficial_url\n')
    for i, url in enumerate(dict.fromkeys(urls), 1):
        h.write(f'cellreports_supplement_{i}\t{url}\n')
PY
  local asset url safe target bytes sha
  while IFS=$'\t' read -r asset url; do
    [[ "$asset" == "asset" || -z "$asset" ]] && continue
    safe="${asset//[^A-Za-z0-9_.-]/_}"
    target="$articles/${safe}.bin"
    if curl --fail --location --retry 3 --retry-all-errors --connect-timeout 30 --max-time 300 --output "$target.part" "$url"; then
      mv -f "$target.part" "$target"
      bytes="$(stat -c '%s' "$target")"; sha="$(sha256sum "$target" | awk '{print $1}')"
      printf '%s\t%s\t%s\t%s\tdownloaded\n' "$asset" "$url" "$bytes" "$sha" >> "$manifest"
    else
      rm -f "$target.part"
      printf '%s\t%s\tNA\tNA\tnot_retrieved\n' "$asset" "$url" >> "$manifest"
    fi
  done < "$discovered"
}

inspect_archives_and_write_audit() {
  python3 - "$PROJECT_ROOT" "$RAW_DIR" "$META_DIR" "$RUN_ID" <<'PY'
import csv, gzip, os, re, subprocess, sys, tarfile, zipfile
from pathlib import Path

root, raw, meta, run_id = map(Path, sys.argv[1:])
audit_path = meta / 'DATASET_AUDIT.tsv'
columns = [
    'dataset_id','priority','sample_or_section_id','donor_id','section_count_role',
    'pathology_or_stage','pathology_region','anatomy','material','platform',
    'reference_build','gene_id_status','processed_matrix_status','spot_coordinates_status',
    'he_image_status','patient_overlap_status','inclusion_status','evidence_source','notes'
]
rows = []

def add(**kwargs):
    rows.append({c: kwargs.get(c, 'NA') for c in columns})

def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as h:
        return list(csv.DictReader(h, delimiter='\t'))

def zip_status(dataset):
    p = raw / dataset / f'{dataset}.processed.zip'
    if not p.is_file():
        return ('missing', 'missing', 'missing', 'not_checked')
    try:
        names = zipfile.ZipFile(p).namelist()
    except Exception as exc:
        return ('unreadable', 'unreadable', 'unreadable', f'zip_error:{type(exc).__name__}')
    low = '\n'.join(x.lower() for x in names)
    matrix = 'present' if any(x in low for x in ('filtered_feature_bc_matrix','raw_feature_bc_matrix','matrix.mtx')) else 'not_found_in_archive'
    coords = 'present' if any(x in low for x in ('tissue_positions','tissue_positions_list','spatial/')) else 'not_found_in_archive'
    image = 'present' if any(x.endswith(s) for x in low.split('\n') for s in ('.jpg','.jpeg','.png','.tif','.tiff')) else 'not_found_in_archive'
    gene = 'feature_id_unchecked; inspect matrix feature file at later smoke test'
    return (matrix, coords, image, gene)

def textval(row, key):
    return (row.get(key) or 'NA').strip() or 'NA'

for dataset, priority, material, platform, reference in [
    ('E-GEAD-622', 'primary_early', 'FFPE', '10x Genomics Visium FFPE; MGI DNBSEQ-G400', 'GRCh38 build 2020-A'),
    ('E-GEAD-619', 'primary_early_and_secondary_advanced', 'fresh_frozen', '10x Genomics Visium; Illumina NovaSeq 6000', 'GRCh38 GCF_000001405.39'),
    ('E-GEAD-579', 'secondary', 'fresh_frozen', '10x Genomics Visium; Illumina NovaSeq 6000', 'GRCh38; IDF does not state build')
]:
    sdrf = raw / dataset / f'{dataset}.sdrf.txt'
    matrix, coords, image, gene = zip_status(dataset)
    for r in read_tsv(sdrf):
        sample = textval(r, 'Source Name')
        donor = textval(r, 'Characteristics[isolate]')
        stage = textval(r, 'Characteristics[disease_stage]')
        tissue = textval(r, 'Characteristics[tissue]')
        desc = textval(r, 'Comment[description]')
        title = textval(r, 'Comment[sample_title]')
        region = 'not_recorded_in_SDRF'
        notes = []
        if dataset == 'E-GEAD-622':
            notes.append('pTis; early lesion; FFPE; all four donors are distinct by isolate ID')
            if sample.lower().startswith('rectum') and tissue.lower() == 'colon':
                notes.append('anatomy conflict: source name implies rectum but SDRF tissue field says colon')
        elif dataset == 'E-GEAD-619':
            notes.append(desc)
            if stage.lower() == 'early':
                notes.append('primary early case')
            elif stage.lower() == 'advanced':
                notes.append('secondary advanced case')
        else:
            region = desc
            notes.append('all four capture areas share donor ACC1; not independent patients')
        add(dataset_id=dataset, priority=priority, sample_or_section_id=sample, donor_id=donor,
            section_count_role='one_section_per_SDRF_row', pathology_or_stage=stage,
            pathology_region=region, anatomy=tissue, material=material, platform=platform,
            reference_build=reference, gene_id_status=gene, processed_matrix_status=matrix,
            spot_coordinates_status=coords, he_image_status=image,
            patient_overlap_status='no_explicit_cross_dataset_identifier; see SAMPLE_OVERLAP_AUDIT.md',
            inclusion_status='included_for_audit', evidence_source=f'{dataset}.sdrf.txt; {dataset}.idf.txt',
            notes='; '.join(x for x in notes if x and x != 'NA'))

# Reuse existing GSE226997 metadata; Stage 4B archives remain immutable.
manifest = root / 'metadata' / 'dataset_manifest.tsv'
if manifest.is_file():
    for r in read_tsv(manifest):
        if r.get('accession') == 'GSE226997':
            add(dataset_id='GSE226997', priority='secondary_existing', sample_or_section_id=r.get('sample_id','NA'),
                donor_id=r.get('donor_id','NA'), section_count_role='one_Visium_section_per_public_identifier',
                pathology_or_stage=r.get('condition','NA'), pathology_region='not_recorded_in_GEO_metadata',
                anatomy=r.get('colon_or_rectum','NA'), material='not_recorded_in_GEO_metadata',
                platform=r.get('platform','NA'), reference_build='not_recorded_in_GEO_metadata',
                gene_id_status='requires archive-open smoke test', processed_matrix_status='existing_archive_size_matches_Stage4B_manifest',
                spot_coordinates_status='requires archive-open smoke test', he_image_status='requires archive-open smoke test',
                patient_overlap_status='no_cross-study identifier; not provable from de-identified public metadata',
                inclusion_status='secondary_existing_audit_only', evidence_source='project metadata/dataset_manifest.tsv; Stage4B file manifest',
                notes='cancer-only; cannot test normal-to-adenoma progression')

# GSE200997: parse public GEO sample titles, without loading the expression matrix.
gse_soft = raw / 'GSE200997' / 'GSE200997_family.soft.gz'
sample_id = title = None
if gse_soft.is_file():
    with gzip.open(gse_soft, 'rt', encoding='utf-8', errors='replace') as h:
        for line in h:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE = '):
                if sample_id and title:
                    m = re.search(r'Patient\s+(\d+),\s*(Tumor|Normal)', title, flags=re.I)
                    donor = f'Patient_{m.group(1)}' if m else 'NA'
                    condition = m.group(2).lower() if m else 'NA'
                    anatomy = 'colon' if 'colon' in title.lower() else 'NA'
                    add(dataset_id='GSE200997', priority='reference_candidate', sample_or_section_id=sample_id,
                        donor_id=donor, section_count_role='one_scRNA_library_per_GEO_sample', pathology_or_stage=condition,
                        pathology_region='tumor_or_adjacent_normal; no spatial region labels', anatomy=anatomy,
                        material='fresh_tissue_implied_by_GEO protocol; verify in GEO record',
                        platform='10x Genomics Single Cell 5\' ; Illumina NextSeq 550', reference_build='GRCh38/hg38; Cell Ranger 3.1.0',
                        gene_id_status='gene symbols/IDs require count-file header audit; not inspected in 10B',
                        processed_matrix_status='downloaded; gzip readability checked', spot_coordinates_status='not_applicable_scRNA',
                        he_image_status='not_applicable_scRNA',
                        patient_overlap_status='different public study; no cross-study linkage available',
                        inclusion_status='reference_candidate_only', evidence_source='GSE200997 GEO SOFT metadata',
                        notes='treatment-naive resectable colon cancer; FAP status not reported; unsuitable for adenoma or spatial inference')
                sample_id, title = line.split(' = ', 1)[1], None
            elif line.startswith('!Sample_title = '):
                title = line.split(' = ', 1)[1]
        if sample_id and title:
            m = re.search(r'Patient\s+(\d+),\s*(Tumor|Normal)', title, flags=re.I)
            donor = f'Patient_{m.group(1)}' if m else 'NA'
            condition = m.group(2).lower() if m else 'NA'
            anatomy = 'colon' if 'colon' in title.lower() else 'NA'
            add(dataset_id='GSE200997', priority='reference_candidate', sample_or_section_id=sample_id,
                donor_id=donor, section_count_role='one_scRNA_library_per_GEO_sample', pathology_or_stage=condition,
                pathology_region='tumor_or_adjacent_normal; no spatial region labels', anatomy=anatomy,
                material='fresh_tissue_implied_by_GEO protocol; verify in GEO record', platform='10x Genomics Single Cell 5\' ; Illumina NextSeq 550',
                reference_build='GRCh38/hg38; Cell Ranger 3.1.0', gene_id_status='gene symbols/IDs require count-file header audit; not inspected in 10B',
                processed_matrix_status='downloaded; gzip readability checked', spot_coordinates_status='not_applicable_scRNA', he_image_status='not_applicable_scRNA',
                patient_overlap_status='different public study; no cross-study linkage available', inclusion_status='reference_candidate_only',
                evidence_source='GSE200997 GEO SOFT metadata', notes='treatment-naive resectable colon cancer; FAP status not reported; unsuitable for adenoma or spatial inference')

with open(audit_path, 'w', encoding='utf-8', newline='') as h:
    w = csv.DictWriter(h, fieldnames=columns, delimiter='\t', lineterminator='\n')
    w.writeheader(); w.writerows(rows)

counts = {}
for r in rows:
    counts.setdefault(r['dataset_id'], {'sections': 0, 'donors': set()})
    counts[r['dataset_id']]['sections'] += 1
    counts[r['dataset_id']]['donors'].add(r['donor_id'])

overlap = meta / 'SAMPLE_OVERLAP_AUDIT.md'
with open(overlap, 'w', encoding='utf-8') as h:
    h.write('# Stage 10B Sample-overlap audit\n\n')
    h.write('This is an identifier-based audit of de-identified public metadata; it cannot prove that two separately deposited studies never shared a participant.\n\n')
    h.write('## Within-dataset structure\n\n')
    h.write('- E-GEAD-622: 4 pTis sections from `case1`–`case4`; treat as 4 distinct patients by explicit isolate IDs.\n')
    h.write('- E-GEAD-619: 2 sections from `case5` (early carcinoma in tubulo-villous adenoma) and `case6` (advanced CRC); treat as 2 distinct patients.\n')
    h.write('- E-GEAD-579: 4 capture areas (`A1`–`D1`) all have isolate `ACC1`; treat as 1 patient with 4 nested sections, never as 4 patients.\n')
    h.write('- GSE226997: 4 cancer sections labelled `patient_1`–`patient_4` in the existing project manifest; patient identity is study-specific.\n')
    h.write('- GSE200997: 23 scRNA libraries from 16 CRC patients: 16 tumor libraries and 7 adjacent-normal libraries; seven donors have matched tumor/normal libraries according to GEO titles.\n\n')
    h.write('## Between-dataset assessment\n\n')
    h.write('- E-GEAD-622 and E-GEAD-619 use non-overlapping explicit isolate labels (`case1`–`case4` versus `case5`–`case6`) and are therefore treated as distinct patients.\n')
    h.write('- E-GEAD-579, E-GEAD-619 and E-GEAD-622 share an institutional research program, but the public records provide no cross-accession subject linkage. Cross-dataset overlap is therefore `not demonstrated; not formally excludable`.\n')
    h.write('- GSE226997 and GSE200997 have different public study identifiers and no participant linkage in the available metadata. Overlap is `not demonstrated; not formally excludable`.\n\n')
    h.write('## GSE200997 suitability decision\n\n')
    h.write('`PASS_WITH_LIMITATIONS` as a non-FAP *candidate* external scRNA cancer-versus-adjacent-normal reference: it is an untreated, resectable colon-cancer cohort with processed counts and annotation. However, FAP status is not explicitly reported, seven normal libraries are paired to tumor donors, it contains no adenoma samples and no spatial coordinates/images. It must not be used to infer adenoma progression or treated as proof of non-FAP status.\n')

print(f'rows={len(rows)} audit={audit_path}')
PY

  # Preserve an archive listing as an audit artifact; do not extract or alter
  # any processed archive. The listing verifies that matrix/coordinate/image
  # assets are packaged, while their biological content stays untouched.
  for dataset in E-GEAD-622 E-GEAD-619 E-GEAD-579; do
    zipinfo -1 "$RAW_DIR/$dataset/$dataset.processed.zip" | sed "s#^#$dataset\t#" >> "$META_DIR/PROCESSED_ARCHIVE_CONTENTS.tsv"
  done
  if [[ -s "$META_DIR/articles/PMC11121171.tar.gz" ]]; then
    tar -tzf "$META_DIR/articles/PMC11121171.tar.gz" > "$META_DIR/PMC11121171_PACKAGE_CONTENTS.txt"
  elif [[ -s "$META_DIR/articles/PMC11121171.fulltext.xml" ]]; then
    grep -q '<article' "$META_DIR/articles/PMC11121171.fulltext.xml"
    printf 'Europe PMC full-text XML retrieved; NCBI OA-package endpoint was unavailable (HTTP 404). Package-only supplementary assets remain unavailable.\n' > "$META_DIR/PMC11121171_PACKAGE_CONTENTS.txt"
  else
    printf 'No PMC article support file was retrieved.\n' > "$META_DIR/PMC11121171_PACKAGE_CONTENTS.txt"
  fi
  gzip -t "$RAW_DIR/GSE200997/GSE200997_family.soft.gz"
  gzip -t "$RAW_DIR/GSE200997/GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz"
  gzip -t "$RAW_DIR/GSE200997/GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz"
}

write_gate_and_summary() {
  local now; now="$(timestamp)"
  cat > "$META_DIR/STAGE10B_GATING.md" <<EOF
# Stage 10B Gating: Data Audit and Download

Run ID: $RUN_ID
Completed: $now
Decision: PASS_WITH_LIMITATIONS

## Technical checks

- Official E-GEAD filelists, IDF and SDRF files downloaded.
- Listed processed archives downloaded with source MD5/size validation where supplied.
- Archive listings were inspected without extraction; the audit records matrix, coordinate and image asset presence.
- GSE200997 public processed count matrix, annotation and GEO metadata downloaded; raw reads were not sought.
- Existing GSE226997 archives were reused without copying or re-downloading; Stage 4B manifest size/hash evidence is recorded.
- All new Stage 10B raw files are read-only after completion.
- No FASTQ, raw sequencing .tar file, Space Ranger run or biological analysis was performed.

## Scientific/data limitations

- E-GEAD-622 is FFPE Visium and E-GEAD-619/E-GEAD-579 are fresh-frozen Visium; platform/material effects require cohort-specific handling.
- E-GEAD-579 has four sections from one donor and cannot provide four independent patients.
- The public metadata cannot formally exclude participant overlap across separately de-identified studies.
- GSE200997 is suitable only as a limited cancer-versus-adjacent-normal external scRNA reference; FAP status is unreported and it has no adenoma or spatial data.
- Stage 10 spatial analysis is not authorized by this gate.

## Required outputs

- $META_DIR/DATASET_AUDIT.tsv
- $META_DIR/SAMPLE_OVERLAP_AUDIT.md
- $META_DIR/DOWNLOAD_MANIFEST.tsv
- $META_DIR/DOWNLOAD_MANIFEST.sha256
- $META_DIR/PROCESSED_ARCHIVE_CONTENTS.tsv
EOF
  cp "$META_DIR/STAGE10B_GATING.md" "$REPORTS_DIR/${STAGE}_GATE_DECISION.md"
  {
    printf '# Stage 10B Summary\n\n'
    printf 'Run ID: %s\n' "$RUN_ID"
    printf 'Decision: PASS_WITH_LIMITATIONS\n\n'
    printf 'Processed archives/matrices and metadata were audited and downloaded without raw-read retrieval, archive extraction or biological analysis.\n\n'
    printf 'Output directory: `%s`\n' "$META_DIR"
    printf 'Raw directory (read-only): `%s`\n' "$RAW_DIR"
    printf 'Log: `%s/full.log`\n' "$LOG_DIR"
    printf 'Manifest SHA256: `%s`\n' "$(awk '{print $1}' "$META_DIR/DOWNLOAD_MANIFEST.sha256")"
    printf '\nSee `STAGE10B_GATING.md` for limitations and `SAMPLE_OVERLAP_AUDIT.md` for patient/section nesting.\n'
  } > "$REPORTS_DIR/${STAGE}_SUMMARY.md"
}

write_submission_record() {
  local record="$STATUS_DIR/${STAGE}_SUBMISSION.tsv"
  [[ -n "$RECOVERY_ID" ]] && record="$STATUS_DIR/${STAGE}_RECOVERY_${RECOVERY_ID}_SUBMISSION.tsv"
  cat > "$record" <<EOF
stage\trun_id\tscheduler\tjob_id_or_pid\tlaunch_time\tlaunch_command\twrapper_script\tgit_commit\texpected_outputs\tfull_log\tstatus_check_command\tfailure_recovery_command\tresource_request\tsmoke_test
$STAGE\t$RUN_ID\ttmux\t$JOB_ID_OR_PID\t$(timestamp)\t$LAUNCH_COMMAND\tscripts/10B_spatial_data_audit_download.sh\t$GIT_COMMIT\t$META_DIR/{DATASET_AUDIT.tsv,SAMPLE_OVERLAP_AUDIT.md,DOWNLOAD_MANIFEST.tsv,DOWNLOAD_MANIFEST.sha256,STAGE10B_GATING.md}\t$LOG_DIR/full.log\ttest -f '$STATUS_DIR/$STAGE.SUCCESS' && echo SUCCESS || test -f '$STATUS_DIR/$STAGE.FAILED' && echo FAILED || test -f '$STATUS_DIR/$STAGE.RUNNING' && echo RUNNING || echo NOT_FOUND\tJOB_ID_OR_PID=RETRY GIT_COMMIT=$GIT_COMMIT bash scripts/10B_spatial_data_audit_download.sh --full --resume-failed --root '$PROJECT_ROOT' --run-id '$RUN_ID' --sources '$SOURCES_FILE'\t1 CPU; sequential I/O; no GPU\t$SMOKE_DIR/SMOKE_TEST.tsv
EOF
}

if [[ "$MODE" == "smoke" ]]; then
  mkdir -p "$SMOKE_DIR"
  fetch_metadata_only
  range_file="$SMOKE_DIR/E-GEAD-622.processed.first_1MiB.bin"
  curl --fail --location --retry 3 --retry-all-errors --connect-timeout 30 --max-time 60 --range 0-1048575 \
    --output "$range_file" "https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-622/E-GEAD-622.processed.zip"
  [[ "$(stat -c '%s' "$range_file")" -gt 0 ]]
  [[ -s "$PROJECT_ROOT/data_raw/GSE226997/GSM7089855_Ajou_Visium_P1.tar.gz" ]]
  cat > "$SMOKE_DIR/SMOKE_TEST.tsv" <<EOF
check\tresult\tdetail
metadata_download\tPASS\tOfficial E-GEAD-622 filelist, IDF and SDRF retrieved; remaining metadata is deferred to the full task.
article_package\tDEFERRED\tThe open-access article package is intentionally deferred to the background full download.
range_download\tPASS\tE-GEAD-622 processed archive returned nonempty first 1 MiB.
existing_gse226997\tPASS\tExisting P1 archive is present and nonempty; it was not copied or downloaded.
tools\tPASS\tcurl, checksums, archive listers and Python standard library available.
disk\tPASS\t$(df -h "$PROJECT_ROOT" | tail -1)
scope\tPASS\tNo raw reads, archive extraction, Space Ranger, scoring or biological analysis.
EOF
  printf '[%s] smoke PASS\n' "$(timestamp)" | tee "$LOG_DIR/smoke.log"
  exit 0
fi

if [[ "$RESUME_FAILED" == 1 ]]; then
  [[ ! -e "$STATUS_DIR/$STAGE.RUNNING" && ! -e "$STATUS_DIR/$STAGE.SUCCESS" && -f "$STATUS_DIR/$STAGE.FAILED" ]] || {
    printf 'A failed marker, with no active/success marker, is required for recovery.\n' >&2; exit 3;
  }
  RECOVERY_ID="$(date +%Y%m%dT%H%M%S%z)"
  RECOVERY_EVENT="$STATUS_DIR/${STAGE}_RECOVERY_${RECOVERY_ID}.tsv"
  mv "$STATUS_DIR/$STAGE.FAILED" "$STATUS_DIR/$STAGE.FAILED.$RECOVERY_ID"
  atomic_text "$RECOVERY_EVENT" \
    "stage=$STAGE" "run_id=$RUN_ID" "recovery_id=$RECOVERY_ID" "time=$(timestamp)" \
    "archived_failed_marker=$STATUS_DIR/$STAGE.FAILED.$RECOVERY_ID" \
    "scope=resume_only; preserve_completed_targets_and_part_files" \
    "command=$LAUNCH_COMMAND"
else
  [[ ! -e "$STATUS_DIR/$STAGE.RUNNING" && ! -e "$STATUS_DIR/$STAGE.SUCCESS" && ! -e "$STATUS_DIR/$STAGE.FAILED" ]] || {
    printf 'Existing Stage 10B marker prevents unsafe overwrite.\n' >&2; exit 3;
  }
fi
unlock_stage_raw_dir
atomic_text "$STATUS_DIR/$STAGE.RUNNING" \
  "stage=$STAGE" "run_id=$RUN_ID" "time=$(timestamp)" "job_id_or_pid=$JOB_ID_OR_PID" \
  "command=$LAUNCH_COMMAND" "expected_outputs=$META_DIR" "log=$LOG_DIR/full.log"
write_submission_record
exec > >(tee -a "$LOG_DIR/full.log") 2>&1
printf '[%s] Stage 10B full download started\n' "$(timestamp)"
download_all_planned_assets
discover_cellreports_supplements
write_manifest
inspect_archives_and_write_audit
lock_raw_dir
write_gate_and_summary
atomic_text "$STATUS_DIR/$STAGE.SUCCESS" \
  "stage=$STAGE" "run_id=$RUN_ID" "time=$(timestamp)" \
  "validation=manifest_sha256=$(awk '{print $1}' "$META_DIR/DOWNLOAD_MANIFEST.sha256")" \
  "summary=$REPORTS_DIR/${STAGE}_SUMMARY.md" "gate=$REPORTS_DIR/${STAGE}_GATE_DECISION.md"
rm -f "$STATUS_DIR/$STAGE.RUNNING"
printf '[%s] Stage 10B full download PASS\n' "$(timestamp)"
