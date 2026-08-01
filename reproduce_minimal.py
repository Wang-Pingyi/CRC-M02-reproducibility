#!/usr/bin/env python3
import csv, hashlib, pathlib, sys
root = pathlib.Path(__file__).resolve().parent
def rows(name):
    with (root/name).open(encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))
def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda:f.read(1024*1024), b""): h.update(c)
    return h.hexdigest()
errors=[]
for line in (root/"SHA256SUMS").read_text(encoding="utf-8").splitlines():
    expected, rel = line.split("  ",1); p=root/rel
    if not p.is_file() or digest(p)!=expected: errors.append(f"checksum:{rel}")
required=["README.md","CODEBOOK.md","DATA_ACCESSIONS.tsv","PATIENT_SAMPLE_MAP.tsv","DATASET_INDEPENDENCE.tsv","GENE_ID_MAPPING.tsv","COHORT_COMMON_GENESETS.tsv","SCRIPT_OUTPUT_MAP.tsv","source_data/stage10fg/STAGE10F_SKIPPED.md","source_data/stage10fg/STAGE10G_SKIPPED.md","source_data/NEGATIVE_AND_NOT_ESTIMABLE.tsv"]
for rel in required:
    if not (root/rel).is_file(): errors.append(f"missing:{rel}")
genes=[r for r in rows("COHORT_COMMON_GENESETS.tsv") if r["geneset_id"]=="M02_SCORE_V1" and r["cohort_or_context"]=="CANONICAL"]
if len(genes)!=36 or len({r["gene"] for r in genes})!=36: errors.append("canonical_m02_not_36_unique")
if not rows("DATASET_INDEPENDENCE.tsv"): errors.append("independence_empty")
if not rows("SCRIPT_OUTPUT_MAP.tsv"): errors.append("script_map_empty")
if not rows("source_data/NEGATIVE_AND_NOT_ESTIMABLE.tsv"): errors.append("boundary_results_empty")
print("PASS" if not errors else "FAIL")
for e in errors: print(e)
sys.exit(1 if errors else 0)
