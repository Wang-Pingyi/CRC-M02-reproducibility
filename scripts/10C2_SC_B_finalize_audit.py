#!/usr/bin/env python3
"""Create the frozen Stage 10C2-SC-B manifests and gate documents.

Inputs are the already completed, M02-blind structural count audit and the
historical Stage 10C lock. This script does not read expression matrices.
"""

from __future__ import annotations

import csv
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_BUNDLE = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
ARCHIVE_SHA256 = "944909a666daf223160d8a2f2f4480dc2e45b451031dd07a280b61821adb121f"
ARCHIVE_MD5 = "35e39ded57c2ddd6729c77c75784cc5c"
RAW_PATH = "${CRC_PROJECT_ROOT}/data_raw/stage10c2_sc/figshare_29925404/v1/LST_count_mat.tar.gz"
DERIVED_PATH = "${CRC_PROJECT_ROOT}/data_processed/stage10c2_sc/figshare_29925404/v1/derived"


SAMPLE_METADATA = {
    "P1_L": ("P1", "LST-G", "tubular adenoma", "rectum", "TRUE", "P1_NORMAL_LSTG", "TRUE"),
    "P1_N": ("P1", "normal mucosa", "normal", "rectum", "TRUE", "P1_NORMAL_LSTG", "TRUE"),
    "P2_L": ("P2", "LST-G", "tubulovillous adenoma", "transverse colon", "FALSE", "NA", "NA"),
    "P3_L": ("P3", "LST-G", "tubulovillous adenoma", "caecum", "FALSE", "P3_MULTI_LESION", "NA"),
    "P3_P": ("P3", "protruded adenoma", "tubular adenoma", "sigmoid colon", "FALSE", "P3_MULTI_LESION", "NA"),
    "P4_L": ("P4", "LST-G", "tubular adenoma", "rectum", "FALSE", "NA", "NA"),
    "P5_L": ("P5", "LST-G", "tubular adenoma", "caecum", "TRUE", "P5_NORMAL_LSTG", "TRUE"),
    "P5_N": ("P5", "normal mucosa", "normal", "caecum", "TRUE", "P5_NORMAL_LSTG", "TRUE"),
    "P6_L": ("P6", "LST-G", "tubulovillous adenoma", "sigmoid colon", "FALSE", "NA", "NA"),
    "P7_P": ("P7", "protruded adenoma", "tubular adenoma", "sigmoid colon", "FALSE", "NA", "NA"),
}

ARCHIVE_MEMBERS = {
    "P1_L.csv": 227602925,
    "P1_N.csv": 172637478,
    "P2_L.csv": 97979898,
    "P3_L.csv": 164255756,
    "P3_P.csv": 112034239,
    "P4_L.csv": 132976259,
    "P5_L.csv": 182434109,
    "P5_N.csv": 174185107,
    "P6_L.csv": 127098954,
    "P7_P.csv": 98494672,
}

EPITHELIAL_POSITIVE = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
EPITHELIAL_EXCLUSION = ["PTPRC", "LST1", "TYROBP", "COL1A1", "COL1A2", "COL3A1", "DCN", "PECAM1", "VWF", "KDR"]
STEM_POSITIVE = ["LGR5", "OLFM4", "SMOC2", "PROM1", "LRIG1", "SOX9"]
DIFFERENTIATION_EXCLUSION = ["KRT20", "CA1", "GUCA2A", "FABP1", "MUC2", "TFF3", "CHGA", "POU2F3"]
CYCLING = ["MKI67", "TOP2A", "UBE2C"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows for {path}")
    fields = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_lock() -> tuple[list[dict[str, str]], set[str]]:
    rows = read_tsv(ROOT / "results/stage10c/STAGE10C_LOCK_MANIFEST.tsv")
    hashes = {row["bundle_sha256"] for row in rows}
    genes = {row["gene"] for row in rows}
    if len(rows) != 36 or len(genes) != 36 or hashes != {EXPECTED_BUNDLE}:
        raise RuntimeError("Primary 36-gene lock failed verification")
    independent_markers = set(EPITHELIAL_POSITIVE + EPITHELIAL_EXCLUSION + STEM_POSITIVE + DIFFERENTIATION_EXCLUSION + CYCLING)
    overlap = independent_markers & genes
    if overlap:
        raise RuntimeError(f"Independent QC markers overlap locked M02 genes: {sorted(overlap)}")
    return rows, genes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--audit",
        type=Path,
        default=ROOT / "data_processed/stage10c2_sc/figshare_29925404/v1/general_qc_matrix_audit.tsv",
        help="Structural/general-QC audit produced by 10C2_SC_B_audit_counts.py",
    )
    args = parser.parse_args()
    _, _ = read_lock()
    matrix_rows = read_tsv(args.audit)
    by_sample = {row["sample_id"]: row for row in matrix_rows}
    if set(by_sample) != set(SAMPLE_METADATA):
        raise RuntimeError("Recovered sample set is not the frozen ten-sample set")

    total_cells = sum(int(row["cells"]) for row in matrix_rows)
    low_feature_cells = sum(int(row["cells_nFeature_lt1000"]) for row in matrix_rows)
    written_qc_cells = sum(int(row["cells_meeting_published_written_qc"]) for row in matrix_rows)
    if total_cells != 31213:
        raise RuntimeError(f"Unexpected cell total: {total_cells}")
    if any(row["numeric_type"] != "nonnegative_integer_counts" for row in matrix_rows):
        raise RuntimeError("At least one matrix is not nonnegative integer counts")

    metadata_dir = ROOT / "data/metadata"
    result_dir = ROOT / "results/stage10c2_sc"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    file_manifest = [{
        "repository": "Figshare",
        "article_id": "29925404",
        "file_id": "57235721",
        "doi_base": "10.6084/m9.figshare.29925404",
        "doi_versioned": "10.6084/m9.figshare.29925404.v1",
        "version": "1",
        "direct_file_name": "LST_count_mat.tar.gz",
        "official_url": "https://ndownloader.figshare.com/files/57235721",
        "download_date": "2026-07-31",
        "bytes": 74626211,
        "official_md5": ARCHIVE_MD5,
        "verified_sha256": ARCHIVE_SHA256,
        "license": "CC BY 4.0",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "server_raw_path": RAW_PATH,
        "raw_file_mode": "0444",
        "validation": "bytes_match;md5_match;sha256_recorded;gzip_magic_1F8B",
        "transport_note": "Official Figshare URL resolved over local IPv6; bytes relayed to server because server IPv4 received HTTP 403; server checksum matched before local temporary deletion",
    }]
    write_tsv(metadata_dir / "stage10c2_sc_file_manifest.tsv", file_manifest)

    patient_rows = []
    for sample in SAMPLE_METADATA:
        patient, morphology, histology, site, paired, pair_id, same_site = SAMPLE_METADATA[sample]
        patient_rows.append({
            "dataset_id": "FIGSHARE29925404",
            "patient_id": patient,
            "sample_id": sample,
            "library_id": f"FIG29925404_{sample}_LIB1",
            "biological_sample_type": "normal mucosa" if morphology == "normal mucosa" else "adenoma/precancer lesion",
            "lesion_morphology": morphology,
            "histology": histology,
            "anatomic_site": site,
            "paired_normal_available": paired,
            "pair_id": pair_id,
            "same_anatomic_site_as_matched_normal": same_site,
            "FAP_status": "NA",
            "germline_APC_status": "NA",
            "technical_replicate": "FALSE",
            "deposited_cells": int(by_sample[sample]["cells"]),
            "matrix_file": f"{sample}.csv",
            "mapping_evidence": "article Table 1 plus exact archive member name",
        })
    write_tsv(metadata_dir / "stage10c2_sc_patient_sample_manifest.tsv", patient_rows)

    inventory_rows = []
    for member, size in ARCHIVE_MEMBERS.items():
        inventory_rows.append({
            "archive": "LST_count_mat.tar.gz",
            "archive_sha256": ARCHIVE_SHA256,
            "member_name": member,
            "member_type": "regular_file",
            "archive_member_mode": "0777",
            "archive_member_mtime": "2024-11-01 23:56",
            "uncompressed_bytes": size,
            "derived_file_mode": "0640",
            "sample_id": member.removesuffix(".csv"),
            "content_interpretation": "raw nonnegative integer UMI count matrix; genes by cells",
            "safe_relative_path": "TRUE",
        })
    write_tsv(result_dir / "STAGE10C2_SC_B_ARCHIVE_INVENTORY.tsv", inventory_rows)

    final_matrix_rows = []
    for row in matrix_rows:
        final_matrix_rows.append({
            "sample_id": row["sample_id"],
            "file_name": row["file_name"],
            "file_bytes": row["file_bytes"],
            "genes": row["genes"],
            "cells": row["cells"],
            "matrix_class": "raw_UMI_counts_in_author_QC_cell_set",
            "orientation": row["orientation"],
            "gene_identifier": row["gene_identifier"],
            "barcode_format": "10x_sequence-dash-1",
            "sample_prefix_embedded": row["sample_prefix_embedded"],
            "barcode_to_sample_mapping": "archive member/file name; prefix sample_id before combining",
            "numeric_type": row["numeric_type"],
            "duplicate_genes": row["duplicate_genes"],
            "duplicate_barcodes_within_file": row["duplicate_barcodes"],
            "negative_values": row["negative_values"],
            "noninteger_values": row["noninteger_values"],
            "sparsity_fraction": f"{float(row['sparsity_fraction']):.6f}",
            "nCount_min": row["nCount_min"],
            "nCount_median": row["nCount_median"],
            "nCount_max": row["nCount_max"],
            "nFeature_min": row["nFeature_min"],
            "nFeature_median": row["nFeature_median"],
            "nFeature_max": row["nFeature_max"],
            "mitochondrial_pct_median": f"{float(row['mitochondrial_pct_median']):.4f}",
            "mitochondrial_pct_p95": f"{float(row['mitochondrial_pct_p95']):.4f}",
            "cells_nFeature_lt1000": row["cells_nFeature_lt1000"],
            "cells_mitochondrial_pct_gt50": row["cells_mitochondrial_pct_gt50"],
            "reusable_cell_annotation": "absent_from_Figshare_archive",
            "epithelial_identifiability": "supported_by_independent_marker_panel",
            "audit_note": "Cell count is part of the exact 31,213-cell deposit; cells below written 1,000-feature cutoff are retained in the deposit and flagged for frozen QC sensitivity",
        })
    write_tsv(result_dir / "STAGE10C2_SC_B_MATRIX_AUDIT.tsv", final_matrix_rows)

    gates = [
        ("G01", "10C-ALIGN decision PASS/PASS_WITH_LIMITATIONS", "PASS", "STAGE10C_ALIGN_DECISION.md is PASS", "required precondition"),
        ("G02", "Original 36-gene primary bundle hash unchanged", "PASS", EXPECTED_BUNDLE, "no historical lock change"),
        ("G03", "10C2-SC-A one-time search closed", "PASS", "STAGE10C2_SC_DECISION.md search closure", "no new cohort considered"),
        ("G04", "All official Figshare files downloaded", "PASS", "Article v1 declares one file and one file was retrieved", "complete object"),
        ("G05", "Official size and MD5 verified", "PASS", f"74626211 bytes; MD5 {ARCHIVE_MD5}", "file integrity"),
        ("G06", "Raw file and raw hierarchy read-only", "PASS", "file 0444; stage raw directories 0555", "raw protected"),
        ("G07", "Archive paths safe and all members inventoried", "PASS", "10 top-level CSV members; no traversal entries", "safe extraction"),
        ("G08", "Required ten samples reconstructed", "PASS", ";".join(SAMPLE_METADATA), "exact sample set"),
        ("G09", "Seven patients reconstructed", "PASS", "P1-P7; P3_L and P3_P share P3", "donor is inferential unit"),
        ("G10", "Only P1 and P5 are normal-LST-G pairs", "PASS", "P1 rectum/rectum; P5 caecum/caecum", "same-site pairs"),
        ("G11", "Paper, archive, and matrices reconcile on cell count", "PASS", f"7 patients; 10 samples; {total_cells} cells", "exact count agreement"),
        ("G12", "Matrix values are reliable counts", "PASS", "All ten matrices are nonnegative integers; 23,756 gene symbols each", "raw UMI counts supported"),
        ("G13", "Reusable cell annotation deposited", "LIMITATION", "No annotation or README member is present", "independent annotation required"),
        ("G14", "Written paper QC exactly reproduced by deposit", "LIMITATION", f"{low_feature_cells} deposited cells have nFeature<1000 although all have nCount>=1000; {written_qc_cells}/{total_cells} meet written criteria", "freeze primary and deposit-preserving QC sensitivity"),
        ("G15", "Epithelial lineage can be identified independently", "PASS", "All independent positive/exclusion markers are represented; marker sets exclude all 36 M02 genes", "annotation feasible without M02"),
        ("G16", "License permits reanalysis", "PASS", "CC BY 4.0 with attribution", "reuse allowed"),
        ("G17", "No prior-cohort matrices technically embedded", "PASS", "Archive contains only the ten Table-1 sample matrices and exact study cell total", "no technical double counting detected"),
        ("G18", "Absolute patient non-overlap proven", "LIMITATION", "De-identification prevents proof; risk recorded as low but not zero", "claim remains directional only"),
        ("G19", "M02 outcome remained unopened", "PASS", "Only file structure, count properties, general QC, and independent markers were audited", "preanalysis blinding preserved"),
        ("G20", "Stage terminal decision", "PASS_AS_DIRECTIONAL_SENSITIVITY_ONLY", "Counts/mapping/license/epithelial gate pass with declared limitations", "10C3 eligible but not started"),
    ]
    write_tsv(result_dir / "STAGE10C2_SC_B_QC_GATE.tsv", [
        {"gate_id": a, "criterion": b, "result": c, "evidence": d, "impact": e}
        for a, b, c, d, e in gates
    ])

    metadata_audit = f"""# Stage 10C2-SC-B Metadata and File Audit

Date: 2026-07-31

Scope: Figshare article 29925404 / DOI `10.6084/m9.figshare.29925404` only.

M02 expression inspected or calculated: **No**.

## Decision

`PASS_AS_DIRECTIONAL_SENSITIVITY_ONLY`

The object is acceptable for a locked, directional sensitivity analysis. It is
not confirmatory: only P1 and P5 supply matched normal-LST-G pairs, hereditary
status is unreported, and reusable cell-level annotations are absent.

## Official file and integrity

Figshare version 1 declares one file, `LST_count_mat.tar.gz` (file ID 57235721),
with 74,626,211 bytes and official MD5 `{ARCHIVE_MD5}` under CC BY 4.0. The
retrieved file matches both the official byte count and MD5; its recorded SHA256
is `{ARCHIVE_SHA256}`.

Figshare returned HTTP 403 over IPv4 from both the server and local client. The
official downloader endpoint was therefore resolved over local IPv6, the
official bytes were relayed to the server, and the server MD5/SHA256 were checked
before the local temporary copy was deleted. The server raw file is mode `0444`;
its raw parent hierarchy is mode `0555`.

## Archive inventory and matrix type

The archive contains exactly ten safe top-level CSV members and no README or
cell-annotation file. Uncompressed members total {sum(ARCHIVE_MEMBERS.values()):,}
bytes. Each CSV has 23,756 gene-symbol rows and cells in columns. Values are
nonnegative integers with no nonfinite values, no duplicate genes, and no
duplicate barcodes within a file. Sparsity ranges from
{min(float(r['sparsity_fraction']) for r in matrix_rows):.3%} to
{max(float(r['sparsity_fraction']) for r in matrix_rows):.3%}. These properties
support interpretation as raw UMI counts for the authors' retained cell set,
not a normalized expression matrix.

Cell barcodes use standard 10x `sequence-1` form but do not contain sample
prefixes. The archive member name is the barcode-to-sample mapping; any combined
object must prepend `sample_id` before concatenation.

## Patient and sample reconstruction

Seven patients and ten samples were recovered exactly:

- P1: `P1_N` normal rectum and `P1_L` rectal LST-G tubular adenoma;
- P2: `P2_L` transverse-colon LST-G tubulovillous adenoma;
- P3: `P3_L` caecal LST-G tubulovillous adenoma and `P3_P` sigmoid protruded
  tubular adenoma;
- P4: `P4_L` rectal LST-G tubular adenoma;
- P5: `P5_N` normal caecum and `P5_L` caecal LST-G tubular adenoma;
- P6: `P6_L` sigmoid LST-G tubulovillous adenoma;
- P7: `P7_P` sigmoid protruded tubular adenoma.

P1 and P5 are the only normal-LST-G pairs; each pair is from the same anatomic
site. P3_L and P3_P are two lesions from one patient and can never be counted as
two donors. There is one matrix/library per sample and no evidence of technical
replicate libraries. `FAP_status` and `germline_APC_status` remain `NA` for every
patient.

## Cell-count reconciliation and QC discrepancy

The matrices contain exactly {total_cells:,} cells, matching the paper's post-QC
total. All ten sample IDs and the paper's 7-patient/10-sample structure match.

However, {low_feature_cells:,} deposited cells have fewer than 1,000 detected
genes, although the Methods text states that cells below 1,000 genes were
excluded. All deposited cells have at least 1,000 UMIs, none exceeds 8,000
detected genes, and none exceeds 50% mitochondrial reads. Thus
{written_qc_cells:,}/{total_cells:,} deposited cells meet the written
`1,000-8,000 genes and <=50% mitochondrial` rule. The exact reason cannot be
proven from the deposit; a possible nCount/nFeature wording mismatch is an
inference, not a fact. Stage 10C3 must retain both the written-rule primary QC
and a deposit-preserving sensitivity, frozen before M02 is opened.

## Annotation and epithelial identifiability

No reusable cell annotation is present in the Figshare object. Independent
lineage identification is nevertheless feasible because all prespecified
epithelial positive markers (`{', '.join(EPITHELIAL_POSITIVE)}`), exclusion
markers, and independent Stem/progenitor markers are present in every gene
index. These panels were checked to have zero overlap with the frozen 36-gene
M02 membership.

## Overlap audit

The archive includes only the ten matrices named in the paper's sample table;
their cell total equals this study's reported total. No accession, filename, or
matrix component from GSE161277, GSE201348, GSE261388, or another public cohort
is embedded, so technical double counting was not detected. This does not prove
absolute patient independence: de-identified public records cannot exclude an
unreported shared patient. The residual patient-overlap risk is recorded as
`low but not zero`.

## Limitations controlling downstream claims

- Only two matched normal-LST-G patients are available.
- LST-G is not interchangeable with ordinary sporadic adenoma.
- FAP and germline APC status are unknown.
- P3 contributes two lesions but one donor.
- No reusable cell-level annotation was deposited.
- The written QC rule and deposited cell set are not fully concordant.
- This cohort can supply directional sensitivity evidence only.
"""
    (result_dir / "STAGE10C2_SC_B_METADATA_AUDIT.md").write_text(metadata_audit, encoding="utf-8")

    qc_gate_rows = "\n".join(
        f"- `{name}`" for name in EPITHELIAL_POSITIVE + EPITHELIAL_EXCLUSION + STEM_POSITIVE + DIFFERENTIATION_EXCLUSION + CYCLING
    )
    preanalysis_lock = f"""# Stage 10C3 Preanalysis Lock

Lock date: 2026-07-31

Dataset: Figshare 29925404 v1 only

Authorized role: directional sensitivity only

Stage 10C3 status: not started

M02 outcome viewed before this lock: **No**

## Immutable upstream lock

- Primary module: frozen 36-gene `Stem_progenitor_SB_M02`.
- Score version: `M02_SCORE_V1`.
- Direction: lesion minus matched normal, positive.
- Primary bundle hash: `{EXPECTED_BUNDLE}`.
- No genes, directions, weights, samples, histology labels, or pair definitions
  may change after an M02 result is viewed.

## Frozen inferential structure

- Patient/donor is the biological replicate.
- P1 and P5 are the only matched normal-LST-G pairs and are the primary
  directional comparison.
- P3_L and P3_P are nested in patient P3 and never supply two degrees of
  freedom. Their LST-G-versus-protruded contrast is descriptive only.
- Unpaired lesion samples may support annotation and direction summaries but
  cannot manufacture normal-lesion paired replication.
- With only two matched patients, report both patient effects and direction;
  do not present this cohort as confirmatory or use cell counts as degrees of
  freedom.

## Frozen cell and sample QC

1. Read each CSV independently and prepend `sample_id` to every barcode.
2. Recompute `nCount`, `nFeature`, mitochondrial and ribosomal percentages.
3. Primary cell QC follows the published written rule: `1,000 <= nFeature <=
   8,000` and mitochondrial percentage `<=50%`.
4. A prespecified deposit-preserving QC sensitivity retains the authors'
   deposited cells after removing zero-count, negative/noninteger, or
   mitochondrial-`>50%` cells. It cannot replace the primary result.
5. Run `scDblFinder` independently by sample after primary QC with a recorded
   deterministic seed; do not infer doublets across pooled patients.
6. Ambient RNA is assessed diagnostically from lineage-conflicting marker
   patterns. No unfiltered droplet matrix exists, so an empty-droplet model is
   not claimed.
7. Primary epithelial sample inclusion requires at least **50 epithelial
   cells**. The only cell-count sensitivity thresholds are **25** and **100**
   epithelial cells. These thresholds cannot change after M02 is opened.
8. The secondary Stem/progenitor state requires at least **20 cells per
   sample**; otherwise it is reported as not estimable for that sample.

## Frozen epithelial identification

Perform unsupervised clustering within each sample/cohort using QC-passing
cells. A cluster is epithelial only when it shows coherent enrichment of
`EPCAM`, `KRT8`, `KRT18`, `KRT19`, and/or `CDH1` (at least two members at the
cluster level) and lacks a dominant immune, stromal, or endothelial exclusion
signature (`PTPRC`, `LST1`, `TYROBP`, `COL1A1`, `COL1A2`, `COL3A1`, `DCN`,
`PECAM1`, `VWF`, `KDR`). Mixed-lineage clusters are excluded or marked
uncertain. Published labels, if later obtained, are auxiliary evidence only.

## Frozen independent Stem/progenitor annotation

The secondary state requires coherent expression of at least two of `LGR5`,
`OLFM4`, `SMOC2`, `PROM1`, `LRIG1`, and `SOX9`, together with low dominance of
differentiation markers (`KRT20`, `CA1`, `GUCA2A`, `FABP1`, `MUC2`, `TFF3`,
`CHGA`, `POU2F3`). Cycling cells identified by `MKI67`, `TOP2A`, and `UBE2C`
are labeled separately unless they also satisfy the independent progenitor
panel. If no stable state satisfies these rules, the secondary analysis is
not estimable; no cluster may be relabeled to obtain a favorable M02 result.

All markers above are outside the frozen 36-gene M02 set. Marker/M02 overlap
was programmatically verified as zero.

## Frozen 10C3 score and reporting constraints

- Aggregate raw counts within patient x sample/condition x epithelial state.
- Apply the previously frozen equal-weight TMM-log2CPM pseudobulk score and
  fixed coverage rules; UCell or gene-wise-z scores remain sensitivity-only.
- Report every eligible patient effect, mapped-gene coverage, direction, and
  all non-estimable results.
- Do not select cells, clusters, QC rules, score variants, or pair definitions
  using M02 direction, P values, or effect size.
- This cohort can only strengthen or weaken directional sensitivity evidence;
  it cannot raise the existing claim ceiling.

## Start gate

Stage 10C3 is eligible for a separate authorization because 10C2-SC-B is
`PASS_AS_DIRECTIONAL_SENSITIVITY_ONLY`. This file does not start 10C3.
"""
    (result_dir / "STAGE10C3_PREANALYSIS_LOCK.md").write_text(preanalysis_lock, encoding="utf-8")

    decision = f"""# Stage 10C2-SC-B Decision

Date: 2026-07-31

Decision: **PASS_AS_DIRECTIONAL_SENSITIVITY_ONLY**

## Basis

- The only official Figshare v1 file was retrieved and verified by byte count,
  official MD5, and SHA256.
- Ten raw integer UMI count matrices reconstruct the exact required sample set,
  7 patients, and {total_cells:,} deposited cells.
- P1 and P5 are the only same-site normal-LST-G pairs.
- P3_L and P3_P are locked to one donor.
- FAP and germline APC status remain `NA`.
- General marker panels independent of all 36 M02 genes support epithelial and
  Stem/progenitor identification.
- CC BY 4.0 permits reanalysis with attribution.
- No prior-cohort matrix is technically embedded; de-identified patient overlap
  remains low but not provably zero.

## Mandatory limitations

- This dataset is not confirmatory and may provide directional sensitivity
  evidence only.
- No reusable cell annotation is present.
- {low_feature_cells:,} deposited cells conflict with the paper's written
  1,000-feature lower threshold; the frozen 10C3 plan includes a primary
  written-rule QC and a deposit-preserving sensitivity.
- Only two matched patients are available; cells cannot be treated as
  independent replicates.

## Authorization

Stage 10C3 is **eligible for a separate investigator authorization** under the
frozen `STAGE10C3_PREANALYSIS_LOCK.md`. It has not been started. This decision
does not authorize 10C2-SP or any spatial inference.

Primary bundle hash verified: `{EXPECTED_BUNDLE}`.

Archive SHA256: `{ARCHIVE_SHA256}`.
"""
    (result_dir / "STAGE10C2_SC_B_DECISION.md").write_text(decision, encoding="utf-8")

    print(f"Wrote Stage 10C2-SC-B outputs; decision PASS_AS_DIRECTIONAL_SENSITIVITY_ONLY; cells={total_cells}; low_nFeature={low_feature_cells}")


if __name__ == "__main__":
    main()
