#!/usr/bin/env python3
"""Audit Figshare 29925404 count CSVs without computing the locked M02 score.

The script only inspects matrix structure, integer-count properties, general
single-cell QC metrics, and marker panels that are independent of the frozen
36-gene module. It never selects or evaluates M02 genes.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


SAMPLES = (
    "P1_L", "P1_N", "P2_L", "P3_L", "P3_P",
    "P4_L", "P5_L", "P5_N", "P6_L", "P7_P",
)

# These marker panels were chosen independently of M02 and intentionally do
# not contain any member of the frozen 36-gene signature.
EPITHELIAL_POSITIVE = {"EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"}
EPITHELIAL_EXCLUSION = {
    "PTPRC", "LST1", "TYROBP", "COL1A1", "COL1A2", "COL3A1", "DCN",
    "PECAM1", "VWF", "KDR",
}
STEM_PROGENITOR_POSITIVE = {"LGR5", "OLFM4", "SMOC2", "PROM1", "LRIG1", "SOX9"}
DIFFERENTIATION_EXCLUSION = {
    "KRT20", "CA1", "GUCA2A", "FABP1", "MUC2", "TFF3", "CHGA", "POU2F3",
}
CYCLING_MARKERS = {"MKI67", "TOP2A", "UBE2C"}

BARCODE_RE = re.compile(r"^[ACGTN]+-[0-9]+$")


def percentile(values: np.ndarray, q: float) -> float:
    return float(np.percentile(values, q)) if values.size else math.nan


def audit_one(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader([handle.readline()]))
        if len(header) < 2:
            raise ValueError(f"No cell columns in {path}")
        barcodes = header[1:]
        n_cells = len(barcodes)
        total_counts = np.zeros(n_cells, dtype=np.int64)
        n_features = np.zeros(n_cells, dtype=np.int32)
        mito_counts = np.zeros(n_cells, dtype=np.int64)
        ribo_counts = np.zeros(n_cells, dtype=np.int64)
        marker_counts: Counter[str] = Counter()
        seen_genes: set[str] = set()
        duplicate_genes = 0
        n_genes = 0
        nnz = 0
        negative_values = 0
        noninteger_values = 0
        nonfinite_values = 0
        max_value = 0
        min_value = 0
        ensembl_like = 0

        marker_union = (
            EPITHELIAL_POSITIVE
            | EPITHELIAL_EXCLUSION
            | STEM_PROGENITOR_POSITIVE
            | DIFFERENTIATION_EXCLUSION
            | CYCLING_MARKERS
        )

        for line_number, line in enumerate(handle, start=2):
            line = line.rstrip("\r\n")
            if not line:
                continue
            try:
                gene_field, numeric = line.split(",", 1)
            except ValueError as exc:
                raise ValueError(f"Malformed row {line_number} in {path}") from exc
            gene = gene_field.strip('"')
            values = np.fromstring(numeric, sep=",", dtype=np.float64)
            if values.size != n_cells:
                raise ValueError(
                    f"Row {line_number} in {path} has {values.size} values; expected {n_cells}"
                )
            nonfinite_values += int(np.count_nonzero(~np.isfinite(values)))
            negative_values += int(np.count_nonzero(values < 0))
            noninteger_values += int(np.count_nonzero(values != np.floor(values)))
            if values.size:
                max_value = max(max_value, int(np.nanmax(values)))
                min_value = min(min_value, int(np.nanmin(values)))
            integer_values = values.astype(np.int64, copy=False)
            nonzero = integer_values != 0
            nnz += int(np.count_nonzero(nonzero))
            total_counts += integer_values
            n_features += nonzero
            if gene.startswith("MT-"):
                mito_counts += integer_values
            if gene.startswith("RPS") or gene.startswith("RPL"):
                ribo_counts += integer_values
            if gene in marker_union:
                marker_counts[gene] = int(np.count_nonzero(nonzero))
            if gene in seen_genes:
                duplicate_genes += 1
            else:
                seen_genes.add(gene)
            if gene.startswith("ENSG"):
                ensembl_like += 1
            n_genes += 1

    safe_total = np.maximum(total_counts, 1)
    mito_pct = mito_counts / safe_total * 100.0
    ribo_pct = ribo_counts / safe_total * 100.0
    paper_qc_mask = (n_features >= 1000) & (n_features <= 8000) & (mito_pct <= 50.0)
    total_values = n_genes * n_cells
    barcode_counts = Counter(barcodes)
    duplicate_barcodes = sum(count - 1 for count in barcode_counts.values() if count > 1)
    pattern_matches = sum(bool(BARCODE_RE.fullmatch(item)) for item in barcodes)

    return {
        "sample_id": path.stem,
        "file_name": path.name,
        "file_bytes": path.stat().st_size,
        "genes": n_genes,
        "cells": n_cells,
        "orientation": "genes_by_cells",
        "gene_identifier": "Ensembl" if ensembl_like > n_genes / 2 else "gene_symbol",
        "duplicate_genes": duplicate_genes,
        "duplicate_barcodes": duplicate_barcodes,
        "barcode_pattern_fraction": pattern_matches / n_cells,
        "sample_prefix_embedded": False,
        "numeric_type": "nonnegative_integer_counts"
        if not (negative_values or noninteger_values or nonfinite_values)
        else "other_numeric",
        "negative_values": negative_values,
        "noninteger_values": noninteger_values,
        "nonfinite_values": nonfinite_values,
        "min_value": min_value,
        "max_value": max_value,
        "nonzero_values": nnz,
        "total_values": total_values,
        "sparsity_fraction": 1.0 - nnz / total_values,
        "zero_count_cells": int(np.count_nonzero(total_counts == 0)),
        "zero_feature_cells": int(np.count_nonzero(n_features == 0)),
        "cells_nCount_lt1000": int(np.count_nonzero(total_counts < 1000)),
        "cells_nFeature_lt1000": int(np.count_nonzero(n_features < 1000)),
        "cells_nFeature_gt8000": int(np.count_nonzero(n_features > 8000)),
        "cells_mitochondrial_pct_gt50": int(np.count_nonzero(mito_pct > 50.0)),
        "cells_meeting_published_written_qc": int(np.count_nonzero(paper_qc_mask)),
        "nCount_min": int(total_counts.min()),
        "nCount_median": float(np.median(total_counts)),
        "nCount_p95": percentile(total_counts, 95),
        "nCount_max": int(total_counts.max()),
        "nFeature_min": int(n_features.min()),
        "nFeature_median": float(np.median(n_features)),
        "nFeature_p95": percentile(n_features, 95),
        "nFeature_max": int(n_features.max()),
        "mitochondrial_pct_median": float(np.median(mito_pct)),
        "mitochondrial_pct_p95": percentile(mito_pct, 95),
        "ribosomal_pct_median": float(np.median(ribo_pct)),
        "ribosomal_pct_p95": percentile(ribo_pct, 95),
        "epithelial_positive_markers_present": ";".join(sorted(EPITHELIAL_POSITIVE & seen_genes)),
        "epithelial_exclusion_markers_present": ";".join(sorted(EPITHELIAL_EXCLUSION & seen_genes)),
        "stem_progenitor_markers_present": ";".join(sorted(STEM_PROGENITOR_POSITIVE & seen_genes)),
        "differentiation_exclusion_markers_present": ";".join(sorted(DIFFERENTIATION_EXCLUSION & seen_genes)),
        "cycling_markers_present": ";".join(sorted(CYCLING_MARKERS & seen_genes)),
        "marker_detected_cell_counts_json": json.dumps(dict(sorted(marker_counts.items())), separators=(",", ":")),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    paths = [args.matrix_dir / f"{sample}.csv" for sample in SAMPLES]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing expected matrices: " + ", ".join(missing))

    rows = [audit_one(path) for path in paths]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps({
        "samples": len(rows),
        "patients_represented": 7,
        "cells_total": sum(int(row["cells"]) for row in rows),
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
