#!/usr/bin/env python3
"""Audit a 10x feature/MTX pair globally, without ROI labels or contrasts."""

from __future__ import annotations

import argparse
import csv
import gzip


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with open(args.mapping, encoding="utf-8", newline="") as handle:
        locked = list(csv.DictReader(handle, delimiter="\t"))

    features: list[tuple[str, str]] = []
    with gzip.open(args.features, "rt", encoding="utf-8", newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if len(fields) >= 2:
                features.append((fields[0].split(".", 1)[0], fields[1]))

    selected: dict[str, tuple[int, str, str, str, int]] = {}
    for row in locked:
        gene = row["canonical_gene"]
        symbol_hits = [i for i, (_, symbol) in enumerate(features) if symbol == gene]
        if len(symbol_hits) == 1:
            i = symbol_hits[0]
            selected[gene] = (i + 1, features[i][0], features[i][1], "exact_symbol", 1)
            continue
        expected_ids = {
            item for item in row["ensembl_ids"].split(";") if item and item != "NA"
        }
        id_hits = [i for i, (stable_id, _) in enumerate(features) if stable_id in expected_ids]
        if len(symbol_hits) > 1 or len(id_hits) > 1:
            selected[gene] = (-1, "NA", "NA", "unresolved", max(len(symbol_hits), len(id_hits)))
        elif len(id_hits) == 1:
            i = id_hits[0]
            selected[gene] = (
                i + 1,
                features[i][0],
                features[i][1],
                row["mapping_class"],
                1,
            )
        else:
            selected[gene] = (-1, "NA", "NA", "absent_from_feature_space", 0)

    row_to_gene = {
        values[0]: gene for gene, values in selected.items() if values[0] > 0
    }
    totals = {gene: 0 for gene in row_to_gene.values()}
    spot_columns = "NA"
    with gzip.open(args.matrix, "rt", encoding="ascii") as handle:
        dimensions_seen = False
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if not dimensions_seen:
                if len(fields) != 3:
                    raise RuntimeError("invalid MatrixMarket dimensions")
                spot_columns = fields[1]
                dimensions_seen = True
                continue
            row_number = int(fields[0])
            if row_number in row_to_gene:
                totals[row_to_gene[row_number]] += int(float(fields[2]))

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "cohort_id",
                "sample_id",
                "canonical_gene",
                "gene_order",
                "mapping_status",
                "feature_id",
                "feature_symbol",
                "feature_duplicate_count",
                "present_but_zero",
                "mapping_source",
                "spot_columns",
                "global_count_sum",
            ]
        )
        for row in sorted(locked, key=lambda item: int(item["gene_order"])):
            gene = row["canonical_gene"]
            _, feature_id, feature_symbol, status, duplicate_count = selected[gene]
            total = totals.get(gene)
            writer.writerow(
                [
                    args.cohort,
                    args.sample,
                    gene,
                    row["gene_order"],
                    status,
                    feature_id,
                    feature_symbol,
                    duplicate_count,
                    "NA" if total is None else ("TRUE" if total == 0 else "FALSE"),
                    "deposited_10x_feature_annotation",
                    spot_columns,
                    "NA" if total is None else total,
                ]
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
