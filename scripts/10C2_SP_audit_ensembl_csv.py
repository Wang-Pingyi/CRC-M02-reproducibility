#!/usr/bin/env python3
"""Audit an Ensembl-ID spatial count CSV without using ROI labels or M02 contrasts."""

from __future__ import annotations

import argparse
import csv
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--sample", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with open(args.mapping, encoding="utf-8", newline="") as handle:
        mapping_rows = list(csv.DictReader(handle, delimiter="\t"))

    id_to_rows: dict[str, list[dict[str, str]]] = {}
    for row in mapping_rows:
        for ensembl_id in row["ensembl_ids"].split(";"):
            if ensembl_id and ensembl_id != "NA":
                id_to_rows.setdefault(ensembl_id, []).append(row)

    observed: dict[str, list[tuple[int, int]]] = {}
    reader = csv.reader(sys.stdin)
    header = next(reader, None)
    if not header or len(header) < 2:
        raise RuntimeError("missing or invalid CSV header")
    for fields in reader:
        if not fields:
            continue
        stable_id = fields[0].split(".", 1)[0]
        if stable_id not in id_to_rows:
            continue
        total = 0
        for value in fields[1:]:
            if value:
                total += int(float(value))
        observed.setdefault(stable_id, []).append((len(fields) - 1, total))

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
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
    for row in sorted(mapping_rows, key=lambda item: int(item["gene_order"])):
        candidate_ids = [
            item for item in row["ensembl_ids"].split(";") if item and item != "NA"
        ]
        hits = [
            (stable_id, metrics)
            for stable_id in candidate_ids
            for metrics in observed.get(stable_id, [])
        ]
        if len(candidate_ids) != 1:
            status = "unresolved"
        elif len(hits) == 0:
            status = "absent_from_feature_space"
        elif len(hits) == 1:
            status = row["mapping_class"]
        else:
            status = "unresolved"
        if len(hits) == 1:
            feature_id, (spot_columns, total) = hits[0]
            zero = "TRUE" if total == 0 else "FALSE"
        else:
            feature_id, spot_columns, total, zero = "NA", "NA", "NA", "NA"
        writer.writerow(
            [
                args.cohort,
                args.sample,
                row["canonical_gene"],
                row["gene_order"],
                status,
                feature_id,
                row["canonical_gene"] if len(hits) == 1 else "NA",
                len(hits),
                zero,
                row["mapping_source"],
                spot_columns,
                total,
            ]
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
