#!/usr/bin/env python3
"""Replace non-unique GSE8671 initials with explicit GEO patient-number keys."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    args = parser.parse_args()

    rows = read_tsv(args.manifest)
    pairs = {row["sample_id"]: row for row in read_tsv(args.pairs)}
    target = [row for row in rows if row["accession"] == "GSE8671"]
    if len(target) != 64 or len(pairs) != 64:
        raise RuntimeError("Expected 64 GSE8671 manifest and pairing rows")

    corrected = 0
    for row in rows:
        if row["accession"] != "GSE8671":
            continue
        pair = pairs.get(row["sample_id"])
        if pair is None:
            raise RuntimeError(f"Pairing absent for {row['sample_id']}")
        original = row["donor_id"]
        row["donor_id"] = pair["verified_donor_id"]
        row["paired_group"] = pair["verified_pair_id"]
        row["donor_id_status"] = "verified_explicit_GEO_patient_number"
        audit_note = (
            f" Pairing correction 2026-07-29: GEO title explicitly assigns patient "
            f"#{pair['patient_number']}; original non-unique initials were "
            f"{pair['original_patient_initials']} (previous donor key {original})."
        )
        if "Pairing correction 2026-07-29" not in row["metadata_notes"]:
            row["metadata_notes"] += audit_note
        corrected += 1

    if corrected != 64:
        raise RuntimeError("Not all GSE8671 rows were corrected")
    fields = list(rows[0])
    temporary = args.manifest.with_suffix(".tsv.stage8b_tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(args.manifest)


if __name__ == "__main__":
    main()
