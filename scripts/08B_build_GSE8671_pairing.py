#!/usr/bin/env python3
"""Recover the 32 explicit GSE8671 patient-number pairs from GEO titles."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


PATIENT_RE = re.compile(r"patient\s*#(\d+)", re.IGNORECASE)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    payload = json.loads(args.geo_json.read_text(encoding="utf-8"))
    samples = payload["samples"]
    grouped: dict[int, list[dict]] = {}
    for sample in samples:
        match = PATIENT_RE.search(sample["title"])
        if not match:
            raise RuntimeError(f"No explicit patient number in {sample['sample_id']}")
        number = int(match.group(1))
        tissue = sample["characteristics"]["Tissue"].lower()
        if tissue not in {"normal", "adenoma"}:
            raise RuntimeError(f"Unexpected tissue for {sample['sample_id']}: {tissue}")
        grouped.setdefault(number, []).append(
            {
                "sample_id": sample["sample_id"],
                "tissue": tissue,
                "initials": sample["characteristics"]["Patient ID"],
                "location": sample["characteristics"]["Location"],
                "source_url": sample["source_url"],
                "title": sample["title"],
            }
        )

    if set(grouped) != set(range(1, 33)):
        raise RuntimeError("Expected explicit patient numbers 1 through 32")

    rows = []
    for number in range(1, 33):
        pair = grouped[number]
        if len(pair) != 2 or {x["tissue"] for x in pair} != {"normal", "adenoma"}:
            raise RuntimeError(f"Patient #{number} is not one normal-adenoma pair")
        normal = next(x for x in pair if x["tissue"] == "normal")
        adenoma = next(x for x in pair if x["tissue"] == "adenoma")
        if normal["initials"] != adenoma["initials"]:
            raise RuntimeError(f"Initials disagree within patient #{number}")
        if normal["location"] != adenoma["location"]:
            raise RuntimeError(f"Location disagrees within patient #{number}")
        pair_id = f"GSE8671_P{number:02d}"
        for sample in (normal, adenoma):
            rows.append(
                {
                    "accession": "GSE8671",
                    "sample_id": sample["sample_id"],
                    "verified_donor_id": pair_id,
                    "verified_pair_id": pair_id,
                    "patient_number": number,
                    "original_patient_initials": sample["initials"],
                    "condition": sample["tissue"],
                    "location": sample["location"],
                    "title": sample["title"],
                    "source_url": sample["source_url"],
                    "verification": "explicit_GEO_title_patient_number_and_matching_pair_metadata",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
