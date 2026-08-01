#!/usr/bin/env python3
"""Build a small, auditable GDC file-to-patient/sample metadata table."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


API = "https://api.gdc.cancer.gov/files"


def read_ids(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    out = []
    for row in rows:
        if row.get("accession") != "TCGA-COAD":
            continue
        version = row.get("source_version", "")
        if "GDC_file_uuid:" not in version:
            continue
        out.append((row["file_name"], version.split("GDC_file_uuid:", 1)[1]))
    if not out:
        raise RuntimeError("No TCGA-COAD GDC UUIDs found in file_manifest.tsv")
    return out


def fetch_batch(ids: list[str], retries: int = 5) -> dict:
    fields = ",".join(
        [
            "file_id",
            "file_name",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.sample_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
            "cases.samples.tissue_type",
        ]
    )
    query = urllib.parse.urlencode(
        {
            "filters": json.dumps(
                {"op": "in", "content": {"field": "files.file_id", "value": ids}}
            ),
            "fields": fields,
            "expand": "cases.samples",
            "format": "JSON",
            "size": str(len(ids)),
        }
    )
    request = urllib.request.Request(
        f"{API}?{query}",
        headers={"User-Agent": "CRC-carcinogenesis-stage8B/1.0"},
    )
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except Exception:
            if attempt == retries:
                raise
            time.sleep(min(2**attempt, 20))
    raise AssertionError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--batch-size", type=int, default=100)
    args = parser.parse_args()

    local_uuid = read_ids(args.file_manifest)
    expected = {uuid: name for name, uuid in local_uuid}
    hits: dict[str, dict] = {}
    ids = list(expected)
    for start in range(0, len(ids), args.batch_size):
        payload = fetch_batch(ids[start : start + args.batch_size])
        for hit in payload["data"]["hits"]:
            hits[hit["file_id"]] = hit

    missing = sorted(set(expected) - set(hits))
    if missing:
        raise RuntimeError(f"GDC metadata missing for {len(missing)} UUIDs")

    output_rows = []
    for uuid in ids:
        hit = hits[uuid]
        cases = hit.get("cases") or []
        if len(cases) != 1:
            raise RuntimeError(f"{uuid} has {len(cases)} cases; expected one")
        case = cases[0]
        samples = case.get("samples") or []
        if len(samples) != 1:
            raise RuntimeError(f"{uuid} has {len(samples)} samples; expected one")
        sample = samples[0]
        output_rows.append(
            {
                "local_file_name": expected[uuid],
                "gdc_file_id": uuid,
                "gdc_file_name": hit.get("file_name", ""),
                "case_id": case.get("case_id", ""),
                "patient_id": case.get("submitter_id", ""),
                "sample_id": sample.get("sample_id", ""),
                "sample_barcode": sample.get("submitter_id", ""),
                "sample_type": sample.get("sample_type", ""),
                "tissue_type": sample.get("tissue_type", ""),
            }
        )

    if len(output_rows) != len(local_uuid):
        raise RuntimeError("GDC output row count does not match manifest")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
