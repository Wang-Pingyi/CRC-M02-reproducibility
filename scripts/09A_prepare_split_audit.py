#!/usr/bin/env python3
"""Prepare the GSE99573 split and test-firewall audit without reading CEL data."""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys
import tarfile
from collections import Counter


EXPECTED = {"training": 265, "testing": 65, "not_used": 8}
EXPECTED_INCLUSION = {
    "training": "include_training_only",
    "testing": "reserve_locked_test_set",
    "not_used": "exclude",
}


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: pathlib.Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_dir", type=pathlib.Path)
    parser.add_argument("run_id")
    args = parser.parse_args()

    project = args.project_dir.resolve()
    manifest_path = project / "metadata" / "dataset_manifest.tsv"
    raw_tar = project / "data_raw" / "GSE99573" / "GSE99573_RAW.tar"
    result_dir = project / "results" / "09A_stool_feasibility" / args.run_id
    result_dir.mkdir(parents=True, exist_ok=True)

    rows = [r for r in read_tsv(manifest_path) if r["accession"] == "GSE99573"]
    if len(rows) != 338:
        raise RuntimeError(f"Expected 338 GSE99573 manifest rows, found {len(rows)}")

    sample_ids = [r["sample_id"] for r in rows]
    donor_ids = [r["donor_id"] for r in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise RuntimeError("GSE99573 sample_id is not unique")
    if not all(donor_ids) or any(x == "NA" for x in donor_ids):
        raise RuntimeError("At least one GSE99573 donor_id is missing")
    if len(donor_ids) != len(set(donor_ids)):
        raise RuntimeError("GSE99573 donor_id is not unique; replicate audit required")

    split_counts = Counter(r["validation_split"] for r in rows)
    if dict(split_counts) != EXPECTED:
        raise RuntimeError(f"Split mismatch: observed {dict(split_counts)}, expected {EXPECTED}")
    for row in rows:
        expected_inclusion = EXPECTED_INCLUSION[row["validation_split"]]
        if row["inclusion"] != expected_inclusion:
            raise RuntimeError(
                f"{row['sample_id']}: inclusion={row['inclusion']} does not match "
                f"split={row['validation_split']}"
            )

    with tarfile.open(raw_tar, mode="r") as archive:
        members = [m.name for m in archive.getmembers() if m.isfile()]

    cel_members = [m for m in members if re.search(r"\.CEL\.gz$", m, re.IGNORECASE)]
    if len(cel_members) != 338:
        raise RuntimeError(f"Expected 338 CEL.gz archive members, found {len(cel_members)}")

    member_by_sample: dict[str, str] = {}
    for sample_id in sample_ids:
        matches = [
            member for member in cel_members
            if pathlib.PurePosixPath(member).name.startswith(sample_id + "_")
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"{sample_id}: expected exactly one CEL.gz archive member, found {len(matches)}"
            )
        member_by_sample[sample_id] = matches[0]

    split_summary: list[dict[str, object]] = []
    for split in ("training", "testing", "not_used"):
        split_rows = [r for r in rows if r["validation_split"] == split]
        for condition in ("normal", "adenoma", "cancer", "benign_not_used"):
            count = sum(r["condition"] == condition for r in split_rows)
            if count:
                split_summary.append(
                    {
                        "validation_split": split,
                        "condition": condition,
                        "n_samples": count,
                        "source": "GEO GSM set and disease-status characteristics",
                        "expression_accessed": "FALSE",
                    }
                )
    write_tsv(
        result_dir / "GSE99573_split_audit.tsv",
        split_summary,
        [
            "validation_split",
            "condition",
            "n_samples",
            "source",
            "expression_accessed",
        ],
    )

    inventory_rows: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda x: x["sample_id"]):
        inventory_rows.append(
            {
                "sample_id": row["sample_id"],
                "donor_id": row["donor_id"],
                "validation_split": row["validation_split"],
                "inclusion": row["inclusion"],
                "raw_archive_member": member_by_sample[row["sample_id"]],
                "raw_member_present": "TRUE",
                "technical_replicate": "FALSE",
                "expression_extracted": "FALSE",
                "expression_accessed": "FALSE",
            }
        )
    write_tsv(
        result_dir / "GSE99573_sample_inventory.tsv",
        inventory_rows,
        [
            "sample_id",
            "donor_id",
            "validation_split",
            "inclusion",
            "raw_archive_member",
            "raw_member_present",
            "technical_replicate",
            "expression_extracted",
            "expression_accessed",
        ],
    )

    test_rows = [
        {
            "sample_id": r["sample_id"],
            "validation_split": "testing",
            "raw_archive_member": member_by_sample[r["sample_id"]],
            "inventory_checked": "TRUE",
            "cel_extracted": "FALSE",
            "expression_accessed": "FALSE",
            "permitted_stage_9A_use": "split_and_archive_inventory_only",
        }
        for r in sorted(rows, key=lambda x: x["sample_id"])
        if r["validation_split"] == "testing"
    ]
    write_tsv(
        result_dir / "locked_test_access_audit.tsv",
        test_rows,
        [
            "sample_id",
            "validation_split",
            "raw_archive_member",
            "inventory_checked",
            "cel_extracted",
            "expression_accessed",
            "permitted_stage_9A_use",
        ],
    )

    training_members = [
        member_by_sample[r["sample_id"]]
        for r in sorted(rows, key=lambda x: x["sample_id"])
        if r["validation_split"] == "training"
    ]
    (result_dir / "training_tar_members.txt").write_text(
        "\n".join(training_members) + "\n", encoding="utf-8"
    )
    (result_dir / "training_sample_ids.txt").write_text(
        "\n".join(
            r["sample_id"]
            for r in sorted(rows, key=lambda x: x["sample_id"])
            if r["validation_split"] == "training"
        )
        + "\n",
        encoding="utf-8",
    )

    validations = [
        ("manifest_rows", len(rows) == 338, str(len(rows))),
        ("sample_id_unique", len(sample_ids) == len(set(sample_ids)), str(len(set(sample_ids)))),
        ("donor_id_complete", all(donor_ids), str(sum(bool(x) for x in donor_ids))),
        ("donor_id_unique", len(donor_ids) == len(set(donor_ids)), str(len(set(donor_ids)))),
        ("training_count", split_counts["training"] == 265, str(split_counts["training"])),
        ("testing_count", split_counts["testing"] == 65, str(split_counts["testing"])),
        ("not_used_count", split_counts["not_used"] == 8, str(split_counts["not_used"])),
        ("archive_cel_count", len(cel_members) == 338, str(len(cel_members))),
        ("one_cel_per_sample", len(member_by_sample) == 338, str(len(member_by_sample))),
        ("test_expression_not_accessed", True, "65 locked inventory rows"),
    ]
    write_tsv(
        result_dir / "split_validation_checks.tsv",
        [
            {
                "check": name,
                "passed": str(passed).upper(),
                "observed": observed,
            }
            for name, passed, observed in validations
        ],
        ["check", "passed", "observed"],
    )
    if not all(passed for _, passed, _ in validations):
        return 2

    print(f"STAGE9A_SPLIT_AUDIT_OK run_id={args.run_id}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
