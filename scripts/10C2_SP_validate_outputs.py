#!/usr/bin/env python3
"""Validate the frozen Stage 10C2-SP gate without reading spatial outcomes."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage10c2_sp"
PRIMARY_EXPECTED = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
SENSITIVITY_EXPECTED = "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8"


def sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def read_tsv(name: str) -> list[dict[str, str]]:
    with (OUT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    required = [
        "STAGE10C2_SP_SOURCE_AUDIT.tsv",
        "STAGE10C2_SP_PATIENT_SLIDE_ROI_MANIFEST.tsv",
        "STAGE10C2_SP_GENE_ID_MAPPING.tsv",
        "STAGE10C2_SP_GENE_COVERAGE.tsv",
        "STAGE10C2_SP_COHORT_COMMON_GENESET.tsv",
        "STAGE10C2_SP_ANALYSIS_PLAN.md",
        "STAGE10C2_SP_LOCK_MANIFEST.tsv",
        "STAGE10C2_SP_DECISION.md",
    ]
    for name in required:
        path = OUT / name
        assert path.is_file() and path.stat().st_size > 0, name
        assert b"\r\n" not in path.read_bytes(), f"CRLF not allowed: {name}"

    source = read_tsv(required[0])
    assert len(source) == 4
    assert len({row["cohort_id"] for row in source}) == 4
    e622_source = next(row for row in source if row["cohort_id"] == "E-GEAD-622")
    assert e622_source["primary_eligible_paired_patients"] == "4"
    assert e622_source["roi_independent_of_m02"] == "TRUE"

    roi = read_tsv(required[1])
    for patient in ("case1", "case2", "case3", "case4"):
        patient_regions = {
            row["roi_category"]
            for row in roi
            if row["cohort_id"] == "E-GEAD-622" and row["patient_id"] == patient
        }
        assert {"Normal", "Adenoma"}.issubset(patient_regions)
    acc1 = [row for row in roi if row["cohort_id"] == "E-GEAD-579"]
    assert len(acc1) == 4 and {row["patient_id"] for row in acc1} == {"ACC1"}

    mapping = read_tsv(required[2])
    assert len(mapping) == 4 * 36
    assert len({(row["cohort_id"], row["canonical_gene"]) for row in mapping}) == len(mapping)
    assert all(row["mapping_frozen_before_region_contrast"] == "TRUE" for row in mapping)

    coverage = {row["cohort_id"]: row for row in read_tsv(required[3])}
    expected_counts = {"E-GEAD-622": 30, "E-GEAD-619": 36, "E-GEAD-579": 36, "GSE226997": 30}
    for cohort, count in expected_counts.items():
        assert int(coverage[cohort]["common_mapped_genes"]) == count
        assert coverage[cohort]["minimum_22of36_pass"] == "TRUE"
        assert coverage[cohort]["inpp5d_common_mapped"] == "TRUE"
    assert coverage["E-GEAD-622"]["coverage_tier"] == "HIGH_COVERAGE"
    assert coverage["E-GEAD-622"]["m02_minus_inpp5d_common_genes"] == "29"

    common = read_tsv(required[4])
    by_cohort: dict[str, list[dict[str, str]]] = {}
    for row in common:
        by_cohort.setdefault(row["cohort_id"], []).append(row)
    assert {cohort: len(rows) for cohort, rows in by_cohort.items()} == expected_counts
    for cohort, rows in by_cohort.items():
        assert all(row["equal_weight_fraction"] == f"1/{len(rows)}" for row in rows)

    primary_lock_path = ROOT / "results" / "stage10c" / "STAGE10C_LOCK_MANIFEST.tsv"
    with primary_lock_path.open(encoding="utf-8", newline="") as handle:
        primary_lock = list(csv.DictReader(handle, delimiter="\t"))
    primary_genes = sorted(row["gene"] for row in primary_lock)
    assert len(primary_genes) == 36 and len(set(primary_genes)) == 36
    primary_gene_hash = sha(("\n".join(primary_genes) + "\n").encode("utf-8"))
    primary_weighted = "".join(
        f"M02_SCORE_V1\tStem_progenitor_SB_M02\t{gene}\t0.0277777777777778\tlesion_minus_normal_positive\n"
        for gene in primary_genes
    )
    primary_weighted_hash = sha(primary_weighted.encode("utf-8"))
    primary_spec = (
        "M02_SCORE_V1\tprimary=equal_weight_mean_TMM_log2CPM_prior2"
        "\tfuture_primary=patient_paired_region_pseudobulk_equal_weight_mean_TMM_log2CPM_prior2"
        "\thistorical_stage7=equal_weight_mean_gene_z_TMM_log2CPM_prior2"
        "\tsensitivity=equal_weight_gene_z_or_UCell_only\tcoverage=60pct_and_min8"
        "\tcoverage_sensitivity=80pct\tmissing=canonical_only_equal_renormalization"
        "\tdirection=lesion_minus_normal_positive\n"
    )
    primary_spec_hash = sha(primary_spec.encode("utf-8"))
    primary_source_hash = primary_lock[0]["source_modules_locked_sha256"]
    primary_bundle = (
        f"gene_set_sha256\t{primary_gene_hash}\n"
        f"weighted_signature_sha256\t{primary_weighted_hash}\n"
        f"scoring_spec_sha256\t{primary_spec_hash}\n"
        f"source_modules_locked_sha256\t{primary_source_hash}\n"
    )
    assert sha(primary_bundle.encode("utf-8")) == PRIMARY_EXPECTED

    sensitivity_path = ROOT / "results" / "stage10c" / "M02_MINUS_INPP5D_SENS_V1.tsv"
    with sensitivity_path.open(encoding="utf-8", newline="") as handle:
        sensitivity_lock = list(csv.DictReader(handle, delimiter="\t"))
    sensitivity_genes = sorted(row["gene"] for row in sensitivity_lock)
    assert len(sensitivity_genes) == 35 and "INPP5D" not in sensitivity_genes
    sensitivity_gene_hash = sha(("\n".join(sensitivity_genes) + "\n").encode("utf-8"))
    sensitivity_weighted = "".join(
        f"M02_MINUS_INPP5D_SENS_V1\tStem_progenitor_SB_M02\t{gene}\t0.0285714285714286\tlesion_minus_normal_positive\n"
        for gene in sensitivity_genes
    )
    sensitivity_weighted_hash = sha(sensitivity_weighted.encode("utf-8"))
    sensitivity_spec = (
        "M02_MINUS_INPP5D_SENS_V1\trole=contamination_sensitivity_only"
        "\tparent_score=M02_SCORE_V1\tparent_module=Stem_progenitor_SB_M02"
        "\tremoved_gene=INPP5D\tgene_count=35\tweight=0.0285714285714286"
        "\tnormalization=TMM_log2CPM_prior2\taggregation=patient_by_region_pseudobulk"
        "\teffect=lesion_minus_matched_normal\tdirection=lesion_minus_normal_positive"
        "\tcoverage=60pct_and_min8_21of35\tcoverage_sensitivity=80pct_28of35"
        "\tmissing=canonical_only_equal_renormalization"
        "\tselection=outcome_independent_no_further_gene_deletion\n"
    )
    sensitivity_spec_hash = sha(sensitivity_spec.encode("utf-8"))
    sensitivity_bundle = (
        f"gene_set_sha256\t{sensitivity_gene_hash}\n"
        f"weighted_signature_sha256\t{sensitivity_weighted_hash}\n"
        f"scoring_spec_sha256\t{sensitivity_spec_hash}\n"
        f"parent_primary_bundle_sha256\t{PRIMARY_EXPECTED}\n"
    )
    assert sha(sensitivity_bundle.encode("utf-8")) == SENSITIVITY_EXPECTED

    lock = read_tsv(required[6])
    lock_by_id = {row["lock_id"]: row for row in lock}
    assert lock_by_id["M02_SCORE_V1_36_gene_primary"]["sha256"] == PRIMARY_EXPECTED
    assert lock_by_id["M02_MINUS_INPP5D_SENS_V1"]["sha256"] == SENSITIVITY_EXPECTED
    for name in required[:6]:
        expected = sha((OUT / name).read_bytes())
        assert lock_by_id[name]["sha256"] == expected
    genes_payload = "".join(f"{row['canonical_gene']}\n" for row in sorted(by_cohort["E-GEAD-622"], key=lambda row: row["canonical_gene"]))
    assert lock_by_id["E-GEAD-622"]["sha256"] == sha(genes_payload.encode("utf-8"))

    decision = (OUT / required[7]).read_text(encoding="utf-8")
    assert "Decision: **PASS**" in decision
    assert "Stage 10D-TECH is allowed" in decision
    assert "Stage 10E and Stage 10F are not authorized" in decision
    assert "No lesion-normal M02 score" in decision

    print("validation=PASS")
    print(f"files={len(required)}")
    print("primary_patients=4")
    print("primary_common_coverage=30/36")
    print(f"primary_bundle={PRIMARY_EXPECTED}")
    print(f"sensitivity_bundle={SENSITIVITY_EXPECTED}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
