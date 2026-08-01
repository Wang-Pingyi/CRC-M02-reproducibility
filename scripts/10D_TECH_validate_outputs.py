#!/usr/bin/env python3
"""Validate Stage 10D-TECH artifacts without reading real spatial outcomes."""

from __future__ import annotations

import csv
import hashlib
import math
import sys
from pathlib import Path


ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage10d_tech"
FIG = ROOT / "figures" / "stage10d_tech"
PRIMARY_BUNDLE = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
SENSITIVITY_BUNDLE = "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8"
COMMON_HASH = "4cd62d74b83673a4d2adf6077bedbdfe73d1cbd369a6a77418d124e0b506d482"
ALLOWED_DECISIONS = {"TECHNICALLY_VALID", "VALID_WITH_LIMITATIONS", "NOT_TECHNICALLY_VALID"}


def sha_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def value_map(path: Path) -> dict[str, str]:
    return {row["parameter"]: row["value"] for row in read_tsv(path)}


def finite(text: str) -> float:
    value = float(text)
    assert math.isfinite(value)
    return value


def main() -> int:
    required = [
        OUT / "STAGE10D_TECH_DESIGN.md",
        OUT / "STAGE10D_TECH_MAPPED_GENESET.tsv",
        OUT / "STAGE10D_TECH_PRIMARY_SCORE_RESULTS.tsv",
        OUT / "STAGE10D_TECH_COVERAGE_RESULTS.tsv",
        OUT / "STAGE10D_TECH_SENSITIVITY_RESULTS.tsv",
        OUT / "STAGE10D_TECH_NULL_MODULE_RESULTS.tsv",
        OUT / "STAGE10D_TECH_DECISION.md",
        OUT / "STAGE10D_TECH_SESSIONINFO.txt",
        OUT / "STAGE10D_TECH_SIMULATION_DESIGN.tsv",
        OUT / "STAGE10D_TECH_RUN_MANIFEST.tsv",
        ROOT / "reports" / "STAGE10D_TECH_SUMMARY.md",
        ROOT / "reports" / "STAGE10D_TECH_GATE_DECISION.md",
    ]
    for path in required:
        assert path.is_file() and path.stat().st_size > 0, path
        assert b"\r\n" not in path.read_bytes(), f"CRLF not permitted: {path}"

    # Re-run the upstream lock validator, then independently verify frozen identities.
    lock36 = read_tsv(ROOT / "results" / "stage10c" / "STAGE10C_LOCK_MANIFEST.tsv")
    assert len(lock36) == 36
    assert {row["bundle_sha256"] for row in lock36} == {PRIMARY_BUNDLE}
    lock35 = read_tsv(ROOT / "results" / "stage10c" / "M02_MINUS_INPP5D_SENS_V1.tsv")
    assert len(lock35) == 35 and all(row["gene"] != "INPP5D" for row in lock35)
    if "bundle_sha256" in lock35[0]:
        assert {row["bundle_sha256"] for row in lock35} == {SENSITIVITY_BUNDLE}

    mapped = read_tsv(OUT / "STAGE10D_TECH_MAPPED_GENESET.tsv")
    assert len(mapped) == 30 and len({row["canonical_gene"] for row in mapped}) == 30
    assert all(row["cohort_id"] == "E-GEAD-622" for row in mapped)
    assert all(row["geneset_sha256"] == COMMON_HASH for row in mapped)
    gene_payload = "".join(f"{gene}\n" for gene in sorted(row["canonical_gene"] for row in mapped))
    assert sha_bytes(gene_payload.encode()) == COMMON_HASH

    params = value_map(ROOT / "config" / "stage10d_tech_parameters.tsv")
    assert params["primary_36_bundle"] == PRIMARY_BUNDLE
    assert params["minus_inpp5d_35_bundle"] == SENSITIVITY_BUNDLE
    assert params["primary_cohort_common_geneset"] == COMMON_HASH
    assert int(params["primary_gene_count"]) == 30
    assert int(params["null_modules"]) == 1000
    assert int(params["seed"]) == 20260731

    primary = read_tsv(OUT / "STAGE10D_TECH_PRIMARY_SCORE_RESULTS.tsv")
    primary_by_metric = {row["metric"]: row for row in primary}
    expected_metrics = {
        "pearson", "spearman", "standardized_bias", "standardized_rmse",
        "calibration_intercept", "calibration_slope", "dynamic_range_ratio",
        "epithelial_residual_slope", "estimated_epithelial_residual_slope",
        "failure_rate", "worst_stress_pearson",
    }
    assert expected_metrics <= set(primary_by_metric)
    for metric in expected_metrics:
        finite(primary_by_metric[metric]["estimate"])
        assert primary_by_metric[metric]["inference_unit"] == "simulation_replicate"
        assert int(primary_by_metric[metric]["n_reference_donors"]) == 4

    coverage = read_tsv(OUT / "STAGE10D_TECH_COVERAGE_RESULTS.tsv")
    tiers = {row["coverage_tier"] for row in coverage}
    assert {"ACTUAL_FROZEN", "ORACLE_TECHNICAL_ONLY", "HIGH_COVERAGE_29", "MINIMUM_COVERAGE_22"} <= tiers
    pearson_coverage = [row for row in coverage if row["metric"] == "pearson"]
    assert len([row for row in pearson_coverage if row["coverage_tier"] == "ACTUAL_FROZEN"]) == 1
    assert len([row for row in pearson_coverage if row["coverage_tier"] == "ORACLE_TECHNICAL_ONLY"]) == 1
    assert len([row for row in pearson_coverage if row["coverage_tier"] == "HIGH_COVERAGE_29"]) == 30
    assert len([row for row in pearson_coverage if row["coverage_tier"] == "MINIMUM_COVERAGE_22"]) == 100
    assert all(int(row["gene_count"]) == 29 for row in pearson_coverage if row["coverage_tier"] == "HIGH_COVERAGE_29")
    assert all(int(row["gene_count"]) == 22 for row in pearson_coverage if row["coverage_tier"] == "MINIMUM_COVERAGE_22")

    sensitivity = read_tsv(OUT / "STAGE10D_TECH_SENSITIVITY_RESULTS.tsv")
    score_ids = {row["score_id"] for row in sensitivity}
    assert {"primary_30", "M02_MINUS_INPP5D_29", "gene_z_30", "UCell_30"} <= score_ids
    assert any(row["sensitivity_type"] == "estimated_epithelial_fraction_threshold" for row in sensitivity)

    nulls = read_tsv(OUT / "STAGE10D_TECH_NULL_MODULE_RESULTS.tsv")
    random_nulls = [row for row in nulls if row["module_id"].startswith("NULL_MATCHED_M30_")]
    assert len(random_nulls) == 1000
    assert all(int(row["gene_count"]) == 30 for row in random_nulls)
    assert len([row for row in nulls if row["module_id"] == "PRIMARY_VS_MATCHED_NULL"]) == 1

    decision_text = (OUT / "STAGE10D_TECH_DECISION.md").read_text(encoding="utf-8")
    decisions = [value for value in ALLOWED_DECISIONS if f"Decision: **{value}**" in decision_text]
    assert len(decisions) == 1
    assert "No real spatial lesion-normal M02 result was read or computed" in decision_text
    assert "Sensitivity scores and the 36-gene oracle were not permitted to rescue" in decision_text

    # Recalculate the decision solely from the frozen primary thresholds.
    metric = lambda name: finite(primary_by_metric[name]["estimate"])
    strict = [
        metric("pearson") >= float(params["strict_pearson_min"]),
        metric("spearman") >= float(params["strict_spearman_min"]),
        abs(metric("standardized_bias")) <= float(params["strict_abs_standardized_bias_max"]),
        metric("standardized_rmse") <= float(params["strict_standardized_rmse_max"]),
        float(params["strict_calibration_slope_min"]) <= metric("calibration_slope") <= float(params["strict_calibration_slope_max"]),
        float(params["strict_dynamic_range_ratio_min"]) <= metric("dynamic_range_ratio") <= float(params["strict_dynamic_range_ratio_max"]),
        abs(metric("epithelial_residual_slope")) <= float(params["strict_abs_epithelial_residual_slope_max"]),
        metric("failure_rate") <= float(params["strict_failure_rate_max"]),
        metric("worst_stress_pearson") >= float(params["strict_worst_stress_pearson_min"]),
    ]
    limited = [
        metric("pearson") >= float(params["limited_pearson_min"]),
        metric("spearman") >= float(params["limited_spearman_min"]),
        metric("standardized_rmse") <= float(params["limited_rmse_max"]),
        metric("failure_rate") <= float(params["limited_failure_rate_max"]),
        metric("worst_stress_pearson") > 0,
    ]
    expected_decision = "TECHNICALLY_VALID" if all(strict) else (
        "VALID_WITH_LIMITATIONS" if all(limited) else "NOT_TECHNICALLY_VALID"
    )
    assert decisions[0] == expected_decision

    for stem in (
        "Fig10D_1_primary_recovery", "Fig10D_2_stress_performance",
        "Fig10D_3_coverage_sensitivity", "Fig10D_4_matched_null",
    ):
        png = FIG / f"{stem}.png"
        pdf = FIG / f"{stem}.pdf"
        assert png.is_file() and png.stat().st_size > 1000 and png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert pdf.is_file() and pdf.stat().st_size > 1000 and pdf.read_bytes()[:4] == b"%PDF"
    for name in (
        "fig1_primary_recovery.tsv", "fig2_stress_performance.tsv",
        "fig3_coverage_performance.tsv", "fig4_null_module.tsv",
    ):
        path = FIG / "source_data" / name
        assert path.is_file() and len(read_tsv(path)) > 0

    # The implementation may read only frozen governance files and GSE161277; reject spatial matrices/ROI outcomes.
    script_text = (ROOT / "scripts" / "10D_TECH_pseudospot_validation.R").read_text(encoding="utf-8")
    forbidden = ["data_raw/stage10b_spatial", "data_processed/stage10b_spatial", "lesion_normal_m02", "adenoma_minus_normal_score"]
    assert all(token not in script_text for token in forbidden)

    manifest = read_tsv(OUT / "STAGE10D_TECH_RUN_MANIFEST.tsv")
    assert manifest and len({row["path"] for row in manifest}) == len(manifest)
    for row in manifest:
        artifact = ROOT / row["path"]
        assert artifact.is_file(), artifact
        assert sha_file(artifact) == row["sha256"], artifact
        assert artifact.stat().st_size == int(float(row["bytes"])), artifact

    print("validation=PASS")
    print(f"decision={expected_decision}")
    print("primary_common_coverage=30/36")
    print(f"common_geneset_sha256={COMMON_HASH}")
    print("spatial_outcome_accessed=NO")
    print("next_stage_started=NO")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

