#!/usr/bin/env python3
"""Protocol-level acceptance checks for Stage 10C3 small outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PRIMARY_BUNDLE = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
SENS_BUNDLE = "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_na(value: str) -> bool:
    return value in {"", "NA", "NaN", "nan"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    results = root / "results" / "stage10c3"
    figures = root / "figures" / "stage10c3"
    checks: list[tuple[str, bool, str]] = []

    def add(name: str, passed: bool, detail: str) -> None:
        checks.append((name, bool(passed), detail))

    required = [
        "STAGE10C3_ANALYSIS_PLAN_LOCKED.md", "STAGE10C3_PATIENT_SAMPLE_SCORES.tsv",
        "STAGE10C3_MATCHED_PAIR_EFFECTS.tsv", "STAGE10C3_SECONDARY_DESCRIPTIVE.tsv",
        "STAGE10C3_SENSITIVITY.tsv", "STAGE10C3_QC.md", "STAGE10C3_CLAIM_LIMITS.md",
        "STAGE10C3_DECISION.md",
    ]
    add("required_outputs", all((results / name).stat().st_size > 0 for name in required), ";".join(required))

    decision_text = (results / "STAGE10C3_DECISION.md").read_text(encoding="utf-8")
    allowed = ["LST_DIRECTIONAL_CONCORDANCE", "MIXED_OR_METHOD_DEPENDENT", "NULL_OR_OPPOSITE", "NOT_ESTIMABLE"]
    found = [item for item in allowed if f"Decision: **{item}**" in decision_text]
    add("allowed_decision", len(found) == 1, found[0] if found else "none_or_multiple")
    add("hashes_recorded", PRIMARY_BUNDLE in decision_text and SENS_BUNDLE in decision_text, "primary_and_35gene_bundle")

    lock = root / "results" / "stage10c2_sc" / "STAGE10C3_PREANALYSIS_LOCK.md"
    plan = results / "STAGE10C3_ANALYSIS_PLAN_LOCKED.md"
    score_file = results / "STAGE10C3_PATIENT_SAMPLE_SCORES.tsv"
    add("lock_precedes_score", lock.stat().st_mtime <= plan.stat().st_mtime < score_file.stat().st_mtime,
        f"lock={lock.stat().st_mtime};plan={plan.stat().st_mtime};score={score_file.stat().st_mtime}")

    coverage = read_tsv(results / "STAGE10C3_GENE_COVERAGE.tsv")
    add("primary_coverage_not_estimated", all(row["mapped36"] == "33" and row["primary_36of36"] == "FALSE" for row in coverage), "33/36 in every scenario")
    add("coverage_status_unambiguous", all(row["status"] == "SENSITIVITY_ONLY_PRIMARY_NOT_ESTIMABLE" for row in coverage), "no sensitivity labeled primary")

    scores = read_tsv(score_file)
    primary_rows = [row for row in scores if row["scenario_id"] == "PRIMARY_ALL_EPI_MIN50" and row["score_method"] == "M02_SCORE_V1"]
    add("primary_scores_are_na", len(primary_rows) == 10 and all(is_na(row["score"]) and row["status"] == "NOT_ESTIMABLE" for row in primary_rows), f"rows={len(primary_rows)}")
    add("patient_is_inference_unit", all(row["inferential_unit"] == "patient" for row in scores), "all score rows")

    sensitivity = read_tsv(results / "STAGE10C3_SENSITIVITY.tsv")
    def sens(method: str, scenario: str = "PRIMARY_ALL_EPI_MIN50") -> dict[str, str]:
        rows = [row for row in sensitivity if row["scenario_id"] == scenario and row["score_method"] == method]
        if len(rows) != 1:
            raise AssertionError(f"Expected one row for {scenario}/{method}, got {len(rows)}")
        return rows[0]
    cov29 = sens("M02_29OF36_COVERAGE_SENS")
    sens35 = sens("M02_MINUS_INPP5D_SENS_V1")
    add("coverage_sensitivity_positive_2of2", cov29["k_positive"] == "2" and float(cov29["P1_effect"]) > 0 and float(cov29["P5_effect"]) > 0, f"P1={cov29['P1_effect']};P5={cov29['P5_effect']}")
    add("35gene_sensitivity_no_reverse", sens35["k_positive"] == "2" and float(sens35["P1_effect"]) >= 0 and float(sens35["P5_effect"]) >= 0, f"P1={sens35['P1_effect']};P5={sens35['P5_effect']}")
    complete = [row for row in sensitivity if row["status"] == "ESTIMABLE" and row["n_pairs"] == "2"]
    systemic_opposite = any(float(row["P1_effect"]) < 0 and float(row["P5_effect"]) < 0 for row in complete)
    add("no_systematic_opposite_sensitivity", not systemic_opposite, f"complete_sensitivities={len(complete)}")
    published = [row for row in sensitivity if row["scenario_id"] == "PUBLISHED_ANNOTATION_SENS"]
    add("published_annotation_visible_not_estimable", len(published) == 1 and published[0]["status"] == "NOT_ESTIMABLE", "archive lacked cell-level labels")

    effects = read_tsv(results / "STAGE10C3_MATCHED_PAIR_EFFECTS.tsv")
    individual = [row for row in effects if row["record_type"] == "individual_patient_effect"]
    add("only_p1_p5_matched_patients", {row["patient_id"] for row in individual} == {"P1", "P5"}, "P3 never counted as independent matched donor")
    inferential_fields = ("ci_lower", "ci_upper", "p_value", "fdr")
    add("no_n2_population_inference", all(is_na(row[field]) for row in effects for field in inferential_fields), "CI/P/FDR all NA")

    secondary = read_tsv(results / "STAGE10C3_SECONDARY_DESCRIPTIVE.tsv")
    lesion_patients = {row["patient_id"] for row in secondary if row["description_type"] == "lesion_sample"}
    add("all_seven_lesion_patients_described", lesion_patients == {f"P{i}" for i in range(1, 8)}, ",".join(sorted(lesion_patients)))
    p3 = [row for row in secondary if row["description_type"] == "within_patient_lesion_contrast"]
    add("p3_nested_contrast", len(p3) == 1 and p3[0]["patient_id"] == "P3", "P3_L-minus-P3_P descriptive only")

    qc = read_tsv(results / "STAGE10C3_QC_SOURCE.tsv")
    primary_qc = [row for row in qc if row["qc_mode"] == "primary_qc"]
    add("all_samples_qc_and_epithelial_gate", len(primary_qc) == 10 and all(row["epithelial_min50"] == "TRUE" for row in primary_qc), f"samples={len(primary_qc)}")
    pair_stem = {row["sample_id"]: row["stem_min20"] for row in primary_qc if row["sample_id"] in {"P1_L", "P1_N", "P5_L", "P5_N"}}
    add("stem_pair_not_rescued", set(pair_stem) == {"P1_L", "P1_N", "P5_L", "P5_N"} and not all(value == "TRUE" for value in pair_stem.values()), str(pair_stem))

    pngs = sorted(figures.glob("*.png"))
    source_tables = sorted(figures.glob("*_source_data.tsv"))
    add("figures_and_source_data", len(pngs) == 3 and len(source_tables) == 3 and all(path.stat().st_size > 1000 for path in pngs + source_tables), f"png={len(pngs)};source={len(source_tables)}")

    output = results / "STAGE10C3_VALIDATION.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "status", "detail"])
        for name, passed, detail in checks:
            writer.writerow([name, "PASS" if passed else "FAIL", detail])

    failed = [name for name, passed, _ in checks if not passed]
    if failed:
        raise SystemExit("Validation failed: " + ", ".join(failed))
    print(f"PASS: {len(checks)} Stage 10C3 acceptance checks")


if __name__ == "__main__":
    main()
