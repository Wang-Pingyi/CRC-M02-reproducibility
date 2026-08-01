#!/usr/bin/env python3
"""Independent, expression-blind validator for the one formal Stage 10E-R run."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

EXPECTED_HASH = "4cd62d74b83673a4d2adf6077bedbdfe73d1cbd369a6a77418d124e0b506d482"


def rows(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def png_ok(path: Path):
    return path.exists() and path.stat().st_size > 100 and path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--run-id", required=True)
    a = ap.parse_args()
    root = Path(a.root).resolve()
    out = root / "results/stage10e_roi_remediation"
    required = [
        "STAGE10E_R_INPUT_WHITELIST.tsv", "STAGE10E_R_DENYLIST_AUDIT.tsv",
        "STAGE10E_R_FAILURE_DIAGNOSIS.tsv", "STAGE10E_R_REMEDIATION_PLAN_LOCKED.md",
        "STAGE10E_R_COORDINATE_AUDIT.tsv", "STAGE10E_R_TRANSFORM_PARAMETERS.tsv",
        "STAGE10E_R_GEOMETRIC_QC.tsv", "STAGE10E_R_ROI_AGREEMENT.tsv",
        "STAGE10E_R_PATIENT_ELIGIBILITY.tsv", "STAGE10E_R_GENESET_INTEGRITY_CHECK.tsv",
        "STAGE10E_R_ACCEPTANCE.md", "STAGE10E_R_DECISION.md", "STAGE10E_R_RUNTIME.tsv",
    ]
    errors = [f"missing:{name}" for name in required if not (out / name).exists() or (out / name).stat().st_size == 0]
    if errors:
        raise SystemExit(";".join(errors))
    integrity = rows(out / "STAGE10E_R_GENESET_INTEGRITY_CHECK.tsv")
    if len(integrity) != 1 or integrity[0]["recomputed_sha256"] != EXPECTED_HASH or integrity[0]["status"] != "PASS":
        errors.append("frozen_30_of_36_hash_failed")
    deny = rows(out / "STAGE10E_R_DENYLIST_AUDIT.tsv")
    bad_allow = [r for r in deny if r["decision"] == "ALLOW" and (
        any(t in r["path"].lower() for t in ("stage10f", "filtered_feature_bc_matrix", "m02_score", "expression_matrix"))
        or Path(r["path"]).suffix.lower() in (".h5", ".h5ad", ".rds", ".mtx"))]
    if bad_allow:
        errors.append("expression_or_outcome_input_allowed")
    hist = rows(out / "STAGE10E_R_HISTORY_BASELINE_SHA256.tsv")
    for r in hist:
        p = root / r["relative_path"]
        if not p.exists() or sha256(p) != r["sha256_before"]:
            errors.append(f"historical_file_changed:{r['relative_path']}")
    agreements = rows(out / "STAGE10E_R_ROI_AGREEMENT.tsv")
    elig = rows(out / "STAGE10E_R_PATIENT_ELIGIBILITY.tsv")
    geom = rows(out / "STAGE10E_R_GEOMETRIC_QC.tsv")
    transforms = rows(out / "STAGE10E_R_TRANSFORM_PARAMETERS.tsv")
    if not (len(agreements) == len(elig) == len(geom) == 4 and len(transforms) == 8):
        errors.append("unexpected_patient_or_transform_count")
    if any(r["nonlinear_used"] != "FALSE" or r["selection_information"] != "geometry_only" for r in transforms):
        errors.append("transform_rule_violation")
    for r in agreements:
        expected = float(r["agreement"]) >= 0.90 and float(r["cohen_kappa"]) >= 0.85
        if (r["registration_pass"] == "TRUE") != expected:
            errors.append(f"threshold_mismatch:{r['patient_id']}")
    pass_count = sum(r["patient_decision"] == "PASS" for r in elig)
    decision_text = (out / "STAGE10E_R_DECISION.md").read_text(encoding="utf-8")
    if pass_count >= 3 and "PASS_REMEDIATED" not in decision_text:
        errors.append("decision_understates_passing_gate")
    if pass_count < 3 and "PASS_REMEDIATED" in decision_text:
        errors.append("decision_overstates_passing_gate")
    for patient in ("case1", "case2", "case3", "case4"):
        png = root / f"figures/stage10e_roi_remediation/Fig10E_R_{patient}_registration_qc.png"
        src = root / f"figures/stage10e_roi_remediation/source_data/Fig10E_R_{patient}_registration_qc_source_data.tsv"
        if not png_ok(png) or not src.exists() or src.stat().st_size == 0:
            errors.append(f"figure_or_source_data_invalid:{patient}")
    decision = "FAIL" if errors else ("PASS_REMEDIATED" if pass_count >= 3 else (
        "PASS_DESCRIPTIVE_ONLY" if pass_count >= 1 else "FAIL"))
    lines = [
        "# Stage 10E-R independent acceptance", "", f"Decision: **{decision}**", "",
        f"- Formal run ID: `{a.run_id}`", f"- Passing patients: {pass_count}/4.",
        f"- Frozen common-set hash: `{EXPECTED_HASH}`.",
        "- Expression/M02 outcome leakage: NO." if not bad_allow else "- Expression/M02 outcome leakage: DETECTED.",
        "- Historical Stage 6-10E files unchanged: YES." if not any(e.startswith("historical_file_changed") for e in errors) else "- Historical files unchanged: NO.",
        f"- Stage 10F future authorization: {'YES, subject to a separate explicit user instruction; not started here' if decision == 'PASS_REMEDIATED' else 'NO'}.",
        "", "## Validation errors", "",
    ]
    lines.extend([f"- {e}" for e in errors] or ["- None."])
    (out / "STAGE10E_R_ACCEPTANCE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if errors:
        raise SystemExit(";".join(errors))
    print(f"decision={decision}")
    print(f"passing_patients={pass_count}/4")


if __name__ == "__main__":
    main()
