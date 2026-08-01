#!/usr/bin/env python3
"""Independent validation of Stage 10E-DESC small outputs; no expression is read."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

EXPECTED_HASH = "4cd62d74b83673a4d2adf6077bedbdfe73d1cbd369a6a77418d124e0b506d482"
ALLOWED_DECISIONS = {"DESCRIPTIVE_POSITIVE", "DESCRIPTIVE_NEGATIVE", "DESCRIPTIVE_METHOD_DEPENDENT", "NOT_ESTIMABLE"}


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def hash_genes(genes):
    return hashlib.sha256(("\n".join(sorted(genes)) + "\n").encode()).hexdigest()


def direction(x, tol=1e-12):
    if x > tol:
        return "positive"
    if x < -tol:
        return "negative"
    return "zero"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    root = Path(args.root).resolve()
    out = root / "results/stage10e_desc"
    fig = root / "figures/stage10e_desc"
    required = [
        "STAGE10E_DESC_ANALYSIS_PLAN_LOCKED.md", "STAGE10E_DESC_MAPPED_GENESET.tsv",
        "STAGE10E_DESC_CASE4_SCORES.tsv", "STAGE10E_DESC_SENSITIVITY.tsv",
        "STAGE10E_DESC_CLAIM_LIMITS.md", "STAGE10E_DESC_DECISION.md",
        "STAGE10E_DESC_RUN_MANIFEST.tsv", "STAGE10E_DESC_SESSIONINFO.txt",
    ]
    errors = [f"missing:{name}" for name in required if not (out / name).exists() or (out / name).stat().st_size == 0]
    for name in ("STAGE10E_DESC_SUMMARY.md", "STAGE10E_DESC_GATE_DECISION.md"):
        if not (root / "reports" / name).exists():
            errors.append(f"missing_report:{name}")
    if errors:
        raise SystemExit(";".join(errors))

    mapped = read_tsv(out / "STAGE10E_DESC_MAPPED_GENESET.tsv")
    genes = [r["canonical_gene"] for r in mapped]
    if len(mapped) != 30 or len(set(genes)) != 30 or hash_genes(genes) != EXPECTED_HASH:
        errors.append("mapped_geneset_or_hash_invalid")
    if any(r.get("patient_id") != "case4" or r.get("expression_feature_present") != "TRUE" for r in mapped):
        errors.append("mapped_table_patient_or_presence_invalid")

    scores = read_tsv(out / "STAGE10E_DESC_CASE4_SCORES.tsv")
    if len(scores) != 2 or {r["region"] for r in scores} != {"Normal", "Adenoma"}:
        errors.append("primary_score_regions_invalid")
    if any(r["patient_id"] != "case4" or r["patient_n"] != "1" for r in scores):
        errors.append("primary_score_patient_boundary_invalid")
    if any(int(r["spot_count"]) < 30 for r in scores):
        errors.append("primary_region_coverage_below_30")
    for r in scores:
        if r["p_value"] != "NOT_COMPUTED_BY_DESIGN" or r["confidence_interval"] != "NOT_COMPUTED_BY_DESIGN" or r["fdr"] != "NOT_COMPUTED_BY_DESIGN":
            errors.append("inferential_quantity_present")

    sensitivity = read_tsv(out / "STAGE10E_DESC_SENSITIVITY.tsv")
    expected_sens = {"primary_reference", "M02_MINUS_INPP5D", "gene_z", "UCell",
                     "epithelial_proxy_unadjusted", "epithelial_proxy_residualized", "fixed_2x2_tiles"}
    if {r["sensitivity"] for r in sensitivity} != expected_sens:
        errors.append("sensitivity_family_invalid")
    if any(r["patient_id"] != "case4" or r["patient_n"] != "1" for r in sensitivity):
        errors.append("sensitivity_patient_boundary_invalid")
    if any(r["p_value"] != "NOT_COMPUTED_BY_DESIGN" or r["confidence_interval"] != "NOT_COMPUTED_BY_DESIGN" or r["fdr"] != "NOT_COMPUTED_BY_DESIGN" for r in sensitivity):
        errors.append("sensitivity_inference_present")

    primary = next(r for r in sensitivity if r["sensitivity"] == "primary_reference")
    pdelta = float(primary["Adenoma_minus_Normal"])
    pdir = direction(pdelta)
    mandatory = [r for r in sensitivity if r["sensitivity"] != "primary_reference"]
    any_na = any(r["direction"] == "not_estimable" for r in mandatory)
    reversal = any(r["direction"] != pdir for r in mandatory if r["direction"] != "not_estimable")
    expected_decision = "NOT_ESTIMABLE" if any_na else (
        "DESCRIPTIVE_METHOD_DEPENDENT" if pdir == "zero" or reversal else
        ("DESCRIPTIVE_POSITIVE" if pdir == "positive" else "DESCRIPTIVE_NEGATIVE")
    )
    decision_text = (out / "STAGE10E_DESC_DECISION.md").read_text(encoding="utf-8")
    decisions = [d for d in ALLOWED_DECISIONS if f"**{d}**" in decision_text]
    if decisions != [expected_decision]:
        errors.append(f"decision_mismatch:expected_{expected_decision}")

    for stem in ("Fig10E_DESC_1_case4_spatial", "Fig10E_DESC_2_case4_paired"):
        png, pdf = fig / f"{stem}.png", fig / f"{stem}.pdf"
        source = fig / "source_data" / f"{stem}_source_data.tsv"
        if not png.exists() or png.stat().st_size < 100 or png.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            errors.append(f"invalid_png:{stem}")
        if not pdf.exists() or pdf.stat().st_size < 100 or pdf.read_bytes()[:4] != b"%PDF":
            errors.append(f"invalid_pdf:{stem}")
        if not source.exists() or source.stat().st_size == 0:
            errors.append(f"missing_source_data:{stem}")

    for skip in (root / "results/stage10f/STAGE10F_SKIPPED.md", root / "results/stage10g/STAGE10G_SKIPPED.md"):
        if not skip.exists():
            errors.append(f"missing_skip_marker:{skip.name}")
    if errors:
        raise SystemExit(";".join(errors))
    print(f"decision={expected_decision}")
    print("patient_n=1")
    print(f"common_geneset_sha256={EXPECTED_HASH}")
    print("inference_computed=NO")


if __name__ == "__main__":
    main()
