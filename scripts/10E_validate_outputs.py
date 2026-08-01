#!/usr/bin/env python3
"""Validate Stage 10E blinded QC artifacts without reading M02 outcomes."""
from __future__ import annotations
import csv, hashlib, sys
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage10e"
FIG = ROOT / "figures" / "stage10e"
PRIMARY = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
SENS = "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8"
COMMON = "4cd62d74b83673a4d2adf6077bedbdfe73d1cbd369a6a77418d124e0b506d482"

def rows(path: Path):
    with path.open(encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))

def hash_file(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    required = [
        ROOT / "data/metadata/spatial_patient_slide_roi_manifest.tsv",
        OUT / "STAGE10E_SPATIAL_QC.tsv", OUT / "STAGE10E_EXCLUSION_LOG.tsv",
        OUT / "STAGE10E_DECONVOLUTION_QC.tsv", OUT / "STAGE10E_GENESET_INTEGRITY_CHECK.tsv",
        OUT / "STAGE10E_ANALYSIS_READY_MANIFEST.tsv", OUT / "STAGE10E_RAW_HASHES.tsv",
        OUT / "STAGE10E_DECISION.md", OUT / "STAGE10E_SESSIONINFO.txt", OUT / "STAGE10E_RUN_MANIFEST.tsv",
        ROOT / "reports/STAGE10E_SUMMARY.md", ROOT / "reports/STAGE10E_GATE_DECISION.md",
    ]
    for path in required:
        assert path.is_file() and path.stat().st_size > 0, path
        assert b"\r\n" not in path.read_bytes(), path
    lock36 = rows(ROOT / "results/stage10c/STAGE10C_LOCK_MANIFEST.tsv")
    lock35 = rows(ROOT / "results/stage10c/M02_MINUS_INPP5D_SENS_V1.tsv")
    lock35_hash_text = (ROOT / "results/stage10c/M02_MINUS_INPP5D_SENS_V1_SHA256.txt").read_text(encoding="utf-8")
    assert len(lock36) == 36 and {x["bundle_sha256"] for x in lock36} == {PRIMARY}
    assert len(lock35) == 35 and "INPP5D" not in {x["gene"] for x in lock35} and SENS in lock35_hash_text
    integrity = rows(OUT / "STAGE10E_GENESET_INTEGRITY_CHECK.tsv")
    assert len(integrity) == 120 and {x["frozen_common_geneset_sha256"] for x in integrity} == {COMMON}
    assert all(x["present_in_feature_space"] == "TRUE" for x in integrity)
    assert all(x["expression_read_or_used_for_outcome"] == "FALSE" for x in integrity)
    qc = rows(OUT / "STAGE10E_SPATIAL_QC.tsv")
    deconv = rows(OUT / "STAGE10E_DECONVOLUTION_QC.tsv")
    assert len(qc) > 0 and len(qc) == len(deconv)
    assert {x["patient_id"] for x in qc} == {"case1", "case2", "case3", "case4"}
    assert all(x["marker_panel_excludes_all_36_M02_genes"] == "TRUE" for x in deconv)
    assert all(x["roi_assignment_expression_blind"] == "TRUE" for x in qc)
    decision = (OUT / "STAGE10E_DECISION.md").read_text(encoding="utf-8")
    allowed = [x for x in ("PASS", "PASS_WITH_LIMITATIONS", "FAIL") if f"Decision: **{x}**" in decision]
    assert len(allowed) == 1
    assert "no M02 expression score or region contrast was computed" in decision
    for stem in ("Fig10E_1_spot_qc", "Fig10E_2_deconvolution", "Fig10E_3_readiness"):
        assert (FIG / f"{stem}.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert (FIG / f"{stem}.pdf").read_bytes()[:4] == b"%PDF"
    for name in ("Fig10E_1_spot_qc_source_data.tsv", "Fig10E_2_reference_deconvolution_source_data.tsv", "Fig10E_3_analysis_ready_source_data.tsv"):
        assert len(rows(FIG / "source_data" / name)) > 0
    manifest = rows(OUT / "STAGE10E_RUN_MANIFEST.tsv")
    for item in manifest:
        p = ROOT / item["path"]
        assert p.is_file() and hash_file(p) == item["sha256"] and p.stat().st_size == int(float(item["bytes"])), p
    print("validation=PASS")
    print(f"decision={allowed[0]}")
    print("m02_lesion_normal_outcome_accessed=NO")

if __name__ == "__main__":
    main()
