#!/usr/bin/env python3
"""Build the Stage 11A immutable artifact inventory without rerunning analyses."""

from __future__ import annotations

import csv
import hashlib
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript_freeze" / "RESULTS_FREEZE_MANIFEST.tsv"
EVIDENCE_COMMIT = "e0f18f35aae957cf1b837ac4bc709a2a8d07cdc0"
REMOTE_ROOT = "${CRC_PROJECT_ROOT}"

REMOTE_ROWS = [
    ("Stage6A", "numeric_result", "main", "CLM03;CLM14", "results/06A_pseudobulk/pseudobulk_results.tsv", 36314540, "d217f8a76153ca2829a74c47f45fe97414c6e81e1189b4ba535f8effb81b7696", "Full donor-level pseudobulk result; large server-only table"),
    ("Stage6A", "environment", "reproducibility", "CLM03;CLM14", "results/06A_pseudobulk/software_versions.tsv", 119, "c849575d85f6f317d3b5d3306d8675249ac9aa52f0daff655e75ebd6bdb1d2ed", "Executed Stage 6A software versions"),
    ("Stage10C", "source_data", "main", "CLM04;CLM16", "results/10C_fap_confounding/20260731_110021/GSE201348_module_donor_effects.tsv", 5473, "914c154b608814c8a2fe2530a3d4db83597f5dda9cfe8214cccbb7d12c7243a1", "Four-donor paired values for Figure 3"),
    ("Stage10C", "numeric_result", "main", "CLM04;CLM16", "results/10C_fap_confounding/20260731_110021/stage10C_locked_module_results.tsv", 13195, "3d18119a4c82d397b12dbc7d7a4390e118ef9bba01289e2436a4233955da95f8", "Locked module effects, confidence intervals, P values and FDR"),
    ("Stage10C", "environment", "reproducibility", "CLM04;CLM16", "results/10C_fap_confounding/20260731_110021/stage10C_sessionInfo.txt", 3623, "c480844ac5e6ccc45d5e322c442771fa0657f6885105b0af0fc58d273407c3ea", "Executed Stage 10C R session"),
    ("Stage10C", "environment", "reproducibility", "CLM04;CLM16", "results/10C_fap_confounding/20260731_110021/stage10C_software_versions.tsv", 197, "0641b25a6f20f0ec7988f33ed8657327c6c1175fbaff54052b56492a71617b73", "Executed Stage 10C software versions"),
]

INCLUDE_DIRS = (
    "scripts",
    "config",
    "protocol",
    "governance",
    "metadata",
    "data/metadata",
    "logs_summary",
    "results_final",
    "figures_final",
    "figures",
    "reports",
    "results/05B_full_qc_integration",
    "results/05C_annotation",
    "results/stage10c",
    "results/stage10c2_sc",
    "results/stage10c2_sp",
    "results/stage10c3",
    "results/stage10d_tech",
    "results/stage10e",
    "results/stage10e_roi_remediation",
    "results/stage10e_desc",
    "results/stage10f",
    "results/stage10g",
    "results/stage10fg",
    "results/stage10h",
    "audit/stage10i",
)

INCLUDE_FILES = (
    "AGENTS.md",
    "PROJECT_CONSTITUTION.md",
    "SERVER_EXECUTION_CONTRACT.md",
    "STAGE10_PROTOCOL_ADDENDUM.md",
    "modules_locked.tsv",
    "modules_locked.tsv.sha256",
    "result_registry.tsv",
    "results/05B_qc_threshold_preflight.tsv",
    "STATUS.md",
    "manuscript_freeze/STAGE10I_R_NARRATIVE_ARCHITECTURE.md",
    "manuscript_freeze/STAGE10I_R_SECTION_CLAIM_MAP.tsv",
    "manuscript_freeze/STAGE10I_R_EVIDENCE_PLACEMENT.tsv",
    "manuscript_freeze/STAGE10I_R_DECISION.md",
    "manuscript_freeze/FINAL_CLAIM_LIST.tsv",
    "manuscript_freeze/FIGURE_BLUEPRINT.md",
    "manuscript_freeze/TABLE_BLUEPRINT.md",
    "manuscript_freeze/SUPPLEMENT_BLUEPRINT.md",
    "manuscript_freeze/SECTION_EVIDENCE_MAP.tsv",
    "scripts/11A_build_freeze_manifest.py",
)

CLAIM_MAP = {
    "logs_summary/stage_6A_key_metrics.tsv": "CLM03;CLM14",
    "logs_summary/stage_6A_candidate_attrition_audit.tsv": "CLM03;CLM14",
    "results_final/stage_6A_exploratory_candidate_modules.tsv": "CLM15",
    "results_final/stage_6A_stage_blind_module_membership.tsv": "CLM15",
    "results/stage10c/STAGE10C_LOCK_MANIFEST.tsv": "CLM12;CLM15",
    "results/stage10c/STAGE10C_STATISTICAL_AUDIT.tsv": "CLM04;CLM16",
    "results_final/stage_7_paired_donor_module_differences.tsv": "CLM05;CLM17",
    "results_final/stage_7_replication_effects_all.tsv": "CLM05;CLM17",
    "results_final/stage_7_replication_summary.tsv": "CLM20",
    "results_final/stage_8B_bulk_validation_summary.tsv": "CLM06;CLM19;CLM20",
    "results_final/stage_8B_meta_analysis_results.tsv": "CLM06;CLM19",
    "figures_final/stage_8B_early_transition_forest_source_data.tsv": "CLM06;CLM19",
    "results_final/stage_9C_stool_test_results.tsv": "CLM07;CLM21",
    "results_final/stage_9C_stool_test_predictions_source_data.tsv": "CLM21",
    "results/stage10c3/STAGE10C3_MATCHED_PAIR_EFFECTS.tsv": "CLM18",
    "results/stage10c3/STAGE10C3_SENSITIVITY.tsv": "CLM18",
    "results/stage10e_desc/STAGE10E_DESC_CASE4_SCORES.tsv": "CLM22",
    "results/stage10e_desc/STAGE10E_DESC_SENSITIVITY.tsv": "CLM22",
    "results/stage10fg/STAGE10FG_BRANCH_CLOSURE.tsv": "CLM07;CLM22",
    "results/stage10fg/STAGE10FG_DECISION.md": "CLM07;CLM22",
    "results/stage10h/STAGE10H_EVIDENCE_MATRIX.tsv": "CLM08;CLM23",
    "results/stage10h/STAGE10H_DATASET_INDEPENDENCE.tsv": "CLM11;CLM23",
    "results/stage10h/STAGE10H_SPATIAL_LIMITATION_MATRIX.tsv": "CLM22;CLM26",
    "audit/stage10i/STAGE10I_NUMERIC_RECONCILIATION.tsv": "CLM03;CLM04;CLM05;CLM06;CLM22;CLM23",
    "audit/stage10i/STAGE10I_FINAL_CLAIM_AUDIT.tsv": "CLM03;CLM04;CLM05;CLM06;CLM07;CLM08;CLM14;CLM15;CLM16;CLM17;CLM18;CLM19;CLM20;CLM21;CLM22;CLM23",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def stage_for(rel: str) -> str:
    low = rel.lower()
    for token in ("10i", "10h", "10fg", "10g", "10f", "10e", "10d", "10c3", "10c2", "10c", "9c", "9b", "9a", "8b", "8a", "7", "6b", "6a", "5c", "5b", "5a", "4b", "4a"):
        filename = Path(rel).name.lower()
        numbered_token = token if token.startswith("10") else "0" + token
        if f"stage{token}" in low or f"stage_{token}" in low or f"/{token}" in low or filename.startswith((token, numbered_token)):
            return f"Stage{token.upper()}"
    if rel.startswith(("governance/", "protocol/")):
        return "Governance"
    if rel.startswith("metadata/"):
        return "Metadata"
    return "Project"


def artifact_class(rel: str) -> str:
    low = rel.lower()
    suffix = Path(rel).suffix.lower()
    if "sessioninfo" in low or "software_versions" in low or "environment" in low:
        return "environment"
    if rel.startswith("scripts/"):
        return "analysis_script"
    if rel.startswith("config/"):
        return "parameter_file"
    if suffix in {".png", ".pdf"} and (rel.startswith("figures") or "/figures" in rel):
        return "figure_snapshot"
    if "source_data" in low or "source-data" in low:
        return "source_data"
    if suffix == ".tsv":
        return "tabular_result_or_audit"
    if rel.startswith(("governance/", "protocol/", "manuscript_freeze/")) or "decision" in low or "claim" in low:
        return "governance_or_claim_file"
    return "report_or_provenance"


def placement(rel: str, cls: str) -> str:
    low = rel.lower()
    if cls in {"analysis_script", "parameter_file", "environment"}:
        return "reproducibility"
    if any(x in low for x in ("stage10c3", "stage10d", "stage10e", "stage10f", "stage10g", "stage_9c", "stage9c", "stage_9a", "stage9a", "stage_9b", "stage9b")):
        return "supplement_or_limit"
    if any(x in low for x in ("stage_6a", "stage6a", "stage_7", "stage7", "stage_8b", "stage8b", "stage10c/", "stage10h", "stage10i")):
        return "main_or_governance"
    return "supporting_reproducibility"


def tracked_files() -> set[str]:
    out = subprocess.check_output(["git", "ls-files"], cwd=ROOT, text=True, encoding="utf-8")
    return {line.strip().replace("\\", "/") for line in out.splitlines() if line.strip()}


def selected(rel: str) -> bool:
    return rel in INCLUDE_FILES or any(rel == d or rel.startswith(d + "/") for d in INCLUDE_DIRS)


def main() -> None:
    tracked = tracked_files()
    rows: list[dict[str, object]] = []
    selected_paths = {rel for rel in tracked if selected(rel)}
    # Freeze small local analysis outputs even when .gitignore intentionally
    # keeps them out of Git. Raw-data and large-object directories are not in
    # INCLUDE_DIRS and therefore cannot enter this inventory accidentally.
    for directory in INCLUDE_DIRS:
        base = ROOT / directory
        if base.is_dir():
            for path in base.rglob("*"):
                if path.is_file():
                    selected_paths.add(path.relative_to(ROOT).as_posix())
    # Include newly created Stage 11A files before they are committed.
    for rel in INCLUDE_FILES:
        if (ROOT / rel).is_file():
            selected_paths.add(rel)
    for idx, rel in enumerate(sorted(selected_paths), start=1):
        path = ROOT / rel
        if not path.is_file():
            raise FileNotFoundError(rel)
        cls = artifact_class(rel)
        rows.append({
            "freeze_id": f"F{idx:04d}",
            "stage": stage_for(rel),
            "artifact_class": cls,
            "manuscript_placement": placement(rel, cls),
            "claim_ids": CLAIM_MAP.get(rel, "NA"),
            "storage": "local_git_worktree" if rel in tracked else "local_analysis_output_untracked_by_design",
            "path": rel,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
            "evidence_commit": EVIDENCE_COMMIT,
            "immutable_status": "FROZEN_NO_SILENT_OVERWRITE",
            "notes": "Inclusion freezes provenance; it does not by itself authorize a manuscript claim.",
        })
    start = len(rows) + 1
    for offset, (stage, cls, place, claims, rel, size, digest, note) in enumerate(REMOTE_ROWS):
        rows.append({
            "freeze_id": f"F{start + offset:04d}",
            "stage": stage,
            "artifact_class": cls,
            "manuscript_placement": place,
            "claim_ids": claims,
            "storage": "remote_server_read_only_reference",
            "path": f"{REMOTE_ROOT}/{rel}",
            "size_bytes": size,
            "sha256": digest,
            "evidence_commit": EVIDENCE_COMMIT,
            "immutable_status": "FROZEN_HASH_VERIFIED_2026-08-01",
            "notes": note,
        })
    fieldnames = [
        "freeze_id", "stage", "artifact_class", "manuscript_placement",
        "claim_ids", "storage", "path", "size_bytes", "sha256",
        "evidence_commit", "immutable_status", "notes",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()
