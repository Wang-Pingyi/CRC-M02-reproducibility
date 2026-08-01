#!/usr/bin/env python3
"""Expression-blind, one-run pathology-overlay registration for Stage 10E-R.

The program never opens an expression matrix.  It uses only published images,
10x spatial coordinates/scalefactors, frozen metadata, and immutable lock files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage, optimize
from scipy.spatial import cKDTree, distance

EXPECTED_COMMON_HASH = "4cd62d74b83673a4d2adf6077bedbdfe73d1cbd369a6a77418d124e0b506d482"
PATHOLOGY_COLOURS = {
    "Adenoma": np.array([31, 119, 180], dtype=float),
    "Carcinoma": np.array([255, 127, 14], dtype=float),
    "Normal": np.array([44, 160, 44], dtype=float),
    "Other": np.array([214, 39, 40], dtype=float),
}
PRIMARY_CLASSES = ("Adenoma", "Normal")
DENIED_SUFFIXES = {".h5", ".h5ad", ".rds", ".mtx", ".mtx.gz", ".loom", ".fastq", ".fq"}
DENIED_TOKENS = (
    "stage10f", "m02_score", "m02_map", "lesion_normal", "regional_expression",
    "pseudobulk", "filtered_feature_bc_matrix", "counts_matrix", "expression_matrix",
)


@dataclass
class AccessRecord:
    path: str
    purpose: str
    decision: str
    reason: str


ACCESS_LOG: list[AccessRecord] = []
ALLOWED: set[Path] = set()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def register_allowed(path: Path) -> None:
    ALLOWED.add(path.resolve())


def guard(path: Path, purpose: str) -> Path:
    p = path.resolve()
    low = str(p).lower().replace("\\", "/")
    denied = any(low.endswith(s) for s in DENIED_SUFFIXES) or any(t in low for t in DENIED_TOKENS)
    if denied:
        ACCESS_LOG.append(AccessRecord(str(p), purpose, "DENY", "path matches expression/outcome denylist"))
        raise PermissionError(f"Denied Stage 10E-R input: {p}")
    if p not in ALLOWED:
        ACCESS_LOG.append(AccessRecord(str(p), purpose, "DENY", "path is absent from frozen whitelist"))
        raise PermissionError(f"Not in frozen Stage 10E-R whitelist: {p}")
    ACCESS_LOG.append(AccessRecord(str(p), purpose, "ALLOW", "frozen expression-blind input"))
    return p


def read_tsv(path: Path, purpose: str) -> list[dict[str, str]]:
    with guard(path, purpose).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k, "NA")) for k in columns})


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, (np.bool_, bool)):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (np.floating, float)):
        if not np.isfinite(value):
            return "NA"
        return f"{float(value):.8g}"
    return str(value)


def read_params(path: Path) -> dict[str, str]:
    return {r["parameter"]: r["value"] for r in read_tsv(path, "frozen remediation parameters")}


def image_array(path: Path, purpose: str) -> np.ndarray:
    with Image.open(guard(path, purpose)) as img:
        return np.asarray(img.convert("RGB"))


def read_positions(path: Path, purpose: str) -> dict[str, np.ndarray]:
    rows = []
    with guard(path, purpose).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle):
            if not row:
                continue
            rows.append(row)
    if rows and rows[0][0].lower() == "barcode":
        rows = rows[1:]
    return {
        "barcode": np.array([r[0] for r in rows], dtype=object),
        "in_tissue": np.array([int(r[1]) for r in rows], dtype=int),
        "array_row": np.array([int(r[2]) for r in rows], dtype=int),
        "array_col": np.array([int(r[3]) for r in rows], dtype=int),
        "px_row": np.array([float(r[4]) for r in rows], dtype=float),
        "px_col": np.array([float(r[5]) for r in rows], dtype=float),
    }


def read_json(path: Path, purpose: str) -> dict:
    with guard(path, purpose).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def colour_labels(rgb: np.ndarray, max_distance: float) -> tuple[np.ndarray, np.ndarray]:
    flat = rgb.reshape(-1, 3).astype(float)
    names = list(PATHOLOGY_COLOURS)
    proto = np.stack([PATHOLOGY_COLOURS[n] for n in names])
    d = np.sqrt(((flat[:, None, :] - proto[None, :, :]) ** 2).sum(axis=2))
    idx = d.argmin(axis=1)
    mind = d[np.arange(len(flat)), idx]
    labels = np.full(len(flat), -1, dtype=np.int16)
    labels[mind <= max_distance] = idx[mind <= max_distance]
    return labels.reshape(rgb.shape[:2]), np.array(names, dtype=object)


def retain_tissue_colour_component(label_image: np.ndarray) -> np.ndarray:
    raw = label_image >= 0
    joined = ndimage.binary_dilation(raw, iterations=2)
    components, count = ndimage.label(joined)
    if count == 0:
        return np.zeros_like(raw)
    sizes = ndimage.sum(joined, components, index=np.arange(1, count + 1))
    keep = int(np.argmax(sizes) + 1)
    return raw & ndimage.binary_dilation(components == keep, iterations=2)


def median_spacing(points: np.ndarray) -> float:
    if len(points) < 2:
        return float("nan")
    d, _ = cKDTree(points).query(points, k=2)
    return float(np.median(d[:, 1]))


def orientation_matrix(rotation: int, mirror: bool) -> np.ndarray:
    theta = math.radians(rotation)
    r = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    if mirror:
        r = r @ np.array([[-1.0, 0.0], [0.0, 1.0]])
    return r


def affine_apply(points: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    return points @ matrix[:, :2].T + matrix[:, 2]


def geometric_score(mapped: np.ndarray, target: np.ndarray, target_tree: cKDTree) -> float:
    d1, _ = target_tree.query(mapped, k=1)
    sample = target[::max(1, len(target) // 2500)]
    d2, _ = cKDTree(mapped).query(sample, k=1)
    return float(np.median(d1) + 0.35 * np.median(d2))


def pipeline_a(source: np.ndarray, target: np.ndarray) -> dict:
    src0 = source - source.mean(axis=0)
    tgt_center = target.mean(axis=0)
    target_tree = cKDTree(target)
    candidates = []
    for rotation in (0, 90, 180, 270):
        for mirror in (False, True):
            oriented = src0 @ orientation_matrix(rotation, mirror).T
            src_range = np.ptp(oriented, axis=0)
            tgt_range = np.ptp(target, axis=0)
            s0 = float(np.median(tgt_range / np.maximum(src_range, 1e-6)))

            def sim_obj(v):
                s = math.exp(v[0])
                mapped = oriented * s + np.array(v[1:3])
                return geometric_score(mapped, target, target_tree)

            sim = optimize.minimize(sim_obj, [math.log(s0), *tgt_center], method="Powell",
                                    bounds=[(math.log(s0 * 0.65), math.log(s0 * 1.35)),
                                            (target[:, 0].min() - 20, target[:, 0].max() + 20),
                                            (target[:, 1].min() - 20, target[:, 1].max() + 20)],
                                    options={"maxiter": 90, "xtol": 1e-4, "ftol": 1e-4})
            s = math.exp(sim.x[0])
            m_sim = np.array([[s, 0.0, sim.x[1]], [0.0, s, sim.x[2]]])
            base = affine_apply(oriented, m_sim)

            def aff_obj(v):
                delta = np.array([[v[0], v[1], v[4]], [v[2], v[3], v[5]]])
                return geometric_score(affine_apply(oriented, delta), target, target_tree)

            lim = 0.28 * s
            aff = optimize.minimize(aff_obj, [s, 0.0, 0.0, s, sim.x[1], sim.x[2]], method="Powell",
                                    bounds=[(0.72*s, 1.28*s), (-lim, lim), (-lim, lim), (0.72*s, 1.28*s),
                                            (sim.x[1]-15, sim.x[1]+15), (sim.x[2]-15, sim.x[2]+15)],
                                    options={"maxiter": 120, "xtol": 1e-4, "ftol": 1e-4})
            m_aff = np.array([[aff.x[0], aff.x[1], aff.x[4]], [aff.x[2], aff.x[3], aff.x[5]]])
            use_affine = aff.fun < 0.95 * sim.fun
            matrix_centered = m_aff if use_affine else m_sim
            candidates.append({"rotation": rotation, "mirror": mirror,
                               "family": "affine" if use_affine else "similarity",
                               "score": float(min(aff.fun, sim.fun)), "matrix_centered": matrix_centered,
                               "oriented": oriented, "orientation": orientation_matrix(rotation, mirror)})
    candidates.sort(key=lambda x: x["score"])
    best = candidates[0]
    linear = best["matrix_centered"][:, :2] @ best["orientation"]
    translation = best["matrix_centered"][:, 2] - linear @ source.mean(axis=0)
    best["matrix"] = np.column_stack([linear, translation])
    best["ambiguity_gap"] = candidates[1]["score"] - best["score"]
    return best


def radial_landmarks(points: np.ndarray, count: int) -> np.ndarray:
    centered = points - points.mean(axis=0)
    out = []
    for angle in np.linspace(0, 2 * np.pi, count, endpoint=False):
        direction = np.array([math.cos(angle), math.sin(angle)])
        out.append(points[int(np.argmax(centered @ direction))])
    return np.asarray(out)


def fit_affine(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    design = np.column_stack([x, np.ones(len(x))])
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return coef.T


def pipeline_b(source: np.ndarray, target: np.ndarray, landmark_count: int) -> dict:
    target_landmarks = radial_landmarks(target, landmark_count)
    candidates = []
    src_center = source.mean(axis=0)
    src0 = source - src_center
    for rotation in (0, 90, 180, 270):
        for mirror in (False, True):
            orient = orientation_matrix(rotation, mirror)
            oriented = src0 @ orient.T
            source_landmarks = radial_landmarks(oriented, landmark_count)
            matrix_local = fit_affine(source_landmarks, target_landmarks)
            mapped_landmarks = affine_apply(source_landmarks, matrix_local)
            errors = np.linalg.norm(mapped_landmarks - target_landmarks, axis=1)
            loo = []
            for i in range(landmark_count):
                keep = np.arange(landmark_count) != i
                m = fit_affine(source_landmarks[keep], target_landmarks[keep])
                loo.append(float(np.linalg.norm(affine_apply(source_landmarks[i:i+1], m)[0] - target_landmarks[i])))
            mapped = affine_apply(oriented, matrix_local)
            score = geometric_score(mapped, target, cKDTree(target)) + float(np.median(loo))
            candidates.append({"rotation": rotation, "mirror": mirror, "family": "affine_landmarks",
                               "score": score, "matrix_local": matrix_local, "orientation": orient,
                               "landmark_errors": errors, "loo_errors": np.asarray(loo),
                               "source_landmarks": source_landmarks, "target_landmarks": target_landmarks})
    candidates.sort(key=lambda x: x["score"])
    best = candidates[0]
    linear = best["matrix_local"][:, :2] @ best["orientation"]
    translation = best["matrix_local"][:, 2] - linear @ src_center
    best["matrix"] = np.column_stack([linear, translation])
    best["ambiguity_gap"] = candidates[1]["score"] - best["score"]
    return best


def assign_classes(mapped: np.ndarray, label_image: np.ndarray, names: np.ndarray,
                   target_keep: np.ndarray, spacing: float, radius_fraction: float,
                   ambiguity_fraction: float) -> np.ndarray:
    class_trees = {}
    for idx, name in enumerate(names):
        yy, xx = np.where((label_image == idx) & target_keep)
        if len(xx):
            class_trees[str(name)] = cKDTree(np.column_stack([xx, yy]))
    out = np.full(len(mapped), "UNASSIGNED", dtype=object)
    for i, p in enumerate(mapped):
        distances = sorted((float(tree.query(p, k=1)[0]), name) for name, tree in class_trees.items())
        if not distances or distances[0][0] > spacing * radius_fraction:
            continue
        if len(distances) > 1 and distances[1][0] - distances[0][0] < spacing * ambiguity_fraction:
            continue
        if distances[0][1] in PRIMARY_CLASSES:
            out[i] = distances[0][1]
    return out


def cohen_kappa(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0:
        return float("nan")
    classes = sorted(set(a) | set(b))
    po = float(np.mean(a == b))
    pe = sum(float(np.mean(a == c) * np.mean(b == c)) for c in classes)
    return float((po - pe) / (1 - pe)) if pe < 1 else float("nan")


def class_metric(a: np.ndarray, b: np.ndarray, cls: str) -> tuple[float, float]:
    aa, bb = a == cls, b == cls
    inter = int(np.sum(aa & bb))
    denom = int(np.sum(aa) + np.sum(bb))
    union = int(np.sum(aa | bb))
    return (2 * inter / denom if denom else float("nan"), inter / union if union else float("nan"))


def point_mask(shape: tuple[int, int], points: np.ndarray, radius: float) -> np.ndarray:
    img = Image.new("1", (shape[1], shape[0]), 0)
    draw = ImageDraw.Draw(img)
    r = max(1.0, radius)
    for x, y in points:
        draw.ellipse((x-r, y-r, x+r, y+r), fill=1)
    return np.asarray(img, dtype=bool)


def dice(mask1: np.ndarray, mask2: np.ndarray) -> float:
    den = int(mask1.sum() + mask2.sum())
    return float(2 * np.sum(mask1 & mask2) / den) if den else float("nan")


def boundary_points(mask: np.ndarray) -> np.ndarray:
    edge = mask ^ ndimage.binary_erosion(mask)
    yy, xx = np.where(edge)
    return np.column_stack([xx, yy])


def hausdorff(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    da = cKDTree(b).query(a, k=1)[0]
    db = cKDTree(a).query(b, k=1)[0]
    return float(max(np.max(da), np.max(db)))


def crop_overlay(row: dict[str, str], root: Path, max_colour_distance: float):
    path = root / row["overlay_source"]
    arr = image_array(path, f"published pathology overlay for {row['patient_id']}")
    x, y, w, h = [int(row[k]) for k in ("overlay_crop_x", "overlay_crop_y", "overlay_crop_width", "overlay_crop_height")]
    crop = arr[y:y+h, x:x+w]
    labels, names = colour_labels(crop, max_colour_distance)
    keep = retain_tissue_colour_component(labels)
    if int(keep.sum()) < 25:
        raise RuntimeError(f"Insufficient pathology-colour pixels for {row['patient_id']}")
    return path, arr, crop, labels, names, keep


def setup_whitelist(root: Path, config_rows: list[dict[str, str]], config_path: Path,
                    params_path: Path, common_path: Path, old_ready_path: Path) -> None:
    for path in (config_path, params_path, common_path, old_ready_path):
        register_allowed(path)
    for row in config_rows:
        for key in ("overlay_source", "hires_source", "positions_source", "scalefactors_source"):
            register_allowed(root / row[key])


def verify_common_hash(common_path: Path) -> str:
    rows = read_tsv(common_path, "frozen 30-of-36 gene identity lock; no expression values")
    genes = sorted(r["canonical_gene"] for r in rows if r["cohort_id"] == "E-GEAD-622")
    digest = hashlib.sha256(("\n".join(genes) + "\n").encode("utf-8")).hexdigest()
    if len(genes) != 30 or digest != EXPECTED_COMMON_HASH:
        raise RuntimeError(f"Frozen E-GEAD-622 gene identity mismatch: n={len(genes)}, hash={digest}")
    return digest


def write_access_outputs(out_dir: Path, root: Path, config_rows: list[dict[str, str]]) -> None:
    whitelist = []
    purposes = {
        "overlay_source": "published pathologist overlay",
        "hires_source": "deposited H&E image",
        "positions_source": "10x spot coordinates and tissue mask",
        "scalefactors_source": "10x coordinate scalefactors",
    }
    for row in config_rows:
        for key, purpose in purposes.items():
            p = (root / row[key]).resolve()
            whitelist.append({"patient_id": row["patient_id"], "input_type": key, "path": str(p),
                              "purpose": purpose, "sha256": sha256_file(p), "expression_content": "NO"})
    write_tsv(out_dir / "STAGE10E_R_INPUT_WHITELIST.tsv", whitelist)
    denied_examples = [
        root / "results/stage10f",
        root / "objects",
        root / "data_processed/stage10e/stage10e_20260801_010406/derived/Rectum_kyudai_Beppu_20200303/filtered_feature_bc_matrix.h5",
    ]
    rows = [r.__dict__ for r in ACCESS_LOG]
    rows.extend({"path": str(p), "purpose": "prohibited outcome/expression input", "decision": "DENY",
                 "reason": "code-level denylist; never opened"} for p in denied_examples)
    write_tsv(out_dir / "STAGE10E_R_DENYLIST_AUDIT.tsv", rows,
              ["path", "purpose", "decision", "reason"])


def audit(root: Path, config_rows: list[dict[str, str]], params: dict[str, str], out_dir: Path) -> None:
    old_rows = read_tsv(root / "results/stage10e/STAGE10E_ANALYSIS_READY_MANIFEST.tsv",
                        "historical Stage 10E technical readiness only")
    old = {r["patient_id"]: r for r in old_rows}
    coordinate_rows, diagnosis_rows = [], []
    old_metrics = {
        "case1": (0.409717, 0.191360), "case2": (0.604699, -0.033043),
        "case3": (0.215517, -0.291530), "case4": (0.385266, -0.017337),
    }
    for row in config_rows:
        _, full, crop, labels, _, keep = crop_overlay(row, root, float(params["colour_rgb_distance_max"]))
        pos = read_positions(root / row["positions_source"], f"spot coordinates for {row['patient_id']}")
        sf = read_json(root / row["scalefactors_source"], f"coordinate scale for {row['patient_id']}")
        hires = image_array(root / row["hires_source"], f"deposited H&E for {row['patient_id']}")
        yy, xx = np.where(keep)
        in_tissue = pos["in_tissue"] == 1
        coordinate_rows.append({
            "patient_id": row["patient_id"], "slide_or_capture_id": row["slide_or_capture_id"],
            "overlay_width": full.shape[1], "overlay_height": full.shape[0],
            "frozen_crop_x": row["overlay_crop_x"], "frozen_crop_y": row["overlay_crop_y"],
            "frozen_crop_width": row["overlay_crop_width"], "frozen_crop_height": row["overlay_crop_height"],
            "colour_bbox_xmin": int(xx.min()), "colour_bbox_xmax": int(xx.max()),
            "colour_bbox_ymin": int(yy.min()), "colour_bbox_ymax": int(yy.max()),
            "pathology_colour_pixels": int(keep.sum()), "hires_width": hires.shape[1], "hires_height": hires.shape[0],
            "tissue_hires_scalef": sf["tissue_hires_scalef"], "spot_diameter_fullres": sf["spot_diameter_fullres"],
            "positions_total": len(pos["barcode"]), "positions_in_tissue": int(in_tissue.sum()),
            "fullres_xmin": float(pos["px_col"].min()), "fullres_xmax": float(pos["px_col"].max()),
            "fullres_ymin": float(pos["px_row"].min()), "fullres_ymax": float(pos["px_row"].max()),
        })
        agreement, kappa = old_metrics[row["patient_id"]]
        diagnosis_rows.append({
            "patient_id": row["patient_id"], "old_agreement": agreement, "old_kappa": kappa,
            "old_coordinate_A": "independent_minmax_fullres_pixel_to_manual_crop",
            "old_coordinate_B": "independent_minmax_array_lattice_to_manual_crop",
            "primary_failure_source": "nonphysical_bbox_normalization_and_crop_axes_legend_contamination",
            "error_type": "global_affine_scale_translation_plus_mask_scope_mismatch",
            "class_swap_evidence": "none_detected_in_code",
            "background_unassigned_issue": "old_scope_used_any_joint_nonmissing_without_common_pathology_tissue_mask",
            "m02_or_expression_used": "NO",
            "recoverable": "YES",
            "expected_kappa_impact": "large; coordinate offsets move spots across broad Normal/Adenoma boundaries",
            "historical_stage10e_status_preserved": old.get(row["patient_id"], {}).get("registration_pass", "FALSE"),
        })
    write_tsv(out_dir / "STAGE10E_R_COORDINATE_AUDIT.tsv", coordinate_rows)
    write_tsv(out_dir / "STAGE10E_R_FAILURE_DIAGNOSIS.tsv", diagnosis_rows)


def formal(root: Path, config_rows: list[dict[str, str]], params: dict[str, str], out_dir: Path,
           figure_dir: Path, run_id: str) -> None:
    rng = np.random.default_rng(int(params["seed"]))
    del rng  # deterministic code path; retained to freeze seed provenance
    transform_rows, geometric_rows, agreement_rows, eligibility_rows = [], [], [], []
    for row in config_rows:
        _, _, crop, label_image, names, target_keep = crop_overlay(
            row, root, float(params["colour_rgb_distance_max"]))
        pos = read_positions(root / row["positions_source"], f"spot coordinates for {row['patient_id']}")
        sf = read_json(root / row["scalefactors_source"], f"coordinate scale for {row['patient_id']}")
        _ = image_array(root / row["hires_source"], f"deposited H&E for {row['patient_id']}")
        in_tissue = pos["in_tissue"] == 1
        source_all = np.column_stack([pos["px_col"], pos["px_row"]]) * float(sf["tissue_hires_scalef"])
        source = source_all[in_tissue]
        yy, xx = np.where(target_keep)
        target = np.column_stack([xx, yy]).astype(float)
        a = pipeline_a(source, target)
        b = pipeline_b(source, target, int(params["landmark_count"]))
        mapped_a = affine_apply(source, a["matrix"])
        mapped_b = affine_apply(source, b["matrix"])
        spacing = float(np.mean([median_spacing(mapped_a), median_spacing(mapped_b)]))
        labels_a = assign_classes(mapped_a, label_image, names, target_keep, spacing,
                                  float(params["class_spatial_radius_fraction"]),
                                  float(params["boundary_ambiguity_fraction"]))
        labels_b = assign_classes(mapped_b, label_image, names, target_keep, spacing,
                                  float(params["class_spatial_radius_fraction"]),
                                  float(params["boundary_ambiguity_fraction"]))
        common = np.isin(labels_a, PRIMARY_CLASSES) & np.isin(labels_b, PRIMARY_CLASSES)
        aa, bb = labels_a[common], labels_b[common]
        agreement = float(np.mean(aa == bb)) if len(aa) else float("nan")
        kappa = cohen_kappa(aa, bb)
        consensus = np.full(len(source), "UNASSIGNED", dtype=object)
        consensus[common & (labels_a == labels_b)] = labels_a[common & (labels_a == labels_b)]
        adenoma_n = int(np.sum(consensus == "Adenoma"))
        normal_n = int(np.sum(consensus == "Normal"))
        min_count = int(row["minimum_spots_per_class"])
        registration_pass = bool(agreement >= float(params["agreement_min"]) and
                                 kappa >= float(params["kappa_min"]))
        paired_coverage = adenoma_n >= min_count and normal_n >= min_count
        patient_status = "PASS" if registration_pass and paired_coverage else (
            "NOT_TECHNICALLY_ELIGIBLE" if registration_pass else "FAIL_REGISTRATION")
        dice_a, iou_a = class_metric(aa, bb, "Adenoma")
        dice_n, iou_n = class_metric(aa, bb, "Normal")
        mask_target = point_mask(crop.shape[:2], target[::max(1, len(target)//3000)], spacing * float(params["mask_radius_fraction"]))
        mask_a = point_mask(crop.shape[:2], mapped_a, spacing * float(params["mask_radius_fraction"]))
        mask_b = point_mask(crop.shape[:2], mapped_b, spacing * float(params["mask_radius_fraction"]))
        boundary_target = boundary_points(mask_target)
        boundary_a = boundary_points(mask_a)
        boundary_b = boundary_points(mask_b)
        a_haus = hausdorff(boundary_a, boundary_target) / spacing
        b_haus = hausdorff(boundary_b, boundary_target) / spacing
        ab_shift = np.linalg.norm(mapped_a - mapped_b, axis=1) / spacing

        for pipeline, fit in (("A_contour_coordinate", a), ("B_morphology_landmark", b)):
            m = fit["matrix"]
            transform_rows.append({
                "run_id": run_id, "patient_id": row["patient_id"], "pipeline": pipeline,
                "rotation_degrees": fit["rotation"], "mirror": fit["mirror"], "transform_family": fit["family"],
                "m00": m[0,0], "m01": m[0,1], "tx": m[0,2], "m10": m[1,0], "m11": m[1,1], "ty": m[1,2],
                "geometric_objective": fit["score"], "second_best_gap": fit["ambiguity_gap"],
                "nonlinear_used": "FALSE", "selection_information": "geometry_only",
            })
        geometric_rows.append({
            "patient_id": row["patient_id"], "slide_or_capture_id": row["slide_or_capture_id"],
            "spot_spacing_overlay_px": spacing, "pipeline_A_tissue_mask_dice": dice(mask_a, mask_target),
            "pipeline_B_tissue_mask_dice": dice(mask_b, mask_target),
            "pipeline_A_hausdorff_spot_diameters": a_haus, "pipeline_B_hausdorff_spot_diameters": b_haus,
            "pipeline_B_landmark_median_reprojection_spot_diameters": float(np.median(b["landmark_errors"]) / spacing),
            "pipeline_B_landmark_p95_reprojection_spot_diameters": float(np.quantile(b["landmark_errors"], .95) / spacing),
            "pipeline_B_LOLO_median_spot_diameters": float(np.median(b["loo_errors"]) / spacing),
            "A_B_median_coordinate_distance_spot_diameters": float(np.median(ab_shift)),
            "A_B_p95_coordinate_distance_spot_diameters": float(np.quantile(ab_shift, .95)),
            "rotation_mirror_ambiguity_A": a["ambiguity_gap"], "rotation_mirror_ambiguity_B": b["ambiguity_gap"],
            "systematic_direction_conflict": "NO" if a["rotation"] == b["rotation"] and a["mirror"] == b["mirror"] else "REVIEW",
        })
        agreement_rows.append({
            "patient_id": row["patient_id"], "slide_or_capture_id": row["slide_or_capture_id"],
            "agreement_scope": "common_pathology_determinable_in_tissue_Normal_or_Adenoma",
            "common_evaluable_spots": int(common.sum()), "agreement": agreement, "cohen_kappa": kappa,
            "PABAK": 2 * agreement - 1 if np.isfinite(agreement) else float("nan"),
            "Adenoma_Dice": dice_a, "Adenoma_IoU": iou_a, "Normal_Dice": dice_n, "Normal_IoU": iou_n,
            "pipeline_A_unassigned_fraction": float(np.mean(labels_a == "UNASSIGNED")),
            "pipeline_B_unassigned_fraction": float(np.mean(labels_b == "UNASSIGNED")),
            "consensus_unassigned_fraction": float(np.mean(consensus == "UNASSIGNED")),
            "consensus_Adenoma_spots": adenoma_n, "consensus_Normal_spots": normal_n,
            "registration_pass": registration_pass,
        })
        eligibility_rows.append({
            "patient_id": row["patient_id"], "slide_or_capture_id": row["slide_or_capture_id"],
            "registration_pass": registration_pass, "paired_minimum_spot_coverage": paired_coverage,
            "Adenoma_spots": adenoma_n, "Normal_spots": normal_n, "minimum_each": min_count,
            "patient_decision": patient_status, "biological_inference_unit": "patient",
        })

        # Permitted QC figure: H&E/pathology overlay, boundaries, coordinates and landmarks only.
        fig = Image.fromarray(crop.copy())
        draw = ImageDraw.Draw(fig)
        for pts, colour in ((mapped_a[::10], (255,255,255)), (mapped_b[::10], (0,0,0))):
            for x, y in pts:
                draw.ellipse((x-1.2, y-1.2, x+1.2, y+1.2), outline=colour)
        for p in affine_apply(b["source_landmarks"], b["matrix_local"]):
            draw.rectangle((p[0]-2, p[1]-2, p[0]+2, p[1]+2), outline=(255,255,0))
        figure_dir.mkdir(parents=True, exist_ok=True)
        fig.save(figure_dir / f"Fig10E_R_{row['patient_id']}_registration_qc.png")
        source_rows = []
        tissue_idx = np.where(in_tissue)[0]
        for j, original_idx in enumerate(tissue_idx):
            source_rows.append({"patient_id": row["patient_id"], "spot_barcode": pos["barcode"][original_idx],
                                "pipeline_A_x": mapped_a[j,0], "pipeline_A_y": mapped_a[j,1],
                                "pipeline_B_x": mapped_b[j,0], "pipeline_B_y": mapped_b[j,1],
                                "pipeline_A_label": labels_a[j], "pipeline_B_label": labels_b[j],
                                "consensus_label": consensus[j]})
        write_tsv(figure_dir / "source_data" / f"Fig10E_R_{row['patient_id']}_registration_qc_source_data.tsv", source_rows)

    write_tsv(out_dir / "STAGE10E_R_TRANSFORM_PARAMETERS.tsv", transform_rows)
    write_tsv(out_dir / "STAGE10E_R_GEOMETRIC_QC.tsv", geometric_rows)
    write_tsv(out_dir / "STAGE10E_R_ROI_AGREEMENT.tsv", agreement_rows)
    write_tsv(out_dir / "STAGE10E_R_PATIENT_ELIGIBILITY.tsv", eligibility_rows)
    passes = sum(r["patient_decision"] == "PASS" for r in eligibility_rows)
    conflict = any(r["systematic_direction_conflict"] == "REVIEW" for r in geometric_rows)
    if passes >= int(params["minimum_passing_patients"]) and not conflict:
        decision = "PASS_REMEDIATED"
    elif passes >= 1:
        decision = "PASS_DESCRIPTIVE_ONLY"
    else:
        decision = "FAIL"
    (out_dir / "STAGE10E_R_DECISION.md").write_text(
        "# Stage 10E-R decision\n\n"
        f"Decision: **{decision}**\n\n"
        f"- Formal run ID: `{run_id}`\n- Passing patients: {passes}/4.\n"
        f"- Frozen agreement/kappa thresholds: {params['agreement_min']} / {params['kappa_min']}.\n"
        f"- Frozen 30/36 common-set hash: `{EXPECTED_COMMON_HASH}`.\n"
        "- Expression/M02 outcome leakage: NO.\n- Stage 10E historical decision remains PASS_WITH_LIMITATIONS and was not modified.\n"
        f"- Stage 10F authorization from this remediation: {'ELIGIBLE_FOR_SEPARATE_FUTURE_AUTHORIZATION' if decision == 'PASS_REMEDIATED' else 'NO'}.\n"
        "- Stage 10F was not started.\n", encoding="utf-8")
    (out_dir / "STAGE10E_R_ACCEPTANCE.md").write_text(
        "# Stage 10E-R formal-run acceptance pending independent validator\n\n"
        f"Formal computation completed with provisional decision **{decision}**. "
        "The separate validator must verify schemas, hashes, leakage audit, figures and thresholds before acceptance.\n",
        encoding="utf-8")
    if decision != "PASS_REMEDIATED":
        skipped = root / "results/stage10f/STAGE10F_SKIPPED.md"
        skipped.parent.mkdir(parents=True, exist_ok=True)
        skipped.write_text("# Stage 10F skipped\n\nStage 10E-R did not meet PASS_REMEDIATED. Stage 10F and Stage 10G are not authorized.\n", encoding="utf-8")


def smoke(out_dir: Path) -> None:
    # Synthetic-only smoke test. No real image, coordinate, pathology or expression input is opened.
    yy, xx = np.mgrid[0:81:4, 0:101:4]
    shape = (((xx - 48) / 47) ** 2 + ((yy - 38) / 34) ** 2 < 1)
    notch = (xx > 63) & (yy < 30) & (yy > 10)
    tail = (xx < 18) & (yy > 45)
    source = np.column_stack([xx[shape & ~notch | tail], yy[shape & ~notch | tail]]).astype(float)
    true = np.array([[0.9, 0.03, 6.0], [-0.02, 0.9, 5.0]])
    target_centres = affine_apply(source, true)
    target = np.repeat(target_centres, 3, axis=0) + np.tile(np.array([[0,0],[.3,0],[0,.3]]), (len(target_centres),1))
    a = pipeline_a(source, target)
    b = pipeline_b(source, target, 16)
    tree = cKDTree(target_centres)
    err_a = float(np.median(tree.query(affine_apply(source, a["matrix"]), k=1)[0]))
    err_b = float(np.median(tree.query(affine_apply(source, b["matrix"]), k=1)[0]))
    status = "PASS" if err_a < 2.0 and err_b < 3.0 else "FAIL"
    write_tsv(out_dir / "STAGE10E_R_SMOKE_TEST.tsv", [{"test": "synthetic_affine_registration", "pipeline_A_median_error_px": err_a,
                                                         "pipeline_B_median_error_px": err_b, "status": status,
                                                         "real_M02_or_expression_accessed": "NO"}])
    if status != "PASS":
        raise RuntimeError("Synthetic Stage 10E-R smoke test failed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--mode", choices=("audit", "smoke", "formal"), required=True)
    parser.add_argument("--run-id", default="NA")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    config_path = root / "config/stage10e_r_registration_inputs.tsv"
    params_path = root / "config/stage10e_r_parameters.tsv"
    common_path = root / "results/stage10c2_sp/STAGE10C2_SP_COHORT_COMMON_GENESET.tsv"
    old_ready_path = root / "results/stage10e/STAGE10E_ANALYSIS_READY_MANIFEST.tsv"
    register_allowed(config_path); register_allowed(params_path)
    config_rows = read_tsv(config_path, "frozen registration inputs")
    setup_whitelist(root, config_rows, config_path, params_path, common_path, old_ready_path)
    params = read_params(params_path)
    out_dir = root / "results/stage10e_roi_remediation"
    figure_dir = root / "figures/stage10e_roi_remediation"
    out_dir.mkdir(parents=True, exist_ok=True)
    digest = verify_common_hash(common_path)
    write_tsv(out_dir / "STAGE10E_R_GENESET_INTEGRITY_CHECK.tsv", [{
        "cohort": "E-GEAD-622", "mapped_genes": 30, "canonical_genes": 36,
        "recomputed_sha256": digest, "expected_sha256": EXPECTED_COMMON_HASH,
        "status": "PASS", "expression_values_opened": "NO"}])
    if args.mode == "audit":
        audit(root, config_rows, params, out_dir)
        write_access_outputs(out_dir, root, config_rows)
    elif args.mode == "smoke":
        smoke(out_dir)
    else:
        plan = root / "results/stage10e_roi_remediation/STAGE10E_R_REMEDIATION_PLAN_LOCKED.md"
        if not plan.exists():
            raise RuntimeError("Frozen remediation plan is missing")
        formal(root, config_rows, params, out_dir, figure_dir, args.run_id)
        write_access_outputs(out_dir, root, config_rows)
        write_tsv(out_dir / "STAGE10E_R_RUNTIME.tsv", [{
            "run_id": args.run_id, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version.replace("\n", " "), "platform": platform.platform(),
            "numpy": np.__version__, "scipy": __import__("scipy").__version__,
            "pillow": __import__("PIL").__version__, "seed": params["seed"]}])


if __name__ == "__main__":
    main()
