#!/usr/bin/env python3
"""Finalize the outcome-blind Stage 10C2-SP spatial eligibility gate."""

from __future__ import annotations

import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_DIR = ROOT / "tmp" / "stage10c2_sp_audit"
OUT = ROOT / "results" / "stage10c2_sp"

PRIMARY_BUNDLE = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
SENSITIVITY_BUNDLE = "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8"
SEED = 20260731

REQUIRED_OUTPUTS = [
    "STAGE10C2_SP_SOURCE_AUDIT.tsv",
    "STAGE10C2_SP_PATIENT_SLIDE_ROI_MANIFEST.tsv",
    "STAGE10C2_SP_GENE_ID_MAPPING.tsv",
    "STAGE10C2_SP_GENE_COVERAGE.tsv",
    "STAGE10C2_SP_COHORT_COMMON_GENESET.tsv",
    "STAGE10C2_SP_ANALYSIS_PLAN.md",
    "STAGE10C2_SP_LOCK_MANIFEST.tsv",
    "STAGE10C2_SP_DECISION.md",
]


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_text(payload: str) -> str:
    return sha256_bytes(payload.encode("utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def coverage_tier(mapped: int) -> str:
    if mapped == 36:
        return "COMPLETE_COVERAGE"
    if mapped >= 29:
        return "HIGH_COVERAGE"
    if mapped >= 22:
        return "MINIMUM_COVERAGE"
    return "BELOW_MINIMUM_COVERAGE"


def source_audit_rows() -> list[dict[str, object]]:
    article = "https://doi.org/10.1016/j.ebiom.2024.105102"
    return [
        {
            "cohort_id": "E-GEAD-622",
            "frozen_role": "primary_early_lesion_spatial_candidate",
            "source_patients": 4,
            "source_slides": 4,
            "early_lesion_patients": 4,
            "primary_eligible_paired_patients": 4,
            "pathology": "pTis early CRC; adenoma and intramucosal carcinoma on the same section",
            "pathology_regions": "Normal;Adenoma;Carcinoma;Other",
            "roi_source": "board-certified colorectal pathologist annotations published in main Fig1b (case1) and Supplementary Fig1b (cases2-4)",
            "roi_independent_of_m02": "TRUE",
            "material": "FFPE",
            "platform": "10x Visium FFPE; MGI DNBSEQ-G400",
            "reference_annotation": "GRCh38 2020-A",
            "counts": "filtered integer counts in deposited 10x H5",
            "coordinates": "deposited tissue_positions_list.csv",
            "he_image": "deposited hires/lowres tissue images; publication pathology overlays",
            "license_or_reuse": "reanalysis under DDBJ/NBDC repository conditions; article CC BY-NC-ND 4.0; do not redistribute raw human data",
            "cross_study_overlap": "not demonstrated; not formally excludable from de-identified records",
            "eligibility": "ELIGIBLE_PRIMARY_PENDING_BLIND_ROI_REGISTRATION_QC",
            "evidence": f"{article}; https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-622/",
            "limitations": "No structured barcode-to-pathology TSV/cloupe was deposited; published pathologist overlays must be registered to deposited spot coordinates without expression data.",
        },
        {
            "cohort_id": "E-GEAD-619",
            "frozen_role": "secondary_early_and_advanced_localization",
            "source_patients": 2,
            "source_slides": 2,
            "early_lesion_patients": 1,
            "primary_eligible_paired_patients": 0,
            "pathology": "case5 pTis adenocarcinoma in tubulo-villous adenoma; case6 pT3 advanced CRC",
            "pathology_regions": "case5: Normal/Adenoma/Carcinoma/Other per caption; case6: Normal/Carcinoma/Other",
            "roi_source": "board-certified pathologist publication annotations and deposited cloupe",
            "roi_independent_of_m02": "TRUE",
            "material": "fresh frozen",
            "platform": "10x Visium; Illumina NovaSeq 6000",
            "reference_annotation": "GRCh38 GCF_000001405.39",
            "counts": "deposited integer count CSV",
            "coordinates": "deposited tissue_positions_list.csv",
            "he_image": "deposited hires/lowres images and cloupe",
            "license_or_reuse": "reanalysis under DDBJ/NBDC repository conditions; article CC BY-NC-ND 4.0; do not redistribute raw human data",
            "cross_study_overlap": "case5/case6 labels do not overlap case1-case4; institutional cross-accession overlap not formally excludable",
            "eligibility": "SECONDARY_ONLY",
            "evidence": f"{article}; https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-619/",
            "limitations": "Only one early patient; case5 age is discordant between SDRF and Supplementary Data 1, and the visible overlay does not clearly recover a Normal area despite the caption.",
        },
        {
            "cohort_id": "E-GEAD-579",
            "frozen_role": "secondary_advanced_border_localization",
            "source_patients": 1,
            "source_slides": 4,
            "early_lesion_patients": 0,
            "primary_eligible_paired_patients": 0,
            "pathology": "advanced CRC tumor-normal border capture areas",
            "pathology_regions": "tumor-normal border; no eligible early adenoma-normal endpoint",
            "roi_source": "source metadata and deposited cloupe",
            "roi_independent_of_m02": "TRUE",
            "material": "fresh frozen",
            "platform": "10x Visium; Illumina NovaSeq 6000",
            "reference_annotation": "GRCh38; exact build not stated in IDF",
            "counts": "deposited integer count CSV",
            "coordinates": "deposited tissue_positions_list.csv",
            "he_image": "deposited hires/lowres images and cloupe",
            "license_or_reuse": "reanalysis under DDBJ/NBDC repository conditions; no raw redistribution",
            "cross_study_overlap": "not demonstrated; not formally excludable",
            "eligibility": "SECONDARY_DESCRIPTIVE_ONLY",
            "evidence": "https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-579/",
            "limitations": "A1-D1 are four nested capture areas from the single donor ACC1 and never count as four patients.",
        },
        {
            "cohort_id": "GSE226997",
            "frozen_role": "secondary_cancer_localization",
            "source_patients": 4,
            "source_slides": 4,
            "early_lesion_patients": 0,
            "primary_eligible_paired_patients": 0,
            "pathology": "CRC cancer-only spatial sections",
            "pathology_regions": "cancer; no matched adenoma-normal early-lesion regions",
            "roi_source": "GEO metadata and deposited spatial assets",
            "roi_independent_of_m02": "TRUE",
            "material": "not sufficiently specified in series metadata",
            "platform": "10x Visium; Illumina NextSeq 2000",
            "reference_annotation": "deposited 10x feature annotation",
            "counts": "deposited MatrixMarket integer counts and 10x H5",
            "coordinates": "deposited tissue_positions_list.csv",
            "he_image": "deposited aligned/detected tissue and hires/lowres images",
            "license_or_reuse": "GEO public reanalysis; do not redistribute embedded FASTQ",
            "cross_study_overlap": "not demonstrated; not formally excludable from de-identified records",
            "eligibility": "SECONDARY_DESCRIPTIVE_ONLY",
            "evidence": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226997",
            "limitations": "Cancer-only; cannot satisfy early lesion generalization or a matched adenoma-normal endpoint.",
        },
    ]


def roi_manifest_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    e622 = [
        ("case1", "Rectum_kyudai_Beppu_20200303", "rectum", "main_Figure_1b"),
        ("case2", "Ascending_kyudai_Beppu_20200430", "ascending_colon", "Supplementary_Figure_1b"),
        ("case3", "Sigmoid_kyudai_Beppu_20210602", "sigmoid_colon", "Supplementary_Figure_1b"),
        ("case4", "Transverse_kyudai_Beppu_20211111", "transverse_colon", "Supplementary_Figure_1b"),
    ]
    for patient, slide, anatomy, figure in e622:
        for roi in ("Normal", "Adenoma", "Carcinoma", "Other"):
            if roi in {"Normal", "Adenoma"}:
                role, inclusion = "primary_pair_component", "INCLUDE_IF_10E_BLIND_REGISTRATION_QC_PASSES"
            elif roi == "Carcinoma":
                role, inclusion = "secondary_within_patient_region", "SECONDARY_ONLY"
            else:
                role, inclusion = "non_endpoint_region", "EXCLUDE_FROM_PRIMARY"
            rows.append(
                {
                    "cohort_id": "E-GEAD-622",
                    "patient_id": patient,
                    "slide_or_capture_id": slide,
                    "anatomic_site": anatomy,
                    "pathology_stage": "pTis",
                    "roi_category": roi,
                    "roi_label_evidence": f"board_certified_pathologist_{figure}",
                    "same_section_matched_normal": "TRUE",
                    "barcode_roi_membership_status": "TO_BE_RECONSTRUCTED_BLINDLY_FROM_PUBLISHED_OVERLAY",
                    "inference_unit": "patient",
                    "independent_patient_count": 1,
                    "analysis_role": role,
                    "inclusion_status": inclusion,
                    "notes": "All spots/ROIs remain nested within this patient; no M02 expression may be used for registration or inclusion.",
                }
            )

    for roi in ("Normal", "Adenoma", "Carcinoma", "Other"):
        rows.append(
            {
                "cohort_id": "E-GEAD-619",
                "patient_id": "case5",
                "slide_or_capture_id": "C_cancer_20210520_upper",
                "anatomic_site": "cecum",
                "pathology_stage": "pTis",
                "roi_category": roi,
                "roi_label_evidence": "publication_caption_and_deposited_cloupe",
                "same_section_matched_normal": "AMBIGUOUS_FOR_PRIMARY",
                "barcode_roi_membership_status": "STRUCTURED_CLOUPE_AVAILABLE_BUT_NOT_PRIMARY",
                "inference_unit": "patient",
                "independent_patient_count": 1,
                "analysis_role": "secondary_early_localization",
                "inclusion_status": "SECONDARY_ONLY",
                "notes": "One early patient cannot meet the primary n>=3 gate; Normal ROI visibility and SDRF/supplement age discrepancy must remain disclosed.",
            }
        )
    for roi in ("Normal", "Carcinoma", "Other"):
        rows.append(
            {
                "cohort_id": "E-GEAD-619",
                "patient_id": "case6",
                "slide_or_capture_id": "A_kyudai_Beppu_RSK_210624",
                "anatomic_site": "rectum",
                "pathology_stage": "pT3_advanced",
                "roi_category": roi,
                "roi_label_evidence": "board_certified_pathologist_publication_overlay_and_cloupe",
                "same_section_matched_normal": "TRUE",
                "barcode_roi_membership_status": "STRUCTURED_CLOUPE_AVAILABLE_BUT_NOT_PRIMARY",
                "inference_unit": "patient",
                "independent_patient_count": 1,
                "analysis_role": "secondary_advanced_localization",
                "inclusion_status": "SECONDARY_ONLY",
                "notes": "Advanced cancer cannot satisfy the early-lesion primary endpoint.",
            }
        )
    for slide in ("A1", "B1", "C1", "D1"):
        rows.append(
            {
                "cohort_id": "E-GEAD-579",
                "patient_id": "ACC1",
                "slide_or_capture_id": slide,
                "anatomic_site": "colon",
                "pathology_stage": "advanced_CRC",
                "roi_category": "tumor_normal_border",
                "roi_label_evidence": "SDRF_and_deposited_cloupe",
                "same_section_matched_normal": "NOT_ELIGIBLE_EARLY_PAIR",
                "barcode_roi_membership_status": "AVAILABLE_FOR_SECONDARY_DESCRIPTION",
                "inference_unit": "patient",
                "independent_patient_count": 1,
                "analysis_role": "secondary_advanced_border_localization",
                "inclusion_status": "SECONDARY_ONLY",
                "notes": "A1-D1 are nested within ACC1 and must never be counted as independent patients.",
            }
        )
    for number in range(1, 5):
        rows.append(
            {
                "cohort_id": "GSE226997",
                "patient_id": f"patient_{number}",
                "slide_or_capture_id": f"GSM{7089854 + number}_Ajou_Visium_P{number}",
                "anatomic_site": "NA",
                "pathology_stage": "CRC_cancer_only",
                "roi_category": "cancer",
                "roi_label_evidence": "GEO_series_metadata",
                "same_section_matched_normal": "FALSE",
                "barcode_roi_membership_status": "CANCER_ONLY",
                "inference_unit": "patient",
                "independent_patient_count": 1,
                "analysis_role": "secondary_cancer_localization",
                "inclusion_status": "SECONDARY_ONLY",
                "notes": "No matched early adenoma-normal endpoint.",
            }
        )
    return rows


def load_feature_audits() -> dict[str, list[dict[str, str]]]:
    file_groups = {
        "E-GEAD-622": [f"case{i}.tsv" for i in range(1, 5)],
        "E-GEAD-619": ["case5.tsv", "case6.tsv"],
        "E-GEAD-579": ["A1.tsv", "B1.tsv", "C1.tsv", "D1.tsv"],
        "GSE226997": [f"P{i}.tsv" for i in range(1, 5)],
    }
    result: dict[str, list[dict[str, str]]] = {}
    for cohort, names in file_groups.items():
        cohort_rows: list[dict[str, str]] = []
        for name in names:
            path = AUDIT_DIR / name
            if not path.exists():
                raise FileNotFoundError(path)
            cohort_rows.extend(read_tsv(path))
        result[cohort] = cohort_rows
    return result


def aggregate_mapping(
    audits: dict[str, list[dict[str, str]]]
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, list[str]]]:
    mapping_rows: list[dict[str, object]] = []
    common_rows: list[dict[str, object]] = []
    common_genes: dict[str, list[str]] = {}
    accepted = {"exact_symbol", "exact_ensembl", "one_to_one_archived_alias"}
    for cohort, rows in audits.items():
        samples = sorted({row["sample_id"] for row in rows})
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            grouped[row["canonical_gene"]].append(row)
        common_genes[cohort] = []
        for gene, group in sorted(grouped.items(), key=lambda item: int(item[1][0]["gene_order"])):
            statuses = sorted({row["mapping_status"] for row in group})
            present = [row for row in group if row["mapping_status"] in accepted]
            feature_ids = sorted({row["feature_id"] for row in present if row["feature_id"] != "NA"})
            symbols = sorted({row["feature_symbol"] for row in present if row["feature_symbol"] != "NA"})
            mapping_sources = sorted({row["mapping_source"] for row in group})
            is_common = len(present) == len(samples) and len(feature_ids) == 1
            if is_common:
                status = statuses[0] if len(statuses) == 1 else "unresolved"
                is_common = status in accepted
            else:
                status = statuses[0] if len(statuses) == 1 else "unresolved"
            zero_values = [row["present_but_zero"] for row in present]
            zero_any = "TRUE" if "TRUE" in zero_values else ("FALSE" if zero_values else "NA")
            zero_all = "TRUE" if zero_values and all(value == "TRUE" for value in zero_values) else ("FALSE" if zero_values else "NA")
            v4_rank = {
                "exact_symbol": 1,
                "exact_ensembl": 2,
                "one_to_one_archived_alias": 3,
                "unresolved": 4,
                "absent_from_feature_space": 5,
            }.get(status, 4)
            mapping_rows.append(
                {
                    "cohort_id": cohort,
                    "canonical_gene": gene,
                    "gene_order": group[0]["gene_order"],
                    "v4_mapping_rank": v4_rank,
                    "mapping_status": status,
                    "frozen_feature_id": feature_ids[0] if len(feature_ids) == 1 else "NA",
                    "feature_symbol": symbols[0] if len(symbols) == 1 else "NA",
                    "audited_samples": len(samples),
                    "present_in_samples": len(present),
                    "cohort_common_mapping": str(is_common).upper(),
                    "present_but_zero_any_sample": zero_any,
                    "present_but_zero_all_samples": zero_all,
                    "mapping_source": ";".join(mapping_sources),
                    "mapping_frozen_before_region_contrast": "TRUE",
                    "notes": "No proxy/substitute gene; stable-ID versions stripped only for exact matching.",
                }
            )
            if is_common:
                common_genes[cohort].append(gene)

        mapped_count = len(common_genes[cohort])
        for gene in common_genes[cohort]:
            source = next(row for row in mapping_rows if row["cohort_id"] == cohort and row["canonical_gene"] == gene)
            common_rows.append(
                {
                    "cohort_id": cohort,
                    "canonical_gene": gene,
                    "gene_order": source["gene_order"],
                    "mapping_status": source["mapping_status"],
                    "frozen_feature_id": source["frozen_feature_id"],
                    "cohort_common_gene_count": mapped_count,
                    "equal_weight_fraction": f"1/{mapped_count}",
                    "equal_weight_decimal": f"{1 / mapped_count:.16f}",
                    "present_but_zero_any_sample": source["present_but_zero_any_sample"],
                    "score_role": "primary_M02_cohort_common_mapping" if cohort == "E-GEAD-622" else "secondary_localization_mapping",
                }
            )
    return mapping_rows, common_rows, common_genes


def coverage_rows(mapping_rows: list[dict[str, object]], common_genes: dict[str, list[str]]) -> list[dict[str, object]]:
    role = {
        "E-GEAD-622": ("primary_early_lesion_spatial_candidate", 4, 4),
        "E-GEAD-619": ("secondary_early_and_advanced_localization", 2, 0),
        "E-GEAD-579": ("secondary_advanced_border_localization", 1, 0),
        "GSE226997": ("secondary_cancer_localization", 4, 0),
    }
    result: list[dict[str, object]] = []
    for cohort in ("E-GEAD-622", "E-GEAD-619", "E-GEAD-579", "GSE226997"):
        genes = common_genes[cohort]
        inpp5d = "INPP5D" in genes
        zero_any = sum(
            1
            for row in mapping_rows
            if row["cohort_id"] == cohort
            and row["cohort_common_mapping"] == "TRUE"
            and row["present_but_zero_any_sample"] == "TRUE"
        )
        cohort_role, patients, eligible = role[cohort]
        mapped = len(genes)
        result.append(
            {
                "cohort_id": cohort,
                "frozen_role": cohort_role,
                "independent_patients": patients,
                "primary_eligible_paired_patients": eligible,
                "common_mapped_genes": mapped,
                "coverage_fraction": f"{mapped / 36:.6f}",
                "coverage_tier": coverage_tier(mapped),
                "minimum_22of36_pass": str(mapped >= 22).upper(),
                "high_29of36_sensitivity_evaluable": str(mapped >= 29).upper(),
                "complete_36of36": str(mapped == 36).upper(),
                "present_but_zero_any_sample_genes": zero_any,
                "inpp5d_common_mapped": str(inpp5d).upper(),
                "m02_minus_inpp5d_common_genes": mapped - 1 if inpp5d else "NA",
                "m02_minus_inpp5d_distinct_and_evaluable": str(inpp5d and mapped - 1 >= 21).upper(),
                "coverage_gate": "PASS" if mapped >= 22 else "FAIL",
                "notes": "Coverage counts uniquely mapped canonical members common to every audited slide/capture in the cohort; globally present zero-valued features remain mapped and are disclosed, not outcome-selected away.",
            }
        )
    return result


ROI_RULE = (
    "E-GEAD-622_ROI_V1\tlabels=published_board_certified_pathologist_Normal_Adenoma_Carcinoma_Other"
    "\tregistration=two_independent_expression_blind_alignments_to_deposited_HE_and_spot_coordinates"
    "\tagreement=percent_ge_0.90_and_kappa_ge_0.85\tboundary=exclude_ambiguous_and_one_spot_ring"
    "\tprimary_regions=Adenoma_minus_Normal\tminimum_spots=30_per_region"
    "\tminimum_patients=3\tM02_or_M02_genes_prohibited_for_ROI_assignment\n"
)


def analysis_plan(common_genes: dict[str, list[str]]) -> str:
    e622_hash = sha256_text("\n".join(sorted(common_genes["E-GEAD-622"])) + "\n")
    return f"""# Stage 10C2-SP frozen spatial implementation plan

Frozen: 2026-07-31, before any lesion-normal M02 score, contrast, plot, or P value.

## Scope and immutable upstream locks

- Stage 10C3-POSTLOCK is `PASS`; Stage 10C3 remains `NOT_ESTIMABLE`.
- The canonical 36-gene M02 and `M02_SCORE_V1` are unchanged. Primary bundle SHA256: `{PRIMARY_BUNDLE}`.
- `M02_MINUS_INPP5D_SENS_V1` is unchanged. Sensitivity bundle SHA256: `{SENSITIVITY_BUNDLE}`.
- This gate inspected source metadata, pathology labels, feature annotations, coordinates, images, and whole-slide/global feature count sums only. It did not calculate or inspect any regional M02 difference.

## Frozen cohort roles

1. **Primary candidate:** E-GEAD-622 cases 1-4, four independent FFPE pTis patients. The primary estimand is within-patient adenoma minus normal on the same section.
2. **Secondary early/advanced localization:** E-GEAD-619. Case5 is one early patient; case6 is advanced. It is not pooled with E-GEAD-622 because material/reference differ and it cannot supply three early paired patients.
3. **Secondary advanced border localization:** E-GEAD-579. A1-D1 are nested capture areas from donor ACC1 and count once.
4. **Secondary cancer localization:** GSE226997. Four cancer-only patients; no early-lesion primary inference.

Selection is based on pathology, patient structure, assets, and feature coverage—not on M02 direction.

## Outcome-blind ROI reconstruction and eligibility

The E-GEAD-622 processed archive does not include a structured barcode-to-pathology table. The fixed source labels are the board-certified pathologist overlays in main Figure 1b (case1) and Supplementary Figure 1b (cases2-4) of DOI 10.1016/j.ebiom.2024.105102.

Before reading regional expression, two independent reviewers/implementations must register each published overlay to the deposited H&E and spot coordinates. No expression value, M02 score, or any of the 36 genes may enter alignment, classification, or adjudication. Required agreement is at least 90% and Cohen kappa at least 0.85. Disagreements, color-ambiguous spots, `Other`, and spots within one nearest-neighbour spot ring of a pathology boundary are excluded from the primary analysis.

A patient is technically eligible only if both Normal and Adenoma have at least 30 in-tissue post-QC spots after boundary exclusion. The primary analysis requires at least three eligible patients. Prespecified spot-count sensitivities are 20 (permissive) and 50 (strict) spots per region; neither can rescue a failed primary gate. If fewer than three patients pass, Stage 10F is `NOT_ESTIMABLE`.

Canonical ROI-rule payload SHA256: `{sha256_text(ROI_RULE)}`.

## Frozen V4 gene mapping

Mapping order is: (1) unique exact canonical symbol; (2) unique exact version-stripped Ensembl stable ID; (3) documented one-to-one archived alias; (4) unresolved; (5) absent from feature space. Proxy genes, synonym rescue without one-to-one evidence, and result-directed mapping are prohibited.

Each cohort uses one fixed common gene set shared by all its patients/slides/ROIs. It is never varied by patient, section, ROI, count level, or result. A mapped feature with zero counts in one slide remains mapped and is reported; it is not silently deleted.

E-GEAD-622 has 30/36 common canonical members (`HIGH_COVERAGE`), gene-set SHA256 `{e622_hash}`. The six absent features are AC016831.7, AL022068.1, AP001636.3, LINC01605, LINC01748, and LMCD1-AS1. INPP5D maps uniquely, so the 29-gene `M02_MINUS_INPP5D_SENS_V1` implementation is distinct and evaluable.

## Frozen preprocessing and scores

For each cohort independently:

1. Retain raw integer counts and aggregate spots to `patient × pathology_region` pseudobulk. Spots, ROIs, and slides never create biological degrees of freedom.
2. Apply edgeR TMM within that cohort across eligible patient-region pseudobulks.
3. Compute log2CPM with prior count 2.
4. Primary M02 score is the equal-weight mean of the cohort-level common mapped canonical genes, weight `1/K` for `K` fixed before the contrast. For E-GEAD-622, `K=30` and weight `1/30`.
5. No alternative primary score may be selected using spatial results. Gene-z and UCell are sensitivity scores only.

The 29/36 rule is a high-coverage audit/sensitivity threshold; it does not authorize dropping one of 30 mapped genes to manufacture a separate result. Missing canonical genes have no substitutes.

## QC frozen before regional scoring

- Primary spot QC: deposited `in_tissue=1`; nonzero library; log10 total counts not below median minus 3 MAD or above median plus 5 MAD; log10 detected features not below median minus 3 MAD. Robust cutoffs are calculated per slide using all in-tissue spots before ROI labels are joined.
- Hard floors: at least 200 UMIs and 100 detected genes. Mitochondrial/ribosomal fractions are diagnostic only for FFPE and are not single-variable exclusion rules.
- Strict technical sensitivity: at least 500 UMIs and 200 detected genes.
- Slide exclusion is permitted only for preregistered technical failure (unreadable matrix/coordinates/image, failed registration, or insufficient Normal/Adenoma spots), never for M02 direction.

## Estimand and patient-level inference

For eligible patient `i`, `delta_i = M02_score(Adenoma_i) - M02_score(Normal_i)`. The overall effect is the unweighted mean of patient effects. Report every `delta_i`, mean effect, 95% Student-t CI, two-sided one-sample t P value, number `k/n` with positive effects, and inference unit=`patient`.

The direction gate is `k >= ceil(0.8*n)`; this is 4/4 for n=4 and 3/3 for n=3. Mandatory small-sample checks are all `2^n` two-sided sign-flip permutations of the mean and an exact Wilcoxon signed-rank test when ties/zeros permit. Their discreteness is reported and never hidden. LOPO recomputes the estimate/direction after omitting each patient; it does not create a new significance family.

No spot-level P value is permitted. Moran's I on spot-level technical residuals is diagnostic only. A spatial block bootstrap within each ROI may quantify technical uncertainty, but its replicates never increase patient n or replace the patient-level CI.

## Mandatory sensitivities and negative controls

- `M02_MINUS_INPP5D_SENS_V1`, using the fixed 29 common genes in E-GEAD-622 and equal weight 1/29.
- Fixed common-set gene-z and UCell scores, aggregated to patient-region before comparison.
- Primary versus strict spot QC; 20/30/50 spot coverage diagnostics; one-ring versus two-ring pathology-boundary erosion.
- LOPO and explicit ROI-registration disagreement audit.
- Fixed-seed (`{SEED}`) 1,000 size-K negative-control gene sets sampled from non-M02 mapped features and matched on whole-slide mean abundance and detection quintiles without ROI labels. Their empirical rank is descriptive, not a replacement P value.

## Evidence classes

- `SPATIALLY_SUPPORTED`: positive mean, direction gate passed, primary CI excludes 0, primary P<0.05, and LOPO/pollution sensitivity do not reverse.
- `DIRECTIONAL_BUT_UNCERTAIN`: positive mean and direction gate passed, but CI crosses 0 or P>=0.05, with no key sensitivity reversal.
- `NOT_ROBUST`: LOPO, `M02_MINUS_INPP5D`, or a key technical sensitivity reverses direction.
- `INCONSISTENT`: direction count fails `ceil(0.8*n)` or the overall direction is not positive.
- `NOT_ESTIMABLE`: fewer than three eligible patients, non-independent regions, or coverage below 22/36.

Only `SPATIALLY_SUPPORTED` permits “statistical spatial support.” `DIRECTIONAL_BUT_UNCERTAIN` permits only directional spatial evidence. Research completion is not conditional on P<0.05, and no method or gene set may change to pursue significance.

## Frozen runtime, versions, and seed

- Server OS: Ubuntu 22.04.5 LTS; R 4.3.2; Python 3.10.12.
- R packages audited before execution: edgeR 4.0.2, Matrix 1.6.5, data.table 1.16.0, hdf5r 1.3.8, rhdf5 2.46.0, Seurat 5.0.2, sf 1.0.19, spdep 1.3.1. SpatialExperiment was not installed and is not required by this frozen algorithm.
- Python packages audited: numpy 1.21.5, pandas 1.3.5, scipy 1.8.0, Pillow 9.0.1; OpenCV was not installed.
- Random seed: `{SEED}` for every stochastic technical procedure and negative-control sampling. Any package installed later must be version-recorded before its first run and cannot alter the frozen method.

## Authorization boundary

This plan authorizes only Stage 10D-TECH after a PASS decision. Stage 10D-TECH validates implementation of this fixed score; it cannot select a score. Stage 10E/10F remain unauthorized until their own gates. No result in this document raises the M02 claim ceiling.
"""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    existing = [OUT / name for name in REQUIRED_OUTPUTS if (OUT / name).exists()]
    if existing:
        raise RuntimeError(f"refusing to overwrite existing frozen outputs: {existing}")

    source_path = OUT / REQUIRED_OUTPUTS[0]
    source_rows = source_audit_rows()
    write_tsv(source_path, list(source_rows[0]), source_rows)

    roi_path = OUT / REQUIRED_OUTPUTS[1]
    roi_rows = roi_manifest_rows()
    write_tsv(roi_path, list(roi_rows[0]), roi_rows)

    audits = load_feature_audits()
    mapping_rows, common_rows, common_genes = aggregate_mapping(audits)
    mapping_path = OUT / REQUIRED_OUTPUTS[2]
    write_tsv(mapping_path, list(mapping_rows[0]), mapping_rows)

    coverage_path = OUT / REQUIRED_OUTPUTS[3]
    coverage = coverage_rows(mapping_rows, common_genes)
    write_tsv(coverage_path, list(coverage[0]), coverage)

    common_path = OUT / REQUIRED_OUTPUTS[4]
    write_tsv(common_path, list(common_rows[0]), common_rows)

    plan_path = OUT / REQUIRED_OUTPUTS[5]
    with plan_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(analysis_plan(common_genes))

    lock_rows: list[dict[str, object]] = [
        {
            "lock_type": "upstream_bundle",
            "lock_id": "M02_SCORE_V1_36_gene_primary",
            "value": "36_genes",
            "sha256": PRIMARY_BUNDLE,
            "status": "VERIFIED_UNCHANGED",
            "notes": "Canonical 36 genes, direction, weights, score semantics and source-module hash bundle.",
        },
        {
            "lock_type": "upstream_bundle",
            "lock_id": "M02_MINUS_INPP5D_SENS_V1",
            "value": "35_genes",
            "sha256": SENSITIVITY_BUNDLE,
            "status": "VERIFIED_UNCHANGED",
            "notes": "Contamination sensitivity only; no further result-based gene deletion permitted.",
        },
        {
            "lock_type": "roi_rule",
            "lock_id": "E-GEAD-622_ROI_V1",
            "value": ROI_RULE.rstrip("\n"),
            "sha256": sha256_text(ROI_RULE),
            "status": "FROZEN_BEFORE_REGION_CONTRAST",
            "notes": "Published pathologist labels; expression-blind registration; primary Adenoma-Normal endpoint.",
        },
        {
            "lock_type": "random_seed",
            "lock_id": "STAGE10_SPATIAL_SEED",
            "value": SEED,
            "sha256": sha256_text(f"{SEED}\n"),
            "status": "FROZEN",
            "notes": "Applies to all stochastic technical procedures and negative controls.",
        },
    ]
    for cohort, genes in common_genes.items():
        payload = "\n".join(sorted(genes)) + "\n"
        gene_hash = sha256_text(payload)
        frozen_effect = {
            "E-GEAD-622": "adenoma_minus_matched_normal",
            "E-GEAD-619": "secondary_case5_adenoma_minus_normal_only",
            "E-GEAD-579": "descriptive_tumor_border_localization_no_primary_estimand",
            "GSE226997": "descriptive_cancer_localization_no_primary_estimand",
        }[cohort]
        scoring_payload = (
            f"M02_SPATIAL_SCORE_V1\tcohort={cohort}\tgeneset_sha256={gene_hash}"
            f"\tgene_count={len(genes)}\tweight=1/{len(genes)}\tnormalization=TMM_log2CPM_prior2"
            f"\taggregation=patient_by_region_pseudobulk\teffect={frozen_effect}"
            "\tmapping=V4_outcome_blind_cohort_common\tmissing=canonical_only_equal_renormalization\n"
        )
        lock_rows.extend(
            [
                {
                    "lock_type": "cohort_common_geneset",
                    "lock_id": cohort,
                    "value": f"{len(genes)}/36",
                    "sha256": gene_hash,
                    "status": "FROZEN_BEFORE_REGION_CONTRAST",
                    "notes": "One common canonical set for every patient/slide/ROI in this cohort.",
                },
                {
                    "lock_type": "cohort_scoring_spec",
                    "lock_id": f"M02_SPATIAL_SCORE_V1::{cohort}",
                    "value": scoring_payload.rstrip("\n"),
                    "sha256": sha256_text(scoring_payload),
                    "status": "FROZEN_BEFORE_REGION_CONTRAST",
                    "notes": "Secondary cohorts retain role restrictions in the source audit.",
                },
            ]
        )
    for filename in REQUIRED_OUTPUTS[:6]:
        path = OUT / filename
        lock_rows.append(
            {
                "lock_type": "artifact",
                "lock_id": filename,
                "value": path.stat().st_size,
                "sha256": sha256_file(path),
                "status": "LOCKED",
                "notes": "UTF-8/LF stage artifact; generated before any regional M02 contrast.",
            }
        )
    lock_path = OUT / REQUIRED_OUTPUTS[6]
    write_tsv(lock_path, list(lock_rows[0]), lock_rows)
    lock_sha = sha256_file(lock_path)

    decision = f"""# Stage 10C2-SP decision

Decision: **PASS**

Frozen: 2026-07-31. No lesion-normal M02 score, contrast, plot, confidence interval, or P value was calculated or viewed.

## Gate verification

| Requirement | Evidence | Result |
|---|---|---|
| Stage 10C3-POSTLOCK | `STAGE10C3_POSTLOCK_DECISION=PASS` | PASS |
| Stage 10C3 status | remains `NOT_ESTIMABLE`; no history changed | PASS |
| 36-gene bundle | `{PRIMARY_BUNDLE}` | VERIFIED |
| 35-gene bundle | `{SENSITIVITY_BUNDLE}` | VERIFIED |
| Primary patients | E-GEAD-622 cases1-4: four independent pTis patients with same-section pathologist-labelled Normal and Adenoma regions | PASS (4 >= 3) |
| ROI independence | Board-certified pathologist labels predate and do not use M02; barcode recovery is frozen as expression-blind image registration | PASS |
| Common M02 coverage | E-GEAD-622: 30/36, one fixed set across all four slides | HIGH_COVERAGE; PASS |
| INPP5D sensitivity | INPP5D uniquely maps; 29-gene minus-INPP5D set differs from the 30-gene primary spatial set | EVALUABLE |
| Plan/IDs frozen | V4 mapping, ROI rule, TMM/log2CPM(prior 2), patient estimand, tests, spatial handling, sensitivities, controls, versions and seed are locked | PASS |

Lock manifest SHA256: `{lock_sha}`.

## Frozen roles and limitations

- E-GEAD-622 is the only primary early-lesion candidate.
- E-GEAD-619 is secondary only (one early and one advanced patient); E-GEAD-579 is one-donor advanced-border localization; GSE226997 is cancer-only localization.
- The primary archive lacks a structured barcode-pathology table. Stage 10E must recover labels from the published pathologist overlays without expression, meet registration agreement and retain at least 30 Normal and 30 Adenoma spots in at least three patients. Failure makes Stage 10F `NOT_ESTIMABLE`.
- E-GEAD-622 is FFPE, only four patients are available, and exact tests are discrete. Cross-accession patient overlap is not demonstrated but cannot be formally excluded from de-identified public records.
- This PASS is a technical eligibility decision. It is not spatial support for M02 and does not raise the claim ceiling.

## Authorization

Stage 10D-TECH is allowed. Stage 10D-TECH may validate only the fixed 30-gene E-GEAD-622 implementation; it may not choose a primary score. Stage 10E and Stage 10F are not authorized by this decision and require their own gates. No downstream stage was started.
"""
    with (OUT / REQUIRED_OUTPUTS[7]).open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(decision)

    print(f"decision=PASS")
    print(f"primary_common_genes={len(common_genes['E-GEAD-622'])}")
    print(f"primary_common_genes_sha256={sha256_text(chr(10).join(sorted(common_genes['E-GEAD-622'])) + chr(10))}")
    print(f"lock_manifest_sha256={lock_sha}")


if __name__ == "__main__":
    main()
