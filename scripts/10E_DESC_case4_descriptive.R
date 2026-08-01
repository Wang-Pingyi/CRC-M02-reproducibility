#!/usr/bin/env Rscript

# Analysis: Stage 10E-DESC case4-only descriptive spatial localization
# Date: 2026-08-01
# Random seed: 20260731 (frozen Stage 10 spatial seed)
# Inference: none; one patient, descriptive only

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: 10E_DESC_case4_descriptive.R <root> <run_id> <smoke|formal>")
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
run_id <- args[[2L]]
mode <- args[[3L]]
if (!mode %in% c("smoke", "formal")) stop("mode must be smoke or formal")

project_lib <- file.path(root, "environment", "R-library")
if (dir.exists(project_lib)) .libPaths(unique(c(project_lib, .libPaths())))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(edgeR)
  library(UCell)
  library(ggplot2)
  library(digest)
})

options(stringsAsFactors = FALSE)
params <- fread(file.path(root, "config", "stage10e_desc_parameters.tsv"))
get_param <- function(x) {
  z <- params[parameter == x, value]
  if (length(z) != 1L) stop("Missing/duplicate parameter: ", x)
  z
}
set.seed(as.integer(get_param("seed")))

expected_hash <- get_param("common_geneset_hash")
expected_patient <- get_param("patient_id")
expected_slide <- get_param("slide_or_capture_id")
if (expected_patient != "case4") stop("Only case4 is authorized")

required <- c(
  "results/stage10e_roi_remediation/STAGE10E_R_ACCEPTANCE.md",
  "results/stage10e_roi_remediation/STAGE10E_R_DECISION.md",
  "results/stage10e_roi_remediation/STAGE10E_R_PATIENT_ELIGIBILITY.tsv",
  "results/stage10e_roi_remediation/STAGE10E_R_ROI_AGREEMENT.tsv",
  "results/stage10e_roi_remediation/STAGE10E_R_GENESET_INTEGRITY_CHECK.tsv",
  "results/stage10c2_sp/STAGE10C2_SP_COHORT_COMMON_GENESET.tsv",
  "results/stage10c2_sp/STAGE10C2_SP_LOCK_MANIFEST.tsv",
  "results/stage10d_tech/STAGE10D_TECH_DECISION.md",
  "results/stage10e/STAGE10E_ANALYSIS_READY_MANIFEST.tsv",
  "results/stage10e_desc/STAGE10E_DESC_ANALYSIS_PLAN_LOCKED.md"
)
missing_required <- required[!file.exists(file.path(root, required))]
if (length(missing_required)) stop("Missing required input: ", paste(missing_required, collapse = ", "))

acceptance_text <- paste(readLines(file.path(root, required[[1L]]), warn = FALSE), collapse = "\n")
decision_text <- paste(readLines(file.path(root, required[[2L]]), warn = FALSE), collapse = "\n")
tech_text <- paste(readLines(file.path(root, "results/stage10d_tech/STAGE10D_TECH_DECISION.md"), warn = FALSE), collapse = "\n")
if (!grepl("PASS_DESCRIPTIVE_ONLY", acceptance_text, fixed = TRUE) ||
    !grepl("PASS_DESCRIPTIVE_ONLY", decision_text, fixed = TRUE)) stop("Stage 10E-R gate mismatch")
if (!grepl("TECHNICALLY_VALID", tech_text, fixed = TRUE)) stop("Stage 10D-TECH is not technically valid")

eligibility <- fread(file.path(root, "results/stage10e_roi_remediation/STAGE10E_R_PATIENT_ELIGIBILITY.tsv"))
if (nrow(eligibility[patient_decision == "PASS"]) != 1L ||
    eligibility[patient_decision == "PASS", patient_id] != "case4") stop("case4 is not the sole frozen eligible patient")
roi_agreement_gate <- fread(file.path(root, "results/stage10e_roi_remediation/STAGE10E_R_ROI_AGREEMENT.tsv"))
if (nrow(roi_agreement_gate[patient_id == "case4" & registration_pass == TRUE]) != 1L) stop("case4 registration gate mismatch")
remediation_integrity <- fread(file.path(root, "results/stage10e_roi_remediation/STAGE10E_R_GENESET_INTEGRITY_CHECK.tsv"))
if (!all(remediation_integrity$status == "PASS") || !all(remediation_integrity$recomputed_sha256 == expected_hash)) {
  stop("Stage 10E-R geneset integrity mismatch")
}
sp_lock <- fread(file.path(root, "results/stage10c2_sp/STAGE10C2_SP_LOCK_MANIFEST.tsv"))
if (nrow(sp_lock[lock_type == "cohort_common_geneset" & lock_id == "E-GEAD-622" & sha256 == expected_hash]) != 1L) {
  stop("Stage 10C2-SP cohort lock mismatch")
}
historical_ready <- fread(file.path(root, "results/stage10e/STAGE10E_ANALYSIS_READY_MANIFEST.tsv"))
if (!all(c("case1", "case2", "case3", "case4") %in% historical_ready$patient_id)) stop("Historical Stage 10E manifest incomplete")

common_all <- fread(file.path(root, "results/stage10c2_sp/STAGE10C2_SP_COHORT_COMMON_GENESET.tsv"))
common30_table <- common_all[cohort_id == "E-GEAD-622"]
common30 <- sort(unique(common30_table$canonical_gene))
recomputed_hash <- digest(paste0(paste(common30, collapse = "\n"), "\n"), algo = "sha256", serialize = FALSE)
if (length(common30) != 30L || recomputed_hash != expected_hash) stop("Frozen 30/36 mapping hash mismatch")

slide_dir <- file.path(root, "data_processed", "stage10e", "stage10e_20260801_010406", "derived", expected_slide)
h5_path <- file.path(slide_dir, "filtered_feature_bc_matrix.h5")
positions_path <- file.path(slide_dir, "tissue_positions_list.csv")
roi_path <- file.path(root, "figures", "stage10e_roi_remediation", "source_data", "Fig10E_R_case4_registration_qc_source_data.tsv")
qc_path <- file.path(root, "results", "stage10e", "STAGE10E_SPATIAL_QC.tsv")
if (!all(file.exists(c(h5_path, positions_path, roi_path, qc_path)))) stop("case4 input missing")

read_positions <- function(path) {
  z <- fread(path, header = FALSE)
  if (ncol(z) != 6L) stop("Unexpected positions schema")
  setnames(z, c("spot_barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  z
}

counts <- Seurat::Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE)
if (is.list(counts)) counts <- counts[[1L]]
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
if (!all(common30 %in% rownames(counts))) stop("Not all frozen 30 genes are present in case4")

roi <- fread(roi_path)
if (!identical(unique(roi$patient_id), "case4")) stop("ROI source contains a non-case4 patient")
roi <- roi[, .(patient_id, spot_barcode, consensus_label, pipeline_A_x, pipeline_A_y, pipeline_B_x, pipeline_B_y)]
if (anyDuplicated(roi$spot_barcode)) stop("Duplicate barcode in frozen ROI")

qc <- fread(qc_path, select = c("patient_id", "slide_or_capture_id", "spot_barcode", "final_qc_pass",
                                "epithelial_reference_fraction_proxy"))
qc <- qc[patient_id == "case4" & slide_or_capture_id == expected_slide]
if (!nrow(qc) || anyDuplicated(qc$spot_barcode)) stop("case4 QC input invalid")
positions <- read_positions(positions_path)

meta <- merge(roi, qc, by = c("patient_id", "spot_barcode"), all.x = TRUE, sort = FALSE)
meta <- merge(meta, positions, by = "spot_barcode", all.x = TRUE, sort = FALSE)
meta[, final_qc_pass := as.logical(final_qc_pass)]
eligible <- meta[final_qc_pass & consensus_label %in% c("Normal", "Adenoma") & spot_barcode %in% colnames(counts)]
if (any(eligible$patient_id != "case4")) stop("Non-case4 patient entered scoring")
spot_counts <- eligible[, .N, by = consensus_label]
if (!all(c("Normal", "Adenoma") %in% spot_counts$consensus_label) ||
    any(spot_counts$N < as.integer(get_param("minimum_spots_per_region")))) stop("case4 frozen ROI coverage is insufficient")

if (mode == "smoke") {
  smoke_dir <- file.path(root, "cache", "stage10e_desc_smoke", run_id)
  dir.create(smoke_dir, recursive = TRUE, showWarnings = FALSE)
  nonmodule <- setdiff(rownames(counts), common30)[seq_len(20L)]
  smoke_spots <- eligible$spot_barcode[seq_len(min(40L, nrow(eligible)))]
  smoke_counts <- counts[nonmodule, smoke_spots, drop = FALSE]
  lib <- Matrix::colSums(counts[, smoke_spots, drop = FALSE])
  smoke_score <- colMeans(edgeR::cpm(smoke_counts, lib.size = lib, log = TRUE, prior.count = 2))
  smoke <- data.table(
    test = c("required_inputs", "sole_patient_gate", "common_hash", "case4_barcode_join", "nonmodule_score_schema"),
    status = "PASS",
    detail = c(length(required), "case4_only", recomputed_hash,
               paste0(nrow(eligible), " eligible spots"), paste0(length(smoke_score), " non-M02 smoke values")),
    M02_biological_result_computed = FALSE
  )
  fwrite(smoke, file.path(smoke_dir, "STAGE10E_DESC_SMOKE_TEST.tsv"), sep = "\t")
  writeLines(capture.output(sessionInfo()), file.path(smoke_dir, "sessionInfo.txt"))
  quit(save = "no", status = 0L)
}

out_dir <- file.path(root, "results", "stage10e_desc")
fig_dir <- file.path(root, "figures", "stage10e_desc")
source_dir <- file.path(fig_dir, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# Primary patient x region raw-count aggregation.
regions <- c("Normal", "Adenoma")
pb <- sapply(regions, function(region) {
  barcodes <- eligible[consensus_label == region, spot_barcode]
  Matrix::rowSums(counts[, barcodes, drop = FALSE])
})
rownames(pb) <- rownames(counts); colnames(pb) <- regions
dge <- edgeR::DGEList(counts = pb)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = as.numeric(get_param("prior_count")))
primary_score <- colMeans(logcpm[common30, regions, drop = FALSE])
primary_delta <- unname(primary_score[["Adenoma"]] - primary_score[["Normal"]])

mapped_out <- copy(common30_table)
mapped_out[, `:=`(patient_id = "case4", score_weight = 1/30,
                  verified_common_hash = recomputed_hash, expression_feature_present = canonical_gene %in% rownames(counts))]
fwrite(mapped_out, file.path(out_dir, "STAGE10E_DESC_MAPPED_GENESET.tsv"), sep = "\t")

scores <- data.table(
  patient_id = "case4", patient_n = 1L, region = regions,
  score_method = "primary_30gene_equal_weight_TMM_log2CPM_prior2",
  score = as.numeric(primary_score[regions]),
  spot_count = vapply(regions, function(x) eligible[consensus_label == x, .N], integer(1)),
  mapped_genes = 30L, canonical_genes = 36L, normalization = "TMM_log2CPM_prior2",
  aggregation = "patient_by_frozen_ROI_pseudobulk", biological_inference = "none_descriptive_only",
  p_value = "NOT_COMPUTED_BY_DESIGN", confidence_interval = "NOT_COMPUTED_BY_DESIGN", fdr = "NOT_COMPUTED_BY_DESIGN"
)
fwrite(scores, file.path(out_dir, "STAGE10E_DESC_CASE4_SCORES.tsv"), sep = "\t")

# Sensitivity 1: fixed minus-INPP5D set.
minus_genes <- setdiff(common30, "INPP5D")
minus_score <- colMeans(logcpm[minus_genes, regions, drop = FALSE])

# Sensitivity 2: gene-wise z on the two region pseudobulks.
gene_z <- t(scale(t(logcpm[common30, regions, drop = FALSE])))
gene_z[!is.finite(gene_z)] <- 0
gene_z_score <- colMeans(gene_z)

# Sensitivity 3: UCell on the same two region pseudobulk vectors.
ucell <- UCell::ScoreSignatures_UCell(as(pb, "dgCMatrix"), features = list(M02 = common30),
                                      maxRank = min(as.integer(get_param("ucell_max_rank")), nrow(pb)), ncores = 1L)
ucell_col <- grep("M02.*UCell|M02", colnames(ucell), value = TRUE)[1L]
if (!length(ucell_col) || is.na(ucell_col)) stop("UCell score column missing")
ucell_score <- setNames(as.numeric(ucell[regions, ucell_col]), regions)

# Sensitivity 4/5 and map source: unsmoothed spot score with composition residualization.
eligible_barcodes <- eligible$spot_barcode
lib_sizes <- Matrix::colSums(counts[, eligible_barcodes, drop = FALSE])
spot_logcpm <- edgeR::cpm(counts[common30, eligible_barcodes, drop = FALSE], lib.size = lib_sizes,
                          log = TRUE, prior.count = as.numeric(get_param("prior_count")))
spot_score <- colMeans(spot_logcpm)
eligible[, spot_m02_score := as.numeric(spot_score[spot_barcode])]
raw_spot_region <- eligible[, .(score = mean(spot_m02_score)), by = consensus_label]
raw_spot_score <- setNames(raw_spot_region$score, raw_spot_region$consensus_label)

epi_ok <- is.finite(eligible$epithelial_reference_fraction_proxy) & is.finite(eligible$spot_m02_score)
if (sum(epi_ok) >= 10L && sd(eligible$epithelial_reference_fraction_proxy[epi_ok]) > 0) {
  composition_fit <- lm(spot_m02_score ~ epithelial_reference_fraction_proxy, data = eligible[epi_ok])
  eligible[, epithelial_residual := NA_real_]
  eligible[epi_ok, epithelial_residual := residuals(composition_fit)]
  residual_region <- eligible[is.finite(epithelial_residual), .(score = mean(epithelial_residual)), by = consensus_label]
  residual_score <- setNames(residual_region$score, residual_region$consensus_label)
} else {
  residual_score <- setNames(c(NA_real_, NA_real_), regions)
}

# Sensitivity 6: fixed 2x2 tiles. Tiles remain technical nested observations.
eligible[, tile_row := ifelse(array_row <= median(array_row), 1L, 2L)]
eligible[, tile_col := ifelse(array_col <= median(array_col), 1L, 2L)]
eligible[, tile_id := paste0("R", tile_row, "C", tile_col)]
tile_groups <- unique(eligible[, .(consensus_label, tile_id)])
tile_pb <- sapply(seq_len(nrow(tile_groups)), function(i) {
  g <- tile_groups[i]
  bc <- eligible[consensus_label == g$consensus_label & tile_id == g$tile_id, spot_barcode]
  Matrix::rowSums(counts[, bc, drop = FALSE])
})
if (is.null(dim(tile_pb))) tile_pb <- matrix(tile_pb, ncol = 1L)
rownames(tile_pb) <- rownames(counts)
colnames(tile_pb) <- paste(tile_groups$consensus_label, tile_groups$tile_id, sep = "__")
tile_dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = tile_pb), method = "TMM")
tile_logcpm <- edgeR::cpm(tile_dge, log = TRUE, prior.count = as.numeric(get_param("prior_count")))
tile_groups[, tile_score := colMeans(tile_logcpm[common30, , drop = FALSE])]
tile_region <- tile_groups[, .(score = mean(tile_score)), by = consensus_label]
tile_score <- setNames(tile_region$score, tile_region$consensus_label)

make_sensitivity <- function(name, role, normal, adenoma, genes, aggregation, note) {
  delta <- adenoma - normal
  primary_direction <- if (primary_delta > 0) "positive" else if (primary_delta < 0) "negative" else "zero"
  direction <- if (!is.finite(delta)) "not_estimable" else if (delta > as.numeric(get_param("decision_tolerance"))) "positive" else
    if (delta < -as.numeric(get_param("decision_tolerance"))) "negative" else "zero"
  data.table(patient_id = "case4", patient_n = 1L, sensitivity = name, role = role,
             Normal_score = normal, Adenoma_score = adenoma, Adenoma_minus_Normal = delta,
             direction = direction, direction_vs_primary = ifelse(direction == "not_estimable", "not_estimable",
                ifelse(direction == primary_direction, "concordant", "opposite_or_zero")),
             mapped_genes = genes, aggregation = aggregation, note = note,
             p_value = "NOT_COMPUTED_BY_DESIGN", confidence_interval = "NOT_COMPUTED_BY_DESIGN", fdr = "NOT_COMPUTED_BY_DESIGN")
}

sensitivity <- rbindlist(list(
  make_sensitivity("primary_reference", "primary_descriptive", primary_score[["Normal"]], primary_score[["Adenoma"]], 30L,
                   "patient_by_frozen_ROI_pseudobulk", "Reference direction only"),
  make_sensitivity("M02_MINUS_INPP5D", "contamination_sensitivity", minus_score[["Normal"]], minus_score[["Adenoma"]], 29L,
                   "patient_by_frozen_ROI_pseudobulk", "INPP5D removed only in prespecified sensitivity"),
  make_sensitivity("gene_z", "score_sensitivity", gene_z_score[["Normal"]], gene_z_score[["Adenoma"]], 30L,
                   "two_region_pseudobulk_gene_z", "Two-region z scale is descriptive and coarse"),
  make_sensitivity("UCell", "score_sensitivity", ucell_score[["Normal"]], ucell_score[["Adenoma"]], 30L,
                   "two_region_pseudobulk_UCell", "Rank score sensitivity"),
  make_sensitivity("epithelial_proxy_unadjusted", "composition_sensitivity", raw_spot_score[["Normal"]], raw_spot_score[["Adenoma"]], 30L,
                   "spot_score_then_patient_region_mean", "Spots are nested technical observations"),
  make_sensitivity("epithelial_proxy_residualized", "composition_sensitivity", residual_score[["Normal"]], residual_score[["Adenoma"]], 30L,
                   "spot_score_residualized_on_M02_excluded_epithelial_proxy_then_region_mean", "No region term and no inference"),
  make_sensitivity("fixed_2x2_tiles", "aggregation_scale_sensitivity", tile_score[["Normal"]], tile_score[["Adenoma"]], 30L,
                   "fixed_2x2_tiles_then_equal_tile_mean", "Tiles do not increase patient n")
), fill = TRUE)
fwrite(sensitivity, file.path(out_dir, "STAGE10E_DESC_SENSITIVITY.tsv"), sep = "\t", na = "NA")

primary_direction <- if (primary_delta > as.numeric(get_param("decision_tolerance"))) "positive" else
  if (primary_delta < -as.numeric(get_param("decision_tolerance"))) "negative" else "zero"
mandatory <- sensitivity[sensitivity != "primary_reference"]
opposite <- mandatory[direction %in% c("positive", "negative", "zero") & direction != primary_direction, .N]
not_estimable <- mandatory[direction == "not_estimable", .N]
decision <- if (!is.finite(primary_delta) || not_estimable > 0L) "NOT_ESTIMABLE" else if (primary_direction == "zero" || opposite > 0L) {
  "DESCRIPTIVE_METHOD_DEPENDENT"
} else if (primary_direction == "positive") "DESCRIPTIVE_POSITIVE" else "DESCRIPTIVE_NEGATIVE"

# Source data and fixed-scale figures.
map_source <- eligible[, .(patient_id = "case4", patient_n = 1L, spot_barcode, consensus_label,
                           pxl_col, pxl_row, array_col, array_row, final_qc_pass,
                           epithelial_reference_fraction_proxy, spot_m02_score,
                           mapped_genes = 30L, canonical_genes = 36L, interpretation = "descriptive only")]
fwrite(map_source, file.path(source_dir, "Fig10E_DESC_1_case4_spatial_source_data.tsv"), sep = "\t")
fwrite(scores[, .(patient_id, patient_n, region, score, spot_count, mapped_genes, canonical_genes,
                  interpretation = "descriptive only")],
       file.path(source_dir, "Fig10E_DESC_2_case4_paired_source_data.tsv"), sep = "\t")

p1 <- ggplot(map_source, aes(x = pxl_col, y = -pxl_row, colour = spot_m02_score)) +
  geom_point(size = 0.55) +
  scale_colour_viridis_c(limits = c(as.numeric(get_param("spatial_colour_min")), as.numeric(get_param("spatial_colour_max"))),
                         oob = scales::squish, name = "30-gene score") +
  coord_equal() + theme_void(base_size = 9) +
  labs(title = "case4 frozen ROI: M02 spatial localization",
       subtitle = "n=1 patient; 30/36 mapped genes; descriptive only",
       caption = "Unsmoothed spot display; spots are not biological replicates")

p2 <- ggplot(scores, aes(x = factor(region, levels = regions), y = score, group = patient_id)) +
  geom_line(linewidth = 0.5, colour = "#666666") + geom_point(size = 2.6, colour = "#0072B2") +
  theme_classic(base_size = 9) + xlab(NULL) + ylab("Equal-weight TMM log2CPM score") +
  labs(title = "case4 Normal and Adenoma ROI scores",
       subtitle = "n=1 patient; 30/36 mapped genes; descriptive only",
       caption = "No P value, confidence interval or population inference")

for (item in list(list(p = p1, name = "Fig10E_DESC_1_case4_spatial"),
                  list(p = p2, name = "Fig10E_DESC_2_case4_paired"))) {
  ggsave(file.path(fig_dir, paste0(item$name, ".png")), item$p, width = 7, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, paste0(item$name, ".pdf")), item$p, width = 7, height = 5, device = cairo_pdf)
}

claim_limits <- c(
  "# Stage 10E-DESC claim limits", "",
  "This stage contains one patient and is descriptive only.", "",
  "Allowed: case4 Normal and Adenoma regional scores, the within-case4 descriptive difference, unsmoothed frozen-ROI localization, and directions of prespecified sensitivities.", "",
  "Prohibited: spatial validation; spatially supported; confirmed; replicated; significant; robust across patients; generalizable; sporadic adenoma validation; mechanistic evidence; biomarker evidence.", "",
  "Spots, tiles and ROIs are nested technical observations and do not increase n beyond one patient. Stage 10F and Stage 10G remain skipped regardless of direction."
)
writeLines(claim_limits, file.path(out_dir, "STAGE10E_DESC_CLAIM_LIMITS.md"))

decision_lines <- c(
  "# Stage 10E-DESC decision", "", paste0("Decision: **", decision, "**"), "",
  "- Patient: case4 only; n=1 patient.",
  paste0("- Frozen mapped set: 30/36; SHA256 `", recomputed_hash, "`."),
  paste0("- Primary descriptive Adenoma minus Normal difference: ", format(primary_delta, digits = 8), "."),
  paste0("- Mandatory sensitivity reversals: ", opposite, "."),
  "- P values, confidence intervals and FDR: not computed by design.",
  "- Stage 10F: remains skipped.", "- Stage 10G: remains skipped.",
  "- Additional spatial analysis: not authorized.", "- Next governance step: Stage 10FG-CLOSE, not started."
)
writeLines(decision_lines, file.path(out_dir, "STAGE10E_DESC_DECISION.md"))

summary_lines <- c(
  "# Stage 10E-DESC summary", "", paste0("Run ID: `", run_id, "`"), "",
  "- Scope: case4-only descriptive localization; n=1 patient.",
  paste0("- Eligible frozen-ROI spots: Normal=", scores[region == "Normal", spot_count],
         "; Adenoma=", scores[region == "Adenoma", spot_count], "."),
  paste0("- Mapped genes: 30/36; hash `", recomputed_hash, "`."),
  paste0("- Normal primary score: ", format(primary_score[["Normal"]], digits = 8), "."),
  paste0("- Adenoma primary score: ", format(primary_score[["Adenoma"]], digits = 8), "."),
  paste0("- Adenoma minus Normal: ", format(primary_delta, digits = 8), "."),
  paste0("- Decision: `", decision, "`."),
  "- Statistical inference: none; no P value, CI or FDR.",
  "- Stage 10F/10G: skipped; no later stage started."
)
writeLines(summary_lines, file.path(root, "reports", "STAGE10E_DESC_SUMMARY.md"))
writeLines(c("# Stage 10E-DESC gate decision", "", paste0("Decision: **", decision, "**"), "",
             "This is a single-patient descriptive decision only. It does not authorize Stage 10F or Stage 10G.",
             "Next stage is not authorized by this file."),
           file.path(root, "reports", "STAGE10E_DESC_GATE_DECISION.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "STAGE10E_DESC_SESSIONINFO.txt"))

manifest <- data.table(
  run_id = run_id, patient_id = "case4", slide_or_capture_id = expected_slide,
  h5_sha256 = digest(file = h5_path, algo = "sha256"), roi_sha256 = digest(file = roi_path, algo = "sha256"),
  qc_sha256 = digest(file = qc_path, algo = "sha256"), common_geneset_sha256 = recomputed_hash,
  plan_sha256 = digest(file = file.path(out_dir, "STAGE10E_DESC_ANALYSIS_PLAN_LOCKED.md"), algo = "sha256"),
  seed = get_param("seed"), biological_patient_n = 1L, stage10f_started = FALSE, stage10g_started = FALSE
)
fwrite(manifest, file.path(out_dir, "STAGE10E_DESC_RUN_MANIFEST.tsv"), sep = "\t")

message("Stage 10E-DESC completed: ", decision)
