#!/usr/bin/env Rscript

# Stage 10E: spatial preprocessing and blinded technical QC.
# This script deliberately never computes an M02 score, lesion-normal contrast,
# M02 heatmap, effect estimate, confidence interval, or P value.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args)) stop("Missing value for ", flag)
  args[hit[length(hit)] + 1L]
}

root <- normalizePath(arg_value("--root", "."), mustWork = TRUE)
run_id <- arg_value("--run-id")
input_dir <- normalizePath(arg_value("--input-dir"), mustWork = TRUE)
smoke <- "--smoke" %in% args
if (is.null(run_id)) stop("--run-id is required")

project_lib <- file.path(root, "environment", "R-library")
if (dir.exists(project_lib)) .libPaths(unique(c(.libPaths(), project_lib)))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(digest)
  library(jpeg)
  library(png)
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260731L)

expected <- list(
  primary_bundle = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30",
  sensitivity_bundle = "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8",
  common_geneset = "4cd62d74b83673a4d2adf6077bedbdfe73d1cbd369a6a77418d124e0b506d482"
)

hash_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)
hash_gene_vector <- function(x) digest(paste0(paste(sort(unique(x)), collapse = "\n"), "\n"),
                                       algo = "sha256", serialize = FALSE)
must_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = "; "))
}

param_path <- file.path(root, "config", "stage10e_parameters.tsv")
panel_path <- file.path(root, "config", "stage10e_roi_overlay_panels.tsv")
decision_sp_path <- file.path(root, "results", "stage10c2_sp", "STAGE10C2_SP_DECISION.md")
decision_tech_path <- file.path(root, "results", "stage10d_tech", "STAGE10D_TECH_DECISION.md")
lock36_path <- file.path(root, "results", "stage10c", "STAGE10C_LOCK_MANIFEST.tsv")
lock35_path <- file.path(root, "results", "stage10c", "M02_MINUS_INPP5D_SENS_V1.tsv")
lock35_sha_path <- file.path(root, "results", "stage10c", "M02_MINUS_INPP5D_SENS_V1_SHA256.txt")
common_path <- file.path(root, "results", "stage10c2_sp", "STAGE10C2_SP_COHORT_COMMON_GENESET.tsv")
patient_path <- file.path(root, "results", "stage10c2_sp", "STAGE10C2_SP_PATIENT_SLIDE_ROI_MANIFEST.tsv")
reference_path <- file.path(root, "objects", "GSE161277_stage7_annotated.rds")
must_exist(c(param_path, panel_path, decision_sp_path, decision_tech_path, lock36_path, lock35_path, lock35_sha_path,
             common_path, patient_path, reference_path))

params <- fread(param_path)
get_param <- function(name) {
  x <- params[parameter == name, value]
  if (length(x) != 1L) stop("Expected one parameter: ", name)
  x
}

sp_decision <- paste(readLines(decision_sp_path, warn = FALSE), collapse = "\n")
tech_decision <- paste(readLines(decision_tech_path, warn = FALSE), collapse = "\n")
if (!grepl("Decision: \\*\\*PASS\\*\\*", sp_decision)) stop("Stage 10C2-SP gate is not PASS")
if (!grepl("Decision: \\*\\*(TECHNICALLY_VALID|VALID_WITH_LIMITATIONS)\\*\\*", tech_decision)) {
  stop("Stage 10D-TECH is not technically valid")
}

lock36 <- fread(lock36_path)
lock35 <- fread(lock35_path)
if (nrow(lock36) != 36L || length(unique(lock36$gene)) != 36L ||
    unique(lock36$bundle_sha256) != expected$primary_bundle) stop("36-gene lock mismatch")
if (nrow(lock35) != 35L || "INPP5D" %in% lock35$gene ||
    !grepl(expected$sensitivity_bundle, paste(readLines(lock35_sha_path, warn = FALSE), collapse = " "), fixed = TRUE)) stop("35-gene lock mismatch")
canonical36 <- lock36$gene[order(lock36$gene_order)]
common_all <- fread(common_path)
common30 <- common_all[cohort_id == "E-GEAD-622", canonical_gene]
if (length(common30) != 30L || hash_gene_vector(common30) != expected$common_geneset) {
  stop("Frozen E-GEAD-622 common gene set mismatch")
}

panels <- fread(panel_path, na.strings = c("NA", ""))
slides <- c("Rectum_kyudai_Beppu_20200303", "Ascending_kyudai_Beppu_20200430",
            "Sigmoid_kyudai_Beppu_20210602", "Transverse_kyudai_Beppu_20211111")
if (!setequal(panels$slide_or_capture_id, slides) || uniqueN(panels$patient_id) != 4L) {
  stop("Frozen four-patient overlay panel configuration mismatch")
}
if (smoke) slides <- slides[1L]

out_root <- if (smoke) file.path(root, "cache", "stage10e_smoke", run_id) else file.path(root, "results", "stage10e")
derived_root <- file.path(root, "data_processed", "stage10e", run_id, "derived")
fig_root <- if (smoke) file.path(out_root, "figures") else file.path(root, "figures", "stage10e")
fig_source <- file.path(fig_root, "source_data")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_source, recursive = TRUE, showWarnings = FALSE)

read_positions <- function(path) {
  x <- fread(path, header = FALSE)
  if (ncol(x) != 6L) stop("Expected six-column tissue_positions_list.csv: ", path)
  setnames(x, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  x[, `:=`(barcode = as.character(barcode), in_tissue = as.integer(in_tissue),
           array_row = as.numeric(array_row), array_col = as.numeric(array_col),
           pxl_row = as.numeric(pxl_row), pxl_col = as.numeric(pxl_col))]
  x
}

robust_bounds <- function(x, floor_value, upper = FALSE, multiplier = 3) {
  med <- median(x, na.rm = TRUE); md <- mad(x, center = med, constant = 1, na.rm = TRUE)
  if (!is.finite(md) || md == 0) md <- max(1, abs(med) * 0.05)
  if (upper) min(floor_value, med + multiplier * md) else max(floor_value, med - multiplier * md)
}

marker_sets <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "ELF3", "CEACAM5", "CEACAM6", "GUCA2A", "MUC2"),
  immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "CD79A", "MS4A1", "LYZ", "FCER1G", "TYROBP", "LST1"),
  fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "PDGFRA"),
  endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "ENG", "ESAM", "CLDN5")
)
if (length(intersect(unique(unlist(marker_sets)), canonical36))) stop("M02 genes in reference marker panel")

message("Building independent broad-cell reference profiles with M02-excluded markers")
reference <- readRDS(reference_path)
ref_counts <- tryCatch(GetAssayData(reference, assay = "RNA", layer = "counts"),
                       error = function(e) GetAssayData(reference, assay = "RNA", slot = "counts"))
ref_meta <- reference[[]]
if (!all(c("major_cell_type", "donor_id", "annotation_status") %in% names(ref_meta))) stop("Reference metadata incomplete")
reference_group <- rep(NA_character_, nrow(ref_meta))
reference_group[ref_meta$major_cell_type == "Epithelial"] <- "epithelial"
reference_group[ref_meta$major_cell_type %in% c("T_NK", "B_cell", "Plasma_cell", "Myeloid")] <- "immune"
reference_group[ref_meta$major_cell_type == "Fibroblast"] <- "fibroblast"
reference_group[ref_meta$major_cell_type == "Endothelial"] <- "endothelial"
ref_keep <- !is.na(reference_group) & ref_meta$annotation_status == "canonical_marker_cluster_annotation"
if (sum(ref_keep) < as.integer(get_param("reference_min_cells")) ||
    length(unique(ref_meta$donor_id[ref_keep])) < as.integer(get_param("reference_min_donors"))) {
  stop("Independent reference coverage below frozen requirement")
}
ref_markers <- intersect(unique(unlist(marker_sets)), rownames(ref_counts))
if (length(ref_markers) < 12L) stop("Reference marker coverage inadequate")
ref_profile <- sapply(names(marker_sets), function(group) {
  cells <- which(ref_keep & reference_group == group)
  if (!length(cells)) stop("No independent reference cells for ", group)
  as.numeric(Matrix::rowSums(ref_counts[ref_markers, cells, drop = FALSE])) + 0.5
})
rownames(ref_profile) <- ref_markers
ref_profile <- sweep(ref_profile, 2L, colSums(ref_profile), "/")
rm(reference, ref_counts, ref_meta); gc(verbose = FALSE)

estimate_reference_proportions <- function(counts, genes) {
  shared <- intersect(rownames(ref_profile), genes)
  if (length(shared) < 12L) stop("Spatial marker coverage inadequate for reference proportions")
  a <- ref_profile[shared, , drop = FALSE]
  a <- sweep(a, 2L, colSums(a), "/")
  y <- counts[shared, , drop = FALSE]
  y <- sweep(y, 2L, pmax(colSums(y), 1), "/")
  ridge <- diag(1e-8, ncol(a))
  ans <- vapply(seq_len(ncol(y)), function(i) {
    b <- tryCatch(solve(crossprod(a) + ridge, crossprod(a, as.numeric(y[, i]))),
                  error = function(e) rep(NA_real_, ncol(a)))
    b <- pmax(as.numeric(b), 0)
    if (!all(is.finite(b)) || sum(b) == 0) return(rep(NA_real_, ncol(a)))
    b / sum(b)
  }, numeric(ncol(a)))
  rownames(ans) <- colnames(a); colnames(ans) <- colnames(y)
  ans
}

read_overlay <- function(path) {
  if (grepl("\\.jpg$|\\.jpeg$", path, ignore.case = TRUE)) readJPEG(path) else readPNG(path)
}
prototypes <- rbind(Adenoma = c(0.1216, 0.4667, 0.7059), Carcinoma = c(1, 0.4980, 0.0549),
                    Normal = c(0.1725, 0.6275, 0.1725), Other = c(0.8392, 0.1529, 0.1569))
classify_colour <- function(image, x, y, radius = 2L) {
  x <- as.integer(round(x)); y <- as.integer(round(y))
  x0 <- max(1L, x - radius); x1 <- min(dim(image)[2L], x + radius)
  y0 <- max(1L, y - radius); y1 <- min(dim(image)[1L], y + radius)
  pix <- matrix(image[y0:y1, x0:x1, 1:3, drop = FALSE], ncol = 3)
  distance <- vapply(seq_len(nrow(prototypes)), function(i) {
    sqrt(rowSums((pix - matrix(prototypes[i, ], nrow(pix), 3, byrow = TRUE))^2))
  }, numeric(nrow(pix)))
  if (is.null(dim(distance))) distance <- matrix(distance, ncol = 1L)
  nearest_by_prototype <- apply(distance, 2L, min)
  best <- which.min(nearest_by_prototype); d <- nearest_by_prototype[best]
  if (!is.finite(d) || d > 0.34) return(NA_character_)
  rownames(prototypes)[best]
}
map_coordinates <- function(pos, panel, method) {
  if (method == "fullres_pixel_bbox_affine") {
    x <- pos$pxl_col; y <- pos$pxl_row
  } else if (method == "array_lattice_bbox_affine") {
    x <- pos$array_col; y <- pos$array_row
  } else stop("Unknown frozen ROI registration method")
  scale01 <- function(z) (z - min(z, na.rm = TRUE)) / (max(z, na.rm = TRUE) - min(z, na.rm = TRUE))
  data.table(x = panel$crop_x + scale01(x) * (panel$crop_width - 1),
             y = panel$crop_y + scale01(y) * (panel$crop_height - 1))
}
cohen_kappa <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 2L) return(NA_real_)
  x <- x[ok]; y <- y[ok]; lev <- union(unique(x), unique(y))
  po <- mean(x == y)
  pe <- sum(vapply(lev, function(z) mean(x == z) * mean(y == z), numeric(1)))
  if (abs(1 - pe) < 1e-12) return(NA_real_)
  (po - pe) / (1 - pe)
}
boundary_ring <- function(pos, label) {
  n <- nrow(pos); out <- rep(FALSE, n)
  eligible <- !is.na(label)
  idx <- which(eligible)
  # Technical diagnostic only; all final endpoint exclusions remain frozen before M02 access.
  for (i in idx) {
    d2 <- (pos$pxl_row - pos$pxl_row[i])^2 + (pos$pxl_col - pos$pxl_col[i])^2
    nearest <- order(d2)[seq_len(min(7L, n))]
    nearest <- setdiff(nearest, i)
    out[i] <- any(!is.na(label[nearest]) & label[nearest] != label[i])
  }
  out
}

qc_tables <- list(); deconv_tables <- list(); map_tables <- list(); gene_tables <- list(); hash_tables <- list()
for (slide in slides) {
  message("Blind QC: ", slide)
  h5_path <- file.path(input_dir, slide, "filtered_feature_bc_matrix.h5")
  pos_path <- file.path(input_dir, slide, "tissue_positions_list.csv")
  must_exist(c(h5_path, pos_path))
  counts <- Seurat::Read10X_h5(h5_path, use.names = TRUE, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[[1L]]
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  pos <- read_positions(pos_path)
  pos <- pos[barcode %in% colnames(counts)]
  if (!nrow(pos) || nrow(pos) != ncol(counts)) stop("Barcode mismatch for ", slide)
  setkey(pos, barcode); counts <- counts[, pos$barcode, drop = FALSE]
  patient <- panels[slide_or_capture_id == slide]
  if (nrow(patient) != 1L) stop("No unique frozen panel for ", slide)

  present_common <- common30 %in% rownames(counts)
  gene_tables[[slide]] <- data.table(cohort_id = "E-GEAD-622", patient_id = patient$patient_id,
    slide_or_capture_id = slide, canonical_gene = common30, present_in_feature_space = present_common,
    frozen_common_geneset_sha256 = expected$common_geneset,
    expression_read_or_used_for_outcome = FALSE,
    integrity_status = ifelse(present_common, "PRESENT_FROZEN_SET_RETAINED", "ANNOTATION_ERROR_STOP_REQUIRED"))
  if (!all(present_common)) stop("Frozen common gene set feature integrity failure for ", slide)

  ncount <- Matrix::colSums(counts); nfeature <- Matrix::colSums(counts > 0)
  mt <- grep("^MT-", rownames(counts), value = TRUE)
  mt_percent <- if (length(mt)) Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(ncount, 1) * 100 else rep(0, length(ncount))
  tissue <- pos$in_tissue == 1L
  low_count <- robust_bounds(ncount[tissue], as.numeric(get_param("minimum_counts")))
  low_feature <- robust_bounds(nfeature[tissue], as.numeric(get_param("minimum_features")))
  high_mt <- robust_bounds(mt_percent[tissue], as.numeric(get_param("maximum_mt_percent")), upper = TRUE,
                              multiplier = as.numeric(get_param("mt_mad_multiplier")))
  high_count_flag <- ncount > (median(ncount[tissue]) + as.numeric(get_param("high_complexity_mad_multiplier")) *
                               max(mad(ncount[tissue], constant = 1), 1))
  q <- as.numeric(get_param("edge_quantile"))
  xr <- rank(pos$pxl_col, ties.method = "average") / nrow(pos); yr <- rank(pos$pxl_row, ties.method = "average") / nrow(pos)
  edge_flag <- xr <= q | xr >= 1 - q | yr <= q | yr >= 1 - q
  final_qc <- tissue & ncount >= low_count & nfeature >= low_feature & mt_percent <= high_mt
  exclusion_reason <- ifelse(!tissue, "outside_tissue", ifelse(ncount < low_count, "low_counts",
                     ifelse(nfeature < low_feature, "low_features", ifelse(mt_percent > high_mt, "high_mt", "retained"))))

  overlay_path <- if (startsWith(patient$overlay_source, "derived/")) {
    file.path(input_dir, sub("^derived/", "", patient$overlay_source))
  } else {
    file.path(root, patient$overlay_source)
  }
  must_exist(overlay_path)
  overlay <- read_overlay(overlay_path)
  xy_a <- map_coordinates(pos, patient, patient$coordinate_method_A)
  xy_b <- map_coordinates(pos, patient, patient$coordinate_method_B)
  label_a <- vapply(seq_len(nrow(pos)), function(i) classify_colour(overlay, xy_a$x[i], xy_a$y[i]), character(1))
  label_b <- vapply(seq_len(nrow(pos)), function(i) classify_colour(overlay, xy_b$x[i], xy_b$y[i]), character(1))
  agreement <- mean(label_a[!is.na(label_a) & !is.na(label_b)] == label_b[!is.na(label_a) & !is.na(label_b)])
  kappa <- cohen_kappa(label_a, label_b)
  consensus <- ifelse(!is.na(label_a) & label_a == label_b, label_a, NA_character_)
  boundary <- boundary_ring(pos, consensus)
  map_tables[[slide]] <- data.table(patient_id = patient$patient_id, slide_or_capture_id = slide,
    spot_barcode = pos$barcode, roi_method_A = label_a, roi_method_B = label_b, roi_consensus = consensus,
    roi_agreement = agreement, roi_kappa = kappa, boundary_ring_flag = boundary,
    roi_assignment_expression_blind = TRUE, pathology_source = patient$pathology_source)

  proportions <- estimate_reference_proportions(counts, rownames(counts))
  deconv_tables[[slide]] <- data.table(patient_id = patient$patient_id, slide_or_capture_id = slide,
    spot_barcode = colnames(counts), epithelial_reference_fraction_proxy = proportions["epithelial", ],
    immune_reference_fraction_proxy = proportions["immune", ], fibroblast_reference_fraction_proxy = proportions["fibroblast", ],
    endothelial_reference_fraction_proxy = proportions["endothelial", ], reference = "GSE161277_independent_broad_profiles",
    marker_panel_excludes_all_36_M02_genes = TRUE)
  qc_tables[[slide]] <- data.table(patient_id = patient$patient_id, slide_or_capture_id = slide,
    spot_barcode = pos$barcode, in_tissue = tissue, nCount = as.numeric(ncount), nFeature = as.numeric(nfeature),
    mitochondrial_percent = as.numeric(mt_percent), tissue_coverage = mean(tissue), edge_flag = edge_flag,
    high_complexity_diagnostic_flag = high_count_flag, low_count_threshold = low_count,
    low_feature_threshold = low_feature, high_mt_threshold = high_mt, final_qc_pass = final_qc,
    exclusion_reason = exclusion_reason, run_id = run_id)
  hash_tables[[slide]] <- data.table(artifact_type = "raw_or_source_input", patient_id = patient$patient_id,
    slide_or_capture_id = slide, path = c(h5_path, pos_path, overlay_path),
    sha256 = vapply(c(h5_path, pos_path, overlay_path), hash_file, character(1)),
    bytes = file.info(c(h5_path, pos_path, overlay_path))$size)
  rm(counts, overlay); gc(verbose = FALSE)
}

qc <- rbindlist(qc_tables); deconv <- rbindlist(deconv_tables); roi <- rbindlist(map_tables); geneset <- rbindlist(gene_tables); hashes <- rbindlist(hash_tables)
setkey(qc, patient_id, slide_or_capture_id, spot_barcode); setkey(roi, patient_id, slide_or_capture_id, spot_barcode)
qc <- roi[qc]
setkey(deconv, patient_id, slide_or_capture_id, spot_barcode); qc <- deconv[qc]
qc[, primary_endpoint_spot_eligible := final_qc_pass & !is.na(roi_consensus) & !boundary_ring_flag & roi_consensus %in% c("Normal", "Adenoma")]

roi_summary <- qc[, .(spots_total = .N, spots_in_tissue = sum(in_tissue), spots_final_qc = sum(final_qc_pass),
  spots_primary_endpoint_eligible = sum(primary_endpoint_spot_eligible),
  median_nCount = median(nCount), median_nFeature = median(nFeature), median_mt_percent = median(mitochondrial_percent),
  roi_agreement = unique(roi_agreement), roi_kappa = unique(roi_kappa),
  registration_pass = unique(roi_agreement) >= as.numeric(get_param("roi_agreement_min")) &
                      unique(roi_kappa) >= as.numeric(get_param("roi_kappa_min"))),
  by = .(patient_id, slide_or_capture_id, roi_category = roi_consensus)]
roi_summary <- roi_summary[!is.na(roi_category)]
patient_pair <- roi_summary[roi_category %in% c("Normal", "Adenoma"),
  .(technical_spot_coverage = sum(spots_primary_endpoint_eligible), registration_pass = all(registration_pass)),
  by = .(patient_id, slide_or_capture_id, roi_category)]
patient_pair[, region_minimum_met := technical_spot_coverage >= as.integer(get_param("minimum_spots_per_region"))]
pair_summary <- dcast(patient_pair, patient_id + slide_or_capture_id + registration_pass ~ roi_category,
                      value.var = "region_minimum_met", fill = FALSE)
if (!"Normal" %in% names(pair_summary)) pair_summary[, Normal := FALSE]
if (!"Adenoma" %in% names(pair_summary)) pair_summary[, Adenoma := FALSE]
pair_summary[, paired_technical_eligible := registration_pass & Normal & Adenoma]

all_qc_readable <- nrow(qc) > 0 && all(geneset$present_in_feature_space)
n_pairs <- sum(pair_summary$paired_technical_eligible)
decision <- if (!all_qc_readable) "FAIL" else if (n_pairs >= as.integer(get_param("minimum_patients"))) "PASS" else "PASS_WITH_LIMITATIONS"

if (!smoke) {
  out_meta <- file.path(root, "data", "metadata", "spatial_patient_slide_roi_manifest.tsv")
  dir.create(dirname(out_meta), recursive = TRUE, showWarnings = FALSE)
  fwrite(roi_summary, out_meta, sep = "\t", na = "NA")
  fwrite(qc, file.path(out_root, "STAGE10E_SPATIAL_QC.tsv"), sep = "\t", na = "NA")
  fwrite(qc[!final_qc_pass | is.na(roi_consensus) | boundary_ring_flag,
    .(patient_id, slide_or_capture_id, spot_barcode, final_qc_pass, roi_method_A, roi_method_B, roi_consensus,
      boundary_ring_flag, exclusion_reason, exclusion_stage = "pre_primary_endpoint_expression_blind")],
    file.path(out_root, "STAGE10E_EXCLUSION_LOG.tsv"), sep = "\t", na = "NA")
  fwrite(deconv, file.path(out_root, "STAGE10E_DECONVOLUTION_QC.tsv"), sep = "\t", na = "NA")
  fwrite(geneset, file.path(out_root, "STAGE10E_GENESET_INTEGRITY_CHECK.tsv"), sep = "\t", na = "NA")
  fwrite(pair_summary, file.path(out_root, "STAGE10E_ANALYSIS_READY_MANIFEST.tsv"), sep = "\t", na = "NA")
  fwrite(hashes, file.path(out_root, "STAGE10E_RAW_HASHES.tsv"), sep = "\t", na = "NA")
  fwrite(qc, file.path(fig_source, "Fig10E_1_spot_qc_source_data.tsv"), sep = "\t", na = "NA")
  fwrite(deconv, file.path(fig_source, "Fig10E_2_reference_deconvolution_source_data.tsv"), sep = "\t", na = "NA")
  fwrite(pair_summary, file.path(fig_source, "Fig10E_3_analysis_ready_source_data.tsv"), sep = "\t", na = "NA")
  theme_qc <- theme_bw(base_size = 10) + theme(panel.grid = element_blank())
  p1 <- ggplot(qc, aes(slide_or_capture_id, nCount, colour = final_qc_pass)) + geom_jitter(width = .16, size = .45) +
    scale_y_log10() + theme_qc + theme(axis.text.x = element_text(angle = 35, hjust = 1)) + labs(x = NULL, y = "Spot counts (log10)", colour = "Final QC")
  p2 <- ggplot(deconv, aes(slide_or_capture_id, epithelial_reference_fraction_proxy)) + geom_boxplot(outlier.size = .3) +
    theme_qc + theme(axis.text.x = element_text(angle = 35, hjust = 1)) + labs(x = NULL, y = "Reference-derived epithelial fraction proxy")
  p3 <- ggplot(pair_summary, aes(patient_id, as.numeric(paired_technical_eligible), fill = registration_pass)) + geom_col() +
    scale_y_continuous(breaks = c(0, 1), labels = c("not eligible", "eligible")) + theme_qc + labs(x = "Patient", y = "Paired technical eligibility", fill = "Registration pass")
  for (x in list(list(p = p1, n = "Fig10E_1_spot_qc"), list(p = p2, n = "Fig10E_2_deconvolution"), list(p = p3, n = "Fig10E_3_readiness"))) {
    ggsave(file.path(fig_root, paste0(x$n, ".png")), x$p, width = 7.2, height = 4.8, dpi = 300)
    ggsave(file.path(fig_root, paste0(x$n, ".pdf")), x$p, width = 7.2, height = 4.8, device = cairo_pdf)
  }
  decision_lines <- c("# Stage 10E decision", "", paste0("Decision: **", decision, "**"), "",
    "## Expression-blind QC gate", "",
    paste0("- Primary cohort: E-GEAD-622; slides processed: ", length(slides), "; patients: ", uniqueN(qc$patient_id), "."),
    paste0("- Frozen common M02 mapping: 30/36; hash ", expected$common_geneset, ". Feature identity was checked only; no M02 expression score or region contrast was computed."),
    paste0("- Broad-cell deconvolution reference: GSE161277; all reference marker panels excluded every locked M02 gene."),
    paste0("- Patients with paired Normal/Adenoma technical coverage after blind QC and frozen boundary exclusion: ", n_pairs, "."),
    "- Spot, ROI and slide observations remain nested within patient; no spatial spot was treated as an independent biological replicate.",
    "- This stage did not calculate, plot, inspect or infer a lesion-normal M02 difference.", "",
    "## Decision interpretation", "",
    if (decision == "PASS") "All Stage 10E technical and blinded-registration conditions were met; Stage 10F still requires separate authorization." else if (decision == "PASS_WITH_LIMITATIONS")
      "Raw data, QC and frozen gene identities are technically usable, but fewer than three patients met the fully paired, blinded ROI-registration technical gate. Stage 10F primary inference is not authorized from this result." else
      "Input readability or frozen geneset integrity failed. Stop and return to governance review; do not repair using M02 results.")
  writeLines(decision_lines, file.path(out_root, "STAGE10E_DECISION.md"), useBytes = TRUE)
  writeLines(c("# Stage 10E server summary", "", paste0("- Run ID: ", run_id), paste0("- Decision: ", decision),
    paste0("- Slides/patients: ", length(slides), "/", uniqueN(qc$patient_id)), paste0("- Eligible paired patients: ", n_pairs),
    "- M02 lesion-normal outcome accessed: NO", "- Next stage started: NO"),
    file.path(root, "reports", "STAGE10E_SUMMARY.md"), useBytes = TRUE)
  writeLines(decision_lines, file.path(root, "reports", "STAGE10E_GATE_DECISION.md"), useBytes = TRUE)
  session <- capture.output(sessionInfo())
  writeLines(session, file.path(out_root, "STAGE10E_SESSIONINFO.txt"), useBytes = TRUE)
  manifest_paths <- c(param_path, panel_path, decision_sp_path, decision_tech_path, lock36_path, lock35_path, lock35_sha_path,
                      common_path, patient_path, reference_path, file.path(root, "scripts", "10E_spatial_blind_qc.R"),
                      file.path(out_root, c("STAGE10E_SPATIAL_QC.tsv", "STAGE10E_EXCLUSION_LOG.tsv", "STAGE10E_DECONVOLUTION_QC.tsv",
                        "STAGE10E_GENESET_INTEGRITY_CHECK.tsv", "STAGE10E_ANALYSIS_READY_MANIFEST.tsv", "STAGE10E_RAW_HASHES.tsv", "STAGE10E_DECISION.md", "STAGE10E_SESSIONINFO.txt")),
                      file.path(root, "data", "metadata", "spatial_patient_slide_roi_manifest.tsv"),
                      list.files(fig_root, recursive = TRUE, full.names = TRUE), file.path(root, "reports", c("STAGE10E_SUMMARY.md", "STAGE10E_GATE_DECISION.md")))
  manifest_paths <- unique(manifest_paths[file.exists(manifest_paths)])
  run_manifest <- data.table(path = normalizePath(manifest_paths), sha256 = vapply(manifest_paths, hash_file, character(1)), bytes = file.info(manifest_paths)$size)
  root_norm <- normalizePath(root)
  run_manifest[, path := sub("^[\\\\/]", "", substring(path, nchar(root_norm) + 1L))]
  fwrite(run_manifest, file.path(out_root, "STAGE10E_RUN_MANIFEST.tsv"), sep = "\t", na = "NA")
}
cat(if (smoke) "STAGE10E_SMOKE_OK\n" else paste0("STAGE10E_COMPLETE decision=", decision, "\n"))
