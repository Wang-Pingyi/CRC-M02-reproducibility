#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)
get_param <- read_stage7_parameters(project_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(ggplot2)
})

min_count <- get_param("qc", "min_nCount", TRUE)
min_feature <- get_param("qc", "min_nFeature", TRUE)
max_mt <- get_param("qc", "max_percent_mt", TRUE)
upper_k <- get_param("qc", "upper_mad_k", TRUE)

raw_dir <- file.path(project_dir, "data_raw", "GSE132465")
count_file <- file.path(
  raw_dir, "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
)
annotation_file <- file.path(
  raw_dir, "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
)

message("Reading official GSE132465 cell annotation")
annotation <- data.table::fread(annotation_file, data.table = FALSE)
required_columns <- c("Index", "Patient", "Class", "Sample", "Cell_type", "Cell_subtype")
if (!all(required_columns %in% colnames(annotation))) {
  stop("GSE132465 annotation schema differs from the locked preflight")
}
if (anyDuplicated(annotation$Index)) stop("GSE132465 annotation cell IDs are duplicated")

message("Reading official dense GSE132465 UMI matrix")
dt <- data.table::fread(
  count_file, data.table = TRUE, check.names = FALSE, showProgress = TRUE
)
gene_column <- colnames(dt)[[1L]]
genes <- as.character(dt[[gene_column]])
if (anyDuplicated(genes)) {
  warning("Duplicate GSE132465 gene symbols made unique without changing source rows")
  genes <- make.unique(genes)
}
matrix_cells <- colnames(dt)[-1L]
if (!setequal(matrix_cells, annotation$Index)) {
  stop("GSE132465 matrix and annotation cell IDs differ")
}
annotation <- annotation[match(matrix_cells, annotation$Index), , drop = FALSE]
if (!identical(matrix_cells, annotation$Index)) stop("GSE132465 cell order mismatch")

counts <- as.matrix(dt[, -1L, with = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- genes
rm(dt)
invisible(gc())

n_count <- colSums(counts)
n_feature <- colSums(counts > 0)
mt_idx <- grepl("^MT-", rownames(counts))
ribo_idx <- grepl("^RP[SL]", rownames(counts))
percent_mt <- if (any(mt_idx)) {
  100 * colSums(counts[mt_idx, , drop = FALSE]) / pmax(n_count, 1)
} else rep(0, ncol(counts))
percent_ribo <- if (any(ribo_idx)) {
  100 * colSums(counts[ribo_idx, , drop = FALSE]) / pmax(n_count, 1)
} else rep(0, ncol(counts))

high_complexity <- rep(FALSE, ncol(counts))
adaptive_max_feature <- rep(NA_real_, ncol(counts))
adaptive_max_count <- rep(NA_real_, ncol(counts))
for (sample_id in unique(annotation$Sample)) {
  idx <- annotation$Sample == sample_id
  fb <- robust_qc_bounds(log10(n_feature[idx] + 1), 3, upper_k)
  cb <- robust_qc_bounds(log10(n_count[idx] + 1), 3, upper_k)
  max_f <- ceiling(10^fb[["upper"]] - 1)
  max_c <- ceiling(10^cb[["upper"]] - 1)
  adaptive_max_feature[idx] <- max_f
  adaptive_max_count[idx] <- max_c
  high_complexity[idx] <- n_feature[idx] > max_f | n_count[idx] > max_c
}
low_quality <- n_count < min_count | n_feature < min_feature | percent_mt > max_mt
retained <- !(low_quality | high_complexity)
reason <- ifelse(
  low_quality, "residual_low_quality",
  ifelse(high_complexity, "residual_high_complexity", "retained")
)

major_map <- c(
  "Epithelial cells" = "Epithelial",
  "T cells" = "T_NK",
  "B cells" = "B_cell",
  "Myeloids" = "Myeloid",
  "Stromal cells" = "Stromal",
  "Mast cells" = "Mast_cell"
)
major_cell_type <- unname(major_map[annotation$Cell_type])
major_cell_type[is.na(major_cell_type)] <- "Other_source_annotation"
condition <- ifelse(annotation$Class == "Tumor", "cancer", "normal")

cell_metadata <- data.frame(
  cell_id = annotation$Index,
  donor_id = annotation$Patient,
  condition = condition,
  sample_id = annotation$Sample,
  source_cell_type = annotation$Cell_type,
  source_cell_subtype = annotation$Cell_subtype,
  major_cell_type = major_cell_type,
  nCount = n_count,
  nFeature = n_feature,
  percent_mt = percent_mt,
  percent_ribo = percent_ribo,
  adaptive_max_nCount = adaptive_max_count,
  adaptive_max_nFeature = adaptive_max_feature,
  exclusion_reason = reason,
  source_doublet_status = "official_processed_matrix_high_gene_outliers_removed",
  stringsAsFactors = FALSE
)
write_stage7_tsv(
  cell_metadata,
  file.path(paths$source, "GSE132465_qc_annotation_cell_source.tsv.gz")
)

qc_summary <- aggregate(
  cbind(
    cells_input = rep(1L, nrow(cell_metadata)),
    cells_retained = as.integer(retained),
    residual_low_quality = as.integer(low_quality),
    residual_high_complexity = as.integer(high_complexity & !low_quality)
  ) ~ sample_id + donor_id + condition,
  data = cell_metadata, FUN = sum
)
qc_summary$doublet_action <- "not_rerun_official_processed_matrix"
qc_summary$ribosomal_action <- "diagnostic_only"
write_stage7_tsv(
  qc_summary, file.path(paths$result, "GSE132465_qc_summary.tsv")
)

locked <- read_locked_stage7_modules(project_dir)
candidate_genes <- unique(locked$membership$gene)
signatures <- major_marker_signatures(candidate_genes)
marker_genes <- unique(unlist(signatures))
marker_genes <- intersect(marker_genes, rownames(counts))
marker_log <- log1p(
  t(t(counts[marker_genes, retained, drop = FALSE]) /
      pmax(n_count[retained], 1)) * 10000
)
marker_scores <- score_signatures(marker_log, signatures)$scores
audit_meta <- cell_metadata[retained, c("cell_id", "major_cell_type"), drop = FALSE]
marker_audit <- do.call(rbind, lapply(unique(audit_meta$major_cell_type), function(type) {
  idx <- audit_meta$major_cell_type == type
  data.frame(
    source_major_cell_type = type,
    signature = colnames(marker_scores),
    mean_signature_score = colMeans(marker_scores[idx, , drop = FALSE], na.rm = TRUE),
    n_cells = sum(idx),
    stringsAsFactors = FALSE
  )
}))
marker_audit$highest_signature_for_source_type <- ave(
  marker_audit$mean_signature_score, marker_audit$source_major_cell_type,
  FUN = function(x) x == max(x, na.rm = TRUE)
)
write_stage7_tsv(
  marker_audit,
  file.path(paths$result, "GSE132465_annotation_evidence.tsv")
)
rm(marker_log, marker_scores)
invisible(gc())

kept_meta <- cell_metadata[retained, , drop = FALSE]
kept_counts <- counts[, retained, drop = FALSE]
group <- paste(
  kept_meta$donor_id, kept_meta$condition, kept_meta$major_cell_type, sep = "|"
)
message("Aggregating GSE132465 raw counts to donor-condition-cell-type")
pseudobulk_counts <- aggregate_counts(kept_counts, group)
group_meta <- unique(
  kept_meta[, c("donor_id", "condition", "major_cell_type"), drop = FALSE]
)
group_meta$pseudobulk_id <- paste(
  group_meta$donor_id, group_meta$condition, group_meta$major_cell_type, sep = "|"
)
group_meta <- group_meta[match(colnames(pseudobulk_counts), group_meta$pseudobulk_id), ]
group_meta$n_cells <- as.integer(table(group)[group_meta$pseudobulk_id])
group_meta$total_umi <- as.numeric(colSums(pseudobulk_counts))
group_meta$cohort <- "GSE132465"
group_meta$matched_primary_donor <- group_meta$donor_id %in% sprintf("SMC%02d", 1:10)

saveRDS(
  list(counts = pseudobulk_counts, metadata = group_meta),
  file.path(paths$processed, "GSE132465_stage7_pseudobulk_raw_counts.rds"),
  compress = "gzip"
)
saveRDS(
  list(
    retained_cell_metadata = kept_meta,
    pseudobulk_metadata = group_meta,
    source_annotation_policy =
      "official annotation retained; canonical-marker concordance independently audited",
    doublet_policy =
      "official processed matrix already high-gene-outlier filtered; not rerun",
    candidate_reselection = FALSE,
    seed = stage7_seed
  ),
  file.path(paths$object, "GSE132465_stage7_processed_summary.rds"),
  compress = "gzip"
)

composition <- aggregate(
  cell_id ~ donor_id + condition + major_cell_type,
  data = kept_meta, FUN = length
)
names(composition)[names(composition) == "cell_id"] <- "n_cells"
write_stage7_tsv(
  composition, file.path(paths$result, "GSE132465_cell_composition.tsv")
)

p_qc <- ggplot(qc_summary, aes(condition, cells_retained, fill = condition)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 1) +
  labs(title = "GSE132465 cells retained after residual QC", y = "Cells per sample") +
  theme_bw(base_size = 9)
p_comp <- ggplot(
  composition, aes(condition, n_cells, fill = major_cell_type)
) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ major_cell_type, scales = "free_y") +
  labs(title = "GSE132465 independent source-annotation composition", y = "Cells") +
  theme_bw(base_size = 8) + theme(legend.position = "none")
ggsave(
  file.path(paths$figure, "GSE132465_qc_retention.pdf"),
  p_qc, width = 6, height = 4
)
ggsave(
  file.path(paths$figure, "GSE132465_qc_retention.png"),
  p_qc, width = 6, height = 4, dpi = 300
)
ggsave(
  file.path(paths$figure, "GSE132465_cell_composition.pdf"),
  p_comp, width = 9, height = 6
)
ggsave(
  file.path(paths$figure, "GSE132465_cell_composition.png"),
  p_comp, width = 9, height = 6, dpi = 300
)
writeLines(
  capture.output(sessionInfo()),
  file.path(paths$result, "GSE132465_sessionInfo.txt")
)
cat(
  "GSE132465_PREPARATION_OK\tretained=", sum(retained),
  "\tepithelial=", sum(retained & major_cell_type == "Epithelial"), "\n", sep = ""
)
