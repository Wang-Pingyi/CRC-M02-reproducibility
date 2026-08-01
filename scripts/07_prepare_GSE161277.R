#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
shared_library <- file.path(project_dir, "environment", "R-library")
.libPaths(c(shared_library, .libPaths()))
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)
get_param <- read_stage7_parameters(project_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(ggplot2)
  library(patchwork)
})

min_count <- get_param("qc", "min_nCount", TRUE)
min_feature <- get_param("qc", "min_nFeature", TRUE)
max_mt <- get_param("qc", "max_percent_mt", TRUE)
lower_k <- get_param("qc", "lower_mad_k", TRUE)
upper_k <- get_param("qc", "upper_mad_k", TRUE)
mt_k <- get_param("qc", "mt_mad_k", TRUE)
annotation_margin <- get_param("annotation", "min_cluster_score_margin", TRUE)
resolution <- get_param("annotation", "resolution_GSE161277", TRUE)
dims_use <- seq_len(as.integer(get_param("annotation", "pca_dimensions", TRUE)))

raw_dir <- file.path(project_dir, "data_raw", "GSE161277")
matrix_files <- sort(list.files(
  raw_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE
))
if (length(matrix_files) != 13L) stop("Expected 13 GSE161277 matrix files")

parse_sample <- function(matrix_file) {
  name <- basename(matrix_file)
  sample_id <- sub("_matrix\\.mtx\\.gz$", "", name)
  gsm <- sub("_.*$", "", sample_id)
  label <- sub("^[^_]+_", "", sample_id)
  donor <- sub("^(Patient[0-9]+).*", "\\1", label)
  condition <- sub("^Patient[0-9]+_", "", label)
  condition <- sub("^adenoma_[12]$", "adenoma", condition)
  condition <- sub("^para-cancer$", "para_cancer", condition)
  data.frame(
    sample_id = sample_id, gsm = gsm, donor_id = donor,
    condition = condition, stringsAsFactors = FALSE
  )
}

retained_counts <- list()
retained_metadata <- list()
qc_summaries <- list()
qc_cell_source <- list()

for (matrix_file in matrix_files) {
  info <- parse_sample(matrix_file)
  prefix <- sub("_matrix\\.mtx\\.gz$", "", matrix_file)
  message("Reading and independently QCing ", info$sample_id)
  counts <- Seurat::ReadMtx(
    mtx = matrix_file,
    features = paste0(prefix, "_features.tsv.gz"),
    cells = paste0(prefix, "_barcodes.tsv.gz"),
    feature.column = 2L,
    unique.features = TRUE
  )
  original_barcodes <- colnames(counts)
  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  percent_mt <- percent_feature_set(counts, "^MT-")
  percent_ribo <- percent_feature_set(counts, "^RP[SL]")

  feature_bounds <- robust_qc_bounds(log10(n_feature + 1), lower_k, upper_k)
  count_bounds <- robust_qc_bounds(log10(n_count + 1), lower_k, upper_k)
  mt_bounds <- robust_qc_bounds(percent_mt, lower_k, mt_k)
  adaptive_min_feature <- max(
    min_feature, floor(10^feature_bounds[["lower"]] - 1)
  )
  adaptive_min_count <- max(min_count, floor(10^count_bounds[["lower"]] - 1))
  adaptive_max_feature <- ceiling(10^feature_bounds[["upper"]] - 1)
  adaptive_max_count <- ceiling(10^count_bounds[["upper"]] - 1)
  adaptive_max_mt <- min(max_mt, mt_bounds[["upper"]])
  if (!is.finite(adaptive_max_mt)) adaptive_max_mt <- max_mt

  low_quality <- n_count < adaptive_min_count |
    n_feature < adaptive_min_feature |
    percent_mt > adaptive_max_mt
  high_complexity <- n_count > adaptive_max_count |
    n_feature > adaptive_max_feature
  doublet_input <- !(low_quality | high_complexity)
  doublet <- rep(FALSE, ncol(counts))
  doublet_score <- rep(NA_real_, ncol(counts))
  doublet_class <- rep("not_evaluated_primary_qc", ncol(counts))
  if (sum(doublet_input) >= 100L && info$condition != "blood") {
    sce <- SingleCellExperiment(list(
      counts = counts[, doublet_input, drop = FALSE]
    ))
    colData(sce)$capture <- info$sample_id
    expected_rate <- min(0.08, 0.08 * ncol(sce) / 10000)
    sce <- scDblFinder(
      sce, samples = "capture", dbr = expected_rate, verbose = FALSE
    )
    doublet[doublet_input] <- colData(sce)$scDblFinder.class == "doublet"
    doublet_score[doublet_input] <- colData(sce)$scDblFinder.score
    doublet_class[doublet_input] <- as.character(
      colData(sce)$scDblFinder.class
    )
    rm(sce)
  }
  retained <- doublet_input & !doublet & info$condition != "blood"
  reason <- ifelse(
    info$condition == "blood", "blood_out_of_scope",
    ifelse(
      low_quality, "low_quality",
      ifelse(high_complexity, "high_complexity",
             ifelse(doublet, "doublet", "retained"))
    )
  )

  unique_cells <- paste(info$sample_id, original_barcodes, sep = "__")
  cell_source <- data.frame(
    cell_id = unique_cells,
    original_barcode = original_barcodes,
    sample_id = info$sample_id,
    gsm = info$gsm,
    donor_id = info$donor_id,
    condition = info$condition,
    nCount = as.numeric(n_count),
    nFeature = as.numeric(n_feature),
    percent_mt = percent_mt,
    percent_ribo = percent_ribo,
    scDblFinder_score = doublet_score,
    scDblFinder_class = doublet_class,
    exclusion_reason = reason,
    stringsAsFactors = FALSE
  )
  qc_cell_source[[info$sample_id]] <- cell_source
  qc_summaries[[info$sample_id]] <- data.frame(
    sample_id = info$sample_id,
    gsm = info$gsm,
    donor_id = info$donor_id,
    condition = info$condition,
    cells_input = ncol(counts),
    cells_primary_qc = sum(doublet_input),
    cells_retained = sum(retained),
    low_quality = sum(low_quality),
    high_complexity = sum(high_complexity & !low_quality),
    doublets = sum(doublet),
    blood_out_of_scope = if (info$condition == "blood") ncol(counts) else 0L,
    adaptive_min_nCount = adaptive_min_count,
    adaptive_min_nFeature = adaptive_min_feature,
    adaptive_max_nCount = adaptive_max_count,
    adaptive_max_nFeature = adaptive_max_feature,
    adaptive_max_percent_mt = adaptive_max_mt,
    ribosomal_action = "diagnostic_only",
    doublet_action = if (info$condition == "blood") {
      "not_run_out_of_scope"
    } else "scDblFinder_per_capture",
    stringsAsFactors = FALSE
  )
  if (any(retained)) {
    kept <- counts[, retained, drop = FALSE]
    colnames(kept) <- unique_cells[retained]
    retained_counts[[info$sample_id]] <- kept
    retained_metadata[[info$sample_id]] <- cell_source[
      retained,
      c("cell_id", "sample_id", "gsm", "donor_id", "condition",
        "nCount", "nFeature", "percent_mt", "percent_ribo",
        "scDblFinder_score"),
      drop = FALSE
    ]
  }
  rm(counts)
  invisible(gc())
}

qc_summary <- do.call(rbind, qc_summaries)
qc_cells <- do.call(rbind, qc_cell_source)
write_stage7_tsv(
  qc_summary, file.path(paths$result, "GSE161277_qc_summary.tsv")
)
write_stage7_tsv(
  qc_cells, file.path(paths$source, "GSE161277_qc_cell_source.tsv.gz")
)

message("Combining retained GSE161277 tissue cells within this cohort only")
combined_counts <- do.call(cbind, retained_counts)
combined_meta <- do.call(rbind, retained_metadata)
rownames(combined_meta) <- combined_meta$cell_id
if (!identical(colnames(combined_counts), rownames(combined_meta))) {
  stop("GSE161277 count/metadata order mismatch")
}

obj <- CreateSeuratObject(
  counts = combined_counts, meta.data = combined_meta,
  project = "GSE161277_stage7", min.cells = 0, min.features = 0
)
locked <- read_locked_stage7_modules(project_dir)
candidate_genes <- unique(locked$membership$gene)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj), candidate_genes)
if (length(VariableFeatures(obj)) < 1000L) {
  stop("Too few non-candidate variable genes remain for independent annotation")
}
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = max(dims_use), verbose = FALSE)
obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = resolution, random.seed = stage7_seed, verbose = FALSE)
obj <- RunUMAP(
  obj, dims = dims_use, seed.use = stage7_seed,
  reduction.name = "gse161277_umap", verbose = FALSE
)

signatures <- major_marker_signatures(candidate_genes)
log_expression <- GetAssayData(obj, assay = "RNA", layer = "data")
scored <- score_signatures(log_expression, signatures)
cluster_assignment <- assign_cluster_annotations(
  scored$scores, obj$seurat_clusters, margin = annotation_margin
)
obj$major_cell_type <- unname(
  cluster_assignment$labels[as.character(obj$seurat_clusters)]
)
obj$annotation_status <- ifelse(
  obj$major_cell_type == "Uncertain",
  "uncertain_excluded_primary", "canonical_marker_cluster_annotation"
)

cluster_evidence <- as.data.frame(cluster_assignment$cluster_scores)
cluster_evidence$cluster <- rownames(cluster_evidence)
cluster_evidence$assigned_major_cell_type <- unname(
  cluster_assignment$labels[cluster_evidence$cluster]
)
cluster_evidence <- cluster_evidence[
  , c("cluster", "assigned_major_cell_type", names(signatures))
]
signature_audit <- data.frame(
  major_cell_type = names(scored$genes),
  genes_used = vapply(scored$genes, paste, collapse = ";", character(1)),
  n_genes_used = lengths(scored$genes),
  candidate_genes_removed = vapply(
    major_marker_signatures(character()), function(x) {
      paste(intersect(x, candidate_genes), collapse = ";")
    }, character(1)
  ),
  annotation_use = "independent_canonical_marker_evidence",
  stringsAsFactors = FALSE
)
write_stage7_tsv(
  cluster_evidence,
  file.path(paths$result, "GSE161277_annotation_evidence.tsv")
)
write_stage7_tsv(
  signature_audit,
  file.path(paths$result, "GSE161277_annotation_signature_audit.tsv")
)

meta <- obj[[]]
meta$cell_id <- rownames(meta)
composition <- aggregate(
  cell_id ~ donor_id + condition + sample_id + major_cell_type,
  data = meta, FUN = length
)
names(composition)[names(composition) == "cell_id"] <- "n_cells"
write_stage7_tsv(
  composition, file.path(paths$result, "GSE161277_cell_composition.tsv")
)

group <- paste(meta$donor_id, meta$condition, meta$major_cell_type, sep = "|")
pseudobulk_counts <- aggregate_counts(
  GetAssayData(obj, assay = "RNA", layer = "counts"), group
)
group_meta <- unique(meta[, c("donor_id", "condition", "major_cell_type")])
group_meta$pseudobulk_id <- paste(
  group_meta$donor_id, group_meta$condition, group_meta$major_cell_type, sep = "|"
)
group_meta <- group_meta[match(colnames(pseudobulk_counts), group_meta$pseudobulk_id), ]
group_meta$n_cells <- as.integer(table(group)[group_meta$pseudobulk_id])
group_meta$total_umi <- as.numeric(Matrix::colSums(pseudobulk_counts))
group_meta$cohort <- "GSE161277"
group_meta$primary_contrast_eligible <- group_meta$condition %in%
  c("normal", "adenoma", "carcinoma") &
  group_meta$donor_id %in% c("Patient1", "Patient2", "Patient3")
group_meta$condition[group_meta$condition == "carcinoma"] <- "cancer"

saveRDS(
  list(counts = pseudobulk_counts, metadata = group_meta),
  file.path(paths$processed, "GSE161277_stage7_pseudobulk_raw_counts.rds"),
  compress = "gzip"
)
obj@misc$stage7 <- list(
  seed = stage7_seed,
  cohort = "GSE161277",
  candidate_reselection = FALSE,
  module_genes_removed_from_annotation_signatures = TRUE,
  primary_unit = "donor",
  doublet_method = "scDblFinder independently per capture",
  ribosomal_action = "diagnostic_only"
)
saveRDS(
  obj, file.path(paths$object, "GSE161277_stage7_annotated.rds"),
  compress = "gzip"
)

umap <- Embeddings(obj, "gse161277_umap")
umap_source <- data.frame(
  cell_id = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
  donor_id = obj$donor_id, condition = obj$condition,
  major_cell_type = obj$major_cell_type, stringsAsFactors = FALSE
)
write_stage7_tsv(
  umap_source, file.path(paths$source, "GSE161277_umap_source.tsv.gz")
)
p_umap <- ggplot(
  umap_source, aes(UMAP_1, UMAP_2, color = major_cell_type)
) +
  geom_point(size = 0.08, alpha = 0.45) +
  labs(title = "GSE161277 independent major-cell annotation", color = "Cell type") +
  theme_bw(base_size = 9) + theme(panel.grid = element_blank())
p_qc <- ggplot(
  qc_summary, aes(condition, cells_retained, fill = condition)
) +
  geom_col() + facet_wrap(~ donor_id, scales = "free_x") +
  labs(title = "GSE161277 retained cells after per-capture QC", y = "Cells") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(paths$figure, "GSE161277_annotation_umap.pdf"),
  p_umap, width = 8, height = 6
)
ggsave(
  file.path(paths$figure, "GSE161277_annotation_umap.png"),
  p_umap, width = 8, height = 6, dpi = 300
)
ggsave(
  file.path(paths$figure, "GSE161277_qc_retention.pdf"),
  p_qc, width = 8, height = 5
)
ggsave(
  file.path(paths$figure, "GSE161277_qc_retention.png"),
  p_qc, width = 8, height = 5, dpi = 300
)
writeLines(
  capture.output(sessionInfo()),
  file.path(paths$result, "GSE161277_sessionInfo.txt")
)
cat(
  "GSE161277_PREPARATION_OK\tcells=", ncol(obj),
  "\tepithelial=", sum(obj$major_cell_type == "Epithelial"), "\n", sep = ""
)
