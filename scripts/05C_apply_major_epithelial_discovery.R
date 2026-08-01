#!/usr/bin/env Rscript

# Analysis: Stage 5C apply blinded major annotation and discover epithelial states
# Date: 2026-07-27
# Scope: annotation-marker evidence only. No lesion-stage differential analysis,
# trajectory analysis, cell-cell communication analysis, or machine learning.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(harmony)
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}

param_file <- file.path(project_root, "config", "annotation_parameters.tsv")
evidence_file <- file.path(project_root, "metadata", "annotation_evidence.tsv")
mapping_file <- file.path(project_root, "config", "major_cluster_annotation.tsv")
input_file <- file.path(
  project_root, "objects", "GSE201348_5C_major_clustered_unannotated.rds"
)
result_dir <- file.path(project_root, "results", "05C_annotation")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_root, "figures", "05C_annotation")
object_dir <- file.path(project_root, "objects")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

params <- fread(param_file)
get_param <- function(section_name, parameter_name, numeric = FALSE) {
  value <- params[
    section == section_name & parameter == parameter_name,
    value
  ]
  if (length(value) != 1L) {
    stop("Expected one parameter: ", section_name, "/", parameter_name)
  }
  if (numeric) as.numeric(value) else value
}

seed <- as.integer(get_param("global", "random_seed", TRUE))
set.seed(seed)

obj <- readRDS(input_file)
mapping <- fread(mapping_file, colClasses = "character")
required_mapping <- c(
  "cluster", "major_cell_type", "confidence", "evidence_summary",
  "uncertainty_note", "reviewer_context_blinded"
)
if (!all(required_mapping %in% colnames(mapping))) {
  stop(
    "Major annotation map is missing: ",
    paste(setdiff(required_mapping, colnames(mapping)), collapse = ";")
  )
}
if (anyDuplicated(mapping$cluster)) stop("Major annotation map has duplicate clusters")
observed_clusters <- sort(unique(as.character(obj$major_cluster)))
if (!setequal(observed_clusters, mapping$cluster)) {
  stop(
    "Major annotation map does not match object clusters. Missing map: ",
    paste(setdiff(observed_clusters, mapping$cluster), collapse = ";"),
    "; extra map: ",
    paste(setdiff(mapping$cluster, observed_clusters), collapse = ";")
  )
}
allowed_major <- fread(evidence_file)[level == "major", unique(annotation)]
if (any(!mapping$major_cell_type %chin% allowed_major)) {
  stop(
    "Unregistered major annotation: ",
    paste(unique(mapping[!major_cell_type %chin% allowed_major, major_cell_type]),
          collapse = ";")
  )
}
if (any(mapping$reviewer_context_blinded != "yes")) {
  stop("All major-cluster decisions must be recorded as context-blinded")
}

type_map <- setNames(mapping$major_cell_type, mapping$cluster)
confidence_map <- setNames(mapping$confidence, mapping$cluster)
uncertainty_map <- setNames(mapping$uncertainty_note, mapping$cluster)
obj$major_cell_type <- unname(type_map[as.character(obj$major_cluster)])
obj$major_annotation_confidence <- unname(
  confidence_map[as.character(obj$major_cluster)]
)
obj$major_annotation_uncertainty <- unname(
  uncertainty_map[as.character(obj$major_cluster)]
)
if (anyNA(obj$major_cell_type)) stop("Major annotation produced missing labels")

fwrite(
  data.table(
    cell_id = colnames(obj),
    major_cluster = as.character(obj$major_cluster),
    major_cell_type = obj$major_cell_type,
    confidence = obj$major_annotation_confidence,
    uncertainty_note = obj$major_annotation_uncertainty
  ),
  file.path(result_dir, "major_cell_annotations.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

major_embedding <- Embeddings(obj, "major_umap")
major_umap_source <- data.table(
  cell_id = rownames(major_embedding),
  UMAP_1 = major_embedding[, 1],
  UMAP_2 = major_embedding[, 2],
  major_cluster = as.character(obj$major_cluster),
  major_cell_type = obj$major_cell_type,
  donor_id = obj$donor_id,
  lesion_stage = obj$lesion_stage
)
fwrite(
  major_umap_source,
  file.path(source_dir, "major_annotated_umap_source.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

major_levels <- sort(unique(major_umap_source$major_cell_type))
major_palette <- setNames(
  grDevices::hcl.colors(length(major_levels), "Dark 3"), major_levels
)
p_major <- ggplot(
  major_umap_source,
  aes(UMAP_1, UMAP_2, color = major_cell_type)
) +
  geom_point(size = 0.05, alpha = 0.35) +
  scale_color_manual(values = major_palette) +
  labs(
    title = "GSE201348 major cell-type annotation",
    color = "Major cell type"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())
ggsave(
  file.path(figure_dir, "major_cell_type_annotation.pdf"),
  p_major, width = 7.5, height = 6
)
ggsave(
  file.path(figure_dir, "major_cell_type_annotation.png"),
  p_major, width = 7.5, height = 6, dpi = 300
)

obj@misc$stage_5C_major_annotation <- list(
  seed = seed,
  mapping_file = "config/major_cluster_annotation.tsv",
  evidence_file = "metadata/annotation_evidence.tsv",
  context_blinded = TRUE
)
saveRDS(
  obj,
  file.path(object_dir, "GSE201348_5C_major_annotated.rds"),
  compress = "gzip"
)

epithelial <- subset(obj, subset = major_cell_type == "Epithelial")
if (ncol(epithelial) < 1000L) {
  stop("Fewer than 1,000 epithelial cells; review major annotation")
}
DefaultAssay(epithelial) <- "RNA"
epithelial <- DietSeurat(
  epithelial,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = NULL,
  graphs = NULL
)
epithelial <- NormalizeData(
  epithelial,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
epithelial <- FindVariableFeatures(
  epithelial,
  selection.method = "vst",
  nfeatures = as.integer(get_param("epithelial", "variable_features", TRUE)),
  verbose = FALSE
)
epithelial <- ScaleData(
  epithelial,
  features = VariableFeatures(epithelial),
  verbose = FALSE
)
epithelial <- RunPCA(
  epithelial,
  features = VariableFeatures(epithelial),
  npcs = as.integer(get_param("epithelial", "pca_dimensions", TRUE)),
  seed.use = seed,
  verbose = FALSE
)
epithelial_dims <- seq_len(as.integer(
  get_param("epithelial", "harmony_dimensions", TRUE)
))
epithelial <- harmony::RunHarmony(
  epithelial,
  group.by.vars = "sample_id",
  reduction.use = "pca",
  dims.use = epithelial_dims,
  theta = get_param("epithelial", "harmony_theta", TRUE),
  lambda = get_param("epithelial", "harmony_lambda", TRUE),
  reduction.save = "epithelial_harmony",
  assay.use = "RNA",
  project.dim = FALSE,
  plot_convergence = FALSE,
  verbose = FALSE
)
epithelial <- FindNeighbors(
  epithelial,
  reduction = "epithelial_harmony",
  dims = epithelial_dims,
  graph.name = c("epithelial_nn", "epithelial_snn"),
  verbose = FALSE
)
epithelial <- FindClusters(
  epithelial,
  graph.name = "epithelial_snn",
  resolution = get_param("epithelial", "cluster_resolution", TRUE),
  random.seed = seed,
  verbose = FALSE
)
epithelial$epithelial_cluster <- as.character(Idents(epithelial))
epithelial <- RunUMAP(
  epithelial,
  reduction = "epithelial_harmony",
  dims = epithelial_dims,
  reduction.name = "epithelial_umap",
  reduction.key = "epithelialUMAP_",
  n.neighbors = as.integer(get_param("epithelial", "umap_neighbors", TRUE)),
  min.dist = get_param("epithelial", "umap_min_dist", TRUE),
  seed.use = seed,
  verbose = FALSE
)
Idents(epithelial) <- "epithelial_cluster"

epithelial_markers <- FindAllMarkers(
  epithelial,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = get_param("epithelial", "marker_min_pct", TRUE),
  logfc.threshold = get_param("epithelial", "marker_logfc", TRUE),
  max.cells.per.ident = as.integer(
    get_param("epithelial", "marker_max_cells_per_cluster", TRUE)
  ),
  random.seed = seed,
  verbose = FALSE
)
epithelial_markers <- as.data.table(epithelial_markers)
fwrite(
  epithelial_markers,
  file.path(result_dir, "epithelial_cluster_markers.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

evidence <- fread(evidence_file)[
  level == "epithelial" & positive_markers != "NA"
]
split_markers <- function(x) {
  unique(trimws(unlist(strsplit(x, ";", fixed = TRUE))))
}
evidence_genes <- intersect(
  unique(c(
    unlist(lapply(evidence$positive_markers, split_markers)),
    unlist(lapply(evidence$exclusion_markers, split_markers))
  )),
  rownames(epithelial)
)
if (length(evidence_genes) < 25L) {
  stop("Too few epithelial-state evidence genes are present")
}

normalized <- GetAssayData(
  epithelial, assay = "RNA", slot = "data"
)[evidence_genes, , drop = FALSE]
groups <- as.character(epithelial$epithelial_cluster)
gene_stats <- rbindlist(lapply(sort(unique(groups)), function(cluster_id) {
  idx <- which(groups == cluster_id)
  data.table(
    cluster = cluster_id,
    gene = evidence_genes,
    avg_expression = as.numeric(
      Matrix::rowMeans(normalized[, idx, drop = FALSE])
    ),
    pct_expressing = as.numeric(
      Matrix::rowMeans(normalized[, idx, drop = FALSE] > 0)
    ),
    cells = length(idx)
  )
}))
fwrite(
  gene_stats,
  file.path(source_dir, "epithelial_marker_dotplot_source.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

score_one_signature <- function(annotation_name, positive, negative, stats) {
  pos <- intersect(split_markers(positive), unique(stats$gene))
  neg <- intersect(split_markers(negative), unique(stats$gene))
  rbindlist(lapply(sort(unique(stats$cluster)), function(cluster_id) {
    x <- stats[cluster == cluster_id]
    data.table(
      cluster = cluster_id,
      annotation = annotation_name,
      positive_genes_present = paste(pos, collapse = ";"),
      exclusion_genes_present = paste(neg, collapse = ";"),
      positive_mean_expression = mean(x[gene %chin% pos, avg_expression]),
      exclusion_mean_expression = mean(x[gene %chin% neg, avg_expression]),
      positive_mean_detection = mean(x[gene %chin% pos, pct_expressing]),
      exclusion_mean_detection = mean(x[gene %chin% neg, pct_expressing])
    )
  }))
}
signature_scores <- rbindlist(lapply(seq_len(nrow(evidence)), function(i) {
  score_one_signature(
    evidence$annotation[i],
    evidence$positive_markers[i],
    evidence$exclusion_markers[i],
    gene_stats
  )
}))
signature_scores[
  ,
  positive_expression_z := as.numeric(scale(positive_mean_expression)),
  by = annotation
]
signature_scores[
  ,
  exclusion_expression_z := as.numeric(scale(exclusion_mean_expression)),
  by = annotation
]
signature_scores[
  ,
  evidence_score := positive_expression_z - 0.5 * exclusion_expression_z
]
signature_scores[
  !is.finite(evidence_score),
  evidence_score := positive_mean_expression - 0.5 * exclusion_mean_expression
]
fwrite(
  signature_scores,
  file.path(result_dir, "epithelial_signature_scores.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

epithelial_meta <- as.data.table(epithelial@meta.data)
cluster_summary <- epithelial_meta[
  ,
  .(
    cells = .N,
    libraries = uniqueN(sample_id),
    biological_samples = uniqueN(biological_sample_id),
    donors = uniqueN(donor_id)
  ),
  by = .(cluster = epithelial_cluster)
]
composition <- rbind(
  epithelial_meta[
    ,
    .N,
    by = .(cluster = epithelial_cluster, level = donor_id)
  ][, variable := "donor_id"],
  epithelial_meta[
    ,
    .N,
    by = .(
      cluster = epithelial_cluster,
      level = lesion_stage
    )
  ][, variable := "lesion_stage"],
  epithelial_meta[
    ,
    .N,
    by = .(
      cluster = epithelial_cluster,
      level = sporadic_or_FAP
    )
  ][, variable := "sporadic_or_FAP"]
)
setcolorder(composition, c("cluster", "variable", "level", "N"))
composition[, fraction := N / sum(N), by = .(cluster, variable)]
fwrite(
  composition,
  file.path(result_dir, "epithelial_cluster_composition.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"), colnames(epithelial_markers)
)
if (length(fc_column) != 1L) stop("Cannot identify marker fold-change column")
setorderv(epithelial_markers, c("cluster", fc_column), c(1L, -1L))
top_markers <- epithelial_markers[
  ,
  .(top_markers = paste(head(gene, 20L), collapse = ";")),
  by = .(cluster = as.character(cluster))
]
top_signatures <- signature_scores[order(cluster, -evidence_score)][
  ,
  .(
    suggested_annotation = annotation[1],
    top_signature_scores = paste(
      paste0(
        head(annotation, 4L), "=",
        sprintf("%.3f", head(evidence_score, 4L))
      ),
      collapse = ";"
    )
  ),
  by = cluster
]
blinded_review <- Reduce(
  function(x, y) merge(x, y, by = "cluster", all = TRUE),
  list(cluster_summary, top_markers, top_signatures)
)
blinded_review[, review_status := "awaiting_marker_blinded_manual_review"]
fwrite(
  blinded_review,
  file.path(result_dir, "epithelial_cluster_review_blinded.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

embedding <- Embeddings(epithelial, "epithelial_umap")
umap_source <- data.table(
  cell_id = rownames(embedding),
  UMAP_1 = embedding[, 1],
  UMAP_2 = embedding[, 2],
  epithelial_cluster = as.character(epithelial$epithelial_cluster),
  donor_id = epithelial$donor_id,
  biological_sample_id = epithelial$biological_sample_id,
  sample_id = epithelial$sample_id,
  lesion_stage = epithelial$lesion_stage,
  sporadic_or_FAP = epithelial$sporadic_or_FAP
)
fwrite(
  umap_source,
  file.path(source_dir, "epithelial_umap_source.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

cluster_levels <- sort(unique(umap_source$epithelial_cluster))
cluster_palette <- setNames(
  grDevices::hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels
)
p_cluster <- ggplot(
  umap_source,
  aes(UMAP_1, UMAP_2, color = epithelial_cluster)
) +
  geom_point(size = 0.08, alpha = 0.4) +
  scale_color_manual(values = cluster_palette) +
  labs(
    title = "GSE201348 epithelial discovery clusters",
    color = "Cluster"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())
ggsave(
  file.path(figure_dir, "epithelial_discovery_clusters.pdf"),
  p_cluster, width = 7, height = 6
)
ggsave(
  file.path(figure_dir, "epithelial_discovery_clusters.png"),
  p_cluster, width = 7, height = 6, dpi = 300
)

dot_genes <- unique(unlist(lapply(
  evidence$positive_markers,
  function(x) head(split_markers(x), 5L)
)))
dot_source <- gene_stats[gene %chin% dot_genes]
dot_source[
  ,
  scaled_expression := as.numeric(scale(avg_expression)),
  by = gene
]
p_dot <- ggplot(
  dot_source,
  aes(
    x = cluster, y = gene, size = pct_expressing, color = scaled_expression
  )
) +
  geom_point() +
  scale_size(range = c(0.2, 4), limits = c(0, 1)) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
  ) +
  labs(
    title = "Curated epithelial-state marker evidence",
    x = "Cluster",
    y = "Marker",
    size = "Fraction expressed",
    color = "Scaled mean"
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(figure_dir, "epithelial_marker_dotplot.pdf"),
  p_dot, width = 9, height = 10
)
ggsave(
  file.path(figure_dir, "epithelial_marker_dotplot.png"),
  p_dot, width = 9, height = 10, dpi = 300
)

epithelial@misc$stage_5C_epithelial_discovery <- list(
  seed = seed,
  parameter_file = "config/annotation_parameters.tsv",
  evidence_file = "metadata/annotation_evidence.tsv",
  annotation_status = "marker_blinded_review_pending",
  context_variables_hidden_during_marker_review = c(
    "donor_id", "sample_id", "lesion_stage", "sporadic_or_FAP"
  )
)
saveRDS(
  epithelial,
  file.path(
    object_dir, "GSE201348_5C_epithelial_clustered_unannotated.rds"
  ),
  compress = "gzip"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "epithelial_discovery_sessionInfo.txt")
)
writeLines(
  c(
    paste0("seed=", seed),
    paste0("major_input_object=", input_file),
    paste0("epithelial_cells=", ncol(epithelial)),
    paste0("epithelial_features=", nrow(epithelial)),
    paste0("epithelial_clusters=", uniqueN(epithelial$epithelial_cluster)),
    "lesion_stage_used_for_annotation=no",
    "donor_identity_used_for_annotation=no",
    "differential_expression_by_stage=not_performed",
    "trajectory_analysis=not_performed",
    "machine_learning=not_performed"
  ),
  file.path(result_dir, "epithelial_discovery_provenance.txt")
)
message(
  "Major annotation applied and epithelial discovery completed: ",
  ncol(epithelial), " epithelial cells; ",
  uniqueN(epithelial$epithelial_cluster), " epithelial clusters"
)
