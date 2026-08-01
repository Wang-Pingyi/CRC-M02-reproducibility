#!/usr/bin/env Rscript

# Analysis: Stage 5C major-compartment annotation discovery
# Date: 2026-07-27
# Random seed: read from config/annotation_parameters.tsv
# Input: validated Stage 5B Harmony-integrated Seurat object
# Scope: clustering and annotation evidence only; no differential, trajectory,
# communication, or machine-learning analysis

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
input_file <- file.path(
  project_root, "objects", "GSE201348_harmony_integrated.rds"
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
    section == section_name & parameter == parameter_name, value
  ]
  if (length(value) != 1L) {
    stop("Expected one parameter: ", section_name, "/", parameter_name)
  }
  if (numeric) as.numeric(value) else value
}

seed <- as.integer(get_param("global", "random_seed", TRUE))
set.seed(seed)

obj <- readRDS(input_file)
if (!inherits(obj, "Seurat")) stop("Stage 5B input is not a Seurat object")
if (!"harmony" %in% names(obj@reductions)) {
  stop("Validated Harmony reduction is absent")
}
required_meta <- c(
  "sample_id", "donor_id", "biological_sample_id", "lesion_stage",
  "sporadic_or_FAP"
)
if (!all(required_meta %in% colnames(obj@meta.data))) {
  stop(
    "Missing metadata: ",
    paste(setdiff(required_meta, colnames(obj@meta.data)), collapse = ";")
  )
}

DefaultAssay(obj) <- "RNA"
dims <- seq_len(as.integer(get_param("major", "harmony_dimensions", TRUE)))
obj <- FindNeighbors(
  obj,
  reduction = "harmony",
  dims = dims,
  graph.name = c("major_nn", "major_snn"),
  verbose = FALSE
)
obj <- FindClusters(
  obj,
  graph.name = "major_snn",
  resolution = get_param("major", "cluster_resolution", TRUE),
  random.seed = seed,
  verbose = FALSE
)
obj$major_cluster <- as.character(Idents(obj))
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = dims,
  reduction.name = "major_umap",
  reduction.key = "majorUMAP_",
  n.neighbors = as.integer(get_param("major", "umap_neighbors", TRUE)),
  min.dist = get_param("major", "umap_min_dist", TRUE),
  seed.use = seed,
  verbose = FALSE
)
Idents(obj) <- "major_cluster"

major_markers <- FindAllMarkers(
  obj,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = get_param("major", "marker_min_pct", TRUE),
  logfc.threshold = get_param("major", "marker_logfc", TRUE),
  max.cells.per.ident = as.integer(
    get_param("major", "marker_max_cells_per_cluster", TRUE)
  ),
  random.seed = seed,
  verbose = FALSE
)
major_markers <- as.data.table(major_markers)
fwrite(
  major_markers,
  file.path(result_dir, "major_cluster_markers.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

evidence <- fread(evidence_file)
major_evidence <- evidence[level == "major"]
split_markers <- function(x) {
  unique(trimws(unlist(strsplit(x, ";", fixed = TRUE))))
}
positive_genes <- unique(unlist(lapply(
  major_evidence$positive_markers, split_markers
)))
exclusion_genes <- unique(unlist(lapply(
  major_evidence$exclusion_markers, split_markers
)))
evidence_genes <- intersect(
  unique(c(positive_genes, exclusion_genes)),
  rownames(obj)
)
if (length(evidence_genes) < 20L) {
  stop("Too few major-compartment evidence genes are present")
}

cluster_gene_stats <- function(object, genes, cluster_field) {
  mat <- GetAssayData(object, assay = "RNA", slot = "data")[genes, , drop = FALSE]
  groups <- as.character(object[[cluster_field]][, 1])
  out <- rbindlist(lapply(sort(unique(groups)), function(cluster_id) {
    idx <- which(groups == cluster_id)
    data.table(
      cluster = cluster_id,
      gene = genes,
      avg_expression = as.numeric(Matrix::rowMeans(mat[, idx, drop = FALSE])),
      pct_expressing = as.numeric(Matrix::rowMeans(mat[, idx, drop = FALSE] > 0)),
      cells = length(idx)
    )
  }))
  out
}

major_gene_stats <- cluster_gene_stats(obj, evidence_genes, "major_cluster")
fwrite(
  major_gene_stats,
  file.path(source_dir, "major_marker_dotplot_source.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

score_one_signature <- function(annotation_name, positive, negative, stats) {
  pos <- intersect(split_markers(positive), unique(stats$gene))
  neg <- intersect(split_markers(negative), unique(stats$gene))
  clusters <- sort(unique(stats$cluster))
  rbindlist(lapply(clusters, function(cluster_id) {
    x <- stats[cluster == cluster_id]
    pos_score <- if (length(pos)) mean(x[gene %chin% pos, avg_expression]) else NA_real_
    neg_score <- if (length(neg)) mean(x[gene %chin% neg, avg_expression]) else NA_real_
    pos_detect <- if (length(pos)) mean(x[gene %chin% pos, pct_expressing]) else NA_real_
    neg_detect <- if (length(neg)) mean(x[gene %chin% neg, pct_expressing]) else NA_real_
    data.table(
      cluster = cluster_id,
      annotation = annotation_name,
      positive_genes_present = paste(pos, collapse = ";"),
      exclusion_genes_present = paste(neg, collapse = ";"),
      positive_mean_expression = pos_score,
      exclusion_mean_expression = neg_score,
      positive_mean_detection = pos_detect,
      exclusion_mean_detection = neg_detect
    )
  }))
}

major_signature_scores <- rbindlist(lapply(
  seq_len(nrow(major_evidence)),
  function(i) score_one_signature(
    major_evidence$annotation[i],
    major_evidence$positive_markers[i],
    major_evidence$exclusion_markers[i],
    major_gene_stats
  )
))
major_signature_scores[
  ,
  positive_expression_z := as.numeric(scale(positive_mean_expression)),
  by = annotation
]
major_signature_scores[
  ,
  exclusion_expression_z := as.numeric(scale(exclusion_mean_expression)),
  by = annotation
]
major_signature_scores[
  ,
  evidence_score := positive_expression_z - 0.5 * exclusion_expression_z
]
major_signature_scores[
  !is.finite(evidence_score),
  evidence_score := positive_mean_expression - 0.5 * exclusion_mean_expression
]
fwrite(
  major_signature_scores,
  file.path(result_dir, "major_signature_scores.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

cluster_summary <- as.data.table(obj@meta.data)[
  ,
  .(
    cells = .N,
    libraries = uniqueN(sample_id),
    biological_samples = uniqueN(biological_sample_id),
    donors = uniqueN(donor_id)
  ),
  by = major_cluster
][order(as.integer(as.character(major_cluster)))]
setnames(cluster_summary, "major_cluster", "cluster")

composition <- rbind(
  as.data.table(obj@meta.data)[
    ,
    .N,
    by = .(cluster = major_cluster, level = donor_id)
  ][, variable := "donor_id"],
  as.data.table(obj@meta.data)[
    ,
    .N,
    by = .(cluster = major_cluster, level = sample_id)
  ][, variable := "sample_id"],
  as.data.table(obj@meta.data)[
    ,
    .N,
    by = .(cluster = major_cluster, level = lesion_stage)
  ][, variable := "lesion_stage"],
  as.data.table(obj@meta.data)[
    ,
    .N,
    by = .(cluster = major_cluster, level = sporadic_or_FAP)
  ][, variable := "sporadic_or_FAP"]
)
setcolorder(composition, c("cluster", "variable", "level", "N"))
composition[, fraction := N / sum(N), by = .(cluster, variable)]
fwrite(
  composition,
  file.path(result_dir, "major_cluster_composition.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

fc_column <- intersect(
  c("avg_log2FC", "avg_logFC"),
  colnames(major_markers)
)
if (length(fc_column) != 1L) stop("Cannot identify marker fold-change column")
setorderv(major_markers, c("cluster", fc_column), c(1L, -1L))
top_markers <- major_markers[
  ,
  .(top_markers = paste(head(gene, 15L), collapse = ";")),
  by = .(cluster = as.character(cluster))
]
ranked_scores <- major_signature_scores[order(cluster, -evidence_score)]
top_signatures <- ranked_scores[
  ,
  .(
    suggested_annotation = annotation[1],
    top_signature_scores = paste(
      paste0(head(annotation, 3L), "=", sprintf("%.3f", head(evidence_score, 3L))),
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
  file.path(result_dir, "major_cluster_review_blinded.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

embedding <- Embeddings(obj, "major_umap")
umap_source <- data.table(
  cell_id = rownames(embedding),
  UMAP_1 = embedding[, 1],
  UMAP_2 = embedding[, 2],
  major_cluster = as.character(obj$major_cluster),
  sample_id = obj$sample_id,
  donor_id = obj$donor_id,
  biological_sample_id = obj$biological_sample_id,
  lesion_stage = obj$lesion_stage,
  sporadic_or_FAP = obj$sporadic_or_FAP
)
fwrite(
  umap_source,
  file.path(source_dir, "major_umap_source.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

cluster_levels <- sort(unique(umap_source$major_cluster))
cluster_palette <- setNames(
  grDevices::hcl.colors(length(cluster_levels), "Dynamic"),
  cluster_levels
)
p_cluster <- ggplot(
  umap_source,
  aes(UMAP_1, UMAP_2, color = major_cluster)
) +
  geom_point(size = 0.05, alpha = 0.35) +
  scale_color_manual(values = cluster_palette) +
  labs(
    title = "GSE201348 major-compartment discovery clusters",
    color = "Cluster"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    legend.key.height = grid::unit(3, "mm")
  )
ggsave(
  file.path(figure_dir, "major_discovery_clusters.pdf"),
  p_cluster, width = 7, height = 6
)
ggsave(
  file.path(figure_dir, "major_discovery_clusters.png"),
  p_cluster, width = 7, height = 6, dpi = 300
)

dot_genes <- unique(unlist(lapply(
  major_evidence$positive_markers,
  function(x) head(split_markers(x), 4L)
)))
dot_source <- major_gene_stats[gene %chin% dot_genes]
dot_source[
  ,
  scaled_expression := as.numeric(scale(avg_expression)),
  by = gene
]
p_dot <- ggplot(
  dot_source,
  aes(x = cluster, y = gene, size = pct_expressing, color = scaled_expression)
) +
  geom_point() +
  scale_size(range = c(0.2, 4), limits = c(0, 1)) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
  ) +
  labs(
    title = "Curated marker evidence across discovery clusters",
    x = "Cluster",
    y = "Marker",
    size = "Fraction expressed",
    color = "Scaled mean"
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(figure_dir, "major_marker_dotplot.pdf"),
  p_dot, width = 9, height = 9
)
ggsave(
  file.path(figure_dir, "major_marker_dotplot.png"),
  p_dot, width = 9, height = 9, dpi = 300
)

obj@misc$stage_5C_major_discovery <- list(
  seed = seed,
  parameter_file = "config/annotation_parameters.tsv",
  evidence_file = "metadata/annotation_evidence.tsv",
  annotation_status = "marker_blinded_review_pending",
  differential_expression = "not_performed",
  trajectory_analysis = "not_performed",
  machine_learning = "not_performed"
)
saveRDS(
  obj,
  file.path(object_dir, "GSE201348_5C_major_clustered_unannotated.rds"),
  compress = "gzip"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "major_discovery_sessionInfo.txt")
)
writeLines(
  c(
    paste0("seed=", seed),
    paste0("input_object=", input_file),
    paste0("cells=", ncol(obj)),
    paste0("features=", nrow(obj)),
    paste0("clusters=", uniqueN(obj$major_cluster)),
    "cell_annotation=marker_blinded_review_pending",
    "differential_expression=not_performed",
    "trajectory_analysis=not_performed",
    "machine_learning=not_performed"
  ),
  file.path(result_dir, "major_discovery_provenance.txt")
)
message(
  "Stage 5C major annotation discovery completed: ",
  ncol(obj), " cells; ", uniqueN(obj$major_cluster), " clusters"
)
