#!/usr/bin/env Rscript

# Resume-only reconstruction of GSE161277 figures after all analytical
# GSE161277 outputs were saved successfully.

set.seed(20260728)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)

suppressPackageStartupMessages(library(ggplot2))

required <- c(
  qc = file.path(paths$result, "GSE161277_qc_summary.tsv"),
  umap = file.path(paths$source, "GSE161277_umap_source.tsv.gz"),
  object = file.path(paths$object, "GSE161277_stage7_annotated.rds"),
  pseudobulk = file.path(
    paths$processed, "GSE161277_stage7_pseudobulk_raw_counts.rds"
  )
)
if (any(!file.exists(required)) || any(file.info(required)$size <= 0)) {
  stop("Cannot resume: a completed GSE161277 analytical output is missing")
}

qc_summary <- utils::read.delim(
  required[["qc"]], check.names = FALSE, stringsAsFactors = FALSE
)
umap_source <- utils::read.delim(
  gzfile(required[["umap"]]), check.names = FALSE, stringsAsFactors = FALSE
)
required_umap <- c(
  "cell_id", "UMAP_1", "UMAP_2", "donor_id", "condition", "major_cell_type"
)
if (!all(required_umap %in% colnames(umap_source))) {
  stop("Saved GSE161277 UMAP source-data schema is incomplete")
}

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
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
cat("GSE161277_FIGURE_RESUME_OK\n")

