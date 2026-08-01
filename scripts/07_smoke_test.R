#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
})

locked <- read_locked_stage7_modules(project_dir)
genes <- unique(locked$membership$gene)
donors <- paste0("D", 1:4)
conditions <- c("normal", "adenoma", "cancer")
metadata <- expand.grid(
  donor_id = donors, condition = conditions, major_cell_type = "Epithelial",
  stringsAsFactors = FALSE
)
metadata$pseudobulk_id <- paste(
  metadata$donor_id, metadata$condition, metadata$major_cell_type, sep = "|"
)
counts <- matrix(
  stats::rpois(length(genes) * nrow(metadata), lambda = 30),
  nrow = length(genes),
  dimnames = list(genes, metadata$pseudobulk_id)
)
first_module <- unique(
  locked$membership$gene[
    locked$membership$module_id == locked$candidates$module_id[[1L]]
  ]
)
counts[first_module, metadata$condition == "adenoma"] <-
  counts[first_module, metadata$condition == "adenoma"] + 20
counts[first_module, metadata$condition == "cancer"] <-
  counts[first_module, metadata$condition == "cancer"] + 40
metadata$n_cells <- 100L
metadata$total_umi <- colSums(counts)

scored <- module_scores_from_counts(
  Matrix(counts, sparse = TRUE), locked$membership, 0.60, 8L
)
if (nrow(scored$scores) != 6L * nrow(metadata)) {
  stop("Smoke-test module-score dimensions are incorrect")
}
fit <- paired_effect(
  scored$scores, metadata, locked$candidates$module_id[[1L]],
  "normal", "adenoma", donors
)
if (fit$summary$n_donors != 4L ||
    !is.finite(fit$summary$effect) ||
    !is.finite(fit$summary$p_value) ||
    !is.finite(fit$summary$lodo_sign_stability)) {
  stop("Smoke-test paired effect failed")
}

specificity_meta <- rbind(
  transform(metadata, major_cell_type = "Epithelial"),
  transform(metadata, major_cell_type = "Myeloid")
)
specificity_meta$pseudobulk_id <- paste(
  specificity_meta$donor_id, specificity_meta$condition,
  specificity_meta$major_cell_type, sep = "|"
)
specificity_scores <- transform(
  specificity_meta,
  module_id = locked$candidates$module_id[[1L]],
  module_score = ifelse(major_cell_type == "Epithelial", 1, -1),
  locked_genes = length(first_module),
  represented_genes = length(first_module),
  gene_coverage = 1,
  evaluable = TRUE
)[, c(
  "pseudobulk_id", "module_id", "module_score", "locked_genes",
  "represented_genes", "gene_coverage", "evaluable"
)]
spec <- specificity_effect(
  specificity_scores, specificity_meta, locked$candidates$module_id[[1L]]
)
if (spec$summary$n_donors != 4L || spec$summary$effect <= 0) {
  stop("Smoke-test epithelial specificity failed")
}
cat("STAGE_7_SMOKE_TEST_OK\n")
