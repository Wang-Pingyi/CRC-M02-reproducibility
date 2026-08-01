#!/usr/bin/env Rscript

# Independent Stage 5C epithelial PCA/Harmony audit.

suppressPackageStartupMessages({
  library(data.table)
  library(RANN)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}
param_file <- file.path(project_root, "config", "annotation_parameters.tsv")
object_file <- file.path(
  project_root, "objects",
  "GSE201348_5C_epithelial_clustered_unannotated.rds"
)
result_dir <- file.path(project_root, "results", "05C_annotation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

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

obj <- readRDS(object_file)
pca <- Embeddings(obj, "pca")
harmony <- Embeddings(obj, "epithelial_harmony")
if (nrow(pca) != ncol(obj) || nrow(harmony) != ncol(obj)) {
  stop("Embedding rows do not match epithelial cells")
}
if (ncol(pca) < 30L || ncol(harmony) < 30L) {
  stop("PCA or Harmony has fewer than 30 dimensions")
}
if (any(!is.finite(pca)) || any(!is.finite(harmony))) {
  stop("PCA or Harmony contains non-finite values")
}

meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
if (!identical(rownames(pca), meta$cell_id) ||
    !identical(rownames(harmony), meta$cell_id)) {
  stop("Embedding and metadata cell order mismatch")
}
cap <- as.integer(get_param(
  "integration_audit", "diagnostic_cell_cap", TRUE
))
diagnostic_index <- if (nrow(meta) <= cap) {
  seq_len(nrow(meta))
} else {
  set.seed(seed)
  meta[
    ,
    .I[sample(
      .N,
      min(
        .N,
        max(50L, floor(as.numeric(cap) * .N / nrow(meta)))
      )
    )],
    by = lesion_stage
  ]$V1
}
diagnostic_index <- sort(unique(diagnostic_index))
if (length(diagnostic_index) > cap) {
  set.seed(seed)
  diagnostic_index <- sort(sample(diagnostic_index, cap))
}
diag_meta <- meta[diagnostic_index]
diag_pca <- pca[diagnostic_index, seq_len(30L), drop = FALSE]
diag_harmony <- harmony[diagnostic_index, seq_len(30L), drop = FALSE]

eta_squared <- function(embedding, group) {
  overall <- colMeans(embedding)
  total_ss <- colSums(
    (embedding - rep(overall, each = nrow(embedding)))^2
  )
  between_ss <- rep(0, ncol(embedding))
  for (level in unique(group)) {
    idx <- which(group == level)
    center <- colMeans(embedding[idx, , drop = FALSE])
    between_ss <- between_ss + length(idx) * (center - overall)^2
  }
  mean(
    ifelse(total_ss > 0, between_ss / total_ss, NA_real_),
    na.rm = TRUE
  )
}

audit_vars <- c(
  "sample_id", "donor_id", "lesion_stage", "sporadic_or_FAP"
)
variance_audit <- rbindlist(lapply(audit_vars, function(variable_name) {
  data.table(
    variable = variable_name,
    representation = c("pca", "epithelial_harmony"),
    mean_eta_squared = c(
      eta_squared(diag_pca, diag_meta[[variable_name]]),
      eta_squared(diag_harmony, diag_meta[[variable_name]])
    )
  )
}))
fwrite(
  variance_audit,
  file.path(result_dir, "epithelial_integration_variance_audit.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

neighbor_audit_one <- function(
  embedding, metadata, representation, k
) {
  nn <- RANN::nn2(
    embedding, k = min(k + 1L, nrow(embedding))
  )$nn.idx[, -1, drop = FALSE]
  rbindlist(lapply(audit_vars, function(variable_name) {
    labels <- as.character(metadata[[variable_name]])
    same <- rowMeans(
      matrix(labels[nn], nrow = nrow(nn)) == labels
    )
    frequencies <- prop.table(table(labels))
    expected <- sum(frequencies^2)
    data.table(
      representation = representation,
      variable = variable_name,
      mean_same_neighbor_fraction = mean(same),
      random_expected_same_fraction = expected,
      excess_same_fraction = mean(same) - expected
    )
  }))
}
k <- as.integer(get_param("integration_audit", "knn_k", TRUE))
neighbor_audit <- rbind(
  neighbor_audit_one(diag_pca, diag_meta, "pca", k),
  neighbor_audit_one(
    diag_harmony, diag_meta, "epithelial_harmony", k
  )
)
fwrite(
  neighbor_audit,
  file.path(result_dir, "epithelial_integration_neighbor_audit.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

sample_pre <- variance_audit[
  variable == "sample_id" & representation == "pca",
  mean_eta_squared
]
sample_post <- variance_audit[
  variable == "sample_id" & representation == "epithelial_harmony",
  mean_eta_squared
]
stage_pre <- variance_audit[
  variable == "lesion_stage" & representation == "pca",
  mean_eta_squared
]
stage_post <- variance_audit[
  variable == "lesion_stage" &
    representation == "epithelial_harmony",
  mean_eta_squared
]
stage_retention <- if (stage_pre > 0) stage_post / stage_pre else NA_real_
stage_neighbor_excess <- neighbor_audit[
  variable == "lesion_stage" &
    representation == "epithelial_harmony",
  excess_same_fraction
]
sample_effect_reduced <- is.finite(sample_pre) &&
  is.finite(sample_post) && sample_post < sample_pre
overcorrection_guard <- is.finite(stage_retention) &&
  stage_retention >= get_param(
    "integration_audit", "min_stage_eta2_retention", TRUE
  ) &&
  stage_neighbor_excess >= get_param(
    "integration_audit", "min_stage_neighbor_excess", TRUE
  )
guard <- data.table(
  metric = c(
    "sample_eta2_pre", "sample_eta2_post",
    "sample_effect_reduced", "stage_eta2_pre", "stage_eta2_post",
    "stage_eta2_retention", "stage_neighbor_excess_post",
    "overcorrection_guard", "finite_30d_embeddings"
  ),
  value = c(
    sample_pre, sample_post, as.numeric(sample_effect_reduced),
    stage_pre, stage_post, stage_retention, stage_neighbor_excess,
    as.numeric(overcorrection_guard), 1
  )
)
fwrite(
  guard,
  file.path(result_dir, "epithelial_integration_guard.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

writeLines(
  capture.output(str(obj[["epithelial_harmony"]]@misc)),
  file.path(result_dir, "epithelial_harmony_misc_structure.txt")
)
writeLines(
  c(
    paste0("diagnostic_cells=", length(diagnostic_index)),
    paste0("sample_effect_reduced=", sample_effect_reduced),
    paste0("stage_eta2_retention=", stage_retention),
    paste0("stage_neighbor_excess_post=", stage_neighbor_excess),
    paste0("overcorrection_guard=", overcorrection_guard)
  ),
  file.path(result_dir, "epithelial_integration_audit_summary.txt")
)
message(
  "Epithelial integration audit: sample_effect_reduced=",
  sample_effect_reduced,
  "; stage_retention=", signif(stage_retention, 4),
  "; stage_neighbor_excess=", signif(stage_neighbor_excess, 4),
  "; overcorrection_guard=", overcorrection_guard
)
if (!sample_effect_reduced || !overcorrection_guard) {
  quit(save = "no", status = 1L)
}
