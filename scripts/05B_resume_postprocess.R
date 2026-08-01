#!/usr/bin/env Rscript

# Resume only the small-table, integration-audit and figure steps after the
# full QC and Harmony objects have already been written successfully.

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(RANN)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else getwd()
config_file <- file.path(project_dir, "config", "qc_parameters.tsv")
result_dir <- file.path(project_dir, "results", "05B_full_qc_integration")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_dir, "figures", "05B_full_qc_integration")
object_dir <- file.path(project_dir, "objects")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

params <- fread(config_file, na.strings = c("NA", ""))
get_global <- function(name, numeric = FALSE) {
  value <- params[scope == "global" & parameter == name, value]
  if (length(value) != 1L) stop("Expected one global parameter: ", name)
  if (numeric) as.numeric(value) else value
}
seed <- as.integer(get_global("random_seed", TRUE))
set.seed(seed)
n_pcs <- as.integer(get_global("pca_dimensions", TRUE))

unintegrated_file <- file.path(object_dir, "GSE201348_unintegrated_normalized.rds")
integrated_file <- file.path(object_dir, "GSE201348_harmony_integrated.rds")
if (!file.exists(unintegrated_file) || !file.exists(integrated_file)) {
  stop("Both completed Stage 5B objects are required for resume")
}
message("Loading completed Stage 5B objects")
unintegrated <- readRDS(unintegrated_file)
integrated <- readRDS(integrated_file)
qc_summary_file <- file.path(result_dir, "library_qc_summary.tsv")
if (!file.exists(qc_summary_file)) {
  stop("Library QC summary is required to validate resumed object dimensions")
}
expected_cells <- sum(fread(qc_summary_file)$cells_after)
if (ncol(unintegrated) != expected_cells || ncol(integrated) != expected_cells) {
  stop("Unexpected cell count in resumed objects")
}
if (!"pca" %in% Reductions(unintegrated) || !"harmony" %in% Reductions(integrated)) {
  stop("Required PCA or Harmony reduction is absent")
}

patient_summary <- as.data.table(integrated@meta.data)[
  , .(
    cells = .N,
    libraries = uniqueN(sample_id),
    biological_tissues = uniqueN(biological_sample_id),
    stages = paste(sort(unique(lesion_stage)), collapse = ";"),
    fap_statuses = paste(sort(unique(sporadic_or_FAP)), collapse = ";"),
    median_nCount = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    median_percent_mt = median(percent_mt),
    median_percent_ribo = median(percent_ribo)
  ),
  by = donor_id
][order(donor_id)]
fwrite(patient_summary, file.path(result_dir, "donor_cell_qc_summary.tsv"), sep = "\t", na = "NA")

biological_sample_summary <- as.data.table(integrated@meta.data)[
  , .(
    donor_id = unique(donor_id),
    cells = .N,
    libraries = uniqueN(sample_id),
    condition = paste(sort(unique(condition)), collapse = ";"),
    lesion_stage = paste(sort(unique(lesion_stage)), collapse = ";"),
    sporadic_or_FAP = paste(sort(unique(sporadic_or_FAP)), collapse = ";"),
    technical_replicate_group = unique(technical_replicate_group),
    qc_risk = paste(sort(unique(qc_risk)), collapse = ";")
  ),
  by = biological_sample_id
][order(donor_id, biological_sample_id)]
fwrite(biological_sample_summary, file.path(result_dir, "biological_sample_cell_summary.tsv"),
       sep = "\t", na = "NA")

distribution_vars <- c("sample_id", "donor_id", "platform", "lesion_stage", "sporadic_or_FAP")
distribution_table <- rbindlist(lapply(distribution_vars, function(v) {
  x <- as.data.table(integrated@meta.data)[, .N, by = v]
  setnames(x, v, "level")
  x[, variable := v]
  rbind(
    copy(x)[, representation := "unintegrated"],
    copy(x)[, representation := "harmony"]
  )
}))
setcolorder(distribution_table, c("representation", "variable", "level", "N"))
fwrite(distribution_table, file.path(result_dir, "pre_post_distribution_audit.tsv"),
       sep = "\t", na = "NA")

meta <- as.data.table(integrated@meta.data, keep.rownames = "cell_id")
pca <- Embeddings(unintegrated, "pca")[, seq_len(n_pcs), drop = FALSE]
harm <- Embeddings(integrated, "harmony")[, seq_len(n_pcs), drop = FALSE]
if (!identical(rownames(pca), meta$cell_id) || !identical(rownames(harm), meta$cell_id)) {
  stop("Embedding and metadata cell order mismatch")
}
if (any(!is.finite(pca)) || any(!is.finite(harm))) {
  stop("Non-finite PCA or Harmony coordinates detected")
}

cap <- as.integer(get_global("diagnostic_cell_cap", TRUE))
diagnostic_index <- if (nrow(meta) <= cap) seq_len(nrow(meta)) else {
  meta[, .I[sample(.N, min(.N, max(50L, floor(as.numeric(cap) * .N / nrow(meta)))))],
       by = lesion_stage]$V1
}
diagnostic_index <- sort(unique(diagnostic_index))
if (length(diagnostic_index) > cap) {
  set.seed(seed)
  diagnostic_index <- sort(sample(diagnostic_index, cap))
}
diag_meta <- meta[diagnostic_index]
diag_pca <- pca[diagnostic_index, , drop = FALSE]
diag_harm <- harm[diagnostic_index, , drop = FALSE]

eta_squared <- function(embedding, group) {
  overall <- colMeans(embedding)
  total_ss <- colSums((embedding - rep(overall, each = nrow(embedding)))^2)
  between_ss <- rep(0, ncol(embedding))
  for (lev in unique(group)) {
    idx <- which(group == lev)
    center <- colMeans(embedding[idx, , drop = FALSE])
    between_ss <- between_ss + length(idx) * (center - overall)^2
  }
  mean(ifelse(total_ss > 0, between_ss / total_ss, NA_real_), na.rm = TRUE)
}

audit_vars <- c("sample_id", "donor_id", "platform", "lesion_stage", "sporadic_or_FAP")
variance_audit <- rbindlist(lapply(audit_vars, function(v) {
  data.table(
    variable = v,
    representation = c("unintegrated", "harmony"),
    mean_eta_squared = c(
      eta_squared(diag_pca, diag_meta[[v]]),
      eta_squared(diag_harm, diag_meta[[v]])
    )
  )
}))
fwrite(variance_audit, file.path(result_dir, "embedding_variance_audit.tsv"),
       sep = "\t", na = "NA")

neighbor_audit_one <- function(embedding, metadata, representation, k) {
  nn <- RANN::nn2(embedding, k = min(k + 1L, nrow(embedding)))$nn.idx[, -1, drop = FALSE]
  rbindlist(lapply(audit_vars, function(v) {
    labels <- as.character(metadata[[v]])
    same <- rowMeans(matrix(labels[nn], nrow = nrow(nn)) == labels)
    frequencies <- prop.table(table(labels))
    expected <- sum(frequencies^2)
    data.table(
      representation = representation,
      variable = v,
      mean_same_neighbor_fraction = mean(same),
      random_expected_same_fraction = expected,
      excess_same_fraction = mean(same) - expected
    )
  }))
}
k <- as.integer(get_global("knn_k", TRUE))
neighbor_audit <- rbind(
  neighbor_audit_one(diag_pca, diag_meta, "unintegrated", k),
  neighbor_audit_one(diag_harm, diag_meta, "harmony", k)
)
fwrite(neighbor_audit, file.path(result_dir, "neighbor_mixing_audit.tsv"),
       sep = "\t", na = "NA")

stage_eta_pre <- variance_audit[variable == "lesion_stage" & representation == "unintegrated",
                                mean_eta_squared]
stage_eta_post <- variance_audit[variable == "lesion_stage" & representation == "harmony",
                                 mean_eta_squared]
stage_excess_post <- neighbor_audit[
  variable == "lesion_stage" & representation == "harmony", excess_same_fraction
]
stage_eta_retention <- if (stage_eta_pre > 0) stage_eta_post / stage_eta_pre else NA_real_
guard_pass <- is.finite(stage_eta_retention) &&
  stage_eta_retention >= get_global("min_stage_eta2_retention", TRUE) &&
  stage_excess_post >= get_global("min_stage_neighbor_excess", TRUE)
overcorrection <- data.table(
  metric = c("stage_eta2_pre", "stage_eta2_post", "stage_eta2_retention",
             "stage_neighbor_excess_post", "overcorrection_guard"),
  value = c(stage_eta_pre, stage_eta_post, stage_eta_retention,
            stage_excess_post, as.numeric(guard_pass)),
  interpretation = c(
    "Mean lesion-stage variance fraction in PCA",
    "Mean lesion-stage variance fraction after Harmony",
    "Fraction of pre-integration lesion-stage signal retained",
    "Same-stage neighbor fraction above random expectation after Harmony",
    "1=guard passed; 0=possible overcorrection"
  )
)
fwrite(overcorrection, file.path(result_dir, "overcorrection_guard.tsv"),
       sep = "\t", na = "NA")

lodo <- meta[, .(
  cells = .N,
  libraries = uniqueN(sample_id),
  biological_tissues = uniqueN(biological_sample_id),
  lesion_stages = paste(sort(unique(lesion_stage)), collapse = ";"),
  fap_statuses = paste(sort(unique(sporadic_or_FAP)), collapse = ";")
), by = .(test_donor = donor_id)][order(test_donor)]
lodo[, fold_id := sprintf("LODO_%02d", seq_len(.N))]
setcolorder(lodo, c("fold_id", "test_donor", "cells", "libraries",
                    "biological_tissues", "lesion_stages", "fap_statuses"))
fwrite(lodo, file.path(result_dir, "leave_one_donor_out_folds.tsv"),
       sep = "\t", na = "NA")

plot_data <- data.table(
  cell_id = diag_meta$cell_id,
  sample_id = diag_meta$sample_id,
  donor_id = diag_meta$donor_id,
  lesion_stage = diag_meta$lesion_stage,
  sporadic_or_FAP = diag_meta$sporadic_or_FAP,
  PCA_1 = diag_pca[, 1],
  PCA_2 = diag_pca[, 2],
  Harmony_1 = diag_harm[, 1],
  Harmony_2 = diag_harm[, 2]
)
fwrite(plot_data, file.path(source_dir, "integration_embedding_plot_source.tsv.gz"),
       sep = "\t", na = "NA")

stage_palette <- c(normal = "#0072B2", adenoma_polyp = "#E69F00", cancer = "#D55E00")
p_pre <- ggplot(plot_data, aes(PCA_1, PCA_2, color = lesion_stage)) +
  geom_point(size = 0.25, alpha = 0.35) +
  scale_color_manual(values = stage_palette) +
  labs(title = "Before integration", color = "Lesion stage") +
  theme_bw(base_size = 10)
p_post <- ggplot(plot_data, aes(Harmony_1, Harmony_2, color = lesion_stage)) +
  geom_point(size = 0.25, alpha = 0.35) +
  scale_color_manual(values = stage_palette) +
  labs(title = "After Harmony", color = "Lesion stage") +
  theme_bw(base_size = 10)
p_stage <- p_pre + p_post + plot_annotation(
  title = "GSE201348 lesion-stage structure before and after integration"
)
ggsave(file.path(figure_dir, "pre_post_integration_by_stage.pdf"), p_stage,
       width = 12, height = 5.5)
ggsave(file.path(figure_dir, "pre_post_integration_by_stage.png"), p_stage,
       width = 12, height = 5.5, dpi = 300)

p_fap_pre <- ggplot(plot_data, aes(PCA_1, PCA_2, color = sporadic_or_FAP)) +
  geom_point(size = 0.25, alpha = 0.35) +
  labs(title = "Before integration", color = "FAP status") +
  theme_bw(base_size = 10)
p_fap_post <- ggplot(plot_data, aes(Harmony_1, Harmony_2, color = sporadic_or_FAP)) +
  geom_point(size = 0.25, alpha = 0.35) +
  labs(title = "After Harmony", color = "FAP status") +
  theme_bw(base_size = 10)
p_fap <- p_fap_pre + p_fap_post + plot_annotation(
  title = "GSE201348 FAP/sporadic structure before and after integration"
)
ggsave(file.path(figure_dir, "pre_post_integration_by_fap.pdf"), p_fap,
       width = 12, height = 5.5)
ggsave(file.path(figure_dir, "pre_post_integration_by_fap.png"), p_fap,
       width = 12, height = 5.5, dpi = 300)

object_manifest <- data.table(
  object = c("unintegrated_normalized", "harmony_integrated"),
  path = c(unintegrated_file, integrated_file),
  bytes = file.info(c(unintegrated_file, integrated_file))$size,
  cells = c(ncol(unintegrated), ncol(integrated)),
  features = c(nrow(unintegrated), nrow(integrated)),
  raw_counts_retained = TRUE,
  annotation_status = "not_performed"
)
fwrite(object_manifest, file.path(result_dir, "object_manifest.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
writeLines(c(
  paste0("seed=", seed),
  paste0("resume_finished_utc=", format(Sys.time(), tz = "UTC")),
  "resume_reason=integer overflow in diagnostic sampling after completed object writes",
  "cell_annotation=not_performed",
  "differential_expression=not_performed",
  "trajectory_analysis=not_performed"
), file.path(result_dir, "run_provenance.txt"))

message("Stage 5B postprocessing completed. overcorrection_guard=", as.integer(guard_pass))
