#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(Seurat)
  library(harmony)
  library(RANN)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else getwd()
raw_dir <- file.path(project_dir, "data_raw", "GSE201348")
config_file <- file.path(project_dir, "config", "qc_parameters.tsv")
metadata_file <- file.path(project_dir, "metadata", "dataset_manifest.tsv")
result_dir <- file.path(project_dir, "results", "05B_full_qc_integration")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_dir, "figures", "05B_full_qc_integration")
object_dir <- file.path(project_dir, "objects")
log_dir <- file.path(project_dir, "logs", "05B_full_qc_integration")
invisible(lapply(c(result_dir, source_dir, figure_dir, object_dir, log_dir),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

params <- fread(config_file, na.strings = c("NA", ""))
get_global <- function(name, numeric = FALSE) {
  value <- params[scope == "global" & parameter == name, value]
  if (length(value) != 1L) stop("Expected one global parameter: ", name)
  if (numeric) as.numeric(value) else value
}
seed <- as.integer(get_global("random_seed", TRUE))
set.seed(seed)

manifest <- fread(metadata_file, na.strings = c("NA", ""))
manifest <- manifest[
  accession == "GSE201348" &
    inclusion %chin% c("include", "include_for_technical_replicate_resolution")
]
if (nrow(manifest) != 72L || uniqueN(manifest$sample_id) != 72L) {
  stop("Expected 72 unique included GSE201348 libraries")
}
if (uniqueN(manifest$biological_sample_id) != 70L) {
  stop("Expected 70 biological tissues after technical-replicate mapping")
}
setorder(manifest, sample_id)

manifest[, lesion_stage := fcase(
  condition %chin% c("normal_mucosa", "unaffected_mucosa"), "normal",
  condition == "polyp_or_adenoma", "adenoma_polyp",
  condition == "cancer", "cancer",
  default = "other"
)]
manifest[, technical_replicate_group := fifelse(
  !is.na(technical_replicate_of), technical_replicate_of, biological_sample_id
)]
manifest[, pathology_qc_flag := fifelse(
  sample_id == "GSM6061645", "discordant_microscopic_fields", "none"
)]
manifest[, qc_risk := fcase(
  sample_id == "GSM6061701", "low_depth_high_risk",
  sample_id %chin% c("GSM6061674", "GSM6061675", "GSM6061708"),
  "excluded_too_few_primary_qc_cells",
  default = "standard"
)]

fwrite(
  manifest[, .(
    sample_id, donor_id, biological_sample_id, sample_role,
    technical_replicate_group, condition, lesion_stage, histology,
    colon_or_rectum, tumor_location, sporadic_or_FAP, platform,
    pathology_qc_flag, qc_risk
  )],
  file.path(source_dir, "included_library_metadata.tsv"),
  sep = "\t", na = "NA"
)

locate_one <- function(sample_id, suffix) {
  hits <- list.files(raw_dir, pattern = paste0("^", sample_id, ".*_", suffix, "$"),
                     full.names = TRUE)
  if (length(hits) != 1L) {
    stop("Expected one ", suffix, " for ", sample_id, "; found ", length(hits))
  }
  hits[[1]]
}

read_10x_sample <- function(sample_id) {
  matrix_file <- locate_one(sample_id, "matrix\\.mtx\\.gz")
  feature_file <- locate_one(sample_id, "features\\.tsv\\.gz")
  barcode_file <- locate_one(sample_id, "barcodes\\.tsv\\.gz")
  features <- fread(feature_file, header = FALSE, data.table = FALSE)
  barcodes <- fread(barcode_file, header = FALSE, data.table = FALSE)
  counts <- as(readMM(gzfile(matrix_file)), "CsparseMatrix")
  if (nrow(counts) != nrow(features) || ncol(counts) != nrow(barcodes)) {
    stop("Dimension mismatch for ", sample_id)
  }
  rownames(counts) <- make.unique(as.character(features[[min(2L, ncol(features))]]))
  colnames(counts) <- paste(sample_id, barcodes[[1]], sep = "_")
  counts
}

robust_bounds <- function(x, lower_k, upper_k, nonnegative = TRUE) {
  med <- median(x, na.rm = TRUE)
  scale <- mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) scale <- IQR(x, na.rm = TRUE) / 1.349
  if (!is.finite(scale) || scale == 0) scale <- max(1, med * 0.1)
  lower <- med - lower_k * scale
  upper <- med + upper_k * scale
  if (nonnegative) lower <- max(0, lower)
  c(lower = lower, upper = upper, median = med, robust_sd = scale)
}

count_lower_k <- get_global("count_lower_mad", TRUE)
count_upper_k <- get_global("count_upper_mad", TRUE)
feature_lower_k <- get_global("feature_lower_mad", TRUE)
feature_upper_k <- get_global("feature_upper_mad", TRUE)
mt_upper_k <- get_global("mitochondrial_upper_mad", TRUE)
ribo_upper_k <- get_global("ribosomal_upper_mad", TRUE)
primary_min_feature <- get_global("primary_min_feature", TRUE)
primary_max_feature <- get_global("primary_max_feature", TRUE)
primary_max_count <- get_global("primary_max_count", TRUE)
primary_max_percent_mt <- get_global("primary_max_percent_mt", TRUE)
doublet_rate_scale <- get_global("doublet_rate_scale", TRUE)
doublet_rate_reference_cells <- get_global("doublet_rate_reference_cells", TRUE)
min_cells_per_library <- as.integer(get_global("min_cells_per_library", TRUE))

counts_list <- vector("list", nrow(manifest))
cell_metadata_list <- vector("list", nrow(manifest))
summary_list <- vector("list", nrow(manifest))
threshold_list <- vector("list", nrow(manifest))
reason_list <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  info <- manifest[i]
  sample_id <- info$sample_id
  message(sprintf("[%02d/%02d] QC %s", i, nrow(manifest), sample_id))
  started <- Sys.time()
  counts <- read_10x_sample(sample_id)
  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt_idx <- grepl("^MT-", rownames(counts), ignore.case = TRUE)
  ribo_idx <- grepl("^RP[SL][0-9A-Z-]+$", rownames(counts), ignore.case = TRUE)
  percent_mt <- if (any(mt_idx)) {
    100 * Matrix::colSums(counts[mt_idx, , drop = FALSE]) / pmax(n_count, 1)
  } else rep(0, length(n_count))
  percent_ribo <- if (any(ribo_idx)) {
    100 * Matrix::colSums(counts[ribo_idx, , drop = FALSE]) / pmax(n_count, 1)
  } else rep(0, length(n_count))

  count_b <- robust_bounds(log10(n_count + 1), count_lower_k, count_upper_k)
  feature_b <- robust_bounds(log10(n_feature + 1), feature_lower_k, feature_upper_k)
  mt_b <- robust_bounds(percent_mt, Inf, mt_upper_k)
  ribo_b <- robust_bounds(percent_ribo, Inf, ribo_upper_k)

  adaptive_low_count <- log10(n_count + 1) < count_b[["lower"]]
  adaptive_low_feature <- log10(n_feature + 1) < feature_b[["lower"]]
  adaptive_high_count <- log10(n_count + 1) > count_b[["upper"]]
  adaptive_high_feature <- log10(n_feature + 1) > feature_b[["upper"]]
  adaptive_high_mt <- percent_mt > mt_b[["upper"]]
  high_ribo <- percent_ribo > ribo_b[["upper"]]
  published_low_feature <- n_feature <= primary_min_feature
  published_high_feature <- n_feature >= primary_max_feature
  published_high_count <- n_count >= primary_max_count
  published_high_mt <- percent_mt >= primary_max_percent_mt
  low_quality <- published_low_feature
  high_complexity <- published_high_feature | published_high_count
  high_mt <- published_high_mt
  doublet_input <- !(low_quality | high_complexity | high_mt)
  library_below_min_cells <- sum(doublet_input) < min_cells_per_library
  library_gate <- rep(library_below_min_cells, ncol(counts))
  if (library_below_min_cells) {
    message(sample_id, " excluded at library gate: ", sum(doublet_input),
            " primary-QC cells")
  }

  expected_doublet_rate <- if (library_below_min_cells) NA_real_ else min(
    doublet_rate_scale * sum(doublet_input) / doublet_rate_reference_cells,
    get_global("max_expected_doublet_rate", TRUE)
  )
  doublet <- rep(FALSE, ncol(counts))
  doublet_score <- rep(NA_real_, ncol(counts))
  doublet_class <- rep("not_evaluated_primary_qc", ncol(counts))
  sce <- NULL
  if (!library_below_min_cells) {
    set.seed(seed + i)
    sce <- SingleCellExperiment(list(counts = counts[, doublet_input, drop = FALSE]))
    colData(sce)$sample_id <- sample_id
    sce <- scDblFinder(
      sce,
      samples = "sample_id",
      dbr = expected_doublet_rate,
      verbose = FALSE
    )
    doublet[doublet_input] <- colData(sce)$scDblFinder.class == "doublet"
    doublet_score[doublet_input] <- colData(sce)$scDblFinder.score
    doublet_class[doublet_input] <- as.character(colData(sce)$scDblFinder.class)
  }

  exclusion_reason <- ifelse(
    low_quality, "low_feature_published_gate",
    ifelse(high_complexity, "high_complexity_published_gate",
           ifelse(high_mt, "high_mitochondrial_published_gate",
                  ifelse(library_gate, "library_below_min_cells",
                         ifelse(doublet, "doublet", "retained"))))
  )
  retained <- exclusion_reason == "retained"
  expected_retained <- if (library_below_min_cells) {
    0L
  } else {
    sum(doublet_input & !doublet)
  }
  if (length(exclusion_reason) != ncol(counts)) {
    stop(sample_id, ": exclusion_reason length does not match cell count")
  }
  if (!library_below_min_cells &&
      sum(exclusion_reason == "doublet") != sum(doublet)) {
    stop(sample_id, ": identified and excluded doublet counts differ")
  }
  if (sum(retained) != expected_retained) {
    stop(sample_id, ": retained-cell reconciliation failed")
  }

  metrics <- data.table(
    sample_id = sample_id,
    cell_id = colnames(counts),
    nCount = as.numeric(n_count),
    nFeature = as.numeric(n_feature),
    percent_mt = as.numeric(percent_mt),
    percent_ribo = as.numeric(percent_ribo),
    scDblFinder_score = doublet_score,
    scDblFinder_class = doublet_class,
    flag_adaptive_low_count = adaptive_low_count,
    flag_adaptive_low_feature = adaptive_low_feature,
    flag_adaptive_high_count = adaptive_high_count,
    flag_adaptive_high_feature = adaptive_high_feature,
    flag_adaptive_high_mt = adaptive_high_mt,
    flag_published_low_feature = published_low_feature,
    flag_published_high_feature = published_high_feature,
    flag_published_high_count = published_high_count,
    flag_published_high_mt = published_high_mt,
    flag_library_below_min_cells = library_gate,
    flag_high_ribo = high_ribo,
    exclusion_reason = exclusion_reason,
    retained = retained
  )
  fwrite(metrics, file.path(source_dir, paste0(sample_id, "_cell_qc_metrics.tsv.gz")),
         sep = "\t", na = "NA")

  threshold_list[[i]] <- data.table(
    sample_id = sample_id,
    metric = c("log10_nCount", "log10_nFeature", "percent_mt", "percent_ribo"),
    adaptive_lower = c(count_b[["lower"]], feature_b[["lower"]], NA, NA),
    adaptive_upper = c(count_b[["upper"]], feature_b[["upper"]],
                       mt_b[["upper"]], ribo_b[["upper"]]),
    primary_lower = c(NA, primary_min_feature, NA, NA),
    primary_upper = c(primary_max_count, primary_max_feature,
                      primary_max_percent_mt, NA),
    median = c(count_b[["median"]], feature_b[["median"]], mt_b[["median"]], ribo_b[["median"]]),
    robust_sd = c(count_b[["robust_sd"]], feature_b[["robust_sd"]],
                  mt_b[["robust_sd"]], ribo_b[["robust_sd"]]),
    filtering_role = c("published_upper_gate", "published_two_sided_gate",
                       "published_upper_gate", "diagnostic_only")
  )
  reason_list[[i]] <- metrics[, .N, by = .(sample_id, exclusion_reason)]
  summary_list[[i]] <- data.table(
    sample_id = sample_id,
    donor_id = info$donor_id,
    biological_sample_id = info$biological_sample_id,
    condition = info$condition,
    lesion_stage = info$lesion_stage,
    histology = info$histology,
    sporadic_or_FAP = info$sporadic_or_FAP,
    platform = info$platform,
    qc_risk = info$qc_risk,
    cells_before = ncol(counts),
    cells_entering_doublet_detection = sum(doublet_input),
    library_below_min_cells = library_below_min_cells,
    library_gate_excluded = sum(exclusion_reason == "library_below_min_cells"),
    low_quality_excluded = sum(low_quality),
    high_complexity_excluded = sum(exclusion_reason == "high_complexity_published_gate"),
    high_mitochondrial_excluded = sum(exclusion_reason == "high_mitochondrial_published_gate"),
    high_ribosomal_flagged = sum(high_ribo),
    ribosomal_only_excluded = 0L,
    expected_doublet_rate = expected_doublet_rate,
    doublets_identified = sum(doublet),
    doublet_final_exclusion = sum(exclusion_reason == "doublet"),
    cells_after = sum(retained),
    retained_percent = 100 * mean(retained),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )

  kept_metrics <- metrics[retained]
  kept_counts <- counts[, retained, drop = FALSE]
  cell_metadata_list[[i]] <- data.table(
    cell_id = kept_metrics$cell_id,
    sample_id = sample_id,
    donor_id = info$donor_id,
    biological_sample_id = info$biological_sample_id,
    sample_role = info$sample_role,
    technical_replicate_group = info$technical_replicate_group,
    condition = info$condition,
    lesion_stage = info$lesion_stage,
    histology = info$histology,
    colon_or_rectum = info$colon_or_rectum,
    tumor_location = info$tumor_location,
    sporadic_or_FAP = info$sporadic_or_FAP,
    platform = info$platform,
    pathology_qc_flag = info$pathology_qc_flag,
    qc_risk = info$qc_risk,
    nCount_RNA = kept_metrics$nCount,
    nFeature_RNA = kept_metrics$nFeature,
    percent_mt = kept_metrics$percent_mt,
    percent_ribo = kept_metrics$percent_ribo,
    scDblFinder_score = kept_metrics$scDblFinder_score
  )
  counts_list[[i]] <- kept_counts
  rm(counts, kept_counts, sce, metrics, kept_metrics)
  gc()
}

qc_summary <- rbindlist(summary_list, fill = TRUE)
thresholds <- rbindlist(threshold_list, fill = TRUE)
exclusion_counts <- rbindlist(reason_list, fill = TRUE)
fwrite(qc_summary, file.path(result_dir, "library_qc_summary.tsv"), sep = "\t", na = "NA")
fwrite(thresholds, file.path(source_dir, "sample_adaptive_thresholds.tsv"), sep = "\t", na = "NA")
fwrite(exclusion_counts, file.path(source_dir, "exclusion_reason_counts.tsv"), sep = "\t", na = "NA")

reference_genes <- rownames(counts_list[[1]])
same_gene_order <- vapply(counts_list, function(x) identical(rownames(x), reference_genes), logical(1))
if (!all(same_gene_order)) stop("Gene order differs among GSE201348 matrices")
all_counts <- do.call(cbind, counts_list)
cell_metadata <- rbindlist(cell_metadata_list, use.names = TRUE)
if (!identical(colnames(all_counts), cell_metadata$cell_id)) {
  stop("Cell metadata order does not match merged count columns")
}
cell_metadata_df <- as.data.frame(cell_metadata)
rownames(cell_metadata_df) <- cell_metadata_df$cell_id

message("Creating unintegrated Seurat object with ", ncol(all_counts), " cells")
unintegrated <- CreateSeuratObject(
  counts = all_counts,
  assay = "RNA",
  project = "GSE201348",
  meta.data = cell_metadata_df,
  min.cells = 0,
  min.features = 0
)
unintegrated <- NormalizeData(
  unintegrated,
  normalization.method = "LogNormalize",
  scale.factor = get_global("normalization_scale_factor", TRUE),
  verbose = FALSE
)
unintegrated <- FindVariableFeatures(
  unintegrated,
  selection.method = "vst",
  nfeatures = as.integer(get_global("variable_features", TRUE)),
  verbose = FALSE
)
unintegrated <- ScaleData(
  unintegrated,
  features = VariableFeatures(unintegrated),
  verbose = FALSE
)
n_pcs <- as.integer(get_global("pca_dimensions", TRUE))
unintegrated <- RunPCA(
  unintegrated,
  features = VariableFeatures(unintegrated),
  npcs = n_pcs,
  seed.use = seed,
  verbose = FALSE
)
unintegrated@misc$stage_5B <- list(
  seed = seed,
  qc_config = "config/qc_parameters.tsv",
  raw_counts_retained = TRUE,
  integration_status = "not_integrated"
)
unintegrated_file <- file.path(object_dir, "GSE201348_unintegrated_normalized.rds")
saveRDS(unintegrated, unintegrated_file, compress = get_global("object_compression"))

message("Running conservative Harmony integration")
integrated <- harmony::RunHarmony(
  unintegrated,
  group.by.vars = get_global("harmony_group"),
  reduction.use = "pca",
  dims.use = seq_len(n_pcs),
  reduction.save = "harmony",
  theta = get_global("harmony_theta", TRUE),
  lambda = get_global("harmony_lambda", TRUE),
  project.dim = FALSE,
  verbose = TRUE
)
integrated@misc$stage_5B <- list(
  seed = seed,
  qc_config = "config/qc_parameters.tsv",
  raw_counts_retained = TRUE,
  integration_method = "Harmony",
  harmony_group = get_global("harmony_group"),
  integration_status = "integrated_no_annotation"
)
integrated_file <- file.path(object_dir, "GSE201348_harmony_integrated.rds")
saveRDS(integrated, integrated_file, compress = get_global("object_compression"))

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

cap <- as.integer(get_global("diagnostic_cell_cap", TRUE))
set.seed(seed)
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
  levels <- unique(group)
  between_ss <- rep(0, ncol(embedding))
  for (lev in levels) {
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
overcorrection <- data.table(
  metric = c("stage_eta2_pre", "stage_eta2_post", "stage_eta2_retention",
             "stage_neighbor_excess_post", "overcorrection_guard"),
  value = c(stage_eta_pre, stage_eta_post, stage_eta_retention,
            stage_excess_post,
            as.numeric(
              is.finite(stage_eta_retention) &&
                stage_eta_retention >= get_global("min_stage_eta2_retention", TRUE) &&
                stage_excess_post >= get_global("min_stage_neighbor_excess", TRUE)
            )),
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
  paste0("run_finished_utc=", format(Sys.time(), tz = "UTC")),
  "cell_annotation=not_performed",
  "differential_expression=not_performed",
  "trajectory_analysis=not_performed"
), file.path(result_dir, "run_provenance.txt"))

message("Stage 5B full QC and Harmony integration completed.")
