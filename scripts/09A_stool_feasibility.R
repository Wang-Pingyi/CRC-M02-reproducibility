#!/usr/bin/env Rscript

set.seed(20260729)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09A_stool_feasibility.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
lib_dir <- file.path(project_dir, "environment", "Rlib_stage9A")
.libPaths(c(lib_dir, .libPaths()))

required <- c(
  "oligo", "pd.hta.2.0", "hta20transcriptcluster.db",
  "AnnotationDbi", "Biobase", "matrixStats"
)
missing_packages <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

result_dir <- file.path(project_dir, "results", "09A_stool_feasibility", run_id)
work_dir <- file.path(project_dir, "data_processed", "09A_stool_feasibility", run_id)
cel_dir <- file.path(work_dir, "training_cel")
figure_dir <- file.path(project_dir, "figures", "09A_stool_feasibility", run_id)
object_dir <- file.path(project_dir, "objects")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
robust_z <- function(x) {
  center <- stats::median(x, na.rm = TRUE)
  scale <- stats::mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) return(rep(0, length(x)))
  (x - center) / scale
}

inventory <- utils::read.delim(
  file.path(result_dir, "GSE99573_sample_inventory.tsv"),
  check.names = FALSE
)
training_inventory <- inventory[inventory$validation_split == "training", , drop = FALSE]
if (nrow(training_inventory) != 265L) stop("Training inventory is not 265 samples")
if (any(training_inventory$expression_extracted != "FALSE")) {
  stop("Split audit unexpectedly marks expression as extracted")
}

cel_files <- list.files(
  cel_dir, pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE
)
cel_sample_ids <- sub("_.*$", "", basename(cel_files))
if (length(cel_files) != 265L || anyDuplicated(cel_sample_ids)) {
  stop("Expected 265 uniquely identified uncompressed training CEL files")
}
if (!setequal(cel_sample_ids, training_inventory$sample_id)) {
  stop("Extracted CEL sample IDs do not equal the frozen training set")
}
test_ids <- inventory$sample_id[inventory$validation_split == "testing"]
if (any(test_ids %in% cel_sample_ids)) {
  stop("Test-set CEL detected in the training work directory")
}
cel_files <- cel_files[match(training_inventory$sample_id, cel_sample_ids)]

message("Reading 265 training CEL files; test CEL files were not extracted")
raw_data <- oligo::read.celfiles(
  cel_files, pkgname = "pd.hta.2.0", verbose = TRUE
)
message("Running training-only core transcript-cluster RMA")
normalized <- oligo::rma(
  raw_data, target = "core", background = TRUE, normalize = TRUE
)
expression <- Biobase::exprs(normalized)
colnames(expression) <- training_inventory$sample_id
if (ncol(expression) != 265L || any(!is.finite(expression))) {
  stop("Invalid normalized training matrix")
}
rm(raw_data)
invisible(gc())

training_object_path <- file.path(
  object_dir, paste0("GSE99573_9A_training_RMA_", run_id, ".rds")
)
saveRDS(
  list(
    expression = expression,
    sample_id = colnames(expression),
    platform = "GPL17586",
    preprocessing = "oligo::rma target=core, training arrays only",
    test_expression_accessed = FALSE,
    seed = 20260729L
  ),
  training_object_path,
  compress = "xz"
)

candidate_modules <- utils::read.delim(
  file.path(
    project_dir, "results_final",
    "stage_6A_exploratory_candidate_modules.tsv"
  ),
  check.names = FALSE
)
membership_all <- utils::read.delim(
  file.path(
    project_dir, "results_final",
    "stage_6A_stage_blind_module_membership.tsv"
  ),
  check.names = FALSE
)
membership <- membership_all[
  membership_all$module_id %in% candidate_modules$module_id,
  c("module_id", "epithelial_state", "gene"),
  drop = FALSE
]
if (nrow(membership) != 747L ||
    length(unique(membership$gene)) != 632L ||
    length(unique(membership$module_id)) != 6L) {
  stop("Locked candidate membership no longer matches the Stage 6A freeze")
}
write_tsv(membership, file.path(result_dir, "locked_module_membership.tsv"))

annotation <- AnnotationDbi::select(
  hta20transcriptcluster.db::hta20transcriptcluster.db,
  keys = rownames(expression),
  keytype = "PROBEID",
  columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
)
annotation <- unique(annotation)
annotation <- annotation[
  !is.na(annotation$SYMBOL) & nzchar(annotation$SYMBOL),
  ,
  drop = FALSE
]
annotation$PROBEID <- as.character(annotation$PROBEID)
annotation <- annotation[annotation$PROBEID %in% rownames(expression), , drop = FALSE]
write_tsv(annotation, file.path(result_dir, "platform_core_annotation.tsv"))

probe_medians <- matrixStats::rowMedians(expression, na.rm = TRUE)
names(probe_medians) <- rownames(expression)
annotated_probe_ids <- intersect(unique(annotation$PROBEID), rownames(expression))
detectability_floor <- as.numeric(
  stats::quantile(
    probe_medians[annotated_probe_ids],
    probs = 0.20,
    names = FALSE,
    na.rm = TRUE,
    type = 8
  )
)

locked_genes <- sort(unique(membership$gene))
mapping_rows <- vector("list", length(locked_genes))
gene_expression <- matrix(
  NA_real_,
  nrow = length(locked_genes),
  ncol = ncol(expression),
  dimnames = list(locked_genes, colnames(expression))
)
for (i in seq_along(locked_genes)) {
  gene <- locked_genes[[i]]
  probe_ids <- intersect(
    unique(annotation$PROBEID[annotation$SYMBOL == gene]),
    rownames(expression)
  )
  if (length(probe_ids)) {
    values <- expression[probe_ids, , drop = FALSE]
    gene_values <- if (nrow(values) == 1L) {
      as.numeric(values[1L, ])
    } else {
      matrixStats::colMedians(values, na.rm = TRUE)
    }
    gene_expression[i, ] <- gene_values
    median_expression <- stats::median(gene_values, na.rm = TRUE)
    fraction_above_floor <- mean(gene_values > detectability_floor, na.rm = TRUE)
    detectable <- median_expression > detectability_floor &&
      fraction_above_floor >= 0.20
  } else {
    median_expression <- NA_real_
    fraction_above_floor <- NA_real_
    detectable <- FALSE
  }
  mapping_rows[[i]] <- data.frame(
    gene = gene,
    n_core_transcript_clusters = length(probe_ids),
    transcript_cluster_ids = if (length(probe_ids)) {
      paste(probe_ids, collapse = ";")
    } else {
      "NA"
    },
    mapping_status = if (length(probe_ids)) "mapped" else "missing_probe",
    training_median_log2_RMA = median_expression,
    training_fraction_above_floor = fraction_above_floor,
    detectability_floor_log2_RMA = detectability_floor,
    detectability_status = if (!length(probe_ids)) {
      "missing_probe"
    } else if (detectable) {
      "detectable"
    } else {
      "low_expression"
    },
    expression_data_scope = if (length(probe_ids)) "training_only" else "none",
    stringsAsFactors = FALSE
  )
}
mapping <- do.call(rbind, mapping_rows)
write_tsv(mapping, file.path(result_dir, "locked_gene_probe_mapping.tsv"))
write_tsv(
  mapping[
    ,
    c(
      "gene", "training_median_log2_RMA", "training_fraction_above_floor",
      "detectability_floor_log2_RMA", "detectability_status",
      "expression_data_scope"
    )
  ],
  file.path(result_dir, "locked_gene_detectability_training.tsv")
)

module_summary <- do.call(
  rbind,
  lapply(unique(membership$module_id), function(module_id) {
    genes <- unique(membership$gene[membership$module_id == module_id])
    z <- mapping[match(genes, mapping$gene), , drop = FALSE]
    data.frame(
      module_id = module_id,
      locked_genes = length(genes),
      mapped_genes = sum(z$mapping_status == "mapped"),
      detectable_genes = sum(z$detectability_status == "detectable"),
      low_expression_genes = sum(z$detectability_status == "low_expression"),
      missing_probe_genes = sum(z$detectability_status == "missing_probe"),
      mapped_fraction = mean(z$mapping_status == "mapped"),
      detectable_fraction = mean(z$detectability_status == "detectable"),
      stringsAsFactors = FALSE
    )
  })
)
write_tsv(
  module_summary,
  file.path(result_dir, "locked_module_detectability_summary.tsv")
)

row_medians <- matrixStats::rowMedians(expression, na.rm = TRUE)
rle <- sweep(expression, 1L, row_medians, FUN = "-")
sample_median <- matrixStats::colMedians(expression, na.rm = TRUE)
sample_iqr <- apply(expression, 2L, stats::IQR, na.rm = TRUE)
rle_median <- matrixStats::colMedians(rle, na.rm = TRUE)
rle_iqr <- apply(rle, 2L, stats::IQR, na.rm = TRUE)
row_variances <- matrixStats::rowVars(expression, na.rm = TRUE)
top_n <- min(5000L, length(row_variances))
top_rows <- order(row_variances, decreasing = TRUE)[seq_len(top_n)]
qc_matrix <- expression[top_rows, , drop = FALSE]
sample_cor <- stats::cor(qc_matrix, method = "pearson", use = "pairwise.complete.obs")
mean_interarray_cor <- (rowSums(sample_cor) - 1) / (ncol(sample_cor) - 1)
pca <- stats::prcomp(t(qc_matrix), center = TRUE, scale. = FALSE, rank. = 5)
pc_scores <- pca$x[, seq_len(min(5L, ncol(pca$x))), drop = FALSE]
pc_z <- apply(pc_scores, 2L, robust_z)
if (is.null(dim(pc_z))) pc_z <- matrix(pc_z, ncol = 1L)
robust_pc_distance <- sqrt(rowSums(pc_z^2))
pc_flag_threshold <- sqrt(stats::qchisq(0.999, df = ncol(pc_z)))

qc <- data.frame(
  sample_id = colnames(expression),
  array_median = sample_median,
  array_IQR = sample_iqr,
  RLE_median = rle_median,
  RLE_IQR = rle_iqr,
  mean_interarray_correlation = mean_interarray_cor[colnames(expression)],
  robust_PC_distance = robust_pc_distance,
  stringsAsFactors = FALSE
)
qc$array_median_robust_z <- robust_z(qc$array_median)
qc$array_IQR_robust_z <- robust_z(qc$array_IQR)
qc$RLE_median_robust_z <- robust_z(qc$RLE_median)
qc$RLE_IQR_robust_z <- robust_z(qc$RLE_IQR)
qc$correlation_robust_z <- robust_z(qc$mean_interarray_correlation)
qc$distribution_flag <- abs(qc$array_median_robust_z) > 5 |
  abs(qc$array_IQR_robust_z) > 5 |
  abs(qc$RLE_median_robust_z) > 5 |
  abs(qc$RLE_IQR_robust_z) > 5
qc$correlation_flag <- qc$correlation_robust_z < -5
qc$PCA_flag <- qc$robust_PC_distance > pc_flag_threshold
qc$qc_review_flag <- qc$distribution_flag | qc$correlation_flag | qc$PCA_flag
qc$automatic_exclusion <- FALSE
qc$outcome_labels_used <- FALSE
write_tsv(qc, file.path(result_dir, "training_sample_qc_metrics.tsv"))
write_tsv(
  qc[qc$qc_review_flag, , drop = FALSE],
  file.path(result_dir, "training_outlier_flags.tsv")
)

grDevices::pdf(
  file.path(figure_dir, "training_array_QC_outcome_blind.pdf"),
  width = 11, height = 8.5
)
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
graphics::boxplot(
  expression, outline = FALSE, las = 2, xaxt = "n",
  main = "Training-only RMA distributions", ylab = "log2 RMA"
)
graphics::boxplot(
  rle, outline = FALSE, las = 2, xaxt = "n",
  main = "Training-only relative log expression", ylab = "RLE"
)
graphics::plot(
  pca$x[, 1], pca$x[, 2], pch = 19,
  col = ifelse(qc$qc_review_flag, "#D55E00", "#0072B2"),
  xlab = "PC1", ylab = "PC2", main = "Outcome-blind training PCA"
)
graphics::abline(h = 0, v = 0, col = "grey85")
graphics::hist(
  qc$mean_interarray_correlation, breaks = 30,
  main = "Mean inter-array correlation", xlab = "Correlation"
)
grDevices::dev.off()

software <- data.frame(
  package = required,
  version = vapply(
    required, function(x) as.character(utils::packageVersion(x)), character(1)
  ),
  stringsAsFactors = FALSE
)
software <- rbind(
  data.frame(
    package = "R",
    version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  ),
  software
)
write_tsv(software, file.path(result_dir, "software_versions.tsv"))

checks <- data.frame(
  check = c(
    "training_samples_265",
    "test_samples_absent_from_cel_directory",
    "normalized_matrix_finite",
    "six_locked_modules",
    "locked_membership_747_rows",
    "locked_unique_genes_632",
    "all_locked_genes_reported",
    "no_automatic_sample_exclusions",
    "outcome_labels_not_used",
    "test_expression_not_accessed"
  ),
  passed = c(
    ncol(expression) == 265L,
    !any(test_ids %in% cel_sample_ids),
    all(is.finite(expression)),
    length(unique(membership$module_id)) == 6L,
    nrow(membership) == 747L,
    length(unique(membership$gene)) == 632L,
    nrow(mapping) == 632L,
    !any(qc$automatic_exclusion),
    !any(qc$outcome_labels_used),
    TRUE
  ),
  observed = c(
    ncol(expression),
    sum(test_ids %in% cel_sample_ids),
    sum(is.finite(expression)),
    length(unique(membership$module_id)),
    nrow(membership),
    length(unique(membership$gene)),
    nrow(mapping),
    sum(qc$automatic_exclusion),
    sum(qc$outcome_labels_used),
    "test CEL inventory only"
  ),
  stringsAsFactors = FALSE
)
write_tsv(checks, file.path(result_dir, "analysis_validation_checks.tsv"))
if (!all(checks$passed)) stop("At least one Stage 9A analysis validation failed")

message("STAGE9A_FEASIBILITY_ANALYSIS_OK run_id=", run_id)
