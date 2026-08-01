#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
shared_library <- file.path(project_dir, "environment", "R-library")
.libPaths(c(shared_library, .libPaths()))
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

gse161_dir <- file.path(project_dir, "data_raw", "GSE161277")
matrix_file <- sort(list.files(
  gse161_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE
))[[1L]]
prefix <- sub("_matrix\\.mtx\\.gz$", "", matrix_file)
pilot_counts <- Seurat::ReadMtx(
  mtx = matrix_file,
  features = paste0(prefix, "_features.tsv.gz"),
  cells = paste0(prefix, "_barcodes.tsv.gz"),
  feature.column = 2L,
  unique.features = TRUE
)
if (nrow(pilot_counts) < 30000L || ncol(pilot_counts) < 100L) {
  stop("GSE161277 real-data pilot has implausible dimensions")
}
pilot_idx <- seq_len(min(500L, ncol(pilot_counts)))
pilot <- pilot_counts[, pilot_idx, drop = FALSE]
pilot_metrics <- data.frame(
  cells = ncol(pilot),
  genes = nrow(pilot),
  median_nCount = stats::median(Matrix::colSums(pilot)),
  median_nFeature = stats::median(Matrix::colSums(pilot > 0)),
  median_percent_mt = stats::median(percent_feature_set(pilot, "^MT-")),
  stringsAsFactors = FALSE
)
rm(pilot_counts, pilot)
invisible(gc())

annotation_file <- file.path(
  project_dir, "data_raw", "GSE132465",
  "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
)
annotation <- utils::read.delim(
  gzfile(annotation_file), nrows = 100L, check.names = FALSE,
  stringsAsFactors = FALSE
)
required_annotation <- c("Index", "Patient", "Class", "Sample", "Cell_type")
if (!all(required_annotation %in% colnames(annotation))) {
  stop("GSE132465 annotation pilot lacks required columns")
}

count_file <- file.path(
  project_dir, "data_raw", "GSE132465",
  "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
)
con <- gzfile(count_file, open = "rt")
header <- readLines(con, n = 1L, warn = FALSE)
first_row <- readLines(con, n = 1L, warn = FALSE)
close(con)
header_fields <- length(strsplit(header, "\t", fixed = TRUE)[[1L]])
row_fields <- length(strsplit(first_row, "\t", fixed = TRUE)[[1L]])
if (header_fields != 63690L || row_fields != header_fields) {
  stop("GSE132465 matrix pilot found unexpected column count")
}

pilot_audit <- rbind(
  data.frame(
    cohort = "GSE161277", test = "10x_real_matrix",
    result = paste0(pilot_metrics$cells, "_cells_piloted"), passed = TRUE
  ),
  data.frame(
    cohort = "GSE132465", test = "annotation_schema",
    result = paste(required_annotation, collapse = ";"), passed = TRUE
  ),
  data.frame(
    cohort = "GSE132465", test = "matrix_column_count",
    result = as.character(header_fields - 1L), passed = TRUE
  )
)
write_stage7_tsv(
  pilot_metrics,
  file.path(paths$result, "preflight", "GSE161277_real_data_pilot_metrics.tsv")
)
write_stage7_tsv(
  pilot_audit,
  file.path(paths$result, "preflight", "real_data_pilot_audit.tsv")
)
cat("STAGE_7_REAL_DATA_PILOT_OK\n")

