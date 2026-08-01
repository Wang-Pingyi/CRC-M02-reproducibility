#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else getwd()
raw_dir <- file.path(project_dir, "data_raw", "GSE201348")
metadata_file <- file.path(project_dir, "metadata", "dataset_manifest.tsv")
config_file <- file.path(project_dir, "config", "qc_parameters.tsv")
output_file <- file.path(project_dir, "results", "05B_qc_threshold_preflight.tsv")

params <- fread(config_file, na.strings = c("NA", ""))
get_global <- function(name) {
  value <- params[scope == "global" & parameter == name, value]
  if (length(value) != 1L) stop("Expected one global parameter: ", name)
  as.numeric(value)
}

manifest <- fread(metadata_file, na.strings = c("NA", ""))
manifest <- manifest[
  accession == "GSE201348" &
    inclusion %chin% c("include", "include_for_technical_replicate_resolution")
]
setorder(manifest, sample_id)
stopifnot(nrow(manifest) == 72L, uniqueN(manifest$sample_id) == 72L)

locate_one <- function(sample_id, suffix) {
  hits <- list.files(raw_dir, pattern = paste0("^", sample_id, ".*_", suffix, "$"),
                     full.names = TRUE)
  if (length(hits) != 1L) stop("Expected one ", suffix, " for ", sample_id)
  hits[[1]]
}

read_metrics <- function(sample_id) {
  matrix_file <- locate_one(sample_id, "matrix\\.mtx\\.gz")
  feature_file <- locate_one(sample_id, "features\\.tsv\\.gz")
  features <- fread(feature_file, header = FALSE, data.table = FALSE)
  counts <- as(readMM(gzfile(matrix_file)), "CsparseMatrix")
  rownames(counts) <- make.unique(as.character(features[[min(2L, ncol(features))]]))
  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt_idx <- grepl("^MT-", rownames(counts), ignore.case = TRUE)
  percent_mt <- if (any(mt_idx)) {
    100 * Matrix::colSums(counts[mt_idx, , drop = FALSE]) / pmax(n_count, 1)
  } else rep(0, length(n_count))
  data.table(n_count, n_feature, percent_mt)
}

min_feature <- get_global("primary_min_feature")
max_feature <- get_global("primary_max_feature")
max_count <- get_global("primary_max_count")
max_mt <- get_global("primary_max_percent_mt")

result <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  info <- manifest[i]
  message(sprintf("[%02d/%02d] %s", i, nrow(manifest), info$sample_id))
  x <- read_metrics(info$sample_id)
  pass <- x$n_feature > min_feature &
    x$n_feature < max_feature &
    x$n_count < max_count &
    x$percent_mt < max_mt
  data.table(
    sample_id = info$sample_id,
    donor_id = info$donor_id,
    biological_sample_id = info$biological_sample_id,
    condition = info$condition,
    sporadic_or_FAP = info$sporadic_or_FAP,
    cells_input = nrow(x),
    cells_primary_gate = sum(pass),
    retention_before_doublets_percent = 100 * mean(pass),
    low_feature = sum(x$n_feature <= min_feature),
    high_feature = sum(x$n_feature >= max_feature),
    high_count = sum(x$n_count >= max_count),
    high_mt = sum(x$percent_mt >= max_mt)
  )
}))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(result, output_file, sep = "\t", na = "NA")
message("Wrote ", output_file)
message("Primary-gate cells: ", sum(result$cells_primary_gate),
        " / ", sum(result$cells_input))
message("Libraries below 100 cells: ", sum(result$cells_primary_gate < 100L))
