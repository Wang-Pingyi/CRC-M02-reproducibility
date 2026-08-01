#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)

required <- c(
  candidate_primary = file.path(
    project_dir, "results", "06A_pseudobulk", "candidate_programs.tsv"
  ),
  candidate_exploratory = file.path(
    project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
  ),
  membership = file.path(
    project_dir, "results", "06A_amendment", "stage_blind_module_membership.tsv"
  ),
  gse132_counts = file.path(
    project_dir, "data_raw", "GSE132465",
    "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
  ),
  gse132_annotation = file.path(
    project_dir, "data_raw", "GSE132465",
    "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
  ),
  manifest = file.path(project_dir, "metadata", "dataset_manifest.tsv"),
  stage6a_report = file.path(
    project_dir, "reports", "stage_6A_pseudobulk_discovery.md"
  ),
  stage6a_amendment = file.path(
    project_dir, "reports", "stage_6A_exploratory_amendment.md"
  )
)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing Stage 7 input(s): ", paste(missing, collapse = "; "))

primary <- utils::read.delim(
  required[["candidate_primary"]], check.names = FALSE, stringsAsFactors = FALSE
)
if (nrow(primary) != 0L) {
  stop("Primary Stage 6A candidate table is no longer empty; evidence hierarchy changed")
}
locked <- read_locked_stage7_modules(project_dir)

gse161_dir <- file.path(project_dir, "data_raw", "GSE161277")
matrix_files <- list.files(gse161_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
feature_files <- list.files(gse161_dir, pattern = "_features\\.tsv\\.gz$", full.names = TRUE)
barcode_files <- list.files(gse161_dir, pattern = "_barcodes\\.tsv\\.gz$", full.names = TRUE)
if (length(matrix_files) != 13L ||
    length(feature_files) != 13L ||
    length(barcode_files) != 13L) {
  stop("GSE161277 must contain 13 complete 10x triplets")
}

all_raw <- c(matrix_files, feature_files, barcode_files, required[c(
  "gse132_counts", "gse132_annotation"
)])
permissions <- file.info(all_raw)$mode
raw_audit <- data.frame(
  file = normalizePath(all_raw, winslash = "/", mustWork = TRUE),
  bytes = file.info(all_raw)$size,
  mode = as.character(permissions),
  owner_writable = bitwAnd(as.integer(permissions), 128L) != 0L,
  stringsAsFactors = FALSE
)
if (any(raw_audit$bytes <= 0)) stop("An input raw-data file is empty")

write_stage7_tsv(
  locked$candidates,
  file.path(paths$result, "preflight", "locked_exploratory_modules.tsv")
)
write_stage7_tsv(
  locked$membership,
  file.path(paths$result, "preflight", "locked_module_membership.tsv")
)
write_stage7_tsv(
  raw_audit,
  file.path(paths$result, "preflight", "raw_input_audit.tsv")
)
writeLines(
  c(
    "primary_candidate_rows=0",
    paste0("locked_exploratory_modules=", nrow(locked$candidates)),
    paste0("locked_membership_rows=", nrow(locked$membership)),
    "validation_reselection=forbidden",
    "cohort_merge=forbidden"
  ),
  file.path(paths$result, "preflight", "evidence_hierarchy.txt")
)
cat(
  "STAGE_7_PREFLIGHT_OK\tmodules=", nrow(locked$candidates),
  "\tgse161_triplets=", length(matrix_files), "\n", sep = ""
)

