#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 08A_rebuild_mapping_audits.R <project_dir> <run_id>")
}

project_dir <- normalizePath(args[[1L]], mustWork = TRUE)
run_id <- args[[2L]]
source(file.path(project_dir, "scripts", "08A_helpers.R"))
paths <- stage8_paths(project_dir, run_id)

for (accession in c("GSE41657", "GSE100179", "GSE8671")) {
  cohort_dir <- file.path(paths$result, accession)
  mapping_file <- file.path(cohort_dir, "locked_module_gene_mapping.tsv")
  summary_file <- file.path(cohort_dir, "locked_module_mapping_summary.tsv")
  if (!file.exists(mapping_file) || !file.exists(summary_file)) {
    stop("Existing mapping outputs are incomplete for ", accession)
  }
  mapping <- utils::read.delim(mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
  previous <- utils::read.delim(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
  method <- unique(previous$normalization_method)
  if (length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("Normalization provenance is ambiguous for ", accession)
  }
  stage8_write_mapping_audit(mapping, method, cohort_dir)
}

source(file.path(project_dir, "scripts", "08A_finalize_report.R"), local = new.env(parent = globalenv()))
cat("STAGE8A_MAPPING_AUDITS_REBUILT\n")
