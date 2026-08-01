#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 08B_preflight.R PROJECT_DIR")
project <- normalizePath(args[1], mustWork = TRUE)

r_files <- setdiff(Sys.glob(file.path(project, "scripts", "08B_*.R")),
                   file.path(project, "scripts", "08B_preflight.R"))
for (f in r_files) {
  parse(file = f)
  cat("PARSE_OK", basename(f), "\n")
}

required <- c("data.table", "edgeR", "sandwich", "lmtest", "survival", "metafor", "ggplot2")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing server R packages: ", paste(missing, collapse = ", "))

source_dir <- file.path(project, "environment", "sources", "08B")
gdc <- read.delim(file.path(source_dir, "gdc_tcga_coad_file_sample_metadata.tsv"), check.names = FALSE)
raw_files <- file.path(project, "data_raw", "TCGA-COAD", gdc$local_file_name)
if (nrow(gdc) != 524L || anyDuplicated(gdc$local_file_name) || any(!file.exists(raw_files))) {
  stop("GDC-to-raw-file reconciliation failed")
}

sample_clin <- read.delim(file.path(source_dir, "cbioportal_coadread_pancan_sample.tsv"),
                          comment.char = "#", check.names = FALSE)
patient_clin <- read.delim(file.path(source_dir, "cbioportal_coadread_pancan_patient.tsv"),
                           comment.char = "#", check.names = FALSE)
if (!all(c("PATIENT_ID", "SAMPLE_ID", "MSI_SCORE_MANTIS") %in% names(sample_clin))) {
  stop("cBioPortal sample clinical columns are incomplete")
}
if (!all(c("PATIENT_ID", "CANCER_TYPE_ACRONYM", "AJCC_PATHOLOGIC_TUMOR_STAGE",
           "OS_STATUS", "OS_MONTHS") %in% names(patient_clin))) {
  stop("cBioPortal patient clinical columns are incomplete")
}

cat("PREFLIGHT_OK GDC_ROWS=524 RAW_MISSING=0 COAD_PATIENT_ROWS=",
    sum(patient_clin$CANCER_TYPE_ACRONYM == "COAD"), "\n", sep = "")
