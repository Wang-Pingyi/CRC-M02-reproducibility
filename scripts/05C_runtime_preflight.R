#!/usr/bin/env Rscript

# Parse Stage 5C scripts and inspect the accepted Stage 5B object without
# modifying it.

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}

script_files <- list.files(
  file.path(project_root, "scripts"),
  pattern = "^05C_.*[.]R$",
  full.names = TRUE
)
if (!length(script_files)) stop("No Stage 5C R scripts found")
for (script_file in script_files) {
  parse(file = script_file)
  message("PARSE_OK ", basename(script_file))
}

parameter_file <- file.path(
  project_root, "config", "annotation_parameters.tsv"
)
evidence_file <- file.path(
  project_root, "metadata", "annotation_evidence.tsv"
)
parameters <- fread(parameter_file)
evidence <- fread(evidence_file)
if (anyDuplicated(parameters[, paste(section, parameter, sep = "/")])) {
  stop("Duplicate annotation parameters")
}
if (anyDuplicated(evidence[, paste(level, annotation, sep = "/")])) {
  stop("Duplicate annotation evidence rows")
}

object_file <- file.path(
  project_root, "objects", "GSE201348_harmony_integrated.rds"
)
obj <- readRDS(object_file)
required_meta <- c(
  "sample_id", "donor_id", "biological_sample_id", "lesion_stage",
  "sporadic_or_FAP"
)
if (!all(required_meta %in% colnames(obj@meta.data))) {
  stop(
    "Missing metadata: ",
    paste(setdiff(required_meta, colnames(obj@meta.data)), collapse = ";")
  )
}
if (!"RNA" %in% names(obj@assays)) stop("RNA assay is missing")
if (!"harmony" %in% names(obj@reductions)) {
  stop("Harmony reduction is missing")
}
if (ncol(Embeddings(obj, "harmony")) < 30L) {
  stop("Harmony has fewer than 30 dimensions")
}

counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
data <- GetAssayData(obj, assay = "RNA", slot = "data")
if (!identical(dim(counts), dim(data))) {
  stop("Counts and normalized-data dimensions differ")
}

result_dir <- file.path(project_root, "results", "05C_annotation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    paste0("R=", getRversion()),
    paste0("Seurat=", packageVersion("Seurat")),
    paste0("SeuratObject=", packageVersion("SeuratObject")),
    paste0("harmony=", packageVersion("harmony")),
    paste0("object_class=", paste(class(obj), collapse = ";")),
    paste0("RNA_assay_class=", paste(class(obj[["RNA"]]), collapse = ";")),
    paste0("features=", nrow(obj)),
    paste0("cells=", ncol(obj)),
    paste0("raw_count_nonzero_entries=", length(counts@x)),
    paste0("reductions=", paste(names(obj@reductions), collapse = ";")),
    paste0(
      "required_metadata_present=",
      paste(required_meta, collapse = ";")
    ),
    paste0("parsed_R_scripts=", length(script_files))
  ),
  file.path(result_dir, "runtime_preflight.txt")
)
message(
  "RUNTIME_PREFLIGHT_OK cells=", ncol(obj),
  " features=", nrow(obj),
  " Seurat=", packageVersion("Seurat")
)
