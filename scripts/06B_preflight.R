#!/usr/bin/env Rscript

# Analysis: Stage 6B trajectory, regulation and communication preflight
# Date: 2026-07-28
# Random seed: 20260728

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
if (dir.exists(private_library)) .libPaths(c(private_library, .libPaths()))

required_files <- c(
  major_object = file.path(project_dir, "objects", "GSE201348_5C_annotated_final.rds"),
  epithelial_object = file.path(
    project_dir, "objects", "GSE201348_5C_epithelial_annotated_CNV.rds"
  ),
  candidates = file.path(
    project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
  ),
  membership = file.path(
    project_dir, "results", "06A_amendment", "stage_blind_module_membership.tsv"
  )
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required Stage 6B input(s): ", paste(missing_files, collapse = "; "))
}

packages <- c(
  "Seurat", "SeuratObject", "SingleCellExperiment", "slingshot", "tradeSeq",
  "Matrix", "edgeR", "limma", "decoupleR", "dorothea", "msigdbr",
  "fgsea", "liana", "OmnipathR", "igraph", "ggplot2"
)
package_audit <- data.frame(
  package = packages,
  installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(
    packages,
    function(x) {
      if (!requireNamespace(x, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(x))
    },
    character(1)
  ),
  library_path = vapply(
    packages,
    function(x) {
      if (!requireNamespace(x, quietly = TRUE)) return(NA_character_)
      normalizePath(find.package(x), winslash = "/", mustWork = TRUE)
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

output_dir <- file.path(project_dir, "results", "06B_regulatory_inference", "preflight")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name) {
  utils::write.table(
    x, file.path(output_dir, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
}
write_tsv(package_audit, "package_audit.tsv")

candidates <- utils::read.delim(required_files[["candidates"]], check.names = FALSE)
membership <- utils::read.delim(required_files[["membership"]], check.names = FALSE)

required_candidate_columns <- c(
  "module_id", "epithelial_state", "passes_LODO", "exploratory_candidate"
)
if (!all(required_candidate_columns %in% colnames(candidates))) {
  stop("Candidate table is missing required columns")
}
candidates <- candidates[
  candidates$passes_LODO & candidates$exploratory_candidate,
  ,
  drop = FALSE
]
if (nrow(candidates) < 1L || anyDuplicated(candidates$module_id)) {
  stop("Locked candidate module set is empty or duplicated")
}
candidate_membership <- membership[membership$module_id %in% candidates$module_id, ]
if (!setequal(unique(candidate_membership$module_id), candidates$module_id)) {
  stop("Not every locked candidate module has a membership definition")
}

find_column <- function(metadata, choices) {
  hit <- choices[choices %in% colnames(metadata)]
  if (length(hit)) hit[[1L]] else NA_character_
}

cat("Reading epithelial object for sequential preflight\n")
epithelial <- readRDS(required_files[["epithelial_object"]])
epithelial_meta <- epithelial[[]]
epithelial_field_values <- c(
  find_column(epithelial_meta, c("donor_id", "patient_id", "donor")),
  find_column(epithelial_meta, c("sample_id", "biological_tissue_id", "orig.ident")),
  find_column(epithelial_meta, c("lesion_stage", "stage", "condition")),
  find_column(epithelial_meta, c("epithelial_state", "epithelial_annotation")),
  find_column(epithelial_meta, c("epithelial_cluster", "seurat_clusters")),
  find_column(epithelial_meta, c("major_cell_type", "major_annotation"))
)
epithelial_audit <- data.frame(
  object = "epithelial",
  cells = ncol(epithelial),
  genes = nrow(epithelial),
  assays = paste(Assays(epithelial), collapse = ";"),
  reductions = paste(Reductions(epithelial), collapse = ";"),
  stringsAsFactors = FALSE
)
rm(epithelial, epithelial_meta)
invisible(gc())

cat("Reading major object for sequential preflight\n")
major <- readRDS(required_files[["major_object"]])
major_meta <- major[[]]
major_field_values <- c(
  find_column(major_meta, c("donor_id", "patient_id", "donor")),
  find_column(major_meta, c("sample_id", "biological_tissue_id", "orig.ident")),
  find_column(major_meta, c("lesion_stage", "stage", "condition")),
  find_column(major_meta, c("epithelial_state", "epithelial_annotation")),
  find_column(major_meta, c("epithelial_cluster", "seurat_clusters")),
  find_column(major_meta, c("major_cell_type", "major_annotation"))
)
major_audit <- data.frame(
  object = "major",
  cells = ncol(major),
  genes = nrow(major),
  assays = paste(Assays(major), collapse = ";"),
  reductions = paste(Reductions(major), collapse = ";"),
  stringsAsFactors = FALSE
)
rm(major, major_meta)
invisible(gc())

field_audit <- data.frame(
  concept = c(
    "donor", "sample", "stage", "epithelial_state", "epithelial_cluster",
    "major_cell_type"
  ),
  epithelial_column = epithelial_field_values,
  major_object_column = major_field_values,
  stringsAsFactors = FALSE
)

required_concepts <- c("donor", "stage", "epithelial_state", "major_cell_type")
if (any(is.na(field_audit$epithelial_column[field_audit$concept %in% required_concepts[1:3]]))) {
  stop("Epithelial object lacks donor, stage or epithelial-state metadata")
}
if (any(is.na(field_audit$major_object_column[field_audit$concept %in% required_concepts[c(1, 2, 4)]]))) {
  stop("Major object lacks donor, stage or major-cell-type metadata")
}

object_audit <- rbind(epithelial_audit, major_audit)
write_tsv(field_audit, "metadata_field_audit.tsv")
write_tsv(object_audit, "object_audit.tsv")
write_tsv(candidates, "locked_candidate_modules.tsv")
write_tsv(candidate_membership, "locked_candidate_module_membership.tsv")

cat("STAGE_6B_PREFLIGHT_FILES_OK\n")
print(package_audit, row.names = FALSE)
print(field_audit, row.names = FALSE)
print(object_audit, row.names = FALSE)
cat(
  "LOCKED_MODULES\t", nrow(candidates),
  "\tMEMBERSHIP_ROWS\t", nrow(candidate_membership), "\n",
  sep = ""
)
