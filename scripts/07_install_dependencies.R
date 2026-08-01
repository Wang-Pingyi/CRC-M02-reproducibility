#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
shared_library <- file.path(project_dir, "environment", "R-library")
stage_library <- file.path(project_dir, "environment", "R", "7-library")
dir.create(stage_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(stage_library, shared_library, .libPaths()))

required <- c(
  "Matrix", "data.table", "Seurat", "SeuratObject",
  "SingleCellExperiment", "scDblFinder", "edgeR", "limma",
  "ggplot2", "patchwork"
)
audit <- data.frame(
  package = required,
  installed = vapply(required, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(required, function(x) {
    if (!requireNamespace(x, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(x))
  }, character(1)),
  library_path = vapply(required, function(x) {
    if (!requireNamespace(x, quietly = TRUE)) return(NA_character_)
    normalizePath(find.package(x), winslash = "/", mustWork = TRUE)
  }, character(1)),
  stringsAsFactors = FALSE
)

out <- file.path(project_dir, "results", "07_singlecell_replication", "preflight")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  audit, file.path(out, "package_audit.tsv"), sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
if (any(!audit$installed)) {
  stop(
    "Missing required Stage 7 package(s); no unpinned network installation attempted: ",
    paste(audit$package[!audit$installed], collapse = "; ")
  )
}
cat("STAGE_7_DEPENDENCIES_OK\n")
print(audit, row.names = FALSE)

