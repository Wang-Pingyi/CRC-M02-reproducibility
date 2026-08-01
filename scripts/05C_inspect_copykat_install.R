#!/usr/bin/env Rscript

# Diagnose the pinned CopyKAT package data required for gene annotation.

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(args[[1]], mustWork = TRUE)
library_dir <- file.path(project_root, "environment", "R_library")
.libPaths(c(library_dir, .libPaths()))

suppressPackageStartupMessages(library(copykat))

data_index <- utils::data(package = "copykat")$results
cat("copykat_version=", as.character(packageVersion("copykat")), "\n")
cat("copykat_library=", find.package("copykat"), "\n")
cat("registered_data_objects:\n")
if (is.null(data_index)) {
  cat("NONE\n")
} else {
  print(data_index[, c("Item", "Title"), drop = FALSE])
}
cat(
  "full.anno_in_namespace=",
  exists("full.anno", envir = asNamespace("copykat"), inherits = FALSE),
  "\n"
)
load_env <- new.env(parent = emptyenv())
try(utils::data("full.anno", package = "copykat", envir = load_env))
cat(
  "full.anno_loadable=",
  exists("full.anno", envir = load_env, inherits = FALSE),
  "\n"
)
sysdata_env <- new.env(parent = emptyenv())
try(utils::data("sysdata", package = "copykat", envir = sysdata_env))
cat(
  "sysdata_objects=",
  paste(sort(ls(sysdata_env)), collapse = ";"),
  "\n"
)
cat("installed_data_files:\n")
print(list.files(
  file.path(find.package("copykat"), "data"),
  recursive = TRUE
))
cat("annotate_genes_definition:\n")
print(getAnywhere("annotate.genes"))
cat("annotation_namespace_objects:\n")
print(grep(
  "anno|gene",
  ls(asNamespace("copykat"), all.names = TRUE),
  value = TRUE,
  ignore.case = TRUE
))
cat("annotateGenes.hg20_definition:\n")
print(get("annotateGenes.hg20", envir = asNamespace("copykat")))
namespace_names <- ls(asNamespace("copykat"), all.names = TRUE)
sysdata_users <- vapply(
  namespace_names,
  function(object_name) {
    object <- get(object_name, envir = asNamespace("copykat"))
    is.function(object) && grepl(
      "full\\.anno|DNA\\.hg20|cyclegenes",
      paste(deparse(body(object)), collapse = "\n")
    )
  },
  logical(1)
)
cat(
  "functions_using_sysdata=",
  paste(namespace_names[sysdata_users], collapse = ";"),
  "\n"
)
