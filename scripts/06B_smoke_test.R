#!/usr/bin/env Rscript

# Technical smoke test: Stage 6B APIs and locked inputs
# Date: 2026-07-28

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
if (dir.exists(private_library)) .libPaths(c(private_library, .libPaths()))

packages <- c(
  "slingshot", "tradeSeq", "decoupleR", "dorothea", "liana",
  "msigdbr", "edgeR", "limma"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Smoke-test packages missing: ", paste(missing, collapse = ", "))

required_functions <- list(
  slingshot = c(
    "slingshot", "slingPseudotime", "slingCurveWeights", "slingLineages"
  ),
  tradeSeq = c("fitGAM", "associationTest"),
  decoupleR = c("run_ulm"),
  liana = c("select_resource")
)
for (package in names(required_functions)) {
  namespace <- asNamespace(package)
  absent <- required_functions[[package]][
    !vapply(required_functions[[package]], exists, logical(1), envir = namespace)
  ]
  if (length(absent)) {
    stop("Missing API in ", package, ": ", paste(absent, collapse = ", "))
  }
}

candidate_path <- file.path(
  project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
)
membership_path <- file.path(
  project_dir, "results", "06A_amendment", "stage_blind_module_membership.tsv"
)
candidates <- utils::read.delim(candidate_path, check.names = FALSE)
candidates <- candidates[candidates$exploratory_candidate & candidates$passes_LODO, ]
membership <- utils::read.delim(membership_path, check.names = FALSE)
membership <- membership[membership$module_id %in% candidates$module_id, ]
if (nrow(candidates) != 6L || anyDuplicated(candidates$module_id)) {
  stop("Locked module set is not exactly six unique modules")
}
if (!setequal(unique(membership$module_id), candidates$module_id)) {
  stop("Locked module membership is incomplete")
}

data("dorothea_hs", package = "dorothea", envir = environment())
network <- get("dorothea_hs", envir = environment())
if (!all(c("tf", "target", "mor", "confidence") %in% colnames(network))) {
  stop("DoRothEA schema is incompatible")
}
resource <- liana::select_resource("Consensus")
if (is.list(resource) && !is.data.frame(resource) && length(resource) == 1L) {
  resource <- resource[[1L]]
}
ligand_ok <- any(c("ligand", "source_genesymbol", "source", "ligand_complex") %in%
  colnames(resource))
receptor_ok <- any(c("receptor", "target_genesymbol", "target", "receptor_complex") %in%
  colnames(resource))
if (!ligand_ok || !receptor_ok || !nrow(resource)) {
  stop("LIANA consensus schema is incompatible")
}

cat(
  "STAGE_6B_SMOKE_TEST_OK",
  paste0("modules=", nrow(candidates)),
  paste0("module_genes=", length(unique(membership$gene))),
  paste0("dorothea_edges=", nrow(network)),
  paste0("liana_interactions=", nrow(resource)),
  sep = "\t"
)
cat("\n")
