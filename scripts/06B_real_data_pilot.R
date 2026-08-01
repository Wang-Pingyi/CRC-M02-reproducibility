#!/usr/bin/env Rscript

# Technical pilot: Stage 6B Slingshot/tradeSeq on a small real-data subset
# Date: 2026-07-28
# Random seed: 20260728

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
if (dir.exists(private_library)) .libPaths(c(private_library, .libPaths()))

required <- c("SeuratObject", "slingshot", "tradeSeq", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Real-data pilot packages missing: ", paste(missing, collapse = ", "))

parameter_path <- file.path(project_dir, "config", "06B_regulatory_parameters.tsv")
object_path <- file.path(
  project_dir, "objects", "GSE201348_5C_epithelial_annotated_CNV.rds"
)
membership_path <- file.path(
  project_dir, "results", "06A_amendment", "stage_blind_module_membership.tsv"
)
candidate_path <- file.path(
  project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
)
parameters <- utils::read.delim(parameter_path, check.names = FALSE)
param <- setNames(parameters$value, parameters$parameter)
p_chr <- function(name) as.character(param[[name]])
p_vec <- function(name) strsplit(p_chr(name), ";", fixed = TRUE)[[1L]]

candidates <- utils::read.delim(candidate_path, check.names = FALSE)
candidates <- candidates[candidates$exploratory_candidate & candidates$passes_LODO, ]
membership <- utils::read.delim(membership_path, check.names = FALSE)
membership <- membership[membership$module_id %in% candidates$module_id, ]
membership <- membership[order(-membership$gene_MAD), ]

object <- readRDS(object_path)
meta <- object[[]]
stage_map <- c(normal = "normal", adenoma_polyp = "adenoma", cancer = "cancer")
meta$stage <- unname(stage_map[meta$lesion_stage])
keep <- meta$epithelial_state %in% p_vec("states")
pilot_meta <- meta[keep, , drop = FALSE]
pilot_meta$stratum <- paste(
  pilot_meta$donor_id, pilot_meta$stage, pilot_meta$epithelial_state, sep = "||"
)
split_index <- split(seq_len(nrow(pilot_meta)), pilot_meta$stratum)
selected <- unlist(lapply(sort(names(split_index)), function(group) {
  idx <- split_index[[group]]
  if (length(idx) < 10L) return(integer())
  set.seed(20260728 + match(group, sort(names(split_index))))
  sort(sample(idx, min(20L, length(idx))))
}), use.names = FALSE)
pilot_meta <- pilot_meta[selected, , drop = FALSE]
pilot_cells <- rownames(pilot_meta)

embedding <- SeuratObject::Embeddings(object, p_chr("reduction"))[
  pilot_cells, seq_len(30L), drop = FALSE
]
sds <- slingshot::slingshot(
  embedding,
  clusterLabels = pilot_meta$epithelial_state,
  start.clus = p_chr("root_state"),
  end.clus = p_vec("terminal_states"),
  stretch = 0
)
pseudotime <- slingshot::slingPseudotime(sds, na = FALSE)
weights <- slingshot::slingCurveWeights(sds)
lineage_definitions <- slingshot::slingLineages(sds)
lineage_terminals <- vapply(
  lineage_definitions, function(x) tail(x, 1L), character(1)
)
keep_lineages <- lineage_terminals %in% p_vec("terminal_states")
pseudotime <- pseudotime[, keep_lineages, drop = FALSE]
weights <- weights[, keep_lineages, drop = FALSE]
lineage_terminals <- lineage_terminals[keep_lineages]
pilot_cells_before_lineage_filter <- length(pilot_cells)
keep_cells <- rowSums(weights) > 0
pseudotime <- pseudotime[keep_cells, , drop = FALSE]
weights <- weights[keep_cells, , drop = FALSE]
pilot_meta <- pilot_meta[keep_cells, , drop = FALSE]
pilot_cells <- rownames(pilot_meta)
if (!ncol(pseudotime) || !identical(dim(pseudotime), dim(weights))) {
  stop("Real-data Slingshot pilot failed")
}
if (!setequal(lineage_terminals, p_vec("terminal_states"))) {
  stop("Real-data Slingshot pilot did not recover every prespecified endpoint")
}

pilot_genes <- head(
  intersect(unique(membership$gene), rownames(object)),
  12L
)
full_counts <- SeuratObject::LayerData(
  object, assay = "RNA", layer = "counts"
)[, pilot_cells, drop = FALSE]
trade_offset <- log(pmax(Matrix::colSums(full_counts), 1))
counts <- full_counts[pilot_genes, , drop = FALSE]
fit <- tradeSeq::fitGAM(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = weights,
  offset = trade_offset,
  nknots = 4L,
  verbose = FALSE,
  parallel = FALSE,
  sce = TRUE
)
association <- tradeSeq::associationTest(fit)
if (
  nrow(association) != length(pilot_genes) ||
    !"pvalue" %in% colnames(association)
) {
  stop("Real-data tradeSeq pilot returned an incompatible result")
}

audit <- data.frame(
  metric = c(
    "pilot_cells_before_lineage_filter", "pilot_cells", "pilot_donors", "pilot_states",
    "slingshot_lineages", "tradeSeq_genes", "finite_tradeSeq_pvalues"
  ),
  value = c(
    pilot_cells_before_lineage_filter, nrow(pilot_meta),
    length(unique(pilot_meta$donor_id)),
    length(unique(pilot_meta$epithelial_state)), ncol(pseudotime),
    nrow(association), sum(is.finite(association$pvalue))
  ),
  stringsAsFactors = FALSE
)
output_dir <- file.path(project_dir, "results", "06B_regulatory_inference", "preflight")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  audit,
  file.path(output_dir, "real_data_pilot_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(audit, row.names = FALSE)
cat("STAGE_6B_REAL_DATA_PILOT_OK\n")
