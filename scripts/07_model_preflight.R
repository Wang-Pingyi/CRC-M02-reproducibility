#!/usr/bin/env Rscript

# Real-data model-input audit. This performs no hypothesis classification and
# must pass before the replication models are run.

set.seed(20260728)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)
get_param <- read_stage7_parameters(project_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
})

locked <- read_locked_stage7_modules(project_dir)
min_fraction <- get_param("global", "module_min_gene_fraction", TRUE)
min_genes <- as.integer(get_param("global", "module_min_genes", TRUE))
min_cells <- as.integer(get_param("global", "min_cells_pseudobulk", TRUE))
min_umi <- get_param("global", "min_umi_pseudobulk", TRUE)

load_eligible <- function(cohort) {
  x <- readRDS(file.path(
    paths$processed, paste0(cohort, "_stage7_pseudobulk_raw_counts.rds")
  ))
  eligible_pseudobulk(x$counts, x$metadata, min_cells, min_umi)
}

g161 <- load_eligible("GSE161277")
g132_raw <- readRDS(file.path(
  paths$processed, "GSE132465_stage7_pseudobulk_raw_counts.rds"
))
g132 <- eligible_pseudobulk(
  g132_raw$counts, g132_raw$metadata, min_cells, min_umi
)

matched161 <- c("Patient1", "Patient2", "Patient3")
m161 <- g161$metadata[
  g161$metadata$major_cell_type == "Epithelial" &
    g161$metadata$donor_id %in% matched161 &
    g161$metadata$condition %in% c("normal", "adenoma", "cancer"),
  ,
  drop = FALSE
]
s161 <- module_scores_from_counts(
  g161$counts[, m161$pseudobulk_id, drop = FALSE],
  locked$membership, min_fraction, min_genes
)
audit <- list()
contrasts <- list(
  normal_to_adenoma = c("normal", "adenoma"),
  adenoma_to_cancer = c("adenoma", "cancer"),
  normal_to_cancer = c("normal", "cancer")
)
for (contrast in names(contrasts)) {
  for (module_id in locked$candidates$module_id) {
    fit <- paired_effect(
      s161$scores, m161, module_id,
      contrasts[[contrast]][1L], contrasts[[contrast]][2L], matched161
    )
    module_scores <- s161$scores[s161$scores$module_id == module_id, ]
    audit[[paste("161", contrast, module_id)]] <- data.frame(
      cohort = "GSE161277", contrast = contrast, module_id = module_id,
      evaluable = all(module_scores$evaluable),
      represented_genes = unique(module_scores$represented_genes),
      locked_genes = unique(module_scores$locked_genes),
      gene_coverage = unique(module_scores$gene_coverage),
      paired_donors = fit$summary$n_donors,
      expected_paired_donors = 3L,
      passed = all(module_scores$evaluable) && fit$summary$n_donors == 3L,
      stringsAsFactors = FALSE
    )
  }
}

matched132 <- sprintf("SMC%02d", 1:10)
eligibility132 <- g132_raw$metadata[
  g132_raw$metadata$major_cell_type == "Epithelial" &
    g132_raw$metadata$donor_id %in% matched132 &
    g132_raw$metadata$condition %in% c("normal", "cancer"),
  c(
    "pseudobulk_id", "donor_id", "condition", "n_cells", "total_umi"
  ),
  drop = FALSE
]
eligibility132$cell_gate_pass <- eligibility132$n_cells >= min_cells
eligibility132$umi_gate_pass <- eligibility132$total_umi >= min_umi
eligibility132$eligible <- eligibility132$cell_gate_pass &
  eligibility132$umi_gate_pass
pair_ok <- vapply(matched132, function(donor) {
  z <- eligibility132[eligibility132$donor_id == donor, ]
  setequal(z$condition[z$eligible], c("normal", "cancer"))
}, logical(1))
eligible_matched132 <- matched132[pair_ok]
eligibility132$paired_analysis_eligible <- eligibility132$donor_id %in%
  eligible_matched132
eligibility132$exclusion_reason <- ifelse(
  eligibility132$paired_analysis_eligible, "included",
  ifelse(
    !eligibility132$cell_gate_pass,
    paste0("fewer_than_", min_cells, "_epithelial_cells"),
    ifelse(
      !eligibility132$umi_gate_pass,
      paste0("fewer_than_", min_umi, "_raw_UMIs"),
      "paired_counterpart_failed_gate"
    )
  )
)
write_stage7_tsv(
  eligibility132,
  file.path(paths$result, "preflight", "GSE132465_pair_eligibility.tsv")
)
m132 <- g132$metadata[
  g132$metadata$major_cell_type == "Epithelial" &
    g132$metadata$condition %in% c("normal", "cancer"),
  ,
  drop = FALSE
]
s132 <- module_scores_from_counts(
  g132$counts[, m132$pseudobulk_id, drop = FALSE],
  locked$membership, min_fraction, min_genes
)
for (module_id in locked$candidates$module_id) {
  fit <- paired_effect(
    s132$scores, m132, module_id, "normal", "cancer", eligible_matched132
  )
  module_scores <- s132$scores[s132$scores$module_id == module_id, ]
  audit[[paste("132", module_id)]] <- data.frame(
    cohort = "GSE132465", contrast = "normal_to_cancer", module_id = module_id,
    evaluable = all(module_scores$evaluable),
    represented_genes = unique(module_scores$represented_genes),
    locked_genes = unique(module_scores$locked_genes),
    gene_coverage = unique(module_scores$gene_coverage),
    paired_donors = fit$summary$n_donors,
    expected_paired_donors = length(eligible_matched132),
    passed = all(module_scores$evaluable) &&
      fit$summary$n_donors == length(eligible_matched132),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, audit)
write_stage7_tsv(
  audit, file.path(paths$result, "preflight", "model_input_audit.tsv")
)
if (any(!audit$passed)) {
  print(audit[!audit$passed, ], row.names = FALSE)
  stop("Stage 7 real-data model preflight failed")
}
cat(
  "STAGE_7_MODEL_PREFLIGHT_OK\trows=", nrow(audit),
  "\tmin_coverage=", min(audit$gene_coverage), "\n", sep = ""
)
