#!/usr/bin/env Rscript

# Register accepted Stage 6B headline values in the project result registry.
# Date: 2026-07-28

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
result_dir <- file.path(project_dir, "results", "06B_regulatory_inference")
registry_path <- file.path(project_dir, "result_registry.tsv")

metrics <- utils::read.delim(
  file.path(result_dir, "stage_6B_key_metrics.tsv"), check.names = FALSE
)
axes <- utils::read.delim(
  file.path(result_dir, "prioritized_ligand_receptor.tsv"), check.names = FALSE
)
registry <- utils::read.delim(registry_path, check.names = FALSE)
registry <- registry[registry$stage != "6B", , drop = FALSE]

unit_map <- c(
  sampled_trajectory_cells = "cells",
  slingshot_lineages = "lineages",
  locked_candidate_modules = "modules",
  locked_candidate_genes = "genes",
  cross_donor_dynamic_gene_lineage_rows = "gene-lineage rows",
  prioritized_regulator_module_pairs = "regulator-module pairs",
  nonredundant_pathway_results = "pathways",
  donor_level_LR_associations_tested = "associations",
  cross_donor_stable_LR_interactions = "associations",
  integrated_candidate_axes = "axes"
)
metric_rows <- data.frame(
  result_id = paste0("6B_", toupper(metrics$metric)),
  stage = "6B",
  metric = metrics$metric,
  value = metrics$value,
  unit = unname(unit_map[metrics$metric]),
  scope = "GSE201348; locked Stage 6A exploratory modules",
  interpretation = paste(
    "Secondary hypothesis-generating computational inference;",
    "does not modify the frozen primary Stage 6A negative result"
  ),
  source_file = "results/06B_regulatory_inference/stage_6B_key_metrics.tsv",
  source_fields = "metric;value",
  generated_by = "scripts/06B_finalize_report.R",
  git_ref = "stage-6B",
  stringsAsFactors = FALSE
)

axis_rows <- do.call(rbind, lapply(seq_len(nrow(axes)), function(i) {
  axis <- axes[i, ]
  data.frame(
    result_id = paste0(
      "6B_", toupper(axis$axis_id),
      c("_RHO", "_FDR", "_DONORS", "_MATCHED_OBSERVATIONS")
    ),
    stage = "6B",
    metric = c(
      "donor_level_spearman_rho", "communication_FDR",
      "unique_donors", "matched_donor_stage_observations"
    ),
    value = c(
      axis$spearman_rho, axis$FDR_communication,
      axis$n_unique_donors, axis$n_matched_observations
    ),
    unit = c("correlation", "probability", "donors", "donor-stage observations"),
    scope = paste(axis$axis_id, axis$module_id, axis$interaction_id, sep = ";"),
    interpretation = paste(
      axis$association_direction,
      "matched donor-level candidate association; not causal"
    ),
    source_file = "results/06B_regulatory_inference/prioritized_ligand_receptor.tsv",
    source_fields = paste(
      "axis_id;module_id;interaction_id;spearman_rho;FDR_communication;",
      "n_unique_donors;n_matched_observations"
    ),
    generated_by = "scripts/06B_liana_donor_communication.R",
    git_ref = "stage-6B",
    stringsAsFactors = FALSE
  )
}))

updated <- rbind(registry, metric_rows, axis_rows)
if (anyDuplicated(updated$result_id)) stop("Duplicate result_id after Stage 6B update")
utils::write.table(
  updated, registry_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
cat("Registered ", nrow(metric_rows) + nrow(axis_rows), " Stage 6B values\n", sep = "")
