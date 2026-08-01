#!/usr/bin/env Rscript

# Reporting: Stage 6B regulatory inference
# Date: 2026-07-28
# Random seed: 20260728

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
if (dir.exists(private_library)) .libPaths(c(private_library, .libPaths()))

result_dir <- file.path(project_dir, "results", "06B_regulatory_inference")
figure_dir <- file.path(project_dir, "figures", "06B_regulatory_inference")
report_dir <- file.path(project_dir, "reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(name) {
  path <- file.path(result_dir, name)
  if (!file.exists(path)) stop("Missing Stage 6B output: ", path)
  utils::read.delim(path, check.names = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

trajectory_audit <- read_tsv("trajectory_audit.tsv")
dynamic <- read_tsv("trajectory_dynamic_genes.tsv")
activity <- read_tsv("regulator_activity.tsv")
regulators <- read_tsv("regulator_module_associations.tsv")
pathways <- read_tsv("pathway_enrichment_pruned.tsv")
communication_audit <- read_tsv("communication_audit.tsv")
prioritized_lr <- read_tsv("prioritized_ligand_receptor.tsv")

metric_value <- function(table, metric) {
  value <- table$value[table$metric == metric]
  if (!length(value)) NA_character_ else as.character(value[[1L]])
}
dynamic_count <- if ("cross_donor_dynamic" %in% colnames(dynamic)) {
  sum(dynamic$cross_donor_dynamic %in% TRUE, na.rm = TRUE)
} else {
  0L
}
regulator_count <- if ("prioritized_regulator" %in% colnames(regulators)) {
  sum(regulators$prioritized_regulator %in% TRUE, na.rm = TRUE)
} else {
  0L
}

key_metrics <- data.frame(
  metric = c(
    "sampled_trajectory_cells", "slingshot_lineages",
    "locked_candidate_modules", "locked_candidate_genes",
    "cross_donor_dynamic_gene_lineage_rows", "prioritized_regulator_module_pairs",
    "nonredundant_pathway_results", "donor_level_LR_associations_tested",
    "cross_donor_stable_LR_interactions", "integrated_candidate_axes"
  ),
  value = c(
    metric_value(trajectory_audit, "sampled_trajectory_cells"),
    metric_value(trajectory_audit, "slingshot_lineages"),
    metric_value(trajectory_audit, "locked_candidate_modules"),
    metric_value(trajectory_audit, "locked_candidate_genes"),
    dynamic_count,
    regulator_count,
    nrow(pathways),
    metric_value(communication_audit, "donor_level_associations_tested"),
    metric_value(communication_audit, "cross_donor_stable_interactions"),
    nrow(prioritized_lr)
  ),
  stringsAsFactors = FALSE
)
write_tsv(key_metrics, file.path(result_dir, "stage_6B_key_metrics.tsv"))

packages <- c(
  "R", "Seurat", "SeuratObject", "slingshot", "tradeSeq", "edgeR",
  "limma", "decoupleR", "dorothea", "msigdbr", "liana", "ggplot2"
)
versions <- data.frame(
  software = packages,
  version = vapply(packages, function(package) {
    if (package == "R") return(as.character(getRversion()))
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
  }, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(versions, file.path(result_dir, "software_versions.tsv"))

axis_lines <- if (nrow(prioritized_lr)) {
  paste0(
    "- ", prioritized_lr$axis_id, ": ", prioritized_lr$axis_label,
    " (predicted association only)."
  )
} else {
  "- No integrated ligand-receptor-regulator-module axis met all locked criteria."
}

report_lines <- c(
  "# Stage 6B regulatory inference",
  "",
  "## Status",
  "",
  "Server analysis completed; independent Codex acceptance remains pending.",
  "The frozen primary Stage 6A negative result was not modified.",
  "",
  "## Locked scope",
  "",
  paste0("- Candidate modules: ", metric_value(trajectory_audit, "locked_candidate_modules"), "."),
  paste0("- Candidate-module genes represented: ", metric_value(trajectory_audit, "locked_candidate_genes"), "."),
  "- Only previously approved epithelial-state labels were used.",
  "- No new program was selected from Stage 6B P values.",
  "",
  "## Slingshot and tradeSeq",
  "",
  paste0("- Deterministically sampled cells: ", metric_value(trajectory_audit, "sampled_trajectory_cells"), "."),
  paste0("- Inferred Slingshot lineages: ", metric_value(trajectory_audit, "slingshot_lineages"), "."),
  paste0("- Donors represented: ", metric_value(trajectory_audit, "donors"), "."),
  paste0("- Candidate gene-lineage rows with tradeSeq support, donor-bin FDR support and cross-donor direction stability: ", dynamic_count, "."),
  "",
  "Stem/progenitor was fixed as the root because it represents the crypt",
  "renewal compartment. Absorptive, Goblet/secretory and Adenoma-like states",
  "were candidate endpoints. The inferred continua describe expression",
  "geometry and must not be interpreted as chronological malignant evolution.",
  "",
  "tradeSeq cell-level tests were auxiliary. Dynamic-gene acceptance required",
  "a donor-by-stage-by-lineage-by-pseudotime-bin raw-count pseudobulk model,",
  "95% confidence intervals, BH FDR and cross-donor sign stability.",
  "",
  "## Predicted transcriptional regulation and pathways",
  "",
  paste0("- Donor-level DoRothEA activity rows: ", nrow(activity), "."),
  paste0("- Prioritized regulator-module associations: ", regulator_count, "."),
  paste0("- Nonredundant Hallmark/Reactome enrichments: ", nrow(pathways), "."),
  "",
  "DoRothEA/decoupleR activities and pathway enrichments are computational",
  "associations. Broad and redundant Reactome terms were removed before",
  "interpretation.",
  "",
  "## Candidate ligand-receptor associations",
  "",
  paste0("- Donor-level ligand-receptor associations tested: ", metric_value(communication_audit, "donor_level_associations_tested"), "."),
  paste0("- Cross-donor stable ligand-receptor associations: ", metric_value(communication_audit, "cross_donor_stable_interactions"), "."),
  paste0("- Integrated candidate axes retained: ", nrow(prioritized_lr), "."),
  "",
  axis_lines,
  "",
  "LIANA supplied the consensus interaction resource. Only interactions with",
  "explicit ligand/cell-surface-ligand and receptor categories were retained.",
  "Statistical filtering",
  "used matched donor-by-stage sender expression, receiver receptor",
  "expression and candidate-module scores. These axes are inferred",
  "associations and do not establish signaling direction or causality.",
  "",
  "## Reproducibility and boundaries",
  "",
  "- Random seed: 20260728.",
  "- Figures were generated as PDF and 300-DPI PNG with source-data tables.",
  "- Large trajectory and tradeSeq objects remain on the server and are excluded from Git.",
  "- No external single-cell, tissue, stool or spatial validation was performed.",
  "- Stage 6B does not authorize the next validation stage automatically.",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(
  report_lines,
  file.path(report_dir, "stage_6B_regulatory_inference.md"),
  useBytes = TRUE
)

output_files <- c(
  "trajectory_source_data.tsv",
  "tradeSeq_candidate_gene_dynamics.tsv",
  "trajectory_dynamic_genes.tsv",
  "regulator_activity.tsv",
  "regulator_module_associations.tsv",
  "pathway_enrichment_all.tsv",
  "pathway_enrichment_pruned.tsv",
  "ligand_receptor_donor_associations.tsv",
  "prioritized_ligand_receptor.tsv",
  "stage_6B_key_metrics.tsv",
  "software_versions.tsv"
)
figure_files <- if (dir.exists(figure_dir)) {
  sort(basename(list.files(figure_dir, full.names = TRUE)))
} else {
  character()
}
manifest <- c(
  "# Stage 6B analysis outputs",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  "## Result tables",
  paste0("- `", output_files, "`"),
  "",
  "## Figures",
  if (length(figure_files)) paste0("- `", figure_files, "`") else "- None generated.",
  "",
  "## Interpretation",
  "- All trajectory, regulator and communication findings are computational inferences.",
  "- The donor is the biological replicate for inferential filtering."
)
writeLines(
  manifest,
  file.path(result_dir, "_analysis_outputs.md"),
  useBytes = TRUE
)
cat("Stage 6B report completed\n")
