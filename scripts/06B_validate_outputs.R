#!/usr/bin/env Rscript

# Validation: Stage 6B regulatory inference outputs
# Date: 2026-07-28

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
result_dir <- file.path(project_dir, "results", "06B_regulatory_inference")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_dir, "figures", "06B_regulatory_inference")
report_path <- file.path(project_dir, "reports", "stage_6B_regulatory_inference.md")
parameter_path <- file.path(project_dir, "config", "06B_regulatory_parameters.tsv")

parameters <- utils::read.delim(parameter_path, check.names = FALSE)
param <- setNames(parameters$value, parameters$parameter)
p_num <- function(name) as.numeric(param[[name]])
p_chr <- function(name) as.character(param[[name]])
p_vec <- function(name) strsplit(p_chr(name), ";", fixed = TRUE)[[1L]]

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing validation input: ", path)
  utils::read.delim(path, check.names = FALSE)
}
result_file <- function(name) file.path(result_dir, name)
source_file <- function(name) file.path(source_dir, name)

required_files <- c(
  result_file("trajectory_source_data.tsv"),
  result_file("tradeSeq_candidate_gene_dynamics.tsv"),
  result_file("trajectory_dynamic_genes.tsv"),
  result_file("trajectory_pseudobulk_manifest.tsv"),
  result_file("regulator_activity.tsv"),
  result_file("regulator_module_associations.tsv"),
  result_file("pathway_enrichment_all.tsv"),
  result_file("pathway_enrichment_pruned.tsv"),
  result_file("communication_group_manifest.tsv"),
  result_file("ligand_receptor_donor_associations.tsv"),
  result_file("prioritized_ligand_receptor.tsv"),
  result_file("stage_6B_key_metrics.tsv"),
  result_file("software_versions.tsv"),
  result_file("_analysis_outputs.md"),
  source_file("trajectory_embedding.tsv"),
  source_file("trajectory_module_scores.tsv"),
  source_file("liana_consensus_monomeric_resource.tsv"),
  source_file("prioritized_ligand_receptor_source_data.tsv"),
  report_path
)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}
add_check(
  "required_outputs_exist",
  all(file.exists(required_files)),
  paste0(sum(file.exists(required_files)), "/", length(required_files))
)
if (!all(file.exists(required_files))) {
  write.table(
    do.call(rbind, checks),
    result_file("validation_checks.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  stop("Required Stage 6B outputs are missing")
}

trajectory <- read_tsv(result_file("trajectory_source_data.tsv"))
trade <- read_tsv(result_file("tradeSeq_candidate_gene_dynamics.tsv"))
dynamic <- read_tsv(result_file("trajectory_dynamic_genes.tsv"))
trajectory_manifest <- read_tsv(result_file("trajectory_pseudobulk_manifest.tsv"))
activity <- read_tsv(result_file("regulator_activity.tsv"))
regulators <- read_tsv(result_file("regulator_module_associations.tsv"))
pathways <- read_tsv(result_file("pathway_enrichment_pruned.tsv"))
communication <- read_tsv(result_file("ligand_receptor_donor_associations.tsv"))
prioritized <- read_tsv(result_file("prioritized_ligand_receptor.tsv"))
communication_resource <- read_tsv(
  source_file("liana_consensus_monomeric_resource.tsv")
)
metrics <- read_tsv(result_file("stage_6B_key_metrics.tsv"))
versions <- read_tsv(result_file("software_versions.tsv"))
locked <- read_tsv(file.path(
  project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
))
locked <- locked[locked$exploratory_candidate & locked$passes_LODO, ]

add_check(
  "locked_candidate_count",
  nrow(locked) == 6L && !anyDuplicated(locked$module_id),
  paste0("modules=", nrow(locked))
)
add_check(
  "trajectory_cells_unique",
  nrow(trajectory) > 0L && !anyDuplicated(trajectory$cell_id),
  paste0("cells=", nrow(trajectory))
)
add_check(
  "trajectory_states_locked",
  all(trajectory$epithelial_state %in% p_vec("states")),
  paste(sort(unique(trajectory$epithelial_state)), collapse = ";")
)
add_check(
  "trajectory_root_present",
  p_chr("root_state") %in% trajectory$epithelial_state,
  p_chr("root_state")
)
add_check(
  "trajectory_has_cross_donor_support",
  length(unique(trajectory$donor_id)) >= p_num("min_donors_dynamic_test"),
  paste0("donors=", length(unique(trajectory$donor_id)))
)
add_check(
  "trajectory_values_valid",
  all(is.finite(trajectory$assigned_weight)) &&
    all(trajectory$assigned_weight >= 0 & trajectory$assigned_weight <= 1) &&
    all(is.finite(trajectory$assigned_pseudotime)),
  "finite pseudotime and weights in [0,1]"
)
add_check(
  "trajectory_pseudobulk_unit_valid",
  !anyDuplicated(trajectory_manifest$group_id) &&
    all(c("donor_id", "stage", "lineage", "pseudotime_bin") %in%
      colnames(trajectory_manifest)),
  paste0("pseudobulks=", nrow(trajectory_manifest))
)
add_check(
  "tradeSeq_statistics_valid",
  nrow(trade) > 0L &&
    all(is.finite(trade$pvalue) & trade$pvalue >= 0 & trade$pvalue <= 1) &&
    all(is.finite(trade$tradeSeq_FDR) &
      trade$tradeSeq_FDR >= 0 & trade$tradeSeq_FDR <= 1),
  paste0("genes=", nrow(trade))
)
required_dynamic <- c(
  "gene", "lineage", "standardized_pseudotime_effect", "CI95_low",
  "CI95_high", "p_value", "FDR", "evaluable_donors",
  "donor_sign_stability", "cross_donor_dynamic"
)
add_check(
  "donor_dynamic_schema",
  all(required_dynamic %in% colnames(dynamic)),
  paste(required_dynamic, collapse = ";")
)
if (nrow(dynamic)) {
  add_check(
    "donor_dynamic_statistics_valid",
    all(is.finite(dynamic$p_value) & dynamic$p_value >= 0 & dynamic$p_value <= 1) &&
      all(is.finite(dynamic$FDR) & dynamic$FDR >= 0 & dynamic$FDR <= 1) &&
      all(dynamic$CI95_low <= dynamic$standardized_pseudotime_effect &
        dynamic$standardized_pseudotime_effect <= dynamic$CI95_high),
    paste0("rows=", nrow(dynamic))
  )
  accepted_dynamic <- dynamic[dynamic$cross_donor_dynamic %in% TRUE, ]
  add_check(
    "accepted_dynamic_rules",
    !nrow(accepted_dynamic) || all(
      accepted_dynamic$FDR < p_num("donor_bin_FDR") &
        accepted_dynamic$tradeSeq_FDR < p_num("tradeseq_FDR") &
        accepted_dynamic$evaluable_donors >= p_num("min_donors_dynamic_test") &
        accepted_dynamic$donor_sign_stability >= p_num("donor_sign_stability")
    ),
    paste0("accepted=", nrow(accepted_dynamic))
  )
}

add_check(
  "regulator_activity_donor_level",
  nrow(activity) > 0L &&
    all(c("pseudobulk_id", "donor_id", "stage", "epithelial_state", "regulator", "score") %in%
      colnames(activity)),
  paste0("rows=", nrow(activity), ";donors=", length(unique(activity$donor_id)))
)
add_check(
  "regulator_association_statistics_valid",
  nrow(regulators) > 0L &&
    all(is.finite(regulators$p_value) & regulators$p_value >= 0 &
      regulators$p_value <= 1) &&
    all(is.finite(regulators$FDR) & regulators$FDR >= 0 & regulators$FDR <= 1) &&
    all(regulators$CI95_low <= regulators$standardized_effect &
      regulators$standardized_effect <= regulators$CI95_high),
  paste0("rows=", nrow(regulators))
)
accepted_regulators <- regulators[regulators$prioritized_regulator %in% TRUE, ]
add_check(
  "prioritized_regulator_rules",
  !nrow(accepted_regulators) || all(
    accepted_regulators$FDR < p_num("regulator_association_FDR") &
      abs(accepted_regulators$standardized_effect) >=
        p_num("min_abs_standardized_effect") &
      accepted_regulators$module_target_count >= p_num("min_module_targets") &
      accepted_regulators$LODO_sign_stability >=
        p_num("regulator_LODO_sign_stability")
  ),
  paste0("prioritized=", nrow(accepted_regulators))
)

add_check(
  "pathway_results_locked_and_nonredundant",
  !nrow(pathways) || all(
    pathways$module_id %in% locked$module_id &
      pathways$overlap_count >= p_num("min_overlap") &
      pathways$FDR < p_num("pathway_FDR") &
      pathways$pathway_genes_in_background >= p_num("min_geneset_size") &
      pathways$pathway_genes_in_background <= p_num("max_geneset_size")
  ),
  paste0("retained=", nrow(pathways))
)

if (nrow(communication)) {
  add_check(
    "communication_resource_categories_valid",
    nrow(communication_resource) > 0L &&
      all(communication_resource$ligand_category %in% p_vec("ligand_categories")) &&
      all(communication_resource$receptor_category %in% p_vec("receptor_categories")),
    paste0("category-valid interactions=", nrow(communication_resource))
  )
  add_check(
    "communication_statistics_valid",
    all(is.finite(communication$p_value) &
      communication$p_value >= 0 & communication$p_value <= 1) &&
      all(is.finite(communication$FDR) &
        communication$FDR >= 0 & communication$FDR <= 1),
    paste0("tested=", nrow(communication))
  )
  stable_communication <- communication[communication$cross_donor_stable %in% TRUE, ]
  add_check(
    "stable_communication_rules",
    !nrow(stable_communication) || all(
      stable_communication$FDR < p_num("communication_association_FDR") &
        abs(stable_communication$spearman_rho) >= p_num("min_abs_spearman_rho") &
        stable_communication$n_unique_donors >= p_num("min_unique_donors") &
        stable_communication$LODO_sign_stability >=
          p_num("communication_LODO_sign_stability")
    ),
    paste0("stable=", nrow(stable_communication))
  )
}
add_check(
  "candidate_axis_limit",
  nrow(prioritized) <= p_num("max_candidate_axes"),
  paste0("axes=", nrow(prioritized))
)
if (nrow(prioritized)) {
  add_check(
    "candidate_axes_integrated_and_stable",
    all(
      prioritized$cross_donor_stable %in% TRUE &
        prioritized$prioritized_regulator %in% TRUE &
        prioritized$module_id %in% locked$module_id &
        prioritized$ligand_category %in% p_vec("ligand_categories") &
        prioritized$receptor_category %in% p_vec("receptor_categories")
    ),
    paste(prioritized$axis_id, collapse = ";")
  )
}

figure_pairs <- c(
  "trajectory_locked_states", "trajectory_pseudotime",
  "regulator_module_associations"
)
figure_ok <- all(vapply(figure_pairs, function(stem) {
  file.exists(file.path(figure_dir, paste0(stem, ".pdf"))) &&
    file.exists(file.path(figure_dir, paste0(stem, ".png")))
}, logical(1)))
if (nrow(prioritized)) {
  figure_ok <- figure_ok &&
    file.exists(file.path(figure_dir, "prioritized_ligand_receptor.pdf")) &&
    file.exists(file.path(figure_dir, "prioritized_ligand_receptor.png"))
}
add_check("figure_pairs_present", figure_ok, paste(figure_pairs, collapse = ";"))
add_check(
  "figure_source_data_present",
  all(file.exists(c(
    source_file("trajectory_embedding.tsv"),
    source_file("top_regulator_effects.tsv"),
    source_file("prioritized_ligand_receptor_source_data.tsv")
  ))),
  "trajectory, regulator and communication source data"
)

report_text <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
add_check(
  "interpretation_boundaries_present",
  grepl("not be interpreted as chronological", report_text, fixed = TRUE) &&
    grepl("do not establish signaling direction or causality", report_text, fixed = TRUE),
  "trajectory and communication limitations"
)
add_check(
  "software_versions_recorded",
  nrow(versions) >= 10L && all(nzchar(versions$version)),
  paste0("rows=", nrow(versions))
)
add_check(
  "key_metrics_unique",
  nrow(metrics) >= 10L && !anyDuplicated(metrics$metric),
  paste0("metrics=", nrow(metrics))
)

validation <- do.call(rbind, checks)
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
write_tsv(validation, result_file("validation_checks.tsv"))
cat(
  "Stage 6B validation:", sum(validation$passed), "/", nrow(validation), "\n"
)
if (!all(validation$passed)) {
  failed <- validation$check[!validation$passed]
  stop("Stage 6B validation failed: ", paste(failed, collapse = "; "))
}
