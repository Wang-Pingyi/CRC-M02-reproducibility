#!/usr/bin/env Rscript

# Independent acceptance audit: Stage 6B regulatory inference
# Date: 2026-07-28
# Statistical unit: donor; no external validation is performed here.

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
result_dir <- file.path(project_dir, "results", "06B_regulatory_inference")
source_dir <- file.path(result_dir, "source_data")

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing acceptance input: ", path)
  utils::read.delim(path, check.names = FALSE)
}
config <- read_tsv(file.path(
  project_dir, "config", "06B_regulatory_parameters.tsv"
))
param <- setNames(config$value, config$parameter)
p_num <- function(name) as.numeric(param[[name]])
p_vec <- function(name) strsplit(as.character(param[[name]]), ";", fixed = TRUE)[[1L]]

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

builtin <- read_tsv(file.path(result_dir, "validation_checks.tsv"))
trajectory <- read_tsv(file.path(result_dir, "trajectory_source_data.tsv"))
dynamic <- read_tsv(file.path(result_dir, "trajectory_dynamic_genes.tsv"))
regulators <- read_tsv(file.path(result_dir, "regulator_module_associations.tsv"))
communication <- read_tsv(file.path(
  result_dir, "ligand_receptor_donor_associations.tsv"
))
prioritized <- read_tsv(file.path(result_dir, "prioritized_ligand_receptor.tsv"))
resource <- read_tsv(file.path(
  source_dir, "liana_consensus_monomeric_resource.tsv"
))
axis_source <- read_tsv(file.path(
  source_dir, "prioritized_ligand_receptor_source_data.tsv"
))
report_path <- file.path(
  project_dir, "reports", "stage_6B_regulatory_inference.md"
)
report_text <- paste(readLines(report_path, warn = FALSE), collapse = "\n")

add_check(
  "builtin_validation_passed",
  nrow(builtin) > 0L && all(builtin$passed %in% TRUE),
  paste0(sum(builtin$passed %in% TRUE), "/", nrow(builtin))
)

before_hash <- readLines(file.path(
  result_dir, "stage_6A_locked_inputs.before.sha256"
), warn = FALSE)
after_hash <- readLines(file.path(
  result_dir, "stage_6A_locked_inputs.after.sha256"
), warn = FALSE)
add_check(
  "stage_6A_inputs_unchanged",
  identical(before_hash, after_hash),
  paste0("locked files=", length(before_hash))
)

state_donor_counts <- table(trajectory$epithelial_state, trajectory$donor_id)
state_donors <- rowSums(state_donor_counts > 0)
state_top_fraction <- apply(state_donor_counts, 1L, max) /
  rowSums(state_donor_counts)
add_check(
  "trajectory_not_single_donor_dominated",
  all(state_donors >= p_num("min_donors_dynamic_test")) &&
    all(state_top_fraction <= 0.35),
  paste(
    paste0(
      names(state_donors), ":donors=", state_donors,
      ",top=", sprintf("%.3f", state_top_fraction)
    ),
    collapse = ";"
  )
)

accepted_dynamic <- dynamic[dynamic$cross_donor_dynamic %in% TRUE, ]
add_check(
  "accepted_dynamics_reconfirm_thresholds",
  !nrow(accepted_dynamic) || all(
    accepted_dynamic$FDR < p_num("donor_bin_FDR") &
      accepted_dynamic$tradeSeq_FDR < p_num("tradeseq_FDR") &
      accepted_dynamic$evaluable_donors >= p_num("min_donors_dynamic_test") &
      accepted_dynamic$donor_sign_stability >= p_num("donor_sign_stability") &
      sign(accepted_dynamic$CI95_low) == sign(accepted_dynamic$CI95_high)
  ),
  paste0("accepted gene-lineage rows=", nrow(accepted_dynamic))
)

accepted_regulators <- regulators[regulators$prioritized_regulator %in% TRUE, ]
add_check(
  "accepted_regulators_reconfirm_thresholds",
  !nrow(accepted_regulators) || all(
    accepted_regulators$FDR < p_num("regulator_association_FDR") &
      abs(accepted_regulators$standardized_effect) >=
        p_num("min_abs_standardized_effect") &
      accepted_regulators$module_target_count >= p_num("min_module_targets") &
      accepted_regulators$LODO_sign_stability >=
        p_num("regulator_LODO_sign_stability")
  ),
  paste0("prioritized regulator-module rows=", nrow(accepted_regulators))
)

add_check(
  "resource_requires_true_ligand_receptor_categories",
  nrow(resource) > 0L &&
    all(resource$ligand_category %in% p_vec("ligand_categories")) &&
    all(resource$receptor_category %in% p_vec("receptor_categories")) &&
    !anyDuplicated(resource[, c("ligand", "receptor")]),
  paste0("category-valid monomeric pairs=", nrow(resource))
)

stable <- communication[communication$cross_donor_stable %in% TRUE, ]
add_check(
  "stable_communication_reconfirm_thresholds",
  !nrow(stable) || all(
    stable$FDR < p_num("communication_association_FDR") &
      abs(stable$spearman_rho) >= p_num("min_abs_spearman_rho") &
      stable$n_unique_donors >= p_num("min_unique_donors") &
      stable$median_ligand_fraction >= p_num("min_expression_fraction") &
      stable$median_receptor_fraction >= p_num("min_expression_fraction") &
      stable$median_ligand_logCPM >= p_num("min_logCPM") &
      stable$median_receptor_logCPM >= p_num("min_logCPM") &
      stable$LODO_sign_stability >=
        p_num("communication_LODO_sign_stability")
  ),
  paste0("stable category-valid associations=", nrow(stable))
)

axis_ok <- nrow(prioritized) <= p_num("max_candidate_axes")
if (nrow(prioritized)) {
  axis_ok <- axis_ok && all(
    prioritized$ligand_category %in% p_vec("ligand_categories") &
      prioritized$receptor_category %in% p_vec("receptor_categories") &
      prioritized$cross_donor_stable %in% TRUE &
      prioritized$prioritized_regulator %in% TRUE &
      prioritized$n_unique_donors >= p_num("min_unique_donors") &
      prioritized$n_evaluable_LODO >= p_num("min_unique_donors") &
      prioritized$LODO_sign_stability_communication >=
        p_num("communication_LODO_sign_stability") &
      grepl("rho=", prioritized$axis_label, fixed = TRUE)
  )
}
add_check(
  "final_axes_are_category_valid_and_robust",
  axis_ok,
  paste0("axes=", nrow(prioritized))
)

source_ok <- if (!nrow(prioritized)) {
  nrow(axis_source) == 0L
} else {
  all(prioritized$row_id %in% axis_source$row_id) &&
    all(vapply(prioritized$row_id, function(id) {
      nrow(axis_source[axis_source$row_id == id, ]) ==
        prioritized$n_matched_observations[prioritized$row_id == id]
    }, logical(1)))
}
add_check(
  "axis_source_rows_match_reported_observations",
  source_ok,
  paste0("source rows=", nrow(axis_source))
)

forbidden_claims <- c(
  "proves", "demonstrates the mechanism", "drives malignant",
  "establishes signaling direction"
)
add_check(
  "interpretation_boundaries_preserved",
  grepl("not be interpreted as chronological", report_text, fixed = TRUE) &&
    grepl("do not establish signaling direction or causality", report_text, fixed = TRUE) &&
    grepl("No external single-cell, tissue, stool or spatial validation", report_text, fixed = TRUE) &&
    !any(vapply(
      forbidden_claims, grepl, logical(1),
      x = tolower(report_text), fixed = TRUE
    )),
  "pseudotime and communication remain explicitly noncausal; Stage 7 absent"
)

output <- do.call(rbind, checks)
output_path <- file.path(result_dir, "independent_acceptance_checks.tsv")
utils::write.table(
  output, output_path, sep = "\t", quote = FALSE, row.names = FALSE
)
print(output, row.names = FALSE)
if (!all(output$passed)) {
  stop("Independent Stage 6B acceptance failed: ",
       paste(output$check[!output$passed], collapse = ", "))
}
cat("STAGE_6B_INDEPENDENT_ACCEPTANCE_OK\n")
