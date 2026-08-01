#!/usr/bin/env Rscript

# Validation for the Stage 6A exploratory amendment.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1L]] else getwd()
result_dir <- file.path(project_dir, "results", "06A_amendment")
figure_dir <- file.path(project_dir, "figures", "06A_amendment")

read_tsv <- function(name) {
  path <- file.path(result_dir, name)
  if (!file.exists(path)) stop("Missing output: ", path)
  utils::read.delim(path, check.names = FALSE)
}

checks <- list()
add_check <- function(check, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

paired_audit <- read_tsv("paired_FAP_state_audit.tsv")
paired_results <- read_tsv("paired_FAP_gene_results.tsv")
paired_hits <- read_tsv("paired_FAP_gene_hits.tsv")
paired_lodo <- read_tsv("paired_FAP_gene_LODO_summary.tsv")
membership <- read_tsv("stage_blind_module_membership.tsv")
module_scores <- read_tsv("source_data/stage_blind_module_scores.tsv")
module_results <- read_tsv("stage_blind_module_results.tsv")
paired_module_results <- read_tsv("paired_FAP_module_results.tsv")
module_audit <- read_tsv("stage_blind_module_audit.tsv")
candidate_modules <- read_tsv("exploratory_candidate_modules.tsv")
metrics <- read_tsv("stage_6A_amendment_key_metrics.tsv")

valid_probability <- function(x) {
  all(is.finite(x) & x >= 0 & x <= 1)
}
valid_ci <- function(x) {
  all(
    is.finite(x$CI95_low) & is.finite(x$CI95_high) &
      x$CI95_low <= x$effect & x$effect <= x$CI95_high
  )
}

add_check(
  "paired_state_audit_unique",
  !anyDuplicated(paired_audit$epithelial_state),
  nrow(paired_audit)
)
add_check(
  "paired_results_nonempty",
  nrow(paired_results) > 0L,
  nrow(paired_results)
)
add_check(
  "paired_results_use_at_least_three_pairs",
  all(paired_results$n_paired_donors >= 3L),
  paste(sort(unique(paired_results$n_paired_donors)), collapse = ";")
)
add_check(
  "paired_results_valid_statistics",
  valid_probability(paired_results$p_value) &
    valid_probability(paired_results$FDR) &
    valid_ci(paired_results),
  "P, FDR and 95% CI"
)
add_check(
  "paired_hits_meet_locked_thresholds",
  nrow(paired_hits) == 0L ||
    all(paired_hits$FDR <= 0.05 & abs(paired_hits$effect) >= 0.25),
  nrow(paired_hits)
)
add_check(
  "paired_LODO_only_for_hits",
  nrow(paired_lodo) == 0L ||
    all(paired_lodo$gene %in% paired_hits$gene),
  nrow(paired_lodo)
)
add_check(
  "module_construction_stage_blind",
  all(!module_audit$stage_labels_used_for_construction) &
    all(!membership$stage_labels_used_for_construction),
  "stage labels absent from construction"
)
add_check(
  "module_membership_unique",
  !anyDuplicated(membership[, c("module_id", "gene")]),
  nrow(membership)
)
testable_membership <- membership[membership$testable, , drop = FALSE]
add_check(
  "testable_modules_meet_size_and_coherence",
  nrow(testable_membership) > 0L &&
    all(testable_membership$module_size >= 10L) &
    all(testable_membership$median_within_module_correlation >= 0.20),
  length(unique(testable_membership$module_id))
)
add_check(
  "module_scores_unique",
  !anyDuplicated(module_scores[, c("module_id", "pseudobulk_id")]),
  nrow(module_scores)
)
add_check(
  "module_results_nonempty",
  nrow(module_results) > 0L,
  nrow(module_results)
)
add_check(
  "module_results_valid_statistics",
  valid_probability(module_results$p_value) &
    valid_probability(module_results$FDR) &
    valid_ci(module_results),
  "P, FDR and 95% CI"
)
add_check(
  "module_results_donor_level",
  all(
    module_results$n_normal_donors >= 3L &
      module_results$n_adenoma_donors >= 3L
  ),
  "minimum donor rule"
)
add_check(
  "paired_module_results_valid",
  nrow(paired_module_results) == 0L ||
    (
      valid_probability(paired_module_results$p_value) &
        valid_probability(paired_module_results$FDR) &
        valid_ci(paired_module_results) &
        all(paired_module_results$n_normal_donors >= 3L)
    ),
  nrow(paired_module_results)
)
add_check(
  "candidate_modules_meet_locked_rules",
  nrow(candidate_modules) == 0L ||
    all(
      candidate_modules$early_FDR <= 0.05 &
        abs(candidate_modules$early_effect) >= 0.50 &
        candidate_modules$paired_direction_concordant &
        candidate_modules$passes_LODO &
        candidate_modules$LODO_sign_stability >= 0.75
    ),
  sum(candidate_modules$exploratory_candidate, na.rm = TRUE)
)
add_check(
  "primary_freeze_hashes_match",
  file.exists(file.path(result_dir, "primary_freeze_sha256.before.tsv")) &&
    file.exists(file.path(result_dir, "primary_freeze_sha256.after.tsv")) &&
    identical(
      readLines(file.path(result_dir, "primary_freeze_sha256.before.tsv")),
      readLines(file.path(result_dir, "primary_freeze_sha256.after.tsv"))
    ),
  "primary Stage 6A outputs unchanged"
)
add_check(
  "report_exists",
  file.exists(file.path(project_dir, "reports", "stage_6A_exploratory_amendment.md")),
  "amendment report"
)
add_check(
  "software_versions_recorded",
  file.exists(file.path(result_dir, "software_versions.tsv")),
  "software_versions.tsv"
)
add_check(
  "key_metrics_complete",
  all(c(
    "paired_states_evaluable", "paired_gene_FDR_hits",
    "stage_blind_modules_testable",
    "exploratory_progressive_modules_final"
  ) %in% metrics$metric),
  paste(metrics$metric, collapse = ";")
)

figure_requirements <- list(
  paired_FAP_gene_hits = file.path(
    result_dir, "source_data", "paired_FAP_hit_counts.tsv"
  ),
  top_stage_blind_module_effects = file.path(
    result_dir, "source_data", "top_stage_blind_module_effects.tsv"
  )
)
for (figure_name in names(figure_requirements)) {
  pdf_path <- file.path(figure_dir, paste0(figure_name, ".pdf"))
  png_path <- file.path(figure_dir, paste0(figure_name, ".png"))
  source_path <- figure_requirements[[figure_name]]
  add_check(
    paste0("figure_", figure_name),
    file.exists(pdf_path) && file.info(pdf_path)$size > 0L &&
      file.exists(png_path) && file.info(png_path)$size > 0L &&
      file.exists(source_path) && file.info(source_path)$size > 0L,
    "PDF, PNG and source data"
  )
}

checks <- do.call(rbind, checks)
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
write_tsv(checks, file.path(result_dir, "validation_checks.tsv"))
cat("Stage 6A amendment validation:", sum(checks$passed), "/", nrow(checks), "\n")
if (!all(checks$passed)) {
  print(checks[!checks$passed, , drop = FALSE])
  quit(status = 1L)
}
