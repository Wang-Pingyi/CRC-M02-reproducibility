#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08B_validate_outputs.R PROJECT_DIR RUN_ID")
project <- normalizePath(args[1], mustWork = TRUE)
run_id <- args[2]
source(file.path(project, "scripts", "08B_helpers.R"))
paths <- stage8b_paths(project, run_id)
modules <- stage8b_locked(paths)$candidates$module_id
checks <- list()
add <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(pass), detail = detail)
}

required <- c("bulk_module_scores.tsv", "module_mapping_coverage.tsv",
  "bulk_cohort_effects.tsv", "bulk_validation_summary.tsv", "meta_analysis_results.tsv",
  "meta_analysis_source_data.tsv", "GSE8671_leave_one_pair_out.tsv",
  "GSE8671_leave_one_pair_out_summary.tsv", "TCGA_sample_selection_audit.tsv",
  "TCGA_module_scores_clinical.tsv", "TCGA_CMScaller_calls.tsv",
  "TCGA_auxiliary_results.tsv", "tissue_validation_gate.tsv")
for (f in required) add(paste0("file_", f), file.exists(file.path(paths$result, f)), f)

if (all(vapply(required, function(f) file.exists(file.path(paths$result, f)), logical(1)))) {
  effects <- read.delim(file.path(paths$result, "bulk_cohort_effects.tsv"), check.names = FALSE)
  meta <- read.delim(file.path(paths$result, "meta_analysis_results.tsv"), check.names = FALSE)
  scores <- read.delim(file.path(paths$result, "bulk_module_scores.tsv"), check.names = FALSE)
  gate <- read.delim(file.path(paths$result, "tissue_validation_gate.tsv"), check.names = FALSE)
  coverage <- read.delim(file.path(paths$result, "module_mapping_coverage.tsv"), check.names = FALSE)
  add("six_locked_modules_scores", setequal(unique(scores$module_id), modules), paste(unique(scores$module_id), collapse = ","))
  add("three_independent_tissue_cohorts", setequal(unique(scores$accession), c("GSE41657", "GSE100179", "GSE8671")), paste(unique(scores$accession), collapse = ","))
  add("unique_sample_module_rows", !anyDuplicated(scores[c("accession", "sample_id", "module_id")]), nrow(scores))
  add("finite_module_scores", all(is.finite(scores$module_score)), "all scores finite")
  add("mapping_coverage_nonzero", all(coverage$mapped_genes > 0 & coverage$coverage <= 1), paste(range(coverage$coverage), collapse = "-"))
  add("primary_early_all_modules_all_cohorts",
      nrow(effects[effects$analysis_set == "primary_all_samples" & effects$endpoint == "adenoma_vs_normal", ]) == 18,
      "expected 6 modules x 3 cohorts")
  add("gse8671_paired_model_only",
      all(effects$model[effects$accession == "GSE8671"] == "verified_pair_difference"),
      "verified donor-pair differences")
  g87 <- scores[scores$accession == "GSE8671", ]
  g87_key <- paste(g87$module_id, g87$donor_id, g87$analysis_group)
  add("gse8671_exactly_32_explicit_pairs",
      length(unique(g87$donor_id)) == 32L &&
        !anyDuplicated(g87_key) &&
        all(table(g87$module_id, g87$analysis_group) == 32L),
      paste0("unique_pairs=", length(unique(g87$donor_id)), "; duplicate module-pair-condition rows=", sum(duplicated(g87_key))))
  add("effects_finite", all(is.finite(effects$effect) & is.finite(effects$standard_error) & effects$standard_error > 0), nrow(effects))
  add("effects_ci_order", all(effects$ci_low <= effects$effect & effects$effect <= effects$ci_high), "effect within CI")
  add("effects_fdr_valid", all(effects$fdr >= 0 & effects$fdr <= 1), "BH FDR")
  add("meta_all_four_endpoints", setequal(unique(meta$endpoint), c("adenoma_vs_normal", "cancer_vs_adenoma", "cancer_vs_normal", "ordered_trend")), paste(unique(meta$endpoint), collapse = ","))
  add("meta_random_effects", all(meta$method == "REML_Knapp-Hartung"), "REML/Knapp-Hartung")
  add("heterogeneity_reported", all(is.finite(meta$tau2) & is.finite(meta$I2) & is.finite(meta$Q) & is.finite(meta$Q_p_value)), "tau2 I2 Q Qp")
  add("gate_all_modules", setequal(gate$module_id, modules), paste(gate$classification, collapse = ","))
  add("gate_valid_classes", all(gate$classification %in% c("replicated", "directionally_consistent_but_underpowered", "not_replicated", "contradictory")), "prespecified classes")
}

figures <- c("stage_8B_early_transition_forest.png", "stage_8B_early_transition_forest.pdf",
             "stage_8B_early_transition_forest_source_data.tsv")
for (f in figures) add(paste0("figure_", f), file.exists(file.path(paths$figure, f)), f)
for (f in figures) add(paste0("canonical_figure_", f), file.exists(file.path(project, "figures_final", f)), f)
add("report_exists", file.exists(file.path(project, "reports", "stage_8B_bulk_validation.md")), "stage report")
add("status_updated", any(grepl(run_id, readLines(file.path(project, "STATUS.md"), warn = FALSE), fixed = TRUE)), run_id)

checks <- do.call(rbind, checks)
stage8b_write_tsv(checks, file.path(paths$result, "validation_checks.tsv"))
stage8b_write_tsv(checks, file.path(project, "logs_summary", "stage_8B_validation_checks.tsv"))
if (!all(checks$pass)) {
  writeLines(paste("Stage 8B validation failed:", paste(checks$check[!checks$pass], collapse = ", ")),
             file.path(paths$log, "NEEDS_CODEX_ATTENTION"))
  quit(status = 1)
}
writeLines(paste("Stage 8B server run ready for independent Codex QC", run_id),
           file.path(paths$log, "READY_FOR_CODEX_QC"))
