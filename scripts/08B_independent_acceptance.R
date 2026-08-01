#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08B_independent_acceptance.R PROJECT_DIR RUN_ID")
project <- normalizePath(args[1], mustWork = TRUE)
run_id <- args[2]
result_dir <- file.path(project, "results", "08B_bulk_validation", run_id)
out_file <- file.path(project, "logs_summary", "stage_8B_independent_acceptance.tsv")

read_tsv <- function(name) read.delim(file.path(result_dir, name), check.names = FALSE)
checks <- list()
add <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(pass),
                                               detail = as.character(detail))
}

scores <- read_tsv("bulk_module_scores.tsv")
effects <- read_tsv("bulk_cohort_effects.tsv")
meta <- read_tsv("meta_analysis_results.tsv")
meta_source <- read_tsv("meta_analysis_source_data.tsv")
gate <- read_tsv("tissue_validation_gate.tsv")
lodo <- read_tsv("GSE8671_leave_one_pair_out.tsv")
coverage <- read_tsv("module_mapping_coverage.tsv")
tcga_selection <- read_tsv("TCGA_sample_selection_audit.tsv")
cms <- read_tsv("TCGA_CMScaller_calls.tsv")
tcga_results <- read_tsv("TCGA_auxiliary_results.tsv")
candidates <- read.delim(file.path(project, "results_final",
                                   "stage_6A_exploratory_candidate_modules.tsv"),
                         check.names = FALSE)
modules <- candidates$module_id[candidates$exploratory_candidate %in% c(TRUE, "TRUE")]

# Independent GSE8671 pair reconstruction using explicit joins, never reshape().
g87 <- scores[scores$accession == "GSE8671", ]
add("GSE8671_32_unique_pair_keys", length(unique(g87$donor_id)) == 32L,
    length(unique(g87$donor_id)))
add("GSE8671_pair_key_pattern", all(grepl("^GSE8671_P(0[1-9]|[12][0-9]|3[0-2])$", g87$donor_id)),
    "P01-P32")
manual <- list()
for (module in modules) {
  z <- g87[g87$module_id == module, ]
  normal <- z[z$analysis_group == "normal", c("donor_id", "module_score")]
  adenoma <- z[z$analysis_group == "adenoma", c("donor_id", "module_score")]
  names(normal)[2] <- "normal"
  names(adenoma)[2] <- "adenoma"
  wide <- merge(normal, adenoma, by = "donor_id", all = TRUE)
  add(paste0("GSE8671_complete_pairs_", module),
      nrow(wide) == 32L && all(complete.cases(wide)) &&
        !anyDuplicated(normal$donor_id) && !anyDuplicated(adenoma$donor_id),
      paste0("n=", nrow(wide)))
  difference <- wide$adenoma - wide$normal
  fit <- lm(difference ~ 1)
  vc <- sandwich::vcovHC(fit, type = "HC3")
  tab <- lmtest::coeftest(fit, vcov. = vc)
  estimate <- unname(tab["(Intercept)", 1])
  se <- unname(tab["(Intercept)", 2])
  reported <- effects[effects$accession == "GSE8671" &
                        effects$module_id == module &
                        effects$endpoint == "adenoma_vs_normal" &
                        effects$analysis_set == "primary_all_samples", ]
  manual[[module]] <- data.frame(module_id = module, estimate = estimate, se = se,
                                  n_pairs = length(difference))
  add(paste0("GSE8671_effect_recomputed_", module),
      nrow(reported) == 1L && abs(reported$effect - estimate) < 1e-10 &&
        abs(reported$standard_error - se) < 1e-10 && reported$n_units == 32L,
      paste0("delta_effect=", signif(reported$effect - estimate, 4),
             "; delta_se=", signif(reported$standard_error - se, 4)))
}

# Cross-check every Meta input row against the frozen cohort-effect table.
effect_key <- paste(effects$accession, effects$module_id, effects$endpoint, effects$analysis_set)
source_key <- paste(meta_source$accession, meta_source$module_id,
                    meta_source$endpoint, meta_source$analysis_set)
idx <- match(source_key, effect_key)
add("meta_source_keys_match_cohort_results", !anyNA(idx) && !anyDuplicated(source_key),
    paste0("rows=", nrow(meta_source)))
add("meta_source_values_match_cohort_results",
    !anyNA(idx) &&
      max(abs(meta_source$effect - effects$effect[idx])) < 1e-12 &&
      max(abs(meta_source$standard_error - effects$standard_error[idx])) < 1e-12,
    "effect and SE exact")
add("TCGA_excluded_from_progression_meta", !any(meta_source$accession == "TCGA-COAD"),
    paste(unique(meta_source$accession), collapse = ","))

# Independently refit all reported random-effects models.
meta_delta <- c()
for (i in seq_len(nrow(meta))) {
  row <- meta[i, ]
  z <- meta_source[meta_source$module_id == row$module_id &
                     meta_source$endpoint == row$endpoint &
                     meta_source$analysis_set == row$analysis_set, ]
  fit <- metafor::rma.uni(yi = z$effect, sei = z$standard_error,
                          method = "REML", test = "knha")
  meta_delta <- c(meta_delta, abs(as.numeric(fit$b[1]) - row$pooled_effect),
                  abs(fit$ci.lb - row$ci_low), abs(fit$ci.ub - row$ci_high),
                  abs(fit$I2 - row$I2), abs(fit$tau2 - row$tau2))
}
add("all_meta_models_independently_recomputed", max(meta_delta, na.rm = TRUE) < 1e-9,
    paste0("max_abs_delta=", signif(max(meta_delta, na.rm = TRUE), 4)))
meta$fdr_recomputed <- ave(meta$p_value, interaction(meta$endpoint, meta$analysis_set),
                           FUN = function(p) p.adjust(p, "BH"))
add("meta_FDR_recomputed", max(abs(meta$fdr - meta$fdr_recomputed)) < 1e-12,
    "BH within endpoint and analysis set")

# Rebuild the prespecified early-transition gate.
early <- effects[effects$analysis_set == "primary_all_samples" &
                   effects$endpoint == "adenoma_vs_normal", ]
early <- merge(early, candidates[, c("module_id", "early_effect")], by = "module_id")
expected_class <- sapply(split(early, early$module_id), function(z) {
  concordant <- sign(z$effect) == sign(z$early_effect)
  opposite <- !concordant & z$fdr < 0.05 &
    ((z$ci_low > 0 & z$early_effect < 0) | (z$ci_high < 0 & z$early_effect > 0))
  m <- meta[meta$module_id == z$module_id[1] & meta$endpoint == "adenoma_vs_normal" &
              meta$analysis_set == "primary_all_samples", ]
  significant <- (nrow(m) == 1L && m$fdr < 0.05) || any(z$fdr < 0.05 & concordant)
  if (any(opposite)) "contradictory" else if (sum(concordant) >= 2 && significant) {
    "replicated"
  } else if (sum(concordant) >= 2) {
    "directionally_consistent_but_underpowered"
  } else {
    "not_replicated"
  }
})
expected_class <- expected_class[match(gate$module_id, names(expected_class))]
add("gate_classification_independently_rebuilt",
    identical(unname(expected_class), gate$classification),
    paste(gate$classification, collapse = ","))
add("no_strong_opposite_primary_effects", !any(gate$strong_opposite %in% c(TRUE, "TRUE")),
    "all six modules")

# Sensitivity and leave-one-pair-out audits.
early_all <- effects[effects$endpoint == "adenoma_vs_normal", ]
early_all <- merge(early_all, candidates[, c("module_id", "early_effect")], by = "module_id")
strong_opposite_sensitivity <- with(early_all,
  sign(effect) != sign(early_effect) & fdr < 0.05 &
    ((ci_low > 0 & early_effect < 0) | (ci_high < 0 & early_effect > 0)))
add("no_strong_opposite_sensitivity_effects", !any(strong_opposite_sensitivity),
    paste0("strong_opposite=", sum(strong_opposite_sensitivity)))
add("all_three_sensitivity_sets_present",
    setequal(unique(effects$analysis_set),
             c("primary_all_samples", "exclude_high_leverage", "exclude_all_qc_flags")),
    paste(unique(effects$analysis_set), collapse = ","))
add("GSE8671_LODO_32_per_module", all(table(lodo$module_id) == 32L),
    paste(range(table(lodo$module_id)), collapse = "-"))
lodo_expected <- candidates$early_effect[match(lodo$module_id, candidates$module_id)]
add("GSE8671_LODO_direction_stable", all(sign(lodo$effect) == sign(lodo_expected)),
    paste0(sum(sign(lodo$effect) == sign(lodo_expected)), "/", nrow(lodo)))

# Mapping and TCGA auxiliary boundaries.
add("module_mapping_coverage_reported",
    nrow(coverage) == 18L && all(coverage$mapped_genes > 0) && all(coverage$coverage <= 1),
    paste0("range=", paste(signif(range(coverage$coverage), 3), collapse = "-")))
selected <- tcga_selection[tcga_selection$selected %in% c(TRUE, "TRUE"), ]
add("TCGA_one_selected_file_per_patient_sample_type",
    !anyDuplicated(paste(selected$patient_id, selected$sample_type)),
    paste0("selected=", nrow(selected)))
add("TCGA_results_labeled_auxiliary",
    all(tcga_results$accession == "TCGA-COAD") &&
      all(tcga_results$analysis_set == "auxiliary"),
    "TCGA auxiliary table is separate")
add("CMS_unclassified_not_forced",
    all(is.na(cms$prediction[cms$FDR > 0.05])),
    paste0("unclassified=", sum(is.na(cms$prediction)), "/", nrow(cms)))

# Artifact and stage-boundary checks.
required_final <- c(
  file.path(project, "results_final", "stage_8B_bulk_validation_summary.tsv"),
  file.path(project, "results_final", "stage_8B_meta_analysis_results.tsv"),
  file.path(project, "figures_final", "stage_8B_early_transition_forest.png"),
  file.path(project, "figures_final", "stage_8B_early_transition_forest.pdf"),
  file.path(project, "figures_final", "stage_8B_early_transition_forest_source_data.tsv"),
  file.path(project, "reports", "stage_8B_bulk_validation.md")
)
add("canonical_artifacts_present", all(file.exists(required_final)),
    paste0(sum(file.exists(required_final)), "/", length(required_final)))
add("stage9_not_started",
    !dir.exists(file.path(project, "results", "09_stool_validation")) &&
      !dir.exists(file.path(project, "results", "09_stool_RNA_validation")),
    "no Stage 9 result directory")

checks <- do.call(rbind, checks)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.table(checks, out_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
if (!all(checks$pass)) {
  stop("Independent Stage 8B acceptance failed: ",
       paste(checks$check[!checks$pass], collapse = ", "))
}
cat("INDEPENDENT_ACCEPTANCE_PASS", sum(checks$pass), "/", nrow(checks), "\n")
