#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
project_lib <- file.path(root, "environment", "R-library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

results_dir <- file.path(root, "results", "stage10c3")
figures_dir <- file.path(root, "figures", "stage10c3")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
required <- file.path(results_dir, c("STAGE10C3_PATIENT_SAMPLE_SCORES.tsv", "STAGE10C3_MATCHED_PAIR_EFFECTS.tsv", "STAGE10C3_SENSITIVITY.tsv", "STAGE10C3_GENE_COVERAGE.tsv", "QC_PREM02.SUCCESS"))
if (!all(file.exists(required))) stop("Cannot finalize: required completed numeric tables are missing")
if (file.exists(file.path(results_dir, "STAGE10C3_DECISION.md"))) stop("Refusing to overwrite an existing Stage 10C3 decision")

scores <- fread(required[1], na.strings = "NA")
effects <- fread(required[2], na.strings = "NA")
sensitivity <- fread(required[3], na.strings = "NA")
coverage <- fread(required[4], na.strings = "NA")
coverage[, status := ifelse(primary_36of36, "PRIMARY_ESTIMABLE", "SENSITIVITY_ONLY_PRIMARY_NOT_ESTIMABLE")]
fwrite(coverage, required[4], sep = "\t", na = "NA")

expected_primary_bundle <- "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
expected_sens_bundle <- "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8"
primary_indiv <- effects[record_type == "individual_patient_effect" & scenario_id == "PRIMARY_ALL_EPI_MIN50" & score_method == "M02_SCORE_V1"]
sens35_indiv <- effects[record_type == "individual_patient_effect" & scenario_id == "PRIMARY_ALL_EPI_MIN50" & score_method == "M02_MINUS_INPP5D_SENS_V1"]
primary_estimable <- nrow(primary_indiv[status == "ESTIMABLE"]) == 2L
primary_both_positive <- primary_estimable && all(primary_indiv$effect_lesion_minus_normal > 0)
sens35_no_reverse <- nrow(sens35_indiv[status == "ESTIMABLE"]) == 2L && all(sens35_indiv$effect_lesion_minus_normal >= 0)
systemic_opposite <- any(sensitivity[status == "ESTIMABLE" & n_pairs == 2L,
                                      is.finite(P1_effect) & is.finite(P5_effect) & P1_effect < 0 & P5_effect < 0])
if (!primary_estimable) {
  decision <- "NOT_ESTIMABLE"
} else if (primary_both_positive && sens35_no_reverse && !systemic_opposite) {
  decision <- "LST_DIRECTIONAL_CONCORDANCE"
} else if (all(primary_indiv$effect_lesion_minus_normal <= 0) && !any(sensitivity[status == "ESTIMABLE", k_positive == 2L])) {
  decision <- "NULL_OR_OPPOSITE"
} else {
  decision <- "MIXED_OR_METHOD_DEPENDENT"
}

display_scores <- scores[scenario_id == "PRIMARY_ALL_EPI_MIN50" & score_method == "M02_29OF36_COVERAGE_SENS"]
secondary <- copy(display_scores)
secondary[, description_type := ifelse(biological_sample_type == "normal mucosa", "normal_sample", "lesion_sample")]
p3_l <- display_scores[sample_id == "P3_L", score]
p3_p <- display_scores[sample_id == "P3_P", score]
if (length(p3_l) == 1L && length(p3_p) == 1L) {
  secondary <- rbind(secondary, data.table(
    sample_id = "P3_L-minus-P3_P", scenario_id = "PRIMARY_ALL_EPI_MIN50", qc_mode = "primary_qc",
    annotation_source = "independent_broad_annotation", compartment = "all_epithelial", min_cells = 50L,
    scenario_role = "secondary_descriptive", score_method = "M02_29OF36_COVERAGE_SENS", score_role = "coverage_sensitivity_descriptive",
    n_cells = NA_integer_, mapped_genes = 33L, required_genes = 29L, score = p3_l - p3_p,
    normalization = "within-patient descriptive contrast", inferential_unit = "patient", status = "ESTIMABLE",
    patient_id = "P3", biological_sample_type = "two lesions", lesion_morphology = "LST-G minus protruded adenoma",
    histology = "within-patient descriptive", anatomic_site = "different sites", pair_id = "P3_MULTI_LESION",
    FAP_status = NA_character_, germline_APC_status = NA_character_, description_type = "within_patient_lesion_contrast"
  ), fill = TRUE)
}
fwrite(secondary, file.path(results_dir, "STAGE10C3_SECONDARY_DESCRIPTIVE.tsv"), sep = "\t", na = "NA")

fig1 <- display_scores[sample_id %in% c("P1_N", "P1_L", "P5_N", "P5_L")]
fig1[, tissue_order := factor(ifelse(biological_sample_type == "normal mucosa", "Normal", "LST-G"), levels = c("Normal", "LST-G"))]
fwrite(fig1, file.path(figures_dir, "FIGURE10C3_1_source_data.tsv"), sep = "\t", na = "NA")
p1 <- ggplot(fig1, aes(tissue_order, score, group = patient_id, color = patient_id)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  labs(x = NULL, y = "29/36 M02 coverage-sensitivity score", color = "Patient", title = "Matched 29/36 coverage sensitivity\n(primary 36/36 score not estimable)") +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(figures_dir, "FIGURE10C3_1_matched_pairs.png"), p1, width = 5.2, height = 4.2, dpi = 300)
ggsave(file.path(figures_dir, "FIGURE10C3_1_matched_pairs.pdf"), p1, width = 5.2, height = 4.2)

fig2 <- copy(display_scores)
fig2[, display_group := fifelse(biological_sample_type == "normal mucosa", "Normal", lesion_morphology)]
fig2[, display_group := factor(display_group, levels = c("Normal", "LST-G", "protruded adenoma"))]
fwrite(fig2, file.path(figures_dir, "FIGURE10C3_2_source_data.tsv"), sep = "\t", na = "NA")
p2 <- ggplot(fig2, aes(display_group, score, color = patient_id, label = sample_id)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.8) +
  geom_text(position = position_jitter(width = 0.08, height = 0), vjust = -0.8, size = 2.7, show.legend = FALSE) +
  labs(x = NULL, y = "29/36 M02 coverage-sensitivity score", color = "Patient", title = "Patient-sample descriptive scores (sensitivity only)") +
  theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
ggsave(file.path(figures_dir, "FIGURE10C3_2_patient_scores.png"), p2, width = 6.5, height = 4.5, dpi = 300)
ggsave(file.path(figures_dir, "FIGURE10C3_2_patient_scores.pdf"), p2, width = 6.5, height = 4.5)

fig3 <- effects[record_type == "individual_patient_effect" & status == "ESTIMABLE"]
fig3[, display_method := paste(scenario_id, score_method, compartment, sep = " | ")]
fwrite(fig3, file.path(figures_dir, "FIGURE10C3_S1_source_data.tsv"), sep = "\t", na = "NA")
p3 <- ggplot(fig3, aes(effect_lesion_minus_normal, display_method, color = patient_id)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") + geom_point(size = 2) +
  labs(x = "Lesion minus normal effect", y = NULL, color = "Patient", title = "Prespecified sensitivity effects") +
  theme_classic(base_size = 9) + theme(legend.position = "top")
height3 <- max(5, 0.20 * uniqueN(fig3$display_method))
ggsave(file.path(figures_dir, "FIGURE10C3_S1_sensitivity.png"), p3, width = 8.5, height = height3, dpi = 300, limitsize = FALSE)
ggsave(file.path(figures_dir, "FIGURE10C3_S1_sensitivity.pdf"), p3, width = 8.5, height = height3, limitsize = FALSE)

primary_text <- if (primary_estimable) paste(sprintf("- %s: effect `%0.6f`", primary_indiv$patient_id, primary_indiv$effect_lesion_minus_normal), collapse = "\n") else "- P1: `NOT_ESTIMABLE` (33/36 genes mapped)\n- P5: `NOT_ESTIMABLE` (33/36 genes mapped)"
sens35_text <- paste(sprintf("- %s: effect `%0.6f`", sens35_indiv$patient_id, sens35_indiv$effect_lesion_minus_normal), collapse = "\n")

writeLines(paste0(
  "# Stage 10C3 Claim Limits\n\nThis stage is an LST-G directional sensitivity analysis only. It is not confirmatory validation, sporadic adenoma validation, generalization, mechanism analysis or biomarker research.\n\n",
  "The patient is the inferential unit. Cells, samples and libraries are nested. P3_L and P3_P are one patient's two lesions.\n\n",
  "FAP status and germline APC status are unknown. The two matched donors cannot establish a population effect.\n\n",
  "Because only 33/36 primary genes mapped, no 29/36 or 35-gene sensitivity may be relabeled as the primary result.\n\n",
  "The strongest wording is permitted only after `LST_DIRECTIONAL_CONCORDANCE`; it is not permitted under the present `NOT_ESTIMABLE` decision.\n\n",
  "No result permits changing M02, restarting dataset search, claiming clinical performance or inferring causality.\n"
), file.path(results_dir, "STAGE10C3_CLAIM_LIMITS.md"), useBytes = TRUE)

writeLines(paste0(
  "# Stage 10C3 Decision\n\nDate: ", format(Sys.Date()), "\n\nDecision: **", decision, "**\n\nRole: LST-G directional sensitivity analysis only.\n\n",
  "## Primary matched effects\n\n", primary_text, "\n\nThe primary 36/36 score has no `k/2`, confidence interval, P value or FDR because three locked genes are absent.\n\n",
  "## Prespecified sensitivities (not substitutes for primary)\n\nThe 29/36 coverage sensitivity was positive in `", sensitivity[scenario_id == "PRIMARY_ALL_EPI_MIN50" & score_method == "M02_29OF36_COVERAGE_SENS", k_positive], "/2` donors.\n\n",
  "Frozen 35-gene contamination sensitivity:\n\n", sens35_text, "\n\nNo complete prespecified sensitivity produced two systematically opposite patient effects: `", !systemic_opposite, "`.\n\n",
  "Primary bundle SHA256: `", expected_primary_bundle, "`\n\n35-gene sensitivity bundle SHA256: `", expected_sens_bundle, "`\n\n",
  "## Interpretation\n\nThe locked primary M02 score is not estimable in this matrix under the prespecified 36/36 requirement. Positive sensitivity behavior is reported descriptively and cannot be called concordant primary validation.\n\n",
  "The next separately authorized stage is Stage 10C2-SP. This file does not start it.\n"
), file.path(results_dir, "STAGE10C3_DECISION.md"), useBytes = TRUE)

writeLines(capture.output(sessionInfo()), file.path(results_dir, "STAGE10C3_SESSIONINFO.txt"), useBytes = TRUE)
message("Stage 10C3 finalization complete: ", decision)
