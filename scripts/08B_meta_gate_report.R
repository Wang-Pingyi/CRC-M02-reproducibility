#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08B_meta_gate_report.R PROJECT_DIR RUN_ID")
project <- normalizePath(args[1], mustWork = TRUE)
run_id <- args[2]
source(file.path(project, "scripts", "08B_helpers.R"))
set.seed(20260729)
paths <- stage8b_paths(project, run_id)
if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")

effects <- read.delim(file.path(paths$result, "bulk_cohort_effects.tsv"), check.names = FALSE)
tcga <- read.delim(file.path(paths$result, "TCGA_auxiliary_results.tsv"), check.names = FALSE)
candidates <- stage8b_locked(paths)$candidates
effects$evidence_role <- "independent_tissue_validation"
tcga$evidence_role <- "TCGA_auxiliary_only"
summary_all <- rbind(effects, tcga[, names(effects)])
stage8b_write_tsv(summary_all, file.path(paths$result, "bulk_validation_summary.tsv"))

comparable <- effects[effects$endpoint %in% c("adenoma_vs_normal", "cancer_vs_adenoma",
                                               "cancer_vs_normal", "ordered_trend"), ]
meta_rows <- list()
meta_input <- list()
for (set_name in unique(comparable$analysis_set)) {
  for (endpoint in unique(comparable$endpoint)) {
    for (module in unique(comparable$module_id)) {
      z <- comparable[comparable$analysis_set == set_name & comparable$endpoint == endpoint &
                        comparable$module_id == module & is.finite(comparable$effect) &
                        is.finite(comparable$standard_error) & comparable$standard_error > 0, ]
      if (nrow(z) < 2L) next
      fit <- metafor::rma.uni(yi = z$effect, sei = z$standard_error, method = "REML", test = "knha")
      pred <- try(metafor::predict.rma(fit), silent = TRUE)
      pi_lb <- if (inherits(pred, "try-error")) NA_real_ else as.numeric(pred$pi.lb)
      pi_ub <- if (inherits(pred, "try-error")) NA_real_ else as.numeric(pred$pi.ub)
      meta_rows[[length(meta_rows) + 1L]] <- data.frame(
        module_id = module, endpoint = endpoint, analysis_set = set_name,
        k_cohorts = nrow(z), pooled_effect = as.numeric(fit$b[1]),
        standard_error = fit$se, ci_low = fit$ci.lb, ci_high = fit$ci.ub,
        p_value = fit$pval, prediction_low = pi_lb, prediction_high = pi_ub,
        tau2 = fit$tau2, I2 = fit$I2, Q = fit$QE, Q_p_value = fit$QEp,
        method = "REML_Knapp-Hartung", caution = ifelse(nrow(z) == 2L, "k=2; heterogeneity imprecise", "NA"),
        stringsAsFactors = FALSE
      )
      z$meta_analysis_set <- set_name
      meta_input[[length(meta_input) + 1L]] <- z
    }
  }
}
meta <- do.call(rbind, meta_rows)
meta$fdr <- NA_real_
key <- interaction(meta$endpoint, meta$analysis_set, drop = TRUE)
for (k in levels(key)) meta$fdr[key == k] <- p.adjust(meta$p_value[key == k], "BH")
stage8b_write_tsv(meta, file.path(paths$result, "meta_analysis_results.tsv"))
meta_source <- do.call(rbind, meta_input)
stage8b_write_tsv(meta_source, file.path(paths$result, "meta_analysis_source_data.tsv"))

# Prespecified tissue gate on the early normal-to-adenoma endpoint.
early <- effects[effects$analysis_set == "primary_all_samples" & effects$endpoint == "adenoma_vs_normal", ]
early <- merge(early, candidates[, c("module_id", "early_effect")], by = "module_id", all.x = TRUE)
early$concordant <- sign(early$effect) == sign(early$early_effect)
early$strong_opposite <- !early$concordant & early$fdr < 0.05 &
  ((early$ci_low > 0 & early$early_effect < 0) | (early$ci_high < 0 & early$early_effect > 0))
early_meta <- meta[meta$analysis_set == "primary_all_samples" & meta$endpoint == "adenoma_vs_normal",
                   c("module_id", "pooled_effect", "ci_low", "ci_high", "p_value", "fdr")]
names(early_meta)[-1] <- paste0("meta_", names(early_meta)[-1])
gate <- do.call(rbind, lapply(split(early, early$module_id), function(z) {
  m <- early_meta[early_meta$module_id == z$module_id[1], , drop = FALSE]
  n_concordant <- sum(z$concordant)
  strong_opposite <- any(z$strong_opposite)
  significant <- (!is.na(m$meta_fdr) && m$meta_fdr < 0.05) || any(z$fdr < 0.05 & z$concordant)
  classification <- if (strong_opposite) "contradictory" else if (n_concordant >= 2 && significant) {
    "replicated"
  } else if (n_concordant >= 2) {
    "directionally_consistent_but_underpowered"
  } else {
    "not_replicated"
  }
  data.frame(module_id = z$module_id[1], expected_early_direction = sign(z$early_effect[1]),
    concordant_cohorts = n_concordant, evaluated_cohorts = nrow(z),
    strong_opposite = strong_opposite, classification = classification,
    meta_effect = if (nrow(m)) m$meta_pooled_effect else NA_real_,
    meta_ci_low = if (nrow(m)) m$meta_ci_low else NA_real_,
    meta_ci_high = if (nrow(m)) m$meta_ci_high else NA_real_,
    meta_fdr = if (nrow(m)) m$meta_fdr else NA_real_, stringsAsFactors = FALSE)
}))
advancing <- gate$classification %in% c("replicated", "directionally_consistent_but_underpowered")
stage_pass <- any(advancing) && !any(gate$strong_opposite[advancing])
gate$stage_gate <- ifelse(stage_pass, "PASS", "FAIL")
stage8b_write_tsv(gate, file.path(paths$result, "tissue_validation_gate.tsv"))

# Forest plot and exact source data.
plot_data <- early
pooled <- meta[meta$analysis_set == "primary_all_samples" & meta$endpoint == "adenoma_vs_normal", ]
pooled_plot <- data.frame(
  module_id = pooled$module_id, accession = "Random-effects meta",
  effect = pooled$pooled_effect, ci_low = pooled$ci_low, ci_high = pooled$ci_high,
  stringsAsFactors = FALSE
)
plot_data <- rbind(plot_data[, c("module_id", "accession", "effect", "ci_low", "ci_high")], pooled_plot)
plot_data$accession <- factor(plot_data$accession,
  levels = c("GSE41657", "GSE100179", "GSE8671", "Random-effects meta"))
stage8b_write_tsv(plot_data, file.path(paths$figure, "stage_8B_early_transition_forest_source_data.tsv"))
p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect, y = accession, color = accession)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high), height = 0.16) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~module_id, scales = "free_x", ncol = 2) +
  ggplot2::labs(x = "Standardized normal-to-adenoma effect (95% CI)", y = NULL,
                title = "Locked epithelial modules in independent tissue cohorts") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(paths$figure, "stage_8B_early_transition_forest.png"), p,
                width = 10, height = 8, dpi = 300)
ggplot2::ggsave(file.path(paths$figure, "stage_8B_early_transition_forest.pdf"), p,
                width = 10, height = 8)

# Small canonical copies; large/intermediate objects remain outside Git.
file.copy(file.path(paths$result, "bulk_validation_summary.tsv"),
          file.path(project, "results_final", "stage_8B_bulk_validation_summary.tsv"), overwrite = TRUE)
file.copy(file.path(paths$result, "meta_analysis_results.tsv"),
          file.path(project, "results_final", "stage_8B_meta_analysis_results.tsv"), overwrite = TRUE)
file.copy(file.path(paths$result, "tissue_validation_gate.tsv"),
          file.path(project, "logs_summary", "stage_8B_tissue_validation_gate.tsv"), overwrite = TRUE)
file.copy(file.path(paths$figure, "stage_8B_early_transition_forest_source_data.tsv"),
          file.path(project, "figures_final", "stage_8B_early_transition_forest_source_data.tsv"), overwrite = TRUE)
file.copy(file.path(paths$figure, "stage_8B_early_transition_forest.png"),
          file.path(project, "figures_final", "stage_8B_early_transition_forest.png"), overwrite = TRUE)
file.copy(file.path(paths$figure, "stage_8B_early_transition_forest.pdf"),
          file.path(project, "figures_final", "stage_8B_early_transition_forest.pdf"), overwrite = TRUE)

class_counts <- as.data.frame(table(gate$classification), stringsAsFactors = FALSE)
class_text <- paste(paste0(class_counts$Freq, " ", class_counts$Var1), collapse = "; ")
report <- c(
  "# Stage 8B tissue validation and cross-cohort Meta-analysis",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Server-computed tissue gate: **", ifelse(stage_pass, "PASS", "FAIL"), "**"),
  paste0("- Module classifications: ", class_text),
  "- Evidence hierarchy: secondary/exploratory; the frozen primary Stage 6A result remains negative.",
  "",
  "## Methods",
  "",
  "GSE41657, GSE100179 and GSE8671 were analyzed independently. Module scores",
  "were constructed from frozen genes after within-cohort probe collapse and",
  "gene standardization. Unpaired models used HC3 robust standard errors;",
  "GSE8671 used one normal-to-adenoma difference per verified donor pair.",
  "The 32 GSE8671 pairs use the explicit patient numbers in GEO sample titles;",
  "non-unique two-letter patient initials are retained only as provenance.",
  "Comparable standardized effects were synthesized using REML random-effects",
  "models with Knapp-Hartung intervals. Cross-platform matrices were not merged.",
  "",
  "TCGA-COAD was analyzed separately as auxiliary cancer-normal, stage, MSI, CMS",
  "and clinical evidence and did not contribute to the adenoma-sequence Meta-analysis.",
  "",
  "## Gate interpretation",
  "",
  "The server-side classification is a prespecified computational gate and awaits",
  "independent Codex quality-control review. It must not be reported as accepted",
  "until input integrity, model coefficients, sensitivities and source tables are audited.",
  "",
  "## Output tables",
  "",
  paste0("- `results/08B_bulk_validation/", run_id, "/bulk_validation_summary.tsv`"),
  paste0("- `results/08B_bulk_validation/", run_id, "/meta_analysis_results.tsv`"),
  paste0("- `results/08B_bulk_validation/", run_id, "/tissue_validation_gate.tsv`"),
  paste0("- `results/08B_bulk_validation/", run_id, "/TCGA_auxiliary_results.tsv`"),
  "",
  "## Stage boundary",
  "",
  "Stage 8B stopped before stool RNA validation. Stage 9 is not authorized."
)
writeLines(report, file.path(project, "reports", "stage_8B_bulk_validation.md"), useBytes = TRUE)

status_path <- file.path(project, "STATUS.md")
status <- readLines(status_path, warn = FALSE)
marker <- paste0("- Stage 8B run ID: `", run_id, "`")
if (!marker %in% status) {
  status <- c(status, "", "- Stage 8B authorization: granted 2026-07-29",
    "- Stage 8B scope: independent tissue-cohort trend models, random-effects Meta-analysis and auxiliary TCGA-COAD analyses",
    marker,
    paste0("- Stage 8B server-computed gate pending independent Codex QC: ", ifelse(stage_pass, "PASS", "FAIL")),
    "- Stage 8B completion marker: `logs/08B_bulk_validation/READY_FOR_CODEX_QC`",
    "- Stage 8B stage boundary: stopped before Stage 9")
  writeLines(status, status_path, useBytes = TRUE)
}
