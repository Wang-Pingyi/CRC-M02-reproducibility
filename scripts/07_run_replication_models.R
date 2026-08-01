#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)
get_param <- read_stage7_parameters(project_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

locked <- read_locked_stage7_modules(project_dir)
candidates <- locked$candidates
membership <- locked$membership
min_fraction <- get_param("global", "module_min_gene_fraction", TRUE)
min_genes <- as.integer(get_param("global", "module_min_genes", TRUE))
min_cells <- as.integer(get_param("global", "min_cells_pseudobulk", TRUE))
min_umi <- get_param("global", "min_umi_pseudobulk", TRUE)
fdr_threshold <- get_param("global", "fdr_threshold", TRUE)
strong_sd <- get_param("global", "strong_effect_sd", TRUE)
lodo_min <- get_param("global", "lodo_min_sign_stability", TRUE)

read_pb <- function(name) {
  readRDS(file.path(paths$processed, paste0(name, "_stage7_pseudobulk_raw_counts.rds")))
}

g161 <- read_pb("GSE161277")
g132 <- read_pb("GSE132465")

prepare_eligible <- function(x) {
  eligible_pseudobulk(x$counts, x$metadata, min_cells, min_umi)
}
g161 <- prepare_eligible(g161)
g132 <- prepare_eligible(g132)

matched_161 <- c("Patient1", "Patient2", "Patient3")
g161_primary_meta <- g161$metadata[
  g161$metadata$major_cell_type == "Epithelial" &
    g161$metadata$donor_id %in% matched_161 &
    g161$metadata$condition %in% c("normal", "adenoma", "cancer"),
  ,
  drop = FALSE
]
g161_primary_counts <- g161$counts[
  , g161_primary_meta$pseudobulk_id, drop = FALSE
]
g161_scored <- module_scores_from_counts(
  g161_primary_counts, membership, min_fraction, min_genes
)

contrasts_161 <- list(
  normal_to_adenoma = c("normal", "adenoma"),
  adenoma_to_cancer = c("adenoma", "cancer"),
  normal_to_cancer = c("normal", "cancer")
)
effects_161 <- list()
details_161 <- list()
for (contrast_name in names(contrasts_161)) {
  conditions <- contrasts_161[[contrast_name]]
  for (module_id in candidates$module_id) {
    fit <- paired_effect(
      g161_scored$scores, g161_primary_meta, module_id,
      conditions[[1L]], conditions[[2L]], matched_161
    )
    fit$summary$cohort <- "GSE161277"
    fit$summary$contrast <- contrast_name
    fit$details$cohort <- "GSE161277"
    fit$details$contrast <- contrast_name
    effects_161[[paste(contrast_name, module_id)]] <- fit$summary
    details_161[[paste(contrast_name, module_id)]] <- fit$details
  }
}
effects_161 <- do.call(rbind, effects_161)
details_161 <- do.call(rbind, details_161)

matched_132 <- sprintf("SMC%02d", 1:10)
g132_epi_meta <- g132$metadata[
  g132$metadata$major_cell_type == "Epithelial" &
    g132$metadata$condition %in% c("normal", "cancer"),
  ,
  drop = FALSE
]
g132_epi_counts <- g132$counts[, g132_epi_meta$pseudobulk_id, drop = FALSE]
g132_epi_scored <- module_scores_from_counts(
  g132_epi_counts, membership, min_fraction, min_genes
)
effects_132 <- list()
details_132 <- list()
for (module_id in candidates$module_id) {
  fit <- paired_effect(
    g132_epi_scored$scores, g132_epi_meta, module_id,
    "normal", "cancer", matched_132
  )
  fit$summary$cohort <- "GSE132465"
  fit$summary$contrast <- "normal_to_cancer"
  fit$details$cohort <- "GSE132465"
  fit$details$contrast <- "normal_to_cancer"
  effects_132[[module_id]] <- fit$summary
  details_132[[module_id]] <- fit$details
}
effects_132 <- do.call(rbind, effects_132)
details_132 <- do.call(rbind, details_132)

g132_all_scored <- module_scores_from_counts(
  g132$counts, membership, min_fraction, min_genes
)
specificity <- list()
specificity_details <- list()
for (module_id in candidates$module_id) {
  fit <- specificity_effect(
    g132_all_scored$scores, g132$metadata, module_id
  )
  fit$summary$cohort <- "GSE132465"
  fit$summary$contrast <- "epithelial_specificity"
  fit$details$cohort <- "GSE132465"
  fit$details$contrast <- "epithelial_specificity"
  specificity[[module_id]] <- fit$summary
  specificity_details[[module_id]] <- fit$details
}
specificity <- do.call(rbind, specificity)
specificity_details <- do.call(rbind, specificity_details)

effects <- add_fdr_by_contrast(rbind(effects_161, effects_132, specificity))
effects$expected_discovery_effect <- NA_real_
effects$expected_direction <- NA_integer_
for (i in seq_len(nrow(effects))) {
  module_row <- candidates[candidates$module_id == effects$module_id[i], ]
  expected <- switch(
    effects$contrast[i],
    normal_to_adenoma = module_row$early_effect,
    adenoma_to_cancer = module_row$cancer_vs_adenoma_effect,
    normal_to_cancer = module_row$cancer_vs_normal_effect,
    epithelial_specificity = 1
  )
  effects$expected_discovery_effect[i] <- expected
  effects$expected_direction[i] <- sign(expected)
}
effects$direction_concordant <- sign(effects$effect) == effects$expected_direction
effects$strong_opposite <- !effects$direction_concordant &
  effects$fdr < fdr_threshold &
  abs(effects$standardized_effect) >= strong_sd
effects$donor_robust <- effects$lodo_sign_stability >= lodo_min
effects$donor_robust[
  effects$cohort == "GSE161277" & effects$n_donors == 3L
] <- effects$lodo_sign_stability[
  effects$cohort == "GSE161277" & effects$n_donors == 3L
] == 1

summary_rows <- lapply(candidates$module_id, function(module_id) {
  m <- effects[effects$module_id == module_id, ]
  early <- m[m$cohort == "GSE161277" & m$contrast == "normal_to_adenoma", ]
  later161 <- m[m$cohort == "GSE161277" & m$contrast == "normal_to_cancer", ]
  later132 <- m[m$cohort == "GSE132465" & m$contrast == "normal_to_cancer", ]
  spec <- m[m$contrast == "epithelial_specificity", ]
  any_strong_opposite <- any(m$strong_opposite %in% TRUE, na.rm = TRUE)
  early_consistent <- isTRUE(early$direction_concordant)
  cancer_consistent <- isTRUE(later132$direction_concordant)
  at_least_one_cohort <- early_consistent || cancer_consistent ||
    isTRUE(later161$direction_concordant)
  relevant_significant <- any(
    c(early$fdr, later161$fdr, later132$fdr) < fdr_threshold, na.rm = TRUE
  )
  donor_robust <- isTRUE(early$donor_robust) || isTRUE(later132$donor_robust)
  specificity_positive <- isTRUE(spec$effect > 0)
  specificity_robust <- isTRUE(spec$donor_robust)

  classification <- if (any_strong_opposite) {
    "contradictory"
  } else if (
    early_consistent && cancer_consistent && relevant_significant &&
      donor_robust && specificity_positive && specificity_robust
  ) {
    "replicated"
  } else if (at_least_one_cohort && donor_robust && specificity_positive) {
    "directionally consistent but underpowered"
  } else {
    "not replicated"
  }
  data.frame(
    module_id = module_id,
    discovery_epithelial_state =
      candidates$epithelial_state[candidates$module_id == module_id],
    evidence_hierarchy = "secondary_exploratory_locked_module",
    GSE161277_normal_to_adenoma_effect = early$effect,
    GSE161277_normal_to_adenoma_ci_low = early$ci_low,
    GSE161277_normal_to_adenoma_ci_high = early$ci_high,
    GSE161277_normal_to_adenoma_fdr = early$fdr,
    GSE161277_normal_to_adenoma_direction_concordant = early_consistent,
    GSE161277_normal_to_adenoma_lodo = early$lodo_sign_stability,
    GSE132465_normal_to_cancer_effect = later132$effect,
    GSE132465_normal_to_cancer_ci_low = later132$ci_low,
    GSE132465_normal_to_cancer_ci_high = later132$ci_high,
    GSE132465_normal_to_cancer_fdr = later132$fdr,
    GSE132465_normal_to_cancer_direction_concordant = cancer_consistent,
    GSE132465_normal_to_cancer_lodo = later132$lodo_sign_stability,
    epithelial_specificity_effect = spec$effect,
    epithelial_specificity_ci_low = spec$ci_low,
    epithelial_specificity_ci_high = spec$ci_high,
    epithelial_specificity_fdr = spec$fdr,
    epithelial_specificity_positive = specificity_positive,
    epithelial_specificity_lodo = spec$lodo_sign_stability,
    strong_opposite_effect = any_strong_opposite,
    donor_robust = donor_robust,
    replication_class = classification,
    validation_gene_reselection = FALSE,
    stringsAsFactors = FALSE
  )
})
replication_summary <- do.call(rbind, summary_rows)

gate_eligible <- replication_summary$replication_class %in% c(
  "replicated", "directionally consistent but underpowered"
) &
  !replication_summary$strong_opposite_effect &
  replication_summary$donor_robust &
  replication_summary$epithelial_specificity_positive
gate_pass <- any(gate_eligible)
gate <- data.frame(
  criterion = c(
    "at_least_one_independent_cohort_direction_consistent",
    "no_strong_opposite_for_gate_eligible_module",
    "not_single_donor_driven",
    "epithelial_specificity_positive",
    "overall_replication_gate"
  ),
  passed = c(
    any(replication_summary$replication_class %in% c(
      "replicated", "directionally consistent but underpowered"
    )),
    any(gate_eligible),
    any(gate_eligible),
    any(gate_eligible),
    gate_pass
  ),
  note = c(
    "Evaluated without changing locked module membership",
    "Strong opposite is FDR<0.10 and |standardized effect|>=0.5",
    "LODO>=0.75; three-donor GSE161277 promotion requires 1.00",
    "Positive epithelial-minus-non-epithelial donor-level effect",
    if (gate_pass) {
      "Gate passed for at least one exploratory module; remain module-level"
    } else {
      "Gate failed; return to annotation/model audit and do not add machine learning"
    }
  ),
  stringsAsFactors = FALSE
)

write_stage7_tsv(
  replication_summary, file.path(paths$result, "replication_summary.tsv")
)
write_stage7_tsv(
  effects, file.path(paths$result, "replication_effects_all.tsv")
)
write_stage7_tsv(
  gate, file.path(paths$result, "replication_gate.tsv")
)
write_stage7_tsv(
  rbind(details_161, details_132),
  file.path(paths$source, "paired_donor_module_differences.tsv")
)
write_stage7_tsv(
  specificity_details,
  file.path(paths$source, "epithelial_specificity_donor_differences.tsv")
)
write_stage7_tsv(
  g161_scored$scores,
  file.path(paths$source, "GSE161277_module_scores.tsv")
)
write_stage7_tsv(
  g132_epi_scored$scores,
  file.path(paths$source, "GSE132465_epithelial_module_scores.tsv")
)
write_stage7_tsv(
  g132_all_scored$scores,
  file.path(paths$source, "GSE132465_all_celltype_module_scores.tsv")
)
write_stage7_tsv(
  g161_scored$normalization,
  file.path(paths$source, "GSE161277_pseudobulk_normalization.tsv")
)
write_stage7_tsv(
  g132_all_scored$normalization,
  file.path(paths$source, "GSE132465_pseudobulk_normalization.tsv")
)

plot_data <- effects[
  effects$contrast %in% c(
    "normal_to_adenoma", "normal_to_cancer", "epithelial_specificity"
  ),
]
plot_data$panel <- paste(plot_data$cohort, plot_data$contrast, sep = ": ")
p <- ggplot(plot_data, aes(effect, module_id, color = direction_concordant)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15) +
  geom_point(size = 2) +
  facet_wrap(~ panel, scales = "free_x") +
  scale_color_manual(values = c(`TRUE` = "#1B7837", `FALSE` = "#B2182B")) +
  labs(
    title = "Locked module effects in independent single-cell cohorts",
    x = "Donor-level module-score effect (95% CI)",
    y = NULL, color = "Direction agrees"
  ) +
  theme_bw(base_size = 9)
ggsave(
  file.path(paths$figure, "stage_7_module_effects.pdf"),
  p, width = 11, height = 7
)
ggsave(
  file.path(paths$figure, "stage_7_module_effects.png"),
  p, width = 11, height = 7, dpi = 300
)
write_stage7_tsv(
  plot_data, file.path(paths$source, "stage_7_module_effects_source.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(paths$result, "stage_7_model_sessionInfo.txt")
)
cat(
  "STAGE_7_MODELS_OK\tgate=", gate_pass,
  "\treplicated=", sum(replication_summary$replication_class == "replicated"),
  "\n", sep = ""
)
