#!/usr/bin/env Rscript

# Stage 10I independent red-team numeric reconciliation
# Date: 2026-08-01
# Random seed: 42
# Execution: Rscript --vanilla scripts/10I_red_team_recompute.R <project_root>
# This script is read-only with respect to all historical Stage 6-10 outputs.

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
out <- file.path(root, "audit", "stage10i", "recompute")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "\"", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("NA", "NaN", ""))
}

write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

clean_num <- function(x) as.numeric(as.character(x))

# 1. Stage 6A prespecified primary normal-to-adenoma results.
pb_path <- file.path(root, "results", "06A_pseudobulk", "pseudobulk_results.tsv")
pb <- read_tsv(pb_path)
primary <- pb[pb$contrast == "adenoma_vs_normal" &
                (is.na(pb$omitted_donor) | pb$omitted_donor == ""), , drop = FALSE]
primary$FDR <- clean_num(primary$FDR)
state_split <- split(primary, primary$epithelial_state)
stage6 <- do.call(rbind, lapply(names(state_split), function(state) {
  z <- state_split[[state]]
  data.frame(
    epithelial_state = state,
    tested_genes = nrow(z),
    fdr_below_0_05 = sum(z$FDR < 0.05, na.rm = TRUE),
    minimum_fdr = min(z$FDR, na.rm = TRUE),
    n_normal_donors = paste(sort(unique(z$n_normal_donors)), collapse = ";"),
    n_adenoma_donors = paste(sort(unique(z$n_adenoma_donors)), collapse = ";"),
    n_cancer_donors = paste(sort(unique(z$n_cancer_donors)), collapse = ";"),
    model_formula = paste(sort(unique(z$model_formula)), collapse = ";"),
    inferential_unit = "donor_pseudobulk",
    stringsAsFactors = FALSE
  )
}))
stage6$primary_discovery_status <- ifelse(stage6$fdr_below_0_05 == 0,
                                          "NEGATIVE_NO_FDR_LT_0.05", "NONZERO_HITS")
write_tsv(stage6, "STAGE6A_PRIMARY_RECOMPUTE.tsv")

# 2. Stage 10C internal four-patient FAP module results.
fap_path <- file.path(root, "results", "10C_fap_confounding", "20260731_110021",
                      "GSE201348_module_donor_effects.tsv")
fap <- read_tsv(fap_path)
fap$difference <- clean_num(fap$difference)
fap_groups <- split(fap, interaction(fap$scope, fap$module_id, drop = TRUE))
fap_re <- do.call(rbind, lapply(fap_groups, function(z) {
  tt <- t.test(z$difference, mu = 0)
  lodo <- vapply(seq_len(nrow(z)), function(i) mean(z$difference[-i]), numeric(1))
  data.frame(
    dataset = unique(z$dataset), scope = unique(z$scope), module_id = unique(z$module_id),
    n_donors = nrow(z), donor_ids = paste(z$donor_id, collapse = ";"),
    effect = mean(z$difference), standard_error = sd(z$difference) / sqrt(nrow(z)),
    ci_low = unname(tt$conf.int[1]), ci_high = unname(tt$conf.int[2]),
    p_value = unname(tt$p.value), lodo_positive_fraction = mean(lodo > 0),
    min_lodo_effect = min(lodo), max_lodo_effect = max(lodo),
    inferential_unit = "donor", stringsAsFactors = FALSE
  )
}))
fap_re$fdr <- ave(fap_re$p_value, fap_re$scope,
                  FUN = function(x) p.adjust(x, method = "BH"))
fap_m02 <- fap_re[fap_re$module_id == "Stem_progenitor_SB_M02", , drop = FALSE]
write_tsv(fap_re, "STAGE10C_FAP_ALL_MODULES_RECOMPUTE.tsv")
write_tsv(fap_m02, "STAGE10C_FAP_M02_RECOMPUTE.tsv")

# 3. GSE161277 three-patient directional evidence, counted once.
g161_path <- file.path(root, "results", "10C_fap_confounding", "20260731_110021",
                       "GSE161277_three_paired_donor_effects.tsv")
g161 <- read_tsv(g161_path)
g161$difference <- clean_num(g161$difference)
g161_groups <- split(g161, g161$module_id)
g161_re <- do.call(rbind, lapply(g161_groups, function(z) {
  tt <- t.test(z$difference, mu = 0)
  lodo <- vapply(seq_len(nrow(z)), function(i) mean(z$difference[-i]), numeric(1))
  data.frame(
    dataset = "GSE161277", module_id = unique(z$module_id), n_donors = nrow(z),
    donor_ids = paste(z$donor_id, collapse = ";"), effect = mean(z$difference),
    standard_error = sd(z$difference) / sqrt(nrow(z)),
    ci_low = unname(tt$conf.int[1]), ci_high = unname(tt$conf.int[2]),
    p_value = unname(tt$p.value), k_positive = sum(z$difference > 0),
    lodo_positive_fraction = mean(lodo > 0), inferential_unit = "donor",
    stringsAsFactors = FALSE
  )
}))
g161_re$fdr <- p.adjust(g161_re$p_value, method = "BH")
g161_m02 <- g161_re[g161_re$module_id == "Stem_progenitor_SB_M02", , drop = FALSE]
write_tsv(g161_re, "GSE161277_ALL_MODULES_RECOMPUTE.tsv")
write_tsv(g161_m02, "GSE161277_M02_RECOMPUTE.tsv")

# 4. Three-accession bulk random-effects meta-analysis, using independent cohort effects.
if (!requireNamespace("metafor", quietly = TRUE)) {
  stop("Package 'metafor' is required for the independent bulk recomputation")
}
bulk_path <- file.path(root, "results_final", "stage_8B_bulk_validation_summary.tsv")
bulk <- read_tsv(bulk_path)
bulk_primary <- bulk[bulk$endpoint == "adenoma_vs_normal" &
                       bulk$analysis_set == "primary_all_samples" &
                       bulk$accession %in% c("GSE41657", "GSE100179", "GSE8671"), , drop = FALSE]
bulk_primary$effect <- clean_num(bulk_primary$effect)
bulk_primary$standard_error <- clean_num(bulk_primary$standard_error)
bulk_groups <- split(bulk_primary, bulk_primary$module_id)
meta_re <- do.call(rbind, lapply(bulk_groups, function(z) {
  z <- z[order(z$accession), , drop = FALSE]
  fit <- metafor::rma.uni(yi = z$effect, sei = z$standard_error,
                          method = "REML", test = "knha")
  data.frame(
    module_id = unique(z$module_id), endpoint = "adenoma_vs_normal",
    k_cohorts = nrow(z), accessions = paste(z$accession, collapse = ";"),
    pooled_effect = as.numeric(fit$b[1]), standard_error = as.numeric(fit$se[1]),
    ci_low = as.numeric(fit$ci.lb), ci_high = as.numeric(fit$ci.ub),
    p_value = as.numeric(fit$pval), tau2 = as.numeric(fit$tau2),
    I2 = as.numeric(fit$I2), Q = as.numeric(fit$QE), Q_p_value = as.numeric(fit$QEp),
    model = "REML_Knapp-Hartung", inferential_unit = "cohort_effect",
    stringsAsFactors = FALSE
  )
}))
meta_re$fdr <- p.adjust(meta_re$p_value, method = "BH")
meta_m02 <- meta_re[meta_re$module_id == "Stem_progenitor_SB_M02", , drop = FALSE]
write_tsv(meta_re, "STAGE8B_META_ALL_MODULES_RECOMPUTE.tsv")
write_tsv(meta_m02, "STAGE8B_META_M02_RECOMPUTE.tsv")

# 5. Evidence-matrix structural and non-duplication checks.
ev_path <- file.path(root, "results", "stage10h", "STAGE10H_EVIDENCE_MATRIX.tsv")
if (!file.exists(ev_path)) {
  ev_path <- file.path(root, "audit", "stage10i", "recompute_inputs",
                       "STAGE10H_EVIDENCE_MATRIX.tsv")
}
ev <- read_tsv(ev_path)
evidence_checks <- data.frame(
  check_id = c(
    "unique_evidence_ids", "gse161277_counted_once", "figshare_counted_once",
    "stage10c3_primary_not_estimable", "case4_single_patient",
    "spatial_patient_level_not_estimable", "no_level_a_or_b_status"
  ),
  observed = c(
    as.character(!anyDuplicated(ev$evidence_id)),
    as.character(sum(grepl("GSE161277", ev$dataset, fixed = TRUE))),
    as.character(sum(grepl("Figshare_29925404", ev$dataset, fixed = TRUE))),
    paste(ev$primary_result_status[ev$stage == "Stage_10C3"], collapse = ";"),
    paste(ev$patients_n[ev$stage == "Stage_10E_DESC"], collapse = ";"),
    paste(ev$primary_result_status[ev$stage == "Stage_10F"], collapse = ";"),
    as.character(!any(grepl("LEVEL_A|LEVEL_B", ev$primary_result_status)))
  ),
  expected = c("TRUE", "1", "1", "NOT_ESTIMABLE", "1_patient",
               "SKIPPED_NOT_ESTIMABLE", "TRUE"),
  stringsAsFactors = FALSE
)
evidence_checks$pass <- evidence_checks$observed == evidence_checks$expected
write_tsv(evidence_checks, "STAGE10H_EVIDENCE_MATRIX_RECOMPUTE.tsv")

# 6. Case4 descriptive score, preserving n=1 and no inferential statistics.
case_path <- file.path(root, "results", "stage10e_desc", "STAGE10E_DESC_CASE4_SCORES.tsv")
case <- read_tsv(case_path)
normal <- clean_num(case$score[case$region == "Normal"])
adenoma <- clean_num(case$score[case$region == "Adenoma"])
case_re <- data.frame(
  patient_id = "case4", patient_n = 1L,
  normal_score = normal, adenoma_score = adenoma,
  adenoma_minus_normal = adenoma - normal,
  mapped_genes = unique(case$mapped_genes), canonical_genes = unique(case$canonical_genes),
  inferential_statistics = "NOT_COMPUTED_BY_DESIGN",
  interpretation = "single_patient_descriptive_only",
  stringsAsFactors = FALSE
)
write_tsv(case_re, "STAGE10E_DESC_CASE4_RECOMPUTE.tsv")

# 7. Reconcile and independently render one core source-data figure.
forest_path <- file.path(root, "figures_final", "stage_8B_early_transition_forest_source_data.tsv")
forest <- read_tsv(forest_path)
forest_m02 <- forest[forest$module_id == "Stem_progenitor_SB_M02", , drop = FALSE]
bulk_m02 <- bulk_primary[bulk_primary$module_id == "Stem_progenitor_SB_M02",
                         c("accession", "effect", "ci_low", "ci_high"), drop = FALSE]
core <- merge(forest_m02, bulk_m02, by = "accession", suffixes = c("_figure_source", "_bulk_table"))
for (v in c("effect", "ci_low", "ci_high")) {
  core[[paste0(v, "_absolute_difference")]] <-
    abs(clean_num(core[[paste0(v, "_figure_source")]]) -
          clean_num(core[[paste0(v, "_bulk_table")]]))
}
core$pass_1e_12 <- apply(core[grepl("absolute_difference$", names(core))], 1,
                             function(x) all(clean_num(x) <= 1e-12))
write_tsv(core, "CORE_FIGURE_SOURCE_RECOMPUTE.tsv")

plot_core <- function(device, path) {
  device(path, width = 7, height = 4.8)
  op <- par(mar = c(5, 7, 3, 1))
  on.exit({par(op); dev.off()}, add = TRUE)
  y <- seq_len(nrow(core))
  x <- clean_num(core$effect_bulk_table)
  lo <- clean_num(core$ci_low_bulk_table)
  hi <- clean_num(core$ci_high_bulk_table)
  plot(x, y, xlim = range(c(0, lo, hi)), ylim = c(0.5, nrow(core) + 0.5),
       yaxt = "n", ylab = "", xlab = "Standardized adenoma-minus-normal effect (95% CI)",
       pch = 19, main = "Stage 10I independent source-data reconciliation")
  axis(2, at = y, labels = core$accession, las = 1)
  abline(v = 0, lty = 2, col = "grey50")
  segments(lo, y, hi, y, lwd = 2)
  points(x, y, pch = 19)
  mtext("Three accessions analyzed separately; no matrix pooling", side = 3, line = 0.2, cex = 0.8)
}
plot_core(function(path, width, height) png(path, width = width, height = height,
                                            units = "in", res = 300),
          file.path(out, "CORE_FIGURE_RECOMPUTED.png"))
plot_core(function(path, width, height) pdf(path, width = width, height = height,
                                            useDingbats = FALSE),
          file.path(out, "CORE_FIGURE_RECOMPUTED.pdf"))

# 8. Compact numeric reconciliation against the frozen historical tables.
locked10c <- read_tsv(file.path(root, "results", "10C_fap_confounding", "20260731_110021",
                                "stage10C_locked_module_results.tsv"))
locked_m02 <- locked10c[locked10c$module_id == "Stem_progenitor_SB_M02" &
                         locked10c$dataset %in% c("GSE201348", "GSE161277"), , drop = FALSE]
existing_meta <- read_tsv(file.path(root, "results_final", "stage_8B_meta_analysis_results.tsv"))
existing_meta <- existing_meta[existing_meta$module_id == "Stem_progenitor_SB_M02" &
                                 existing_meta$endpoint == "adenoma_vs_normal" &
                                 existing_meta$analysis_set == "primary_all_samples", , drop = FALSE]
existing_case <- read_tsv(file.path(root, "results", "stage10e_desc", "STAGE10E_DESC_SENSITIVITY.tsv"))
existing_case <- existing_case[existing_case$sensitivity == "primary_reference", , drop = FALSE]

scalar_numeric <- function(x, label) {
  if (length(x) != 1L) stop(sprintf("Expected one historical value for %s; found %d", label, length(x)))
  clean_num(x)
}
locked_value <- function(dataset, scope, contrast, field) {
  idx <- locked_m02$dataset == dataset & locked_m02$scope == scope & locked_m02$contrast == contrast
  scalar_numeric(locked_m02[[field]][idx], paste(dataset, scope, contrast, field, sep = "/"))
}

recon <- data.frame(
  item = c(
    "Stage6A_primary_total_FDR_hits",
    "Stage10C_FAP_all_epithelial_effect", "Stage10C_FAP_all_epithelial_p",
    "Stage10C_FAP_all_epithelial_fdr", "Stage10C_FAP_stem_effect",
    "Stage10C_FAP_stem_p", "Stage10C_FAP_stem_fdr",
    "GSE161277_M02_effect", "GSE161277_M02_p", "GSE161277_M02_fdr",
    "Stage8B_M02_meta_effect", "Stage8B_M02_meta_ci_low", "Stage8B_M02_meta_ci_high",
    "Stage8B_M02_meta_p", "Stage8B_M02_meta_fdr", "Stage8B_M02_meta_I2",
    "case4_adenoma_minus_normal"
  ),
  recomputed = c(
    sum(stage6$fdr_below_0_05),
    fap_m02$effect[fap_m02$scope == "All_epithelial"],
    fap_m02$p_value[fap_m02$scope == "All_epithelial"],
    fap_m02$fdr[fap_m02$scope == "All_epithelial"],
    fap_m02$effect[fap_m02$scope == "Stem_progenitor"],
    fap_m02$p_value[fap_m02$scope == "Stem_progenitor"],
    fap_m02$fdr[fap_m02$scope == "Stem_progenitor"],
    g161_m02$effect, g161_m02$p_value, g161_m02$fdr,
    meta_m02$pooled_effect, meta_m02$ci_low, meta_m02$ci_high,
    meta_m02$p_value, meta_m02$fdr, meta_m02$I2,
    case_re$adenoma_minus_normal
  ),
  historical = c(
    0,
    locked_value("GSE201348", "All_epithelial", "FAP_adenoma_vs_normal", "effect"),
    locked_value("GSE201348", "All_epithelial", "FAP_adenoma_vs_normal", "p_value"),
    locked_value("GSE201348", "All_epithelial", "FAP_adenoma_vs_normal", "FDR"),
    locked_value("GSE201348", "Stem_progenitor", "FAP_adenoma_vs_normal", "effect"),
    locked_value("GSE201348", "Stem_progenitor", "FAP_adenoma_vs_normal", "p_value"),
    locked_value("GSE201348", "Stem_progenitor", "FAP_adenoma_vs_normal", "FDR"),
    locked_value("GSE161277", "Epithelial", "paired_normal_to_adenoma", "effect"),
    locked_value("GSE161277", "Epithelial", "paired_normal_to_adenoma", "p_value"),
    locked_value("GSE161277", "Epithelial", "paired_normal_to_adenoma", "FDR"),
    scalar_numeric(existing_meta$pooled_effect, "meta/effect"),
    scalar_numeric(existing_meta$ci_low, "meta/ci_low"),
    scalar_numeric(existing_meta$ci_high, "meta/ci_high"),
    scalar_numeric(existing_meta$p_value, "meta/p"),
    scalar_numeric(existing_meta$fdr, "meta/fdr"),
    scalar_numeric(existing_meta$I2, "meta/I2"),
    scalar_numeric(existing_case$Adenoma_minus_Normal, "case4/difference")
  ),
  tolerance = c(0, rep(1e-10, 16)),
  stringsAsFactors = FALSE
)
recon$absolute_difference <- abs(recon$recomputed - recon$historical)
recon$status <- ifelse(recon$absolute_difference <= recon$tolerance, "MATCH", "MISMATCH")
write_tsv(recon, "STAGE10I_NUMERIC_RECOMPUTE_RAW.tsv")

session_lines <- sub("[[:space:]]+$", "", capture.output(sessionInfo()))
writeLines(session_lines, file.path(out, "STAGE10I_RECOMPUTE_SESSIONINFO.txt"))
writeLines(c(
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "seed=42", "execution=Rscript --vanilla", "historical_outputs_modified=FALSE",
  paste0("all_numeric_matches=", all(recon$status == "MATCH")),
  paste0("all_evidence_structure_checks_pass=", all(evidence_checks$pass)),
  paste0("all_core_figure_source_checks_pass=", all(core$pass_1e_12))
), file.path(out, "RECOMPUTE.SUCCESS"))
