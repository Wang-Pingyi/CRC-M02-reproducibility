#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args0 == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args0)) stop("Missing value for ", flag)
  args0[hit[length(hit)] + 1L]
}

root <- normalizePath(arg_value("--root", "."), mustWork = TRUE)
resume_after_report_failure <- "--resume-after-report-failure" %in% args0
project_lib <- file.path(root, "environment", "R-library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(UCell)
  library(ggplot2)
})

set.seed(42)
options(stringsAsFactors = FALSE)

results_dir <- file.path(root, "results", "stage10c3")
figures_dir <- file.path(root, "figures", "stage10c3")
objects_dir <- file.path(root, "objects", "stage10c3")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

plan_path <- file.path(results_dir, "STAGE10C3_ANALYSIS_PLAN_LOCKED.md")
qc_success <- file.path(results_dir, "QC_PREM02.SUCCESS")
decision_b <- file.path(root, "results", "stage10c2_sc", "STAGE10C2_SC_B_DECISION.md")
lock_path <- file.path(root, "results", "stage10c", "STAGE10C_LOCK_MANIFEST.tsv")
sens_lock_path <- file.path(root, "results", "stage10c", "M02_MINUS_INPP5D_SENS_V1.tsv")
manifest_path <- file.path(root, "data", "metadata", "stage10c2_sc_patient_sample_manifest.tsv")

required <- c(plan_path, qc_success, decision_b, lock_path, sens_lock_path, manifest_path)
if (!all(file.exists(required))) stop("Missing required frozen input: ", paste(required[!file.exists(required)], collapse = ", "))
if (!grepl("PASS_AS_DIRECTIONAL_SENSITIVITY_ONLY", paste(readLines(decision_b, warn = FALSE), collapse = "\n"), fixed = TRUE)) stop("10C2-SC-B decision gate failed")
if (file.mtime(plan_path) > file.mtime(qc_success)) stop("Analysis plan does not predate QC completion")

final_score_path <- file.path(results_dir, "STAGE10C3_PATIENT_SAMPLE_SCORES.tsv")
if (file.exists(final_score_path)) {
  failed_marker <- file.path(root, "status", "STAGE10C3.FAILED")
  failed_history <- list.files(file.path(root, "status", "history"), pattern = "^STAGE10C3\\..*\\.FAILED$", full.names = TRUE)
  decision_path <- file.path(results_dir, "STAGE10C3_DECISION.md")
  allowed_resume <- resume_after_report_failure && (file.exists(failed_marker) || length(failed_history) > 0L) && !file.exists(decision_path)
  if (!allowed_resume) stop("Existing Stage 10C3 score output found; refusing silent overwrite")
  message("Explicitly resuming after report-generation failure; frozen QC and annotation objects are reused")
}

lock <- fread(lock_path)
sens_lock <- fread(sens_lock_path)
manifest <- fread(manifest_path, na.strings = c("NA", ""))
expected_primary_bundle <- "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30"
expected_sens_bundle <- "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8"
if (nrow(lock) != 36L || unique(lock$bundle_sha256) != expected_primary_bundle) stop("Primary bundle lock mismatch")
if (nrow(sens_lock) != 35L || any(sens_lock$gene == "INPP5D") || unique(sens_lock$parent_primary_bundle_sha256) != expected_primary_bundle) stop("35-gene sensitivity lock mismatch")

m02 <- lock$gene[order(lock$gene, method = "radix")]
m02_35 <- sens_lock$gene[order(sens_lock$gene, method = "radix")]

scenario_defs <- data.table(
  scenario_id = c("PRIMARY_ALL_EPI_MIN50", "CELL_MIN25_SENS", "CELL_MIN100_SENS", "DEPOSIT_QC_SENS", "STEM_PROGENITOR_SECONDARY"),
  qc_mode = c("primary_qc", "primary_qc", "primary_qc", "deposit_preserving", "primary_qc"),
  annotation_source = "independent_broad_annotation",
  compartment = c("all_epithelial", "all_epithelial", "all_epithelial", "all_epithelial", "Stem/progenitor"),
  min_cells = c(50L, 25L, 100L, 50L, 20L),
  analysis_role = c("primary", "cell_count_sensitivity", "cell_count_sensitivity", "qc_sensitivity", "secondary")
)

object_cache <- new.env(parent = emptyenv())
load_object <- function(sample_id, qc_mode) {
  key <- paste(sample_id, qc_mode, sep = "__")
  if (!exists(key, envir = object_cache, inherits = FALSE)) {
    path <- file.path(objects_dir, paste0(key, ".rds"))
    if (!file.exists(path)) stop("Missing M02-blind object: ", path)
    assign(key, readRDS(path), envir = object_cache)
  }
  get(key, envir = object_cache, inherits = FALSE)
}

get_cells <- function(obj, compartment) {
  if (compartment == "all_epithelial") return(obj$cell_meta[broad_type == "epithelial", cell_id])
  obj$cell_meta[grepl("^Stem/progenitor", epithelial_state), cell_id]
}

ucell_mean <- function(counts, cells, genes) {
  if (!length(cells)) return(NA_real_)
  mat <- counts[, cells, drop = FALSE]
  score <- UCell::ScoreSignatures_UCell(mat, features = list(M02 = genes), maxRank = min(1500L, nrow(mat)))
  score_col <- grep("M02_UCell$", colnames(score), value = TRUE)
  if (length(score_col) != 1L) stop("Unexpected UCell output")
  mean(as.numeric(score[, score_col]))
}

score_scenario <- function(def) {
  pb <- list()
  n_cells <- integer()
  ucell <- numeric()
  eligible_samples <- character()
  all_genes <- NULL
  for (sample_id in manifest$sample_id) {
    obj <- load_object(sample_id, def$qc_mode)
    cells <- get_cells(obj, def$compartment)
    if (length(cells) < def$min_cells) next
    if (is.null(all_genes)) all_genes <- rownames(obj$counts)
    if (!identical(all_genes, rownames(obj$counts))) stop("Gene order mismatch across sample objects")
    pb[[sample_id]] <- Matrix::rowSums(obj$counts[, cells, drop = FALSE])
    n_cells[sample_id] <- length(cells)
    ucell[sample_id] <- ucell_mean(obj$counts, cells, m02)
    eligible_samples <- c(eligible_samples, sample_id)
  }
  if (!length(pb)) return(list(scores = data.table(), coverage = data.table()))
  count_matrix <- do.call(cbind, pb)
  rownames(count_matrix) <- all_genes
  colnames(count_matrix) <- names(pb)
  dge <- DGEList(counts = count_matrix)
  dge <- calcNormFactors(dge, method = "TMM")
  logcpm <- cpm(dge, log = TRUE, prior.count = 2)
  mapped36 <- intersect(m02, rownames(logcpm))
  mapped35 <- intersect(m02_35, rownames(logcpm))
  if (length(mapped36) < 29L || length(mapped35) < 28L) return(list(scores = data.table(), coverage = data.table(
    scenario_id = def$scenario_id, mapped36 = length(mapped36), mapped35 = length(mapped35), status = "NOT_ESTIMABLE")))
  primary_estimable <- length(mapped36) == 36L
  gene_z <- t(scale(t(logcpm[mapped36, , drop = FALSE])))
  gene_z[!is.finite(gene_z)] <- 0
  values <- list(
    M02_SCORE_V1 = if (primary_estimable) colMeans(logcpm[mapped36, , drop = FALSE]) else rep(NA_real_, ncol(logcpm)),
    M02_29OF36_COVERAGE_SENS = colMeans(logcpm[mapped36, , drop = FALSE]),
    M02_MINUS_INPP5D_SENS_V1 = colMeans(logcpm[mapped35, , drop = FALSE]),
    M02_GENE_Z_SENS = colMeans(gene_z),
    M02_UCELL_SENS = unname(ucell[colnames(logcpm)])
  )
  roles <- c(M02_SCORE_V1 = "primary", M02_29OF36_COVERAGE_SENS = "coverage_sensitivity",
             M02_MINUS_INPP5D_SENS_V1 = "contamination_sensitivity", M02_GENE_Z_SENS = "score_sensitivity",
             M02_UCELL_SENS = "score_sensitivity")
  rows <- rbindlist(lapply(names(values), function(method) data.table(
    scenario_id = def$scenario_id, qc_mode = def$qc_mode, annotation_source = def$annotation_source,
    compartment = def$compartment, min_cells = def$min_cells, scenario_role = def$analysis_role,
    score_method = method, score_role = roles[[method]], sample_id = colnames(logcpm),
    n_cells = unname(n_cells[colnames(logcpm)]), mapped_genes = if (method == "M02_MINUS_INPP5D_SENS_V1") length(mapped35) else length(mapped36),
    required_genes = if (method == "M02_SCORE_V1") 36L else if (method == "M02_MINUS_INPP5D_SENS_V1") 35L else 29L,
    score = as.numeric(values[[method]]), normalization = if (method == "M02_UCELL_SENS") "UCell_then_patient_sample_mean" else if (method == "M02_GENE_Z_SENS") "TMM_log2CPM_prior2_then_gene_z" else "TMM_log2CPM_prior2",
    inferential_unit = "patient", status = ifelse(is.finite(values[[method]]), "ESTIMABLE", "NOT_ESTIMABLE")
  )), fill = TRUE)
  rows <- merge(rows, manifest[, .(patient_id, sample_id, biological_sample_type, lesion_morphology, histology, anatomic_site, pair_id, FAP_status, germline_APC_status)], by = "sample_id", all.x = TRUE, sort = FALSE)
  coverage <- data.table(scenario_id = def$scenario_id, mapped36 = length(mapped36), mapped35 = length(mapped35),
                         primary_36of36 = primary_estimable, sensitivity_29of36 = length(mapped36) >= 29L,
                         sensitivity_35gene_80pct = length(mapped35) >= 28L,
                         status = if (primary_estimable) "PRIMARY_ESTIMABLE" else "SENSITIVITY_ONLY_PRIMARY_NOT_ESTIMABLE")
  list(scores = rows, coverage = coverage)
}

scenario_results <- lapply(seq_len(nrow(scenario_defs)), function(i) score_scenario(scenario_defs[i]))
scores <- rbindlist(lapply(scenario_results, `[[`, "scores"), fill = TRUE)
coverage <- rbindlist(lapply(scenario_results, `[[`, "coverage"), fill = TRUE)
if (!nrow(scores)) stop("No estimable score scenario")

pair_spec <- data.table(patient_id = c("P1", "P5"), lesion_sample = c("P1_L", "P5_L"), normal_sample = c("P1_N", "P5_N"))
effect_rows <- list()
combo_cols <- c("scenario_id", "qc_mode", "annotation_source", "compartment", "min_cells", "scenario_role", "score_method", "score_role")
combos <- unique(scores[, ..combo_cols])
for (j in seq_len(nrow(combos))) {
  combo <- combos[j]
  sub <- scores[combo, on = combo_cols]
  indiv <- list()
  for (i in seq_len(nrow(pair_spec))) {
    p <- pair_spec[i]
    lesion <- sub[sample_id == p$lesion_sample, score]
    normal <- sub[sample_id == p$normal_sample, score]
    estimable <- length(lesion) == 1L && length(normal) == 1L && is.finite(lesion) && is.finite(normal)
    eff <- if (estimable) lesion - normal else NA_real_
    indiv[[i]] <- data.table(
      record_type = "individual_patient_effect", scenario_id = combo$scenario_id, qc_mode = combo$qc_mode,
      annotation_source = combo$annotation_source, compartment = combo$compartment, min_cells = combo$min_cells,
      scenario_role = combo$scenario_role, score_method = combo$score_method, score_role = combo$score_role,
      patient_id = p$patient_id, lesion_sample = p$lesion_sample, normal_sample = p$normal_sample,
      lesion_score = if (length(lesion)) lesion[1] else NA_real_, normal_score = if (length(normal)) normal[1] else NA_real_,
      effect_lesion_minus_normal = eff, direction_positive = if (estimable) eff > 0 else NA,
      k_positive = NA_integer_, n_pairs = 1L, descriptive_mean_effect = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_, fdr = NA_real_,
      inferential_unit = "patient", status = if (estimable) "ESTIMABLE" else "NOT_ESTIMABLE",
      reason = if (estimable) "individual matched-patient directional effect; no population inference" else "required matched sample score unavailable"
    )
  }
  indiv <- rbindlist(indiv)
  evaluable <- indiv[status == "ESTIMABLE"]
  summary_row <- data.table(
    record_type = "direction_summary", scenario_id = combo$scenario_id, qc_mode = combo$qc_mode,
    annotation_source = combo$annotation_source, compartment = combo$compartment, min_cells = combo$min_cells,
    scenario_role = combo$scenario_role, score_method = combo$score_method, score_role = combo$score_role,
    patient_id = "P1+P5", lesion_sample = "P1_L;P5_L", normal_sample = "P1_N;P5_N",
    lesion_score = NA_real_, normal_score = NA_real_, effect_lesion_minus_normal = NA_real_, direction_positive = NA,
    k_positive = sum(evaluable$effect_lesion_minus_normal > 0), n_pairs = nrow(evaluable),
    descriptive_mean_effect = if (nrow(evaluable)) mean(evaluable$effect_lesion_minus_normal) else NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_, fdr = NA_real_, inferential_unit = "patient",
    status = if (nrow(evaluable) == 2L) "ESTIMABLE" else "NOT_ESTIMABLE",
    reason = "two directional patient effects; no population inference"
  )
  effect_rows[[length(effect_rows) + 1L]] <- rbind(indiv, summary_row, fill = TRUE)
}
effects <- rbindlist(effect_rows, fill = TRUE)

sensitivity <- effects[record_type == "direction_summary", .(
  scenario_id, qc_mode, annotation_source, compartment, min_cells, scenario_role, score_method, score_role,
  P1_effect = effects[record_type == "individual_patient_effect" & patient_id == "P1"][.SD, on = combo_cols, x.effect_lesion_minus_normal],
  P5_effect = effects[record_type == "individual_patient_effect" & patient_id == "P5"][.SD, on = combo_cols, x.effect_lesion_minus_normal],
  k_positive, n_pairs, descriptive_mean_effect, status, reason
)]
sensitivity <- rbind(sensitivity, data.table(
  scenario_id = "PUBLISHED_ANNOTATION_SENS", qc_mode = "primary_qc", annotation_source = "published_cell_level_annotation",
  compartment = "all_epithelial", min_cells = 50L, scenario_role = "annotation_sensitivity", score_method = "M02_SCORE_V1",
  score_role = "primary_under_annotation_sensitivity", P1_effect = NA_real_, P5_effect = NA_real_, k_positive = NA_integer_, n_pairs = 0L,
  descriptive_mean_effect = NA_real_, status = "NOT_ESTIMABLE",
  reason = "official Figshare archive contains no deposited cell-level published annotation; no reconstructed label substituted"
), fill = TRUE)

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

primary_scores <- scores[scenario_id == "PRIMARY_ALL_EPI_MIN50" & score_method == "M02_SCORE_V1"]
display_scores <- scores[scenario_id == "PRIMARY_ALL_EPI_MIN50" & score_method == "M02_29OF36_COVERAGE_SENS"]
secondary <- copy(display_scores)
secondary[, description_type := ifelse(biological_sample_type == "normal mucosa", "normal_sample", "lesion_sample")]
p3_l <- primary_scores[sample_id == "P3_L", score]
p3_p <- primary_scores[sample_id == "P3_P", score]
if (length(p3_l) == 1L && length(p3_p) == 1L) {
  secondary <- rbind(secondary, data.table(
    sample_id = "P3_L-minus-P3_P", scenario_id = "PRIMARY_ALL_EPI_MIN50", qc_mode = "primary_qc",
    annotation_source = "independent_broad_annotation", compartment = "all_epithelial", min_cells = 50L,
    scenario_role = "secondary_descriptive", score_method = "M02_SCORE_V1", score_role = "primary_descriptive",
    n_cells = NA_integer_, mapped_genes = 36L, required_genes = 36L, score = p3_l - p3_p,
    normalization = "within-patient descriptive contrast", inferential_unit = "patient", status = "ESTIMABLE",
    patient_id = "P3", biological_sample_type = "two lesions", lesion_morphology = "LST-G minus protruded adenoma",
    histology = "within-patient descriptive", anatomic_site = "different sites", pair_id = "P3_MULTI_LESION",
    FAP_status = NA_character_, germline_APC_status = NA_character_, description_type = "within_patient_lesion_contrast"
  ), fill = TRUE)
}

fwrite(scores, final_score_path, sep = "\t", na = "NA")
fwrite(effects, file.path(results_dir, "STAGE10C3_MATCHED_PAIR_EFFECTS.tsv"), sep = "\t", na = "NA")
fwrite(secondary, file.path(results_dir, "STAGE10C3_SECONDARY_DESCRIPTIVE.tsv"), sep = "\t", na = "NA")
fwrite(sensitivity, file.path(results_dir, "STAGE10C3_SENSITIVITY.tsv"), sep = "\t", na = "NA")
fwrite(coverage, file.path(results_dir, "STAGE10C3_GENE_COVERAGE.tsv"), sep = "\t", na = "NA")

fig1 <- display_scores[sample_id %in% c("P1_N", "P1_L", "P5_N", "P5_L")]
fig1[, tissue_order := factor(ifelse(biological_sample_type == "normal mucosa", "Normal", "LST-G"), levels = c("Normal", "LST-G"))]
fwrite(fig1, file.path(figures_dir, "FIGURE10C3_1_source_data.tsv"), sep = "\t", na = "NA")
p1 <- ggplot(fig1, aes(tissue_order, score, group = patient_id, color = patient_id)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  labs(x = NULL, y = "29/36 M02 coverage-sensitivity score", color = "Patient", title = "Matched 29/36 coverage sensitivity\n(primary 36/36 score not estimable)") +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(figures_dir, "FIGURE10C3_1_matched_pairs.png"), p1, width = 5.2, height = 4.2, dpi = 300)
ggsave(file.path(figures_dir, "FIGURE10C3_1_matched_pairs.pdf"), p1, width = 5.2, height = 4.2)

fig2 <- display_scores
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
ggsave(file.path(figures_dir, "FIGURE10C3_S1_sensitivity.png"), p3, width = 8.5, height = max(5, 0.20 * uniqueN(fig3$display_method)), dpi = 300, limitsize = FALSE)
ggsave(file.path(figures_dir, "FIGURE10C3_S1_sensitivity.pdf"), p3, width = 8.5, height = max(5, 0.20 * uniqueN(fig3$display_method)), limitsize = FALSE)

primary_effect_text <- if (primary_estimable) {
  paste(sprintf("- %s: effect `%0.6f` (%s)", primary_indiv$patient_id, primary_indiv$effect_lesion_minus_normal,
                ifelse(primary_indiv$effect_lesion_minus_normal > 0, "positive", "non-positive")), collapse = "\n")
} else {
  "- P1: `NOT_ESTIMABLE` (33/36 genes mapped)\n- P5: `NOT_ESTIMABLE` (33/36 genes mapped)"
}
primary_k <- if (primary_estimable) sum(primary_indiv$effect_lesion_minus_normal > 0) else NA_integer_
sens35_text <- paste(sprintf("- %s: effect `%0.6f`", sens35_indiv$patient_id, sens35_indiv$effect_lesion_minus_normal), collapse = "\n")

claim_md <- paste0(
  "# Stage 10C3 Claim Limits\n\n",
  "This stage is an LST-G directional sensitivity analysis only. It is not confirmatory validation, ",
  "sporadic adenoma validation, generalization, mechanism analysis or biomarker research.\n\n",
  "The patient is the inferential unit. Cells, samples and libraries are nested. P3_L and P3_P are one patient's two lesions.\n\n",
  "FAP status and germline APC status are unknown for every participant. The two matched donors cannot establish a population effect.\n\n",
  "The strongest permitted wording, and only if the decision is `LST_DIRECTIONAL_CONCORDANCE`, is:\n\n",
  "> the locked M02 showed concordant directional behavior in two matched normal-LST-G donors, with unknown hereditary status.\n\n",
  "No result permits changing M02, restarting dataset search, claiming clinical performance or inferring causality.\n"
)
writeLines(claim_md, file.path(results_dir, "STAGE10C3_CLAIM_LIMITS.md"), useBytes = TRUE)

decision_md <- paste0(
  "# Stage 10C3 Decision\n\n",
  "Date: ", format(Sys.Date()), "\n\n",
  "Decision: **", decision, "**\n\n",
  "Role: LST-G directional sensitivity analysis only.\n\n",
  "## Primary matched effects\n\n", primary_effect_text, "\n\n",
  "Direction summary: `", ifelse(is.na(primary_k), "NA", primary_k), "/2` for the primary score because 36/36 coverage was unavailable. ",
  "No population confidence interval, P value or FDR was estimated from two donors.\n\n",
  "## Frozen contamination sensitivity\n\n", sens35_text, "\n\n",
  "No complete prespecified sensitivity produced two systematically opposite patient effects: `", !systemic_opposite, "`.\n\n",
  "Primary bundle SHA256: `", expected_primary_bundle, "`\n\n",
  "35-gene sensitivity bundle SHA256: `", expected_sens_bundle, "`\n\n",
  "## Interpretation ceiling\n\n",
  if (decision == "LST_DIRECTIONAL_CONCORDANCE") "the locked M02 showed concordant directional behavior in two matched normal-LST-G donors, with unknown hereditary status.\n\n" else "The locked M02 did not meet the prespecified concordant directional criterion.\n\n",
  "The next separately authorized stage is Stage 10C2-SP. This file does not start it.\n"
)
writeLines(decision_md, file.path(results_dir, "STAGE10C3_DECISION.md"), useBytes = TRUE)

writeLines(capture.output(sessionInfo()), file.path(results_dir, "STAGE10C3_SESSIONINFO.txt"), useBytes = TRUE)
message("Stage 10C3 scoring complete: ", decision)
