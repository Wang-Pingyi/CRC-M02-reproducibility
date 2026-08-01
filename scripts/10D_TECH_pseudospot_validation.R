#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args)) stop("Missing value for ", flag)
  args[hit[length(hit)] + 1L]
}

root <- normalizePath(arg_value("--root", "."), mustWork = TRUE)
smoke <- "--smoke" %in% args
project_lib <- file.path(root, "environment", "R-library")
if (dir.exists(project_lib)) .libPaths(unique(c(.libPaths(), project_lib)))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(edgeR)
  library(UCell)
  library(ggplot2)
  library(digest)
})

if (packageVersion("Seurat") < "5.0.0" || packageVersion("SeuratObject") < "5.0.0") {
  stop("Stage 10D-TECH requires Seurat/SeuratObject >=5 to read the frozen reference object")
}

options(stringsAsFactors = FALSE, scipen = 999)

expected <- list(
  primary_bundle = "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30",
  sensitivity_bundle = "78f7bec53e00fd3226ac89151872e9bcc8eedbab031dcfc44c41ad32ae40b8c8",
  common_geneset = "4cd62d74b83673a4d2adf6077bedbdfe73d1cbd369a6a77418d124e0b506d482"
)

parameter_path <- file.path(root, "config", "stage10d_tech_parameters.tsv")
design_path <- file.path(root, "results", "stage10d_tech", "STAGE10D_TECH_DESIGN.md")
mapped_path <- file.path(root, "results", "stage10d_tech", "STAGE10D_TECH_MAPPED_GENESET.tsv")
lock36_path <- file.path(root, "results", "stage10c", "STAGE10C_LOCK_MANIFEST.tsv")
lock35_path <- file.path(root, "results", "stage10c", "M02_MINUS_INPP5D_SENS_V1.tsv")
common_path <- file.path(root, "results", "stage10c2_sp", "STAGE10C2_SP_COHORT_COMMON_GENESET.tsv")
decision_sp_path <- file.path(root, "results", "stage10c2_sp", "STAGE10C2_SP_DECISION.md")
reference_path <- file.path(root, "objects", "GSE161277_stage7_annotated.rds")

required_inputs <- c(parameter_path, design_path, mapped_path, lock36_path, lock35_path,
                     common_path, decision_sp_path, reference_path)
if (!all(file.exists(required_inputs))) {
  stop("Missing Stage 10D-TECH input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = "; "))
}

params <- fread(parameter_path, na.strings = c("NA", ""))
get_param <- function(name) {
  value <- params[parameter == name, value]
  if (length(value) != 1L) stop("Expected exactly one parameter: ", name)
  value
}
as_num_levels <- function(name) as.numeric(strsplit(get_param(name), ";", fixed = TRUE)[[1L]])

seed <- as.integer(get_param("seed"))
set.seed(seed)

hash_gene_vector <- function(x) digest(paste0(paste(sort(unique(x)), collapse = "\n"), "\n"),
                                       algo = "sha256", serialize = FALSE)
hash_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)

decision_sp <- paste(readLines(decision_sp_path, warn = FALSE), collapse = "\n")
if (!grepl("Decision: \\*\\*PASS\\*\\*", decision_sp)) stop("Stage 10C2-SP is not PASS")

lock36 <- fread(lock36_path)
lock35 <- fread(lock35_path)
if (nrow(lock36) != 36L || length(unique(lock36$gene)) != 36L ||
    unique(lock36$bundle_sha256) != expected$primary_bundle) stop("Primary 36-gene lock mismatch")
if (nrow(lock35) != 35L || "INPP5D" %in% lock35$gene) stop("35-gene sensitivity lock mismatch")
if ("bundle_sha256" %in% names(lock35) && unique(lock35$bundle_sha256) != expected$sensitivity_bundle) {
  stop("35-gene sensitivity bundle mismatch")
}

canonical36 <- lock36$gene[order(lock36$gene_order)]
common_all <- fread(common_path)
common30 <- common_all[cohort_id == "E-GEAD-622", canonical_gene]
if (length(common30) != 30L || hash_gene_vector(common30) != expected$common_geneset) {
  stop("Frozen E-GEAD-622 cohort-common gene set mismatch")
}
mapped_out <- fread(mapped_path)
if (nrow(mapped_out) != 30L || !setequal(mapped_out$canonical_gene, common30) ||
    unique(mapped_out$geneset_sha256) != expected$common_geneset) {
  stop("Stage 10D mapped-gene output does not reproduce the frozen set")
}

marker_sets <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "ELF3", "CEACAM5", "CEACAM6", "GUCA2A", "MUC2"),
  immune = c("PTPRC", "CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "CD79A", "MS4A1", "LYZ", "FCER1G", "TYROBP", "LST1"),
  stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "PDGFRA", "PECAM1", "VWF", "EMCN", "KDR")
)
if (length(intersect(unique(unlist(marker_sets)), canonical36))) {
  stop("M02 genes are prohibited from reference markers and proportion model")
}

message("Reading independent GSE161277 cell-level reference")
reference <- readRDS(reference_path)
counts_all <- tryCatch(
  GetAssayData(reference, assay = "RNA", layer = "counts"),
  error = function(e) GetAssayData(reference, assay = "RNA", slot = "counts")
)
meta <- reference[[]]
if (!all(c("major_cell_type", "donor_id", "annotation_status") %in% names(meta))) {
  stop("Reference metadata lacks frozen broad annotation fields")
}
if (!all(canonical36 %in% rownames(counts_all))) stop("36-gene oracle not available in reference")

cell_group <- rep(NA_character_, nrow(meta))
cell_group[meta$major_cell_type == "Epithelial"] <- "epithelial"
cell_group[meta$major_cell_type %in% c("B_cell", "Plasma_cell", "Myeloid", "T_NK")] <- "immune"
cell_group[meta$major_cell_type %in% c("Fibroblast", "Endothelial")] <- "stromal"
qualified <- !is.na(cell_group) & meta$annotation_status == "canonical_marker_cluster_annotation"
if (sum(qualified) < 1000L || length(unique(meta$donor_id[qualified])) < 3L) {
  stop("Independent reference does not retain sufficient labelled cells/donors")
}
counts_ref <- counts_all[, qualified, drop = FALSE]
cell_group <- cell_group[qualified]
reference_donors <- unique(meta$donor_id[qualified])
rm(reference, counts_all, meta)
gc(verbose = FALSE)

if (!inherits(counts_ref, "dgCMatrix")) counts_ref <- as(counts_ref, "dgCMatrix")
detection <- tabulate(counts_ref@i + 1L, nbins = nrow(counts_ref))
names(detection) <- rownames(counts_ref)
totals <- Matrix::rowSums(counts_ref)
names(totals) <- rownames(counts_ref)

background_n <- as.integer(get_param("background_gene_count"))
eligible_background <- setdiff(
  rownames(counts_ref)[detection >= 10L & totals > 0 &
                         !grepl("^MT-|^RP[SL]", rownames(counts_ref))],
  canonical36
)
ord <- order(detection[eligible_background], totals[eligible_background], decreasing = TRUE)
background <- eligible_background[ord][seq_len(min(background_n, length(eligible_background)))]
marker_genes <- intersect(unique(unlist(marker_sets)), rownames(counts_ref))
simulation_genes <- unique(c(canonical36, marker_genes, background))
if (length(background) < 1000L || any(lengths(lapply(marker_sets, intersect, y = simulation_genes)) < 4L)) {
  stop("Reference background or marker coverage is insufficient")
}

counts_selected <- counts_ref[simulation_genes, , drop = FALSE]
profile <- vapply(c("epithelial", "immune", "stromal"), function(group) {
  as.numeric(Matrix::rowSums(counts_selected[, cell_group == group, drop = FALSE])) + 0.5
}, numeric(length(simulation_genes)))
rownames(profile) <- simulation_genes
profile <- sweep(profile, 2L, colSums(profile), "/")

marker_use <- lapply(marker_sets, intersect, y = simulation_genes)
marker_union <- unique(unlist(marker_use))
reference_signature <- profile[marker_union, , drop = FALSE]
reference_signature <- sweep(reference_signature, 2L, colSums(reference_signature), "/")

estimate_proportions <- function(count_matrix) {
  y <- count_matrix[marker_union, , drop = FALSE]
  y <- sweep(y, 2L, pmax(colSums(y), 1), "/")
  ridge <- diag(1e-8, ncol(reference_signature))
  solve_one <- function(v) {
    b <- tryCatch(solve(crossprod(reference_signature) + ridge,
                        crossprod(reference_signature, v)), error = function(e) rep(NA_real_, 3L))
    b <- pmax(as.numeric(b), 0)
    if (!all(is.finite(b)) || sum(b) <= 0) return(rep(NA_real_, 3L))
    b / sum(b)
  }
  out <- vapply(seq_len(ncol(y)), function(i) solve_one(y[, i]), numeric(3L))
  rownames(out) <- colnames(reference_signature)
  colnames(out) <- colnames(y)
  out
}

smooth_field <- function(nspots, rho) {
  side <- as.integer(round(sqrt(nspots)))
  if (side * side != nspots) stop("Region spot count must be a square")
  raw <- matrix(rnorm(nspots), nrow = side)
  sm <- raw
  for (iteration in seq_len(3L)) {
    up <- rbind(sm[1L, , drop = FALSE], sm[-side, , drop = FALSE])
    down <- rbind(sm[-1L, , drop = FALSE], sm[side, , drop = FALSE])
    left <- cbind(sm[, 1L, drop = FALSE], sm[, -side, drop = FALSE])
    right <- cbind(sm[, -1L, drop = FALSE], sm[, side, drop = FALSE])
    sm <- (sm + up + down + left + right) / 5
  }
  raw <- as.numeric(scale(as.numeric(raw)))
  sm <- as.numeric(scale(as.numeric(sm)))
  out <- sqrt(max(0, 1 - rho^2)) * raw + rho * sm
  as.numeric(scale(out))
}

level_values <- list(
  umi_depth = as_num_levels("umi_depth_levels"),
  dropout = as_num_levels("dropout_levels"),
  epithelial_fraction = as_num_levels("epithelial_fraction_levels"),
  immune_contamination = as_num_levels("immune_contamination_levels"),
  region_spots = as_num_levels("region_spot_levels"),
  batch_log_sd = as_num_levels("batch_log_sd_levels"),
  spatial_rho = as_num_levels("spatial_rho_levels")
)
n_sim <- if (smoke) 3L else as.integer(get_param("n_simulation_replicates"))
n_regions <- if (smoke) 4L else as.integer(get_param("n_regions_per_replicate"))
design <- data.table(sim_id = sprintf("SIM%03d", seq_len(n_sim)))
set.seed(seed)
for (nm in names(level_values)) {
  values <- rep(level_values[[nm]], length.out = n_sim)
  design[, (nm) := sample(values, length(values), replace = FALSE)]
}
if (smoke) design[, region_spots := 16]
design_out_dir <- if (smoke) file.path(root, "cache", "stage10d_tech_smoke") else
  file.path(root, "results", "stage10d_tech")
dir.create(design_out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(design, file.path(design_out_dir, "STAGE10D_TECH_SIMULATION_DESIGN.tsv"), sep = "\t", na = "NA")

program_scale <- as.numeric(get_param("program_log2fc_per_latent_unit"))
all_region_ids <- unlist(lapply(design$sim_id, function(x) paste0(x, "_R", seq_len(n_regions))))
region_counts <- matrix(0L, nrow = length(simulation_genes), ncol = length(all_region_ids),
                        dimnames = list(simulation_genes, all_region_ids))
region_meta <- vector("list", length(all_region_ids))
region_cursor <- 0L

message("Generating blinded pseudo-spots")
for (s in seq_len(nrow(design))) {
  d <- design[s]
  gene_shift <- list(
    A = rep(1, length(simulation_genes)),
    B = exp(rnorm(length(simulation_genes), mean = 0, sd = d$batch_log_sd))
  )
  targets <- seq(-1.75, 1.75, length.out = n_regions)
  for (r in seq_len(n_regions)) {
    region_cursor <- region_cursor + 1L
    field <- smooth_field(as.integer(d$region_spots), d$spatial_rho)
    epi_fraction <- plogis(qlogis(d$epithelial_fraction) + 0.45 * field + rnorm(length(field), 0, 0.12))
    immune_share <- plogis(qlogis(d$immune_contamination) + 0.30 * field + rnorm(length(field), 0, 0.10))
    latent <- targets[r] + 0.40 * field
    batch <- if (r %% 2L) "A" else "B"
    aggregate_counts <- integer(length(simulation_genes))
    for (spot in seq_along(field)) {
      epithelial_profile <- profile[, "epithelial"]
      epithelial_profile[canonical36] <- epithelial_profile[canonical36] *
        2^(program_scale * latent[spot])
      epithelial_profile <- epithelial_profile / sum(epithelial_profile)
      non_epithelial <- immune_share[spot] * profile[, "immune"] +
        (1 - immune_share[spot]) * profile[, "stromal"]
      probability <- epi_fraction[spot] * epithelial_profile +
        (1 - epi_fraction[spot]) * non_epithelial
      probability <- probability * gene_shift[[batch]]
      if (d$dropout > 0) {
        relative_abundance <- sqrt(probability / max(probability))
        drop_probability <- d$dropout * (1 - relative_abundance)
        probability[runif(length(probability)) < drop_probability] <- 0
      }
      if (sum(probability) <= 0) stop("All genes dropped in pseudo-spot")
      probability <- probability / sum(probability)
      umi <- max(200L, rpois(1L, lambda = d$umi_depth))
      aggregate_counts <- aggregate_counts + as.integer(rmultinom(1L, size = umi, prob = probability)[, 1L])
    }
    region_counts[, region_cursor] <- aggregate_counts
    region_meta[[region_cursor]] <- data.table(
      sim_id = d$sim_id, region_id = colnames(region_counts)[region_cursor], region_index = r,
      batch = batch, latent_truth = mean(latent), true_epithelial_fraction = mean(epi_fraction),
      true_immune_share_non_epi = mean(immune_share), total_umi = sum(aggregate_counts),
      umi_depth = d$umi_depth, dropout = d$dropout,
      epithelial_fraction_setting = d$epithelial_fraction,
      immune_contamination_setting = d$immune_contamination,
      region_spots = d$region_spots, batch_log_sd = d$batch_log_sd, spatial_rho = d$spatial_rho
    )
  }
  if (s %% 12L == 0L || s == nrow(design)) message("simulations_completed=", s, "/", nrow(design))
}
region_meta <- rbindlist(region_meta)

message("Applying within-replicate TMM and fixed scores")
logcpm <- matrix(NA_real_, nrow = nrow(region_counts), ncol = ncol(region_counts),
                 dimnames = dimnames(region_counts))
for (sim in unique(region_meta$sim_id)) {
  idx <- which(region_meta$sim_id == sim)
  y <- DGEList(counts = region_counts[, idx, drop = FALSE])
  y <- calcNormFactors(y, method = "TMM")
  logcpm[, idx] <- cpm(y, log = TRUE, prior.count = 2)
}

proportions <- estimate_proportions(region_counts)
region_meta[, estimated_epithelial_fraction := proportions["epithelial", region_id]]
region_meta[, primary_score := colMeans(logcpm[common30, region_id, drop = FALSE])]
region_meta[, oracle36_score := colMeans(logcpm[canonical36, region_id, drop = FALSE])]
minus_inpp5d <- setdiff(common30, "INPP5D")
if (length(minus_inpp5d) != 29L) stop("M02_MINUS_INPP5D is not distinguishable from primary")
region_meta[, minus_inpp5d_score := colMeans(logcpm[minus_inpp5d, region_id, drop = FALSE])]

gene_z_score <- rep(NA_real_, nrow(region_meta))
for (sim in unique(region_meta$sim_id)) {
  idx <- which(region_meta$sim_id == sim)
  z <- t(scale(t(logcpm[common30, region_meta$region_id[idx], drop = FALSE])))
  z[!is.finite(z)] <- 0
  gene_z_score[idx] <- colMeans(z)
}
region_meta[, gene_z_score := gene_z_score]

ucell_result <- UCell::ScoreSignatures_UCell(
  as(region_counts, "dgCMatrix"), features = list(M02 = common30),
  maxRank = as.integer(get_param("ucell_max_rank")), ncores = 1L
)
ucell_column <- grep("M02.*UCell|M02", colnames(ucell_result), value = TRUE)[1L]
if (!length(ucell_column) || is.na(ucell_column)) stop("UCell output column not found")
region_meta[, ucell_score := as.numeric(ucell_result[region_id, ucell_column])]

z_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

replicate_metrics <- function(meta_dt, score_column) {
  rbindlist(lapply(unique(meta_dt$sim_id), function(sim) {
    z <- meta_dt[sim_id == sim]
    ok <- is.finite(z[[score_column]]) & is.finite(z$oracle36_score) & is.finite(z$latent_truth)
    z <- z[ok]
    if (nrow(z) < 3L || sd(z[[score_column]]) <= 0 || sd(z$oracle36_score) <= 0) {
      return(data.table(sim_id = sim, pearson = NA_real_, spearman = NA_real_,
                        standardized_bias = NA_real_, standardized_rmse = NA_real_,
                        calibration_intercept = NA_real_, calibration_slope = NA_real_,
                        dynamic_range_ratio = NA_real_, epithelial_residual_slope = NA_real_,
                        estimated_epithelial_residual_slope = NA_real_, failed = TRUE))
    }
    oracle_sd <- sd(z$oracle36_score)
    obs_std <- (z[[score_column]] - mean(z[[score_column]])) / oracle_sd
    oracle_std <- (z$oracle36_score - mean(z$oracle36_score)) / oracle_sd
    truth_z <- z_safe(z$latent_truth)
    effect_obs <- coef(lm(z[[score_column]] ~ z$latent_truth))[2L]
    effect_oracle <- coef(lm(z$oracle36_score ~ z$latent_truth))[2L]
    effect_ratio <- if (is.finite(effect_oracle) && abs(effect_oracle) > 1e-12) effect_obs / effect_oracle else NA_real_
    calibration <- coef(lm(obs_std ~ oracle_std))
    error <- obs_std - oracle_std
    epi_slope <- coef(lm(error ~ z_safe(z$true_epithelial_fraction)))[2L]
    estimated_epi_slope <- coef(lm(error ~ z_safe(z$estimated_epithelial_fraction)))[2L]
    pearson <- suppressWarnings(cor(obs_std, truth_z, method = "pearson"))
    spearman <- suppressWarnings(cor(obs_std, truth_z, method = "spearman"))
    dynamic <- IQR(obs_std) / IQR(oracle_std)
    data.table(
      sim_id = sim, pearson = pearson, spearman = spearman,
      standardized_bias = effect_ratio - 1,
      standardized_rmse = sqrt(mean(error^2)),
      calibration_intercept = unname(calibration[1L]),
      calibration_slope = unname(calibration[2L]),
      dynamic_range_ratio = dynamic,
      epithelial_residual_slope = unname(epi_slope),
      estimated_epithelial_residual_slope = unname(estimated_epi_slope),
      failed = !is.finite(pearson) || pearson <= 0
    )
  }))
}

metric_estimators <- c(
  pearson = "median", spearman = "median", standardized_bias = "mean",
  standardized_rmse = "median", calibration_intercept = "median",
  calibration_slope = "median", dynamic_range_ratio = "median",
  epithelial_residual_slope = "median", estimated_epithelial_residual_slope = "median",
  failure_rate = "mean"
)

summarize_rep_metrics <- function(rep_dt, score_id, B = 500L) {
  if (!nrow(rep_dt) || !"failed" %in% names(rep_dt)) {
    return(data.table(
      score_id = score_id, metric = names(metric_estimators), estimate = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, fdr = NA_real_,
      n_reference_donors = length(reference_donors), n_simulation_replicates = 0L,
      inference_unit = "simulation_replicate",
      interval_method = "NOT_ESTIMABLE: no region passed the pre-specified filter"
    ))
  }
  values <- copy(rep_dt)
  values[, failure_rate := as.numeric(failed)]
  summarize_one <- function(x, estimator) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(NA_real_, NA_real_, NA_real_))
    fun <- if (estimator == "mean") mean else median
    estimate <- fun(x)
    if (length(x) < 3L || B <= 0L) return(c(estimate, NA_real_, NA_real_))
    boots <- replicate(B, fun(sample(x, length(x), replace = TRUE)))
    c(estimate, unname(quantile(boots, 0.025, na.rm = TRUE)),
      unname(quantile(boots, 0.975, na.rm = TRUE)))
  }
  out <- rbindlist(lapply(names(metric_estimators), function(metric) {
    estimate <- summarize_one(values[[metric]], metric_estimators[[metric]])
    data.table(score_id = score_id, metric = metric, estimate = estimate[1L],
               ci_low = estimate[2L], ci_high = estimate[3L])
  }))
  out[, `:=`(
    p_value = NA_real_, fdr = NA_real_, n_reference_donors = length(reference_donors),
    n_simulation_replicates = uniqueN(values$sim_id), inference_unit = "simulation_replicate",
    interval_method = paste0("percentile cluster bootstrap; B=", B)
  )]
  out
}

primary_rep <- replicate_metrics(region_meta, "primary_score")
bootstrap_B <- if (smoke) 20L else as.integer(get_param("bootstrap_replicates"))
primary_results <- summarize_rep_metrics(primary_rep, "M02_SPATIAL_SCORE_V1_EGEAD622_30", bootstrap_B)

stress_factors <- c("umi_depth", "dropout", "epithelial_fraction_setting",
                    "immune_contamination_setting", "region_spots", "batch_log_sd", "spatial_rho")
stress_results <- rbindlist(lapply(stress_factors, function(factor_name) {
  rbindlist(lapply(sort(unique(region_meta[[factor_name]])), function(level) {
    z <- region_meta[get(factor_name) == level]
    centered <- z[, .(
      latent_z = z_safe(latent_truth),
      primary_centered = primary_score - mean(primary_score)
    ), by = sim_id]
    data.table(factor = factor_name, level = as.character(level),
               pearson = suppressWarnings(cor(centered$primary_centered, centered$latent_z,
                                               use = "complete.obs")),
               n_regions = nrow(centered), n_simulation_replicates = uniqueN(centered$sim_id))
  }))
}))
worst_stress <- min(stress_results$pearson, na.rm = TRUE)
primary_results <- rbind(
  primary_results,
  data.table(score_id = "M02_SPATIAL_SCORE_V1_EGEAD622_30", metric = "worst_stress_pearson",
             estimate = worst_stress, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
             fdr = NA_real_, n_reference_donors = length(reference_donors),
             n_simulation_replicates = n_sim, inference_unit = "simulation_replicate",
             interval_method = "minimum across pre-specified factor levels"), fill = TRUE
)

score_vector_metrics <- function(score, score_id, B = if (smoke) 10L else 100L, subset = rep(TRUE, nrow(region_meta))) {
  temp <- copy(region_meta[subset])
  temp[, candidate_score := score[subset]]
  summarize_rep_metrics(replicate_metrics(temp, "candidate_score"), score_id, B)
}

message("Evaluating fixed coverage panels")
coverage_results <- list()
coverage_results[[1L]] <- score_vector_metrics(region_meta$primary_score, "actual_frozen_30of36", bootstrap_B)
coverage_results[[1L]][, `:=`(coverage_tier = "ACTUAL_FROZEN", panel_id = "actual30", gene_count = 30L,
                              genes = paste(sort(common30), collapse = ";"))]
coverage_results[[2L]] <- score_vector_metrics(region_meta$oracle36_score, "oracle_36of36", bootstrap_B)
coverage_results[[2L]][, `:=`(coverage_tier = "ORACLE_TECHNICAL_ONLY", panel_id = "oracle36", gene_count = 36L,
                              genes = paste(sort(canonical36), collapse = ";"))]

coverage_index <- 2L
for (drop_gene in common30) {
  coverage_index <- coverage_index + 1L
  genes <- setdiff(common30, drop_gene)
  score <- colMeans(logcpm[genes, region_meta$region_id, drop = FALSE])
  row <- score_vector_metrics(score, paste0("29of36_drop_", drop_gene))
  row[, `:=`(coverage_tier = "HIGH_COVERAGE_29", panel_id = paste0("drop_", drop_gene),
             gene_count = 29L, genes = paste(sort(genes), collapse = ";"))]
  coverage_results[[coverage_index]] <- row
}

set.seed(seed + 22L)
minimum_panels <- replicate(if (smoke) 3L else as.integer(get_param("minimum_coverage_panels")),
                            sort(sample(common30, 22L, replace = FALSE)), simplify = FALSE)
for (i in seq_along(minimum_panels)) {
  coverage_index <- coverage_index + 1L
  genes <- minimum_panels[[i]]
  score <- colMeans(logcpm[genes, region_meta$region_id, drop = FALSE])
  row <- score_vector_metrics(score, sprintf("22of36_seeded_%03d", i))
  row[, `:=`(coverage_tier = "MINIMUM_COVERAGE_22", panel_id = sprintf("seeded_%03d", i),
             gene_count = 22L, genes = paste(genes, collapse = ";"))]
  coverage_results[[coverage_index]] <- row
}
coverage_results <- rbindlist(coverage_results, fill = TRUE)
setcolorder(coverage_results, c("coverage_tier", "panel_id", "score_id", "gene_count", "genes", "metric",
                               "estimate", "ci_low", "ci_high", "p_value", "fdr", "n_reference_donors",
                               "n_simulation_replicates", "inference_unit", "interval_method"))

message("Evaluating pre-specified score and epithelial-fraction sensitivities")
sensitivity_results <- list(
  score_vector_metrics(region_meta$primary_score, "primary_30", bootstrap_B),
  score_vector_metrics(region_meta$minus_inpp5d_score, "M02_MINUS_INPP5D_29", bootstrap_B),
  score_vector_metrics(region_meta$gene_z_score, "gene_z_30", bootstrap_B),
  score_vector_metrics(region_meta$ucell_score, "UCell_30", bootstrap_B)
)
sensitivity_results[[1L]][, `:=`(sensitivity_type = "score", threshold = NA_real_, gene_count = 30L)]
sensitivity_results[[2L]][, `:=`(sensitivity_type = "contamination", threshold = NA_real_, gene_count = 29L)]
sensitivity_results[[3L]][, `:=`(sensitivity_type = "score", threshold = NA_real_, gene_count = 30L)]
sensitivity_results[[4L]][, `:=`(sensitivity_type = "score", threshold = NA_real_, gene_count = 30L)]
for (threshold in as_num_levels("epithelial_thresholds")) {
  keep <- is.finite(region_meta$estimated_epithelial_fraction) &
    region_meta$estimated_epithelial_fraction >= threshold
  row <- score_vector_metrics(region_meta$primary_score,
                              paste0("primary30_estimated_epi_ge_", threshold), bootstrap_B, keep)
  row[, `:=`(sensitivity_type = "estimated_epithelial_fraction_threshold", threshold = threshold,
             gene_count = 30L, n_regions_retained = sum(keep))]
  sensitivity_results[[length(sensitivity_results) + 1L]] <- row
}
sensitivity_results <- rbindlist(sensitivity_results, fill = TRUE)

message("Generating expression/detection-matched random modules of actual size m=30")
gene_mean <- rowMeans(logcpm)
gene_detection <- rowMeans(region_counts > 0)
null_candidates <- setdiff(background, canonical36)
rank_mean <- rank(gene_mean, ties.method = "average") / length(gene_mean)
rank_detection <- rank(gene_detection, ties.method = "average") / length(gene_detection)
names(rank_mean) <- names(gene_mean)
names(rank_detection) <- names(gene_detection)
match_pools <- lapply(common30, function(gene) {
  distance <- (rank_mean[null_candidates] - rank_mean[gene])^2 +
    (rank_detection[null_candidates] - rank_detection[gene])^2
  names(sort(distance))[seq_len(min(80L, length(distance)))]
})
names(match_pools) <- common30

draw_matched_module <- function() {
  selected <- character()
  for (gene in common30) {
    pool <- setdiff(match_pools[[gene]], selected)
    if (!length(pool)) pool <- setdiff(null_candidates, selected)
    selected <- c(selected, sample(pool, 1L))
  }
  selected
}
set.seed(seed + 1000L)
n_null <- if (smoke) 10L else as.integer(get_param("null_modules"))
null_results <- rbindlist(lapply(seq_len(n_null), function(i) {
  genes <- draw_matched_module()
  score <- colMeans(logcpm[genes, region_meta$region_id, drop = FALSE])
  null_meta <- copy(region_meta)
  null_meta[, null_score := score]
  rep_metrics <- replicate_metrics(null_meta, "null_score")
  data.table(
    module_id = sprintf("NULL_MATCHED_M30_%04d", i), gene_count = length(genes),
    genes = paste(sort(genes), collapse = ";"),
    pearson = median(rep_metrics$pearson, na.rm = TRUE),
    spearman = median(rep_metrics$spearman, na.rm = TRUE),
    standardized_rmse = median(rep_metrics$standardized_rmse, na.rm = TRUE),
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, fdr = NA_real_,
    n_reference_donors = length(reference_donors), n_simulation_replicates = n_sim,
    inference_unit = "simulation_replicate",
    interval_method = "empirical null draw; no per-module interval"
  )
}))
primary_pearson <- primary_results[metric == "pearson", estimate]
null_empirical_p <- (1 + sum(null_results$pearson >= primary_pearson, na.rm = TRUE)) / (1 + nrow(null_results))
null_results <- rbind(
  null_results,
  data.table(module_id = "PRIMARY_VS_MATCHED_NULL", gene_count = 30L, genes = paste(sort(common30), collapse = ";"),
             pearson = primary_pearson, spearman = primary_results[metric == "spearman", estimate],
             standardized_rmse = primary_results[metric == "standardized_rmse", estimate],
             ci_low = NA_real_, ci_high = NA_real_, p_value = null_empirical_p, fdr = NA_real_,
             n_reference_donors = length(reference_donors), n_simulation_replicates = n_sim,
             inference_unit = "simulation_replicate",
             interval_method = paste0("empirical one-sided comparison to ", n_null, " matched modules")), fill = TRUE
)

metric_value <- function(name) primary_results[metric == name, estimate][1L]
strict_checks <- c(
  pearson = metric_value("pearson") >= as.numeric(get_param("strict_pearson_min")),
  spearman = metric_value("spearman") >= as.numeric(get_param("strict_spearman_min")),
  bias = abs(metric_value("standardized_bias")) <= as.numeric(get_param("strict_abs_standardized_bias_max")),
  rmse = metric_value("standardized_rmse") <= as.numeric(get_param("strict_standardized_rmse_max")),
  calibration_slope = metric_value("calibration_slope") >= as.numeric(get_param("strict_calibration_slope_min")) &&
    metric_value("calibration_slope") <= as.numeric(get_param("strict_calibration_slope_max")),
  dynamic_range = metric_value("dynamic_range_ratio") >= as.numeric(get_param("strict_dynamic_range_ratio_min")) &&
    metric_value("dynamic_range_ratio") <= as.numeric(get_param("strict_dynamic_range_ratio_max")),
  epithelial_dependence = abs(metric_value("epithelial_residual_slope")) <=
    as.numeric(get_param("strict_abs_epithelial_residual_slope_max")),
  failure_rate = metric_value("failure_rate") <= as.numeric(get_param("strict_failure_rate_max")),
  stress_recovery = metric_value("worst_stress_pearson") >= as.numeric(get_param("strict_worst_stress_pearson_min"))
)
limited_checks <- c(
  pearson = metric_value("pearson") >= as.numeric(get_param("limited_pearson_min")),
  spearman = metric_value("spearman") >= as.numeric(get_param("limited_spearman_min")),
  rmse = metric_value("standardized_rmse") <= as.numeric(get_param("limited_rmse_max")),
  failure_rate = metric_value("failure_rate") <= as.numeric(get_param("limited_failure_rate_max")),
  no_negative_stress = metric_value("worst_stress_pearson") > 0
)
decision <- if (all(strict_checks)) "TECHNICALLY_VALID" else if (all(limited_checks)) {
  "VALID_WITH_LIMITATIONS"
} else "NOT_TECHNICALLY_VALID"

out_dir <- if (smoke) file.path(root, "cache", "stage10d_tech_smoke", "results") else
  file.path(root, "results", "stage10d_tech")
fig_dir <- if (smoke) file.path(root, "cache", "stage10d_tech_smoke", "figures") else
  file.path(root, "figures", "stage10d_tech")
source_dir <- file.path(fig_dir, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(primary_results, file.path(out_dir, "STAGE10D_TECH_PRIMARY_SCORE_RESULTS.tsv"), sep = "\t", na = "NA")
fwrite(coverage_results, file.path(out_dir, "STAGE10D_TECH_COVERAGE_RESULTS.tsv"), sep = "\t", na = "NA")
fwrite(sensitivity_results, file.path(out_dir, "STAGE10D_TECH_SENSITIVITY_RESULTS.tsv"), sep = "\t", na = "NA")
fwrite(null_results, file.path(out_dir, "STAGE10D_TECH_NULL_MODULE_RESULTS.tsv"), sep = "\t", na = "NA")
fwrite(stress_results, file.path(source_dir, "fig2_stress_performance.tsv"), sep = "\t", na = "NA")
fwrite(region_meta, file.path(source_dir, "fig1_primary_recovery.tsv"), sep = "\t", na = "NA")
fwrite(coverage_results, file.path(source_dir, "fig3_coverage_performance.tsv"), sep = "\t", na = "NA")
fwrite(null_results, file.path(source_dir, "fig4_null_module.tsv"), sep = "\t", na = "NA")

theme_stage10d <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
p1 <- ggplot(region_meta, aes(latent_truth, primary_score, color = true_epithelial_fraction)) +
  geom_point(alpha = 0.55, size = 1.3) + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
  scale_color_viridis_c(option = "C") + theme_stage10d +
  labs(x = "Known synthetic M02 amplitude", y = "Frozen 30-gene TMM-log2CPM score",
       color = "True epithelial\nfraction", title = "Blinded pseudo-spot recovery")
p2 <- ggplot(stress_results, aes(level, factor, fill = pearson)) +
  geom_tile(color = "white") + scale_fill_gradient2(limits = c(-1, 1), low = "#B2182B", mid = "white", high = "#2166AC") +
  theme_stage10d + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Pre-specified level", y = NULL, fill = "Pearson r", title = "Recovery across stress factors")
p3 <- ggplot(coverage_results[metric == "pearson" & coverage_tier %in% c("HIGH_COVERAGE_29", "MINIMUM_COVERAGE_22")],
             aes(coverage_tier, estimate, fill = coverage_tier)) +
  geom_violin(trim = FALSE, alpha = 0.6) + geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = primary_pearson, linetype = 2) + theme_stage10d + theme(legend.position = "none") +
  labs(x = NULL, y = "Median replicate Pearson r", title = "Coverage sensitivity (no panel selection)")
p4 <- ggplot(null_results[module_id != "PRIMARY_VS_MATCHED_NULL"], aes(pearson)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(xintercept = primary_pearson, color = "#B2182B", linewidth = 0.8) + theme_stage10d +
  labs(x = "Median replicate Pearson r", y = "Matched random modules",
       title = "Frozen score versus size-30 matched null modules")
for (spec in list(list(plot = p1, name = "Fig10D_1_primary_recovery"),
                  list(plot = p2, name = "Fig10D_2_stress_performance"),
                  list(plot = p3, name = "Fig10D_3_coverage_sensitivity"),
                  list(plot = p4, name = "Fig10D_4_matched_null"))) {
  ggsave(file.path(fig_dir, paste0(spec$name, ".png")), spec$plot, width = 7.2, height = 5.0, dpi = 300)
  ggsave(file.path(fig_dir, paste0(spec$name, ".pdf")), spec$plot, width = 7.2, height = 5.0, device = cairo_pdf)
}

failed_strict <- names(strict_checks)[!strict_checks]
decision_lines <- c(
  "# Stage 10D-TECH decision",
  "",
  paste0("Decision: **", decision, "**"),
  "",
  "This decision uses only the frozen E-GEAD-622 30-gene technical object and pre-specified pseudo-spot simulations. No real spatial lesion-normal M02 result was read or computed.",
  "",
  "## Primary technical metrics",
  "",
  paste0("- Pearson correlation: ", format(metric_value("pearson"), digits = 4)),
  paste0("- Spearman correlation: ", format(metric_value("spearman"), digits = 4)),
  paste0("- Standardized effect bias: ", format(metric_value("standardized_bias"), digits = 4)),
  paste0("- Standardized RMSE: ", format(metric_value("standardized_rmse"), digits = 4)),
  paste0("- Calibration slope: ", format(metric_value("calibration_slope"), digits = 4)),
  paste0("- Dynamic-range ratio: ", format(metric_value("dynamic_range_ratio"), digits = 4)),
  paste0("- Epithelial residual slope: ", format(metric_value("epithelial_residual_slope"), digits = 4)),
  paste0("- Failure rate: ", format(metric_value("failure_rate"), digits = 4)),
  paste0("- Worst pre-specified stress-level Pearson correlation: ", format(metric_value("worst_stress_pearson"), digits = 4)),
  "",
  "## Gate interpretation",
  "",
  if (length(failed_strict)) paste0("Strict criteria not met: ", paste(failed_strict, collapse = ", "), ".") else
    "All strict criteria were met.",
  "Sensitivity scores and the 36-gene oracle were not permitted to rescue the primary score.",
  "This is a technical decision only and provides no spatial or biological evidence.",
  "",
  if (decision == "NOT_TECHNICALLY_VALID")
    "Stage 10E/10F must not use an alternative score as rescue; the spatial primary branch stops unless governance is prospectively revised without viewing outcomes." else
    "Stage 10E may be considered only after independent acceptance of this completed stage; Stage 10F remains unauthorized here."
)
writeLines(decision_lines, file.path(out_dir, "STAGE10D_TECH_DECISION.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(out_dir, "STAGE10D_TECH_SESSIONINFO.txt"))

summary_lines <- c(
  "# Stage 10D-TECH server summary",
  "",
  paste0("- Decision: ", decision),
  paste0("- Reference: GSE161277, ", length(reference_donors), " donors; pre-existing M02-excluded broad annotations"),
  paste0("- Frozen primary set: 30/36; SHA256 ", expected$common_geneset),
  paste0("- Simulations: ", n_sim, " replicates x ", n_regions, " regions"),
  paste0("- Primary Pearson/Spearman: ", format(metric_value("pearson"), digits = 4), " / ",
         format(metric_value("spearman"), digits = 4)),
  paste0("- RMSE/failure rate: ", format(metric_value("standardized_rmse"), digits = 4), " / ",
         format(metric_value("failure_rate"), digits = 4)),
  paste0("- Failed strict criteria: ", if (length(failed_strict)) paste(failed_strict, collapse = ", ") else "none"),
  "- Spatial lesion-normal data accessed: NO",
  "- Next: independent output/hash acceptance; do not start Stage 10E or 10F automatically."
)
if (!smoke) {
  dir.create(file.path(root, "reports"), recursive = TRUE, showWarnings = FALSE)
  writeLines(summary_lines, file.path(root, "reports", "STAGE10D_TECH_SUMMARY.md"), useBytes = TRUE)
  writeLines(decision_lines, file.path(root, "reports", "STAGE10D_TECH_GATE_DECISION.md"), useBytes = TRUE)
}

script_path <- file.path(root, "scripts", "10D_TECH_pseudospot_validation.R")
manifest_inputs <- data.table(
  artifact_type = "input",
  path = c(parameter_path, design_path, mapped_path, lock36_path, lock35_path, common_path,
           decision_sp_path, reference_path, script_path),
  sha256 = NA_character_, bytes = NA_real_
)
manifest_inputs[, path := normalizePath(path, mustWork = TRUE)]
manifest_inputs[, `:=`(sha256 = vapply(path, hash_file, character(1)), bytes = file.info(path)$size)]

if (!smoke) {
  output_paths <- c(
    file.path(out_dir, c("STAGE10D_TECH_PRIMARY_SCORE_RESULTS.tsv", "STAGE10D_TECH_COVERAGE_RESULTS.tsv",
                         "STAGE10D_TECH_SENSITIVITY_RESULTS.tsv", "STAGE10D_TECH_NULL_MODULE_RESULTS.tsv",
                         "STAGE10D_TECH_DECISION.md", "STAGE10D_TECH_SESSIONINFO.txt",
                         "STAGE10D_TECH_SIMULATION_DESIGN.tsv")),
    list.files(fig_dir, recursive = TRUE, full.names = TRUE),
    file.path(root, "reports", c("STAGE10D_TECH_SUMMARY.md", "STAGE10D_TECH_GATE_DECISION.md"))
  )
  output_paths <- unique(output_paths[file.exists(output_paths)])
  manifest_outputs <- data.table(artifact_type = "output", path = normalizePath(output_paths),
                                 sha256 = vapply(output_paths, hash_file, character(1)),
                                 bytes = file.info(output_paths)$size)
  manifest <- rbind(manifest_inputs, manifest_outputs, fill = TRUE)
  root_norm <- normalizePath(root)
  manifest[, path := vapply(path, function(p) {
    if (startsWith(p, root_norm)) sub("^[\\\\/]", "", substring(p, nchar(root_norm) + 1L)) else p
  }, character(1))]
  fwrite(manifest, file.path(out_dir, "STAGE10D_TECH_RUN_MANIFEST.tsv"), sep = "\t", na = "NA")
}

cat(if (smoke) "STAGE10D_TECH_SMOKE_OK\n" else paste0("STAGE10D_TECH_COMPLETE decision=", decision, "\n"))
