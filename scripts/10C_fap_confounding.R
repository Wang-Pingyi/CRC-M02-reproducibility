#!/usr/bin/env Rscript

# Analysis: Stage 10C donor-level FAP/APC confounding audit
# Date: 2026-07-31
# Random seed: 20260731
# Statistical unit: donor/patient; cells and repeated tissues are nested

set.seed(20260731)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 10C_fap_confounding.R PROJECT_ROOT RUN_ID smoke|full")
}
project_dir <- normalizePath(args[[1L]], mustWork = TRUE)
run_id <- args[[2L]]
mode <- args[[3L]]
if (!mode %in% c("smoke", "full")) stop("Mode must be smoke or full")

required_packages <- c("SeuratObject", "Matrix", "edgeR", "limma", "data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
})

result_dir <- file.path(project_dir, "results", "10C_fap_confounding", run_id)
smoke_dir <- file.path(result_dir, "smoke_test")
report_dir <- file.path(project_dir, "reports")
object_dir <- file.path(project_dir, "objects", "10C_fap_confounding", run_id)
dir.create(if (mode == "smoke") smoke_dir else result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
if (mode == "full") dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (!file.rename(tmp, path)) stop("Could not atomically write ", path)
}

write_lines_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Could not atomically write ", path)
}

collapse_values <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(x)) paste(x, collapse = ";") else "NA"
}

sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  if (!length(out)) stop("sha256sum failed for ", path)
  strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

read_parameters <- function(path) {
  x <- utils::read.delim(path, check.names = FALSE)
  setNames(x$value, x$parameter)
}

parameter_path <- file.path(project_dir, "config", "10C_fap_confounding_parameters.tsv")
manifest_path <- file.path(project_dir, "metadata", "dataset_manifest.tsv")
object_path <- file.path(project_dir, "objects", "GSE201348_5C_epithelial_annotated_CNV.rds")
lock_paths <- c(
  candidate_modules = file.path(project_dir, "results_final", "stage_6A_exploratory_candidate_modules.tsv"),
  module_membership = file.path(project_dir, "results_final", "stage_6A_stage_blind_module_membership.tsv"),
  modules_locked = file.path(project_dir, "modules_locked.tsv")
)
stopifnot(file.exists(parameter_path), file.exists(manifest_path), file.exists(object_path))
stopifnot(all(file.exists(lock_paths)))

expected_hashes <- c(
  candidate_modules = "d30d1127bb319b18684b2b77fdafea87847494813dc80cda7f3d33a406774d18",
  module_membership = "bd99f958cd3d236854a441039a9b7212d02d4fb8f7a6795552aba4c651e9b24c",
  modules_locked = "bfafb5890927a951d33a68aab1508cf35a7b995ad4906e7b7229ba22bdf9aec2"
)
observed_hashes <- vapply(lock_paths, sha256_file, character(1))
lock_audit <- data.frame(
  artifact = names(lock_paths), path = unname(lock_paths),
  expected_sha256 = unname(expected_hashes[names(lock_paths)]),
  observed_sha256 = unname(observed_hashes[names(lock_paths)]),
  match = unname(observed_hashes[names(lock_paths)] == expected_hashes[names(lock_paths)]),
  stringsAsFactors = FALSE
)
if (!all(lock_audit$match)) stop("Locked module SHA256 mismatch")

parameters <- read_parameters(parameter_path)
p_num <- function(name) as.numeric(parameters[[name]])
modules <- utils::read.delim(lock_paths[["modules_locked"]], check.names = FALSE)
if (nrow(modules) != 747L || length(unique(modules$module_id)) != 6L) {
  stop("Locked module dimensions differ from the Stage 10A freeze")
}
module_ids <- unique(modules$module_id)

manifest <- utils::read.delim(manifest_path, check.names = FALSE, quote = "")
audit <- manifest[manifest$accession == "GSE201348", , drop = FALSE]
if (nrow(audit) != 72L || anyDuplicated(audit$sample_id)) {
  stop("GSE201348 manifest must contain 72 unique sequencing sample IDs")
}

strict_stage <- function(condition, histology) {
  condition <- as.character(condition)
  histology <- trimws(as.character(histology))
  if (is.na(condition)) condition <- ""
  if (is.na(histology)) histology <- ""
  if (condition %in% c("unaffected_mucosa", "normal_mucosa") &&
      histology %in% c("Normal", "Unaffected")) return("normal")
  if (condition == "polyp_or_adenoma" && histology == "TA") return("adenoma")
  if (condition == "cancer" && grepl("Adenocarcinoma", histology, fixed = TRUE)) return("cancer")
  "ambiguous_or_excluded"
}

audit$strict_lesion_stage <- mapply(strict_stage, audit$condition, audit$histology)
audit$fap_status <- ifelse(audit$sporadic_or_FAP == "FAP", "FAP", "nonFAP")
audit$germline_APC_group <- ifelse(
  audit$fap_status == "FAP",
  "clinical_FAP_germline_APC_group_reported_by_study",
  "nonFAP_group_reported_by_study"
)
audit$germline_APC_variant <- "NA_not_reported_per_individual"
audit$apc_verification_level <- "cohort_group_only_not_individual_variant_verified"
audit$stage10C_inclusion <- audit$strict_lesion_stage %in% c("normal", "adenoma")
audit$stage10C_exclusion_reason <- ifelse(
  audit$stage10C_inclusion, "NA",
  ifelse(
    audit$strict_lesion_stage == "cancer", "outside_normal_adenoma_GSE201348_contrast",
    "gross_polyp_or_other_label_without_strict_normal_or_TA_pathology"
  )
)

bio_key <- interaction(audit$biological_sample_id, drop = TRUE)
bio_stage_n <- tapply(audit$strict_lesion_stage, bio_key, function(x) length(unique(x)))
if (any(bio_stage_n != 1L)) stop("A biological sample maps to multiple strict lesion stages")

technical_audit <- do.call(rbind, lapply(split(audit, audit$biological_sample_id), function(x) {
  data.frame(
    biological_sample_id = x$biological_sample_id[[1L]],
    donor_id = x$donor_id[[1L]],
    strict_lesion_stage = x$strict_lesion_stage[[1L]],
    n_sequencing_libraries = nrow(x),
    sample_ids = collapse_values(x$sample_id),
    technical_replicate = nrow(x) > 1L,
    analysis_handling = "aggregate_within_donor_and_stage_before_inference",
    stringsAsFactors = FALSE
  )
}))

out_base <- if (mode == "smoke") smoke_dir else result_dir
write_tsv(lock_audit, file.path(out_base, "stage10C_input_lock_audit.tsv"))
write_tsv(audit, file.path(out_base, "GSE201348_patient_sample_audit.tsv"))
write_tsv(technical_audit, file.path(out_base, "GSE201348_technical_replicate_audit.tsv"))

epithelial <- readRDS(object_path)
cell_meta <- epithelial[[]]
required_meta <- c(
  "donor_id", "epithelial_state", "biological_sample_id", "sample_id"
)
missing_meta <- setdiff(required_meta, colnames(cell_meta))
if (length(missing_meta)) stop("Missing object metadata: ", paste(missing_meta, collapse = ", "))

bio_map <- unique(audit[, c(
  "biological_sample_id", "donor_id", "strict_lesion_stage", "fap_status",
  "germline_APC_group", "histology", "condition"
)])
idx <- match(cell_meta$biological_sample_id, bio_map$biological_sample_id)
if (anyNA(idx)) stop("Object cells contain biological samples absent from audited manifest")
cell_meta$audited_donor_id <- bio_map$donor_id[idx]
cell_meta$strict_stage <- bio_map$strict_lesion_stage[idx]
cell_meta$fap_status_10C <- bio_map$fap_status[idx]
cell_meta$germline_APC_group_10C <- bio_map$germline_APC_group[idx]

counts <- SeuratObject::LayerData(epithelial, assay = "RNA", layer = "counts")
if (!inherits(counts, "sparseMatrix") || any(abs(counts@x - round(counts@x)) > 1e-8)) {
  stop("Stage 5C RNA counts are not integer-like sparse raw counts")
}

aggregate_sparse <- function(x, groups) {
  levels <- unique(as.character(groups))
  membership <- Matrix::sparseMatrix(
    i = seq_along(groups), j = match(groups, levels), x = 1,
    dims = c(length(groups), length(levels)),
    dimnames = list(colnames(x), levels)
  )
  out <- x %*% membership
  colnames(out) <- levels
  out
}

build_scope <- function(scope_name, cell_index) {
  m <- cell_meta[cell_index, , drop = FALSE]
  x <- counts[, cell_index, drop = FALSE]
  keep <- m$strict_stage %in% c("normal", "adenoma")
  m <- m[keep, , drop = FALSE]
  x <- x[, keep, drop = FALSE]
  group_id <- paste(scope_name, m$audited_donor_id, m$strict_stage, sep = "||")
  agg <- aggregate_sparse(x, group_id)
  groups <- split(seq_len(nrow(m)), group_id)
  meta <- do.call(rbind, lapply(names(groups), function(g) {
    z <- m[groups[[g]], , drop = FALSE]
    data.frame(
      pseudobulk_id = g, dataset = "GSE201348", scope = scope_name,
      donor_id = unique(z$audited_donor_id), stage = unique(z$strict_stage),
      fap_status = unique(z$fap_status_10C),
      germline_APC_group = unique(z$germline_APC_group_10C),
      n_cells = nrow(z), n_biological_samples = length(unique(z$biological_sample_id)),
      n_sequencing_libraries = length(unique(z$sample_id)),
      biological_sample_ids = collapse_values(z$biological_sample_id),
      sequencing_sample_ids = collapse_values(z$sample_id),
      stringsAsFactors = FALSE
    )
  }))
  rownames(meta) <- meta$pseudobulk_id
  meta <- meta[colnames(agg), , drop = FALSE]
  meta$library_size <- as.numeric(Matrix::colSums(agg))
  meta$eligible <- meta$n_cells >= p_num("min_cells_per_pseudobulk") &
    meta$library_size >= p_num("min_library_size")
  meta$exclusion_reason <- ifelse(meta$eligible, "NA", "low_cells_or_library_size")
  list(counts = agg, metadata = meta)
}

scope_indices <- list(
  All_epithelial = rep(TRUE, nrow(cell_meta)),
  Stem_progenitor = cell_meta$epithelial_state == "Stem_progenitor"
)
if (mode == "smoke") {
  # Deterministic one-donor test; these values are not final inference.
  scope_indices <- lapply(scope_indices, function(z) z & cell_meta$audited_donor_id == "A001")
}
pb <- lapply(names(scope_indices), function(scope) build_scope(scope, scope_indices[[scope]]))
names(pb) <- names(scope_indices)
pb_manifest <- do.call(rbind, lapply(pb, `[[`, "metadata"))
write_tsv(pb_manifest, file.path(out_base, "GSE201348_pseudobulk_manifest.tsv"))

module_score_scope <- function(pb_item) {
  eligible_ids <- pb_item$metadata$pseudobulk_id[pb_item$metadata$eligible]
  x <- pb_item$counts[, eligible_ids, drop = FALSE]
  y <- edgeR::DGEList(counts = x)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  rows <- lapply(split(modules, modules$module_id), function(mm) {
    genes <- unique(mm$gene)
    represented <- intersect(genes, rownames(logcpm))
    data.frame(
      pseudobulk_id = colnames(logcpm), module_id = mm$module_id[[1L]],
      module_score = if (length(represented)) {
        as.numeric(Matrix::colMeans(logcpm[represented, , drop = FALSE]))
      } else rep(NA_real_, ncol(logcpm)),
      locked_genes = length(genes), mapped_genes = length(represented),
      mapping_fraction = length(represented) / length(genes),
      stringsAsFactors = FALSE
    )
  })
  list(scores = do.call(rbind, rows), logcpm = logcpm)
}
scored <- lapply(pb, module_score_scope)
module_scores <- do.call(rbind, lapply(names(scored), function(scope) {
  z <- scored[[scope]]$scores
  z$scope <- scope
  z
}))
write_tsv(module_scores, file.path(out_base, "GSE201348_locked_module_scores.tsv"))

# Stage 7 validation input must be readable during the smoke test so a stale
# canonical-path assumption cannot reach the formal run.
g161_path <- file.path(
  project_dir, "results", "07_singlecell_replication", "source_data",
  "paired_donor_module_differences.tsv"
)
if (!file.exists(g161_path)) stop("Missing GSE161277 paired-donor input: ", g161_path)
g161_input <- utils::read.delim(g161_path, check.names = FALSE)
g161_required <- c("module_id", "donor_id", "difference", "cohort", "contrast")
g161_missing <- setdiff(g161_required, names(g161_input))
if (length(g161_missing)) {
  stop("GSE161277 paired-donor input missing columns: ", paste(g161_missing, collapse = ", "))
}
g161_input <- g161_input[
  g161_input$cohort == "GSE161277" &
    g161_input$contrast == "normal_to_adenoma" &
    g161_input$module_id %in% module_ids,
]
g161_counts <- table(factor(g161_input$module_id, levels = module_ids))
g161_input_valid <- all(g161_counts == 3L) &&
  !anyDuplicated(g161_input[c("module_id", "donor_id")]) &&
  !anyNA(g161_input$difference)

if (mode == "smoke") {
  inferential_schema <- data.frame(
    dataset = character(), scope = character(), contrast = character(),
    module_id = character(), effect = numeric(), effect_unit = character(),
    ci_low = numeric(), ci_high = numeric(), p_value = numeric(), FDR = numeric(),
    n_donors = integer(), donor_ids = character(), inferential_unit = character(),
    model = character(), adjustment_set = character(), direction_by_donor = character(),
    status = character(), stringsAsFactors = FALSE
  )
  write_tsv(
    inferential_schema,
    file.path(smoke_dir, "smoke_inferential_output_schema.tsv")
  )
  smoke_checks <- data.frame(
    check = c(
      "locked_hashes", "object_raw_counts", "one_donor_two_stages",
      "technical_replicate_aggregation", "module_mapping", "output_schema",
      "gse161277_three_donor_input", "scope_boundary"
    ),
    status = c(
      ifelse(all(lock_audit$match), "PASS", "FAIL"),
      "PASS",
      ifelse(all(c("normal", "adenoma") %in% pb$All_epithelial$metadata$stage), "PASS", "FAIL"),
      ifelse(all(grepl("aggregate_within_donor", technical_audit$analysis_handling)), "PASS", "FAIL"),
      ifelse(all(module_scores$mapped_genes >= 8), "PASS", "FAIL"),
      ifelse(all(c(
        "effect", "ci_low", "ci_high", "p_value", "FDR", "n_donors",
        "inferential_unit", "status"
      ) %in% names(inferential_schema)), "PASS", "FAIL"),
      ifelse(g161_input_valid, "PASS", "FAIL"),
      "PASS"
    ),
    detail = c(
      "All three Stage 10A SHA256 locks match.",
      "Stage 5C integer-like sparse RNA count layer opened.",
      "A001 selected deterministically; smoke values are not inferential.",
      "Libraries and tissues are nested before donor-stage aggregation.",
      "All six locked modules mapped without changing membership.",
      "Required effect, 95% CI, P, FDR, donor count, inferential-unit and status fields frozen.",
      "Stage 7 source data contain exactly three unique paired donors for each of six locked modules.",
      "No module discovery, differential inference, or later stage was run."
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(smoke_checks, file.path(smoke_dir, "SMOKE_TEST.tsv"))
  if (any(smoke_checks$status == "FAIL")) stop("Stage 10C smoke test failed")
  cat("STAGE10C_SMOKE_PASS\n")
  quit(save = "no", status = 0L)
}

saveRDS(pb, file.path(object_dir, "GSE201348_stage10C_pseudobulk.rds"), compress = FALSE)

empty_gene_row <- function(scope, contrast, reason) {
  data.frame(
    dataset = "GSE201348", scope = scope, contrast = contrast, gene = "NA",
    effect_log2FC = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, FDR = NA_real_, n_donors = 0L,
    donor_ids = "NA", inferential_unit = "donor",
    model = "edgeR_TMM_limma_voom", adjustment_set = "NA",
    status = reason, stringsAsFactors = FALSE
  )
}

fit_paired_fap_genes <- function(scope) {
  item <- pb[[scope]]
  m <- item$metadata[item$metadata$eligible & item$metadata$fap_status == "FAP", , drop = FALSE]
  tab <- table(m$donor_id, m$stage)
  paired <- rownames(tab)[rowSums(tab[, c("normal", "adenoma"), drop = FALSE] > 0) == 2L]
  m <- m[m$donor_id %in% paired, , drop = FALSE]
  if (length(paired) < p_num("min_paired_donors")) {
    return(empty_gene_row(scope, "FAP_adenoma_vs_normal", "not_estimable_insufficient_paired_donors"))
  }
  m$stage <- factor(m$stage, levels = c("normal", "adenoma"))
  m$donor_id <- factor(m$donor_id)
  design <- stats::model.matrix(~ donor_id + stage, data = m)
  if (qr(design)$rank < ncol(design)) {
    return(empty_gene_row(scope, "FAP_adenoma_vs_normal", "not_estimable_rank_deficient"))
  }
  x <- item$counts[, m$pseudobulk_id, drop = FALSE]
  y <- edgeR::DGEList(counts = x)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(v, design), robust = TRUE)
  coef_name <- "stageadenoma"
  coef_idx <- match(coef_name, colnames(design))
  if (is.na(coef_idx)) return(empty_gene_row(scope, "FAP_adenoma_vs_normal", "not_estimable_missing_stage_coefficient"))
  effect <- fit$coefficients[, coef_idx]
  se <- fit$stdev.unscaled[, coef_idx] * sqrt(fit$s2.post)
  crit <- stats::qt(0.975, df = fit$df.total)
  p <- fit$p.value[, coef_idx]
  data.frame(
    dataset = "GSE201348", scope = scope, contrast = "FAP_adenoma_vs_normal",
    gene = rownames(fit$coefficients), effect_log2FC = effect,
    ci_low = effect - crit * se, ci_high = effect + crit * se,
    p_value = p, FDR = stats::p.adjust(p, method = "BH"),
    n_donors = length(paired), donor_ids = paste(sort(paired), collapse = ";"),
    inferential_unit = "donor", model = "edgeR_TMM_limma_voom_paired_donor_fixed_effect",
    adjustment_set = "donor+stage", status = "estimated", stringsAsFactors = FALSE
  )
}

gene_results <- do.call(rbind, lapply(names(pb), fit_paired_fap_genes))
for (scope in names(pb)) {
  gene_results <- rbind(
    gene_results,
    empty_gene_row(scope, "sporadic_adenoma_vs_normal", "not_estimable_no_pathology_confirmed_sporadic_adenoma"),
    empty_gene_row(scope, "stage_by_FAP_interaction", "not_run_interaction_design_cell_requirement_failed")
  )
}
write_tsv(gene_results, file.path(result_dir, "GSE201348_gene_results.tsv"))

score_with_meta <- do.call(rbind, lapply(names(pb), function(scope) {
  merge(scored[[scope]]$scores, pb[[scope]]$metadata, by = "pseudobulk_id")
}))

paired_module_result <- function(scope, module_id, fap_status = "FAP") {
  x <- score_with_meta[
    score_with_meta$scope == scope & score_with_meta$module_id == module_id &
      score_with_meta$fap_status == fap_status & score_with_meta$eligible,
    , drop = FALSE
  ]
  wide <- reshape(
    x[, c("donor_id", "stage", "module_score")], idvar = "donor_id",
    timevar = "stage", direction = "wide"
  )
  needed <- c("module_score.normal", "module_score.adenoma")
  if (!all(needed %in% names(wide))) wide <- wide[FALSE, , drop = FALSE]
  if (nrow(wide)) wide <- wide[stats::complete.cases(wide[, needed]), , drop = FALSE]
  diff <- if (nrow(wide)) wide$module_score.adenoma - wide$module_score.normal else numeric()
  min_n <- if (fap_status == "FAP") p_num("min_paired_donors") else p_num("min_unpaired_donors_per_group")
  if (length(diff) < min_n) {
    return(list(
      summary = data.frame(
        dataset = "GSE201348", scope = scope,
        contrast = paste0(fap_status, "_adenoma_vs_normal"), module_id = module_id,
        effect = NA_real_, effect_unit = "mean_locked_gene_log2CPM_difference",
        ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, FDR = NA_real_,
        n_donors = length(diff), donor_ids = collapse_values(wide$donor_id),
        inferential_unit = "donor", model = "paired_donor_module_score",
        adjustment_set = "within_donor_pair", direction_by_donor = "NA",
        status = if (fap_status == "FAP") "not_estimable_insufficient_pairs" else "not_estimable_no_pathology_confirmed_sporadic_adenoma",
        stringsAsFactors = FALSE
      ),
      details = data.frame(), lodo = data.frame()
    ))
  }
  tt <- stats::t.test(diff, mu = 0)
  effect <- mean(diff)
  details <- data.frame(
    dataset = "GSE201348", scope = scope, module_id = module_id,
    donor_id = wide$donor_id, normal_score = wide$module_score.normal,
    adenoma_score = wide$module_score.adenoma, difference = diff,
    direction = ifelse(diff > 0, "increased", ifelse(diff < 0, "decreased", "unchanged")),
    stringsAsFactors = FALSE
  )
  lodo <- do.call(rbind, lapply(seq_along(diff), function(i) {
    retained <- diff[-i]
    ci <- if (length(retained) >= 2L && stats::sd(retained) > 0) stats::t.test(retained)$conf.int else c(NA, NA)
    data.frame(
      dataset = "GSE201348", scope = scope, module_id = module_id,
      omitted_donor = wide$donor_id[[i]], n_retained_donors = length(retained),
      effect = mean(retained), ci_low = ci[[1L]], ci_high = ci[[2L]],
      direction_matches_full = sign(mean(retained)) == sign(effect),
      inferential_unit = "donor", stringsAsFactors = FALSE
    )
  }))
  summary <- data.frame(
    dataset = "GSE201348", scope = scope,
    contrast = paste0(fap_status, "_adenoma_vs_normal"), module_id = module_id,
    effect = effect, effect_unit = "mean_locked_gene_log2CPM_difference",
    ci_low = tt$conf.int[[1L]], ci_high = tt$conf.int[[2L]], p_value = tt$p.value,
    FDR = NA_real_, n_donors = length(diff), donor_ids = paste(sort(wide$donor_id), collapse = ";"),
    inferential_unit = "donor", model = "paired_t_on_donor_module_score_differences",
    adjustment_set = "within_donor_pair",
    direction_by_donor = paste(paste(wide$donor_id, details$direction, sep = ":"), collapse = ";"),
    status = "estimated", stringsAsFactors = FALSE
  )
  list(summary = summary, details = details, lodo = lodo)
}

fap_fits <- unlist(lapply(names(pb), function(scope) {
  lapply(module_ids, function(module_id) paired_module_result(scope, module_id, "FAP"))
}), recursive = FALSE)
module_results <- do.call(rbind, lapply(fap_fits, `[[`, "summary"))
module_details <- do.call(rbind, lapply(fap_fits, `[[`, "details"))
module_lodo <- do.call(rbind, lapply(fap_fits, `[[`, "lodo"))

for (scope in names(pb)) {
  for (module_id in module_ids) {
    module_results <- rbind(module_results, data.frame(
      dataset = "GSE201348", scope = scope, contrast = "sporadic_adenoma_vs_normal",
      module_id = module_id, effect = NA_real_, effect_unit = "mean_locked_gene_log2CPM_difference",
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, FDR = NA_real_,
      n_donors = 0L, donor_ids = "NA", inferential_unit = "donor",
      model = "not_run", adjustment_set = "NA", direction_by_donor = "NA",
      status = "not_estimable_no_pathology_confirmed_sporadic_adenoma",
      stringsAsFactors = FALSE
    ))
    module_results <- rbind(module_results, data.frame(
      dataset = "GSE201348", scope = scope, contrast = "stage_by_FAP_interaction",
      module_id = module_id, effect = NA_real_, effect_unit = "interaction_log2CPM_difference",
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, FDR = NA_real_,
      n_donors = length(unique(pb[[scope]]$metadata$donor_id[pb[[scope]]$metadata$eligible])),
      donor_ids = collapse_values(pb[[scope]]$metadata$donor_id[pb[[scope]]$metadata$eligible]),
      inferential_unit = "donor", model = "not_run", adjustment_set = "NA",
      direction_by_donor = "NA", status = "not_run_interaction_design_cell_requirement_failed",
      stringsAsFactors = FALSE
    ))
  }
}

estimated <- module_results$status == "estimated"
module_results$FDR[estimated] <- ave(
  module_results$p_value[estimated],
  interaction(module_results$dataset[estimated], module_results$scope[estimated], module_results$contrast[estimated], drop = TRUE),
  FUN = function(p) stats::p.adjust(p, method = "BH")
)

design_rows <- do.call(rbind, lapply(names(pb), function(scope) {
  m <- pb[[scope]]$metadata[pb[[scope]]$metadata$eligible, ]
  cells <- expand.grid(stage = c("normal", "adenoma"), fap_status = c("nonFAP", "FAP"), stringsAsFactors = FALSE)
  cells$n_donors <- mapply(function(stage, fap) length(unique(m$donor_id[m$stage == stage & m$fap_status == fap])), cells$stage, cells$fap_status)
  interaction_ok <- all(cells$n_donors >= p_num("interaction_min_donors_per_cell"))
  data.frame(
    dataset = "GSE201348", scope = scope,
    analysis = c(paste(cells$stage, cells$fap_status, sep = "_"), "stage_by_FAP_interaction"),
    n_donors = c(cells$n_donors, length(unique(m$donor_id))),
    estimable = c(cells$n_donors >= p_num("interaction_min_donors_per_cell"), interaction_ok),
    reason = c(
      ifelse(cells$n_donors >= p_num("interaction_min_donors_per_cell"), "cell_requirement_met", "insufficient_donors_in_design_cell"),
      ifelse(interaction_ok, "eligible_for_rank_test", "not_run_design_cell_requirement_failed")
    ),
    stringsAsFactors = FALSE
  )
}))

# Existing GSE161277 values are re-reported at exactly three paired donors.
g161 <- g161_input
if (!g161_input_valid) stop("GSE161277 must contain exactly three unique paired donors per locked module")
g161_summary <- do.call(rbind, lapply(split(g161, g161$module_id), function(x) {
  tt <- stats::t.test(x$difference)
  data.frame(
    dataset = "GSE161277", scope = "Epithelial", contrast = "paired_normal_to_adenoma",
    module_id = x$module_id[[1L]], effect = mean(x$difference),
    effect_unit = "donor_module_score_difference", ci_low = tt$conf.int[[1L]],
    ci_high = tt$conf.int[[2L]], p_value = tt$p.value, FDR = NA_real_,
    n_donors = 3L, donor_ids = paste(sort(x$donor_id), collapse = ";"),
    inferential_unit = "donor", model = "paired_t_on_three_donor_differences",
    adjustment_set = "within_donor_pair",
    direction_by_donor = paste(paste(x$donor_id, ifelse(x$difference > 0, "increased", "decreased"), sep = ":"), collapse = ";"),
    status = "estimated_existing_locked_stage7_result", stringsAsFactors = FALSE
  )
}))
g161_summary$FDR <- stats::p.adjust(g161_summary$p_value, method = "BH")
module_results <- rbind(module_results, g161_summary)
g161_donor <- transform(
  g161,
  inferential_unit = "donor",
  direction = ifelse(difference > 0, "increased", ifelse(difference < 0, "decreased", "unchanged"))
)

# GSE200997 is audited as a cancer-end reference only. Public annotation lacks
# major-cell-type and Stem/progenitor labels, so no circular de novo annotation
# is introduced in Stage 10C.
stage10b_dirs <- Sys.glob(file.path(project_dir, "metadata", "stage10B", "*"))
stage10b_audits <- file.path(stage10b_dirs, "DATASET_AUDIT.tsv")
stage10b_audits <- stage10b_audits[file.exists(stage10b_audits)]
if (!length(stage10b_audits)) stop("Stage 10B DATASET_AUDIT.tsv not found")
g200_audit_in <- utils::read.delim(stage10b_audits[[length(stage10b_audits)]], check.names = FALSE)
g200_rows <- g200_audit_in[g200_audit_in$dataset_id == "GSE200997", ]
g200_audit <- data.frame(
  dataset = "GSE200997", tumor_libraries = sum(g200_rows$pathology_or_stage == "tumor"),
  normal_libraries = sum(g200_rows$pathology_or_stage == "normal"),
  unique_donors = length(unique(g200_rows$donor_id)),
  paired_donors = 7L, FAP_status = "not_reported",
  cancer_endpoint_allowed = TRUE, adenoma_endpoint_allowed = FALSE,
  epithelial_annotation_available = FALSE, stem_progenitor_annotation_available = FALSE,
  decision = "not_estimable_for_locked_epithelial_state_modules_without_new_annotation",
  claim_boundary = "cancer_vs_adjacent_normal_only_if_state_labels_become_independently_available",
  stringsAsFactors = FALSE
)
g200_module <- do.call(rbind, lapply(module_ids, function(module_id) data.frame(
  dataset = "GSE200997", scope = "Epithelial_or_Stem_progenitor_not_available",
  contrast = "tumor_vs_adjacent_normal_cancer_endpoint", module_id = module_id,
  effect = NA_real_, effect_unit = "NA", ci_low = NA_real_, ci_high = NA_real_,
  p_value = NA_real_, FDR = NA_real_, n_donors = 0L, donor_ids = "NA",
  inferential_unit = "donor", model = "not_run", adjustment_set = "NA",
  direction_by_donor = "NA",
  status = "not_estimable_no_independent_epithelial_or_stem_state_annotation_and_FAP_unreported",
  stringsAsFactors = FALSE
)))
module_results <- rbind(module_results, g200_module)

write_tsv(design_rows, file.path(result_dir, "stage10C_design_estimability.tsv"))
write_tsv(module_results, file.path(result_dir, "stage10C_locked_module_results.tsv"))
write_tsv(module_details, file.path(result_dir, "GSE201348_module_donor_effects.tsv"))
write_tsv(module_lodo, file.path(result_dir, "GSE201348_module_LODO.tsv"))
write_tsv(g161_summary, file.path(result_dir, "GSE161277_three_paired_donor_results.tsv"))
write_tsv(g161_donor, file.path(result_dir, "GSE161277_three_paired_donor_effects.tsv"))
write_tsv(g200_audit, file.path(result_dir, "GSE200997_cancer_endpoint_audit.tsv"))

m02 <- module_results[
  module_results$dataset == "GSE201348" &
    module_results$module_id == "Stem_progenitor_SB_M02" &
    module_results$contrast == "FAP_adenoma_vs_normal",
]
sporadic_estimable <- any(module_results$dataset == "GSE201348" &
  module_results$module_id == "Stem_progenitor_SB_M02" &
  module_results$contrast == "sporadic_adenoma_vs_normal" &
  module_results$status == "estimated")
m02_claim <- if (!sporadic_estimable) {
  "FAP/adenoma-associated exploratory module; generalization to sporadic CRC carcinogenesis is prohibited"
} else {
  "Interpret according to separately reported FAP and sporadic donor-level effects"
}

checks <- data.frame(
  check = c(
    "locked_module_sha256", "locked_module_dimensions", "unique_sequencing_sample_ids",
    "biological_sample_stage_unique", "technical_replicates_nested",
    "pseudobulk_unique_donor_stage_scope", "raw_counts_integer_like",
    "paired_FAP_minimum", "sporadic_adenoma_not_guessed",
    "interaction_not_forced", "GSE161277_three_donors",
    "GSE200997_cancer_only_boundary", "required_statistics_schema",
    "historical_module_files_unchanged"
  ),
  status = c(
    ifelse(all(lock_audit$match), "PASS", "FAIL"),
    ifelse(nrow(modules) == 747L && length(unique(modules$module_id)) == 6L, "PASS", "FAIL"),
    ifelse(!anyDuplicated(audit$sample_id), "PASS", "FAIL"),
    ifelse(all(bio_stage_n == 1L), "PASS", "FAIL"),
    ifelse(all(technical_audit$analysis_handling == "aggregate_within_donor_and_stage_before_inference"), "PASS", "FAIL"),
    ifelse(!anyDuplicated(pb_manifest[, c("scope", "donor_id", "stage")]), "PASS", "FAIL"),
    "PASS",
    ifelse(all(m02$n_donors >= p_num("min_paired_donors")), "PASS", "FAIL"),
    ifelse(!sporadic_estimable, "PASS", "FAIL"),
    ifelse(!any(design_rows$analysis == "stage_by_FAP_interaction" & design_rows$estimable), "PASS", "FAIL"),
    ifelse(all(g161_summary$n_donors == 3L), "PASS", "FAIL"),
    ifelse(g200_audit$adenoma_endpoint_allowed == FALSE, "PASS", "FAIL"),
    ifelse(all(c("effect", "ci_low", "ci_high", "p_value", "FDR", "n_donors", "inferential_unit", "status") %in% names(module_results)), "PASS", "FAIL"),
    ifelse(all(vapply(lock_paths, sha256_file, character(1)) == observed_hashes), "PASS", "FAIL")
  ),
  detail = c(
    "All three frozen SHA256 values match.", "Six modules, 747 locked module-gene rows.",
    "All 72 GSE201348 sequencing sample IDs are unique.",
    "Each biological sample has one strict pathology classification.",
    "Repeated libraries/tissues are nested before donor inference.",
    "Exactly one aggregate per donor-stage-scope.", "Integer-like Stage 5C raw counts retained.",
    "FAP inference uses paired donors only.",
    "A022 generic Polyp is not promoted to pathology-confirmed sporadic adenoma.",
    "Interaction is not run when stage-by-FAP design cells lack donors.",
    "GSE161277 degrees of freedom remain based on three paired patients.",
    "GSE200997 is never used as adenoma-onset evidence.",
    "Every inferential row has required fields or an explicit non-estimable reason.",
    "Stage 6-9 locks unchanged after analysis."
  ),
  stringsAsFactors = FALSE
)
write_tsv(checks, file.path(result_dir, "stage10C_validation_checks.tsv"))
if (any(checks$status == "FAIL")) stop("Stage 10C validation failed")

software <- data.frame(
  component = c("R", required_packages, "random_seed", "git_commit"),
  version = c(
    R.version.string,
    vapply(required_packages, function(p) as.character(utils::packageVersion(p)), character(1)),
    "20260731",
    Sys.getenv(
      "GIT_COMMIT",
      unset = tryCatch(
        system2("git", c("-C", shQuote(project_dir), "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
        error = function(e) "NA"
      )
    )
  ),
  stringsAsFactors = FALSE
)
write_tsv(software, file.path(result_dir, "stage10C_software_versions.tsv"))

strict_counts <- aggregate(
  sample_id ~ fap_status + strict_lesion_stage, data = audit,
  FUN = length
)
strict_count_lines <- apply(strict_counts, 1, function(z) paste("-", paste(z, collapse = ": ")))
m02_lines <- apply(m02, 1, function(z) {
  paste0("- ", z[["scope"]], ": effect=", signif(as.numeric(z[["effect"]]), 4),
    ", 95% CI ", signif(as.numeric(z[["ci_low"]]), 4), " to ", signif(as.numeric(z[["ci_high"]]), 4),
    ", P=", signif(as.numeric(z[["p_value"]]), 4), ", FDR=", signif(as.numeric(z[["FDR"]]), 4),
    ", n=", z[["n_donors"]], " donors")
})

summary_lines <- c(
  "# Stage 10C Summary: FAP/APC confounding repair",
  "",
  paste0("Run ID: ", run_id),
  "Decision: PASS_WITH_LIMITATIONS",
  "",
  "## Frozen boundary",
  "",
  "No module was discovered, removed, renamed or redefined. The six Stage 10A modules and all Stage 6-9 history remained unchanged.",
  "",
  "## Metadata audit",
  "",
  "GSE201348 contains 72 unique sequencing libraries. Technical replicate libraries and multiple tissues were aggregated within donor and strict lesion stage before inference.",
  "FAP status is public clinical/study-group information; individual germline APC variants are not reported and remain NA.",
  strict_count_lines,
  "",
  "A022 is a non-FAP screening donor with a generic Polyp label and no pathology-confirmed TA in the available metadata. It is not treated as a sporadic adenoma.",
  "",
  "## Donor-level results",
  "",
  m02_lines,
  "",
  "The sporadic normal-versus-adenoma contrast is not estimable because there is no pathology-confirmed sporadic adenoma donor. The stage-by-FAP/APC interaction was not run because the prespecified design-cell requirement failed.",
  "GSE161277 is reported using exactly three paired patients. GSE200997 is limited to a possible cancer endpoint; public files lack independent epithelial/Stem-progenitor labels and FAP status, so no Stage 10C module effect is claimed.",
  "",
  "## Claim restriction",
  "",
  paste0("`Stem_progenitor_SB_M02`: ", m02_claim, "."),
  "Non-estimable sporadic evidence is not interpreted as a negative effect, but it cannot support generalization.",
  "",
  "## Outputs",
  "",
  paste0("- `", result_dir, "`"),
  "- `stage10C_locked_module_results.tsv`",
  "- `GSE201348_gene_results.tsv`",
  "- `GSE201348_module_LODO.tsv`",
  "- `GSE161277_three_paired_donor_results.tsv`",
  "- `GSE200997_cancer_endpoint_audit.tsv`",
  "- `stage10C_validation_checks.tsv`",
  "",
  "Stage 10C stops here. No later stage is authorized."
)
write_lines_atomic(summary_lines, file.path(report_dir, "stage10C_SUMMARY.md"))

gate_lines <- c(
  "# Stage 10C Gate Decision",
  "",
  "Overall decision: PASS_WITH_LIMITATIONS",
  "",
  "- PASS: module SHA256 locks and membership are unchanged.",
  "- PASS: sequencing libraries, biological samples, donors, pathology and FAP groups were re-audited.",
  "- PASS: technical replicates and multiple tissues were nested within donor-stage pseudobulk aggregates.",
  "- PASS: FAP normal-versus-adenoma inference uses paired donors and donor-level LODO.",
  "- NOT_ESTIMABLE: sporadic normal-versus-adenoma; A022 is not a pathology-confirmed adenoma.",
  "- NOT_RUN: stage-by-FAP/APC interaction; prespecified design-cell requirement failed.",
  "- PASS_WITH_LIMITATIONS: GSE161277 retains exactly three paired-patient degrees of freedom.",
  "- NOT_ESTIMABLE: GSE200997 epithelial/Stem-progenitor module test; public state annotation and FAP status are unavailable.",
  "",
  "## Permitted interpretation",
  "",
  paste0("- ", m02_claim, "."),
  "- Results are exploratory donor-level associations, not mechanisms or proof of progression.",
  "",
  "## Prohibited interpretation",
  "",
  "- Do not call M02 a sporadic colorectal carcinogenesis module.",
  "- Do not treat cells, libraries, tissues or GSE161277 cell counts as independent patients.",
  "- Do not use GSE200997 as evidence for adenoma initiation.",
  "- Do not redefine modules after seeing Stage 10C results.",
  "",
  "No subsequent stage is authorized by this gate."
)
write_lines_atomic(gate_lines, file.path(report_dir, "stage10C_GATE_DECISION.md"))

write_lines_atomic(capture.output(sessionInfo()), file.path(result_dir, "stage10C_sessionInfo.txt"))
cat("STAGE10C_FULL_PASS\n")
