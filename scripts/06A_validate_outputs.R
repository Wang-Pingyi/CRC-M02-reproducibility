#!/usr/bin/env Rscript

# Independent output validation for Stage 6A.
# Random seed: 20260728

set.seed(20260728)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
result_dir <- file.path(project_dir, "results", "06A_pseudobulk")
figure_dir <- file.path(project_dir, "figures", "06A_pseudobulk")

read_tsv <- function(name) {
  path <- file.path(result_dir, name)
  if (!file.exists(path)) stop("Missing output: ", path)
  utils::read.delim(path, check.names = FALSE)
}

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

results <- read_tsv("pseudobulk_results.tsv")
manifest <- read_tsv("pseudobulk_sample_manifest.tsv")
fit_audit <- read_tsv("model_fit_audit.tsv")
contrasts <- read_tsv("contrast_definitions.tsv")
candidate_programs <- read_tsv("candidate_programs.tsv")
lodo <- read_tsv("leave_one_donor_out_summary.tsv")
covariates <- read_tsv("model_covariate_audit.tsv")
metrics <- read_tsv("stage_6A_key_metrics.tsv")

required_result_columns <- c(
  "gene", "epithelial_state", "contrast", "log2FC", "CI95_low", "CI95_high",
  "p_value", "FDR", "n_normal_donors", "n_adenoma_donors", "n_cancer_donors",
  "model_formula"
)
add_check(
  "required_result_columns",
  all(required_result_columns %in% colnames(results)),
  paste(setdiff(required_result_columns, colnames(results)), collapse = ";")
)
add_check("results_nonempty", nrow(results) > 0L, nrow(results))
add_check(
  "three_primary_contrasts",
  setequal(
    unique(results$contrast),
    c("adenoma_vs_normal", "cancer_vs_adenoma", "cancer_vs_normal")
  ),
  paste(sort(unique(results$contrast)), collapse = ";")
)
add_check(
  "p_values_valid",
  all(is.finite(results$p_value) & results$p_value >= 0 & results$p_value <= 1),
  "all exact P values are finite and within [0,1]"
)
add_check(
  "fdr_valid",
  all(is.finite(results$FDR) & results$FDR >= 0 & results$FDR <= 1),
  "all BH FDR values are finite and within [0,1]"
)
add_check(
  "confidence_intervals_ordered",
  all(results$CI95_low <= results$log2FC & results$log2FC <= results$CI95_high),
  "all log2FC estimates lie within their 95% CI"
)
add_check(
  "minimum_donors",
  all(results$n_normal_donors >= 3 & results$n_adenoma_donors >= 3 & results$n_cancer_donors >= 3),
  "all reported contrasts have at least three donors per stage"
)
add_check(
  "unique_pseudobulk_ids",
  !anyDuplicated(manifest$pseudobulk_id),
  nrow(manifest)
)
add_check(
  "aggregation_unit_columns",
  all(c("donor_id", "stage", "epithelial_state") %in% colnames(manifest)),
  "donor × stage × epithelial state retained"
)
add_check(
  "raw_counts_object_exists",
  file.exists(file.path(project_dir, "objects", "GSE201348_6A_epithelial_pseudobulk.rds")),
  "pseudobulk count object"
)
add_check(
  "covariate_audit_complete",
  all(c("donor_id", "FAP_status", "colon_or_rectum", "tumor_location") %in% covariates$variable),
  paste(covariates$variable, collapse = ";")
)
add_check(
  "fit_audit_has_success",
  any(fit_audit$status == "ok"),
  paste(fit_audit$epithelial_state[fit_audit$status == "ok"], collapse = ";")
)
add_check(
  "contrast_direction_documented",
  all(nzchar(contrasts$interpretation)),
  "positive log2FC direction documented"
)
add_check(
  "report_exists",
  file.exists(file.path(project_dir, "reports", "stage_6A_pseudobulk_discovery.md")),
  "stage report generated"
)
add_check(
  "software_versions_exist",
  file.exists(file.path(result_dir, "software_versions.tsv")),
  "R and package versions recorded"
)
add_check(
  "candidate_program_schema",
  nrow(candidate_programs) == 0L ||
    all(c("program_id", "gene", "epithelial_state", "direction", "passes_lodo") %in%
      colnames(candidate_programs)),
  nrow(candidate_programs)
)
add_check(
  "candidate_programs_lodo_stable",
  nrow(candidate_programs) == 0L || all(candidate_programs$passes_lodo),
  "no LODO-unstable gene promoted"
)
add_check(
  "lodo_audit_written",
  file.exists(file.path(result_dir, "leave_one_donor_out_results.tsv")),
  nrow(lodo)
)
add_check(
  "no_single_cell_primary_test",
  !any(grepl("single.cell|wilcox|FindMarkers", results$model_formula, ignore.case = TRUE)),
  "primary model formulas are pseudobulk limma-voom"
)
add_check(
  "no_copykat_selection",
  !any(grepl("copykat|malignan", colnames(manifest), ignore.case = TRUE)),
  "malignancy labels absent from aggregation manifest"
)
add_check(
  "metrics_written",
  all(c(
    "input_epithelial_cells", "donors", "differential_result_rows",
    "candidate_programs"
  ) %in% metrics$metric),
  paste(metrics$metric, collapse = ";")
)

figure_bases <- c("candidate_gene_counts_by_state", "candidate_program_scores")
for (base in figure_bases) {
  pdf_path <- file.path(figure_dir, paste0(base, ".pdf"))
  png_path <- file.path(figure_dir, paste0(base, ".png"))
  source_name <- if (base == "candidate_gene_counts_by_state") {
    "candidate_gene_counts_by_state.tsv"
  } else {
    "candidate_program_scores.tsv"
  }
  source_path <- file.path(result_dir, "source_data", source_name)
  figures_expected <- if (base == "candidate_gene_counts_by_state") {
    as.numeric(metrics$value[metrics$metric == "lodo_stable_candidate_genes"]) > 0
  } else {
    as.numeric(metrics$value[metrics$metric == "candidate_programs"]) > 0
  }
  add_check(
    paste0(base, "_files"),
    !figures_expected ||
      (file.exists(pdf_path) && file.info(pdf_path)$size > 0 &&
        file.exists(png_path) && file.info(png_path)$size > 0 &&
        file.exists(source_path) && file.info(source_path)$size > 0),
    if (figures_expected) "PDF, PNG and source data required" else "not expected: zero candidates"
  )
}

checks <- do.call(rbind, checks)
validation_path <- file.path(result_dir, "validation_checks.tsv")
utils::write.table(
  checks, validation_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
cat("Stage 6A validation:", sum(checks$passed), "/", nrow(checks), "checks passed\n")
if (!all(checks$passed)) {
  print(checks[!checks$passed, , drop = FALSE])
  quit(status = 1L)
}
quit(status = 0L)
