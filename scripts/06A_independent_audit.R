#!/usr/bin/env Rscript

# Independent Stage 6A statistical audit.
# This script does not refit models or change thresholds.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1L]] else getwd()
result_dir <- file.path(project_dir, "results", "06A_pseudobulk")

results <- utils::read.delim(
  file.path(result_dir, "pseudobulk_results.tsv"),
  check.names = FALSE
)
parameters <- utils::read.delim(
  file.path(project_dir, "config", "06A_pseudobulk_parameters.tsv"),
  check.names = FALSE
)
param <- setNames(parameters$value, parameters$parameter)
p_num <- function(name) as.numeric(param[[name]])

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

states <- sort(unique(results$epithelial_state))
attrition_rows <- list()
top_rows <- list()

for (state_name in states) {
  state_results <- results[results$epithelial_state == state_name, ]
  wide <- reshape(
    state_results[, c("gene", "contrast", "log2FC", "CI95_low", "CI95_high", "p_value", "FDR")],
    idvar = "gene",
    timevar = "contrast",
    direction = "wide"
  )
  required <- c(
    "log2FC.adenoma_vs_normal", "FDR.adenoma_vs_normal",
    "log2FC.cancer_vs_adenoma", "FDR.cancer_vs_adenoma",
    "log2FC.cancer_vs_normal", "FDR.cancer_vs_normal"
  )
  wide <- wide[stats::complete.cases(wide[, required]), ]
  direction <- sign(wide$log2FC.adenoma_vs_normal)

  step_1 <- wide$FDR.adenoma_vs_normal <= p_num("fdr_threshold")
  step_2 <- step_1 &
    abs(wide$log2FC.adenoma_vs_normal) >= p_num("min_abs_log2fc")
  step_3 <- step_2 &
    sign(wide$log2FC.cancer_vs_normal) == direction
  step_4 <- step_3 &
    abs(wide$log2FC.cancer_vs_normal) >=
      p_num("sustained_fraction") * abs(wide$log2FC.adenoma_vs_normal)
  step_5 <- step_4 &
    (
      sign(wide$log2FC.cancer_vs_adenoma) == direction |
        abs(wide$log2FC.cancer_vs_adenoma) <= p_num("plateau_tolerance_log2fc")
    )

  attrition_rows[[state_name]] <- data.frame(
    epithelial_state = state_name,
    genes_tested = nrow(wide),
    min_p_adenoma_vs_normal = min(wide$p_value.adenoma_vs_normal),
    min_FDR_adenoma_vs_normal = min(wide$FDR.adenoma_vs_normal),
    adenoma_vs_normal_nominal_p_0_05 = sum(wide$p_value.adenoma_vs_normal <= 0.05),
    adenoma_vs_normal_FDR = sum(step_1),
    adenoma_vs_normal_FDR_0_10 = sum(wide$FDR.adenoma_vs_normal <= 0.10),
    plus_abs_log2FC = sum(step_2),
    plus_cancer_vs_normal_same_direction = sum(step_3),
    plus_cancer_vs_normal_sustained_magnitude = sum(step_4),
    plus_cancer_vs_adenoma_nonopposing_or_plateau = sum(step_5),
    cancer_vs_adenoma_FDR = sum(wide$FDR.cancer_vs_adenoma <= p_num("fdr_threshold")),
    cancer_vs_normal_FDR = sum(wide$FDR.cancer_vs_normal <= p_num("fdr_threshold")),
    stringsAsFactors = FALSE
  )

  ranked <- wide[order(wide$FDR.adenoma_vs_normal, -abs(wide$log2FC.adenoma_vs_normal)), ]
  ranked <- head(ranked, 20L)
  ranked$epithelial_state <- state_name
  ranked$passes_adenoma_FDR_and_effect <- with(
    ranked,
    FDR.adenoma_vs_normal <= p_num("fdr_threshold") &
      abs(log2FC.adenoma_vs_normal) >= p_num("min_abs_log2fc")
  )
  ranked$cancer_vs_normal_same_direction <- with(
    ranked,
    sign(log2FC.cancer_vs_normal) == sign(log2FC.adenoma_vs_normal)
  )
  ranked$cancer_vs_adenoma_nonopposing_or_plateau <- with(
    ranked,
    sign(log2FC.cancer_vs_adenoma) == sign(log2FC.adenoma_vs_normal) |
      abs(log2FC.cancer_vs_adenoma) <= p_num("plateau_tolerance_log2fc")
  )
  top_rows[[state_name]] <- ranked
}

attrition <- do.call(rbind, attrition_rows)
top_genes <- do.call(rbind, top_rows)
write_tsv(attrition, file.path(result_dir, "independent_candidate_attrition_audit.tsv"))
write_tsv(top_genes, file.path(result_dir, "independent_top_adenoma_changes.tsv"))

formula_audit <- unique(results[, c(
  "epithelial_state", "n_normal_donors", "n_adenoma_donors",
  "n_cancer_donors", "donor_correlation", "model_formula"
)])
formula_audit <- formula_audit[order(formula_audit$epithelial_state), ]
write_tsv(formula_audit, file.path(result_dir, "independent_model_formula_audit.tsv"))

checks <- data.frame(
  check = c(
    "all_three_contrasts_present_per_state",
    "unique_gene_state_contrast_rows",
    "confidence_intervals_finite_and_ordered",
    "p_values_and_fdr_valid",
    "minimum_three_donors_per_stage",
    "candidate_attrition_reproduces_zero"
  ),
  passed = c(
    all(vapply(
      split(results$contrast, results$epithelial_state),
      function(x) setequal(unique(x), c(
        "adenoma_vs_normal", "cancer_vs_adenoma", "cancer_vs_normal"
      )),
      logical(1)
    )),
    !anyDuplicated(results[, c("gene", "epithelial_state", "contrast")]),
    all(
      is.finite(results$CI95_low) & is.finite(results$CI95_high) &
        results$CI95_low <= results$log2FC &
        results$log2FC <= results$CI95_high
    ),
    all(
      is.finite(results$p_value) & results$p_value >= 0 & results$p_value <= 1 &
        is.finite(results$FDR) & results$FDR >= 0 & results$FDR <= 1
    ),
    all(
      results$n_normal_donors >= 3 &
        results$n_adenoma_donors >= 3 &
        results$n_cancer_donors >= 3
    ),
    sum(attrition$plus_cancer_vs_adenoma_nonopposing_or_plateau) == 0
  ),
  stringsAsFactors = FALSE
)

pseudobulk_object <- readRDS(file.path(
  project_dir, "objects", "GSE201348_6A_epithelial_pseudobulk.rds"
))
epithelial_object <- readRDS(file.path(
  project_dir, "objects", "GSE201348_5C_epithelial_annotated_CNV.rds"
))
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop("SeuratObject is required for raw-count conservation audit")
}
original_counts <- SeuratObject::LayerData(
  epithelial_object, assay = "RNA", layer = "counts"
)
pseudobulk_counts <- pseudobulk_object$counts
pseudobulk_meta <- pseudobulk_object$metadata

conservation_checks <- data.frame(
  check = c(
    "pseudobulk_cell_count_conserved",
    "pseudobulk_UMI_count_conserved",
    "pseudobulk_gene_order_conserved",
    "pseudobulk_columns_match_metadata",
    "pseudobulk_counts_nonnegative_integer_like"
  ),
  passed = c(
    sum(pseudobulk_meta$n_cells) == ncol(epithelial_object),
    isTRUE(all.equal(
      as.numeric(sum(pseudobulk_counts)),
      as.numeric(sum(original_counts)),
      tolerance = 0
    )),
    identical(rownames(pseudobulk_counts), rownames(original_counts)),
    identical(colnames(pseudobulk_counts), rownames(pseudobulk_meta)),
    all(pseudobulk_counts@x >= 0) &
      all(abs(pseudobulk_counts@x - round(pseudobulk_counts@x)) < 1e-8)
  ),
  stringsAsFactors = FALSE
)
checks <- rbind(checks, conservation_checks)
write_tsv(checks, file.path(result_dir, "independent_validation_checks.tsv"))

cat("Independent Stage 6A audit\n")
print(attrition, row.names = FALSE)
cat("Checks:", sum(checks$passed), "/", nrow(checks), "passed\n")
if (!all(checks$passed)) quit(status = 1L)
