#!/usr/bin/env Rscript

# Register accepted Stage 7 headline values in the project result registry.

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
result_dir <- file.path(project_dir, "results", "07_singlecell_replication")
registry_path <- file.path(project_dir, "result_registry.tsv")

summary <- utils::read.delim(
  file.path(result_dir, "replication_summary.tsv"), check.names = FALSE
)
gate <- utils::read.delim(
  file.path(result_dir, "replication_gate.tsv"), check.names = FALSE
)
qc161 <- utils::read.delim(
  file.path(result_dir, "GSE161277_qc_summary.tsv"), check.names = FALSE
)
qc132 <- utils::read.delim(
  file.path(result_dir, "GSE132465_qc_summary.tsv"), check.names = FALSE
)
eligibility <- utils::read.delim(
  file.path(result_dir, "preflight", "GSE132465_pair_eligibility.tsv"),
  check.names = FALSE
)
validation <- utils::read.delim(
  file.path(result_dir, "validation_checks.tsv"), check.names = FALSE
)
acceptance <- utils::read.delim(
  file.path(result_dir, "acceptance_qc_checks.tsv"), check.names = FALSE
)
registry <- utils::read.delim(registry_path, check.names = FALSE)
registry <- registry[registry$stage != "7", , drop = FALSE]

row <- function(
  result_id, metric, value, unit, scope, interpretation, source_file,
  source_fields, generated_by
) {
  data.frame(
    result_id = result_id, stage = "7", metric = metric, value = value,
    unit = unit, scope = scope, interpretation = interpretation,
    source_file = source_file, source_fields = source_fields,
    generated_by = generated_by, git_ref = "stage-7",
    stringsAsFactors = FALSE
  )
}

headline <- rbind(
  row(
    "7_GSE161277_INPUT_CELLS", "input_cells", sum(qc161$cells_input),
    "cells", "GSE161277", "All deposited cells before per-capture QC",
    "results/07_singlecell_replication/GSE161277_qc_summary.tsv",
    "cells_input", "scripts/07_prepare_GSE161277.R"
  ),
  row(
    "7_GSE161277_RETAINED_CELLS", "retained_tissue_cells",
    sum(qc161$cells_retained), "cells", "GSE161277",
    "Tissue cells after adaptive QC and per-capture doublet detection",
    "results/07_singlecell_replication/GSE161277_qc_summary.tsv",
    "cells_retained", "scripts/07_prepare_GSE161277.R"
  ),
  row(
    "7_GSE161277_MATCHED_DONORS", "matched_donors", 3, "donors",
    "GSE161277 normal-adenoma-cancer", "Patient-level biological replicates",
    "results/07_singlecell_replication/replication_effects_all.tsv",
    "cohort;contrast;n_donors", "scripts/07_run_replication_models.R"
  ),
  row(
    "7_GSE132465_INPUT_CELLS", "input_cells", sum(qc132$cells_input),
    "cells", "GSE132465", "Deposited processed cells before residual QC",
    "results/07_singlecell_replication/GSE132465_qc_summary.tsv",
    "cells_input", "scripts/07_prepare_GSE132465.R"
  ),
  row(
    "7_GSE132465_RETAINED_CELLS", "retained_cells",
    sum(qc132$cells_retained), "cells", "GSE132465",
    "Cells retained after conservative residual QC",
    "results/07_singlecell_replication/GSE132465_qc_summary.tsv",
    "cells_retained", "scripts/07_prepare_GSE132465.R"
  ),
  row(
    "7_GSE132465_EVALUABLE_PAIRS", "evaluable_epithelial_pairs",
    length(unique(
      eligibility$donor_id[eligibility$paired_analysis_eligible]
    )),
    "donor pairs", "GSE132465 matched tumor-normal",
    "Eight of ten source-matched donors passed both-condition epithelial gates",
    "results/07_singlecell_replication/preflight/GSE132465_pair_eligibility.tsv",
    "donor_id;paired_analysis_eligible;exclusion_reason",
    "scripts/07_model_preflight.R"
  ),
  row(
    "7_REPLICATED_MODULES", "replicated_locked_modules",
    sum(summary$replication_class == "replicated"), "modules",
    "GSE161277 and GSE132465",
    "Secondary exploratory module-level replication; not primary discovery",
    "results/07_singlecell_replication/replication_summary.tsv",
    "module_id;replication_class", "scripts/07_run_replication_models.R"
  ),
  row(
    "7_UNDERPOWERED_MODULES", "directionally_consistent_underpowered_modules",
    sum(summary$replication_class ==
          "directionally consistent but underpowered"),
    "modules", "GSE161277 and GSE132465",
    "Directionally supportive without full replication criterion",
    "results/07_singlecell_replication/replication_summary.tsv",
    "module_id;replication_class", "scripts/07_run_replication_models.R"
  ),
  row(
    "7_REPLICATION_GATE", "replication_gate_passed",
    gate$passed[gate$criterion == "overall_replication_gate"],
    "boolean", "Stage 7", "Exploratory module-level gate only",
    "results/07_singlecell_replication/replication_gate.tsv",
    "criterion;passed;note", "scripts/07_run_replication_models.R"
  ),
  row(
    "7_SERVER_VALIDATION", "server_validation_checks_passed",
    sum(validation$passed), "checks", "Stage 7",
    "All technical server checks passed",
    "results/07_singlecell_replication/validation_checks.tsv",
    "check;passed", "scripts/07_validate_outputs.R"
  ),
  row(
    "7_ACCEPTANCE_VALIDATION", "independent_acceptance_checks_passed",
    sum(acceptance$passed), "checks", "Stage 7",
    "All donor-level effects and classifications independently recomputed",
    "results/07_singlecell_replication/acceptance_qc_checks.tsv",
    "check;passed", "scripts/07_acceptance_qc.R"
  )
)

module_rows <- do.call(rbind, lapply(seq_len(nrow(summary)), function(i) {
  x <- summary[i, ]
  scope <- paste(x$module_id, x$discovery_epithelial_state, sep = ";")
  source <- "results/07_singlecell_replication/replication_summary.tsv"
  rbind(
    row(
      paste0("7_", toupper(x$module_id), "_GSE161277_EARLY_EFFECT"),
      "GSE161277_adenoma_minus_normal_effect",
      x$GSE161277_normal_to_adenoma_effect, "standardized module score",
      scope, "Three matched donors; secondary exploratory replication",
      source,
      "module_id;GSE161277_normal_to_adenoma_effect",
      "scripts/07_run_replication_models.R"
    ),
    row(
      paste0("7_", toupper(x$module_id), "_GSE161277_EARLY_FDR"),
      "GSE161277_adenoma_minus_normal_FDR",
      x$GSE161277_normal_to_adenoma_fdr, "probability",
      scope, "BH FDR across six locked modules",
      source,
      "module_id;GSE161277_normal_to_adenoma_fdr",
      "scripts/07_run_replication_models.R"
    ),
    row(
      paste0("7_", toupper(x$module_id), "_GSE132465_CANCER_EFFECT"),
      "GSE132465_cancer_minus_normal_effect",
      x$GSE132465_normal_to_cancer_effect, "standardized module score",
      scope, "Eight evaluable matched donors",
      source,
      "module_id;GSE132465_normal_to_cancer_effect",
      "scripts/07_run_replication_models.R"
    ),
    row(
      paste0("7_", toupper(x$module_id), "_GSE132465_CANCER_FDR"),
      "GSE132465_cancer_minus_normal_FDR",
      x$GSE132465_normal_to_cancer_fdr, "probability",
      scope, "BH FDR across six locked modules",
      source,
      "module_id;GSE132465_normal_to_cancer_fdr",
      "scripts/07_run_replication_models.R"
    ),
    row(
      paste0("7_", toupper(x$module_id), "_SPECIFICITY_EFFECT"),
      "epithelial_specificity_effect",
      x$epithelial_specificity_effect, "standardized module score",
      scope, "Epithelial minus median non-epithelial donor-level effect",
      source,
      "module_id;epithelial_specificity_effect",
      "scripts/07_run_replication_models.R"
    ),
    row(
      paste0("7_", toupper(x$module_id), "_CLASS"),
      "replication_class", x$replication_class, "category",
      scope, "Prespecified locked-module classification",
      source,
      "module_id;replication_class",
      "scripts/07_run_replication_models.R"
    )
  )
}))

updated <- rbind(registry, headline, module_rows)
if (anyDuplicated(updated$result_id)) stop("Duplicate result_id after Stage 7 update")
utils::write.table(
  updated, registry_path, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
cat(
  "Registered ", nrow(headline) + nrow(module_rows),
  " accepted Stage 7 values\n", sep = ""
)

