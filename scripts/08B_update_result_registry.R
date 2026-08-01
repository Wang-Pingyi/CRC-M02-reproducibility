#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08B_update_result_registry.R PROJECT_DIR RUN_ID")
project <- normalizePath(args[1], mustWork = TRUE)
run_id <- args[2]
registry_file <- file.path(project, "result_registry.tsv")
result_dir <- file.path(project, "results", "08B_bulk_validation", run_id)
registry <- read.delim(registry_file, check.names = FALSE)
meta <- read.delim(file.path(result_dir, "meta_analysis_results.tsv"), check.names = FALSE)
gate <- read.delim(file.path(result_dir, "tissue_validation_gate.tsv"), check.names = FALSE)
scores <- read.delim(file.path(result_dir, "bulk_module_scores.tsv"), check.names = FALSE)

registry <- registry[registry$stage != "8B", , drop = FALSE]
rows <- list()
add <- function(id, metric, value, unit, scope, interpretation, source, fields) {
  rows[[length(rows) + 1L]] <<- data.frame(
    result_id = id, stage = "8B", metric = metric, value = as.character(value),
    unit = unit, scope = scope, interpretation = interpretation,
    source_file = source, source_fields = fields,
    generated_by = "scripts/08B_meta_gate_report.R;scripts/08B_independent_acceptance.R",
    git_ref = "stage-8B", stringsAsFactors = FALSE
  )
}

add("8B_GSE8671_VERIFIED_PAIRS", "verified_donor_pairs", 32, "pairs", "GSE8671",
    "Explicit GEO patient #1-#32; independent reconstruction passed",
    "metadata/GSE8671_verified_pairs.tsv", "verified_pair_id;sample_id;condition")
add("8B_TISSUE_GATE", "tissue_validation_gate", unique(gate$stage_gate), "category",
    "GSE41657;GSE100179;GSE8671",
    "Prespecified exploratory-module gate; primary Stage 6A result remains negative",
    paste0("results/08B_bulk_validation/", run_id, "/tissue_validation_gate.tsv"),
    "stage_gate")

early <- meta[meta$analysis_set == "primary_all_samples" &
                meta$endpoint == "adenoma_vs_normal", ]
for (i in seq_len(nrow(early))) {
  x <- early[i, ]
  prefix <- paste0("8B_", toupper(x$module_id), "_EARLY_META_")
  scope <- paste0(x$module_id, "; three independent tissue cohorts")
  source <- paste0("results/08B_bulk_validation/", run_id, "/meta_analysis_results.tsv")
  add(paste0(prefix, "EFFECT"), "normal_to_adenoma_random_effect", x$pooled_effect,
      "cohort-standardized module score", scope, "REML Knapp-Hartung",
      source, "module_id;endpoint;analysis_set;pooled_effect")
  add(paste0(prefix, "CI_LOW"), "normal_to_adenoma_95CI_low", x$ci_low,
      "cohort-standardized module score", scope, "REML Knapp-Hartung",
      source, "module_id;endpoint;analysis_set;ci_low")
  add(paste0(prefix, "CI_HIGH"), "normal_to_adenoma_95CI_high", x$ci_high,
      "cohort-standardized module score", scope, "REML Knapp-Hartung",
      source, "module_id;endpoint;analysis_set;ci_high")
  add(paste0(prefix, "FDR"), "normal_to_adenoma_meta_FDR", x$fdr,
      "probability", scope, "BH across six locked modules",
      source, "module_id;endpoint;analysis_set;fdr")
  add(paste0(prefix, "I2"), "normal_to_adenoma_I2", x$I2,
      "percent", scope, "Between-cohort heterogeneity",
      source, "module_id;endpoint;analysis_set;I2")
  classification <- gate$classification[match(x$module_id, gate$module_id)]
  add(paste0(prefix, "CLASS"), "tissue_replication_class", classification,
      "category", scope, "Prespecified gate classification",
      paste0("results/08B_bulk_validation/", run_id, "/tissue_validation_gate.tsv"),
      "module_id;classification")
}

new_rows <- do.call(rbind, rows)
registry <- rbind(registry, new_rows[, names(registry), drop = FALSE])
write.table(registry, registry_file, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "NA")
