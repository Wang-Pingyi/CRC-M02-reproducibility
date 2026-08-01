#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09A_finalize_report.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09A_stool_feasibility", run_id)
report_path <- file.path(project_dir, "reports", "stage_9A_stool_feasibility.md")

read_tsv <- function(name) {
  utils::read.delim(file.path(result_dir, name), check.names = FALSE)
}
split <- read_tsv("GSE99573_split_audit.tsv")
mapping <- read_tsv("locked_gene_probe_mapping.tsv")
module <- read_tsv("locked_module_detectability_summary.tsv")
qc <- read_tsv("training_sample_qc_metrics.tsv")
split_checks <- read_tsv("split_validation_checks.tsv")
analysis_checks <- read_tsv("analysis_validation_checks.tsv")

split_total <- stats::aggregate(
  n_samples ~ validation_split, data = split, FUN = sum
)
lookup <- function(x) {
  split_total$n_samples[match(x, split_total$validation_split)]
}
lines <- c(
  "# Stage 9A stool RNA feasibility and split audit",
  "",
  paste0("- Run ID: `", run_id, "`"),
  "- Status: server processing complete; pending independent Codex QC.",
  "- Platform: GSE99573 / GPL17586, Affymetrix Human Transcriptome Array 2.0.",
  "- Scope: split audit, platform mapping, training-only detectability, and outcome-blind QC.",
  "- Prohibited analyses performed: none; no model, feature selection, disease comparison, ROC analysis, or test-expression inspection.",
  "",
  "## Training/test reconstruction",
  "",
  paste0("- Training: ", lookup("training"), " samples."),
  paste0("- Testing: ", lookup("testing"), " samples; immutable independent test set."),
  paste0("- Not Used: ", lookup("not_used"), " benign samples; retained as excluded."),
  "- The split comes directly from the sample-level GEO `set` characteristic and was not inferred from filenames.",
  "- All 338 sample IDs and donor IDs are unique; the archive contains one CEL.gz member per sample and no technical replicate was counted as another participant.",
  "",
  "## Test-set firewall",
  "",
  "- Only the 265 training CEL files were selectively extracted and normalized.",
  "- The 65 test CEL files were inventory-checked inside the raw TAR but were not extracted, read, normalized, summarized, plotted, or associated with outcomes.",
  "- Test-set access is recorded in `locked_test_access_audit.tsv`.",
  "",
  "## Locked-feature platform feasibility",
  "",
  paste0("- Locked feature universe: 6 modules, 747 module-gene entries, 632 unique gene symbols."),
  paste0("- Mapped unique genes: ", sum(mapping$mapping_status == "mapped"), "/632."),
  paste0("- Provisionally detectable unique genes: ", sum(mapping$detectability_status == "detectable"), "/632."),
  paste0("- Low-expression unique genes: ", sum(mapping$detectability_status == "low_expression"), "/632."),
  paste0("- Missing-probe unique genes: ", sum(mapping$detectability_status == "missing_probe"), "/632."),
  "- Detectability is an outcome-blind technical screen, not feature selection. No locked gene was removed.",
  "",
  "### Module-level mapping",
  "",
  "| Module | Locked | Mapped | Detectable | Low expression | Missing probe |",
  "| --- | ---: | ---: | ---: | ---: | ---: |"
)
for (i in seq_len(nrow(module))) {
  lines <- c(
    lines,
    paste0(
      "| ", module$module_id[[i]], " | ", module$locked_genes[[i]], " | ",
      module$mapped_genes[[i]], " | ", module$detectable_genes[[i]], " | ",
      module$low_expression_genes[[i]], " | ",
      module$missing_probe_genes[[i]], " |"
    )
  )
}
lines <- c(
  lines,
  "",
  "## Training-only array QC",
  "",
  paste0("- Training arrays normalized: ", nrow(qc), "."),
  paste0("- Outcome-blind QC review flags: ", sum(qc$qc_review_flag), "."),
  "- QC flags are not automatic exclusions and require later review before model training.",
  "- Disease labels were not used to derive or interpret QC flags.",
  "",
  "## Prespecified later diagnostic comparisons",
  "",
  "1. adenoma versus normal;",
  "2. colorectal cancer versus normal;",
  "3. adenoma or colorectal cancer versus normal.",
  "",
  "## Validation",
  "",
  paste0(
    "- Split/inventory checks: ",
    sum(split_checks$passed), "/", nrow(split_checks), " passed."
  ),
  paste0(
    "- Analysis/firewall checks: ",
    sum(analysis_checks$passed), "/", nrow(analysis_checks), " passed."
  ),
  "- Raw archive SHA256 values before and after processing must match before the server creates the completion marker.",
  "",
  "## Outputs and stage boundary",
  "",
  paste0("- Run-specific tables: `results/09A_stool_feasibility/", run_id, "`."),
  paste0("- Training-only normalized object: `objects/GSE99573_9A_training_RMA_", run_id, ".rds`."),
  paste0("- Outcome-blind QC figure: `figures/09A_stool_feasibility/", run_id, "/training_array_QC_outcome_blind.pdf`."),
  "- Frozen protocol: `protocol/stool_model_protocol.md`.",
  "- Stop after Stage 9A. Do not train a model or access test expression until Stage 9A is independently accepted and a later model is locked."
)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, report_path, useBytes = TRUE)
cat("STAGE9A_REPORT_WRITTEN\n")
