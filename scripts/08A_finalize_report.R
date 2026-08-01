#!/usr/bin/env Rscript
set.seed(20260728)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08A_finalize_report.R <project_dir> <run_id>")
source(file.path(args[[1]], "scripts", "08A_helpers.R"))
project_dir <- normalizePath(args[[1]], mustWork = TRUE); run_id <- args[[2]]; paths <- stage8_paths(project_dir, run_id)
cohorts <- c("GSE41657", "GSE100179", "GSE8671")
qcs <- lapply(cohorts, function(a) utils::read.delim(file.path(paths$result, a, "qc_sample_metrics.tsv"), check.names = FALSE))
maps <- lapply(cohorts, function(a) utils::read.delim(file.path(paths$result, a, "locked_module_mapping_summary.tsv"), check.names = FALSE))
qc <- do.call(rbind, qcs); map <- do.call(rbind, maps)
stage8_write_tsv(qc, file.path(paths$result, "all_cohort_qc_sample_metrics.tsv"))
stage8_write_tsv(map, file.path(paths$result, "all_cohort_locked_module_mapping_summary.tsv"))
gse8671_provenance <- utils::read.delim(
  file.path(paths$result, "GSE8671", "input_processing_provenance.tsv"),
  check.names = FALSE
)
gse8671_method_line <- if (identical(gse8671_provenance$processing_status[[1]], "raw_RMA_completed")) {
  "- GSE8671: Affymetrix U133 Plus 2.0 raw CEL; `affy::rma` preprocessing with the official `hgu133plus2cdf` environment."
} else {
  "- GSE8671: raw CEL RMA could not be run; the explicitly audited fallback is the official GEO series processed matrix. See `input_processing_provenance.tsv`."
}
versions <- data.frame(package = names(sessionInfo()$otherPkgs), version = vapply(sessionInfo()$otherPkgs, function(x) as.character(x$Version), character(1)), stringsAsFactors = FALSE)
stage8_write_tsv(versions, file.path(paths$result, "software_versions.tsv"))
lines <- c(
  "# Stage 8A bulk tissue preprocessing (server-generated; pending Codex QC)", "",
  paste0("- Run ID: `", run_id, "`"),
  "- Scope: preprocessing and QC only. No ordered-trend test, differential expression, cross-platform merge, or candidate re-selection was performed.",
  "- GSE41258: not activated because it is outside the confirmed Stage 4B download inventory and the prespecified core tissue cohorts.",
  "- Raw archives: read-only; archives were extracted only to a run-specific `data_processed/` work directory.", "",
  "## Cohort methods", "",
  "- GSE41657: Agilent one-color raw TXT; `limma` normexp background correction (offset 50), log2 transformation and quantile normalization.",
  "- GSE100179: official GEO HTA 2.0 gene-level RMA-Sketch matrix; its labels were reconciled one-to-one against raw-CEL filenames before use.",
  gse8671_method_line, "",
  "## QC review flags", ""
 )
for (a in cohorts) {
  z <- qc[qc$accession == a, , drop = FALSE]
  lines <- c(lines, paste0("- ", a, ": ", nrow(z), " samples; ", sum(z$qc_flag_review), " robust-MAD QC review flag(s). Flags are not automatic exclusions."))
 }
lines <- c(lines, "", "## Locked exploratory module mapping", "")
for (a in cohorts) {
  z <- map[map$accession == a, , drop = FALSE]
  lines <- c(lines, paste0("- ", a, ": ", sum(z$mapped), "/", sum(z$locked_genes), " locked module-gene entries mapped; unmapped entries are saved explicitly in the cohort result directory."))
}
lines <- c(lines, "", "## Outputs", "", paste0("- Run-specific normalized matrices and extracted working copies: `", paths$root, "`"), paste0("- QC, probe annotation, duplicate-probe preservation, mapping audit, source data and software versions: `", paths$result, "`"), paste0("- QC figures: `", paths$figure, "`"), "", "## Next action", "", "Server processing is complete and requires Codex QC/acceptance. Do not enter Stage 8B until the investigator approves Stage 8A.")
writeLines(lines, paths$report)
cat("STAGE_8A_REPORT_WRITTEN\n")
