#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08A_validate_outputs.R <project_dir> <run_id>")
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "08A_bulk_preprocessing", run_id)
matrix_dir <- file.path(project_dir, "data_processed", "stage_8A_bulk_preprocessing", run_id)
figure_dir <- file.path(project_dir, "figures", "08A_bulk_preprocessing", run_id)
cohorts <- c(GSE41657 = 88L, GSE100179 = 60L, GSE8671 = 64L)
checks <- list()
metrics <- list()

add_check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, pass = isTRUE(pass), detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

manifest <- utils::read.delim(file.path(project_dir, "metadata", "dataset_manifest.tsv"), check.names = FALSE)
membership <- utils::read.delim(
  file.path(project_dir, "results", "07_singlecell_replication", "preflight", "locked_module_membership.tsv"),
  check.names = FALSE
)
locked_modules <- unique(membership$module_id)

for (accession in names(cohorts)) {
  expected_n <- cohorts[[accession]]
  matrix_file <- file.path(matrix_dir, paste0(accession, "_normalized_probe_matrix.rds"))
  sample_file <- file.path(result_dir, accession, "sample_metadata.tsv")
  qc_file <- file.path(result_dir, accession, "qc_sample_metrics.tsv")
  feature_file <- file.path(result_dir, accession, "qc_feature_filter.tsv")
  map_file <- file.path(result_dir, accession, "locked_module_mapping_summary.tsv")
  unmapped_file <- file.path(result_dir, accession, "unmapped_locked_genes.tsv")
  required <- c(matrix_file, sample_file, qc_file, feature_file, map_file, unmapped_file)
  add_check(paste0(accession, "_required_outputs"), all(file.exists(required)), paste(basename(required[!file.exists(required)]), collapse = ","))
  if (!all(file.exists(required))) next
  mat <- readRDS(matrix_file)
  meta <- utils::read.delim(sample_file, check.names = FALSE)
  qc <- utils::read.delim(qc_file, check.names = FALSE)
  fmap <- utils::read.delim(map_file, check.names = FALSE)
  unmapped <- utils::read.delim(unmapped_file, check.names = FALSE)
  feature <- utils::read.delim(feature_file, check.names = FALSE)
  manifest_ids <- manifest$sample_id[manifest$accession == accession & manifest$inclusion == "include"]
  add_check(paste0(accession, "_sample_count"), ncol(mat) == expected_n && nrow(meta) == expected_n && nrow(qc) == expected_n, paste(ncol(mat), nrow(meta), nrow(qc), sep = "/"))
  add_check(paste0(accession, "_sample_identity"), identical(sort(colnames(mat)), sort(manifest_ids)), "matrix columns versus included manifest sample_id")
  add_check(paste0(accession, "_unique_samples"), !anyDuplicated(colnames(mat)) && !anyDuplicated(meta$sample_id), "no duplicate sample IDs")
  add_check(paste0(accession, "_locked_modules"), setequal(fmap$module_id, locked_modules), paste(length(unique(fmap$module_id)), "modules"))
  mapping_ok <- all(fmap$mapped >= 0L) &&
    all(fmap$unmapped_genes >= 0L) &&
    all(fmap$mapped <= fmap$locked_genes) &&
    all(fmap$mapped + fmap$unmapped_genes == fmap$locked_genes)
  add_check(paste0(accession, "_mapping_accounted"), mapping_ok, "unique mapped genes plus unique unmapped genes equals locked")
  unmapped_ok <- !any(unmapped$mapped) &&
    !anyDuplicated(unmapped[c("module_id", "epithelial_state", "gene")]) &&
    nrow(unmapped) == sum(fmap$unmapped_genes)
  add_check(paste0(accession, "_unmapped_audit"), unmapped_ok, paste(nrow(unmapped), "unique unmapped module-gene entries"))
  add_check(paste0(accession, "_qc_feature_audit"), feature$full_matrix_preserved[[1]] && feature$complete_features_used_for_pca[[1]] >= 10, paste(feature$complete_features_used_for_pca[[1]], "complete features"))
  add_check(paste0(accession, "_figures"), all(file.exists(file.path(figure_dir, paste0(accession, "_qc.", c("png", "pdf"))))), "PNG and PDF")
  finite <- mat[is.finite(mat)]
  q <- stats::quantile(finite, c(0, 0.01, 0.5, 0.99, 1), na.rm = TRUE)
  plausible_log2 <- q[[4]] < 30 && q[[1]] > -20
  add_check(paste0(accession, "_plausible_log2_scale"), plausible_log2, paste(names(q), signif(q, 6), collapse = "; "))
  metrics[[accession]] <- data.frame(
    accession = accession, features = nrow(mat), samples = ncol(mat),
    nonfinite_values = sum(!is.finite(mat)), qc_review_flags = sum(qc$qc_flag_review),
    min = q[[1]], p01 = q[[2]], median = q[[3]], p99 = q[[4]], max = q[[5]],
    stringsAsFactors = FALSE
  )
}

before <- file.path(result_dir, "..", paste0("raw_inputs_", run_id, ".before.sha256"))
after <- file.path(result_dir, "..", paste0("raw_inputs_", run_id, ".after.sha256"))
before <- normalizePath(before, mustWork = FALSE)
after <- normalizePath(after, mustWork = FALSE)
hash_ok <- file.exists(before) && file.exists(after) && identical(readLines(before), readLines(after))
add_check("raw_input_sha256_unchanged", hash_ok, "before and after manifests identical")

prov <- file.path(result_dir, "GSE8671", "input_processing_provenance.tsv")
prov_ok <- file.exists(prov)
prov_detail <- "missing provenance file"
if (prov_ok) {
  p <- utils::read.delim(prov, check.names = FALSE)
  required <- c("selected_input", "processing_status", "reason", "source_url", "source_sha256")
  shape_ok <- nrow(p) == 1L && all(required %in% names(p))
  if (shape_ok) {
    raw_path <- identical(p$selected_input[[1L]], "raw_CEL") &&
      identical(p$processing_status[[1L]], "raw_RMA_completed")
    fallback_path <- identical(p$selected_input[[1L]], "official_GEO_series_matrix") &&
      identical(p$processing_status[[1L]], "raw_RMA_not_run")
    source_ok <- !is.na(p$source_url[[1L]]) && nzchar(p$source_url[[1L]]) &&
      !is.na(p$source_sha256[[1L]]) && nchar(p$source_sha256[[1L]]) == 64L
    reason_ok <- !is.na(p$reason[[1L]]) && nzchar(p$reason[[1L]])
    prov_ok <- (raw_path || fallback_path) && source_ok && reason_ok
    prov_detail <- paste(p$selected_input[[1L]], p$processing_status[[1L]], "with URL/reason/SHA256")
  } else {
    prov_ok <- FALSE
    prov_detail <- "provenance columns or row count invalid"
  }
}
add_check("GSE8671_input_provenance", prov_ok, prov_detail)

report <- file.path(project_dir, "reports", "stage_8A_bulk_preprocessing.md")
report_text <- if (file.exists(report)) paste(readLines(report, warn = FALSE), collapse = "\n") else ""
boundary_ok <- grepl("No ordered-trend test", report_text, fixed = TRUE) &&
  grepl("do not enter Stage 8B", report_text, ignore.case = TRUE)
add_check("stage_boundary", boundary_ok, "preprocessing only; no Stage 8B")

check_df <- do.call(rbind, checks)
metric_df <- do.call(rbind, metrics)
utils::write.table(check_df, file.path(result_dir, "validation_checks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(metric_df, file.path(result_dir, "matrix_qc_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("STAGE8A_VALIDATION", sum(check_df$pass), "/", nrow(check_df), "\n")
print(metric_df)
if (!all(check_df$pass)) {
  print(check_df[!check_df$pass, , drop = FALSE])
  quit(status = 1L)
}
