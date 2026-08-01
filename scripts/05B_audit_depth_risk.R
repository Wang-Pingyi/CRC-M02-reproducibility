#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else getwd()
result_dir <- file.path(project_root, "results", "05B_full_qc_integration")
source_dir <- file.path(result_dir, "source_data")

metric_files <- sort(list.files(
  source_dir,
  pattern = "^GSM[0-9]+_cell_qc_metrics\\.tsv\\.gz$",
  full.names = TRUE
))
if (length(metric_files) != 72L) {
  stop("Expected 72 per-library QC metric files; found ", length(metric_files))
}

audit_one <- function(path) {
  x <- fread(
    path,
    select = c("sample_id", "nCount", "nFeature", "percent_mt", "exclusion_reason")
  )
  retained <- x[exclusion_reason == "retained"]
  if (nrow(retained) == 0L) {
    return(data.table(
      sample_id = unique(x$sample_id),
      cells_retained = 0L,
      nCount_p10 = NA_real_,
      nCount_median = NA_real_,
      nCount_p90 = NA_real_,
      nFeature_p10 = NA_real_,
      nFeature_median = NA_real_,
      nFeature_p90 = NA_real_,
      retained_nCount_lt_100 = 0L,
      retained_nCount_lt_200 = 0L,
      retained_nFeature_lt_100 = 0L,
      retained_nFeature_le_400 = 0L,
      retained_nFeature_ge_4000 = 0L,
      retained_nCount_ge_10000 = 0L,
      retained_percent_mt_ge_5 = 0L,
      retained_any_original_qc_violation = 0L
    ))
  }
  data.table(
    sample_id = unique(x$sample_id),
    cells_retained = nrow(retained),
    nCount_p10 = as.numeric(quantile(retained$nCount, 0.10, names = FALSE)),
    nCount_median = median(retained$nCount),
    nCount_p90 = as.numeric(quantile(retained$nCount, 0.90, names = FALSE)),
    nFeature_p10 = as.numeric(quantile(retained$nFeature, 0.10, names = FALSE)),
    nFeature_median = median(retained$nFeature),
    nFeature_p90 = as.numeric(quantile(retained$nFeature, 0.90, names = FALSE)),
    retained_nCount_lt_100 = sum(retained$nCount < 100),
    retained_nCount_lt_200 = sum(retained$nCount < 200),
    retained_nFeature_lt_100 = sum(retained$nFeature < 100),
    retained_nFeature_le_400 = sum(retained$nFeature <= 400),
    retained_nFeature_ge_4000 = sum(retained$nFeature >= 4000),
    retained_nCount_ge_10000 = sum(retained$nCount >= 10000),
    retained_percent_mt_ge_5 = sum(retained$percent_mt >= 5),
    retained_any_original_qc_violation = sum(
      retained$nFeature <= 400 |
        retained$nFeature >= 4000 |
        retained$nCount >= 10000 |
        retained$percent_mt >= 5
    )
  )
}

depth_audit <- rbindlist(lapply(metric_files, audit_one))
depth_audit[, original_qc_violation_percent :=
              100 * retained_any_original_qc_violation / cells_retained]
depth_audit[, depth_risk := fifelse(
  cells_retained == 0L,
  "excluded_library",
  fifelse(
  nCount_median < 100 | nFeature_median < 100,
  "critical_median_below_100",
  fifelse(
    nCount_median < 200 | nFeature_median < 200,
    "high_median_below_200",
    fifelse(
      retained_any_original_qc_violation > 0,
      "original_qc_violation_present",
      "none"
    )
  )
  )
)]
setorder(depth_audit, nFeature_median, nCount_median, sample_id)
fwrite(
  depth_audit,
  file.path(result_dir, "library_depth_risk_audit.tsv"),
  sep = "\t",
  quote = TRUE
)

summary <- data.table(
  metric = c(
    "cells_retained",
    "retained_nCount_lt_100",
    "retained_nCount_lt_200",
    "retained_nFeature_lt_100",
    "retained_nFeature_le_400",
    "retained_nFeature_ge_4000",
    "retained_nCount_ge_10000",
    "retained_percent_mt_ge_5",
    "retained_any_original_qc_violation",
    "libraries_critical_median_below_100",
    "libraries_high_median_below_200"
  ),
  value = c(
    sum(depth_audit$cells_retained),
    sum(depth_audit$retained_nCount_lt_100),
    sum(depth_audit$retained_nCount_lt_200),
    sum(depth_audit$retained_nFeature_lt_100),
    sum(depth_audit$retained_nFeature_le_400),
    sum(depth_audit$retained_nFeature_ge_4000),
    sum(depth_audit$retained_nCount_ge_10000),
    sum(depth_audit$retained_percent_mt_ge_5),
    sum(depth_audit$retained_any_original_qc_violation),
    sum(depth_audit$depth_risk == "critical_median_below_100"),
    sum(depth_audit$depth_risk == "high_median_below_200")
  )
)
fwrite(
  summary,
  file.path(result_dir, "original_qc_threshold_comparison.tsv"),
  sep = "\t",
  quote = TRUE
)

message(
  "Depth-risk audit completed: ",
  summary[metric == "retained_any_original_qc_violation", value],
  " retained cells violate at least one original-study QC threshold"
)
