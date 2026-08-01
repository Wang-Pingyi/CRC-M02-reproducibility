#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09C_finalize_report.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09C_external_test", run_id)
report_path <- file.path(project_dir, "reports", "stage_9C_external_test.md")

read_tsv <- function(name) {
  utils::read.delim(
    file.path(result_dir, name),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
results <- read_tsv("stool_test_results.tsv")
fixed <- read_tsv("fixed_specificity_results.tsv")
checks <- read_tsv("stage_9C_validation_checks.tsv")
baseline <- read_tsv("clinical_baseline_comparison.tsv")
provenance <- read_tsv("stage_9C_provenance.tsv")
model_sha <- provenance$value[
  provenance$field == "model_sha256_after"
]

endpoint_labels <- c(
  adenoma_vs_normal = "Adenoma versus normal (primary)",
  cancer_vs_normal = "Colorectal cancer versus normal (secondary)",
  neoplasia_vs_normal =
    "Adenoma or colorectal cancer versus normal (secondary composite)"
)
metric_row <- function(endpoint, metric) {
  results[
    results$endpoint == endpoint & results$metric == metric,
    ,
    drop = FALSE
  ]
}
format_ci <- function(endpoint, metric, digits = 3L) {
  z <- metric_row(endpoint, metric)
  if (!nrow(z) || !all(is.finite(c(z$estimate, z$ci_low, z$ci_high)))) {
    return("Not estimable")
  }
  paste0(
    formatC(z$estimate, format = "f", digits = digits),
    " (95% CI ",
    formatC(z$ci_low, format = "f", digits = digits),
    " to ",
    formatC(z$ci_high, format = "f", digits = digits),
    ")"
  )
}
auc_rows <- results[results$metric == "AUC", , drop = FALSE]
weak_result <- all(
  auc_rows$ci_low <= 0.5 |
    auc_rows$estimate < 0.70
)
interpretation <- if (weak_result) {
  paste(
    "The one-time test result is an exploratory translational result;",
    "it does not establish a clinically useful stool classifier."
  )
} else {
  paste(
    "The one-time test result shows potentially useful discrimination,",
    "but remains an exploratory translational result requiring independent",
    "prospective confirmation."
  )
}

lines <- c(
  "# Stage 9C one-time independent stool test",
  "",
  paste0("- Run ID: `", run_id, "`"),
  "- Dataset: GSE99573 immutable 65-participant testing split.",
  "- Test inventory: 22 normal, 20 adenoma and 23 colorectal cancer.",
  "- Analysis unit: participant.",
  "- Preprocessing: independent 65-array core RMA batch; no outcome-informed preprocessing or exclusions.",
  "- Model application: stored original-scale coefficients, unchanged feature mapping and locked threshold 0.5.",
  paste0("- Locked model SHA256 before and after testing: `", model_sha, "`."),
  "- Model retrained: no.",
  "- Test evaluation count: one.",
  "",
  "## Independent-test performance at the locked threshold",
  "",
  "| Endpoint | AUC | Sensitivity | Specificity | PPV | NPV | Brier score |",
  "| --- | --- | --- | --- | --- | --- | --- |"
)
for (endpoint in names(endpoint_labels)) {
  lines <- c(
    lines,
    paste0(
      "| ", endpoint_labels[[endpoint]], " | ",
      format_ci(endpoint, "AUC"), " | ",
      format_ci(endpoint, "sensitivity"), " | ",
      format_ci(endpoint, "specificity"), " | ",
      format_ci(endpoint, "PPV"), " | ",
      format_ci(endpoint, "NPV"), " | ",
      format_ci(endpoint, "Brier"), " |"
    )
  )
}
lines <- c(
  lines,
  "",
  "AUC confidence intervals use DeLong's method except that a completely constant score has an exact conditional AUC and CI of 0.5. Classification-proportion confidence intervals use the Wilson score method. Undefined predictive values remain not estimable.",
  "",
  "## Training-locked 90% specificity operating point",
  "",
  "| Endpoint | Frozen threshold | Training specificity | Test specificity | Test sensitivity |",
  "| --- | ---: | ---: | --- | --- |"
)
for (endpoint in names(endpoint_labels)) {
  z <- fixed[fixed$endpoint == endpoint, , drop = FALSE]
  threshold_text <- if (is.infinite(z$training_locked_threshold)) {
    "Inf"
  } else {
    formatC(z$training_locked_threshold, format = "f", digits = 6)
  }
  lines <- c(
    lines,
    paste0(
      "| ", endpoint_labels[[endpoint]], " | ", threshold_text, " | ",
      formatC(z$training_specificity, format = "f", digits = 3), " | ",
      formatC(z$test_specificity, format = "f", digits = 3),
      " (95% CI ",
      formatC(z$test_specificity_ci_low, format = "f", digits = 3),
      " to ",
      formatC(z$test_specificity_ci_high, format = "f", digits = 3),
      ") | ",
      formatC(z$test_sensitivity, format = "f", digits = 3),
      " (95% CI ",
      formatC(z$test_sensitivity_ci_low, format = "f", digits = 3),
      " to ",
      formatC(z$test_sensitivity_ci_high, format = "f", digits = 3),
      ") |"
    )
  )
}
lines <- c(
  lines,
  "",
  "The high-specificity thresholds were derived from Stage 9B training-only out-of-fold predictions before test access and were not recalibrated in the test set.",
  "",
  "## Calibration",
  "",
  "| Endpoint | Calibration intercept | Calibration slope |",
  "| --- | --- | --- |"
)
for (endpoint in names(endpoint_labels)) {
  lines <- c(
    lines,
    paste0(
      "| ", endpoint_labels[[endpoint]], " | ",
      format_ci(endpoint, "calibration_intercept"), " | ",
      format_ci(endpoint, "calibration_slope"), " |"
    )
  )
}
lines <- c(
  lines,
  "",
  "Calibration estimates are computational performance summaries. Constant or nearly constant predictions can make calibration slopes unstable or not estimable.",
  "",
  "## Clinical baseline comparison",
  "",
  paste0("- Status: `", baseline$status, "`."),
  paste0("- Reason: ", baseline$reason),
  "- Disease category and cancer stage were not used as predictors because they define or reveal the endpoint.",
  "- The published SVM was not reconstructed because a frozen model and participant-level predictions were not available to this project.",
  "",
  "## Interpretation",
  "",
  interpretation,
  "All negative, weak or undefined results are retained. No endpoint, sample, feature, coefficient, threshold or algorithm was changed after viewing the test results.",
  "",
  "## Reproducibility and validation",
  "",
  paste0("- Server validation: ", sum(checks$passed), "/", nrow(checks), " checks passed."),
  "- Source tables: `stool_test_predictions_source_data.tsv`, `ROC_source_data.tsv`, `calibration_source_data.tsv`, `fixed_specificity_results.tsv` and `test_sample_qc_metrics.tsv`.",
  "- Figures: test ROC and calibration in PDF and 300-dpi PNG.",
  "",
  "## Stage boundary",
  "",
  "Stage 9C stops here. The independent test set must not be reused for model selection, threshold adjustment, feature replacement or retraining."
)
writeLines(lines, report_path, useBytes = TRUE)

manifest_lines <- c(
  "# Stage 9C analysis outputs",
  "",
  paste0("Generated by run `", run_id, "`."),
  "",
  "## Tables",
  "",
  "- `stool_test_results.tsv` — complete endpoint performance.",
  "- `fixed_specificity_results.tsv` — training-locked high-specificity operating points.",
  "- `clinical_baseline_comparison.tsv` — clinical baseline availability audit.",
  "- `stage_9C_validation_checks.tsv` — server validation.",
  "",
  "## Source data",
  "",
  "- `stool_test_predictions_source_data.tsv` — participant-level locked predictions.",
  "- `ROC_source_data.tsv` — ROC curve coordinates.",
  "- `calibration_source_data.tsv` — calibration plot bins.",
  "- `test_sample_qc_metrics.tsv` — outcome-blind array diagnostics.",
  "- `test_feature_reconstruction_audit.tsv` — frozen feature mapping audit.",
  "",
  "## Figures",
  "",
  "- `stool_test_ROC.pdf` / `stool_test_ROC.png`.",
  "- `stool_test_calibration.pdf` / `stool_test_calibration.png`."
)
writeLines(
  manifest_lines,
  file.path(result_dir, "_analysis_outputs.md"),
  useBytes = TRUE
)
cat("STAGE9C_REPORT_WRITTEN\n")
