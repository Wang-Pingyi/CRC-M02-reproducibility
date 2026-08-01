#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09B_finalize_report.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09B_model_training", run_id)
report_path <- file.path(project_dir, "reports", "stage_9B_model_training.md")

read_tsv <- function(name) {
  utils::read.delim(
    file.path(result_dir, name),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
metrics <- read_tsv("nested_cv_performance.tsv")
summary <- read_tsv("model_lock_summary.tsv")
checks <- read_tsv("stage_9B_validation_checks.tsv")
lock_manifest <- read_tsv("model_lock_manifest.tsv")
coefficients <- read_tsv("locked_model_coefficients.tsv")
estimability <- read_tsv("metric_estimability_audit.tsv")
sha <- lock_manifest$value[lock_manifest$field == "model_sha256"]

metric_value <- function(endpoint, metric, field = "estimate") {
  z <- metrics[metrics$endpoint == endpoint & metrics$metric == metric, ]
  z[[field]][[1L]]
}
format_ci <- function(endpoint, metric, digits = 3L) {
  estimate <- metric_value(endpoint, metric, "estimate")
  low <- metric_value(endpoint, metric, "ci_low")
  high <- metric_value(endpoint, metric, "ci_high")
  if (!all(is.finite(c(estimate, low, high)))) {
    reason <- estimability$nonestimable_reason[
      estimability$endpoint == endpoint & estimability$metric == metric
    ]
    return(paste0("Not estimable (", reason[[1L]], ")"))
  }
  paste0(
    formatC(estimate, format = "f", digits = digits),
    " (95% CI ",
    formatC(low, format = "f", digits = digits),
    " to ",
    formatC(high, format = "f", digits = digits),
    ")"
  )
}

endpoint_labels <- c(
  adenoma_vs_normal = "Adenoma versus normal (primary)",
  cancer_vs_normal = "Colorectal cancer versus normal (secondary)",
  neoplasia_vs_normal =
    "Adenoma or colorectal cancer versus normal (secondary composite)"
)
selected <- coefficients[
  coefficients$selected & coefficients$feature != "(Intercept)",
  ,
  drop = FALSE
]
selected_description <- if (nrow(selected) == 0L) {
  "None"
} else {
  paste0(
    selected$endpoint,
    ": ",
    selected$feature,
    " (coefficient ",
    formatC(selected$coefficient, format = "f", digits = 3),
    ")",
    collapse = "; "
  )
}
AUCs <- vapply(
  names(endpoint_labels),
  function(endpoint) metric_value(endpoint, "AUC"),
  numeric(1)
)
AUC_lows <- vapply(
  names(endpoint_labels),
  function(endpoint) metric_value(endpoint, "AUC", "ci_low"),
  numeric(1)
)
AUC_highs <- vapply(
  names(endpoint_labels),
  function(endpoint) metric_value(endpoint, "AUC", "ci_high"),
  numeric(1)
)
lines <- c(
  "# Stage 9B training and locked stool models",
  "",
  paste0("- Run ID: `", run_id, "`"),
  "- Status: server processing and validation complete; independent acceptance is recorded separately in `reports/stage_9B_acceptance_audit.md`.",
  "- Training data: GSE99573 GEO training partition only.",
  "- Independent test data: not extracted, read, normalized, summarized, or modeled.",
  "- Algorithm: one prespecified LASSO-logistic pipeline for all endpoints; no algorithm search.",
  "- Selection: 5-fold outer / 5-fold inner nested CV; final `lambda.1se` from deterministic 10-fold training CV.",
  "- Classification threshold: 0.5, fixed before fitting.",
  "- Feature universe: 521 previously locked and detectable genes; no univariable screening.",
  "- QC handling: all 265 training arrays retained; nine review flags and zero automatic exclusions.",
  "",
  "## Locked endpoint models",
  "",
  "| Endpoint | N | Cases | Controls | Selected genes | Lambda | Threshold |",
  "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"
)
for (endpoint in names(endpoint_labels)) {
  z <- summary[summary$endpoint == endpoint, , drop = FALSE]
  lines <- c(
    lines,
    paste0(
      "| ", endpoint_labels[[endpoint]], " | ", z$training_n, " | ",
      z$training_cases, " | ", z$training_controls, " | ",
      z$selected_features, " | ",
      formatC(z$lambda, format = "e", digits = 3), " | ",
      formatC(z$threshold, format = "f", digits = 1), " |"
    )
  )
}
lines <- c(
  lines,
  "",
  "## Nested-CV training performance",
  "",
  "These values are internal out-of-fold training estimates and are not independent validation.",
  "",
  "| Endpoint | AUC | Sensitivity | Specificity | Brier score | Calibration intercept | Calibration slope |",
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
      format_ci(endpoint, "Brier"), " | ",
      format_ci(endpoint, "calibration_intercept"), " | ",
      format_ci(endpoint, "calibration_slope"), " |"
    )
  )
}
lines <- c(
  lines,
  "",
  "## Prespecified-threshold behavior",
  "",
  "| Endpoint | PPV | NPV |",
  "| --- | --- | --- |"
)
for (endpoint in names(endpoint_labels)) {
  lines <- c(
    lines,
    paste0(
      "| ", endpoint_labels[[endpoint]], " | ",
      format_ci(endpoint, "PPV"), " | ",
      format_ci(endpoint, "NPV"), " |"
    )
  )
}
lines <- c(
  lines,
  "",
  "At the frozen threshold of 0.5, the primary model predicted all samples as normal, while the cancer and composite models predicted all samples as cases. The corresponding PPV or NPV is therefore mathematically not estimable and is explicitly retained as `NA`; it is not imputed.",
  "",
  "## Scientific result",
  "",
  "- The prespecified `lambda.1se` rule selected zero genes for adenoma versus normal and cancer versus normal; these are valid intercept-only locked models.",
  paste0("- The composite model selected one gene: ", selected_description, "."),
  paste0(
    "- All three nested-CV AUC estimates were near chance (range ",
    formatC(min(AUCs), format = "f", digits = 3),
    " to ",
    formatC(max(AUCs), format = "f", digits = 3),
    "), and every 95% CI included 0.5 (",
    ifelse(all(AUC_lows <= 0.5 & AUC_highs >= 0.5), "confirmed", "not confirmed"),
    ")."
  ),
  "- The training data therefore do not support useful discrimination by these frozen candidate features under the prespecified pipeline. This is a negative scientific result, not a reason to change the algorithm, lambda rule, threshold or feature set.",
  "",
  "## Model lock",
  "",
  "- Canonical artifact: `objects/locked_stool_model.rds`.",
  paste0("- SHA256: `", sha, "`."),
  "- The model file is read-only and contains all endpoint formulas, original-scale coefficients, mappings, lambdas, thresholds, fold rules, random seed and software versions.",
  "- The desktop workflow creates the `stool-model-locked` Git tag only after independent Codex QC commits this checksum and report.",
  "",
  "## Validation and interpretation",
  "",
  paste0("- Server validation: ", sum(checks$passed), "/", nrow(checks), " checks passed."),
  "- Penalized coefficients are predictive parameters, not causal or inferential gene effects.",
  "- The six source modules remain exploratory; training performance cannot overturn the frozen negative primary Stage 6A result.",
  "- No feature or model was changed to improve AUC after viewing results.",
  "",
  "## Stage boundary",
  "",
  "Stop after Stage 9B. Do not run, normalize, inspect, or evaluate the 65-sample independent test set until the investigator separately authorizes the locked one-time evaluation."
)
writeLines(lines, report_path, useBytes = TRUE)
cat("STAGE9B_REPORT_WRITTEN\n")
