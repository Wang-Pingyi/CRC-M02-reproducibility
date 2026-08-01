#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09C_independent_acceptance.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09C_external_test", run_id)
model_path <- file.path(project_dir, "objects", "locked_stool_model.rds")
test_object_path <- file.path(
  project_dir, "objects", paste0("GSE99573_9C_test_RMA_", run_id, ".rds")
)
expected_model_sha <-
  "3d4cd825b81fada4d0a3a92907d766dba4cdfea8e3ef7ce926b6790339c861fb"
expected_raw_sha <-
  "80cc17d2317d26da0dd739fe76f292b3c2212f5b6ced4c85611f94f3fe2dcdff"

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
result_tsv <- function(name) read_tsv(file.path(result_dir, name))
checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    observed = paste(observed, collapse = ";"),
    expected = paste(expected, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
auc_rank <- function(y, p) {
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

model <- readRDS(model_path)
test_object <- readRDS(test_object_path)
predictions <- result_tsv("stool_test_predictions_source_data.tsv")
results <- result_tsv("stool_test_results.tsv")
fixed <- result_tsv("fixed_specificity_results.tsv")
server_checks <- result_tsv("stage_9C_validation_checks.tsv")
provenance <- result_tsv("stage_9C_provenance.tsv")
baseline <- result_tsv("clinical_baseline_comparison.tsv")
access <- result_tsv("test_access_transition.tsv")
qc <- result_tsv("test_sample_qc_metrics.tsv")
manifest <- read_tsv(file.path(project_dir, "metadata", "dataset_manifest.tsv"))
test_ids <- sort(manifest$sample_id[
  manifest$accession == "GSE99573" &
    manifest$validation_split == "testing"
])
forbidden_ids <- manifest$sample_id[
  manifest$accession == "GSE99573" &
    manifest$validation_split != "testing"
]

add_check(
  "server_validation_passed",
  nrow(server_checks) == 21L && all(server_checks$passed),
  paste(sum(server_checks$passed), nrow(server_checks), sep = "/"),
  "21/21"
)
model_sha <- digest::digest(model_path, algo = "sha256", file = TRUE)
add_check(
  "locked_model_hash_anchored_after_test",
  identical(model_sha, expected_model_sha),
  model_sha,
  expected_model_sha
)
raw_path <- file.path(
  project_dir, "data_raw", "GSE99573", "GSE99573_RAW.tar"
)
raw_sha <- digest::digest(raw_path, algo = "sha256", file = TRUE)
add_check(
  "raw_archive_hash_unchanged",
  identical(raw_sha, expected_raw_sha),
  raw_sha,
  expected_raw_sha
)
add_check(
  "model_identity_and_lock_unchanged",
  identical(model$artifact, "CRC_carcinogenesis_locked_stool_model") &&
    identical(model$run_id, "20260730_001838") &&
    identical(model$test_expression_accessed, FALSE) &&
    bitwAnd(
      as.integer(file.info(model_path)$mode),
      as.integer(strtoi("222", base = 8))
    ) == 0L,
  c(model$artifact, model$run_id, model$test_expression_accessed),
  c(
    "CRC_carcinogenesis_locked_stool_model",
    "20260730_001838", FALSE
  )
)
add_check(
  "one_time_marker_and_completion_consistent",
  file.exists(file.path(
    project_dir, "logs", "09C_external_test", "ONE_TIME_TEST_STARTED"
  )) &&
    file.exists(file.path(
      project_dir, "logs", "09C_external_test", "READY_FOR_CODEX_QC"
    )) &&
    !file.exists(file.path(
      project_dir, "logs", "09C_external_test", "NEEDS_CODEX_ATTENTION"
    )),
  "started=yes;ready=yes;failure=no",
  "started=yes;ready=yes;failure=no"
)
add_check(
  "test_object_dimensions_and_finiteness",
  identical(dim(test_object$expression), c(70523L, 65L)) &&
    identical(dim(test_object$gene_expression), c(521L, 65L)) &&
    all(is.finite(test_object$expression)) &&
    all(is.finite(test_object$gene_expression)),
  c(dim(test_object$expression), dim(test_object$gene_expression)),
  c(70523L, 65L, 521L, 65L)
)
add_check(
  "test_participants_exact_and_forbidden_absent",
  setequal(test_object$sample_id, test_ids) &&
    !any(test_object$sample_id %in% forbidden_ids) &&
    !any(predictions$sample_id %in% forbidden_ids),
  c(
    length(unique(test_object$sample_id)),
    sum(test_object$sample_id %in% forbidden_ids),
    sum(predictions$sample_id %in% forbidden_ids)
  ),
  c(65L, 0L, 0L)
)
add_check(
  "all_test_samples_retained_without_label_QC",
  nrow(qc) == 65L &&
    all(is.finite(qc$mean_interarray_correlation)) &&
    !any(qc$automatic_exclusion) &&
    !any(qc$outcome_labels_used_for_QC),
  c(
    nrow(qc), sum(qc$automatic_exclusion),
    sum(qc$outcome_labels_used_for_QC)
  ),
  c(65L, 0L, 0L)
)

endpoint_order <- c(
  "adenoma_vs_normal", "cancer_vs_normal", "neoplasia_vs_normal"
)
max_probability_error <- 0
metric_audit <- list()
for (endpoint in endpoint_order) {
  z <- predictions[predictions$endpoint == endpoint, , drop = FALSE]
  endpoint_model <- model$models[[endpoint]]
  x <- t(test_object$gene_expression[, z$sample_id, drop = FALSE])
  beta <- endpoint_model$coefficients[colnames(x)]
  p <- stats::plogis(
    as.numeric(endpoint_model$intercept + x %*% beta)
  )
  max_probability_error <- max(
    max_probability_error,
    max(abs(p - z$predicted_probability))
  )
  y <- as.integer(z$outcome)
  predicted <- as.integer(p >= 0.5)
  tp <- sum(predicted == 1L & y == 1L)
  tn <- sum(predicted == 0L & y == 0L)
  fp <- sum(predicted == 1L & y == 0L)
  fn <- sum(predicted == 0L & y == 1L)
  expected <- c(
    AUC = auc_rank(y, p),
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    PPV = if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_,
    NPV = if ((tn + fn) > 0L) tn / (tn + fn) else NA_real_,
    accuracy = (tp + tn) / length(y),
    Brier = mean((p - y)^2)
  )
  observed <- vapply(
    names(expected),
    function(metric) {
      results$estimate[
        results$endpoint == endpoint & results$metric == metric
      ]
    },
    numeric(1)
  )
  difference <- abs(expected - observed)
  difference[is.na(expected) & is.na(observed)] <- 0
  metric_audit[[endpoint]] <- data.frame(
    endpoint = endpoint,
    metric = names(expected),
    independently_recalculated = as.numeric(expected),
    reported = as.numeric(observed),
    absolute_difference = as.numeric(difference),
    stringsAsFactors = FALSE
  )
}
metric_audit <- do.call(rbind, metric_audit)
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
write_tsv(
  metric_audit,
  file.path(
    project_dir, "results_final",
    "stage_9C_independent_metric_recalculation.tsv"
  )
)
add_check(
  "locked_probabilities_independently_reproduce",
  is.finite(max_probability_error) && max_probability_error < 1e-12,
  format(max_probability_error, scientific = TRUE),
  "<1e-12"
)
add_check(
  "key_metrics_independently_reproduce",
  all(is.finite(metric_audit$absolute_difference)) &&
    max(metric_audit$absolute_difference) < 1e-12,
  format(max(metric_audit$absolute_difference), scientific = TRUE),
  "<1e-12"
)
auc <- results[results$metric == "AUC", , drop = FALSE]
add_check(
  "AUC_inventory_complete_and_CI_ordered",
  nrow(auc) == 3L &&
    all(is.finite(auc$estimate)) &&
    all(is.finite(auc$ci_low)) &&
    all(is.finite(auc$ci_high)) &&
    all(auc$ci_low <= auc$estimate & auc$estimate <= auc$ci_high),
  sprintf(
    "%s=%.4f[%.4f,%.4f]",
    auc$endpoint, auc$estimate, auc$ci_low, auc$ci_high
  ),
  "three finite ordered AUC estimates and CIs"
)
add_check(
  "fixed_specificity_operating_points_not_test_optimized",
  nrow(fixed) == 3L &&
    all(fixed$threshold_source == "Stage_9B_training_outer_OOF") &&
    all(fixed$target_specificity == 0.9),
  c(
    nrow(fixed), paste(fixed$threshold_source, collapse = "/"),
    paste(fixed$target_specificity, collapse = "/")
  ),
  "3;training-only thresholds;0.9/0.9/0.9"
)
add_check(
  "provenance_prohibitions_preserved",
  identical(
    provenance$value[provenance$field == "model_retrained"], "FALSE"
  ) &&
    identical(
      provenance$value[provenance$field == "features_changed"], "FALSE"
    ) &&
    identical(
      provenance$value[provenance$field == "threshold_changed"], "FALSE"
    ) &&
    identical(
      provenance$value[provenance$field == "test_samples"], "65"
    ),
  paste(
    provenance$value[provenance$field %in% c(
      "model_retrained", "features_changed", "threshold_changed",
      "test_samples"
    )],
    collapse = "/"
  ),
  "FALSE/FALSE/FALSE/65"
)
add_check(
  "access_transition_once_for_all_65",
  nrow(access) == 65L &&
    setequal(access$sample_id, test_ids) &&
    all(access$stage_9C_prediction_generated) &&
    !any(access$stage_9C_excluded),
  c(
    nrow(access), sum(access$stage_9C_prediction_generated),
    sum(access$stage_9C_excluded)
  ),
  c(65L, 65L, 0L)
)
add_check(
  "clinical_baseline_limitation_honest",
  nrow(baseline) == 1L &&
    !baseline$clinical_baseline_available &&
    identical(baseline$status, "not_evaluable") &&
    !baseline$test_outcomes_used_to_build_baseline,
  c(
    baseline$clinical_baseline_available, baseline$status,
    baseline$test_outcomes_used_to_build_baseline
  ),
  c(FALSE, "not_evaluable", FALSE)
)
report_path <- file.path(
  project_dir, "reports", "stage_9C_external_test.md"
)
report_text <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
add_check(
  "report_discloses_exploratory_negative_result",
  grepl(
    "does not establish a clinically useful stool classifier",
    report_text, fixed = TRUE
  ) &&
    grepl("No endpoint, sample, feature", report_text, fixed = TRUE) &&
    grepl("must not be reused", report_text, fixed = TRUE),
  "required interpretation and stopping statements searched",
  "exploratory limitation, no post-test changes, no reuse"
)
figure_paths <- file.path(
  project_dir, "figures", "09C_external_test", run_id,
  c(
    "stool_test_ROC.pdf", "stool_test_ROC.png",
    "stool_test_calibration.pdf", "stool_test_calibration.png"
  )
)
source_paths <- file.path(
  result_dir,
  c(
    "stool_test_predictions_source_data.tsv",
    "ROC_source_data.tsv", "calibration_source_data.tsv",
    "fixed_specificity_results.tsv", "test_sample_qc_metrics.tsv"
  )
)
add_check(
  "source_data_and_figures_present",
  all(file.exists(figure_paths)) &&
    all(file.info(figure_paths)$size > 1000L) &&
    all(file.exists(source_paths)) &&
    all(file.info(source_paths)$size > 0L),
  c(sum(file.exists(figure_paths)), sum(file.exists(source_paths))),
  c(4L, 5L)
)

checks_df <- do.call(rbind, checks)
write_tsv(
  checks_df,
  file.path(
    project_dir, "results_final",
    "stage_9C_independent_acceptance_checks.tsv"
  )
)
if (!all(checks_df$passed)) {
  print(checks_df[!checks_df$passed, , drop = FALSE])
  stop("Stage 9C independent acceptance failed")
}
cat(
  "STAGE9C_INDEPENDENT_ACCEPTANCE_PASS ",
  sum(checks_df$passed), "/", nrow(checks_df), "\n",
  sep = ""
)
