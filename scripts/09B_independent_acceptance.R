#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09B_independent_acceptance.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09B_model_training", run_id)
model_path <- file.path(project_dir, "objects", "locked_stool_model.rds")
output_path <- file.path(
  project_dir,
  "results_final",
  "stage_9B_independent_acceptance_checks.tsv"
)
expected_sha256 <-
  "3d4cd825b81fada4d0a3a92907d766dba4cdfea8e3ef7ce926b6790339c861fb"

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

model <- readRDS(model_path)
manifest <- read_tsv(file.path(project_dir, "metadata", "dataset_manifest.tsv"))
manifest <- manifest[manifest$accession == "GSE99573", , drop = FALSE]
training_ids <- sort(manifest$sample_id[manifest$validation_split == "training"])
test_ids <- sort(manifest$sample_id[manifest$validation_split == "testing"])
not_used_ids <- sort(
  manifest$sample_id[manifest$validation_split == "not_used"]
)
predictions <- result_tsv("nested_cv_predictions_training_only.tsv")
performance <- result_tsv("nested_cv_performance.tsv")
coefficients <- result_tsv("locked_model_coefficients.tsv")
summary <- result_tsv("model_lock_summary.tsv")
folds <- result_tsv("fold_assignments.tsv")
server_checks <- result_tsv("stage_9B_validation_checks.tsv")
estimability <- result_tsv("metric_estimability_audit.tsv")
lock_manifest <- result_tsv("model_lock_manifest.tsv")
model_sha <- digest::digest(model_path, algo = "sha256", file = TRUE)

add_check(
  "server_validation_all_passed",
  nrow(server_checks) == 26L && all(server_checks$passed),
  paste(sum(server_checks$passed), nrow(server_checks), sep = "/"),
  "26/26"
)
add_check(
  "model_sha256_anchored",
  identical(model_sha, expected_sha256),
  model_sha,
  expected_sha256
)
add_check(
  "model_sha256_matches_lock_manifest",
  identical(
    model_sha,
    lock_manifest$value[lock_manifest$field == "model_sha256"]
  ),
  model_sha,
  lock_manifest$value[lock_manifest$field == "model_sha256"]
)
mode <- file.info(model_path)$mode
add_check(
  "model_is_read_only",
  bitwAnd(as.integer(mode), as.integer(strtoi("222", base = 8))) == 0L,
  mode,
  "no write bits"
)
add_check(
  "model_run_and_identity",
  identical(model$run_id, run_id) &&
    identical(model$artifact, "CRC_carcinogenesis_locked_stool_model") &&
    identical(model$artifact_version, "1.0.0"),
  c(model$run_id, model$artifact, model$artifact_version),
  c(run_id, "CRC_carcinogenesis_locked_stool_model", "1.0.0")
)
add_check(
  "test_firewall_in_model",
  identical(model$test_expression_accessed, FALSE) &&
    identical(model$training_split, "GEO set=Training"),
  c(model$test_expression_accessed, model$training_split),
  c(FALSE, "GEO set=Training")
)
add_check(
  "manifest_split_reconstructed",
  length(training_ids) == 265L &&
    length(test_ids) == 65L &&
    length(not_used_ids) == 8L,
  c(length(training_ids), length(test_ids), length(not_used_ids)),
  c(265L, 65L, 8L)
)
analysis_ids <- unique(c(predictions$sample_id, folds$sample_id))
add_check(
  "analysis_outputs_training_only",
  all(analysis_ids %in% training_ids) &&
    !any(analysis_ids %in% c(test_ids, not_used_ids)),
  c(
    sum(analysis_ids %in% training_ids),
    sum(analysis_ids %in% test_ids),
    sum(analysis_ids %in% not_used_ids)
  ),
  c(length(analysis_ids), 0L, 0L)
)
expected_endpoints <- c(
  "adenoma_vs_normal", "cancer_vs_normal", "neoplasia_vs_normal"
)
add_check(
  "endpoint_hierarchy_preserved",
  identical(names(model$models), expected_endpoints) &&
    identical(
      model$endpoint_hierarchy$endpoint[
        model$endpoint_hierarchy$role == "primary"
      ],
      "adenoma_vs_normal"
    ),
  c(
    names(model$models),
    model$endpoint_hierarchy$endpoint[
      model$endpoint_hierarchy$role == "primary"
    ]
  ),
  c(expected_endpoints, "adenoma_vs_normal")
)
add_check(
  "frozen_algorithm_preserved",
  identical(model$random_seed, 42L) &&
    identical(model$algorithm$family, "binomial") &&
    identical(model$algorithm$alpha, 1) &&
    identical(model$algorithm$lambda_rule, "lambda.1se") &&
    identical(model$algorithm$threshold, 0.5) &&
    identical(model$algorithm$outer_folds, 5L) &&
    identical(model$algorithm$inner_folds, 5L) &&
    identical(model$algorithm$final_folds, 10L),
  c(
    model$random_seed, model$algorithm$family, model$algorithm$alpha,
    model$algorithm$lambda_rule, model$algorithm$threshold,
    model$algorithm$outer_folds, model$algorithm$inner_folds,
    model$algorithm$final_folds
  ),
  c(42L, "binomial", 1, "lambda.1se", 0.5, 5L, 5L, 10L)
)
add_check(
  "frozen_feature_universe_preserved",
  nrow(model$feature_universe) == 521L &&
    !anyDuplicated(model$feature_universe$gene),
  c(nrow(model$feature_universe), length(unique(model$feature_universe$gene))),
  c(521L, 521L)
)
endpoint_counts <- table(
  factor(predictions$endpoint, levels = expected_endpoints)
)
add_check(
  "oof_predictions_complete_and_unique",
  identical(as.integer(endpoint_counts), c(171L, 183L, 265L)) &&
    !anyDuplicated(predictions[, c("endpoint", "sample_id")]),
  as.integer(endpoint_counts),
  c(171L, 183L, 265L)
)
add_check(
  "oof_probabilities_and_threshold_valid",
  all(is.finite(predictions$predicted_probability)) &&
    all(predictions$predicted_probability > 0) &&
    all(predictions$predicted_probability < 1) &&
    all(predictions$threshold == 0.5),
  c(
    min(predictions$predicted_probability),
    max(predictions$predicted_probability),
    unique(predictions$threshold)
  ),
  "finite probabilities strictly between 0 and 1; threshold 0.5"
)
predicted_positive <- vapply(
  expected_endpoints,
  function(endpoint) {
    sum(
      predictions$predicted_class[predictions$endpoint == endpoint] == 1L
    )
  },
  integer(1)
)
add_check(
  "frozen_threshold_degenerate_predictions_explicit",
  all(predicted_positive == c(0L, 183L, 265L)),
  predicted_positive,
  c(0L, 183L, 265L)
)
selected <- coefficients[
  coefficients$selected & coefficients$feature != "(Intercept)",
  ,
  drop = FALSE
]
selected_counts <- vapply(
  expected_endpoints,
  function(endpoint) sum(selected$endpoint == endpoint),
  integer(1)
)
add_check(
  "selected_feature_counts_reconciled",
  all(selected_counts == c(0L, 0L, 1L)) &&
    all(as.integer(summary$selected_features) == selected_counts),
  selected_counts,
  c(0L, 0L, 1L)
)
add_check(
  "composite_selected_feature_is_frozen_B4GALT7",
  nrow(selected) == 1L &&
    identical(selected$endpoint, "neoplasia_vs_normal") &&
    identical(selected$feature, "B4GALT7") &&
    selected$feature %in% model$feature_universe$gene,
  c(selected$endpoint, selected$feature),
  c("neoplasia_vs_normal", "B4GALT7")
)
nonestimable <- estimability[!estimability$observed_estimable, , drop = FALSE]
expected_nonestimable <- c(
  "adenoma_vs_normal:PPV",
  "cancer_vs_normal:NPV",
  "neoplasia_vs_normal:NPV"
)
observed_nonestimable <- sort(
  paste(nonestimable$endpoint, nonestimable$metric, sep = ":")
)
add_check(
  "undefined_predictive_values_are_mathematically_expected",
  identical(observed_nonestimable, sort(expected_nonestimable)) &&
    all(!nonestimable$expected_estimable),
  observed_nonestimable,
  sort(expected_nonestimable)
)
auc <- performance[performance$metric == "AUC", , drop = FALSE]
add_check(
  "nested_cv_auc_near_chance_and_CI_includes_half",
  nrow(auc) == 3L &&
    all(auc$estimate >= 0.45 & auc$estimate <= 0.55) &&
    all(auc$ci_low <= 0.5 & auc$ci_high >= 0.5),
  sprintf(
    "%s=%.4f[%.4f,%.4f]",
    auc$endpoint, auc$estimate, auc$ci_low, auc$ci_high
  ),
  "all AUC 0.45-0.55 and all 95% CIs include 0.5"
)
report_path <- file.path(project_dir, "reports", "stage_9B_model_training.md")
report_text <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
add_check(
  "report_discloses_negative_result_and_no_test_access",
  grepl("negative scientific result", report_text, fixed = TRUE) &&
    grepl("Independent test data: not extracted", report_text, fixed = TRUE) &&
    grepl("do not support useful discrimination", report_text, fixed = TRUE),
  "required disclosures searched",
  "negative result and test firewall disclosures present"
)

checks_df <- do.call(rbind, checks)
utils::write.table(
  checks_df,
  output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
if (!all(checks_df$passed)) {
  print(checks_df[!checks_df$passed, , drop = FALSE])
  stop("Stage 9B independent acceptance failed")
}
cat(
  "STAGE9B_INDEPENDENT_ACCEPTANCE_PASS ",
  sum(checks_df$passed), "/", nrow(checks_df), "\n",
  sep = ""
)
