#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09B_validate_locked_model.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09B_model_training", run_id)
model_path <- file.path(project_dir, "objects", "locked_stool_model.rds")

read_tsv <- function(name) {
  utils::read.delim(
    file.path(result_dir, name),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected),
    stringsAsFactors = FALSE
  )
}

model <- readRDS(model_path)
manifest <- utils::read.delim(
  file.path(project_dir, "metadata", "dataset_manifest.tsv"),
  check.names = FALSE
)
manifest <- manifest[manifest$accession == "GSE99573", , drop = FALSE]
training_ids <- sort(
  manifest$sample_id[manifest$validation_split == "training"]
)
test_ids <- sort(
  manifest$sample_id[manifest$validation_split == "testing"]
)
predictions <- read_tsv("nested_cv_predictions_training_only.tsv")
performance <- read_tsv("nested_cv_performance.tsv")
coefficients <- read_tsv("locked_model_coefficients.tsv")
features <- read_tsv("locked_feature_manifest.tsv")
summary <- read_tsv("model_lock_summary.tsv")
folds <- read_tsv("fold_assignments.tsv")
lock_manifest <- read_tsv("model_lock_manifest.tsv")

add_check(
  "model_file_exists_nonempty",
  file.exists(model_path) && file.info(model_path)$size > 1000,
  file.info(model_path)$size,
  ">1000"
)
mode <- file.info(model_path)$mode
add_check(
  "model_file_read_only",
  bitwAnd(as.integer(mode), as.integer(strtoi("222", base = 8))) == 0L,
  as.character(mode),
  "no write bits"
)
add_check(
  "artifact_identity",
  identical(model$artifact, "CRC_carcinogenesis_locked_stool_model") &&
    identical(model$artifact_version, "1.0.0"),
  paste(model$artifact, model$artifact_version),
  "CRC_carcinogenesis_locked_stool_model 1.0.0"
)
add_check(
  "test_expression_not_accessed",
  identical(model$test_expression_accessed, FALSE),
  model$test_expression_accessed,
  FALSE
)
add_check(
  "training_sample_count_265",
  identical(model$training_sample_count, 265L),
  model$training_sample_count,
  265L
)
add_check(
  "feature_universe_521_unique",
  nrow(model$feature_universe) == 521L &&
    !anyDuplicated(model$feature_universe$gene),
  paste(
    nrow(model$feature_universe),
    length(unique(model$feature_universe$gene)),
    sep = "/"
  ),
  "521/521"
)
add_check(
  "three_endpoint_models",
  length(model$models) == 3L &&
    identical(
      names(model$models),
      c("adenoma_vs_normal", "cancer_vs_normal", "neoplasia_vs_normal")
    ),
  paste(names(model$models), collapse = ";"),
  "three prespecified endpoints"
)
add_check(
  "primary_endpoint_frozen",
  identical(
    model$endpoint_hierarchy$role[
      model$endpoint_hierarchy$endpoint == "adenoma_vs_normal"
    ],
    "primary"
  ),
  model$endpoint_hierarchy$endpoint[
    model$endpoint_hierarchy$role == "primary"
  ],
  "adenoma_vs_normal"
)
add_check(
  "algorithm_frozen",
  identical(model$algorithm$alpha, 1) &&
    identical(model$algorithm$lambda_rule, "lambda.1se") &&
    identical(model$algorithm$threshold, 0.5),
  paste(
    model$algorithm$alpha,
    model$algorithm$lambda_rule,
    model$algorithm$threshold,
    sep = "/"
  ),
  "1/lambda.1se/0.5"
)
add_check(
  "fold_rules_frozen",
  identical(model$algorithm$outer_folds, 5L) &&
    identical(model$algorithm$inner_folds, 5L) &&
    identical(model$algorithm$final_folds, 10L),
  paste(
    model$algorithm$outer_folds,
    model$algorithm$inner_folds,
    model$algorithm$final_folds,
    sep = "/"
  ),
  "5/5/10"
)
add_check(
  "prediction_row_count",
  nrow(predictions) == 619L,
  nrow(predictions),
  619L
)
add_check(
  "prediction_ids_training_only",
  all(predictions$sample_id %in% training_ids) &&
    !any(predictions$sample_id %in% test_ids),
  paste(
    sum(predictions$sample_id %in% training_ids),
    sum(predictions$sample_id %in% test_ids),
    sep = "/"
  ),
  "619/0"
)
add_check(
  "one_oof_prediction_per_endpoint_sample",
  !anyDuplicated(predictions[, c("endpoint", "sample_id")]),
  nrow(unique(predictions[, c("endpoint", "sample_id")])),
  619L
)
add_check(
  "prediction_probabilities_valid",
  all(is.finite(predictions$predicted_probability)) &&
    all(predictions$predicted_probability > 0) &&
    all(predictions$predicted_probability < 1),
  sprintf(
    "%.6f to %.6f",
    min(predictions$predicted_probability),
    max(predictions$predicted_probability)
  ),
  "strictly between 0 and 1"
)
add_check(
  "thresholds_fixed_0_5",
  all(predictions$threshold == 0.5) &&
    all(summary$threshold == 0.5),
  paste(unique(c(predictions$threshold, summary$threshold)), collapse = ";"),
  "0.5"
)
add_check(
  "performance_rows_complete_unique",
  nrow(performance) == 33L &&
    !anyDuplicated(performance[, c("endpoint", "metric")]),
  paste(
    nrow(performance),
    nrow(unique(performance[, c("endpoint", "metric")])),
    sep = "/"
  ),
  "33/33"
)

# PPV is mathematically undefined when no sample is predicted positive, and
# NPV is undefined when no sample is predicted negative. These are legitimate
# frozen-threshold results, not missing analysis outputs.
prediction_support <- do.call(
  rbind,
  lapply(split(predictions, predictions$endpoint), function(z) {
    data.frame(
      endpoint = z$endpoint[[1L]],
      predicted_negative = sum(z$predicted_class == 0L),
      predicted_positive = sum(z$predicted_class == 1L),
      stringsAsFactors = FALSE
    )
  })
)
rownames(prediction_support) <- NULL
performance$expected_estimable <- TRUE
for (i in seq_len(nrow(performance))) {
  support <- prediction_support[
    prediction_support$endpoint == performance$endpoint[[i]],
    ,
    drop = FALSE
  ]
  if (performance$metric[[i]] == "PPV") {
    performance$expected_estimable[[i]] <- support$predicted_positive > 0L
  } else if (performance$metric[[i]] == "NPV") {
    performance$expected_estimable[[i]] <- support$predicted_negative > 0L
  }
}
performance$observed_estimable <-
  is.finite(performance$estimate) &
  is.finite(performance$ci_low) &
  is.finite(performance$ci_high)
performance$nonestimable_reason <- ifelse(
  performance$expected_estimable,
  "",
  ifelse(
    performance$metric == "PPV",
    "No predicted-positive samples at the prespecified threshold 0.5",
    "No predicted-negative samples at the prespecified threshold 0.5"
  )
)
estimability_audit <- merge(
  performance[
    ,
    c(
      "endpoint", "metric", "expected_estimable", "observed_estimable",
      "nonestimable_reason"
    )
  ],
  prediction_support,
  by = "endpoint",
  all.x = TRUE,
  sort = FALSE
)
write_tsv(
  estimability_audit,
  file.path(result_dir, "metric_estimability_audit.tsv")
)
add_check(
  "performance_defined_metrics_finite",
  all(performance$observed_estimable[performance$expected_estimable]),
  sum(
    performance$observed_estimable[performance$expected_estimable]
  ),
  sum(performance$expected_estimable)
)
add_check(
  "performance_estimability_matches_prediction_support",
  identical(
    performance$observed_estimable,
    performance$expected_estimable
  ),
  paste(
    sum(!performance$observed_estimable),
    paste(
      paste0(
        performance$endpoint[!performance$observed_estimable],
        ":",
        performance$metric[!performance$observed_estimable]
      ),
      collapse = ";"
    ),
    sep = "/"
  ),
  "3/adenoma_vs_normal:PPV;cancer_vs_normal:NPV;neoplasia_vs_normal:NPV"
)
add_check(
  "feature_coefficients_only_from_frozen_universe",
  all(
    coefficients$feature[coefficients$feature != "(Intercept)"] %in%
      model$feature_universe$gene
  ),
  length(unique(
    coefficients$feature[coefficients$feature != "(Intercept)"]
  )),
  521L
)
selected <- coefficients[
  coefficients$feature != "(Intercept)" & coefficients$selected,
  ,
  drop = FALSE
]
selected_counts <- vapply(
  summary$endpoint,
  function(endpoint) sum(selected$endpoint == endpoint),
  integer(1)
)
add_check(
  "selected_features_match_locked_summary_and_universe",
  length(selected_counts) == nrow(summary) &&
    all(selected_counts == as.integer(summary$selected_features)) &&
    all(selected$feature %in% model$feature_universe$gene),
  paste(selected_counts, collapse = "/"),
  paste(summary$selected_features, collapse = "/")
)
add_check(
  "endpoint_sample_counts",
  identical(summary$training_n, c(171L, 183L, 265L)),
  paste(summary$training_n, collapse = "/"),
  "171/183/265"
)
add_check(
  "all_training_arrays_retained_policy",
  grepl("All 265 training arrays retained", model$qc_policy, fixed = TRUE),
  model$qc_policy,
  "all 265 retained"
)
add_check(
  "fold_assignments_complete",
  nrow(folds) == 619L &&
    all(folds$outer_fold %in% seq_len(5L)) &&
    all(folds$final_10fold_id %in% seq_len(10L)),
  paste(
    nrow(folds),
    paste(sort(unique(folds$outer_fold)), collapse = ","),
    paste(sort(unique(folds$final_10fold_id)), collapse = ","),
    sep = "/"
  ),
  "619/1-5/1-10"
)
add_check(
  "lock_manifest_test_flag_false",
  identical(
    lock_manifest$value[
      lock_manifest$field == "test_expression_accessed"
    ],
    "FALSE"
  ),
  lock_manifest$value[
    lock_manifest$field == "test_expression_accessed"
  ],
  "FALSE"
)
add_check(
  "model_sha256_matches_manifest",
  identical(
    digest::digest(model_path, algo = "sha256", file = TRUE),
    lock_manifest$value[lock_manifest$field == "model_sha256"]
  ),
  digest::digest(model_path, algo = "sha256", file = TRUE),
  lock_manifest$value[lock_manifest$field == "model_sha256"]
)
add_check(
  "no_test_or_not_used_output_scope",
  !any(test_ids %in% unique(c(predictions$sample_id, folds$sample_id))) &&
    !any(
      manifest$sample_id[manifest$validation_split == "not_used"] %in%
        unique(c(predictions$sample_id, folds$sample_id))
    ),
  "0 forbidden sample IDs",
  "0"
)

checks_df <- do.call(rbind, checks)
write_tsv(
  checks_df,
  file.path(result_dir, "stage_9B_validation_checks.tsv")
)
if (!all(checks_df$passed)) {
  print(checks_df[!checks_df$passed, , drop = FALSE])
  stop("Stage 9B validation failed")
}
cat(
  "STAGE9B_VALIDATION_PASS ",
  sum(checks_df$passed), "/", nrow(checks_df), "\n",
  sep = ""
)
