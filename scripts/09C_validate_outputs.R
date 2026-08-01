#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09C_validate_outputs.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09C_external_test", run_id)
object_path <- file.path(
  project_dir, "objects", paste0("GSE99573_9C_test_RMA_", run_id, ".rds")
)
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
    observed = paste(observed, collapse = ";"),
    expected = paste(expected, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

lock_config <- utils::read.delim(
  file.path(project_dir, "config", "stage_9C_lock_verification.tsv"),
  check.names = FALSE
)
expected_sha <- lock_config$value[lock_config$field == "model_sha256"]
model_sha <- digest::digest(model_path, algo = "sha256", file = TRUE)
locked <- readRDS(model_path)
test_object <- readRDS(object_path)
predictions <- read_tsv("stool_test_predictions_source_data.tsv")
results <- read_tsv("stool_test_results.tsv")
fixed_specificity <- read_tsv("fixed_specificity_results.tsv")
feature_audit <- read_tsv("test_feature_reconstruction_audit.tsv")
qc <- read_tsv("test_sample_qc_metrics.tsv")
baseline <- read_tsv("clinical_baseline_comparison.tsv")
provenance <- read_tsv("stage_9C_provenance.tsv")
access <- read_tsv("test_access_transition.tsv")
operating_points <- utils::read.delim(
  file.path(project_dir, "config", "stool_test_operating_points.tsv"),
  check.names = FALSE
)
operating_points$threshold <- as.numeric(operating_points$threshold)
manifest <- utils::read.delim(
  file.path(project_dir, "metadata", "dataset_manifest.tsv"),
  check.names = FALSE
)
test_ids <- sort(manifest$sample_id[
  manifest$accession == "GSE99573" &
    manifest$validation_split == "testing"
])
training_ids <- manifest$sample_id[
  manifest$accession == "GSE99573" &
    manifest$validation_split == "training"
]
not_used_ids <- manifest$sample_id[
  manifest$accession == "GSE99573" &
    manifest$validation_split == "not_used"
]

add_check(
  "locked_model_sha256_unchanged",
  identical(model_sha, expected_sha),
  model_sha,
  expected_sha
)
mode <- file.info(model_path)$mode
add_check(
  "locked_model_still_read_only",
  bitwAnd(as.integer(mode), as.integer(strtoi("222", base = 8))) == 0L,
  mode,
  "no write bits"
)
add_check(
  "locked_model_identity_unchanged",
  identical(locked$artifact, "CRC_carcinogenesis_locked_stool_model") &&
    identical(locked$run_id, "20260730_001838") &&
    identical(locked$test_expression_accessed, FALSE),
  c(locked$artifact, locked$run_id, locked$test_expression_accessed),
  c(
    "CRC_carcinogenesis_locked_stool_model",
    "20260730_001838", FALSE
  )
)
add_check(
  "test_object_complete",
  identical(dim(test_object$expression), c(70523L, 65L)) &&
    identical(dim(test_object$gene_expression), c(521L, 65L)) &&
    all(is.finite(test_object$expression)) &&
    all(is.finite(test_object$gene_expression)),
  c(dim(test_object$expression), dim(test_object$gene_expression)),
  c(70523L, 65L, 521L, 65L)
)
add_check(
  "test_object_ids_exactly_frozen_test",
  setequal(test_object$sample_id, test_ids) &&
    !any(test_object$sample_id %in% c(training_ids, not_used_ids)),
  c(
    sum(test_object$sample_id %in% test_ids),
    sum(test_object$sample_id %in% training_ids),
    sum(test_object$sample_id %in% not_used_ids)
  ),
  c(65L, 0L, 0L)
)
add_check(
  "test_preprocessing_outcome_blind",
  identical(test_object$test_outcomes_used_during_preprocessing, FALSE) &&
    identical(test_object$automatic_exclusions, 0L),
  c(
    test_object$test_outcomes_used_during_preprocessing,
    test_object$automatic_exclusions
  ),
  c(FALSE, 0L)
)
endpoint_order <- c(
  "adenoma_vs_normal", "cancer_vs_normal", "neoplasia_vs_normal"
)
prediction_counts <- table(
  factor(predictions$endpoint, levels = endpoint_order)
)
add_check(
  "prediction_rows_complete",
  nrow(predictions) == 152L &&
    identical(as.integer(prediction_counts), c(42L, 45L, 65L)) &&
    !anyDuplicated(predictions[, c("endpoint", "sample_id")]),
  c(nrow(predictions), as.integer(prediction_counts)),
  c(152L, 42L, 45L, 65L)
)
add_check(
  "prediction_ids_test_only",
  all(predictions$sample_id %in% test_ids) &&
    !any(predictions$sample_id %in% c(training_ids, not_used_ids)),
  c(
    sum(predictions$sample_id %in% test_ids),
    sum(predictions$sample_id %in% training_ids),
    sum(predictions$sample_id %in% not_used_ids)
  ),
  c(152L, 0L, 0L)
)
add_check(
  "endpoint_case_control_counts",
  all(c(
    sum(
      predictions$endpoint == "adenoma_vs_normal" &
        predictions$outcome == 1L
    ) == 20L,
    sum(
      predictions$endpoint == "adenoma_vs_normal" &
        predictions$outcome == 0L
    ) == 22L,
    sum(
      predictions$endpoint == "cancer_vs_normal" &
        predictions$outcome == 1L
    ) == 23L,
    sum(
      predictions$endpoint == "cancer_vs_normal" &
        predictions$outcome == 0L
    ) == 22L,
    sum(
      predictions$endpoint == "neoplasia_vs_normal" &
        predictions$outcome == 1L
    ) == 43L,
    sum(
      predictions$endpoint == "neoplasia_vs_normal" &
        predictions$outcome == 0L
    ) == 22L
  )),
  "20/22;23/22;43/22",
  "20/22;23/22;43/22"
)
add_check(
  "probabilities_valid_and_threshold_frozen",
  all(is.finite(predictions$predicted_probability)) &&
    all(predictions$predicted_probability > 0) &&
    all(predictions$predicted_probability < 1) &&
    all(predictions$locked_threshold == 0.5),
  c(
    min(predictions$predicted_probability),
    max(predictions$predicted_probability),
    unique(predictions$locked_threshold)
  ),
  "finite probabilities strictly between 0 and 1; threshold 0.5"
)

# Recalculate every locked probability from the stored test gene matrix and
# locked coefficients; this validation does not fit or tune any model.
max_abs_error <- 0
for (endpoint in endpoint_order) {
  z <- predictions[predictions$endpoint == endpoint, , drop = FALSE]
  model <- locked$models[[endpoint]]
  x <- t(test_object$gene_expression[, z$sample_id, drop = FALSE])
  beta <- model$coefficients[colnames(x)]
  recalculated <- stats::plogis(
    as.numeric(model$intercept + x %*% beta)
  )
  max_abs_error <- max(
    max_abs_error,
    max(abs(recalculated - z$predicted_probability))
  )
}
add_check(
  "predictions_reproduce_locked_coefficients",
  is.finite(max_abs_error) && max_abs_error < 1e-12,
  format(max_abs_error, scientific = TRUE),
  "<1e-12"
)
add_check(
  "test_result_metric_inventory",
  nrow(results) == 27L &&
    all(table(results$endpoint) == 9L) &&
    setequal(
      unique(results$metric),
      c(
        "AUC", "sensitivity", "specificity", "PPV", "NPV", "accuracy",
        "Brier", "calibration_intercept", "calibration_slope"
      )
    ),
  c(nrow(results), paste(table(results$endpoint), collapse = "/")),
  "27;9/9/9"
)
auc <- results[results$metric == "AUC", , drop = FALSE]
add_check(
  "AUC_and_95CI_reported_for_all_endpoints",
  nrow(auc) == 3L &&
    all(is.finite(auc$estimate)) &&
    all(is.finite(auc$ci_low)) &&
    all(is.finite(auc$ci_high)) &&
    all(auc$ci_low <= auc$estimate & auc$estimate <= auc$ci_high),
  sprintf(
    "%s=%.4f[%.4f,%.4f]",
    auc$endpoint, auc$estimate, auc$ci_low, auc$ci_high
  ),
  "three finite ordered AUC estimates and 95% CIs"
)
classification <- results[
  results$metric %in% c(
    "sensitivity", "specificity", "PPV", "NPV", "accuracy"
  ),
  ,
  drop = FALSE
]
classification_estimable <- classification$denominator > 0L
add_check(
  "classification_CI_estimability_correct",
  all(is.finite(classification$estimate[classification_estimable])) &&
    all(is.finite(classification$ci_low[classification_estimable])) &&
    all(is.finite(classification$ci_high[classification_estimable])) &&
    all(!is.finite(
      classification$estimate[!classification_estimable]
    )),
  c(
    sum(classification_estimable),
    sum(is.finite(classification$estimate))
  ),
  "finite exactly when denominator is nonzero"
)
add_check(
  "fixed_specificity_rows_and_thresholds_frozen",
  nrow(fixed_specificity) == 3L &&
    all(fixed_specificity$threshold_source ==
      "Stage_9B_training_outer_OOF") &&
    all(
      fixed_specificity$training_locked_threshold ==
        operating_points$threshold[
          match(
            paste(
              fixed_specificity$endpoint,
              "training_locked_specificity_90"
            ),
            paste(
              operating_points$endpoint,
              operating_points$operating_point
            )
          )
        ] |
        (
          is.infinite(fixed_specificity$training_locked_threshold) &
          is.infinite(operating_points$threshold[
            match(
              paste(
                fixed_specificity$endpoint,
                "training_locked_specificity_90"
              ),
              paste(
                operating_points$endpoint,
                operating_points$operating_point
              )
            )
          ])
        )
    ),
  paste(fixed_specificity$training_locked_threshold, collapse = "/"),
  "Inf/Inf/0.67534779307876"
)
add_check(
  "all_521_features_reconstructed_without_substitution",
  nrow(feature_audit) == 521L &&
    all(feature_audit$available_transcript_clusters > 0L) &&
    all(
      feature_audit$reconstruction_status ==
        "reconstructed_without_substitution"
    ),
  c(
    nrow(feature_audit),
    min(feature_audit$available_transcript_clusters)
  ),
  "521;minimum available clusters >=1"
)
add_check(
  "all_65_test_arrays_retained",
  nrow(qc) == 65L &&
    !any(qc$automatic_exclusion) &&
    !any(qc$outcome_labels_used_for_QC),
  c(
    nrow(qc), sum(qc$automatic_exclusion),
    sum(qc$outcome_labels_used_for_QC)
  ),
  c(65L, 0L, 0L)
)
add_check(
  "clinical_baseline_unavailability_explicit",
  nrow(baseline) == 1L &&
    !baseline$clinical_baseline_available &&
    identical(baseline$status, "not_evaluable") &&
    !baseline$test_outcomes_used_to_build_baseline,
  c(
    baseline$clinical_baseline_available,
    baseline$status,
    baseline$test_outcomes_used_to_build_baseline
  ),
  c(FALSE, "not_evaluable", FALSE)
)
add_check(
  "one_time_access_transition_complete",
  nrow(access) == 65L &&
    setequal(access$sample_id, test_ids) &&
    all(access$stage_9C_CEL_extracted) &&
    all(access$stage_9C_expression_normalized) &&
    all(access$stage_9C_prediction_generated) &&
    !any(access$stage_9C_excluded),
  c(
    nrow(access),
    sum(access$stage_9C_prediction_generated),
    sum(access$stage_9C_excluded)
  ),
  c(65L, 65L, 0L)
)
add_check(
  "no_retraining_or_lock_change",
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
      provenance$value[provenance$field == "model_sha256_before"],
      provenance$value[provenance$field == "model_sha256_after"]
    ),
  paste(
    provenance$value[provenance$field %in% c(
      "model_retrained", "features_changed", "threshold_changed"
    )],
    collapse = "/"
  ),
  "FALSE/FALSE/FALSE"
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
    "fixed_specificity_results.tsv", "test_sample_qc_metrics.tsv",
    "test_feature_reconstruction_audit.tsv",
    "clinical_baseline_comparison.tsv", "test_access_transition.tsv"
  )
)
add_check(
  "figures_and_source_data_complete",
  all(file.exists(figure_paths)) &&
    all(file.info(figure_paths)$size > 1000L) &&
    all(file.exists(source_paths)) &&
    all(file.info(source_paths)$size > 0L),
  c(sum(file.exists(figure_paths)), sum(file.exists(source_paths))),
  c(4L, 8L)
)

checks_df <- do.call(rbind, checks)
write_tsv(
  checks_df,
  file.path(result_dir, "stage_9C_validation_checks.tsv")
)
if (!all(checks_df$passed)) {
  print(checks_df[!checks_df$passed, , drop = FALSE])
  stop("Stage 9C validation failed")
}
cat(
  "STAGE9C_VALIDATION_PASS ",
  sum(checks_df$passed), "/", nrow(checks_df), "\n",
  sep = ""
)
