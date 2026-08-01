#!/usr/bin/env Rscript

# Analysis: one-time independent evaluation of locked GSE99573 stool models
# Date: 2026-07-30
# Random seed: 42
# Unit: participant

set.seed(42)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 09C_test_once.R <project_dir> <run_id> <setup_git_commit>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
setup_git_commit <- args[[3]]
lib_dir <- file.path(project_dir, "environment", "Rlib_stage9A")
.libPaths(c(lib_dir, .libPaths()))

required <- c(
  "oligo", "pd.hta.2.0", "Biobase", "matrixStats", "pROC", "ggplot2",
  "digest"
)
missing_packages <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

result_dir <- file.path(project_dir, "results", "09C_external_test", run_id)
figure_dir <- file.path(project_dir, "figures", "09C_external_test", run_id)
work_dir <- file.path(project_dir, "data_processed", "09C_external_test", run_id)
cel_dir <- file.path(work_dir, "test_cel")
object_path <- file.path(
  project_dir, "objects", paste0("GSE99573_9C_test_RMA_", run_id, ".rds")
)
started_marker <- file.path(
  project_dir, "logs", "09C_external_test", "ONE_TIME_TEST_STARTED"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(started_marker)) {
  stop("One-time test start marker is missing")
}
if (file.exists(file.path(result_dir, "stool_test_results.tsv"))) {
  stop("Refusing to overwrite an existing one-time test result")
}
if (file.exists(object_path)) {
  stop("Refusing to overwrite an existing test RMA object")
}

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
clamp_probability <- function(p) {
  pmin(pmax(as.numeric(p), 1e-8), 1 - 1e-8)
}
wilson_ci <- function(successes, denominator, conf = 0.95) {
  if (!is.finite(denominator) || denominator <= 0L) {
    return(c(NA_real_, NA_real_))
  }
  p <- successes / denominator
  z <- stats::qnorm(1 - (1 - conf) / 2)
  denom <- 1 + z^2 / denominator
  center <- (p + z^2 / (2 * denominator)) / denom
  half <- z * sqrt(
    p * (1 - p) / denominator + z^2 / (4 * denominator^2)
  ) / denom
  c(max(0, center - half), min(1, center + half))
}
auc_with_ci <- function(y, p) {
  if (length(unique(p)) == 1L) {
    return(c(
      estimate = 0.5, ci_low = 0.5, ci_high = 0.5,
      method_code = 1
    ))
  }
  roc <- pROC::roc(
    response = y, predictor = p, levels = c(0, 1),
    direction = "<", quiet = TRUE
  )
  ci <- as.numeric(pROC::ci.auc(roc, method = "delong"))
  c(
    estimate = as.numeric(pROC::auc(roc)),
    ci_low = ci[[1L]], ci_high = ci[[3L]], method_code = 0
  )
}
calibration_values <- function(y, p) {
  p <- clamp_probability(p)
  lp <- stats::qlogis(p)
  intercept_fit <- try(
    suppressWarnings(
      stats::glm(y ~ 1L + offset(lp), family = stats::binomial())
    ),
    silent = TRUE
  )
  slope_fit <- try(
    suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())),
    silent = TRUE
  )
  intercept <- if (inherits(intercept_fit, "try-error")) {
    NA_real_
  } else {
    unname(stats::coef(intercept_fit)[[1L]])
  }
  slope <- if (inherits(slope_fit, "try-error")) {
    NA_real_
  } else {
    unname(stats::coef(slope_fit)[["lp"]])
  }
  c(calibration_intercept = intercept, calibration_slope = slope)
}
bootstrap_continuous_metrics <- function(y, p, n_boot, seed) {
  point <- c(
    Brier = mean((p - y)^2),
    calibration_values(y, p)
  )
  set.seed(seed)
  by_class <- split(seq_along(y), y)
  boot <- matrix(
    NA_real_, nrow = n_boot, ncol = length(point),
    dimnames = list(NULL, names(point))
  )
  for (b in seq_len(n_boot)) {
    idx <- unlist(
      lapply(by_class, function(z) sample(z, length(z), replace = TRUE)),
      use.names = FALSE
    )
    boot[b, ] <- c(
      Brier = mean((p[idx] - y[idx])^2),
      calibration_values(y[idx], p[idx])
    )
  }
  do.call(
    rbind,
    lapply(names(point), function(metric) {
      finite <- is.finite(boot[, metric])
      ci <- if (sum(finite) >= 100L) {
        stats::quantile(
          boot[finite, metric], c(0.025, 0.975),
          names = FALSE, type = 8
        )
      } else {
        c(NA_real_, NA_real_)
      }
      data.frame(
        metric = metric,
        estimate = unname(point[[metric]]),
        ci_low = ci[[1L]],
        ci_high = ci[[2L]],
        ci_method = paste0(
          "stratified bootstrap; ", n_boot,
          " replicates; finite=", sum(finite)
        ),
        stringsAsFactors = FALSE
      )
    })
  )
}
classification_metrics <- function(y, p, threshold) {
  predicted <- as.integer(p >= threshold)
  tp <- sum(predicted == 1L & y == 1L)
  tn <- sum(predicted == 0L & y == 0L)
  fp <- sum(predicted == 1L & y == 0L)
  fn <- sum(predicted == 0L & y == 1L)
  definitions <- list(
    sensitivity = c(tp, tp + fn),
    specificity = c(tn, tn + fp),
    PPV = c(tp, tp + fp),
    NPV = c(tn, tn + fn),
    accuracy = c(tp + tn, length(y))
  )
  rows <- lapply(names(definitions), function(metric) {
    z <- definitions[[metric]]
    ci <- wilson_ci(z[[1L]], z[[2L]])
    data.frame(
      metric = metric,
      estimate = if (z[[2L]] > 0L) z[[1L]] / z[[2L]] else NA_real_,
      ci_low = ci[[1L]],
      ci_high = ci[[2L]],
      ci_method = if (z[[2L]] > 0L) {
        "Wilson score"
      } else {
        "not estimable; zero denominator"
      },
      numerator = z[[1L]],
      denominator = z[[2L]],
      stringsAsFactors = FALSE
    )
  })
  list(
    table = do.call(rbind, rows),
    predicted = predicted,
    counts = c(tp = tp, tn = tn, fp = fp, fn = fn)
  )
}

lock_config <- utils::read.delim(
  file.path(project_dir, "config", "stage_9C_lock_verification.tsv"),
  check.names = FALSE
)
expected_model_sha <- lock_config$value[
  lock_config$field == "model_sha256"
]
model_path <- file.path(project_dir, "objects", "locked_stool_model.rds")
model_sha_before <- digest::digest(model_path, algo = "sha256", file = TRUE)
if (!identical(model_sha_before, expected_model_sha)) {
  stop("Locked model SHA256 mismatch before test evaluation")
}
locked <- readRDS(model_path)
if (!identical(locked$test_expression_accessed, FALSE) ||
    !identical(locked$run_id, "20260730_001838")) {
  stop("Locked model identity or test firewall flag is invalid")
}

cel_files <- sort(list.files(
  cel_dir, pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE
))
test_ids_from_files <- sub("_.*$", "", basename(cel_files))
if (length(cel_files) != 65L || anyDuplicated(test_ids_from_files)) {
  stop("Expected exactly 65 unique uncompressed test CEL files")
}

message("Reading the authorized 65-sample independent test batch")
raw_data <- oligo::read.celfiles(
  cel_files, pkgname = "pd.hta.2.0", verbose = TRUE
)
message("Running independent test-batch core RMA")
normalized <- oligo::rma(
  raw_data, target = "core", background = TRUE, normalize = TRUE
)
expression <- Biobase::exprs(normalized)
colnames(expression) <- test_ids_from_files
if (!identical(dim(expression), c(70523L, 65L)) ||
    any(!is.finite(expression))) {
  stop("Invalid independent test RMA matrix")
}
rm(raw_data, normalized)
invisible(gc())

# Outcome-blind test-array diagnostics; no sample is excluded.
row_medians <- matrixStats::rowMedians(expression, na.rm = TRUE)
rle <- sweep(expression, 1L, row_medians, FUN = "-")
row_variances <- matrixStats::rowVars(expression, na.rm = TRUE)
top_rows <- order(row_variances, decreasing = TRUE)[seq_len(5000L)]
sample_cor <- stats::cor(
  expression[top_rows, , drop = FALSE],
  method = "pearson", use = "pairwise.complete.obs"
)
qc <- data.frame(
  sample_id = colnames(expression),
  array_median = matrixStats::colMedians(expression, na.rm = TRUE),
  array_IQR = apply(expression, 2L, stats::IQR, na.rm = TRUE),
  RLE_median = matrixStats::colMedians(rle, na.rm = TRUE),
  RLE_IQR = apply(rle, 2L, stats::IQR, na.rm = TRUE),
  mean_interarray_correlation =
    (rowSums(sample_cor) - 1) / (ncol(sample_cor) - 1),
  automatic_exclusion = FALSE,
  outcome_labels_used_for_QC = FALSE,
  stringsAsFactors = FALSE
)
write_tsv(qc, file.path(result_dir, "test_sample_qc_metrics.tsv"))
rm(rle, sample_cor)
invisible(gc())

feature_map <- locked$feature_universe
if (nrow(feature_map) != 521L || anyDuplicated(feature_map$gene)) {
  stop("Locked feature universe is not 521 unique genes")
}
gene_expression <- matrix(
  NA_real_, nrow = nrow(feature_map), ncol = ncol(expression),
  dimnames = list(feature_map$gene, colnames(expression))
)
mapping_audit <- vector("list", nrow(feature_map))
for (i in seq_len(nrow(feature_map))) {
  probe_ids <- strsplit(
    feature_map$transcript_cluster_ids[[i]], ";", fixed = TRUE
  )[[1L]]
  available <- intersect(probe_ids, rownames(expression))
  if (!length(available)) {
    stop(
      "Frozen test feature cannot be reconstructed: ",
      feature_map$gene[[i]]
    )
  }
  values <- expression[available, , drop = FALSE]
  gene_expression[i, ] <- if (nrow(values) == 1L) {
    as.numeric(values[1L, ])
  } else {
    matrixStats::colMedians(values, na.rm = TRUE)
  }
  mapping_audit[[i]] <- data.frame(
    gene = feature_map$gene[[i]],
    frozen_transcript_clusters = length(probe_ids),
    available_transcript_clusters = length(available),
    reconstruction_status = "reconstructed_without_substitution",
    stringsAsFactors = FALSE
  )
}
mapping_audit <- do.call(rbind, mapping_audit)
if (any(!is.finite(gene_expression))) {
  stop("Nonfinite frozen test gene expression")
}
write_tsv(
  mapping_audit,
  file.path(result_dir, "test_feature_reconstruction_audit.tsv")
)

saveRDS(
  list(
    expression = expression,
    gene_expression = gene_expression,
    sample_id = colnames(expression),
    platform = "GPL17586",
    preprocessing = paste(
      "oligo::rma target=core; background=TRUE; normalize=TRUE;",
      "65-sample independent test batch"
    ),
    model_sha256 = model_sha_before,
    test_expression_accessed = TRUE,
    test_outcomes_used_during_preprocessing = FALSE,
    automatic_exclusions = 0L,
    seed = 42L
  ),
  object_path,
  compress = "xz"
)
rm(expression)
invisible(gc())

manifest <- utils::read.delim(
  file.path(project_dir, "metadata", "dataset_manifest.tsv"),
  check.names = FALSE
)
test_meta <- manifest[
  manifest$accession == "GSE99573" &
    manifest$validation_split == "testing",
  c("sample_id", "donor_id", "condition", "histology"),
  drop = FALSE
]
if (nrow(test_meta) != 65L ||
    anyDuplicated(test_meta$sample_id) ||
    !setequal(test_meta$sample_id, colnames(gene_expression))) {
  stop("Frozen test metadata does not match normalized test arrays")
}
test_meta <- test_meta[
  match(colnames(gene_expression), test_meta$sample_id),
  ,
  drop = FALSE
]
if (!identical(
  as.integer(table(factor(
    test_meta$condition, levels = c("normal", "adenoma", "cancer")
  ))),
  c(22L, 20L, 23L)
)) {
  stop("Unexpected test condition inventory")
}

endpoint_definitions <- list(
  adenoma_vs_normal = list(
    role = "primary", included = c("normal", "adenoma"),
    case = "adenoma", label = "Adenoma versus normal"
  ),
  cancer_vs_normal = list(
    role = "secondary", included = c("normal", "cancer"),
    case = "cancer", label = "Colorectal cancer versus normal"
  ),
  neoplasia_vs_normal = list(
    role = "secondary_screening_composite",
    included = c("normal", "adenoma", "cancer"),
    case = c("adenoma", "cancer"),
    label = "Adenoma or colorectal cancer versus normal"
  )
)
operating_points <- utils::read.delim(
  file.path(project_dir, "config", "stool_test_operating_points.tsv"),
  check.names = FALSE
)
operating_points$threshold <- as.numeric(operating_points$threshold)

all_predictions <- list()
all_results <- list()
all_fixed_specificity <- list()
all_roc <- list()
all_calibration <- list()
for (endpoint_index in seq_along(endpoint_definitions)) {
  endpoint <- names(endpoint_definitions)[[endpoint_index]]
  definition <- endpoint_definitions[[endpoint]]
  endpoint_model <- locked$models[[endpoint]]
  keep <- test_meta$condition %in% definition$included
  meta <- test_meta[keep, , drop = FALSE]
  y <- as.integer(meta$condition %in% definition$case)
  x <- t(gene_expression[, meta$sample_id, drop = FALSE])
  beta <- endpoint_model$coefficients
  if (!setequal(names(beta), rownames(gene_expression))) {
    stop("Coefficient-feature mismatch for ", endpoint)
  }
  beta <- beta[colnames(x)]
  linear_predictor <- as.numeric(
    endpoint_model$intercept + x %*% beta
  )
  probability <- stats::plogis(linear_predictor)
  if (any(!is.finite(probability))) {
    stop("Nonfinite locked prediction for ", endpoint)
  }
  locked_metrics <- classification_metrics(
    y, probability, endpoint_model$threshold
  )
  predictions <- data.frame(
    endpoint = endpoint,
    endpoint_role = definition$role,
    sample_id = meta$sample_id,
    donor_id = meta$donor_id,
    condition = meta$condition,
    outcome = y,
    predicted_probability = probability,
    locked_threshold = endpoint_model$threshold,
    locked_predicted_class = locked_metrics$predicted,
    model_sha256 = model_sha_before,
    evaluation_scope = "one_time_independent_test",
    stringsAsFactors = FALSE
  )

  class_table <- locked_metrics$table
  class_table$endpoint <- endpoint
  class_table$endpoint_role <- definition$role
  class_table$threshold <- endpoint_model$threshold
  class_table$evaluation_scope <- "one_time_independent_test"

  auc <- auc_with_ci(y, probability)
  auc_row <- data.frame(
    metric = "AUC",
    estimate = auc[["estimate"]],
    ci_low = auc[["ci_low"]],
    ci_high = auc[["ci_high"]],
    ci_method = if (auc[["method_code"]] == 1) {
      "constant-score exact conditional AUC"
    } else {
      "DeLong"
    },
    numerator = NA_integer_,
    denominator = length(y),
    endpoint = endpoint,
    endpoint_role = definition$role,
    threshold = NA_real_,
    evaluation_scope = "one_time_independent_test",
    stringsAsFactors = FALSE
  )
  continuous <- bootstrap_continuous_metrics(
    y, probability, 2000L, 920000L + endpoint_index
  )
  continuous$numerator <- NA_integer_
  continuous$denominator <- length(y)
  continuous$endpoint <- endpoint
  continuous$endpoint_role <- definition$role
  continuous$threshold <- NA_real_
  continuous$evaluation_scope <- "one_time_independent_test"
  result_rows <- rbind(
    auc_row,
    class_table[
      ,
      c(
        "metric", "estimate", "ci_low", "ci_high", "ci_method",
        "numerator", "denominator", "endpoint", "endpoint_role",
        "threshold", "evaluation_scope"
      )
    ],
    continuous[
      ,
      c(
        "metric", "estimate", "ci_low", "ci_high", "ci_method",
        "numerator", "denominator", "endpoint", "endpoint_role",
        "threshold", "evaluation_scope"
      )
    ]
  )

  high_spec <- operating_points[
    operating_points$endpoint == endpoint &
      operating_points$operating_point == "training_locked_specificity_90",
    ,
    drop = FALSE
  ]
  if (nrow(high_spec) != 1L) {
    stop("Missing frozen high-specificity operating point for ", endpoint)
  }
  high_spec_metrics <- classification_metrics(
    y, probability, high_spec$threshold[[1L]]
  )
  predictions$specificity90_threshold <- high_spec$threshold[[1L]]
  predictions$specificity90_predicted_class <-
    high_spec_metrics$predicted
  sens <- high_spec_metrics$table[
    high_spec_metrics$table$metric == "sensitivity", , drop = FALSE
  ]
  spec <- high_spec_metrics$table[
    high_spec_metrics$table$metric == "specificity", , drop = FALSE
  ]
  fixed_row <- data.frame(
    endpoint = endpoint,
    endpoint_role = definition$role,
    target_specificity = high_spec$target_specificity,
    training_locked_threshold = high_spec$threshold,
    training_specificity = high_spec$training_specificity,
    training_sensitivity = high_spec$training_sensitivity,
    test_specificity = spec$estimate,
    test_specificity_ci_low = spec$ci_low,
    test_specificity_ci_high = spec$ci_high,
    test_sensitivity = sens$estimate,
    test_sensitivity_ci_low = sens$ci_low,
    test_sensitivity_ci_high = sens$ci_high,
    threshold_source = "Stage_9B_training_outer_OOF",
    stringsAsFactors = FALSE
  )

  if (length(unique(probability)) == 1L) {
    roc_source <- data.frame(
      false_positive_rate = c(0, 1),
      true_positive_rate = c(0, 1),
      threshold = c(Inf, -Inf),
      endpoint = endpoint,
      stringsAsFactors = FALSE
    )
  } else {
    roc <- pROC::roc(
      response = y, predictor = probability, levels = c(0, 1),
      direction = "<", quiet = TRUE
    )
    coords <- pROC::coords(
      roc, x = "all",
      ret = c("threshold", "specificity", "sensitivity"),
      transpose = FALSE
    )
    roc_source <- data.frame(
      false_positive_rate = 1 - coords$specificity,
      true_positive_rate = coords$sensitivity,
      threshold = coords$threshold,
      endpoint = endpoint,
      stringsAsFactors = FALSE
    )
  }

  if (length(unique(probability)) < 3L) {
    calibration <- data.frame(
      bin = "all",
      mean_predicted_probability = mean(probability),
      observed_case_fraction = mean(y),
      n = length(y),
      endpoint = endpoint,
      stringsAsFactors = FALSE
    )
  } else {
    breaks <- unique(as.numeric(stats::quantile(
      probability, probs = seq(0, 1, length.out = 11),
      names = FALSE, type = 8
    )))
    bins <- cut(probability, breaks = breaks, include.lowest = TRUE)
    calibration <- stats::aggregate(
      cbind(probability, y) ~ bins,
      FUN = mean
    )
    n_by_bin <- as.data.frame(table(bins), stringsAsFactors = FALSE)
    calibration$n <- n_by_bin$Freq[
      match(as.character(calibration$bins), as.character(n_by_bin$bins))
    ]
    names(calibration)[names(calibration) == "bins"] <- "bin"
    names(calibration)[names(calibration) == "probability"] <-
      "mean_predicted_probability"
    names(calibration)[names(calibration) == "y"] <-
      "observed_case_fraction"
    calibration$endpoint <- endpoint
  }

  all_predictions[[endpoint]] <- predictions
  all_results[[endpoint]] <- result_rows
  all_fixed_specificity[[endpoint]] <- fixed_row
  all_roc[[endpoint]] <- roc_source
  all_calibration[[endpoint]] <- calibration
}

predictions <- do.call(rbind, all_predictions)
test_results <- do.call(rbind, all_results)
fixed_specificity <- do.call(rbind, all_fixed_specificity)
roc_source <- do.call(rbind, all_roc)
calibration_source <- do.call(rbind, all_calibration)
rownames(predictions) <- NULL
rownames(test_results) <- NULL

write_tsv(
  test_results,
  file.path(result_dir, "stool_test_results.tsv")
)
write_tsv(
  predictions,
  file.path(result_dir, "stool_test_predictions_source_data.tsv")
)
write_tsv(
  fixed_specificity,
  file.path(result_dir, "fixed_specificity_results.tsv")
)
write_tsv(
  roc_source,
  file.path(result_dir, "ROC_source_data.tsv")
)
write_tsv(
  calibration_source,
  file.path(result_dir, "calibration_source_data.tsv")
)
write_tsv(
  data.frame(
    comparison = "locked_molecular_model_vs_clinical_baseline",
    clinical_baseline_available = FALSE,
    status = "not_evaluable",
    reason = paste(
      "Public project metadata provide no non-outcome baseline predictors",
      "(age, sex, FIT or equivalent); disease category and cancer stage",
      "define or reveal the endpoint and were not used as predictors."
    ),
    test_outcomes_used_to_build_baseline = FALSE,
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "clinical_baseline_comparison.tsv")
)
write_tsv(
  data.frame(
    sample_id = test_meta$sample_id,
    prior_stage_9A_expression_accessed = FALSE,
    stage_9C_CEL_extracted = TRUE,
    stage_9C_expression_normalized = TRUE,
    stage_9C_prediction_generated = TRUE,
    stage_9C_excluded = FALSE,
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "test_access_transition.tsv")
)

endpoint_labels <- c(
  adenoma_vs_normal = "Adenoma vs normal",
  cancer_vs_normal = "Cancer vs normal",
  neoplasia_vs_normal = "Adenoma/cancer vs normal"
)
plot_roc <- roc_source
plot_roc$endpoint_label <- factor(
  plot_roc$endpoint,
  levels = names(endpoint_labels),
  labels = endpoint_labels
)
colors <- c("#0072B2", "#D55E00", "#009E73")
roc_plot <- ggplot2::ggplot(
  plot_roc,
  ggplot2::aes(
    false_positive_rate, true_positive_rate, color = endpoint_label
  )
) +
  ggplot2::geom_path(linewidth = 0.8) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", color = "grey55"
  ) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Locked stool models in the independent test set",
    x = "False-positive rate",
    y = "True-positive rate",
    color = "Endpoint"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )
ggplot2::ggsave(
  file.path(figure_dir, "stool_test_ROC.pdf"),
  roc_plot, width = 7, height = 6, units = "in"
)
ggplot2::ggsave(
  file.path(figure_dir, "stool_test_ROC.png"),
  roc_plot, width = 7, height = 6, units = "in", dpi = 300
)

plot_calibration <- calibration_source
plot_calibration$endpoint_label <- factor(
  plot_calibration$endpoint,
  levels = names(endpoint_labels),
  labels = endpoint_labels
)
calibration_plot <- ggplot2::ggplot(
  plot_calibration,
  ggplot2::aes(
    mean_predicted_probability, observed_case_fraction,
    color = endpoint_label
  )
) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", color = "grey55"
  ) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.85) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Independent-test calibration",
    x = "Mean predicted probability",
    y = "Observed case fraction",
    color = "Endpoint",
    size = "Participants"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )
ggplot2::ggsave(
  file.path(figure_dir, "stool_test_calibration.pdf"),
  calibration_plot, width = 7, height = 6, units = "in"
)
ggplot2::ggsave(
  file.path(figure_dir, "stool_test_calibration.png"),
  calibration_plot, width = 7, height = 6, units = "in", dpi = 300
)

software <- data.frame(
  package = c("R", required),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(
      required,
      function(x) as.character(utils::packageVersion(x)),
      character(1)
    )
  ),
  stringsAsFactors = FALSE
)
write_tsv(software, file.path(result_dir, "software_versions.tsv"))

model_sha_after <- digest::digest(model_path, algo = "sha256", file = TRUE)
if (!identical(model_sha_before, model_sha_after)) {
  stop("Locked model changed during independent evaluation")
}
write_tsv(
  data.frame(
    field = c(
      "run_id", "setup_git_commit", "model_sha256_before",
      "model_sha256_after", "test_samples", "automatic_exclusions",
      "model_retrained", "features_changed", "threshold_changed",
      "clinical_baseline_available"
    ),
    value = c(
      run_id, setup_git_commit, model_sha_before, model_sha_after,
      "65", "0", "FALSE", "FALSE", "FALSE", "FALSE"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "stage_9C_provenance.tsv")
)

cat("STAGE9C_ONE_TIME_TEST_ANALYSIS_OK run_id=", run_id, "\n", sep = "")
