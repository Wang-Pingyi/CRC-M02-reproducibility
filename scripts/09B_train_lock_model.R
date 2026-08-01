#!/usr/bin/env Rscript

set.seed(42)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 09B_train_lock_model.R <project_dir> <run_id> <setup_git_commit>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
setup_git_commit <- args[[3]]

required <- c(
  "glmnet", "Matrix", "matrixStats", "pROC", "ggplot2", "digest"
)
missing_packages <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

result_dir <- file.path(project_dir, "results", "09B_model_training", run_id)
figure_dir <- file.path(project_dir, "figures", "09B_model_training", run_id)
object_dir <- file.path(project_dir, "objects")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

model_path <- file.path(object_dir, "locked_stool_model.rds")
if (file.exists(model_path)) {
  stop("Refusing to overwrite existing locked model: ", model_path)
}
temporary_model_path <- file.path(
  object_dir, paste0(".locked_stool_model.", run_id, ".tmp.rds")
)
if (file.exists(temporary_model_path)) {
  stop("Temporary model path already exists: ", temporary_model_path)
}

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  folds <- integer(length(y))
  for (class_value in sort(unique(y))) {
    idx <- which(y == class_value)
    idx <- sample(idx, length(idx), replace = FALSE)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  if (any(folds == 0L)) stop("Failed to assign all fold IDs")
  folds
}

clamp_probability <- function(p) {
  pmin(pmax(as.numeric(p), 1e-6), 1 - 1e-6)
}

auc_rank <- function(y, p) {
  y <- as.integer(y)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
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
    suppressWarnings(
      stats::glm(y ~ lp, family = stats::binomial())
    ),
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

performance_values <- function(y, p, threshold = 0.5) {
  y <- as.integer(y)
  p <- clamp_probability(p)
  predicted <- as.integer(p >= threshold)
  tp <- sum(predicted == 1L & y == 1L)
  tn <- sum(predicted == 0L & y == 0L)
  fp <- sum(predicted == 1L & y == 0L)
  fn <- sum(predicted == 0L & y == 1L)
  calibration <- calibration_values(y, p)
  c(
    n = length(y),
    prevalence = mean(y),
    AUC = auc_rank(y, p),
    sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
    specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
    PPV = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    NPV = if ((tn + fn) > 0) tn / (tn + fn) else NA_real_,
    accuracy = mean(predicted == y),
    Brier = mean((p - y)^2),
    calibration
  )
}

bootstrap_performance <- function(y, p, threshold, n_boot, seed) {
  set.seed(seed)
  point <- performance_values(y, p, threshold)
  metrics <- names(point)
  boot <- matrix(
    NA_real_, nrow = n_boot, ncol = length(metrics),
    dimnames = list(NULL, metrics)
  )
  class_indices <- split(seq_along(y), y)
  for (b in seq_len(n_boot)) {
    idx <- unlist(
      lapply(
        class_indices,
        function(z) sample(z, length(z), replace = TRUE)
      ),
      use.names = FALSE
    )
    boot[b, ] <- performance_values(y[idx], p[idx], threshold)
  }
  do.call(
    rbind,
    lapply(metrics, function(metric) {
      values <- boot[, metric]
      data.frame(
        metric = metric,
        estimate = unname(point[[metric]]),
        ci_low = as.numeric(
          stats::quantile(values, 0.025, na.rm = TRUE, names = FALSE)
        ),
        ci_high = as.numeric(
          stats::quantile(values, 0.975, na.rm = TRUE, names = FALSE)
        ),
        bootstrap_replicates = n_boot,
        stringsAsFactors = FALSE
      )
    })
  )
}

stage9a_object_path <- file.path(
  project_dir, "objects",
  "GSE99573_9A_training_RMA_20260729_213907.rds"
)
stage9a <- readRDS(stage9a_object_path)
if (!identical(stage9a$test_expression_accessed, FALSE)) {
  stop("Stage 9A object does not prove test-expression isolation")
}
expression <- stage9a$expression
if (!identical(dim(expression), c(70523L, 265L))) {
  stop("Unexpected Stage 9A training matrix dimensions")
}
if (any(!is.finite(expression))) stop("Training expression contains nonfinite values")

manifest <- utils::read.delim(
  file.path(project_dir, "metadata", "dataset_manifest.tsv"),
  check.names = FALSE
)
manifest <- manifest[manifest$accession == "GSE99573", , drop = FALSE]
training_meta <- manifest[
  manifest$validation_split == "training",
  c("sample_id", "donor_id", "condition"),
  drop = FALSE
]
test_ids <- manifest$sample_id[manifest$validation_split == "testing"]
not_used_ids <- manifest$sample_id[manifest$validation_split == "not_used"]
if (nrow(training_meta) != 265L ||
    anyDuplicated(training_meta$sample_id) ||
    anyDuplicated(training_meta$donor_id)) {
  stop("Training metadata does not contain 265 independent participants")
}
if (!setequal(colnames(expression), training_meta$sample_id)) {
  stop("Training expression columns do not match the frozen training split")
}
if (any(colnames(expression) %in% c(test_ids, not_used_ids))) {
  stop("Forbidden test or Not Used sample detected in training expression")
}
training_meta <- training_meta[
  match(colnames(expression), training_meta$sample_id),
  ,
  drop = FALSE
]

mapping <- utils::read.delim(
  file.path(
    project_dir, "results_final",
    "stage_9A_locked_gene_probe_mapping.tsv"
  ),
  check.names = FALSE
)
eligible <- mapping[mapping$detectability_status == "detectable", , drop = FALSE]
eligible <- eligible[order(eligible$gene), , drop = FALSE]
if (nrow(eligible) != 521L || anyDuplicated(eligible$gene)) {
  stop("Frozen detectable feature universe must contain 521 unique genes")
}

gene_expression <- matrix(
  NA_real_,
  nrow = nrow(eligible),
  ncol = ncol(expression),
  dimnames = list(eligible$gene, colnames(expression))
)
for (i in seq_len(nrow(eligible))) {
  probe_ids <- strsplit(
    eligible$transcript_cluster_ids[[i]], ";", fixed = TRUE
  )[[1L]]
  probe_ids <- intersect(probe_ids, rownames(expression))
  if (!length(probe_ids)) {
    stop("No transcript cluster available for eligible gene ", eligible$gene[[i]])
  }
  values <- expression[probe_ids, , drop = FALSE]
  gene_expression[i, ] <- if (nrow(values) == 1L) {
    as.numeric(values[1L, ])
  } else {
    matrixStats::colMedians(values, na.rm = TRUE)
  }
}
if (any(!is.finite(gene_expression))) {
  stop("Frozen gene-expression matrix contains nonfinite values")
}

qc <- utils::read.delim(
  file.path(
    project_dir, "results_final",
    "stage_9A_training_sample_qc_metrics.tsv"
  ),
  check.names = FALSE
)
if (nrow(qc) != 265L || sum(qc$qc_review_flag) != 9L) {
  stop("Stage 9A QC inventory no longer matches the accepted result")
}

membership <- utils::read.delim(
  file.path(
    project_dir, "results_final",
    "stage_6A_stage_blind_module_membership.tsv"
  ),
  check.names = FALSE
)
candidate_modules <- utils::read.delim(
  file.path(
    project_dir, "results_final",
    "stage_6A_exploratory_candidate_modules.tsv"
  ),
  check.names = FALSE
)
membership <- membership[
  membership$module_id %in% candidate_modules$module_id &
    membership$gene %in% eligible$gene,
  c("module_id", "epithelial_state", "gene"),
  drop = FALSE
]
gene_modules <- stats::aggregate(
  module_id ~ gene,
  data = membership,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)
eligible$module_id <- gene_modules$module_id[
  match(eligible$gene, gene_modules$gene)
]

endpoint_definitions <- list(
  adenoma_vs_normal = list(
    role = "primary",
    included = c("normal", "adenoma"),
    case = "adenoma",
    label = "Adenoma versus normal"
  ),
  cancer_vs_normal = list(
    role = "secondary",
    included = c("normal", "cancer"),
    case = "cancer",
    label = "Colorectal cancer versus normal"
  ),
  neoplasia_vs_normal = list(
    role = "secondary_screening_composite",
    included = c("normal", "adenoma", "cancer"),
    case = c("adenoma", "cancer"),
    label = "Adenoma or colorectal cancer versus normal"
  )
)

all_predictions <- list()
all_metrics <- list()
all_fold_assignments <- list()
all_coefficients <- list()
all_outer_models <- list()
locked_models <- list()
threshold <- 0.5
outer_k <- 5L
inner_k <- 5L
final_k <- 10L
n_boot <- 2000L

for (endpoint_index in seq_along(endpoint_definitions)) {
  endpoint <- names(endpoint_definitions)[[endpoint_index]]
  definition <- endpoint_definitions[[endpoint]]
  keep <- training_meta$condition %in% definition$included
  meta <- training_meta[keep, , drop = FALSE]
  y <- as.integer(meta$condition %in% definition$case)
  x <- t(gene_expression[, meta$sample_id, drop = FALSE])
  storage.mode(x) <- "double"
  if (length(unique(y)) != 2L) stop("Endpoint is not binary: ", endpoint)
  expected_n <- c(
    adenoma_vs_normal = 171L,
    cancer_vs_normal = 183L,
    neoplasia_vs_normal = 265L
  )[[endpoint]]
  if (length(y) != expected_n) stop("Unexpected endpoint sample count: ", endpoint)

  outer_seed <- 42000L + endpoint_index * 100L
  outer_fold <- make_stratified_folds(y, outer_k, outer_seed)
  oof <- rep(NA_real_, length(y))
  outer_rows <- list()
  for (outer in seq_len(outer_k)) {
    train_idx <- which(outer_fold != outer)
    holdout_idx <- which(outer_fold == outer)
    inner_fold <- make_stratified_folds(
      y[train_idx], inner_k, outer_seed + outer
    )
    fit <- glmnet::cv.glmnet(
      x = x[train_idx, , drop = FALSE],
      y = y[train_idx],
      family = "binomial",
      alpha = 1,
      foldid = inner_fold,
      type.measure = "deviance",
      nlambda = 100,
      standardize = TRUE,
      intercept = TRUE,
      parallel = FALSE,
      keep = FALSE
    )
    oof[holdout_idx] <- as.numeric(
      stats::predict(
        fit,
        newx = x[holdout_idx, , drop = FALSE],
        s = "lambda.1se",
        type = "response"
      )
    )
    coefficients <- as.matrix(stats::coef(fit, s = "lambda.1se"))
    outer_rows[[outer]] <- data.frame(
      endpoint = endpoint,
      outer_fold = outer,
      training_n = length(train_idx),
      holdout_n = length(holdout_idx),
      lambda_min = fit$lambda.min,
      lambda_1se = fit$lambda.1se,
      nonzero_features = sum(coefficients[-1L, 1L] != 0),
      stringsAsFactors = FALSE
    )
  }
  if (any(!is.finite(oof))) stop("Nonfinite outer prediction: ", endpoint)

  predictions <- data.frame(
    endpoint = endpoint,
    endpoint_role = definition$role,
    sample_id = meta$sample_id,
    donor_id = meta$donor_id,
    outcome = y,
    outer_fold = outer_fold,
    predicted_probability = oof,
    threshold = threshold,
    predicted_class = as.integer(oof >= threshold),
    data_scope = "training_outer_out_of_fold_only",
    stringsAsFactors = FALSE
  )
  metrics <- bootstrap_performance(
    y, oof, threshold, n_boot, 420000L + endpoint_index
  )
  metrics$endpoint <- endpoint
  metrics$endpoint_role <- definition$role
  metrics$evaluation_scope <- "nested_CV_training_only"
  metrics <- metrics[
    ,
    c(
      "endpoint", "endpoint_role", "metric", "estimate", "ci_low",
      "ci_high", "bootstrap_replicates", "evaluation_scope"
    )
  ]

  final_fold <- make_stratified_folds(
    y, final_k, 421000L + endpoint_index
  )
  final_cv <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 1,
    foldid = final_fold,
    type.measure = "deviance",
    nlambda = 100,
    standardize = TRUE,
    intercept = TRUE,
    parallel = FALSE,
    keep = FALSE
  )
  coefficient_matrix <- as.matrix(
    stats::coef(final_cv, s = "lambda.1se")
  )
  coefficient_vector <- as.numeric(coefficient_matrix[, 1L])
  names(coefficient_vector) <- rownames(coefficient_matrix)
  intercept <- unname(coefficient_vector[["(Intercept)"]])
  beta <- coefficient_vector[names(coefficient_vector) != "(Intercept)"]
  selected <- names(beta)[beta != 0]

  coefficient_table <- data.frame(
    endpoint = endpoint,
    endpoint_role = definition$role,
    feature = names(beta),
    coefficient = as.numeric(beta),
    selected = as.numeric(beta) != 0,
    module_id = eligible$module_id[match(names(beta), eligible$gene)],
    lambda_rule = "lambda.1se",
    lambda = final_cv$lambda.1se,
    threshold = threshold,
    stringsAsFactors = FALSE
  )
  intercept_row <- data.frame(
    endpoint = endpoint,
    endpoint_role = definition$role,
    feature = "(Intercept)",
    coefficient = intercept,
    selected = TRUE,
    module_id = "NA",
    lambda_rule = "lambda.1se",
    lambda = final_cv$lambda.1se,
    threshold = threshold,
    stringsAsFactors = FALSE
  )
  coefficient_table <- rbind(intercept_row, coefficient_table)

  locked_models[[endpoint]] <- list(
    endpoint = endpoint,
    endpoint_label = definition$label,
    endpoint_role = definition$role,
    case_conditions = definition$case,
    control_condition = "normal",
    training_n = length(y),
    training_cases = sum(y == 1L),
    training_controls = sum(y == 0L),
    training_sample_ids = meta$sample_id,
    alpha = 1,
    lambda_rule = "lambda.1se",
    lambda = final_cv$lambda.1se,
    lambda_min = final_cv$lambda.min,
    threshold = threshold,
    intercept = intercept,
    coefficients = beta,
    selected_features = selected,
    selected_feature_count = length(selected),
    formula_text = paste0(
      "logit(P(case)) = ", format(intercept, digits = 17),
      if (length(selected)) {
        paste0(
          " + sum(beta_j * log2_RMA_gene_j) for ",
          length(selected), " nonzero locked features"
        )
      } else {
        " (intercept-only)"
      }
    ),
    standardization = paste(
      "glmnet internal training standardization; stored coefficients are",
      "returned on the original gene-level log2 RMA scale"
    ),
    final_fold_id = stats::setNames(final_fold, meta$sample_id),
    nested_cv_metrics = metrics
  )

  all_predictions[[endpoint]] <- predictions
  all_metrics[[endpoint]] <- metrics
  all_fold_assignments[[endpoint]] <- data.frame(
    endpoint = endpoint,
    sample_id = meta$sample_id,
    outcome = y,
    outer_fold = outer_fold,
    final_10fold_id = final_fold,
    stringsAsFactors = FALSE
  )
  all_coefficients[[endpoint]] <- coefficient_table
  all_outer_models[[endpoint]] <- do.call(rbind, outer_rows)
}

predictions <- do.call(rbind, all_predictions)
metrics <- do.call(rbind, all_metrics)
fold_assignments <- do.call(rbind, all_fold_assignments)
coefficient_table <- do.call(rbind, all_coefficients)
outer_model_summary <- do.call(rbind, all_outer_models)

write_tsv(
  predictions,
  file.path(result_dir, "nested_cv_predictions_training_only.tsv")
)
write_tsv(
  metrics,
  file.path(result_dir, "nested_cv_performance.tsv")
)
write_tsv(
  fold_assignments,
  file.path(result_dir, "fold_assignments.tsv")
)
write_tsv(
  coefficient_table,
  file.path(result_dir, "locked_model_coefficients.tsv")
)
write_tsv(
  outer_model_summary,
  file.path(result_dir, "outer_fold_model_summary.tsv")
)
write_tsv(
  eligible[
    ,
    c(
      "gene", "transcript_cluster_ids", "module_id",
      "training_median_log2_RMA", "training_fraction_above_floor",
      "detectability_status"
    )
  ],
  file.path(result_dir, "locked_feature_manifest.tsv")
)

endpoint_lock_summary <- do.call(
  rbind,
  lapply(locked_models, function(model) {
    data.frame(
      endpoint = model$endpoint,
      endpoint_role = model$endpoint_role,
      training_n = model$training_n,
      training_cases = model$training_cases,
      training_controls = model$training_controls,
      eligible_features = length(model$coefficients),
      selected_features = model$selected_feature_count,
      alpha = model$alpha,
      lambda_rule = model$lambda_rule,
      lambda = model$lambda,
      threshold = model$threshold,
      stringsAsFactors = FALSE
    )
  })
)
write_tsv(
  endpoint_lock_summary,
  file.path(result_dir, "model_lock_summary.tsv")
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

locked_artifact <- list(
  artifact = "CRC_carcinogenesis_locked_stool_model",
  artifact_version = "1.0.0",
  stage = "9B",
  run_id = run_id,
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  setup_git_commit = setup_git_commit,
  random_seed = 42L,
  training_accession = "GSE99573",
  training_split = "GEO set=Training",
  training_sample_count = 265L,
  test_sample_count = 65L,
  test_expression_accessed = FALSE,
  qc_policy = paste(
    "All 265 training arrays retained; nine Stage 9A review flags recorded;",
    "zero outcome-driven exclusions"
  ),
  feature_universe = eligible[
    ,
    c("gene", "transcript_cluster_ids", "module_id"),
    drop = FALSE
  ],
  feature_generation = list(
    platform = "GPL17586",
    normalization = paste(
      "oligo::rma target=core, background=TRUE, normalize=TRUE;",
      "training batch only"
    ),
    collapse = "median across all frozen mapped core transcript clusters",
    expression_scale = "log2 RMA",
    missing_selected_feature_action = "abort; do not substitute or refit"
  ),
  algorithm = list(
    family = "binomial",
    method = "LASSO logistic regression",
    package = "glmnet",
    alpha = 1,
    lambda_rule = "lambda.1se",
    lambda_metric = "binomial deviance",
    nlambda = 100L,
    threshold = 0.5,
    outer_folds = 5L,
    inner_folds = 5L,
    final_folds = 10L,
    bootstrap_replicates = 2000L
  ),
  endpoint_hierarchy = data.frame(
    endpoint = names(endpoint_definitions),
    role = vapply(endpoint_definitions, `[[`, character(1), "role"),
    label = vapply(endpoint_definitions, `[[`, character(1), "label"),
    stringsAsFactors = FALSE
  ),
  models = locked_models,
  training_object_sha256 = digest::digest(
    stage9a_object_path, algo = "sha256", file = TRUE
  ),
  locked_feature_manifest_sha256 = digest::digest(
    file.path(
      project_dir, "results_final",
      "stage_9A_locked_gene_probe_mapping.tsv"
    ),
    algo = "sha256",
    file = TRUE
  ),
  software_versions = software
)

saveRDS(locked_artifact, temporary_model_path, compress = "xz")
if (!file.rename(temporary_model_path, model_path)) {
  stop("Atomic move to canonical locked model path failed")
}
Sys.chmod(model_path, mode = "0444", use_umask = FALSE)

model_sha256 <- digest::digest(model_path, algo = "sha256", file = TRUE)
writeLines(
  paste(model_sha256, model_path),
  file.path(result_dir, "locked_stool_model.sha256")
)
write_tsv(
  data.frame(
    field = c(
      "model_path", "model_sha256", "run_id", "setup_git_commit",
      "test_expression_accessed", "random_seed", "algorithm",
      "lambda_rule", "threshold", "endpoint_count"
    ),
    value = c(
      model_path, model_sha256, run_id, setup_git_commit, "FALSE", "42",
      "LASSO logistic regression", "lambda.1se", "0.5", "3"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "model_lock_manifest.tsv")
)

plot_data <- predictions
plot_data$endpoint <- factor(
  plot_data$endpoint,
  levels = names(endpoint_definitions),
  labels = vapply(endpoint_definitions, `[[`, character(1), "label")
)
roc_plot <- ggplot2::ggplot()
roc_colors <- c("#0072B2", "#D55E00", "#009E73")
for (i in seq_along(levels(plot_data$endpoint))) {
  endpoint_label <- levels(plot_data$endpoint)[[i]]
  z <- plot_data[plot_data$endpoint == endpoint_label, , drop = FALSE]
  roc <- pROC::roc(
    response = z$outcome,
    predictor = z$predicted_probability,
    direction = "<",
    quiet = TRUE
  )
  coords <- pROC::coords(
    roc, x = "all", ret = c("specificity", "sensitivity"),
    transpose = FALSE
  )
  curve <- data.frame(
    false_positive_rate = 1 - coords$specificity,
    true_positive_rate = coords$sensitivity,
    endpoint = endpoint_label
  )
  roc_plot <- roc_plot + ggplot2::geom_path(
    data = curve,
    ggplot2::aes(
      x = false_positive_rate,
      y = true_positive_rate,
      color = endpoint
    ),
    linewidth = 0.8
  )
}
roc_plot <- roc_plot +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", color = "grey55"
  ) +
  ggplot2::scale_color_manual(values = roc_colors) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Nested cross-validation ROC curves",
    subtitle = "GSE99573 training partition only; not independent test performance",
    x = "False-positive rate",
    y = "True-positive rate",
    color = "Endpoint"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
  filename = file.path(figure_dir, "nested_cv_ROC_training_only.pdf"),
  plot = roc_plot, width = 7, height = 6, units = "in"
)
ggplot2::ggsave(
  filename = file.path(figure_dir, "nested_cv_ROC_training_only.png"),
  plot = roc_plot, width = 7, height = 6, units = "in", dpi = 300
)

calibration_rows <- do.call(
  rbind,
  lapply(split(predictions, predictions$endpoint), function(z) {
    breaks <- unique(
      as.numeric(
        stats::quantile(
          z$predicted_probability,
          probs = seq(0, 1, length.out = 11),
          na.rm = TRUE,
          type = 8
        )
      )
    )
    if (length(breaks) < 3L) {
      breaks <- seq(
        min(z$predicted_probability),
        max(z$predicted_probability),
        length.out = 3L
      )
    }
    bins <- cut(
      z$predicted_probability,
      breaks = breaks,
      include.lowest = TRUE
    )
    out <- stats::aggregate(
      cbind(predicted_probability, outcome) ~ bins,
      data = z,
      FUN = mean
    )
    n <- as.data.frame(table(bins), stringsAsFactors = FALSE)
    out$n <- n$Freq[match(as.character(out$bins), as.character(n$bins))]
    out$endpoint <- z$endpoint[[1L]]
    out
  })
)
write_tsv(
  calibration_rows,
  file.path(result_dir, "calibration_plot_source_data.tsv")
)
calibration_rows$endpoint <- factor(
  calibration_rows$endpoint,
  levels = names(endpoint_definitions),
  labels = vapply(endpoint_definitions, `[[`, character(1), "label")
)
calibration_plot <- ggplot2::ggplot(
  calibration_rows,
  ggplot2::aes(
    x = predicted_probability,
    y = outcome,
    color = endpoint
  )
) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", color = "grey55"
  ) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.85) +
  ggplot2::scale_color_manual(values = roc_colors) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Nested cross-validation calibration",
    subtitle = "Decile summaries from training-only out-of-fold predictions",
    x = "Mean predicted probability",
    y = "Observed case fraction",
    color = "Endpoint",
    size = "Participants"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(
  filename = file.path(figure_dir, "nested_cv_calibration_training_only.pdf"),
  plot = calibration_plot, width = 7, height = 6, units = "in"
)
ggplot2::ggsave(
  filename = file.path(figure_dir, "nested_cv_calibration_training_only.png"),
  plot = calibration_plot, width = 7, height = 6, units = "in", dpi = 300
)

cat("STAGE9B_MODEL_TRAINING_OK run_id=", run_id, "\n", sep = "")
