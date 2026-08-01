#!/usr/bin/env Rscript

# Figure-only renderer for accepted Stage 9C source-data tables.
# This script does not read CEL/RDS files, generate predictions, or fit models.

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09C_render_figures.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09C_external_test", run_id)
figure_dir <- file.path(project_dir, "figures", "09C_external_test", run_id)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Missing package: ggplot2")
}
read_tsv <- function(name) {
  utils::read.delim(
    file.path(result_dir, name),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
roc_source <- read_tsv("ROC_source_data.tsv")
calibration_source <- read_tsv("calibration_source_data.tsv")
endpoint_labels <- c(
  adenoma_vs_normal = "Adenoma vs normal",
  cancer_vs_normal = "Cancer vs normal",
  neoplasia_vs_normal = "Adenoma/cancer vs normal"
)
colors <- c("#0072B2", "#D55E00", "#009E73")

roc_source$endpoint_label <- factor(
  roc_source$endpoint,
  levels = names(endpoint_labels),
  labels = endpoint_labels
)
roc_plot <- ggplot2::ggplot(
  roc_source,
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
  roc_plot, width = 7, height = 6.3, units = "in"
)
ggplot2::ggsave(
  file.path(figure_dir, "stool_test_ROC.png"),
  roc_plot, width = 7, height = 6.3, units = "in", dpi = 300
)

calibration_source$endpoint_label <- factor(
  calibration_source$endpoint,
  levels = names(endpoint_labels),
  labels = endpoint_labels
)
calibration_plot <- ggplot2::ggplot(
  calibration_source,
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
  calibration_plot, width = 7, height = 6.3, units = "in"
)
ggplot2::ggsave(
  file.path(figure_dir, "stool_test_calibration.png"),
  calibration_plot, width = 7, height = 6.3, units = "in", dpi = 300
)
cat("STAGE9C_FIGURES_RENDERED_FROM_SOURCE_DATA_ONLY\n")
