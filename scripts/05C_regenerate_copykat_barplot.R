#!/usr/bin/env Rscript

# Rebuild the CopyKAT tissue-composition figure from its source-data table.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(args[[1]], mustWork = TRUE)
source_file <- file.path(
  project_root, "results", "05C_annotation", "source_data",
  "copykat_direct_prediction_composition.tsv"
)
figure_dir <- file.path(project_root, "figures", "05C_annotation")

source <- fread(source_file)
sample_order <- source[
  order(
    factor(
      lesion_stage,
      levels = c("normal", "adenoma_polyp", "cancer")
    ),
    biological_sample_id
  ),
  unique(biological_sample_id)
]
source[
  ,
  biological_sample_plot := factor(
    biological_sample_id, levels = sample_order
  )
]

plot <- ggplot(
  source,
  aes(
    x = biological_sample_plot,
    y = fraction,
    fill = copykat_prediction
  )
) +
  geom_col(width = 0.85) +
  facet_grid(. ~ lesion_stage, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      aneuploid = "#B2182B",
      diploid = "#2166AC",
      not.defined = "#BDBDBD"
    )
  ) +
  labs(
    title = "Direct CopyKAT predictions by biological tissue",
    x = "Biological tissue sample",
    y = "Sampled epithelial-cell fraction",
    fill = "CopyKAT"
  ) +
  theme_bw(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(figure_dir, "copykat_direct_predictions_by_sample.pdf"),
  plot, width = 14, height = 6
)
ggsave(
  file.path(figure_dir, "copykat_direct_predictions_by_sample.png"),
  plot, width = 14, height = 6, dpi = 300
)
message("CopyKAT direct-prediction figure regenerated")
