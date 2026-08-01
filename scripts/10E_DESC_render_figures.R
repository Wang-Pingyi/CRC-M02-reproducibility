#!/usr/bin/env Rscript

# Post-run rendering-only fix: reads frozen source-data tables, never expression.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 10E_DESC_render_figures.R <root>")
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
project_lib <- file.path(root, "environment", "R-library")
if (dir.exists(project_lib)) .libPaths(unique(c(project_lib, .libPaths())))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

fig_dir <- file.path(root, "figures", "stage10e_desc")
source_dir <- file.path(fig_dir, "source_data")
map_source <- fread(file.path(source_dir, "Fig10E_DESC_1_case4_spatial_source_data.tsv"))
paired_source <- fread(file.path(source_dir, "Fig10E_DESC_2_case4_paired_source_data.tsv"))
if (!identical(unique(map_source$patient_id), "case4") || !identical(unique(paired_source$patient_id), "case4")) {
  stop("Rendering source is not case4-only")
}

white_background <- theme(
  plot.background = element_rect(fill = "white", colour = NA),
  panel.background = element_rect(fill = "white", colour = NA),
  legend.background = element_rect(fill = "white", colour = NA),
  legend.key = element_rect(fill = "white", colour = NA)
)

p1 <- ggplot(map_source, aes(x = pxl_col, y = -pxl_row, colour = spot_m02_score)) +
  geom_point(size = 0.55) +
  scale_colour_viridis_c(limits = c(0, 12), oob = scales::squish, name = "30-gene score") +
  coord_equal() + theme_void(base_size = 9) + white_background +
  labs(title = "case4 frozen ROI: M02 spatial localization",
       subtitle = "n=1 patient; 30/36 mapped genes; descriptive only",
       caption = "Unsmoothed spot display; spots are not biological replicates")

p2 <- ggplot(paired_source, aes(x = factor(region, levels = c("Normal", "Adenoma")), y = score, group = patient_id)) +
  geom_line(linewidth = 0.5, colour = "#666666") + geom_point(size = 2.6, colour = "#0072B2") +
  theme_classic(base_size = 9) + white_background + xlab(NULL) + ylab("Equal-weight TMM log2CPM score") +
  labs(title = "case4 Normal and Adenoma ROI scores",
       subtitle = "n=1 patient; 30/36 mapped genes; descriptive only",
       caption = "No P value, confidence interval or population inference")

for (item in list(list(p = p1, name = "Fig10E_DESC_1_case4_spatial"),
                  list(p = p2, name = "Fig10E_DESC_2_case4_paired"))) {
  ggsave(file.path(fig_dir, paste0(item$name, ".png")), item$p, width = 7, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(item$name, ".pdf")), item$p, width = 7, height = 5, device = cairo_pdf, bg = "white")
}

writeLines(c(
  "stage=STAGE10E_DESC", "action=render_only_white_background_fix",
  "expression_or_H5_read=NO", "source_data_changed=NO", "analysis_values_changed=NO"
), file.path(root, "results", "stage10e_desc", "STAGE10E_DESC_RENDER_FIX.txt"))
