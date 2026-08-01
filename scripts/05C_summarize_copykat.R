#!/usr/bin/env Rscript

# Analysis: Stage 5C summarize CopyKAT, assign conservative CNV-supported
# malignancy evidence, and write final annotated objects.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}

param_file <- file.path(project_root, "config", "annotation_parameters.tsv")
plan_file <- file.path(
  project_root, "results", "05C_annotation", "copykat_sample_plan.tsv"
)
epithelial_file <- file.path(
  project_root, "objects", "GSE201348_5C_epithelial_annotated_preCNV.rds"
)
major_file <- file.path(
  project_root, "objects", "GSE201348_5C_annotated_preCNV.rds"
)
copykat_dir <- file.path(project_root, "cache", "05C_copykat_outputs")
result_dir <- file.path(project_root, "results", "05C_annotation")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_root, "figures", "05C_annotation")
object_dir <- file.path(project_root, "objects")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

params <- fread(param_file)
get_param <- function(section_name, parameter_name, numeric = FALSE) {
  value <- params[
    section == section_name & parameter == parameter_name,
    value
  ]
  if (length(value) != 1L) {
    stop("Expected one parameter: ", section_name, "/", parameter_name)
  }
  if (numeric) as.numeric(value) else value
}
min_group <- as.integer(get_param(
  "cnv", "min_cells_per_donor_cluster", TRUE
))
aneuploid_threshold <- get_param(
  "cnv", "aneuploid_fraction_candidate", TRUE
)
diploid_threshold <- get_param(
  "cnv", "diploid_fraction_reference", TRUE
)

plan <- fread(plan_file, colClasses = "character")
numeric_plan <- c(
  "query_cells_available", "query_cells_selected",
  "same_sample_reference_available", "same_sample_reference_selected",
  "healthy_reference_supplement", "total_reference_selected",
  "input_cells"
)
plan[, (numeric_plan) := lapply(.SD, as.integer), .SDcols = numeric_plan]

status_rows <- lapply(plan$sample_key, function(sample_key) {
  planned_status <- plan$status[
    match(sample_key, plan$sample_key)
  ]
  status_path <- file.path(
    copykat_dir, sample_key, "run_status.tsv"
  )
  if (!file.exists(status_path)) {
    return(data.table(
      sample_key = sample_key,
      run_status = ifelse(
        planned_status == "ready",
        "missing",
        planned_status
      ),
      run_message = "No CopyKAT run-status file",
      elapsed_minutes = NA_real_,
      predictions = NA_integer_
    ))
  }
  x <- fread(status_path)
  data.table(
    sample_key = sample_key,
    run_status = x$status[1],
    run_message = x$message[1],
    elapsed_minutes = x$elapsed_minutes[1],
    predictions = x$predictions[1]
  )
})
run_audit <- merge(
  plan,
  rbindlist(status_rows, fill = TRUE),
  by = "sample_key",
  all.x = TRUE
)
fwrite(
  run_audit,
  file.path(result_dir, "copykat_run_audit.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

completed_keys <- run_audit[run_status == "completed", sample_key]
prediction_list <- lapply(completed_keys, function(sample_key) {
  path <- file.path(
    copykat_dir, sample_key, "copykat_predictions.tsv.gz"
  )
  if (!file.exists(path)) {
    stop("Completed run lacks prediction file: ", sample_key)
  }
  fread(path)
})
predictions <- if (length(prediction_list)) {
  rbindlist(prediction_list, fill = TRUE)
} else {
  data.table()
}
if (!nrow(predictions)) stop("No completed CopyKAT predictions")
query_predictions <- predictions[
  copykat_role == "epithelial_query"
]
query_predictions[
  ,
  epithelial_cluster := as.character(epithelial_cluster)
]
if (anyDuplicated(query_predictions$cell_id)) {
  duplicated_cells <- query_predictions[
    duplicated(cell_id) | duplicated(cell_id, fromLast = TRUE),
    unique(cell_id)
  ]
  stop(
    "Epithelial query cells appeared in multiple CopyKAT runs: ",
    paste(head(duplicated_cells, 10L), collapse = ";")
  )
}
query_predictions[
  ,
  defined_prediction := copykat_prediction %chin% c(
    "aneuploid", "diploid"
  )
]
group_summary <- query_predictions[
  ,
  .(
    sampled_query_cells = .N,
    defined_cells = sum(defined_prediction),
    aneuploid_cells = sum(copykat_prediction == "aneuploid"),
    diploid_cells = sum(copykat_prediction == "diploid"),
    not_defined_cells = sum(!defined_prediction)
  ),
  by = .(biological_sample_id, donor_id, epithelial_cluster)
]
group_summary[
  ,
  `:=`(
    aneuploid_fraction_defined = fifelse(
      defined_cells > 0, aneuploid_cells / defined_cells, NA_real_
    ),
    diploid_fraction_defined = fifelse(
      defined_cells > 0, diploid_cells / defined_cells, NA_real_
    )
  )
]
group_summary[
  ,
  cnv_cluster_consensus := fcase(
    defined_cells < min_group,
    "uncertain_insufficient_defined_cells",
    aneuploid_fraction_defined >= aneuploid_threshold,
    "candidate_malignant_CNV_cluster_supported",
    diploid_fraction_defined >= diploid_threshold,
    "likely_non_malignant_CNV_cluster_supported",
    default = "uncertain_mixed_CNV"
  )
]
fwrite(
  group_summary,
  file.path(result_dir, "copykat_sample_cluster_summary.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)
fwrite(
  query_predictions,
  file.path(result_dir, "copykat_epithelial_predictions_sampled.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

epithelial <- readRDS(epithelial_file)
major <- readRDS(major_file)
epi_row_names <- rownames(epithelial@meta.data)
epi_meta <- as.data.table(epithelial@meta.data)
if ("cell_id" %in% colnames(epi_meta)) {
  if (!identical(as.character(epi_meta$cell_id), epi_row_names)) {
    stop("Existing epithelial cell_id column does not match metadata rows")
  }
} else {
  epi_meta[, cell_id := epi_row_names]
}
if (anyDuplicated(colnames(epi_meta))) {
  stop("Epithelial metadata contains duplicated column names")
}
epi_meta[
  ,
  epithelial_cluster := as.character(epithelial_cluster)
]
direct <- query_predictions[
  ,
  .(cell_id, copykat_prediction_direct = copykat_prediction)
]
epi_calls <- merge(epi_meta, direct, by = "cell_id", all.x = TRUE)
epi_calls <- merge(
  epi_calls,
  group_summary[
    ,
    .(
      biological_sample_id,
      epithelial_cluster,
      sampled_query_cells,
      defined_cells,
      aneuploid_fraction_defined,
      diploid_fraction_defined,
      cnv_cluster_consensus
    )
  ],
  by = c("biological_sample_id", "epithelial_cluster"),
  all.x = TRUE
)
epi_calls[
  ,
  malignancy_evidence := fcase(
    copykat_prediction_direct == "aneuploid",
    "candidate_malignant_CNV_direct",
    copykat_prediction_direct == "diploid",
    "likely_non_malignant_CNV_direct",
    cnv_cluster_consensus ==
      "candidate_malignant_CNV_cluster_supported",
    "candidate_malignant_CNV_cluster_supported",
    cnv_cluster_consensus ==
      "likely_non_malignant_CNV_cluster_supported",
    "likely_non_malignant_CNV_cluster_supported",
    is.na(cnv_cluster_consensus),
    "not_evaluated_CNV",
    default = "uncertain_CNV"
  )
]
epi_calls[
  ,
  malignancy_interpretation := fcase(
    grepl("^candidate_malignant", malignancy_evidence),
    "candidate_malignant_epithelial",
    grepl("^likely_non_malignant", malignancy_evidence),
    "likely_non_malignant_epithelial",
    malignancy_evidence == "not_evaluated_CNV",
    "not_evaluated",
    default = "uncertain"
  )
]
epi_calls[
  ,
  malignancy_call_basis := fifelse(
    !is.na(copykat_prediction_direct) &
      copykat_prediction_direct %chin% c("aneuploid", "diploid"),
    "direct_CopyKAT_cell_call",
    fifelse(
      !is.na(cnv_cluster_consensus),
      "same_tissue_epithelial_cluster_consensus",
      "no_CNV_result"
    )
  )
]
epi_calls[
  ,
  prohibited_origin_only_call := FALSE
]

call_map <- setNames(epi_calls$malignancy_evidence, epi_calls$cell_id)
interpretation_map <- setNames(
  epi_calls$malignancy_interpretation, epi_calls$cell_id
)
direct_map <- setNames(
  epi_calls$copykat_prediction_direct, epi_calls$cell_id
)
basis_map <- setNames(epi_calls$malignancy_call_basis, epi_calls$cell_id)
epithelial$copykat_prediction_direct <- unname(
  direct_map[colnames(epithelial)]
)
epithelial$malignancy_evidence <- unname(
  call_map[colnames(epithelial)]
)
epithelial$malignancy_interpretation <- unname(
  interpretation_map[colnames(epithelial)]
)
epithelial$malignancy_call_basis <- unname(
  basis_map[colnames(epithelial)]
)
if (anyNA(epithelial$malignancy_evidence)) {
  stop("Missing final malignancy evidence in epithelial object")
}

major$copykat_prediction_direct <- NA_character_
major$malignancy_evidence <- NA_character_
major$malignancy_interpretation <- NA_character_
major$malignancy_call_basis <- NA_character_
epi_positions <- match(colnames(epithelial), colnames(major))
if (anyNA(epi_positions)) stop("Epithelial cell absent from major object")
major$copykat_prediction_direct[epi_positions] <-
  epithelial$copykat_prediction_direct
major$malignancy_evidence[epi_positions] <-
  epithelial$malignancy_evidence
major$malignancy_interpretation[epi_positions] <-
  epithelial$malignancy_interpretation
major$malignancy_call_basis[epi_positions] <-
  epithelial$malignancy_call_basis

fwrite(
  epi_calls[
    ,
    .(
      cell_id,
      donor_id,
      sample_id,
      biological_sample_id,
      lesion_stage,
      epithelial_cluster,
      epithelial_state,
      copykat_prediction_direct,
      cnv_cluster_consensus,
      malignancy_evidence,
      malignancy_interpretation,
      malignancy_call_basis
    )
  ],
  file.path(result_dir, "epithelial_cell_annotations_final.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

malignancy_summary <- epi_calls[
  ,
  .N,
  by = .(
    donor_id,
    biological_sample_id,
    lesion_stage,
    epithelial_state,
    malignancy_interpretation
  )
]
malignancy_summary[
  ,
  fraction_within_sample_state := N / sum(N),
  by = .(biological_sample_id, epithelial_state)
]
fwrite(
  malignancy_summary,
  file.path(result_dir, "malignancy_annotation_summary.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

embedding <- Embeddings(epithelial, "epithelial_umap")
umap_source <- data.table(
  cell_id = rownames(embedding),
  UMAP_1 = embedding[, 1],
  UMAP_2 = embedding[, 2],
  epithelial_state = epithelial$epithelial_state,
  malignancy_interpretation = epithelial$malignancy_interpretation,
  malignancy_evidence = epithelial$malignancy_evidence,
  donor_id = epithelial$donor_id,
  biological_sample_id = epithelial$biological_sample_id,
  lesion_stage = epithelial$lesion_stage
)
fwrite(
  umap_source,
  file.path(source_dir, "epithelial_malignancy_umap_source.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

malignancy_palette <- c(
  candidate_malignant_epithelial = "#B2182B",
  likely_non_malignant_epithelial = "#2166AC",
  uncertain = "#FDB863",
  not_evaluated = "#BDBDBD"
)
p_malignancy <- ggplot(
  umap_source,
  aes(UMAP_1, UMAP_2, color = malignancy_interpretation)
) +
  geom_point(size = 0.08, alpha = 0.45) +
  scale_color_manual(
    values = malignancy_palette,
    drop = FALSE
  ) +
  labs(
    title = "CopyKAT-supported epithelial malignancy evidence",
    subtitle = "Aneuploidy is auxiliary evidence, not proof of malignancy",
    color = "CNV interpretation"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())
ggsave(
  file.path(figure_dir, "epithelial_copykat_malignancy.pdf"),
  p_malignancy, width = 8, height = 6
)
ggsave(
  file.path(figure_dir, "epithelial_copykat_malignancy.png"),
  p_malignancy, width = 8, height = 6, dpi = 300
)

direct_bar_source <- query_predictions[
  queried_epithelial == TRUE,
  .N,
  by = .(
    biological_sample_id,
    donor_id,
    lesion_stage,
    copykat_prediction
  )
]
direct_bar_source[
  ,
  fraction := N / sum(N),
  by = biological_sample_id
]
sample_order <- direct_bar_source[
  order(
    factor(
      lesion_stage,
      levels = c("normal", "adenoma_polyp", "cancer")
    ),
    biological_sample_id
  ),
  unique(biological_sample_id)
]
direct_bar_source[
  ,
  biological_sample_plot := factor(
    biological_sample_id, levels = sample_order
  )
]
fwrite(
  direct_bar_source,
  file.path(source_dir, "copykat_direct_prediction_composition.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)
p_bar <- ggplot(
  direct_bar_source,
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
  p_bar, width = 14, height = 6
)
ggsave(
  file.path(figure_dir, "copykat_direct_predictions_by_sample.png"),
  p_bar, width = 14, height = 6, dpi = 300
)

epithelial@misc$stage_5C_CNV <- list(
  method = "CopyKAT",
  unit = "biological_tissue_sample",
  raw_counts = TRUE,
  malignant_interpretation = "candidate_only",
  tumor_origin_alone_used = FALSE,
  EPCAM_alone_used = FALSE,
  direct_prediction_file =
    "results/05C_annotation/copykat_epithelial_predictions_sampled.tsv.gz",
  cluster_summary_file =
    "results/05C_annotation/copykat_sample_cluster_summary.tsv"
)
major@misc$stage_5C_CNV <- epithelial@misc$stage_5C_CNV
saveRDS(
  epithelial,
  file.path(object_dir, "GSE201348_5C_epithelial_annotated_CNV.rds"),
  compress = "gzip"
)
saveRDS(
  major,
  file.path(object_dir, "GSE201348_5C_annotated_final.rds"),
  compress = "gzip"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "copykat_summary_sessionInfo.txt")
)
writeLines(
  c(
    paste0("completed_copykat_tissues=", length(completed_keys)),
    paste0(
      "failed_or_missing_ready_tissues=",
      run_audit[
        status == "ready" & run_status != "completed",
        .N
      ]
    ),
    paste0("sampled_epithelial_cells=", nrow(query_predictions)),
    paste0(
      "direct_defined_epithelial_cells=",
      query_predictions[defined_prediction == TRUE, .N]
    ),
    paste0(
      "candidate_malignant_cells_final=",
      epi_calls[
        malignancy_interpretation == "candidate_malignant_epithelial",
        .N
      ]
    ),
    "malignant_interpretation=candidate_only",
    "tumor_tissue_origin_as_sole_call=no",
    "EPCAM_as_sole_call=no",
    "differential_expression=not_performed",
    "trajectory_analysis=not_performed",
    "machine_learning=not_performed"
  ),
  file.path(result_dir, "copykat_summary_provenance.txt")
)
message(
  "CopyKAT summary completed: ", length(completed_keys),
  " tissues; ", nrow(query_predictions), " sampled epithelial cells"
)
