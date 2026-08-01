#!/usr/bin/env Rscript

# Independent Stage 5C validation. Exits non-zero if any required check fails.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}

result_dir <- file.path(project_root, "results", "05C_annotation")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_root, "figures", "05C_annotation")
object_dir <- file.path(project_root, "objects")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail)
  )
}

required_small_files <- c(
  "config/annotation_parameters.tsv",
  "config/major_cluster_annotation.tsv",
  "config/epithelial_cluster_annotation.tsv",
  "metadata/annotation_evidence.tsv",
  "results/05C_annotation/major_cluster_markers.tsv.gz",
  "results/05C_annotation/major_signature_scores.tsv",
  "results/05C_annotation/major_cluster_composition.tsv",
  "results/05C_annotation/epithelial_cluster_markers.tsv.gz",
  "results/05C_annotation/epithelial_signature_scores.tsv",
  "results/05C_annotation/epithelial_cluster_composition.tsv",
  "results/05C_annotation/epithelial_donor_dominance_audit.tsv",
  "results/05C_annotation/copykat_run_audit.tsv",
  "results/05C_annotation/copykat_sample_cluster_summary.tsv",
  "results/05C_annotation/malignancy_annotation_summary.tsv",
  "results/05C_annotation/epithelial_cell_annotations_final.tsv.gz"
)
for (relative_path in required_small_files) {
  path <- file.path(project_root, relative_path)
  add_check(
    paste0("required_file:", relative_path),
    file.exists(path) && file.info(path)$size > 0,
    if (file.exists(path)) {
      paste0("bytes=", file.info(path)$size)
    } else {
      "missing"
    }
  )
}

figure_source_pairs <- data.table(
  figure = c(
    "major_discovery_clusters",
    "major_marker_dotplot",
    "major_cell_type_annotation",
    "epithelial_discovery_clusters",
    "epithelial_marker_dotplot",
    "epithelial_state_annotation",
    "epithelial_cluster_donor_composition",
    "epithelial_copykat_malignancy",
    "copykat_direct_predictions_by_sample"
  ),
  source = c(
    "major_umap_source.tsv.gz",
    "major_marker_dotplot_source.tsv.gz",
    "major_annotated_umap_source.tsv.gz",
    "epithelial_umap_source.tsv.gz",
    "epithelial_marker_dotplot_source.tsv.gz",
    "epithelial_annotated_umap_source.tsv.gz",
    "epithelial_cluster_donor_composition.tsv",
    "epithelial_malignancy_umap_source.tsv.gz",
    "copykat_direct_prediction_composition.tsv"
  )
)
for (i in seq_len(nrow(figure_source_pairs))) {
  pdf_path <- file.path(
    figure_dir, paste0(figure_source_pairs$figure[i], ".pdf")
  )
  png_path <- file.path(
    figure_dir, paste0(figure_source_pairs$figure[i], ".png")
  )
  source_path <- file.path(source_dir, figure_source_pairs$source[i])
  passed <- all(
    file.exists(c(pdf_path, png_path, source_path))
  ) && all(file.info(c(pdf_path, png_path, source_path))$size > 0)
  add_check(
    paste0("figure_source_pair:", figure_source_pairs$figure[i]),
    passed,
    paste0("source=", figure_source_pairs$source[i])
  )
}

evidence <- fread(file.path(
  project_root, "metadata", "annotation_evidence.tsv"
))
required_evidence_columns <- c(
  "level", "annotation", "positive_markers", "exclusion_markers",
  "primary_interpretation", "auxiliary_source", "decision_rule"
)
add_check(
  "annotation_evidence_columns",
  all(required_evidence_columns %in% colnames(evidence)),
  paste(colnames(evidence), collapse = ";")
)
add_check(
  "annotation_evidence_no_missing_marker_rules",
  all(
    !is.na(evidence$positive_markers) &
      evidence$positive_markers != "" &
      !is.na(evidence$exclusion_markers) &
      evidence$exclusion_markers != ""
  ),
  paste0("rows=", nrow(evidence))
)
required_major_types <- c(
  "Epithelial", "T_NK", "B_cell", "Plasma_cell", "Myeloid",
  "Fibroblast", "Endothelial", "Mast_cell"
)
required_epithelial_states <- c(
  "Stem_progenitor", "Cycling", "Absorptive", "BEST4_positive",
  "Goblet_secretory", "Adenoma_like"
)
add_check(
  "required_major_types_registered",
  all(required_major_types %chin% evidence[level == "major", annotation]),
  paste(
    setdiff(required_major_types, evidence[level == "major", annotation]),
    collapse = ";"
  )
)
add_check(
  "required_epithelial_states_registered",
  all(
    required_epithelial_states %chin%
      evidence[level == "epithelial", annotation]
  ),
  paste(
    setdiff(
      required_epithelial_states,
      evidence[level == "epithelial", annotation]
    ),
    collapse = ";"
  )
)

major_map <- fread(file.path(
  project_root, "config", "major_cluster_annotation.tsv"
))
epithelial_map <- fread(file.path(
  project_root, "config", "epithelial_cluster_annotation.tsv"
))
add_check(
  "major_map_unique_complete",
  !anyDuplicated(major_map$cluster) &&
    all(major_map$reviewer_context_blinded == "yes"),
  paste0("clusters=", nrow(major_map))
)
add_check(
  "epithelial_map_unique_complete",
  !anyDuplicated(epithelial_map$cluster) &&
    all(epithelial_map$reviewer_context_blinded == "yes"),
  paste0("clusters=", nrow(epithelial_map))
)

run_audit <- fread(file.path(
  result_dir, "copykat_run_audit.tsv"
))
unfinished_ready <- run_audit[
  status == "ready" & run_status != "completed"
]
add_check(
  "copykat_all_ready_tissues_completed",
  nrow(unfinished_ready) == 0L,
  paste0("unfinished_ready_tissues=", nrow(unfinished_ready))
)
add_check(
  "copykat_has_completed_tissues",
  run_audit[run_status == "completed", .N] > 0L,
  paste0("completed_tissues=", run_audit[run_status == "completed", .N])
)

dominance <- fread(file.path(
  result_dir, "epithelial_donor_dominance_audit.tsv"
))
add_check(
  "donor_dominance_all_clusters_audited",
  setequal(
    epithelial_map$cluster,
    dominance[level == "epithelial_cluster", as.character(group)]
  ),
  paste0(
    "audited_clusters=",
    dominance[level == "epithelial_cluster", uniqueN(group)]
  )
)

accepted <- readRDS(file.path(
  object_dir, "GSE201348_harmony_integrated.rds"
))
final <- readRDS(file.path(
  object_dir, "GSE201348_5C_annotated_final.rds"
))
epithelial <- readRDS(file.path(
  object_dir, "GSE201348_5C_epithelial_annotated_CNV.rds"
))
add_check(
  "final_dimensions_match_5B",
  identical(dim(final), dim(accepted)),
  paste0(
    "5B=", paste(dim(accepted), collapse = "x"),
    ";5C=", paste(dim(final), collapse = "x")
  )
)
add_check(
  "final_cell_names_match_5B",
  identical(colnames(final), colnames(accepted)),
  paste0("cells=", ncol(final))
)
accepted_nnz <- length(GetAssayData(
  accepted, assay = "RNA", slot = "counts"
)@x)
final_nnz <- length(GetAssayData(
  final, assay = "RNA", slot = "counts"
)@x)
add_check(
  "raw_count_nonzero_entries_preserved",
  accepted_nnz == final_nnz,
  paste0("5B=", accepted_nnz, ";5C=", final_nnz)
)
add_check(
  "all_cells_have_major_annotation",
  all(!is.na(final$major_cell_type) & final$major_cell_type != ""),
  paste0("missing=", sum(is.na(final$major_cell_type)))
)
add_check(
  "all_epithelial_cells_have_state",
  all(
    !is.na(epithelial$epithelial_state) &
      epithelial$epithelial_state != ""
  ),
  paste0("epithelial_cells=", ncol(epithelial))
)
add_check(
  "all_epithelial_cells_have_CNV_interpretation",
  all(
    !is.na(epithelial$malignancy_interpretation) &
      epithelial$malignancy_interpretation != ""
  ),
  paste0(
    "missing=",
    sum(is.na(epithelial$malignancy_interpretation))
  )
)
candidate <- which(
  epithelial$malignancy_interpretation ==
    "candidate_malignant_epithelial"
)
candidate_basis_valid <- if (length(candidate)) {
  all(
    epithelial$malignancy_call_basis[candidate] %chin% c(
      "direct_CopyKAT_cell_call",
      "same_tissue_epithelial_cluster_consensus"
    )
  )
} else {
  TRUE
}
add_check(
  "candidate_malignancy_has_CNV_basis",
  candidate_basis_valid,
  paste0("candidate_cells=", length(candidate))
)
add_check(
  "no_origin_or_EPCAM_only_malignancy_call",
  isFALSE(final@misc$stage_5C_CNV$tumor_origin_alone_used) &&
    isFALSE(final@misc$stage_5C_CNV$EPCAM_alone_used),
  "tumor_origin_alone=no;EPCAM_alone=no"
)

marker_files <- c(
  file.path(result_dir, "major_cluster_markers.tsv.gz"),
  file.path(result_dir, "epithelial_cluster_markers.tsv.gz")
)
for (path in marker_files) {
  marker_table <- fread(path)
  add_check(
    paste0("marker_table_nonempty:", basename(path)),
    nrow(marker_table) > 0L &&
      all(c("gene", "cluster") %in% colnames(marker_table)),
    paste0("rows=", nrow(marker_table))
  )
}

prohibited_result_patterns <- c(
  "trajectory", "pseudotime", "cellchat", "nichenet",
  "machine_learning", "classifier", "stage_DE", "differential_by_stage"
)
result_names <- list.files(
  result_dir, recursive = TRUE, full.names = FALSE
)
prohibited_hits <- result_names[
  Reduce(
    `|`,
    lapply(
      prohibited_result_patterns,
      function(pattern) grepl(
        pattern, result_names, ignore.case = TRUE, fixed = TRUE
      )
    )
  )
]
add_check(
  "stage_scope_no_prohibited_analysis_outputs",
  length(prohibited_hits) == 0L,
  paste(prohibited_hits, collapse = ";")
)

validation <- rbindlist(checks)
fwrite(
  validation,
  file.path(result_dir, "validation_checks.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "validation_sessionInfo.txt")
)
failed <- validation[passed == FALSE]
message(
  "Stage 5C validation: ",
  validation[passed == TRUE, .N], "/", nrow(validation), " passed"
)
if (nrow(failed)) {
  message(
    "Failed checks: ",
    paste(failed$check, collapse = ";")
  )
  quit(save = "no", status = 1L)
}
