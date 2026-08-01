#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else getwd()
result_dir <- file.path(project_root, "results", "05B_full_qc_integration")
object_dir <- file.path(project_root, "objects")

qc <- fread(file.path(result_dir, "library_qc_summary.tsv"))
expected_cells <- sum(qc$cells_after)
expected_features <- 33538L
included_meta <- fread(file.path(result_dir, "source_data", "included_library_metadata.tsv"))
active_samples <- qc[cells_after > 0, sample_id]
expected_libraries <- length(active_samples)
expected_biological_samples <- uniqueN(
  included_meta[sample_id %chin% active_samples, biological_sample_id]
)
expected_donors <- uniqueN(included_meta[sample_id %chin% active_samples, donor_id])
required_metadata <- c(
  "sample_id", "donor_id", "biological_sample_id", "lesion_stage",
  "sporadic_or_FAP", "platform", "nCount_RNA", "nFeature_RNA",
  "percent_mt", "percent_ribo", "scDblFinder_score"
)
forbidden_annotation <- c("seurat_clusters", "cell_type", "celltype", "annotation")

checks <- list()
add_check <- function(object, check, value, pass) {
  checks[[length(checks) + 1L]] <<- data.table(
    object = object,
    check = check,
    value = as.character(value),
    status = if (isTRUE(pass)) "PASS" else "FAIL"
  )
}

included_qc <- qc[library_below_min_cells == FALSE]
excluded_qc <- qc[library_below_min_cells == TRUE]
add_check(
  "qc_summary", "included_doublet_reconciliation",
  sum(included_qc$doublets_identified - included_qc$doublet_final_exclusion),
  all(included_qc$doublets_identified == included_qc$doublet_final_exclusion)
)
add_check(
  "qc_summary", "included_retained_reconciliation",
  sum(included_qc$cells_entering_doublet_detection -
        included_qc$doublet_final_exclusion - included_qc$cells_after),
  all(included_qc$cells_after ==
        included_qc$cells_entering_doublet_detection -
        included_qc$doublet_final_exclusion)
)
add_check(
  "qc_summary", "excluded_library_reconciliation",
  nrow(excluded_qc),
  all(excluded_qc$cells_after == 0L) &&
    all(excluded_qc$library_gate_excluded ==
          excluded_qc$cells_entering_doublet_detection)
)

cell_qc_files <- list.files(
  file.path(result_dir, "source_data"),
  pattern = "^GSM[0-9]+_cell_qc_metrics[.]tsv[.]gz$",
  full.names = TRUE
)
cell_qc_audit <- rbindlist(lapply(cell_qc_files, function(path) {
  x <- fread(
    path,
    select = c(
      "nCount", "nFeature", "percent_mt", "scDblFinder_class", "retained"
    )
  )
  x[retained == TRUE, .(
    retained_cells = .N,
    retained_doublet_cells = sum(scDblFinder_class == "doublet", na.rm = TRUE),
    published_gate_violations = sum(
      nFeature <= 400L |
        nFeature >= 4000L |
        nCount >= 10000L |
        percent_mt >= 5
    )
  )]
}))
add_check(
  "cell_qc_source", "source_file_count",
  length(cell_qc_files), length(cell_qc_files) == nrow(qc)
)
add_check(
  "cell_qc_source", "retained_cell_reconciliation",
  sum(cell_qc_audit$retained_cells),
  sum(cell_qc_audit$retained_cells) == expected_cells
)
add_check(
  "cell_qc_source", "retained_doublet_cells",
  sum(cell_qc_audit$retained_doublet_cells),
  sum(cell_qc_audit$retained_doublet_cells) == 0L
)
add_check(
  "cell_qc_source", "published_gate_violations",
  sum(cell_qc_audit$published_gate_violations),
  sum(cell_qc_audit$published_gate_violations) == 0L
)

validate_object <- function(path, object_name, required_reductions,
                            forbidden_reductions = character()) {
  obj <- readRDS(path)
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  reductions <- names(obj@reductions)
  metadata_names <- colnames(obj@meta.data)

  add_check(object_name, "rds_readable", TRUE, inherits(obj, "Seurat"))
  add_check(object_name, "cell_count", ncol(obj), ncol(obj) == expected_cells)
  add_check(object_name, "feature_count", nrow(obj), nrow(obj) == expected_features)
  add_check(object_name, "raw_counts_dimensions",
            paste(dim(counts), collapse = "x"),
            identical(dim(counts), c(expected_features, expected_cells)))
  add_check(object_name, "raw_counts_nonzero_entries", length(counts@x),
            length(counts@x) > 0L)
  published_gate_violations <- sum(
    obj$nFeature_RNA <= 400L |
      obj$nFeature_RNA >= 4000L |
      obj$nCount_RNA >= 10000L |
      obj$percent_mt >= 5
  )
  add_check(object_name, "published_gate_violations",
            published_gate_violations, published_gate_violations == 0L)
  add_check(object_name, "required_metadata_present",
            paste(setdiff(required_metadata, metadata_names), collapse = ";"),
            all(required_metadata %in% metadata_names))
  add_check(object_name, "unique_libraries",
            uniqueN(obj$sample_id), uniqueN(obj$sample_id) == expected_libraries)
  add_check(object_name, "unique_biological_samples",
            uniqueN(obj$biological_sample_id),
            uniqueN(obj$biological_sample_id) == expected_biological_samples)
  add_check(object_name, "unique_donors",
            uniqueN(obj$donor_id), uniqueN(obj$donor_id) == expected_donors)
  add_check(object_name, "required_reductions",
            paste(reductions, collapse = ";"),
            all(required_reductions %in% reductions))
  add_check(object_name, "forbidden_reductions_absent",
            paste(intersect(forbidden_reductions, reductions), collapse = ";"),
            !any(forbidden_reductions %in% reductions))
  add_check(object_name, "annotation_fields_absent",
            paste(intersect(forbidden_annotation, metadata_names), collapse = ";"),
            !any(forbidden_annotation %in% metadata_names))

  reduction_dims <- vapply(required_reductions, function(x) {
    ncol(Embeddings(obj, reduction = x))
  }, integer(1))
  add_check(object_name, "reduction_dimensions",
            paste(names(reduction_dims), reduction_dims, sep = "=", collapse = ";"),
            all(reduction_dims == 30L))

  rm(counts, obj)
  invisible(gc())
}

validate_object(
  file.path(object_dir, "GSE201348_unintegrated_normalized.rds"),
  "unintegrated_normalized",
  required_reductions = "pca",
  forbidden_reductions = "harmony"
)
validate_object(
  file.path(object_dir, "GSE201348_harmony_integrated.rds"),
  "harmony_integrated",
  required_reductions = c("pca", "harmony")
)

validation <- rbindlist(checks)
fwrite(
  validation,
  file.path(result_dir, "validation_checks.tsv"),
  sep = "\t",
  quote = TRUE
)
if (any(validation$status != "PASS")) {
  print(validation[status != "PASS"])
  stop("Stage 5B object validation failed")
}
message("Stage 5B object validation passed: ", nrow(validation), " checks")
