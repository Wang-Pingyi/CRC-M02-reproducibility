#!/usr/bin/env Rscript

# Stage 6A preflight for GSE201348 donor-level pseudobulk analysis.
# Random seed: 20260728

set.seed(20260728)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1L]] else getwd()
object_path <- file.path(
  project_dir,
  "objects",
  "GSE201348_5C_epithelial_annotated_CNV.rds"
)

stopifnot(file.exists(object_path))

cat("timestamp\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("R_version\t", R.version.string, "\n", sep = "")

packages <- c("Seurat", "SeuratObject", "Matrix", "edgeR", "limma")
for (pkg in packages) {
  available <- requireNamespace(pkg, quietly = TRUE)
  version <- if (available) as.character(utils::packageVersion(pkg)) else "NA"
  cat("package\t", pkg, "\t", available, "\t", version, "\n", sep = "")
}

obj <- readRDS(object_path)
meta <- obj[[]]

cat("object_class\t", paste(class(obj), collapse = ";"), "\n", sep = "")
cat("cells\t", ncol(obj), "\n", sep = "")
cat("features\t", nrow(obj), "\n", sep = "")
cat("metadata_columns\t", paste(colnames(meta), collapse = ";"), "\n", sep = "")
cat("assays\t", paste(names(obj@assays), collapse = ";"), "\n", sep = "")
cat("default_assay\t", SeuratObject::DefaultAssay(obj), "\n", sep = "")

candidate_columns <- c(
  "donor_id", "patient_id", "condition", "lesion_stage", "stage",
  "epithelial_state", "epithelial_state_final", "histology",
  "sporadic_or_FAP", "fap_status", "tumor_location",
  "colon_or_rectum", "sample_id", "biological_sample_id", "library_id"
)

for (column in intersect(candidate_columns, colnames(meta))) {
  values <- as.character(meta[[column]])
  tab <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
  rendered <- paste(names(tab), as.integer(tab), sep = "=", collapse = ";")
  cat("metadata_summary\t", column, "\t", rendered, "\n", sep = "")
}

cross_fields <- intersect(
  c(
    "donor_id", "lesion_stage", "condition", "sporadic_or_FAP",
    "colon_or_rectum", "tumor_location"
  ),
  colnames(meta)
)
crosswalk_input <- meta[, cross_fields, drop = FALSE]
crosswalk_input$cell_count <- 1L
crosswalk <- stats::aggregate(
  cell_count ~ .,
  data = crosswalk_input,
  FUN = sum,
  na.action = stats::na.pass
)
crosswalk <- crosswalk[do.call(order, crosswalk[c("donor_id", "lesion_stage")]), , drop = FALSE]
cat("crosswalk_begin\n")
utils::write.table(
  crosswalk,
  file = stdout(),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
cat("crosswalk_end\n")

counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
cat("counts_class\t", paste(class(counts), collapse = ";"), "\n", sep = "")
cat("counts_dimensions\t", paste(dim(counts), collapse = "x"), "\n", sep = "")
cat("counts_nonzero\t", length(counts@x), "\n", sep = "")
cat("counts_integer_like\t", all(abs(counts@x - round(counts@x)) < 1e-8), "\n", sep = "")

sessionInfo()
