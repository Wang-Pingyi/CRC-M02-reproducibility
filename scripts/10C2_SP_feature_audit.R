#!/usr/bin/env Rscript

# Stage 10C2-SP feature-space audit only. This script must never read ROI
# labels or calculate regional M02 scores.

packages <- c(
  "hdf5r", "rhdf5", "Matrix", "edgeR", "data.table",
  "SpatialExperiment", "Seurat", "sf", "spdep"
)

print_environment <- function() {
  cat(sprintf("R\t%s\n", R.version.string))
  for (package_name in packages) {
    version <- if (requireNamespace(package_name, quietly = TRUE)) {
      as.character(utils::packageVersion(package_name))
    } else {
      "NA"
    }
    cat(sprintf("%s\t%s\n", package_name, version))
  }
}

read_locked_genes <- function(lock_manifest) {
  lock <- utils::read.delim(
    lock_manifest,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(nrow(lock) == 36L, length(unique(lock$gene)) == 36L)
  lock[order(lock$gene_order), c("gene", "gene_order")]
}

audit_tenx_h5 <- function(h5_path, lock_manifest, cohort_id, sample_id) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r is required")
  }
  locked <- read_locked_genes(lock_manifest)
  h5 <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  feature_id <- as.character(h5[["matrix/features/id"]][])
  feature_name <- as.character(h5[["matrix/features/name"]][])
  if (length(feature_id) != length(feature_name)) {
    stop("feature id/name length mismatch")
  }
  target_index <- match(locked$gene, feature_name)
  target_duplicates <- vapply(
    locked$gene,
    function(gene) sum(feature_name == gene),
    integer(1)
  )
  nonzero <- rep(FALSE, nrow(locked))
  names(nonzero) <- locked$gene
  present_rows <- target_index[!is.na(target_index)] - 1L
  if (length(present_rows)) {
    indices <- h5[["matrix/indices"]][]
    data <- h5[["matrix/data"]][]
    nonzero_rows <- unique(indices[data > 0])
    nonzero[!is.na(target_index)] <- present_rows %in% nonzero_rows
  }
  data.frame(
    cohort_id = cohort_id,
    sample_id = sample_id,
    canonical_gene = locked$gene,
    gene_order = locked$gene_order,
    mapping_status = ifelse(
      is.na(target_index),
      "absent_from_feature_space",
      ifelse(target_duplicates == 1L, "exact_symbol", "unresolved")
    ),
    feature_id = ifelse(is.na(target_index), "NA", feature_id[target_index]),
    feature_symbol = ifelse(is.na(target_index), "NA", feature_name[target_index]),
    feature_duplicate_count = target_duplicates,
    present_but_zero = ifelse(is.na(target_index), "NA", ifelse(nonzero, "FALSE", "TRUE")),
    mapping_source = "deposited_10x_feature_annotation",
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || identical(args[[1]], "environment")) {
  print_environment()
} else if (identical(args[[1]], "h5")) {
  if (length(args) != 6L) {
    stop("usage: h5 <h5_path> <lock_manifest> <cohort_id> <sample_id> <output_tsv>")
  }
  result <- audit_tenx_h5(args[[2]], args[[3]], args[[4]], args[[5]])
  utils::write.table(
    result,
    file = args[[6]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
} else {
  stop("unknown command")
}
