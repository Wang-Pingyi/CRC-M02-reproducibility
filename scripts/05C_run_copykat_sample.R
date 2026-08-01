#!/usr/bin/env Rscript

# Run CopyKAT for one biological tissue input prepared by
# 05C_prepare_copykat_inputs.R.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: 05C_run_copykat_sample.R PROJECT_ROOT SAMPLE_KEY")
}
project_root <- normalizePath(args[[1]], mustWork = TRUE)
sample_key <- args[[2]]
library_dir <- file.path(project_root, "environment", "R_library")
.libPaths(c(library_dir, .libPaths()))
if (!requireNamespace("copykat", quietly = TRUE)) {
  stop("CopyKAT is not installed in the project environment")
}

# CopyKAT 1.2.5 stores its internal genome annotations in the public
# `sysdata` dataset, but the four functions below refer to those objects as
# namespace-level variables. The pinned upstream revision does not load that
# dataset into its namespace when installed from source. Load the official
# package dataset into a process-local child environment and rebind only the
# affected function environments. This does not modify the installed package.
copykat_namespace <- asNamespace("copykat")
copykat_sysdata <- new.env(parent = copykat_namespace)
suppressWarnings(utils::data(
  "sysdata", package = "copykat", envir = copykat_sysdata
))
required_sysdata <- c(
  "full.anno", "full.anno.mm10", "DNA.hg20", "cyclegenes"
)
if (!all(vapply(
  required_sysdata,
  exists,
  logical(1),
  envir = copykat_sysdata,
  inherits = FALSE
))) {
  stop("CopyKAT sysdata genome annotations could not be loaded")
}
sysdata_function_names <- c(
  "annotateGenes.hg20",
  "annotateGenes.mm10",
  "convert.all.bins.hg20",
  "copykat"
)
for (function_name in sysdata_function_names) {
  patched_function <- get(
    function_name, envir = copykat_namespace, inherits = FALSE
  )
  environment(patched_function) <- copykat_sysdata
  assignInNamespace(
    function_name, patched_function, ns = "copykat"
  )
}

param_file <- file.path(project_root, "config", "annotation_parameters.tsv")
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

input_file <- file.path(
  project_root, "cache", "05C_copykat_inputs", paste0(sample_key, ".rds")
)
output_dir <- file.path(
  project_root, "cache", "05C_copykat_outputs", sample_key
)
log_dir <- file.path(project_root, "logs", "05C_copykat")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(output_dir, "run_status.tsv")
start_time <- Sys.time()

write_status <- function(status, message_text, predictions = NA_integer_) {
  fwrite(
    data.table(
      sample_key = sample_key,
      status = status,
      message = message_text,
      start_time = format(start_time, "%Y-%m-%dT%H:%M:%S%z"),
      end_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      elapsed_minutes = as.numeric(
        difftime(Sys.time(), start_time, units = "mins")
      ),
      predictions = predictions,
      copykat_version = as.character(packageVersion("copykat")),
      copykat_remote_sha = packageDescription("copykat")$RemoteSha,
      sysdata_runtime_patch = TRUE
    ),
    status_file,
    sep = "\t",
    quote = TRUE,
    na = "NA"
  )
}

tryCatch({
  input <- readRDS(input_file)
  if (!all(input$normal_reference_cell_ids %in% colnames(
    input$raw_counts
  ))) {
    stop("Known normal reference cell names are missing from raw matrix")
  }
  if (length(intersect(
    input$query_cell_ids, input$normal_reference_cell_ids
  ))) {
    stop("Query and reference cells overlap")
  }
  set.seed(input$seed)
  raw_dense <- as.matrix(input$raw_counts)
  storage.mode(raw_dense) <- "numeric"

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  fit <- copykat::copykat(
    rawmat = raw_dense,
    id.type = "S",
    cell.line = "no",
    ngene.chr = as.integer(get_param("cnv", "ngene_chr", TRUE)),
    win.size = as.integer(get_param("cnv", "win_size", TRUE)),
    KS.cut = get_param("cnv", "ks_cut", TRUE),
    sam.name = sample_key,
    distance = "euclidean",
    norm.cell.names = input$normal_reference_cell_ids,
    output.seg = "FALSE",
    plot.genes = "FALSE",
    genome = "hg20",
    n.cores = as.integer(get_param("cnv", "n_cores", TRUE))
  )
  prediction <- as.data.table(fit$prediction)
  cell_column <- intersect(
    c("cell.names", "cell_name", "cell"), colnames(prediction)
  )
  prediction_column <- intersect(
    c("copykat.pred", "copykat_pred", "prediction"), colnames(prediction)
  )
  if (length(cell_column) != 1L || length(prediction_column) != 1L) {
    stop(
      "Unexpected CopyKAT prediction columns: ",
      paste(colnames(prediction), collapse = ";")
    )
  }
  setnames(
    prediction,
    c(cell_column, prediction_column),
    c("cell_id", "copykat_prediction")
  )
  manifest <- copy(input$cell_manifest)
  prediction <- merge(
    manifest,
    prediction[, .(cell_id, copykat_prediction)],
    by = "cell_id",
    all.x = TRUE
  )
  prediction[
    is.na(copykat_prediction),
    copykat_prediction := "not.defined"
  ]
  prediction[, sample_key := input$sample_key]
  prediction[
    ,
    queried_epithelial := copykat_role == "epithelial_query"
  ]
  fwrite(
    prediction,
    file.path(output_dir, "copykat_predictions.tsv.gz"),
    sep = "\t",
    quote = TRUE,
    na = "NA"
  )
  saveRDS(
    fit,
    file.path(output_dir, "copykat_fit.rds"),
    compress = "gzip"
  )
  fwrite(
    prediction[
      ,
      .N,
      by = .(copykat_role, copykat_prediction)
    ][order(copykat_role, copykat_prediction)],
    file.path(output_dir, "prediction_summary.tsv"),
    sep = "\t",
    quote = TRUE,
    na = "NA"
  )
  write_status("completed", "CopyKAT completed", nrow(prediction))
  message(
    sample_key, ": CopyKAT completed with ", nrow(prediction),
    " predictions"
  )
}, error = function(e) {
  write_status("failed", conditionMessage(e))
  message(sample_key, ": CopyKAT failed: ", conditionMessage(e))
  quit(save = "no", status = 1L)
})
