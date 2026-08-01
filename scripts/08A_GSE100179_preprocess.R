#!/usr/bin/env Rscript
set.seed(20260728)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08A_GSE100179_preprocess.R <project_dir> <run_id>")
source(file.path(args[[1]], "scripts", "08A_helpers.R"))
stage8_preprocess_gse100179(normalizePath(args[[1]], mustWork = TRUE), args[[2]])
