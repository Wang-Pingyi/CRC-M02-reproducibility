#!/usr/bin/env Rscript

# Rebuild preprocessCore without pthread support to avoid the documented
# OpenBLAS/pthread_create(EINVAL) failure during RMA on this Linux server.

user_lib <- .libPaths()[1L]
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
source_tarball <- if (length(args)) args[[1L]] else file.path(
  "environment", "sources", "preprocessCore_1.64.0.tar.gz"
)
if (!file.exists(source_tarball)) {
  stop("Official preprocessCore source tarball is missing: ", source_tarball)
}

install.packages(
  source_tarball,
  lib = user_lib,
  repos = NULL,
  type = "source",
  configure.args = "--disable-threading",
  clean = TRUE
)

pkg_path <- find.package("preprocessCore", lib.loc = user_lib)
desc <- read.dcf(file.path(pkg_path, "DESCRIPTION"))
cat(
  "Installed preprocessCore ", desc[1L, "Version"],
  " without threading in ", pkg_path, "\n",
  sep = ""
)
