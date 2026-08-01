#!/usr/bin/env Rscript

# Environment setup: Stage 6B project-private R dependencies
# Date: 2026-07-28
# This script does not modify the system R library.

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  download.file.method = "libcurl",
  timeout = 1200
)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
source_dir <- file.path(project_dir, "environment", "sources")
dir.create(private_library, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(private_library, .libPaths()))

download_with_curl <- function(url, destination, expected_md5) {
  if (file.exists(destination)) {
    observed <- unname(tools::md5sum(destination))
    if (identical(tolower(observed), tolower(expected_md5))) return(invisible(destination))
    unlink(destination)
  }
  status <- system2(
    "curl",
    c(
      "-L", "--fail", "--retry", "8", "--retry-delay", "5",
      "--connect-timeout", "20", "--max-time", "1200",
      "-o", shQuote(destination), shQuote(url)
    )
  )
  if (status != 0L || !file.exists(destination)) {
    stop("Dependency download failed: ", url)
  }
  observed <- unname(tools::md5sum(destination))
  if (!identical(tolower(observed), tolower(expected_md5))) {
    stop(
      "Dependency MD5 mismatch for ", basename(destination),
      ": expected ", expected_md5, ", observed ", observed
    )
  }
  invisible(destination)
}

install_source <- function(package, archive) {
  if (requireNamespace(package, quietly = TRUE)) return(invisible(TRUE))
  utils::install.packages(
    archive,
    repos = NULL,
    type = "source",
    lib = private_library,
    dependencies = FALSE,
    INSTALL_opts = c("--no-multiarch", "--no-test-load")
  )
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Dependency installation failed: ", package)
  }
  invisible(TRUE)
}

mirror <- "https://bioconductor.statistik.tu-dortmund.de/packages/3.18"
sources <- data.frame(
  package = c("bcellViper", "tradeSeq", "dorothea", "logger", "OmnipathR"),
  version = c("1.38.0", "1.16.0", "1.14.1", "0.4.2", "3.10.1"),
  url = c(
    sprintf("%s/data/experiment/src/contrib/bcellViper_1.38.0.tar.gz", mirror),
    sprintf("%s/bioc/src/contrib/tradeSeq_1.16.0.tar.gz", mirror),
    sprintf("%s/data/experiment/src/contrib/dorothea_1.14.1.tar.gz", mirror),
    "https://cloud.r-project.org/src/contrib/logger_0.4.2.tar.gz",
    sprintf("%s/bioc/src/contrib/OmnipathR_3.10.1.tar.gz", mirror)
  ),
  md5 = c(
    "7e93bbaa204826358c77282e2a370074",
    "5de25ecb5472b592837cd3937ff467c7",
    "962242cff51a595005ff65d6dad90121",
    "6f6a3cc1ed3caf7c0dc8d8efc305c841",
    "c51296b40e8d580f14dc3ff6839495f6"
  ),
  stringsAsFactors = FALSE
)
sources$archive <- file.path(
  source_dir,
  sprintf("%s_%s.tar.gz", sources$package, sources$version)
)

for (i in seq_len(nrow(sources))) {
  if (!requireNamespace(sources$package[i], quietly = TRUE)) {
    download_with_curl(sources$url[i], sources$archive[i], sources$md5[i])
    install_source(sources$package[i], sources$archive[i])
  }
}

liana_archive <- file.path(source_dir, "liana_6cab46c54234.tar.gz")
liana_sha256 <- "4c839c52084c1a90aa0e9d6f97b4f3e72e17b1f58ddad1060ed372349b39d168"
if (!file.exists(liana_archive)) {
  stop("Pinned LIANA source archive is missing: ", liana_archive)
}
observed_liana_sha256 <- strsplit(
  system2("sha256sum", shQuote(liana_archive), stdout = TRUE),
  "[[:space:]]+"
)[[1L]][1L]
if (!identical(tolower(observed_liana_sha256), liana_sha256)) {
  stop("Pinned LIANA source SHA256 mismatch")
}

if (!requireNamespace("liana", quietly = TRUE)) {
  liana_imports <- c(
    "reticulate", "stringr", "magrittr", "purrr", "tidyr", "tibble",
    "dplyr", "readr", "rlang", "OmnipathR", "SingleCellExperiment",
    "scran", "scater", "scuttle", "SeuratObject", "tidyselect",
    "ggplot2", "ComplexHeatmap", "RColorBrewer", "basilisk.utils", "basilisk"
  )
  missing_imports <- liana_imports[
    !vapply(liana_imports, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_imports)) {
    stop(
      "Pinned LIANA imports are missing: ",
      paste(missing_imports, collapse = ", ")
    )
  }
  # The pinned upstream archive stores configure with CRLF line endings. Linux
  # cannot execute its /bin/bash shebang until the line endings are normalized.
  # Apply this packaging-only correction to a temporary copy; retain and verify
  # the original source archive unchanged above.
  liana_build_dir <- tempfile("liana-build-")
  dir.create(liana_build_dir, recursive = TRUE)
  utils::untar(liana_archive, exdir = liana_build_dir)
  liana_configure <- file.path(liana_build_dir, "liana", "configure")
  if (!file.exists(liana_configure)) stop("LIANA configure script is missing")
  configure_lines <- readLines(liana_configure, warn = FALSE)
  writeLines(configure_lines, liana_configure, useBytes = TRUE)
  Sys.chmod(liana_configure, mode = "0755")
  corrected_archive <- tempfile("liana-lf-", fileext = ".tar.gz")
  previous_dir <- setwd(liana_build_dir)
  tar_status <- system2(
    "tar", c("-czf", shQuote(corrected_archive), "liana")
  )
  setwd(previous_dir)
  if (tar_status != 0L || !file.exists(corrected_archive)) {
    stop("Could not create the temporary LF-normalized LIANA archive")
  }
  install_source("liana", corrected_archive)
}

required <- c("tradeSeq", "dorothea", "liana")
still_missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing)) {
  stop("Stage 6B dependencies remain unavailable: ", paste(still_missing, collapse = ", "))
}

versions <- data.frame(
  package = required,
  version = vapply(
    required, function(x) as.character(utils::packageVersion(x)), character(1)
  ),
  library_path = vapply(
    required,
    function(x) normalizePath(find.package(x), winslash = "/", mustWork = TRUE),
    character(1)
  ),
  stringsAsFactors = FALSE
)
output_dir <- file.path(project_dir, "results", "06B_regulatory_inference", "preflight")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  versions,
  file.path(output_dir, "installed_dependency_versions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

source_manifest <- rbind(
  transform(
    sources[, c("package", "version", "url", "md5")],
    installation_transform = "none"
  ),
  data.frame(
    package = "liana",
    version = "git-6cab46c54234",
    url = "https://github.com/saezlab/liana/commit/6cab46c54234f861ea176c3de77c4b8aa45ecb3d",
    md5 = paste0("SHA256:", liana_sha256),
    installation_transform = "temporary configure CRLF-to-LF normalization",
    stringsAsFactors = FALSE
  )
)
utils::write.table(
  source_manifest,
  file.path(output_dir, "dependency_source_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(versions, row.names = FALSE)
cat("STAGE_6B_DEPENDENCIES_OK\n")
