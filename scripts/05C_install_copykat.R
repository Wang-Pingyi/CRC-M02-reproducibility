#!/usr/bin/env Rscript

# Install the pinned official CopyKAT revision into the project R library.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}
param_file <- file.path(project_root, "config", "annotation_parameters.tsv")
library_dir <- file.path(project_root, "environment", "R_library")
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

params <- fread(param_file)
pinned_sha <- params[
  section == "cnv" & parameter == "copykat_git_commit",
  value
]
if (length(pinned_sha) != 1L) stop("Missing pinned CopyKAT commit")
get_cnv_parameter <- function(parameter_name) {
  value <- params[
    section == "cnv" & parameter == parameter_name,
    value
  ]
  if (length(value) != 1L) {
    stop("Missing or duplicated cnv parameter: ", parameter_name)
  }
  value
}
transport_version <- get_cnv_parameter("transport_version")
transport_url <- get_cnv_parameter("transport_source_url")
transport_sha256 <- get_cnv_parameter("transport_sha256")

installed_transport <- if (requireNamespace("transport", quietly = TRUE)) {
  as.character(packageVersion("transport"))
} else {
  NA_character_
}
transport_version_matches <- function(installed, pinned) {
  !is.na(installed) &&
    package_version(installed) == package_version(pinned)
}
if (!transport_version_matches(
  installed_transport, transport_version
)) {
  transport_archive <- tempfile(
    pattern = paste0("transport_", transport_version, "_"),
    fileext = ".tar.gz"
  )
  on.exit(unlink(transport_archive), add = TRUE)
  curl_status <- system2(
    "curl",
    c(
      "--fail", "--location",
      "--retry", "8",
      "--retry-delay", "2",
      "--retry-all-errors",
      "--connect-timeout", "30",
      "--max-time", "300",
      "--output", transport_archive,
      transport_url
    )
  )
  if (!identical(curl_status, 0L)) {
    stop("Failed to download pinned transport archive with curl")
  }
  sha_output <- system2(
    "sha256sum",
    transport_archive,
    stdout = TRUE,
    stderr = TRUE
  )
  observed_sha256 <- strsplit(sha_output[[1]], "[[:space:]]+")[[1]][1]
  if (!identical(observed_sha256, transport_sha256)) {
    stop(
      "Pinned transport archive SHA256 mismatch: observed=",
      observed_sha256, "; expected=", transport_sha256
    )
  }
  install.packages(
    transport_archive,
    repos = NULL,
    type = "source",
    lib = library_dir
  )
}
if (
  !requireNamespace("transport", quietly = TRUE) ||
    !transport_version_matches(
      as.character(packageVersion("transport")),
      transport_version
    )
) {
  stop("Pinned transport installation failed")
}

installed_sha <- NA_character_
if (requireNamespace("copykat", quietly = TRUE)) {
  installed_sha <- packageDescription("copykat")$RemoteSha
}
if (is.null(installed_sha)) installed_sha <- NA_character_

if (is.na(installed_sha) || !startsWith(pinned_sha, installed_sha)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages(
      "remotes",
      repos = "https://cloud.r-project.org",
      lib = library_dir
    )
  }
  remotes::install_github(
    paste0("navinlabcode/copykat@", pinned_sha),
    lib = library_dir,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never",
    force = TRUE
  )
}

if (!requireNamespace("copykat", quietly = TRUE)) {
  stop("CopyKAT installation failed")
}
description <- packageDescription("copykat")
installed_sha <- description$RemoteSha
if (is.null(installed_sha) || !startsWith(pinned_sha, installed_sha)) {
  stop(
    "Installed CopyKAT revision does not match pin: installed=",
    installed_sha, "; expected=", pinned_sha
  )
}

result_dir <- file.path(project_root, "results", "05C_annotation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    paste0("package=copykat"),
    paste0("version=", as.character(packageVersion("copykat"))),
    paste0("RemoteSha=", installed_sha),
    paste0("pinned_sha=", pinned_sha),
    paste0(
      "transport_version=",
      as.character(packageVersion("transport"))
    ),
    paste0("transport_source_url=", transport_url),
    paste0("transport_sha256=", transport_sha256),
    paste0(
      "RcppEigen_version=",
      as.character(packageVersion("RcppEigen"))
    ),
    paste0("library=", find.package("copykat"))
  ),
  file.path(result_dir, "copykat_package_provenance.txt")
)
message(
  "CopyKAT ", packageVersion("copykat"),
  " installed at pinned revision ", installed_sha
)
