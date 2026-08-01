#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1L]] else getwd()
targets <- c(
  "scripts/06A_preflight.R",
  "scripts/06A_smoke_test.R",
  "scripts/06A_pseudobulk_discovery.R",
  "scripts/06A_validate_outputs.R",
  "scripts/06A_independent_audit.R",
  "scripts/06A_amendment_smoke_test.R",
  "scripts/06A_amendment_exploratory.R",
  "scripts/06A_validate_amendment.R"
)

for (target in targets) {
  path <- file.path(project_dir, target)
  parse(file = path)
  cat("PARSE_OK\t", target, "\n", sep = "")
}
