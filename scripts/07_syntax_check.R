#!/usr/bin/env Rscript

scripts <- c(
  "scripts/07_helpers.R",
  "scripts/07_install_dependencies.R",
  "scripts/07_smoke_test.R",
  "scripts/07_preflight.R",
  "scripts/07_real_data_pilot.R",
  "scripts/07_prepare_GSE161277.R",
  "scripts/07_resume_GSE161277_figures.R",
  "scripts/07_prepare_GSE132465.R",
  "scripts/07_model_preflight.R",
  "scripts/07_run_replication_models.R",
  "scripts/07_finalize_report.R",
  "scripts/07_validate_outputs.R",
  "scripts/07_acceptance_qc.R",
  "scripts/07_update_result_registry.R"
)
for (script in scripts) {
  if (!file.exists(script)) stop("Missing script: ", script)
  parse(file = script)
  cat("PARSE_OK\t", script, "\n", sep = "")
}
