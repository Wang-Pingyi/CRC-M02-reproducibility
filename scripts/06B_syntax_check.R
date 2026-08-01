#!/usr/bin/env Rscript

# Syntax validation: Stage 6B scripts

scripts <- c(
  "scripts/06B_preflight.R",
  "scripts/06B_install_dependencies.R",
  "scripts/06B_real_data_pilot.R",
  "scripts/06B_trajectory_tradeSeq.R",
  "scripts/06B_regulator_pathway.R",
  "scripts/06B_liana_donor_communication.R",
  "scripts/06B_finalize_report.R",
  "scripts/06B_validate_outputs.R",
  "scripts/06B_acceptance_qc.R",
  "scripts/06B_update_result_registry.R"
)
for (script in scripts) {
  if (!file.exists(script)) stop("Missing script: ", script)
  parse(file = script)
  cat("PARSE_OK\t", script, "\n", sep = "")
}
