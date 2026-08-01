#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)

required <- c(
  file.path(paths$result, "replication_summary.tsv"),
  file.path(paths$result, "replication_effects_all.tsv"),
  file.path(paths$result, "replication_gate.tsv"),
  file.path(paths$result, "GSE161277_qc_summary.tsv"),
  file.path(paths$result, "GSE132465_qc_summary.tsv"),
  file.path(paths$result, "GSE161277_annotation_evidence.tsv"),
  file.path(paths$result, "GSE132465_annotation_evidence.tsv"),
  file.path(
    paths$result, "preflight", "GSE132465_pair_eligibility.tsv"
  ),
  file.path(paths$source, "paired_donor_module_differences.tsv"),
  file.path(paths$source, "epithelial_specificity_donor_differences.tsv"),
  file.path(paths$source, "stage_7_module_effects_source.tsv"),
  file.path(paths$figure, "stage_7_module_effects.pdf"),
  file.path(paths$figure, "stage_7_module_effects.png"),
  file.path(project_dir, "reports", "stage_7_singlecell_replication.md"),
  file.path(project_dir, "STATUS.md"),
  file.path(paths$object, "GSE161277_stage7_annotated.rds"),
  file.path(paths$object, "GSE132465_stage7_processed_summary.rds"),
  file.path(paths$processed, "GSE161277_stage7_pseudobulk_raw_counts.rds"),
  file.path(paths$processed, "GSE132465_stage7_pseudobulk_raw_counts.rds")
)

summary <- utils::read.delim(
  file.path(paths$result, "replication_summary.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
effects <- utils::read.delim(
  file.path(paths$result, "replication_effects_all.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
gate <- utils::read.delim(
  file.path(paths$result, "replication_gate.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
primary <- utils::read.delim(
  file.path(project_dir, "results", "06A_pseudobulk", "candidate_programs.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
allowed <- c(
  "replicated", "directionally consistent but underpowered",
  "not replicated", "contradictory"
)
checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

add_check(
  "required_outputs_exist",
  all(file.exists(required)),
  paste(basename(required[!file.exists(required)]), collapse = ";")
)
add_check(
  "required_outputs_nonempty",
  all(file.info(required)$size > 0),
  paste(basename(required[file.info(required)$size <= 0]), collapse = ";")
)
add_check("primary_candidate_table_still_empty", nrow(primary) == 0L, nrow(primary))
add_check("six_locked_modules_reported", nrow(summary) == 6L, nrow(summary))
add_check(
  "unique_module_rows", !anyDuplicated(summary$module_id),
  paste(summary$module_id, collapse = ";")
)
add_check(
  "classification_vocabulary",
  all(summary$replication_class %in% allowed),
  paste(unique(summary$replication_class), collapse = ";")
)
add_check(
  "no_validation_gene_reselection",
  all(!summary$validation_gene_reselection),
  paste(unique(summary$validation_gene_reselection), collapse = ";")
)
add_check("all_prespecified_effect_rows", nrow(effects) == 30L, nrow(effects))
add_check(
  "unique_cohort_contrast_module_effects",
  !anyDuplicated(paste(effects$cohort, effects$contrast, effects$module_id)),
  "expected unique cohort-contrast-module rows"
)
add_check(
  "valid_confidence_intervals",
  all(
    is.na(effects$ci_low) | is.na(effects$ci_high) |
      effects$ci_low <= effects$effect & effects$effect <= effects$ci_high
  ),
  "effect lies within reported interval"
)
add_check(
  "valid_p_values",
  all(is.na(effects$p_value) | effects$p_value >= 0 & effects$p_value <= 1),
  "P in [0,1]"
)
add_check(
  "valid_fdr_values",
  all(is.na(effects$fdr) | effects$fdr >= 0 & effects$fdr <= 1),
  "FDR in [0,1]"
)
g161 <- effects[
  effects$cohort == "GSE161277" & effects$contrast != "epithelial_specificity",
]
add_check(
  "GSE161277_matched_donor_count",
  all(g161$n_donors == 3L),
  paste(unique(g161$n_donors), collapse = ";")
)
g132 <- effects[
  effects$cohort == "GSE132465" & effects$contrast == "normal_to_cancer",
]
eligibility132 <- utils::read.delim(
  file.path(
    paths$result, "preflight", "GSE132465_pair_eligibility.tsv"
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)
expected_pairs132 <- length(unique(
  eligibility132$donor_id[eligibility132$paired_analysis_eligible]
))
add_check(
  "GSE132465_matched_donor_count",
  all(g132$n_donors == expected_pairs132) && expected_pairs132 == 8L,
  paste0(
    "modeled=", paste(unique(g132$n_donors), collapse = ";"),
    "; eligible=", expected_pairs132
  )
)
add_check(
  "GSE132465_pair_exclusions_traceable",
  setequal(
    unique(eligibility132$donor_id[
      !eligibility132$paired_analysis_eligible
    ]),
    c("SMC05", "SMC08")
  ),
  paste(
    unique(eligibility132$donor_id[
      !eligibility132$paired_analysis_eligible
    ]),
    collapse = ";"
  )
)
add_check(
  "gate_has_all_criteria",
  nrow(gate) == 5L && sum(gate$criterion == "overall_replication_gate") == 1L,
  paste(gate$criterion, collapse = ";")
)
add_check(
  "raw_data_unchanged_by_hash",
  file.exists(file.path(paths$result, "raw_inputs.after.sha256")) &&
    file.exists(file.path(paths$result, "raw_inputs.before.sha256")) &&
    identical(
      readLines(file.path(paths$result, "raw_inputs.before.sha256")),
      readLines(file.path(paths$result, "raw_inputs.after.sha256"))
    ),
  "before/after raw hashes identical"
)
add_check(
  "locked_inputs_unchanged_by_hash",
  file.exists(file.path(paths$result, "locked_inputs.after.sha256")) &&
    file.exists(file.path(paths$result, "locked_inputs.before.sha256")) &&
    identical(
      readLines(file.path(paths$result, "locked_inputs.before.sha256")),
      readLines(file.path(paths$result, "locked_inputs.after.sha256"))
    ),
  "before/after locked-input hashes identical"
)
report <- readLines(
  file.path(project_dir, "reports", "stage_7_singlecell_replication.md"),
  warn = FALSE
)
add_check(
  "report_preserves_exploratory_hierarchy",
  any(grepl("secondary and module-level", report, fixed = TRUE)),
  "primary negative result is not upgraded"
)
add_check(
  "report_stops_before_stage8",
  any(grepl("Stage 8 was not started", report, fixed = TRUE)),
  "stage boundary present"
)

checks <- do.call(rbind, checks)
write_stage7_tsv(checks, file.path(paths$result, "validation_checks.tsv"))
if (any(!checks$passed)) {
  print(checks[!checks$passed, ], row.names = FALSE)
  stop("Stage 7 validation failed")
}
cat("STAGE_7_VALIDATION_OK\t", nrow(checks), "/", nrow(checks), "\n", sep = "")
