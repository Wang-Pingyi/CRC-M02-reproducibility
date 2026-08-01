#!/usr/bin/env Rscript

# Analysis: Stage 5C descriptive annotation audit
# Date: 2026-07-27
# Random seed: not applicable (deterministic tabulation only)
# Unit: cells for annotation inventory; biological tissue/donor retained in
#       all stratified outputs; no inferential testing is performed.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}
result_dir <- file.path(project_root, "results", "05C_annotation")

major <- fread(file.path(
  result_dir, "major_cell_annotations.tsv.gz"
))
epithelial <- fread(file.path(
  result_dir, "epithelial_cell_annotations_final.tsv.gz"
))
run_audit <- fread(file.path(
  result_dir, "copykat_run_audit.tsv"
))
dominance <- fread(file.path(
  result_dir, "epithelial_donor_dominance_audit.tsv"
))

stopifnot(
  nrow(major) == 201649L,
  nrow(epithelial) == 166239L,
  !anyDuplicated(major$cell_id),
  !anyDuplicated(epithelial$cell_id),
  all(epithelial$cell_id %chin% major$cell_id)
)

major_inventory <- major[
  ,
  .(cells = .N),
  by = major_cell_type
][order(-cells)]
major_inventory[, fraction := cells / sum(cells)]
fwrite(
  major_inventory,
  file.path(result_dir, "major_cell_type_inventory.tsv"),
  sep = "\t",
  quote = TRUE
)

epithelial_inventory <- epithelial[
  ,
  .(
    cells = .N,
    donors = uniqueN(donor_id),
    biological_tissues = uniqueN(biological_sample_id)
  ),
  by = epithelial_state
][order(-cells)]
epithelial_inventory[, fraction := cells / sum(cells)]
fwrite(
  epithelial_inventory,
  file.path(result_dir, "epithelial_state_inventory.tsv"),
  sep = "\t",
  quote = TRUE
)

direct <- epithelial[!is.na(copykat_prediction_direct)]
direct_by_stage <- direct[
  ,
  .(cells = .N),
  by = .(lesion_stage, copykat_prediction_direct)
]
direct_by_stage[
  ,
  fraction_of_sampled := cells / sum(cells),
  by = lesion_stage
]
direct_by_stage[
  ,
  defined_denominator := sum(
    cells[copykat_prediction_direct %chin% c("aneuploid", "diploid")]
  ),
  by = lesion_stage
]
direct_by_stage[
  ,
  fraction_among_defined := fifelse(
    copykat_prediction_direct %chin% c("aneuploid", "diploid") &
      defined_denominator > 0,
    cells / defined_denominator,
    NA_real_
  )
]
fwrite(
  direct_by_stage[order(lesion_stage, copykat_prediction_direct)],
  file.path(result_dir, "copykat_direct_predictions_by_stage.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

final_by_stage <- epithelial[
  ,
  .(cells = .N),
  by = .(lesion_stage, malignancy_interpretation)
]
final_by_stage[
  ,
  fraction_within_stage := cells / sum(cells),
  by = lesion_stage
]
fwrite(
  final_by_stage[order(lesion_stage, malignancy_interpretation)],
  file.path(result_dir, "malignancy_interpretation_by_stage.tsv"),
  sep = "\t",
  quote = TRUE
)

final_overall <- epithelial[
  ,
  .(cells = .N),
  by = malignancy_interpretation
][order(-cells)]
final_overall[, fraction := cells / sum(cells)]
fwrite(
  final_overall,
  file.path(result_dir, "malignancy_interpretation_inventory.tsv"),
  sep = "\t",
  quote = TRUE
)

metrics <- data.table(
  metric = c(
    "all_cells",
    "epithelial_cells",
    "epithelial_fraction",
    "major_cell_types",
    "epithelial_states",
    "epithelial_clusters",
    "donors",
    "biological_tissues",
    "copykat_ready_tissues",
    "copykat_completed_tissues",
    "copykat_failed_ready_tissues",
    "copykat_not_evaluable_tissues",
    "copykat_sampled_epithelial_cells",
    "copykat_direct_defined_cells",
    "copykat_direct_defined_fraction",
    "candidate_malignant_cells",
    "candidate_malignant_fraction",
    "donor_dominated_epithelial_clusters"
  ),
  value = c(
    nrow(major),
    nrow(epithelial),
    nrow(epithelial) / nrow(major),
    uniqueN(major$major_cell_type),
    uniqueN(epithelial$epithelial_state),
    uniqueN(epithelial$epithelial_cluster),
    uniqueN(epithelial$donor_id),
    uniqueN(epithelial$biological_sample_id),
    run_audit[status == "ready", .N],
    run_audit[run_status == "completed", .N],
    run_audit[status == "ready" & run_status != "completed", .N],
    run_audit[status != "ready", .N],
    nrow(direct),
    direct[
      copykat_prediction_direct %chin% c("aneuploid", "diploid"),
      .N
    ],
    direct[
      copykat_prediction_direct %chin% c("aneuploid", "diploid"),
      .N
    ] / nrow(direct),
    epithelial[
      malignancy_interpretation == "candidate_malignant_epithelial",
      .N
    ],
    epithelial[
      malignancy_interpretation == "candidate_malignant_epithelial",
      .N
    ] / nrow(epithelial),
    dominance[
      level == "epithelial_cluster" &
        donor_dominance_flag == TRUE,
      .N
    ]
  ),
  note = c(
    "Inventory only",
    "Major annotation Epithelial compartment",
    "Cell-level descriptive fraction; not an abundance estimate",
    "Unique primary major-cell labels",
    "Unique final epithelial-state labels",
    "Unsupervised epithelial clusters audited",
    "Biological donors",
    "Biological tissue samples with retained epithelial cells",
    "Prespecified evaluable tissue runs",
    "Completed CopyKAT tissue runs",
    "Must be zero",
    "Too few epithelial cells; no forced inference",
    "Deterministically sampled query cells",
    "Direct aneuploid or diploid calls",
    "Descriptive technical coverage",
    "CNV-supported candidate label, not proof of malignancy",
    "Cell-level descriptive fraction; not a patient-level endpoint",
    "Flagged at top-donor fraction >= 0.5"
  )
)
fwrite(
  metrics,
  file.path(result_dir, "stage_5C_key_metrics.tsv"),
  sep = "\t",
  quote = TRUE
)

writeLines(
  c(
    "scope=descriptive_annotation_audit",
    "inferential_testing=not_performed",
    "differential_expression=not_performed",
    "trajectory_analysis=not_performed",
    "machine_learning=not_performed",
    "candidate_malignant_interpretation=CNV_supported_candidate_only",
    "cell_level_fractions=not_patient_level_effect_estimates"
  ),
  file.path(result_dir, "stage_5C_descriptive_audit_provenance.txt")
)

message(
  "Stage 5C descriptive audit completed: ",
  nrow(major), " total cells; ",
  nrow(epithelial), " epithelial cells"
)
