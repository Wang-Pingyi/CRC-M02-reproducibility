#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 09A_independent_acceptance.R <project_dir> <run_id>")
}
project_dir <- normalizePath(args[[1]], mustWork = TRUE)
run_id <- args[[2]]
result_dir <- file.path(project_dir, "results", "09A_stool_feasibility", run_id)
cel_dir <- file.path(
  project_dir, "data_processed", "09A_stool_feasibility", run_id,
  "training_cel"
)
object_path <- file.path(
  project_dir, "objects",
  paste0("GSE99573_9A_training_RMA_", run_id, ".rds")
)

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
checks <- list()
add_check <- function(name, passed, observed, expected, severity = "blocker") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected),
    severity = severity,
    stringsAsFactors = FALSE
  )
}

manifest <- read_tsv(file.path(project_dir, "metadata", "dataset_manifest.tsv"))
manifest <- manifest[manifest$accession == "GSE99573", , drop = FALSE]
split <- read_tsv(file.path(result_dir, "GSE99573_split_audit.tsv"))
inventory <- read_tsv(file.path(result_dir, "GSE99573_sample_inventory.tsv"))
test_audit <- read_tsv(file.path(result_dir, "locked_test_access_audit.tsv"))
membership <- read_tsv(file.path(result_dir, "locked_module_membership.tsv"))
mapping <- read_tsv(file.path(result_dir, "locked_gene_probe_mapping.tsv"))
module <- read_tsv(file.path(result_dir, "locked_module_detectability_summary.tsv"))
qc <- read_tsv(file.path(result_dir, "training_sample_qc_metrics.tsv"))
outliers <- read_tsv(file.path(result_dir, "training_outlier_flags.tsv"))
split_checks <- read_tsv(file.path(result_dir, "split_validation_checks.tsv"))
analysis_checks <- read_tsv(
  file.path(result_dir, "analysis_validation_checks.tsv")
)

expected_split <- c(training = 265L, testing = 65L, not_used = 8L)
observed_split <- table(manifest$validation_split)
add_check(
  "manifest_rows_338", nrow(manifest) == 338L, nrow(manifest), 338L
)
add_check(
  "manifest_sample_id_unique",
  !anyDuplicated(manifest$sample_id),
  length(unique(manifest$sample_id)),
  338L
)
add_check(
  "manifest_donor_id_complete_unique",
  all(!is.na(manifest$donor_id) & nzchar(manifest$donor_id)) &&
    !anyDuplicated(manifest$donor_id),
  paste(
    sum(!is.na(manifest$donor_id) & nzchar(manifest$donor_id)),
    length(unique(manifest$donor_id)),
    sep = "/"
  ),
  "338/338"
)
for (partition in names(expected_split)) {
  add_check(
    paste0("manifest_split_", partition),
    observed_split[[partition]] == expected_split[[partition]],
    observed_split[[partition]],
    expected_split[[partition]]
  )
}
split_recomputed <- stats::aggregate(
  sample_id ~ validation_split + condition,
  data = manifest,
  FUN = length
)
names(split_recomputed)[[3L]] <- "n_samples"
split_compare <- merge(
  split_recomputed[, c("validation_split", "condition", "n_samples")],
  split[, c("validation_split", "condition", "n_samples")],
  by = c("validation_split", "condition"),
  all = TRUE,
  suffixes = c("_manifest", "_reported")
)
add_check(
  "reported_split_matches_manifest",
  nrow(split_compare) == 7L &&
    all(split_compare$n_samples_manifest == split_compare$n_samples_reported),
  paste(
    if (nrow(split_compare)) {
      sum(split_compare$n_samples_manifest == split_compare$n_samples_reported)
    } else {
      0L
    },
    nrow(split_compare),
    sep = "/"
  ),
  "7/7"
)

expected_test_ids <- sort(
  manifest$sample_id[manifest$validation_split == "testing"]
)
add_check(
  "test_audit_exactly_65_locked_ids",
  nrow(test_audit) == 65L &&
    identical(sort(test_audit$sample_id), expected_test_ids),
  paste(nrow(test_audit), length(unique(test_audit$sample_id)), sep = "/"),
  "65/65"
)
add_check(
  "test_audit_no_expression_access",
  all(!test_audit$cel_extracted) && all(!test_audit$expression_accessed),
  paste(
    sum(test_audit$cel_extracted),
    sum(test_audit$expression_accessed),
    sep = "/"
  ),
  "0/0"
)

cel_files <- list.files(
  cel_dir, pattern = "\\.CEL$", full.names = FALSE, ignore.case = TRUE
)
cel_gz_files <- list.files(
  cel_dir, pattern = "\\.CEL\\.gz$", full.names = FALSE, ignore.case = TRUE
)
cel_ids <- sub("_.*$", "", cel_files)
expected_training_ids <- sort(
  manifest$sample_id[manifest$validation_split == "training"]
)
add_check(
  "working_cel_exactly_training",
  length(cel_files) == 265L &&
    identical(sort(cel_ids), expected_training_ids),
  paste(length(cel_files), length(unique(cel_ids)), sep = "/"),
  "265/265"
)
add_check(
  "working_cel_gz_exactly_265",
  length(cel_gz_files) == 265L,
  length(cel_gz_files),
  265L
)
add_check(
  "test_ids_absent_from_working_directory",
  !any(expected_test_ids %in% cel_ids),
  sum(expected_test_ids %in% cel_ids),
  0L
)
not_used_ids <- manifest$sample_id[manifest$validation_split == "not_used"]
add_check(
  "not_used_ids_absent_from_working_directory",
  !any(not_used_ids %in% cel_ids),
  sum(not_used_ids %in% cel_ids),
  0L
)

before_hash <- readLines(
  file.path(result_dir, "raw_inputs.before.sha256"), warn = FALSE
)
after_hash <- readLines(
  file.path(result_dir, "raw_inputs.after.sha256"), warn = FALSE
)
add_check(
  "raw_hashes_unchanged",
  identical(before_hash, after_hash),
  paste(length(before_hash), length(after_hash), sep = "/"),
  "2/2"
)
add_check(
  "server_split_checks_all_pass",
  nrow(split_checks) == 10L && all(split_checks$passed),
  paste(sum(split_checks$passed), nrow(split_checks), sep = "/"),
  "10/10"
)
add_check(
  "server_analysis_checks_all_pass",
  nrow(analysis_checks) == 10L && all(analysis_checks$passed),
  paste(sum(analysis_checks$passed), nrow(analysis_checks), sep = "/"),
  "10/10"
)

object <- readRDS(object_path)
expression <- object$expression
add_check(
  "RMA_object_dimensions",
  identical(dim(expression), c(70523L, 265L)),
  paste(dim(expression), collapse = "x"),
  "70523x265"
)
add_check(
  "RMA_object_all_finite",
  all(is.finite(expression)),
  sum(is.finite(expression)),
  length(expression)
)
add_check(
  "RMA_object_training_ids_exact",
  identical(sort(colnames(expression)), expected_training_ids),
  length(unique(colnames(expression))),
  265L
)
add_check(
  "RMA_object_test_firewall_flag",
  identical(object$test_expression_accessed, FALSE),
  object$test_expression_accessed,
  FALSE
)
add_check(
  "RMA_object_plausible_log2_range",
  min(expression) > -1 && max(expression) < 25,
  sprintf("%.4f to %.4f", min(expression), max(expression)),
  "-1 < range < 25",
  "caveat"
)

add_check(
  "locked_membership_747_rows",
  nrow(membership) == 747L,
  nrow(membership),
  747L
)
add_check(
  "locked_membership_six_modules",
  length(unique(membership$module_id)) == 6L,
  length(unique(membership$module_id)),
  6L
)
add_check(
  "locked_mapping_632_unique_genes",
  nrow(mapping) == 632L && !anyDuplicated(mapping$gene),
  paste(nrow(mapping), length(unique(mapping$gene)), sep = "/"),
  "632/632"
)
add_check(
  "mapping_status_partition",
  sum(mapping$mapping_status == "mapped") == 567L &&
    sum(mapping$mapping_status == "missing_probe") == 65L,
  paste(
    sum(mapping$mapping_status == "mapped"),
    sum(mapping$mapping_status == "missing_probe"),
    sep = "/"
  ),
  "567/65"
)
add_check(
  "detectability_status_partition",
  sum(mapping$detectability_status == "detectable") == 521L &&
    sum(mapping$detectability_status == "low_expression") == 46L &&
    sum(mapping$detectability_status == "missing_probe") == 65L,
  paste(
    sum(mapping$detectability_status == "detectable"),
    sum(mapping$detectability_status == "low_expression"),
    sum(mapping$detectability_status == "missing_probe"),
    sep = "/"
  ),
  "521/46/65"
)
add_check(
  "module_summary_six_rows",
  nrow(module) == 6L && !anyDuplicated(module$module_id),
  nrow(module),
  6L
)
add_check(
  "module_accounting_complete",
  all(
    module$mapped_genes + module$missing_probe_genes == module$locked_genes
  ) &&
    all(
      module$detectable_genes + module$low_expression_genes ==
        module$mapped_genes
    ),
  "all module partitions reconcile",
  "TRUE"
)

add_check(
  "QC_265_unique_training_samples",
  nrow(qc) == 265L && !anyDuplicated(qc$sample_id) &&
    identical(sort(qc$sample_id), expected_training_ids),
  paste(nrow(qc), length(unique(qc$sample_id)), sep = "/"),
  "265/265"
)
add_check(
  "QC_no_outcome_labels_or_auto_exclusion",
  !("condition" %in% names(qc)) &&
    !any(qc$outcome_labels_used) &&
    !any(qc$automatic_exclusion),
  paste(
    "condition_column", "condition" %in% names(qc),
    "outcome_used", sum(qc$outcome_labels_used),
    "excluded", sum(qc$automatic_exclusion)
  ),
  "FALSE/0/0"
)
add_check(
  "QC_flag_table_matches",
  nrow(outliers) == sum(qc$qc_review_flag) &&
    identical(
      sort(outliers$sample_id),
      sort(qc$sample_id[qc$qc_review_flag])
    ),
  paste(nrow(outliers), sum(qc$qc_review_flag), sep = "/"),
  "9/9"
)
add_check(
  "QC_metrics_finite",
  all(vapply(
    qc[vapply(qc, is.numeric, logical(1))],
    function(x) all(is.finite(x)),
    logical(1)
  )),
  "all numeric QC columns finite",
  "TRUE"
)

figure_path <- file.path(
  project_dir, "figures", "09A_stool_feasibility", run_id,
  "training_array_QC_outcome_blind.pdf"
)
add_check(
  "QC_figure_present",
  file.exists(figure_path) && file.info(figure_path)$size > 1000,
  if (file.exists(figure_path)) file.info(figure_path)$size else 0,
  ">1000 bytes"
)
report_path <- file.path(
  project_dir, "reports", "stage_9A_stool_feasibility.md"
)
report <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
required_report_terms <- c(
  "Training: 265", "Testing: 65", "Not Used: 8",
  "Mapped unique genes: 567/632",
  "Provisionally detectable unique genes: 521/632",
  "Outcome-blind QC review flags: 9",
  "no model", "Stop after Stage 9A"
)
add_check(
  "report_key_claims_present",
  all(vapply(
    required_report_terms,
    function(term) grepl(term, report, fixed = TRUE),
    logical(1)
  )),
  sum(vapply(
    required_report_terms,
    function(term) grepl(term, report, fixed = TRUE),
    logical(1)
  )),
  length(required_report_terms)
)

checks_df <- do.call(rbind, checks)
write_tsv(
  checks_df,
  file.path(result_dir, "independent_acceptance_checks.tsv")
)
blockers <- checks_df$severity == "blocker" & !checks_df$passed
cat(
  "STAGE9A_INDEPENDENT_ACCEPTANCE ",
  if (any(blockers)) "FAIL" else "PASS",
  " ",
  sum(checks_df$passed),
  "/",
  nrow(checks_df),
  "\n",
  sep = ""
)
if (any(blockers)) {
  print(checks_df[blockers, , drop = FALSE])
  quit(status = 2L)
}
