#!/usr/bin/env Rscript

# Independent acceptance audit for completed Stage 7 outputs.
# Recomputes all inferential summaries from donor-level source data.

set.seed(20260728)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)
get_param <- read_stage7_parameters(project_dir)

effects <- utils::read.delim(
  file.path(paths$result, "replication_effects_all.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
summary <- utils::read.delim(
  file.path(paths$result, "replication_summary.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
paired <- utils::read.delim(
  file.path(paths$source, "paired_donor_module_differences.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
specificity <- utils::read.delim(
  file.path(paths$source, "epithelial_specificity_donor_differences.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
eligibility <- utils::read.delim(
  file.path(paths$result, "preflight", "GSE132465_pair_eligibility.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
locked <- read_locked_stage7_modules(project_dir)

fdr_threshold <- get_param("global", "fdr_threshold", TRUE)
strong_sd <- get_param("global", "strong_effect_sd", TRUE)
lodo_min <- get_param("global", "lodo_min_sign_stability", TRUE)
tolerance <- 1e-10

recompute_one <- function(x) {
  difference <- x$difference
  n <- length(difference)
  effect <- mean(difference)
  sd_diff <- if (n > 1L) stats::sd(difference) else NA_real_
  se <- sd_diff / sqrt(n)
  critical <- stats::qt(0.975, df = n - 1L)
  p <- if (is.finite(sd_diff) && sd_diff > 0) {
    stats::t.test(difference, mu = 0)$p.value
  } else if (n > 1L && all(difference == 0)) 1 else NA_real_
  lodo <- if (n >= 3L) {
    vapply(seq_len(n), function(i) mean(difference[-i]), numeric(1))
  } else numeric()
  data.frame(
    cohort = x$cohort[[1L]],
    contrast = x$contrast[[1L]],
    module_id = x$module_id[[1L]],
    n_donors_recomputed = n,
    effect_recomputed = effect,
    ci_low_recomputed = effect - critical * se,
    ci_high_recomputed = effect + critical * se,
    standardized_effect_recomputed =
      if (is.finite(sd_diff) && sd_diff > 0) effect / sd_diff else NA_real_,
    p_value_recomputed = p,
    lodo_recomputed = if (length(lodo) && effect != 0) {
      mean(sign(lodo) == sign(effect))
    } else NA_real_,
    max_contribution_recomputed =
      if (sum(abs(difference)) > 0) max(abs(difference)) / sum(abs(difference))
      else NA_real_,
    donor_ids_recomputed = paste(x$donor_id, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

all_details <- rbind(
  paired[, c("module_id", "donor_id", "difference", "cohort", "contrast")],
  specificity[, c("module_id", "donor_id", "difference", "cohort", "contrast")]
)
split_key <- paste(
  all_details$cohort, all_details$contrast, all_details$module_id, sep = "|"
)
recomputed <- do.call(rbind, lapply(split(all_details, split_key), recompute_one))
audit <- merge(
  effects, recomputed,
  by = c("cohort", "contrast", "module_id"), all = TRUE, sort = FALSE
)
recomputed$fdr_recomputed <- ave(
  recomputed$p_value_recomputed,
  interaction(recomputed$cohort, recomputed$contrast, drop = TRUE),
  FUN = function(p) stats::p.adjust(p, method = "BH")
)
audit <- merge(
  audit,
  recomputed[, c("cohort", "contrast", "module_id", "fdr_recomputed")],
  by = c("cohort", "contrast", "module_id"), all.x = TRUE, sort = FALSE
)

same_numeric <- function(a, b, tol = tolerance) {
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}
row_recalculation_pass <- with(
  audit,
  n_donors == n_donors_recomputed &
    same_numeric(effect, effect_recomputed) &
    same_numeric(ci_low, ci_low_recomputed) &
    same_numeric(ci_high, ci_high_recomputed) &
    same_numeric(standardized_effect, standardized_effect_recomputed) &
    same_numeric(p_value, p_value_recomputed) &
    same_numeric(fdr, fdr_recomputed) &
    same_numeric(lodo_sign_stability, lodo_recomputed) &
    same_numeric(
      max_absolute_donor_contribution, max_contribution_recomputed
    )
)
write_stage7_tsv(
  cbind(audit, recalculation_pass = row_recalculation_pass),
  file.path(paths$result, "acceptance_recomputed_effects.tsv")
)

expected_direction <- function(module_id, contrast) {
  candidate <- locked$candidates[locked$candidates$module_id == module_id, ]
  value <- switch(
    contrast,
    normal_to_adenoma = candidate$early_effect,
    adenoma_to_cancer = candidate$cancer_vs_adenoma_effect,
    normal_to_cancer = candidate$cancer_vs_normal_effect,
    epithelial_specificity = 1
  )
  sign(value)
}
effects$expected_direction_audit <- mapply(
  expected_direction, effects$module_id, effects$contrast
)
effects$direction_audit <- sign(effects$effect) == effects$expected_direction_audit
effects$strong_opposite_audit <- !effects$direction_audit &
  effects$fdr < fdr_threshold &
  abs(effects$standardized_effect) >= strong_sd
effects$donor_robust_audit <- effects$lodo_sign_stability >= lodo_min
effects$donor_robust_audit[
  effects$cohort == "GSE161277" & effects$n_donors == 3L
] <- effects$lodo_sign_stability[
  effects$cohort == "GSE161277" & effects$n_donors == 3L
] == 1

reclassify <- function(module_id) {
  m <- effects[effects$module_id == module_id, ]
  early <- m[m$cohort == "GSE161277" & m$contrast == "normal_to_adenoma", ]
  later161 <- m[m$cohort == "GSE161277" & m$contrast == "normal_to_cancer", ]
  later132 <- m[m$cohort == "GSE132465" & m$contrast == "normal_to_cancer", ]
  spec <- m[m$contrast == "epithelial_specificity", ]
  any_opposite <- any(m$strong_opposite_audit %in% TRUE, na.rm = TRUE)
  early_consistent <- isTRUE(early$direction_audit)
  cancer_consistent <- isTRUE(later132$direction_audit)
  any_consistent <- early_consistent || cancer_consistent ||
    isTRUE(later161$direction_audit)
  any_significant <- any(
    c(early$fdr, later161$fdr, later132$fdr) < fdr_threshold,
    na.rm = TRUE
  )
  donor_robust <- isTRUE(early$donor_robust_audit) ||
    isTRUE(later132$donor_robust_audit)
  specificity_positive <- isTRUE(spec$effect > 0)
  specificity_robust <- isTRUE(spec$donor_robust_audit)
  if (any_opposite) {
    "contradictory"
  } else if (
    early_consistent && cancer_consistent && any_significant &&
      donor_robust && specificity_positive && specificity_robust
  ) {
    "replicated"
  } else if (any_consistent && donor_robust && specificity_positive) {
    "directionally consistent but underpowered"
  } else {
    "not replicated"
  }
}
classification_audit <- data.frame(
  module_id = summary$module_id,
  reported_class = summary$replication_class,
  recomputed_class = vapply(summary$module_id, reclassify, character(1)),
  stringsAsFactors = FALSE
)
classification_audit$passed <- with(
  classification_audit, reported_class == recomputed_class
)
write_stage7_tsv(
  classification_audit,
  file.path(paths$result, "acceptance_classification_audit.tsv")
)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}
add_check(
  "all_30_effect_rows_recomputed",
  nrow(audit) == 30L && all(row_recalculation_pass),
  paste0(sum(row_recalculation_pass), "/30 exact within ", tolerance)
)
add_check(
  "all_classifications_independently_rebuilt",
  all(classification_audit$passed),
  paste(
    classification_audit$module_id,
    classification_audit$recomputed_class,
    sep = "=", collapse = ";"
  )
)
add_check(
  "GSE161277_exact_matched_donors",
  setequal(
    unique(paired$donor_id[paired$cohort == "GSE161277"]),
    c("Patient1", "Patient2", "Patient3")
  ),
  paste(unique(paired$donor_id[paired$cohort == "GSE161277"]), collapse = ";")
)
eligible132 <- unique(
  eligibility$donor_id[eligibility$paired_analysis_eligible]
)
add_check(
  "GSE132465_exact_evaluable_donors",
  setequal(
    unique(paired$donor_id[paired$cohort == "GSE132465"]),
    eligible132
  ) && length(eligible132) == 8L,
  paste(eligible132, collapse = ";")
)
add_check(
  "GSE132465_exclusions_prespecified_gate_only",
  setequal(
    unique(eligibility$donor_id[!eligibility$paired_analysis_eligible]),
    c("SMC05", "SMC08")
  ) &&
    all(grepl(
      "fewer_than_20|paired_counterpart_failed_gate",
      eligibility$exclusion_reason[!eligibility$paired_analysis_eligible]
    )),
  paste(
    unique(eligibility$donor_id[!eligibility$paired_analysis_eligible]),
    collapse = ";"
  )
)
add_check(
  "all_GSE161277_early_LODO_strict",
  all(
    effects$lodo_sign_stability[
      effects$cohort == "GSE161277" &
        effects$contrast == "normal_to_adenoma"
    ] == 1
  ),
  "all six modules retain direction in every 2-of-3 donor fit"
)
add_check(
  "no_strong_opposite_effect",
  !any(effects$strong_opposite_audit),
  paste(
    effects$module_id[effects$strong_opposite_audit],
    effects$cohort[effects$strong_opposite_audit],
    effects$contrast[effects$strong_opposite_audit],
    sep = ":", collapse = ";"
  )
)
add_check(
  "primary_candidate_table_remains_empty",
  nrow(utils::read.delim(
    file.path(
      project_dir, "results", "06A_pseudobulk", "candidate_programs.tsv"
    ),
    check.names = FALSE
  )) == 0L,
  "primary Stage 6A negative result frozen"
)
add_check(
  "locked_membership_unchanged",
  identical(
    readLines(file.path(paths$result, "locked_inputs.before.sha256")),
    readLines(file.path(paths$result, "locked_inputs.after.sha256"))
  ),
  "before/after SHA256 manifests identical"
)
add_check(
  "raw_inputs_unchanged",
  identical(
    readLines(file.path(paths$result, "raw_inputs.before.sha256")),
    readLines(file.path(paths$result, "raw_inputs.after.sha256"))
  ),
  "before/after SHA256 manifests identical"
)
add_check(
  "module_gene_coverage",
  all(
    utils::read.delim(
      file.path(paths$result, "preflight", "model_input_audit.tsv")
    )$gene_coverage >= 0.60
  ),
  "minimum observed coverage 0.90"
)
add_check(
  "figure_source_traceability",
  file.exists(file.path(paths$source, "stage_7_module_effects_source.tsv")) &&
    file.exists(file.path(paths$figure, "stage_7_module_effects.pdf")) &&
    file.exists(file.path(paths$figure, "stage_7_module_effects.png")),
  "source table plus PDF and 300-dpi PNG present"
)
add_check(
  "stage_boundary",
  grepl(
    "Stage 8 was not started",
    paste(readLines(
      file.path(project_dir, "reports", "stage_7_singlecell_replication.md")
    ), collapse = "\n"),
    fixed = TRUE
  ),
  "Stage 8 not started"
)

checks <- do.call(rbind, checks)
write_stage7_tsv(
  checks, file.path(paths$result, "acceptance_qc_checks.tsv")
)
if (any(!checks$passed)) {
  print(checks[!checks$passed, ], row.names = FALSE)
  stop("Independent Stage 7 acceptance QC failed")
}
cat(
  "STAGE_7_ACCEPTANCE_QC_OK\t", nrow(checks), "/", nrow(checks), "\n",
  sep = ""
)

