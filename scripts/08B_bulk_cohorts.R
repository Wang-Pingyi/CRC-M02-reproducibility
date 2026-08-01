#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08B_bulk_cohorts.R PROJECT_DIR RUN_ID")
project <- normalizePath(args[1], mustWork = TRUE)
run_id <- args[2]
source(file.path(project, "scripts", "08B_helpers.R"))
set.seed(20260729)
paths <- stage8b_paths(project, run_id)

cohorts <- c("GSE41657", "GSE100179", "GSE8671")
score_sets <- list()
coverage_sets <- list()
for (cohort in cohorts) {
  x <- stage8b_module_scores(
    file.path(paths$input_run, paste0(cohort, "_normalized_probe_matrix.rds")),
    file.path(paths$input_results, cohort, "locked_module_gene_mapping.tsv"),
    file.path(paths$input_results, cohort, "sample_metadata.tsv"),
    paths, cohort
  )
  score_sets[[cohort]] <- x$scores
  coverage_sets[[cohort]] <- x$coverage
}
all_scores <- do.call(rbind, score_sets)
all_coverage <- do.call(rbind, coverage_sets)

# Freeze cohort stage labels before fitting.
all_scores$analysis_group <- NA_character_
all_scores$stage_score <- NA_real_
i <- all_scores$accession == "GSE41657"
all_scores$analysis_group[i & all_scores$condition == "normal_mucosa"] <- "normal"
all_scores$analysis_group[i & all_scores$histology == "Low-grade dysplasia"] <- "low_grade_adenoma"
all_scores$analysis_group[i & all_scores$histology == "High-grade dysplasia"] <- "high_grade_adenoma"
all_scores$analysis_group[i & all_scores$condition == "cancer"] <- "cancer"
all_scores$stage_score[i] <- c(normal = 0, low_grade_adenoma = 1, high_grade_adenoma = 2, cancer = 3)[all_scores$analysis_group[i]]

i <- all_scores$accession == "GSE100179"
all_scores$analysis_group[i & all_scores$condition == "normal_mucosa"] <- "normal"
all_scores$analysis_group[i & all_scores$condition == "adenoma"] <- "adenoma"
all_scores$analysis_group[i & all_scores$condition == "cancer"] <- "cancer"
all_scores$stage_score[i] <- c(normal = 0, adenoma = 1, cancer = 2)[all_scores$analysis_group[i]]

i <- all_scores$accession == "GSE8671"
all_scores$analysis_group[i & all_scores$condition == "normal_mucosa"] <- "normal"
all_scores$analysis_group[i & all_scores$condition == "adenoma"] <- "adenoma"
all_scores$stage_score[i] <- c(normal = 0, adenoma = 1)[all_scores$analysis_group[i]]
if (anyNA(all_scores$analysis_group) || anyNA(all_scores$stage_score)) stop("Unmapped Stage 8B group labels")

qc <- do.call(rbind, lapply(cohorts, function(cohort) {
  read.delim(file.path(paths$input_results, cohort, "qc_sample_metrics.tsv"), check.names = FALSE)
}))
all_scores <- merge(all_scores, qc[, c("accession", "sample_id", "qc_flag_review")],
                    by = c("accession", "sample_id"), all.x = TRUE, suffixes = c("", "_from_qc"))
if ("qc_flag_review_from_qc" %in% names(all_scores)) all_scores$qc_flag_review <- all_scores$qc_flag_review_from_qc

stage8b_write_tsv(all_scores, file.path(paths$result, "bulk_module_scores.tsv"))
stage8b_write_tsv(all_coverage, file.path(paths$result, "module_mapping_coverage.tsv"))

high_leverage <- list(
  GSE41657 = c("GSM1021168", "GSM1021177"),
  GSE100179 = character(),
  GSE8671 = "GSM215093"
)

filter_set <- function(x, set_name) {
  if (set_name == "primary_all_samples") return(x)
  remove <- if (set_name == "exclude_high_leverage") x$sample_id %in% high_leverage[[unique(x$accession)]] else x$qc_flag_review %in% c(TRUE, "TRUE")
  if (unique(x$accession) == "GSE8671" && any(remove)) {
    remove_donors <- unique(x$donor_id[remove])
    remove <- x$donor_id %in% remove_donors
  }
  x[!remove, , drop = FALSE]
}

fit_cohort <- function(x, analysis_set) {
  accession <- unique(x$accession)
  output <- list()
  for (module in unique(x$module_id)) {
    z <- filter_set(x[x$module_id == module, ], analysis_set)
    if (accession == "GSE41657") {
      z$adenoma_binary <- ifelse(grepl("adenoma", z$analysis_group), "adenoma", z$analysis_group)
      output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, "ordered_trend", analysis_set,
        stage8b_fit_trend(z), "HC3_robust_ordinal")
      output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, "adenoma_vs_normal", analysis_set,
        stage8b_fit_binary(transform(z, analysis_group = adenoma_binary), "normal", "adenoma"), "HC3_robust_binary")
      output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, "cancer_vs_adenoma", analysis_set,
        stage8b_fit_binary(transform(z, analysis_group = adenoma_binary), "adenoma", "cancer"), "HC3_robust_binary")
      output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, "cancer_vs_normal", analysis_set,
        stage8b_fit_binary(z, "normal", "cancer"), "HC3_robust_binary")
      for (pair in list(c("normal", "low_grade_adenoma"), c("low_grade_adenoma", "high_grade_adenoma"), c("high_grade_adenoma", "cancer"))) {
        endpoint <- paste0(pair[2], "_vs_", pair[1])
        output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, endpoint, analysis_set,
          stage8b_fit_binary(z, pair[1], pair[2]), "HC3_robust_adjacent")
      }
    } else if (accession == "GSE100179") {
      output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, "ordered_trend", analysis_set,
        stage8b_fit_trend(z), "HC3_robust_ordinal")
      for (pair in list(c("normal", "adenoma"), c("adenoma", "cancer"), c("normal", "cancer"))) {
        endpoint <- paste0(pair[2], "_vs_", pair[1])
        output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, endpoint, analysis_set,
          stage8b_fit_binary(z, pair[1], pair[2]), "HC3_robust_binary")
      }
    } else {
      output[[length(output) + 1L]] <- stage8b_effect_row(accession, module, "adenoma_vs_normal", analysis_set,
        stage8b_fit_paired(z, "normal", "adenoma"), "verified_pair_difference",
        "One difference per verified donor pair")
    }
  }
  do.call(rbind, output)
}

effects <- do.call(rbind, unlist(lapply(cohorts, function(cohort) {
  x <- all_scores[all_scores$accession == cohort, ]
  lapply(c("primary_all_samples", "exclude_high_leverage", "exclude_all_qc_flags"),
         function(set_name) fit_cohort(x, set_name))
}), recursive = FALSE))
effects <- stage8b_add_fdr(effects)
stage8b_write_tsv(effects, file.path(paths$result, "bulk_cohort_effects.tsv"))

means <- aggregate(module_score ~ accession + module_id + analysis_group, all_scores,
                   function(x) c(n = length(x), mean = mean(x), sd = sd(x)))
means <- data.frame(means[c("accession", "module_id", "analysis_group")],
                    n = means$module_score[, "n"], mean = means$module_score[, "mean"],
                    sd = means$module_score[, "sd"])
stage8b_write_tsv(means, file.path(paths$result, "bulk_stage_means.tsv"))

# Leave-one-pair-out audit for GSE8671.
g <- all_scores[all_scores$accession == "GSE8671", ]
lodo <- list()
for (module in unique(g$module_id)) {
  z <- g[g$module_id == module, ]
  for (donor in unique(z$donor_id)) {
    fit <- stage8b_fit_paired(z[z$donor_id != donor, ], "normal", "adenoma")
    lodo[[length(lodo) + 1L]] <- data.frame(module_id = module, omitted_donor = donor,
      effect = fit["effect"], ci_low = fit["ci_low"], ci_high = fit["ci_high"],
      p_value = fit["p_value"], n_pairs = fit["n"])
  }
}
lodo <- do.call(rbind, lodo)
stage8b_write_tsv(lodo, file.path(paths$result, "GSE8671_leave_one_pair_out.tsv"))
lodo_summary <- do.call(rbind, lapply(split(lodo, lodo$module_id), function(z) {
  data.frame(module_id = z$module_id[1], n_fits = nrow(z), effect_min = min(z$effect),
    effect_max = max(z$effect), sign_stability = max(mean(z$effect > 0), mean(z$effect < 0)))
}))
stage8b_write_tsv(lodo_summary, file.path(paths$result, "GSE8671_leave_one_pair_out_summary.tsv"))

saveRDS(list(scores = all_scores, effects = effects, coverage = all_coverage),
        file.path(paths$result, "bulk_stage8B_intermediate.rds"))
