#!/usr/bin/env Rscript

# Analysis: Stage 6A exploratory amendment
# Components: paired FAP gene-level sensitivity and stage-blind modules
# Date: 2026-07-28
# Random seed: 20260728

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else getwd()

required_packages <- c("Matrix", "edgeR", "limma", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

input_path <- file.path(
  project_dir, "objects", "GSE201348_6A_epithelial_pseudobulk.rds"
)
primary_result_path <- file.path(
  project_dir, "results", "06A_pseudobulk", "pseudobulk_results.tsv"
)
parameter_path <- file.path(
  project_dir, "config", "06A_amendment_parameters.tsv"
)
stopifnot(file.exists(input_path), file.exists(primary_result_path), file.exists(parameter_path))

result_dir <- file.path(project_dir, "results", "06A_amendment")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_dir, "figures", "06A_amendment")
report_dir <- file.path(project_dir, "reports")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

parameters <- utils::read.delim(parameter_path, check.names = FALSE)
param <- setNames(parameters$value, parameters$parameter)
p_num <- function(name) as.numeric(param[[name]])

write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

pb <- readRDS(input_path)
counts <- pb$counts
meta <- pb$metadata
stopifnot(identical(colnames(counts), rownames(meta)))
meta$stage <- factor(meta$stage, levels = c("normal", "adenoma", "cancer"))

primary_results <- utils::read.delim(primary_result_path, check.names = FALSE)
primary_early <- primary_results[
  primary_results$contrast == "adenoma_vs_normal",
  c("gene", "epithelial_state", "log2FC", "FDR"),
  drop = FALSE
]
colnames(primary_early)[3:4] <- c("primary_log2FC", "primary_FDR")

make_result_table <- function(fit, coefficient, feature_name, state_name, analysis_name) {
  tt <- limma::topTable(
    fit, coef = coefficient, number = Inf, sort.by = "none", adjust.method = "BH"
  )
  se <- fit$stdev.unscaled[, coefficient] * fit$sigma
  df <- fit$df.total
  critical <- stats::qt(0.975, df = df)
  data.frame(
    feature = rownames(tt),
    epithelial_state = state_name,
    analysis = analysis_name,
    effect = tt$logFC,
    CI95_low = tt$logFC - critical * se,
    CI95_high = tt$logFC + critical * se,
    p_value = tt$P.Value,
    FDR = tt$adj.P.Val,
    standard_error = se,
    df_total = df,
    feature_type = feature_name,
    stringsAsFactors = FALSE
  )
}

safe_ebayes <- function(fit, context, robust_preferred = TRUE) {
  if (!robust_preferred) {
    return(limma::eBayes(fit, robust = FALSE))
  }
  tryCatch(
    limma::eBayes(fit, robust = TRUE),
    error = function(e) {
      warning(
        "Robust eBayes failed; using standard eBayes. context=", context,
        "; coefficients=", paste(dim(fit$coefficients), collapse = "x"),
        "; stdev.unscaled=", paste(dim(fit$stdev.unscaled), collapse = "x"),
        "; reason=", conditionMessage(e)
      )
      limma::eBayes(fit, robust = FALSE)
    }
  )
}

fit_paired_gene_state <- function(state_name, omit_donor = NA_character_) {
  state_meta <- meta[
    meta$eligible &
      meta$epithelial_state == state_name &
      meta$fap_binary == "FAP" &
      meta$stage %in% c("normal", "adenoma"),
    ,
    drop = FALSE
  ]
  if (!is.na(omit_donor)) {
    state_meta <- state_meta[state_meta$donor_id != omit_donor, , drop = FALSE]
  }
  normal_donors <- unique(state_meta$donor_id[state_meta$stage == "normal"])
  adenoma_donors <- unique(state_meta$donor_id[state_meta$stage == "adenoma"])
  paired_donors <- sort(intersect(normal_donors, adenoma_donors))
  state_meta <- state_meta[state_meta$donor_id %in% paired_donors, , drop = FALSE]
  state_meta <- state_meta[order(state_meta$donor_id, state_meta$stage), , drop = FALSE]

  if (length(paired_donors) < p_num("min_paired_FAP_donors")) {
    return(list(status = "insufficient_pairs", donors = paired_donors, metadata = state_meta))
  }

  state_counts <- counts[, state_meta$pseudobulk_id, drop = FALSE]
  dge <- edgeR::DGEList(state_counts)
  keep <- edgeR::filterByExpr(
    dge,
    group = state_meta$stage,
    min.count = p_num("gene_min_count"),
    min.total.count = p_num("gene_min_total_count")
  )
  dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  state_meta$donor_id <- factor(state_meta$donor_id)
  state_meta$stage <- factor(state_meta$stage, levels = c("normal", "adenoma"))
  design <- stats::model.matrix(~ donor_id + stage, data = state_meta)
  voom_object <- limma::voom(dge, design = design, plot = FALSE)
  fit <- safe_ebayes(
    limma::lmFit(voom_object, design),
    context = paste0(
      "paired_gene:", state_name, ":omit=",
      ifelse(is.na(omit_donor), "none", omit_donor)
    )
  )
  coefficient <- match("stageadenoma", colnames(fit$coefficients))
  result <- make_result_table(
    fit, coefficient, "gene", state_name, "paired_FAP_adenoma_vs_normal"
  )
  result$gene <- result$feature
  result$n_paired_donors <- length(paired_donors)
  result$paired_donors <- paste(paired_donors, collapse = ";")
  result$omitted_donor <- ifelse(is.na(omit_donor), "", omit_donor)
  list(
    status = "ok",
    donors = paired_donors,
    metadata = state_meta,
    result = result,
    voom = voom_object
  )
}

states <- sort(unique(meta$epithelial_state))
paired_audit_rows <- list()
paired_fit_list <- list()
paired_results_list <- list()

for (state_name in states) {
  cat("Paired FAP fit:", state_name, "\n")
  fitted <- fit_paired_gene_state(state_name)
  paired_fit_list[[state_name]] <- fitted
  paired_audit_rows[[state_name]] <- data.frame(
    epithelial_state = state_name,
    status = fitted$status,
    n_paired_donors = length(fitted$donors),
    paired_donors = paste(fitted$donors, collapse = ";"),
    stringsAsFactors = FALSE
  )
  if (fitted$status == "ok") {
    paired_results_list[[state_name]] <- fitted$result
  }
}

paired_audit <- do.call(rbind, paired_audit_rows)
write_tsv(paired_audit, file.path(result_dir, "paired_FAP_state_audit.tsv"))

empty_paired <- data.frame(
  feature = character(), epithelial_state = character(), analysis = character(),
  effect = numeric(), CI95_low = numeric(), CI95_high = numeric(),
  p_value = numeric(), FDR = numeric(), standard_error = numeric(),
  df_total = numeric(), feature_type = character(), gene = character(),
  n_paired_donors = integer(), paired_donors = character(),
  omitted_donor = character(), primary_log2FC = numeric(),
  primary_FDR = numeric(), direction_matches_primary = logical(),
  stringsAsFactors = FALSE
)
paired_results <- if (length(paired_results_list)) {
  do.call(rbind, paired_results_list)
} else {
  empty_paired
}
if (nrow(paired_results)) {
  paired_results <- merge(
    paired_results,
    primary_early,
    by = c("gene", "epithelial_state"),
    all.x = TRUE
  )
  paired_results$direction_matches_primary <- with(
    paired_results,
    is.na(primary_log2FC) | sign(effect) == sign(primary_log2FC)
  )
  paired_results <- paired_results[
    order(paired_results$epithelial_state, paired_results$FDR),
  ]
}
write_tsv(paired_results, file.path(result_dir, "paired_FAP_gene_results.tsv"))

paired_hits <- paired_results[
  paired_results$FDR <= p_num("exploratory_FDR_threshold") &
    abs(paired_results$effect) >= p_num("paired_gene_min_abs_log2fc"),
  ,
  drop = FALSE
]
write_tsv(paired_hits, file.path(result_dir, "paired_FAP_gene_hits.tsv"))

paired_lodo_rows <- list()
if (nrow(paired_hits)) {
  for (state_name in unique(paired_hits$epithelial_state)) {
    state_hits <- paired_hits[paired_hits$epithelial_state == state_name, ]
    for (omitted in paired_fit_list[[state_name]]$donors) {
      lodo_fit <- fit_paired_gene_state(state_name, omit_donor = omitted)
      if (lodo_fit$status != "ok") next
      selected <- lodo_fit$result[
        lodo_fit$result$gene %in% state_hits$gene,
        c("gene", "epithelial_state", "effect", "omitted_donor"),
        drop = FALSE
      ]
      paired_lodo_rows[[length(paired_lodo_rows) + 1L]] <- selected
    }
  }
}
paired_lodo <- if (length(paired_lodo_rows)) {
  do.call(rbind, paired_lodo_rows)
} else {
  data.frame(
    gene = character(), epithelial_state = character(),
    effect = numeric(), omitted_donor = character()
  )
}
write_tsv(paired_lodo, file.path(result_dir, "paired_FAP_gene_LODO.tsv"))

paired_lodo_summary <- data.frame(
  gene = character(), epithelial_state = character(), full_effect = numeric(),
  n_evaluable_omissions = integer(), sign_stability = numeric(),
  median_LODO_effect = numeric(), passes_LODO = logical()
)
if (nrow(paired_hits)) {
  paired_lodo_summary <- do.call(rbind, lapply(seq_len(nrow(paired_hits)), function(i) {
    hit <- paired_hits[i, ]
    values <- paired_lodo$effect[
      paired_lodo$gene == hit$gene &
        paired_lodo$epithelial_state == hit$epithelial_state
    ]
    data.frame(
      gene = hit$gene,
      epithelial_state = hit$epithelial_state,
      full_effect = hit$effect,
      n_evaluable_omissions = length(values),
      sign_stability = if (length(values)) mean(sign(values) == sign(hit$effect)) else NA_real_,
      median_LODO_effect = if (length(values)) stats::median(values) else NA_real_,
      passes_LODO = if (length(values)) {
        mean(sign(values) == sign(hit$effect)) >= p_num("lodo_sign_stability_threshold")
      } else {
        FALSE
      },
      stringsAsFactors = FALSE
    )
  }))
}
write_tsv(
  paired_lodo_summary,
  file.path(result_dir, "paired_FAP_gene_LODO_summary.tsv")
)

fit_module_scores <- function(
  score_matrix,
  score_meta,
  state_name,
  paired_only = FALSE,
  omitted_donor = NA_character_
) {
  if (!is.na(omitted_donor)) {
    keep <- score_meta$donor_id != omitted_donor
    score_meta <- score_meta[keep, , drop = FALSE]
    score_matrix <- score_matrix[, keep, drop = FALSE]
  }

  if (paired_only) {
    score_meta <- score_meta[
      score_meta$fap_binary == "FAP" &
        score_meta$stage %in% c("normal", "adenoma"),
      ,
      drop = FALSE
    ]
    normal_donors <- unique(score_meta$donor_id[score_meta$stage == "normal"])
    adenoma_donors <- unique(score_meta$donor_id[score_meta$stage == "adenoma"])
    paired_donors <- sort(intersect(normal_donors, adenoma_donors))
    score_meta <- score_meta[score_meta$donor_id %in% paired_donors, , drop = FALSE]
    score_matrix <- score_matrix[, score_meta$pseudobulk_id, drop = FALSE]
    if (length(paired_donors) < p_num("min_paired_FAP_donors")) {
      return(list(status = "insufficient_pairs", results = NULL, donors = paired_donors))
    }
    score_meta$donor_id <- factor(score_meta$donor_id)
    score_meta$stage <- factor(score_meta$stage, levels = c("normal", "adenoma"))
    design <- stats::model.matrix(~ donor_id + stage, data = score_meta)
    # Robust empirical Bayes needs more than one feature. Candidate-restricted
    # LODO fits can legitimately contain a single module, so use standard
    # eBayes for that edge case instead of dropping the sensitivity estimate.
    fit <- safe_ebayes(
      limma::lmFit(score_matrix, design),
      context = paste0(
        "paired_module:", state_name, ":omit=",
        ifelse(is.na(omitted_donor), "none", omitted_donor)
      ),
      robust_preferred = nrow(score_matrix) >= 2L
    )
    coefficient <- match("stageadenoma", colnames(fit$coefficients))
    result <- make_result_table(
      fit, coefficient, "module", state_name, "paired_FAP_module_adenoma_vs_normal"
    )
    result$contrast <- "adenoma_vs_normal"
    result$n_normal_donors <- length(paired_donors)
    result$n_adenoma_donors <- length(paired_donors)
    result$n_cancer_donors <- 0L
    result$donor_correlation <- NA_real_
    result$model_formula <- "~ donor_id + stage"
    return(list(status = "ok", results = result, donors = paired_donors))
  }

  score_meta$stage <- droplevels(score_meta$stage)
  donor_counts <- table(score_meta$stage, score_meta$donor_id) > 0
  donors_per_stage <- rowSums(donor_counts)
  if (!all(c("normal", "adenoma") %in% names(donors_per_stage)) ||
      any(donors_per_stage[c("normal", "adenoma")] < 3L)) {
    return(list(
      status = "insufficient_normal_adenoma_donors",
      results = NULL,
      donors = donors_per_stage
    ))
  }

  has_cancer <- "cancer" %in% names(donors_per_stage) &&
    donors_per_stage[["cancer"]] >= 3L
  keep_stages <- if (has_cancer) c("normal", "adenoma", "cancer") else c("normal", "adenoma")
  keep <- score_meta$stage %in% keep_stages
  score_meta <- score_meta[keep, , drop = FALSE]
  score_matrix <- score_matrix[, keep, drop = FALSE]
  score_meta$stage <- factor(score_meta$stage, levels = keep_stages)
  score_meta$fap_binary <- factor(score_meta$fap_binary, levels = c("nonFAP", "FAP"))

  design_formula <- ~ 0 + stage + fap_binary
  design <- stats::model.matrix(design_formula, data = score_meta)
  if (qr(design)$rank < ncol(design)) {
    design_formula <- ~ 0 + stage
    design <- stats::model.matrix(design_formula, data = score_meta)
  }
  colnames(design) <- sub("^stage", "", colnames(design))
  correlation <- 0
  if (nrow(score_matrix) >= 2L) {
    duplicate <- tryCatch(
      limma::duplicateCorrelation(
        score_matrix, design = design, block = score_meta$donor_id
      ),
      error = function(e) NULL
    )
    if (!is.null(duplicate) && is.finite(duplicate$consensus.correlation)) {
      correlation <- duplicate$consensus.correlation
    }
  }
  fit <- limma::lmFit(
    score_matrix,
    design = design,
    block = score_meta$donor_id,
    correlation = correlation
  )
  contrast_matrix <- if (has_cancer) {
    limma::makeContrasts(
      adenoma_vs_normal = adenoma - normal,
      cancer_vs_adenoma = cancer - adenoma,
      cancer_vs_normal = cancer - normal,
      levels = design
    )
  } else {
    limma::makeContrasts(
      adenoma_vs_normal = adenoma - normal,
      levels = design
    )
  }
  fit <- safe_ebayes(
    limma::contrasts.fit(fit, contrast_matrix),
    context = paste0(
      "all_cohort_module:", state_name, ":omit=",
      ifelse(is.na(omitted_donor), "none", omitted_donor)
    ),
    robust_preferred = nrow(score_matrix) >= 2L
  )
  result_rows <- lapply(colnames(contrast_matrix), function(contrast_name) {
    coefficient <- match(contrast_name, colnames(fit$coefficients))
    result <- make_result_table(
      fit, coefficient, "module", state_name, "stage_blind_module_all_cohort"
    )
    result$contrast <- contrast_name
    result$n_normal_donors <- unname(donors_per_stage["normal"])
    result$n_adenoma_donors <- unname(donors_per_stage["adenoma"])
    result$n_cancer_donors <- if (has_cancer) unname(donors_per_stage["cancer"]) else 0L
    result$donor_correlation <- correlation
    result$model_formula <- paste(deparse(design_formula), collapse = "")
    result
  })
  list(
    status = "ok",
    results = do.call(rbind, result_rows),
    donors = donors_per_stage
  )
}

module_membership_rows <- list()
module_score_rows <- list()
module_result_rows <- list()
paired_module_result_rows <- list()
module_audit_rows <- list()
module_score_objects <- list()

for (state_name in states) {
  cat("Stage-blind modules:", state_name, "\n")
  state_meta <- meta[
    meta$eligible & meta$epithelial_state == state_name,
    ,
    drop = FALSE
  ]
  if (nrow(state_meta) < 6L) {
    module_audit_rows[[state_name]] <- data.frame(
      epithelial_state = state_name,
      status = "insufficient_pseudobulks",
      n_pseudobulks = nrow(state_meta),
      genes_after_expression_filter = 0L,
      genes_clustered = 0L,
      modules_raw = 0L,
      modules_testable = 0L,
      stage_labels_used_for_construction = FALSE
    )
    next
  }

  state_counts <- counts[, state_meta$pseudobulk_id, drop = FALSE]
  dge <- edgeR::DGEList(state_counts)
  keep <- edgeR::filterByExpr(
    dge,
    group = rep(1L, nrow(state_meta)),
    min.count = p_num("gene_min_count"),
    min.total.count = p_num("gene_min_total_count")
  )
  dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  variability <- apply(logcpm, 1, stats::mad)
  variability <- variability[is.finite(variability) & variability > 0]
  selected_genes <- names(sort(variability, decreasing = TRUE))
  selected_genes <- head(
    selected_genes,
    min(length(selected_genes), as.integer(p_num("module_top_variable_genes")))
  )
  if (length(selected_genes) < 20L) {
    module_audit_rows[[state_name]] <- data.frame(
      epithelial_state = state_name,
      status = "insufficient_variable_genes",
      n_pseudobulks = nrow(state_meta),
      genes_after_expression_filter = nrow(logcpm),
      genes_clustered = length(selected_genes),
      modules_raw = 0L,
      modules_testable = 0L,
      stage_labels_used_for_construction = FALSE
    )
    next
  }
  expression <- logcpm[selected_genes, , drop = FALSE]
  scaled <- t(scale(t(expression)))
  scaled[!is.finite(scaled)] <- 0
  gene_correlation <- stats::cor(t(scaled), use = "pairwise.complete.obs")
  gene_correlation[!is.finite(gene_correlation)] <- 0
  diag(gene_correlation) <- 1
  module_count <- min(
    as.integer(p_num("module_max_count")),
    length(selected_genes),
    max(2L, floor(length(selected_genes) / 30L))
  )
  membership <- stats::cutree(
    stats::hclust(stats::as.dist(1 - gene_correlation), method = "average"),
    k = module_count
  )

  state_members <- list()
  state_scores <- list()
  testable_count <- 0L
  for (module_number in sort(unique(membership))) {
    module_genes <- names(membership)[membership == module_number]
    if (length(module_genes) > 1L) {
      within <- gene_correlation[module_genes, module_genes, drop = FALSE]
      median_correlation <- stats::median(within[upper.tri(within)], na.rm = TRUE)
    } else {
      median_correlation <- NA_real_
    }
    testable <- length(module_genes) >= p_num("module_min_genes") &&
      is.finite(median_correlation) &&
      median_correlation >= p_num("module_min_median_correlation")
    module_id <- sprintf(
      "%s_SB_M%02d", gsub("[^A-Za-z0-9]+", "_", state_name), module_number
    )
    member_table <- data.frame(
      module_id = module_id,
      epithelial_state = state_name,
      gene = module_genes,
      gene_MAD = variability[module_genes],
      module_size = length(module_genes),
      median_within_module_correlation = median_correlation,
      testable = testable,
      stage_labels_used_for_construction = FALSE,
      stringsAsFactors = FALSE
    )
    module_membership_rows[[length(module_membership_rows) + 1L]] <- member_table
    if (!testable) next
    testable_count <- testable_count + 1L
    score <- colMeans(scaled[module_genes, , drop = FALSE])
    state_scores[[module_id]] <- score
    module_score_rows[[length(module_score_rows) + 1L]] <- data.frame(
      module_id = module_id,
      epithelial_state = state_name,
      pseudobulk_id = colnames(scaled),
      donor_id = state_meta$donor_id,
      stage = as.character(state_meta$stage),
      fap_binary = state_meta$fap_binary,
      score = score,
      stringsAsFactors = FALSE
    )
  }

  module_audit_rows[[state_name]] <- data.frame(
    epithelial_state = state_name,
    status = if (testable_count) "ok" else "no_coherent_modules",
    n_pseudobulks = nrow(state_meta),
    genes_after_expression_filter = nrow(logcpm),
    genes_clustered = length(selected_genes),
    modules_raw = module_count,
    modules_testable = testable_count,
    stage_labels_used_for_construction = FALSE,
    stringsAsFactors = FALSE
  )
  if (!testable_count) next

  score_matrix <- do.call(rbind, state_scores)
  colnames(score_matrix) <- state_meta$pseudobulk_id
  module_score_objects[[state_name]] <- list(scores = score_matrix, metadata = state_meta)
  all_fit <- fit_module_scores(score_matrix, state_meta, state_name, paired_only = FALSE)
  if (all_fit$status == "ok") {
    module_result_rows[[state_name]] <- all_fit$results
  }
  paired_fit <- fit_module_scores(score_matrix, state_meta, state_name, paired_only = TRUE)
  if (paired_fit$status == "ok") {
    paired_module_result_rows[[state_name]] <- paired_fit$results
  }
}

module_membership <- if (length(module_membership_rows)) {
  do.call(rbind, module_membership_rows)
} else {
  data.frame(
    module_id = character(), epithelial_state = character(), gene = character(),
    gene_MAD = numeric(), module_size = integer(),
    median_within_module_correlation = numeric(), testable = logical(),
    stage_labels_used_for_construction = logical()
  )
}
module_scores <- if (length(module_score_rows)) {
  do.call(rbind, module_score_rows)
} else {
  data.frame(
    module_id = character(), epithelial_state = character(),
    pseudobulk_id = character(), donor_id = character(),
    stage = character(), fap_binary = character(), score = numeric()
  )
}
module_results <- if (length(module_result_rows)) {
  do.call(rbind, module_result_rows)
} else {
  data.frame(
    feature = character(), epithelial_state = character(), analysis = character(),
    effect = numeric(), CI95_low = numeric(), CI95_high = numeric(),
    p_value = numeric(), FDR = numeric(), standard_error = numeric(),
    df_total = numeric(), feature_type = character(), contrast = character(),
    n_normal_donors = integer(), n_adenoma_donors = integer(),
    n_cancer_donors = integer(), donor_correlation = numeric(),
    model_formula = character()
  )
}
paired_module_results <- if (length(paired_module_result_rows)) {
  do.call(rbind, paired_module_result_rows)
} else {
  module_results[0, ]
}
module_audit <- do.call(rbind, module_audit_rows)

write_tsv(module_membership, file.path(result_dir, "stage_blind_module_membership.tsv"))
write_tsv(module_scores, file.path(source_dir, "stage_blind_module_scores.tsv"))
write_tsv(module_results, file.path(result_dir, "stage_blind_module_results.tsv"))
write_tsv(
  paired_module_results,
  file.path(result_dir, "paired_FAP_module_results.tsv")
)
write_tsv(module_audit, file.path(result_dir, "stage_blind_module_audit.tsv"))

exploratory_modules <- data.frame(
  module_id = character(), epithelial_state = character(),
  early_effect = numeric(), early_FDR = numeric(),
  cancer_vs_adenoma_effect = numeric(), cancer_vs_normal_effect = numeric(),
  paired_FAP_effect = numeric(), paired_FAP_FDR = numeric(),
  paired_direction_concordant = logical(), n_evaluable_LODO = integer(),
  LODO_sign_stability = numeric(), passes_LODO = logical(),
  exploratory_candidate = logical()
)
module_lodo <- data.frame(
  module_id = character(), epithelial_state = character(),
  omitted_donor = character(), effect = numeric()
)

if (nrow(module_results)) {
  wide <- reshape(
    module_results[, c("feature", "epithelial_state", "contrast", "effect", "FDR")],
    idvar = c("feature", "epithelial_state"),
    timevar = "contrast",
    direction = "wide"
  )
  needed <- c(
    "effect.adenoma_vs_normal", "FDR.adenoma_vs_normal",
    "effect.cancer_vs_adenoma", "effect.cancer_vs_normal"
  )
  if (all(needed %in% colnames(wide))) {
    wide <- wide[stats::complete.cases(wide[, needed, drop = FALSE]), ]
    wide <- wide[
      wide$FDR.adenoma_vs_normal <= p_num("exploratory_FDR_threshold") &
        abs(wide$effect.adenoma_vs_normal) >= p_num("module_min_abs_effect") &
        sign(wide$effect.cancer_vs_normal) == sign(wide$effect.adenoma_vs_normal) &
        abs(wide$effect.cancer_vs_normal) >=
          p_num("sustained_fraction") * abs(wide$effect.adenoma_vs_normal) &
        (
          sign(wide$effect.cancer_vs_adenoma) == sign(wide$effect.adenoma_vs_normal) |
            abs(wide$effect.cancer_vs_adenoma) <= p_num("module_plateau_tolerance")
        ),
      ,
      drop = FALSE
    ]
  } else {
    wide <- wide[0, , drop = FALSE]
  }

  if (nrow(wide)) {
    paired_lookup <- paired_module_results[
      paired_module_results$contrast == "adenoma_vs_normal",
      c("feature", "epithelial_state", "effect", "FDR"),
      drop = FALSE
    ]
    colnames(paired_lookup)[3:4] <- c("paired_FAP_effect", "paired_FAP_FDR")
    wide <- merge(wide, paired_lookup, by = c("feature", "epithelial_state"), all.x = TRUE)
    wide$paired_direction_concordant <- with(
      wide,
      is.finite(paired_FAP_effect) &
        sign(paired_FAP_effect) == sign(effect.adenoma_vs_normal)
    )

    lodo_rows <- list()
    for (state_name in unique(wide$epithelial_state)) {
      state_candidates <- wide[wide$epithelial_state == state_name, ]
      score_object <- module_score_objects[[state_name]]
      donors <- sort(unique(score_object$metadata$donor_id))
      for (omitted in donors) {
        fitted <- fit_module_scores(
          score_object$scores[state_candidates$feature, , drop = FALSE],
          score_object$metadata,
          state_name,
          paired_only = FALSE,
          omitted_donor = omitted
        )
        if (fitted$status != "ok") next
        selected <- fitted$results[
          fitted$results$contrast == "adenoma_vs_normal",
          c("feature", "epithelial_state", "effect"),
          drop = FALSE
        ]
        selected$omitted_donor <- omitted
        colnames(selected)[1] <- "module_id"
        lodo_rows[[length(lodo_rows) + 1L]] <- selected
      }
    }
    if (length(lodo_rows)) module_lodo <- do.call(rbind, lodo_rows)

    exploratory_modules <- do.call(rbind, lapply(seq_len(nrow(wide)), function(i) {
      candidate <- wide[i, ]
      values <- module_lodo$effect[
        module_lodo$module_id == candidate$feature &
          module_lodo$epithelial_state == candidate$epithelial_state
      ]
      stability <- if (length(values)) {
        mean(sign(values) == sign(candidate$effect.adenoma_vs_normal))
      } else {
        NA_real_
      }
      passes_lodo <- is.finite(stability) &&
        stability >= p_num("lodo_sign_stability_threshold")
      data.frame(
        module_id = candidate$feature,
        epithelial_state = candidate$epithelial_state,
        early_effect = candidate$effect.adenoma_vs_normal,
        early_FDR = candidate$FDR.adenoma_vs_normal,
        cancer_vs_adenoma_effect = candidate$effect.cancer_vs_adenoma,
        cancer_vs_normal_effect = candidate$effect.cancer_vs_normal,
        paired_FAP_effect = candidate$paired_FAP_effect,
        paired_FAP_FDR = candidate$paired_FAP_FDR,
        paired_direction_concordant = candidate$paired_direction_concordant,
        n_evaluable_LODO = length(values),
        LODO_sign_stability = stability,
        passes_LODO = passes_lodo,
        exploratory_candidate = candidate$paired_direction_concordant & passes_lodo,
        stringsAsFactors = FALSE
      )
    }))
  }
}
write_tsv(module_lodo, file.path(result_dir, "stage_blind_module_LODO.tsv"))
write_tsv(
  exploratory_modules,
  file.path(result_dir, "exploratory_candidate_modules.tsv")
)

paired_hit_summary <- if (nrow(paired_hits)) {
  aggregate(
    gene ~ epithelial_state,
    data = paired_hits,
    FUN = length
  )
} else {
  data.frame(epithelial_state = character(), gene = integer())
}
paired_hit_counts <- merge(
  paired_audit,
  paired_hit_summary,
  by = "epithelial_state",
  all.x = TRUE
)
colnames(paired_hit_counts)[colnames(paired_hit_counts) == "gene"] <- "FDR_hit_genes"
paired_hit_counts$FDR_hit_genes[is.na(paired_hit_counts$FDR_hit_genes)] <- 0L
write_tsv(
  paired_hit_counts,
  file.path(source_dir, "paired_FAP_hit_counts.tsv")
)

paired_plot <- ggplot(
  paired_hit_counts,
  aes(x = epithelial_state, y = FDR_hit_genes)
) +
  geom_col(fill = "#0072B2") +
  labs(
    x = "Epithelial state",
    y = "Paired FAP genes at FDR < 0.05"
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(figure_dir, "paired_FAP_gene_hits.pdf"),
  paired_plot, width = 7, height = 4
)
ggsave(
  file.path(figure_dir, "paired_FAP_gene_hits.png"),
  paired_plot, width = 7, height = 4, dpi = 300
)

early_modules <- module_results[
  module_results$contrast == "adenoma_vs_normal",
  ,
  drop = FALSE
]
early_modules <- early_modules[order(early_modules$p_value), ]
module_plot_data <- do.call(
  rbind,
  lapply(split(early_modules, early_modules$epithelial_state), head, n = 5L)
)
write_tsv(
  module_plot_data,
  file.path(source_dir, "top_stage_blind_module_effects.tsv")
)
if (nrow(module_plot_data)) {
  module_plot_data$label <- paste(
    module_plot_data$epithelial_state,
    module_plot_data$feature,
    sep = " | "
  )
  module_plot_data$label <- factor(
    module_plot_data$label,
    levels = rev(module_plot_data$label)
  )
  module_plot <- ggplot(
    module_plot_data,
    aes(x = effect, y = label)
  ) +
    geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.4) +
    geom_errorbarh(
      aes(xmin = CI95_low, xmax = CI95_high),
      height = 0.2, colour = "#555555"
    ) +
    geom_point(aes(colour = FDR < 0.05), size = 1.8) +
    scale_colour_manual(values = c(`FALSE` = "#777777", `TRUE` = "#D55E00")) +
    labs(
      x = "Standardized adenoma-minus-normal effect (95% CI)",
      y = NULL,
      colour = "FDR < 0.05"
    ) +
    theme_bw(base_size = 8)
  ggsave(
    file.path(figure_dir, "top_stage_blind_module_effects.pdf"),
    module_plot, width = 7, height = 8
  )
  ggsave(
    file.path(figure_dir, "top_stage_blind_module_effects.png"),
    module_plot, width = 7, height = 8, dpi = 300
  )
}

software <- data.frame(
  software = c("R", required_packages),
  version = c(
    R.version.string,
    vapply(required_packages, function(x) as.character(utils::packageVersion(x)), character(1))
  ),
  stringsAsFactors = FALSE
)
write_tsv(software, file.path(result_dir, "software_versions.tsv"))

key_metrics <- data.frame(
  metric = c(
    "paired_states_evaluable",
    "paired_gene_result_rows",
    "paired_gene_FDR_hits",
    "paired_gene_LODO_stable_hits",
    "stage_blind_modules_raw",
    "stage_blind_modules_testable",
    "stage_blind_module_result_rows",
    "exploratory_progressive_modules_pre_LODO",
    "exploratory_progressive_modules_final"
  ),
  value = c(
    sum(paired_audit$status == "ok"),
    nrow(paired_results),
    nrow(paired_hits),
    sum(paired_lodo_summary$passes_LODO, na.rm = TRUE),
    sum(module_audit$modules_raw),
    sum(module_audit$modules_testable),
    nrow(module_results),
    nrow(exploratory_modules),
    sum(exploratory_modules$exploratory_candidate, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_tsv(key_metrics, file.path(result_dir, "stage_6A_amendment_key_metrics.tsv"))

report_lines <- c(
  "# Stage 6A exploratory amendment",
  "",
  "## Status",
  "",
  "Server analysis completed; independent Codex acceptance remains pending.",
  "The primary Stage 6A negative result remains frozen and unchanged.",
  "",
  "## Paired FAP sensitivity analysis",
  "",
  sprintf("- Evaluable epithelial states: %d.", sum(paired_audit$status == "ok")),
  sprintf("- Gene-state result rows: %s.", format(nrow(paired_results), big.mark = ",")),
  sprintf("- Genes at exploratory FDR < 0.05 and |log2FC| >= 0.25: %d.", nrow(paired_hits)),
  sprintf("- LODO-stable paired hits: %d.", sum(paired_lodo_summary$passes_LODO, na.rm = TRUE)),
  "",
  "Donor was a fixed blocking factor. Cells and tissues were not treated as",
  "independent biological replicates.",
  "",
  "## Stage-blind module analysis",
  "",
  sprintf("- Raw stage-blind modules: %d.", sum(module_audit$modules_raw)),
  sprintf("- Coherent testable modules: %d.", sum(module_audit$modules_testable)),
  sprintf("- Module result rows: %d.", nrow(module_results)),
  sprintf("- Pre-LODO progressive modules: %d.", nrow(exploratory_modules)),
  sprintf(
    "- Final exploratory modules after paired-direction and LODO checks: %d.",
    sum(exploratory_modules$exploratory_candidate, na.rm = TRUE)
  ),
  "",
  "Lesion-stage labels were not used for expression filtering, variance ranking,",
  "gene correlation, hierarchical clustering or module-score orientation.",
  "",
  "## Interpretation boundary",
  "",
  "All findings are secondary and exploratory. They do not replace the primary",
  "Stage 6A result and do not authorize Stage 6B.",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(report_lines, file.path(report_dir, "stage_6A_exploratory_amendment.md"))

manifest <- c(
  "# Stage 6A exploratory amendment outputs",
  "",
  "- `paired_FAP_gene_results.tsv`",
  "- `paired_FAP_gene_hits.tsv`",
  "- `paired_FAP_gene_LODO_summary.tsv`",
  "- `stage_blind_module_membership.tsv`",
  "- `stage_blind_module_results.tsv`",
  "- `paired_FAP_module_results.tsv`",
  "- `stage_blind_module_LODO.tsv`",
  "- `exploratory_candidate_modules.tsv`",
  "- `source_data/stage_blind_module_scores.tsv`",
  "- `source_data/paired_FAP_hit_counts.tsv`",
  "- `source_data/top_stage_blind_module_effects.tsv`"
)
writeLines(manifest, file.path(result_dir, "_analysis_outputs.md"))

saveRDS(
  list(
    parameters = parameters,
    paired_audit = paired_audit,
    module_audit = module_audit,
    key_metrics = key_metrics
  ),
  file.path(project_dir, "objects", "GSE201348_6A_amendment_audit.rds")
)

cat("Stage 6A exploratory amendment completed\n")
