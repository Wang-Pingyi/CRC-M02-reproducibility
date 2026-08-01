#!/usr/bin/env Rscript

# Analysis: Stage 6B donor-level regulator activity and pathway enrichment
# Date: 2026-07-28
# Random seed: 20260728
# Primary inferential unit: donor-level epithelial pseudobulk

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
if (dir.exists(private_library)) .libPaths(c(private_library, .libPaths()))

required_packages <- c(
  "Matrix", "edgeR", "limma", "decoupleR", "dorothea",
  "msigdbr", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing Stage 6B regulator packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(decoupleR)
  library(dorothea)
  library(msigdbr)
  library(ggplot2)
})

parameter_path <- file.path(project_dir, "config", "06B_regulatory_parameters.tsv")
pseudobulk_path <- file.path(
  project_dir, "objects", "GSE201348_6A_epithelial_pseudobulk.rds"
)
candidate_path <- file.path(
  project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
)
membership_path <- file.path(
  project_dir, "results", "06A_amendment", "stage_blind_module_membership.tsv"
)
module_score_path <- file.path(
  project_dir, "results", "06A_amendment", "source_data",
  "stage_blind_module_scores.tsv"
)
stopifnot(
  file.exists(parameter_path), file.exists(pseudobulk_path),
  file.exists(candidate_path), file.exists(membership_path),
  file.exists(module_score_path)
)

result_dir <- file.path(project_dir, "results", "06B_regulatory_inference")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_dir, "figures", "06B_regulatory_inference")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

parameters <- utils::read.delim(parameter_path, check.names = FALSE)
param <- setNames(parameters$value, parameters$parameter)
p_num <- function(name) as.numeric(param[[name]])
p_chr <- function(name) as.character(param[[name]])
p_vec <- function(name) strsplit(p_chr(name), ";", fixed = TRUE)[[1L]]
set.seed(as.integer(p_num("random_seed")))
write_tsv <- function(x, path) {
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}
safe_ebayes <- function(fit, context) {
  tryCatch(
    limma::eBayes(fit, robust = TRUE),
    error = function(e) {
      warning(
        "Robust eBayes failed; using standard eBayes. context=", context,
        "; reason=", conditionMessage(e)
      )
      limma::eBayes(fit, robust = FALSE)
    }
  )
}

cat("Stage 6B regulator analysis started\n")
pb <- readRDS(pseudobulk_path)
counts <- pb$counts
meta <- pb$metadata
stopifnot(identical(colnames(counts), rownames(meta)))

candidates <- utils::read.delim(candidate_path, check.names = FALSE)
candidates <- candidates[candidates$exploratory_candidate & candidates$passes_LODO, ]
membership_all <- utils::read.delim(membership_path, check.names = FALSE)
membership_locked <- membership_all[
  membership_all$module_id %in% candidates$module_id,
  ,
  drop = FALSE
]
module_scores <- utils::read.delim(module_score_path, check.names = FALSE)
module_scores <- module_scores[module_scores$module_id %in% candidates$module_id, ]

candidate_states <- unique(candidates$epithelial_state)
analysis_meta <- meta[
  meta$eligible & meta$epithelial_state %in% candidate_states,
  ,
  drop = FALSE
]
analysis_counts <- counts[, analysis_meta$pseudobulk_id, drop = FALSE]
dge <- edgeR::DGEList(analysis_counts)
keep <- edgeR::filterByExpr(
  dge,
  group = interaction(analysis_meta$epithelial_state, analysis_meta$stage)
)
dge <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

data("dorothea_hs", package = "dorothea", envir = environment())
network <- get("dorothea_hs", envir = environment())
required_network <- c("tf", "target", "mor", "confidence")
if (length(setdiff(required_network, colnames(network)))) {
  stop("Unexpected DoRothEA network schema")
}
network <- network[
  network$confidence %in% p_vec("dorothea_confidence") &
    network$target %in% rownames(logcpm),
  ,
  drop = FALSE
]
if (!nrow(network)) stop("DoRothEA network has no genes in expression matrix")

activity <- decoupleR::run_ulm(
  mat = logcpm,
  network = network,
  .source = "tf",
  .target = "target",
  .mor = "mor",
  minsize = as.integer(p_num("min_regulon_size"))
)
required_activity <- c("source", "condition", "score", "p_value")
if (length(setdiff(required_activity, colnames(activity)))) {
  stop("Unexpected decoupleR run_ulm output schema")
}
colnames(activity)[colnames(activity) == "source"] <- "regulator"
colnames(activity)[colnames(activity) == "condition"] <- "pseudobulk_id"
activity <- merge(
  activity,
  analysis_meta[, c("pseudobulk_id", "donor_id", "stage", "epithelial_state")],
  by = "pseudobulk_id",
  all.x = TRUE
)
activity <- activity[
  order(activity$epithelial_state, activity$regulator, activity$pseudobulk_id),
]
write_tsv(activity, file.path(result_dir, "regulator_activity.tsv"))

fit_regulator_module <- function(module_id, omitted_donor = NA_character_) {
  module_state <- candidates$epithelial_state[candidates$module_id == module_id]
  scores <- module_scores[
    module_scores$module_id == module_id &
      module_scores$epithelial_state == module_state,
    ,
    drop = FALSE
  ]
  if (!is.na(omitted_donor)) {
    scores <- scores[scores$donor_id != omitted_donor, , drop = FALSE]
  }
  state_activity <- activity[
    activity$epithelial_state == module_state &
      activity$pseudobulk_id %in% scores$pseudobulk_id,
    ,
    drop = FALSE
  ]
  regulators <- sort(unique(state_activity$regulator))
  sample_ids <- scores$pseudobulk_id[
    scores$pseudobulk_id %in% unique(state_activity$pseudobulk_id)
  ]
  scores <- scores[match(sample_ids, scores$pseudobulk_id), , drop = FALSE]
  activity_matrix <- matrix(
    NA_real_,
    nrow = length(regulators),
    ncol = length(sample_ids),
    dimnames = list(regulators, sample_ids)
  )
  activity_matrix[
    cbind(
      match(state_activity$regulator, regulators),
      match(state_activity$pseudobulk_id, sample_ids)
    )
  ] <- state_activity$score
  keep_regulator <- rowSums(is.finite(activity_matrix)) == ncol(activity_matrix)
  activity_matrix <- activity_matrix[keep_regulator, , drop = FALSE]
  if (nrow(activity_matrix) < 2L || ncol(activity_matrix) < 8L) return(NULL)

  scores$module_score_standardized <- as.numeric(scale(scores$score))
  scores$stage <- factor(scores$stage, levels = c("normal", "adenoma", "cancer"))
  design <- stats::model.matrix(
    ~ 0 + stage + module_score_standardized,
    data = scores
  )
  if (qr(design)$rank < ncol(design)) {
    design <- stats::model.matrix(~ module_score_standardized, data = scores)
  }
  duplicate <- tryCatch(
    limma::duplicateCorrelation(
      activity_matrix, design = design, block = scores$donor_id
    ),
    error = function(e) NULL
  )
  correlation <- if (
    !is.null(duplicate) && is.finite(duplicate$consensus.correlation)
  ) duplicate$consensus.correlation else 0
  fit <- limma::lmFit(
    activity_matrix, design = design,
    block = scores$donor_id, correlation = correlation
  )
  fit <- safe_ebayes(fit, paste0("regulator:", module_id))
  coefficient <- match("module_score_standardized", colnames(fit$coefficients))
  tt <- limma::topTable(
    fit, coef = coefficient, number = Inf, sort.by = "none", adjust.method = "BH"
  )
  se <- fit$stdev.unscaled[, coefficient] * fit$sigma
  critical <- stats::qt(0.975, df = fit$df.total)
  result <- data.frame(
    module_id = module_id,
    epithelial_state = module_state,
    regulator = rownames(tt),
    standardized_effect = tt$logFC,
    CI95_low = tt$logFC - critical * se,
    CI95_high = tt$logFC + critical * se,
    p_value = tt$P.Value,
    FDR = tt$adj.P.Val,
    n_donors = length(unique(scores$donor_id)),
    n_pseudobulks = nrow(scores),
    donor_correlation = correlation,
    omitted_donor = ifelse(is.na(omitted_donor), "", omitted_donor),
    stringsAsFactors = FALSE
  )
  module_genes <- membership_locked$gene[membership_locked$module_id == module_id]
  overlap <- lapply(result$regulator, function(tf) {
    tf_network <- network[network$tf == tf & network$target %in% module_genes, ]
    data.frame(
      regulator = tf,
      module_target_count = length(unique(tf_network$target)),
      module_targets = paste(sort(unique(tf_network$target)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  overlap <- do.call(rbind, overlap)
  merge(result, overlap, by = "regulator", all.x = TRUE)
}

association_rows <- lapply(candidates$module_id, fit_regulator_module)
association_rows <- association_rows[!vapply(association_rows, is.null, logical(1))]
regulator_associations <- if (length(association_rows)) {
  do.call(rbind, association_rows)
} else {
  data.frame()
}

if (nrow(regulator_associations)) {
  provisional <- regulator_associations[
    regulator_associations$FDR < p_num("regulator_association_FDR") &
      abs(regulator_associations$standardized_effect) >=
        p_num("min_abs_standardized_effect") &
      regulator_associations$module_target_count >= p_num("min_module_targets"),
    ,
    drop = FALSE
  ]
  regulator_associations$n_evaluable_LODO <- NA_integer_
  regulator_associations$LODO_sign_stability <- NA_real_
  regulator_associations$prioritized_regulator <- FALSE
  lodo_rows <- list()
  if (nrow(provisional)) {
    for (module_id in unique(provisional$module_id)) {
      module_hits <- provisional[provisional$module_id == module_id, ]
      donors <- sort(unique(
        module_scores$donor_id[module_scores$module_id == module_id]
      ))
      for (omitted in donors) {
        omitted_fit <- fit_regulator_module(module_id, omitted)
        if (is.null(omitted_fit)) next
        selected <- omitted_fit[
          omitted_fit$regulator %in% module_hits$regulator,
          c("module_id", "regulator", "standardized_effect"),
          drop = FALSE
        ]
        selected$omitted_donor <- omitted
        lodo_rows[[length(lodo_rows) + 1L]] <- selected
      }
    }
  }
  regulator_lodo <- if (length(lodo_rows)) {
    do.call(rbind, lodo_rows)
  } else {
    data.frame(
      module_id = character(), regulator = character(),
      standardized_effect = numeric(), omitted_donor = character()
    )
  }
  if (nrow(provisional)) {
    stability <- do.call(rbind, lapply(seq_len(nrow(provisional)), function(i) {
      row <- provisional[i, ]
      values <- regulator_lodo[
        regulator_lodo$module_id == row$module_id &
          regulator_lodo$regulator == row$regulator,
        ,
        drop = FALSE
      ]
      data.frame(
        module_id = row$module_id,
        regulator = row$regulator,
        n_evaluable_LODO = nrow(values),
        LODO_sign_stability = if (nrow(values)) {
          mean(sign(values$standardized_effect) == sign(row$standardized_effect))
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }))
    stability_key <- paste(stability$module_id, stability$regulator, sep = "||")
    association_key <- paste(
      regulator_associations$module_id,
      regulator_associations$regulator,
      sep = "||"
    )
    matched <- match(association_key, stability_key)
    regulator_associations$n_evaluable_LODO <- stability$n_evaluable_LODO[matched]
    regulator_associations$LODO_sign_stability <- stability$LODO_sign_stability[matched]
    regulator_associations$prioritized_regulator <- with(
      regulator_associations,
      FDR < p_num("regulator_association_FDR") &
        abs(standardized_effect) >= p_num("min_abs_standardized_effect") &
        module_target_count >= p_num("min_module_targets") &
        LODO_sign_stability >= p_num("regulator_LODO_sign_stability")
    )
    regulator_associations$prioritized_regulator[
      is.na(regulator_associations$prioritized_regulator)
    ] <- FALSE
  }
} else {
  regulator_lodo <- data.frame()
}
write_tsv(
  regulator_associations,
  file.path(result_dir, "regulator_module_associations.tsv")
)
write_tsv(
  regulator_lodo,
  file.path(source_dir, "regulator_module_LODO.tsv")
)

gene_sets <- msigdbr::msigdbr(species = "Homo sapiens")
gene_sets <- gene_sets[
  gene_sets$gs_cat == "H" | gene_sets$gs_subcat == "CP:REACTOME",
  c("gs_name", "gene_symbol"),
  drop = FALSE
]
generic_pattern <- paste(
  c(
    "^REACTOME_SIGNALING_BY_", "^REACTOME_METABOLISM_OF_",
    "^REACTOME_DISEASES_OF_", "^REACTOME_GENE_EXPRESSION",
    "^REACTOME_CELLULAR_RESPONSES_TO_", "^REACTOME_TRANSPORT_OF_SMALL_MOLECULES$"
  ),
  collapse = "|"
)
gene_sets <- gene_sets[!grepl(generic_pattern, gene_sets$gs_name), ]
gene_sets <- split(gene_sets$gene_symbol, gene_sets$gs_name)
gene_sets <- lapply(gene_sets, unique)

enrichment_rows <- list()
for (module_id in candidates$module_id) {
  state <- candidates$epithelial_state[candidates$module_id == module_id]
  module_genes <- unique(membership_locked$gene[membership_locked$module_id == module_id])
  background <- unique(
    membership_all$gene[
      membership_all$epithelial_state == state & membership_all$testable
    ]
  )
  module_genes <- intersect(module_genes, background)
  for (set_name in names(gene_sets)) {
    set_genes <- intersect(gene_sets[[set_name]], background)
    if (
      length(set_genes) < p_num("min_geneset_size") ||
        length(set_genes) > p_num("max_geneset_size")
    ) next
    overlap <- intersect(module_genes, set_genes)
    if (length(overlap) < p_num("min_overlap")) next
    a <- length(overlap)
    b <- length(module_genes) - a
    c_value <- length(set_genes) - a
    d <- length(background) - a - b - c_value
    odds_ratio <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c_value + 0.5))
    enrichment_rows[[length(enrichment_rows) + 1L]] <- data.frame(
      module_id = module_id,
      epithelial_state = state,
      pathway = set_name,
      module_genes = length(module_genes),
      pathway_genes_in_background = length(set_genes),
      overlap_count = a,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      odds_ratio = odds_ratio,
      p_value = stats::phyper(
        a - 1L, length(set_genes),
        length(background) - length(set_genes),
        length(module_genes), lower.tail = FALSE
      ),
      stringsAsFactors = FALSE
    )
  }
}
pathway_enrichment <- if (length(enrichment_rows)) {
  do.call(rbind, enrichment_rows)
} else {
  data.frame(
    module_id = character(), epithelial_state = character(),
    pathway = character(), module_genes = integer(),
    pathway_genes_in_background = integer(), overlap_count = integer(),
    overlap_genes = character(), odds_ratio = numeric(),
    p_value = numeric(), FDR = numeric()
  )
}
if (nrow(pathway_enrichment)) {
  pathway_enrichment$FDR <- ave(
    pathway_enrichment$p_value,
    pathway_enrichment$module_id,
    FUN = function(x) stats::p.adjust(x, method = "BH")
  )
  pathway_enrichment <- pathway_enrichment[
    order(pathway_enrichment$module_id, pathway_enrichment$FDR),
  ]
}
write_tsv(
  pathway_enrichment,
  file.path(result_dir, "pathway_enrichment_all.tsv")
)

pruned_rows <- list()
if (nrow(pathway_enrichment)) {
  for (module_id in unique(pathway_enrichment$module_id)) {
    candidates_pathway <- pathway_enrichment[
      pathway_enrichment$module_id == module_id &
        pathway_enrichment$FDR < p_num("pathway_FDR"),
      ,
      drop = FALSE
    ]
    retained_sets <- list()
    for (i in seq_len(nrow(candidates_pathway))) {
      genes <- strsplit(candidates_pathway$overlap_genes[i], ";", fixed = TRUE)[[1L]]
      redundant <- any(vapply(
        retained_sets,
        function(existing) {
          length(intersect(genes, existing)) / length(union(genes, existing)) >=
            p_num("redundancy_jaccard")
        },
        logical(1)
      ))
      if (!redundant) {
        pruned_rows[[length(pruned_rows) + 1L]] <- candidates_pathway[i, ]
        retained_sets[[length(retained_sets) + 1L]] <- genes
      }
    }
  }
}
pathway_pruned <- if (length(pruned_rows)) {
  do.call(rbind, pruned_rows)
} else {
  pathway_enrichment[0, , drop = FALSE]
}
write_tsv(
  pathway_pruned,
  file.path(result_dir, "pathway_enrichment_pruned.tsv")
)

if (nrow(regulator_associations)) {
  plot_data <- regulator_associations[
    regulator_associations$prioritized_regulator %in% TRUE,
    ,
    drop = FALSE
  ]
  if (!nrow(plot_data)) {
    plot_data <- do.call(
      rbind,
      lapply(
        split(regulator_associations, regulator_associations$module_id),
        function(x) head(x[order(x$FDR), ], 3L)
      )
    )
  }
  write_tsv(plot_data, file.path(source_dir, "top_regulator_effects.tsv"))
  regulator_plot <- ggplot(
    plot_data,
    aes(x = standardized_effect, y = reorder(regulator, standardized_effect))
  ) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbarh(aes(xmin = CI95_low, xmax = CI95_high), height = 0.2) +
    geom_point(aes(color = epithelial_state), size = 2) +
    facet_wrap(~ module_id, scales = "free_y") +
    theme_classic(base_size = 9) +
    labs(
      x = "Standardized activity association (95% CI)",
      y = "Predicted regulator", color = "State",
      title = "Donor-level DoRothEA regulator associations"
    )
  ggsave(
    file.path(figure_dir, "regulator_module_associations.pdf"),
    regulator_plot, width = 7, height = 6
  )
  ggsave(
    file.path(figure_dir, "regulator_module_associations.png"),
    regulator_plot, width = 7, height = 6, dpi = 300
  )
}

cat("Stage 6B regulator and pathway analysis completed\n")
