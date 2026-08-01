#!/usr/bin/env Rscript

# Analysis: GSE201348 donor-level epithelial pseudobulk discovery
# Date: 2026-07-28
# Random seed: 20260728
# Primary method: limma-voom with donor blocking

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else getwd()

required_packages <- c("SeuratObject", "Matrix", "edgeR", "limma", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

input_epithelial <- file.path(
  project_dir, "objects", "GSE201348_5C_epithelial_annotated_CNV.rds"
)
input_full <- file.path(
  project_dir, "objects", "GSE201348_5C_annotated_final.rds"
)
parameter_path <- file.path(project_dir, "config", "06A_pseudobulk_parameters.tsv")
stopifnot(file.exists(input_epithelial), file.exists(input_full), file.exists(parameter_path))

result_dir <- file.path(project_dir, "results", "06A_pseudobulk")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_dir, "figures", "06A_pseudobulk")
object_dir <- file.path(project_dir, "objects")
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
    x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

collapse_values <- function(x) {
  values <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(values) == 0L) "NA" else paste(values, collapse = ";")
}

aggregate_sparse_counts <- function(counts, group_id) {
  group_levels <- unique(group_id)
  membership <- Matrix::sparseMatrix(
    i = seq_along(group_id),
    j = match(group_id, group_levels),
    x = 1,
    dims = c(length(group_id), length(group_levels)),
    dimnames = list(colnames(counts), group_levels)
  )
  aggregated <- counts %*% membership
  colnames(aggregated) <- group_levels
  aggregated
}

cat("Stage 6A started:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n")
cat("Project:", project_dir, "\n")

epithelial <- readRDS(input_epithelial)
meta <- epithelial[[]]
required_meta <- c(
  "donor_id", "lesion_stage", "epithelial_state", "sporadic_or_FAP",
  "colon_or_rectum", "tumor_location", "biological_sample_id", "sample_id"
)
missing_meta <- setdiff(required_meta, colnames(meta))
if (length(missing_meta) > 0L) {
  stop("Missing required epithelial metadata: ", paste(missing_meta, collapse = ", "))
}

stage_map <- c(normal = "normal", adenoma_polyp = "adenoma", cancer = "cancer")
if (!all(unique(meta$lesion_stage) %in% names(stage_map))) {
  stop("Unexpected lesion_stage values: ", paste(setdiff(unique(meta$lesion_stage), names(stage_map)), collapse = ", "))
}
meta$stage <- unname(stage_map[as.character(meta$lesion_stage)])
meta$stage <- factor(meta$stage, levels = c("normal", "adenoma", "cancer"))
meta$fap_binary <- ifelse(meta$sporadic_or_FAP == "FAP", "FAP", "nonFAP")
meta$group_id <- paste(meta$donor_id, meta$stage, meta$epithelial_state, sep = "||")

counts <- SeuratObject::LayerData(epithelial, assay = "RNA", layer = "counts")
if (!inherits(counts, "sparseMatrix") || any(abs(counts@x - round(counts@x)) > 1e-8)) {
  stop("RNA counts are not an integer-like sparse raw-count matrix")
}

pseudobulk_counts <- aggregate_sparse_counts(counts, meta$group_id)
split_index <- split(seq_len(nrow(meta)), meta$group_id)
pseudobulk_meta <- do.call(
  rbind,
  lapply(names(split_index), function(group) {
    idx <- split_index[[group]]
    data.frame(
      pseudobulk_id = group,
      donor_id = unique(meta$donor_id[idx]),
      stage = as.character(unique(meta$stage[idx])),
      epithelial_state = unique(meta$epithelial_state[idx]),
      fap_binary = unique(meta$fap_binary[idx]),
      sporadic_or_FAP = collapse_values(meta$sporadic_or_FAP[idx]),
      colon_or_rectum = collapse_values(meta$colon_or_rectum[idx]),
      tumor_location = collapse_values(meta$tumor_location[idx]),
      biological_samples = collapse_values(meta$biological_sample_id[idx]),
      libraries = collapse_values(meta$sample_id[idx]),
      n_biological_samples = length(unique(meta$biological_sample_id[idx])),
      n_libraries = length(unique(meta$sample_id[idx])),
      n_cells = length(idx),
      library_size = sum(pseudobulk_counts[, group]),
      stringsAsFactors = FALSE
    )
  })
)
rownames(pseudobulk_meta) <- pseudobulk_meta$pseudobulk_id
pseudobulk_meta <- pseudobulk_meta[colnames(pseudobulk_counts), , drop = FALSE]
stopifnot(identical(colnames(pseudobulk_counts), rownames(pseudobulk_meta)))
pseudobulk_meta$eligible <- pseudobulk_meta$n_cells >= p_num("min_cells_per_pseudobulk") &
  pseudobulk_meta$library_size >= p_num("min_library_size")
pseudobulk_meta$exclusion_reason <- ifelse(
  pseudobulk_meta$eligible,
  "",
  ifelse(
    pseudobulk_meta$n_cells < p_num("min_cells_per_pseudobulk") &
      pseudobulk_meta$library_size < p_num("min_library_size"),
    "low_cell_count_and_library_size",
    ifelse(
      pseudobulk_meta$n_cells < p_num("min_cells_per_pseudobulk"),
      "low_cell_count",
      "low_library_size"
    )
  )
)

write_tsv(
  pseudobulk_meta[order(pseudobulk_meta$epithelial_state, pseudobulk_meta$stage, pseudobulk_meta$donor_id), ],
  file.path(result_dir, "pseudobulk_sample_manifest.tsv")
)
saveRDS(
  list(counts = pseudobulk_counts, metadata = pseudobulk_meta, parameters = parameters),
  file.path(object_dir, "GSE201348_6A_epithelial_pseudobulk.rds"),
  compress = FALSE
)

contrast_definitions <- data.frame(
  contrast = c("adenoma_vs_normal", "cancer_vs_adenoma", "cancer_vs_normal"),
  numerator = c("adenoma", "cancer", "cancer"),
  denominator = c("normal", "adenoma", "normal"),
  interpretation = c(
    "positive log2FC means higher in adenoma than normal",
    "positive log2FC means higher in cancer than adenoma",
    "positive log2FC means higher in cancer than normal"
  ),
  stringsAsFactors = FALSE
)
write_tsv(contrast_definitions, file.path(result_dir, "contrast_definitions.tsv"))

fit_state <- function(state_name, omitted_donor = NA_character_) {
  state_meta <- pseudobulk_meta[
    pseudobulk_meta$epithelial_state == state_name & pseudobulk_meta$eligible,
    ,
    drop = FALSE
  ]
  if (!is.na(omitted_donor)) {
    state_meta <- state_meta[state_meta$donor_id != omitted_donor, , drop = FALSE]
  }
  state_meta$stage <- factor(state_meta$stage, levels = c("normal", "adenoma", "cancer"))
  state_meta$fap_binary <- factor(state_meta$fap_binary, levels = c("nonFAP", "FAP"))

  donor_counts <- table(state_meta$stage, state_meta$donor_id) > 0
  donors_per_stage <- rowSums(donor_counts)
  if (any(donors_per_stage < p_num("min_donors_per_contrast_group"))) {
    return(list(
      status = "insufficient_donors",
      donors_per_stage = donors_per_stage,
      metadata = state_meta
    ))
  }

  state_counts <- pseudobulk_counts[, state_meta$pseudobulk_id, drop = FALSE]
  dge <- edgeR::DGEList(counts = state_counts)
  keep <- edgeR::filterByExpr(
    dge,
    group = state_meta$stage,
    min.count = p_num("gene_min_count"),
    min.total.count = p_num("gene_min_total_count")
  )
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  design_formula <- ~ 0 + stage + fap_binary
  design <- stats::model.matrix(design_formula, data = state_meta)
  if (qr(design)$rank < ncol(design)) {
    design_formula <- ~ 0 + stage
    design <- stats::model.matrix(design_formula, data = state_meta)
  }
  colnames(design) <- sub("^stage", "", colnames(design))
  rownames(design) <- state_meta$pseudobulk_id

  voom_first <- limma::voom(dge, design = design, plot = FALSE)
  duplicate <- tryCatch(
    limma::duplicateCorrelation(
      voom_first,
      design = design,
      block = state_meta$donor_id
    ),
    error = function(e) NULL
  )
  consensus <- if (is.null(duplicate) || !is.finite(duplicate$consensus.correlation)) {
    0
  } else {
    duplicate$consensus.correlation
  }
  voom_final <- limma::voom(
    dge,
    design = design,
    plot = FALSE,
    block = state_meta$donor_id,
    correlation = consensus
  )
  fit <- limma::lmFit(
    voom_final,
    design = design,
    block = state_meta$donor_id,
    correlation = consensus
  )

  contrast_matrix <- limma::makeContrasts(
    adenoma_vs_normal = adenoma - normal,
    cancer_vs_adenoma = cancer - adenoma,
    cancer_vs_normal = cancer - normal,
    levels = design
  )
  fit_contrasts <- limma::contrasts.fit(fit, contrast_matrix)
  fit_contrasts <- limma::eBayes(fit_contrasts, robust = TRUE)

  tables <- lapply(colnames(contrast_matrix), function(contrast_name) {
    coefficient <- match(contrast_name, colnames(fit_contrasts$coefficients))
    tt <- limma::topTable(
      fit_contrasts,
      coef = coefficient,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH"
    )
    standard_error <- fit_contrasts$stdev.unscaled[, coefficient] * fit_contrasts$sigma
    degrees_freedom <- fit_contrasts$df.total
    critical_t <- stats::qt(0.975, df = degrees_freedom)
    data.frame(
      gene = rownames(tt),
      epithelial_state = state_name,
      contrast = contrast_name,
      log2FC = tt$logFC,
      CI95_low = tt$logFC - critical_t * standard_error,
      CI95_high = tt$logFC + critical_t * standard_error,
      p_value = tt$P.Value,
      FDR = tt$adj.P.Val,
      average_logCPM = tt$AveExpr,
      t_statistic = tt$t,
      B_statistic = tt$B,
      standard_error = standard_error,
      df_total = degrees_freedom,
      n_normal_donors = unname(donors_per_stage["normal"]),
      n_adenoma_donors = unname(donors_per_stage["adenoma"]),
      n_cancer_donors = unname(donors_per_stage["cancer"]),
      n_pseudobulk_samples = nrow(state_meta),
      donor_correlation = consensus,
      model_formula = paste(deparse(design_formula), collapse = ""),
      omitted_donor = ifelse(is.na(omitted_donor), "", omitted_donor),
      stringsAsFactors = FALSE
    )
  })

  list(
    status = "ok",
    results = do.call(rbind, tables),
    voom = voom_final,
    metadata = state_meta,
    design = design,
    correlation = consensus,
    donors_per_stage = donors_per_stage
  )
}

states <- sort(unique(pseudobulk_meta$epithelial_state))
full_fits <- list()
all_results <- list()
fit_audit <- list()
for (state_name in states) {
  cat("Primary fit:", state_name, "\n")
  fitted <- fit_state(state_name)
  full_fits[[state_name]] <- fitted
  donor_text <- paste(names(fitted$donors_per_stage), fitted$donors_per_stage, sep = "=", collapse = ";")
  fit_audit[[state_name]] <- data.frame(
    epithelial_state = state_name,
    status = fitted$status,
    donors_per_stage = donor_text,
    n_eligible_pseudobulks = nrow(fitted$metadata),
    donor_correlation = ifelse(fitted$status == "ok", fitted$correlation, NA_real_),
    stringsAsFactors = FALSE
  )
  if (fitted$status == "ok") {
    all_results[[state_name]] <- fitted$results
  }
}

fit_audit <- do.call(rbind, fit_audit)
write_tsv(fit_audit, file.path(result_dir, "model_fit_audit.tsv"))
if (length(all_results) == 0L) {
  stop("No epithelial state had enough donors for all three primary contrasts")
}
pseudobulk_results <- do.call(rbind, all_results)
pseudobulk_results <- pseudobulk_results[
  order(pseudobulk_results$epithelial_state, pseudobulk_results$contrast, pseudobulk_results$FDR),
]
write_tsv(pseudobulk_results, file.path(result_dir, "pseudobulk_results.tsv"))

wide_effects <- reshape(
  pseudobulk_results[, c("gene", "epithelial_state", "contrast", "log2FC", "FDR")],
  idvar = c("gene", "epithelial_state"),
  timevar = "contrast",
  direction = "wide"
)
required_effect_columns <- c(
  "log2FC.adenoma_vs_normal", "FDR.adenoma_vs_normal",
  "log2FC.cancer_vs_adenoma", "FDR.cancer_vs_adenoma",
  "log2FC.cancer_vs_normal", "FDR.cancer_vs_normal"
)
wide_effects <- wide_effects[stats::complete.cases(wide_effects[, required_effect_columns]), ]

direction <- sign(wide_effects$log2FC.adenoma_vs_normal)
plateau_or_same_direction <- (
  sign(wide_effects$log2FC.cancer_vs_adenoma) == direction |
    abs(wide_effects$log2FC.cancer_vs_adenoma) <= p_num("plateau_tolerance_log2fc")
)
provisional <- wide_effects[
  wide_effects$FDR.adenoma_vs_normal <= p_num("fdr_threshold") &
    abs(wide_effects$log2FC.adenoma_vs_normal) >= p_num("min_abs_log2fc") &
    sign(wide_effects$log2FC.cancer_vs_normal) == direction &
    abs(wide_effects$log2FC.cancer_vs_normal) >=
      p_num("sustained_fraction") * abs(wide_effects$log2FC.adenoma_vs_normal) &
    plateau_or_same_direction,
  ,
  drop = FALSE
]
provisional$direction <- ifelse(provisional$log2FC.adenoma_vs_normal > 0, "up", "down")
write_tsv(provisional, file.path(result_dir, "provisional_progression_genes.tsv"))

lodo_rows <- list()
if (nrow(provisional) > 0L) {
  for (state_name in unique(provisional$epithelial_state)) {
    state_candidates <- provisional[provisional$epithelial_state == state_name, ]
    donors <- sort(unique(full_fits[[state_name]]$metadata$donor_id))
    for (omitted in donors) {
      cat("LODO fit:", state_name, "without", omitted, "\n")
      lodo_fit <- fit_state(state_name, omitted_donor = omitted)
      if (lodo_fit$status != "ok") {
        lodo_rows[[length(lodo_rows) + 1L]] <- data.frame(
          gene = state_candidates$gene,
          epithelial_state = state_name,
          omitted_donor = omitted,
          contrast = NA_character_,
          log2FC = NA_real_,
          status = lodo_fit$status,
          stringsAsFactors = FALSE
        )
        next
      }
      selected <- lodo_fit$results[
        lodo_fit$results$gene %in% state_candidates$gene,
        c("gene", "epithelial_state", "omitted_donor", "contrast", "log2FC"),
        drop = FALSE
      ]
      selected$status <- "ok"
      lodo_rows[[length(lodo_rows) + 1L]] <- selected
    }
  }
}
lodo_results <- if (length(lodo_rows) > 0L) {
  do.call(rbind, lodo_rows)
} else {
  data.frame(
    gene = character(), epithelial_state = character(), omitted_donor = character(),
    contrast = character(), log2FC = numeric(), status = character()
  )
}
write_tsv(lodo_results, file.path(result_dir, "leave_one_donor_out_results.tsv"))

lodo_summary <- data.frame(
  gene = character(),
  epithelial_state = character(),
  direction = character(),
  n_lodo_adenoma_vs_normal = numeric(),
  sign_stability_adenoma_vs_normal = numeric(),
  median_lodo_log2FC_adenoma_vs_normal = numeric(),
  n_lodo_cancer_vs_adenoma = numeric(),
  sign_stability_cancer_vs_adenoma = numeric(),
  median_lodo_log2FC_cancer_vs_adenoma = numeric(),
  n_lodo_cancer_vs_normal = numeric(),
  sign_stability_cancer_vs_normal = numeric(),
  median_lodo_log2FC_cancer_vs_normal = numeric(),
  stringsAsFactors = FALSE
)
if (nrow(provisional) > 0L) {
  lodo_summary_rows <- lapply(seq_len(nrow(provisional)), function(i) {
    candidate <- provisional[i, ]
    subset_lodo <- lodo_results[
      lodo_results$gene == candidate$gene &
        lodo_results$epithelial_state == candidate$epithelial_state &
        lodo_results$status == "ok",
      ,
      drop = FALSE
    ]
    summarize_contrast <- function(contrast_name, expected_direction) {
      values <- subset_lodo$log2FC[subset_lodo$contrast == contrast_name]
      values <- values[is.finite(values)]
      c(
        n = length(values),
        sign_stability = if (length(values)) mean(sign(values) == expected_direction) else NA_real_,
        median_log2FC = if (length(values)) stats::median(values) else NA_real_,
        min_log2FC = if (length(values)) min(values) else NA_real_,
        max_log2FC = if (length(values)) max(values) else NA_real_
      )
    }
    expected <- ifelse(candidate$direction == "up", 1, -1)
    ad <- summarize_contrast("adenoma_vs_normal", expected)
    ca <- summarize_contrast("cancer_vs_adenoma", expected)
    cn <- summarize_contrast("cancer_vs_normal", expected)
    data.frame(
      gene = candidate$gene,
      epithelial_state = candidate$epithelial_state,
      direction = candidate$direction,
      n_lodo_adenoma_vs_normal = ad["n"],
      sign_stability_adenoma_vs_normal = ad["sign_stability"],
      median_lodo_log2FC_adenoma_vs_normal = ad["median_log2FC"],
      n_lodo_cancer_vs_adenoma = ca["n"],
      sign_stability_cancer_vs_adenoma = ca["sign_stability"],
      median_lodo_log2FC_cancer_vs_adenoma = ca["median_log2FC"],
      n_lodo_cancer_vs_normal = cn["n"],
      sign_stability_cancer_vs_normal = cn["sign_stability"],
      median_lodo_log2FC_cancer_vs_normal = cn["median_log2FC"],
      stringsAsFactors = FALSE
    )
  })
  lodo_summary <- do.call(rbind, lodo_summary_rows)
}
write_tsv(lodo_summary, file.path(result_dir, "leave_one_donor_out_summary.tsv"))

cat("Computing epithelial-relative specificity from the full Stage 5C object\n")
rm(epithelial, counts)
invisible(gc())
full_object <- readRDS(input_full)
full_meta <- full_object[[]]
full_counts <- SeuratObject::LayerData(full_object, assay = "RNA", layer = "counts")
specificity_group <- paste(full_meta$donor_id, full_meta$major_cell_type, sep = "||")
specificity_counts <- aggregate_sparse_counts(full_counts, specificity_group)
specificity_split <- split(seq_len(nrow(full_meta)), specificity_group)
specificity_meta <- do.call(
  rbind,
  lapply(names(specificity_split), function(group) {
    idx <- specificity_split[[group]]
    data.frame(
      group_id = group,
      donor_id = unique(full_meta$donor_id[idx]),
      major_cell_type = unique(full_meta$major_cell_type[idx]),
      n_cells = length(idx),
      stringsAsFactors = FALSE
    )
  })
)
rownames(specificity_meta) <- specificity_meta$group_id
specificity_keep <- specificity_meta$n_cells >= p_num("min_cells_per_pseudobulk")
specificity_counts <- specificity_counts[, specificity_meta$group_id[specificity_keep], drop = FALSE]
specificity_meta <- specificity_meta[specificity_keep, , drop = FALSE]
specificity_dge <- edgeR::DGEList(counts = specificity_counts)
specificity_dge <- edgeR::calcNormFactors(specificity_dge)
specificity_logcpm <- edgeR::cpm(specificity_dge, log = TRUE, prior.count = 1)
cell_types <- sort(unique(specificity_meta$major_cell_type))
cell_type_medians <- sapply(cell_types, function(cell_type) {
  rowMeans(
    specificity_logcpm[, specificity_meta$major_cell_type == cell_type, drop = FALSE],
    na.rm = TRUE
  )
})
if (is.null(dim(cell_type_medians))) {
  cell_type_medians <- matrix(
    cell_type_medians,
    ncol = 1L,
    dimnames = list(rownames(specificity_logcpm), cell_types)
  )
}
non_epithelial_types <- setdiff(colnames(cell_type_medians), "Epithelial")
specificity <- data.frame(
  gene = rownames(cell_type_medians),
  epithelial_mean_logCPM = cell_type_medians[, "Epithelial"],
  strongest_non_epithelial_mean_logCPM = apply(
    cell_type_medians[, non_epithelial_types, drop = FALSE], 1, max, na.rm = TRUE
  ),
  strongest_non_epithelial_type = non_epithelial_types[
    max.col(cell_type_medians[, non_epithelial_types, drop = FALSE], ties.method = "first")
  ],
  stringsAsFactors = FALSE
)
specificity$epithelial_specificity_delta <- specificity$epithelial_mean_logCPM -
  specificity$strongest_non_epithelial_mean_logCPM
specificity$epithelial_specific <- specificity$epithelial_specificity_delta >=
  p_num("epithelial_specificity_log2cpm_delta")
write_tsv(specificity, file.path(result_dir, "epithelial_specificity.tsv"))
rm(full_object, full_counts, specificity_counts, specificity_logcpm)
invisible(gc())

candidate_genes <- provisional
if (nrow(candidate_genes) > 0L) {
  candidate_genes <- merge(
    candidate_genes,
    lodo_summary,
    by = c("gene", "epithelial_state", "direction"),
    all.x = TRUE
  )
  candidate_genes <- merge(candidate_genes, specificity, by = "gene", all.x = TRUE)
  candidate_genes$passes_lodo <- with(
    candidate_genes,
    sign_stability_adenoma_vs_normal >= p_num("lodo_sign_stability_threshold") &
      sign_stability_cancer_vs_normal >= p_num("lodo_sign_stability_threshold")
  )
  candidate_genes$single_donor_driven <- !candidate_genes$passes_lodo
  candidate_genes$priority_candidate <- candidate_genes$passes_lodo &
    candidate_genes$epithelial_specific
  candidate_genes <- candidate_genes[candidate_genes$passes_lodo, , drop = FALSE]
}
write_tsv(candidate_genes, file.path(result_dir, "candidate_genes.tsv"))

candidate_programs <- data.frame(
  program_id = character(),
  gene = character(),
  epithelial_state = character(),
  direction = character(),
  passes_lodo = logical(),
  stringsAsFactors = FALSE
)
program_summary <- data.frame(
  program_id = character(),
  epithelial_state = character(),
  direction = character(),
  n_genes = integer(),
  n_epithelial_specific_genes = integer(),
  median_log2FC_adenoma_vs_normal = numeric(),
  median_log2FC_cancer_vs_adenoma = numeric(),
  median_log2FC_cancer_vs_normal = numeric(),
  min_lodo_sign_stability = numeric(),
  stringsAsFactors = FALSE
)
module_scores <- data.frame(
  program_id = character(),
  epithelial_state = character(),
  direction = character(),
  pseudobulk_id = character(),
  donor_id = character(),
  stage = character(),
  module_score = numeric(),
  stringsAsFactors = FALSE
)
if (nrow(candidate_genes) > 0L) {
  program_rows <- list()
  summary_rows <- list()
  score_rows <- list()
  strata <- unique(candidate_genes[, c("epithelial_state", "direction")])
  for (i in seq_len(nrow(strata))) {
    state_name <- strata$epithelial_state[i]
    direction_name <- strata$direction[i]
    genes <- candidate_genes$gene[
      candidate_genes$epithelial_state == state_name &
        candidate_genes$direction == direction_name
    ]
    fit <- full_fits[[state_name]]
    expression <- fit$voom$E[genes, , drop = FALSE]
    n_modules <- max(1L, ceiling(length(genes) / p_num("max_genes_per_program")))
    if (length(genes) == 1L || n_modules == 1L) {
      membership <- rep(1L, length(genes))
    } else {
      scaled <- t(scale(t(expression)))
      scaled[!is.finite(scaled)] <- 0
      gene_correlation <- stats::cor(t(scaled), use = "pairwise.complete.obs")
      gene_correlation[!is.finite(gene_correlation)] <- 0
      diag(gene_correlation) <- 1
      membership <- stats::cutree(
        stats::hclust(stats::as.dist(1 - gene_correlation)),
        k = min(n_modules, length(genes))
      )
    }
    for (module_number in sort(unique(membership))) {
      module_genes <- genes[membership == module_number]
      program_id <- sprintf(
        "%s_progressive_%s_M%02d",
        gsub("[^A-Za-z0-9]+", "_", state_name),
        direction_name,
        module_number
      )
      members <- candidate_genes[
        candidate_genes$gene %in% module_genes &
          candidate_genes$epithelial_state == state_name,
        ,
        drop = FALSE
      ]
      members$program_id <- program_id
      program_rows[[length(program_rows) + 1L]] <- members

      z_expression <- t(scale(t(expression[module_genes, , drop = FALSE])))
      z_expression[!is.finite(z_expression)] <- 0
      score <- colMeans(z_expression)
      score_meta <- fit$metadata[colnames(expression), , drop = FALSE]
      score_rows[[length(score_rows) + 1L]] <- data.frame(
        program_id = program_id,
        epithelial_state = state_name,
        direction = direction_name,
        pseudobulk_id = colnames(expression),
        donor_id = score_meta$donor_id,
        stage = score_meta$stage,
        module_score = score,
        stringsAsFactors = FALSE
      )
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        program_id = program_id,
        epithelial_state = state_name,
        direction = direction_name,
        n_genes = length(module_genes),
        n_epithelial_specific_genes = sum(members$epithelial_specific, na.rm = TRUE),
        median_log2FC_adenoma_vs_normal = stats::median(members$log2FC.adenoma_vs_normal),
        median_log2FC_cancer_vs_adenoma = stats::median(members$log2FC.cancer_vs_adenoma),
        median_log2FC_cancer_vs_normal = stats::median(members$log2FC.cancer_vs_normal),
        min_lodo_sign_stability = min(
          members$sign_stability_adenoma_vs_normal,
          members$sign_stability_cancer_vs_normal,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  candidate_programs <- do.call(rbind, program_rows)
  program_summary <- do.call(rbind, summary_rows)
  module_scores <- do.call(rbind, score_rows)
}
write_tsv(candidate_programs, file.path(result_dir, "candidate_programs.tsv"))
write_tsv(program_summary, file.path(result_dir, "candidate_program_summary.tsv"))
write_tsv(module_scores, file.path(source_dir, "candidate_program_scores.tsv"))

if (nrow(candidate_genes) > 0L) {
  candidate_counts <- aggregate(
    gene ~ epithelial_state + direction,
    data = candidate_genes,
    FUN = length
  )
  colnames(candidate_counts)[colnames(candidate_counts) == "gene"] <- "candidate_gene_count"
} else {
  candidate_counts <- data.frame(
    epithelial_state = character(),
    direction = character(),
    candidate_gene_count = integer(),
    stringsAsFactors = FALSE
  )
}
write_tsv(candidate_counts, file.path(source_dir, "candidate_gene_counts_by_state.tsv"))

if (nrow(candidate_counts) > 0L) {
  count_plot <- ggplot(
    candidate_counts,
    aes(x = epithelial_state, y = candidate_gene_count, fill = direction)
  ) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c(down = "#0072B2", up = "#D55E00")) +
    labs(
      x = "Epithelial state",
      y = "LODO-stable candidate genes",
      fill = "Direction"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(figure_dir, "candidate_gene_counts_by_state.pdf"), count_plot, width = 7, height = 4)
  ggsave(file.path(figure_dir, "candidate_gene_counts_by_state.png"), count_plot, width = 7, height = 4, dpi = 300)
}

if (nrow(module_scores) > 0L) {
  selected_programs <- head(
    program_summary$program_id[order(-program_summary$n_genes)],
    12L
  )
  score_plot_data <- module_scores[module_scores$program_id %in% selected_programs, ]
  score_plot_data$stage <- factor(score_plot_data$stage, levels = c("normal", "adenoma", "cancer"))
  score_plot <- ggplot(
    score_plot_data,
    aes(x = stage, y = module_score, group = donor_id)
  ) +
    geom_line(alpha = 0.35, colour = "#666666") +
    geom_point(aes(colour = stage), size = 1.5) +
    facet_wrap(~ program_id, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c(normal = "#009E73", adenoma = "#E69F00", cancer = "#D55E00")) +
    labs(x = "Lesion stage", y = "Mean standardized expression", colour = "Stage") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(figure_dir, "candidate_program_scores.pdf"), score_plot, width = 8, height = 9)
  ggsave(file.path(figure_dir, "candidate_program_scores.png"), score_plot, width = 8, height = 9, dpi = 300)
}

model_covariate_audit <- data.frame(
  variable = c("donor_id", "FAP_status", "colon_or_rectum", "tumor_location"),
  handling = c(
    "random block via duplicateCorrelation",
    "fixed covariate when design is full rank",
    "audited but not modeled",
    "audited but not modeled"
  ),
  reason = c(
    "donor is the biological replicate and may contribute multiple stages",
    "explicit metadata; non-FAP and sporadic combined as nonFAP",
    "some donor-stage aggregates mix colon and rectum; sporadic cancer location is unavailable",
    "multiple sites can contribute to one donor-stage aggregate and sporadic cancer location is unavailable"
  ),
  stringsAsFactors = FALSE
)
write_tsv(model_covariate_audit, file.path(result_dir, "model_covariate_audit.tsv"))

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
    "input_epithelial_cells", "input_genes", "donors",
    "pseudobulk_samples_total", "pseudobulk_samples_eligible",
    "states_fit_successfully", "differential_result_rows",
    "provisional_progression_genes", "lodo_stable_candidate_genes",
    "candidate_programs"
  ),
  value = c(
    nrow(meta), nrow(pseudobulk_counts), length(unique(meta$donor_id)),
    nrow(pseudobulk_meta), sum(pseudobulk_meta$eligible),
    sum(fit_audit$status == "ok"), nrow(pseudobulk_results),
    nrow(provisional), nrow(candidate_genes), nrow(program_summary)
  ),
  stringsAsFactors = FALSE
)
write_tsv(key_metrics, file.path(result_dir, "stage_6A_key_metrics.tsv"))

report_lines <- c(
  "# Stage 6A — donor-level epithelial pseudobulk discovery",
  "",
  "## Execution status",
  "",
  "Server-side statistical analysis completed. This is an automatically generated",
  "analysis report pending independent Codex validation and investigator acceptance.",
  "",
  "## Statistical unit and inputs",
  "",
  sprintf("- Biological replicate: donor."),
  sprintf("- Input epithelial cells: %s.", format(nrow(meta), big.mark = ",")),
  sprintf("- Donors: %d.", length(unique(meta$donor_id))),
  sprintf("- Pseudobulk unit: donor × lesion stage × epithelial state."),
  sprintf("- Eligible pseudobulk samples: %d/%d.", sum(pseudobulk_meta$eligible), nrow(pseudobulk_meta)),
  "",
  "## Primary model",
  "",
  "- Raw counts were aggregated before statistical testing.",
  "- edgeR TMM normalization and gene filtering were followed by limma-voom.",
  "- Repeated stages from the same donor were handled using duplicateCorrelation.",
  "- FAP status was included when estimable; non-FAP and sporadic were combined as nonFAP.",
  "- Tumor location and colon/rectum were audited but not forced into the model because",
  "  donor-stage aggregates can contain mixed sites and sporadic cancer location is unavailable.",
  "- Positive log2FC always means higher expression in the later lesion stage.",
  "- Benjamini-Hochberg FDR controlled multiple testing.",
  "",
  "## Outputs",
  "",
  sprintf("- Successfully modeled epithelial states: %d/%d.", sum(fit_audit$status == "ok"), nrow(fit_audit)),
  sprintf("- Differential result rows: %s.", format(nrow(pseudobulk_results), big.mark = ",")),
  sprintf("- Provisional progression genes: %d.", nrow(provisional)),
  sprintf("- LODO-stable candidate genes: %d.", nrow(candidate_genes)),
  sprintf("- Data-driven candidate programs: %d.", nrow(program_summary)),
  "",
  "The complete results include log2FC, 95% confidence intervals, exact P values,",
  "Benjamini-Hochberg FDR, donor counts, model formula, and donor correlation.",
  "",
  "## Sensitivity and interpretation safeguards",
  "",
  "- Leave-one-donor-out models were restricted to prespecified provisional progression genes.",
  "- Genes failing the prespecified sign-stability threshold were not promoted to candidate programs.",
  "- Epithelial specificity was assessed against donor-level major-cell-type pseudobulks.",
  "- Programs were constructed from data-supported direction and coexpression, not fashionable themes.",
  "- No cell-level test was used as primary evidence.",
  "- No prognostic model, trajectory analysis, cell communication analysis, or machine learning was run.",
  "- CopyKAT candidate malignancy labels were not used to select cells.",
  "",
  "## Required next action",
  "",
  "Run the independent Stage 6A validation script, inspect all warnings and source-data",
  "tables, update STATUS.md and result_registry.tsv, then stop for investigator approval.",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
)
writeLines(report_lines, file.path(report_dir, "stage_6A_pseudobulk_discovery.md"))

manifest_lines <- c(
  "# Stage 6A analysis outputs",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "",
  "## Primary tables",
  "",
  "- `pseudobulk_results.tsv` — all donor-level differential results",
  "- `candidate_programs.tsv` — LODO-stable data-driven program membership",
  "- `candidate_program_summary.tsv` — program-level summaries",
  "- `pseudobulk_sample_manifest.tsv` — aggregation and eligibility audit",
  "- `leave_one_donor_out_results.tsv` — donor omission effects",
  "- `epithelial_specificity.tsv` — major-compartment specificity audit",
  "",
  "## Source data",
  "",
  "- `source_data/candidate_program_scores.tsv`",
  "- `source_data/candidate_gene_counts_by_state.tsv`",
  "",
  "## Figures",
  "",
  "- `figures/06A_pseudobulk/candidate_gene_counts_by_state.pdf/.png`",
  "- `figures/06A_pseudobulk/candidate_program_scores.pdf/.png`"
)
writeLines(manifest_lines, file.path(result_dir, "_analysis_outputs.md"))

saveRDS(
  list(
    parameters = parameters,
    fit_audit = fit_audit,
    contrast_definitions = contrast_definitions,
    candidate_program_summary = program_summary
  ),
  file.path(object_dir, "GSE201348_6A_model_audit.rds")
)

cat("Stage 6A analysis completed:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n")
