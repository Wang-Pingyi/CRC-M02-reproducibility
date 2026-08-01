#!/usr/bin/env Rscript

# Analysis: Stage 6B Slingshot trajectory and tradeSeq dynamics
# Date: 2026-07-28
# Random seed: 20260728
# Primary inferential unit: donor-level pseudobulk within pseudotime bins

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
if (dir.exists(private_library)) .libPaths(c(private_library, .libPaths()))

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "slingshot", "tradeSeq",
  "edgeR", "limma", "BiocParallel", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing Stage 6B trajectory packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(slingshot)
  library(tradeSeq)
  library(edgeR)
  library(limma)
  library(BiocParallel)
  library(ggplot2)
})

parameter_path <- file.path(project_dir, "config", "06B_regulatory_parameters.tsv")
epithelial_path <- file.path(
  project_dir, "objects", "GSE201348_5C_epithelial_annotated_CNV.rds"
)
candidate_path <- file.path(
  project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
)
membership_path <- file.path(
  project_dir, "results", "06A_amendment", "stage_blind_module_membership.tsv"
)
stopifnot(
  file.exists(parameter_path), file.exists(epithelial_path),
  file.exists(candidate_path), file.exists(membership_path)
)

result_dir <- file.path(project_dir, "results", "06B_regulatory_inference")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_dir, "figures", "06B_regulatory_inference")
object_dir <- file.path(project_dir, "objects")
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
aggregate_sparse <- function(counts, groups) {
  levels <- unique(groups)
  membership <- Matrix::sparseMatrix(
    i = seq_along(groups),
    j = match(groups, levels),
    x = 1,
    dims = c(length(groups), length(levels)),
    dimnames = list(colnames(counts), levels)
  )
  output <- counts %*% membership
  colnames(output) <- levels
  output
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

cat("Stage 6B trajectory started:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n")
epithelial <- readRDS(epithelial_path)
meta <- epithelial[[]]
required_meta <- c(
  "donor_id", "lesion_stage", "epithelial_state", "biological_sample_id"
)
if (length(setdiff(required_meta, colnames(meta)))) {
  stop("Epithelial object lacks required metadata")
}
stage_map <- c(normal = "normal", adenoma_polyp = "adenoma", cancer = "cancer")
if (!all(meta$lesion_stage %in% names(stage_map))) {
  stop("Unexpected lesion_stage values")
}
meta$stage <- unname(stage_map[meta$lesion_stage])

candidates <- utils::read.delim(candidate_path, check.names = FALSE)
candidates <- candidates[candidates$exploratory_candidate & candidates$passes_LODO, ]
membership <- utils::read.delim(membership_path, check.names = FALSE)
membership <- membership[membership$module_id %in% candidates$module_id, ]
candidate_genes <- sort(intersect(unique(membership$gene), rownames(epithelial)))
if (length(candidate_genes) < 10L) stop("Too few locked candidate genes in object")

trajectory_states <- p_vec("states")
keep <- meta$epithelial_state %in% trajectory_states
meta_keep <- meta[keep, , drop = FALSE]
meta_keep$cell_id <- rownames(meta_keep)
meta_keep$stratum <- paste(
  meta_keep$donor_id, meta_keep$stage, meta_keep$epithelial_state, sep = "||"
)

split_index <- split(seq_len(nrow(meta_keep)), meta_keep$stratum)
min_stratum <- as.integer(p_num("min_cells_per_stratum"))
max_stratum <- as.integer(p_num("max_cells_per_donor_stage_state"))
selected_index <- unlist(
  lapply(sort(names(split_index)), function(group_name) {
    idx <- split_index[[group_name]]
    if (length(idx) < min_stratum) return(integer())
    set.seed(as.integer(p_num("random_seed")) + match(group_name, sort(names(split_index))))
    sort(sample(idx, min(length(idx), max_stratum), replace = FALSE))
  }),
  use.names = FALSE
)
if (length(selected_index) < 1000L) stop("Trajectory downsampling retained too few cells")
selected_cells <- meta_keep$cell_id[selected_index]
trajectory_meta <- meta_keep[selected_index, , drop = FALSE]

reduction_name <- p_chr("reduction")
if (!reduction_name %in% Reductions(epithelial)) {
  stop("Missing locked trajectory reduction: ", reduction_name)
}
embedding <- Embeddings(epithelial, reduction = reduction_name)
embedding <- embedding[selected_cells, seq_len(as.integer(p_num("dimensions"))), drop = FALSE]
if (!identical(rownames(embedding), trajectory_meta$cell_id)) {
  stop("Trajectory metadata and embedding are misaligned")
}

state_counts <- table(trajectory_meta$epithelial_state)
if (!all(c(p_chr("root_state"), p_vec("terminal_states")) %in% names(state_counts))) {
  stop("Root or terminal state absent after deterministic sampling")
}

sds <- slingshot::slingshot(
  embedding,
  clusterLabels = trajectory_meta$epithelial_state,
  start.clus = p_chr("root_state"),
  end.clus = p_vec("terminal_states"),
  stretch = 0
)
pseudotime <- slingshot::slingPseudotime(sds, na = FALSE)
weights <- slingshot::slingCurveWeights(sds)
lineage_definitions <- slingshot::slingLineages(sds)
if (length(lineage_definitions) != ncol(pseudotime)) {
  stop("Slingshot lineage definitions are misaligned with pseudotime")
}
lineage_terminals <- vapply(
  lineage_definitions, function(x) tail(x, 1L), character(1)
)
keep_lineages <- lineage_terminals %in% p_vec("terminal_states")
pseudotime <- pseudotime[, keep_lineages, drop = FALSE]
weights <- weights[, keep_lineages, drop = FALSE]
lineage_definitions <- lineage_definitions[keep_lineages]
lineage_terminals <- lineage_terminals[keep_lineages]
sampled_cells_before_lineage_filter <- nrow(trajectory_meta)
keep_cells <- rowSums(weights) > 0
pseudotime <- pseudotime[keep_cells, , drop = FALSE]
weights <- weights[keep_cells, , drop = FALSE]
trajectory_meta <- trajectory_meta[keep_cells, , drop = FALSE]
selected_cells <- trajectory_meta$cell_id
embedding <- embedding[keep_cells, , drop = FALSE]
if (
  !identical(dim(pseudotime), dim(weights)) ||
    ncol(pseudotime) < 1L ||
    !setequal(lineage_terminals, p_vec("terminal_states"))
) {
  stop("Slingshot returned invalid pseudotime or lineage weights")
}
lineage_names <- make.unique(paste0(
  p_chr("root_state"), "_to_", lineage_terminals
))
colnames(pseudotime) <- lineage_names
colnames(weights) <- lineage_names
lineage_manifest <- do.call(rbind, lapply(seq_along(lineage_names), function(i) {
  data.frame(
    lineage = lineage_names[i],
    root_state = p_chr("root_state"),
    terminal_state = lineage_terminals[i],
    state_path = paste(lineage_definitions[[i]], collapse = " -> "),
    root_basis = paste(
      "Stem/progenitor cells are the renewing crypt compartment and were",
      "prespecified as the computational origin; pseudotime is not real time."
    ),
    stringsAsFactors = FALSE
  )
}))
write_tsv(
  lineage_manifest,
  file.path(result_dir, "trajectory_lineage_manifest.tsv")
)

assigned_index <- max.col(weights, ties.method = "first")
assigned_weight <- weights[cbind(seq_len(nrow(weights)), assigned_index)]
trajectory_meta$assigned_lineage <- lineage_names[assigned_index]
trajectory_meta$assigned_weight <- assigned_weight
trajectory_meta$assigned_pseudotime <- pseudotime[
  cbind(seq_len(nrow(pseudotime)), assigned_index)
]

trajectory_source <- trajectory_meta[, c(
  "cell_id", "donor_id", "biological_sample_id", "stage", "lesion_stage",
  "epithelial_state", "assigned_lineage", "assigned_weight",
  "assigned_pseudotime"
)]
for (j in seq_len(ncol(pseudotime))) {
  trajectory_source[[paste0("pseudotime_", lineage_names[j])]] <- pseudotime[, j]
  trajectory_source[[paste0("weight_", lineage_names[j])]] <- weights[, j]
}
write_tsv(
  trajectory_source,
  file.path(result_dir, "trajectory_source_data.tsv")
)

full_counts_selected <- LayerData(
  epithelial, assay = "RNA", layer = "counts"
)[, selected_cells, drop = FALSE]
trade_offset <- log(pmax(Matrix::colSums(full_counts_selected), 1))
counts_selected <- full_counts_selected[candidate_genes, , drop = FALSE]
if (!identical(colnames(counts_selected), trajectory_meta$cell_id)) {
  stop("Counts and trajectory metadata are misaligned")
}

bp <- BiocParallel::MulticoreParam(
  workers = as.integer(p_num("workers")),
  progressbar = TRUE,
  stop.on.error = TRUE
)
trade_fit <- tradeSeq::fitGAM(
  counts = counts_selected,
  pseudotime = pseudotime,
  cellWeights = weights,
  offset = trade_offset,
  nknots = as.integer(p_num("nknots")),
  verbose = TRUE,
  parallel = TRUE,
  BPPARAM = bp,
  sce = TRUE
)
association <- tradeSeq::associationTest(trade_fit, lineages = TRUE)
association$gene <- rownames(association)
association$tradeSeq_FDR <- stats::p.adjust(association$pvalue, method = "BH")
write_tsv(
  association,
  file.path(result_dir, "tradeSeq_candidate_gene_dynamics.tsv")
)

minimum_weight <- p_num("min_lineage_weight")
bin_count <- as.integer(p_num("pseudotime_bins"))
assigned <- trajectory_meta$assigned_weight >= minimum_weight &
  is.finite(trajectory_meta$assigned_pseudotime)
bin_meta <- trajectory_meta[assigned, , drop = FALSE]
bin_meta$pseudotime_scaled <- NA_real_
bin_meta$pseudotime_bin <- NA_integer_
for (lineage in lineage_names) {
  idx <- which(bin_meta$assigned_lineage == lineage)
  values <- bin_meta$assigned_pseudotime[idx]
  range_value <- range(values, finite = TRUE)
  if (diff(range_value) <= 0) next
  scaled <- (values - range_value[1L]) / diff(range_value)
  bin_meta$pseudotime_scaled[idx] <- scaled
  bin_meta$pseudotime_bin[idx] <- pmin(
    bin_count,
    floor(scaled * bin_count) + 1L
  )
}
bin_meta$group_id <- paste(
  bin_meta$donor_id, bin_meta$stage, bin_meta$assigned_lineage,
  bin_meta$pseudotime_bin, sep = "||"
)
bin_counts <- counts_selected[, rownames(bin_meta), drop = FALSE]
bin_counts <- aggregate_sparse(bin_counts, bin_meta$group_id)

bin_split <- split(seq_len(nrow(bin_meta)), bin_meta$group_id)
bin_manifest <- do.call(rbind, lapply(names(bin_split), function(group_id) {
  idx <- bin_split[[group_id]]
  data.frame(
    group_id = group_id,
    donor_id = unique(bin_meta$donor_id[idx]),
    stage = unique(bin_meta$stage[idx]),
    lineage = unique(bin_meta$assigned_lineage[idx]),
    pseudotime_bin = unique(bin_meta$pseudotime_bin[idx]),
    mean_pseudotime = mean(bin_meta$pseudotime_scaled[idx]),
    n_cells = length(idx),
    stringsAsFactors = FALSE
  )
}))
rownames(bin_manifest) <- bin_manifest$group_id
bin_manifest <- bin_manifest[colnames(bin_counts), , drop = FALSE]
bin_manifest$eligible <- bin_manifest$n_cells >=
  p_num("min_cells_per_donor_lineage_bin")
write_tsv(bin_manifest, file.path(result_dir, "trajectory_pseudobulk_manifest.tsv"))

dynamic_rows <- list()
donor_slope_rows <- list()
module_score_rows <- list()
for (lineage in lineage_names) {
  lineage_meta <- bin_manifest[
    bin_manifest$lineage == lineage & bin_manifest$eligible,
    ,
    drop = FALSE
  ]
  if (length(unique(lineage_meta$donor_id)) < p_num("min_donors_dynamic_test")) next
  lineage_counts <- bin_counts[, lineage_meta$group_id, drop = FALSE]
  dge <- edgeR::DGEList(lineage_counts)
  keep_gene <- edgeR::filterByExpr(dge, group = lineage_meta$pseudotime_bin)
  if (sum(keep_gene) < 5L) next
  dge <- edgeR::calcNormFactors(dge[keep_gene, , keep.lib.sizes = FALSE])
  lineage_meta$stage <- factor(lineage_meta$stage, levels = c("normal", "adenoma", "cancer"))
  lineage_meta$pseudotime_scaled <- as.numeric(scale(lineage_meta$mean_pseudotime))
  design <- stats::model.matrix(~ 0 + stage + pseudotime_scaled, data = lineage_meta)
  if (qr(design)$rank < ncol(design)) {
    design <- stats::model.matrix(~ pseudotime_scaled, data = lineage_meta)
  }
  duplicate <- tryCatch(
    limma::duplicateCorrelation(
      edgeR::cpm(dge, log = TRUE, prior.count = 1),
      design = design,
      block = lineage_meta$donor_id
    ),
    error = function(e) NULL
  )
  correlation <- if (
    !is.null(duplicate) && is.finite(duplicate$consensus.correlation)
  ) duplicate$consensus.correlation else 0
  voom_object <- limma::voom(
    dge, design = design, plot = FALSE,
    block = lineage_meta$donor_id, correlation = correlation
  )
  fit <- limma::lmFit(
    voom_object, design = design,
    block = lineage_meta$donor_id, correlation = correlation
  )
  fit <- safe_ebayes(fit, paste0("trajectory:", lineage))
  coefficient <- match("pseudotime_scaled", colnames(fit$coefficients))
  if (is.na(coefficient)) stop("Missing pseudotime coefficient")
  tt <- limma::topTable(
    fit, coef = coefficient, number = Inf, sort.by = "none", adjust.method = "BH"
  )
  se <- fit$stdev.unscaled[, coefficient] * fit$sigma
  critical <- stats::qt(0.975, df = fit$df.total)
  result <- data.frame(
    gene = rownames(tt),
    lineage = lineage,
    standardized_pseudotime_effect = tt$logFC,
    CI95_low = tt$logFC - critical * se,
    CI95_high = tt$logFC + critical * se,
    p_value = tt$P.Value,
    FDR = tt$adj.P.Val,
    n_donors = length(unique(lineage_meta$donor_id)),
    n_pseudobulks = nrow(lineage_meta),
    donor_correlation = correlation,
    stringsAsFactors = FALSE
  )

  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  candidate_for_stability <- result$gene[
    result$FDR < p_num("donor_bin_FDR")
  ]
  if (length(candidate_for_stability)) {
    for (gene in candidate_for_stability) {
      donor_values <- lapply(
        split(seq_len(nrow(lineage_meta)), lineage_meta$donor_id),
        function(idx) {
          if (length(unique(lineage_meta$pseudotime_bin[idx])) < 3L) return(NULL)
          slope <- stats::coef(stats::lm(
            as.numeric(logcpm[gene, idx]) ~ lineage_meta$mean_pseudotime[idx]
          ))[2L]
          if (!is.finite(slope)) return(NULL)
          data.frame(
            gene = gene,
            lineage = lineage,
            donor_id = lineage_meta$donor_id[idx[1L]],
            slope = unname(slope),
            stringsAsFactors = FALSE
          )
        }
      )
      donor_values <- donor_values[!vapply(donor_values, is.null, logical(1))]
      if (length(donor_values)) {
        donor_slope_rows[[length(donor_slope_rows) + 1L]] <- do.call(rbind, donor_values)
      }
    }
  }
  dynamic_rows[[lineage]] <- result

  scaled_expression <- t(scale(t(logcpm)))
  scaled_expression[!is.finite(scaled_expression)] <- 0
  for (module_id in candidates$module_id) {
    genes <- intersect(
      membership$gene[membership$module_id == module_id],
      rownames(scaled_expression)
    )
    if (!length(genes)) next
    module_score_rows[[length(module_score_rows) + 1L]] <- data.frame(
      module_id = module_id,
      lineage = lineage,
      group_id = lineage_meta$group_id,
      donor_id = lineage_meta$donor_id,
      stage = as.character(lineage_meta$stage),
      pseudotime_bin = lineage_meta$pseudotime_bin,
      mean_pseudotime = lineage_meta$mean_pseudotime,
      n_cells = lineage_meta$n_cells,
      module_score = colMeans(scaled_expression[genes, , drop = FALSE]),
      genes_scored = length(genes),
      stringsAsFactors = FALSE
    )
  }
}

dynamic <- if (length(dynamic_rows)) {
  do.call(rbind, dynamic_rows)
} else {
  data.frame(
    gene = character(), lineage = character(),
    standardized_pseudotime_effect = numeric(),
    CI95_low = numeric(), CI95_high = numeric(),
    p_value = numeric(), FDR = numeric(), n_donors = integer(),
    n_pseudobulks = integer(), donor_correlation = numeric(),
    evaluable_donors = integer(), donor_sign_stability = numeric(),
    tradeSeq_p_value = numeric(), tradeSeq_FDR = numeric(),
    cross_donor_dynamic = logical()
  )
}
donor_slopes <- if (length(donor_slope_rows)) {
  do.call(rbind, donor_slope_rows)
} else {
  data.frame(gene = character(), lineage = character(), donor_id = character(), slope = numeric())
}
if (nrow(dynamic)) {
  stability <- do.call(rbind, lapply(seq_len(nrow(dynamic)), function(i) {
    row <- dynamic[i, ]
    values <- donor_slopes[
      donor_slopes$gene == row$gene & donor_slopes$lineage == row$lineage,
      ,
      drop = FALSE
    ]
    data.frame(
      gene = row$gene,
      lineage = row$lineage,
      evaluable_donors = nrow(values),
      donor_sign_stability = if (nrow(values)) {
        mean(sign(values$slope) == sign(row$standardized_pseudotime_effect))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  dynamic <- merge(dynamic, stability, by = c("gene", "lineage"), all.x = TRUE)
  trade_lookup <- association[, c("gene", "pvalue", "tradeSeq_FDR"), drop = FALSE]
  colnames(trade_lookup)[2L] <- "tradeSeq_p_value"
  dynamic <- merge(dynamic, trade_lookup, by = "gene", all.x = TRUE)
  dynamic$cross_donor_dynamic <- with(
    dynamic,
    FDR < p_num("donor_bin_FDR") &
      tradeSeq_FDR < p_num("tradeseq_FDR") &
      evaluable_donors >= p_num("min_donors_dynamic_test") &
      donor_sign_stability >= p_num("donor_sign_stability")
  )
}
write_tsv(dynamic, file.path(result_dir, "trajectory_dynamic_genes.tsv"))
write_tsv(donor_slopes, file.path(source_dir, "trajectory_gene_donor_slopes.tsv"))

module_scores <- if (length(module_score_rows)) {
  do.call(rbind, module_score_rows)
} else {
  data.frame(
    module_id = character(), lineage = character(), group_id = character(),
    donor_id = character(), stage = character(), pseudotime_bin = integer(),
    mean_pseudotime = numeric(), n_cells = integer(),
    module_score = numeric(), genes_scored = integer()
  )
}
write_tsv(
  module_scores,
  file.path(source_dir, "trajectory_module_scores.tsv")
)

plot_source <- trajectory_source
plot_source$plot_x <- embedding[, 1L]
plot_source$plot_y <- embedding[, 2L]
write_tsv(plot_source, file.path(source_dir, "trajectory_embedding.tsv"))

state_plot <- ggplot(
  plot_source,
  aes(x = plot_x, y = plot_y, color = epithelial_state)
) +
  geom_point(size = 0.15, alpha = 0.35) +
  theme_classic(base_size = 9) +
  labs(
    x = "Epithelial Harmony 1", y = "Epithelial Harmony 2",
    color = "Locked state",
    title = "Slingshot input geometry and locked epithelial states"
  )
ggsave(
  file.path(figure_dir, "trajectory_locked_states.pdf"),
  state_plot, width = 7, height = 5
)
ggsave(
  file.path(figure_dir, "trajectory_locked_states.png"),
  state_plot, width = 7, height = 5, dpi = 300
)

pseudotime_plot <- ggplot(
  plot_source[plot_source$assigned_weight >= minimum_weight, ],
  aes(x = plot_x, y = plot_y, color = assigned_pseudotime)
) +
  geom_point(size = 0.15, alpha = 0.4) +
  scale_color_viridis_c(option = "C") +
  facet_wrap(~ assigned_lineage) +
  theme_classic(base_size = 9) +
  labs(
    x = "Epithelial Harmony 1", y = "Epithelial Harmony 2",
    color = "Pseudotime",
    title = "Inferred expression continua (not chronological time)"
  )
ggsave(
  file.path(figure_dir, "trajectory_pseudotime.pdf"),
  pseudotime_plot, width = 7, height = 5
)
ggsave(
  file.path(figure_dir, "trajectory_pseudotime.png"),
  pseudotime_plot, width = 7, height = 5, dpi = 300
)

saveRDS(
  list(
    slingshot = sds,
    cell_ids = trajectory_meta$cell_id,
    pseudotime = pseudotime,
    weights = weights,
    parameters = parameters
  ),
  file.path(object_dir, "GSE201348_6B_slingshot.rds"),
  compress = FALSE
)
saveRDS(
  trade_fit,
  file.path(object_dir, "GSE201348_6B_tradeSeq_candidate_modules.rds"),
  compress = FALSE
)

trajectory_audit <- data.frame(
  metric = c(
    "input_epithelial_cells", "sampled_cells_before_lineage_filter",
    "sampled_trajectory_cells", "locked_candidate_modules",
    "locked_candidate_genes", "slingshot_lineages", "donors",
    "cross_donor_dynamic_gene_lineage_rows"
  ),
  value = c(
    ncol(epithelial), sampled_cells_before_lineage_filter,
    nrow(trajectory_meta), nrow(candidates),
    length(candidate_genes), ncol(pseudotime),
    length(unique(trajectory_meta$donor_id)),
    if (nrow(dynamic)) sum(dynamic$cross_donor_dynamic, na.rm = TRUE) else 0L
  ),
  stringsAsFactors = FALSE
)
write_tsv(trajectory_audit, file.path(result_dir, "trajectory_audit.tsv"))
cat("Stage 6B trajectory completed\n")
