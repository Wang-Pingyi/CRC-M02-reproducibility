#!/usr/bin/env Rscript

# Analysis: Stage 6B LIANA-resource donor-level communication inference
# Date: 2026-07-28
# Random seed: 20260728
# Primary inferential unit: matched donor-by-stage aggregates

set.seed(20260728)
options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
private_library <- file.path(project_dir, "environment", "R", "6B-library")
if (dir.exists(private_library)) .libPaths(c(private_library, .libPaths()))

required_packages <- c(
  "SeuratObject", "Matrix", "edgeR", "liana", "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing Stage 6B communication packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(liana)
  library(ggplot2)
})

parameter_path <- file.path(project_dir, "config", "06B_regulatory_parameters.tsv")
major_object_path <- file.path(
  project_dir, "objects", "GSE201348_5C_annotated_final.rds"
)
candidate_path <- file.path(
  project_dir, "results", "06A_amendment", "exploratory_candidate_modules.tsv"
)
module_score_path <- file.path(
  project_dir, "results", "06A_amendment", "source_data",
  "stage_blind_module_scores.tsv"
)
regulator_path <- file.path(
  project_dir, "results", "06B_regulatory_inference",
  "regulator_module_associations.tsv"
)
stopifnot(
  file.exists(parameter_path), file.exists(major_object_path),
  file.exists(candidate_path), file.exists(module_score_path),
  file.exists(regulator_path)
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
first_existing <- function(values, choices) {
  hit <- choices[choices %in% values]
  if (length(hit)) hit[[1L]] else NA_character_
}

cat("Stage 6B donor-level communication started\n")
resource <- liana::select_resource("Consensus")
if (is.list(resource) && !is.data.frame(resource) && length(resource) == 1L) {
  resource <- resource[[1L]]
}
if (!is.data.frame(resource)) {
  stop("LIANA consensus resource did not resolve to a data frame")
}
ligand_column <- first_existing(
  colnames(resource),
  c("ligand", "source_genesymbol", "source", "ligand_complex")
)
receptor_column <- first_existing(
  colnames(resource),
  c("receptor", "target_genesymbol", "target", "receptor_complex")
)
ligand_category_column <- first_existing(
  colnames(resource), c("category_intercell_source", "ligand_category")
)
receptor_category_column <- first_existing(
  colnames(resource), c("category_intercell_target", "receptor_category")
)
if (
  is.na(ligand_column) || is.na(receptor_column) ||
    is.na(ligand_category_column) || is.na(receptor_category_column)
) {
  stop("Unable to identify ligand, receptor or intercellular-category columns")
}
resource_compact <- unique(data.frame(
  ligand = as.character(resource[[ligand_column]]),
  receptor = as.character(resource[[receptor_column]]),
  ligand_category = as.character(resource[[ligand_category_column]]),
  receptor_category = as.character(resource[[receptor_category_column]]),
  stringsAsFactors = FALSE
))
resource_compact <- resource_compact[
  !is.na(resource_compact$ligand) & !is.na(resource_compact$receptor) &
    nzchar(resource_compact$ligand) & nzchar(resource_compact$receptor) &
    resource_compact$ligand_category %in% p_vec("ligand_categories") &
    resource_compact$receptor_category %in% p_vec("receptor_categories"),
]
resource_compact <- resource_compact[
  !duplicated(resource_compact[, c("ligand", "receptor")]),
]

major <- readRDS(major_object_path)
meta <- major[[]]
required_meta <- c(
  "donor_id", "lesion_stage", "major_cell_type", "epithelial_state"
)
if (length(setdiff(required_meta, colnames(meta)))) {
  stop("Major-cell object lacks communication metadata")
}
stage_map <- c(normal = "normal", adenoma_polyp = "adenoma", cancer = "cancer")
if (!all(meta$lesion_stage %in% names(stage_map))) stop("Unexpected lesion stage")
meta$stage <- unname(stage_map[meta$lesion_stage])

available_genes <- rownames(major)
resource_compact <- resource_compact[
  resource_compact$ligand %in% available_genes &
    resource_compact$receptor %in% available_genes,
]
resource_compact$interaction_id <- paste(
  resource_compact$ligand, resource_compact$receptor, sep = "->"
)
if (!nrow(resource_compact)) {
  stop("No monomeric LIANA interactions map to the expression matrix")
}
write_tsv(
  resource_compact,
  file.path(source_dir, "liana_consensus_monomeric_resource.tsv")
)

candidates <- utils::read.delim(candidate_path, check.names = FALSE)
candidates <- candidates[candidates$exploratory_candidate & candidates$passes_LODO, ]
module_scores <- utils::read.delim(module_score_path, check.names = FALSE)
module_scores <- module_scores[module_scores$module_id %in% candidates$module_id, ]
regulators <- utils::read.delim(regulator_path, check.names = FALSE)
if (!"prioritized_regulator" %in% colnames(regulators)) {
  regulators$prioritized_regulator <- FALSE
}

sender_types <- p_vec("sender_cell_types")
receiver_states <- p_vec("receiver_states")
sender_index <- which(meta$major_cell_type %in% sender_types)
receiver_index <- which(
  meta$major_cell_type == "Epithelial" &
    meta$epithelial_state %in% receiver_states
)
if (!length(sender_index) || !length(receiver_index)) {
  stop("No configured sender or receiver cells")
}

communication_index <- c(sender_index, receiver_index)
communication_meta <- meta[communication_index, , drop = FALSE]
communication_meta$role <- c(
  rep("sender", length(sender_index)),
  rep("receiver", length(receiver_index))
)
communication_meta$cell_group <- ifelse(
  communication_meta$role == "sender",
  communication_meta$major_cell_type,
  communication_meta$epithelial_state
)
communication_meta$group_id <- paste(
  communication_meta$role, communication_meta$donor_id,
  communication_meta$stage, communication_meta$cell_group,
  sep = "||"
)

resource_genes <- sort(unique(c(
  resource_compact$ligand, resource_compact$receptor
)))
counts <- LayerData(
  major, assay = "RNA", layer = "counts"
)[resource_genes, communication_index, drop = FALSE]
if (!identical(colnames(counts), rownames(communication_meta))) {
  stop("Communication counts and metadata are misaligned")
}
group_counts <- aggregate_sparse(counts, communication_meta$group_id)
group_detected <- aggregate_sparse(counts > 0, communication_meta$group_id)
group_split <- split(seq_len(nrow(communication_meta)), communication_meta$group_id)
group_manifest <- do.call(rbind, lapply(names(group_split), function(group_id) {
  idx <- group_split[[group_id]]
  data.frame(
    group_id = group_id,
    role = unique(communication_meta$role[idx]),
    donor_id = unique(communication_meta$donor_id[idx]),
    stage = unique(communication_meta$stage[idx]),
    cell_group = unique(communication_meta$cell_group[idx]),
    n_cells = length(idx),
    stringsAsFactors = FALSE
  )
}))
rownames(group_manifest) <- group_manifest$group_id
group_manifest <- group_manifest[colnames(group_counts), , drop = FALSE]
group_manifest$eligible <- group_manifest$n_cells >= p_num("min_cells_per_group")
write_tsv(group_manifest, file.path(result_dir, "communication_group_manifest.tsv"))

eligible_groups <- group_manifest$group_id[group_manifest$eligible]
group_counts <- group_counts[, eligible_groups, drop = FALSE]
group_detected <- group_detected[, eligible_groups, drop = FALSE]
eligible_manifest <- group_manifest[eligible_groups, , drop = FALSE]
expression_logcpm <- edgeR::cpm(group_counts, log = TRUE, prior.count = 1)
expression_fraction <- sweep(
  as.matrix(group_detected),
  2L,
  eligible_manifest$n_cells,
  "/"
)
rm(major, counts, group_counts, group_detected)
invisible(gc())

get_expression <- function(gene, manifest_rows) {
  ids <- manifest_rows$group_id
  data.frame(
    group_id = ids,
    donor_id = manifest_rows$donor_id,
    stage = manifest_rows$stage,
    cell_group = manifest_rows$cell_group,
    logCPM = as.numeric(expression_logcpm[gene, ids]),
    expression_fraction = as.numeric(expression_fraction[gene, ids]),
    stringsAsFactors = FALSE
  )
}

association_rows <- list()
source_rows <- list()
for (module_id in candidates$module_id) {
  receiver_state <- candidates$epithelial_state[candidates$module_id == module_id]
  if (!receiver_state %in% receiver_states) next
  score_data <- module_scores[module_scores$module_id == module_id, ]
  receiver_manifest <- eligible_manifest[
    eligible_manifest$role == "receiver" &
      eligible_manifest$cell_group == receiver_state,
    ,
    drop = FALSE
  ]
  score_receiver <- merge(
    score_data[, c("module_id", "donor_id", "stage", "score")],
    receiver_manifest[, c("group_id", "donor_id", "stage")],
    by = c("donor_id", "stage")
  )
  for (sender_type in sender_types) {
    sender_manifest <- eligible_manifest[
      eligible_manifest$role == "sender" &
        eligible_manifest$cell_group == sender_type,
      ,
      drop = FALSE
    ]
    matched_base <- merge(
      score_receiver,
      sender_manifest[, c("group_id", "donor_id", "stage")],
      by = c("donor_id", "stage"),
      suffixes = c("_receiver", "_sender")
    )
    if (
      nrow(matched_base) < p_num("min_matched_observations") ||
        length(unique(matched_base$donor_id)) < p_num("min_unique_donors")
    ) next

    for (interaction_index in seq_len(nrow(resource_compact))) {
      ligand <- resource_compact$ligand[interaction_index]
      receptor <- resource_compact$receptor[interaction_index]
      ligand_log <- expression_logcpm[ligand, matched_base$group_id_sender]
      ligand_fraction <- expression_fraction[ligand, matched_base$group_id_sender]
      receptor_log <- expression_logcpm[receptor, matched_base$group_id_receiver]
      receptor_fraction <- expression_fraction[receptor, matched_base$group_id_receiver]
      expressed <- ligand_log >= p_num("min_logCPM") &
        ligand_fraction >= p_num("min_expression_fraction") &
        receptor_log >= p_num("min_logCPM") &
        receptor_fraction >= p_num("min_expression_fraction")
      if (
        sum(expressed) < p_num("min_matched_observations") ||
          length(unique(matched_base$donor_id[expressed])) <
            p_num("min_unique_donors")
      ) next
      test <- suppressWarnings(stats::cor.test(
        as.numeric(ligand_log[expressed]),
        matched_base$score[expressed],
        method = "spearman",
        exact = FALSE
      ))
      if (!is.finite(test$estimate) || !is.finite(test$p.value)) next
      row_id <- paste(module_id, sender_type, ligand, receptor, sep = "||")
      association_rows[[length(association_rows) + 1L]] <- data.frame(
        row_id = row_id,
        module_id = module_id,
        receiver_state = receiver_state,
        sender_cell_type = sender_type,
        ligand = ligand,
        receptor = receptor,
        interaction_id = resource_compact$interaction_id[interaction_index],
        spearman_rho = unname(test$estimate),
        p_value = test$p.value,
        n_matched_observations = sum(expressed),
        n_unique_donors = length(unique(matched_base$donor_id[expressed])),
        median_ligand_logCPM = stats::median(ligand_log[expressed]),
        median_ligand_fraction = stats::median(ligand_fraction[expressed]),
        median_receptor_logCPM = stats::median(receptor_log[expressed]),
        median_receptor_fraction = stats::median(receptor_fraction[expressed]),
        stringsAsFactors = FALSE
      )
      source_rows[[length(source_rows) + 1L]] <- data.frame(
        row_id = row_id,
        donor_id = matched_base$donor_id[expressed],
        stage = matched_base$stage[expressed],
        module_score = matched_base$score[expressed],
        ligand_logCPM = as.numeric(ligand_log[expressed]),
        ligand_fraction = as.numeric(ligand_fraction[expressed]),
        receptor_logCPM = as.numeric(receptor_log[expressed]),
        receptor_fraction = as.numeric(receptor_fraction[expressed]),
        stringsAsFactors = FALSE
      )
    }
  }
}

associations <- if (length(association_rows)) {
  do.call(rbind, association_rows)
} else {
  data.frame(
    row_id = character(), module_id = character(),
    receiver_state = character(), sender_cell_type = character(),
    ligand = character(), receptor = character(), interaction_id = character(),
    spearman_rho = numeric(), p_value = numeric(),
    n_matched_observations = integer(), n_unique_donors = integer(),
    median_ligand_logCPM = numeric(), median_ligand_fraction = numeric(),
    median_receptor_logCPM = numeric(), median_receptor_fraction = numeric()
  )
}
association_source <- if (length(source_rows)) {
  do.call(rbind, source_rows)
} else {
  data.frame(
    row_id = character(), donor_id = character(), stage = character(),
    module_score = numeric(), ligand_logCPM = numeric(),
    ligand_fraction = numeric(), receptor_logCPM = numeric(),
    receptor_fraction = numeric()
  )
}

if (nrow(associations)) {
  associations$FDR <- ave(
    associations$p_value,
    associations$module_id,
    FUN = function(x) stats::p.adjust(x, method = "BH")
  )
  provisional <- associations[
    associations$FDR < p_num("communication_association_FDR") &
      abs(associations$spearman_rho) >= p_num("min_abs_spearman_rho"),
    ,
    drop = FALSE
  ]
  stability_rows <- lapply(seq_len(nrow(provisional)), function(i) {
    row <- provisional[i, ]
    values <- association_source[association_source$row_id == row$row_id, ]
    donor_ids <- sort(unique(values$donor_id))
    lodo_rho <- vapply(donor_ids, function(omitted) {
      retained <- values$donor_id != omitted
      if (
        sum(retained) < p_num("min_matched_observations") - 1L ||
          length(unique(values$donor_id[retained])) < p_num("min_unique_donors") - 1L
      ) return(NA_real_)
      suppressWarnings(stats::cor(
        values$ligand_logCPM[retained],
        values$module_score[retained],
        method = "spearman",
        use = "complete.obs"
      ))
    }, numeric(1))
    lodo_rho <- lodo_rho[is.finite(lodo_rho)]
    data.frame(
      row_id = row$row_id,
      n_evaluable_LODO = length(lodo_rho),
      LODO_sign_stability = if (length(lodo_rho)) {
        mean(sign(lodo_rho) == sign(row$spearman_rho))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  stability <- if (length(stability_rows)) do.call(rbind, stability_rows) else data.frame()
  associations$n_evaluable_LODO <- NA_integer_
  associations$LODO_sign_stability <- NA_real_
  if (nrow(stability)) {
    matched <- match(associations$row_id, stability$row_id)
    associations$n_evaluable_LODO <- stability$n_evaluable_LODO[matched]
    associations$LODO_sign_stability <- stability$LODO_sign_stability[matched]
  }
  associations$cross_donor_stable <- with(
    associations,
    FDR < p_num("communication_association_FDR") &
      abs(spearman_rho) >= p_num("min_abs_spearman_rho") &
      LODO_sign_stability >= p_num("communication_LODO_sign_stability")
  )
  associations$cross_donor_stable[is.na(associations$cross_donor_stable)] <- FALSE
} else {
  associations$FDR <- numeric()
  associations$n_evaluable_LODO <- integer()
  associations$LODO_sign_stability <- numeric()
  associations$cross_donor_stable <- logical()
}
write_tsv(
  associations,
  file.path(result_dir, "ligand_receptor_donor_associations.tsv")
)

stable_lr <- associations[associations$cross_donor_stable %in% TRUE, , drop = FALSE]
prioritized_regulators <- regulators[
  regulators$prioritized_regulator %in% TRUE,
  ,
  drop = FALSE
]
if (nrow(stable_lr) && nrow(prioritized_regulators)) {
  best_regulator <- do.call(
    rbind,
    lapply(
      split(prioritized_regulators, prioritized_regulators$module_id),
      function(x) x[order(x$FDR, -abs(x$standardized_effect)), ][1L, ]
    )
  )
  integrated <- merge(
    stable_lr,
    best_regulator[, c(
      "module_id", "regulator", "standardized_effect", "FDR",
      "module_target_count", "LODO_sign_stability"
    )],
    by = "module_id",
    suffixes = c("_communication", "_regulator")
  )
  integrated <- integrated[
    order(
      integrated$FDR_communication,
      integrated$FDR_regulator,
      -abs(integrated$spearman_rho)
    ),
  ]
  selected_rows <- integer()
  used_modules <- character()
  for (i in seq_len(nrow(integrated))) {
    if (integrated$module_id[i] %in% used_modules) next
    selected_rows <- c(selected_rows, i)
    used_modules <- c(used_modules, integrated$module_id[i])
    if (length(selected_rows) >= p_num("max_candidate_axes")) break
  }
  prioritized <- integrated[selected_rows, , drop = FALSE]
  prioritized$axis_id <- paste0("Axis_", seq_len(nrow(prioritized)))
  prioritized$association_direction <- ifelse(
    prioritized$spearman_rho >= 0, "positive", "inverse"
  )
  prioritized$axis_label <- paste(
    prioritized$sender_cell_type,
    paste0(prioritized$ligand, "-", prioritized$receptor),
    "association in", prioritized$receiver_state, "; predicted",
    paste0(prioritized$regulator, "-", prioritized$module_id),
    sprintf(
      "(%s donor-level association, rho=%.3f)",
      prioritized$association_direction, prioritized$spearman_rho
    )
  )
  prioritized$interpretation <- paste(
    "Integrated candidate association with explicit LIANA ligand and receptor",
    "categories, matched donor-level expression, module correlation and inferred",
    "TF activity; the notation does not imply activation, signaling direction or causality"
  )
} else {
  prioritized <- data.frame(
    axis_id = character(), axis_label = character(), module_id = character(),
    receiver_state = character(), sender_cell_type = character(),
    ligand = character(), receptor = character(), regulator = character(),
    interpretation = character()
  )
}
write_tsv(
  prioritized,
  file.path(result_dir, "prioritized_ligand_receptor.tsv")
)

if (nrow(prioritized)) {
  selected_source <- association_source[
    association_source$row_id %in% prioritized$row_id,
    ,
    drop = FALSE
  ]
  write_tsv(
    selected_source,
    file.path(source_dir, "prioritized_ligand_receptor_source_data.tsv")
  )
  plot_data <- merge(
    selected_source,
    prioritized[, c("row_id", "axis_id", "axis_label")],
    by = "row_id"
  )
  communication_plot <- ggplot(
    plot_data,
    aes(x = ligand_logCPM, y = module_score, color = stage)
  ) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30") +
    facet_wrap(~ axis_id, scales = "free") +
    theme_classic(base_size = 9) +
    labs(
      x = "Sender ligand logCPM",
      y = "Receiver candidate-module score",
      color = "Lesion stage",
      title = "Matched donor-level candidate ligand-receptor associations"
    )
  ggsave(
    file.path(figure_dir, "prioritized_ligand_receptor.pdf"),
    communication_plot, width = 7, height = 4
  )
  ggsave(
    file.path(figure_dir, "prioritized_ligand_receptor.png"),
    communication_plot, width = 7, height = 4, dpi = 300
  )
} else {
  write_tsv(
    data.frame(
      row_id = character(), donor_id = character(), stage = character(),
      module_score = numeric(), ligand_logCPM = numeric(),
      ligand_fraction = numeric(), receptor_logCPM = numeric(),
      receptor_fraction = numeric()
    ),
    file.path(source_dir, "prioritized_ligand_receptor_source_data.tsv")
  )
}

communication_audit <- data.frame(
  metric = c(
    "LIANA_category_valid_monomeric_interactions",
    "eligible_sender_receiver_groups",
    "donor_level_associations_tested", "cross_donor_stable_interactions",
    "integrated_candidate_axes"
  ),
  value = c(
    nrow(resource_compact), nrow(eligible_manifest), nrow(associations),
    sum(associations$cross_donor_stable, na.rm = TRUE), nrow(prioritized)
  ),
  stringsAsFactors = FALSE
)
write_tsv(communication_audit, file.path(result_dir, "communication_audit.tsv"))
cat("Stage 6B donor-level communication completed\n")
