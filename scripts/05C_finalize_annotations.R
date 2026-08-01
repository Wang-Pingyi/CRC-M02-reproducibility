#!/usr/bin/env Rscript

# Analysis: Stage 5C finalize epithelial annotation and donor-dominance audit
# Date: 2026-07-27
# Scope: descriptive annotation audit only; no donor-level hypothesis testing.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}

param_file <- file.path(project_root, "config", "annotation_parameters.tsv")
evidence_file <- file.path(project_root, "metadata", "annotation_evidence.tsv")
mapping_file <- file.path(
  project_root, "config", "epithelial_cluster_annotation.tsv"
)
epithelial_file <- file.path(
  project_root, "objects",
  "GSE201348_5C_epithelial_clustered_unannotated.rds"
)
major_file <- file.path(
  project_root, "objects", "GSE201348_5C_major_annotated.rds"
)
result_dir <- file.path(project_root, "results", "05C_annotation")
source_dir <- file.path(result_dir, "source_data")
figure_dir <- file.path(project_root, "figures", "05C_annotation")
object_dir <- file.path(project_root, "objects")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

params <- fread(param_file)
get_param <- function(section_name, parameter_name, numeric = FALSE) {
  value <- params[
    section == section_name & parameter == parameter_name,
    value
  ]
  if (length(value) != 1L) {
    stop("Expected one parameter: ", section_name, "/", parameter_name)
  }
  if (numeric) as.numeric(value) else value
}
seed <- as.integer(get_param("global", "random_seed", TRUE))
set.seed(seed)

epithelial <- readRDS(epithelial_file)
major <- readRDS(major_file)
mapping <- fread(mapping_file, colClasses = "character")
required_mapping <- c(
  "cluster", "epithelial_state_blinded", "epithelial_state_final",
  "confidence", "positive_evidence", "exclusion_evidence",
  "uncertainty_note", "reviewer_context_blinded",
  "post_blind_context_decision"
)
if (!all(required_mapping %in% colnames(mapping))) {
  stop(
    "Epithelial annotation map is missing: ",
    paste(setdiff(required_mapping, colnames(mapping)), collapse = ";")
  )
}
if (anyDuplicated(mapping$cluster)) {
  stop("Epithelial annotation map has duplicate clusters")
}
observed_clusters <- sort(unique(as.character(
  epithelial$epithelial_cluster
)))
if (!setequal(observed_clusters, mapping$cluster)) {
  stop(
    "Epithelial annotation map does not match object clusters. Missing map: ",
    paste(setdiff(observed_clusters, mapping$cluster), collapse = ";"),
    "; extra map: ",
    paste(setdiff(mapping$cluster, observed_clusters), collapse = ";")
  )
}
allowed_states <- fread(evidence_file)[
  level == "epithelial", unique(annotation)
]
if (any(!mapping$epithelial_state_final %chin% allowed_states)) {
  stop(
    "Unregistered epithelial state: ",
    paste(
      unique(mapping[
        !epithelial_state_final %chin% allowed_states,
        epithelial_state_final
      ]),
      collapse = ";"
    )
  )
}
if (any(mapping$reviewer_context_blinded != "yes")) {
  stop("All epithelial-cluster reviews must begin context-blinded")
}
adenoma_rows <- mapping[epithelial_state_final == "Adenoma_like"]
if (nrow(adenoma_rows) &&
    any(adenoma_rows$post_blind_context_decision == "not_used")) {
  stop("Adenoma_like requires a documented post-blind multi-donor audit")
}

state_map <- setNames(mapping$epithelial_state_final, mapping$cluster)
blind_map <- setNames(mapping$epithelial_state_blinded, mapping$cluster)
confidence_map <- setNames(mapping$confidence, mapping$cluster)
uncertainty_map <- setNames(mapping$uncertainty_note, mapping$cluster)
epithelial$epithelial_state <- unname(
  state_map[as.character(epithelial$epithelial_cluster)]
)
epithelial$epithelial_state_blinded <- unname(
  blind_map[as.character(epithelial$epithelial_cluster)]
)
epithelial$epithelial_annotation_confidence <- unname(
  confidence_map[as.character(epithelial$epithelial_cluster)]
)
epithelial$epithelial_annotation_uncertainty <- unname(
  uncertainty_map[as.character(epithelial$epithelial_cluster)]
)
if (anyNA(epithelial$epithelial_state)) {
  stop("Epithelial annotation produced missing labels")
}

meta <- as.data.table(epithelial@meta.data, keep.rownames = "cell_id")
cluster_donor <- meta[
  ,
  .N,
  by = .(epithelial_cluster, donor_id)
]
cluster_donor[, fraction := N / sum(N), by = epithelial_cluster]
state_donor <- meta[
  ,
  .N,
  by = .(epithelial_state, donor_id)
]
state_donor[, fraction := N / sum(N), by = epithelial_state]

dominance_summary <- function(counts, group_field) {
  counts[, group_value := get(group_field)]
  counts[
    order(group_value, -fraction),
    {
      fractions <- fraction
      positive <- fractions[fractions > 0]
      entropy <- -sum(positive * log(positive))
      .(
        cells = sum(N),
        donors = uniqueN(donor_id),
        top_donor = donor_id[1],
        top_donor_cells = N[1],
        top_donor_fraction = fraction[1],
        effective_donors = exp(entropy)
      )
    },
    by = group_value
  ]
}

cluster_dominance <- dominance_summary(
  copy(cluster_donor), "epithelial_cluster"
)
setnames(cluster_dominance, "group_value", "group")
cluster_dominance[, level := "epithelial_cluster"]
state_dominance <- dominance_summary(copy(state_donor), "epithelial_state")
setnames(state_dominance, "group_value", "group")
state_dominance[, level := "epithelial_state"]
dominance <- rbind(cluster_dominance, state_dominance, fill = TRUE)
dominance[
  ,
  `:=`(
    donor_dominance_flag = top_donor_fraction >= get_param(
      "dominance", "top_donor_fraction_flag", TRUE
    ),
    insufficient_donor_flag = donors < get_param(
      "dominance", "min_donors_for_reproducible", TRUE
    ),
    small_state_flag = cells < get_param(
      "dominance", "min_cells_for_stable_state", TRUE
    )
  )
]
dominance[
  ,
  reproducibility_status := fifelse(
    donor_dominance_flag | insufficient_donor_flag | small_state_flag,
    "uncertain_or_donor_dominated",
    "multi_donor_supported"
  )
]
fwrite(
  dominance,
  file.path(result_dir, "epithelial_donor_dominance_audit.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)
fwrite(
  cluster_donor,
  file.path(source_dir, "epithelial_cluster_donor_composition.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)
fwrite(
  state_donor,
  file.path(source_dir, "epithelial_state_donor_composition.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

context_audit <- rbind(
  meta[
    ,
    .N,
    by = .(
      annotation = epithelial_cluster,
      context_level = lesion_stage,
      donor_id
    )
  ][
    ,
    `:=`(
      annotation_level = "cluster",
      context_variable = "lesion_stage"
    )
  ],
  meta[
    ,
    .N,
    by = .(
      annotation = epithelial_state,
      context_level = lesion_stage,
      donor_id
    )
  ][
    ,
    `:=`(
      annotation_level = "state",
      context_variable = "lesion_stage"
    )
  ],
  meta[
    ,
    .N,
    by = .(
      annotation = epithelial_state,
      context_level = sporadic_or_FAP,
      donor_id
    )
  ][
    ,
    `:=`(
      annotation_level = "state",
      context_variable = "sporadic_or_FAP"
    )
  ]
)
setcolorder(
  context_audit,
  c(
    "annotation_level", "annotation", "context_variable",
    "context_level", "donor_id", "N"
  )
)
context_audit[
  ,
  fraction_within_donor_annotation := N / sum(N),
  by = .(annotation_level, annotation, context_variable, donor_id)
]
fwrite(
  context_audit,
  file.path(result_dir, "epithelial_post_blind_context_audit.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

embedding <- Embeddings(epithelial, "epithelial_umap")
umap_source <- data.table(
  cell_id = rownames(embedding),
  UMAP_1 = embedding[, 1],
  UMAP_2 = embedding[, 2],
  epithelial_cluster = as.character(epithelial$epithelial_cluster),
  epithelial_state = epithelial$epithelial_state,
  annotation_confidence = epithelial$epithelial_annotation_confidence,
  donor_id = epithelial$donor_id,
  lesion_stage = epithelial$lesion_stage,
  sporadic_or_FAP = epithelial$sporadic_or_FAP
)
fwrite(
  umap_source,
  file.path(source_dir, "epithelial_annotated_umap_source.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

state_levels <- sort(unique(umap_source$epithelial_state))
state_palette <- setNames(
  grDevices::hcl.colors(length(state_levels), "Dark 3"), state_levels
)
p_state <- ggplot(
  umap_source,
  aes(UMAP_1, UMAP_2, color = epithelial_state)
) +
  geom_point(size = 0.08, alpha = 0.4) +
  scale_color_manual(values = state_palette) +
  labs(
    title = "GSE201348 epithelial-state annotation",
    color = "Epithelial state"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())
ggsave(
  file.path(figure_dir, "epithelial_state_annotation.pdf"),
  p_state, width = 8, height = 6
)
ggsave(
  file.path(figure_dir, "epithelial_state_annotation.png"),
  p_state, width = 8, height = 6, dpi = 300
)

p_donor <- ggplot(
  cluster_donor,
  aes(
    x = epithelial_cluster, y = fraction, fill = donor_id
  )
) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Donor contribution to epithelial clusters",
    x = "Epithelial cluster",
    y = "Cell fraction",
    fill = "Donor"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
ggsave(
  file.path(figure_dir, "epithelial_cluster_donor_composition.pdf"),
  p_donor, width = 10, height = 6
)
ggsave(
  file.path(figure_dir, "epithelial_cluster_donor_composition.png"),
  p_donor, width = 10, height = 6, dpi = 300
)

fwrite(
  data.table(
    cell_id = colnames(epithelial),
    donor_id = epithelial$donor_id,
    sample_id = epithelial$sample_id,
    biological_sample_id = epithelial$biological_sample_id,
    epithelial_cluster = epithelial$epithelial_cluster,
    epithelial_state_blinded = epithelial$epithelial_state_blinded,
    epithelial_state = epithelial$epithelial_state,
    confidence = epithelial$epithelial_annotation_confidence,
    uncertainty_note = epithelial$epithelial_annotation_uncertainty
  ),
  file.path(result_dir, "epithelial_cell_annotations_preCNV.tsv.gz"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)

major$epithelial_cluster <- NA_character_
major$epithelial_state <- NA_character_
major$epithelial_annotation_confidence <- NA_character_
epi_cells <- colnames(epithelial)
if (!all(epi_cells %in% colnames(major))) {
  stop("Epithelial cells are not all present in the major object")
}
major$epithelial_cluster[
  match(epi_cells, colnames(major))
] <- as.character(epithelial$epithelial_cluster)
major$epithelial_state[
  match(epi_cells, colnames(major))
] <- epithelial$epithelial_state
major$epithelial_annotation_confidence[
  match(epi_cells, colnames(major))
] <- epithelial$epithelial_annotation_confidence

epithelial@misc$stage_5C_epithelial_annotation <- list(
  seed = seed,
  mapping_file = "config/epithelial_cluster_annotation.tsv",
  evidence_file = "metadata/annotation_evidence.tsv",
  context_blinded_initial_review = TRUE,
  donor_dominance_audit = "results/05C_annotation/epithelial_donor_dominance_audit.tsv"
)
saveRDS(
  epithelial,
  file.path(object_dir, "GSE201348_5C_epithelial_annotated_preCNV.rds"),
  compress = "gzip"
)
saveRDS(
  major,
  file.path(object_dir, "GSE201348_5C_annotated_preCNV.rds"),
  compress = "gzip"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(result_dir, "annotation_finalize_sessionInfo.txt")
)
writeLines(
  c(
    paste0("seed=", seed),
    paste0("epithelial_cells=", ncol(epithelial)),
    paste0("epithelial_clusters=", uniqueN(epithelial$epithelial_cluster)),
    paste0("epithelial_states=", uniqueN(epithelial$epithelial_state)),
    paste0(
      "donor_dominated_clusters=",
      dominance[
        level == "epithelial_cluster" & donor_dominance_flag == TRUE,
        .N
      ]
    ),
    "differential_expression_by_stage=not_performed",
    "trajectory_analysis=not_performed",
    "machine_learning=not_performed"
  ),
  file.path(result_dir, "annotation_finalize_provenance.txt")
)
message(
  "Stage 5C annotation finalized before CNV: ",
  ncol(epithelial), " epithelial cells; ",
  uniqueN(epithelial$epithelial_state), " states"
)
