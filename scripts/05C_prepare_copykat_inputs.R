#!/usr/bin/env Rscript

# Analysis: Stage 5C prepare per-biological-tissue CopyKAT inputs
# Date: 2026-07-27
# CopyKAT is intentionally run one biological tissue at a time to avoid
# cross-sample batch effects. Raw UMI counts are used.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  getwd()
}

param_file <- file.path(project_root, "config", "annotation_parameters.tsv")
input_file <- file.path(
  project_root, "objects", "GSE201348_5C_annotated_preCNV.rds"
)
cache_dir <- file.path(project_root, "cache", "05C_copykat_inputs")
result_dir <- file.path(project_root, "results", "05C_annotation")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

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
max_query <- as.integer(get_param(
  "cnv", "max_query_cells_per_biological_sample", TRUE
))
min_query <- as.integer(get_param(
  "cnv", "min_query_cells_per_biological_sample", TRUE
))
reference_target <- as.integer(get_param(
  "cnv", "normal_reference_cells_per_run", TRUE
))
min_reference <- as.integer(get_param(
  "cnv", "min_reference_cells", TRUE
))
set.seed(seed)

obj <- readRDS(input_file)
required_meta <- c(
  "major_cell_type", "epithelial_cluster", "epithelial_state",
  "sample_id", "biological_sample_id", "donor_id", "lesion_stage"
)
if (!all(required_meta %in% colnames(obj@meta.data))) {
  stop(
    "Input object is missing: ",
    paste(setdiff(required_meta, colnames(obj@meta.data)), collapse = ";")
  )
}
DefaultAssay(obj) <- "RNA"
meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")

sample_stratified <- function(cells, strata, target, local_seed) {
  if (length(cells) <= target) return(cells)
  set.seed(local_seed)
  split_cells <- split(cells, strata[cells])
  selected <- character()
  base_n <- max(1L, floor(target / length(split_cells)))
  for (group_cells in split_cells) {
    selected <- c(
      selected,
      sample(group_cells, min(length(group_cells), base_n))
    )
  }
  if (length(selected) > target) {
    selected <- sample(selected, target)
  }
  if (length(selected) < target) {
    remaining <- setdiff(cells, selected)
    selected <- c(
      selected,
      sample(remaining, min(length(remaining), target - length(selected)))
    )
  }
  selected
}

healthy_reference_pool <- meta[
  donor_id %chin% c("B001", "B004") &
    lesion_stage %chin% c("normal", "Normal") &
    major_cell_type != "Epithelial",
  cell_id
]
if (length(healthy_reference_pool) < min_reference) {
  stop(
    "B001/B004 healthy non-epithelial reference pool has only ",
    length(healthy_reference_pool), " cells"
  )
}

epithelial_meta <- meta[major_cell_type == "Epithelial"]
sample_ids <- sort(unique(epithelial_meta$biological_sample_id))
plan_rows <- vector("list", length(sample_ids))
strata_lookup <- setNames(
  epithelial_meta$epithelial_cluster, epithelial_meta$cell_id
)

for (i in seq_along(sample_ids)) {
  biological_sample <- sample_ids[i]
  query_all <- epithelial_meta[
    biological_sample_id == biological_sample,
    cell_id
  ]
  sample_key <- sprintf(
    "%03d_%s",
    i,
    gsub("[^A-Za-z0-9_.-]", "_", biological_sample)
  )
  if (length(query_all) < min_query) {
    plan_rows[[i]] <- data.table(
      sample_key = sample_key,
      biological_sample_id = biological_sample,
      donor_id = paste(
        unique(epithelial_meta[
          biological_sample_id == biological_sample, donor_id
        ]),
        collapse = ";"
      ),
      lesion_stage = paste(
        unique(epithelial_meta[
          biological_sample_id == biological_sample, lesion_stage
        ]),
        collapse = ";"
      ),
      query_cells_available = length(query_all),
      query_cells_selected = 0L,
      same_sample_reference_available = 0L,
      same_sample_reference_selected = 0L,
      healthy_reference_supplement = 0L,
      total_reference_selected = 0L,
      input_cells = 0L,
      status = "not_evaluable_too_few_epithelial_cells",
      input_rds = NA_character_
    )
    next
  }

  query_selected <- sample_stratified(
    query_all,
    strata_lookup,
    max_query,
    seed + i * 1009L
  )
  same_sample_ref <- meta[
    biological_sample_id == biological_sample &
      major_cell_type != "Epithelial",
    cell_id
  ]
  set.seed(seed + i * 1013L)
  same_sample_selected <- if (length(same_sample_ref)) {
    sample(
      same_sample_ref,
      min(length(same_sample_ref), reference_target)
    )
  } else {
    character()
  }
  supplement_needed <- max(0L, reference_target - length(
    same_sample_selected
  ))
  supplement_pool <- setdiff(
    healthy_reference_pool,
    c(query_selected, same_sample_selected)
  )
  set.seed(seed + i * 1019L)
  supplement <- if (supplement_needed > 0L) {
    sample(
      supplement_pool,
      min(length(supplement_pool), supplement_needed)
    )
  } else {
    character()
  }
  reference_selected <- unique(c(same_sample_selected, supplement))

  if (length(reference_selected) < min_reference) {
    plan_rows[[i]] <- data.table(
      sample_key = sample_key,
      biological_sample_id = biological_sample,
      donor_id = paste(
        unique(epithelial_meta[
          biological_sample_id == biological_sample, donor_id
        ]),
        collapse = ";"
      ),
      lesion_stage = paste(
        unique(epithelial_meta[
          biological_sample_id == biological_sample, lesion_stage
        ]),
        collapse = ";"
      ),
      query_cells_available = length(query_all),
      query_cells_selected = length(query_selected),
      same_sample_reference_available = length(same_sample_ref),
      same_sample_reference_selected = length(same_sample_selected),
      healthy_reference_supplement = length(supplement),
      total_reference_selected = length(reference_selected),
      input_cells = length(query_selected) + length(reference_selected),
      status = "not_evaluable_insufficient_reference_cells",
      input_rds = NA_character_
    )
    next
  }

  selected_cells <- c(query_selected, reference_selected)
  if (anyDuplicated(selected_cells)) stop("Query/reference cell overlap")
  counts <- GetAssayData(
    obj, assay = "RNA", slot = "counts"
  )[, selected_cells, drop = FALSE]
  cell_manifest <- meta[
    match(selected_cells, cell_id),
    .(
      cell_id,
      sample_id,
      biological_sample_id,
      donor_id,
      lesion_stage,
      major_cell_type,
      epithelial_cluster,
      epithelial_state
    )
  ]
  cell_manifest[
    ,
    copykat_role := fifelse(
      cell_id %chin% query_selected,
      "epithelial_query",
      fifelse(
        biological_sample_id == biological_sample,
        "same_tissue_known_normal_reference",
        "B001_B004_normal_reference_supplement"
      )
    )
  ]
  input_path <- file.path(cache_dir, paste0(sample_key, ".rds"))
  saveRDS(
    list(
      raw_counts = counts,
      cell_manifest = cell_manifest,
      query_cell_ids = query_selected,
      normal_reference_cell_ids = reference_selected,
      biological_sample_id = biological_sample,
      sample_key = sample_key,
      seed = seed + i * 1021L
    ),
    input_path,
    compress = "gzip"
  )
  plan_rows[[i]] <- data.table(
    sample_key = sample_key,
    biological_sample_id = biological_sample,
    donor_id = paste(
      unique(epithelial_meta[
        biological_sample_id == biological_sample, donor_id
      ]),
      collapse = ";"
    ),
    lesion_stage = paste(
      unique(epithelial_meta[
        biological_sample_id == biological_sample, lesion_stage
      ]),
      collapse = ";"
    ),
    query_cells_available = length(query_all),
    query_cells_selected = length(query_selected),
    same_sample_reference_available = length(same_sample_ref),
    same_sample_reference_selected = length(same_sample_selected),
    healthy_reference_supplement = length(supplement),
    total_reference_selected = length(reference_selected),
    input_cells = length(selected_cells),
    status = "ready",
    input_rds = file.path(
      "cache", "05C_copykat_inputs", basename(input_path)
    )
  )
}

plan <- rbindlist(plan_rows, fill = TRUE)
fwrite(
  plan,
  file.path(result_dir, "copykat_sample_plan.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)
writeLines(
  c(
    paste0("seed=", seed),
    paste0("biological_samples_with_epithelial_cells=", length(sample_ids)),
    paste0("copykat_ready_samples=", plan[status == "ready", .N]),
    paste0("not_evaluable_samples=", plan[status != "ready", .N]),
    paste0("selected_epithelial_query_cells=", sum(plan$query_cells_selected)),
    paste0("copykat_input_cells=", sum(plan$input_cells)),
    "copykat_unit=biological_tissue_sample",
    "raw_umi_counts=yes",
    "cross_sample_expression_matrix=no"
  ),
  file.path(result_dir, "copykat_input_provenance.txt")
)
message(
  "Prepared ", plan[status == "ready", .N],
  " per-tissue CopyKAT inputs; ",
  plan[status != "ready", .N], " tissues not evaluable"
)
