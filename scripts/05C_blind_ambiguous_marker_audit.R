#!/usr/bin/env Rscript

# Marker-only self-review of ambiguous major clusters. Donor, lesion stage,
# sample identity and FAP status are intentionally not used.

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
input_file <- file.path(
  project_root, "objects", "GSE201348_5C_major_clustered_unannotated.rds"
)
result_dir <- file.path(project_root, "results", "05C_annotation")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_file)
programs <- list(
  `15/Tuft` = c(
    "POU2F3", "TRPM5", "SH2D6", "SH2D7", "LRMP", "IL17RB"
  ),
  `15/Mast_cell` = c(
    "TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT"
  ),
  `17/Enteroendocrine` = c(
    "CHGA", "CHGB", "TPH1", "SCGN", "PAX6", "RFX6"
  ),
  `17/Pericyte_SMC` = c(
    "RGS5", "ACTA2", "MYH11", "PDGFRB", "CSPG4"
  ),
  `18/Fibroblast_FDC` = c(
    "FDCSP", "PAPPA", "VCAM1", "COL1A1", "DCN", "LUM"
  ),
  `18/B_cell` = c(
    "CD79A", "MS4A1", "CD37", "CD74", "CD22"
  ),
  `18/Enteric_glia` = c(
    "PLP1", "S100B", "SOX10", "SLC1A3", "S100A1"
  ),
  `19/Epithelial` = c(
    "EPCAM", "KRT8", "KRT18", "KRT19", "SOX9", "ELF3"
  ),
  `19/Immune` = c(
    "PTPRC", "CD3D", "MS4A1", "LYZ", "FCER1G"
  )
)

normalized <- GetAssayData(obj, assay = "RNA", slot = "data")
summary_rows <- lapply(names(programs), function(program_name) {
  pieces <- strsplit(program_name, "/", fixed = TRUE)[[1]]
  cluster_id <- pieces[1]
  label <- pieces[2]
  cells <- colnames(obj)[as.character(obj$major_cluster) == cluster_id]
  genes <- intersect(programs[[program_name]], rownames(normalized))
  detected <- Matrix::colSums(
    normalized[genes, cells, drop = FALSE] > 0
  )
  data.table(
    cluster = cluster_id,
    candidate_program = label,
    cells = length(cells),
    genes_requested = paste(programs[[program_name]], collapse = ";"),
    genes_present = paste(genes, collapse = ";"),
    mean_markers_detected_per_cell = mean(detected),
    cells_with_at_least_one_marker = sum(detected >= 1),
    fraction_with_at_least_one_marker = mean(detected >= 1),
    cells_with_at_least_two_markers = sum(detected >= 2),
    fraction_with_at_least_two_markers = mean(detected >= 2)
  )
})
audit <- rbindlist(summary_rows)
fwrite(
  audit,
  file.path(result_dir, "major_ambiguous_cluster_blind_audit.tsv"),
  sep = "\t",
  quote = TRUE,
  na = "NA"
)
message("Ambiguous-cluster blind marker audit completed")
