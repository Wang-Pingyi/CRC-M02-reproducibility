#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  args0 <- commandArgs(trailingOnly = TRUE)
})

arg_value <- function(flag, default = NULL) {
  hit <- which(args0 == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args0)) stop("Missing value for ", flag)
  args0[hit[length(hit)] + 1L]
}

root <- normalizePath(arg_value("--root", "."), mustWork = TRUE)
smoke <- "--smoke" %in% args0
project_lib <- file.path(root, "environment", "R-library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

set.seed(42)
options(stringsAsFactors = FALSE)

sample_manifest_path <- file.path(root, "data", "metadata", "stage10c2_sc_patient_sample_manifest.tsv")
lock_path <- file.path(root, "results", "stage10c", "STAGE10C_LOCK_MANIFEST.tsv")
plan_path <- file.path(root, "results", "stage10c3", "STAGE10C3_ANALYSIS_PLAN_LOCKED.md")
matrix_dir <- file.path(root, "data_processed", "stage10c2_sc", "figshare_29925404", "v1", "derived")

stopifnot(file.exists(sample_manifest_path), file.exists(lock_path), file.exists(plan_path), dir.exists(matrix_dir))
manifest <- fread(sample_manifest_path, na.strings = c("NA", ""))
lock <- fread(lock_path)
stopifnot(nrow(lock) == 36L, unique(lock$bundle_sha256) == "d4e34472243b0259650aca3123a5df7e767e9a408a409cef08861937f2360a30")
m02_genes <- lock$gene

panels <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  immune = c("PTPRC", "LST1", "TYROBP", "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "FCER1G"),
  stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1", "COL6A2"),
  endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG"),
  mast = c("TPSAB1", "TPSB2", "MS4A2", "KIT"),
  erythroid = c("HBA1", "HBA2", "HBB")
)
stem_panel <- c("LGR5", "OLFM4", "SMOC2", "PROM1", "LRIG1", "SOX9")
diff_panel <- c("KRT20", "CA1", "GUCA2A", "FABP1", "MUC2", "TFF3", "CHGA", "POU2F3")
cycle_panel <- c("MKI67", "TOP2A", "UBE2C")
annotation_genes <- unique(c(unlist(panels), stem_panel, diff_panel, cycle_panel))
if (length(intersect(annotation_genes, m02_genes))) stop("Annotation panel overlaps frozen M02")

samples <- manifest$sample_id
if (smoke) samples <- "P2_L"
if (!all(samples %in% manifest$sample_id)) stop("Unknown sample requested")

out_results <- if (smoke) file.path(root, "cache", "stage10c3_smoke", "results") else file.path(root, "results", "stage10c3")
out_objects <- if (smoke) file.path(root, "cache", "stage10c3_smoke", "objects") else file.path(root, "objects", "stage10c3")
dir.create(out_results, recursive = TRUE, showWarnings = FALSE)
dir.create(out_objects, recursive = TRUE, showWarnings = FALSE)

get_data <- function(obj) {
  tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
           error = function(e) GetAssayData(obj, assay = "RNA", slot = "data"))
}

read_count_csv <- function(path, sample_id) {
  dt <- fread(path, check.names = FALSE, showProgress = FALSE, nThread = min(8L, parallel::detectCores()))
  genes <- dt[[1L]]
  if (anyDuplicated(genes)) stop("Duplicate genes in ", sample_id)
  barcodes <- names(dt)[-1L]
  dense <- as.matrix(dt[, -1L, with = FALSE])
  storage.mode(dense) <- "numeric"
  if (any(!is.finite(dense)) || any(dense < 0) || any(dense != floor(dense))) stop("Invalid counts in ", sample_id)
  rownames(dense) <- genes
  colnames(dense) <- paste(sample_id, barcodes, sep = "__")
  sparse <- as(dense, "dgCMatrix")
  rm(dt, dense)
  gc(verbose = FALSE)
  sparse
}

qc_metrics <- function(counts) {
  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt <- grepl("^MT-", rownames(counts))
  rb <- grepl("^RP[SL]", rownames(counts))
  pct_mt <- if (any(mt)) 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(n_count, 1) else rep(0, ncol(counts))
  pct_rb <- if (any(rb)) 100 * Matrix::colSums(counts[rb, , drop = FALSE]) / pmax(n_count, 1) else rep(0, ncol(counts))
  data.table(cell_id = colnames(counts), nCount = as.numeric(n_count), nFeature = as.integer(n_feature),
             percent_mt = as.numeric(pct_mt), percent_ribo = as.numeric(pct_rb))
}

run_clusters <- function(counts, resolution = 0.4, seed = 42L) {
  if (ncol(counts) < 30L) return(setNames(rep("0", ncol(counts)), colnames(counts)))
  obj <- CreateSeuratObject(counts = counts, project = "stage10c3", min.cells = 0, min.features = 0)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000L, nrow(obj) - 1L), verbose = FALSE)
  vf <- setdiff(VariableFeatures(obj), m02_genes)
  vf <- vf[vf %in% rownames(obj)]
  if (length(vf) < 50L) stop("Too few non-M02 variable features")
  obj <- ScaleData(obj, features = vf, verbose = FALSE)
  npcs <- min(30L, length(vf) - 1L, ncol(obj) - 1L)
  set.seed(seed)
  obj <- RunPCA(obj, features = vf, npcs = npcs, verbose = FALSE, seed.use = seed)
  ndims <- min(20L, npcs)
  obj <- FindNeighbors(obj, dims = seq_len(ndims), verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, random.seed = seed, verbose = FALSE)
  setNames(as.character(Idents(obj)), colnames(obj))
}

cluster_panel_stats <- function(counts, clusters, panel_list) {
  obj <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  logdata <- get_data(obj)
  out <- list()
  for (cl in sort(unique(clusters))) {
    cells <- names(clusters)[clusters == cl]
    for (pn in names(panel_list)) {
      genes <- intersect(panel_list[[pn]], rownames(counts))
      if (!length(genes)) {
        score <- 0
        coherent <- 0L
      } else {
        score <- mean(as.numeric(Matrix::rowMeans(logdata[genes, cells, drop = FALSE])))
        frac <- as.numeric(Matrix::rowMeans(counts[genes, cells, drop = FALSE] > 0))
        coherent <- sum(frac >= 0.10)
      }
      out[[length(out) + 1L]] <- data.table(cluster = cl, panel = pn, score = score, coherent_genes = coherent, n_cells = length(cells))
    }
  }
  rbindlist(out)
}

annotate_broad <- function(counts, clusters) {
  stats <- cluster_panel_stats(counts, clusters, panels)
  wide_score <- dcast(stats, cluster ~ panel, value.var = "score", fill = 0)
  wide_coh <- dcast(stats, cluster ~ panel, value.var = "coherent_genes", fill = 0)
  labels <- lapply(seq_len(nrow(wide_score)), function(i) {
    cl <- wide_score$cluster[i]
    s <- as.list(wide_score[i, -1L])
    coh <- as.list(wide_coh[cluster == cl, -1L])
    exclusion_max <- max(s$immune, s$stromal, s$endothelial)
    if (coh$epithelial >= 2L && s$epithelial > exclusion_max) return(data.table(cluster = cl, broad_type = "epithelial", ambient_flag = "none"))
    if (coh$epithelial >= 2L && s$epithelial <= exclusion_max) return(data.table(cluster = cl, broad_type = "uncertain_mixed", ambient_flag = "lineage_conflict"))
    candidates <- c(immune = s$immune, stromal = s$stromal, endothelial = s$endothelial, mast = s$mast, erythroid = s$erythroid)
    winner <- names(which.max(candidates))
    required <- if (winner %in% c("mast", "erythroid")) 1L else 2L
    if (coh[[winner]] >= required) data.table(cluster = cl, broad_type = winner, ambient_flag = "none")
    else data.table(cluster = cl, broad_type = "uncertain", ambient_flag = "insufficient_coherence")
  })
  list(labels = rbindlist(labels), stats = stats)
}

annotate_epithelial_state <- function(counts, broad, seed) {
  epi_cells <- names(broad)[broad == "epithelial"]
  state <- setNames(rep(NA_character_, length(broad)), names(broad))
  subcluster <- setNames(rep(NA_character_, length(broad)), names(broad))
  if (length(epi_cells) < 30L) return(list(state = state, subcluster = subcluster, stats = data.table()))
  epi_counts <- counts[, epi_cells, drop = FALSE]
  epi_cl <- run_clusters(epi_counts, resolution = 0.35, seed = seed + 1000L)
  subcluster[names(epi_cl)] <- epi_cl
  epi_panels <- list(stem = stem_panel, differentiation = diff_panel, cycling = cycle_panel)
  stats <- cluster_panel_stats(epi_counts, epi_cl, epi_panels)
  ws <- dcast(stats, cluster ~ panel, value.var = "score", fill = 0)
  wc <- dcast(stats, cluster ~ panel, value.var = "coherent_genes", fill = 0)
  cl_state <- lapply(seq_len(nrow(ws)), function(i) {
    cl <- ws$cluster[i]
    coh <- wc[cluster == cl]
    stem_ok <- coh$stem >= 2L && ws$stem[i] > ws$differentiation[i]
    cycling_ok <- coh$cycling >= 2L
    label <- if (stem_ok && cycling_ok) "Stem/progenitor_cycling" else if (stem_ok) "Stem/progenitor" else if (cycling_ok) "Cycling" else if (ws$differentiation[i] >= ws$stem[i]) "Differentiated" else "Other_epithelial"
    data.table(cluster = cl, epithelial_state = label)
  })
  cl_state <- rbindlist(cl_state)
  state[names(epi_cl)] <- cl_state$epithelial_state[match(epi_cl, cl_state$cluster)]
  list(state = state, subcluster = subcluster, stats = stats)
}

process_mode <- function(counts, sample_id, mode, sample_index) {
  metrics <- qc_metrics(counts)
  primary_mask <- metrics$nCount > 0 & metrics$nFeature >= 1000L & metrics$nFeature <= 8000L & metrics$percent_mt <= 50
  deposit_mask <- metrics$nCount > 0 & metrics$percent_mt <= 50
  keep <- if (mode == "primary_qc") primary_mask else deposit_mask
  selected <- counts[, metrics$cell_id[keep], drop = FALSE]
  predicted_doublets <- 0L
  if (mode == "primary_qc") {
    sce <- SingleCellExperiment(assays = list(counts = selected))
    set.seed(42L + sample_index)
    sce <- scDblFinder(sce, BPPARAM = SerialParam(), verbose = FALSE)
    dbl <- as.character(colData(sce)$scDblFinder.class)
    predicted_doublets <- sum(dbl == "doublet")
    selected <- selected[, colnames(sce)[dbl == "singlet"], drop = FALSE]
  }
  clusters <- run_clusters(selected, resolution = 0.4, seed = 42L + sample_index)
  broad_result <- annotate_broad(selected, clusters)
  broad_map <- setNames(broad_result$labels$broad_type, broad_result$labels$cluster)
  ambient_map <- setNames(broad_result$labels$ambient_flag, broad_result$labels$cluster)
  broad <- setNames(unname(broad_map[clusters]), names(clusters))
  ambient <- setNames(unname(ambient_map[clusters]), names(clusters))
  state_result <- annotate_epithelial_state(selected, broad, 42L + sample_index)
  cell_meta <- data.table(
    cell_id = colnames(selected), sample_id = sample_id, qc_mode = mode,
    broad_cluster = unname(clusters[colnames(selected)]),
    broad_type = unname(broad[colnames(selected)]),
    ambient_flag = unname(ambient[colnames(selected)]),
    epithelial_subcluster = unname(state_result$subcluster[colnames(selected)]),
    epithelial_state = unname(state_result$state[colnames(selected)])
  )
  cell_meta <- merge(cell_meta, metrics, by = "cell_id", all.x = TRUE, sort = FALSE)
  summary <- data.table(
    sample_id = sample_id, qc_mode = mode, deposited_cells = ncol(counts),
    cells_after_general_qc = sum(keep), predicted_doublets = predicted_doublets,
    cells_after_doublet_filter = ncol(selected), epithelial_cells = sum(cell_meta$broad_type == "epithelial"),
    stem_progenitor_cells = sum(grepl("^Stem/progenitor", cell_meta$epithelial_state), na.rm = TRUE),
    uncertain_mixed_cells = sum(cell_meta$broad_type == "uncertain_mixed"),
    epithelial_min25 = sum(cell_meta$broad_type == "epithelial") >= 25L,
    epithelial_min50 = sum(cell_meta$broad_type == "epithelial") >= 50L,
    epithelial_min100 = sum(cell_meta$broad_type == "epithelial") >= 100L,
    stem_min20 = sum(grepl("^Stem/progenitor", cell_meta$epithelial_state), na.rm = TRUE) >= 20L
  )
  list(counts = selected, cell_meta = cell_meta, summary = summary,
       broad_cluster_stats = broad_result$stats, epithelial_cluster_stats = state_result$stats)
}

qc_rows <- list()
cluster_rows <- list()
started <- Sys.time()
modes <- if (smoke) "primary_qc" else c("primary_qc", "deposit_preserving")

for (i in seq_along(samples)) {
  sample_id <- samples[i]
  message(sprintf("[%s] reading %s", format(Sys.time(), "%F %T"), sample_id))
  counts <- read_count_csv(file.path(matrix_dir, paste0(sample_id, ".csv")), sample_id)
  for (mode in modes) {
    message(sprintf("[%s] annotating %s / %s", format(Sys.time(), "%F %T"), sample_id, mode))
    obj <- process_mode(counts, sample_id, mode, match(sample_id, manifest$sample_id))
    saveRDS(obj[c("counts", "cell_meta")], file.path(out_objects, paste0(sample_id, "__", mode, ".rds")), compress = "xz")
    qc_rows[[length(qc_rows) + 1L]] <- obj$summary
    if (nrow(obj$broad_cluster_stats)) {
      tmp <- copy(obj$broad_cluster_stats); tmp[, `:=`(sample_id = sample_id, qc_mode = mode, annotation_level = "broad")]
      cluster_rows[[length(cluster_rows) + 1L]] <- tmp
    }
    if (nrow(obj$epithelial_cluster_stats)) {
      tmp <- copy(obj$epithelial_cluster_stats); tmp[, `:=`(sample_id = sample_id, qc_mode = mode, annotation_level = "epithelial_state")]
      cluster_rows[[length(cluster_rows) + 1L]] <- tmp
    }
    rm(obj); gc(verbose = FALSE)
  }
  rm(counts); gc(verbose = FALSE)
}

qc <- rbindlist(qc_rows, fill = TRUE)
cluster_stats <- rbindlist(cluster_rows, fill = TRUE)
fwrite(qc, file.path(out_results, "STAGE10C3_QC_SOURCE.tsv"), sep = "\t", na = "NA")
fwrite(cluster_stats, file.path(out_results, "STAGE10C3_ANNOTATION_SOURCE.tsv"), sep = "\t", na = "NA")

package_names <- c("R", "data.table", "Matrix", "Seurat", "SingleCellExperiment", "scDblFinder", "BiocParallel")
package_versions <- c(as.character(getRversion()), vapply(package_names[-1L], function(x) as.character(packageVersion(x)), character(1)))
pkg_text <- paste(sprintf("- %s: `%s`", package_names, package_versions), collapse = "\n")
qc_table <- paste(capture.output(print(qc)), collapse = "\n")
qc_md <- paste0(
  "# Stage 10C3 QC and M02-blind Annotation\n\n",
  "Generated: ", format(Sys.time(), "%F %T %Z"), "\n\n",
  "This QC and broad annotation was completed before M02 expression scoring. ",
  "No frozen M02 gene was used as a highly variable feature, marker, cluster name or exclusion criterion.\n\n",
  "## Rules applied\n\n",
  "- Each sample was read and checked independently.\n",
  "- Primary QC: 1,000-8,000 detected genes, mitochondrial percentage <=50%, then per-sample scDblFinder.\n",
  "- Deposit-preserving sensitivity: valid deposited cells with mitochondrial percentage <=50%.\n",
  "- Epithelial annotation: independent fixed positive and exclusion panels with cluster-level coherence.\n",
  "- Stem/progenitor: independent six-gene positive panel against fixed differentiation markers.\n",
  "- Empty-droplet ambient modeling was not claimed because unfiltered droplets were unavailable. Mixed-lineage clusters were diagnosed and excluded.\n",
  "- Published cell-level annotation sensitivity is unavailable because no such file was deposited.\n\n",
  "## Sample summary\n\n```text\n", qc_table, "\n```\n\n",
  "## Software\n\n", pkg_text, "\n\n",
  "Random seed: `42`; per-sample scDblFinder seed: `42 + frozen sample index`.\n"
)
writeLines(qc_md, file.path(out_results, "STAGE10C3_QC.md"), useBytes = TRUE)
writeLines(c(
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("plan_mtime=", format(file.mtime(plan_path), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("samples=", paste(samples, collapse = ",")),
  paste0("runtime_seconds=", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)),
  "m02_expression_scored=FALSE"
), file.path(out_results, "QC_PREM02.SUCCESS"), useBytes = TRUE)

message("QC and M02-blind annotation complete")
