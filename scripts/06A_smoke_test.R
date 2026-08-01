#!/usr/bin/env Rscript

set.seed(20260728)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1L]] else getwd()

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(limma)
})

obj <- readRDS(file.path(
  project_dir, "objects", "GSE201348_5C_epithelial_annotated_CNV.rds"
))
meta <- obj[[]]
selected <- which(meta$epithelial_state == "Absorptive")
meta <- meta[selected, , drop = FALSE]
counts <- LayerData(obj, assay = "RNA", layer = "counts")[, selected, drop = FALSE]
meta$stage <- factor(
  c(normal = "normal", adenoma_polyp = "adenoma", cancer = "cancer")[
    meta$lesion_stage
  ],
  levels = c("normal", "adenoma", "cancer")
)
meta$fap_binary <- factor(
  ifelse(meta$sporadic_or_FAP == "FAP", "FAP", "nonFAP"),
  levels = c("nonFAP", "FAP")
)
group_id <- paste(meta$donor_id, meta$stage, sep = "||")
group_levels <- unique(group_id)
membership <- sparseMatrix(
  i = seq_along(group_id),
  j = match(group_id, group_levels),
  x = 1,
  dims = c(length(group_id), length(group_levels))
)
aggregated <- counts %*% membership
colnames(aggregated) <- group_levels

group_meta <- do.call(rbind, lapply(group_levels, function(group) {
  idx <- which(group_id == group)
  data.frame(
    group_id = group,
    donor_id = unique(meta$donor_id[idx]),
    stage = unique(meta$stage[idx]),
    fap_binary = unique(meta$fap_binary[idx]),
    n_cells = length(idx),
    library_size = sum(aggregated[, group]),
    stringsAsFactors = FALSE
  )
}))
keep_samples <- group_meta$n_cells >= 20 & group_meta$library_size >= 10000
group_meta <- group_meta[keep_samples, , drop = FALSE]
aggregated <- aggregated[, group_meta$group_id, drop = FALSE]
group_meta$stage <- factor(group_meta$stage, levels = c("normal", "adenoma", "cancer"))
group_meta$fap_binary <- factor(group_meta$fap_binary, levels = c("nonFAP", "FAP"))

dge <- DGEList(aggregated)
keep_genes <- filterByExpr(dge, group = group_meta$stage, min.count = 10, min.total.count = 15)
dge <- calcNormFactors(dge[keep_genes, , keep.lib.sizes = FALSE])
design <- model.matrix(~ 0 + stage + fap_binary, group_meta)
if (qr(design)$rank < ncol(design)) {
  design <- model.matrix(~ 0 + stage, group_meta)
}
colnames(design) <- sub("^stage", "", colnames(design))
voom_first <- voom(dge, design, plot = FALSE)
correlation <- duplicateCorrelation(voom_first, design, block = group_meta$donor_id)$consensus.correlation
voom_final <- voom(
  dge, design, plot = FALSE, block = group_meta$donor_id, correlation = correlation
)
fit <- lmFit(
  voom_final, design, block = group_meta$donor_id, correlation = correlation
)
contrast <- makeContrasts(
  adenoma_vs_normal = adenoma - normal,
  cancer_vs_adenoma = cancer - adenoma,
  cancer_vs_normal = cancer - normal,
  levels = design
)
fit <- eBayes(contrasts.fit(fit, contrast), robust = TRUE)
stopifnot(all(is.finite(fit$coefficients)))
cat(
  "SMOKE_TEST_OK",
  "state=Absorptive",
  paste0("pseudobulks=", nrow(group_meta)),
  paste0("genes=", nrow(dge)),
  paste0("correlation=", signif(correlation, 4)),
  sep = "\t"
)
cat("\n")
