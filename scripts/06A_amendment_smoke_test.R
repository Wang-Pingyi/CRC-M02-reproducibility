#!/usr/bin/env Rscript

set.seed(20260728)
args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) args[[1L]] else getwd()

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
})

pb <- readRDS(file.path(
  project_dir, "objects", "GSE201348_6A_epithelial_pseudobulk.rds"
))
counts <- pb$counts
meta <- pb$metadata
stopifnot(identical(colnames(counts), rownames(meta)))

state_meta <- meta[
  meta$eligible &
    meta$epithelial_state == "Absorptive" &
    meta$fap_binary == "FAP" &
    meta$stage %in% c("normal", "adenoma"),
  ,
  drop = FALSE
]
paired <- intersect(
  state_meta$donor_id[state_meta$stage == "normal"],
  state_meta$donor_id[state_meta$stage == "adenoma"]
)
stopifnot(length(unique(paired)) >= 3L)
state_meta <- state_meta[state_meta$donor_id %in% paired, , drop = FALSE]
state_counts <- counts[, state_meta$pseudobulk_id, drop = FALSE]

dge <- DGEList(state_counts)
keep <- filterByExpr(dge, group = state_meta$stage, min.count = 10, min.total.count = 15)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
state_meta$donor_id <- factor(state_meta$donor_id)
state_meta$stage <- factor(state_meta$stage, levels = c("normal", "adenoma"))
paired_design <- model.matrix(~ donor_id + stage, state_meta)
paired_voom <- voom(dge, paired_design, plot = FALSE)
paired_fit <- eBayes(lmFit(paired_voom, paired_design), robust = TRUE)
stopifnot("stageadenoma" %in% colnames(paired_fit$coefficients))

all_meta <- meta[
  meta$eligible & meta$epithelial_state == "Absorptive",
  ,
  drop = FALSE
]
all_counts <- counts[, all_meta$pseudobulk_id, drop = FALSE]
all_dge <- DGEList(all_counts)
all_keep <- filterByExpr(
  all_dge, group = rep(1L, nrow(all_meta)), min.count = 10, min.total.count = 15
)
all_dge <- calcNormFactors(all_dge[all_keep, , keep.lib.sizes = FALSE])
logcpm <- cpm(all_dge, log = TRUE, prior.count = 1)
variability <- sort(apply(logcpm, 1, mad), decreasing = TRUE)
genes <- names(head(variability[is.finite(variability) & variability > 0], 200L))
scaled <- t(scale(t(logcpm[genes, , drop = FALSE])))
scaled[!is.finite(scaled)] <- 0
correlation <- cor(t(scaled), use = "pairwise.complete.obs")
correlation[!is.finite(correlation)] <- 0
diag(correlation) <- 1
membership <- cutree(hclust(as.dist(1 - correlation), method = "average"), k = 5L)
score_matrix <- do.call(rbind, lapply(sort(unique(membership)), function(module) {
  colMeans(scaled[membership == module, , drop = FALSE])
}))
rownames(score_matrix) <- paste0("smoke_module_", seq_len(nrow(score_matrix)))

all_meta$stage <- factor(all_meta$stage, levels = c("normal", "adenoma", "cancer"))
all_meta$fap_binary <- factor(all_meta$fap_binary, levels = c("nonFAP", "FAP"))
module_design <- model.matrix(~ 0 + stage + fap_binary, all_meta)
colnames(module_design) <- sub("^stage", "", colnames(module_design))
duplicate <- duplicateCorrelation(
  score_matrix, module_design, block = all_meta$donor_id
)
module_fit <- lmFit(
  score_matrix,
  module_design,
  block = all_meta$donor_id,
  correlation = duplicate$consensus.correlation
)
module_contrast <- makeContrasts(
  adenoma_vs_normal = adenoma - normal,
  cancer_vs_adenoma = cancer - adenoma,
  cancer_vs_normal = cancer - normal,
  levels = module_design
)
module_fit <- eBayes(contrasts.fit(module_fit, module_contrast), robust = TRUE)
stopifnot(all(is.finite(module_fit$coefficients)))

cat(
  "AMENDMENT_SMOKE_TEST_OK",
  paste0("paired_donors=", length(unique(paired))),
  paste0("paired_genes=", nrow(dge)),
  paste0("modules=", nrow(score_matrix)),
  sep = "\t"
)
cat("\n")
