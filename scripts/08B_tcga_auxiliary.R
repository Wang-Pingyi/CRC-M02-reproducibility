#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 08B_tcga_auxiliary.R PROJECT_DIR RUN_ID")
project <- normalizePath(args[1], mustWork = TRUE)
run_id <- args[2]
source(file.path(project, "scripts", "08B_helpers.R"))
set.seed(20260729)
paths <- stage8b_paths(project, run_id)

required <- c("data.table", "edgeR", "sandwich", "lmtest", "survival", "CMScaller")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

gdc <- read.delim(file.path(paths$source, "gdc_tcga_coad_file_sample_metadata.tsv"), check.names = FALSE)
sample_clin <- read.delim(file.path(paths$source, "cbioportal_coadread_pancan_sample.tsv"),
                          comment.char = "#", check.names = FALSE)
patient_clin <- read.delim(file.path(paths$source, "cbioportal_coadread_pancan_patient.tsv"),
                           comment.char = "#", check.names = FALSE)
patient_clin <- patient_clin[patient_clin$CANCER_TYPE_ACRONYM == "COAD", , drop = FALSE]
sample_clin <- sample_clin[sample_clin$PATIENT_ID %in% patient_clin$PATIENT_ID, , drop = FALSE]

raw_dir <- file.path(project, "data_raw", "TCGA-COAD")
gdc$path <- file.path(raw_dir, gdc$local_file_name)
if (any(!file.exists(gdc$path))) stop("TCGA raw count files absent: ", sum(!file.exists(gdc$path)))
gdc <- gdc[order(gdc$patient_id, gdc$sample_type, gdc$sample_barcode, gdc$gdc_file_id), ]
gdc$selected <- !duplicated(paste(gdc$patient_id, gdc$sample_type))
selection_audit <- gdc
stage8b_write_tsv(selection_audit, file.path(paths$result, "TCGA_sample_selection_audit.tsv"))
gdc <- gdc[gdc$selected, ]

message("Reading ", nrow(gdc), " selected TCGA count files")
first <- data.table::fread(gdc$path[1], skip = 1, select = c("gene_id", "gene_name", "gene_type", "unstranded"))
keep <- !startsWith(first$gene_id, "N_")
gene_info <- as.data.frame(first[keep, c("gene_id", "gene_name", "gene_type")])
counts <- matrix(0, nrow = nrow(gene_info), ncol = nrow(gdc),
                 dimnames = list(gene_info$gene_id, gdc$sample_barcode))
counts[, 1] <- first$unstranded[keep]
if (nrow(gdc) > 1L) {
  for (j in 2:nrow(gdc)) {
    z <- data.table::fread(gdc$path[j], skip = 1, select = c("gene_id", "unstranded"))
    z <- z[!startsWith(gene_id, "N_")]
    if (!identical(z$gene_id, gene_info$gene_id)) stop("Gene order mismatch: ", gdc$path[j])
    counts[, j] <- z$unstranded
    if (j %% 50L == 0L) message("TCGA files read: ", j, "/", nrow(gdc))
  }
}
storage.mode(counts) <- "integer"

y <- edgeR::DGEList(counts = counts)
keep_expr <- edgeR::filterByExpr(y, group = gdc$sample_type)
y <- edgeR::calcNormFactors(y[keep_expr, , keep.lib.sizes = FALSE], method = "TMM")
logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
symbols <- toupper(gene_info$gene_name[match(rownames(logcpm), gene_info$gene_id)])
valid_symbol <- !is.na(symbols) & nzchar(symbols)
symbol_logcpm <- rowsum(logcpm[valid_symbol, , drop = FALSE], group = symbols[valid_symbol], reorder = FALSE)
symbol_n <- as.numeric(table(factor(symbols[valid_symbol], levels = rownames(symbol_logcpm))))
symbol_logcpm <- symbol_logcpm / symbol_n

locked <- stage8b_locked(paths)
score_list <- list()
coverage <- list()
for (module in locked$candidates$module_id) {
  expected <- locked$membership$gene[locked$membership$module_id == module]
  genes <- intersect(expected, rownames(symbol_logcpm))
  gz <- t(scale(t(symbol_logcpm[genes, , drop = FALSE])))
  gz <- gz[apply(gz, 1, function(z) all(is.finite(z))), , drop = FALSE]
  score <- as.numeric(scale(colMeans(gz)))
  score_list[[module]] <- data.frame(
    sample_barcode = colnames(symbol_logcpm), module_id = module,
    module_score = score, mapped_genes = nrow(gz), stringsAsFactors = FALSE
  )
  coverage[[module]] <- data.frame(accession = "TCGA-COAD", module_id = module,
    mapped_genes = nrow(gz), locked_genes = length(expected), coverage = nrow(gz) / length(expected))
}
scores <- do.call(rbind, score_list)
scores <- merge(scores, gdc[, c("sample_barcode", "patient_id", "sample_type", "tissue_type", "gdc_file_id")],
                by = "sample_barcode", all.x = TRUE, sort = FALSE)
stage8b_write_tsv(do.call(rbind, coverage), file.path(paths$result, "TCGA_module_mapping_coverage.tsv"))

# cBioPortal sample identifiers omit the vial character used by GDC.
scores$cbio_sample_id <- substr(scores$sample_barcode, 1, 15)
scores <- merge(scores, sample_clin[, c("PATIENT_ID", "SAMPLE_ID", "MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE")],
                by.x = c("patient_id", "cbio_sample_id"), by.y = c("PATIENT_ID", "SAMPLE_ID"),
                all.x = TRUE, sort = FALSE)
scores <- merge(scores, patient_clin, by.x = "patient_id", by.y = "PATIENT_ID", all.x = TRUE, sort = FALSE)
stage8b_write_tsv(scores, file.path(paths$result, "TCGA_module_scores_clinical.tsv"))

result <- list()
for (module in unique(scores$module_id)) {
  z <- scores[scores$module_id == module, ]
  # Paired primary tumor - solid tissue normal.
  p <- z[z$sample_type %in% c("Primary Tumor", "Solid Tissue Normal"), ]
  p$analysis_group <- ifelse(p$sample_type == "Primary Tumor", "cancer", "normal")
  pair_fit <- try(stage8b_fit_paired(transform(p, donor_id = patient_id), "normal", "cancer"), silent = TRUE)
  if (!inherits(pair_fit, "try-error")) {
    result[[length(result) + 1L]] <- stage8b_effect_row("TCGA-COAD", module, "cancer_vs_normal_paired",
      "auxiliary", pair_fit, "paired_patient_difference", "Auxiliary TCGA comparison")
  }
}

# Deterministic one-primary-tumor-per-patient table.
tumor <- scores[scores$sample_type == "Primary Tumor", ]
tumor <- tumor[!duplicated(paste(tumor$patient_id, tumor$module_id)), ]
stage_roman <- toupper(as.character(tumor$AJCC_PATHOLOGIC_TUMOR_STAGE))
tumor$stage_numeric <- ifelse(grepl("STAGE IV", stage_roman), 4,
  ifelse(grepl("STAGE III", stage_roman), 3,
    ifelse(grepl("STAGE II", stage_roman), 2,
      ifelse(grepl("STAGE I", stage_roman), 1, NA_real_))))
tumor$msi_group <- ifelse(tumor$MSI_SCORE_MANTIS > 0.6, "MSI",
                          ifelse(tumor$MSI_SCORE_MANTIS < 0.4, "MSS", "indeterminate"))

for (module in unique(tumor$module_id)) {
  z <- tumor[tumor$module_id == module, ]
  zs <- z[is.finite(z$stage_numeric), ]
  if (nrow(zs) >= 20L && length(unique(zs$stage_numeric)) >= 3L) {
    fit <- lm(module_score ~ stage_numeric, data = zs)
    cf <- stage8b_robust_coef(fit, "stage_numeric")
    result[[length(result) + 1L]] <- stage8b_effect_row("TCGA-COAD", module, "pathologic_stage_ordinal",
      "auxiliary", cf, "HC3_robust_ordinal", "Primary tumors only")
  }
  zm <- z[z$msi_group %in% c("MSS", "MSI"), ]
  if (sum(zm$msi_group == "MSI") >= 5L && sum(zm$msi_group == "MSS") >= 5L) {
    zm$analysis_group <- zm$msi_group
    cf <- stage8b_fit_binary(zm, "MSS", "MSI")
    result[[length(result) + 1L]] <- stage8b_effect_row("TCGA-COAD", module, "MSI_vs_MSS",
      "auxiliary", cf, "HC3_robust_binary", "MANTIS; indeterminate excluded")
  }
}

# CMS classification from raw tumor counts, once per selected primary tumor.
tumor_samples <- unique(scores$sample_barcode[scores$sample_type == "Primary Tumor"])
cms_counts <- counts[, tumor_samples, drop = FALSE]
cms_symbols <- toupper(gene_info$gene_name)
cms_keep <- !is.na(cms_symbols) & nzchar(cms_symbols)
cms_counts <- rowsum(cms_counts[cms_keep, , drop = FALSE], group = cms_symbols[cms_keep], reorder = FALSE)
cms <- CMScaller::CMScaller(cms_counts, rowNames = "symbol", RNAseq = TRUE,
                            nPerm = 1000, seed = 20260729, FDR = 0.05,
                            doPlot = FALSE, verbose = TRUE)
cms$sample_barcode <- rownames(cms)
stage8b_write_tsv(cms, file.path(paths$result, "TCGA_CMScaller_calls.tsv"))
tumor <- merge(tumor, cms[, c("sample_barcode", "prediction", "p.value", "FDR")],
               by = "sample_barcode", all.x = TRUE, sort = FALSE)
for (module in unique(tumor$module_id)) {
  z <- tumor[tumor$module_id == module & !is.na(tumor$prediction), ]
  z$prediction <- factor(z$prediction)
  if (nlevels(z$prediction) < 2L) next
  fit <- lm(module_score ~ prediction, data = z)
  vc <- sandwich::vcovHC(fit, type = "HC3")
  tab <- lmtest::coeftest(fit, vcov. = vc)
  for (term in setdiff(rownames(tab), "(Intercept)")) {
    cf <- stage8b_robust_coef(fit, term)
    result[[length(result) + 1L]] <- stage8b_effect_row("TCGA-COAD", module,
      paste0("CMS_", sub("^prediction", "", term), "_vs_", levels(z$prediction)[1]),
      "auxiliary", cf, "HC3_robust_CMS", "Confident CMScaller predictions only")
  }
}

# Exploratory overall-survival association.
for (module in unique(tumor$module_id)) {
  z <- tumor[tumor$module_id == module, ]
  z$event <- as.integer(grepl("^1:", z$OS_STATUS))
  z$time <- suppressWarnings(as.numeric(z$OS_MONTHS))
  z$age <- suppressWarnings(as.numeric(z$AGE))
  z$sex <- factor(z$SEX)
  z <- z[is.finite(z$time) & z$time > 0 & !is.na(z$event) & is.finite(z$age) &
           is.finite(z$stage_numeric) & !is.na(z$sex), ]
  if (nrow(z) < 50L || sum(z$event) < 15L) next
  fit <- survival::coxph(survival::Surv(time, event) ~ module_score + age + sex + stage_numeric,
                        data = z, ties = "efron")
  sm <- summary(fit)
  beta <- sm$coefficients["module_score", "coef"]
  se <- sm$coefficients["module_score", "se(coef)"]
  pval <- sm$coefficients["module_score", "Pr(>|z|)"]
  ph <- survival::cox.zph(fit)
  result[[length(result) + 1L]] <- data.frame(
    accession = "TCGA-COAD", module_id = module, endpoint = "overall_survival",
    analysis_set = "auxiliary", model = "Cox_age_sex_stage_adjusted",
    effect = beta, standard_error = se, ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se,
    p_value = pval, n_units = nrow(z), residual_df = NA_real_,
    note = paste0("log-HR per 1-SD; PH_p=", signif(ph$table["module_score", "p"], 4)),
    stringsAsFactors = FALSE
  )
}

tcga_results <- do.call(rbind, result)
tcga_results <- stage8b_add_fdr(tcga_results)
stage8b_write_tsv(tcga_results, file.path(paths$result, "TCGA_auxiliary_results.tsv"))
saveRDS(list(gdc = gdc, counts = counts, logcpm = logcpm, scores = scores, cms = cms),
        file.path(paths$result, "TCGA_stage8B_intermediate.rds"), compress = FALSE)
