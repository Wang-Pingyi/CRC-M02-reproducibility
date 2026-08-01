#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(SingleCellExperiment)
  library(scDblFinder)
})

set.seed(20260726)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else getwd()
raw_dir <- file.path(project_dir, "data_raw", "GSE201348")
result_dir <- file.path(project_dir, "results", "05A_qc_pilot")
figure_dir <- file.path(project_dir, "figures", "05A_qc_pilot")
source_dir <- file.path(result_dir, "source_data")
log_dir <- file.path(project_dir, "logs", "05A_qc_pilot")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pilot <- data.frame(
  sample_id = c("GSM6061709", "GSM6061647", "GSM6061646", "GSM6061701", "GSM6061702"),
  donor_id = c("B001", "A001", "A001", "A022", "CRC1_8810"),
  condition = c("normal_mucosa", "unaffected_mucosa", "polyp_or_adenoma",
                "polyp_or_adenoma", "cancer"),
  histology = c("Normal", "Normal", "TA", "Polyp", "Adenocarcinoma"),
  fap_status = c("non-FAP", "FAP", "FAP", "non-FAP", "sporadic"),
  approved_min_count = c(NA, NA, NA, 100, NA),
  approved_min_feature = c(NA, NA, NA, 100, NA),
  qc_risk = c("standard", "standard", "standard", "low_depth_high_risk", "standard"),
  rationale = c(
    "healthy non-FAP normal mucosa",
    "FAP unaffected mucosa",
    "FAP tubular adenoma",
    "only audited non-FAP polyp donor",
    "sporadic adenocarcinoma with unambiguous pathology"
  ),
  stringsAsFactors = FALSE
)
fwrite(pilot, file.path(source_dir, "pilot_sample_selection.tsv"), sep = "\t", na = "NA")

locate_one <- function(sample_id, suffix) {
  hits <- list.files(raw_dir, pattern = paste0("^", sample_id, ".*_", suffix, "$"),
                     full.names = TRUE)
  if (length(hits) != 1L) {
    stop("Expected exactly one ", suffix, " file for ", sample_id, "; found ", length(hits))
  }
  hits[[1]]
}

read_10x_sample <- function(sample_id) {
  matrix_file <- locate_one(sample_id, "matrix\\.mtx\\.gz")
  feature_file <- locate_one(sample_id, "features\\.tsv\\.gz")
  barcode_file <- locate_one(sample_id, "barcodes\\.tsv\\.gz")
  features <- fread(feature_file, header = FALSE, data.table = FALSE)
  barcodes <- fread(barcode_file, header = FALSE, data.table = FALSE)
  counts <- readMM(gzfile(matrix_file))
  counts <- as(counts, "dgCMatrix")
  if (nrow(counts) != nrow(features) || ncol(counts) != nrow(barcodes)) {
    stop("Dimension mismatch for ", sample_id)
  }
  gene_names <- make.unique(as.character(features[[min(2L, ncol(features))]]))
  rownames(counts) <- gene_names
  colnames(counts) <- paste(sample_id, barcodes[[1]], sep = "_")
  counts
}

robust_bounds <- function(x, lower_k = 3, upper_k = 4, nonnegative = TRUE) {
  med <- median(x, na.rm = TRUE)
  scale <- mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) scale <- IQR(x, na.rm = TRUE) / 1.349
  if (!is.finite(scale) || scale == 0) scale <- max(1, med * 0.1)
  lower <- med - lower_k * scale
  upper <- med + upper_k * scale
  if (nonnegative) lower <- max(0, lower)
  c(lower = lower, upper = upper, median = med, robust_sd = scale)
}

run_one <- function(row) {
  sample_id <- row$sample_id
  message("Processing ", sample_id)
  started <- Sys.time()
  gc(reset = TRUE)
  counts <- read_10x_sample(sample_id)
  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt_idx <- grepl("^MT-", rownames(counts), ignore.case = TRUE)
  ribo_idx <- grepl("^RP[SL][0-9A-Z-]+$", rownames(counts), ignore.case = TRUE)
  pct_mt <- if (any(mt_idx)) 100 * Matrix::colSums(counts[mt_idx, , drop = FALSE]) / pmax(n_count, 1) else 0
  pct_ribo <- if (any(ribo_idx)) 100 * Matrix::colSums(counts[ribo_idx, , drop = FALSE]) / pmax(n_count, 1) else 0

  count_b <- robust_bounds(log10(n_count + 1))
  feature_b <- robust_bounds(log10(n_feature + 1))
  mt_b <- robust_bounds(pct_mt, lower_k = Inf, upper_k = 3)
  ribo_b <- robust_bounds(pct_ribo, lower_k = Inf, upper_k = 4)

  low_count <- log10(n_count + 1) < count_b[["lower"]]
  high_count <- log10(n_count + 1) > count_b[["upper"]]
  low_feature <- log10(n_feature + 1) < feature_b[["lower"]]
  high_feature <- log10(n_feature + 1) > feature_b[["upper"]]
  high_mt <- pct_mt > mt_b[["upper"]]
  high_ribo <- pct_ribo > ribo_b[["upper"]]

  approved_low_count <- if (is.na(row$approved_min_count)) {
    rep(FALSE, length(n_count))
  } else {
    n_count < row$approved_min_count
  }
  approved_low_feature <- if (is.na(row$approved_min_feature)) {
    rep(FALSE, length(n_feature))
  } else {
    n_feature < row$approved_min_feature
  }
  low_quality <- low_count | low_feature | approved_low_count | approved_low_feature
  doublet_input <- !low_quality
  if (sum(doublet_input) < 100L) {
    stop("Too few cells remain for scDblFinder after prefiltering ", sample_id)
  }

  sce <- SingleCellExperiment(list(counts = counts[, doublet_input, drop = FALSE]))
  colData(sce)$sample_id <- sample_id
  sce <- scDblFinder(sce, samples = "sample_id", verbose = FALSE)
  doublet <- rep(FALSE, ncol(counts))
  doublet_score <- rep(NA_real_, ncol(counts))
  doublet_class <- rep("not_evaluated_low_quality", ncol(counts))
  doublet[doublet_input] <- colData(sce)$scDblFinder.class == "doublet"
  doublet_score[doublet_input] <- colData(sce)$scDblFinder.score
  doublet_class[doublet_input] <- as.character(colData(sce)$scDblFinder.class)

  exclusion_reason <- ifelse(
    low_quality, "low_quality",
    ifelse(high_count | high_feature, "abnormally_high_complexity",
           ifelse(high_mt, "high_mitochondrial_fraction",
                  ifelse(doublet, "doublet", "retained")))
  )

  metrics <- data.frame(
    sample_id = sample_id,
    cell_id = colnames(counts),
    nCount = as.numeric(n_count),
    nFeature = as.numeric(n_feature),
    percent_mt = as.numeric(pct_mt),
    percent_ribo = as.numeric(pct_ribo),
    scDblFinder_score = doublet_score,
    scDblFinder_class = doublet_class,
    flag_approved_min_count = approved_low_count,
    flag_approved_min_feature = approved_low_feature,
    flag_low_count = low_count,
    flag_high_count = high_count,
    flag_low_feature = low_feature,
    flag_high_feature = high_feature,
    flag_high_mt = high_mt,
    flag_high_ribo = high_ribo,
    exclusion_reason = exclusion_reason,
    retained = exclusion_reason == "retained",
    stringsAsFactors = FALSE
  )
  fwrite(metrics, file.path(source_dir, paste0(sample_id, "_cell_qc_metrics.tsv.gz")),
         sep = "\t", na = "NA")

  long <- rbind(
    data.frame(sample_id, metric = "log10_nCount", value = log10(n_count + 1)),
    data.frame(sample_id, metric = "log10_nFeature", value = log10(n_feature + 1)),
    data.frame(sample_id, metric = "percent_mt", value = pct_mt),
    data.frame(sample_id, metric = "percent_ribo", value = pct_ribo)
  )
  thresholds <- data.frame(
    sample_id = sample_id,
    metric = c("log10_nCount", "log10_nFeature", "percent_mt", "percent_ribo"),
    lower = c(count_b[["lower"]], feature_b[["lower"]], NA, NA),
    effective_lower = c(
      if (is.na(row$approved_min_count)) count_b[["lower"]] else
        max(count_b[["lower"]], log10(row$approved_min_count + 1)),
      if (is.na(row$approved_min_feature)) feature_b[["lower"]] else
        max(feature_b[["lower"]], log10(row$approved_min_feature + 1)),
      NA,
      NA
    ),
    upper = c(count_b[["upper"]], feature_b[["upper"]], mt_b[["upper"]], ribo_b[["upper"]]),
    median = c(count_b[["median"]], feature_b[["median"]], mt_b[["median"]], ribo_b[["median"]]),
    robust_sd = c(count_b[["robust_sd"]], feature_b[["robust_sd"]], mt_b[["robust_sd"]], ribo_b[["robust_sd"]]),
    filtering_role = c("exclusion", "exclusion", "exclusion", "diagnostic_only")
  )
  fwrite(long, file.path(source_dir, paste0(sample_id, "_distribution_source.tsv.gz")),
         sep = "\t", na = "NA")

  p <- ggplot(long, aes(x = value)) +
    geom_histogram(bins = 60, fill = "#3C78A8", color = "white", linewidth = 0.1) +
    geom_vline(data = thresholds, aes(xintercept = effective_lower), color = "#D55E00",
               linetype = 2, na.rm = TRUE) +
    geom_vline(data = thresholds, aes(xintercept = upper), color = "#D55E00",
               linetype = 2, na.rm = TRUE) +
    facet_wrap(~metric, scales = "free", ncol = 2) +
    labs(title = paste0(sample_id, " pre-QC distributions"),
         subtitle = "Dashed lines: sample-specific robust thresholds",
         x = NULL, y = "Cells") +
    theme_bw(base_size = 10)
  ggsave(file.path(figure_dir, paste0(sample_id, "_qc_distributions.pdf")), p,
         width = 8.5, height = 6.5)
  ggsave(file.path(figure_dir, paste0(sample_id, "_qc_distributions.png")), p,
         width = 8.5, height = 6.5, dpi = 180)

  reason_counts <- as.data.frame(table(metrics$exclusion_reason), stringsAsFactors = FALSE)
  names(reason_counts) <- c("exclusion_reason", "n_cells")
  reason_counts$sample_id <- sample_id
  reason_counts <- reason_counts[, c("sample_id", "exclusion_reason", "n_cells")]

  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  peak_gc_mb <- max(gc()[, "max used"] * c(0.000056, 0.000008))
  summary <- data.frame(
    sample_id = sample_id,
    donor_id = row$donor_id,
    condition = row$condition,
    histology = row$histology,
    fap_status = row$fap_status,
    cells_before = ncol(counts),
    cells_entering_doublet_detection = sum(doublet_input),
    cells_after = sum(metrics$retained),
    retained_percent = 100 * mean(metrics$retained),
    doublets = sum(doublet),
    low_quality_excluded = sum(low_quality),
    high_ribosomal_flagged = sum(high_ribo),
    ribosomal_only_excluded = 0L,
    qc_risk = row$qc_risk,
    elapsed_seconds = elapsed,
    approximate_peak_R_memory_mb = peak_gc_mb,
    input_sparse_object_mb = as.numeric(object.size(counts)) / 1024^2,
    ambient_rna_assessment = "limited: only filtered matrix supplied; no empty droplets",
    stringsAsFactors = FALSE
  )
  rm(sce, counts)
  gc()
  list(summary = summary, thresholds = thresholds, reasons = reason_counts)
}

results <- lapply(seq_len(nrow(pilot)), function(i) run_one(pilot[i, , drop = FALSE]))
summary_table <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
threshold_table <- rbindlist(lapply(results, `[[`, "thresholds"), fill = TRUE)
reason_table <- rbindlist(lapply(results, `[[`, "reasons"), fill = TRUE)
fwrite(summary_table, file.path(result_dir, "qc_pilot_summary.tsv"), sep = "\t", na = "NA")
fwrite(threshold_table, file.path(source_dir, "adaptive_thresholds.tsv"), sep = "\t", na = "NA")
fwrite(reason_table, file.path(source_dir, "exclusion_reason_counts.tsv"), sep = "\t", na = "NA")

output_bytes <- sum(file.info(list.files(c(result_dir, figure_dir), recursive = TRUE,
                                         full.names = TRUE))$size, na.rm = TRUE)
resource <- data.frame(
  run_started = format(min(file.info(list.files(source_dir, full.names = TRUE))$ctime), tz = "UTC"),
  run_finished = format(Sys.time(), tz = "UTC"),
  elapsed_seconds_total = sum(summary_table$elapsed_seconds),
  max_sample_peak_R_memory_mb = max(summary_table$approximate_peak_R_memory_mb),
  output_disk_mb = output_bytes / 1024^2,
  projected_72_sample_cpu_hours = sum(summary_table$elapsed_seconds) / nrow(summary_table) * 72 / 3600,
  projected_72_sample_output_gb = output_bytes / nrow(summary_table) * 72 / 1024^3
)
fwrite(resource, file.path(result_dir, "resource_estimate.tsv"), sep = "\t", na = "NA")

writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
writeLines(c(
  "seed=20260726",
  paste0("git_commit=", tryCatch(system2("git", c("-C", shQuote(project_dir), "rev-parse", "HEAD"),
                                         stdout = TRUE, stderr = FALSE), error = function(e) "NA"))
), file.path(result_dir, "run_provenance.txt"))

message("QC pilot completed successfully.")
