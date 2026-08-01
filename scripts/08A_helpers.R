#!/usr/bin/env Rscript

# Stage 8A shared helpers. This script performs preprocessing/QC only; it
# deliberately contains no ordered-trend, differential-expression, or model-fitting step.

options(stringsAsFactors = FALSE, warn = 1)

stage8_paths <- function(project_dir, run_id) {
  root <- file.path(project_dir, "data_processed", "stage_8A_bulk_preprocessing", run_id)
  out <- list(
    root = root,
    work = file.path(root, "work"),
    result = file.path(project_dir, "results", "08A_bulk_preprocessing", run_id),
    figure = file.path(project_dir, "figures", "08A_bulk_preprocessing", run_id),
    report = file.path(project_dir, "reports", "stage_8A_bulk_preprocessing.md"),
    manifest = file.path(project_dir, "metadata", "dataset_manifest.tsv"),
    file_manifest = file.path(project_dir, "metadata", "file_manifest.tsv"),
    membership = file.path(project_dir, "results", "07_singlecell_replication", "preflight", "locked_module_membership.tsv")
  )
  for (x in c(out$root, out$work, out$result, out$figure)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  out
}

stage8_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

stage8_manifest <- function(paths, accession) {
  x <- utils::read.delim(paths$manifest, check.names = FALSE, stringsAsFactors = FALSE)
  x <- x[x$accession == accession & x$inclusion == "include", , drop = FALSE]
  if (!nrow(x)) stop("No included manifest rows for ", accession)
  if (anyDuplicated(x$sample_id)) stop("Duplicated sample_id in manifest for ", accession)
  x
}

stage8_locked_membership <- function(paths) {
  if (!file.exists(paths$membership)) stop("Locked Stage 7 module membership is absent")
  x <- utils::read.delim(paths$membership, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("module_id", "epithelial_state", "gene")
  if (!all(required %in% names(x))) stop("Locked membership columns are incomplete")
  x$gene <- toupper(trimws(x$gene))
  x <- x[nzchar(x$gene) & !is.na(x$gene), required, drop = FALSE]
  x[!duplicated(x), , drop = FALSE]
}

stage8_extract <- function(raw_tar, out_dir) {
  if (!file.exists(raw_tar)) stop("Missing archive: ", raw_tar)
  if (length(list.files(out_dir, all.files = FALSE, recursive = TRUE))) return(invisible(out_dir))
  utils::untar(raw_tar, exdir = out_dir)
  invisible(out_dir)
}

stage8_map_symbols <- function(probe, symbol, membership, accession, evidence) {
  probe <- as.character(probe)
  symbol <- as.character(symbol)
  spl <- strsplit(toupper(symbol), "\\s*(//|;|,|\\|)\\s*")
  expanded <- do.call(rbind, Filter(Negate(is.null), lapply(seq_along(spl), function(i) {
    genes <- trimws(spl[[i]])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) return(NULL)
    data.frame(probe_id = rep(probe[[i]], length(genes)), gene = genes, stringsAsFactors = FALSE)
  })))
  if (is.null(expanded)) expanded <- data.frame(probe_id = character(), gene = character(), stringsAsFactors = FALSE)
  out <- merge(membership, expanded, by = "gene", all.x = TRUE, sort = FALSE)
  out$accession <- accession
  out$mapping_evidence <- evidence
  out$mapped <- !is.na(out$probe_id)
  out <- out[, c("accession", "module_id", "epithelial_state", "gene", "mapped", "probe_id", "mapping_evidence")]
  out
}

stage8_gene_median <- function(mat, symbols) {
  symbols <- toupper(trimws(as.character(symbols)))
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  split_rows <- split(seq_along(symbols), symbols)
  ans <- vapply(split_rows, function(idx) {
    if (length(idx) == 1L) mat[idx, ] else apply(mat[idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  }, numeric(ncol(mat)))
  if (is.null(dim(ans))) ans <- matrix(ans, nrow = 1L)
  ans <- t(ans)
  rownames(ans) <- names(split_rows)
  colnames(ans) <- colnames(mat)
  ans
}

stage8_qc <- function(mat, sample_meta, accession, paths) {
  qc_mat <- mat
  qc_mat[!is.finite(qc_mat)] <- NA_real_
  complete_features <- rowSums(is.na(qc_mat)) == 0L
  if (sum(complete_features) < 10L) stop(accession, " has fewer than 10 complete features for PCA")
  feature_audit <- data.frame(
    accession = accession,
    total_features = nrow(mat),
    complete_features_used_for_pca = sum(complete_features),
    features_excluded_from_pca_nonfinite = sum(!complete_features),
    full_matrix_preserved = TRUE,
    stringsAsFactors = FALSE
  )
  stage8_write_tsv(feature_audit, file.path(paths$result, accession, "qc_feature_filter.tsv"))
  med <- apply(qc_mat, 2, stats::median, na.rm = TRUE)
  iqr <- apply(qc_mat, 2, stats::IQR, na.rm = TRUE)
  cor_mat <- stats::cor(qc_mat, use = "pairwise.complete.obs")
  connectivity <- colMeans(cor_mat, na.rm = TRUE)
  pca <- stats::prcomp(t(qc_mat[complete_features, , drop = FALSE]), center = TRUE, scale. = FALSE)
  z_mad <- function(x) {
    d <- stats::mad(x, constant = 1, na.rm = TRUE)
    if (is.na(d) || d == 0) rep(0, length(x)) else (x - stats::median(x, na.rm = TRUE)) / d
  }
  qc <- data.frame(
    sample_id = colnames(mat), median_log2_expression = med, iqr_log2_expression = iqr,
    mean_sample_correlation = connectivity, median_z_mad = z_mad(med),
    iqr_z_mad = z_mad(iqr), connectivity_z_mad = z_mad(connectivity),
    stringsAsFactors = FALSE
  )
  qc$qc_flag_review <- abs(qc$median_z_mad) > 3 | abs(qc$iqr_z_mad) > 3 | abs(qc$connectivity_z_mad) > 3
  qc <- merge(sample_meta, qc, by = "sample_id", all.y = TRUE, sort = FALSE)
  stage8_write_tsv(qc, file.path(paths$result, accession, "qc_sample_metrics.tsv"))
  pca_df <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
  pca_df <- merge(sample_meta, pca_df, by = "sample_id", all.y = TRUE, sort = FALSE)
  stage8_write_tsv(pca_df, file.path(paths$result, accession, "qc_pca_coordinates.tsv"))
  stage8_plot_qc(qc, pca_df, accession, paths)
  list(qc = qc, pca = pca_df)
}

stage8_plot_qc <- function(qc, pca, accession, paths) {
  png <- file.path(paths$figure, paste0(accession, "_qc.png"))
  pdf <- file.path(paths$figure, paste0(accession, "_qc.pdf"))
  draw <- function() {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    on.exit(par(op), add = TRUE)
    boxplot(median_log2_expression ~ condition, data = qc, las = 2, ylab = "Median log2 expression", main = paste(accession, "QC"))
    lev <- unique(pca$condition)
    col <- grDevices::hcl.colors(length(lev), "Dark 3")
    plot(pca$PC1, pca$PC2, col = col[match(pca$condition, lev)], pch = 19,
         xlab = "PC1", ylab = "PC2", main = paste(accession, "PCA (QC only)"))
    legend("topright", legend = lev, col = col, pch = 19, cex = 0.75, bty = "n")
  }
  grDevices::png(png, width = 2400, height = 1200, res = 300); draw(); grDevices::dev.off()
  grDevices::pdf(pdf, width = 10, height = 4); draw(); grDevices::dev.off()
}

stage8_write_mapping_audit <- function(mapping, method, d) {
  gene_status <- aggregate(
    mapped ~ accession + module_id + epithelial_state + gene,
    data = mapping,
    FUN = any
  )
  mapped <- aggregate(
    mapped ~ accession + module_id + epithelial_state,
    data = gene_status,
    FUN = sum
  )
  total <- aggregate(
    gene ~ accession + module_id + epithelial_state,
    data = gene_status,
    FUN = length
  )
  names(total)[4] <- "locked_genes"
  audit <- merge(total, mapped, by = c("accession", "module_id", "epithelial_state"), all.x = TRUE)
  audit$unmapped_genes <- audit$locked_genes - audit$mapped
  audit$normalization_method <- method
  stage8_write_tsv(audit, file.path(d, "locked_module_mapping_summary.tsv"))
  stage8_write_tsv(gene_status[!gene_status$mapped, ], file.path(d, "unmapped_locked_genes.tsv"))
  invisible(audit)
}

stage8_save_outputs <- function(mat, probe_annotation, mapping, accession, method, paths) {
  d <- file.path(paths$result, accession)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  saveRDS(mat, file.path(paths$root, paste0(accession, "_normalized_probe_matrix.rds")), compress = "xz")
  stage8_write_tsv(probe_annotation, file.path(d, "probe_annotation.tsv"))
  stage8_write_tsv(mapping, file.path(d, "locked_module_gene_mapping.tsv"))
  stage8_write_mapping_audit(mapping, method, d)
}

stage8_preprocess_gse41657 <- function(project_dir, run_id) {
  if (!requireNamespace("limma", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) stop("GSE41657 requires limma, AnnotationDbi and org.Hs.eg.db")
  paths <- stage8_paths(project_dir, run_id); accession <- "GSE41657"
  raw_dir <- file.path(project_dir, "data_raw", accession)
  ext <- file.path(paths$work, accession, "raw_txt")
  stage8_extract(file.path(raw_dir, "GSE41657_RAW.tar"), ext)
  files <- sort(list.files(ext, pattern = "\\.txt\\.gz$", full.names = TRUE))
  meta <- stage8_manifest(paths, accession)
  sample_id <- sub("_.*$", "", basename(files))
  if (length(files) != nrow(meta) || !setequal(sample_id, meta$sample_id)) stop("GSE41657 raw TXT files do not reconcile to included manifest samples")
  raw <- limma::read.maimages(files, source = "agilent", green.only = TRUE)
  corrected <- limma::backgroundCorrect(raw, method = "normexp", offset = 50)
  log_expression <- log2(pmax(corrected$E, 1))
  mat <- limma::normalizeBetweenArrays(log_expression, method = "quantile")
  colnames(mat) <- sample_id
  probes <- make.unique(as.character(raw$genes$ProbeName))
  rownames(mat) <- probes
  acc <- sub("\\..*$", "", as.character(raw$genes$SystematicName))
  symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = acc, keytype = "REFSEQ", column = "SYMBOL", multiVals = "first")
  symbols <- unname(symbols[acc])
  probe_annotation <- data.frame(probe_id = probes, refseq_accession = acc, gene_symbol = symbols, stringsAsFactors = FALSE)
  membership <- stage8_locked_membership(paths)
  mapping <- stage8_map_symbols(probes, symbols, membership, accession, "Agilent SystematicName RefSeq -> org.Hs.eg.db SYMBOL")
  stage8_save_outputs(mat, probe_annotation, mapping, accession, "Agilent one-color normexp background correction (offset 50), log2 transform, and quantile normalization", paths)
  stage8_write_tsv(stage8_manifest(paths, accession), file.path(paths$result, accession, "sample_metadata.tsv"))
  stage8_qc(mat, stage8_manifest(paths, accession), accession, paths)
}

stage8_preprocess_gse100179 <- function(project_dir, run_id) {
  paths <- stage8_paths(project_dir, run_id); accession <- "GSE100179"
  f <- file.path(project_dir, "data_raw", accession, "GSE100179_2017-07-21-annotated-RMA-SKETCH.RMA-GENE-FULL-Group1.txt.gz")
  if (!file.exists(f)) stop("GSE100179 official RMA-Sketch matrix is absent")
  x <- utils::read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
  signal <- grep("-Signal$", names(x), value = TRUE)
  if (length(signal) != 60L || !"Probe.Set.ID" %in% make.names(names(x))) stop("GSE100179 official matrix has unexpected columns")
  sample_label <- sub("_\\(HTA.*$", "", signal)
  mat <- as.matrix(x[, signal, drop = FALSE]); storage.mode(mat) <- "numeric"
  raw_dir <- file.path(project_dir, "data_raw", accession)
  ext <- file.path(paths$work, accession, "raw_cel_names")
  stage8_extract(file.path(raw_dir, "GSE100179_RAW.tar"), ext)
  cel <- basename(list.files(ext, pattern = "\\.CEL\\.gz$", full.names = TRUE))
  raw_label <- sub("^GSM[0-9]+_", "", cel); raw_label <- sub("_HTA.*$", "", raw_label)
  raw_id <- sub("_.*$", "", cel)
  idx <- match(tolower(sample_label), tolower(raw_label))
  if (anyNA(idx) || anyDuplicated(idx)) stop("GSE100179 processed-matrix labels cannot be uniquely reconciled to raw CEL labels")
  colnames(mat) <- raw_id[idx]
  meta <- stage8_manifest(paths, accession)
  if (!setequal(colnames(mat), meta$sample_id)) stop("GSE100179 RMA matrix does not reconcile to included manifest samples")
  mat <- mat[, meta$sample_id, drop = FALSE]
  symbols <- toupper(trimws(as.character(x[["Gene Symbol"]])))
  probe <- as.character(x[["Probe Set ID"]])
  probe_annotation <- data.frame(probe_id = probe, gene_symbol = symbols, stringsAsFactors = FALSE)
  membership <- stage8_locked_membership(paths)
  mapping <- stage8_map_symbols(probe, symbols, membership, accession, "Official GEO RMA-Sketch gene-symbol annotation")
  stage8_save_outputs(mat, probe_annotation, mapping, accession, "Official GEO HTA 2.0 RMA-Sketch gene-level matrix; sample labels reconciled to raw CEL archive", paths)
  stage8_write_tsv(meta, file.path(paths$result, accession, "sample_metadata.tsv"))
  stage8_qc(mat, meta, accession, paths)
}

stage8_preprocess_gse8671 <- function(project_dir, run_id) {
  if (!requireNamespace("hgu133plus2.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) stop("GSE8671 requires AnnotationDbi and hgu133plus2.db")
  paths <- stage8_paths(project_dir, run_id); accession <- "GSE8671"
  raw_dir <- file.path(project_dir, "data_raw", accession)
  ext <- file.path(paths$work, accession, "raw_cel")
  stage8_extract(file.path(raw_dir, "GSE8671_RAW.tar"), ext)
  cel <- sort(list.files(ext, pattern = "\\.CEL\\.gz$", full.names = TRUE))
  ids <- sub("\\.CEL\\.gz$", "", basename(cel))
  meta <- stage8_manifest(paths, accession)
  if (length(cel) != nrow(meta) || !setequal(ids, meta$sample_id)) stop("GSE8671 CEL files do not reconcile to included manifest samples")
  use_raw_rma <- requireNamespace("affy", quietly = TRUE) && requireNamespace("hgu133plus2cdf", quietly = TRUE)
  if (use_raw_rma) {
    files <- utils::read.delim(paths$file_manifest, check.names = FALSE, stringsAsFactors = FALSE)
    raw_record <- files[
      files$accession == accession &
        files$file_name == "GSE8671_RAW.tar" &
        files$status == "verified",
      , drop = FALSE
    ]
    if (nrow(raw_record) != 1L ||
        !nzchar(raw_record$official_url[[1L]]) ||
        nchar(raw_record$sha256[[1L]]) != 64L) {
      stop("GSE8671 verified raw-CEL provenance is incomplete in file_manifest.tsv")
    }
    raw <- affy::ReadAffy(filenames = cel)
    eset <- affy::rma(raw)
    mat <- Biobase::exprs(eset); colnames(mat) <- ids
    mat <- mat[, meta$sample_id, drop = FALSE]
    method <- "Affymetrix U133 Plus 2.0 RMA from raw CEL files"
    provenance <- data.frame(
      accession = accession, selected_input = "raw_CEL", processing_status = "raw_RMA_completed",
      reason = "Verified official raw CEL archive and hgu133plus2cdf were available",
      source_url = raw_record$official_url[[1L]],
      source_sha256 = raw_record$sha256[[1L]],
      stringsAsFactors = FALSE
    )
  } else {
    geo_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8671/matrix/GSE8671_series_matrix.txt.gz"
    geo_file <- file.path(paths$work, accession, "GSE8671_series_matrix.txt.gz")
    if (!file.exists(geo_file)) utils::download.file(geo_url, geo_file, mode = "wb", method = "libcurl", quiet = TRUE)
    lines <- readLines(gzfile(geo_file), warn = FALSE)
    begin <- which(lines == "!series_matrix_table_begin")
    end <- which(lines == "!series_matrix_table_end")
    if (length(begin) != 1L || length(end) != 1L || end <= begin) stop("GSE8671 GEO series matrix table boundaries are invalid")
    tab <- utils::read.delim(gzfile(geo_file), skip = begin, nrows = end - begin - 1L, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(tab) != 65L) stop("GSE8671 GEO matrix does not contain 64 sample columns")
    mat <- as.matrix(tab[, -1L, drop = FALSE]); storage.mode(mat) <- "numeric"
    rownames(mat) <- as.character(tab[[1L]])
    if (!setequal(colnames(mat), meta$sample_id)) stop("GSE8671 GEO matrix samples do not reconcile to the manifest")
    mat <- mat[, meta$sample_id, drop = FALSE]
    method <- "Official GEO GSE8671 series processed matrix; no secondary normalization"
    provenance <- data.frame(
      accession = accession, selected_input = "official_GEO_series_matrix", processing_status = "raw_RMA_not_run",
      reason = "hgu133plus2cdf absent and Bioconductor repository unreachable; fallback is explicit and auditable",
      source_url = geo_url,
      source_sha256 = sub("\\s+.*$", "", system2("sha256sum", shQuote(geo_file), stdout = TRUE)), stringsAsFactors = FALSE
    )
  }
  probe <- rownames(mat)
  ann <- AnnotationDbi::select(hgu133plus2.db::hgu133plus2.db, keys = probe, keytype = "PROBEID", columns = "SYMBOL")
  ann <- ann[!duplicated(ann$PROBEID), , drop = FALSE]
  symbols <- ann$SYMBOL[match(probe, ann$PROBEID)]
  probe_annotation <- data.frame(probe_id = probe, gene_symbol = symbols, stringsAsFactors = FALSE)
  membership <- stage8_locked_membership(paths)
  mapping <- stage8_map_symbols(probe, symbols, membership, accession, "hgu133plus2.db PROBEID -> SYMBOL")
  stage8_save_outputs(mat, probe_annotation, mapping, accession, method, paths)
  stage8_write_tsv(meta, file.path(paths$result, accession, "sample_metadata.tsv"))
  stage8_write_tsv(provenance, file.path(paths$result, accession, "input_processing_provenance.tsv"))
  stage8_qc(mat, meta, accession, paths)
}
