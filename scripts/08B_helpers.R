#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

stage8b_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

stage8b_paths <- function(project, run_id) {
  x <- list(
    project = project,
    run_id = run_id,
    input_run = file.path(project, "data_processed", "stage_8A_bulk_preprocessing", "20260729_105540"),
    input_results = file.path(project, "results", "08A_bulk_preprocessing", "20260729_105540"),
    result = file.path(project, "results", "08B_bulk_validation", run_id),
    figure = file.path(project, "figures", "08B_bulk_validation", run_id),
    log = file.path(project, "logs", "08B_bulk_validation"),
    source = file.path(project, "environment", "sources", "08B"),
    membership = file.path(project, "results_final", "stage_6A_stage_blind_module_membership.tsv"),
    candidates = file.path(project, "results_final", "stage_6A_exploratory_candidate_modules.tsv")
  )
  for (d in x[c("result", "figure", "log", "source")]) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x
}

stage8b_locked <- function(paths) {
  cand <- read.delim(paths$candidates, check.names = FALSE)
  mem <- read.delim(paths$membership, check.names = FALSE)
  cand <- cand[cand$exploratory_candidate %in% c(TRUE, "TRUE"), , drop = FALSE]
  mem <- mem[mem$module_id %in% cand$module_id, c("module_id", "epithelial_state", "gene")]
  mem$gene <- toupper(trimws(mem$gene))
  mem <- unique(mem[nzchar(mem$gene), ])
  if (length(unique(mem$module_id)) != nrow(cand)) stop("Frozen module membership is incomplete")
  list(candidates = cand, membership = mem)
}

stage8b_module_scores <- function(matrix_file, mapping_file, metadata_file, paths, cohort) {
  expr <- readRDS(matrix_file)
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  meta <- read.delim(metadata_file, check.names = FALSE)
  meta$original_patient_initials <- NA_character_
  if (cohort == "GSE8671") {
    pair_file <- file.path(paths$project, "metadata", "GSE8671_verified_pairs.tsv")
    pairs <- read.delim(pair_file, check.names = FALSE)
    idx <- match(meta$sample_id, pairs$sample_id)
    if (anyNA(idx) || anyDuplicated(pairs$sample_id) || length(unique(pairs$verified_pair_id)) != 32L) {
      stop("GSE8671 explicit patient-number pairing audit failed")
    }
    meta$original_patient_initials <- pairs$original_patient_initials[idx]
    meta$donor_id <- pairs$verified_donor_id[idx]
    meta$paired_group <- pairs$verified_pair_id[idx]
    meta$donor_id_status <- "verified_explicit_GEO_patient_number"
  }
  map <- read.delim(mapping_file, check.names = FALSE)
  if (is.null(rownames(expr))) {
    annotation_file <- file.path(dirname(mapping_file), "probe_annotation.tsv")
    annotation <- read.delim(annotation_file, check.names = FALSE)
    if (!"probe_id" %in% names(annotation) || nrow(annotation) != nrow(expr) ||
        anyDuplicated(annotation$probe_id)) {
      stop(cohort, ": matrix row names are absent and probe annotation cannot restore them")
    }
    rownames(expr) <- annotation$probe_id
  }
  locked <- stage8b_locked(paths)
  map <- map[map$module_id %in% locked$candidates$module_id &
               map$mapped %in% c(TRUE, "TRUE") &
               !is.na(map$probe_id) & nzchar(map$probe_id), ]
  if (!all(meta$sample_id %in% colnames(expr))) stop(cohort, ": metadata samples absent from matrix")
  expr <- expr[, meta$sample_id, drop = FALSE]
  rows <- list()
  coverage <- list()
  for (module in locked$candidates$module_id) {
    mm <- unique(map[map$module_id == module, c("gene", "probe_id")])
    mm$gene <- toupper(mm$gene)
    mm <- mm[mm$probe_id %in% rownames(expr), ]
    genes <- intersect(unique(mm$gene), locked$membership$gene[locked$membership$module_id == module])
    if (!length(genes)) {
      stop(cohort, " / ", module, ": no mapped frozen genes match matrix row names")
    }
    gene_matrix <- t(vapply(genes, function(gene) {
      probes <- unique(mm$probe_id[mm$gene == gene])
      z <- expr[probes, , drop = FALSE]
      if (nrow(z) == 1L) as.numeric(z[1, ]) else apply(z, 2, median, na.rm = TRUE)
    }, numeric(ncol(expr))))
    rownames(gene_matrix) <- genes
    gene_z <- t(scale(t(gene_matrix)))
    gene_z <- gene_z[apply(gene_z, 1, function(z) all(is.finite(z))), , drop = FALSE]
    score_raw <- colMeans(gene_z)
    score <- as.numeric(scale(score_raw))
    rows[[module]] <- data.frame(
      accession = cohort,
      sample_id = colnames(expr),
      module_id = module,
      module_score = score,
      mapped_genes = nrow(gene_z),
      stringsAsFactors = FALSE
    )
    n_locked <- sum(locked$membership$module_id == module)
    coverage[[module]] <- data.frame(
      accession = cohort, module_id = module, mapped_genes = nrow(gene_z),
      locked_genes = n_locked, coverage = nrow(gene_z) / n_locked
    )
  }
  scores <- do.call(rbind, rows)
  scores <- merge(scores, meta, by = c("accession", "sample_id"), all.x = TRUE, sort = FALSE)
  list(scores = scores, coverage = do.call(rbind, coverage))
}

stage8b_robust_coef <- function(fit, term) {
  vc <- sandwich::vcovHC(fit, type = "HC3")
  tab <- lmtest::coeftest(fit, vcov. = vc)
  if (!term %in% rownames(tab)) stop("Model term absent: ", term)
  est <- unname(tab[term, 1])
  se <- unname(tab[term, 2])
  df <- max(stats::df.residual(fit), 1)
  crit <- stats::qt(0.975, df)
  c(effect = est, se = se, ci_low = est - crit * se, ci_high = est + crit * se,
    p_value = unname(tab[term, 4]), n = stats::nobs(fit), df = df)
}

stage8b_fit_trend <- function(x, score_col = "module_score") {
  fit <- stats::lm(stats::as.formula(paste(score_col, "~ stage_score")), data = x)
  stage8b_robust_coef(fit, "stage_score")
}

stage8b_fit_binary <- function(x, early, late, score_col = "module_score") {
  z <- x[x$analysis_group %in% c(early, late), , drop = FALSE]
  z$analysis_group <- factor(z$analysis_group, levels = c(early, late))
  fit <- stats::lm(stats::as.formula(paste(score_col, "~ analysis_group")), data = z)
  stage8b_robust_coef(fit, paste0("analysis_group", late))
}

stage8b_fit_paired <- function(x, early, late, score_col = "module_score") {
  z <- x[x$analysis_group %in% c(early, late), c("donor_id", "analysis_group", score_col)]
  names(z)[3] <- "score"
  wide <- reshape(z, idvar = "donor_id", timevar = "analysis_group", direction = "wide")
  wide <- wide[complete.cases(wide), ]
  wide$difference <- wide[[paste0("score.", late)]] - wide[[paste0("score.", early)]]
  fit <- stats::lm(difference ~ 1, data = wide)
  out <- stage8b_robust_coef(fit, "(Intercept)")
  out["n"] <- nrow(wide)
  out
}

stage8b_effect_row <- function(accession, module_id, endpoint, analysis_set, fit, model, note = NA_character_) {
  data.frame(
    accession = accession, module_id = module_id, endpoint = endpoint,
    analysis_set = analysis_set, model = model, effect = fit["effect"],
    standard_error = fit["se"], ci_low = fit["ci_low"], ci_high = fit["ci_high"],
    p_value = fit["p_value"], n_units = fit["n"], residual_df = fit["df"],
    note = note, stringsAsFactors = FALSE
  )
}

stage8b_add_fdr <- function(x) {
  x$fdr <- NA_real_
  key <- interaction(x$accession, x$endpoint, x$analysis_set, drop = TRUE)
  for (k in levels(key)) {
    idx <- which(key == k)
    x$fdr[idx] <- p.adjust(x$p_value[idx], method = "BH")
  }
  x
}
