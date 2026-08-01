# Shared helpers for Stage 7 independent single-cell replication.

stage7_seed <- 20260728L

stage7_paths <- function(project_dir) {
  list(
    project = normalizePath(project_dir, winslash = "/", mustWork = TRUE),
    result = file.path(project_dir, "results", "07_singlecell_replication"),
    source = file.path(
      project_dir, "results", "07_singlecell_replication", "source_data"
    ),
    figure = file.path(project_dir, "figures", "07_singlecell_replication"),
    object = file.path(project_dir, "objects"),
    processed = file.path(project_dir, "data_processed"),
    log = file.path(project_dir, "logs", "07_singlecell_replication")
  )
}

stage7_init <- function(project_dir) {
  paths <- stage7_paths(project_dir)
  invisible(lapply(
    paths[c("result", "source", "figure", "object", "processed", "log")],
    dir.create, recursive = TRUE, showWarnings = FALSE
  ))
  paths
}

read_stage7_parameters <- function(project_dir) {
  x <- utils::read.delim(
    file.path(project_dir, "config", "07_replication_parameters.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  function(section, parameter, numeric = FALSE) {
    hit <- x$value[x$section == section & x$parameter == parameter]
    if (length(hit) != 1L) {
      stop("Stage 7 parameter is absent or duplicated: ", section, "/", parameter)
    }
    if (numeric) as.numeric(hit) else hit
  }
}

write_stage7_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  target <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else path
  on.exit(if (inherits(target, "connection")) close(target), add = TRUE)
  utils::write.table(
    x, target, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

read_locked_stage7_modules <- function(project_dir) {
  candidates <- utils::read.delim(
    file.path(
      project_dir, "results", "06A_amendment",
      "exploratory_candidate_modules.tsv"
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  membership <- utils::read.delim(
    file.path(
      project_dir, "results", "06A_amendment",
      "stage_blind_module_membership.tsv"
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  candidates <- candidates[
    candidates$passes_LODO & candidates$exploratory_candidate, ,
    drop = FALSE
  ]
  membership <- membership[membership$module_id %in% candidates$module_id, ]
  if (nrow(candidates) != 6L || anyDuplicated(candidates$module_id)) {
    stop("Expected exactly six unique locked exploratory modules")
  }
  if (!setequal(unique(membership$module_id), candidates$module_id)) {
    stop("Locked module membership is incomplete")
  }
  list(candidates = candidates, membership = membership)
}

robust_qc_bounds <- function(x, lower_k = 3, upper_k = 5) {
  med <- stats::median(x, na.rm = TRUE)
  spread <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) {
    spread <- stats::IQR(x, na.rm = TRUE) / 1.349
  }
  if (!is.finite(spread) || spread == 0) spread <- 0
  c(lower = med - lower_k * spread, upper = med + upper_k * spread)
}

percent_feature_set <- function(counts, pattern) {
  idx <- grepl(pattern, rownames(counts), ignore.case = FALSE)
  totals <- Matrix::colSums(counts)
  out <- rep(0, ncol(counts))
  if (any(idx)) {
    out <- 100 * Matrix::colSums(counts[idx, , drop = FALSE]) / pmax(totals, 1)
  }
  as.numeric(out)
}

major_marker_signatures <- function(candidate_genes = character()) {
  signatures <- list(
    Epithelial = c(
      "EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "ELF3", "CEACAM5",
      "CEACAM6", "GUCA2A", "MUC2"
    ),
    T_NK = c(
      "PTPRC", "CD3D", "CD3E", "TRAC", "CD247", "IL7R", "NKG7", "GNLY",
      "KLRD1", "CCL5"
    ),
    B_cell = c("PTPRC", "CD79A", "MS4A1", "CD37", "CD74", "HLA-DRA", "CD22"),
    Plasma_cell = c(
      "PTPRC", "JCHAIN", "MZB1", "SDC1", "IGHG1", "IGHA1", "IGHA2"
    ),
    Myeloid = c(
      "PTPRC", "LYZ", "FCER1G", "TYROBP", "LST1", "AIF1", "CTSS",
      "FCGR3A", "CD68", "CD14"
    ),
    Fibroblast = c(
      "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "PDGFRA"
    ),
    Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "ENG", "ESAM", "CLDN5"),
    Mast_cell = c("KIT", "MS4A2", "TPSAB1", "TPSB2", "CPA3", "HPGDS"),
    Pericyte_SMC = c("RGS5", "CSPG4", "MCAM", "PDGFRB", "ACTA2", "MYH11"),
    Enteric_glia = c("PLP1", "S100B", "SOX10", "SLC1A3", "S100A1")
  )
  lapply(signatures, function(x) setdiff(unique(x), candidate_genes))
}

score_signatures <- function(log_expression, signatures) {
  present <- lapply(signatures, intersect, y = rownames(log_expression))
  if (any(lengths(present) < 2L)) {
    warning(
      "Some annotation signatures have fewer than two represented genes: ",
      paste(names(present)[lengths(present) < 2L], collapse = ";")
    )
  }
  out <- vapply(
    present,
    function(genes) {
      if (!length(genes)) return(rep(NA_real_, ncol(log_expression)))
      Matrix::colMeans(log_expression[genes, , drop = FALSE])
    },
    numeric(ncol(log_expression))
  )
  rownames(out) <- colnames(log_expression)
  list(scores = out, genes = present)
}

assign_cluster_annotations <- function(cell_scores, clusters, margin = 0.05) {
  cluster_levels <- sort(unique(as.character(clusters)))
  cluster_scores <- do.call(rbind, lapply(cluster_levels, function(cl) {
    idx <- as.character(clusters) == cl
    colMeans(cell_scores[idx, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(cluster_scores) <- cluster_levels
  labels <- apply(cluster_scores, 1, function(x) {
    ord <- order(x, decreasing = TRUE, na.last = TRUE)
    if (!length(ord) || !is.finite(x[ord[1L]])) return("Uncertain")
    if (length(ord) > 1L && is.finite(x[ord[2L]]) &&
        x[ord[1L]] - x[ord[2L]] < margin) return("Uncertain")
    colnames(cluster_scores)[ord[1L]]
  })
  list(
    labels = setNames(labels, cluster_levels),
    cluster_scores = cluster_scores
  )
}

aggregate_counts <- function(counts, groups) {
  if (length(groups) != ncol(counts)) stop("Group vector length does not match cells")
  group_levels <- unique(as.character(groups))
  indicator <- Matrix::sparseMatrix(
    i = seq_along(groups),
    j = match(as.character(groups), group_levels),
    x = 1,
    dims = c(length(groups), length(group_levels)),
    dimnames = list(colnames(counts), group_levels)
  )
  out <- counts %*% indicator
  colnames(out) <- group_levels
  out
}

eligible_pseudobulk <- function(counts, metadata, min_cells, min_umi) {
  ok <- metadata$n_cells >= min_cells &
    metadata$total_umi >= min_umi &
    colnames(counts) %in% metadata$pseudobulk_id
  ids <- metadata$pseudobulk_id[ok]
  list(
    counts = counts[, ids, drop = FALSE],
    metadata = metadata[match(ids, metadata$pseudobulk_id), , drop = FALSE]
  )
}

module_scores_from_counts <- function(
  counts, membership, min_fraction = 0.60, min_genes = 8L
) {
  if (ncol(counts) < 2L) stop("At least two pseudobulks are required")
  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  gene_z <- t(scale(t(logcpm)))
  gene_z[!is.finite(gene_z)] <- 0

  modules <- split(membership$gene, membership$module_id)
  score_list <- lapply(names(modules), function(module_id) {
    locked <- unique(modules[[module_id]])
    represented <- intersect(locked, rownames(gene_z))
    coverage <- length(represented) / length(locked)
    evaluable <- length(represented) >= min_genes && coverage >= min_fraction
    score <- if (evaluable) {
      colMeans(gene_z[represented, , drop = FALSE])
    } else {
      rep(NA_real_, ncol(gene_z))
    }
    data.frame(
      pseudobulk_id = colnames(gene_z),
      module_id = module_id,
      module_score = as.numeric(score),
      locked_genes = length(locked),
      represented_genes = length(represented),
      gene_coverage = coverage,
      evaluable = evaluable,
      stringsAsFactors = FALSE
    )
  })
  list(
    scores = do.call(rbind, score_list),
    logcpm = logcpm,
    normalization = data.frame(
      pseudobulk_id = colnames(y),
      library_size = y$samples$lib.size,
      norm_factor = y$samples$norm.factors,
      stringsAsFactors = FALSE
    )
  )
}

paired_effect <- function(
  scores, metadata, module_id, condition_early, condition_late,
  donor_subset = NULL
) {
  x <- merge(scores[scores$module_id == module_id, ], metadata, by = "pseudobulk_id")
  x <- x[x$condition %in% c(condition_early, condition_late), ]
  if (!is.null(donor_subset)) x <- x[x$donor_id %in% donor_subset, ]
  wide <- reshape(
    x[, c("donor_id", "condition", "module_score")],
    idvar = "donor_id", timevar = "condition", direction = "wide"
  )
  early_col <- paste0("module_score.", condition_early)
  late_col <- paste0("module_score.", condition_late)
  wide <- wide[stats::complete.cases(wide[, c(early_col, late_col)]), ]
  diff <- wide[[late_col]] - wide[[early_col]]
  n <- length(diff)
  effect <- if (n) mean(diff) else NA_real_
  sd_diff <- if (n > 1L) stats::sd(diff) else NA_real_
  se <- if (n > 1L) sd_diff / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(0.975, df = n - 1L) else NA_real_
  p <- if (n > 1L && is.finite(sd_diff) && sd_diff > 0) {
    stats::t.test(diff, mu = 0)$p.value
  } else if (n > 1L && all(diff == 0)) {
    1
  } else {
    NA_real_
  }
  standardized <- if (is.finite(sd_diff) && sd_diff > 0) effect / sd_diff else NA_real_
  lodo <- if (n >= 3L) vapply(seq_len(n), function(i) mean(diff[-i]), numeric(1)) else numeric()
  sign_stability <- if (length(lodo) && is.finite(effect) && effect != 0) {
    mean(sign(lodo) == sign(effect))
  } else {
    NA_real_
  }
  max_contribution <- if (n && sum(abs(diff)) > 0) max(abs(diff)) / sum(abs(diff)) else NA_real_
  details <- if (n) {
    data.frame(
      module_id = rep(module_id, n),
      donor_id = wide$donor_id,
      early_score = wide[[early_col]],
      late_score = wide[[late_col]],
      difference = diff,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      module_id = character(), donor_id = character(),
      early_score = numeric(), late_score = numeric(), difference = numeric(),
      stringsAsFactors = FALSE
    )
  }
  summary <- data.frame(
    module_id = module_id,
    early_condition = condition_early,
    late_condition = condition_late,
    n_donors = n,
    effect = effect,
    ci_low = effect - critical * se,
    ci_high = effect + critical * se,
    standardized_effect = standardized,
    p_value = p,
    lodo_sign_stability = sign_stability,
    max_absolute_donor_contribution = max_contribution,
    donor_ids = paste(wide$donor_id, collapse = ";"),
    stringsAsFactors = FALSE
  )
  list(summary = summary, details = details)
}

specificity_effect <- function(scores, metadata, module_id) {
  x <- merge(scores[scores$module_id == module_id, ], metadata, by = "pseudobulk_id")
  x <- x[is.finite(x$module_score), ]
  donors <- unique(x$donor_id)
  details <- do.call(rbind, lapply(donors, function(donor) {
    d <- x[x$donor_id == donor, ]
    conditions <- unique(d$condition)
    condition_diffs <- vapply(conditions, function(cond) {
      z <- d[d$condition == cond, ]
      epithelial <- z$module_score[z$major_cell_type == "Epithelial"]
      other <- z$module_score[z$major_cell_type != "Epithelial"]
      if (length(epithelial) != 1L || !length(other)) return(NA_real_)
      epithelial - stats::median(other)
    }, numeric(1))
    condition_diffs <- condition_diffs[is.finite(condition_diffs)]
    if (!length(condition_diffs)) return(NULL)
    data.frame(
      module_id = module_id, donor_id = donor,
      difference = mean(condition_diffs),
      n_conditions = length(condition_diffs),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(details) || !nrow(details)) {
    details <- data.frame(
      module_id = character(), donor_id = character(),
      difference = numeric(), n_conditions = integer()
    )
  }
  diff <- details$difference
  n <- length(diff)
  effect <- if (n) mean(diff) else NA_real_
  sd_diff <- if (n > 1L) stats::sd(diff) else NA_real_
  se <- if (n > 1L) sd_diff / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt(0.975, n - 1L) else NA_real_
  p <- if (n > 1L && is.finite(sd_diff) && sd_diff > 0) {
    stats::t.test(diff)$p.value
  } else {
    NA_real_
  }
  lodo <- if (n >= 3L) vapply(seq_len(n), function(i) mean(diff[-i]), numeric(1)) else numeric()
  summary <- data.frame(
    module_id = module_id,
    early_condition = "median_non_epithelial",
    late_condition = "epithelial",
    n_donors = n,
    effect = effect,
    ci_low = effect - critical * se,
    ci_high = effect + critical * se,
    standardized_effect = if (is.finite(sd_diff) && sd_diff > 0) effect / sd_diff else NA_real_,
    p_value = p,
    lodo_sign_stability = if (length(lodo) && is.finite(effect) && effect != 0) {
      mean(sign(lodo) == sign(effect))
    } else NA_real_,
    max_absolute_donor_contribution = if (n && sum(abs(diff)) > 0) {
      max(abs(diff)) / sum(abs(diff))
    } else NA_real_,
    donor_ids = paste(details$donor_id, collapse = ";"),
    stringsAsFactors = FALSE
  )
  list(summary = summary, details = details)
}

add_fdr_by_contrast <- function(x) {
  x$fdr <- ave(
    x$p_value, interaction(x$cohort, x$contrast, drop = TRUE),
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
  x
}
