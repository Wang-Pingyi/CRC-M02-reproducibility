#!/usr/bin/env Rscript

set.seed(20260728)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else getwd()
source(file.path(project_dir, "scripts", "07_helpers.R"))
paths <- stage7_init(project_dir)

summary <- utils::read.delim(
  file.path(paths$result, "replication_summary.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
effects <- utils::read.delim(
  file.path(paths$result, "replication_effects_all.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
gate <- utils::read.delim(
  file.path(paths$result, "replication_gate.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
qc161 <- utils::read.delim(
  file.path(paths$result, "GSE161277_qc_summary.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
qc132 <- utils::read.delim(
  file.path(paths$result, "GSE132465_qc_summary.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
comp161 <- utils::read.delim(
  file.path(paths$result, "GSE161277_cell_composition.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
comp132 <- utils::read.delim(
  file.path(paths$result, "GSE132465_cell_composition.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
eligibility132 <- utils::read.delim(
  file.path(
    paths$result, "preflight", "GSE132465_pair_eligibility.tsv"
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)
evaluable132 <- unique(
  eligibility132$donor_id[eligibility132$paired_analysis_eligible]
)
excluded132 <- unique(
  eligibility132$donor_id[!eligibility132$paired_analysis_eligible]
)

overall <- gate$passed[gate$criterion == "overall_replication_gate"]
class_counts <- table(factor(
  summary$replication_class,
  levels = c(
    "replicated", "directionally consistent but underpowered",
    "not replicated", "contradictory"
  )
))
fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

module_lines <- vapply(seq_len(nrow(summary)), function(i) {
  x <- summary[i, ]
  paste0(
    "| ", x$module_id, " | ", x$replication_class, " | ",
    fmt(x$GSE161277_normal_to_adenoma_effect), " (",
    fmt(x$GSE161277_normal_to_adenoma_ci_low), ", ",
    fmt(x$GSE161277_normal_to_adenoma_ci_high), "); FDR ",
    fmt(x$GSE161277_normal_to_adenoma_fdr), " | ",
    fmt(x$GSE132465_normal_to_cancer_effect), " (",
    fmt(x$GSE132465_normal_to_cancer_ci_low), ", ",
    fmt(x$GSE132465_normal_to_cancer_ci_high), "); FDR ",
    fmt(x$GSE132465_normal_to_cancer_fdr), " | ",
    fmt(x$epithelial_specificity_effect), "; FDR ",
    fmt(x$epithelial_specificity_fdr), " |"
  )
}, character(1))

report <- c(
  "# Stage 7 independent single-cell replication",
  "",
  "## Stage decision",
  "",
  paste0(
    "**Replication gate: ", if (isTRUE(overall)) "passed" else "not passed",
    ".**"
  ),
  "",
  paste0(
    "Classification counts: replicated ", class_counts[["replicated"]],
    "; directionally consistent but underpowered ",
    class_counts[["directionally consistent but underpowered"]],
    "; not replicated ", class_counts[["not replicated"]],
    "; contradictory ", class_counts[["contradictory"]], "."
  ),
  "",
  "The primary Stage 6A candidate table remained empty. This stage tested only",
  "the six previously locked exploratory stage-blind modules. Therefore all",
  "evidence remains secondary and module-level even if the gate passed.",
  "",
  "## Cohort processing and independence",
  "",
  paste0(
    "- GSE161277: ", sum(qc161$cells_input), " input cells across ",
    nrow(qc161), " captures; ", sum(qc161$cells_retained),
    " tissue cells retained after per-capture adaptive QC and scDblFinder."
  ),
  paste0(
    "- GSE161277 epithelial cells retained for annotation/pseudobulk inventory: ",
    sum(comp161$n_cells[comp161$major_cell_type == "Epithelial"]), "."
  ),
  "- GSE161277 blood was audited but excluded; para-cancer tissue stayed",
  "  separate from normal mucosa and was excluded from primary stage contrasts.",
  "- Patient3 adenoma 1 and adenoma 2 remained traceable as separate tissues but",
  "  were pooled at donor-stage raw-count aggregation.",
  paste0(
    "- GSE132465: ", sum(qc132$cells_input), " deposited processed cells; ",
    sum(qc132$cells_retained), " retained after conservative residual QC."
  ),
  paste0(
    "- GSE132465 epithelial cells retained: ",
    sum(comp132$n_cells[comp132$major_cell_type == "Epithelial"]), "."
  ),
  "- GSE132465 official major-cell annotations were retained and independently",
  "  audited with canonical markers. A new doublet caller was not applied because",
  "  the official processed matrix had already removed high-gene outliers.",
  "- The two validation cohorts were normalized, annotated and modeled",
  "  independently. No expression matrix was merged with GSE201348 or across",
  "  validation cohorts.",
  "",
  "## Statistical analysis",
  "",
  "- Biological replicate: donor/patient.",
  "- Raw counts were summed to donor × condition × major-cell-type pseudobulks.",
  "- TMM logCPM values were used for frozen-module scoring.",
  "- GSE161277 primary contrasts used the three fully matched donors",
  "  Patient1–Patient3.",
  paste0(
    "- GSE132465 contained 10 source-matched donors; ", length(evaluable132),
    " passed both-condition epithelial pseudobulk gates and entered the paired",
    " analysis (", paste(evaluable132, collapse = ", "), ")."
  ),
  paste0(
    "- Excluded source-matched donors: ", paste(excluded132, collapse = ", "),
    " because at least one condition had fewer than 20 epithelial cells."
  ),
  "- Effects are paired donor differences with 95% t confidence intervals,",
  "  exact P values and BH FDR across the six locked modules per contrast.",
  "- Leave-one-donor-out direction stability was evaluated; no cell-level",
  "  significance was used.",
  "",
  "## Locked-module replication results",
  "",
  "| Module | Classification | GSE161277 adenoma−normal effect (95% CI); FDR | GSE132465 cancer−normal effect (95% CI); FDR | Epithelial specificity effect; FDR |",
  "|---|---|---:|---:|---:|",
  module_lines,
  "",
  "Full effects, confidence intervals, FDR values, LODO statistics and donor",
  "differences are in `results/07_singlecell_replication/` and its",
  "`source_data/` directory. Both significant and nonsignificant validation",
  "results are retained.",
  "",
  "## Replication gate audit",
  "",
  paste0(
    "- ", gate$criterion, ": ", ifelse(gate$passed, "PASS", "FAIL"),
    " — ", gate$note
  ),
  "",
  "A gate pass permits only continuation to the separately authorized tissue",
  "validation stage; it does not convert the exploratory modules into primary",
  "discoveries. A gate failure requires annotation/model review and forbids",
  "candidate reselection or machine-learning escalation.",
  "",
  "## Reproducibility and boundaries",
  "",
  "- Random seed: 20260728.",
  "- Module membership was hashed before and after the run.",
  "- Candidate genes were removed from GSE161277 annotation signatures to avoid",
  "  circular annotation.",
  "- Every figure has a saved source-data table.",
  "- Raw files were read only and not decompressed or overwritten.",
  "- No trajectory, cell-communication, prognostic model or machine learning was",
  "  performed.",
  "- Stage 8 was not started.",
  "",
  "## Limitations",
  "",
  "GSE161277 has only three fully matched donors, so confidence intervals are",
  "necessarily wide and LODO is stringent. GSE132465 provides stronger paired",
  "tumor–normal and cell-specificity evidence but does not contain adenomas.",
  "Source annotations and platform/population differences introduce residual",
  "heterogeneity. Results must be described as replication or non-replication",
  "of locked computational modules, not as causal mechanisms."
)
writeLines(
  report, file.path(project_dir, "reports", "stage_7_singlecell_replication.md")
)

status_path <- file.path(project_dir, "STATUS.md")
status <- readLines(status_path, warn = FALSE)
stage_block <- c(
  "",
  "## Stage 7",
  "",
  "- Stage 7 authorization: granted 2026-07-28",
  "- Stage 7 execution: complete on server; awaiting independent Codex QC",
  paste0("- Stage 7 replication gate: ", if (isTRUE(overall)) "passed" else "not passed"),
  paste0("- Stage 7 replicated modules: ", class_counts[["replicated"]]),
  paste0(
    "- Stage 7 directionally consistent but underpowered modules: ",
    class_counts[["directionally consistent but underpowered"]]
  ),
  paste0("- Stage 7 contradictory modules: ", class_counts[["contradictory"]]),
  "- Evidence hierarchy: secondary/exploratory; primary Stage 6A negative result frozen",
  "- Report: `reports/stage_7_singlecell_replication.md`",
  "- Completion marker: `logs/07_singlecell_replication/READY_FOR_CODEX_QC`",
  "- Stage 8: not started"
)
if (any(status == "## Stage 7")) {
  status <- status[seq_len(match("## Stage 7", status) - 1L)]
}
writeLines(c(status, stage_block), status_path)
cat("STAGE_7_REPORT_OK\n")
