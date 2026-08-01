#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT="${CRC_PROJECT_ROOT}"
RUN_ID="20260729_172640"
RESULT="$PROJECT/results/08B_bulk_validation/$RUN_ID"
AUDIT="$RESULT/audit_pre_GSE8671_pairing_correction"
LOG="$PROJECT/logs/08B_bulk_validation/stage_8B_${RUN_ID}_pairing_correction.log"
READY="$PROJECT/logs/08B_bulk_validation/READY_FOR_CODEX_QC"
ATTENTION="$PROJECT/logs/08B_bulk_validation/NEEDS_CODEX_ATTENTION"

mkdir -p "$AUDIT" "$PROJECT/figures_final" "$PROJECT/results_final" "$PROJECT/logs_summary"
rm -f "$READY" "$ATTENTION"
for file in bulk_module_scores.tsv bulk_cohort_effects.tsv bulk_validation_summary.tsv \
  meta_analysis_results.tsv meta_analysis_source_data.tsv \
  GSE8671_leave_one_pair_out.tsv GSE8671_leave_one_pair_out_summary.tsv \
  tissue_validation_gate.tsv validation_checks.tsv; do
  if [[ -f "$RESULT/$file" && ! -f "$AUDIT/$file" ]]; then
    cp -p "$RESULT/$file" "$AUDIT/$file"
  fi
done
if [[ -f "$PROJECT/reports/stage_8B_bulk_validation.md" &&
      ! -f "$AUDIT/stage_8B_bulk_validation.md" ]]; then
  cp -p "$PROJECT/reports/stage_8B_bulk_validation.md" "$AUDIT/"
fi

on_error() {
  local status=$?
  printf 'Stage 8B GSE8671 pairing correction failed (exit %s) at %s. See %s\n' \
    "$status" "$(date --iso-8601=seconds)" "$LOG" > "$ATTENTION"
  exit "$status"
}
trap on_error ERR
exec > >(tee -a "$LOG") 2>&1

echo "Pairing correction start: $(date --iso-8601=seconds)"
Rscript "$PROJECT/scripts/08B_bulk_cohorts.R" "$PROJECT" "$RUN_ID"
Rscript "$PROJECT/scripts/08B_meta_gate_report.R" "$PROJECT" "$RUN_ID"
Rscript "$PROJECT/scripts/08B_validate_outputs.R" "$PROJECT" "$RUN_ID"
echo "Pairing correction complete: $(date --iso-8601=seconds)"
trap - ERR
