#!/usr/bin/env bash
# Stage 10C background wrapper. Launch once, write terminal marker, and stop.
set -Eeuo pipefail

if (($# != 2)); then
  printf 'Usage: 10C_run_server.sh PROJECT_ROOT RUN_ID\n' >&2
  exit 2
fi

project_root="$1"
run_id="$2"
stage='stage10C'
status_dir="$project_root/status"
result_dir="$project_root/results/10C_fap_confounding/$run_id"
log_dir="$project_root/logs/10C_fap_confounding/$run_id"
report_dir="$project_root/reports"
job_id="${JOB_ID_OR_PID:-not_recorded}"
git_commit="${GIT_COMMIT:-not_recorded}"
launch_command="${LAUNCH_COMMAND:-not_recorded}"

mkdir -p "$status_dir" "$result_dir" "$log_dir" "$report_dir"

atomic_lines() {
  local target="$1"; shift
  local tmp="${target}.tmp.$$"
  printf '%s\n' "$@" > "$tmp"
  mv "$tmp" "$target"
}

lock_hashes() {
  sha256sum \
    "$project_root/results_final/stage_6A_exploratory_candidate_modules.tsv" \
    "$project_root/results_final/stage_6A_stage_blind_module_membership.tsv" \
    "$project_root/modules_locked.tsv"
}

on_error() {
  local code="$?" line="$1"
  set +e
  atomic_lines "$status_dir/$stage.FAILED" \
    "stage=$stage" "run_id=$run_id" "time=$(date --iso-8601=seconds)" \
    "exit_code=$code" "line=$line" "job_id_or_pid=$job_id" \
    "log=$log_dir/full.log" "preserved_outputs=$result_dir" \
    "recovery=NEW_RUN_ID=YYYYMMDD_HHMMSS JOB_ID_OR_PID=RETRY GIT_COMMIT=$git_commit bash scripts/10C_run_server.sh '$project_root' \"\$NEW_RUN_ID\""
  rm -f "$status_dir/$stage.RUNNING"
  exit "$code"
}
trap 'on_error $LINENO' ERR

for tool in Rscript sha256sum git; do
  command -v "$tool" >/dev/null 2>&1 || { printf 'Missing tool: %s\n' "$tool" >&2; exit 127; }
done

[[ ! -e "$status_dir/$stage.RUNNING" && ! -e "$status_dir/$stage.SUCCESS" && ! -e "$status_dir/$stage.FAILED" ]] || {
  printf 'Existing Stage 10C marker prevents unsafe overwrite.\n' >&2
  exit 3
}

lock_text="$(lock_hashes)"
atomic_lines "$status_dir/$stage.RUNNING" \
  "stage=$stage" "run_id=$run_id" "time=$(date --iso-8601=seconds)" \
  "scheduler=tmux" "job_id_or_pid=$job_id" "git_commit=$git_commit" \
  "command=$launch_command" "expected_outputs=$result_dir" "log=$log_dir/full.log"

cat > "$status_dir/${stage}_SUBMISSION.tsv.tmp.$$" <<EOF
stage	run_id	scheduler	job_id_or_pid	launch_time	launch_command	wrapper_script	git_commit	input_lock_sha256	expected_outputs	full_log	status_check_command	failure_recovery_command	resource_request	smoke_test
$stage	$run_id	tmux	$job_id	$(date --iso-8601=seconds)	$launch_command	scripts/10C_run_server.sh	$git_commit	$(printf '%s' "$lock_text" | tr '\n' ';')	$result_dir	$log_dir/full.log	test -f '$status_dir/$stage.SUCCESS' && echo SUCCESS || test -f '$status_dir/$stage.FAILED' && echo FAILED || test -f '$status_dir/$stage.RUNNING' && echo RUNNING || echo NOT_FOUND	NEW_RUN_ID=YYYYMMDD_HHMMSS bash scripts/10C_run_server.sh '$project_root' \"\$NEW_RUN_ID\"	4 CPU; <=32 GiB RAM; no GPU	$result_dir/smoke_test/SMOKE_TEST.tsv
EOF
mv "$status_dir/${stage}_SUBMISSION.tsv.tmp.$$" "$status_dir/${stage}_SUBMISSION.tsv"

exec > >(tee -a "$log_dir/full.log") 2>&1
printf '[%s] Stage 10C started\n' "$(date --iso-8601=seconds)"
printf '%s\n' "$lock_text" > "$result_dir/stage10C_input_lock_sha256.txt"

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

Rscript "$project_root/scripts/10C_fap_confounding.R" "$project_root" "$run_id" full

required=(
  GSE201348_patient_sample_audit.tsv
  GSE201348_pseudobulk_manifest.tsv
  GSE201348_gene_results.tsv
  stage10C_locked_module_results.tsv
  GSE201348_module_LODO.tsv
  GSE161277_three_paired_donor_results.tsv
  GSE200997_cancer_endpoint_audit.tsv
  stage10C_validation_checks.tsv
  stage10C_software_versions.tsv
)
for file in "${required[@]}"; do
  [[ -s "$result_dir/$file" ]] || { printf 'Missing required output: %s\n' "$file" >&2; exit 4; }
done
[[ -s "$report_dir/${stage}_SUMMARY.md" && -s "$report_dir/${stage}_GATE_DECISION.md" ]]
[[ "$(wc -l < "$report_dir/${stage}_SUMMARY.md")" -le 200 ]]
[[ "$(wc -l < "$report_dir/${stage}_GATE_DECISION.md")" -le 200 ]]
awk -F '\t' 'NR>1 && $2!="PASS"{bad=1} END{exit bad}' "$result_dir/stage10C_validation_checks.tsv"
lock_hashes > "$result_dir/stage10C_output_lock_sha256.txt"
cmp "$result_dir/stage10C_input_lock_sha256.txt" "$result_dir/stage10C_output_lock_sha256.txt"

atomic_lines "$status_dir/$stage.SUCCESS" \
  "stage=$stage" "run_id=$run_id" "time=$(date --iso-8601=seconds)" \
  "job_id_or_pid=$job_id" "git_commit=$git_commit" \
  "validation=all_stage10C_checks_PASS" \
  "summary=$report_dir/${stage}_SUMMARY.md" \
  "gate=$report_dir/${stage}_GATE_DECISION.md" \
  "results=$result_dir"
rm -f "$status_dir/$stage.RUNNING"
printf '[%s] Stage 10C SUCCESS\n' "$(date --iso-8601=seconds)"
