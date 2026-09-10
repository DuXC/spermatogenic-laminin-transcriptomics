#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_python="${SPERM_PYTHON:-python3}"
task_rscript="${SPERM_RSCRIPT:-Rscript}"
protein_dir="$task_root/12_higher_tier_evidence_plan_20260910"
extension_dir="$task_root/13_public_data_submission_20260910"
"$task_rscript" "$protein_dir/05_scripts/05_analyze_ecm_proteomics.R" "$protein_dir"
"$task_python" "$protein_dir/05_scripts/06_ecm_figure_and_qa.py"
"$task_python" "$extension_dir/05_scripts/01_common_member_sensitivity.py"
"$task_rscript" "$extension_dir/05_scripts/02_verify_common_member.R" "$extension_dir"
"$task_python" "$extension_dir/05_scripts/03_plot_common_member.py"
