#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
SPERM_PYTHON=${SPERM_PYTHON:-python3}
SPERM_RSCRIPT=${SPERM_RSCRIPT:-Rscript}
mkdir -p logs env/R_library
case "${1:-verify}" in
  verify) "$SPERM_PYTHON" scripts/verify_release_data.py ;;
  bulk) "$SPERM_PYTHON" scripts/01_prepare.py; "$SPERM_RSCRIPT" scripts/02_bulk_analysis.R ;;
  pseudobulk) "$SPERM_RSCRIPT" 07_single_cell/scripts/03_pseudobulk.R ;;
  sources) "$SPERM_PYTHON" 07_single_cell/scripts/02_annotate_aggregate.py ;;
  figures)
    "$SPERM_PYTHON" scripts/03_figures.py
    "$SPERM_PYTHON" 07_single_cell/scripts/04_figures.py
    "$SPERM_PYTHON" 08_manuscript_refinement_20260909/scripts/02_tables.py
    "$SPERM_PYTHON" 08_manuscript_refinement_20260909/scripts/03a_source_figures.py
    "$SPERM_PYTHON" 08_manuscript_refinement_20260909/scripts/03b_evidence_figures.py ;;
  *) echo 'Usage: bash scripts/run_reproduction.sh verify|bulk|pseudobulk|sources|figures'; exit 2 ;;
esac
