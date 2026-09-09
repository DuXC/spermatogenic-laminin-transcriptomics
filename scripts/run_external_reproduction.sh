#!/usr/bin/env bash
set -euo pipefail
SPERM_ROOT=$(cd "$(dirname "$0")/.." && pwd)
SPERM_PYTHON=${SPERM_PYTHON:-python3}
SPERM_RSCRIPT=${SPERM_RSCRIPT:-Rscript}
cd "$SPERM_ROOT/10_external_validation_20260910"
mkdir -p logs 04_results/bulk 04_results/single_cell 05_figures ../env/R_library
case "${1:-verify}" in
 prepare) "$SPERM_PYTHON" scripts/02_prepare_bulk.py ;;
 bulk) "$SPERM_RSCRIPT" scripts/03_bulk_models.R ;;
 aggregate) "$SPERM_PYTHON" scripts/04_prepare_single_cell.py; "$SPERM_PYTHON" scripts/05_join_single_cell.py ;;
 pseudobulk) "$SPERM_RSCRIPT" scripts/06_single_cell_models.R ;;
 verify) "$SPERM_PYTHON" scripts/07_verify_numerical.py ;;
 figures) "$SPERM_PYTHON" scripts/08_figures.py ;;
 all) "$SPERM_PYTHON" scripts/02_prepare_bulk.py; "$SPERM_RSCRIPT" scripts/03_bulk_models.R; "$SPERM_RSCRIPT" scripts/06_single_cell_models.R; "$SPERM_PYTHON" scripts/07_verify_numerical.py; "$SPERM_PYTHON" scripts/08_figures.py ;;
 *) echo "Usage: $0 {prepare|bulk|aggregate|pseudobulk|verify|figures|all}" >&2; exit 2 ;;
esac
