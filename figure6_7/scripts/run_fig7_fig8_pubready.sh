#!/usr/bin/env bash
# Locate the public repository root without relying on a machine-specific path.
_gms_find_repo_root() {
  local d
  d=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
  while [ "$d" != "/" ]; do
    if [ -f "$d/README.md" ] && [ -f "$d/FIGURE_MAP.tsv" ]; then
      printf '%s\n' "$d"
      return 0
    fi
    d=$(dirname "$d")
  done
  return 1
}
REPO_ROOT=${REPO_ROOT:-$(_gms_find_repo_root)}
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${MODULE_DIR}/.." && pwd)"

COUNTS_CSV="${1:-${COUNTS_CSV:-${REPO_ROOT}/data/processed/raw_counts_rsemgenes.tsv}}"
SAMPLE_SHEET_CSV="${2:-${SAMPLE_SHEET_CSV:-${REPO_ROOT}/config/sample_sheet.csv}}"
OUT_DIR="${3:-${OUT_DIR:-${MODULE_DIR}/outputs_pubready}}"
COUNT_INPUT_MODE="${COUNT_INPUT_MODE:-rsem_expected_counts}"
FGSEA_NPERM="${FGSEA_NPERM:-10000}"
PLOT_DPI="${PLOT_DPI:-450}"

mkdir -p "${OUT_DIR}" "${MODULE_DIR}/logs"
LOG="${MODULE_DIR}/logs/run_fig78_pubready_$(date +%Y%m%d_%H%M%S).log"

echo "[Fig7/8 pubready] repo_root=${REPO_ROOT}" | tee "${LOG}"
echo "[Fig7/8 pubready] counts=${COUNTS_CSV}" | tee -a "${LOG}"
echo "[Fig7/8 pubready] sample_sheet=${SAMPLE_SHEET_CSV}" | tee -a "${LOG}"
echo "[Fig7/8 pubready] out_dir=${OUT_DIR}" | tee -a "${LOG}"
echo "[Fig7/8 pubready] count_input_mode=${COUNT_INPUT_MODE}" | tee -a "${LOG}"
echo "[Fig7/8 pubready] fgsea_nperm=${FGSEA_NPERM}" | tee -a "${LOG}"
echo "[Fig7/8 pubready] plot_dpi=${PLOT_DPI}" | tee -a "${LOG}"

Rscript "${SCRIPT_DIR}/run.R" \
  --counts_csv="${COUNTS_CSV}" \
  --sample_sheet_csv="${SAMPLE_SHEET_CSV}" \
  --count_input_mode="${COUNT_INPUT_MODE}" \
  --out_dir="${OUT_DIR}" \
  --svg \
  --nperm="${FGSEA_NPERM}" \
  --plot_dpi="${PLOT_DPI}" \
  --fig7_composite_width=23.5 \
  --fig7_composite_height=24.5 \
  --fig8_composite_width=19.2 \
  --fig8_composite_height=13.8 \
  "${@:4}" 2>&1 | tee -a "${LOG}"

echo "[Fig7/8 pubready] complete: ${OUT_DIR}" | tee -a "${LOG}"
