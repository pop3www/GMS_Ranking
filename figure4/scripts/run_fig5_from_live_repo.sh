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

# Run from repo root:
#   ${REPO_ROOT}
#
# Primary Figure 5 input is validated RSEM expected counts normalized by
# median-ratio size factors and transformed to log2(normalized count + 1), matching the manuscript Methods.
# TPM can be supplied explicitly for sensitivity by setting FIG5_EXPR_INPUT and
# FIG5_MATRIX_TYPE=tpm FIG5_NORMALIZATION=none.
#
# Raw validation inputs live in data/raw/, not a nested raw RSEM subdirectory.

EXPR_INPUT="${FIG5_EXPR_INPUT:-data/processed/raw_counts_rsemgenes.tsv}"
SAMPLE_SHEET="${FIG5_SAMPLE_SHEET:-config/sample_sheet.csv}"
CANONICAL_EXPR="${FIG5_CANONICAL_EXPR:-core/outputs/tables/rank_geometry_expression_gene_by_sample.csv}"
CONFIG="${FIG5_CONFIG:-config/rank_geometry.yaml}"
OUTDIR="${FIG5_OUTDIR:-figure4/outputs}"
SOURCE_DATA_DIR="${FIG5_SOURCE_DATA_DIR:-figure4/source_data}"
MATRIX_TYPE="${FIG5_MATRIX_TYPE:-raw_counts}"
NORMALIZATION="${FIG5_NORMALIZATION:-median_ratio}"
TRANSFORM="${FIG5_TRANSFORM:-auto}"
MIN_TOTAL_COUNT="${FIG5_MIN_TOTAL_COUNT:-80}"

if [[ ! -f "${EXPR_INPUT}" ]]; then
  cat >&2 <<EOF
[ERROR] Figure 4 primary expression input not found:
  ${EXPR_INPUT}

Expected primary input after validation:
  data/processed/raw_counts_rsemgenes.tsv

To regenerate standardized processed files from validated raw inputs, run:
  python core/scripts/validate_expression_inputs.py \
    --counts data/raw/RawCountFile_rsemgenes.txt \
    --tpm data/raw/TPMCountFile_rsemgenes.txt \
    --sample-sheet config/sample_sheet.csv \
    --out-dir data/processed \
    --write-standardized

Important: the raw files are in data/raw/, not a nested raw RSEM subdirectory.
For a TPM sensitivity run, set:
  FIG5_EXPR_INPUT=data/processed/TPMCountFile_rsemgenes.csv \
  FIG5_MATRIX_TYPE=tpm \
  FIG5_NORMALIZATION=none \
  FIG5_TRANSFORM=none
EOF
  exit 2
fi

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
  echo "[ERROR] Sample sheet not found: ${SAMPLE_SHEET}" >&2
  exit 2
fi

if [[ ! -f "${CONFIG}" ]]; then
  echo "[ERROR] Figure 5 config not found: ${CONFIG}" >&2
  exit 2
fi

python core/scripts/prepare_rank_geometry_expression.py \
  --input "${EXPR_INPUT}" \
  --sample-sheet "${SAMPLE_SHEET}" \
  --output "${CANONICAL_EXPR}" \
  --matrix-type "${MATRIX_TYPE}" \
  --normalization "${NORMALIZATION}" \
  --transform "${TRANSFORM}" \
  --min-total-count "${MIN_TOTAL_COUNT}"

python figure4/scripts/run.py \
  --expr "${CANONICAL_EXPR}" \
  --sample-sheet "${SAMPLE_SHEET}" \
  --config "${CONFIG}" \
  --outdir "${OUTDIR}" \
  --source-data-dir "${SOURCE_DATA_DIR}"
