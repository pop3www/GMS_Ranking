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
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_FILE_DEFAULT="${SCRIPT_DIR}/fig6_env.sh"

if [[ -f "${ENV_FILE_DEFAULT}" ]]; then
  # shellcheck disable=SC1090
  source "${ENV_FILE_DEFAULT}"
fi

PYTHON_BIN="${PYTHON_BIN:-python}"

usage() {
  cat >&2 <<EOF
Usage: $0 [matrix.csv/tsv] [output_dir] [sample_sheet.csv]

Default Figure 6 primary input is count-derived expression:
  data/processed/log2_sf_norm_rsem_expected_counts.tsv
built from:
  data/processed/raw_counts_rsemgenes.tsv
when FIG6_AUTO_BUILD_COUNT_MATRIX=1.

TPM sensitivity example:
  MATRIX_PATH=\$GMS_ROOT/data/processed/TPMCountFile_rsemgenes.csv \\
  FIG6_MATRIX_SCALE=tpm_log1p FIG6_APPLY_LOG1P=1 bash ./run_fig6_bundle.sh

Resolution order:
  1) explicit CLI args
  2) exported env vars from ./fig6_env.sh
  3) repo-native curated sample sheet: ${REPO_ROOT}/config/sample_sheet.csv

Optional env vars:
  GMS_ROOT, FIG6_DIR, FIG6_COUNTS_PATH, FIG6_COUNT_DERIVED_MATRIX,
  MATRIX_PATH, FIG6_MATRIX_SCALE, FIG6_APPLY_LOG1P, FIG6_AUTO_BUILD_COUNT_MATRIX,
  FIG6_OUT, SAMPLE_SHEET, DDR_GENE_LIST, FIG6_MODELS, FIG6_FIGURE_MODEL,
  FIG6_MAX_GENES, FIG6_N_PERMUTATIONS, FORCE_LEGACY_SAMPLE_SHEET=1,
  SKIP_GENE_SCREEN=1, SKIP_MANUSCRIPT_ASSETS=1, SKIP_ENV_PREFLIGHT=1
EOF
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

MATRIX_PATH="${1:-${MATRIX_PATH:-}}"
if [[ -z "${MATRIX_PATH}" ]]; then
  usage
  exit 1
fi

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="${2:-${FIG6_OUT:-${OUTPUT_DIR:-${SCRIPT_DIR}/outputs/${RUN_STAMP}}}}"
SAMPLE_SHEET_ARG="${3:-${SAMPLE_SHEET:-}}"
DEFAULT_CURATED_SAMPLE_SHEET="${REPO_ROOT}/config/sample_sheet.csv"
DEFAULT_DDR_GENE_LIST="${SCRIPT_DIR}/source_data/ddr_genes_frozen.txt"
FORCE_LEGACY_SAMPLE_SHEET="${FORCE_LEGACY_SAMPLE_SHEET:-0}"
SKIP_GENE_SCREEN="${SKIP_GENE_SCREEN:-0}"
SKIP_MANUSCRIPT_ASSETS="${SKIP_MANUSCRIPT_ASSETS:-0}"
SKIP_ENV_PREFLIGHT="${SKIP_ENV_PREFLIGHT:-0}"
FIG6_FIGURE_MODEL="${FIG6_FIGURE_MODEL:-ridge_logistic}"
FIG6_MAX_GENES="${FIG6_MAX_GENES:-500}"
FIG6_N_PERMUTATIONS="${FIG6_N_PERMUTATIONS:-200}"
FIG6_MANUSCRIPT_CV_SCHEME="${FIG6_MANUSCRIPT_CV_SCHEME:-leave_one_sample_out}"
FIG6_MATRIX_SCALE="${FIG6_MATRIX_SCALE:-unspecified}"
FIG6_APPLY_LOG1P="${FIG6_APPLY_LOG1P:-0}"
FIG6_AUTO_BUILD_COUNT_MATRIX="${FIG6_AUTO_BUILD_COUNT_MATRIX:-1}"
FIG6_COUNTS_PATH="${FIG6_COUNTS_PATH:-${REPO_ROOT}/data/processed/raw_counts_rsemgenes.tsv}"
FIG6_COUNT_DERIVED_MATRIX="${FIG6_COUNT_DERIVED_MATRIX:-${REPO_ROOT}/data/processed/log2_sf_norm_rsem_expected_counts.tsv}"
FIG6_SIZE_FACTOR_METHOD="${FIG6_SIZE_FACTOR_METHOD:-median_ratio}"
FIG6_MIN_TOTAL_EXPECTED_COUNT="${FIG6_MIN_TOTAL_EXPECTED_COUNT:-80}"
FIG6_PSEUDOCOUNT="${FIG6_PSEUDOCOUNT:-1}"
IFS=' ' read -r -a FIG6_MODELS_ARRAY <<< "${FIG6_MODELS:-ridge_logistic gbt_svd}"

mkdir -p "${OUTPUT_DIR}" "${SCRIPT_DIR}/logs" "${SCRIPT_DIR}/outputs" "${SCRIPT_DIR}/source_data"

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  echo "[ERROR] Python executable not found: ${PYTHON_BIN}" >&2
  exit 1
fi

if [[ -n "${SAMPLE_SHEET_ARG}" ]]; then
  RESOLVED_SAMPLE_SHEET="${SAMPLE_SHEET_ARG}"
elif [[ -f "${DEFAULT_CURATED_SAMPLE_SHEET}" ]]; then
  RESOLVED_SAMPLE_SHEET="${DEFAULT_CURATED_SAMPLE_SHEET}"
elif [[ "${FORCE_LEGACY_SAMPLE_SHEET}" == "1" ]]; then
  RESOLVED_SAMPLE_SHEET="${OUTPUT_DIR}/legacy_sample_sheet.csv"
else
  echo "[ERROR] Curated sample sheet not found: ${DEFAULT_CURATED_SAMPLE_SHEET}" >&2
  echo "[ERROR] Refusing unsafe column-order metadata inference." >&2
  echo "[ERROR] To force compatibility-mode behavior, set FORCE_LEGACY_SAMPLE_SHEET=1 explicitly." >&2
  exit 1
fi

if [[ ! -f "${RESOLVED_SAMPLE_SHEET}" && "${RESOLVED_SAMPLE_SHEET}" != "${OUTPUT_DIR}/legacy_sample_sheet.csv" ]]; then
  echo "[ERROR] Sample sheet not found: ${RESOLVED_SAMPLE_SHEET}" >&2
  exit 1
fi

if [[ ! -f "${MATRIX_PATH}" ]]; then
  if [[ "${FIG6_AUTO_BUILD_COUNT_MATRIX}" == "1" && "${MATRIX_PATH}" == "${FIG6_COUNT_DERIVED_MATRIX}" ]]; then
    if [[ ! -f "${FIG6_COUNTS_PATH}" ]]; then
      echo "[ERROR] Count-derived matrix is missing and count source was not found: ${FIG6_COUNTS_PATH}" >&2
      exit 1
    fi
    echo "[INFO] Building count-derived Figure 6 matrix: ${MATRIX_PATH}"
    "${PYTHON_BIN}" "${SCRIPT_DIR}/build_count_expression_matrix.py" \
      --counts "${FIG6_COUNTS_PATH}" \
      --sample-sheet "${RESOLVED_SAMPLE_SHEET}" \
      --output "${MATRIX_PATH}" \
      --min-total-expected-count "${FIG6_MIN_TOTAL_EXPECTED_COUNT}" \
      --pseudocount "${FIG6_PSEUDOCOUNT}" \
      --size-factor-method "${FIG6_SIZE_FACTOR_METHOD}"
  else
    echo "[ERROR] Matrix file not found: ${MATRIX_PATH}" >&2
    exit 1
  fi
fi

if [[ "${RESOLVED_SAMPLE_SHEET}" == "${OUTPUT_DIR}/legacy_sample_sheet.csv" ]]; then
  echo "[WARN] FORCE_LEGACY_SAMPLE_SHEET=1: building compatibility sample sheet from matrix names."
  "${PYTHON_BIN}" "${SCRIPT_DIR}/build_legacy_sample_sheet.py" \
    --matrix "${MATRIX_PATH}" \
    --output "${RESOLVED_SAMPLE_SHEET}"
else
  echo "[INFO] Using curated sample sheet: ${RESOLVED_SAMPLE_SHEET}"
fi

if [[ -n "${DDR_GENE_LIST:-}" ]]; then
  RESOLVED_DDR_GENE_LIST="${DDR_GENE_LIST}"
elif [[ -f "${DEFAULT_DDR_GENE_LIST}" ]]; then
  RESOLVED_DDR_GENE_LIST="${DEFAULT_DDR_GENE_LIST}"
else
  RESOLVED_DDR_GENE_LIST=""
fi

if [[ -n "${RESOLVED_DDR_GENE_LIST}" && ! -f "${RESOLVED_DDR_GENE_LIST}" ]]; then
  echo "[WARN] DDR gene list requested but not found: ${RESOLVED_DDR_GENE_LIST}" >&2
  echo "[WARN] Falling back to built-in DEFAULT_DDR_GENES inside fig6_ml_core.py" >&2
  RESOLVED_DDR_GENE_LIST=""
fi

MATRIX_EXTRA_ARGS=(--matrix-scale "${FIG6_MATRIX_SCALE}")
if [[ "${FIG6_APPLY_LOG1P}" == "1" ]]; then
  MATRIX_EXTRA_ARGS+=(--log1p)
fi

if [[ ( "${FIG6_MATRIX_SCALE}" == "log2_sf_norm_rsem_expected_counts" || "${FIG6_MATRIX_SCALE}" == "log2_size_factor_normalized_rsem_expected_counts" ) && "${FIG6_APPLY_LOG1P}" == "1" ]]; then
  echo "[ERROR] FIG6_APPLY_LOG1P=1 is invalid for pre-transformed count-derived matrix scale." >&2
  exit 1
fi

echo "[INFO] Repo root        : ${REPO_ROOT}"
echo "[INFO] Figure6 dir      : ${SCRIPT_DIR}"
echo "[INFO] Matrix           : ${MATRIX_PATH}"
echo "[INFO] Matrix scale     : ${FIG6_MATRIX_SCALE}"
echo "[INFO] Apply log1p      : ${FIG6_APPLY_LOG1P}"
echo "[INFO] Count source     : ${FIG6_COUNTS_PATH}"
echo "[INFO] Sample sheet     : ${RESOLVED_SAMPLE_SHEET}"
echo "[INFO] DDR gene list    : ${RESOLVED_DDR_GENE_LIST:-<built-in default>}"
echo "[INFO] Output dir       : ${OUTPUT_DIR}"
echo "[INFO] Models           : ${FIG6_MODELS_ARRAY[*]}"
echo "[INFO] Figure model     : ${FIG6_FIGURE_MODEL}"
echo "[INFO] Max genes        : ${FIG6_MAX_GENES}"
echo "[INFO] N permutations   : ${FIG6_N_PERMUTATIONS}"

if [[ "${SKIP_ENV_PREFLIGHT}" != "1" ]]; then
  echo "[INFO] Running environment preflight"
  "${PYTHON_BIN}" "${SCRIPT_DIR}/preflight_fig6_env.py"
else
  echo "[INFO] SKIP_ENV_PREFLIGHT=1, skipping environment preflight"
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/run_fig6_sample_level.py" \
  --matrix "${MATRIX_PATH}" \
  --sample-sheet "${RESOLVED_SAMPLE_SHEET}" \
  --output-dir "${OUTPUT_DIR}" \
  "${MATRIX_EXTRA_ARGS[@]}" \
  --models "${FIG6_MODELS_ARRAY[@]}" \
  --figure-model "${FIG6_FIGURE_MODEL}" \
  --max-genes "${FIG6_MAX_GENES}" \
  --n-permutations "${FIG6_N_PERMUTATIONS}"

if [[ "${SKIP_GENE_SCREEN}" != "1" ]]; then
  GENE_SCREEN_ARGS=(
    --matrix "${MATRIX_PATH}"
    --sample-sheet "${RESOLVED_SAMPLE_SHEET}"
    --output-dir "${OUTPUT_DIR}/gene_screen"
    "${MATRIX_EXTRA_ARGS[@]}"
  )
  if [[ -n "${RESOLVED_DDR_GENE_LIST}" ]]; then
    GENE_SCREEN_ARGS+=( --ddr-genes "${RESOLVED_DDR_GENE_LIST}" )
  fi
  "${PYTHON_BIN}" "${SCRIPT_DIR}/fig6_gene_screen.py" "${GENE_SCREEN_ARGS[@]}"
else
  echo "[INFO] SKIP_GENE_SCREEN=1, skipping fig6_gene_screen.py"
fi

if [[ "${SKIP_MANUSCRIPT_ASSETS}" != "1" ]]; then
  "${PYTHON_BIN}" "${SCRIPT_DIR}/build_fig6_manuscript_assets.py" \
    --run-dir "${OUTPUT_DIR}" \
    --model "${FIG6_FIGURE_MODEL}" \
    --cv-scheme "${FIG6_MANUSCRIPT_CV_SCHEME}" \
    --include-control-summary
else
  echo "[INFO] SKIP_MANUSCRIPT_ASSETS=1, skipping build_fig6_manuscript_assets.py"
fi

echo "[OK] Figure 6 run finished: ${OUTPUT_DIR}"
