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
# Source this file before running Figure 6 locally or via SLURM.
# Existing exported variables are respected; unset variables receive repo-native defaults.

export GMS_ROOT="${GMS_ROOT:-${REPO_ROOT}}"
export FIG6_DIR="${FIG6_DIR:-$GMS_ROOT/figure5}"

# Primary expression contract for the revised manuscript:
# count-derived expression, not TPM, is the default for Figure 6.
export SAMPLE_SHEET="${SAMPLE_SHEET:-$GMS_ROOT/config/sample_sheet.csv}"
export FIG6_COUNTS_PATH="${FIG6_COUNTS_PATH:-$GMS_ROOT/data/processed/raw_counts_rsemgenes.tsv}"
export FIG6_COUNT_DERIVED_MATRIX="${FIG6_COUNT_DERIVED_MATRIX:-$GMS_ROOT/data/processed/log2_sf_norm_rsem_expected_counts.tsv}"
export FIG6_TPM_MATRIX="${FIG6_TPM_MATRIX:-$GMS_ROOT/data/processed/TPMCountFile_rsemgenes.csv}"
export FIG6_SIZE_FACTOR_PATH="${FIG6_SIZE_FACTOR_PATH:-$GMS_ROOT/data/processed/fig6_rsem_expected_count_size_factors.tsv}"
export FIG6_MIN_TOTAL_EXPECTED_COUNT="${FIG6_MIN_TOTAL_EXPECTED_COUNT:-80}"
export FIG6_PSEUDOCOUNT="${FIG6_PSEUDOCOUNT:-1}"
export FIG6_SIZE_FACTOR_METHOD="${FIG6_SIZE_FACTOR_METHOD:-median_ratio}"
export FIG6_MATRIX_SCALE="${FIG6_MATRIX_SCALE:-log2_size_factor_normalized_rsem_expected_counts}"
export FIG6_AUTO_BUILD_COUNT_MATRIX="${FIG6_AUTO_BUILD_COUNT_MATRIX:-1}"

# Unless explicitly overridden, use the count-derived primary matrix without an extra log1p.
case "${FIG6_MATRIX_SCALE}" in
  log2_size_factor_normalized_rsem_expected_counts|log2_sf_norm_rsem_expected_counts)
    export MATRIX_PATH="${MATRIX_PATH:-$FIG6_COUNT_DERIVED_MATRIX}"
    export FIG6_APPLY_LOG1P="${FIG6_APPLY_LOG1P:-0}"
    ;;
  tpm_log1p|log1p_tpm_sensitivity)
    export MATRIX_PATH="${MATRIX_PATH:-$FIG6_TPM_MATRIX}"
    export FIG6_APPLY_LOG1P="${FIG6_APPLY_LOG1P:-1}"
    ;;
  custom_pretransformed)
    : "${MATRIX_PATH:?Set MATRIX_PATH for custom_pretransformed mode}"
    export FIG6_APPLY_LOG1P="${FIG6_APPLY_LOG1P:-0}"
    ;;
  custom_log1p)
    : "${MATRIX_PATH:?Set MATRIX_PATH for custom_log1p mode}"
    export FIG6_APPLY_LOG1P="${FIG6_APPLY_LOG1P:-1}"
    ;;
  *)
    echo "[WARN] Unknown FIG6_MATRIX_SCALE=${FIG6_MATRIX_SCALE}; using MATRIX_PATH as supplied." >&2
    export MATRIX_PATH="${MATRIX_PATH:-$FIG6_COUNT_DERIVED_MATRIX}"
    export FIG6_APPLY_LOG1P="${FIG6_APPLY_LOG1P:-0}"
    ;;
esac

export DDR_GENE_LIST="${DDR_GENE_LIST:-$FIG6_DIR/source_data/ddr_genes_frozen.txt}"

# Runtime / outputs
export FIG6_ENV_NAME="${FIG6_ENV_NAME:-fig6_repro}"
export PYTHON_BIN="${PYTHON_BIN:-python}"
export FIG6_OUT="${FIG6_OUT:-$FIG6_DIR/outputs/$(date +%Y%m%d_%H%M%S)}"

# Reproducibility / plotting / thread control
export MPLBACKEND="${MPLBACKEND:-Agg}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-4}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-4}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-4}"

# Optional pipeline overrides
export FIG6_MODELS="${FIG6_MODELS:-ridge_logistic gbt_svd}"
export FIG6_FIGURE_MODEL="${FIG6_FIGURE_MODEL:-ridge_logistic}"
export FIG6_MAX_GENES="${FIG6_MAX_GENES:-500}"
export FIG6_N_PERMUTATIONS="${FIG6_N_PERMUTATIONS:-200}"
export FIG6_MANUSCRIPT_CV_SCHEME="${FIG6_MANUSCRIPT_CV_SCHEME:-leave_one_sample_out}"

mkdir -p "$FIG6_DIR/outputs" "$FIG6_DIR/logs" "$FIG6_DIR/source_data" "$GMS_ROOT/data/processed"
