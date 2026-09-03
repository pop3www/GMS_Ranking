#!/usr/bin/env bash
set -euo pipefail

# Install and validate the count/TPM sheets from the current directory or a supplied source directory.
# Run from repo root. Example:
#   bash core/scripts/install_expression_inputs.sh /path/to/downloaded/files

SRC_DIR="${1:-$(pwd)}"
RAW_DIR="data/raw/rsem"
PROC_DIR="data/processed"
mkdir -p "$RAW_DIR" "$PROC_DIR"

copy_if_exists() {
  local src="$1"
  local dst="$2"
  if [[ -f "$src" ]]; then
    cp -p "$src" "$dst"
    echo "Copied $src -> $dst"
  else
    echo "Missing optional/input file: $src" >&2
  fi
}

copy_if_exists "$SRC_DIR/RawCountFile_rsemgenes.txt" "$RAW_DIR/RawCountFile_rsemgenes.txt"
copy_if_exists "$SRC_DIR/TPMCountFile_rsemgenes.txt" "$RAW_DIR/TPMCountFile_rsemgenes.txt"
copy_if_exists "$SRC_DIR/TPMCountFile_rsemgenes.csv" "$RAW_DIR/TPMCountFile_rsemgenes.csv"
copy_if_exists "$SRC_DIR/David_Levens_CS033646_32RNA_012523_HT3VMDRX2.xlsx" "$RAW_DIR/David_Levens_CS033646_32RNA_012523_HT3VMDRX2.xlsx"

TPM_INPUT="$RAW_DIR/TPMCountFile_rsemgenes.txt"
if [[ ! -f "$TPM_INPUT" && -f "$RAW_DIR/TPMCountFile_rsemgenes.csv" ]]; then
  TPM_INPUT="$RAW_DIR/TPMCountFile_rsemgenes.csv"
fi

python core/scripts/validate_expression_inputs.py \
  --counts "$RAW_DIR/RawCountFile_rsemgenes.txt" \
  --tpm "$TPM_INPUT" \
  --sample-sheet config/sample_sheet.csv \
  --out-dir "$PROC_DIR" \
  --write-standardized

echo
cat "$PROC_DIR/expression_input_validation.md"
