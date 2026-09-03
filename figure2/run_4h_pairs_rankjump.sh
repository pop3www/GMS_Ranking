#!/usr/bin/env bash
# 4h pairwise runner (rankjump + motif summaries)

set -euo pipefail

SCRIPT_VERSION="2026-01-09 v10 (env fallbacks for step10; stable defaults)"

TPM=${1:-}
OUTROOT=${2:-}
if [[ -z "${TPM}" || -z "${OUTROOT}" ]]; then
  echo "Usage: $0 <TPMCountFile_rsemgenes.txt> <outroot>" >&2
  exit 2
fi

mkdir -p "${OUTROOT}"

# ---- configurable knobs (via env) ----
: "${DE_METHOD:=limma}"          # limma | none
: "${FDR_THRESH:=0.05}"
: "${LFC_THRESH:=0.32}"
: "${PVAL_MAX:=1e-6}"            # FIMO p-value cutoff
: "${QVAL_MAX:=1}"               # FIMO q-value cutoff
: "${HEAD_FRAC:=0.10}"
: "${TAIL_FRAC:=0.10}"
: "${JUMP_FRAC:=0.01}"           # top fraction for rankjump
: "${MIN_BASELINE_TPM:=1}"

: "${FIMO_EBOX_DIR:=out/fimo_ebox}"
: "${FIMO_TETO_DIR:=""}"  # recommend setting explicitly
: "${TETO_MOTIF:=""}"     # used only by 08_build_supp_fimo_figs.R
: "${POS_WEIGHTS:=""}"    # used only by 08_build_supp_fimo_figs.R

echo "[run_4h_pairs_rankjump] script version: ${SCRIPT_VERSION}" >&2

check_fimo_dir() {
  local d="$1"
  if [[ -z "${d}" ]]; then
    return 1
  fi
  if [[ -f "${d}/fimo.tsv" || -f "${d}/fimo.txt" ]]; then
    return 0
  fi
  return 1
}

if ! check_fimo_dir "${FIMO_EBOX_DIR}"; then
  echo "[run_4h_pairs_rankjump] WARNING: no fimo.tsv or fimo.txt in: ${FIMO_EBOX_DIR}" >&2
fi
if [[ -n "${FIMO_TETO_DIR}" ]] && ! check_fimo_dir "${FIMO_TETO_DIR}"; then
  echo "[run_4h_pairs_rankjump] WARNING: no fimo.tsv or fimo.txt in: ${FIMO_TETO_DIR}" >&2
fi

run_one() {
  local name="$1"
  local A_rx="$2"
  local B_rx="$3"

  local OUT="${OUTROOT}/fimo_summary_4h_${name}"
  mkdir -p "${OUT}"

  echo "==> [4h] ${name} -> ${OUT}"

  # 06: bands + DE calls
  Rscript 06_make_rank_bands_from_TPM.R \
    --tpm "${TPM}" \
    --id-col gene_id \
    --baseline-pattern "${B_rx}" \
    --groupA-pattern "${A_rx}" \
    --groupB-pattern "${B_rx}" \
    --head-frac "${HEAD_FRAC}" \
    --tail-frac "${TAIL_FRAC}" \
    --de-method "${DE_METHOD}" \
    --fdr-thresh "${FDR_THRESH}" \
    --lfc-thresh "${LFC_THRESH}" \
    --out "${OUT}"

  # 07: FIMO enrichment summaries
  local args07=(
    --bands "${OUT}/rank_bands.csv"
    --id-col gene
    --pval-max "${PVAL_MAX}"
    --qval-max "${QVAL_MAX}"
    --out "${OUT}"
  )

  if check_fimo_dir "${FIMO_EBOX_DIR}"; then
    args07+=( --fimo-ebox "${FIMO_EBOX_DIR}" )
  fi
  if check_fimo_dir "${FIMO_TETO_DIR}"; then
    args07+=( --fimo-teto "${FIMO_TETO_DIR}" )
  fi

  Rscript 07_summarize_fimo.R "${args07[@]}"

  # 08: supplemental plots (logo + weights + density + ORs)
  # Only run if motif/weights paths are provided and exist.
  if [[ -n "${TETO_MOTIF}" && -f "${TETO_MOTIF}" && -n "${POS_WEIGHTS}" && -f "${POS_WEIGHTS}" ]]; then
    Rscript 08_build_supp_fimo_figs.R \
      --summary-dir "${OUT}" \
      --motif-meme "${TETO_MOTIF}" \
      --pos-weights "${POS_WEIGHTS}" \
      --out "${OUT}"
  else
    echo "[08] motif/weights not provided; skipping 08_build_supp_fimo_figs.R" >&2
  fi

  # 10: rank-jump motif enrichment (robust: pass both CLI args and env fallbacks)
  TPM="${TPM}" ID_COL="gene_id" \
    GROUPA_PATTERN="${A_rx}" GROUPB_PATTERN="${B_rx}" \
    BANDS="${OUT}/rank_bands.csv" FIMO_HITS="${OUT}/fimo_hits_with_bands.csv" \
    LABEL="${name}" TOP_FRAC="${JUMP_FRAC}" MIN_BASELINE_TPM="${MIN_BASELINE_TPM}" OUTDIR="${OUT}" \
    Rscript 10_rankjump_motif_analysis.R \
      --tpm "${TPM}" \
      --id-col gene_id \
      --groupA-pattern "${A_rx}" \
      --groupB-pattern "${B_rx}" \
      --bands "${OUT}/rank_bands.csv" \
      --fimo-hits "${OUT}/fimo_hits_with_bands.csv" \
      --label "${name}" \
      --top-frac "${JUMP_FRAC}" \
      --min-baseline-tpm "${MIN_BASELINE_TPM}" \
      --out "${OUT}"
}

# --------- comparisons (4h) ---------
run_one "DT_vs_D"  '(^|_)4DT_RP[0-9]+($|_)' '(^|_)4D_RP[0-9]+($|_)'
