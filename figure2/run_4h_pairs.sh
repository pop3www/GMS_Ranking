#!/usr/bin/env bash
set -euo pipefail

# Run 4h pairwise contrasts and summarize motif enrichment (Up vs Down) per band.
#
# Usage:
#   bash run_4h_pairs.sh <prepared_expression_matrix.tsv> out
#
# Environment overrides (optional):
#   FIMO_EBOX_DIR   default: out/fimo_ebox
#   FIMO_TETO_DIR   default: out/fimo_teto19_weighted
#   HEAD_FRAC       default: 0.10
#   TAIL_FRAC       default: 0.10
#   LFC_THRESH      default: 0.32
#   DE_METHOD       default: none   (or limma)
#   FDR_THRESH      default: 0.05   (if DE_METHOD=limma)
#   PVAL_MAX        optional: stricter FIMO p-value filter in summarization (e.g. 1e-6)
#   QVAL_MAX        optional: stricter FIMO q-value filter in summarization
#   TSS_WINDOW      optional: restrict hits to +/- bp around promoter center in summarization (requires FIMO start/stop)
#   PROMOTERS_BED   optional: BED file to estimate promoter lengths for TSS_WINDOW (default: ref/promoters_2kb.bed if exists)
#   TETO_MOTIF      default: motifs/tetO_TRE3Gs19_weighted.meme
#   POS_WEIGHTS     default: motifs/tetO_positional_weights_19.tsv

EXPR="${1:?Prepared expression matrix file required}"
OUTROOT="${2:-out}"

FIMO_EBOX_DIR="${FIMO_EBOX_DIR:-out/fimo_ebox}"
FIMO_TETO_DIR="${FIMO_TETO_DIR:-out/fimo_teto19_weighted}"

HEAD_FRAC="${HEAD_FRAC:-0.10}"
TAIL_FRAC="${TAIL_FRAC:-0.10}"
LFC_THRESH="${LFC_THRESH:-0.32}"

DE_METHOD="${DE_METHOD:-none}"
FDR_THRESH="${FDR_THRESH:-0.05}"
EXPR_SCALE="${EXPR_SCALE:-${FIG2_EXPR_SCALE:-normready}}"

PVAL_MAX="${PVAL_MAX:-}"
QVAL_MAX="${QVAL_MAX:-}"
TSS_WINDOW="${TSS_WINDOW:-}"

PROMOTERS_BED="${PROMOTERS_BED:-ref/promoters_2kb.bed}"
if [[ ! -f "$PROMOTERS_BED" ]]; then
  PROMOTERS_BED=""
fi

TETO_MOTIF="${TETO_MOTIF:-motifs/tetO_TRE3Gs19_weighted.meme}"
POS_WEIGHTS="${POS_WEIGHTS:-motifs/tetO_positional_weights_19.tsv}"

mkdir -p "$OUTROOT"

run_one () {
  local name="$1"
  local baseline_rx="$2"
  local A_rx="$3"
  local B_rx="$4"

  local OUT="$OUTROOT/fimo_summary_4h_${name}"
  mkdir -p "$OUT"

  echo "==> [4h] ${name} -> ${OUT}"

  # 06
  Rscript 06_make_rank_bands_from_TPM.R \
    --tpm "$EXPR" --id-col gene_id \
    --baseline-pattern "$baseline_rx" \
    --groupA-pattern   "$A_rx" \
    --groupB-pattern   "$B_rx" \
    --head-frac "$HEAD_FRAC" --tail-frac "$TAIL_FRAC" --lfc-thresh "$LFC_THRESH" \
    --de-method "$DE_METHOD" --fdr-thresh "$FDR_THRESH" \
    --input-scale "$EXPR_SCALE" \
    --out "$OUT"

  # 07
  args07=()
  [[ -n "$PVAL_MAX" ]] && args07+=( --pval-max "$PVAL_MAX" )
  [[ -n "$QVAL_MAX"  ]] && args07+=( --qval-max "$QVAL_MAX" )
  [[ -n "$TSS_WINDOW" ]] && args07+=( --tss-window "$TSS_WINDOW" )
  [[ -n "$PROMOTERS_BED" ]] && args07+=( --promoters "$PROMOTERS_BED" )

  Rscript 07_summarize_fimo.R \
    --fimo-ebox "$FIMO_EBOX_DIR" \
    --fimo-teto "$FIMO_TETO_DIR" \
    --bands "$OUT/rank_bands.csv" \
    --id-col gene \
    --out "$OUT" \
    "${args07[@]}"

  # 08 figures (optional but useful)
  Rscript 08_build_supp_fimo_figs.R \
    --summary-dir "$OUT" \
    --motif-meme  "$TETO_MOTIF" \
    --pos-weights "$POS_WEIGHTS" \
    --out         "$OUT"
}

# Regex patterns for 4h only
RX_4D='(^|_)4D_RP[0-9]+($|_)'
RX_4DT='(^|_)4DT_RP[0-9]+($|_)'
RX_4TAM='(^|_)4_Tam_RP[0-9]+($|_)'
RX_4CTRL='(^|_)4_ctrl_RP[0-9]+($|_)'

run_one "DT_vs_D"     "$RX_4D"    "$RX_4DT"  "$RX_4D"
run_one "DT_vs_Tam"   "$RX_4TAM"  "$RX_4DT"  "$RX_4TAM"
run_one "DT_vs_Ctrl"  "$RX_4CTRL" "$RX_4DT"  "$RX_4CTRL"
run_one "Tam_vs_Ctrl" "$RX_4CTRL" "$RX_4TAM" "$RX_4CTRL"

# Combine plots across contrasts
Rscript combine_pairs_OR.R "$OUTROOT"

echo "[DONE] Pairwise summaries under: $OUTROOT/fimo_summary_4h_*"
