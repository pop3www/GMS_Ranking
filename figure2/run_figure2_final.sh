#!/usr/bin/env bash
set -euo pipefail

# Figure 2 final build.
# Input is RSEM expected counts (fractional raw counts), NOT TPM.
# The script normalizes counts once by median-ratio size factors and uses that
# normalized matrix for all expression-derived Figure 2 panels.

COUNTS=${1:-data/RawCountFile_rsemgenes.tsv}
OUT=${2:-results_fig2_teto19_weighted_FINAL}

: "${FIMO_EBOX_DIR:=out/fimo_ebox}"
: "${FIMO_TETO_DIR:=out/fimo_teto19_weighted}"
: "${TETO_MOTIF:=motifs/tetO_TRE3Gs19_weighted.meme}"
: "${POS_WEIGHTS:=motifs/tetO_positional_weights_19.tsv}"
: "${FIG2A_BOOTSTRAP:=5000}"
: "${FIG2_MIN_TOTAL:=80}"

export FIMO_EBOX_DIR FIMO_TETO_DIR TETO_MOTIF POS_WEIGHTS

mkdir -p "$OUT"
EXPR_DIR="$OUT/FIGURE2_EXPRESSION_INPUTS"
EXPR="$EXPR_DIR/Fig2_expression_normcounts.tsv"

printf '[run_figure2_final v40]\n'
printf '  COUNTS=%s\n' "$COUNTS"
printf '  OUT=%s\n' "$OUT"
printf '  FIMO_EBOX_DIR=%s\n' "$FIMO_EBOX_DIR"
printf '  FIMO_TETO_DIR=%s\n' "$FIMO_TETO_DIR"
printf '  FIG2A_BOOTSTRAP=%s\n' "$FIG2A_BOOTSTRAP"
printf '  FIG2_MIN_TOTAL=%s\n' "$FIG2_MIN_TOTAL"

# Remove only generated Figure-2 outputs to avoid schema mixing; keep raw FIMO intact.
rm -rf "$OUT"/fimo_summary_4h_* "$OUT"/priming_compare_DTvsD_vs_TamvsCtrl "$OUT"/FIGURE2_FINAL "$OUT"/FIGURE2_AUXILIARY "$EXPR_DIR"
mkdir -p "$OUT" "$EXPR_DIR"

# Prepare normalized expected-count matrix for all expression-derived Figure 2 panels.
Rscript 05_prepare_fig2_expression_matrix.R \
  --counts "$COUNTS" \
  --id-col gene_id \
  --min-total "$FIG2_MIN_TOTAL" \
  --out "$EXPR_DIR"

if [[ ! -s "$EXPR" ]]; then
  echo "[ERROR] normalized Figure 2 expression matrix was not created: $EXPR" >&2
  exit 2
fi

# A: delta-density controls on normalized RSEM expected-count scale.
Rscript 12_fig2A_delta_density_controls.R \
  --expr "$EXPR" \
  --id-col gene_id \
  --input-scale normready \
  --bootstrap "$FIG2A_BOOTSTRAP" \
  --out "$OUT"

# B-F source summaries. run_4h_pairs.sh takes a generic expression matrix even though legacy flag names say TPM.
EXPR_SCALE=normcounts bash run_4h_pairs.sh "$EXPR" "$OUT"

# B/C priming comparison.
Rscript 11_priming_amplification_compare.R \
  --rankA "$OUT/fimo_summary_4h_DT_vs_D/rank_bands.csv" \
  --rankB "$OUT/fimo_summary_4h_Tam_vs_Ctrl/rank_bands.csv" \
  --labelA DT_vs_D \
  --labelB Tam_vs_Ctrl \
  --band-source B \
  --motif-hits "$OUT/fimo_summary_4h_DT_vs_D/fimo_hits_with_bands.csv" \
  --out "$OUT/priming_compare_DTvsD_vs_TamvsCtrl"

# final copy/deposit.
bash export_fig2_outputs.sh "$OUT"

# assemble publication-ready composite.
Rscript assemble_figure2_final.R --figure-dir "$OUT/FIGURE2_FINAL"

echo "[run_figure2_final] DONE"
