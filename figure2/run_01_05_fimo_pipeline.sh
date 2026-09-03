#!/usr/bin/env bash
set -euo pipefail

SCRIPT_VERSION="2026-01-09 v12 (rebuild FIMO out/ for EBOX + tetO19; safe skip/force; no symlink traps)"

echo "[run_01_05_fimo_pipeline] version: ${SCRIPT_VERSION}" >&2

# Inputs (override via env)
PROMOTERS_FA="${PROMOTERS_FA:-data/promoters_2kb.fa}"
EBOX_MOTIF="${EBOX_MOTIF:-motifs/MA0147.3_MYCe_meme.txt}"
TETO19_MOTIF="${TETO19_MOTIF:-motifs/tetO_TRE3Gs19.meme}"
TRE7_FA="${TRE7_FA:-motifs/TRE3Gs_7x.fa}"

# Outputs (override via env)
OUTDIR="${OUTDIR:-out}"
EBOX_OUT="${EBOX_OUT:-${OUTDIR}/fimo_ebox}"
TETO19_OUT="${TETO19_OUT:-${OUTDIR}/fimo_teto19}"
POSCTRL19_OUT="${POSCTRL19_OUT:-${OUTDIR}/fimo_poscontrol_teto19}"

# FIMO thresholds (scan-time). You can still filter later in 07_summarize_fimo.R via PVAL_MAX / QVAL_MAX.
THRESH_EBOX="${THRESH_EBOX:-1e-4}"
THRESH_TETO19="${THRESH_TETO19:-1e-4}"
THRESH_POSCTRL19="${THRESH_POSCTRL19:-1e-4}"

# Controls
RUN_EBOX="${RUN_EBOX:-1}"
RUN_TETO19="${RUN_TETO19:-1}"
RUN_POSCTRL19="${RUN_POSCTRL19:-1}"

# If FORCE=1, existing output dirs are moved aside before re-running.
FORCE="${FORCE:-0}"

# ------------------------------------------------------------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "[run_01_05_fimo_pipeline] ERROR: '$1' not found in PATH" >&2; exit 2; }; }
need_file() { [[ -s "$1" ]] || { echo "[run_01_05_fimo_pipeline] ERROR: missing/empty file: $1" >&2; exit 2; }; }

need fimo
need_file "$PROMOTERS_FA"

if [[ "$RUN_EBOX" == "1" ]]; then
  need_file "$EBOX_MOTIF"
fi
if [[ "$RUN_TETO19" == "1" ]]; then
  need_file "$TETO19_MOTIF"
fi
if [[ "$RUN_POSCTRL19" == "1" ]]; then
  need_file "$TRE7_FA"
fi

mkdir -p "$OUTDIR"

# small sanity check: tetO19 motif length vs positional weights (if present)
if [[ -f "motifs/tetO_positional_weights_19.tsv" ]]; then
  w_meme=$(awk 'toupper($1)=="W"{print $2; exit}' "$TETO19_MOTIF" 2>/dev/null || true)
  n_wts=$(awk 'NF>=2 && $1 ~ /^[0-9]+$/ {n++} END{print n+0}' motifs/tetO_positional_weights_19.tsv)
  if [[ -n "$w_meme" && "$w_meme" != "0" && "$n_wts" != "0" && "$w_meme" != "$n_wts" ]]; then
    echo "[run_01_05_fimo_pipeline] WARNING: tetO19 motif width (w=$w_meme) != weights rows ($n_wts)." >&2
    echo "[run_01_05_fimo_pipeline]          If you changed weights, regenerate the .meme to match before running FIMO." >&2
  fi
fi

stamp() { date +%Y%m%d_%H%M%S; }

run_one() {
  local label="$1" motif="$2" fasta="$3" outdir="$4" thresh="$5"

  echo "[run_01_05_fimo_pipeline] ==> ${label}" >&2
  echo "  motif:  $motif" >&2
  echo "  fasta:  $fasta" >&2
  echo "  out:    $outdir" >&2
  echo "  thresh: $thresh" >&2

  if [[ -s "$outdir/fimo.tsv" && "$FORCE" != "1" ]]; then
    echo "[run_01_05_fimo_pipeline]   found existing fimo.tsv; skipping (set FORCE=1 to rebuild)" >&2
    return 0
  fi

  if [[ -e "$outdir" && "$FORCE" == "1" ]]; then
    local bak="${outdir}.bak_$(stamp)"
    echo "[run_01_05_fimo_pipeline]   moving existing '$outdir' -> '$bak'" >&2
    mv "$outdir" "$bak"
  fi

  mkdir -p "$outdir"

  # MEME suite: --thresh sets p-value threshold
  fimo --oc "$outdir" --thresh "$thresh" "$motif" "$fasta"

  if [[ ! -s "$outdir/fimo.tsv" && ! -s "$outdir/fimo.txt" ]]; then
    echo "[run_01_05_fimo_pipeline] ERROR: FIMO did not produce fimo.tsv or fimo.txt in: $outdir" >&2
    exit 3
  fi

  # quick stats
  local hits=0
  if [[ -s "$outdir/fimo.tsv" ]]; then
    hits=$(awk 'NR>1{n++} END{print n+0}' "$outdir/fimo.tsv")
  elif [[ -s "$outdir/fimo.txt" ]]; then
    hits=$(grep -vc '^#' "$outdir/fimo.txt" || true)
  fi
  echo "[run_01_05_fimo_pipeline]   done. hits: $hits" >&2
}

if [[ "$RUN_EBOX" == "1" ]]; then
  run_one "EBOX (MYC E-box)" "$EBOX_MOTIF" "$PROMOTERS_FA" "$EBOX_OUT" "$THRESH_EBOX"
fi

if [[ "$RUN_TETO19" == "1" ]]; then
  run_one "TETO19" "$TETO19_MOTIF" "$PROMOTERS_FA" "$TETO19_OUT" "$THRESH_TETO19"
fi

if [[ "$RUN_POSCTRL19" == "1" ]]; then
  run_one "POSCONTROL TETO19 (TRE3Gs_7x)" "$TETO19_MOTIF" "$TRE7_FA" "$POSCTRL19_OUT" "$THRESH_POSCTRL19"
fi

echo "[run_01_05_fimo_pipeline] All done." >&2
