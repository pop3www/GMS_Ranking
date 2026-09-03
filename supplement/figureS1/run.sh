#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TRACE="$ROOT/supplement/figureS1/source_data/GEO_TRACE.tsv"
SOURCE_NOTE="$ROOT/supplement/figureS1/source_data/SOURCE_TRACE.txt"
DISPLAY="$ROOT/final_outputs/FigureS1.pdf"

[[ -s "$TRACE" ]] || { echo "ERROR: missing $TRACE" >&2; exit 1; }
[[ -s "$DISPLAY" ]] || { echo "ERROR: missing $DISPLAY" >&2; exit 1; }
grep -q 'GSE318271' "$TRACE" || { echo "ERROR: GSE318271 absent from GEO trace" >&2; exit 1; }
grep -q 'GSE318584' "$TRACE" || { echo "ERROR: GSE318584 absent from GEO trace" >&2; exit 1; }

echo "Figure S1 provenance validation passed."
echo "  ChIP-seq source: GSE318271"
echo "  RNA-seq rank source: GSE318584"
echo "  Archived display: final_outputs/FigureS1.pdf"
if [[ -s "$SOURCE_NOTE" ]]; then
  echo "  Additional provenance: supplement/figureS1/source_data/SOURCE_TRACE.txt"
fi
echo "Full reconstruction of input-normalized promoter profiles begins from the GEO records."
