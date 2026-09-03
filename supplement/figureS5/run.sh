#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$ROOT"; exec python3 figure3/differential_expression/scripts/export_figureS5.py --jaccard figure3/differential_expression/source_data/figure3_panelC_jaccard.tsv --spearman figure3/differential_expression/source_data/figure3_panelD_spearman.tsv --output-dir final_outputs --stem FigureS5
