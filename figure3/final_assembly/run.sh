#!/usr/bin/env bash
set -Eeuo pipefail
cd "$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
Rscript 'figure3/final_assembly/90_composite_figures/figure_3_4_merged/scripts/build_figure_3_4_merged.R'
