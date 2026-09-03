#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
Rscript 90_composite_figures/figure_3_4_merged/scripts/build_figure_3_4_merged.R "$@"
