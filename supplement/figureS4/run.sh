#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
exec Rscript supplement/figureS4/scripts/run_figS4_gene_length_sensitivity.R "$@"
