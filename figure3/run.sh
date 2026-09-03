#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; python3 figure3/differential_expression/scripts/run.py "$@"; python3 figure3/scalar_metrics/scripts/run.py "$@"; exec bash figure3/final_assembly/run.sh
