#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
MODE="${1:-list}"
commands=(
  "python3 core/scripts/validate_expression_inputs.py"
  "bash figure1/run.sh"
  "bash figure2/run.sh"
  "bash supplement/figureS1/run.sh"
  "python3 figure3/differential_expression/scripts/run.py"
  "python3 figure3/scalar_metrics/scripts/run.py"
  "bash figure3/final_assembly/run.sh"
  "bash figure4/run.sh"
  "bash figure5/run.sh"
  "bash figure6_7/run.sh"
  "bash supplement/figureS4/run.sh"
  "bash supplement/figureS5/run.sh"
)
case "$MODE" in
  list) printf '%s\n' "${commands[@]}" ;;
  execute)
    echo "Running canonical workflows; substantial time and memory may be required."
    for cmd in "${commands[@]}"; do
      echo "+ $cmd"
      bash -lc "$cmd"
    done
    ;;
  *) echo "Usage: $0 [list|execute]" >&2; exit 2 ;;
esac
