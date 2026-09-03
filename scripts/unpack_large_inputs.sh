#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)

unpack_one() {
  local gz=$1
  local out=${gz%.gz}
  if [ -f "$out" ]; then
    printf 'Already present: %s\n' "${out#$ROOT/}"
  elif [ -f "$gz" ]; then
    gzip -dk "$gz"
    printf 'Unpacked: %s\n' "${out#$ROOT/}"
  else
    printf 'Missing compressed input: %s\n' "${gz#$ROOT/}" >&2
    return 1
  fi
}

unpack_one "$ROOT/figure2/data/promoters_2kb.fa.gz"
unpack_one "$ROOT/figure3/differential_expression/source_data/figure3_unified_de_results_used.tsv.gz"
