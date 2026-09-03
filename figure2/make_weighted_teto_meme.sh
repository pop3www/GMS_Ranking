#!/usr/bin/env bash
set -euo pipefail
out=${1:-motifs/tetO_custom.meme}
seq=${2:-TCTATCATTGATAGG}
weights_csv=${3:-"0.70,0.92,0.98,0.95,0.95,0.90,0.80,0.25,0.80,0.90,0.95,0.95,0.98,0.92,0.70"}
mkdir -p "$(dirname "$out")"
L=${#seq}
IFS=',' read -r -a W <<< "$weights_csv"
if [ "${#W[@]}" -ne "$L" ]; then
  echo "[ERR] weights length (${#W[@]}) != motif length ($L)"; exit 1
fi

{
  echo "MEME version 4"
  echo
  echo "ALPHABET= ACGT"
  echo "strands: + -"
  echo
  echo "Background letter frequencies:"
  echo "A 0.25 C 0.25 G 0.25 T 0.25"
  echo
  echo "MOTIF tetO_weighted_${L}"
  echo "letter-probability matrix: alength= 4 w= ${L} nsites= 7 E= 0"
  for ((i=0;i<L;i++)); do
    base=${seq:i:1}
    w=${W[$i]}
    other=$(awk -v w="$w" 'BEGIN{printf "%.6f",(1.0-w)/3.0}')
    pa=$other; pc=$other; pg=$other; pt=$other
    case $base in
      A) pa=$w;;
      C) pc=$w;;
      G) pg=$w;;
      T) pt=$w;;
      *) echo "[ERR] non-ACGT base at pos $((i+1)): $base" >&2; exit 1;;
    esac
    printf "%.6f %.6f %.6f %.6f\n" "$pa" "$pc" "$pg" "$pt"
  done
} > "$out"
echo "[OK] wrote $out"