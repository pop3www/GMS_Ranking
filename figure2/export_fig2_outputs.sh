#!/usr/bin/env bash
set -euo pipefail

RES=${1:-}
if [[ -z "${RES}" ]]; then
  echo "Usage: export_fig2_outputs.sh <results_root>"
  echo "Example: export_fig2_outputs.sh results_fig2_teto19_weighted_FINAL"
  exit 1
fi

FIGDIR="${RES}/FIGURE2_FINAL"
mkdir -p "${FIGDIR}"

copy_if_exists () {
  local src="$1"; local dst="$2"
  if [[ -f "$src" ]]; then
    cp -f "$src" "$dst"
    echo "[copy] $dst"
  else
    echo "[missing] $src" >&2
  fi
}

# Fig2A: delta-density controls, expected at results root.
for ext in svg png pdf; do
  copy_if_exists "${RES}/Fig2A_delta_density_controls.${ext}" "${FIGDIR}/Fig2A_delta_density_controls.${ext}"
done

# Fig2B/C: priming comparison.
PRIM="${RES}/priming_compare_DTvsD_vs_TamvsCtrl"
for ext in svg png pdf; do
  copy_if_exists "${PRIM}/Fig2B_priming_dependent_component.${ext}" "${FIGDIR}/Fig2B_priming_dependent_component.${ext}"
  copy_if_exists "${PRIM}/Fig2C_DTvsD_vs_TamvsCtrl_scatter.${ext}" "${FIGDIR}/Fig2C_DTvsD_vs_TamvsCtrl_scatter.${ext}"
done

# Fig2D-F: motif scan panels from DT_vs_D summary.
SUM="${RES}/fimo_summary_4h_DT_vs_D"
for ext in svg png pdf; do
  copy_if_exists "${SUM}/Fig2D_tetO_logo_weights.${ext}" "${FIGDIR}/Fig2D_tetO_logo_weights.${ext}"
  copy_if_exists "${SUM}/Fig2E_hit_density_by_band.${ext}" "${FIGDIR}/Fig2E_hit_density_by_band.${ext}"
  copy_if_exists "${SUM}/Fig2F_enrichment_OR.${ext}" "${FIGDIR}/Fig2F_enrichment_OR.${ext}"
done

# Source/supplemental auxiliary figures (optional but useful for traceability).
AUX="${RES}/FIGURE2_AUXILIARY"
mkdir -p "${AUX}"
for f in \
  "${SUM}/Supp_tetO_logo_with_weights.svg" "${SUM}/Supp_tetO_logo_with_weights.png" \
  "${SUM}/Supp_tetO_logo.svg" "${SUM}/Supp_tetO_logo.png" \
  "${SUM}/Supp_tetO_weights.svg" "${SUM}/Supp_tetO_weights.png" \
  "${SUM}/Supp_tetO_IC.svg" "${SUM}/Supp_tetO_IC.png" \
  "${SUM}/Supp_hit_density_by_band.svg" "${SUM}/Supp_hit_density_by_band.png" \
  "${SUM}/Supp_gene_fraction_by_band.svg" "${SUM}/Supp_gene_fraction_by_band.png" \
  "${SUM}/Supp_enrichment_OR.svg" "${SUM}/Supp_enrichment_OR.png" \
  "${SUM}/Supp_enrichment_OR_up_vs_down.svg" "${SUM}/Supp_enrichment_OR_up_vs_down.png" \
  "${SUM}/Supp_enrichment_OR_by_band.svg" "${SUM}/Supp_enrichment_OR_by_band.png"; do
  [[ -f "$f" ]] && cp -f "$f" "${AUX}/$(basename "$f")"
done

# Manifest
{
  echo "Figure 2 final export for: ${RES}"
  date
  echo
  echo "Expected final panels:"
  find "${FIGDIR}" -maxdepth 1 -type f -exec basename {} \; | sort | sed "s/^/- /"
  echo
  echo "Auxiliary/source panels:"
  find "${AUX}" -maxdepth 1 -type f -exec basename {} \; | sort | sed "s/^/- /"
} > "${RES}/FIGURE2_FINAL_MANIFEST.txt"



if compgen -G "${FIGDIR}/Figure2_combined.*" > /dev/null; then
  echo >> "${RES}/FIGURE2_FINAL_MANIFEST.txt"
  echo "Combined figure outputs:" >> "${RES}/FIGURE2_FINAL_MANIFEST.txt"
  for f in "${FIGDIR}"/Figure2_combined.*; do [[ -f "$f" ]] && echo "- $(basename "$f")"; done >> "${RES}/FIGURE2_FINAL_MANIFEST.txt"
fi

echo "[OK] Final Figure 2 panels in: ${FIGDIR}"
ls -1 "${FIGDIR}" | sed 's/^/  - /'
