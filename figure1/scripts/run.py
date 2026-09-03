#!/usr/bin/env python3
"""
Production Figure 1: transcriptome-wide density and delta-density summaries.

Current manuscript-aligned default:
  - expression scale: size-factor-normalized RSEM expected counts
  - low-count filter: total expected counts >= 80 across all samples
  - delta-density contrasts: config/contrasts.csv matched numerator/denominator pairs
  - uncertainty: gene-resampling bootstrap bands (not replicate-level CIs)

The script also supports TPM and raw-count sensitivity output via command-line options.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, ConnectionPatch
from scipy.ndimage import gaussian_filter1d


CONDITION_ORDER = [
    "Ctrl",
    "Tam",
    "D",
    "DT",
    "D_L_CPT",
    "DT_L_CPT",
    "D_H_CPT",
    "DT_H_CPT",
]

# Okabe-Ito-like palette, keeping the legacy color intent where possible.
CONDITION_COLORS = {
    "Ctrl": "#000000",
    "Tam": "#E69F00",
    "D": "#56B4E9",
    "DT": "#009E73",
    "D_L_CPT": "#F0E442",
    "DT_L_CPT": "#0072B2",
    "D_H_CPT": "#D55E00",
    "DT_H_CPT": "#CC79A7",
}

TIME_ORDER = ["4h", "24h"]
TIME_LINESTYLES = {"4h": (0, (4, 2)), "24h": "solid"}


@dataclass(frozen=True)
class ExpressionData:
    log_expr: pd.DataFrame          # rows genes, columns condition keys like "4h|DT"
    values_by_group: Dict[Tuple[str, str], np.ndarray]
    gene_ids: List[str]
    sample_size_factors: Optional[pd.Series]
    scale_name: str
    axis_label: str
    matrix_path: str
    n_genes_before_filter: int
    n_genes_after_filter: int


def _repo_path(root: Path, rel_or_abs: str) -> Path:
    p = Path(rel_or_abs)
    return p if p.is_absolute() else root / p


def _first_existing(root: Path, candidates: Sequence[str]) -> Path:
    for c in candidates:
        p = _repo_path(root, c)
        if p.exists():
            return p
    raise FileNotFoundError("None of these candidate files exists: " + ", ".join(candidates))


def _read_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing matrix: {path}")
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    if df.shape[1] < 2:
        raise ValueError(f"Matrix has fewer than two columns: {path}")
    first = df.columns[0]
    if first.lower() not in {"gene_id", "gene", "genes", "name", "id"}:
        # Still treat first column as gene IDs; legacy files sometimes use nonstandard names.
        pass
    df = df.rename(columns={first: "gene_id"})
    if df["gene_id"].duplicated().any():
        dup = df.loc[df["gene_id"].duplicated(), "gene_id"].head().tolist()
        raise ValueError(f"Duplicated gene IDs in {path}: {dup}")
    df = df.set_index("gene_id")
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.isna().any().any():
        n_bad = int(df.isna().sum().sum())
        raise ValueError(f"Matrix contains {n_bad} non-numeric/missing values after parsing: {path}")
    if (df < 0).any().any():
        raise ValueError(f"Matrix contains negative values: {path}")
    return df.astype(float)


def _condition_from_group_label(group_label: str) -> str:
    label = str(group_label)
    label = re.sub(r"^(4|24)_?", "", label)
    if label.lower() == "ctrl":
        return "Ctrl"
    return label


def _read_sample_sheet(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing sample sheet: {path}")
    ss = pd.read_csv(path)
    required = {"sample_id", "group_label", "time_h"}
    missing = sorted(required - set(ss.columns))
    if missing:
        raise ValueError(f"Sample sheet missing columns {missing}: {path}")
    ss = ss.copy()
    ss["condition"] = ss["group_label"].map(_condition_from_group_label)
    ss["timepoint"] = ss["time_h"].map(lambda x: f"{int(x)}h")
    return ss


def _size_factors_median_ratio(counts: pd.DataFrame) -> pd.Series:
    """DESeq2-style median-ratio size factors using genes positive in all samples."""
    arr = counts.to_numpy(dtype=float)
    valid = np.all(arr > 0, axis=1)
    if valid.sum() < 100:
        raise ValueError(
            f"Only {valid.sum()} genes are positive in all samples; cannot estimate stable size factors."
        )
    gm = np.exp(np.mean(np.log(arr[valid, :]), axis=1))
    ratios = arr[valid, :] / gm[:, None]
    sf = np.median(ratios, axis=0)
    if not np.all(np.isfinite(sf)) or np.any(sf <= 0):
        raise ValueError("Invalid size factors computed from count matrix.")
    # Center to geometric mean 1 for interpretability; DESeq2 factors are relative, so this is harmless.
    sf = sf / np.exp(np.mean(np.log(sf)))
    return pd.Series(sf, index=counts.columns, name="median_ratio_size_factor")


def _load_expression(
    *,
    root: Path,
    expression_source: str,
    counts_path: Path,
    tpm_path: Path,
    sample_sheet_path: Path,
    min_total_count: Optional[float],
) -> Tuple[ExpressionData, pd.DataFrame]:
    ss = _read_sample_sheet(sample_sheet_path)

    if expression_source in {"norm_counts", "raw_counts"}:
        mat = _read_matrix(counts_path)
        matrix_path = str(counts_path)
    elif expression_source == "tpm":
        mat = _read_matrix(tpm_path)
        matrix_path = str(tpm_path)
    else:
        raise ValueError(f"Unsupported expression source: {expression_source}")

    sample_ids = ss["sample_id"].tolist()
    missing_samples = [s for s in sample_ids if s not in mat.columns]
    if missing_samples:
        raise ValueError(f"Expression matrix missing {len(missing_samples)} sample columns, e.g. {missing_samples[:5]}")
    mat = mat.loc[:, sample_ids]

    n_before = int(mat.shape[0])
    if min_total_count is not None and min_total_count > 0:
        # Filter using counts when available. If expression_source is TPM, still use counts if present; otherwise
        # filtering by TPM totals is avoided because it has a different meaning.
        if counts_path.exists():
            counts_for_filter = _read_matrix(counts_path).loc[:, sample_ids]
            keep = counts_for_filter.sum(axis=1) >= float(min_total_count)
            mat = mat.loc[mat.index.intersection(keep[keep].index), :]
        elif expression_source != "tpm":
            keep = mat.sum(axis=1) >= float(min_total_count)
            mat = mat.loc[keep, :]
        else:
            raise FileNotFoundError(
                "--min-total-count was requested for TPM input, but no count matrix was available for filtering."
            )
    n_after = int(mat.shape[0])

    size_factors: Optional[pd.Series] = None
    if expression_source == "norm_counts":
        # Always estimate size factors from the filtered count matrix used for the figure.
        size_factors = _size_factors_median_ratio(mat)
        expr = mat.divide(size_factors, axis=1)
        scale_name = "median-ratio normalized RSEM expected counts"
        axis_label = "log2(normalized count + 1)"
    elif expression_source == "raw_counts":
        expr = mat
        scale_name = "unnormalized RSEM expected counts"
        axis_label = "log2(RSEM expected count + 1)"
    else:
        expr = mat
        size_factors = None
        scale_name = "TPM"
        axis_label = "log2(TPM + 1)"

    values_by_group: Dict[Tuple[str, str], np.ndarray] = {}
    log_expr_cols = {}
    for timepoint in TIME_ORDER:
        for condition in CONDITION_ORDER:
            sample_ids_group = ss.loc[
                (ss["timepoint"] == timepoint) & (ss["condition"] == condition), "sample_id"
            ].tolist()
            if not sample_ids_group:
                continue
            mean_expr = expr.loc[:, sample_ids_group].mean(axis=1)
            log_values = np.log2(mean_expr.to_numpy(dtype=float) + 1.0)
            key = f"{timepoint}|{condition}"
            log_expr_cols[key] = log_values
            values_by_group[(timepoint, condition)] = log_values

    log_expr = pd.DataFrame(log_expr_cols, index=mat.index)
    return (
        ExpressionData(
            log_expr=log_expr,
            values_by_group=values_by_group,
            gene_ids=mat.index.astype(str).tolist(),
            sample_size_factors=size_factors,
            scale_name=scale_name,
            axis_label=axis_label,
            matrix_path=matrix_path,
            n_genes_before_filter=n_before,
            n_genes_after_filter=n_after,
        ),
        ss,
    )


def _grid_from_values(values: Iterable[np.ndarray], grid_n: int, x_min: Optional[float], x_max: Optional[float]) -> np.ndarray:
    all_values = np.concatenate([np.asarray(v, dtype=float) for v in values])
    all_values = all_values[np.isfinite(all_values)]
    if all_values.size == 0:
        raise ValueError("No finite expression values available for density estimation.")
    lo = math.floor(float(np.min(all_values))) if x_min is None else float(x_min)
    hi = math.ceil(float(np.max(all_values))) if x_max is None else float(x_max)
    if hi <= lo:
        hi = lo + 1.0
    return np.linspace(lo, hi, int(grid_n))


def _density(values: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.full(grid.shape, np.nan)
    dx = float(grid[1] - grid[0])
    edges = np.concatenate(([grid[0] - dx / 2.0], (grid[:-1] + grid[1:]) / 2.0, [grid[-1] + dx / 2.0]))
    hist, _ = np.histogram(values, bins=edges)
    den = hist.astype(float) / max(values.size * dx, 1.0)
    sigma = max(float(bandwidth) / dx, 0.01)
    den = gaussian_filter1d(den, sigma=sigma, mode="constant")
    area = np.trapezoid(den, grid)
    if area > 0 and np.isfinite(area):
        den = den / area
    return den


def _bootstrap_density_delta_ci(
    num_values: np.ndarray,
    den_values: np.ndarray,
    grid: np.ndarray,
    bandwidth: float,
    n_boot: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Approximate gene-resampling CI using binned multinomial bootstrap + Gaussian smoothing."""
    if n_boot <= 0:
        return np.full(grid.shape, np.nan), np.full(grid.shape, np.nan)

    num_values = np.asarray(num_values, dtype=float)
    den_values = np.asarray(den_values, dtype=float)
    num_values = num_values[np.isfinite(num_values)]
    den_values = den_values[np.isfinite(den_values)]
    if min(num_values.size, den_values.size) < 10:
        return np.full(grid.shape, np.nan), np.full(grid.shape, np.nan)

    dx = float(grid[1] - grid[0])
    edges = np.concatenate(([grid[0] - dx / 2.0], (grid[:-1] + grid[1:]) / 2.0, [grid[-1] + dx / 2.0]))

    h_num, _ = np.histogram(num_values, bins=edges)
    h_den, _ = np.histogram(den_values, bins=edges)
    if h_num.sum() == 0 or h_den.sum() == 0:
        return np.full(grid.shape, np.nan), np.full(grid.shape, np.nan)

    p_num = h_num.astype(float) / h_num.sum()
    p_den = h_den.astype(float) / h_den.sum()
    sigma = max(float(bandwidth) / dx, 0.01)

    # Work in blocks to avoid holding very large matrices when n_boot=5000.
    block = min(500, n_boot)
    out = []
    done = 0
    while done < n_boot:
        b = min(block, n_boot - done)
        c_num = rng.multinomial(num_values.size, p_num, size=b).astype(float) / (num_values.size * dx)
        c_den = rng.multinomial(den_values.size, p_den, size=b).astype(float) / (den_values.size * dx)
        d_num = gaussian_filter1d(c_num, sigma=sigma, axis=1, mode="constant")
        d_den = gaussian_filter1d(c_den, sigma=sigma, axis=1, mode="constant")
        # Normalize each bootstrap density to unit integral.
        area_num = np.trapezoid(d_num, grid, axis=1)
        area_den = np.trapezoid(d_den, grid, axis=1)
        d_num = d_num / area_num[:, None]
        d_den = d_den / area_den[:, None]
        out.append(d_num - d_den)
        done += b
    mat = np.vstack(out)
    return np.nanpercentile(mat, 2.5, axis=0), np.nanpercentile(mat, 97.5, axis=0)


def _build_densities(expr: ExpressionData, grid: np.ndarray, bandwidth: float) -> pd.DataFrame:
    rows = []
    for (timepoint, condition), vals in expr.values_by_group.items():
        y = _density(vals, grid, bandwidth)
        rows.append(pd.DataFrame({
            "timepoint": timepoint,
            "condition": condition,
            "x": grid,
            "density": y,
        }))
    return pd.concat(rows, ignore_index=True)


def _build_delta_curves(
    dens: pd.DataFrame,
    contrasts: pd.DataFrame,
    delta_mode: str,
) -> pd.DataFrame:
    rows = []
    dens_key = {(timepoint, condition): sub["density"].to_numpy() for (timepoint, condition), sub in dens.groupby(["timepoint", "condition"])}
    grid = dens["x"].drop_duplicates().to_numpy()

    if delta_mode == "contrasts":
        for _, r in contrasts.iterrows():
            timepoint = f"{int(r['time_h'])}h"
            num = str(r["numerator"])
            den = str(r["denominator"])
            if (timepoint, num) not in dens_key or (timepoint, den) not in dens_key:
                continue
            rows.append(pd.DataFrame({
                "contrast_id": r["contrast_id"],
                "timepoint": timepoint,
                "numerator": num,
                "denominator": den,
                "x": grid,
                "delta_density": dens_key[(timepoint, num)] - dens_key[(timepoint, den)],
            }))
    elif delta_mode == "ctrl":
        for timepoint in TIME_ORDER:
            if (timepoint, "Ctrl") not in dens_key:
                continue
            for cond in CONDITION_ORDER:
                if cond == "Ctrl" or (timepoint, cond) not in dens_key:
                    continue
                rows.append(pd.DataFrame({
                    "contrast_id": f"{cond}_vs_Ctrl_{timepoint}",
                    "timepoint": timepoint,
                    "numerator": cond,
                    "denominator": "Ctrl",
                    "x": grid,
                    "delta_density": dens_key[(timepoint, cond)] - dens_key[(timepoint, "Ctrl")],
                }))
    else:
        raise ValueError(f"Unsupported delta mode: {delta_mode}")
    if not rows:
        raise ValueError("No delta curves could be constructed. Check sample labels and contrasts.")
    return pd.concat(rows, ignore_index=True)


def _build_bootstrap_bands(
    expr: ExpressionData,
    delta: pd.DataFrame,
    grid: np.ndarray,
    bandwidth: float,
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    pairs = delta[["contrast_id", "timepoint", "numerator", "denominator"]].drop_duplicates()
    for _, r in pairs.iterrows():
        key_num = (r["timepoint"], r["numerator"])
        key_den = (r["timepoint"], r["denominator"])
        if key_num not in expr.values_by_group or key_den not in expr.values_by_group:
            continue
        ymin, ymax = _bootstrap_density_delta_ci(
            expr.values_by_group[key_num],
            expr.values_by_group[key_den],
            grid,
            bandwidth,
            n_boot,
            rng,
        )
        rows.append(pd.DataFrame({
            "contrast_id": r["contrast_id"],
            "timepoint": r["timepoint"],
            "numerator": r["numerator"],
            "denominator": r["denominator"],
            "x": grid,
            "ci_low": ymin,
            "ci_high": ymax,
        }))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _condition_sort_key(c: str) -> int:
    try:
        return CONDITION_ORDER.index(c)
    except ValueError:
        return len(CONDITION_ORDER)


def _write_source_data(
    *,
    expr: ExpressionData,
    sample_sheet: pd.DataFrame,
    contrasts: pd.DataFrame,
    dens: pd.DataFrame,
    delta: pd.DataFrame,
    bands: pd.DataFrame,
    source_dir: Path,
    args: argparse.Namespace,
) -> None:
    source_dir.mkdir(parents=True, exist_ok=True)

    # Long condition-level expression table.
    long_expr = expr.log_expr.copy()
    long_expr.insert(0, "gene_id", expr.log_expr.index)
    long_expr = long_expr.melt(id_vars="gene_id", var_name="time_condition", value_name="log_expression")
    long_expr[["timepoint", "condition"]] = long_expr["time_condition"].str.split("|", regex=False, expand=True)
    long_expr = long_expr.drop(columns=["time_condition"])
    long_expr.to_csv(source_dir / "figure1_condition_mean_log_expression.tsv", sep="\t", index=False)

    dens.to_csv(source_dir / "figure1_raw_density_curves.tsv", sep="\t", index=False)
    delta.to_csv(source_dir / "figure1_delta_density_curves.tsv", sep="\t", index=False)
    if not bands.empty:
        bands.to_csv(source_dir / "figure1_delta_density_gene_bootstrap.tsv", sep="\t", index=False)

    if expr.sample_size_factors is not None:
        expr.sample_size_factors.rename_axis("sample_id").reset_index().to_csv(
            source_dir / "figure1_median_ratio_size_factors.tsv", sep="\t", index=False
        )

    contrasts.to_csv(source_dir / "figure1_contrast_table_used.tsv", sep="\t", index=False)

    manifest = {
        "figure": "Figure 1",
        "script": "figure1/scripts/run.py",
        "expression_source_argument": args.expression_source,
        "expression_scale_name": expr.scale_name,
        "axis_label": expr.axis_label,
        "matrix_path": expr.matrix_path,
        "sample_sheet": str(args.sample_sheet),
        "contrasts": str(args.contrasts),
        "delta_mode": args.delta_mode,
        "low_count_filter_total_counts": args.min_total_count,
        "n_genes_before_filter": expr.n_genes_before_filter,
        "n_genes_after_filter": expr.n_genes_after_filter,
        "density_bandwidth_log2_units": args.bandwidth,
        "grid_n": args.grid_n,
        "bootstrap_resamples": args.bootstrap_resamples,
        "bootstrap_interpretation": "gene-resampling density bands; not biological-replicate confidence intervals",
        "panel_c_linetype": {"4h": "dashed", "24h": "solid"},
        "condition_order": CONDITION_ORDER,
    }
    with open(source_dir / "figure1_manifest.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    with open(source_dir / "figure1_manifest.txt", "w", encoding="utf-8") as fh:
        for k, v in manifest.items():
            fh.write(f"{k}: {v}\n")


def _plot_figure(
    *,
    expr: ExpressionData,
    dens: pd.DataFrame,
    delta: pd.DataFrame,
    bands: pd.DataFrame,
    out_dir: Path,
    stem: str,
    dpi: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "legend.title_fontsize": 9.5,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig = plt.figure(figsize=(8.6, 8.2), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=3,
        ncols=3,
        width_ratios=[1.0, 1.0, 0.42],
        height_ratios=[1.0, 1.0, 0.95],
        left=0.10,
        right=0.98,
        bottom=0.08,
        top=0.96,
        wspace=0.32,
        hspace=0.44,
    )
    ax_a = {"4h": fig.add_subplot(gs[0, 0]), "24h": fig.add_subplot(gs[0, 1])}
    ax_b = {"4h": fig.add_subplot(gs[1, 0]), "24h": fig.add_subplot(gs[1, 1])}
    ax_c_zoom = fig.add_subplot(gs[2, 0])
    ax_c_full = fig.add_subplot(gs[2, 1])
    ax_leg = fig.add_subplot(gs[:, 2])
    ax_leg.axis("off")

    # Panel A and B.
    for timepoint, ax in ax_a.items():
        ax.axhline(0, color="0.35", lw=0.8)
        dd = delta[delta["timepoint"] == timepoint]
        bb = bands[bands["timepoint"] == timepoint] if not bands.empty else pd.DataFrame()
        for _, rr in dd[["contrast_id", "numerator", "denominator"]].drop_duplicates().iterrows():
            cid = rr["contrast_id"]
            num = rr["numerator"]
            sub = dd[dd["contrast_id"] == cid]
            color = CONDITION_COLORS.get(num, "0.4")
            if not bb.empty:
                bsub = bb[bb["contrast_id"] == cid]
                if not bsub.empty and bsub["ci_low"].notna().any():
                    ax.fill_between(bsub["x"], bsub["ci_low"], bsub["ci_high"], color=color, alpha=0.06, lw=0)
            ax.plot(sub["x"], sub["delta_density"], color=color, lw=1.0)
        ax.set_title(timepoint, fontweight="bold", pad=4)
        ax.set_xlim(float(delta["x"].min()), float(delta["x"].max()))
        ax.set_ylim(-0.10, 0.06)
        ax.grid(True, color="0.90", lw=0.8)
        ax.set_xlabel(expr.axis_label)
    ax_a["4h"].set_ylabel("Δ density")
    ax_a["24h"].set_ylabel("")

    for timepoint, ax in ax_b.items():
        ax.axhline(0, color="0.35", lw=0.8)
        dd = delta[delta["timepoint"] == timepoint]
        bb = bands[bands["timepoint"] == timepoint] if not bands.empty else pd.DataFrame()
        for _, rr in dd[["contrast_id", "numerator", "denominator"]].drop_duplicates().iterrows():
            cid = rr["contrast_id"]
            num = rr["numerator"]
            sub = dd[dd["contrast_id"] == cid]
            color = CONDITION_COLORS.get(num, "0.4")
            if not bb.empty:
                bsub = bb[bb["contrast_id"] == cid]
                if not bsub.empty and bsub["ci_low"].notna().any():
                    ax.fill_between(bsub["x"], bsub["ci_low"], bsub["ci_high"], color=color, alpha=0.08, lw=0)
            ax.plot(sub["x"], sub["delta_density"], color=color, lw=1.05)
        ax.set_title(timepoint, fontweight="bold", pad=4)
        ax.set_xlim(5.0, 15.0)
        ax.set_ylim(-0.035, 0.035)
        ax.grid(True, color="0.90", lw=0.8)
        ax.set_xlabel(expr.axis_label)
    ax_b["4h"].set_ylabel("Δ density")
    ax_b["24h"].set_ylabel("")

    # Panel C raw density overlays.
    for timepoint in TIME_ORDER:
        for condition in CONDITION_ORDER:
            sub = dens[(dens["timepoint"] == timepoint) & (dens["condition"] == condition)]
            if sub.empty:
                continue
            color = CONDITION_COLORS.get(condition, "0.4")
            ls = TIME_LINESTYLES[timepoint]
            ax_c_zoom.plot(sub["x"], sub["density"], color=color, lw=0.95, ls=ls, alpha=0.98)
            ax_c_full.plot(sub["x"], sub["density"], color=color, lw=0.95, ls=ls, alpha=0.98)

    ax_c_zoom.set_xlim(5.0, 15.0)
    # Adaptive zoom y-limit but keep close to the legacy window.
    max_zoom = dens[(dens["x"] >= 5.0) & (dens["x"] <= 15.0)]["density"].max()
    ax_c_zoom.set_ylim(0, max(0.09, max_zoom * 1.12))
    ax_c_full.set_xlim(float(dens["x"].min()), float(dens["x"].max()))
    ax_c_full.set_ylim(0, max(0.20, dens["density"].max() * 1.12))
    for ax in (ax_c_zoom, ax_c_full):
        ax.grid(True, color="0.90", lw=0.8)
        ax.set_xlabel(expr.axis_label)
    ax_c_zoom.set_ylabel("Density")
    ax_c_full.set_ylabel("")

    # Zoom region indicator on full panel plus connector lines.
    y0, y1 = ax_c_zoom.get_ylim()
    rect_h = min(y1, ax_c_full.get_ylim()[1] * 0.18)
    rect = Rectangle((5.0, 0), 10.0, rect_h, facecolor="0.88", edgecolor="0.2", lw=0.8, alpha=0.5)
    ax_c_full.add_patch(rect)
    con1 = ConnectionPatch(
        xyA=(15.0, y1), coordsA=ax_c_zoom.transData,
        xyB=(5.0, rect_h), coordsB=ax_c_full.transData,
        color="0.25", lw=0.8,
    )
    con2 = ConnectionPatch(
        xyA=(15.0, y0), coordsA=ax_c_zoom.transData,
        xyB=(5.0, 0), coordsB=ax_c_full.transData,
        color="0.25", lw=0.8,
    )
    fig.add_artist(con1)
    fig.add_artist(con2)

    # Panel labels: uppercase, matching figure legend style.
    fig.text(0.015, 0.955, "A", fontsize=17, fontweight="bold", va="top", ha="left")
    fig.text(0.015, 0.620, "B", fontsize=17, fontweight="bold", va="top", ha="left")
    fig.text(0.015, 0.300, "C", fontsize=17, fontweight="bold", va="top", ha="left")

    condition_handles = [
        Line2D([0], [0], color=CONDITION_COLORS[c], lw=2.0, label=c)
        for c in CONDITION_ORDER
        if c in set(dens["condition"])
    ]
    time_handles = [
        Line2D([0], [0], color="0.25", lw=1.8, ls=TIME_LINESTYLES[t], label=t)
        for t in TIME_ORDER
    ]
    leg1 = ax_leg.legend(
        handles=condition_handles,
        title="Condition",
        loc="upper left",
        frameon=False,
        borderaxespad=0,
        handlelength=2.3,
    )
    ax_leg.add_artist(leg1)
    ax_leg.legend(
        handles=time_handles,
        title="Panel C time",
        loc="upper left",
        bbox_to_anchor=(0.0, 0.56),
        frameon=False,
        borderaxespad=0,
        handlelength=2.3,
    )

    for ext in ["pdf", "svg", "png"]:
        out = out_dir / f"{stem}.{ext}"
        if ext == "png":
            fig.savefig(out, dpi=dpi)
        else:
            fig.savefig(out)
    plt.close(fig)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root; default current working directory.")
    parser.add_argument(
        "--expression-source",
        choices=["norm_counts", "tpm", "raw_counts"],
        default="norm_counts",
        help="Expression scale for Figure 1. Default is normalized RSEM expected counts.",
    )
    parser.add_argument("--counts", default=None, help="Count matrix path. Defaults to processed then raw candidates.")
    parser.add_argument("--tpm", default=None, help="TPM matrix path. Defaults to processed candidates.")
    parser.add_argument("--sample-sheet", default="config/sample_sheet.csv")
    parser.add_argument("--contrasts", default="config/contrasts.csv")
    parser.add_argument(
        "--delta-mode",
        choices=["contrasts", "ctrl"],
        default="contrasts",
        help="Default uses config/contrasts.csv matched numerator/denominator pairs. 'ctrl' reproduces legacy Ctrl-subtraction.",
    )
    parser.add_argument("--min-total-count", type=float, default=80.0, help="Low-count filter based on total expected counts. Use 0 to disable.")
    parser.add_argument("--grid-n", type=int, default=1024)
    parser.add_argument("--bandwidth", type=float, default=0.25, help="Gaussian KDE bandwidth in log2-expression units.")
    parser.add_argument("--bootstrap-resamples", type=int, default=5000)
    parser.add_argument("--quick", action="store_true", help="Quick QA run with 200 bootstrap resamples.")
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--out-dir", default="figure1/outputs")
    parser.add_argument("--source-dir", default="figure1/source_data")
    parser.add_argument("--output-stem", default="figure1_main")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float, default=None)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    root = Path(args.repo_root).resolve()
    counts_path = _repo_path(root, args.counts) if args.counts else _first_existing(root, [
        "data/processed/raw_counts_rsemgenes.tsv",
        "data/processed/RawCountFile_rsemgenes.tsv",
        "data/raw/RawCountFile_rsemgenes.txt",
        "data/raw/rsem/RawCountFile_rsemgenes.txt",
    ])
    tpm_path = _repo_path(root, args.tpm) if args.tpm else _first_existing(root, [
        "data/processed/TPMCountFile_rsemgenes.csv",
        "data/processed/TPMCountFile_rsemgenes.tsv",
        "data/processed/tpm_rsemgenes.tsv",
        "data/raw/TPMCountFile_rsemgenes.csv",
        "data/raw/TPMCountFile_rsemgenes.txt",
        "data/raw/rsem/TPMCountFile_rsemgenes.txt",
    ])
    sample_sheet_path = _repo_path(root, args.sample_sheet)
    contrasts_path = _repo_path(root, args.contrasts)
    out_dir = _repo_path(root, args.out_dir)
    source_dir = _repo_path(root, args.source_dir)

    if args.quick:
        args.bootstrap_resamples = min(args.bootstrap_resamples, 200)

    min_total_count = None if args.min_total_count <= 0 else float(args.min_total_count)
    expr, sample_sheet = _load_expression(
        root=root,
        expression_source=args.expression_source,
        counts_path=counts_path,
        tpm_path=tpm_path,
        sample_sheet_path=sample_sheet_path,
        min_total_count=min_total_count,
    )
    contrasts = pd.read_csv(contrasts_path)
    dens = _build_densities(expr, _grid_from_values(expr.values_by_group.values(), args.grid_n, args.x_min, args.x_max), args.bandwidth)
    grid = dens["x"].drop_duplicates().to_numpy()
    delta = _build_delta_curves(dens, contrasts, args.delta_mode)
    bands = _build_bootstrap_bands(expr, delta, grid, args.bandwidth, args.bootstrap_resamples, args.seed)

    _write_source_data(
        expr=expr,
        sample_sheet=sample_sheet,
        contrasts=contrasts,
        dens=dens,
        delta=delta,
        bands=bands,
        source_dir=source_dir,
        args=args,
    )
    _plot_figure(expr=expr, dens=dens, delta=delta, bands=bands, out_dir=out_dir, stem=args.output_stem, dpi=args.dpi)

    done = out_dir / f"{args.output_stem}.done"
    done.write_text(
        f"Figure 1 complete\n"
        f"expression_source={args.expression_source}\n"
        f"expression_scale={expr.scale_name}\n"
        f"delta_mode={args.delta_mode}\n"
        f"n_genes_after_filter={expr.n_genes_after_filter}\n"
        f"bootstrap_resamples={args.bootstrap_resamples}\n",
        encoding="utf-8",
    )
    print(done)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
