from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

import matplotlib

matplotlib.use("Agg")
# Keep text editable in vector outputs.  Matplotlib's default SVG backend
# often converts text to paths, which makes final figure editing painful in
# Illustrator/Inkscape.  These settings preserve SVG text nodes and embed
# TrueType text in PDF/PS outputs where possible.
matplotlib.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
})

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd

from utils import condition_color_map, ensure_dir


def _get_cfg(config: Mapping[str, Any], *keys: str, default: Any = None) -> Any:
    node: Any = config
    for key in keys:
        if not isinstance(node, Mapping) or key not in node:
            return default
        node = node[key]
    return node


def _style_cfg(config: Mapping[str, Any], key: str, default: Any) -> Any:
    return _get_cfg(config, "figure", "style", key, default=default)


def _fmt_thousands(x: float, pos: int) -> str:
    return f"{int(x):,}"


def _panel_baseline_mode(config: Mapping[str, Any]) -> str:
    return str(
        _get_cfg(config, "figure", "baseline_mode", default=_get_cfg(config, "baseline", "default_mode", default="matched"))
    )


def _window_defaults(config: Mapping[str, Any]) -> tuple[int, list[int]]:
    windows = [int(w) for w in _get_cfg(config, "rolling", "windows", default=[25, 100, 250, 500])]
    default_window = int(_get_cfg(config, "rolling", "default_window", default=windows[len(windows) // 2]))
    profile_windows = [int(w) for w in _get_cfg(config, "figure", "profile_windows", default=[windows[0], windows[-1]])]
    if len(profile_windows) == 1:
        profile_windows = [profile_windows[0], profile_windows[0]]
    return default_window, profile_windows


def _color_for_condition(condition_label: str, config: Mapping[str, Any]) -> Optional[str]:
    palette = condition_color_map(config)
    return palette.get(condition_label)


def _display_condition_label(condition_label: str, config: Mapping[str, Any]) -> str:
    """Return manuscript-facing label for legends without changing data selectors."""
    label_map = _get_cfg(config, "figure", "condition_display_labels", default={})
    if isinstance(label_map, Mapping):
        return str(label_map.get(condition_label, condition_label))
    return condition_label


def _contrast_state_note(condition_label: str, baseline_label: str, config: Mapping[str, Any]) -> Optional[str]:
    """Return concise Dox/MYC-ER state note for a condition-baseline contrast."""
    notes = _get_cfg(config, "figure", "contrast_state_notes", default={})
    if not isinstance(notes, Mapping):
        return None
    keys = (
        f"{condition_label}|{baseline_label}",
        f"{condition_label} vs {baseline_label}",
        f"{condition_label}::{baseline_label}",
    )
    for key in keys:
        if key in notes:
            value = str(notes[key]).strip()
            return value or None
    return None


def _rank_shift_axis_label(label: str) -> str:
    """Append rank-shift wording without burying it after multiline state notes."""
    text = str(label)
    if "rank shift" in text.lower():
        return text
    if "\n" in text:
        first, rest = text.split("\n", 1)
        return f"{first} rank shift\n{rest}"
    return f"{text} rank shift"


def _shading_limits(n_genes: int, head_pct: float, tail_pct: float) -> tuple[int, int]:
    head_end = max(1, min(int(math.ceil(n_genes * (head_pct / 100.0))), n_genes))
    tail_start = max(1, min(int(math.floor(n_genes * (1 - tail_pct / 100.0))) + 1, n_genes))
    tail_start = max(tail_start, head_end + 1)
    return head_end, tail_start


def _shade_head_tail(ax: plt.Axes, n_genes: int, config: Mapping[str, Any]) -> tuple[int, int]:
    head_pct = float(_get_cfg(config, "analysis", "band_percentiles", "head", default=10))
    tail_pct = float(_get_cfg(config, "analysis", "band_percentiles", "tail", default=10))
    head_color = _get_cfg(config, "figure", "head_color", default="#7FC8F8")
    tail_color = _get_cfg(config, "figure", "tail_color", default="#F4A261")
    shade_alpha = float(_get_cfg(config, "figure", "shade_alpha", default=0.12))
    head_end, tail_start = _shading_limits(n_genes, head_pct=head_pct, tail_pct=tail_pct)
    ax.axvspan(1, head_end, color=head_color, alpha=shade_alpha, linewidth=0, zorder=0)
    ax.axvspan(tail_start, n_genes, color=tail_color, alpha=shade_alpha, linewidth=0, zorder=0)
    return head_end, tail_start


def _apply_tick_style(ax: plt.Axes, config: Mapping[str, Any]) -> None:
    ax.tick_params(axis="both", labelsize=float(_style_cfg(config, "tick_label_fontsize", 11)))


def _robust_abs_limit(values: pd.Series | np.ndarray, percentile: Any, min_abs: float, pad: float = 1.05) -> Optional[float]:
    """Return a symmetric robust display limit, or None when disabled.

    The source-data tables remain unmodified; this is only a plotting control to
    prevent a few unstable tail ranks from flattening the central rank geometry.
    """
    if percentile is None or str(percentile).lower() in {"none", "null", "false", "off"}:
        return None
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    pct = float(percentile)
    if pct <= 0 or pct > 100:
        raise ValueError(f"Display-limit percentile must be in (0, 100], got {pct}.")
    lim = float(np.nanpercentile(np.abs(arr), pct)) * float(pad)
    return max(float(min_abs), lim, 1.0)


def _set_symmetric_ylim(ax: plt.Axes, lim: Optional[float]) -> None:
    if lim is not None and np.isfinite(lim) and lim > 0:
        ax.set_ylim(-lim, lim)


def _set_common_axis_style(ax: plt.Axes, n_genes: int, config: Mapping[str, Any]) -> None:
    ax.set_xlim(1, n_genes)
    ax.axhline(
        0,
        linestyle=str(_style_cfg(config, "zero_line_style", "--")),
        color=str(_style_cfg(config, "zero_line_color", "#6B6B6B")),
        linewidth=float(_style_cfg(config, "zero_line_width", 1.0)),
        alpha=float(_style_cfg(config, "zero_line_alpha", 0.9)),
        zorder=1,
    )
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_thousands))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    _apply_tick_style(ax, config)


def _add_panel_letter(ax: plt.Axes, letter: str, config: Mapping[str, Any]) -> None:
    ax.text(
        float(_style_cfg(config, "panel_letter_x", -0.14)),
        float(_style_cfg(config, "panel_letter_y", 1.08)),
        letter.upper(),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=float(_style_cfg(config, "panel_letter_fontsize", 18)),
        fontweight=str(_style_cfg(config, "panel_letter_weight", "bold")),
        clip_on=False,
    )


def _select_comparison(
    long_df: pd.DataFrame,
    *,
    baseline_mode: str,
    condition_label: str,
    baseline_label: Optional[str] = None,
) -> pd.DataFrame:
    mask = (long_df["baseline_mode"] == baseline_mode) & (long_df["condition_label"] == condition_label)
    if baseline_label is not None:
        mask &= long_df["baseline_label"] == baseline_label
    subset = long_df.loc[mask].copy().sort_values("baseline_rank", kind="mergesort")
    if subset.empty:
        raise ValueError(
            f"No comparison rows found for baseline_mode='{baseline_mode}', "
            f"condition_label='{condition_label}', baseline_label='{baseline_label}'."
        )
    if subset["comparison_id"].nunique() != 1:
        raise ValueError(
            "Selector resolved to multiple comparisons. Provide an explicit baseline_label in config if needed."
        )
    return subset


def _select_rolling(
    rolling_df: pd.DataFrame,
    *,
    baseline_mode: str,
    condition_label: str,
    window: int,
    baseline_label: Optional[str] = None,
) -> pd.DataFrame:
    mask = (
        (rolling_df["baseline_mode"] == baseline_mode)
        & (rolling_df["condition_label"] == condition_label)
        & (rolling_df["window"] == int(window))
    )
    if baseline_label is not None:
        mask &= rolling_df["baseline_label"] == baseline_label
    subset = rolling_df.loc[mask].copy().sort_values("baseline_rank", kind="mergesort")
    if subset.empty:
        raise ValueError(
            f"No rolling summary found for baseline_mode='{baseline_mode}', "
            f"condition_label='{condition_label}', baseline_label='{baseline_label}', window={window}."
        )
    if subset["comparison_id"].nunique() != 1:
        raise ValueError(
            "Rolling selector resolved to multiple comparisons. Provide an explicit baseline_label in config if needed."
        )
    return subset


def prepare_scatter_panel_data(
    long_df: pd.DataFrame,
    rolling_df: pd.DataFrame,
    *,
    condition_label: str,
    baseline_mode: str,
    window: int,
    baseline_label: Optional[str] = None,
) -> pd.DataFrame:
    scatter_df = _select_comparison(
        long_df,
        baseline_mode=baseline_mode,
        condition_label=condition_label,
        baseline_label=baseline_label,
    )
    summary_df = _select_rolling(
        rolling_df,
        baseline_mode=baseline_mode,
        condition_label=condition_label,
        baseline_label=baseline_label,
        window=window,
    )[["baseline_rank", "median", "q25", "q75"]].rename(
        columns={"median": "rolling_median", "q25": "rolling_q25", "q75": "rolling_q75"}
    )
    return scatter_df.merge(summary_df, on="baseline_rank", how="left")


def prepare_profile_panel_data(
    rolling_df: pd.DataFrame,
    *,
    condition_labels: list[str],
    baseline_mode: str,
    window: int,
    baseline_labels: Optional[list[str]] = None,
) -> pd.DataFrame:
    mask = (rolling_df["baseline_mode"] == baseline_mode) & (rolling_df["window"] == int(window))
    subset = rolling_df.loc[mask & rolling_df["condition_label"].isin(condition_labels)].copy()
    if baseline_labels is not None:
        subset = subset.loc[subset["baseline_label"].isin(baseline_labels)].copy()
    if subset.empty:
        raise ValueError(
            f"No positional profile data found for baseline_mode='{baseline_mode}', "
            f"conditions={condition_labels}, window={window}."
        )
    # Preserve the manuscript/config order in profile panels.  A lexical sort would
    # place "24 h ..." before "4 h ...", which reverses the intended early-to-late
    # visual reading and legend order.
    order_lookup = {str(label): idx for idx, label in enumerate(condition_labels)}
    subset["_condition_order"] = subset["condition_label"].map(order_lookup).fillna(len(order_lookup)).astype(int)
    subset = subset.sort_values(["_condition_order", "baseline_rank"], kind="mergesort").drop(columns=["_condition_order"])
    return subset.reset_index(drop=True)


def prepare_hexbin_panel_data(
    long_df: pd.DataFrame,
    *,
    x_selector: Mapping[str, Any],
    y_selector: Mapping[str, Any],
) -> pd.DataFrame:
    x_df = _select_comparison(
        long_df,
        baseline_mode=str(x_selector["baseline_mode"]),
        condition_label=str(x_selector["condition_label"]),
        baseline_label=x_selector.get("baseline_label"),
    )
    y_df = _select_comparison(
        long_df,
        baseline_mode=str(y_selector["baseline_mode"]),
        condition_label=str(y_selector["condition_label"]),
        baseline_label=y_selector.get("baseline_label"),
    )

    x_cols = [
        "gene_id",
        "comparison_id",
        "baseline_mode",
        "baseline_label",
        "condition_label",
        "baseline_rank",
        "condition_rank",
        "rank_shift",
        "band",
    ]
    y_cols = list(x_cols)
    merged = x_df[x_cols].merge(
        y_df[y_cols],
        on="gene_id",
        how="inner",
        suffixes=("_x", "_y"),
    )
    merged = merged.rename(
        columns={
            "rank_shift_x": "x_rank_shift",
            "rank_shift_y": "y_rank_shift",
            "baseline_rank_x": "x_baseline_rank",
            "baseline_rank_y": "y_baseline_rank",
            "condition_rank_x": "x_condition_rank",
            "condition_rank_y": "y_condition_rank",
            "condition_label_x": "x_condition_label",
            "condition_label_y": "y_condition_label",
            "baseline_mode_x": "x_baseline_mode",
            "baseline_mode_y": "y_baseline_mode",
            "baseline_label_x": "x_baseline_label",
            "baseline_label_y": "y_baseline_label",
            "band_x": "x_band",
            "band_y": "y_band",
            "comparison_id_x": "x_comparison_id",
            "comparison_id_y": "y_comparison_id",
        }
    )
    merged["x_axis_label"] = str(x_selector.get("axis_label") or f"{merged['x_condition_label'].iloc[0]} on {merged['x_baseline_label'].iloc[0]}")
    merged["y_axis_label"] = str(y_selector.get("axis_label") or f"{merged['y_condition_label'].iloc[0]} on {merged['y_baseline_label'].iloc[0]}")
    return merged


def build_panel_source_data(
    long_df: pd.DataFrame,
    rolling_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> Dict[str, pd.DataFrame]:
    baseline_mode = _panel_baseline_mode(config)
    default_window, profile_windows = _window_defaults(config)

    panel_a_condition = str(_get_cfg(config, "figure", "panel_a", "condition_label", default="4 h DT"))
    panel_b_condition = str(_get_cfg(config, "figure", "panel_b", "condition_label", default="24 h DT"))
    profile_conditions = list(_get_cfg(config, "figure", "profile_conditions", default=[]))
    if not profile_conditions:
        profile_conditions = list(_get_cfg(config, "analysis", "condition_order", default=[]))
    profile_baselines = _get_cfg(config, "figure", "profile_baseline_labels", default=None)

    panel_data = {
        "A": prepare_scatter_panel_data(
            long_df,
            rolling_df,
            condition_label=panel_a_condition,
            baseline_mode=baseline_mode,
            window=default_window,
            baseline_label=_get_cfg(config, "figure", "panel_a", "baseline_label", default=None),
        ),
        "B": prepare_scatter_panel_data(
            long_df,
            rolling_df,
            condition_label=panel_b_condition,
            baseline_mode=baseline_mode,
            window=default_window,
            baseline_label=_get_cfg(config, "figure", "panel_b", "baseline_label", default=None),
        ),
        "C": prepare_profile_panel_data(
            rolling_df,
            condition_labels=profile_conditions,
            baseline_mode=baseline_mode,
            window=profile_windows[0],
            baseline_labels=profile_baselines,
        ),
        "D": prepare_profile_panel_data(
            rolling_df,
            condition_labels=profile_conditions,
            baseline_mode=baseline_mode,
            window=profile_windows[1],
            baseline_labels=profile_baselines,
        ),
    }

    robustness_panels = list(_get_cfg(config, "robustness", "panels", default=[]))
    if len(robustness_panels) >= 1:
        panel_data["E"] = prepare_hexbin_panel_data(
            long_df,
            x_selector=robustness_panels[0]["x"],
            y_selector=robustness_panels[0]["y"],
        )
    if len(robustness_panels) >= 2:
        panel_data["F"] = prepare_hexbin_panel_data(
            long_df,
            x_selector=robustness_panels[1]["x"],
            y_selector=robustness_panels[1]["y"],
        )

    return panel_data


def _plot_scatter_panel(
    ax: plt.Axes,
    panel_df: pd.DataFrame,
    config: Mapping[str, Any],
    title: str,
) -> None:
    n_genes = int(panel_df["n_genes"].iloc[0]) if "n_genes" in panel_df.columns else int(panel_df["baseline_rank"].max())
    condition_label = str(panel_df["condition_label"].iloc[0])
    baseline_label = str(panel_df["baseline_label"].iloc[0])
    color = _color_for_condition(condition_label, config) or "#1f77b4"

    ax.scatter(
        panel_df["baseline_rank"],
        panel_df["rank_shift"],
        s=float(_get_cfg(config, "figure", "scatter_size", default=4)),
        alpha=float(_get_cfg(config, "figure", "scatter_alpha", default=0.12)),
        color=color,
        linewidths=0,
        zorder=2,
    )
    ax.fill_between(
        panel_df["baseline_rank"],
        panel_df["rolling_q25"],
        panel_df["rolling_q75"],
        color=color,
        alpha=float(_get_cfg(config, "figure", "ribbon_alpha", default=0.20)),
        linewidth=0,
        zorder=3,
    )
    ax.plot(panel_df["baseline_rank"], panel_df["rolling_median"], color="black", linewidth=2.2, zorder=4)

    head_end, tail_start = _shade_head_tail(ax, n_genes, config)
    head_median = panel_df.loc[panel_df["baseline_rank"] <= head_end, "rank_shift"].median()
    tail_median = panel_df.loc[panel_df["baseline_rank"] >= tail_start, "rank_shift"].median()

    state_note = _contrast_state_note(condition_label, baseline_label, config)
    if state_note:
        annotation_text = (
            f"Contrast: {condition_label} vs {baseline_label}\n"
            f"{state_note}\n"
            f"Head med: {head_median:.0f} | Tail med: {tail_median:.0f}"
        )
    else:
        annotation_text = (
            f"Condition: {condition_label}\n"
            f"Baseline: {baseline_label}\n"
            f"Head med: {head_median:.0f} | Tail med: {tail_median:.0f}"
        )

    ax.text(
        0.02,
        0.98,
        annotation_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=float(_style_cfg(config, "annotation_fontsize", 10)),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.92, edgecolor="0.7"),
        zorder=5,
    )

    ax.set_title(title, fontsize=float(_style_cfg(config, "title_fontsize", 14)))
    ax.set_xlabel("Baseline rank (head → tail)", fontsize=float(_style_cfg(config, "axis_label_fontsize", 13)))
    ax.set_ylabel("Rank shift (baseline − condition)", fontsize=float(_style_cfg(config, "axis_label_fontsize", 13)))
    _set_common_axis_style(ax, n_genes, config)
    scatter_lim = _robust_abs_limit(
        panel_df["rank_shift"],
        _get_cfg(config, "figure", "scatter_y_limit_percentile", default=None),
        float(_get_cfg(config, "figure", "scatter_y_limit_min_abs", default=750)),
        float(_get_cfg(config, "figure", "scatter_y_limit_pad", default=1.05)),
    )
    _set_symmetric_ylim(ax, scatter_lim)



def _plot_profile_panel(
    ax: plt.Axes,
    panel_df: pd.DataFrame,
    config: Mapping[str, Any],
    title: str,
    show_iqr: bool = False,
) -> None:
    n_genes = int(panel_df["n_genes"].max())
    _shade_head_tail(ax, n_genes, config)

    for condition_label, group in panel_df.groupby("condition_label", sort=False):
        group = group.sort_values("baseline_rank", kind="mergesort")
        color = _color_for_condition(str(condition_label), config)
        display_label = _display_condition_label(str(condition_label), config)
        ax.plot(group["baseline_rank"], group["median"], linewidth=2.2, label=display_label, color=color, zorder=3)
        if show_iqr:
            ax.fill_between(group["baseline_rank"], group["q25"], group["q75"], alpha=0.12, color=color, zorder=2)

    ax.set_title(title, fontsize=float(_style_cfg(config, "title_fontsize", 14)))
    ax.set_xlabel("Baseline rank (head → tail)", fontsize=float(_style_cfg(config, "axis_label_fontsize", 13)))
    ax.set_ylabel("Rolling median of rank shift", fontsize=float(_style_cfg(config, "axis_label_fontsize", 13)))
    _set_common_axis_style(ax, n_genes, config)
    profile_lim = _robust_abs_limit(
        panel_df["median"],
        _get_cfg(config, "figure", "profile_y_limit_percentile", default=None),
        float(_get_cfg(config, "figure", "profile_y_limit_min_abs", default=250)),
        float(_get_cfg(config, "figure", "profile_y_limit_pad", default=1.08)),
    )
    _set_symmetric_ylim(ax, profile_lim)
    legend = ax.legend(
        title=str(_get_cfg(config, "figure", "profile_legend_title", default="Activation contrast")),
        fontsize=float(_get_cfg(config, "figure", "profile_legend_fontsize", default=_style_cfg(config, "legend_fontsize", 10))),
        frameon=True,
        loc=str(_get_cfg(config, "figure", "profile_legend_loc", default="upper left")),
        title_fontsize=float(_get_cfg(config, "figure", "profile_legend_title_fontsize", default=_style_cfg(config, "legend_title_fontsize", 10))),
        labelspacing=float(_get_cfg(config, "figure", "profile_legend_labelspacing", default=0.55)),
        borderpad=float(_get_cfg(config, "figure", "profile_legend_borderpad", default=0.45)),
        handlelength=float(_get_cfg(config, "figure", "profile_legend_handlelength", default=2.0)),
        handletextpad=float(_get_cfg(config, "figure", "profile_legend_handletextpad", default=0.7)),
    )
    if legend is not None:
        legend.get_frame().set_alpha(0.95)



def _smooth_histogram2d(hist: np.ndarray, passes: int = 2) -> np.ndarray:
    """Lightweight histogram smoothing without adding a SciPy dependency."""
    arr = np.asarray(hist, dtype=float)
    if arr.size == 0 or passes <= 0:
        return arr
    kernel = np.asarray([1.0, 4.0, 6.0, 4.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    out = arr.copy()
    for _ in range(int(passes)):
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), axis=0, arr=out)
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), axis=1, arr=out)
    return out


def _density_enclosed_thresholds(hist: np.ndarray, fractions: list[float]) -> list[float]:
    """Return density thresholds enclosing requested cumulative mass fractions.

    These are soft density envelopes, not hard biological boundaries.  For
    example, 0.50 marks the contour around the densest half of the displayed gene
    density; 0.90 marks a broader envelope.
    """
    flat = np.asarray(hist, dtype=float).ravel()
    flat = flat[np.isfinite(flat) & (flat > 0)]
    if flat.size == 0:
        return []
    order = np.sort(flat)[::-1]
    csum = np.cumsum(order)
    total = float(csum[-1])
    levels: list[float] = []
    for frac in fractions:
        frac = float(frac)
        if not (0 < frac < 1):
            continue
        idx = int(np.searchsorted(csum, frac * total, side="left"))
        idx = max(0, min(idx, order.size - 1))
        levels.append(float(order[idx]))
    return sorted(set(levels))


def _finite_xy_for_hexbin(panel_df: pd.DataFrame, lim: float) -> tuple[np.ndarray, np.ndarray]:
    x = pd.to_numeric(panel_df["x_rank_shift"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(panel_df["y_rank_shift"], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y) & (np.abs(x) <= lim) & (np.abs(y) <= lim)
    return x[finite], y[finite]


def _histogram2d_for_hex_guides(
    panel_df: pd.DataFrame,
    lim: float,
    bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x, y = _finite_xy_for_hexbin(panel_df, lim)
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins, range=[[-lim, lim], [-lim, lim]])
    return x, y, hist, xedges, yedges


def _add_density_contours(
    ax: plt.Axes,
    x: pd.Series,
    y: pd.Series,
    lim: float,
    config: Mapping[str, Any],
) -> None:
    """Add soft density envelopes to E/F without hard class boundaries."""
    if not bool(_get_cfg(config, "figure", "hexbin_show_density_contours", default=True)):
        return

    bins = int(_get_cfg(config, "figure", "hexbin_contour_bins", default=60))
    hist, xedges, yedges = np.histogram2d(x.to_numpy(), y.to_numpy(), bins=bins, range=[[-lim, lim], [-lim, lim]])
    if np.count_nonzero(hist > 0) < 4:
        return

    smooth_passes = int(_get_cfg(config, "figure", "hexbin_contour_smoothing_passes", default=2))
    hist_for_contour = _smooth_histogram2d(hist, passes=smooth_passes)

    # Enclosed-mass contours behave like soft boundaries of the dense core rather
    # than arbitrary fractions of the maximum bin height.
    enclosed = _get_cfg(
        config,
        "figure",
        "hexbin_contour_cumulative_fractions",
        default=_get_cfg(config, "figure", "hexbin_contour_enclosed_fractions", default=[0.50, 0.80]),
    )
    levels = _density_enclosed_thresholds(hist_for_contour, [float(v) for v in enclosed])

    # Backward-compatible fallback: allow direct density levels in config.
    if not levels:
        raw_levels = _get_cfg(config, "figure", "hexbin_contour_levels", default=[0.35, 0.60, 0.80])
        max_count = float(np.nanmax(hist_for_contour))
        levels = []
        for value in raw_levels:
            level = float(value)
            if 0 < level < 1:
                level = level * max_count
            if 0 < level <= max_count:
                levels.append(level)
        levels = sorted(set(levels))
    if not levels:
        return

    xcenters = (xedges[:-1] + xedges[1:]) / 2.0
    ycenters = (yedges[:-1] + yedges[1:]) / 2.0
    xx, yy = np.meshgrid(xcenters, ycenters)

    # Draw a subtle dark underlay plus a white overlay. This keeps contours
    # visible over both the dense core and sparse white regions.
    underlay_width = float(_get_cfg(config, "figure", "hexbin_contour_underlay_linewidth", default=2.3))
    if underlay_width > 0:
        ax.contour(
            xx,
            yy,
            hist_for_contour.T,
            levels=levels,
            colors=str(_get_cfg(config, "figure", "hexbin_contour_underlay_color", default="#202020")),
            linewidths=underlay_width,
            alpha=float(_get_cfg(config, "figure", "hexbin_contour_underlay_alpha", default=0.38)),
            zorder=4,
        )
    ax.contour(
        xx,
        yy,
        hist_for_contour.T,
        levels=levels,
        colors=str(_get_cfg(config, "figure", "hexbin_contour_color", default="#FFFFFF")),
        linewidths=float(_get_cfg(config, "figure", "hexbin_contour_linewidth", default=1.2)),
        alpha=float(_get_cfg(config, "figure", "hexbin_contour_alpha", default=0.92)),
        zorder=5,
    )


def _add_hexbin_zero_axes(ax: plt.Axes, config: Mapping[str, Any]) -> None:
    if not bool(_get_cfg(config, "figure", "hexbin_show_zero_axes", default=True)):
        return
    color = str(_get_cfg(config, "figure", "hexbin_zero_axis_color", default="#6B6B6B"))
    alpha = float(_get_cfg(config, "figure", "hexbin_zero_axis_alpha", default=0.35))
    width = float(_get_cfg(config, "figure", "hexbin_zero_axis_width", default=0.8))
    ax.axhline(0, color=color, linewidth=width, alpha=alpha, zorder=3)
    ax.axvline(0, color=color, linewidth=width, alpha=alpha, zorder=3)


def _density_core_center(
    panel_df: pd.DataFrame,
    lim: float,
    config: Mapping[str, Any],
) -> Optional[tuple[float, float, float]]:
    """Return a density-weighted center of the densest E/F cloud and y-x offset."""
    method = str(_get_cfg(config, "figure", "hexbin_core_method", default="density_peak")).lower()
    x, y = _finite_xy_for_hexbin(panel_df, lim)
    if x.size < 10:
        return None
    if method in {"median", "coordinate_median"}:
        core_x = float(np.nanmedian(x))
        core_y = float(np.nanmedian(y))
        return (core_x, core_y, core_y - core_x) if np.isfinite(core_x) and np.isfinite(core_y) else None

    bins = int(_get_cfg(config, "figure", "hexbin_core_bins", default=_get_cfg(config, "figure", "hexbin_contour_bins", default=60)))
    _, _, hist, xedges, yedges = _histogram2d_for_hex_guides(panel_df, lim, bins)
    if hist.max() <= 0:
        return None
    smooth_passes = int(_get_cfg(config, "figure", "hexbin_core_smoothing_passes", default=_get_cfg(config, "figure", "hexbin_contour_smoothing_passes", default=2)))
    core_hist = _smooth_histogram2d(hist, passes=smooth_passes)
    positive = core_hist[core_hist > 0]
    fraction = float(_get_cfg(config, "figure", "hexbin_core_density_fraction", default=0.55))
    fraction = min(max(fraction, 0.05), 1.0)
    threshold = max(float(positive.max()) * fraction, float(np.nanpercentile(positive, 90.0)))
    mask = core_hist >= threshold
    if not np.any(mask):
        mask = core_hist == core_hist.max()

    xcenters = (xedges[:-1] + xedges[1:]) / 2.0
    ycenters = (yedges[:-1] + yedges[1:]) / 2.0
    xx, yy = np.meshgrid(xcenters, ycenters, indexing="ij")
    weights = core_hist * mask
    total = float(weights.sum())
    if total <= 0:
        return None
    core_x = float((xx * weights).sum() / total)
    core_y = float((yy * weights).sum() / total)
    if not (np.isfinite(core_x) and np.isfinite(core_y)):
        return None
    return core_x, core_y, core_y - core_x


def _plot_parallel_offset_line(ax: plt.Axes, lim: float, offset: float, config: Mapping[str, Any]) -> None:
    """Draw y=x+offset as a soft guide through the displaced dense core."""
    if not bool(_get_cfg(config, "figure", "hexbin_show_core_parallel_line", default=True)):
        return
    if not np.isfinite(offset):
        return
    x0 = max(-lim, -lim - offset)
    x1 = min(lim, lim - offset)
    if x0 >= x1:
        return
    xs = np.array([x0, x1], dtype=float)
    ys = xs + float(offset)

    ax.plot(
        xs,
        ys,
        linestyle=str(_get_cfg(config, "figure", "hexbin_core_parallel_line_style", default="-")),
        color=str(_get_cfg(config, "figure", "hexbin_core_parallel_line_underlay_color", default="#FFFFFF")),
        linewidth=float(_get_cfg(config, "figure", "hexbin_core_parallel_line_underlay_width", default=3.2)),
        alpha=float(_get_cfg(config, "figure", "hexbin_core_parallel_line_underlay_alpha", default=0.88)),
        zorder=5.5,
    )
    ax.plot(
        xs,
        ys,
        linestyle=str(_get_cfg(config, "figure", "hexbin_core_parallel_line_style", default="-")),
        color=str(_get_cfg(config, "figure", "hexbin_core_parallel_line_color", default="#202020")),
        linewidth=float(_get_cfg(config, "figure", "hexbin_core_parallel_line_width", default=1.15)),
        alpha=float(_get_cfg(config, "figure", "hexbin_core_parallel_line_alpha", default=0.86)),
        zorder=5.6,
    )


def _add_hexbin_core_marker(ax: plt.Axes, panel_df: pd.DataFrame, lim: float, config: Mapping[str, Any]) -> None:
    """Mark the dense core and its modest offset from the identity line.

    This is intentionally a soft guide, not a hard boundary. It helps the reader
    see the manuscript claim that Dox ON and Dox OFF rank-shift maps share the
    same overall geometry while the densest portion is modestly displaced.
    """
    if not bool(_get_cfg(config, "figure", "hexbin_show_core_marker", default=True)):
        return

    core = _density_core_center(panel_df, lim, config)
    if core is None:
        return
    core_x, core_y, offset = core

    _plot_parallel_offset_line(ax, lim, offset, config)

    identity_x = (core_x + core_y) / 2.0
    identity_y = identity_x
    if bool(_get_cfg(config, "figure", "hexbin_show_core_projection", default=True)):
        # Perpendicular projection from the dense core to the identity line.
        ax.plot(
            [core_x, identity_x],
            [core_y, identity_y],
            color=str(_get_cfg(config, "figure", "hexbin_core_projection_color", default="#202020")),
            linewidth=float(_get_cfg(config, "figure", "hexbin_core_projection_width", default=1.1)),
            alpha=float(_get_cfg(config, "figure", "hexbin_core_projection_alpha", default=0.90)),
            zorder=7,
        )

    if bool(_get_cfg(config, "figure", "hexbin_show_identity_foot_marker", default=True)):
        ax.scatter(
            [identity_x],
            [identity_y],
            s=float(_get_cfg(config, "figure", "hexbin_identity_foot_marker_size", default=32)),
            marker=str(_get_cfg(config, "figure", "hexbin_identity_foot_marker", default="o")),
            facecolors=str(_get_cfg(config, "figure", "hexbin_identity_foot_marker_facecolor", default="none")),
            edgecolors=str(_get_cfg(config, "figure", "hexbin_identity_foot_marker_edgecolor", default="#202020")),
            linewidths=float(_get_cfg(config, "figure", "hexbin_identity_foot_marker_linewidth", default=1.0)),
            zorder=7,
        )

    ax.scatter(
        [core_x],
        [core_y],
        s=float(_get_cfg(config, "figure", "hexbin_core_marker_size", default=76)),
        marker=str(_get_cfg(config, "figure", "hexbin_core_marker", default="o")),
        facecolors=str(_get_cfg(config, "figure", "hexbin_core_marker_facecolor", default="#FFFFFF")),
        edgecolors=str(_get_cfg(config, "figure", "hexbin_core_marker_edgecolor", default="#202020")),
        linewidths=float(_get_cfg(config, "figure", "hexbin_core_marker_linewidth", default=1.4)),
        zorder=8,
    )

    if bool(_get_cfg(config, "figure", "hexbin_annotate_core_offset", default=True)):
        ax.text(
            0.03,
            0.97,
            f"Dense-core guide\nΔ(y−x) = {offset:+.0f} ranks",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=float(_style_cfg(config, "hexbin_annotation_fontsize", _style_cfg(config, "annotation_fontsize", 10))),
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.90, edgecolor="0.7"),
            zorder=8,
        )


def _plot_hexbin_panel(ax: plt.Axes, panel_df: pd.DataFrame, config: Mapping[str, Any], title: str) -> None:
    gridsize = int(_get_cfg(config, "figure", "hexbin_gridsize", default=60))
    values = panel_df[["x_rank_shift", "y_rank_shift"]].to_numpy()
    limit_pct = float(_get_cfg(config, "figure", "hexbin_limit_percentile", default=99.0))
    lim = float(np.nanpercentile(np.abs(values), limit_pct))
    lim = max(1.0, lim)
    hb = ax.hexbin(
        panel_df["x_rank_shift"],
        panel_df["y_rank_shift"],
        gridsize=gridsize,
        mincnt=int(_get_cfg(config, "figure", "hexbin_mincnt", default=1)),
        bins=_get_cfg(config, "figure", "hexbin_bins", default="log"),
        cmap=_get_cfg(config, "figure", "hexbin_cmap", default="viridis"),
        extent=(-lim, lim, -lim, lim),
        linewidths=0,
    )
    _add_density_contours(ax, panel_df["x_rank_shift"], panel_df["y_rank_shift"], lim, config)
    _add_hexbin_zero_axes(ax, config)
    ax.plot(
        [-lim, lim],
        [-lim, lim],
        linestyle=str(_style_cfg(config, "identity_line_style", "--")),
        color=str(_style_cfg(config, "identity_line_color", "#2F2F2F")),
        linewidth=float(_style_cfg(config, "identity_line_width", 1.5)),
        alpha=float(_style_cfg(config, "identity_line_alpha", 0.95)),
        zorder=5,
    )
    _add_hexbin_core_marker(ax, panel_df, lim, config)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    if bool(_get_cfg(config, "figure", "hexbin_equal_aspect", default=True)):
        ax.set_aspect("equal", adjustable="box")
    x_label = str(panel_df["x_axis_label"].iloc[0])
    y_label = str(panel_df["y_axis_label"].iloc[0])
    ax.set_xlabel(_rank_shift_axis_label(x_label), fontsize=float(_style_cfg(config, "axis_label_fontsize", 13)))
    ax.set_ylabel(_rank_shift_axis_label(y_label), fontsize=float(_style_cfg(config, "axis_label_fontsize", 13)))
    ax.set_title(title, fontsize=float(_style_cfg(config, "title_fontsize", 14)))
    _apply_tick_style(ax, config)
    cbar = plt.colorbar(hb, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(str(_get_cfg(config, "figure", "hexbin_colorbar_label", default="Local gene density")), fontsize=float(_style_cfg(config, "colorbar_label_fontsize", 10)))
    cbar.ax.tick_params(labelsize=float(_style_cfg(config, "colorbar_tick_fontsize", 9)))



def plot_figure5(
    long_df: pd.DataFrame,
    rolling_df: pd.DataFrame,
    config: Mapping[str, Any],
    outdir: Path | str,
) -> Path:
    outdir = ensure_dir(Path(outdir))
    figures_dir = ensure_dir(outdir / "figures")

    panel_data = build_panel_source_data(long_df, rolling_df, config)

    panel_a_title = str(_get_cfg(config, "figure", "panel_a", "title", default="Christmas tree 4 h"))
    panel_b_title = str(_get_cfg(config, "figure", "panel_b", "title", default="Christmas tree 24 h"))
    _, profile_windows = _window_defaults(config)

    panel_c_title = str(
        _get_cfg(
            config,
            "figure",
            "panel_c_title",
            default=f"Positional profile ({profile_windows[0]}-gene rolling window)",
        )
    )
    panel_d_title = str(
        _get_cfg(
            config,
            "figure",
            "panel_d_title",
            default=f"Positional profile ({profile_windows[1]}-gene rolling window)",
        )
    )

    robustness_panels = list(_get_cfg(config, "robustness", "panels", default=[]))
    panel_e_title = robustness_panels[0]["title"] if len(robustness_panels) >= 1 else "Robustness comparison E"
    panel_f_title = robustness_panels[1]["title"] if len(robustness_panels) >= 2 else "Robustness comparison F"

    fig_size = _get_cfg(config, "figure", "figure_size", default=[15, 13])
    fig, axes = plt.subplots(3, 2, figsize=(float(fig_size[0]), float(fig_size[1])))

    _plot_scatter_panel(axes[0, 0], panel_data["A"], config, panel_a_title)
    _add_panel_letter(axes[0, 0], "A", config)

    _plot_scatter_panel(axes[0, 1], panel_data["B"], config, panel_b_title)
    _add_panel_letter(axes[0, 1], "B", config)

    _plot_profile_panel(
        axes[1, 0],
        panel_data["C"],
        config,
        panel_c_title,
        show_iqr=bool(_get_cfg(config, "figure", "profile_show_iqr", default=False)),
    )
    _add_panel_letter(axes[1, 0], "C", config)

    _plot_profile_panel(
        axes[1, 1],
        panel_data["D"],
        config,
        panel_d_title,
        show_iqr=bool(_get_cfg(config, "figure", "profile_show_iqr", default=False)),
    )
    _add_panel_letter(axes[1, 1], "D", config)

    if "E" in panel_data:
        _plot_hexbin_panel(axes[2, 0], panel_data["E"], config, panel_e_title)
        _add_panel_letter(axes[2, 0], "E", config)
    else:
        axes[2, 0].axis("off")

    if "F" in panel_data:
        _plot_hexbin_panel(axes[2, 1], panel_data["F"], config, panel_f_title)
        _add_panel_letter(axes[2, 1], "F", config)
    else:
        axes[2, 1].axis("off")

    plt.tight_layout(rect=(0.06, 0.03, 0.99, 0.98))

    prefix = str(_get_cfg(config, "outputs", "figure_prefix", default="Figure5_rank_geometry"))
    formats = _get_cfg(config, "outputs", "formats", default=["png", "svg"])
    formats = [str(fmt).lower().lstrip(".") for fmt in formats]
    if "png" not in formats:
        formats.insert(0, "png")

    output_paths: dict[str, Path] = {}
    for fmt in formats:
        out_path = figures_dir / f"{prefix}.{fmt}"
        if fmt == "png":
            plt.savefig(out_path, dpi=int(_get_cfg(config, "figure", "dpi", default=300)))
        else:
            plt.savefig(out_path)
        output_paths[fmt] = out_path

    plt.close(fig)

    return output_paths["png"]
