#!/usr/bin/env python3
"""Compose publication-ready Figure 6 and Supplementary Fig. 2.

This script reads a completed figure5 run directory and writes
PNG, SVG, and PDF outputs. It intentionally keeps the main figure A-D and the
support figure A-D. Extra gene diagnostics (AUC>=0.80, MYC/DOX overlap, DDR
overlap) are not included by default; create them only with --make-s3.

Default manuscript-aligned choices:
  - no figure-level title inside artwork (caption supplies title)
  - no p-values in metric panels unless --show-pvalues is set
  - score panels show individual replicates and thin replicate connectors
  - default support output is the local two-panel Supplementary Fig. 2
    (Dox-control metrics + HCPT-local heatmap); the full gene-screen S2 is optional
  - full-mode zero-significant-gene tasks are shown explicitly as undefined in region-composition panel
"""
from __future__ import annotations

import argparse
import os
import re
import textwrap
from pathlib import Path
from typing import Iterable, Optional, Sequence

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.image as mpimg
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

MAIN_TASKS = ["MYC@4h", "MYC@24h"]
CONTROL_TASKS = ["DOX@4h", "DOX@24h"]
TASK_ORDER = MAIN_TASKS + CONTROL_TASKS
CONTEXT_ORDER = ["None_CPT", "LCPT", "HCPT"]
CONTEXT_LABEL = {"None_CPT": "no CPT", "LCPT": "LCPT", "HCPT": "HCPT", "None_MYC": "MYC OFF", "MYC": "MYC ON"}
CONTEXT_MARKERS = {"None_CPT": "o", "LCPT": "s", "HCPT": "^", "None_MYC": "o", "MYC": "s"}
LABEL_COLORS = {0: "#E69F00", 1: "#2CA02C"}
METRIC_COLORS = {"balanced_accuracy": "#1F77B4", "roc_auc": "#FF7F0E"}
METRIC_MARKERS = {"balanced_accuracy": "o", "roc_auc": "s"}
REGION_COLORS = {"head": "#1f77b4", "mid": "#ff7f0e", "tail": "#2ca02c"}
CV_LABELS = {"leave_one_sample_out": "Leave-one-sample-out", "leave_one_context_out": "Leave-one-context-out"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, help="Completed Figure 6 run directory")
    p.add_argument("--output-dir", default=None, help="Output directory; default <run-dir>/publication_figures")
    p.add_argument("--matrix", default=None, help="Optional expression matrix for vector HCPT heatmap")
    p.add_argument("--model", default="ridge_logistic")
    p.add_argument("--main-cv", default="leave_one_sample_out")
    p.add_argument("--robust-cv", default="leave_one_context_out")
    p.add_argument("--formats", nargs="+", default=["png", "svg", "pdf"], choices=["png", "svg", "pdf"])
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--with-title", action="store_true", help="Add figure-level titles inside artwork; default off")
    p.add_argument("--show-pvalues", action="store_true", help="Show permutation p-values; default off")
    p.add_argument("--show-metric-values", action="store_true", help="Annotate metric values in metric panels; default off for publication cleanliness")
    p.add_argument("--no-replicate-lines", action="store_true", help="Hide thin replicate connectors")
    p.add_argument("--highlight-samples", nargs="*", default=["29_24D_H_CPT_RP1", "31_24D_H_CPT_RP2"])
    p.add_argument("--support-mode", choices=["full", "local", "none"], default="local",
                   help="local=publication default Dox metrics + HCPT local only; full=optional A-D Dox/gene screen support; none=skip support")
    p.add_argument("--gene-count-scale", choices=["linear", "log"], default="linear")
    p.add_argument("--make-s3", action="store_true", help="Also compose optional extra gene diagnostics; off by default")
    return p.parse_args()


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)



def write_publication_manifest(out_dir: Path, args: argparse.Namespace, written: list[Path]) -> None:
    """Write a small manifest documenting publication-art settings."""
    import json
    manifest = {
        "script": "compose_publication_figure6.py",
        "run_dir": str(Path(args.run_dir).resolve()),
        "support_mode": args.support_mode,
        "model": args.model,
        "main_cv": args.main_cv,
        "robust_cv": args.robust_cv,
        "formats": list(args.formats),
        "show_pvalues": bool(args.show_pvalues),
        "show_metric_values": bool(args.show_metric_values),
        "with_title": bool(args.with_title),
        "matrix": str(Path(args.matrix).resolve()) if args.matrix else None,
        "outputs": [str(p) for p in written],
        "notes": [
            "Main Figure 6 is metric-first and sample-level.",
            "Metric p-values and numeric labels are hidden by default.",
            "Support-mode local is publication default; full gene-screen support is optional.",
            "No t-SNE/embedding panels are included in publication-ready Figure 6 by default."
        ],
    }
    with open(out_dir / "publication_figure6_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)

def save_all(fig: plt.Figure, out_dir: Path, stem: str, formats: Sequence[str], dpi: int) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        path = out_dir / f"{stem}.{fmt}"
        if fmt == "png":
            fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
        written.append(path)
    return written


def panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.065, 1.035, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top", ha="left")


def clean_sample_id(sample_id: str) -> str:
    s = str(sample_id)
    m = re.match(r"^\d+_(.*)_RP(\d+)$", s)
    if m:
        return f"{m.group(1)} R{m.group(2)}"
    return s.replace("_RP", " R")


def short_rep_label(row: pd.Series) -> str:
    label = "ON" if int(row.get("label", 0)) == 1 else "OFF"
    rep = row.get("replicate", None)
    try:
        rep_s = f"R{int(rep)}"
    except Exception:
        m = re.search(r"RP(\d+)$", str(row.get("sample_id", "")))
        rep_s = f"R{m.group(1)}" if m else ""
    return f"{label} {rep_s}"


def symbol_from_gene(g: str) -> str:
    s = str(g)
    return s.split("_")[-1] if "_" in s else s


def detect_sep(path: Path) -> str:
    return "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","


def plot_task_panel(ax: plt.Axes, task_df: pd.DataFrame) -> None:
    """Compact task-design grid for publication artwork.

    Shows the two MYC tasks, the three matched Dox+ contexts, and the
    no-averaging replicate rule without the sparse text-box layout.
    """
    ax.axis("off")
    sub = task_df[task_df["task_name"].isin(MAIN_TASKS)].copy()
    ax.set_title("Matched MYC classification tasks", fontsize=9, pad=7)

    # Layout in axes coordinates.
    left = 0.02
    task_w = 0.15
    cell_w = 0.205
    total_w = 0.16
    top = 0.80
    header_h = 0.12
    row_h = 0.18
    gap = 0.018

    cols = ["no CPT", "LCPT", "HCPT"]
    x_task = left
    x_cells = [left + task_w + i * cell_w for i in range(3)]
    x_total = left + task_w + 3 * cell_w

    def box(x, y, w, h, text, face="white", weight="normal", fs=6.8, color="black"):
        rect = plt.Rectangle((x, y), w, h, transform=ax.transAxes,
                             facecolor=face, edgecolor="0.55", linewidth=0.55)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, transform=ax.transAxes, ha="center", va="center",
                fontsize=fs, fontweight=weight, color=color)

    # Header row.
    y_header = top
    box(x_task, y_header, task_w, header_h, "task", face="0.92", weight="bold")
    for x, c in zip(x_cells, cols):
        box(x, y_header, cell_w, header_h, c, face="0.92", weight="bold")
    box(x_total, y_header, total_w, header_h, "samples", face="0.92", weight="bold")

    rows = [
        ("MYC@4h", ["D → DT", "D_L → DT_L", "D_H → DT_H"]),
        ("MYC@24h", ["24D → 24DT", "24D_L → 24DT_L", "24D_H → 24DT_H"]),
    ]
    for r, (task, contrasts) in enumerate(rows):
        y = y_header - (r + 1) * row_h
        if task in set(sub["task_name"]):
            row = sub[sub["task_name"] == task].iloc[0]
            n = int(row.get("n_samples", 12))
            npos = int(row.get("n_positive", 6))
            nneg = int(row.get("n_negative", 6))
        else:
            n, npos, nneg = 12, 6, 6
        box(x_task, y, task_w, row_h, task, weight="bold", fs=6.9)
        for x, ctext in zip(x_cells, contrasts):
            box(x, y, cell_w, row_h, ctext, fs=6.2)
            # Small OFF/ON cue underneath each contrast.
            ax.text(x + cell_w/2, y + 0.035, "2 OFF + 2 ON", transform=ax.transAxes,
                    ha="center", va="center", fontsize=5.7, color="0.35")
        box(x_total, y, total_w, row_h, f"n={n}\n{nneg} OFF / {npos} ON", fs=6.3)

    # Footnote / methods cue.
    ax.text(left, 0.20, "Primary CV: leave-one-sample-out (LOSO); robustness CV: leave-one-context-out (LOCO).",
            transform=ax.transAxes, ha="left", va="center", fontsize=7.0)
    ax.text(left, 0.12, "Individual biological replicates are shown in score panels; no averaging.",
            transform=ax.transAxes, ha="left", va="center", fontsize=7.0)

def plot_metric_panel(ax: plt.Axes, metrics: pd.DataFrame, tasks: list[str], title: str, main_cv: str, robust_cv: str,
                      show_pvalues: bool = False, show_metric_values: bool = False) -> None:
    sub = metrics[metrics["task_name"].isin(tasks)].copy()
    xbase = {task: i for i, task in enumerate(tasks)}
    metric_off = {"balanced_accuracy": -0.055, "roc_auc": 0.055}
    cv_off = {main_cv: -0.10, robust_cv: 0.10}
    for _, row in sub.iterrows():
        cv = row["cv_scheme"]
        if cv not in {main_cv, robust_cv}:
            continue
        for met in ["balanced_accuracy", "roc_auc"]:
            if met not in row or pd.isna(row[met]):
                continue
            x = xbase[row["task_name"]] + metric_off[met] + cv_off.get(cv, 0.0)
            val = float(row[met])
            filled = cv == main_cv
            ax.scatter([x], [val], s=50, marker=METRIC_MARKERS[met],
                       facecolors=METRIC_COLORS[met] if filled else "white",
                       edgecolors=METRIC_COLORS[met], linewidths=1.4, zorder=3)
            if show_metric_values:
                yoff = 0.035 if met == "balanced_accuracy" else 0.070
                ax.text(x, min(val + yoff, 1.065), f"{val:.2f}", ha="center", va="bottom", fontsize=6.8)
            if show_pvalues and cv == main_cv and "permutation_pvalue" in row and pd.notna(row["permutation_pvalue"]):
                ax.text(x, min(val + 0.095, 1.08), f"p={row['permutation_pvalue']:.3g}", ha="center", va="bottom", fontsize=6.2)
    ax.axhline(0.5, linestyle="--", color="0.45", linewidth=0.8)
    ax.set_xlim(-0.45, len(tasks) - 0.55)
    ax.set_ylim(0, 1.08)
    ax.set_xticks(range(len(tasks)))
    ax.set_xticklabels(tasks)
    ax.set_ylabel("Performance")
    ax.set_title(title, pad=6)
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    handles = [
        Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=METRIC_COLORS["balanced_accuracy"], markeredgecolor=METRIC_COLORS["balanced_accuracy"], label="Balanced accuracy"),
        Line2D([0], [0], marker="s", linestyle="None", markerfacecolor=METRIC_COLORS["roc_auc"], markeredgecolor=METRIC_COLORS["roc_auc"], label="ROC-AUC"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor="black", label="LOSO"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor="white", label="LOCO"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, ncol=2, columnspacing=1.4, handletextpad=0.5)


def load_predictions(run_dir: Path, task: str, model: str, cv: str) -> pd.DataFrame:
    return read_csv(run_dir / "predictions" / f"oof_predictions_{task}_{model}_{cv}.csv")


def plot_oof_scores(ax: plt.Axes, df: pd.DataFrame, title: str, highlight_samples: list[str], draw_lines: bool = True) -> None:
    df = df.copy()
    df["replicate"] = df.get("replicate", 1).astype(int)
    df["context_sort"] = df["context_value"].map({"None_CPT": 0, "LCPT": 1, "HCPT": 2}).fillna(9)
    df = df.sort_values(["context_sort", "label", "replicate", "sample_id"])
    xpos, ticks, x = {}, [], 0.0
    centers, boundaries = [], []
    for ctx in CONTEXT_ORDER:
        ctx_df = df[df["context_value"] == ctx]
        if ctx_df.empty:
            continue
        start = x
        for lab in [0, 1]:
            lab_df = ctx_df[ctx_df["label"] == lab].sort_values("replicate")
            for _, row in lab_df.iterrows():
                xpos[row["sample_id"]] = x
                ticks.append((x, short_rep_label(row)))
                x += 1.0
        centers.append(((start + x - 1) / 2, CONTEXT_LABEL.get(ctx, ctx)))
        boundaries.append(x - 0.5)
        x += 0.65
    if draw_lines:
        for (_, ctx, lab), sub in df.groupby(["task_name", "context_value", "label"]):
            sub = sub.sort_values("replicate")
            xs = [xpos[s] for s in sub["sample_id"] if s in xpos]
            ys = sub["pred_score"].tolist()
            if len(xs) >= 2:
                ax.plot(xs, ys, color="0.72", linewidth=0.9, zorder=1)
    for _, row in df.iterrows():
        sid = row["sample_id"]
        ctx = row["context_value"]
        lab = int(row["label"])
        edge = "black"; lw = 0.5; size = 52
        if sid in highlight_samples:
            edge = "#7f0000"; lw = 1.5; size = 68
        ax.scatter([xpos[sid]], [row["pred_score"]], marker=CONTEXT_MARKERS.get(ctx, "o"),
                   color=LABEL_COLORS[lab], edgecolor=edge, linewidth=lw, s=size, zorder=3)
        if sid in highlight_samples:
            rep_txt = "R1" if str(sid).endswith("RP1") else "R2"
            # Keep the callout short; the x-axis already identifies OFF/ON and context.
            dx = -0.26 if rep_txt == "R1" else 0.10
            dy = 0.055 if rep_txt == "R2" else 0.025
            ax.text(xpos[sid] + dx, row["pred_score"] + dy, rep_txt, rotation=0,
                    ha="center", va="bottom", fontsize=6.8, fontweight="bold", color="#7f0000")
    ax.axhline(0.5, linestyle="--", color="0.45", linewidth=0.8)
    for b in boundaries[:-1]:
        ax.axvline(b + 0.325, color="0.82", linewidth=0.8)
    for c, lab in centers:
        ax.text(c, 1.05, lab, ha="center", va="bottom", fontsize=7)
    ax.set_ylim(-0.02, 1.08)
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], rotation=60, ha="right")
    ax.set_ylabel("Out-of-fold P(MYC ON)")
    ax.set_title(title, pad=6)
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)


def plot_gene_counts(ax: plt.Axes, gene_dir: Path, scale: str) -> None:
    df = read_csv(gene_dir / "gene_screen_task_summary.csv")
    df["task_name"] = pd.Categorical(df["task_name"], TASK_ORDER, ordered=True)
    df = df.sort_values("task_name")
    vals = df["n_sig_fdr"].astype(int).to_numpy()
    labels = df["task_name"].astype(str).tolist()
    ax.bar(range(len(vals)), vals, color="0.35")
    if scale == "log":
        ax.set_yscale("symlog", linthresh=10)
    ymax = max(vals.max() if len(vals) else 1, 1)
    for i, v in enumerate(vals):
        ax.text(i, v + ymax * 0.025, f"{v}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Significant genes")
    ax.set_title("Gene-screen yield (FDR <= 0.25)")
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)


def plot_region(ax: plt.Axes, gene_dir: Path) -> None:
    """Plot baseline-region composition among significant genes.

    Tasks with zero significant genes have undefined composition.  Earlier
    versions silently omitted those bars, which made Supplementary Fig. 2C
    look incomplete.  We now draw an explicit neutral placeholder labelled
    "0 sig." / "not defined" so the absence is interpretable.
    """
    task_summary_path = gene_dir / "gene_screen_task_summary.csv"
    task_summary = read_csv(task_summary_path) if task_summary_path.exists() else pd.DataFrame()
    n_sig = {}
    if not task_summary.empty and {"task_name", "n_sig_fdr"}.issubset(task_summary.columns):
        n_sig = dict(zip(task_summary["task_name"].astype(str), task_summary["n_sig_fdr"].fillna(0).astype(int)))

    p = gene_dir / "gene_screen_region_summary.csv"
    if p.exists():
        df = read_csv(p)
        if {"task_name", "region", "fraction_of_sig"}.issubset(df.columns):
            wide = df.pivot_table(index="task_name", columns="region", values="fraction_of_sig", aggfunc="sum").fillna(0)
        else:
            wide = pd.DataFrame(index=TASK_ORDER, columns=["head", "mid", "tail"]).fillna(0)
    elif not task_summary.empty and {"task_name", "head_frac_sig", "mid_frac_sig", "tail_frac_sig"}.issubset(task_summary.columns):
        wide = task_summary.set_index("task_name")[["head_frac_sig", "mid_frac_sig", "tail_frac_sig"]]
        wide.columns = ["head", "mid", "tail"]
        wide = wide.fillna(0)
    else:
        wide = pd.DataFrame(index=TASK_ORDER, columns=["head", "mid", "tail"]).fillna(0)

    wide = wide.reindex(TASK_ORDER).fillna(0)
    xs = np.arange(len(wide))
    bottom = np.zeros(len(wide))
    positive_mask = np.array([n_sig.get(task, int(wide.loc[task].sum() > 0)) > 0 for task in wide.index])

    for region in ["head", "mid", "tail"]:
        vals = wide[region].to_numpy(dtype=float) if region in wide.columns else np.zeros(len(wide))
        vals = np.where(positive_mask, vals, 0.0)
        ax.bar(xs, vals, bottom=bottom, color=REGION_COLORS[region], label=region, linewidth=0.0)
        bottom += vals

    # Explicit placeholders for tasks with no significant genes.
    for i, task in enumerate(wide.index):
        if not positive_mask[i]:
            ax.bar(i, 1.0, width=0.72, bottom=0.0, facecolor="white", edgecolor="0.55",
                   linewidth=0.8, linestyle=(0, (3, 2)), zorder=1)
            ax.text(i, 0.50, "0 sig.\nnot defined", ha="center", va="center",
                    fontsize=6.3, color="0.30", linespacing=0.95)

    ax.set_xticks(xs)
    ax.set_xticklabels(wide.index.tolist())
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction of significant genes")
    ax.set_title("Baseline-region composition\namong significant genes")
    handles = [Patch(facecolor=REGION_COLORS[r], label=r) for r in ["head", "mid", "tail"]]
    if any(~positive_mask):
        handles.append(Patch(facecolor="white", edgecolor="0.55", linestyle=(0, (3, 2)), label="not defined"))
    ax.legend(handles=handles, frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)


def find_matrix(args: argparse.Namespace) -> Optional[Path]:
    if args.matrix and Path(args.matrix).exists():
        return Path(args.matrix)
    for key in ["FIG6_COUNT_DERIVED_MATRIX", "MATRIX_PATH"]:
        val = os.environ.get(key)
        if val and Path(val).exists():
            return Path(val)
    return None


def plot_hcpt_heatmap_vector(ax: plt.Axes, run_dir: Path, matrix_path: Optional[Path]) -> bool:
    if matrix_path is None or not matrix_path.exists():
        return False
    p = run_dir / "gene_screen" / "per_gene_auc_MYC@24h.csv"
    if not p.exists():
        return False
    try:
        markers = read_csv(p)
        if "mean_delta" in markers.columns:
            markers = markers[markers["mean_delta"] > 0].copy()
        markers["fdr_sort"] = markers.get("fdr", pd.Series([1.0] * len(markers))).fillna(1.0)
        markers["max_auc_sort"] = markers.get("max_auc", markers.get("auc", pd.Series([0.5] * len(markers))))
        markers["abs_delta"] = markers.get("mean_delta", pd.Series([0.0] * len(markers))).abs()
        markers = markers.sort_values(["max_auc_sort", "fdr_sort", "abs_delta"], ascending=[False, True, False]).head(20)
        genes = markers["gene"].astype(str).tolist()
        mat = pd.read_csv(matrix_path, sep="\t" if matrix_path.suffix.lower() in {".tsv", ".txt"} else ",")
        gene_col = mat.columns[0]
        samples = ["29_24D_H_CPT_RP1", "31_24D_H_CPT_RP2", "30_24DT_H_CPT_RP1", "32_24DT_H_CPT_RP2"]
        if any(s not in mat.columns for s in samples):
            return False
        sub = mat[mat[gene_col].astype(str).isin(genes)].copy()
        if sub.empty:
            return False
        sub["_order"] = sub[gene_col].astype(str).map({g: i for i, g in enumerate(genes)})
        sub = sub.sort_values("_order")
        vals = sub[samples].astype(float).to_numpy()
        sd = vals.std(axis=1, keepdims=True); sd[sd == 0] = 1
        z = (vals - vals.mean(axis=1, keepdims=True)) / sd
        labels = [symbol_from_gene(x) for x in sub[gene_col].astype(str)]
        im = ax.imshow(z, aspect="auto", cmap="viridis", vmin=-1.6, vmax=1.6)
        ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=5.8)
        ax.set_xticks(range(len(samples))); ax.set_xticklabels([clean_sample_id(s) for s in samples], rotation=60, ha="right")
        ax.axvline(1.5, color="black", linewidth=0.8)
        ax.set_title("24 h HCPT focus")
        cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label("Row z-score", fontsize=7)
        return True
    except Exception as e:
        print(f"[WARN] vector heatmap failed: {e}")
        return False


def plot_hcpt_local(ax: plt.Axes, run_dir: Path, matrix_path: Optional[Path]) -> None:
    if plot_hcpt_heatmap_vector(ax, run_dir, matrix_path):
        return
    img = run_dir / "gene_screen" / "figures" / "gene_screen_hcpt24h_top_myc_heatmap.png"
    if img.exists():
        ax.imshow(mpimg.imread(img)); ax.axis("off"); ax.set_title("24 h HCPT focus")
    else:
        ax.axis("off"); ax.text(0.5, 0.5, "HCPT focus panel unavailable", ha="center", va="center")


def compose_main(run_dir: Path, out_dir: Path, args: argparse.Namespace) -> list[Path]:
    task_summary = read_csv(run_dir / "summary" / "task_membership_summary.csv")
    metrics = read_csv(run_dir / "summary" / "metrics_all_models.csv")
    metrics = metrics[(metrics["model_name"] == args.model) & (metrics["task_name"].isin(MAIN_TASKS))]
    pred4 = load_predictions(run_dir, "MYC@4h", args.model, args.main_cv)
    pred24 = load_predictions(run_dir, "MYC@24h", args.model, args.main_cv)

    # Manual axes avoid constrained-layout failures caused by long labels and rotated ticks.
    fig = plt.figure(figsize=(9.0, 7.0))
    axes = [
        fig.add_axes([0.06, 0.58, 0.42, 0.34]),
        fig.add_axes([0.58, 0.58, 0.36, 0.34]),
        fig.add_axes([0.06, 0.10, 0.42, 0.40]),
        fig.add_axes([0.58, 0.10, 0.36, 0.40]),
    ]
    plot_task_panel(axes[0], task_summary)
    plot_metric_panel(axes[1], metrics, MAIN_TASKS, "MYC task metrics", args.main_cv, args.robust_cv, args.show_pvalues, args.show_metric_values)
    plot_oof_scores(axes[2], pred4, "MYC@4h - held-out sample scores", [], not args.no_replicate_lines)
    plot_oof_scores(axes[3], pred24, "MYC@24h - held-out sample scores", args.highlight_samples, not args.no_replicate_lines)
    for ax, letter in zip(axes, "ABCD"):
        panel_label(ax, letter)
    if args.with_title:
        fig.suptitle("Figure 6 | Sample-level classifier readouts provide a secondary view of MYC-state separation", fontsize=10.5, y=0.985)
    return save_all(fig, out_dir, "Figure6_main_publication_ready", args.formats, args.dpi)


def compose_support(run_dir: Path, out_dir: Path, args: argparse.Namespace, matrix_path: Optional[Path]) -> list[Path]:
    metrics = read_csv(run_dir / "summary" / "metrics_all_models.csv")
    metrics = metrics[(metrics["model_name"] == args.model) & (metrics["task_name"].isin(CONTROL_TASKS))]
    if args.support_mode == "local":
        fig = plt.figure(figsize=(9.0, 3.8))
        axes = [fig.add_axes([0.06, 0.18, 0.42, 0.70]), fig.add_axes([0.58, 0.18, 0.36, 0.70])]
        plot_metric_panel(axes[0], metrics, CONTROL_TASKS, "Dox-only control metrics", args.main_cv, args.robust_cv, args.show_pvalues, args.show_metric_values)
        plot_hcpt_local(axes[1], run_dir, matrix_path)
        letters = "AB"; stem = "FigureS2_local_support_publication_ready"
    else:
        fig = plt.figure(figsize=(9.0, 7.0))
        axes = [
            fig.add_axes([0.06, 0.58, 0.42, 0.34]),
            fig.add_axes([0.58, 0.58, 0.36, 0.34]),
            fig.add_axes([0.06, 0.10, 0.42, 0.38]),
            fig.add_axes([0.58, 0.10, 0.36, 0.38]),
        ]
        plot_metric_panel(axes[0], metrics, CONTROL_TASKS, "Dox-only control metrics", args.main_cv, args.robust_cv, args.show_pvalues, args.show_metric_values)
        plot_gene_counts(axes[1], run_dir / "gene_screen", args.gene_count_scale)
        plot_region(axes[2], run_dir / "gene_screen")
        plot_hcpt_local(axes[3], run_dir, matrix_path)
        letters = "ABCD"; stem = "FigureS2_dox_gene_support_publication_ready"
    for ax, letter in zip(axes, letters):
        panel_label(ax, letter)
    if args.with_title:
        fig.suptitle("Supplementary Fig. 2 | Doxycycline priming and gene-level support for Fig. 6", fontsize=10.5, y=0.985)
    return save_all(fig, out_dir, stem, args.formats, args.dpi)


def compose_s3(run_dir: Path, out_dir: Path, args: argparse.Namespace) -> list[Path]:
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.5), constrained_layout=True)
    imgs = [
        ("A", run_dir / "gene_screen" / "figures" / "gene_screen_auc80_counts.png"),
        ("B", run_dir / "gene_screen" / "figures" / "gene_screen_myc24_dox24_sig_overlap.png"),
        ("C", run_dir / "gene_screen" / "figures" / "gene_screen_ddr_overlap_summary.png"),
    ]
    for ax, (lab, path) in zip(axes, imgs):
        if path.exists():
            ax.imshow(mpimg.imread(path)); ax.axis("off")
        else:
            ax.axis("off"); ax.text(0.5, 0.5, f"Missing {path.name}", ha="center", va="center")
        panel_label(ax, lab)
    return save_all(fig, out_dir, "FigureS3_optional_gene_diagnostics_not_default", args.formats, args.dpi)


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    out_dir = Path(args.output_dir).resolve() if args.output_dir else run_dir / "publication_figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = Path(args.matrix).resolve() if args.matrix and Path(args.matrix).exists() else None
    if matrix_path is None:
        for key in ["FIG6_COUNT_DERIVED_MATRIX", "MATRIX_PATH"]:
            val = os.environ.get(key)
            if val and Path(val).exists():
                matrix_path = Path(val).resolve(); break
    if matrix_path:
        print(f"[INFO] Using matrix for vector HCPT heatmap: {matrix_path}")
    else:
        print("[INFO] No matrix supplied/found; support heatmap will use existing PNG if available.")
    written = []
    written += compose_main(run_dir, out_dir, args)
    if args.support_mode != "none":
        written += compose_support(run_dir, out_dir, args, matrix_path)
    if args.make_s3:
        written += compose_s3(run_dir, out_dir, args)
    write_publication_manifest(out_dir, args, written)
    for w in written:
        print(f"[OK] Wrote {w}")
    print(f"[OK] Wrote {out_dir / 'publication_figure6_manifest.json'}")
    print("[DONE] Publication-ready figures composed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
