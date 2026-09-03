#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

METHODS = ["DESeq2", "edgeR", "limma/voom", "RankProd", "PenDA", "RankCompV3"]
KINDS = ["All", "Up", "Down"]


def matrix_from_long(df: pd.DataFrame, value_col: str) -> np.ndarray:
    pivot = df.pivot(index="method2", columns="method1", values=value_col)
    return pivot.reindex(index=METHODS[::-1], columns=METHODS).to_numpy(dtype=float)


def annotate(ax, matrix: np.ndarray, mode: str) -> None:
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if not np.isfinite(value):
                continue
            if mode == "jaccard":
                label = "<0.01" if 0 < value < 0.01 else f"{value:.2f}"
                color = "white" if value > 0.55 else "black"
            else:
                label = f"{value:.2f}"
                color = "white" if value > 0.55 or value < -0.55 else "black"
            ax.text(j, i, label, ha="center", va="center", fontsize=7.4, color=color)


def format_heatmap(ax, title: str, show_y: bool = True) -> None:
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=7)
    ax.set_xticks(np.arange(len(METHODS)))
    ax.set_xticklabels(METHODS, rotation=45, ha="right", fontsize=7.7)
    ax.set_yticks(np.arange(len(METHODS)))
    ax.set_yticklabels(METHODS[::-1] if show_y else [], fontsize=7.7)
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)


def main() -> int:
    parser = argparse.ArgumentParser(description="Render Supplementary Figure S5 from frozen source tables.")
    parser.add_argument("--jaccard", required=True, type=Path)
    parser.add_argument("--spearman", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--stem", default="FigureS5")
    args = parser.parse_args()

    jac = pd.read_csv(args.jaccard, sep="\t")
    cor = pd.read_csv(args.spearman, sep="\t")
    required_jac = {"method1", "method2", "value", "kind"}
    required_cor = {"method1", "method2", "value"}
    if not required_jac.issubset(jac.columns):
        raise SystemExit(f"Jaccard table missing columns: {sorted(required_jac - set(jac.columns))}")
    if not required_cor.issubset(cor.columns):
        raise SystemExit(f"Spearman table missing columns: {sorted(required_cor - set(cor.columns))}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(15.0, 5.7), constrained_layout=False)
    outer = GridSpec(1, 2, figure=fig, width_ratios=[3.75, 1.25], wspace=0.32,
                     left=0.045, right=0.965, bottom=0.20, top=0.84)
    left = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0], wspace=0.18)

    jaccard_image = None
    left_axes = []
    for idx, kind in enumerate(KINDS):
        ax = fig.add_subplot(left[idx])
        left_axes.append(ax)
        mat = matrix_from_long(jac[jac["kind"] == kind], "value")
        jaccard_image = ax.imshow(mat, vmin=0, vmax=1, cmap="Blues", aspect="equal")
        annotate(ax, mat, "jaccard")
        format_heatmap(ax, kind, show_y=(idx == 0))

    ax_cor = fig.add_subplot(outer[1])
    cor_mat = matrix_from_long(cor, "value")
    cor_image = ax_cor.imshow(cor_mat, vmin=-1, vmax=1, cmap="RdBu", aspect="equal")
    annotate(ax_cor, cor_mat, "spearman")
    format_heatmap(ax_cor, "Pairwise Spearman concordance", show_y=True)

    fig.text(0.045, 0.925, "A", fontsize=18, fontweight="bold", va="top")
    fig.text(0.073, 0.915, "Pairwise Jaccard overlap", fontsize=13.5, fontweight="bold", va="top")
    fig.text(0.765, 0.925, "B", fontsize=18, fontweight="bold", va="top")

    cb1 = fig.colorbar(jaccard_image, ax=left_axes, fraction=0.018, pad=0.018, shrink=0.90)
    cb1.set_label("Jaccard", fontsize=9)
    cb1.ax.tick_params(labelsize=8)

    cb2 = fig.colorbar(cor_image, ax=ax_cor, fraction=0.050, pad=0.035, shrink=0.90)
    cb2.set_label(r"Spearman $\rho$", fontsize=9)
    cb2.ax.tick_params(labelsize=8)

    for ext in ("pdf", "svg", "png"):
        path = args.output_dir / f"{args.stem}.{ext}"
        kwargs = {"bbox_inches": "tight", "facecolor": "white"}
        if ext == "png":
            kwargs["dpi"] = 600
        fig.savefig(path, **kwargs)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
