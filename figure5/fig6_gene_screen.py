from __future__ import annotations

# SCRIPT_DIR_BOOTSTRAP
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import argparse
from typing import Dict, Set, Tuple

import re

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.metrics import roc_auc_score

from fig6_ml_core import (
    DEFAULT_DDR_GENES,
    benjamini_hochberg,
    build_manifest,
    build_task_table,
    compute_baseline_expression_regions,
    default_task_specs,
    expression_to_sample_by_gene,
    read_sample_sheet,
    read_wide_matrix,
    save_dataframe,
    validate_and_align_inputs,
    write_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Separate gene-level support analysis for Figure 6. "
            "This keeps per-gene AUROC screens distinct from the sample-level classifier pipeline."
        )
    )
    parser.add_argument("--matrix", required=True, help="Wide gene-by-sample matrix (CSV/TSV).")
    parser.add_argument("--sample-sheet", required=True, help="Sample sheet CSV/TSV matching the matrix columns.")
    parser.add_argument("--output-dir", default="fig6_gene_screen", help="Output directory.")
    parser.add_argument("--id-col", default=None, help="Gene identifier column. Defaults to the first column.")
    parser.add_argument("--log1p", action="store_true", help="Apply log1p to the matrix after loading. Use only for TPM/raw abundance matrices, not pre-transformed count-derived matrices.")
    parser.add_argument("--matrix-scale", "--input-scale", dest="matrix_scale", default="unspecified", help="Audit label for matrix scale, e.g. log2_size_factor_normalized_rsem_expected_counts or tpm_log1p. Does not change numeric values except for --log1p.")
    parser.add_argument(
        "--duplicate-policy",
        default="error",
        choices=["error", "mean"],
        help="How to handle duplicated gene identifiers in the input matrix.",
    )
    parser.add_argument(
        "--fill-missing-value",
        type=float,
        default=None,
        help="Optional numeric fill value for missing matrix entries after parsing.",
    )
    parser.add_argument(
        "--min-per-class",
        type=int,
        default=2,
        help="Minimum samples per class required for per-gene AUROC / p-value computation.",
    )
    parser.add_argument(
        "--candidate-low-task",
        default="DOX@4h",
        help="First task for the distortion-candidate comparison.",
    )
    parser.add_argument(
        "--candidate-high-task",
        default="DOX@24h",
        help="Second task for the distortion-candidate comparison.",
    )
    parser.add_argument("--low-auc-threshold", type=float, default=0.60)
    parser.add_argument("--high-auc-threshold", type=float, default=0.70)
    parser.add_argument("--head-pct", type=float, default=10.0)
    parser.add_argument("--tail-pct", type=float, default=10.0)
    parser.add_argument(
        "--ddr-genes",
        default=None,
        help="Optional newline-delimited gene-set file for frozen DDR overlap tests.",
    )
    parser.add_argument(
        "--skip-distortion-candidates",
        action="store_true",
        help="Only write the per-gene AUROC tables and skip pairwise candidate logic.",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.25,
        help="Task-level gene-screen significance threshold used by summary figures and text.",
    )
    parser.add_argument(
        "--strong-auc-threshold",
        type=float,
        default=0.80,
        help="Threshold on max(AUC, 1-AUC) for strong univariate separation summaries.",
    )
    parser.add_argument(
        "--focus-task",
        default="MYC@24h",
        help="Task used for the focused context heatmap.",
    )
    parser.add_argument(
        "--focus-context",
        default="HCPT",
        help="Context value used for the focused context heatmap.",
    )
    parser.add_argument(
        "--focus-top-n",
        type=int,
        default=20,
        help="Number of top genes to show in the focused context heatmap.",
    )
    parser.add_argument(
        "--skip-summary-figures",
        action="store_true",
        help="Skip writing the automatic gene-screen summary figures and text.",
    )
    return parser.parse_args()


TASK_ORDER = ["MYC@4h", "MYC@24h", "DOX@4h", "DOX@24h"]
REGION_ORDER = ["head", "mid", "tail"]
TASK_LABELS = {
    "MYC@4h": "MYC@4h",
    "MYC@24h": "MYC@24h",
    "DOX@4h": "DOX@4h",
    "DOX@24h": "DOX@24h",
}


def load_gene_set(path: str | None) -> Set[str]:
    if path is None:
        return {g.upper() for g in DEFAULT_DDR_GENES}
    genes = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.add(gene.upper())
    return genes


GENE_SYMBOL_SPLIT_PATTERN = re.compile(r"[|_]")


def normalize_gene_symbol(gene_id: str) -> str:
    """Extract a gene-symbol-like token from composite IDs such as ENSG..._TP53."""
    token = str(gene_id).strip()
    if not token:
        return token
    parts = [p for p in GENE_SYMBOL_SPLIT_PATTERN.split(token) if p]
    if len(parts) >= 2:
        tail = parts[-1].strip()
        if tail and not tail.upper().startswith("ENSG"):
            return tail
    return token


def _task_sort_key(task_name: str) -> int:
    return TASK_ORDER.index(task_name) if task_name in TASK_ORDER else len(TASK_ORDER)


def _max_auc(series: pd.Series) -> pd.Series:
    arr = series.to_numpy(dtype=float)
    return pd.Series(np.maximum(arr, 1.0 - arr), index=series.index)


def per_gene_screen(
    X_task: pd.DataFrame,
    task_meta: pd.DataFrame,
    positive_label: str,
    min_per_class: int,
) -> pd.DataFrame:
    y = (task_meta["label_text"] == positive_label).to_numpy(dtype=bool)
    n_pos = int(y.sum())
    n_neg = int((~y).sum())
    if min(n_pos, n_neg) < int(min_per_class):
        raise ValueError(
            f"Not enough samples per class for per-gene screen: pos={n_pos}, neg={n_neg}."
        )

    rows = []
    for gene in X_task.columns:
        values = X_task[gene].to_numpy(dtype=float)
        pos = values[y]
        neg = values[~y]
        try:
            auc = roc_auc_score(y.astype(int), values)
        except ValueError:
            continue
        try:
            pval = mannwhitneyu(pos, neg, alternative="two-sided").pvalue
        except ValueError:
            pval = 1.0
        rows.append(
            {
                "gene": str(gene),
                "gene_symbol": normalize_gene_symbol(str(gene)),
                "auc": float(auc),
                "max_auc": float(max(auc, 1.0 - auc)),
                "mean_positive": float(np.mean(pos)),
                "mean_negative": float(np.mean(neg)),
                "mean_delta": float(np.mean(pos) - np.mean(neg)),
                "pvalue": float(pval),
                "n_positive": n_pos,
                "n_negative": n_neg,
            }
        )

    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["fdr"] = benjamini_hochberg(res["pvalue"].to_numpy(dtype=float))
    res = res.sort_values(["fdr", "max_auc", "pvalue"], ascending=[True, False, True]).reset_index(drop=True)
    return res


def summarize_task_level_counts(
    per_task_results: Dict[str, pd.DataFrame],
    baseline_regions: pd.DataFrame,
    fdr_threshold: float,
    strong_auc_threshold: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    base = baseline_regions[["gene", "region"]].copy()
    task_rows = []
    region_rows = []
    for task_name, res in per_task_results.items():
        merged = res.merge(base, on="gene", how="left") if not res.empty else res.copy()
        sig = merged.loc[merged["fdr"] <= float(fdr_threshold)].copy() if not merged.empty else merged.copy()
        strong = merged.loc[merged["max_auc"] >= float(strong_auc_threshold)].copy() if not merged.empty else merged.copy()
        task_rows.append(
            {
                "task_name": task_name,
                "n_genes_tested": int(len(merged)),
                "n_sig_fdr": int(len(sig)),
                "n_strong_auc": int(len(strong)),
                "n_auc_ge_0p8": int((merged["auc"] >= 0.8).sum()) if not merged.empty else 0,
                "n_auc_le_0p2": int((merged["auc"] <= 0.2).sum()) if not merged.empty else 0,
                "head_frac_sig": float((sig["region"] == "head").mean()) if len(sig) else np.nan,
                "mid_frac_sig": float((sig["region"] == "mid").mean()) if len(sig) else np.nan,
                "tail_frac_sig": float((sig["region"] == "tail").mean()) if len(sig) else np.nan,
            }
        )
        for region in REGION_ORDER:
            n_region = int((sig["region"] == region).sum()) if len(sig) else 0
            region_rows.append(
                {
                    "task_name": task_name,
                    "region": region,
                    "n_sig_fdr": n_region,
                    "fraction_of_sig": float(n_region / len(sig)) if len(sig) else 0.0,
                    "fdr_threshold": float(fdr_threshold),
                }
            )
    task_df = pd.DataFrame(task_rows).sort_values("task_name", key=lambda s: s.map(_task_sort_key)).reset_index(drop=True)
    region_df = pd.DataFrame(region_rows).sort_values(["task_name", "region"], key=lambda s: s.map(lambda x: _task_sort_key(x) if x in TASK_ORDER else REGION_ORDER.index(x) if x in REGION_ORDER else 99)).reset_index(drop=True)
    return task_df, region_df


def _make_bar_plot(ax, df: pd.DataFrame, value_col: str, title: str, ylabel: str, annotate_fmt: str = "{:.0f}") -> None:
    task_names = df["task_name"].tolist()
    x = np.arange(len(task_names))
    vals = df[value_col].to_numpy(dtype=float)
    bars = ax.bar(x, vals)
    ax.set_xticks(x)
    ax.set_xticklabels(task_names)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.25)
    ymax = max(vals.max() if len(vals) else 0.0, 1.0)
    ax.set_ylim(0, ymax * 1.18)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2.0, b.get_height() + ymax * 0.03, annotate_fmt.format(v), ha="center", va="bottom", fontsize=9)


def plot_task_counts(task_summary: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    _make_bar_plot(
        ax,
        task_summary,
        "n_sig_fdr",
        "Gene-screen counts at FDR threshold",
        "Significant genes",
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_strong_auc_counts(task_summary: pd.DataFrame, output_path: Path, auc_threshold: float) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    _make_bar_plot(
        ax,
        task_summary,
        "n_strong_auc",
        f"Strong per-gene separation (max AUC ≥ {auc_threshold:.2f})",
        "Genes",
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_region_composition(region_summary: pd.DataFrame, output_path: Path) -> None:
    """Plot baseline-region composition and explicitly mark tasks with zero significant genes.

    When a task has no FDR-significant genes, region composition is undefined rather
    than zero.  Earlier diagnostic plots silently omitted such tasks, which made
    DOX@4h look like a plotting error.  This function keeps the x-axis task
    present and draws a neutral dashed placeholder labeled "0 sig.; not defined".
    """
    tasks = [t for t in TASK_ORDER if t in set(region_summary["task_name"].astype(str))]
    extra = [t for t in region_summary["task_name"].drop_duplicates().astype(str).tolist() if t not in tasks]
    tasks.extend(extra)
    x = np.arange(len(tasks))
    bottom = np.zeros(len(tasks), dtype=float)
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)

    total_sig_by_task = {
        task: int(region_summary.loc[region_summary["task_name"].astype(str) == task, "n_sig_fdr"].sum())
        for task in tasks
    }

    for region in REGION_ORDER:
        vals = []
        for task in tasks:
            row = region_summary.loc[
                (region_summary["task_name"].astype(str) == task)
                & (region_summary["region"].astype(str) == region)
            ]
            vals.append(float(row["fraction_of_sig"].iloc[0]) if not row.empty else 0.0)
        vals = np.asarray(vals, dtype=float)
        ax.bar(x, vals, bottom=bottom, label=region)
        bottom += vals

    # Explicitly show undefined composition for zero-hit tasks.
    for i, task in enumerate(tasks):
        if total_sig_by_task.get(task, 0) == 0:
            ax.bar(
                i,
                1.0,
                bottom=0.0,
                width=0.70,
                facecolor="none",
                edgecolor="0.55",
                linestyle="--",
                linewidth=1.1,
                label="not defined" if i == next((j for j, t in enumerate(tasks) if total_sig_by_task.get(t, 0) == 0), i) else None,
            )
            ax.text(
                i,
                0.50,
                "0 sig.\nnot defined",
                ha="center",
                va="center",
                fontsize=8,
                color="0.35",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(tasks)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction of significant genes")
    ax.set_title("Baseline-region composition among significant genes")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_focus_heatmap(
    X_task: pd.DataFrame,
    task_meta: pd.DataFrame,
    per_gene_df: pd.DataFrame,
    focus_task: str,
    focus_context: str,
    output_path: Path,
    summary_path: Path,
    top_n: int,
) -> None:
    focus_meta = task_meta.loc[task_meta["context_value"].astype(str) == str(focus_context)].copy()
    if focus_meta.empty:
        return
    ranked = per_gene_df.copy()
    ranked = ranked.loc[ranked["auc"] > 0.5].copy()
    ranked = ranked.sort_values(["fdr", "auc", "max_auc"], ascending=[True, False, False]).reset_index(drop=True)
    if ranked.empty:
        return
    gene_list = ranked["gene"].head(int(top_n)).tolist()
    gene_symbols = ranked.set_index("gene").loc[gene_list, "gene_symbol"].tolist()
    matrix = X_task.loc[focus_meta["sample_id"], gene_list].transpose().copy()
    if matrix.empty:
        return
    values = matrix.to_numpy(dtype=float)
    row_means = values.mean(axis=1, keepdims=True)
    row_sds = values.std(axis=1, keepdims=True)
    row_sds[row_sds == 0] = 1.0
    z = (values - row_means) / row_sds

    ordered_meta = focus_meta.sort_values(["label", "sample_id"], ascending=[True, True]).reset_index(drop=True)
    ordered_ids = ordered_meta["sample_id"].tolist()
    z_df = pd.DataFrame(z, index=gene_symbols, columns=matrix.columns)
    z_df = z_df[ordered_ids]

    fig, ax = plt.subplots(figsize=(max(6, len(ordered_ids) * 0.8), max(6, len(gene_symbols) * 0.32)), constrained_layout=True)
    im = ax.imshow(z_df.to_numpy(dtype=float), aspect="auto", interpolation="nearest")
    ax.set_title(f"{focus_task} {focus_context} top-marker heatmap")
    ax.set_xticks(np.arange(len(ordered_ids)))
    ax.set_xticklabels(ordered_ids, rotation=60, ha="right")
    ax.set_yticks(np.arange(len(gene_symbols)))
    ax.set_yticklabels(gene_symbols, fontsize=8)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Top per-gene markers")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Row z-score")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary_rows = []
    for sid in ordered_ids:
        row = ordered_meta.loc[ordered_meta["sample_id"] == sid].iloc[0]
        sample_vals = X_task.loc[sid, gene_list].to_numpy(dtype=float)
        summary_rows.append(
            {
                "sample_id": sid,
                "group_label": row.get("group_label", sid),
                "label_text": row.get("label_text", ""),
                "context_value": row.get("context_value", ""),
                "mean_top_marker_expression": float(np.mean(sample_vals)),
            }
        )
    save_dataframe(pd.DataFrame(summary_rows), summary_path)


def overlap_summary(sig_a: pd.DataFrame, sig_b: pd.DataFrame, task_a: str, task_b: str) -> pd.DataFrame:
    set_a = set(sig_a["gene_symbol"].str.upper())
    set_b = set(sig_b["gene_symbol"].str.upper())
    shared = set_a & set_b
    return pd.DataFrame([
        {
            "task_a": task_a,
            "task_b": task_b,
            "n_task_a_sig": len(set_a),
            "n_task_b_sig": len(set_b),
            "n_shared_sig": len(shared),
            "n_task_a_only": len(set_a - set_b),
            "n_task_b_only": len(set_b - set_a),
            "jaccard": float(len(shared) / len(set_a | set_b)) if (set_a | set_b) else np.nan,
        }
    ])


def plot_overlap_summary(overlap_df: pd.DataFrame, output_path: Path) -> None:
    row = overlap_df.iloc[0]
    labels = [f"{row['task_a']} only", "shared", f"{row['task_b']} only"]
    values = [int(row['n_task_a_only']), int(row['n_shared_sig']), int(row['n_task_b_only'])]
    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    bars = ax.bar(np.arange(3), values)
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("Significant genes")
    ax.set_title(f"Significant-hit overlap: {row['task_a']} vs {row['task_b']}")
    ax.grid(axis="y", alpha=0.25)
    ymax = max(values + [1])
    ax.set_ylim(0, ymax * 1.18)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width()/2.0, b.get_height() + ymax*0.03, str(v), ha='center', va='bottom', fontsize=9)
    ax.text(0.98, 0.98, f"Jaccard = {float(row['jaccard']):.2f}", ha="right", va="top", transform=ax.transAxes)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_ddr_summary(overlap_df: pd.DataFrame, output_path: Path) -> None:
    row = overlap_df.iloc[0]
    fig, ax = plt.subplots(figsize=(6.5, 3.8), constrained_layout=True)
    values = [int(row["n_candidate_genes"]), int(row["n_ddr_genes_in_candidates"])]
    labels = ["candidates", "DDR in candidates"]
    bars = ax.bar(np.arange(2), values)
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Genes")
    ax.set_title(f"DDR overlap: {row['low_task']} → {row['high_task']}")
    ax.grid(axis="y", alpha=0.25)
    ymax = max(values + [1])
    ax.set_ylim(0, ymax * 1.25)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width()/2.0, b.get_height() + ymax*0.04, str(v), ha='center', va='bottom', fontsize=9)
    ax.text(
        0.98,
        0.98,
        f"OR = {float(row['fisher_odds_ratio']):.2f}\np = {float(row['fisher_pvalue']):.3g}",
        ha="right",
        va="top",
        transform=ax.transAxes,
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_gene_screen_text(
    task_summary: pd.DataFrame,
    region_summary: pd.DataFrame,
    fdr_threshold: float,
    strong_auc_threshold: float,
) -> Tuple[str, str]:
    by_task = task_summary.set_index("task_name")
    base = (
        f"Univariate per-gene screens were analysed separately from the sample-level classifiers. At FDR ≤ {fdr_threshold:.2f}, "
        f"MYC@4h yielded {int(by_task.loc['MYC@4h', 'n_sig_fdr']) if 'MYC@4h' in by_task.index else 0} significant genes, "
        f"MYC@24h yielded {int(by_task.loc['MYC@24h', 'n_sig_fdr']) if 'MYC@24h' in by_task.index else 0}, "
        f"DOX@4h yielded {int(by_task.loc['DOX@4h', 'n_sig_fdr']) if 'DOX@4h' in by_task.index else 0}, and "
        f"DOX@24h yielded {int(by_task.loc['DOX@24h', 'n_sig_fdr']) if 'DOX@24h' in by_task.index else 0}. "
        f"Using a strong-separation threshold of max(AUC, 1−AUC) ≥ {strong_auc_threshold:.2f}, the corresponding counts were "
        f"{int(by_task.loc['MYC@4h', 'n_strong_auc']) if 'MYC@4h' in by_task.index else 0}, "
        f"{int(by_task.loc['MYC@24h', 'n_strong_auc']) if 'MYC@24h' in by_task.index else 0}, "
        f"{int(by_task.loc['DOX@4h', 'n_strong_auc']) if 'DOX@4h' in by_task.index else 0}, and "
        f"{int(by_task.loc['DOX@24h', 'n_strong_auc']) if 'DOX@24h' in by_task.index else 0}."
    )
    dox_bridge = ""
    if {'MYC@24h', 'DOX@24h'}.issubset(by_task.index):
        myc_sig = int(by_task.loc['MYC@24h', 'n_sig_fdr'])
        dox_sig = int(by_task.loc['DOX@24h', 'n_sig_fdr'])
        myc_head = float(by_task.loc['MYC@24h', 'head_frac_sig']) if pd.notna(by_task.loc['MYC@24h', 'head_frac_sig']) else np.nan
        dox_head = float(by_task.loc['DOX@24h', 'head_frac_sig']) if pd.notna(by_task.loc['DOX@24h', 'head_frac_sig']) else np.nan
        bridge = (
            f"At 24 h, the Dox-only task remained transcriptionally detectable ({dox_sig} genes at FDR ≤ {fdr_threshold:.2f}), "
            f"confirming that priming is not chemically silent. However, the late MYC task remained more head-concentrated "
            f"({myc_head:.0%} of significant MYC@24h genes in the baseline head versus {dox_head:.0%} for DOX@24h)"
            if np.isfinite(myc_head) and np.isfinite(dox_head)
            else "At 24 h, the Dox-only task remained transcriptionally detectable, confirming that priming is not chemically silent."
        )
        if np.isfinite(myc_head) and np.isfinite(dox_head):
            if dox_head < myc_head:
                bridge += ", consistent with a measurable but differently structured Dox-only response that does not reproduce the head-focused late MYC pattern."
            else:
                bridge += ", but the late MYC and late Dox tasks differed more in structure than in raw detectability."
        elif dox_sig >= myc_sig:
            bridge += " The Dox-only signature was broad, but breadth alone did not recreate the late MYC distortion geometry."
        dox_bridge = bridge
    return base, dox_bridge


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    text_dir = output_dir / "text"
    figures_dir.mkdir(parents=True, exist_ok=True)
    text_dir.mkdir(parents=True, exist_ok=True)

    expr = read_wide_matrix(
        matrix_path=args.matrix,
        id_col=args.id_col,
        apply_log1p=args.log1p,
        duplicate_policy=args.duplicate_policy,
        fill_missing_value=args.fill_missing_value,
    )
    meta = read_sample_sheet(args.sample_sheet)
    expr, meta = validate_and_align_inputs(expr, meta)
    X_samples = expression_to_sample_by_gene(expr)

    task_specs = default_task_specs(include_dox_controls=True)
    task_lookup = {task.name: task for task in task_specs}

    per_task_results: Dict[str, pd.DataFrame] = {}
    task_tables: Dict[str, Tuple[pd.DataFrame, pd.DataFrame]] = {}
    task_membership_rows = []
    for task_name, task in task_lookup.items():
        X_task, task_meta = build_task_table(X_samples, meta, task)
        task_tables[task_name] = (X_task, task_meta)
        save_dataframe(task_meta, output_dir / f"task_samples_{task_name}.csv")
        task_membership_rows.append(
            {
                "task_name": task.name,
                "time_h": task.time_h,
                "label_col": task.label_col,
                "positive_label": task.positive_label,
                "negative_label": task.negative_label,
                "context_col": task.context_col,
                "n_samples": int(task_meta.shape[0]),
                "n_positive": int(task_meta["label"].sum()),
                "n_negative": int((1 - task_meta["label"]).sum()),
            }
        )
        res = per_gene_screen(
            X_task=X_task,
            task_meta=task_meta,
            positive_label=task.positive_label,
            min_per_class=args.min_per_class,
        )
        per_task_results[task_name] = res
        save_dataframe(res, output_dir / f"per_gene_auc_{task_name}.csv")

    membership_df = pd.DataFrame(task_membership_rows).sort_values("task_name", key=lambda s: s.map(_task_sort_key)).reset_index(drop=True)
    save_dataframe(membership_df, output_dir / "task_membership_summary.csv")

    baseline_regions = compute_baseline_expression_regions(
        X_samples=X_samples,
        meta=meta,
        head_pct=args.head_pct,
        tail_pct=args.tail_pct,
    )
    save_dataframe(baseline_regions, output_dir / "baseline_expression_regions.csv")

    task_summary_df, region_summary_df = summarize_task_level_counts(
        per_task_results=per_task_results,
        baseline_regions=baseline_regions,
        fdr_threshold=args.fdr_threshold,
        strong_auc_threshold=args.strong_auc_threshold,
    )
    save_dataframe(task_summary_df, output_dir / "gene_screen_task_summary.csv")
    save_dataframe(region_summary_df, output_dir / "gene_screen_region_summary.csv")

    text_main, text_dox = make_gene_screen_text(
        task_summary_df,
        region_summary_df,
        fdr_threshold=args.fdr_threshold,
        strong_auc_threshold=args.strong_auc_threshold,
    )
    (text_dir / "gene_screen_results_paragraph.txt").write_text(text_main + "\n", encoding="utf-8")
    if text_dox:
        (text_dir / "gene_screen_dox_bridge_paragraph.txt").write_text(text_dox + "\n", encoding="utf-8")

    if not args.skip_summary_figures:
        plot_task_counts(task_summary_df, figures_dir / "gene_screen_significant_counts.png")
        plot_strong_auc_counts(task_summary_df, figures_dir / "gene_screen_auc80_counts.png", args.strong_auc_threshold)
        plot_region_composition(region_summary_df, figures_dir / "gene_screen_region_composition_sig.png")

        if args.focus_task in task_tables and args.focus_task in per_task_results:
            X_focus, meta_focus = task_tables[args.focus_task]
            plot_focus_heatmap(
                X_task=X_focus,
                task_meta=meta_focus,
                per_gene_df=per_task_results[args.focus_task],
                focus_task=args.focus_task,
                focus_context=args.focus_context,
                output_path=figures_dir / "gene_screen_hcpt24h_top_myc_heatmap.png",
                summary_path=output_dir / "gene_screen_hcpt24h_focus_summary.csv",
                top_n=args.focus_top_n,
            )

        if "MYC@24h" in per_task_results and "DOX@24h" in per_task_results:
            sig_a = per_task_results["MYC@24h"].loc[per_task_results["MYC@24h"]["fdr"] <= args.fdr_threshold].copy()
            sig_b = per_task_results["DOX@24h"].loc[per_task_results["DOX@24h"]["fdr"] <= args.fdr_threshold].copy()
            ov = overlap_summary(sig_a, sig_b, "MYC@24h", "DOX@24h")
            save_dataframe(ov, output_dir / "gene_screen_myc24_dox24_sig_overlap.csv")
            plot_overlap_summary(ov, figures_dir / "gene_screen_myc24_dox24_sig_overlap.png")

    ddr_genes = load_gene_set(args.ddr_genes)
    save_dataframe(
        pd.DataFrame({"gene_symbol": sorted(ddr_genes)}),
        output_dir / "ddr_gene_set_used.csv",
    )

    if args.skip_distortion_candidates:
        manifest = build_manifest(
            arguments=vars(args),
            input_files={"matrix": args.matrix, "sample_sheet": args.sample_sheet},
            script_files=[Path(__file__), Path(__file__).with_name("fig6_ml_core.py")],
            task_specs=task_specs,
            extra={"outputs_dir": str(output_dir.resolve()), "matrix_scale": args.matrix_scale, "matrix_transform": {"log1p_applied_in_loader": bool(args.log1p)}},
        )
        write_json(manifest, output_dir / "run_manifest_gene_screen.json")
        print(f"[OK] Wrote per-gene AUROC tables to {output_dir.resolve()}")
        return 0

    low_task = args.candidate_low_task
    high_task = args.candidate_high_task
    if low_task not in per_task_results or high_task not in per_task_results:
        raise ValueError(
            f"Candidate tasks not found. Available tasks: {sorted(per_task_results.keys())}"
        )

    low_df = per_task_results[low_task]
    high_df = per_task_results[high_task]
    merged = high_df[["gene", "gene_symbol", "auc", "fdr"]].rename(
        columns={"auc": f"auc_{high_task}", "fdr": f"fdr_{high_task}"}
    ).merge(
        low_df[["gene", "auc", "fdr"]].rename(columns={"auc": f"auc_{low_task}", "fdr": f"fdr_{low_task}"}),
        on="gene",
        how="inner",
    )

    candidates = merged.loc[
        (merged[f"auc_{low_task}"] <= args.low_auc_threshold)
        & (merged[f"auc_{high_task}"] >= args.high_auc_threshold)
    ].copy()
    candidates = candidates.merge(baseline_regions, on="gene", how="left")
    candidates["is_DDR"] = candidates["gene_symbol"].str.upper().isin(ddr_genes)
    save_dataframe(
        candidates.sort_values([f"auc_{high_task}", f"auc_{low_task}"], ascending=[False, True]),
        output_dir / f"distortion_candidates_{low_task}_to_{high_task}.csv",
    )

    background_genes = set(merged["gene_symbol"].str.upper())
    candidate_genes = set(candidates["gene_symbol"].str.upper())
    a = len(candidate_genes & ddr_genes)
    b = len(candidate_genes - ddr_genes)
    c = len((background_genes & ddr_genes) - candidate_genes)
    d = len(background_genes - ddr_genes - candidate_genes)
    table = np.array([[a, b], [c, d]], dtype=int)
    try:
        odds_ratio, pvalue = fisher_exact(table, alternative="greater")
    except ValueError:
        odds_ratio, pvalue = np.nan, np.nan

    overlap_summary_df = pd.DataFrame(
        [
            {
                "low_task": low_task,
                "high_task": high_task,
                "n_background_genes": len(background_genes),
                "n_candidate_genes": len(candidate_genes),
                "n_ddr_genes_in_background": len(background_genes & ddr_genes),
                "n_ddr_genes_in_candidates": a,
                "fisher_odds_ratio": odds_ratio,
                "fisher_pvalue": pvalue,
            }
        ]
    )
    save_dataframe(overlap_summary_df, output_dir / f"ddr_overlap_{low_task}_to_{high_task}.csv")
    if not args.skip_summary_figures:
        plot_ddr_summary(overlap_summary_df, figures_dir / "gene_screen_ddr_overlap_summary.png")

    manifest = build_manifest(
        arguments=vars(args),
        input_files={"matrix": args.matrix, "sample_sheet": args.sample_sheet},
        script_files=[Path(__file__), Path(__file__).with_name("fig6_ml_core.py")],
        task_specs=task_specs,
        extra={
            "outputs_dir": str(output_dir.resolve()),
            "candidate_comparison": {"low_task": low_task, "high_task": high_task},
            "ddr_gene_set_source": args.ddr_genes if args.ddr_genes else "DEFAULT_DDR_GENES",
            "summary_figures_dir": str(figures_dir.resolve()),
            "matrix_scale": args.matrix_scale,
            "matrix_transform": {"log1p_applied_in_loader": bool(args.log1p)},
        },
    )
    write_json(manifest, output_dir / "run_manifest_gene_screen.json")

    print(f"[OK] Wrote gene-screen outputs to {output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
