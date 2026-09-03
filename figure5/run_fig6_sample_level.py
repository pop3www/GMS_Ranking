from __future__ import annotations

# SCRIPT_DIR_BOOTSTRAP
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import argparse
from typing import Dict, List

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig6_ml_core import (
    InputContractError,
    build_cv_schemes,
    build_estimator,
    build_manifest,
    build_task_table,
    compute_display_embedding,
    cross_validate_binary_model,
    default_task_specs,
    expression_to_sample_by_gene,
    extract_linear_coefficients,
    fit_final_model,
    permutation_test_balanced_accuracy,
    read_sample_sheet,
    read_wide_matrix,
    save_array_as_csv,
    save_dataframe,
    validate_and_align_inputs,
    write_json,
)

TASK_ORDER = ["MYC@4h", "MYC@24h", "DOX@4h", "DOX@24h"]
MODEL_ORDER = ["ridge_logistic", "gbt_svd"]
CONTEXT_ORDER = {
    "cpt_level": ["None_CPT", "LCPT", "HCPT"],
    "myc": ["None_MYC", "MYC"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Production rewrite of Figure 6: sample-level transcriptome classifiers, "
            "out-of-fold metrics, and sample-level embeddings."
        )
    )
    parser.add_argument("--matrix", required=True, help="Wide gene-by-sample matrix (CSV/TSV).")
    parser.add_argument("--sample-sheet", required=True, help="Sample sheet CSV/TSV matching the matrix columns.")
    parser.add_argument("--output-dir", default="fig6_outputs", help="Output directory.")
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
        "--models",
        nargs="+",
        default=["ridge_logistic", "gbt_svd"],
        choices=["ridge_logistic", "gbt_svd"],
        help="One or more classifier backends to evaluate.",
    )
    parser.add_argument(
        "--figure-model",
        default="ridge_logistic",
        choices=["ridge_logistic", "gbt_svd"],
        help="Model used for the main panel-style plots and embeddings.",
    )
    parser.add_argument("--max-genes", type=int, default=500, help="Top variable genes kept within each training fold.")
    parser.add_argument(
        "--max-svd-components",
        type=int,
        default=6,
        help="Upper bound on fold-wise SVD dimensionality for the GBT model.",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=0,
        help=(
            "Optional number of label permutations for empirical balanced-accuracy p-values. "
            "Applied to leave-one-sample-out only."
        ),
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--skip-dox-controls",
        action="store_true",
        help="Run only the two MYC tasks and skip DOX control tasks.",
    )
    return parser.parse_args()


def _task_sort_key(task_name: str) -> int:
    return TASK_ORDER.index(task_name) if task_name in TASK_ORDER else len(TASK_ORDER)


def _model_sort_key(model_name: str) -> int:
    return MODEL_ORDER.index(model_name) if model_name in MODEL_ORDER else len(MODEL_ORDER)


def _context_marker_map(context_values: List[str]) -> Dict[str, str]:
    markers = ["o", "s", "^", "D", "P", "X", "v", "*"]
    return {ctx: markers[i % len(markers)] for i, ctx in enumerate(context_values)}


def plot_metric_summary(metrics_df: pd.DataFrame, out_path: Path, cv_scheme: str):
    sub = metrics_df.loc[metrics_df["cv_scheme"] == cv_scheme].copy()
    if sub.empty:
        return

    sub["task_order"] = sub["task_name"].map(_task_sort_key)
    sub["model_order"] = sub["model_name"].map(_model_sort_key)
    sub = sub.sort_values(["task_order", "model_order"]).reset_index(drop=True)

    tasks = sub["task_name"].drop_duplicates().tolist()
    models = sub["model_name"].drop_duplicates().tolist()
    x = np.arange(len(tasks), dtype=float)
    offsets = np.linspace(-0.15, 0.15, num=max(len(models), 1))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    metric_fields = [
        ("balanced_accuracy", "Balanced accuracy"),
        ("roc_auc", "ROC-AUC"),
    ]
    for ax, (metric, title) in zip(axes, metric_fields):
        for offset, model_name in zip(offsets, models):
            model_sub = sub.loc[sub["model_name"] == model_name].set_index("task_name")
            y = [model_sub.loc[t, metric] if t in model_sub.index else np.nan for t in tasks]
            ax.scatter(x + offset, y, s=60, label=model_name)
            for xi, yi in zip(x + offset, y):
                if pd.notna(yi):
                    ax.text(xi, yi + 0.015, f"{yi:.2f}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(tasks, rotation=30, ha="right")
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel(title)
        ax.set_title(f"{title} ({cv_scheme})")
        ax.axhline(0.5, linestyle="--", linewidth=1)
        ax.grid(axis="y", alpha=0.25)
    axes[1].legend(frameon=False, loc="lower right")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_oof_scores(pred_df: pd.DataFrame, task_label: str, positive_label: str, out_path: Path):
    if pred_df.empty:
        return

    context_order = CONTEXT_ORDER.get(
        pred_df["context_col"].iloc[0] if "context_col" in pred_df.columns else "",
        sorted(pred_df["context_value"].unique().tolist()),
    )
    order_lookup = {ctx: i for i, ctx in enumerate(context_order)}

    plot_df = pred_df.copy()
    plot_df["context_order"] = plot_df["context_value"].map(lambda x: order_lookup.get(x, 999))
    plot_df = plot_df.sort_values(
        ["context_order", "label", "group_label", "replicate", "sample_id"]
    ).reset_index(drop=True)
    plot_df["x"] = np.arange(len(plot_df), dtype=float)

    label_colors = {
        int(plot_df["label"].min()): "tab:orange",
        int(plot_df["label"].max()): "tab:green",
    }
    marker_map = _context_marker_map(plot_df["context_value"].drop_duplicates().tolist())

    fig, ax = plt.subplots(figsize=(max(7.5, 0.6 * len(plot_df)), 4.5), constrained_layout=True)
    for _, row in plot_df.iterrows():
        ax.scatter(
            row["x"],
            row["pred_score"],
            s=80,
            color=label_colors[int(row["label"])],
            marker=marker_map[row["context_value"]],
            edgecolor="black",
            linewidth=0.6,
        )
        ax.text(
            row["x"],
            row["pred_score"] + 0.03,
            f"{row['group_label']}_r{row['replicate']}",
            rotation=60,
            fontsize=7,
            ha="left",
            va="bottom",
        )

    prev_context = None
    for _, row in plot_df.iterrows():
        if prev_context is not None and row["context_value"] != prev_context:
            ax.axvline(row["x"] - 0.5, color="grey", linewidth=1, alpha=0.4)
        prev_context = row["context_value"]

    ax.axhline(0.5, linestyle="--", linewidth=1)
    ax.set_ylim(-0.02, 1.05)
    ax.set_ylabel(f"Out-of-fold P({positive_label})")
    ax.set_xlabel("Held-out samples")
    ax.set_title(f"{task_label} — sample-level out-of-fold scores")
    ax.set_xticks(plot_df["x"])
    ax.set_xticklabels(plot_df["sample_id"], rotation=70, ha="right", fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)



def plot_embedding(embedding_df: pd.DataFrame, task_label: str, out_path: Path):
    if embedding_df.empty:
        return

    marker_map = _context_marker_map(embedding_df["context_value"].drop_duplicates().tolist())
    label_colors = {
        int(embedding_df["label"].min()): "tab:orange",
        int(embedding_df["label"].max()): "tab:green",
    }

    fig, ax = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
    for _, row in embedding_df.iterrows():
        ax.scatter(
            row["dim1"],
            row["dim2"],
            s=95,
            color=label_colors[int(row["label"])],
            marker=marker_map[row["context_value"]],
            edgecolor="black",
            linewidth=0.6,
        )
        ax.text(
            row["dim1"],
            row["dim2"],
            f" {row['group_label']}_r{row['replicate']}",
            fontsize=8,
            ha="left",
            va="center",
        )

    ax.set_xlabel("Embedding dimension 1")
    ax.set_ylabel("Embedding dimension 2")
    ax.set_title(f"{task_label} — sample-level embedding")
    ax.grid(alpha=0.25)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)



def main() -> int:
    args = parse_args()
    if args.figure_model not in args.models:
        raise InputContractError("--figure-model must also be listed in --models.")

    output_dir = Path(args.output_dir)
    summary_dir = output_dir / "summary"
    pred_dir = output_dir / "predictions"
    emb_dir = output_dir / "embeddings"
    coef_dir = output_dir / "coefficients"
    fig_dir = output_dir / "figures"

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

    tasks = default_task_specs(include_dox_controls=not args.skip_dox_controls)
    all_metrics = []
    task_membership_rows = []

    for task in tasks:
        X_task, task_meta = build_task_table(X_samples, meta, task)
        task_meta["context_col"] = task.context_col
        save_dataframe(task_meta, summary_dir / f"task_samples_{task.name}.csv")
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
                "context_values": ";".join(task_meta["context_value"].drop_duplicates().tolist()),
            }
        )
        schemes = build_cv_schemes(task_meta)

        for model_name in args.models:
            estimator = build_estimator(
                model_name=model_name,
                random_state=args.seed,
                max_genes=args.max_genes,
                max_svd_components=args.max_svd_components,
            )

            for scheme_name, splits in schemes.items():
                cv_res = cross_validate_binary_model(
                    estimator=estimator,
                    X=X_task,
                    task_meta=task_meta,
                    splits=splits,
                    task_name=task.name,
                    model_name=model_name,
                    cv_scheme=scheme_name,
                )

                pred_out = pred_dir / f"oof_predictions_{task.name}_{model_name}_{scheme_name}.csv"
                save_dataframe(cv_res.predictions, pred_out)
                conf_out = pred_dir / f"confusion_{task.name}_{model_name}_{scheme_name}.csv"
                save_array_as_csv(
                    cv_res.confusion,
                    row_labels=[task.negative_label, task.positive_label],
                    col_labels=[task.negative_label, task.positive_label],
                    out_path=conf_out,
                )

                metric_row = {
                    "task_name": task.name,
                    "task_description": task.description,
                    "model_name": model_name,
                    "cv_scheme": scheme_name,
                    "positive_label": task.positive_label,
                    "negative_label": task.negative_label,
                    "include_in_main": task.include_in_main,
                }
                metric_row.update(cv_res.metrics)

                if scheme_name == "leave_one_sample_out":
                    perm_stats = permutation_test_balanced_accuracy(
                        estimator=estimator,
                        X=X_task,
                        task_meta=task_meta,
                        splits=splits,
                        n_permutations=args.n_permutations,
                        random_state=args.seed,
                        permute_within_context=True,
                    )
                else:
                    perm_stats = {
                        "null_mean_balanced_accuracy": np.nan,
                        "null_sd_balanced_accuracy": np.nan,
                        "permutation_pvalue": np.nan,
                        "n_valid_permutations": 0,
                    }
                metric_row.update(perm_stats)
                all_metrics.append(metric_row)

                if model_name == args.figure_model and scheme_name == "leave_one_sample_out":
                    plot_oof_scores(
                        pred_df=cv_res.predictions,
                        task_label=task.name,
                        positive_label=task.positive_label,
                        out_path=fig_dir / f"oof_scores_{task.name}_{model_name}.png",
                    )

            final_estimator = fit_final_model(
                estimator=estimator,
                X=X_task,
                task_meta=task_meta,
                task_name=task.name,
                model_name=model_name,
            )
            embedding = compute_display_embedding(final_estimator, X_task)
            embedding = pd.concat([task_meta.reset_index(drop=True), embedding], axis=1)
            emb_out = emb_dir / f"embedding_{task.name}_{model_name}.csv"
            save_dataframe(embedding, emb_out)

            coef_df = extract_linear_coefficients(final_estimator)
            if coef_df is not None:
                save_dataframe(coef_df, coef_dir / f"coefficients_{task.name}_{model_name}.csv")

            if model_name == args.figure_model:
                plot_embedding(
                    embedding_df=embedding,
                    task_label=task.name,
                    out_path=fig_dir / f"embedding_{task.name}_{model_name}.png",
                )

    metrics_df = pd.DataFrame(all_metrics)
    if metrics_df.empty:
        raise RuntimeError("No Figure 6 outputs were generated. Check the matrix/sample sheet contract.")

    metrics_df["task_order"] = metrics_df["task_name"].map(_task_sort_key)
    metrics_df["model_order"] = metrics_df["model_name"].map(_model_sort_key)
    metrics_df = metrics_df.sort_values(["task_order", "model_order", "cv_scheme"]).reset_index(drop=True)

    save_dataframe(metrics_df, summary_dir / "metrics_all_models.csv")
    save_dataframe(metrics_df.loc[metrics_df["include_in_main"]], summary_dir / "metrics_main_text.csv")
    save_dataframe(pd.DataFrame(task_membership_rows), summary_dir / "task_membership_summary.csv")
    save_dataframe(pd.DataFrame([task.__dict__ for task in tasks]), summary_dir / "task_definitions.csv")

    for scheme_name in metrics_df["cv_scheme"].drop_duplicates().tolist():
        plot_metric_summary(
            metrics_df=metrics_df,
            out_path=fig_dir / f"metric_summary_{scheme_name}.png",
            cv_scheme=scheme_name,
        )

    manifest = build_manifest(
        arguments=vars(args),
        input_files={
            "matrix": args.matrix,
            "sample_sheet": args.sample_sheet,
        },
        script_files=[
            Path(__file__),
            Path(__file__).with_name("fig6_ml_core.py"),
        ],
        task_specs=tasks,
        extra={
            "outputs": {
                "summary_dir": str(summary_dir.resolve()),
                "predictions_dir": str(pred_dir.resolve()),
                "embeddings_dir": str(emb_dir.resolve()),
                "coefficients_dir": str(coef_dir.resolve()),
                "figures_dir": str(fig_dir.resolve()),
            },
            "task_membership_summary_file": str((summary_dir / "task_membership_summary.csv").resolve()),
            "metrics_file": str((summary_dir / "metrics_all_models.csv").resolve()),
            "matrix_scale": args.matrix_scale,
            "matrix_transform": {"log1p_applied_in_loader": bool(args.log1p)},
        },
    )
    write_json(manifest, summary_dir / "run_manifest_fig6.json")

    print("[OK] Figure 6 rewrite completed.")
    print(f"[OK] Matrix scale: {args.matrix_scale}; apply_log1p={args.log1p}")
    print(f"[OK] Outputs written to: {output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
