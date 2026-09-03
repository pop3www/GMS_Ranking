from __future__ import annotations

# SCRIPT_DIR_BOOTSTRAP
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import argparse
import json
from typing import Dict, Iterable, List

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig6_ml_core import InputContractError, build_manifest, save_dataframe, write_json

TASK_ORDER = ["MYC@4h", "MYC@24h", "DOX@4h", "DOX@24h"]
MODEL_LABELS = {
    "ridge_logistic": "ridge-logistic",
    "gbt_svd": "gradient-boosting/SVD",
}
CONTEXT_ORDER = {
    "cpt_level": ["None_CPT", "LCPT", "HCPT"],
    "myc": ["None_MYC", "MYC"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build manuscript-facing Figure 6 assets from an existing sample-level run: "
            "metric-first composite panel, supplementary embeddings, and text snippets."
        )
    )
    parser.add_argument("--run-dir", required=True, help="Output directory created by run_fig6_sample_level.py")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for manuscript assets. Defaults to <run-dir>/manuscript_assets.",
    )
    parser.add_argument(
        "--model",
        default="ridge_logistic",
        choices=["ridge_logistic", "gbt_svd"],
        help="Classifier model used for the main-text Figure 6 assets.",
    )
    parser.add_argument(
        "--cv-scheme",
        default="leave_one_sample_out",
        help="Cross-validation scheme used for the main metric panel.",
    )
    parser.add_argument(
        "--input-scale",
        default=None,
        help="Audit label for the expression matrix scale. If omitted, read from run_manifest_fig6.json when available.",
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        default=["MYC@4h", "MYC@24h"],
        help="Main-text task order. Defaults to the two MYC tasks.",
    )
    parser.add_argument(
        "--include-control-summary",
        action="store_true",
        help="Also write a small control-task metric summary if DOX control tasks are present.",
    )
    parser.add_argument(
        "--show-pvalues",
        action="store_true",
        help="Show empirical permutation p-values in metric panels. Default is off to match the frozen manuscript figure guidance.",
    )
    return parser.parse_args()



def _task_sort_key(task_name: str) -> int:
    return TASK_ORDER.index(task_name) if task_name in TASK_ORDER else len(TASK_ORDER)



def _marker_map(values: Iterable[str]) -> Dict[str, str]:
    markers = ["o", "s", "^", "D", "P", "X", "v", "*"]
    values = list(dict.fromkeys(values))
    return {val: markers[i % len(markers)] for i, val in enumerate(values)}



def _load_required_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required {label} file not found: {path}")
    return pd.read_csv(path)



def _format_context_values(text: str) -> str:
    raw_values = [v for v in str(text).split(";") if v]
    order = {"None_CPT": 0, "LCPT": 1, "HCPT": 2, "None_MYC": 0, "MYC": 1}
    raw_values = sorted(raw_values, key=lambda x: (order.get(x, 999), x))
    if raw_values == ["None_CPT", "LCPT", "HCPT"]:
        return "3 matched backgrounds"
    if raw_values == ["None_MYC", "MYC"]:
        return "MYC OFF/ON control"
    values = [v.replace("None_CPT", "no CPT") for v in raw_values]
    return ", ".join(values)



def _plot_task_overview(ax, task_summary: pd.DataFrame, cv_scheme: str):
    ax.axis("off")
    show = task_summary.copy()
    show["backgrounds"] = show["context_values"].map(_format_context_values)
    show["n"] = show["n_samples"].astype(int)
    show["MYC ON"] = show["n_positive"].astype(int)
    show["MYC OFF"] = show["n_negative"].astype(int)
    cv_label = {
        "leave_one_sample_out": "LOO by sample",
        "leave_one_context_out": "LOO by context",
    }.get(cv_scheme, cv_scheme.replace("_", " "))
    show["CV"] = cv_label
    show = show[["task_name", "n", "MYC ON", "MYC OFF", "backgrounds", "CV"]]
    show = show.rename(columns={"task_name": "task"})
    table = ax.table(
        cellText=show.values,
        colLabels=show.columns,
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.5)
    ax.set_title("Task definition and sample counts", pad=10)



def _plot_metric_panel(ax, metrics_df: pd.DataFrame, tasks: List[str], model_label: str, show_pvalues: bool = False):
    metric_labels = {
        "balanced_accuracy": ("Balanced accuracy", "o"),
        "roc_auc": ("ROC–AUC", "s"),
    }
    xpos = np.arange(len(tasks), dtype=float)
    offsets = {"balanced_accuracy": -0.08, "roc_auc": 0.08}

    for metric, (label, marker) in metric_labels.items():
        vals = []
        for task in tasks:
            row = metrics_df.loc[metrics_df["task_name"] == task]
            vals.append(float(row[metric].iloc[0]) if not row.empty else np.nan)
        ax.scatter(xpos + offsets[metric], vals, s=90, marker=marker, label=label)
        for xi, yi, task in zip(xpos + offsets[metric], vals, tasks):
            if pd.notna(yi):
                row = metrics_df.loc[metrics_df["task_name"] == task].iloc[0]
                pval = row.get("permutation_pvalue", np.nan)
                ann = f"{yi:.2f}"
                if show_pvalues and pd.notna(pval):
                    ann += f"\np={pval:.3g}"
                ax.text(xi, yi + 0.03, ann, ha="center", va="bottom", fontsize=8)

    ax.axhline(0.5, linestyle="--", linewidth=1)
    ax.set_ylim(0.0, 1.05)
    ax.set_xticks(xpos)
    ax.set_xticklabels(tasks)
    ax.set_ylabel("Cross-validated performance")
    ax.set_title(f"Primary metrics ({model_label})")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="lower right")



def _plot_oof_scores(ax, pred_df: pd.DataFrame, task_label: str, positive_label: str):
    if pred_df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No predictions found for {task_label}", ha="center", va="center")
        return

    context_col = pred_df["context_col"].iloc[0] if "context_col" in pred_df.columns else "context_value"
    order_values = CONTEXT_ORDER.get(str(context_col), list(dict.fromkeys(pred_df["context_value"].tolist())))
    order_lookup = {ctx: i for i, ctx in enumerate(order_values)}

    plot_df = pred_df.copy()
    plot_df["context_order"] = plot_df["context_value"].map(lambda x: order_lookup.get(x, 999))
    plot_df = plot_df.sort_values(["context_order", "label", "group_label", "replicate", "sample_id"]).reset_index(drop=True)
    plot_df["x"] = np.arange(len(plot_df), dtype=float)

    label_colors = {
        int(plot_df["label"].min()): "tab:orange",
        int(plot_df["label"].max()): "tab:green",
    }
    marker_map = _marker_map(plot_df["context_value"].tolist())

    for _, row in plot_df.iterrows():
        ax.scatter(
            row["x"],
            row["pred_score"],
            s=75,
            color=label_colors[int(row["label"])],
            marker=marker_map[row["context_value"]],
            edgecolor="black",
            linewidth=0.5,
        )
    prev_context = None
    for _, row in plot_df.iterrows():
        if prev_context is not None and row["context_value"] != prev_context:
            ax.axvline(row["x"] - 0.5, color="grey", linewidth=1, alpha=0.4)
        prev_context = row["context_value"]

    ax.axhline(0.5, linestyle="--", linewidth=1)
    ax.set_ylim(-0.02, 1.05)
    ax.set_ylabel(f"OOF P({positive_label})")
    ax.set_xlabel("Held-out samples")
    ax.set_title(task_label)
    ax.set_xticks(plot_df["x"])
    ax.set_xticklabels(plot_df["sample_id"], rotation=70, ha="right", fontsize=7)
    ax.grid(axis="y", alpha=0.25)



def _plot_embedding(ax, emb_df: pd.DataFrame, task_label: str):
    if emb_df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No embedding found for {task_label}", ha="center", va="center")
        return

    marker_map = _marker_map(emb_df["context_value"].tolist())
    label_colors = {
        int(emb_df["label"].min()): "tab:orange",
        int(emb_df["label"].max()): "tab:green",
    }
    for _, row in emb_df.iterrows():
        ax.scatter(
            row["dim1"],
            row["dim2"],
            s=85,
            color=label_colors[int(row["label"])],
            marker=marker_map[row["context_value"]],
            edgecolor="black",
            linewidth=0.5,
        )
        ax.text(row["dim1"], row["dim2"], f" {row['group_label']}_r{row['replicate']}", fontsize=7)
    ax.set_title(task_label)
    ax.set_xlabel("Dimension 1")
    ax.set_ylabel("Dimension 2")
    ax.grid(alpha=0.25)



def _trend_sentence(row4: pd.Series, row24: pd.Series) -> str:
    ba4 = float(row4["balanced_accuracy"])
    ba24 = float(row24["balanced_accuracy"])
    auc4 = float(row4["roc_auc"])
    auc24 = float(row24["roc_auc"])
    if (ba4 > ba24 and auc4 >= auc24) or (ba4 >= ba24 and auc4 > auc24):
        return (
            "Read conservatively, these results are more consistent with a clearer MYC-state boundary at 4 h than at "
            "24 h, rather than with a complete loss of late-state signal."
        )
    if (ba24 > ba4 and auc24 >= auc4) or (ba24 >= ba4 and auc24 > auc4):
        return (
            "Read conservatively, these results do not support a simpler early-cleaner versus late-weaker ordering, so "
            "the text should describe Figure 6 as a changed sample-level readout over time rather than as a monotonic "
            "decline in separation."
        )
    return (
        "The two primary metrics do not establish a uniform early-versus-late ordering, so Figure 6 should be framed as "
        "a secondary sample-level readout whose quantitative structure changes across time rather than as a simple "
        "monotonic decline in separation."
    )



def _results_paragraph(
    metrics_df: pd.DataFrame,
    task_summary: pd.DataFrame,
    model: str,
    cv_scheme: str,
) -> str:
    metric_lookup = {(row["task_name"], row["cv_scheme"]): row for _, row in metrics_df.iterrows()}
    summary_lookup = {row["task_name"]: row for _, row in task_summary.iterrows()}

    required = [("MYC@4h", cv_scheme), ("MYC@24h", cv_scheme)]
    missing = [f"{task}/{scheme}" for task, scheme in required if (task, scheme) not in metric_lookup]
    if missing:
        raise InputContractError("The results paragraph generator is missing required rows: " + ", ".join(missing))

    r4 = metric_lookup[("MYC@4h", cv_scheme)]
    r24 = metric_lookup[("MYC@24h", cv_scheme)]
    s4 = summary_lookup["MYC@4h"]
    s24 = summary_lookup["MYC@24h"]
    trend = _trend_sentence(r4, r24)
    model_label = MODEL_LABELS.get(model, model)

    paragraph = (
        "To ask whether the MYC ON/OFF boundary remained learnable at the sample level, we trained leakage-controlled "
        "classifiers on sample-level transcriptomes restricted to matched Dox-positive backgrounds at 4 h and 24 h. "
        "Figure 6 is presented as a secondary quantitative readout: the primary evidence is the cross-validated metric "
        f"panel and the held-out sample scores, not the two-dimensional embedding. Using the {model_label} model under "
        f"{cv_scheme.replace('_', '-')} cross-validation, MYC@4h yielded balanced accuracy {float(r4['balanced_accuracy']):.2f} "
        f"and AUROC {float(r4['roc_auc']):.2f} across {int(s4['n_samples'])} samples ({int(s4['n_positive'])} MYC ON, "
        f"{int(s4['n_negative'])} MYC OFF), whereas MYC@24h yielded balanced accuracy {float(r24['balanced_accuracy']):.2f} "
        f"and AUROC {float(r24['roc_auc']):.2f} across {int(s24['n_samples'])} samples ({int(s24['n_positive'])} MYC ON, "
        f"{int(s24['n_negative'])} MYC OFF)."
    )

    if ("MYC@4h", "leave_one_context_out") in metric_lookup and ("MYC@24h", "leave_one_context_out") in metric_lookup:
        c4 = metric_lookup[("MYC@4h", "leave_one_context_out")]
        c24 = metric_lookup[("MYC@24h", "leave_one_context_out")]
        paragraph += (
            f" Under the more stringent leave-one-context-out evaluation, performance changed from balanced accuracy "
            f"{float(c4['balanced_accuracy']):.2f} and AUROC {float(c4['roc_auc']):.2f} at 4 h to balanced accuracy "
            f"{float(c24['balanced_accuracy']):.2f} and AUROC {float(c24['roc_auc']):.2f} at 24 h, indicating that the "
            "later state generalized less cleanly across perturbation backgrounds."
        )

    paragraph += (
        f" {trend} Held-out score panels display individual biological replicates rather than replicate averages, "
        "so replicate-specific instability, including the late HCPT MYC-OFF ambiguity, remains visible. "
        "Doxycycline-control tasks were retained only as secondary controls outside the main figure."
    )
    return paragraph


def _dox_bridge_paragraph(metrics_df: pd.DataFrame, model: str, cv_scheme: str) -> str:
    metric_lookup = {(row["task_name"], row["cv_scheme"]): row for _, row in metrics_df.iterrows()}
    if ("DOX@4h", cv_scheme) not in metric_lookup or ("DOX@24h", cv_scheme) not in metric_lookup:
        return ""
    r4 = metric_lookup[("DOX@4h", cv_scheme)]
    r24 = metric_lookup[("DOX@24h", cv_scheme)]
    model_label = MODEL_LABELS.get(model, model)
    return (
        f"Doxycycline-only control tasks were retained as secondary controls. Using the {model_label} model under "
        f"{cv_scheme.replace('_', '-')} cross-validation, DOX@4h yielded balanced accuracy {float(r4['balanced_accuracy']):.2f} "
        f"and AUROC {float(r4['roc_auc']):.2f}, whereas DOX@24h yielded balanced accuracy {float(r24['balanced_accuracy']):.2f} "
        f"and AUROC {float(r24['roc_auc']):.2f}. These control readouts confirm that priming is not chemically silent "
        "at the sample level, but they are interpreted separately from the main MYC readout and do not by themselves "
        "define the distortion regime."
    )



def _legend_text() -> str:
    return (
        "Figure 6 | Sample-level classifier readouts provide a secondary view of MYC-state separation. "
        "(A) Task definition and sample counts for the two main-text MYC tasks, each restricted to Dox-positive "
        "backgrounds with matched no-CPT, LCPT, and HCPT contexts. (B) Leave-one-sample-out balanced accuracy and "
        "ROC–AUC for the primary ridge-logistic classifier; metrics, not the embedding, are the primary evidence in "
        "this figure. (C–D) Out-of-fold P(MYC ON) for held-out samples at 4 h (C) and 24 h (D), with marker shape "
        "denoting CPT background and colour denoting true class. Points are individual biological replicates rather than "
        "replicate averages; in the publication figure, thin grey lines connect replicate observations within the same "
        "class and background to make replicate spread explicit. A secondary two-dimensional display of the same "
        "sample-level structure can be shown separately, but is not interpreted independently of the cross-validated "
        "metric panel."
    )



def _load_run_scale(run_dir: Path) -> tuple[str, bool | None]:
    manifest_path = run_dir / "summary" / "run_manifest_fig6.json"
    if not manifest_path.exists():
        return "unspecified", None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception:
        return "unspecified", None
    args = manifest.get("arguments", {}) if isinstance(manifest, dict) else {}
    extra = manifest.get("extra", {}) if isinstance(manifest, dict) else {}
    scale = args.get("input_scale") or args.get("matrix_scale") or extra.get("matrix_scale") or "unspecified"
    log1p = args.get("log1p")
    if isinstance(log1p, str):
        log1p = log1p.lower() in {"1", "true", "yes"}
    elif log1p is not None:
        log1p = bool(log1p)
    return str(scale), log1p


def _scale_methods_sentence(matrix_scale: str, log1p: bool | None) -> str:
    if matrix_scale in {"log2_sf_norm_rsem_expected_counts", "log2_size_factor_normalized_rsem_expected_counts"}:
        return (
            "The primary expression input was a count-derived matrix, computed as "
            "log2(size-factor-normalized RSEM expected count + 1); no additional log transform was applied in the "
            "classifier workflow. "
        )
    if matrix_scale in {"tpm_log1p", "TPM_log1p"}:
        return (
            "This run used log1p-transformed TPM expression and should be interpreted as a TPM-scale sensitivity "
            "analysis rather than the primary count-derived run. "
        )
    if log1p is True:
        return "The input matrix was log1p-transformed inside the workflow before modeling. "
    if log1p is False:
        return "The input matrix was used on its supplied scale without an additional log1p transform. "
    return "The input expression scale is recorded in the run manifest. "


def _methods_blurb(matrix_scale: str = "unspecified", log1p: bool | None = None) -> str:
    scale_sentence = _scale_methods_sentence(matrix_scale, log1p)
    return (
        "Classifier analyses were used as secondary sample-level readouts rather than as stand-alone proof of mechanism. "
        + scale_sentence
        + "For each task, one row represented one sample and genes were used as predictors; sample metadata were used only "
        "to define task membership, labels, and grouped cross-validation splits. The main-text tasks were MYC@4h and "
        "MYC@24h, each restricted to Dox-positive backgrounds with matched no-CPT, LCPT, and HCPT contexts. Within each "
        "training fold, the top-variance genes were selected and standardized. A ridge-penalized logistic regression model "
        "served as the primary classifier, with a shallow gradient-boosting classifier on truncated-SVD features used as a "
        "sensitivity analysis. Primary performance metrics were balanced accuracy and ROC–AUC under leave-one-sample-out "
        "cross-validation, with leave-one-context-out analyses used as robustness checks when each context contained both "
        "classes. Out-of-fold sample scores were displayed directly as individual biological replicates, without replicate "
        "averaging, and two-dimensional sample-level embeddings were shown "
        "only as secondary qualitative summaries of the model-transformed feature space. Univariate per-gene AUROC screens "
        "were run in a separate workflow and were not interpreted as model-native feature importance."
    )



def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    output_dir = Path(args.output_dir) if args.output_dir is not None else run_dir / "manuscript_assets"
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_dir = run_dir / "summary"
    pred_dir = run_dir / "predictions"
    emb_dir = run_dir / "embeddings"

    metrics = _load_required_csv(summary_dir / "metrics_all_models.csv", "metrics")
    task_summary = _load_required_csv(summary_dir / "task_membership_summary.csv", "task summary")

    metrics = metrics.copy()
    metrics["task_order"] = metrics["task_name"].map(_task_sort_key)
    task_summary["task_order"] = task_summary["task_name"].map(_task_sort_key)

    main_metrics = metrics.loc[
        (metrics["model_name"] == args.model) & (metrics["cv_scheme"] == args.cv_scheme) & (metrics["task_name"].isin(args.tasks))
    ].copy()
    main_metrics = main_metrics.sort_values("task_order").reset_index(drop=True)
    main_task_summary = task_summary.loc[task_summary["task_name"].isin(args.tasks)].copy()
    main_task_summary = main_task_summary.sort_values("task_order").reset_index(drop=True)

    if main_metrics.empty:
        raise InputContractError(
            f"No metric rows found for model={args.model!r}, cv_scheme={args.cv_scheme!r}, tasks={args.tasks}."
        )
    missing_tasks = [t for t in args.tasks if t not in set(main_metrics["task_name"])]
    if missing_tasks:
        raise InputContractError("Missing metric rows for tasks: " + ", ".join(missing_tasks))

    figures_dir = output_dir / "figures"
    text_dir = output_dir / "text"
    figures_dir.mkdir(parents=True, exist_ok=True)
    text_dir.mkdir(parents=True, exist_ok=True)

    save_dataframe(main_metrics, output_dir / "figure6_main_metrics.csv")
    save_dataframe(main_task_summary, output_dir / "figure6_task_summary.csv")

    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 1.15])
    ax_table = fig.add_subplot(gs[0, 0])
    ax_metrics = fig.add_subplot(gs[0, 1])
    ax_oof_4h = fig.add_subplot(gs[1, 0])
    ax_oof_24h = fig.add_subplot(gs[1, 1])

    _plot_task_overview(ax_table, main_task_summary, args.cv_scheme)
    _plot_metric_panel(ax_metrics, main_metrics, args.tasks, MODEL_LABELS.get(args.model, args.model), show_pvalues=args.show_pvalues)

    for ax, task in zip([ax_oof_4h, ax_oof_24h], args.tasks[:2]):
        pred_path = pred_dir / f"oof_predictions_{task}_{args.model}_{args.cv_scheme}.csv"
        pred_df = _load_required_csv(pred_path, f"predictions for {task}")
        pos_label = str(main_metrics.loc[main_metrics["task_name"] == task, "positive_label"].iloc[0])
        _plot_oof_scores(ax, pred_df, task, pos_label)

    fig.suptitle("Figure 6 (main-text layout): metric-first sample-level classifier readout", fontsize=14)
    fig.savefig(figures_dir / "figure6_metric_first_main.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    emb_tasks = [t for t in args.tasks if (emb_dir / f"embedding_{t}_{args.model}.csv").exists()]
    if emb_tasks:
        fig2, axes = plt.subplots(1, len(emb_tasks), figsize=(6 * len(emb_tasks), 5), constrained_layout=True)
        if len(emb_tasks) == 1:
            axes = [axes]
        for ax, task in zip(axes, emb_tasks):
            emb_df = _load_required_csv(emb_dir / f"embedding_{task}_{args.model}.csv", f"embedding for {task}")
            _plot_embedding(ax, emb_df, task)
        fig2.suptitle("Figure 6 secondary sample-level embeddings", fontsize=14)
        fig2.savefig(figures_dir / "figure6_secondary_embeddings.png", dpi=300, bbox_inches="tight")
        plt.close(fig2)

    if args.include_control_summary:
        control_metrics = metrics.loc[
            (metrics["model_name"] == args.model) & (metrics["cv_scheme"] == args.cv_scheme) & (~metrics["task_name"].isin(args.tasks))
        ].copy()
        if not control_metrics.empty:
            control_metrics = control_metrics.sort_values("task_order").reset_index(drop=True)
            save_dataframe(control_metrics, output_dir / "figure6_control_metrics.csv")
            fig3, ax3 = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
            _plot_metric_panel(ax3, control_metrics, control_metrics["task_name"].tolist(), MODEL_LABELS.get(args.model, args.model), show_pvalues=args.show_pvalues)
            ax3.set_title("Control-task metrics (supplementary)")
            fig3.savefig(figures_dir / "figure6_control_metrics.png", dpi=300, bbox_inches="tight")
            plt.close(fig3)

    model_metrics = metrics.loc[metrics["model_name"] == args.model].copy()
    results_paragraph = _results_paragraph(model_metrics, task_summary, args.model, args.cv_scheme)
    legend_text = _legend_text()
    matrix_scale, matrix_log1p = _load_run_scale(run_dir)
    if args.input_scale is not None:
        matrix_scale = args.input_scale
    methods_blurb = _methods_blurb(matrix_scale, matrix_log1p)
    (text_dir / "figure6_results_paragraph.txt").write_text(results_paragraph + "\n", encoding="utf-8")
    (text_dir / "figure6_legend.txt").write_text(legend_text + "\n", encoding="utf-8")
    (text_dir / "figure6_methods_blurb.txt").write_text(methods_blurb + "\n", encoding="utf-8")
    dox_bridge = _dox_bridge_paragraph(model_metrics, args.model, args.cv_scheme)
    if dox_bridge:
        (text_dir / "figure6_dox_bridge_paragraph.txt").write_text(dox_bridge + "\n", encoding="utf-8")

    manifest = build_manifest(
        arguments=vars(args),
        input_files={
            "metrics_all_models": summary_dir / "metrics_all_models.csv",
            "task_membership_summary": summary_dir / "task_membership_summary.csv",
        },
        script_files=[Path(__file__), Path(__file__).with_name("fig6_ml_core.py")],
        extra={
            "run_dir": str(run_dir.resolve()),
            "output_dir": str(output_dir.resolve()),
            "main_metric_rows": int(len(main_metrics)),
            "main_figure": str((figures_dir / "figure6_metric_first_main.png").resolve()),
        },
    )
    write_json(manifest, output_dir / "run_manifest_manuscript_assets.json")

    print("[OK] Wrote manuscript-facing Figure 6 assets.")
    print(f"[OK] Output directory: {output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
