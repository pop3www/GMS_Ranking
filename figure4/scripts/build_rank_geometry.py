from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

from utils import (
    add_baseline_bands,
    canonicalize_label,
    canonicalize_mapping,
    compute_rolling_summary,
    load_yaml,
    read_table,
    stable_desc_rank,
)


REQUIRED_LONG_COLUMNS = [
    "gene_id",
    "baseline_mode",
    "comparison_id",
    "baseline_label",
    "condition_label",
    "baseline_value",
    "condition_value",
    "baseline_rank",
    "condition_rank",
    "rank_shift",
    "baseline_percentile",
    "band",
    "time_h",
    "dox",
    "tam",
    "cpt_level",
]


def _get_cfg(config: Mapping[str, Any], *keys: str, default: Any = None) -> Any:
    node: Any = config
    for key in keys:
        if not isinstance(node, Mapping) or key not in node:
            return default
        node = node[key]
    return node


def _normalize_sample_sheet(sample_sheet: pd.DataFrame, config: Mapping[str, Any]) -> pd.DataFrame:
    sscfg = dict(_get_cfg(config, "sample_sheet", default={}))
    required = {
        "sample_id_column": sscfg.get("sample_id_column", "sample_id"),
        "condition_label_column": sscfg.get("condition_label_column", "condition_label"),
        "time_h_column": sscfg.get("time_h_column", "time_h"),
        "dox_column": sscfg.get("dox_column", "dox"),
        "tam_column": sscfg.get("tam_column", "tam"),
        "cpt_level_column": sscfg.get("cpt_level_column", "cpt_level"),
    }
    missing = [col for col in required.values() if col not in sample_sheet.columns]
    if missing:
        raise ValueError(
            "Sample sheet is missing required columns: "
            + ", ".join(missing)
            + ". Update config/sample_sheet column names or the sample sheet."
        )

    label_map = dict(_get_cfg(config, "analysis", "label_map", default={}))
    out = sample_sheet.copy()
    out[required["sample_id_column"]] = out[required["sample_id_column"]].astype(str)
    out = out.rename(
        columns={
            required["sample_id_column"]: "sample_id",
            required["condition_label_column"]: "condition_label_raw",
            required["time_h_column"]: "time_h",
            required["dox_column"]: "dox",
            required["tam_column"]: "tam",
            required["cpt_level_column"]: "cpt_level",
        }
    )
    out["condition_label"] = out["condition_label_raw"].map(lambda x: canonicalize_label(x, label_map))
    return out


def aggregate_expression_by_condition(
    expr_df: pd.DataFrame,
    sample_sheet: pd.DataFrame,
    config: Mapping[str, Any],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate a gene-by-sample matrix to a gene-by-condition matrix."""
    exprcfg = dict(_get_cfg(config, "expression", default={}))
    gene_id_col = exprcfg.get("gene_id_column", "gene_id")
    aggregate = str(_get_cfg(config, "analysis", "aggregate", default="mean")).lower()

    if gene_id_col not in expr_df.columns:
        raise ValueError(f"Expression matrix is missing the configured gene_id column: {gene_id_col}")

    out = expr_df.copy()
    out.columns = [str(c) for c in out.columns]
    out[gene_id_col] = out[gene_id_col].astype(str)

    if out[gene_id_col].duplicated().any():
        dupes = out.loc[out[gene_id_col].duplicated(), gene_id_col].head(10).tolist()
        raise ValueError(
            "Expression matrix gene_id column contains duplicates. "
            f"Examples: {dupes}. Deduplicate upstream before rank-geometry analysis."
        )

    sample_ids = sample_sheet["sample_id"].astype(str).tolist()
    missing_samples = [sid for sid in sample_ids if sid not in out.columns]
    if missing_samples:
        raise ValueError(
            "The following sample IDs are present in the sample sheet but missing from the expression matrix: "
            + ", ".join(missing_samples[:20])
        )

    grouped = sample_sheet.groupby("condition_label", sort=False)

    metadata_rows: List[Dict[str, Any]] = []
    agg_df = pd.DataFrame({"gene_id": out[gene_id_col]})

    for condition_label, group in grouped:
        for meta_col in ["time_h", "dox", "tam", "cpt_level"]:
            unique_values = pd.unique(group[meta_col].dropna())
            if len(unique_values) > 1:
                raise ValueError(
                    f"Condition '{condition_label}' has inconsistent '{meta_col}' values in the sample sheet: "
                    f"{unique_values.tolist()}"
                )
        metadata_rows.append(
            {
                "condition_label": condition_label,
                "time_h": group["time_h"].dropna().iloc[0] if group["time_h"].notna().any() else None,
                "dox": group["dox"].dropna().iloc[0] if group["dox"].notna().any() else None,
                "tam": group["tam"].dropna().iloc[0] if group["tam"].notna().any() else None,
                "cpt_level": group["cpt_level"].dropna().iloc[0] if group["cpt_level"].notna().any() else None,
                "n_samples": int(len(group)),
                "sample_ids": ";".join(group["sample_id"].astype(str).tolist()),
            }
        )

        samples = group["sample_id"].astype(str).tolist()
        values = out[samples]
        if aggregate == "mean":
            agg_df[condition_label] = values.mean(axis=1, skipna=True)
        elif aggregate == "median":
            agg_df[condition_label] = values.median(axis=1, skipna=True)
        else:
            raise ValueError(f"Unsupported aggregation method '{aggregate}'. Use 'mean' or 'median'.")

    metadata_df = pd.DataFrame(metadata_rows).set_index("condition_label", drop=False)
    return agg_df, metadata_df


def _resolve_condition_order(
    available_conditions: Sequence[str],
    config: Mapping[str, Any],
) -> List[str]:
    label_map = dict(_get_cfg(config, "analysis", "label_map", default={}))
    configured_order = _get_cfg(config, "analysis", "condition_order", default=None)
    if configured_order is None:
        return list(available_conditions)

    resolved = [canonicalize_label(label, label_map) for label in configured_order]
    return [label for label in resolved if label in set(available_conditions)]


def _mapping_for_mode(mode: str, config: Mapping[str, Any], label_map: Mapping[str, str]) -> Optional[Dict[str, str]]:
    if mode == "matched":
        raw_mapping = _get_cfg(config, "baseline", "matched", "mapping", default={})
    elif mode == "global":
        return None
    else:
        raw_mapping = _get_cfg(config, "baseline", "named_modes", mode, "mapping", default={})
        if not raw_mapping:
            raw_mapping = _get_cfg(config, "baseline", mode, "mapping", default={})
    if not raw_mapping:
        return None
    return canonicalize_mapping(dict(raw_mapping), label_map=label_map)


def build_baseline_assignments(
    available_conditions: Sequence[str],
    config: Mapping[str, Any],
) -> List[Dict[str, str]]:
    """Build comparison assignments for enabled baseline modes."""
    label_map = dict(_get_cfg(config, "analysis", "label_map", default={}))
    condition_order = _resolve_condition_order(available_conditions, config)
    available = set(available_conditions)

    baseline_cfg = dict(_get_cfg(config, "baseline", default={}))
    enabled_modes = list(baseline_cfg.get("enabled_modes", []))
    if not enabled_modes:
        enabled_modes = [str(baseline_cfg.get("default_mode", "matched"))]

    assignments: List[Dict[str, str]] = []

    for mode in enabled_modes:
        if mode == "global":
            global_label = canonicalize_label(_get_cfg(config, "baseline", "global", "label", default=""), label_map)
            if global_label not in available:
                raise ValueError(
                    f"Global baseline '{global_label}' is not available among aggregated conditions: {sorted(available)}"
                )
            targets = [cond for cond in condition_order if cond != global_label]
            for condition_label in targets:
                assignments.append(
                    {
                        "baseline_mode": "global",
                        "baseline_label": global_label,
                        "condition_label": condition_label,
                    }
                )
            continue

        mapping = _mapping_for_mode(mode, config, label_map)
        if not mapping:
            raise ValueError(
                f"Baseline mode '{mode}' requires a mapping under baseline.matched.mapping or baseline.named_modes.{mode}.mapping"
            )

        for condition_label in condition_order:
            if condition_label not in mapping:
                continue
            baseline_label = mapping[condition_label]
            if baseline_label not in available:
                raise ValueError(
                    f"Baseline '{baseline_label}' for condition '{condition_label}' in mode '{mode}' "
                    "is not present in the aggregated condition table."
                )
            assignments.append(
                {
                    "baseline_mode": mode,
                    "baseline_label": baseline_label,
                    "condition_label": condition_label,
                }
            )

    if not assignments:
        raise ValueError("No baseline assignments were generated. Check condition_order and baseline config.")

    return assignments


def build_rank_geometry_from_frames(
    expr_df: pd.DataFrame,
    sample_sheet_df: pd.DataFrame,
    config: Mapping[str, Any],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build the source-of-truth long dataframe and rolling summaries."""
    sample_sheet = _normalize_sample_sheet(sample_sheet_df, config)
    condition_values_df, condition_metadata_df = aggregate_expression_by_condition(expr_df, sample_sheet, config)
    available_conditions = [c for c in condition_values_df.columns if c != "gene_id"]
    assignments = build_baseline_assignments(available_conditions, config)

    band_cfg = dict(_get_cfg(config, "analysis", "band_percentiles", default={}))
    head_pct = float(band_cfg.get("head", 10))
    mid_pct = float(band_cfg.get("mid", 50))
    tail_pct = float(band_cfg.get("tail", 10))

    frames: List[pd.DataFrame] = []

    for assignment in assignments:
        baseline_label = assignment["baseline_label"]
        condition_label = assignment["condition_label"]
        baseline_mode = assignment["baseline_mode"]

        baseline_df = stable_desc_rank(
            condition_values_df[["gene_id", baseline_label]].rename(columns={baseline_label: "baseline_value"}),
            gene_col="gene_id",
            value_col="baseline_value",
            rank_col="baseline_rank",
        )
        condition_df = stable_desc_rank(
            condition_values_df[["gene_id", condition_label]].rename(columns={condition_label: "condition_value"}),
            gene_col="gene_id",
            value_col="condition_value",
            rank_col="condition_rank",
        )

        merged = baseline_df.merge(condition_df, on="gene_id", how="inner")
        merged["baseline_label"] = baseline_label
        merged["condition_label"] = condition_label
        merged["baseline_mode"] = baseline_mode
        merged["comparison_id"] = f"{baseline_mode}::{baseline_label}__vs__{condition_label}"
        merged["rank_shift"] = merged["baseline_rank"] - merged["condition_rank"]
        merged["order"] = merged["baseline_rank"]
        merged["n_genes"] = int(len(merged))

        condition_meta = condition_metadata_df.loc[condition_label]
        merged["time_h"] = condition_meta["time_h"]
        merged["dox"] = condition_meta["dox"]
        merged["tam"] = condition_meta["tam"]
        merged["cpt_level"] = condition_meta["cpt_level"]

        merged = add_baseline_bands(
            merged,
            baseline_rank_col="baseline_rank",
            head_pct=head_pct,
            mid_pct=mid_pct,
            tail_pct=tail_pct,
        )
        frames.append(merged)

    long_df = pd.concat(frames, ignore_index=True)
    long_df = long_df.sort_values(
        by=["baseline_mode", "baseline_label", "condition_label", "baseline_rank", "gene_id"],
        kind="mergesort",
    ).reset_index(drop=True)

    rolling_cfg = dict(_get_cfg(config, "rolling", default={}))
    windows = rolling_cfg.get("windows", [25, 100, 250, 500])
    rolling_df = compute_rolling_summary(
        long_df,
        windows=windows,
        center=bool(rolling_cfg.get("center", True)),
        min_periods=rolling_cfg.get("min_periods", None),
    ).sort_values(
        by=["baseline_mode", "baseline_label", "condition_label", "window", "baseline_rank"],
        kind="mergesort",
    ).reset_index(drop=True)

    missing_required = [col for col in REQUIRED_LONG_COLUMNS if col not in long_df.columns]
    if missing_required:
        raise RuntimeError("Long rank-geometry dataframe is missing required columns: " + ", ".join(missing_required))

    ordered_cols = REQUIRED_LONG_COLUMNS + [col for col in long_df.columns if col not in set(REQUIRED_LONG_COLUMNS)]
    long_df = long_df[ordered_cols]

    condition_values_df = condition_values_df.copy()
    condition_values_df = condition_values_df.sort_values("gene_id", kind="mergesort").reset_index(drop=True)

    return long_df, rolling_df, condition_values_df, condition_metadata_df.reset_index(drop=True)


def build_rank_geometry_from_files(
    expr_path: Path | str,
    sample_sheet_path: Path | str,
    config_path: Path | str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, Any]]:
    """File-backed wrapper used by run.py."""
    config = load_yaml(config_path)
    expr_df = read_table(expr_path)
    sample_sheet_df = read_table(sample_sheet_path)
    long_df, rolling_df, condition_values_df, condition_metadata_df = build_rank_geometry_from_frames(
        expr_df=expr_df,
        sample_sheet_df=sample_sheet_df,
        config=config,
    )
    return long_df, rolling_df, condition_values_df, condition_metadata_df, config
