from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
import yaml


def ensure_dir(path: Path | str) -> Path:
    """Create a directory if needed and return it as a Path."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _clean_column_name(name: Any) -> str:
    """Normalize input column names without changing biological/sample labels."""
    return str(name).replace("\ufeff", "").strip()


def read_table(path: Path | str) -> pd.DataFrame:
    """Read a table with light format inference and BOM-safe column handling.

    Supports .csv, .tsv/.txt/.tab, .parquet, and gzip-compressed text tables
    such as .csv.gz or .tsv.gz. The first column in the live TPM file carries
    a UTF-8 BOM in some exports, so column names are cleaned after import.
    """
    path = Path(path)
    name = path.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        df = pd.read_csv(path, encoding="utf-8-sig")
    elif name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz", ".tab", ".tab.gz")):
        df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    elif name.endswith(".parquet"):
        df = pd.read_parquet(path)
    else:
        raise ValueError(
            f"Unsupported table format for {path}. Use .csv, .csv.gz, .tsv/.txt/.tab, .tsv.gz/.txt.gz/.tab.gz, or .parquet."
        )
    df.columns = [_clean_column_name(col) for col in df.columns]
    return df


def write_table(df: pd.DataFrame, path: Path | str) -> Path:
    """Write a dataframe to CSV and return the output path."""
    path = Path(path)
    ensure_dir(path.parent)
    df.to_csv(path, index=False)
    return path


def load_yaml(path: Path | str) -> Dict[str, Any]:
    """Load a YAML file into a dict."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return data or {}


def dump_yaml(data: Mapping[str, Any], path: Path | str) -> Path:
    """Write a mapping to YAML and return the output path."""
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(dict(data), handle, sort_keys=False)
    return path


def sha256_file(path: Path | str, chunk_size: int = 1024 * 1024) -> str:
    """Compute a SHA256 hash for a file."""
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def canonicalize_label(label: Any, label_map: Optional[Mapping[str, str]] = None) -> str:
    """Map a raw label to a manuscript label when a mapping is supplied."""
    if pd.isna(label):
        return ""
    label_str = str(label)
    if label_map is None:
        return label_str
    return str(label_map.get(label_str, label_str))


def canonicalize_mapping(
    mapping: Mapping[str, str],
    label_map: Optional[Mapping[str, str]] = None,
) -> Dict[str, str]:
    """Apply label canonicalization to both keys and values of a mapping."""
    out: Dict[str, str] = {}
    for key, value in mapping.items():
        out[canonicalize_label(key, label_map)] = canonicalize_label(value, label_map)
    return out


def stable_desc_rank(
    df: pd.DataFrame,
    gene_col: str,
    value_col: str,
    rank_col: str,
) -> pd.DataFrame:
    """
    Deterministically rank genes in descending expression order.

    Ties are broken by gene_id ascending after a stable sort.
    Rank 1 corresponds to the expression head.
    """
    out = df[[gene_col, value_col]].copy()
    out[gene_col] = out[gene_col].astype(str)
    out = out.sort_values(
        by=[value_col, gene_col],
        ascending=[False, True],
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)
    out[rank_col] = np.arange(1, len(out) + 1, dtype=int)
    return out


def add_baseline_bands(
    df: pd.DataFrame,
    baseline_rank_col: str = "baseline_rank",
    head_pct: float = 10.0,
    mid_pct: float = 50.0,
    tail_pct: float = 10.0,
) -> pd.DataFrame:
    """
    Add baseline percentiles and manuscript bands using the fixed baseline ranking.

    baseline_percentile is measured from the head:
    - low percentile = high expression
    - head = top head_pct
    - mid  = middle mid_pct centered on the median
    - tail = bottom tail_pct
    - everything else = other
    """
    out = df.copy()
    n = len(out)
    if n == 0:
        out["baseline_percentile"] = pd.Series(dtype=float)
        out["band"] = pd.Series(dtype=object)
        return out

    out["baseline_percentile"] = ((out[baseline_rank_col].astype(float) - 0.5) / float(n)) * 100.0

    mid_lower = (100.0 - float(mid_pct)) / 2.0
    mid_upper = 100.0 - mid_lower

    conditions = [
        out["baseline_percentile"] <= float(head_pct),
        (out["baseline_percentile"] >= mid_lower) & (out["baseline_percentile"] <= mid_upper),
        out["baseline_percentile"] >= (100.0 - float(tail_pct)),
    ]
    choices = ["head", "mid", "tail"]
    out["band"] = np.select(conditions, choices, default="other")
    return out


def compute_rolling_summary(
    long_df: pd.DataFrame,
    windows: Sequence[int],
    *,
    center: bool = True,
    min_periods: Optional[int] = None,
    group_cols: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Compute rolling median / q25 / q75 summaries on baseline rank for each comparison."""
    if group_cols is None:
        group_cols = ["baseline_mode", "baseline_label", "condition_label", "comparison_id"]

    windows = [int(w) for w in windows]
    frames: List[pd.DataFrame] = []

    if long_df.empty:
        return pd.DataFrame(
            columns=list(group_cols)
            + ["window", "baseline_rank", "median", "q25", "q75", "n_genes", "time_h", "dox", "tam", "cpt_level"]
        )

    meta_cols = [col for col in ["time_h", "dox", "tam", "cpt_level"] if col in long_df.columns]

    for keys, group in long_df.groupby(list(group_cols), sort=False):
        group = group.sort_values("baseline_rank").reset_index(drop=True)
        meta = {col: group[col].iloc[0] for col in meta_cols}
        for window in windows:
            mp = int(min_periods) if min_periods is not None else int(window)
            roller = group["rank_shift"].rolling(window=window, center=center, min_periods=mp)
            block = pd.DataFrame(
                {
                    "baseline_rank": group["baseline_rank"].to_numpy(),
                    "window": int(window),
                    "median": roller.median().to_numpy(),
                    "q25": roller.quantile(0.25).to_numpy(),
                    "q75": roller.quantile(0.75).to_numpy(),
                    "n_genes": len(group),
                }
            )
            if isinstance(keys, tuple):
                for col, value in zip(group_cols, keys):
                    block[col] = value
            else:
                block[group_cols[0]] = keys
            for col, value in meta.items():
                block[col] = value
            frames.append(block)

    order = list(group_cols) + ["window", "baseline_rank", "median", "q25", "q75", "n_genes"] + meta_cols
    return pd.concat(frames, ignore_index=True)[order]


def condition_color_map(config: Mapping[str, Any]) -> Dict[str, str]:
    """Get a condition palette from config when present."""
    return dict(config.get("figure", {}).get("condition_palette", {}))


def json_dump(data: Mapping[str, Any], path: Path | str) -> Path:
    """Write JSON with pretty formatting."""
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(dict(data), handle, indent=2, sort_keys=False)
    return path
