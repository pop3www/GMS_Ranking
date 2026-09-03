from __future__ import annotations

# SCRIPT_DIR_BOOTSTRAP
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import argparse
import json
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import pandas as pd

from fig6_ml_core import InputContractError, read_sample_sheet, sha256_file


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build a count-derived expression matrix for Figure 6 from RSEM expected counts. "
            "The output is log2(size-factor-normalized expected count + pseudocount) and is "
            "intended to be used without an additional --log1p transform."
        )
    )
    p.add_argument("--counts", required=True, help="Gene-by-sample RSEM expected-count matrix, CSV/TSV.")
    p.add_argument("--sample-sheet", required=True, help="Sample sheet used to align and validate sample columns.")
    p.add_argument("--output", required=True, help="Output matrix path, CSV or TSV.")
    p.add_argument("--id-col", default=None, help="Gene identifier column. Defaults to first column.")
    p.add_argument("--min-total-expected-count", type=float, default=80.0, help="Keep genes with total expected count >= this threshold. Default: 80.")
    p.add_argument("--pseudocount", type=float, default=1.0, help="Pseudocount added after size-factor normalization before log2. Default: 1.0.")
    p.add_argument("--duplicate-policy", choices=["error", "mean"], default="error", help="How to handle duplicated gene identifiers. Default: error.")
    p.add_argument("--size-factor-method", choices=["median_ratio", "library_size"], default="median_ratio", help="Primary size-factor estimator. Default: median_ratio; falls back to library_size if needed.")
    return p.parse_args()


def _sep_for(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _read_matrix(path: Path, id_col: Optional[str], duplicate_policy: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Counts matrix not found: {path}")
    df = pd.read_csv(path, sep=_sep_for(path))
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 2:
        raise InputContractError("Counts matrix must contain one gene column and at least one sample column.")
    if id_col is None:
        id_col = df.columns[0]
    if id_col not in df.columns:
        raise InputContractError(f"Gene identifier column not found: {id_col}")
    df = df.dropna(subset=[id_col]).copy()
    df[id_col] = df[id_col].astype(str).str.strip()
    sample_cols = [c for c in df.columns if c != id_col]
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if df[sample_cols].isna().any().any():
        n = int(df[sample_cols].isna().sum().sum())
        raise InputContractError(f"Counts matrix contains {n} missing/non-numeric values after parsing.")
    if (df[sample_cols] < 0).any().any():
        raise InputContractError("Counts matrix contains negative values, which are invalid for count-derived normalization.")
    duplicated = df[id_col].duplicated(keep=False)
    if duplicated.any():
        if duplicate_policy == "error":
            raise InputContractError(
                f"Counts matrix contains {int(duplicated.sum())} duplicated gene identifiers. Deduplicate upstream or use --duplicate-policy mean."
            )
        df = df.groupby(id_col, as_index=False)[sample_cols].mean()
    return df.set_index(id_col)


def _library_size_factors(counts: pd.DataFrame) -> pd.Series:
    lib = counts.sum(axis=0).astype(float)
    if (lib <= 0).any():
        bad = lib.index[lib <= 0].tolist()
        raise InputContractError("At least one sample has non-positive total expected count: " + ", ".join(bad))
    return lib / float(np.median(lib))


def _median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    positive = counts.gt(0).all(axis=1)
    usable = counts.loc[positive]
    if usable.shape[0] < 10:
        raise InputContractError("Fewer than 10 genes are positive in all samples; median-ratio size factors are unstable.")
    log_counts = np.log(usable.to_numpy(dtype=float))
    geo_means = np.exp(log_counts.mean(axis=1))
    ratios = usable.divide(geo_means, axis=0)
    sf = ratios.median(axis=0).astype(float)
    if sf.isna().any() or (sf <= 0).any():
        raise InputContractError("Median-ratio size factor estimation produced invalid values.")
    return sf / float(np.median(sf))


def estimate_size_factors(counts: pd.DataFrame, method: str) -> tuple[pd.Series, str]:
    if method == "library_size":
        return _library_size_factors(counts), "library_size"
    try:
        return _median_ratio_size_factors(counts), "median_ratio"
    except Exception as exc:
        sf = _library_size_factors(counts)
        sf.attrs["fallback_reason"] = f"median_ratio_failed: {type(exc).__name__}: {exc}"
        return sf, "library_size_fallback"


def main() -> int:
    args = parse_args()
    counts_path = Path(args.counts)
    out_path = Path(args.output)
    sample_sheet_path = Path(args.sample_sheet)

    counts = _read_matrix(counts_path, args.id_col, args.duplicate_policy)
    meta = read_sample_sheet(sample_sheet_path)
    sample_ids = meta["sample_id"].astype(str).str.strip().tolist()

    missing_in_counts = [s for s in sample_ids if s not in counts.columns]
    extra_in_counts = [s for s in counts.columns.astype(str) if s not in set(sample_ids)]
    if missing_in_counts or extra_in_counts:
        msg = []
        if missing_in_counts:
            msg.append("sample-sheet rows absent from counts: " + ", ".join(missing_in_counts[:10]))
        if extra_in_counts:
            msg.append("counts columns absent from sample sheet: " + ", ".join(extra_in_counts[:10]))
        raise InputContractError("Sample identifiers do not match. " + " | ".join(msg))

    counts = counts.loc[:, sample_ids]
    n_genes_before = int(counts.shape[0])
    total_counts = counts.sum(axis=1)
    keep = total_counts >= float(args.min_total_expected_count)
    counts_filt = counts.loc[keep].copy()
    if counts_filt.empty:
        raise InputContractError("No genes remain after min-total-expected-count filtering.")

    size_factors, sf_method_used = estimate_size_factors(counts_filt, args.size_factor_method)
    norm = counts_filt.divide(size_factors, axis=1)
    transformed = np.log2(norm + float(args.pseudocount))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    gene_col = counts.index.name or "gene"
    out_df = transformed.reset_index().rename(columns={transformed.index.name or "index": gene_col})
    if out_df.columns[0] != gene_col:
        out_df = out_df.rename(columns={out_df.columns[0]: gene_col})
    out_df.to_csv(out_path, sep=_sep_for(out_path), index=False)

    sf_df = pd.DataFrame({"sample_id": size_factors.index, "size_factor": size_factors.values})
    sf_path = out_path.with_suffix(out_path.suffix + ".size_factors.csv")
    sf_df.to_csv(sf_path, index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_counts": str(counts_path.resolve()),
        "source_counts_sha256": sha256_file(counts_path),
        "sample_sheet": str(sample_sheet_path.resolve()),
        "sample_sheet_sha256": sha256_file(sample_sheet_path),
        "output_matrix": str(out_path.resolve()),
        "output_matrix_sha256": sha256_file(out_path),
        "size_factors_file": str(sf_path.resolve()),
        "size_factors_sha256": sha256_file(sf_path),
        "n_genes_before_filter": n_genes_before,
        "n_genes_after_filter": int(counts_filt.shape[0]),
        "n_samples": int(counts_filt.shape[1]),
        "min_total_expected_count": float(args.min_total_expected_count),
        "pseudocount": float(args.pseudocount),
        "size_factor_method_requested": args.size_factor_method,
        "size_factor_method_used": sf_method_used,
        "size_factor_fallback_reason": size_factors.attrs.get("fallback_reason", ""),
        "matrix_scale": "log2_size_factor_normalized_rsem_expected_counts",
        "notes": [
            "Input counts are allowed to be fractional RSEM expected counts.",
            "Output is already log2-transformed; downstream Figure 6 runs should not pass --log1p.",
        ],
    }
    manifest_path = out_path.with_suffix(out_path.suffix + ".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"[OK] Wrote count-derived expression matrix: {out_path}")
    print(f"[OK] Wrote size factors: {sf_path}")
    print(f"[OK] Wrote manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
