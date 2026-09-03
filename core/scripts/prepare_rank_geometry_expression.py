#!/usr/bin/env python3
"""Prepare the canonical Figure 5 gene-by-sample expression matrix.

Figure 5 should be built from the same count-scale expression convention used
in the manuscript: filtered RSEM expected counts normalized by
median-ratio size factors and represented on the manuscript scale
log2(normalized expected count + 1). TPM remains useful as an optional
sensitivity input, but it is not the default primary matrix for the
rank-geometry panels.

Primary live-repo convention:
    data/processed/raw_counts_rsemgenes.tsv

Validated raw inputs, if processed files must be regenerated:
    data/raw/RawCountFile_rsemgenes.txt
    data/raw/TPMCountFile_rsemgenes.txt

Example, from repo root:
    python core/scripts/prepare_rank_geometry_expression.py \
      --input data/processed/raw_counts_rsemgenes.tsv \
      --sample-sheet config/sample_sheet.csv \
      --output core/outputs/tables/rank_geometry_expression_gene_by_sample.csv \
      --matrix-type raw_counts \
      --normalization median_ratio \
      --transform log2p1 \
      --min-total-count 80
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd


VALIDATION_HELP = """
Expected primary Figure 5 input after validation is now:
  data/processed/raw_counts_rsemgenes.tsv

If it is missing, regenerate standardized processed expression files from the
validated raw files with:
  python core/scripts/validate_expression_inputs.py \\
    --counts data/raw/RawCountFile_rsemgenes.txt \\
    --tpm data/raw/TPMCountFile_rsemgenes.txt \\
    --sample-sheet config/sample_sheet.csv \\
    --out-dir data/processed \\
    --write-standardized

Important: the raw files are in data/raw/, not a nested raw RSEM subdirectory.
TPM can be supplied as an explicit sensitivity input with --matrix-type tpm
--normalization none, but the manuscript uses log2(size-factor-normalized
RSEM expected count + 1) for the primary distribution/rank-space workflow.
""".strip()


MatrixType = Literal["auto", "raw_counts", "normalized_counts", "tpm"]
Normalization = Literal["auto", "none", "median_ratio"]
Transform = Literal["auto", "none", "log2p1"]


def _clean_column_name(name: Any) -> str:
    return str(name).replace("\ufeff", "").strip()


def read_table_bom_safe(path: Path) -> pd.DataFrame:
    """Read CSV/TSV/TXT tables with BOM-safe column cleanup."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input table not found: {path}\n\n{VALIDATION_HELP}")

    name = path.name.lower()
    if name.endswith((".csv", ".csv.gz")):
        df = pd.read_csv(path, encoding="utf-8-sig")
    elif name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz", ".tab", ".tab.gz")):
        df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    elif name.endswith(".parquet"):
        df = pd.read_parquet(path)
    else:
        raise ValueError(
            f"Unsupported table format for {path}. Use .csv, .tsv/.txt/.tab, gzip-compressed text, or .parquet."
        )
    df.columns = [_clean_column_name(col) for col in df.columns]
    return df


def infer_matrix_type(path: Path, explicit: MatrixType) -> str:
    if explicit != "auto":
        return explicit
    name = path.name.lower()
    if "tpm" in name:
        return "tpm"
    if "raw" in name or "count" in name or "expected" in name:
        return "raw_counts"
    return "normalized_counts"


def resolve_normalization(matrix_type: str, explicit: Normalization) -> str:
    if explicit != "auto":
        return explicit
    if matrix_type == "raw_counts":
        return "median_ratio"
    return "none"


def resolve_transform(matrix_type: str, normalization: str, explicit: Transform) -> str:
    """Resolve the output transform for the canonical matrix.

    Rank ordering is invariant to log2(x + 1), but the manuscript
    describes the rank/distribution input as log2(size-factor-normalized RSEM
    expected count + 1).  The default therefore writes this manuscript-scale
    matrix for primary raw-count inputs, while keeping TPM/other sensitivity
    matrices unchanged unless explicitly requested.
    """
    if explicit != "auto":
        return explicit
    if matrix_type == "raw_counts" and normalization == "median_ratio":
        return "log2p1"
    return "none"


def coerce_numeric_matrix(df: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame:
    values = df[sample_ids].apply(pd.to_numeric, errors="coerce")
    bad = values.isna().sum().sum()
    if bad:
        raise ValueError(f"Expression matrix contains {bad} missing or non-numeric values in sample columns.")
    if (values < 0).any().any():
        raise ValueError("Expression matrix contains negative values in sample columns.")
    return values.astype(float)


def median_ratio_size_factors(values: pd.DataFrame) -> tuple[pd.Series, str]:
    """Estimate DESeq2-style median-ratio size factors.

    The primary path uses genes with positive counts in all samples, matching the
    usual median-ratio idea. If no such genes exist, the function falls back to a
    positive-count geometric mean and sample-wise positive-ratio medians. The
    size factors are normalized to geometric mean 1 for stable reporting.
    """
    arr = values.to_numpy(dtype=float)
    sample_ids = list(values.columns)

    if arr.shape[0] == 0 or arr.shape[1] == 0:
        raise ValueError("Cannot estimate size factors from an empty matrix.")

    positive_all = np.all(arr > 0, axis=1)
    method = "median_ratio_all_positive_genes"

    if positive_all.any():
        used = arr[positive_all, :]
        geom = np.exp(np.mean(np.log(used), axis=1))
        ratios = used / geom[:, None]
        size_factors = np.median(ratios, axis=0)
    else:
        # Rare fallback for extremely sparse matrices.
        method = "median_ratio_poscounts_fallback"
        positive = arr > 0
        n_pos = positive.sum(axis=1)
        use = n_pos > 0
        if not use.any():
            raise ValueError("Cannot estimate size factors: all expression values are zero.")
        logs = np.where(positive[use, :], np.log(np.where(positive[use, :], arr[use, :], 1.0)), np.nan)
        geom = np.exp(np.nanmean(logs, axis=1))
        ratios = arr[use, :] / geom[:, None]
        ratios[ratios <= 0] = np.nan
        size_factors = np.nanmedian(ratios, axis=0)

    if not np.all(np.isfinite(size_factors)) or np.any(size_factors <= 0):
        method = method + "+library_size_fallback"
        totals = arr.sum(axis=0)
        if np.any(totals <= 0):
            raise ValueError("Cannot estimate fallback library-size factors: at least one sample has total count <= 0.")
        size_factors = totals / np.exp(np.mean(np.log(totals)))
    else:
        size_factors = size_factors / np.exp(np.mean(np.log(size_factors)))

    return pd.Series(size_factors, index=sample_ids, name="size_factor"), method


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare canonical Figure 5 expression matrix.")
    parser.add_argument(
        "--input",
        default="data/processed/raw_counts_rsemgenes.tsv",
        help=(
            "Input gene-by-sample expression table. Default uses validated processed RSEM expected counts; "
            "TPM should be supplied only as an explicit sensitivity input."
        ),
    )
    parser.add_argument("--sample-sheet", default="config/sample_sheet.csv", help="Sample sheet with sample_id column.")
    parser.add_argument(
        "--output",
        default="core/outputs/tables/rank_geometry_expression_gene_by_sample.csv",
        help="Output clean expression CSV.",
    )
    parser.add_argument("--gene-id-column", default="gene_id", help="Gene ID column after BOM cleanup.")
    parser.add_argument("--sample-id-column", default="sample_id", help="Sample ID column in the sample sheet.")
    parser.add_argument(
        "--matrix-type",
        choices=["auto", "raw_counts", "normalized_counts", "tpm"],
        default="auto",
        help="Input scale. auto infers from filename. Primary Figure 5 should use raw_counts.",
    )
    parser.add_argument(
        "--normalization",
        choices=["auto", "none", "median_ratio"],
        default="auto",
        help="Normalization applied before writing canonical matrix. auto applies median_ratio to raw_counts only.",
    )
    parser.add_argument(
        "--min-total-count",
        type=float,
        default=80.0,
        help="Filter genes with total input counts below this value when matrix-type is raw_counts. Use 0 to disable.",
    )
    parser.add_argument(
        "--transform",
        choices=["auto", "none", "log2p1"],
        default="auto",
        help=(
            "Transform written expression values after filtering/normalization. "
            "auto writes log2(normalized count + 1) for the primary raw-count path, "
            "matching the frozen manuscript, and leaves non-count sensitivity inputs unchanged."
        ),
    )
    parser.add_argument(
        "--write-size-factors",
        default="core/outputs/tables/rank_geometry_size_factors.csv",
        help="Where to write size factors/normalization metadata. Set to empty string to skip.",
    )
    parser.add_argument(
        "--write-manifest",
        default="core/outputs/tables/rank_geometry_expression_manifest.json",
        help="Where to write expression-preparation manifest. Set to empty string to skip.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    sample_sheet_path = Path(args.sample_sheet)
    output_path = Path(args.output)

    expr = read_table_bom_safe(input_path)
    sample_sheet = read_table_bom_safe(sample_sheet_path)

    if args.gene_id_column not in expr.columns:
        raise ValueError(
            f"Expression matrix is missing gene ID column '{args.gene_id_column}'. "
            f"Columns start: {expr.columns[:10].tolist()}"
        )
    if args.sample_id_column not in sample_sheet.columns:
        raise ValueError(f"Sample sheet is missing sample ID column '{args.sample_id_column}'.")

    expr[args.gene_id_column] = expr[args.gene_id_column].astype(str)
    if expr[args.gene_id_column].duplicated().any():
        examples = expr.loc[expr[args.gene_id_column].duplicated(), args.gene_id_column].head(10).tolist()
        raise ValueError(f"Duplicate gene IDs found in expression matrix; examples: {examples}")

    sample_ids = sample_sheet[args.sample_id_column].astype(str).tolist()
    missing = [sid for sid in sample_ids if sid not in expr.columns]
    if missing:
        raise ValueError("Sample IDs in sample sheet but missing from expression matrix: " + ", ".join(missing[:30]))

    matrix_type = infer_matrix_type(input_path, args.matrix_type)
    normalization = resolve_normalization(matrix_type, args.normalization)
    transform = resolve_transform(matrix_type, normalization, args.transform)

    values = coerce_numeric_matrix(expr, sample_ids)
    n_genes_initial = int(expr.shape[0])

    filter_applied = False
    if matrix_type == "raw_counts" and args.min_total_count and args.min_total_count > 0:
        total = values.sum(axis=1)
        keep = total >= float(args.min_total_count)
        expr = expr.loc[keep].reset_index(drop=True)
        values = values.loc[keep].reset_index(drop=True)
        filter_applied = True

    size_factor_method = "none"
    size_factors = pd.Series(1.0, index=sample_ids, name="size_factor")
    if normalization == "median_ratio":
        size_factors, size_factor_method = median_ratio_size_factors(values)
        values = values.divide(size_factors, axis=1)
    elif normalization == "none":
        pass
    else:
        raise ValueError(f"Unsupported normalization: {normalization}")

    if transform == "log2p1":
        if (values < 0).any().any():
            raise ValueError("Cannot apply log2p1 transform to negative values.")
        values = np.log2(values + 1.0)
    elif transform == "none":
        pass
    else:
        raise ValueError(f"Unsupported transform: {transform}")

    out = pd.concat([expr[[args.gene_id_column]].reset_index(drop=True), values.reset_index(drop=True)], axis=1)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)

    if args.write_size_factors:
        sf_path = Path(args.write_size_factors)
        sf_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            {
                "sample_id": sample_ids,
                "size_factor": [float(size_factors[sid]) for sid in sample_ids],
                "normalization": normalization,
                "size_factor_method": size_factor_method,
            }
        ).to_csv(sf_path, index=False)
    else:
        sf_path = None

    manifest = {
        "input": str(input_path),
        "sample_sheet": str(sample_sheet_path),
        "output": str(output_path),
        "matrix_type": matrix_type,
        "normalization": normalization,
        "size_factor_method": size_factor_method,
        "min_total_count": float(args.min_total_count),
        "filter_applied": bool(filter_applied),
        "transform": transform,
        "n_genes_initial": n_genes_initial,
        "n_genes_written": int(out.shape[0]),
        "n_samples": len(sample_ids),
        "scale_label": (
            "log2(size-factor-normalized RSEM expected count + 1)"
            if matrix_type == "raw_counts" and normalization == "median_ratio" and transform == "log2p1"
            else (
                "median-ratio-size-factor-normalized RSEM expected count"
                if matrix_type == "raw_counts" and normalization == "median_ratio"
                else matrix_type
            )
        ),
        "size_factors_file": str(sf_path) if sf_path else None,
    }
    if args.write_manifest:
        manifest_path = Path(args.write_manifest)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"Input: {input_path}")
    print(f"Matrix type: {matrix_type}")
    print(f"Normalization: {normalization} ({size_factor_method})")
    print(f"Transform: {transform}")
    print(f"Wrote: {output_path}")
    print(f"Genes initial: {n_genes_initial}")
    print(f"Genes written: {out.shape[0]}")
    print(f"Samples: {len(sample_ids)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
