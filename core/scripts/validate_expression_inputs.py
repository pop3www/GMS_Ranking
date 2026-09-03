#!/usr/bin/env python3
"""Validate and standardize RSEM gene-level count/TPM matrices for GMS_Ranking.

Typical use from repo root:

  python core/scripts/validate_expression_inputs.py \
    --counts data/raw/rsem/RawCountFile_rsemgenes.txt \
    --tpm data/raw/rsem/TPMCountFile_rsemgenes.txt \
    --sample-sheet config/sample_sheet.csv \
    --out-dir data/processed \
    --write-standardized

The script checks sample-name consistency against config/sample_sheet.csv, gene uniqueness,
non-negative numeric values, count/TPM row/column agreement, and whether the count matrix
contains fractional values. Fractional values are common for RSEM expected counts; they are
valid inputs for edgeR/limma workflows but should be handled explicitly for DESeq2.

By default, simple repeated non-ENSG IDs are normalized in standardized outputs and validation
(e.g., LYT2_LYT2 -> LYT2). Original source files in data/raw/rsem are not modified.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def _normalize_simple_repeated_id(x: str) -> str:
    x = str(x)
    if x.startswith("ENSG"):
        return x
    if "_" not in x:
        return x
    left, right = x.rsplit("_", 1)
    if left and left == right:
        return left
    return x


def _read_matrix(path: Path, normalize_repeated_ids: bool = True) -> Tuple[pd.DataFrame, Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    suffix = path.suffix.lower()
    sep = "," if suffix == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    if df.shape[1] < 2:
        # Retry with automatic delimiter if extension was misleading.
        df = pd.read_csv(path, sep=None, engine="python")
    first = df.columns[0]
    if first != "gene_id":
        df = df.rename(columns={first: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str)

    id_changes: Dict[str, str] = {}
    if normalize_repeated_ids:
        new_ids = []
        for old in df["gene_id"].tolist():
            new = _normalize_simple_repeated_id(old)
            if new != old:
                id_changes[old] = new
            new_ids.append(new)
        df["gene_id"] = new_ids
    return df, id_changes


def _numeric_block(df: pd.DataFrame, label: str) -> pd.DataFrame:
    sample_cols = [c for c in df.columns if c != "gene_id"]
    out = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    if out.isna().any().any():
        bad = int(out.isna().sum().sum())
        raise ValueError(f"{label}: found {bad} non-numeric/NA values in sample columns")
    return out


def _check_matrix(
    df: pd.DataFrame,
    label: str,
    expected_samples: List[str],
) -> Tuple[Dict[str, object], pd.DataFrame]:
    samples = [c for c in df.columns if c != "gene_id"]
    nums = _numeric_block(df, label)
    missing = sorted(set(expected_samples) - set(samples))
    extra = sorted(set(samples) - set(expected_samples))
    duplicated_genes = int(df["gene_id"].duplicated().sum())
    duplicated_samples = int(pd.Index(samples).duplicated().sum())
    finite = bool(np.isfinite(nums.to_numpy(dtype=float)).all())
    min_value = float(nums.min().min())
    max_value = float(nums.max().max())
    negative_values = int((nums < 0).sum().sum())
    fractional_values = int((~np.isclose(nums.to_numpy(dtype=float), np.rint(nums.to_numpy(dtype=float)))).sum())

    summary = {
        "label": label,
        "n_genes": int(df.shape[0]),
        "n_samples": int(len(samples)),
        "sample_set_matches_sample_sheet": len(missing) == 0 and len(extra) == 0,
        "missing_samples": missing,
        "extra_samples": extra,
        "duplicated_genes": duplicated_genes,
        "duplicated_sample_columns": duplicated_samples,
        "all_values_finite": finite,
        "min_value": min_value,
        "max_value": max_value,
        "negative_values": negative_values,
        "fractional_values": fractional_values,
    }
    return summary, nums


def _write_standardized(df: pd.DataFrame, expected_samples: List[str], out_path: Path, sep: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = df[["gene_id"] + expected_samples].copy()
    ordered.to_csv(out_path, sep=sep, index=False)


def _status_from_summary(summary: Dict[str, object], require_integer_counts: bool = False) -> str:
    fail = []
    if not summary["sample_set_matches_sample_sheet"]:
        fail.append("sample-set mismatch")
    if summary["duplicated_genes"]:
        fail.append("duplicated genes")
    if summary["duplicated_sample_columns"]:
        fail.append("duplicated sample columns")
    if not summary["all_values_finite"]:
        fail.append("non-finite values")
    if summary["negative_values"]:
        fail.append("negative values")
    if require_integer_counts and summary["fractional_values"]:
        fail.append("fractional count-like values")
    return "PASS" if not fail else "CHECK: " + "; ".join(fail)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True, help="Raw/RSEM count matrix, genes x samples")
    parser.add_argument("--tpm", type=Path, default=None, help="TPM matrix, genes x samples")
    parser.add_argument("--sample-sheet", type=Path, default=Path("config/sample_sheet.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--write-standardized", action="store_true", help="Write matrices reordered to sample_sheet order")
    parser.add_argument("--require-integer-counts", action="store_true", help="Fail count check if count matrix has fractional values")
    parser.add_argument("--no-normalize-repeated-ids", action="store_true", help="Do not normalize simple repeated non-ENSG IDs such as LYT2_LYT2")
    args = parser.parse_args()

    normalize_ids = not args.no_normalize_repeated_ids
    sample_sheet = pd.read_csv(args.sample_sheet)
    if "sample_id" not in sample_sheet.columns:
        raise ValueError(f"sample sheet lacks sample_id column: {args.sample_sheet}")
    expected_samples = sample_sheet["sample_id"].astype(str).tolist()

    counts, count_id_changes = _read_matrix(args.counts, normalize_repeated_ids=normalize_ids)
    counts_summary, _ = _check_matrix(counts, "counts", expected_samples)
    counts_status = _status_from_summary(counts_summary, require_integer_counts=args.require_integer_counts)

    report: Dict[str, object] = {
        "sample_sheet": str(args.sample_sheet),
        "n_expected_samples": len(expected_samples),
        "normalize_simple_repeated_nonensg_ids": normalize_ids,
        "counts_path": str(args.counts),
        "counts_id_changes": count_id_changes,
        "counts": counts_summary,
        "counts_status": counts_status,
        "notes": [],
    }

    if counts_summary["fractional_values"]:
        report["notes"].append(
            "Count matrix contains fractional values. This is compatible with RSEM expected counts; "
            "for DESeq2, handle integer conversion explicitly in the DE script and document it."
        )
    if count_id_changes:
        report["notes"].append(f"Normalized count gene IDs: {count_id_changes}")

    tpm_summary: Optional[Dict[str, object]] = None
    tpm_status: Optional[str] = None
    if args.tpm is not None:
        tpm, tpm_id_changes = _read_matrix(args.tpm, normalize_repeated_ids=normalize_ids)
        tpm_summary, _ = _check_matrix(tpm, "tpm", expected_samples)
        tpm_status = _status_from_summary(tpm_summary, require_integer_counts=False)
        counts_ids = counts["gene_id"].tolist()
        tpm_ids = tpm["gene_id"].tolist()
        counts_set = set(counts_ids)
        tpm_set = set(tpm_ids)
        report["tpm_path"] = str(args.tpm)
        report["tpm_id_changes"] = tpm_id_changes
        report["tpm"] = tpm_summary
        report["tpm_status"] = tpm_status
        report["counts_tpm_same_gene_order"] = bool(counts_ids == tpm_ids)
        report["counts_tpm_same_gene_set"] = bool(counts_set == tpm_set)
        report["counts_tpm_common_genes"] = int(len(counts_set & tpm_set))
        report["genes_in_counts_not_tpm"] = sorted(counts_set - tpm_set)[:50]
        report["genes_in_tpm_not_counts"] = sorted(tpm_set - counts_set)[:50]
        report["counts_tpm_same_sample_set"] = bool(set(counts.columns) == set(tpm.columns))
        if tpm_id_changes:
            report["notes"].append(f"Normalized TPM gene IDs: {tpm_id_changes}")
        if not report["counts_tpm_same_gene_order"]:
            report["notes"].append("Counts and TPM matrices do not have identical gene order after ID normalization.")
    else:
        tpm = None

    # Overall pass: matrix-level checks must pass. Counts/TPM gene mismatch is reported, not used
    # as a hard failure, because Figure 1/4 can use TPM and DE can use counts independently.
    overall_ok = counts_status == "PASS" and (args.tpm is None or tpm_status == "PASS")
    report["overall_status"] = "PASS" if overall_ok else "CHECK"

    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / "expression_input_validation.json"
    md_path = args.out_dir / "expression_input_validation.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    md_lines = [
        "# Expression input validation",
        "",
        f"Overall status: **{report['overall_status']}**",
        "",
        f"Sample sheet: `{args.sample_sheet}`",
        f"Expected samples: {len(expected_samples)}",
        f"Normalize simple repeated non-ENSG IDs: {normalize_ids}",
        "",
        "## Counts",
        f"Path: `{args.counts}`",
        f"Status: **{counts_status}**",
        f"Genes: {counts_summary['n_genes']}",
        f"Samples: {counts_summary['n_samples']}",
        f"Sample set matches sample sheet: {counts_summary['sample_set_matches_sample_sheet']}",
        f"Duplicated genes: {counts_summary['duplicated_genes']}",
        f"Negative values: {counts_summary['negative_values']}",
        f"Fractional values: {counts_summary['fractional_values']}",
        "",
    ]
    if args.tpm is not None and tpm_summary is not None:
        md_lines += [
            "## TPM",
            f"Path: `{args.tpm}`",
            f"Status: **{tpm_status}**",
            f"Genes: {tpm_summary['n_genes']}",
            f"Samples: {tpm_summary['n_samples']}",
            f"Sample set matches sample sheet: {tpm_summary['sample_set_matches_sample_sheet']}",
            f"Duplicated genes: {tpm_summary['duplicated_genes']}",
            f"Negative values: {tpm_summary['negative_values']}",
            f"Counts/TPM same gene order: {report['counts_tpm_same_gene_order']}",
            f"Counts/TPM same gene set: {report['counts_tpm_same_gene_set']}",
            f"Counts/TPM common genes: {report['counts_tpm_common_genes']}",
            "",
        ]
    if report["notes"]:
        md_lines += ["## Notes"] + [f"- {x}" for x in report["notes"]] + [""]
    md_path.write_text("\n".join(md_lines) + "\n")

    if args.write_standardized:
        _write_standardized(counts, expected_samples, args.out_dir / "RawCountFile_rsemgenes.tsv", sep="\t")
        _write_standardized(counts, expected_samples, args.out_dir / "raw_counts_rsemgenes.tsv", sep="\t")
        if args.tpm is not None and tpm is not None:
            _write_standardized(tpm, expected_samples, args.out_dir / "TPMCountFile_rsemgenes.tsv", sep="\t")
            _write_standardized(tpm, expected_samples, args.out_dir / "TPMCountFile_rsemgenes.csv", sep=",")
            _write_standardized(tpm, expected_samples, args.out_dir / "tpm_rsemgenes.tsv", sep="\t")

    print(f"Validation report: {md_path}")
    print(f"Machine-readable report: {json_path}")
    print(f"Overall status: {report['overall_status']}")
    if counts_summary["fractional_values"]:
        print("Note: count matrix has fractional values; treat as RSEM expected counts.")
    return 0 if report["overall_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
