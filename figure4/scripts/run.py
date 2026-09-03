from __future__ import annotations

import argparse
import sys
from pathlib import Path


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from build_rank_geometry import build_rank_geometry_from_files
from export_source_data import write_source_data_bundle
from plot_figure5 import plot_figure5
from utils import ensure_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and plot Figure 5 rank geometry.")
    parser.add_argument(
        "--expr",
        required=True,
        help="Path to the canonical gene-by-sample expression matrix (.csv, .tsv/.txt/.tab, or .parquet).",
    )
    parser.add_argument(
        "--sample-sheet",
        required=True,
        help="Path to config/sample_sheet.csv.",
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to config/rank_geometry.yaml.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for figures, tables, and manifests.",
    )
    parser.add_argument(
        "--source-data-dir",
        default=None,
        help="Optional directory for panel-wise source-data CSVs. Defaults to sibling source_data/ when outdir ends with outputs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    outdir = ensure_dir(Path(args.outdir))

    long_df, rolling_df, condition_values_df, condition_metadata_df, config = build_rank_geometry_from_files(
        expr_path=args.expr,
        sample_sheet_path=args.sample_sheet,
        config_path=args.config,
    )

    figure_path = plot_figure5(
        long_df=long_df,
        rolling_df=rolling_df,
        config=config,
        outdir=outdir,
    )

    bundle_paths = write_source_data_bundle(
        long_df=long_df,
        rolling_df=rolling_df,
        condition_values_df=condition_values_df,
        condition_metadata_df=condition_metadata_df,
        config=config,
        outdir=outdir,
        expr_path=args.expr,
        sample_sheet_path=args.sample_sheet,
        config_path=args.config,
        figure_path=figure_path,
        source_data_dir=args.source_data_dir,
    )

    print("Figure 5 rank-geometry bundle complete.")
    print(f"Figure: {figure_path}")
    print(f"Long table: {bundle_paths['long_table']}")
    print(f"Rolling summary: {bundle_paths['rolling_table']}")
    print(f"Manifest: {bundle_paths['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
