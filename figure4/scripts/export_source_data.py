from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

import pandas as pd

from plot_figure5 import build_panel_source_data
from utils import dump_yaml, ensure_dir, json_dump, sha256_file, write_table


def write_source_data_bundle(
    *,
    long_df: pd.DataFrame,
    rolling_df: pd.DataFrame,
    condition_values_df: pd.DataFrame,
    condition_metadata_df: pd.DataFrame,
    config: Mapping[str, Any],
    outdir: Path | str,
    expr_path: Path | str,
    sample_sheet_path: Path | str,
    config_path: Path | str,
    figure_path: Optional[Path | str] = None,
    source_data_dir: Optional[Path | str] = None,
) -> Dict[str, Path]:
    outdir = ensure_dir(Path(outdir))
    tables_dir = ensure_dir(outdir / "tables")
    manifests_dir = ensure_dir(outdir / "manifests")

    if source_data_dir is None:
        if outdir.name == "outputs":
            source_data_dir = outdir.parent / "source_data"
        else:
            source_data_dir = outdir / "source_data"
    source_data_dir = ensure_dir(Path(source_data_dir))

    paths: Dict[str, Path] = {}
    paths["long_table"] = write_table(long_df, tables_dir / "rank_geometry_long.csv")
    paths["rolling_table"] = write_table(rolling_df, tables_dir / "rank_geometry_rolling_summary.csv")
    paths["condition_values"] = write_table(condition_values_df, tables_dir / "condition_aggregates.csv")
    paths["condition_metadata"] = write_table(condition_metadata_df, tables_dir / "condition_metadata.csv")

    panel_data = build_panel_source_data(long_df, rolling_df, config)
    for panel, df in panel_data.items():
        paths[f"panel_{panel}"] = write_table(df, Path(source_data_dir) / f"Figure5_panel_{panel}_source_data.csv")

    config_snapshot_path = manifests_dir / "rank_geometry_config.snapshot.yaml"
    dump_yaml(config, config_snapshot_path)
    paths["config_snapshot"] = config_snapshot_path

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "expression_matrix": {
                "path": str(Path(expr_path).resolve()),
                "sha256": sha256_file(expr_path),
            },
            "sample_sheet": {
                "path": str(Path(sample_sheet_path).resolve()),
                "sha256": sha256_file(sample_sheet_path),
            },
            "config": {
                "path": str(Path(config_path).resolve()),
                "sha256": sha256_file(config_path),
            },
        },
        "outputs": {
            "figure": str(Path(figure_path).resolve()) if figure_path else None,
            "long_table": str(paths["long_table"].resolve()),
            "rolling_table": str(paths["rolling_table"].resolve()),
            "condition_values": str(paths["condition_values"].resolve()),
            "condition_metadata": str(paths["condition_metadata"].resolve()),
            "source_data_dir": str(Path(source_data_dir).resolve()),
            "source_data_files": {
                panel: str(paths[f"panel_{panel}"].resolve()) for panel in panel_data
            },
        },
        "summary": {
            "n_long_rows": int(len(long_df)),
            "n_rolling_rows": int(len(rolling_df)),
            "baseline_modes": sorted(long_df["baseline_mode"].dropna().unique().tolist()),
            "conditions": sorted(long_df["condition_label"].dropna().unique().tolist()),
        },
        "config_snapshot_inline": dict(config),
    }
    manifest_path = manifests_dir / "run_manifest.json"
    json_dump(manifest, manifest_path)
    paths["manifest"] = manifest_path

    return paths
