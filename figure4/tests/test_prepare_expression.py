from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PREP_PATH = ROOT / ".." / "core" / "scripts" / "prepare_rank_geometry_expression.py"
# In installed repo layout tests run from repo root, so use direct repo-relative fallback.
if not PREP_PATH.resolve().exists():
    PREP_PATH = Path("core/scripts/prepare_rank_geometry_expression.py")

spec = importlib.util.spec_from_file_location("prepare_rank_geometry_expression", PREP_PATH)
prep = importlib.util.module_from_spec(spec)
assert spec is not None and spec.loader is not None
spec.loader.exec_module(prep)


def test_median_ratio_size_factors_recover_known_library_scaling() -> None:
    counts = pd.DataFrame(
        {
            "s1": [100.0, 200.0, 300.0, 400.0],
            "s2": [200.0, 400.0, 600.0, 800.0],
        }
    )
    sf, method = prep.median_ratio_size_factors(counts)
    assert method == "median_ratio_all_positive_genes"
    ratio = sf["s2"] / sf["s1"]
    assert np.isclose(ratio, 2.0)


def test_prepare_script_filters_counts_and_writes_normalized_matrix(tmp_path, monkeypatch) -> None:
    expr = pd.DataFrame(
        {
            "gene_id": ["low", "g1", "g2"],
            "s1": [10.0, 100.0, 200.0],
            "s2": [10.0, 200.0, 400.0],
        }
    )
    samples = pd.DataFrame({"sample_id": ["s1", "s2"], "group_label": ["A", "A"]})
    expr_path = tmp_path / "raw_counts.tsv"
    sample_path = tmp_path / "sample_sheet.csv"
    out_path = tmp_path / "canonical.csv"
    sf_path = tmp_path / "sf.csv"
    manifest_path = tmp_path / "manifest.json"
    expr.to_csv(expr_path, sep="\t", index=False)
    samples.to_csv(sample_path, index=False)

    argv = [
        "prepare_rank_geometry_expression.py",
        "--input",
        str(expr_path),
        "--sample-sheet",
        str(sample_path),
        "--output",
        str(out_path),
        "--matrix-type",
        "raw_counts",
        "--normalization",
        "median_ratio",
        "--transform",
        "log2p1",
        "--min-total-count",
        "80",
        "--write-size-factors",
        str(sf_path),
        "--write-manifest",
        str(manifest_path),
    ]
    monkeypatch.setattr("sys.argv", argv)
    assert prep.main() == 0

    out = pd.read_csv(out_path)
    assert out["gene_id"].tolist() == ["g1", "g2"]
    assert sf_path.exists()
    assert manifest_path.exists()
    # s2 is exactly double s1 before normalization; after median-ratio normalization
    # and the manuscript-scale log2(x + 1) transform, retained genes should match
    # between samples.
    assert np.allclose(out["s1"].to_numpy(), out["s2"].to_numpy())
    assert np.all(out[["s1", "s2"]].to_numpy() < 10.0)  # log2-normalized scale, not raw counts.
    manifest = json.loads(manifest_path.read_text())
    assert manifest["transform"] == "log2p1"
    assert manifest["scale_label"] == "log2(size-factor-normalized RSEM expected count + 1)"
    assert (out[["s1", "s2"]].to_numpy() < 20).all()
    manifest = json.loads(manifest_path.read_text())
    assert manifest["transform"] == "log2p1"
    assert manifest["scale_label"] == "log2(size-factor-normalized RSEM expected count + 1)"
