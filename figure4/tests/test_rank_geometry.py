from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


THIS_DIR = Path(__file__).resolve().parent
MODULE_DIR = THIS_DIR.parent / "scripts"
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

from build_rank_geometry import build_rank_geometry_from_frames
from utils import add_baseline_bands, stable_desc_rank


def _toy_config() -> dict:
    return {
        "expression": {
            "gene_id_column": "gene_id",
        },
        "sample_sheet": {
            "sample_id_column": "sample_id",
            "condition_label_column": "condition_label",
            "time_h_column": "time_h",
            "dox_column": "dox",
            "tam_column": "tam",
            "cpt_level_column": "cpt_level",
        },
        "analysis": {
            "aggregate": "mean",
            "condition_order": ["Early", "Late"],
            "label_map": {},
            "band_percentiles": {
                "head": 10,
                "mid": 50,
                "tail": 10,
            },
        },
        "baseline": {
            "enabled_modes": ["global"],
            "default_mode": "global",
            "global": {
                "label": "Ctrl",
            },
        },
        "rolling": {
            "windows": [4],
            "default_window": 4,
            "center": True,
            "min_periods": 4,
        },
        "figure": {
            "baseline_mode": "global",
            "panel_a": {"condition_label": "Early", "title": "Early"},
            "panel_b": {"condition_label": "Late", "title": "Late"},
            "profile_conditions": ["Early", "Late"],
            "profile_windows": [4, 4],
        },
        "robustness": {"panels": []},
        "outputs": {"figure_prefix": "Figure5_rank_geometry"},
    }


def _toy_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = [f"g{i:02d}" for i in range(1, 21)]
    ctrl = [200 - 10 * (i - 1) for i in range(1, 21)]
    early = ctrl.copy()
    late = ctrl.copy()

    # Preserve early ordering with mild value changes.
    early[0] = 199
    early[1] = 188
    early[2] = 182
    early[3] = 171

    # Introduce late apex clipping at the head.
    late[0] = 145
    late[1] = 144
    late[2] = 143
    late[3] = 142
    late[4] = 170
    late[5] = 169
    late[6] = 168
    late[7] = 167

    expr = pd.DataFrame(
        {
            "gene_id": genes,
            "ctrl_r1": ctrl,
            "ctrl_r2": ctrl,
            "early_r1": early,
            "early_r2": early,
            "late_r1": late,
            "late_r2": late,
        }
    )

    sample_sheet = pd.DataFrame(
        {
            "sample_id": ["ctrl_r1", "ctrl_r2", "early_r1", "early_r2", "late_r1", "late_r2"],
            "condition_label": ["Ctrl", "Ctrl", "Early", "Early", "Late", "Late"],
            "time_h": [0, 0, 4, 4, 24, 24],
            "dox": [0, 0, 0, 0, 0, 0],
            "tam": [0, 0, 1, 1, 1, 1],
            "cpt_level": ["none", "none", "none", "none", "none", "none"],
        }
    )
    return expr, sample_sheet


def test_stable_desc_rank_breaks_ties_by_gene_id() -> None:
    df = pd.DataFrame(
        {
            "gene_id": ["B", "A", "C"],
            "value": [10, 10, 9],
        }
    )
    ranked = stable_desc_rank(df, gene_col="gene_id", value_col="value", rank_col="rank")
    assert ranked["gene_id"].tolist() == ["A", "B", "C"]
    assert ranked["rank"].tolist() == [1, 2, 3]


def test_add_baseline_bands_is_fixed_on_baseline_rank() -> None:
    df = pd.DataFrame({"baseline_rank": list(range(1, 21))})
    banded = add_baseline_bands(df, baseline_rank_col="baseline_rank", head_pct=10, mid_pct=50, tail_pct=10)

    head_genes = banded.loc[banded["band"] == "head", "baseline_rank"].tolist()
    mid_genes = banded.loc[banded["band"] == "mid", "baseline_rank"].tolist()
    tail_genes = banded.loc[banded["band"] == "tail", "baseline_rank"].tolist()

    assert head_genes == [1, 2]
    assert mid_genes == list(range(6, 16))
    assert tail_genes == [19, 20]


def test_rank_shift_sign_convention_is_baseline_minus_condition() -> None:
    expr = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "ctrl_r1": [10, 9, 8],
            "ctrl_r2": [10, 9, 8],
            "early_r1": [8, 10, 9],
            "early_r2": [8, 10, 9],
        }
    )
    sample_sheet = pd.DataFrame(
        {
            "sample_id": ["ctrl_r1", "ctrl_r2", "early_r1", "early_r2"],
            "condition_label": ["Ctrl", "Ctrl", "Early", "Early"],
            "time_h": [0, 0, 4, 4],
            "dox": [0, 0, 0, 0],
            "tam": [0, 0, 1, 1],
            "cpt_level": ["none", "none", "none", "none"],
        }
    )
    config = _toy_config()
    config["analysis"]["condition_order"] = ["Early"]
    long_df, _, _, _ = build_rank_geometry_from_frames(expr, sample_sheet, config)

    shifts = long_df.set_index("gene_id")["rank_shift"].to_dict()
    assert shifts["g1"] == -2
    assert shifts["g2"] == 1
    assert shifts["g3"] == 1


def test_late_apex_clipping_is_more_negative_in_the_head_than_early_state() -> None:
    expr, sample_sheet = _toy_frames()
    config = _toy_config()

    long_df, rolling_df, _, _ = build_rank_geometry_from_frames(expr, sample_sheet, config)

    early_head = long_df.loc[(long_df["condition_label"] == "Early") & (long_df["band"] == "head"), "rank_shift"].median()
    late_head = long_df.loc[(long_df["condition_label"] == "Late") & (long_df["band"] == "head"), "rank_shift"].median()

    assert early_head >= late_head
    assert late_head < 0

    early_roll = rolling_df.loc[(rolling_df["condition_label"] == "Early") & (rolling_df["baseline_rank"] <= 4), "median"].dropna().mean()
    late_roll = rolling_df.loc[(rolling_df["condition_label"] == "Late") & (rolling_df["baseline_rank"] <= 4), "median"].dropna().mean()

    assert late_roll < early_roll


def test_named_baseline_mode_builds_additional_assignments() -> None:
    expr, sample_sheet = _toy_frames()
    config = _toy_config()
    config["analysis"]["condition_order"] = ["Early"]
    config["baseline"] = {
        "enabled_modes": ["global", "unprimed_matched"],
        "default_mode": "global",
        "global": {"label": "Ctrl"},
        "named_modes": {
            "unprimed_matched": {
                "mapping": {
                    "Early": "Ctrl",
                }
            }
        },
    }
    long_df, _, _, _ = build_rank_geometry_from_frames(expr, sample_sheet, config)
    assert sorted(long_df["baseline_mode"].unique().tolist()) == ["global", "unprimed_matched"]


def test_with_dox_without_dox_panel_selects_matched_contrasts() -> None:
    from plot_figure5 import build_panel_source_data

    expr = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "ctrl_r1": [10, 8, 6, 4],
            "d_r1": [9, 8, 7, 4],
            "tam_r1": [8, 10, 6, 4],
            "dt_r1": [7, 9, 8, 4],
        }
    )
    sample_sheet = pd.DataFrame(
        {
            "sample_id": ["ctrl_r1", "d_r1", "tam_r1", "dt_r1"],
            "condition_label": ["Ctrl", "D", "Tam", "DT"],
            "time_h": [4, 4, 4, 4],
            "dox": [0, 1, 0, 1],
            "tam": [0, 0, 1, 1],
            "cpt_level": ["none", "none", "none", "none"],
        }
    )
    config = _toy_config()
    config["analysis"]["condition_order"] = ["Tam", "DT"]
    config["baseline"] = {
        "enabled_modes": ["matched"],
        "default_mode": "matched",
        "matched": {"mapping": {"Tam": "Ctrl", "DT": "D"}},
    }
    config["figure"] = {
        "baseline_mode": "matched",
        "panel_a": {"condition_label": "Tam", "baseline_label": "Ctrl", "title": "Tam"},
        "panel_b": {"condition_label": "DT", "baseline_label": "D", "title": "DT"},
        "profile_conditions": ["Tam", "DT"],
        "profile_windows": [4, 4],
    }
    config["robustness"] = {
        "panels": [
            {
                "panel": "E",
                "x": {
                    "condition_label": "Tam",
                    "baseline_mode": "matched",
                    "baseline_label": "Ctrl",
                    "axis_label": "without Dox (Tam vs Ctrl)",
                },
                "y": {
                    "condition_label": "DT",
                    "baseline_mode": "matched",
                    "baseline_label": "D",
                    "axis_label": "with Dox (DT vs D)",
                },
            }
        ]
    }

    long_df, rolling_df, _, _ = build_rank_geometry_from_frames(expr, sample_sheet, config)
    panel_data = build_panel_source_data(long_df, rolling_df, config)
    panel_e = panel_data["E"]

    assert set(panel_e["x_baseline_label"].unique()) == {"Ctrl"}
    assert set(panel_e["x_condition_label"].unique()) == {"Tam"}
    assert set(panel_e["y_baseline_label"].unique()) == {"D"}
    assert set(panel_e["y_condition_label"].unique()) == {"DT"}
    assert panel_e["x_axis_label"].iloc[0] == "without Dox (Tam vs Ctrl)"
    assert panel_e["y_axis_label"].iloc[0] == "with Dox (DT vs D)"


def test_live_repo_sample_sheet_columns_can_be_mapped_to_manuscript_labels() -> None:
    expr = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "1_4_ctrl_RP1": [10, 8, 6, 4],
            "2_4_Tam_RP1": [8, 10, 6, 4],
            "5_4D_RP1": [9, 8, 7, 4],
            "6_4DT_RP1": [7, 9, 8, 4],
        }
    )
    sample_sheet = pd.DataFrame(
        {
            "sample_id": ["1_4_ctrl_RP1", "2_4_Tam_RP1", "5_4D_RP1", "6_4DT_RP1"],
            "group_label": ["4_ctrl", "4_Tam", "4D", "4DT"],
            "time_h": [4, 4, 4, 4],
            "dox": ["None_Dox", "None_Dox", "Dox", "Dox"],
            "myc": ["None_MYC", "MYC", "None_MYC", "MYC"],
            "cpt_level": ["None_CPT", "None_CPT", "None_CPT", "None_CPT"],
        }
    )
    config = _toy_config()
    config["sample_sheet"] = {
        "sample_id_column": "sample_id",
        "condition_label_column": "group_label",
        "time_h_column": "time_h",
        "dox_column": "dox",
        "tam_column": "myc",
        "cpt_level_column": "cpt_level",
    }
    config["analysis"]["label_map"] = {
        "4_ctrl": "4 h Ctrl",
        "4_Tam": "4 h Tam",
        "4D": "4 h D",
        "4DT": "4 h DT",
    }
    config["analysis"]["condition_order"] = ["4 h Tam", "4 h DT"]
    config["baseline"] = {
        "enabled_modes": ["matched"],
        "default_mode": "matched",
        "matched": {"mapping": {"4 h Tam": "4 h Ctrl", "4 h DT": "4 h D"}},
    }
    long_df, _, condition_values, metadata = build_rank_geometry_from_frames(expr, sample_sheet, config)
    assert {"4 h Ctrl", "4 h D", "4 h Tam", "4 h DT"}.issubset(set(condition_values.columns))
    assert set(long_df["condition_label"].unique()) == {"4 h Tam", "4 h DT"}
    assert set(long_df["baseline_label"].unique()) == {"4 h Ctrl", "4 h D"}
    assert metadata.loc[metadata["condition_label"] == "4 h DT", "tam"].iloc[0] == "MYC"


def test_read_table_strips_utf8_bom_from_gene_id_header(tmp_path) -> None:
    from utils import read_table

    p = tmp_path / "matrix.csv"
    p.write_bytes("\ufeffgene_id,s1\ng1,1\n".encode("utf-8"))
    df = read_table(p)
    assert df.columns.tolist() == ["gene_id", "s1"]
