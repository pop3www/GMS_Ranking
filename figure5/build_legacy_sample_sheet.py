from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd

from fig6_ml_core import build_manifest, write_json

LEGACY_GROUPS = [
    "4DT_LCPT",
    "4D_LCPT",
    "4DT_LCPT",
    "4D_HCPT",
    "4DT_HCPT",
    "4_ctrl",
    "4D_HCPT",
    "4DT_HCPT",
    "24_ctrl",
    "24_Tam",
    "24_ctrl",
    "24_Tam",
    "24D",
    "24DT",
    "24D",
    "24DT",
    "4_Tam",
    "24D_LCPT",
    "24DT_LCPT",
    "24D_LCPT",
    "24DT_LCPT",
    "24D_HCPT",
    "24DT_HCPT",
    "24D_HCPT",
    "24DT_HCPT",
    "4_ctrl",
    "4_Tam",
    "4D",
    "4DT",
    "4D",
    "4DT",
    "4D_LCPT",
]

LEGACY_TIME = [
    4, 4, 4, 4, 4, 4, 4, 4,
    24, 24, 24, 24, 24, 24, 24, 24,
    4,
    24, 24, 24, 24, 24, 24, 24, 24,
    4, 4, 4, 4, 4, 4, 4,
]

LEGACY_DOX = [
    "Dox", "Dox", "Dox", "Dox", "Dox", "None_Dox", "Dox", "Dox",
    "None_Dox", "None_Dox", "None_Dox", "None_Dox", "Dox", "Dox", "Dox", "Dox",
    "None_Dox",
    "Dox", "Dox", "Dox", "Dox", "Dox", "Dox", "Dox", "Dox",
    "None_Dox", "None_Dox", "Dox", "Dox", "Dox", "Dox", "Dox",
]

LEGACY_MYC = [
    "MYC", "None_MYC", "MYC", "None_MYC", "MYC", "None_MYC", "None_MYC", "MYC",
    "None_MYC", "MYC", "None_MYC", "MYC", "None_MYC", "MYC", "None_MYC", "MYC",
    "MYC",
    "None_MYC", "MYC", "None_MYC", "MYC", "None_MYC", "MYC", "None_MYC", "MYC",
    "None_MYC", "MYC", "None_MYC", "MYC", "None_MYC", "MYC", "None_MYC",
]

LEGACY_LCPT = [
    "LCPT", "LCPT", "LCPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT",
    "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT",
    "None_CPT",
    "LCPT", "LCPT", "LCPT", "LCPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT",
    "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "LCPT",
]

LEGACY_HCPT = [
    "None_CPT", "None_CPT", "None_CPT", "HCPT", "HCPT", "None_CPT", "HCPT", "HCPT",
    "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT",
    "None_CPT",
    "None_CPT", "None_CPT", "None_CPT", "None_CPT", "HCPT", "HCPT", "HCPT", "HCPT",
    "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT", "None_CPT",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct a compatibility sample sheet from the archived 32-sample column order."
        )
    )
    parser.add_argument(
        "--matrix",
        default=None,
        help=(
            "Optional wide matrix used only to read the sample column names. "
            "If omitted, placeholder sample IDs are generated."
        ),
    )
    parser.add_argument(
        "--output",
        default="legacy_sample_sheet.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--id-col",
        default=None,
        help="Optional gene identifier column if the matrix header is available and the first column should not be inferred.",
    )
    parser.add_argument(
        "--allow-length-mismatch",
        action="store_true",
        help="Allow matrices whose number of sample columns does not match the archived 32-column contract.",
    )
    return parser.parse_args()



def sample_ids_from_matrix(matrix_path: str, id_col: str | None = None) -> List[str]:
    path = Path(matrix_path)
    if not path.exists():
        raise FileNotFoundError(f"Matrix not found: {path}")
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    header = pd.read_csv(path, sep=sep, nrows=0)
    if header.shape[1] < 2:
        raise ValueError("The matrix must contain at least one gene column and one sample column.")
    if id_col is None:
        id_col = header.columns[0]
    if id_col not in header.columns:
        raise ValueError(f"Gene identifier column '{id_col}' not present in matrix header.")
    return [str(c) for c in header.columns if c != id_col]



def main() -> int:
    args = parse_args()
    if args.matrix is not None:
        sample_ids = sample_ids_from_matrix(args.matrix, args.id_col)
    else:
        sample_ids = [f"sample_{i:02d}" for i in range(1, len(LEGACY_GROUPS) + 1)]

    if len(sample_ids) != len(LEGACY_GROUPS) and not args.allow_length_mismatch:
        raise ValueError(
            f"Legacy metadata expects {len(LEGACY_GROUPS)} sample columns, "
            f"but found {len(sample_ids)}. Use --allow-length-mismatch only if you know the column order has changed."
        )

    n = min(len(sample_ids), len(LEGACY_GROUPS))
    rows = []
    replicate_counts = {}
    for i in range(n):
        group_label = LEGACY_GROUPS[i]
        replicate_counts[group_label] = replicate_counts.get(group_label, 0) + 1
        cpt_level = "HCPT" if LEGACY_HCPT[i] == "HCPT" else LEGACY_LCPT[i]
        if cpt_level not in {"LCPT", "HCPT"}:
            cpt_level = "None_CPT"

        row = {
            "legacy_position": i + 1,
            "sample_id": sample_ids[i],
            "group_label": group_label,
            "time_h": LEGACY_TIME[i],
            "dox": LEGACY_DOX[i],
            "myc": LEGACY_MYC[i],
            "cpt_level": cpt_level,
            "replicate": replicate_counts[group_label],
            "context_myc_task": cpt_level,
            "context_dox_task": LEGACY_MYC[i],
        }
        rows.append(row)

    out_df = pd.DataFrame(rows)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)

    manifest = build_manifest(
        arguments=vars(args),
        input_files=({"matrix": args.matrix} if args.matrix is not None else {}),
        script_files=[Path(__file__), Path(__file__).with_name("fig6_ml_core.py")],
        extra={
            "n_rows_written": int(len(out_df)),
            "output_csv": str(out_path.resolve()),
            "legacy_contract_n_samples": int(len(LEGACY_GROUPS)),
        },
    )
    write_json(manifest, out_path.with_suffix(out_path.suffix + ".manifest.json"))

    print(f"[OK] Wrote {len(out_df)} rows to {out_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
