from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import warnings

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut, LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


REQUIRED_SAMPLE_SHEET_COLUMNS = [
    "sample_id",
    "group_label",
    "time_h",
    "dox",
    "myc",
    "cpt_level",
]

DEFAULT_DDR_GENES = {
    "TP53",
    "CDKN1A",
    "GADD45A",
    "BRCA1",
    "BRCA2",
    "RAD51",
    "ATM",
    "ATR",
    "CHEK1",
    "CHEK2",
    "MDM2",
    "PCNA",
}


class InputContractError(ValueError):
    """Raised when the matrix / sample-sheet contract is violated."""


@dataclass(frozen=True)
class TaskSpec:
    name: str
    time_h: int
    subset_mode: str  # "myc_task" or "dox_task"
    label_col: str
    positive_label: str
    negative_label: str
    context_col: str
    description: str
    include_in_main: bool = True


def default_task_specs(include_dox_controls: bool = True) -> List[TaskSpec]:
    tasks = [
        TaskSpec(
            name="MYC@4h",
            time_h=4,
            subset_mode="myc_task",
            label_col="myc",
            positive_label="MYC",
            negative_label="None_MYC",
            context_col="cpt_level",
            description=(
                "Sample-level MYC classification at 4 h within Dox-positive, "
                "matched CPT backgrounds."
            ),
            include_in_main=True,
        ),
        TaskSpec(
            name="MYC@24h",
            time_h=24,
            subset_mode="myc_task",
            label_col="myc",
            positive_label="MYC",
            negative_label="None_MYC",
            context_col="cpt_level",
            description=(
                "Sample-level MYC classification at 24 h within Dox-positive, "
                "matched CPT backgrounds."
            ),
            include_in_main=True,
        ),
    ]
    if include_dox_controls:
        tasks.extend(
            [
                TaskSpec(
                    name="DOX@4h",
                    time_h=4,
                    subset_mode="dox_task",
                    label_col="dox",
                    positive_label="Dox",
                    negative_label="None_Dox",
                    context_col="myc",
                    description=(
                        "Sample-level Dox classification at 4 h in the no-CPT arm, "
                        "used as a control task."
                    ),
                    include_in_main=False,
                ),
                TaskSpec(
                    name="DOX@24h",
                    time_h=24,
                    subset_mode="dox_task",
                    label_col="dox",
                    positive_label="Dox",
                    negative_label="None_Dox",
                    context_col="myc",
                    description=(
                        "Sample-level Dox classification at 24 h in the no-CPT arm, "
                        "used as a control task."
                    ),
                    include_in_main=False,
                ),
            ]
        )
    return tasks


class TopVarianceGeneSelector(BaseEstimator, TransformerMixin):
    """Select the highest-variance genes on the training fold only."""

    def __init__(self, max_genes: int = 500, min_genes: int = 2):
        self.max_genes = int(max_genes)
        self.min_genes = int(min_genes)

    def fit(self, X, y=None):
        X_arr = _to_numpy(X)
        if X_arr.ndim != 2:
            raise ValueError("TopVarianceGeneSelector expects a 2D matrix.")
        n_features = X_arr.shape[1]
        if n_features == 0:
            raise ValueError("The input matrix has zero features.")
        if isinstance(X, pd.DataFrame):
            self.feature_names_in_ = np.asarray(X.columns.astype(str), dtype=object)
        else:
            self.feature_names_in_ = np.asarray([f"feature_{i}" for i in range(n_features)], dtype=object)
        variances = np.nanvar(X_arr, axis=0)
        variances = np.where(np.isnan(variances), -np.inf, variances)
        keep_n = min(max(self.min_genes, 1), n_features) if self.max_genes <= 0 else min(max(self.min_genes, min(self.max_genes, n_features)), n_features)
        order = np.lexsort((self.feature_names_in_.astype(str), -variances))
        self.support_indices_ = np.sort(order[:keep_n])
        self.feature_names_out_ = self.feature_names_in_[self.support_indices_]
        return self

    def transform(self, X):
        X_arr = _to_numpy(X)
        return X_arr[:, self.support_indices_]

    def get_feature_names_out(self, input_features=None):
        return np.asarray(self.feature_names_out_, dtype=object)


class AdaptiveTruncatedSVD(BaseEstimator, TransformerMixin):
    """Choose a valid SVD dimensionality from the training fold automatically."""

    def __init__(self, max_components: int = 6, random_state: int = 42):
        self.max_components = int(max_components)
        self.random_state = int(random_state)

    def fit(self, X, y=None):
        X_arr = _to_numpy(X)
        n_samples, n_features = X_arr.shape
        max_valid = min(self.max_components, max(n_samples - 1, 1), max(n_features - 1, 1))
        if max_valid < 2:
            self.model_ = None
            self.n_components_ = X_arr.shape[1]
            return self
        self.model_ = TruncatedSVD(n_components=max_valid, random_state=self.random_state)
        self.model_.fit(X_arr)
        self.n_components_ = max_valid
        return self

    def transform(self, X):
        X_arr = _to_numpy(X)
        if getattr(self, "model_", None) is None:
            return X_arr
        return self.model_.transform(X_arr)


def _to_numpy(X) -> np.ndarray:
    if isinstance(X, pd.DataFrame):
        return X.to_numpy(dtype=float, copy=False)
    if isinstance(X, pd.Series):
        return X.to_frame().to_numpy(dtype=float, copy=False)
    return np.asarray(X, dtype=float)


def _read_delimited_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    sep = "," if suffix == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def read_wide_matrix(
    matrix_path: str | Path,
    id_col: Optional[str] = None,
    apply_log1p: bool = False,
    duplicate_policy: str = "error",
    fill_missing_value: Optional[float] = None,
) -> pd.DataFrame:
    """Read a gene-by-sample matrix and return a DataFrame indexed by gene."""

    path = Path(matrix_path)
    if not path.exists():
        raise FileNotFoundError(f"Matrix file not found: {path}")

    df = _read_delimited_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 2:
        raise InputContractError(
            "The expression matrix must have one gene identifier column and at least one sample column."
        )

    if id_col is None:
        id_col = df.columns[0]
    if id_col not in df.columns:
        raise InputContractError(f"Gene identifier column '{id_col}' not found in {path}.")

    df = df.dropna(subset=[id_col]).copy()
    df[id_col] = df[id_col].astype(str).str.strip()

    sample_cols = [c for c in df.columns if c != id_col]
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if df[sample_cols].isna().any().any():
        missing_total = int(df[sample_cols].isna().sum().sum())
        if fill_missing_value is None:
            raise InputContractError(
                f"The matrix contains {missing_total} missing numeric values after parsing. "
                "Provide a clean matrix or set a fill value explicitly."
            )
        warnings.warn(
            f"Filling {missing_total} missing matrix values with {fill_missing_value}.",
            RuntimeWarning,
        )
        df[sample_cols] = df[sample_cols].fillna(float(fill_missing_value))

    duplicated = df[id_col].duplicated(keep=False)
    if duplicated.any():
        dup_count = int(duplicated.sum())
        if duplicate_policy == "error":
            raise InputContractError(
                f"The matrix contains {dup_count} duplicated gene identifiers. "
                "Deduplicate upstream or use duplicate_policy='mean'."
            )
        if duplicate_policy != "mean":
            raise InputContractError(f"Unsupported duplicate_policy: {duplicate_policy}")
        df = df.groupby(id_col, as_index=False)[sample_cols].mean()

    if apply_log1p:
        if (df[sample_cols] < 0).any().any():
            raise InputContractError("Cannot apply log1p because the matrix contains negative values.")
        df[sample_cols] = np.log1p(df[sample_cols])

    expr = df.set_index(id_col)
    if expr.index.duplicated().any():
        raise InputContractError("Gene identifiers remain duplicated after processing.")
    return expr


_NORMALIZE_DOX = {
    "dox": "Dox",
    "yes": "Dox",
    "1": "Dox",
    "true": "Dox",
    "none_dox": "None_Dox",
    "no_dox": "None_Dox",
    "none": "None_Dox",
    "no": "None_Dox",
    "0": "None_Dox",
    "false": "None_Dox",
}

_NORMALIZE_MYC = {
    "myc": "MYC",
    "on": "MYC",
    "yes": "MYC",
    "1": "MYC",
    "true": "MYC",
    "none_myc": "None_MYC",
    "off": "None_MYC",
    "no": "None_MYC",
    "0": "None_MYC",
    "false": "None_MYC",
}

_NORMALIZE_CPT = {
    "lcpt": "LCPT",
    "low": "LCPT",
    "low_cpt": "LCPT",
    "hcpt": "HCPT",
    "high": "HCPT",
    "high_cpt": "HCPT",
    "none": "None_CPT",
    "none_cpt": "None_CPT",
    "no_cpt": "None_CPT",
    "0": "None_CPT",
    "": "None_CPT",
}


def _normalize_text_label(value: object, mapping: Dict[str, str], field_name: str) -> str:
    text = str(value).strip()
    key = text.lower().replace("-", "_").replace(" ", "_")
    if key in mapping:
        return mapping[key]
    raise InputContractError(
        f"Unrecognized value '{value}' in sample sheet field '{field_name}'."
    )


def _normalize_time_h(value: object) -> int:
    text = str(value).strip().lower().replace("hours", "").replace("hour", "").replace("h", "")
    if text in {"4", "24"}:
        return int(text)
    raise InputContractError(
        f"Unrecognized time_h value '{value}'. Expected 4/24 or equivalent text labels."
    )


def read_sample_sheet(sample_sheet_path: str | Path) -> pd.DataFrame:
    path = Path(sample_sheet_path)
    if not path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {path}")

    meta = _read_delimited_table(path)
    meta.columns = [str(c).strip() for c in meta.columns]
    missing_cols = [c for c in REQUIRED_SAMPLE_SHEET_COLUMNS if c not in meta.columns]
    if missing_cols:
        raise InputContractError(
            "The sample sheet is missing required columns: " + ", ".join(missing_cols)
        )

    meta = meta.copy()
    if meta["sample_id"].duplicated().any():
        dup = meta.loc[meta["sample_id"].duplicated(), "sample_id"].tolist()
        raise InputContractError(
            "The sample sheet contains duplicated sample_id entries: " + ", ".join(map(str, dup[:10]))
        )

    meta["sample_id"] = meta["sample_id"].astype(str).str.strip()
    meta["group_label"] = meta["group_label"].astype(str).str.strip()
    meta["time_h"] = meta["time_h"].map(_normalize_time_h)
    meta["dox"] = meta["dox"].map(lambda x: _normalize_text_label(x, _NORMALIZE_DOX, "dox"))
    meta["myc"] = meta["myc"].map(lambda x: _normalize_text_label(x, _NORMALIZE_MYC, "myc"))
    meta["cpt_level"] = meta["cpt_level"].map(lambda x: _normalize_text_label(x, _NORMALIZE_CPT, "cpt_level"))

    if "replicate" not in meta.columns:
        meta["replicate"] = meta.groupby("group_label").cumcount() + 1
    else:
        meta["replicate"] = pd.to_numeric(meta["replicate"], errors="coerce")
        if meta["replicate"].isna().any():
            raise InputContractError("The sample sheet column 'replicate' must be numeric if provided.")
        meta["replicate"] = meta["replicate"].astype(int)

    if "background_context" not in meta.columns:
        meta["background_context"] = np.where(
            meta["cpt_level"] != "None_CPT",
            meta["cpt_level"],
            np.where(meta["myc"] == "MYC", "MYC_background", "baseline_background"),
        )

    return meta


def validate_and_align_inputs(expr: pd.DataFrame, meta: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sample_cols = [str(c).strip() for c in expr.columns]
    meta_ids = meta["sample_id"].astype(str).tolist()
    missing_in_meta = [s for s in sample_cols if s not in set(meta_ids)]
    missing_in_matrix = [s for s in meta_ids if s not in set(sample_cols)]
    if missing_in_meta or missing_in_matrix:
        lines = []
        if missing_in_meta:
            lines.append("Missing metadata for matrix columns: " + ", ".join(missing_in_meta[:10]))
        if missing_in_matrix:
            lines.append("Sample-sheet rows absent from matrix: " + ", ".join(missing_in_matrix[:10]))
        raise InputContractError("Sample identifiers do not match between matrix and sample sheet. " + " | ".join(lines))

    meta_aligned = meta.set_index("sample_id").loc[sample_cols].reset_index()
    return expr.loc[:, sample_cols], meta_aligned


def expression_to_sample_by_gene(expr: pd.DataFrame) -> pd.DataFrame:
    X = expr.T.copy()
    X.index.name = "sample_id"
    X.columns = X.columns.astype(str)
    return X


def build_task_table(
    X_samples: pd.DataFrame,
    meta: pd.DataFrame,
    task: TaskSpec,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    meta = meta.copy()
    if task.subset_mode == "myc_task":
        mask = (meta["time_h"] == task.time_h) & (meta["dox"] == "Dox")
    elif task.subset_mode == "dox_task":
        mask = (meta["time_h"] == task.time_h) & (meta["cpt_level"] == "None_CPT")
    else:
        raise InputContractError(f"Unknown subset_mode: {task.subset_mode}")

    task_meta = meta.loc[mask].copy()
    if task_meta.empty:
        raise InputContractError(f"Task {task.name} selected zero samples.")

    task_meta["label_text"] = np.where(
        task_meta[task.label_col] == task.positive_label,
        task.positive_label,
        task.negative_label,
    )
    label_counts = task_meta["label_text"].value_counts()
    if len(label_counts) != 2:
        raise InputContractError(
            f"Task {task.name} does not have two labels after filtering. Counts: {label_counts.to_dict()}"
        )

    task_meta["label"] = (task_meta[task.label_col] == task.positive_label).astype(int)
    task_meta["task_name"] = task.name
    task_meta["context_value"] = task_meta[task.context_col].astype(str)
    task_meta = task_meta.sort_values(
        ["context_value", "label", "group_label", "replicate", "sample_id"],
        kind="mergesort",
    ).reset_index(drop=True)

    X_task = X_samples.loc[task_meta["sample_id"]].copy()
    if X_task.shape[0] != task_meta.shape[0]:
        raise InputContractError(f"Failed to align expression rows for task {task.name}.")
    return X_task, task_meta.reset_index(drop=True)


def build_cv_schemes(task_meta: pd.DataFrame) -> Dict[str, List[Tuple[np.ndarray, np.ndarray]]]:
    y = task_meta["label"].to_numpy()
    schemes: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = {}
    loo = list(LeaveOneOut().split(np.zeros(len(task_meta)), y))
    schemes["leave_one_sample_out"] = loo

    groups = task_meta["context_value"].to_numpy()
    group_ok = True
    for group_name, sub in task_meta.groupby("context_value"):
        if sub["label"].nunique() < 2:
            group_ok = False
            break
    if group_ok and len(np.unique(groups)) >= 2:
        logo = list(LeaveOneGroupOut().split(np.zeros(len(task_meta)), y, groups=groups))
        if logo:
            schemes["leave_one_context_out"] = logo
    return schemes


def build_estimator(
    model_name: str,
    random_state: int = 42,
    max_genes: int = 500,
    max_svd_components: int = 6,
):
    model_name = str(model_name).strip().lower()
    selector = TopVarianceGeneSelector(max_genes=max_genes)
    scaler = StandardScaler(with_mean=True, with_std=True)

    if model_name == "ridge_logistic":
        clf = LogisticRegression(
            C=1.0,
            class_weight="balanced",
            solver="liblinear",
            max_iter=10000,
            random_state=random_state,
        )
        return Pipeline(
            steps=[
                ("select", selector),
                ("scale", scaler),
                ("clf", clf),
            ]
        )

    if model_name == "gbt_svd":
        clf = GradientBoostingClassifier(
            learning_rate=0.05,
            n_estimators=150,
            max_depth=1,
            subsample=0.9,
            random_state=random_state,
        )
        return Pipeline(
            steps=[
                ("select", selector),
                ("scale", scaler),
                ("svd", AdaptiveTruncatedSVD(max_components=max_svd_components, random_state=random_state)),
                ("clf", clf),
            ]
        )

    raise InputContractError(
        f"Unsupported model_name '{model_name}'. Use 'ridge_logistic' or 'gbt_svd'."
    )


@dataclass
class CVResult:
    metrics: Dict[str, float]
    predictions: pd.DataFrame
    confusion: np.ndarray


def cross_validate_binary_model(
    estimator,
    X: pd.DataFrame,
    task_meta: pd.DataFrame,
    splits: Sequence[Tuple[np.ndarray, np.ndarray]],
    task_name: str,
    model_name: str,
    cv_scheme: str,
) -> CVResult:
    y = task_meta["label"].to_numpy(dtype=int)
    pred_label = np.full(len(y), -1, dtype=int)
    pred_score = np.full(len(y), np.nan, dtype=float)
    fold_id = np.full(len(y), -1, dtype=int)

    for fold_idx, (train_idx, test_idx) in enumerate(splits, start=1):
        if len(np.unique(y[train_idx])) < 2:
            raise InputContractError(
                f"Fold {fold_idx} of {task_name}/{model_name}/{cv_scheme} has only one training class."
            )
        est = clone(estimator)
        est.fit(X.iloc[train_idx], y[train_idx])

        if hasattr(est, "predict_proba"):
            score = est.predict_proba(X.iloc[test_idx])[:, 1]
            pred = (score >= 0.5).astype(int)
        elif hasattr(est, "decision_function"):
            score = est.decision_function(X.iloc[test_idx])
            pred = (score >= 0.0).astype(int)
        else:
            pred = est.predict(X.iloc[test_idx]).astype(int)
            score = pred.astype(float)

        pred_label[test_idx] = pred
        pred_score[test_idx] = score
        fold_id[test_idx] = fold_idx

    if (pred_label < 0).any():
        raise RuntimeError(
            f"Internal error: missing held-out predictions for {task_name}/{model_name}/{cv_scheme}."
        )

    bal_acc = balanced_accuracy_score(y, pred_label)
    try:
        auc = roc_auc_score(y, pred_score)
    except ValueError:
        auc = np.nan
    cm = confusion_matrix(y, pred_label, labels=[0, 1])

    pos_label_text = str(task_meta.loc[task_meta["label"] == 1, "label_text"].iloc[0])
    neg_label_text = str(task_meta.loc[task_meta["label"] == 0, "label_text"].iloc[0])

    pred_df = task_meta.copy()
    pred_df["pred_label"] = pred_label
    pred_df["pred_label_text"] = np.where(pred_label == 1, pos_label_text, neg_label_text)
    pred_df["pred_score"] = pred_score
    pred_df["fold_id"] = fold_id
    pred_df["correct"] = pred_df["pred_label"].eq(pred_df["label"])
    pred_df["task_name"] = task_name
    pred_df["model_name"] = model_name
    pred_df["cv_scheme"] = cv_scheme

    metrics = {
        "n_samples": float(len(y)),
        "n_positive": float(int(y.sum())),
        "n_negative": float(int((1 - y).sum())),
        "balanced_accuracy": float(bal_acc),
        "roc_auc": float(auc) if not np.isnan(auc) else np.nan,
    }
    return CVResult(metrics=metrics, predictions=pred_df, confusion=cm)


def fit_final_model(
    estimator,
    X: pd.DataFrame,
    task_meta: pd.DataFrame,
    task_name: str,
    model_name: str,
):
    y = task_meta["label"].to_numpy(dtype=int)
    est = clone(estimator)
    est.fit(X, y)
    return est


def compute_display_embedding(
    fitted_estimator,
    X: pd.DataFrame,
) -> pd.DataFrame:
    steps = dict(fitted_estimator.named_steps)
    if "select" not in steps or "scale" not in steps:
        raise InputContractError("Expected 'select' and 'scale' steps in the fitted pipeline.")

    X_trans = steps["scale"].transform(steps["select"].transform(X))

    if "svd" in steps:
        X_embed = steps["svd"].transform(X_trans)
    else:
        n_components = min(2, X_trans.shape[0] - 1, X_trans.shape[1])
        if n_components < 2:
            pad = np.zeros((X_trans.shape[0], 2), dtype=float)
            if X_trans.shape[1] >= 1:
                pad[:, 0] = X_trans[:, 0]
            return pd.DataFrame({"dim1": pad[:, 0], "dim2": pad[:, 1]})
        svd = TruncatedSVD(n_components=2, random_state=42)
        X_embed = svd.fit_transform(X_trans)

    if X_embed.shape[1] == 1:
        X_embed = np.column_stack([X_embed[:, 0], np.zeros(X_embed.shape[0], dtype=float)])
    elif X_embed.shape[1] < 2:
        X_embed = np.zeros((X_embed.shape[0], 2), dtype=float)

    return pd.DataFrame({"dim1": X_embed[:, 0], "dim2": X_embed[:, 1]})


def extract_linear_coefficients(fitted_estimator) -> Optional[pd.DataFrame]:
    steps = dict(fitted_estimator.named_steps)
    clf = steps.get("clf")
    selector = steps.get("select")
    if clf is None or selector is None or not hasattr(clf, "coef_"):
        return None
    coef = np.asarray(clf.coef_).ravel()
    genes = selector.get_feature_names_out()
    coef_df = pd.DataFrame(
        {
            "gene": genes,
            "coefficient": coef,
            "abs_coefficient": np.abs(coef),
        }
    ).sort_values(["abs_coefficient", "coefficient"], ascending=[False, False])
    return coef_df.reset_index(drop=True)


def permutation_test_balanced_accuracy(
    estimator,
    X: pd.DataFrame,
    task_meta: pd.DataFrame,
    splits: Sequence[Tuple[np.ndarray, np.ndarray]],
    n_permutations: int = 0,
    random_state: int = 42,
    permute_within_context: bool = True,
) -> Dict[str, float]:
    if n_permutations <= 0:
        return {
            "null_mean_balanced_accuracy": np.nan,
            "null_sd_balanced_accuracy": np.nan,
            "permutation_pvalue": np.nan,
            "n_valid_permutations": 0,
        }

    rng = np.random.default_rng(random_state)
    y_true = task_meta["label"].to_numpy(dtype=int)
    groups = task_meta["context_value"].to_numpy()

    observed = cross_validate_binary_model(
        estimator=estimator,
        X=X,
        task_meta=task_meta,
        splits=splits,
        task_name=str(task_meta["task_name"].iloc[0]),
        model_name="_perm_observed",
        cv_scheme="_perm_observed",
    ).metrics["balanced_accuracy"]

    null_scores: List[float] = []
    for _ in range(int(n_permutations)):
        y_perm = y_true.copy()
        if permute_within_context:
            for group in np.unique(groups):
                idx = np.flatnonzero(groups == group)
                shuffled = y_true[idx].copy()
                rng.shuffle(shuffled)
                y_perm[idx] = shuffled
        else:
            rng.shuffle(y_perm)

        perm_meta = task_meta.copy()
        perm_meta["label"] = y_perm
        try:
            perm_result = cross_validate_binary_model(
                estimator=estimator,
                X=X,
                task_meta=perm_meta,
                splits=splits,
                task_name=str(task_meta["task_name"].iloc[0]),
                model_name="_perm",
                cv_scheme="_perm",
            )
        except InputContractError:
            continue
        null_scores.append(perm_result.metrics["balanced_accuracy"])

    if not null_scores:
        return {
            "null_mean_balanced_accuracy": np.nan,
            "null_sd_balanced_accuracy": np.nan,
            "permutation_pvalue": np.nan,
            "n_valid_permutations": 0,
        }

    null_arr = np.asarray(null_scores, dtype=float)
    pval = (1.0 + float(np.sum(null_arr >= observed))) / (len(null_arr) + 1.0)
    return {
        "null_mean_balanced_accuracy": float(np.nanmean(null_arr)),
        "null_sd_balanced_accuracy": float(np.nanstd(null_arr, ddof=1)) if len(null_arr) > 1 else 0.0,
        "permutation_pvalue": float(pval),
        "n_valid_permutations": int(len(null_arr)),
    }


def benjamini_hochberg(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return np.asarray([], dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty_like(adj)
    out[order] = adj
    return out


def compute_per_gene_auc(
    X_task: pd.DataFrame,
    task_meta: pd.DataFrame,
    positive_label: str,
    min_per_class: int = 2,
) -> pd.DataFrame:
    y = (task_meta["label_text"] == positive_label).to_numpy(dtype=bool)
    n_pos = int(y.sum())
    n_neg = int((~y).sum())
    if min(n_pos, n_neg) < int(min_per_class):
        raise InputContractError(
            f"Not enough samples per class for per-gene AUROC: pos={n_pos}, neg={n_neg}."
        )

    rows = []
    for gene in X_task.columns:
        x = X_task[gene].to_numpy(dtype=float)
        try:
            auc = roc_auc_score(y.astype(int), x)
        except ValueError:
            continue
        pos = x[y]
        neg = x[~y]
        # Use a Mann–Whitney-compatible AUC-only screen; p-value is approximated later upstream if needed.
        rows.append((gene, auc, float(np.mean(pos)), float(np.mean(neg)), float(np.mean(pos) - np.mean(neg))))

    res = pd.DataFrame(rows, columns=["gene", "auc", "mean_positive", "mean_negative", "mean_delta"])
    if res.empty:
        return res
    res = res.sort_values("auc", ascending=False).reset_index(drop=True)
    return res


def compute_baseline_expression_regions(
    X_samples: pd.DataFrame,
    meta: pd.DataFrame,
    head_pct: float = 10.0,
    tail_pct: float = 10.0,
) -> pd.DataFrame:
    base_mask = (meta["time_h"] == 4) & (meta["dox"] == "None_Dox")
    base_samples = meta.loc[base_mask, "sample_id"].tolist()
    if not base_samples:
        raise InputContractError("Could not find 4 h / None_Dox samples for baseline head-tail ranking.")

    baseline_mean = X_samples.loc[base_samples].mean(axis=0).sort_values(ascending=False)
    n = len(baseline_mean)
    head_n = max(1, int(np.ceil(n * head_pct / 100.0)))
    tail_n = max(1, int(np.ceil(n * tail_pct / 100.0)))

    region = pd.Series("mid", index=baseline_mean.index, dtype=object)
    region.iloc[:head_n] = "head"
    region.iloc[-tail_n:] = "tail"
    return pd.DataFrame(
        {
            "gene": baseline_mean.index.astype(str),
            "baseline_mean": baseline_mean.to_numpy(dtype=float),
            "region": region.to_numpy(dtype=object),
        }
    )


def save_dataframe(df: pd.DataFrame, out_path: str | Path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False, float_format="%.10g")


def save_array_as_csv(arr: np.ndarray, row_labels: Sequence[str], col_labels: Sequence[str], out_path: str | Path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(arr, index=row_labels, columns=col_labels)
    df.to_csv(out_path, float_format="%.10g")


def write_json(data: Dict[str, object], out_path: str | Path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def sha256_file(path: str | Path) -> str:
    path = Path(path)
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def runtime_versions() -> Dict[str, str]:
    import matplotlib
    import scipy
    import sklearn

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
    }


def build_manifest(
    arguments: Dict[str, object],
    input_files: Dict[str, str | Path],
    script_files: Optional[Sequence[str | Path]] = None,
    task_specs: Optional[Sequence[TaskSpec]] = None,
    extra: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    manifest: Dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "python_executable": sys.executable,
        "working_directory": str(Path.cwd().resolve()),
        "software_versions": runtime_versions(),
        "arguments": {k: (str(v) if isinstance(v, Path) else v) for k, v in arguments.items()},
        "input_files": {},
        "script_files": {},
    }

    for label, file_path in input_files.items():
        p = Path(file_path).resolve()
        manifest["input_files"][label] = {
            "path": str(p),
            "sha256": sha256_file(p),
            "size_bytes": int(p.stat().st_size),
        }

    for file_path in script_files or []:
        p = Path(file_path).resolve()
        manifest["script_files"][p.name] = {
            "path": str(p),
            "sha256": sha256_file(p),
            "size_bytes": int(p.stat().st_size),
        }

    if task_specs is not None:
        manifest["task_specs"] = [asdict(task) for task in task_specs]
    if extra:
        manifest["extra"] = extra
    return manifest
