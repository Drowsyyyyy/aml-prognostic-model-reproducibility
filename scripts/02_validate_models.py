"""Repeated nested cross-validation for TCGA-LAML survival models."""

from __future__ import annotations

import json
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
OUTPUT = ROOT / "outputs"
TABLES = OUTPUT / "tables"
FIGURES = OUTPUT / "figures"
RANDOM_SEED = 20260724
OUTER_SPLITS = 5
OUTER_REPEATS = 5
INNER_SPLITS = 3
TOP_GENES = 250
L1_RATIO = 0.9
COXNET_ALPHAS = np.geomspace(0.5, 0.01, 18)
COXPH_ALPHAS = np.array([0.01, 0.1, 1.0, 10.0])
N_PERMUTATIONS = 20
CLINICAL_PENALTY_WEIGHT = 0.05

NUMERIC_CLINICAL = ["age", "log1p_wbc", "bm_blast_pct"]
CATEGORICAL_CLINICAL = ["sex", "risk_cyto"]


def make_survival(clinical: pd.DataFrame) -> np.ndarray:
    return Surv.from_arrays(
        event=clinical["event"].astype(bool).to_numpy(),
        time=clinical["os_months"].astype(float).to_numpy(),
    )


def cindex(y: np.ndarray, risk: np.ndarray) -> float:
    return float(
        concordance_index_censored(y["event"], y["time"], risk)[0]
    )


def clinical_transformer() -> ColumnTransformer:
    numeric = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
        ]
    )
    categorical = Pipeline(
        [
            ("impute", SimpleImputer(strategy="most_frequent")),
            (
                "onehot",
                OneHotEncoder(
                    handle_unknown="ignore", drop="first", sparse_output=False
                ),
            ),
        ]
    )
    return ColumnTransformer(
        [
            ("numeric", numeric, NUMERIC_CLINICAL),
            ("categorical", categorical, CATEGORICAL_CLINICAL),
        ],
        remainder="drop",
        verbose_feature_names_out=False,
    )


def fit_expression_transform(
    train: pd.DataFrame, other: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    variances = train.var(axis=0, skipna=True).fillna(0.0)
    genes = variances.nlargest(min(TOP_GENES, len(variances))).index.tolist()
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    train_array = scaler.fit_transform(imputer.fit_transform(train[genes]))
    other_array = scaler.transform(imputer.transform(other[genes]))
    return train_array, other_array, genes


def fit_clinical_transform(
    train: pd.DataFrame, other: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    transformer = clinical_transformer()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        train_array = transformer.fit_transform(train)
        other_array = transformer.transform(other)
    names = transformer.get_feature_names_out().tolist()
    return train_array, other_array, names


def tune_clinical(
    clinical: pd.DataFrame,
    y: np.ndarray,
    event: np.ndarray,
    seed: int,
) -> float:
    scores = defaultdict(list)
    splitter = StratifiedKFold(
        n_splits=INNER_SPLITS, shuffle=True, random_state=seed
    )
    for train_idx, valid_idx in splitter.split(clinical, event):
        x_train, x_valid, _ = fit_clinical_transform(
            clinical.iloc[train_idx], clinical.iloc[valid_idx]
        )
        for alpha in COXPH_ALPHAS:
            model = CoxPHSurvivalAnalysis(alpha=float(alpha))
            model.fit(x_train, y[train_idx])
            scores[float(alpha)].append(
                cindex(y[valid_idx], model.predict(x_valid))
            )
    ranked = sorted(
        ((np.mean(values), alpha) for alpha, values in scores.items()),
        key=lambda item: (-item[0], -item[1]),
    )
    return float(ranked[0][1])


def tune_coxnet(
    clinical: pd.DataFrame,
    expression: pd.DataFrame,
    y: np.ndarray,
    event: np.ndarray,
    seed: int,
    combined: bool,
) -> float:
    scores = defaultdict(list)
    splitter = StratifiedKFold(
        n_splits=INNER_SPLITS, shuffle=True, random_state=seed
    )
    for train_idx, valid_idx in splitter.split(expression, event):
        xg_train, xg_valid, _ = fit_expression_transform(
            expression.iloc[train_idx], expression.iloc[valid_idx]
        )
        if combined:
            xc_train, xc_valid, _ = fit_clinical_transform(
                clinical.iloc[train_idx], clinical.iloc[valid_idx]
            )
            x_train = np.column_stack([xc_train, xg_train])
            x_valid = np.column_stack([xc_valid, xg_valid])
            penalty_factor = np.concatenate(
                [
                    np.full(xc_train.shape[1], CLINICAL_PENALTY_WEIGHT),
                    np.ones(xg_train.shape[1]),
                ]
            )
        else:
            x_train, x_valid = xg_train, xg_valid
            penalty_factor = None

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = CoxnetSurvivalAnalysis(
                alphas=COXNET_ALPHAS,
                l1_ratio=L1_RATIO,
                penalty_factor=penalty_factor,
                max_iter=100000,
                tol=1e-7,
            ).fit(x_train, y[train_idx])
        for alpha in model.alphas_:
            scores[float(alpha)].append(
                cindex(y[valid_idx], model.predict(x_valid, alpha=alpha))
            )

    complete = [
        (np.mean(values), alpha)
        for alpha, values in scores.items()
        if len(values) == INNER_SPLITS
    ]
    if not complete:
        raise RuntimeError("No Coxnet alpha was fitted in every inner fold.")
    ranked = sorted(complete, key=lambda item: (-item[0], -item[1]))
    return float(ranked[0][1])


def fit_predict_outer(
    clinical_train: pd.DataFrame,
    clinical_test: pd.DataFrame,
    expression_train: pd.DataFrame,
    expression_test: pd.DataFrame,
    y_train: np.ndarray,
    model_name: str,
    alpha: float,
) -> tuple[np.ndarray, dict[str, float], list[str]]:
    if model_name == "Clinical":
        x_train, x_test, names = fit_clinical_transform(
            clinical_train, clinical_test
        )
        model = CoxPHSurvivalAnalysis(alpha=alpha).fit(x_train, y_train)
        coefficients = dict(zip(names, model.coef_, strict=True))
        return model.predict(x_test), coefficients, names

    xg_train, xg_test, genes = fit_expression_transform(
        expression_train, expression_test
    )
    if model_name == "Combined":
        xc_train, xc_test, clinical_names = fit_clinical_transform(
            clinical_train, clinical_test
        )
        x_train = np.column_stack([xc_train, xg_train])
        x_test = np.column_stack([xc_test, xg_test])
        names = clinical_names + [f"gene:{gene}" for gene in genes]
        penalty_factor = np.concatenate(
            [
                np.full(xc_train.shape[1], CLINICAL_PENALTY_WEIGHT),
                np.ones(xg_train.shape[1]),
            ]
        )
    else:
        x_train, x_test = xg_train, xg_test
        names = [f"gene:{gene}" for gene in genes]
        penalty_factor = None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = CoxnetSurvivalAnalysis(
            alphas=COXNET_ALPHAS,
            l1_ratio=L1_RATIO,
            penalty_factor=penalty_factor,
            max_iter=100000,
            tol=1e-7,
        ).fit(x_train, y_train)
    fitted_alpha = float(model.alphas_[np.argmin(np.abs(model.alphas_ - alpha))])
    position = int(np.argmin(np.abs(model.alphas_ - fitted_alpha)))
    coefficients = dict(zip(names, model.coef_[:, position], strict=True))
    return (
        model.predict(x_test, alpha=fitted_alpha),
        coefficients,
        names,
    )


def zscore(values: np.ndarray) -> np.ndarray:
    std = float(np.std(values))
    if std == 0.0:
        return np.zeros_like(values, dtype=float)
    return (values - float(np.mean(values))) / std


def summarize_performance(folds: pd.DataFrame) -> pd.DataFrame:
    return (
        folds.groupby("model")["c_index"]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reset_index()
        .sort_values("mean", ascending=False)
    )


def permutation_control(
    expression: pd.DataFrame,
    y: np.ndarray,
    event: np.ndarray,
    fixed_alpha: float,
) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED + 9000)
    rows = []
    for permutation in range(1, N_PERMUTATIONS + 1):
        order = rng.permutation(len(y))
        y_permuted = y[order]
        splitter = StratifiedKFold(
            n_splits=OUTER_SPLITS,
            shuffle=True,
            random_state=RANDOM_SEED + permutation,
        )
        fold_scores = []
        for fold, (train_idx, test_idx) in enumerate(
            splitter.split(expression, y_permuted["event"]), start=1
        ):
            x_train, x_test, _ = fit_expression_transform(
                expression.iloc[train_idx], expression.iloc[test_idx]
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = CoxnetSurvivalAnalysis(
                    alphas=COXNET_ALPHAS,
                    l1_ratio=L1_RATIO,
                    max_iter=100000,
                    tol=1e-7,
                ).fit(x_train, y_permuted[train_idx])
            fitted_alpha = float(
                model.alphas_[np.argmin(np.abs(model.alphas_ - fixed_alpha))]
            )
            score = cindex(
                y_permuted[test_idx],
                model.predict(x_test, alpha=fitted_alpha),
            )
            fold_scores.append(score)
            rows.append(
                {
                    "permutation": permutation,
                    "fold": fold,
                    "c_index": score,
                }
            )
        rows.append(
            {
                "permutation": permutation,
                "fold": "mean",
                "c_index": float(np.mean(fold_scores)),
            }
        )
    return pd.DataFrame(rows)


def save_figures(
    fold_results: pd.DataFrame,
    gene_stability: pd.DataFrame,
    permutation: pd.DataFrame,
    clinical: pd.DataFrame,
    oof: pd.DataFrame,
) -> dict[str, float]:
    sns.set_theme(style="whitegrid", context="talk")

    plt.figure(figsize=(9, 6))
    order = (
        fold_results.groupby("model")["c_index"]
        .mean()
        .sort_values(ascending=False)
        .index
    )
    sns.boxplot(
        data=fold_results,
        x="model",
        y="c_index",
        order=order,
        color="#7aa6c2",
        width=0.55,
        showfliers=False,
    )
    sns.stripplot(
        data=fold_results,
        x="model",
        y="c_index",
        order=order,
        color="#183b56",
        alpha=0.55,
        size=4,
    )
    plt.axhline(0.5, color="#a33", linestyle="--", linewidth=1.5, label="Chance")
    plt.ylim(0.35, 0.85)
    plt.ylabel("Harrell C-index")
    plt.xlabel("")
    plt.title("Repeated outer-fold performance")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(FIGURES / "model_cindex.png", dpi=180)
    plt.close()

    stable = gene_stability.head(20).sort_values("selection_rate")
    plt.figure(figsize=(9, 7))
    colors = ["#c45a47" if value < 0 else "#3478a8" for value in stable["median_coef"]]
    plt.barh(stable["gene"], stable["selection_rate"], color=colors)
    plt.xlim(0, 1)
    plt.xlabel("Fraction of outer folds with non-zero coefficient")
    plt.title("Expression-model feature stability")
    plt.tight_layout()
    plt.savefig(FIGURES / "gene_stability.png", dpi=180)
    plt.close()

    permutation_means = permutation[permutation["fold"].eq("mean")].copy()
    actual_mean = float(
        fold_results.loc[
            fold_results["model"].eq("Expression"), "c_index"
        ].mean()
    )
    plt.figure(figsize=(9, 6))
    sns.histplot(
        permutation_means["c_index"],
        bins=10,
        color="#999999",
        edgecolor="white",
    )
    plt.axvline(actual_mean, color="#b33", linewidth=2.5, label="Observed mean")
    plt.axvline(0.5, color="#183b56", linestyle="--", label="Chance")
    plt.xlabel("Mean five-fold C-index")
    plt.title("Outcome-label permutation negative control")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(FIGURES / "permutation_control.png", dpi=180)
    plt.close()

    average_oof = (
        oof.groupby(["patient_id", "model"], as_index=False)["risk_z"].mean()
        .pivot(index="patient_id", columns="model", values="risk_z")
    )
    risk = average_oof["Combined"].reindex(clinical.index)
    cutoff = float(risk.median())
    high = risk >= cutoff
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(9, 6))
    for mask, label, color in [
        (~high, "Lower predicted risk", "#3478a8"),
        (high, "Higher predicted risk", "#c45a47"),
    ]:
        kmf.fit(
            clinical.loc[mask, "os_months"],
            clinical.loc[mask, "event"],
            label=label,
        )
        kmf.plot_survival_function(ci_show=True, color=color)
    test = logrank_test(
        clinical.loc[high, "os_months"],
        clinical.loc[~high, "os_months"],
        event_observed_A=clinical.loc[high, "event"],
        event_observed_B=clinical.loc[~high, "event"],
    )
    plt.xlabel("Months")
    plt.ylabel("Estimated survival probability")
    plt.title("Combined-model repeated OOF risk groups")
    plt.tight_layout()
    plt.savefig(FIGURES / "oof_risk_groups_km.png", dpi=180)
    plt.close()
    return {
        "combined_oof_median_cutoff": cutoff,
        "combined_oof_logrank_p": float(test.p_value),
    }


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    clinical = pd.read_csv(
        PROCESSED / "clinical_survival.csv", index_col="patient_id"
    )
    expression = pd.read_csv(
        PROCESSED / "expression_log2.csv", index_col="patient_id"
    )
    expression = expression.loc[clinical.index]
    y = make_survival(clinical)
    event = clinical["event"].astype(bool).to_numpy()

    outer = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS,
        n_repeats=OUTER_REPEATS,
        random_state=RANDOM_SEED,
    )
    fold_rows = []
    oof_rows = []
    coefficient_rows = []
    eligibility = Counter()

    for outer_number, (train_idx, test_idx) in enumerate(
        outer.split(expression, event), start=1
    ):
        repeat = (outer_number - 1) // OUTER_SPLITS + 1
        fold = (outer_number - 1) % OUTER_SPLITS + 1
        seed = RANDOM_SEED + outer_number
        clinical_train = clinical.iloc[train_idx]
        clinical_test = clinical.iloc[test_idx]
        expression_train = expression.iloc[train_idx]
        expression_test = expression.iloc[test_idx]

        tuned = {
            "Clinical": tune_clinical(
                clinical_train, y[train_idx], event[train_idx], seed
            ),
            "Expression": tune_coxnet(
                clinical_train,
                expression_train,
                y[train_idx],
                event[train_idx],
                seed,
                combined=False,
            ),
            "Combined": tune_coxnet(
                clinical_train,
                expression_train,
                y[train_idx],
                event[train_idx],
                seed,
                combined=True,
            ),
        }

        for model_name in ["Clinical", "Expression", "Combined"]:
            predictions, coefficients, fitted_names = fit_predict_outer(
                clinical_train,
                clinical_test,
                expression_train,
                expression_test,
                y[train_idx],
                model_name,
                tuned[model_name],
            )
            score = cindex(y[test_idx], predictions)
            fold_rows.append(
                {
                    "repeat": repeat,
                    "fold": fold,
                    "outer_number": outer_number,
                    "model": model_name,
                    "c_index": score,
                    "selected_alpha": tuned[model_name],
                    "test_n": len(test_idx),
                    "test_events": int(event[test_idx].sum()),
                }
            )
            for patient_id, risk_value in zip(
                clinical.index[test_idx], zscore(predictions), strict=True
            ):
                oof_rows.append(
                    {
                        "patient_id": patient_id,
                        "repeat": repeat,
                        "fold": fold,
                        "model": model_name,
                        "risk_z": float(risk_value),
                    }
                )
            if model_name == "Expression":
                for name in fitted_names:
                    eligibility[name.removeprefix("gene:")] += 1
                for name, value in coefficients.items():
                    gene = name.removeprefix("gene:")
                    if abs(value) > 1e-8:
                        coefficient_rows.append(
                            {
                                "outer_number": outer_number,
                                "gene": gene,
                                "coefficient": float(value),
                            }
                        )
        print(
            f"Completed outer fold {outer_number}/{OUTER_SPLITS * OUTER_REPEATS}"
        )

    folds = pd.DataFrame(fold_rows)
    oof = pd.DataFrame(oof_rows)
    coefficients = pd.DataFrame(coefficient_rows)
    summary = summarize_performance(folds)

    selected_counts = coefficients.groupby("gene").size().rename("selected_folds")
    medians = coefficients.groupby("gene")["coefficient"].median().rename("median_coef")
    positive = (
        coefficients.assign(positive=lambda frame: frame["coefficient"] > 0)
        .groupby("gene")["positive"]
        .mean()
        .rename("positive_fraction")
    )
    gene_stability = pd.concat([selected_counts, medians, positive], axis=1)
    gene_stability["eligible_folds"] = [
        eligibility[gene] for gene in gene_stability.index
    ]
    gene_stability["selection_rate"] = (
        gene_stability["selected_folds"] / (OUTER_SPLITS * OUTER_REPEATS)
    )
    gene_stability["selection_rate_if_eligible"] = (
        gene_stability["selected_folds"] / gene_stability["eligible_folds"]
    )
    gene_stability = (
        gene_stability.reset_index()
        .sort_values(
            ["selection_rate", "selection_rate_if_eligible"],
            ascending=False,
        )
    )

    expression_alpha = float(
        folds.loc[folds["model"].eq("Expression"), "selected_alpha"].median()
    )
    permutation = permutation_control(expression, y, event, expression_alpha)

    folds.to_csv(TABLES / "fold_performance.csv", index=False)
    summary.to_csv(TABLES / "performance_summary.csv", index=False)
    oof.to_csv(TABLES / "oof_predictions.csv", index=False)
    coefficients.to_csv(TABLES / "expression_nonzero_coefficients.csv", index=False)
    gene_stability.to_csv(TABLES / "gene_stability.csv", index=False)
    permutation.to_csv(TABLES / "permutation_control.csv", index=False)

    plot_stats = save_figures(
        folds, gene_stability, permutation, clinical, oof
    )
    analysis_summary = {
        "random_seed": RANDOM_SEED,
        "participants": int(len(clinical)),
        "events": int(event.sum()),
        "outer_folds": int(OUTER_SPLITS * OUTER_REPEATS),
        "top_genes_per_training_fold": TOP_GENES,
        "l1_ratio": L1_RATIO,
        "clinical_penalty_weight_in_combined": CLINICAL_PENALTY_WEIGHT,
        "model_mean_cindex": {
            row["model"]: float(row["mean"])
            for row in summary.to_dict(orient="records")
        },
        "permutations": N_PERMUTATIONS,
        "expression_alpha_for_permutation": expression_alpha,
        "permutation_mean_cindex_mean": float(
            permutation.loc[
                permutation["fold"].eq("mean"), "c_index"
            ].mean()
        ),
        "permutation_mean_cindex_min": float(
            permutation.loc[
                permutation["fold"].eq("mean"), "c_index"
            ].min()
        ),
        "permutation_mean_cindex_max": float(
            permutation.loc[
                permutation["fold"].eq("mean"), "c_index"
            ].max()
        ),
        **plot_stats,
    }
    (OUTPUT / "analysis_summary.json").write_text(
        json.dumps(analysis_summary, indent=2), encoding="utf-8"
    )
    print("\nPerformance summary")
    print(summary.to_string(index=False))
    print("\nAnalysis summary")
    print(json.dumps(analysis_summary, indent=2))


if __name__ == "__main__":
    main()
