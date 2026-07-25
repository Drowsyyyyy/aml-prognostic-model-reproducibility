"""Create compact descriptive and model-comparison result tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from lifelines import KaplanMeierFitter
from sksurv.metrics import concordance_index_censored


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "outputs" / "tables"


def median_iqr(series: pd.Series) -> str:
    return (
        f"{series.median():.1f} "
        f"({series.quantile(0.25):.1f}–{series.quantile(0.75):.1f})"
    )


def main() -> None:
    clinical = pd.read_csv(
        PROCESSED / "clinical_survival.csv", index_col="patient_id"
    )
    folds = pd.read_csv(TABLES / "fold_performance.csv")
    oof = pd.read_csv(TABLES / "oof_predictions.csv")
    permutation = pd.read_csv(TABLES / "permutation_control.csv")

    km = KaplanMeierFitter().fit(clinical["os_months"], clinical["event"])
    rows = [
        ("Participants", str(len(clinical))),
        ("Deaths", f"{int(clinical['event'].sum())} ({100 * clinical['event'].mean():.1f}%)"),
        ("Censored", f"{int((~clinical['event']).sum())} ({100 * (~clinical['event']).mean():.1f}%)"),
        ("Age, years", median_iqr(clinical["age"])),
        ("Male", f"{int((clinical['sex'] == 'Male').sum())} ({100 * (clinical['sex'] == 'Male').mean():.1f}%)"),
        ("Female", f"{int((clinical['sex'] == 'Female').sum())} ({100 * (clinical['sex'] == 'Female').mean():.1f}%)"),
        ("WBC", median_iqr(clinical["wbc"])),
        ("Bone-marrow blasts, %", median_iqr(clinical["bm_blast_pct"])),
        ("Cytogenetic risk: good", f"{int((clinical['risk_cyto'] == 'Good').sum())}"),
        ("Cytogenetic risk: intermediate", f"{int((clinical['risk_cyto'] == 'Intermediate').sum())}"),
        ("Cytogenetic risk: poor", f"{int((clinical['risk_cyto'] == 'Poor').sum())}"),
        ("Cytogenetic risk: unknown", f"{int((clinical['risk_cyto'] == 'Unknown').sum())}"),
        ("Kaplan–Meier median survival, months", f"{km.median_survival_time_:.1f}"),
        ("Estimated 12-month survival", f"{float(km.predict(12)):.3f}"),
        ("Estimated 24-month survival", f"{float(km.predict(24)):.3f}"),
    ]
    pd.DataFrame(rows, columns=["characteristic", "value"]).to_csv(
        TABLES / "cohort_summary.csv", index=False
    )

    pivot = folds.pivot(index="outer_number", columns="model", values="c_index")
    comparisons = []
    for first, second in [
        ("Combined", "Clinical"),
        ("Expression", "Clinical"),
        ("Combined", "Expression"),
    ]:
        difference = pivot[first] - pivot[second]
        comparisons.append(
            {
                "comparison": f"{first} minus {second}",
                "mean_cindex_difference": float(difference.mean()),
                "median_cindex_difference": float(difference.median()),
                "folds_first_better": int((difference > 0).sum()),
                "total_folds": int(len(difference)),
            }
        )
    pd.DataFrame(comparisons).to_csv(
        TABLES / "paired_model_comparisons.csv", index=False
    )

    average_oof = (
        oof.groupby(["patient_id", "model"])["risk_z"].mean().unstack()
    )
    pooled_rows = []
    for model in ["Clinical", "Combined", "Expression"]:
        aligned = clinical.loc[average_oof.index]
        score = concordance_index_censored(
            aligned["event"].astype(bool),
            aligned["os_months"],
            average_oof[model],
        )[0]
        pooled_rows.append({"model": model, "pooled_oof_c_index": float(score)})
    pd.DataFrame(pooled_rows).to_csv(
        TABLES / "pooled_oof_performance.csv", index=False
    )

    observed = float(
        folds.loc[folds["model"].eq("Expression"), "c_index"].mean()
    )
    permutation_means = permutation.loc[
        permutation["fold"].astype(str).eq("mean"), "c_index"
    ]
    exceedances = int((permutation_means >= observed).sum())
    negative_control = pd.DataFrame(
        [
            {
                "observed_expression_mean_cindex": observed,
                "permutation_mean_cindex": float(permutation_means.mean()),
                "permutation_min": float(permutation_means.min()),
                "permutation_max": float(permutation_means.max()),
                "permutations_at_least_observed": exceedances,
                "permutations": int(len(permutation_means)),
                "descriptive_plus_one_p": float(
                    (exceedances + 1) / (len(permutation_means) + 1)
                ),
            }
        ]
    )
    negative_control.to_csv(
        TABLES / "negative_control_summary.csv", index=False
    )

    print("Created cohort_summary.csv")
    print("Created paired_model_comparisons.csv")
    print("Created pooled_oof_performance.csv")
    print("Created negative_control_summary.csv")


if __name__ == "__main__":
    main()
