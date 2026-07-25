"""Fail-fast integrity checks for generated analysis artifacts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"


def main() -> None:
    clinical = pd.read_csv(
        PROCESSED / "clinical_survival.csv", index_col="patient_id"
    )
    expression = pd.read_csv(
        PROCESSED / "expression_log2.csv", index_col="patient_id"
    )
    folds = pd.read_csv(TABLES / "fold_performance.csv")
    oof = pd.read_csv(TABLES / "oof_predictions.csv")
    performance = pd.read_csv(TABLES / "performance_summary.csv")

    assert len(clinical) == 173
    assert clinical.index.is_unique and expression.index.is_unique
    assert clinical.index.equals(expression.index)
    assert clinical[["os_months", "event"]].notna().all().all()
    assert (clinical["os_months"] >= 0).all()
    assert len(folds) == 25 * 3
    assert set(folds["model"]) == {"Clinical", "Expression", "Combined"}
    assert folds["c_index"].between(0, 1).all()
    assert len(performance) == 3
    assert not performance.isna().any().any()

    counts = oof.groupby(["patient_id", "model"]).size()
    assert counts.eq(5).all()
    assert oof["risk_z"].notna().all()

    required_figures = {
        "model_cindex.png",
        "gene_stability.png",
        "permutation_control.png",
        "oof_risk_groups_km.png",
    }
    found = {path.name for path in FIGURES.glob("*.png") if path.stat().st_size > 0}
    assert required_figures.issubset(found)

    print("All quality checks passed.")
    print(f"Clinical-expression rows aligned: {len(clinical)}")
    print(f"Outer model evaluations: {len(folds)}")
    print(f"OOF predictions: {len(oof)}")


if __name__ == "__main__":
    main()

