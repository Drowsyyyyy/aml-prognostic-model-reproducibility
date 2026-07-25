"""Prepare a linked clinical-expression TCGA-LAML analysis cohort."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "laml_tcga_pub"
PROCESSED = ROOT / "data" / "processed"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_cbioportal_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", low_memory=False)


def main() -> None:
    required = [
        RAW / "data_clinical_patient.txt",
        RAW / "data_clinical_sample.txt",
        RAW / "data_mrna_seq_v2_rsem.txt",
        RAW / "LICENSE",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required source files:\n" + "\n".join(missing))

    PROCESSED.mkdir(parents=True, exist_ok=True)

    patients = read_cbioportal_table(required[0])
    samples = read_cbioportal_table(required[1])
    expression_raw = read_cbioportal_table(required[2])

    sample_to_patient = (
        samples[["SAMPLE_ID", "PATIENT_ID"]]
        .drop_duplicates()
        .set_index("SAMPLE_ID")["PATIENT_ID"]
    )
    expression_samples = [
        column
        for column in expression_raw.columns[2:]
        if column in sample_to_patient.index
    ]

    expression_raw["Hugo_Symbol"] = (
        expression_raw["Hugo_Symbol"].astype("string").str.strip()
    )
    expression_raw = expression_raw[
        expression_raw["Hugo_Symbol"].notna()
        & expression_raw["Hugo_Symbol"].ne("")
    ].copy()
    numeric_expression = expression_raw[expression_samples].apply(
        pd.to_numeric, errors="coerce"
    )
    numeric_expression.insert(0, "Hugo_Symbol", expression_raw["Hugo_Symbol"])
    expression_by_gene = numeric_expression.groupby("Hugo_Symbol").mean()

    expression = np.log2(expression_by_gene.T + 1.0)
    expression.index = [sample_to_patient.loc[sample] for sample in expression.index]
    expression.index.name = "patient_id"
    if expression.index.duplicated().any():
        expression = expression.groupby(level=0).mean()

    clinical = patients.set_index("PATIENT_ID").copy()
    selected = pd.DataFrame(index=clinical.index)
    selected["os_months"] = pd.to_numeric(clinical["OS_MONTHS"], errors="coerce")
    selected["event"] = clinical["OS_STATUS"].astype("string").str.startswith("1")
    selected.loc[clinical["OS_STATUS"].isna(), "event"] = pd.NA
    selected["age"] = pd.to_numeric(clinical["AGE"], errors="coerce")
    selected["sex"] = clinical["SEX"].astype("string").str.strip()
    selected["wbc"] = pd.to_numeric(clinical["WBC"], errors="coerce")
    selected["log1p_wbc"] = np.log1p(selected["wbc"])
    selected["bm_blast_pct"] = pd.to_numeric(
        clinical["BM_BLAST_PERCENTAGE"], errors="coerce"
    )
    selected["risk_cyto"] = (
        clinical["RISK_CYTO"]
        .astype("string")
        .str.strip()
        .replace({"N.D.": "Unknown", "": pd.NA})
    )

    endpoint_complete = selected["os_months"].notna() & selected["event"].notna()
    eligible_ids = expression.index.intersection(selected.index[endpoint_complete])
    cohort = selected.loc[eligible_ids].copy()
    cohort.index.name = "patient_id"
    cohort["event"] = cohort["event"].astype(bool)
    expression = expression.loc[eligible_ids]
    expression.index.name = "patient_id"

    if not cohort.index.equals(expression.index):
        raise RuntimeError("Clinical and expression row order does not match.")
    if cohort.index.duplicated().any() or expression.index.duplicated().any():
        raise RuntimeError("Analytic cohort contains duplicated participant IDs.")

    cohort.to_csv(PROCESSED / "clinical_survival.csv")
    expression.to_csv(PROCESSED / "expression_log2.csv")

    missingness = pd.DataFrame(
        {
            "variable": cohort.columns,
            "missing_n": [int(cohort[c].isna().sum()) for c in cohort.columns],
            "missing_pct": [
                float(100.0 * cohort[c].isna().mean()) for c in cohort.columns
            ],
        }
    )
    missingness.to_csv(PROCESSED / "missingness.csv", index=False)

    data_dictionary = pd.DataFrame(
        [
            ("patient_id", "Public de-identified TCGA participant identifier"),
            ("os_months", "Overall-survival or censoring time in months"),
            ("event", "True if death was observed; False if censored"),
            ("age", "Age at diagnosis in years"),
            ("sex", "Sex recorded in the source file"),
            ("wbc", "Baseline white-blood-cell count from the source file"),
            ("log1p_wbc", "Natural log of one plus WBC"),
            ("bm_blast_pct", "Baseline bone-marrow blast percentage"),
            ("risk_cyto", "Cytogenetic risk; N.D. relabeled Unknown"),
            ("expression_log2", "RSEM abundance transformed as log2(x + 1)"),
        ],
        columns=["field", "definition"],
    )
    data_dictionary.to_csv(PROCESSED / "data_dictionary.csv", index=False)

    manifest = {
        "source": "cBioPortal DataHub laml_tcga_pub",
        "download_url": "https://datahub.assets.cbioportal.org/laml_tcga_pub.tar.gz",
        "source_files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in required
        },
        "source_patient_rows": int(len(patients)),
        "source_sample_rows": int(len(samples)),
        "expression_sample_columns": int(len(expression_samples)),
        "unique_expression_genes": int(expression.shape[1]),
        "analytic_participants": int(len(cohort)),
        "deaths": int(cohort["event"].sum()),
        "censored": int((~cohort["event"]).sum()),
        "median_os_months": float(cohort["os_months"].median()),
    }
    (PROCESSED / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
