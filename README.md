# AML Prognostic-Model Reproducibility Pilot

> A leakage-resistant comparison of clinical and gene-expression survival
> models using public TCGA-LAML data.

**Status:** Internal-validation pilot completed · External validation pending  
**Intended use:** Reproducible research and portfolio demonstration  
**Not intended for:** Patient-level prediction or clinical decision-making

## Research question

Does baseline RNA-seq gene expression improve overall-survival risk ranking
beyond a compact diagnosis-time clinical model in acute myeloid leukemia?

## Main finding

Gene expression showed some internal prognostic signal, but did **not** provide a
consistent improvement over the clinical baseline.

| Model | Mean outer-fold C-index | SD | Median | Range |
|---|---:|---:|---:|---:|
| Clinical | **0.695** | 0.058 | 0.698 | 0.602–0.834 |
| Clinical + expression | **0.689** | 0.046 | 0.689 | 0.598–0.805 |
| Expression only | **0.621** | 0.055 | 0.614 | 0.515–0.758 |

The combined model exceeded the clinical model in 10 of 25 held-out folds. Its
mean paired C-index difference was -0.006.

![Repeated outer-fold performance](outputs/figures/model_cindex.png)

## Why this project matters

High-dimensional molecular models can look convincing when only a favorable
split, a pooled risk-group plot, or selected genes are reported. This project
instead emphasizes:

- training-fold-only preprocessing;
- repeated nested cross-validation;
- outcome-label permutation controls;
- feature-selection stability;
- transparent reporting of negative results and analysis-stage decisions.

The contribution is methodological rather than clinical: the workflow can show
when added molecular complexity does not become reliable incremental
prediction.

## Cohort

The analytic cohort contained 173 participants from the public, de-identified
TCGA-LAML study:

- 173 RNA-seq profiles linked to complete overall-survival outcomes;
- 114 observed deaths;
- median age 58 years;
- Kaplan–Meier median survival 18.1 months.

The clinical baseline used age, sex, `log1p(WBC)`, bone-marrow blast percentage,
and cytogenetic risk. Treatment, response, transplantation, and other
post-baseline fields were excluded to reduce temporal leakage.

## Validation design

```text
173 participants
    │
    ├── five-fold outer cross-validation × five repeats
    │       └── 25 held-out fold evaluations
    │
    └── within every outer training set
            ├── three-fold inner cross-validation
            ├── fold-local imputation and encoding
            ├── fold-local gene selection and scaling
            └── regularization tuning
```

Three Cox models were evaluated:

1. a ridge-penalized clinical Cox model;
2. an elastic-net expression Cox model (`l1_ratio = 0.9`);
3. a combined model with clinical predictors receiving 5% of the relative gene
   penalty weight.

The primary metric was Harrell's C-index. It measures risk ranking, not
calibration, treatment benefit, or clinical utility.

## Negative control

The paired survival outcome was permuted 20 times. Each permutation received one
five-fold expression-model evaluation with a fixed regularization value derived
from the observed-outcome analysis.

- observed expression mean C-index: 0.621;
- permutation mean: 0.519;
- permutation range: 0.442–0.626;
- permutations at least as large as observed: 1/20;
- descriptive plus-one proportion: 2/21 = 0.095.

![Outcome-label permutation control](outputs/figures/permutation_control.png)

This is a descriptive negative control, not a fully nested or highly resolved
permutation test.

## Feature stability

IRX1, DOCK1, and AREG were the most frequently retained expression features.
Selection frequency in one cohort is not biomarker validation.

![Expression-feature stability](outputs/figures/gene_stability.png)

## Repository map

| Path | Purpose |
|---|---|
| `research_plan.md` | Prespecified question, predictors, validation, and interpretation rule |
| `research_log.md` | Data decisions, failed attempts, and numerical-stability revision |
| `report.md` | Full English report |
| `report_ko.md` | Korean explanatory report |
| `scripts/00_download_data.py` | Download and verify the source study archive |
| `scripts/01_prepare_data.py` | Link clinical and expression data |
| `scripts/02_validate_models.py` | Nested validation and negative control |
| `scripts/03_summarize_results.py` | Aggregate result tables |
| `scripts/99_quality_checks.py` | Artifact-integrity checks |
| `outputs/tables/` | Non-patient-level result tables |
| `outputs/figures/` | Main result figures |

## Reproduce

### 1. Create the environment

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2. Download the source study

```powershell
python scripts\00_download_data.py
```

The script downloads the public cBioPortal `laml_tcga_pub` archive, verifies its
recorded SHA-256 hash, and extracts it under `data/raw/`.

### 3. Run the workflow

```powershell
.\run_all.ps1
```

Or run each analysis stage explicitly:

```powershell
python scripts\01_prepare_data.py
python scripts\02_validate_models.py
python scripts\03_summarize_results.py
python scripts\99_quality_checks.py
```

## Data availability

Patient-level source and derived files are not redistributed in this repository.
The repository contains the source URL, file hashes, preparation code, aggregate
tables, and figures.

Source archive:

```text
https://datahub.assets.cbioportal.org/laml_tcga_pub.tar.gz
```

The results shown here are in whole or part based upon data generated by the
TCGA Research Network: <https://www.cancer.gov/tcga>.

## Interpretation

The justified conclusion is:

> Some internal molecular prognostic signal was detectable, but incremental
> value beyond the clinical baseline was not demonstrated.

The results do not establish a deployable model, a validated biomarker, an
actionable threshold, or treatment benefit.

## Limitations and next step

- single retrospective cohort;
- 173 participants and 114 events;
- internal validation only;
- compact rather than mutation-aware clinical baseline;
- no calibration or decision-curve analysis;
- descriptive 20-permutation negative control;
- no assay harmonization or locked external validation.

The next phase is external validation in an independent AML expression cohort,
with platform harmonization and no outcome-based external retuning.

## Portfolio and AI disclosure

This repository is a retrospective portfolio publication of a completed
AI-assisted pilot. Git was not used prospectively throughout the original
analysis, and the commit history begins with the portfolio publication. AI
assistance supported analysis design, code generation, debugging,
documentation, and interpretation. See `PORTFOLIO_DISCLOSURE.md`.

The repository does not imply faculty supervision. Any academic description
must reflect the work, understanding, and supervision that actually occurred.

## License and source terms

No license is granted for third-party TCGA or cBioPortal data. Source data remain
subject to their original terms. In the absence of a repository software
license, default copyright applies to original repository content.

