# Prespecified research plan

## Question

In the public TCGA-LAML cohort, how stable is the apparent prognostic signal
from baseline gene expression when evaluated with leakage-resistant repeated
cross-validation, and does it improve risk ranking beyond a compact clinical
baseline?

## Cohort and endpoint

- Inclusion: a TCGA-LAML participant with RNA-seq expression, overall-survival
  time, and overall-survival status.
- Primary endpoint: overall survival in months.
- Event: death (`OS_STATUS = 1:DECEASED`).
- Time zero: diagnosis/baseline as represented by the source clinical file.

## Predictors

Clinical baseline:

- age;
- sex;
- white-blood-cell count, transformed as `log1p(WBC)`;
- bone-marrow blast percentage;
- cytogenetic risk.

Expression predictors:

- log2-transformed RSEM abundance (`log2(x + 1)`);
- blank gene symbols removed;
- duplicated gene symbols collapsed by their mean;
- the 250 most variable genes selected within each training fold.

Treatment, transplantation, follow-up, and post-baseline fields are excluded to
reduce temporal leakage.

## Missing data

- Missing endpoint time or event: exclude from the analytic cohort.
- Numeric clinical predictor: training-fold median imputation.
- Categorical clinical predictor: training-fold most-frequent imputation,
  followed by one-hot encoding.
- Expression zero is treated as an observed value, not as missing.
- Any expression missingness is imputed from the training-fold gene median.

In the downloaded version of the source data, the selected clinical fields and
survival endpoint are complete among expression-profiled participants, but the
pipeline retains fold-local imputation so the rule is explicit and reusable.

## Models

1. Penalized Cox proportional-hazards model with clinical predictors.
2. Elastic-net Cox model with gene-expression predictors.
3. Elastic-net Cox model with clinical and gene-expression predictors. Clinical
   coefficients receive 5% of the gene-expression penalty weight. This preserves
   the clinical baseline relative to the genes while avoiding numerical
   separation in small training folds.

For the expression models, `l1_ratio = 0.9`. The regularization strength is
selected by three-fold inner cross-validation.

## Validation

- Outer validation: five-fold cross-validation repeated five times.
- Inner tuning: three-fold cross-validation within each outer training set.
- Stratification: death/censoring indicator.
- Primary metric: Harrell concordance index (C-index).
- Report: mean, standard deviation, median, and range across outer folds.
- Negative control: permute the paired survival outcome and repeat five-fold
  expression-model validation 20 times using a fixed regularization value
  derived from the real-outcome analysis.
- Feature stability: frequency and sign consistency of non-zero expression
  coefficients across outer folds.

## Interpretation rule

The analysis will not call a model clinically useful based on a C-index alone.
Any apparent improvement is described as internal predictive evidence that
requires independent-cohort validation, calibration assessment, comparison with
current clinical standards, and prospective evaluation.
