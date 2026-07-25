# Research log

## 2026-07-24 — Data acquisition and cohort definition

- Downloaded the public `laml_tcga_pub` study archive from the cBioPortal
  DataHub.
- Retained the distributor-provided license text.
- Recorded SHA-256 hashes of the clinical, sample, expression, and license files.
- Found 200 clinical participants and 173 RNA-seq profiles.
- Linked all 173 expression profiles to unique participants.
- Confirmed complete overall-survival time and status in the 173-person analytic
  cohort.
- Selected diagnosis-time clinical variables only; excluded treatment,
  transplantation, response, and follow-up variables.

## 2026-07-24 — Prespecified analysis

- Defined overall survival as the endpoint.
- Defined three comparisons: clinical, expression, and combined.
- Selected repeated nested cross-validation and Harrell C-index.
- Required fold-local preprocessing and gene selection.
- Added a permuted-outcome negative control and coefficient-stability analysis.

## 2026-07-24 — Implementation notes

- The original table export preserved inconsistent index labels. The preparation
  script was corrected so both derived matrices use `patient_id`.
- The initial combined model penalized clinical and expression coefficients
  equally. This was recognized as an unfair incremental-value comparison because
  the clinical baseline could be shrunk together with thousands of gene
  candidates.
- A second specification attempted to leave clinical coefficients unpenalized.
  One inner-training fold produced numerical separation and Coxnet failed.
- Final specification: clinical coefficients receive 5% of the expression
  penalty weight. This retained a largely protected baseline while resolving
  the numerical instability.
- The permutation splitter was aligned to the permuted event indicator rather
  than the original event indicator.
- All final results were regenerated after these changes.

## 2026-07-24 — Final results

- Clinical mean outer-fold C-index: 0.695.
- Combined mean outer-fold C-index: 0.689.
- Expression mean outer-fold C-index: 0.621.
- Combined model exceeded clinical performance in 10 of 25 folds.
- Expression model exceeded clinical performance in 4 of 25 folds.
- Negative-control permutation mean: 0.519; range 0.442–0.626.
- One of 20 permutation means equaled or exceeded the observed expression-model
  mean.
- Final automated integrity checks passed.

## Interpretation decision

The project reports no demonstrated incremental value from adding expression to
the clinical baseline. It does not label the frequently selected genes as
validated biomarkers and does not describe internal cross-validation as
external or clinical validation.

## Work required before legitimate academic use

- Reproduce the workflow personally from the raw archive.
- Read and be able to explain every analysis decision and code block.
- Have a real supervisor review the clinical question, methods, and
  interpretation.
- Record the supervisor's actual contributions prospectively.
- Add external validation before making a generalizability claim.
- Disclose AI assistance according to the relevant institution's policy.

