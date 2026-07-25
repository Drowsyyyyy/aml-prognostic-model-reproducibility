# Data directory

Expected input:

```text
data/raw/laml_tcga_pub/
  data_clinical_patient.txt
  data_clinical_sample.txt
  data_mrna_seq_v2_rsem.txt
  LICENSE
```

Download archive:

`https://datahub.assets.cbioportal.org/laml_tcga_pub.tar.gz`

Recorded archive SHA-256:

`be0b4bb8a0481acd6e4127f4e0a15889de9eb7fd4bc6e663cc90f65db419b0b9`

The public repository includes `scripts/00_download_data.py` to download,
verify, and extract this recorded source version.

The preparation script records file hashes and cohort counts in
`data/processed/manifest.json`. Raw and processed patient-level data are excluded
from version control.
