# MIMIC-IV data access and safe research-code release

## Credentialing and database

MIMIC-IV v3.1 is PhysioNet **credentialed health data** and is not included
in this repository. Access requires training, credentialing and the applicable
data use agreement (DUA).

- Dataset: https://physionet.org/content/mimiciv/3.1/
- DUA: https://physionet.org/content/mimiciv/view-dua/3.1/
- Original PostgreSQL build scripts:
  https://github.com/MIT-LCP/mimic-code/tree/main/mimic-iv/buildmimic/postgres
- Derived-concept PostgreSQL SQL:
  https://github.com/MIT-LCP/mimic-code/tree/main/mimic-iv/concepts_postgres

Local DB must contain `mimiciv_icu`, `mimiciv_hosp`, and the derived
`mimiciv_derived` schema with at least `sepsis3` and `vitalsign`.

**For exact reproduction, record the commit/version of the official
derived SQL used to build `mimiciv_derived`, as upstream definitions
can change.** The source's hard cohort-count gates may legitimately
fail for a different source dataset or derived-concept version.

## What is safe to publish

This release includes frozen **source code only**, plus an
author-supplied, non-patient-level reference result summary. It intentionally
excludes all data extracts, MIMIC tables, per-stay records, patient IDs,
patient-level derivatives, private outputs, database credentials and `.env`.

- Do **NOT** upload raw MIMIC-IV data or individual-level derivatives to
  public GitHub or other public services.
- Treat any MIMIC-derived datasets and models as sensitive; sharing
  derived datasets should follow PhysioNet's restricted-access pathway.
- Review even aggregate reports before any further public dissemination.
- Never commit `outputs/`, database files, credentials or logs.
- Do not paste MIMIC restricted records into external AI tools or services.
- The original analysis uses subject identifiers **internally** for
  clustering, but does not export those identifiers.

PhysioNet DUA clause 9 calls for publicly contributing analysis code when
openly disseminating results. That is why the code is released separately
from the credentialed dataset.

## Code licensing

No software license is automatically granted here on the coauthors' behalf.
Before making the GitHub repository public, the authors should jointly choose
an appropriate code license (e.g., MIT or Apache-2.0) and add `LICENSE`.
