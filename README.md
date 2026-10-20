# NE-first versus AVP-first initial de-escalation: target trial emulation

**Reproducibility code for a MIMIC-IV v3.1 observational
clone–censor–weight (CCW) target trial emulation.**

Authors: **Zheng Yao**, **Zhaotao Fang**, **Shijia Zhu**.

> This repository contains **analysis code, documentation, offline
> synthetic tests, and published aggregate reference metrics only**.
> It does **not** include MIMIC-IV data or MIMIC-derived patient-level
> datasets. Access to the underlying data is restricted by PhysioNet.

## Research question

Among adult patients with Sepsis-3 sepsis receiving concurrent
norepinephrine (NE) and vasopressin (AVP) who have reached a stable
pre-weaning decision point, does reducing **NE first** versus
**AVP first** change the risk of persistent hypotension in the next
**6 hours**?

This analysis concerns the **first clinically meaningful reduction**,
**not the order of complete drug discontinuation**.

## Reference primary results

| Quantity | Frozen G10 analysis |
|---|---:|
| Eligible ICU stays | 2,260 |
| Unique patients | 2,168 |
| NE-first estimated 6-hour risk | 19.74% |
| AVP-first estimated 6-hour risk | 18.45% |
| Risk difference (NE minus AVP) | +1.29 percentage points |
| RD 95% bootstrap CI | −8.07 to +10.38 percentage points |
| Risk ratio | 1.07 |
| RR 95% bootstrap CI | 0.68 to 1.79 |
| Subject-level bootstrap | 500/500 successful; seed 20261015 |

**No clear advantage was observed; this does not demonstrate equivalence.**
The observed 204 NE-first and 139 AVP-first qualifying actions are measures
of empirical adherence/support within the 120-minute grace period, **not**
the emulated arm sizes: each of the 2,260 eligible stays was cloned to
both strategies.

## Repo contents

```text
analysis/
  G10_v11_MIMIC_NE_AVP_final_CCW_subject_bootstrap.py   # frozen main analysis
  G11_MIMIC_NE_AVP_manuscript_Table1.py               # frozen Table 1
run_pipeline.py                  # convenience runner; analysis code unchanged
preflight.py                     # read-only PostgreSQL schema checks
verify_results.py               # checks private outputs against reference
check_integrity.py              # SHA-256 check of frozen analysis scripts
expected_results.json
source_sha256.json
requirements.txt
tests/
docs/
CITATION.cff
```

## Prerequisites

- **Authorized, credentialed** access to **MIMIC-IV v3.1** under
  PhysioNet's training and data use requirements.
- Local PostgreSQL containing raw schemas `mimiciv_icu`, `mimiciv_hosp`
  and derived concepts `mimiciv_derived.sepsis3` and
  `mimiciv_derived.vitalsign`.
- Python **3.11+** (3.11–3.13 is intended; inspect the captured environment
  metadata for your run).
- Sufficient local compute/time for 500 subject-bootstrap fits.

For PostgreSQL installation and derived SQL, use the official
[MIT-LCP MIMIC Code Repository](https://github.com/MIT-LCP/mimic-code).
The MIMIC-IV data files are not distributed here. For exact cohort
reproduction, the official **derived-concept SQL commit/version** matters
and should be recorded by the reproducing investigator.

## Quick start

From the repository root (bash/macOS/Linux; analogous commands
work in Windows PowerShell):

**1. Create a private Python environment:**

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

**2. Confirm unchanged source and run fully offline synthetic tests:**

```bash
python check_integrity.py
python run_pipeline.py --self-test-only
python -m unittest discover -s tests -v
```

These tests **do not access MIMIC data or connect to PostgreSQL**.

**3. Check the authorized database (read only):**

```bash
python preflight.py \
  --db-name mimiciv \
  --db-host /var/run/postgresql \
  --db-port 5432 \
  --db-user "$(whoami)"
```

The defaults match the authors' local Unix-socket setup. With a local
network PostgreSQL server, use `--db-host 127.0.0.1` and your own
`--db-user`. For secure password authentication use PostgreSQL's
usual `PGPASSWORD` or a local `.pgpass` **outside the repository**;
never commit a password.

**4. Run the complete manuscript analysis and Table 1:**

```bash
python run_pipeline.py \
  --step all \
  --db-name mimiciv \
  --db-host /var/run/postgresql \
  --n-bootstrap 500 \
  --seed 20261015
```

By default, results are written to:
`outputs/<timestamp>/g10/` and `outputs/<timestamp>/g11/`.
These folders are **private local files**, excluded by `.gitignore`.
The pipeline checks expected output automatically and writes
`run_environment.json` with dependency versions and script hashes.

To run one step only:

```bash
python run_pipeline.py --step g10  # 500-bootstrap primary inference
python run_pipeline.py --step g11  # baseline Table 1, no bootstrap
```

**5. Verify existing results without rerunning the database analysis:**

```bash
python verify_results.py \
  --g10-dir outputs/<timestamp>/g10 \
  --g11-dir outputs/<timestamp>/g11
```

You can also point the verifier at **your pre-existing local**
G10 v1.1 and G11 output folders. Do **not** upload those folders to
the public repository.

**Jupyter alternative**, from the repository root:

```python
import subprocess, sys

subprocess.run(
    [sys.executable, "run_pipeline.py", "--self-test-only"],
    check=True,
)
# Run the real database analysis when local prerequisites are satisfied:
# subprocess.run([sys.executable, "run_pipeline.py", "--step", "all"], check=True)
```

## Scientific methods and disclosure

See [frozen protocol](docs/FROZEN_PROTOCOL.md) for exact eligibility,
common time zero, action thresholds, censoring, model and outcome
definitions. See [reproducibility notes](docs/REPRODUCIBILITY_NOTES.md)
for limitations of cross-environment replay.

This observational analysis was refined in a series of design-stage
and pilot estimator audits before the **final** G10 bootstrap analysis
was frozen. It was **not prospectively registered**. G11 uses the same
common eligible cohort to export manuscript Table 1.

Since this is a historical frozen pipeline, some intermediate diagnostic
constants and comments remain inside the unchanged source for provenance.
The manuscript estimand is the frozen G10 **primary DOVSS-scale 120-min
grace / untruncated marginal IPCW weighted-KM / 6-h outcome**.

## Data use / GitHub safety

**Publish the code, not the data.** See
[data access and privacy](docs/DATA_ACCESS_AND_PRIVACY.md).

PhysioNet:
- [MIMIC-IV v3.1](https://physionet.org/content/mimiciv/3.1/)
- [Credentialed Health Data Use Agreement](https://physionet.org/content/mimiciv/view-dua/3.1/)

MIMIC-IV citation:
Johnson A, Bulgarelli L, Pollard T, et al. **MIMIC-IV (version 3.1)**.
PhysioNet, 2024. DOI: [10.13026/kpb9-mt58](https://doi.org/10.13026/kpb9-mt58).

Before releasing publicly, review
[RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) and jointly select a
suitable software license. No data-sharing permission is granted by
this code repository.
