# Quality-assurance record for public release

**Release date:** 2026-10-20

## Completed checks (no MIMIC patient records placed in repository)

| Check | Result |
|---|---|
| G10 v1.1 SHA-256 matches frozen original | PASS |
| G11 Table 1 SHA-256 matches frozen original | PASS |
| Python source compile for all repository scripts | PASS |
| `run_pipeline.py --self-test-only`, G10 | PASS |
| `run_pipeline.py --self-test-only`, G11 | PASS |
| Repository's offline unit tests | 5/5 PASS |
| Verifier against authors' previously completed G10 aggregate output | PASS |
| Verifier against authors' previously completed G11 aggregate output | PASS |
| Synthetic verifier fixture flags deliberately corrupted cohort count | PASS |

The private G10 reference run reported 2260 eligible ICU stays, 2168 unique
patients, 204/139 unambiguous qualifying first actions, and 500 of
500 successful subject-cluster bootstrap replicates. The reference primary
risks were 0.197381 vs 0.184511; RD=+0.012870, RR=1.069751.

The private G11 reference run reported the same 2260 ICU stays/2168
subjects and its Table 1 export was found and parsed.

This was a **verification of archived aggregate reference outputs**.
The full database analysis was **not re-executed inside the packaging
environment**, which does not contain credentialed MIMIC-IV v3.1 data
or a local PostgreSQL installation with the source database.

## What downstream researchers should provide

1. Proof of authorized MIMIC-IV v3.1 access (not copied into repo).
2. PostgreSQL/database version and MIT-LCP derived-concept SQL commit.
3. `run_environment.json` generated locally (do not publish without review).
4. Pass/fail verification summary (carefully screened before sharing).
5. Scientific limitations noted in `docs/FROZEN_PROTOCOL.md`.

Do not upload restricted data or patient-level derived records as "proof"
of reproducibility.
