# Reproducibility, provenance, and interpretation notes

## Source-code provenance

The copies in `analysis/` are byte-for-byte identical to the frozen local
G10 v1.1 and G11 scripts that generated the manuscript analyses. The
wrapper `run_pipeline.py` controls the DB connection and output path by
setting *module-level configuration values* without editing any analysis
function, treatment rule, outcome definition, estimator, or bootstrap logic.

`source_sha256.json` is the integrity record. Verify it with:

```bash
python check_integrity.py
```

## Scope of verification

In an environment without credentialed local MIMIC-IV data, a user can:

- Compile every Python file.
- Run built-in G10/G11 synthetic edge-case tests.
- Validate the expected-reference consistency.
- Test output verification on aggregate synthetic fixtures.

The **live database extraction**, frozen cohort reconstruction, model fitting,
and all 500 subject-cluster bootstrap fits require a credentialed local copy
of the real dataset and cannot be honestly claimed as retested without it.

## What is expected to reproduce

- G10 2260 eligible ICU stays; 2168 unique subjects.
- Observed qualifying actions within 120 min: 204 NE-first, 139 AVP-first,
  10 exact simultaneous, 1907 no action.
- Primary untruncated 6-h risks near 0.197381 (NE) and 0.184511 (AVP).
- RD near +0.012870; RR near 1.069751.
- Bootstrap 500/500 successful with seed 20261015.

G11 redoes the same frozen cohort reconstruction for the manuscript
Table 1; it is not a second intervention analysis.

**Tolerance policy**: exact counts are required. Point-estimate tolerances
allow for software/derived-concept variation; bootstrap confidence limits
have looser tolerances. Passing checks is evidence of consistent output,
not proof of exact bitwise computational identity.

## Interpretive caveats

- The original analysis was developed iteratively. It was frozen before
  final bootstrap inference, but it was **not prospectively preregistered**.
- A G10 cohort mismatch is not something to circumvent by weakening
  the hardcoded gates: first verify MIMIC 3.1 and derived concepts.
- 204/139 counts are observed strategy *support*, not sample sizes of
  randomized or naturally observed arms in the CCW analysis.
- No superiority or equivalence was established.
- 12/24-h analyses are exploratory due to natural censoring.
