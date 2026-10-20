# Frozen analysis specification

This document describes **the reported analysis**, not a newly designed target
trial. The design was frozen for final G10 bootstrap inference **after earlier
methodological/pilot analyses**; it was not prospectively preregistered.

## Eligibility and common time zero

Adult ICU stays with Sepsis-3 sepsis onset by the candidate time, continuously
overlapping norepinephrine (NE) and vasopressin (AVP) for at least 120 minutes.
Select the first 15-minute-grid candidate meeting the following *pre-action*
criteria:

- MAP stability in the previous 120 min: at least 3 measurements, spanning
  at least 90 min, most recent within 30 min, all MAP >=65 mmHg.
- Pre-action left-limit NE <=0.30 mcg/kg/min and AVP <=3.60 U/h.
- No meaningful escalation in the previous 120 min.
- No exact-t0 NE/AVP rate transition.
- No single DOVSS-scale NE/AVP rate transition during the previous 120 min.
- No other non-study vasopressor active at t0 or during previous 120 min.

Operational decision-point name:
`dovss_stable120_no_otherpressor2h`.

> The 120-minute DOVSS-scale stability rule is an **operational
> analogue**, not identical to DOVSS trial eligibility: it permits
> subthreshold infusion microtitrations rather than demanding a fully
> constant infusion for the entire 2 hours.

> Clinical limitation: Sepsis-3-derived sepsis plus concurrent NE/AVP is an
> operational septic-shock recovery phenotype, not confirmation that every
> patient meets the entire bedside Sepsis-3 shock definition including adequate
> fluid resuscitation.

## Strategies

Each eligible ICU stay is cloned to both:

1. **NE-first:** first meaningful observed decrease from t0 NE dose is
   >=0.10 mcg/kg/min, or complete NE discontinuation.
2. **AVP-first:** first meaningful observed decrease from t0 AVP dose is
   >=0.60 U/h (0.01 U/min), or complete AVP discontinuation.

Primary grace window = 120 min, no future persistence criterion.
An exact simultaneous first qualifying NE and AVP transition is unresolved.
Before divergence, both strategies remain compatible with observed care;
an outcome during compatibility can contribute to both clone arms.
Once assigned first action is completed, no further *artificial* censoring
is applied: the intervention is first action, **not complete drug cessation**.

## Outcome

First occurrence of **two MAP values <65 mmHg separated by 15–60 min**
within 6 h of the common time zero.

Twelve- and 24-hour windows are exploratory because natural ICU discharge/death
censoring becomes more consequential and is not additionally IPC weighted.

## Estimation

- Artificial censoring in 15-min intervals; strategy-specific logistic models.
- Primary stabilized numerator: interval time only.
- Denominator: baseline and time-varying covariates as implemented in G10.
- Probability clip = [0.01, 0.999].
- Primary: untruncated stabilized IPCW + weighted Kaplan–Meier.
- Contrasts: NE-first minus AVP-first RD; NE-first / AVP-first RR.
- 500 bootstrap replicates at `subject_id` level, seed `20261015`, refit
  censor models per replicate; include all eligible ICU stays and both clones
  for a resampled subject.

## Sensitivity analyses

- 1st/99th percentile cumulative weight truncation.
- DOVSS action, 240-min grace.
- Immediate cumulative >=30% relative reduction, 120-min grace.
- Baseline-stabilized weighted outcome model + standardization, point estimate.

The reference 6-h outcome risks were 19.74% for NE-first and 18.45% for
AVP-first; RD +1.287 percentage points (95% bootstrap CI -8.071 to +10.378).
**These data do not demonstrate strategy equivalence.**
