# G10 MIMIC-IV NE/AVP final frozen-protocol CCW + subject-cluster bootstrap
# Version 1.1 — 2026-10-15
#
# READ-ONLY. Patient identifiers are used internally for clustering but are NEVER
# exported. All saved outputs are aggregate or bootstrap-replicate summaries.
#
# FROZEN PRIMARY TARGET TRIAL
# ---------------------------
# Population / decision point:
#   First time at which ALL are true:
#     - adult
#     - Sepsis-3 onset already present
#     - continuous NE + AVP for >=2 h
#     - MAP stable during prior 2 h
#     - LEFT-LIMIT (pre-action) NE <=0.30 mcg/kg/min
#     - LEFT-LIMIT AVP <=3.60 U/h
#     - no clinically meaningful escalation in prior 2 h
#     - no NE/AVP rate transition exactly at t0
#     - no single DOVSS-scale NE/AVP change during prior 2 h
#         NE absolute transition >=0.10 mcg/kg/min
#         AVP absolute transition >=0.60 U/h (=0.01 U/min)
#     - no non-study pressor exposure during prior 2 h and none active at t0
#
# This is the G9 scenario:
#   dovss_stable120_no_otherpressor2h
#
# Primary strategies:
#   NE-first: first observed cumulative fall from t0 NE dose >=0.10 mcg/kg/min
#             (or complete NE stop) within 120 min.
#   AVP-first: first observed cumulative fall from t0 AVP dose >=0.60 U/h
#              (or complete AVP stop) within 120 min.
#   No future persistence requirement is used. Only exactly simultaneous
#   threshold crossings are unresolved.
#
# Primary endpoint:
#   first episode of two MAP values <65 mmHg separated by 15-60 min,
#   within 6 h of t0.
#
# Primary estimator:
#   marginal/time-only stabilized inverse-probability-of-artificial-censoring
#   weighting + weighted Kaplan-Meier.
#
# Primary weights:
#   cumulative stabilized weights without percentile truncation, after the
#   pre-specified per-interval probability clip [0.01, 0.999].
#
# Sensitivities:
#   - same t0, same DOVSS-scale strategy, grace 240 min
#   - same t0, immediate cumulative >=30% reduction, grace 120 min
#   - 1st/99th percentile truncation of cumulative stabilized weights
#   - baseline-stabilized weighted outcome-model standardization (point estimate)
#
# Uncertainty:
#   subject-level cluster bootstrap. A sampled subject contributes ALL eligible
#   ICU stays and both clones; censor models are refit in every replicate.
#   Default bootstrap replicates = 500, configurable with environment variable
#   G10_N_BOOT (e.g. G10_N_BOOT=1000).
#
# Primary inference is for the 6-h endpoint. 12-h and 24-h estimates are reported
# as secondary/exploratory because natural ICU discharge/death censoring is more
# consequential at longer horizons and is not additionally IPC-weighted here.
#
# This script is intended to be the first inferential frozen-protocol analysis.
# It still reports extensive diagnostics and will not export patient-level data.

from __future__ import annotations

import csv
import getpass
import json
import math
import os
import statistics
import time
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
DBNAME = "mimiciv"
PORT = 5432
SOCKET_DIR = "/var/run/postgresql"
USER = getpass.getuser()
TIMEOUT_MS = 300000

OUT = (
    Path.cwd()
    / "MIMIC_NE_AVP_CCW"
    / "G10_final_CCW_bootstrap"
    / datetime.now().strftime("%Y%m%d_%H%M%S")
)

NE_ITEMID = 221906
AVP_ITEMID = 222315

STEP_MIN = 15
SEARCH_MAX_DAYS = 7
DUAL_LOOKBACK_MIN = 120

MAP_LOOKBACK_MIN = 120
MAP_MIN_READINGS = 3
MAP_MIN_SPAN_MIN = 90
MAP_MAX_AGE_MIN = 30
MAP_THRESHOLD = 65.0

NO_ESCALATION_LOOKBACK_MIN = 120
LEGACY_TAPER_REL = 0.10
LEGACY_NE_TAPER_ABS = 0.01
LEGACY_AVP_TAPER_ABS = 0.12
LEGACY_TAPER_PERSIST_MIN = 15

PRIMARY_NE_MAX = 0.30
PRIMARY_AVP_MAX = 3.60

ACTION_CONFIRM_MIN = 30
SIMULTANEOUS_MIN = 0

DOVSS_NE_ABS_DROP = 0.10
DOVSS_AVP_ABS_DROP = 0.60

PRIMARY_GRACE_MIN = 120
OUTCOME_WINDOWS_H = (6, 12, 24)

SPECS = [
    dict(
        name="primary_dovss120_immediate",
        mode="dovss_scale_immediate",
        grace_minutes=120,
    ),
    dict(
        name="sensitivity_dovss240_immediate",
        mode="dovss_scale_immediate",
        grace_minutes=240,
    ),
    dict(
        name="sensitivity_cum30_120_immediate",
        mode="cumulative30_immediate",
        grace_minutes=120,
    ),
]

NE_HARD_MAX = 10.0
AVP_HARD_MAX = 30.0

# Require exact G5 replication before CCW proceeds.
EXPECTED_G5_PRIMARY_N = 2744
EXPECTED_FINAL_ELIGIBLE_N = 2260
EXPECTED_FINAL_PRIMARY_ACTIONS = {
    "NE_first": 204,
    "AVP_first": 139,
    "near_simultaneous": 10,
    "no_clear_action": 1907,
}
BOOTSTRAP_REPS = int(os.environ.get("G10_N_BOOT", "500"))
BOOTSTRAP_SEED = int(os.environ.get("G10_BOOT_SEED", "20261015"))
BOOTSTRAP_MIN_SUCCESS_FRACTION = 0.90
EXPECTED_G6_CCW_N = 2094
EXPECTED_G6_PRIMARY_ACTIONS = {
    "NE_first": 175,
    "AVP_first": 130,
    "near_simultaneous": 21,
    "no_clear_action": 1768,
}

# Logistic model safeguards.
RIDGE_LAMBDA = 1e-3
LOGIT_MAX_ITER = 100
PROB_CLIP_LO = 0.01
PROB_CLIP_HI = 0.999

# Primary stabilized-weight truncation.
WEIGHT_LO_Q = 0.01
WEIGHT_HI_Q = 0.99

LOG: List[str] = []
ERRORS: List[dict] = []
WARNINGS: List[str] = []


def log(*parts) -> None:
    s = " ".join(str(x) for x in parts)
    print(s)
    LOG.append(s)


def warn(*parts) -> None:
    s = " ".join(str(x) for x in parts)
    print("WARNING:", s)
    LOG.append("WARNING: " + s)
    WARNINGS.append(s)


def save_csv(
    name: str,
    rows: Iterable[dict],
    columns: Optional[Sequence[str]] = None,
) -> Path:
    rows = list(rows)
    path = OUT / name
    if columns is None:
        columns = list(rows[0].keys()) if rows else ["message"]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=list(columns),
            extrasaction="ignore",
        )
        w.writeheader()
        if rows:
            w.writerows(rows)
    return path


def median_or_none(xs: Iterable[Optional[float]]) -> Optional[float]:
    vals = []
    for x in xs:
        if x is None:
            continue
        x = float(x)
        if math.isfinite(x):
            vals.append(x)
    return statistics.median(vals) if vals else None


def qtile_or_none(xs: Iterable[Optional[float]], q: float) -> Optional[float]:
    vals = sorted(
        float(x)
        for x in xs
        if x is not None and math.isfinite(float(x))
    )
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def ceil_dt(dt: datetime, minutes: int) -> datetime:
    epoch = datetime(1970, 1, 1, tzinfo=dt.tzinfo)
    total = (dt - epoch).total_seconds()
    step = minutes * 60
    return epoch + timedelta(seconds=math.ceil(total / step) * step)


def min_time(*times: Optional[datetime]) -> Optional[datetime]:
    vals = [x for x in times if x is not None]
    return min(vals) if vals else None


# ---------------------------------------------------------------------
# Timeline reconstruction
# ---------------------------------------------------------------------
@dataclass(frozen=True)
class Segment:
    start: datetime
    end: datetime
    drug: str
    rate: float


@dataclass(frozen=True)
class Atomic:
    start: datetime
    end: datetime
    ne: Optional[float]
    avp: Optional[float]
    ne_conflict: bool = False
    avp_conflict: bool = False


@dataclass(frozen=True)
class ActionEvent:
    time: datetime
    drug: str
    kind: str
    from_rate: float
    to_rate: float


class Timeline:
    def __init__(self, atoms: Sequence[Atomic]):
        self.atoms = list(atoms)
        self.starts = [x.start for x in self.atoms]

    def atom_at(self, t: datetime) -> Optional[Atomic]:
        i = bisect_right(self.starts, t) - 1
        if i < 0:
            return None
        a = self.atoms[i]
        if a.start <= t < a.end:
            return a
        return None

    def rate_at(self, t: datetime, drug: str) -> Optional[float]:
        a = self.atom_at(t)
        if a is None:
            return 0.0
        if drug == "NE":
            return None if a.ne_conflict else float(a.ne or 0.0)
        return None if a.avp_conflict else float(a.avp or 0.0)

    def first_dual_time(self) -> Optional[datetime]:
        for a in self.atoms:
            if (
                not a.ne_conflict
                and not a.avp_conflict
                and a.ne is not None
                and a.avp is not None
                and a.ne > 0
                and a.avp > 0
            ):
                return a.start
        return None

    def continuous_dual(self, start: datetime, end: datetime) -> bool:
        if end <= start:
            return False

        cursor = start
        for a in self.atoms:
            if a.end <= start:
                continue
            if a.start >= end:
                break

            left = max(a.start, start)
            right = min(a.end, end)

            if right <= left:
                continue

            if left > cursor + timedelta(seconds=1):
                return False

            if (
                a.ne_conflict
                or a.avp_conflict
                or a.ne is None
                or a.avp is None
                or a.ne <= 0
                or a.avp <= 0
            ):
                return False

            cursor = max(cursor, right)
            if cursor >= end:
                return True

        return cursor >= end

    def transitions(self, drug: str) -> List[Tuple[datetime, float, float, float]]:
        out = []
        prev_rate = 0.0
        prev_conflict = False
        prev_end = None

        for a in self.atoms:
            conflict = a.ne_conflict if drug == "NE" else a.avp_conflict
            raw = a.ne if drug == "NE" else a.avp
            rate = 0.0 if raw is None else float(raw)

            if prev_end is not None and a.start > prev_end:
                prev_rate = 0.0
                prev_conflict = False

            if (
                not conflict
                and not prev_conflict
                and abs(rate - prev_rate) > 1e-12
            ):
                persist = max(
                    0.0,
                    (a.end - a.start).total_seconds() / 60.0,
                )
                out.append((a.start, prev_rate, rate, persist))

            prev_rate = rate
            prev_conflict = conflict
            prev_end = a.end

        return out


def normalize_rate(
    itemid,
    rate,
    rateuom,
    patientweight,
) -> Tuple[Optional[str], Optional[float], str]:
    if rate is None:
        return None, None, "missing_rate"

    try:
        r = float(rate)
    except Exception:
        return None, None, "bad_rate"

    if not math.isfinite(r) or r <= 0:
        return None, None, "nonpositive_rate"

    u = (rateuom or "").strip().lower()

    if int(itemid) == NE_ITEMID:
        drug = "NE"

        if u == "mcg/kg/min":
            rr = r
        elif u == "mg/kg/min":
            if patientweight is not None and abs(float(patientweight) - 1.0) < 1e-9:
                rr = r
            else:
                rr = r * 1000.0
        else:
            return drug, None, "unexpected_unit"

        if rr > NE_HARD_MAX:
            return drug, None, "hard_implausible"

        return drug, rr, "accepted"

    if int(itemid) == AVP_ITEMID:
        drug = "AVP"

        if u == "units/hour":
            rr = r
        elif u == "units/min":
            rr = r * 60.0
        else:
            return drug, None, "unexpected_unit"

        if rr > AVP_HARD_MAX:
            return drug, None, "hard_implausible"

        return drug, rr, "accepted"

    return None, None, "wrong_itemid"


def segments_from_rows(rows: Sequence[dict]) -> Tuple[List[Segment], Counter]:
    segs: List[Segment] = []
    qc = Counter()

    for r in rows:
        status = str(r.get("statusdescription") or "").strip().lower()
        ordercat = str(r.get("ordercategorydescription") or "").strip().lower()

        if status in {"rewritten", "flushed"}:
            qc["rewritten_or_flushed"] += 1
            continue

        if int(r["itemid"]) == NE_ITEMID and (
            "bolus" in status or "bolus" in ordercat
        ):
            qc["ne_bolus"] += 1
            continue

        st = r.get("starttime")
        en = r.get("endtime")

        if st is None or en is None or en <= st:
            qc["invalid_interval"] += 1
            continue

        drug, rr, reason = normalize_rate(
            r["itemid"],
            r.get("rate"),
            r.get("rateuom"),
            r.get("patientweight"),
        )

        qc[f"{drug or 'UNK'}_{reason}"] += 1

        if rr is None:
            continue

        segs.append(
            Segment(st, en, drug, rr)
        )

    return segs, qc


def build_atomic_timeline(
    segments: Sequence[Segment],
) -> Tuple[Timeline, Counter]:
    qc = Counter()

    if not segments:
        return Timeline([]), qc

    bps = sorted(
        {s.start for s in segments}
        | {s.end for s in segments}
    )

    atoms: List[Atomic] = []

    for left, right in zip(bps[:-1], bps[1:]):
        if right <= left:
            continue

        active = {"NE": [], "AVP": []}

        for s in segments:
            if s.start <= left and s.end >= right:
                active[s.drug].append(float(s.rate))

        values = {}
        conflicts = {}

        for drug in ("NE", "AVP"):
            rs = active[drug]

            if not rs:
                values[drug] = 0.0
                conflicts[drug] = False
                continue

            lo = min(rs)
            hi = max(rs)
            tol = max(
                1e-9,
                0.01 * max(abs(lo), abs(hi)),
            )
            conflict = (hi - lo) > tol

            conflicts[drug] = conflict

            if conflict:
                qc[f"{drug}_conflict_atomic"] += 1
                values[drug] = None
            else:
                values[drug] = statistics.median(rs)

        cur = Atomic(
            left,
            right,
            values["NE"],
            values["AVP"],
            conflicts["NE"],
            conflicts["AVP"],
        )

        if atoms:
            p = atoms[-1]
            same = (
                p.ne == cur.ne
                and p.avp == cur.avp
                and p.ne_conflict == cur.ne_conflict
                and p.avp_conflict == cur.avp_conflict
                and p.end == cur.start
            )

            if same:
                atoms[-1] = Atomic(
                    p.start,
                    cur.end,
                    p.ne,
                    p.avp,
                    p.ne_conflict,
                    p.avp_conflict,
                )
                continue

        atoms.append(cur)

    return Timeline(atoms), qc


# ---------------------------------------------------------------------
# Fixed primary t0 from G5
# ---------------------------------------------------------------------
def legacy_meaningful_reduction(
    from_rate: float,
    to_rate: float,
    drug: str,
) -> bool:
    if from_rate <= 0 or to_rate >= from_rate:
        return False

    abs_thr = (
        LEGACY_NE_TAPER_ABS
        if drug == "NE"
        else LEGACY_AVP_TAPER_ABS
    )

    delta = from_rate - to_rate
    eps = 1e-12

    return (
        delta + eps >= abs_thr
        and delta / from_rate + eps >= LEGACY_TAPER_REL
    )


def has_recent_escalation(
    tl: Timeline,
    t: datetime,
    drug: str,
) -> bool:
    start = t - timedelta(
        minutes=NO_ESCALATION_LOOKBACK_MIN
    )
    abs_thr = (
        LEGACY_NE_TAPER_ABS
        if drug == "NE"
        else LEGACY_AVP_TAPER_ABS
    )

    for tt, fr, to, persist in tl.transitions(drug):
        if tt <= start or tt > t:
            continue

        if (
            fr <= 0
            or to <= fr
            or persist < LEGACY_TAPER_PERSIST_MIN
        ):
            continue

        delta = to - fr

        if (
            delta + 1e-12 >= abs_thr
            and delta / fr + 1e-12 >= LEGACY_TAPER_REL
        ):
            return True

    return False


def map_stable(
    map_rows: Sequence[Tuple[datetime, float]],
    t: datetime,
) -> bool:
    lo = t - timedelta(minutes=MAP_LOOKBACK_MIN)

    vals = [
        (tt, float(v))
        for tt, v in map_rows
        if lo <= tt <= t and v is not None
    ]

    if len(vals) < MAP_MIN_READINGS:
        return False

    times = [x[0] for x in vals]

    if (
        max(times) - min(times)
    ).total_seconds() < MAP_MIN_SPAN_MIN * 60:
        return False

    if (
        t - max(times)
    ).total_seconds() > MAP_MAX_AGE_MIN * 60:
        return False

    return all(v >= MAP_THRESHOLD for _, v in vals)


def find_primary_t0(
    tl: Timeline,
    map_rows: Sequence[Tuple[datetime, float]],
    meta: dict,
) -> Tuple[Optional[datetime], dict]:
    first_dual = tl.first_dual_time()

    stage = dict(
        two_hour_dual=False,
        sepsis_before_t=False,
        adult_ok=False,
        map_stable=False,
        dose_ok=False,
        no_escalation=False,
        eligible=False,
    )

    if first_dual is None:
        return None, stage

    start = ceil_dt(
        first_dual + timedelta(minutes=DUAL_LOOKBACK_MIN),
        STEP_MIN,
    )
    stop = min(
        meta["outtime"],
        first_dual + timedelta(days=SEARCH_MAX_DAYS),
    )

    t = start

    while t <= stop:
        if not tl.continuous_dual(
            t - timedelta(minutes=DUAL_LOOKBACK_MIN),
            t,
        ):
            t += timedelta(minutes=STEP_MIN)
            continue

        stage["two_hour_dual"] = True

        sus = meta.get("suspected_infection_time")
        sofa = meta.get("sofa_time")

        if (
            sus is None
            or sofa is None
            or sus > t
            or sofa > t
        ):
            t += timedelta(minutes=STEP_MIN)
            continue

        stage["sepsis_before_t"] = True

        age = meta.get("age")
        adult_ok = (
            age is not None
            and float(age) >= 18.0
        )

        if not adult_ok:
            t += timedelta(minutes=STEP_MIN)
            continue

        stage["adult_ok"] = True

        if not map_stable(map_rows, t):
            t += timedelta(minutes=STEP_MIN)
            continue

        stage["map_stable"] = True

        ne = tl.rate_at(t, "NE")
        avp = tl.rate_at(t, "AVP")

        if (
            ne is None
            or avp is None
            or ne <= 0
            or avp <= 0
            or ne > PRIMARY_NE_MAX + 1e-12
            or avp > PRIMARY_AVP_MAX + 1e-12
        ):
            t += timedelta(minutes=STEP_MIN)
            continue

        stage["dose_ok"] = True

        noesc = not (
            has_recent_escalation(tl, t, "NE")
            or has_recent_escalation(tl, t, "AVP")
        )

        if not noesc:
            t += timedelta(minutes=STEP_MIN)
            continue

        stage["no_escalation"] = True
        stage["eligible"] = True

        return t, stage

    return None, stage


# ---------------------------------------------------------------------
# Clinically meaningful action definitions
# ---------------------------------------------------------------------
def _drug_conflict(a: Atomic, drug: str) -> bool:
    return (
        a.ne_conflict
        if drug == "NE"
        else a.avp_conflict
    )


def _drug_value(a: Atomic, drug: str) -> float:
    raw = a.ne if drug == "NE" else a.avp
    return 0.0 if raw is None else float(raw)


def stays_at_or_below(
    tl: Timeline,
    drug: str,
    start: datetime,
    minutes: int,
    threshold: float,
    observation_end: datetime,
) -> bool:
    end = start + timedelta(minutes=minutes)

    if end > observation_end:
        return False

    for a in tl.atoms:
        if a.end <= start:
            continue
        if a.start >= end:
            break

        if _drug_conflict(a, drug):
            return False

        if _drug_value(a, drug) > threshold + 1e-12:
            return False

    return True


def cumulative_event_for_drug(
    tl: Timeline,
    t0: datetime,
    grace_min: int,
    drug: str,
    baseline: float,
    mode: str,
    observation_end: datetime,
) -> Optional[ActionEvent]:
    if baseline <= 0:
        return None

    if mode == "cumulative30_30m":
        threshold = baseline * 0.70
    elif mode == "dovss_scale_30m":
        abs_drop = (
            DOVSS_NE_ABS_DROP
            if drug == "NE"
            else DOVSS_AVP_ABS_DROP
        )
        threshold = max(
            0.0,
            baseline - abs_drop,
        )
    else:
        raise ValueError(
            f"Unknown action mode: {mode}"
        )

    grace_end = t0 + timedelta(
        minutes=grace_min
    )

    for tt, fr, to, _persist in tl.transitions(drug):
        if (
            tt <= t0
            or tt > grace_end
            or fr <= 0
        ):
            continue

        qualifies = (
            to <= 0
            or to <= threshold + 1e-12
        )

        if not qualifies:
            continue

        if not stays_at_or_below(
            tl,
            drug,
            tt,
            ACTION_CONFIRM_MIN,
            threshold,
            observation_end,
        ):
            continue

        return ActionEvent(
            time=tt,
            drug=drug,
            kind=("stop" if to <= 0 else "taper"),
            from_rate=baseline,
            to_rate=to,
        )

    return None


def get_spec_events(
    tl: Timeline,
    rec: dict,
    spec: dict,
) -> Tuple[Optional[ActionEvent], Optional[ActionEvent]]:
    ne = cumulative_event_for_drug(
        tl,
        rec["t0"],
        spec["grace_minutes"],
        "NE",
        rec["ne0"],
        spec["mode"],
        rec["outtime"],
    )

    avp = cumulative_event_for_drug(
        tl,
        rec["t0"],
        spec["grace_minutes"],
        "AVP",
        rec["avp0"],
        spec["mode"],
        rec["outtime"],
    )

    return ne, avp


def classify_events(
    ne: Optional[ActionEvent],
    avp: Optional[ActionEvent],
) -> str:
    if ne is None and avp is None:
        return "no_clear_action"

    if ne is not None and avp is not None:
        diff = abs(
            (ne.time - avp.time).total_seconds()
        ) / 60.0

        if diff <= SIMULTANEOUS_MIN:
            return "near_simultaneous"

        return (
            "NE_first"
            if ne.time < avp.time
            else "AVP_first"
        )

    return (
        "NE_first"
        if ne is not None
        else "AVP_first"
    )


def clone_protocol_times(
    rec: dict,
    tl: Timeline,
    spec: dict,
    arm: str,
) -> dict:
    """
    Returns protocol completion and artificial censoring times.

    Once the assigned drug is the first qualifying reduction, the strategy is
    completed and no further artificial censoring is imposed.

    If the opposite drug qualifies first, the clone is censored at that action.
    Near-simultaneous qualifying reductions are treated as protocol deviation
    for BOTH arms because the intervention "which drug first" is not resolved.
    No qualifying action by grace end -> artificial censor at grace end.
    """
    ne, avp = get_spec_events(
        tl,
        rec,
        spec,
    )

    observed = classify_events(ne, avp)
    grace_end = rec["t0"] + timedelta(
        minutes=spec["grace_minutes"]
    )

    if observed == "no_clear_action":
        return dict(
            observed_strategy=observed,
            completion_time=None,
            censor_time=grace_end,
            censor_reason="no_action_by_grace",
        )

    if observed == "near_simultaneous":
        first_time = min(
            ne.time,
            avp.time,
        )
        return dict(
            observed_strategy=observed,
            completion_time=None,
            censor_time=first_time,
            censor_reason="near_simultaneous",
        )

    first_event = (
        ne if observed == "NE_first" else avp
    )

    if arm == observed:
        return dict(
            observed_strategy=observed,
            completion_time=first_event.time,
            censor_time=None,
            censor_reason=None,
        )

    return dict(
        observed_strategy=observed,
        completion_time=None,
        censor_time=first_event.time,
        censor_reason="opposite_action_first",
    )


# ---------------------------------------------------------------------
# Measurement summaries
# ---------------------------------------------------------------------
def window_values(
    rows: Sequence[Tuple[datetime, float]],
    t: datetime,
    lookback_min: int,
) -> List[Tuple[datetime, float]]:
    lo = t - timedelta(minutes=lookback_min)
    return [
        (tt, float(v))
        for tt, v in rows
        if lo <= tt < t and v is not None
    ]


def mean_or_none(vals: Sequence[float]) -> Optional[float]:
    return (
        float(statistics.mean(vals))
        if vals
        else None
    )


def min_or_none(vals: Sequence[float]) -> Optional[float]:
    return (
        float(min(vals))
        if vals
        else None
    )


def latest_before(
    rows: Sequence[Tuple[datetime, float]],
    t: datetime,
    lookback_min: int,
) -> Optional[float]:
    vals = window_values(
        rows,
        t,
        lookback_min,
    )
    if not vals:
        return None
    vals.sort(key=lambda x: x[0])
    return float(vals[-1][1])


def recent_summary(
    rows: Sequence[Tuple[datetime, float]],
    t: datetime,
    lookback_min: int,
) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    vals = window_values(
        rows,
        t,
        lookback_min,
    )
    ys = [v for _, v in vals]

    return (
        mean_or_none(ys),
        min_or_none(ys),
        (
            sorted(vals, key=lambda x: x[0])[-1][1]
            if vals
            else None
        ),
    )


def latest_lactate(
    rows: Sequence[Tuple[datetime, float]],
    t: datetime,
    lookback_h: int = 24,
) -> Optional[float]:
    lo = t - timedelta(hours=lookback_h)

    vals = [
        (tt, float(v))
        for tt, v in rows
        if lo <= tt < t and v is not None
    ]

    if not vals:
        return None

    vals.sort(key=lambda x: x[0])
    return float(vals[-1][1])


def first_pair_lt65_time(
    map_rows: Sequence[Tuple[datetime, float]],
    t0: datetime,
    horizon_h: int,
) -> Optional[datetime]:
    end = t0 + timedelta(hours=horizon_h)

    lows = sorted(
        [
            (tt, float(v))
            for tt, v in map_rows
            if (
                t0 < tt <= end
                and v is not None
                and float(v) < 65.0
            )
        ],
        key=lambda x: x[0],
    )

    for j in range(len(lows)):
        t2 = lows[j][0]

        for i in range(j - 1, -1, -1):
            t1 = lows[i][0]
            delta = (
                t2 - t1
            ).total_seconds() / 60.0

            if delta > 60.0:
                break

            if 15.0 <= delta <= 60.0:
                return t2

    return None


def first_lt60_time(
    map_rows: Sequence[Tuple[datetime, float]],
    t0: datetime,
    horizon_h: int,
) -> Optional[datetime]:
    end = t0 + timedelta(hours=horizon_h)

    vals = sorted(
        [
            tt
            for tt, v in map_rows
            if (
                t0 < tt <= end
                and v is not None
                and float(v) < 60.0
            )
        ]
    )

    return vals[0] if vals else None


# ---------------------------------------------------------------------
# Other vasoactive / steroid exposure
# ---------------------------------------------------------------------
def classify_vasoactive_label(label: str) -> str:
    s = (label or "").lower()

    if "norepinephrine" in s or "vasopressin" in s:
        return "study_drug"

    pressor_terms = [
        "epinephrine",
        "phenylephrine",
        "dopamine",
        "angiotensin",
        "methylene blue",
        "isoproterenol",
        "isoprenaline",
    ]

    inotrope_terms = [
        "dobutamine",
        "milrinone",
    ]

    if any(term in s for term in pressor_terms):
        return "other_pressor"

    if any(term in s for term in inotrope_terms):
        return "inotrope"

    return "other"


def interval_active(
    rows: Sequence[dict],
    t: datetime,
    classes: Sequence[str],
    label_classes: Dict[int, str],
) -> bool:
    classes = set(classes)

    for r in rows:
        itemid = int(r["itemid"])

        if label_classes.get(
            itemid,
            "other",
        ) not in classes:
            continue

        st = r.get("starttime")
        en = r.get("endtime")

        if (
            st is not None
            and en is not None
            and st <= t < en
        ):
            return True

    return False


def interval_overlap(
    rows: Sequence[dict],
    start: datetime,
    end: datetime,
    classes: Sequence[str],
    label_classes: Dict[int, str],
) -> bool:
    classes = set(classes)

    for r in rows:
        itemid = int(r["itemid"])

        if label_classes.get(
            itemid,
            "other",
        ) not in classes:
            continue

        st = r.get("starttime")
        en = r.get("endtime")

        if (
            st is not None
            and en is not None
            and st < end
            and en > start
        ):
            return True

    return False


def steroid_active(
    rows: Sequence[dict],
    t: datetime,
) -> bool:
    for r in rows:
        st = r.get("starttime")
        en = r.get("stoptime")

        if st is None:
            continue

        if st <= t and (
            en is None or t <= en
        ):
            return True

    return False


def steroid_prior24h(
    rows: Sequence[dict],
    t: datetime,
) -> bool:
    lo = t - timedelta(hours=24)

    for r in rows:
        st = r.get("starttime")
        en = r.get("stoptime")

        if st is None:
            continue

        effective_end = en or t

        if (
            st < t
            and effective_end > lo
        ):
            return True

    return False


# ---------------------------------------------------------------------
# Categorical collapsing
# ---------------------------------------------------------------------
def race_group(x) -> str:
    s = str(x or "").upper()

    if "WHITE" in s:
        return "White"
    if "BLACK" in s:
        return "Black"
    if "ASIAN" in s:
        return "Asian"
    if "HISPANIC" in s or "LATINO" in s:
        return "Hispanic"
    if not s or s == "NAN":
        return "Missing"
    return "Other"


def collapse_rare(
    values: Sequence[str],
    min_count: int = 30,
) -> Dict[str, str]:
    c = Counter(values)

    return {
        v: (
            v
            if c[v] >= min_count
            else "Other"
        )
        for v in c
    }


# ---------------------------------------------------------------------
# Ridge logistic regression (NumPy only)
# ---------------------------------------------------------------------
@dataclass
class DesignSpec:
    columns: List[str]
    continuous_columns: List[str]
    categorical_columns: List[str]
    medians: Dict[str, float]
    means: Dict[str, float]
    sds: Dict[str, float]
    dummy_columns: List[str]


@dataclass
class LogisticModel:
    beta: np.ndarray
    design: DesignSpec
    ridge_lambda: float
    converged: bool
    iterations: int


def _prepare_design_fit(
    df: pd.DataFrame,
    continuous: Sequence[str],
    categorical: Sequence[str],
) -> Tuple[np.ndarray, DesignSpec]:
    work = pd.DataFrame(index=df.index)

    medians = {}
    means = {}
    sds = {}

    for col in continuous:
        x = pd.to_numeric(
            df[col],
            errors="coerce",
        )

        med = float(
            x.median()
            if x.notna().any()
            else 0.0
        )

        miss = x.isna().astype(float)
        x = x.fillna(med)

        mu = float(x.mean())
        sd = float(x.std(ddof=0))

        if not math.isfinite(sd) or sd < 1e-8:
            sd = 1.0

        work[col] = (
            x - mu
        ) / sd

        work[
            f"{col}__missing"
        ] = miss

        medians[col] = med
        means[col] = mu
        sds[col] = sd

    cat_df = pd.DataFrame(index=df.index)

    for col in categorical:
        cat_df[col] = (
            df[col]
            .astype("string")
            .fillna("Missing")
        )

    if categorical:
        dummies = pd.get_dummies(
            cat_df,
            columns=list(categorical),
            drop_first=False,
            dtype=float,
        )
    else:
        dummies = pd.DataFrame(
            index=df.index
        )

    dummy_columns = list(
        dummies.columns
    )

    full = pd.concat(
        [work, dummies],
        axis=1,
    )

    full.insert(
        0,
        "intercept",
        1.0,
    )

    columns = list(full.columns)

    spec = DesignSpec(
        columns=columns,
        continuous_columns=list(
            continuous
        ),
        categorical_columns=list(
            categorical
        ),
        medians=medians,
        means=means,
        sds=sds,
        dummy_columns=dummy_columns,
    )

    return full.to_numpy(dtype=float), spec


def _prepare_design_apply(
    df: pd.DataFrame,
    spec: DesignSpec,
) -> np.ndarray:
    work = pd.DataFrame(index=df.index)

    for col in spec.continuous_columns:
        x = pd.to_numeric(
            df[col],
            errors="coerce",
        )

        miss = x.isna().astype(float)

        x = x.fillna(
            spec.medians[col]
        )

        work[col] = (
            x - spec.means[col]
        ) / spec.sds[col]

        work[
            f"{col}__missing"
        ] = miss

    cat_df = pd.DataFrame(index=df.index)

    for col in spec.categorical_columns:
        cat_df[col] = (
            df[col]
            .astype("string")
            .fillna("Missing")
        )

    if spec.categorical_columns:
        dummies = pd.get_dummies(
            cat_df,
            columns=spec.categorical_columns,
            drop_first=False,
            dtype=float,
        )
    else:
        dummies = pd.DataFrame(
            index=df.index
        )

    for col in spec.dummy_columns:
        if col not in dummies:
            dummies[col] = 0.0

    dummies = dummies[
        spec.dummy_columns
    ]

    full = pd.concat(
        [work, dummies],
        axis=1,
    )

    full.insert(
        0,
        "intercept",
        1.0,
    )

    for col in spec.columns:
        if col not in full:
            full[col] = 0.0

    full = full[
        spec.columns
    ]

    return full.to_numpy(dtype=float)


def sigmoid(z: np.ndarray) -> np.ndarray:
    z = np.clip(
        z,
        -35.0,
        35.0,
    )
    return 1.0 / (
        1.0 + np.exp(-z)
    )


def fit_logistic_ridge(
    df: pd.DataFrame,
    y_col: str,
    continuous: Sequence[str],
    categorical: Sequence[str],
    ridge_lambda: float = RIDGE_LAMBDA,
) -> LogisticModel:
    X, design = _prepare_design_fit(
        df,
        continuous,
        categorical,
    )

    y = (
        pd.to_numeric(
            df[y_col],
            errors="coerce",
        )
        .fillna(0)
        .to_numpy(dtype=float)
    )

    beta = np.zeros(
        X.shape[1],
        dtype=float,
    )

    penalty = np.eye(
        X.shape[1],
        dtype=float,
    )
    penalty[0, 0] = 0.0

    converged = False

    for iteration in range(
        1,
        LOGIT_MAX_ITER + 1,
    ):
        eta = X @ beta
        p = sigmoid(eta)

        w = np.clip(
            p * (1.0 - p),
            1e-6,
            None,
        )

        grad = (
            X.T @ (y - p)
            - ridge_lambda
            * penalty
            @ beta
        )

        hess_pos = (
            X.T
            @ (X * w[:, None])
            + ridge_lambda
            * penalty
            + 1e-8
            * np.eye(X.shape[1])
        )

        try:
            step = np.linalg.solve(
                hess_pos,
                grad,
            )
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(
                hess_pos
            ) @ grad

        beta_new = beta + step

        if np.max(
            np.abs(beta_new - beta)
        ) < 1e-7:
            beta = beta_new
            converged = True
            break

        beta = beta_new

    return LogisticModel(
        beta=beta,
        design=design,
        ridge_lambda=ridge_lambda,
        converged=converged,
        iterations=iteration,
    )


def predict_logistic(
    model: LogisticModel,
    df: pd.DataFrame,
) -> np.ndarray:
    X = _prepare_design_apply(
        df,
        model.design,
    )
    return sigmoid(
        X @ model.beta
    )


# ---------------------------------------------------------------------
# Weight-model / clone-row construction
# ---------------------------------------------------------------------
BASELINE_CONT = [
    "age",
    "sofa_score",
    "ne0",
    "avp0",
    "delay_h",
    "baseline_lactate",
    "baseline_map_mean",
    "baseline_map_min",
    "baseline_hr_mean",
    "steroid_prior24h",
    "inotrope_at_t0",
]

BASELINE_CAT = [
    "gender",
    "race_group",
    "careunit_group",
]

NUMERATOR_CONT = BASELINE_CONT + [
    "elapsed_h",
    "elapsed_h2",
]

NUMERATOR_CAT = BASELINE_CAT

DENOMINATOR_CONT = NUMERATOR_CONT + [
    "current_ne",
    "current_avp",
    "recent_map_mean",
    "recent_map_min",
    "recent_map_latest",
    "recent_hr_mean",
    "recent_lactate",
    "other_pressor_active",
    "inotrope_active",
    "steroid_active",
]

DENOMINATOR_CAT = BASELINE_CAT


def interval_features(
    rec: dict,
    tl: Timeline,
    t: datetime,
    by_stay_map: Dict[int, list],
    by_stay_hr: Dict[int, list],
    by_hadm_lac: Dict[int, list],
    by_stay_other: Dict[int, list],
    label_classes: Dict[int, str],
    by_hadm_steroid: Dict[int, list],
) -> dict:
    sid = rec["sid"]
    hadm_id = rec["hadm_id"]

    # Use values just BEFORE interval start to avoid conditioning on the
    # treatment action occurring at the interval boundary.
    teval = t - timedelta(
        microseconds=1
    )

    current_ne = tl.rate_at(
        teval,
        "NE",
    )
    current_avp = tl.rate_at(
        teval,
        "AVP",
    )

    map_mean, map_minimum, map_latest = recent_summary(
        by_stay_map.get(
            sid,
            [],
        ),
        teval,
        60,
    )

    hr_mean, _hr_min, _hr_latest = recent_summary(
        by_stay_hr.get(
            sid,
            [],
        ),
        teval,
        60,
    )

    lac = latest_lactate(
        by_hadm_lac.get(
            hadm_id,
            [],
        ),
        teval,
        24,
    )

    other_rows = by_stay_other.get(
        sid,
        [],
    )

    steroid_rows = by_hadm_steroid.get(
        hadm_id,
        [],
    )

    elapsed_h = (
        t - rec["t0"]
    ).total_seconds() / 3600.0

    return dict(
        elapsed_h=elapsed_h,
        elapsed_h2=elapsed_h ** 2,
        current_ne=current_ne,
        current_avp=current_avp,
        recent_map_mean=map_mean,
        recent_map_min=map_minimum,
        recent_map_latest=map_latest,
        recent_hr_mean=hr_mean,
        recent_lactate=lac,
        other_pressor_active=int(
            interval_active(
                other_rows,
                teval,
                ["other_pressor"],
                label_classes,
            )
        ),
        inotrope_active=int(
            interval_active(
                other_rows,
                teval,
                ["inotrope"],
                label_classes,
            )
        ),
        steroid_active=int(
            steroid_active(
                steroid_rows,
                teval,
            )
        ),
    )


def build_weight_rows_for_spec(
    records: Sequence[dict],
    timelines: Dict[int, Timeline],
    spec: dict,
    by_stay_map: Dict[int, list],
    by_stay_hr: Dict[int, list],
    by_hadm_lac: Dict[int, list],
    by_stay_other: Dict[int, list],
    label_classes: Dict[int, str],
    by_hadm_steroid: Dict[int, list],
    outcome_h: int,
) -> Tuple[pd.DataFrame, List[dict]]:
    rows = []
    clone_protocol = []

    for rec in records:
        sid = rec["sid"]
        tl = timelines[sid]

        outcome_time = first_pair_lt65_time(
            by_stay_map.get(
                sid,
                [],
            ),
            rec["t0"],
            outcome_h,
        )

        for arm in (
            "NE_first",
            "AVP_first",
        ):
            proto = clone_protocol_times(
                rec,
                tl,
                spec,
                arm,
            )

            clone_id = (
                f"{sid}|{spec['name']}|{arm}"
            )

            proto_row = dict(
                clone_id=clone_id,
                sid=sid,
                subject_id=rec["subject_id"],
                arm=arm,
                spec_name=spec["name"],
                observed_strategy=proto[
                    "observed_strategy"
                ],
                completion_time=proto[
                    "completion_time"
                ],
                censor_time=proto[
                    "censor_time"
                ],
                censor_reason=proto[
                    "censor_reason"
                ],
                outcome_time=outcome_time,
                natural_end=rec["outtime"],
            )
            clone_protocol.append(
                proto_row
            )

            grace_end = rec[
                "t0"
            ] + timedelta(
                minutes=spec[
                    "grace_minutes"
                ]
            )

            # Weight model is needed only until outcome, strategy completion,
            # artificial censoring, natural censoring, or grace end.
            stop_for_weight = min_time(
                outcome_time,
                proto["completion_time"],
                proto["censor_time"],
                rec["outtime"],
                grace_end,
            )

            n_intervals = int(
                math.ceil(
                    spec["grace_minutes"]
                    / STEP_MIN
                )
            )

            for k in range(
                n_intervals
            ):
                start = rec[
                    "t0"
                ] + timedelta(
                    minutes=k * STEP_MIN
                )
                end = min(
                    start + timedelta(
                        minutes=STEP_MIN
                    ),
                    grace_end,
                )

                if (
                    stop_for_weight is not None
                    and start >= stop_for_weight
                ):
                    break

                row = dict(
                    clone_id=clone_id,
                    sid=sid,
                    subject_id=rec["subject_id"],
                    arm=arm,
                    spec_name=spec[
                        "name"
                    ],
                    interval_index=k,
                    interval_start=start,
                    interval_end=end,
                    age=rec["age"],
                    sofa_score=rec[
                        "sofa_score"
                    ],
                    ne0=rec["ne0"],
                    avp0=rec["avp0"],
                    delay_h=rec[
                        "delay_h"
                    ],
                    baseline_lactate=rec[
                        "baseline_lactate"
                    ],
                    baseline_map_mean=rec[
                        "baseline_map_mean"
                    ],
                    baseline_map_min=rec[
                        "baseline_map_min"
                    ],
                    baseline_hr_mean=rec[
                        "baseline_hr_mean"
                    ],
                    steroid_prior24h=rec[
                        "steroid_prior24h"
                    ],
                    inotrope_at_t0=rec[
                        "inotrope_at_t0"
                    ],
                    gender=rec[
                        "gender"
                    ],
                    race_group=rec[
                        "race_group"
                    ],
                    careunit_group=rec[
                        "careunit_group"
                    ],
                )

                row.update(
                    interval_features(
                        rec,
                        tl,
                        start,
                        by_stay_map,
                        by_stay_hr,
                        by_hadm_lac,
                        by_stay_other,
                        label_classes,
                        by_hadm_steroid,
                    )
                )

                censor_time = proto[
                    "censor_time"
                ]
                completion_time = proto[
                    "completion_time"
                ]

                outcome_before_end = (
                    outcome_time is not None
                    and start
                    < outcome_time
                    <= end
                )

                censor_before_end = (
                    censor_time is not None
                    and start
                    < censor_time
                    <= end
                )

                completion_before_end = (
                    completion_time is not None
                    and start
                    < completion_time
                    <= end
                )

                natural_before_end = (
                    rec["outtime"] is not None
                    and start
                    < rec["outtime"]
                    <= end
                )

                # Exact ordering within the interval.
                # Priority rule for exact ties:
                # an observed outcome at the same timestamp as an artificial
                # censor is retained as an outcome (not discarded by censoring).
                candidate_times = [
                    x for x in (
                        outcome_time if outcome_before_end else None,
                        censor_time if censor_before_end else None,
                        completion_time if completion_before_end else None,
                        rec["outtime"] if natural_before_end else None,
                    )
                    if x is not None
                ]
                earliest = min(candidate_times) if candidate_times else None

                outcome_is_first = bool(
                    outcome_before_end
                    and earliest == outcome_time
                )

                artificial_censor_this_interval = int(
                    censor_before_end
                    and earliest == censor_time
                    and not outcome_is_first
                )

                strategy_completed_this_interval = int(
                    completion_before_end
                    and earliest == completion_time
                    and not outcome_is_first
                    and not artificial_censor_this_interval
                )

                row[
                    "remain_uncensored"
                ] = (
                    0
                    if artificial_censor_this_interval
                    else 1
                )

                row[
                    "artificial_censor"
                ] = (
                    artificial_censor_this_interval
                )

                row[
                    "strategy_completed"
                ] = (
                    strategy_completed_this_interval
                )

                row[
                    "outcome_this_interval"
                ] = int(
                    outcome_is_first
                )

                rows.append(row)

                if (
                    artificial_censor_this_interval
                    or row[
                        "strategy_completed"
                    ]
                    or row[
                        "outcome_this_interval"
                    ]
                    or (
                        natural_before_end
                        and earliest
                        == rec["outtime"]
                    )
                ):
                    break

    return pd.DataFrame(rows), clone_protocol


def fit_weights(
    weight_df: pd.DataFrame,
    spec_name: str,
) -> Tuple[pd.DataFrame, List[dict]]:
    out = weight_df.copy()
    diagnostics = []

    out[
        "p_num_stay"
    ] = np.nan
    out[
        "p_den_stay"
    ] = np.nan

    for arm in (
        "NE_first",
        "AVP_first",
    ):
        mask = (
            out["arm"] == arm
        )
        d = out.loc[
            mask
        ].copy()

        y = d[
            "remain_uncensored"
        ].astype(int)

        if y.nunique() < 2:
            raise RuntimeError(
                f"{spec_name}/{arm}: "
                "censor model outcome has only one class."
            )

        num_model = fit_logistic_ridge(
            d,
            "remain_uncensored",
            NUMERATOR_CONT,
            NUMERATOR_CAT,
        )

        den_model = fit_logistic_ridge(
            d,
            "remain_uncensored",
            DENOMINATOR_CONT,
            DENOMINATOR_CAT,
        )

        p_num = predict_logistic(
            num_model,
            d,
        )
        p_den = predict_logistic(
            den_model,
            d,
        )

        out.loc[
            mask,
            "p_num_stay",
        ] = p_num

        out.loc[
            mask,
            "p_den_stay",
        ] = p_den

        diagnostics.append(
            dict(
                spec_name=spec_name,
                arm=arm,
                n_weight_rows=len(d),
                censor_events=int(
                    (
                        d[
                            "remain_uncensored"
                        ]
                        == 0
                    ).sum()
                ),
                numerator_converged=int(
                    num_model.converged
                ),
                denominator_converged=int(
                    den_model.converged
                ),
                numerator_iterations=num_model.iterations,
                denominator_iterations=den_model.iterations,
                min_raw_p_den=float(
                    np.min(p_den)
                ),
                p01_raw_p_den=float(
                    np.quantile(
                        p_den,
                        0.01,
                    )
                ),
                p05_raw_p_den=float(
                    np.quantile(
                        p_den,
                        0.05,
                    )
                ),
                median_raw_p_den=float(
                    np.median(p_den)
                ),
                max_raw_p_den=float(
                    np.max(p_den)
                ),
            )
        )

    out[
        "p_num_clip"
    ] = out[
        "p_num_stay"
    ].clip(
        PROB_CLIP_LO,
        PROB_CLIP_HI,
    )

    out[
        "p_den_clip"
    ] = out[
        "p_den_stay"
    ].clip(
        PROB_CLIP_LO,
        PROB_CLIP_HI,
    )

    # Start-of-interval weights.
    out[
        "sw_start"
    ] = np.nan

    out[
        "sw_end"
    ] = np.nan

    for clone_id, idx in out.groupby(
        "clone_id",
        sort=False,
    ).groups.items():
        indices = list(idx)
        indices.sort(
            key=lambda i: out.at[
                i,
                "interval_index",
            ]
        )

        sw = 1.0

        for i in indices:
            out.at[
                i,
                "sw_start",
            ] = sw

            # If artificially censored or outcome occurs in this interval,
            # no post-interval weight is needed.
            if int(
                out.at[
                    i,
                    "artificial_censor",
                ]
            ) == 1:
                out.at[
                    i,
                    "sw_end",
                ] = sw
                continue

            if int(
                out.at[
                    i,
                    "outcome_this_interval",
                ]
            ) == 1:
                out.at[
                    i,
                    "sw_end",
                ] = sw
                continue

            ratio = (
                float(
                    out.at[
                        i,
                        "p_num_clip",
                    ]
                )
                / float(
                    out.at[
                        i,
                        "p_den_clip",
                    ]
                )
            )

            sw = sw * ratio

            out.at[
                i,
                "sw_end",
            ] = sw

    return out, diagnostics


def clone_weight_paths(
    weight_df: pd.DataFrame,
) -> Dict[str, dict]:
    paths = {}

    for clone_id, g in weight_df.groupby(
        "clone_id",
        sort=False,
    ):
        g = g.sort_values(
            "interval_index"
        )

        rows = []

        for _, r in g.iterrows():
            rows.append(
                dict(
                    start=r[
                        "interval_start"
                    ],
                    end=r[
                        "interval_end"
                    ],
                    sw_start=float(
                        r[
                            "sw_start"
                        ]
                    ),
                    sw_end=float(
                        r[
                            "sw_end"
                        ]
                    ),
                    artificial_censor=int(
                        r[
                            "artificial_censor"
                        ]
                    ),
                    strategy_completed=int(
                        r[
                            "strategy_completed"
                        ]
                    ),
                    outcome_this_interval=int(
                        r[
                            "outcome_this_interval"
                        ]
                    ),
                )
            )

        final_sw = (
            rows[-1][
                "sw_end"
            ]
            if rows
            else 1.0
        )

        paths[clone_id] = dict(
            rows=rows,
            final_sw=float(
                final_sw
            ),
        )

    return paths


def weight_at_time(
    path: dict,
    t: datetime,
) -> float:
    rows = path[
        "rows"
    ]

    if not rows:
        return 1.0

    for row in rows:
        if (
            row["start"]
            <= t
            < row["end"]
        ):
            return float(
                row["sw_start"]
            )

    if t >= rows[-1]["end"]:
        return float(
            path["final_sw"]
        )

    return 1.0


# ---------------------------------------------------------------------
# Weighted survival risk
# ---------------------------------------------------------------------
def build_analysis_rows(
    records: Sequence[dict],
    timelines: Dict[int, Timeline],
    spec: dict,
    clone_protocol: Sequence[dict],
    weight_paths: Dict[str, dict],
    by_stay_map: Dict[int, list],
    horizon_h: int,
) -> pd.DataFrame:
    protocol_lookup = {
        x["clone_id"]: x
        for x in clone_protocol
    }

    rec_lookup = {
        r["sid"]: r
        for r in records
    }

    rows = []

    max_intervals = int(
        math.ceil(
            horizon_h
            * 60
            / STEP_MIN
        )
    )

    for clone_id, proto in protocol_lookup.items():
        rec = rec_lookup[
            proto["sid"]
        ]

        outcome_time = first_pair_lt65_time(
            by_stay_map.get(
                rec["sid"],
                [],
            ),
            rec["t0"],
            horizon_h,
        )

        artificial_censor = proto[
            "censor_time"
        ]

        natural_end = min(
            rec["outtime"],
            rec["t0"]
            + timedelta(
                hours=horizon_h
            ),
        )

        for k in range(
            max_intervals
        ):
            start = rec[
                "t0"
            ] + timedelta(
                minutes=k * STEP_MIN
            )
            end = min(
                start + timedelta(
                    minutes=STEP_MIN
                ),
                rec["t0"]
                + timedelta(
                    hours=horizon_h
                ),
            )

            if start >= natural_end:
                break

            if (
                artificial_censor is not None
                and start >= artificial_censor
            ):
                break

            if (
                outcome_time is not None
                and start >= outcome_time
            ):
                break

            event_this = int(
                outcome_time is not None
                and start
                < outcome_time
                <= end
                and (
                    artificial_censor is None
                    or outcome_time
                    <= artificial_censor
                )
                and outcome_time
                <= natural_end
            )

            censor_this = int(
                artificial_censor is not None
                and start
                < artificial_censor
                <= end
                and (
                    outcome_time is None
                    or artificial_censor
                    < outcome_time
                )
                and artificial_censor
                <= natural_end
            )

            w = weight_at_time(
                weight_paths.get(
                    clone_id,
                    dict(
                        rows=[],
                        final_sw=1.0,
                    ),
                ),
                start,
            )

            rows.append(
                dict(
                    clone_id=clone_id,
                    sid=rec["sid"],
                    subject_id=rec["subject_id"],
                    arm=proto[
                        "arm"
                    ],
                    spec_name=spec[
                        "name"
                    ],
                    horizon_h=horizon_h,
                    interval_index=k,
                    interval_start=start,
                    event=event_this,
                    artificial_censor=censor_this,
                    weight_raw=w,
                )
            )

            if event_this or censor_this:
                break

    return pd.DataFrame(
        rows
    )


def truncate_weights(
    df: pd.DataFrame,
    weight_col: str,
) -> pd.Series:
    out = pd.Series(
        index=df.index,
        dtype=float,
    )

    for arm in (
        "NE_first",
        "AVP_first",
    ):
        mask = (
            df["arm"] == arm
        )

        vals = df.loc[
            mask,
            weight_col,
        ].astype(float)

        if vals.empty:
            continue

        lo = float(
            vals.quantile(
                WEIGHT_LO_Q
            )
        )
        hi = float(
            vals.quantile(
                WEIGHT_HI_Q
            )
        )

        out.loc[
            mask
        ] = vals.clip(
            lo,
            hi,
        )

    return out


def weighted_km_risk(
    df: pd.DataFrame,
    weight_col: str,
) -> List[dict]:
    results = []

    for arm in (
        "NE_first",
        "AVP_first",
    ):
        d = df[
            df["arm"] == arm
        ].copy()

        survival = 1.0
        event_sum = 0

        for k in sorted(
            d[
                "interval_index"
            ].unique()
        ):
            x = d[
                d[
                    "interval_index"
                ]
                == k
            ]

            denom = float(
                x[
                    weight_col
                ].sum()
            )

            numer = float(
                (
                    x[
                        weight_col
                    ]
                    * x[
                        "event"
                    ]
                ).sum()
            )

            event_sum += int(
                x[
                    "event"
                ].sum()
            )

            if denom > 0:
                hazard = min(
                    max(
                        numer / denom,
                        0.0,
                    ),
                    1.0,
                )
                survival *= (
                    1.0 - hazard
                )

        results.append(
            dict(
                arm=arm,
                risk=1.0 - survival,
                observed_events=event_sum,
                person_period_rows=len(d),
            )
        )

    return results


# ---------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------
def effective_sample_size(
    weights: Sequence[float],
) -> Optional[float]:
    w = np.asarray(
        list(weights),
        dtype=float,
    )

    w = w[
        np.isfinite(w)
        & (w > 0)
    ]

    if len(w) == 0:
        return None

    denom = float(
        np.sum(w ** 2)
    )

    if denom <= 0:
        return None

    return float(
        np.sum(w) ** 2
        / denom
    )


def weight_distribution_rows(
    weight_df: pd.DataFrame,
    spec_name: str,
) -> List[dict]:
    rows = []

    for arm in (
        "NE_first",
        "AVP_first",
    ):
        d = weight_df[
            weight_df[
                "arm"
            ]
            == arm
        ]

        vals = (
            d[
                "sw_start"
            ]
            .dropna()
            .astype(float)
        )

        if vals.empty:
            continue

        rows.append(
            dict(
                spec_name=spec_name,
                arm=arm,
                n_person_periods=len(vals),
                mean_weight=float(
                    vals.mean()
                ),
                median_weight=float(
                    vals.median()
                ),
                p95_weight=float(
                    vals.quantile(
                        0.95
                    )
                ),
                p99_weight=float(
                    vals.quantile(
                        0.99
                    )
                ),
                max_weight=float(
                    vals.max()
                ),
                ess_person_period=effective_sample_size(
                    vals
                ),
            )
        )

    return rows


def smd_cont(
    x1,
    x0,
    w1=None,
    w0=None,
) -> Optional[float]:
    x1 = np.asarray(
        x1,
        dtype=float,
    )
    x0 = np.asarray(
        x0,
        dtype=float,
    )

    if w1 is None:
        w1 = np.ones_like(
            x1
        )
    else:
        w1 = np.asarray(
            w1,
            dtype=float,
        )

    if w0 is None:
        w0 = np.ones_like(
            x0
        )
    else:
        w0 = np.asarray(
            w0,
            dtype=float,
        )

    mask1 = (
        np.isfinite(x1)
        & np.isfinite(w1)
        & (w1 > 0)
    )
    mask0 = (
        np.isfinite(x0)
        & np.isfinite(w0)
        & (w0 > 0)
    )

    x1 = x1[
        mask1
    ]
    x0 = x0[
        mask0
    ]
    w1 = w1[
        mask1
    ]
    w0 = w0[
        mask0
    ]

    if (
        len(x1) < 2
        or len(x0) < 2
    ):
        return None

    m1 = float(
        np.average(
            x1,
            weights=w1,
        )
    )
    m0 = float(
        np.average(
            x0,
            weights=w0,
        )
    )

    v1 = float(
        np.average(
            (x1 - m1) ** 2,
            weights=w1,
        )
    )
    v0 = float(
        np.average(
            (x0 - m0) ** 2,
            weights=w0,
        )
    )

    pooled = math.sqrt(
        max(
            (v1 + v0) / 2.0,
            1e-12,
        )
    )

    return (
        (m1 - m0)
        / pooled
    )


def balance_at_completion(
    records: Sequence[dict],
    clone_protocol: Sequence[dict],
    weight_paths: Dict[str, dict],
    spec_name: str,
) -> List[dict]:
    rec_lookup = {
        r["sid"]: r
        for r in records
    }

    rows = []

    for p in clone_protocol:
        if p[
            "completion_time"
        ] is None:
            continue

        rec = rec_lookup[
            p["sid"]
        ]

        # Use the stabilized adherence weight accumulated through the
        # strategy-completion interval.
        w = float(
            weight_paths[
                p[
                    "clone_id"
                ]
            ][
                "final_sw"
            ]
        )

        rows.append(
            dict(
                arm=p[
                    "arm"
                ],
                weight=w,
                age=rec[
                    "age"
                ],
                sofa_score=rec[
                    "sofa_score"
                ],
                ne0=rec[
                    "ne0"
                ],
                avp0=rec[
                    "avp0"
                ],
                delay_h=rec[
                    "delay_h"
                ],
                baseline_lactate=rec[
                    "baseline_lactate"
                ],
                baseline_map_mean=rec[
                    "baseline_map_mean"
                ],
                baseline_hr_mean=rec[
                    "baseline_hr_mean"
                ],
            )
        )

    if not rows:
        return []

    d = pd.DataFrame(
        rows
    )

    out = []

    for var in [
        "age",
        "sofa_score",
        "ne0",
        "avp0",
        "delay_h",
        "baseline_lactate",
        "baseline_map_mean",
        "baseline_hr_mean",
    ]:
        a = d[
            d["arm"]
            == "NE_first"
        ]
        b = d[
            d["arm"]
            == "AVP_first"
        ]

        unweighted = smd_cont(
            a[var],
            b[var],
        )

        weighted = smd_cont(
            a[var],
            b[var],
            a["weight"],
            b["weight"],
        )

        out.append(
            dict(
                spec_name=spec_name,
                variable=var,
                n_NE=len(a),
                n_AVP=len(b),
                unweighted_SMD=unweighted,
                weighted_SMD=weighted,
                abs_weighted_SMD=(
                    abs(weighted)
                    if weighted
                    is not None
                    else None
                ),
            )
        )

    return out



# ---------------------------------------------------------------------
# G7 estimator repair
# ---------------------------------------------------------------------
G7_BASELINE_NUM_CONT = list(BASELINE_CONT)
G7_BASELINE_NUM_CAT = list(BASELINE_CAT) + ["interval_cat"]

G7_MARGINAL_NUM_CONT: List[str] = []
G7_MARGINAL_NUM_CAT = ["interval_cat"]

G7_DEN_CONT = list(BASELINE_CONT) + [
    "current_ne",
    "current_avp",
    "ne_fraction_t0",
    "avp_fraction_t0",
    "recent_map_mean",
    "recent_map_min",
    "recent_map_latest",
    "recent_hr_mean",
    "recent_lactate",
    "other_pressor_active",
    "inotrope_active",
    "steroid_active",
    "ne0_sq",
    "avp0_sq",
    "current_ne_sq",
    "current_avp_sq",
    "ne0_x_avp0",
    "current_ne_x_avp",
]

G7_DEN_CAT = list(BASELINE_CAT) + [
    "interval_cat",
    "ne0_bin",
    "avp0_bin",
    "current_ne_bin",
    "current_avp_bin",
    "recent_map_state",
    "recent_lactate_state",
]

G7_BALANCE_VARS = [
    "age",
    "sofa_score",
    "ne0",
    "avp0",
    "delay_h",
    "baseline_lactate",
    "baseline_map_mean",
    "baseline_map_min",
    "baseline_hr_mean",
]

G7_OUTCOME_CONT = list(BASELINE_CONT) + [
    "time_h",
    "time_h2",
    "time_h3",
    "arm_ne",
    "arm_time_h",
    "arm_time_h2",
    "arm_time_h3",
]
G7_OUTCOME_CAT = list(BASELINE_CAT)


def ne_dose_bin_g7(v) -> str:
    if v is None or not math.isfinite(float(v)):
        return "missing"
    x = float(v)
    if x <= 0.05 + 1e-12:
        return "<=0.05"
    if x <= 0.10 + 1e-12:
        return "0.05-0.10"
    if x <= 0.20 + 1e-12:
        return "0.10-0.20"
    return ">0.20"


def avp_dose_bin_g7(v) -> str:
    if v is None or not math.isfinite(float(v)):
        return "missing"
    x = float(v)
    if x <= 1.20 + 1e-12:
        return "<=1.2"
    if x <= 1.80 + 1e-12:
        return "1.2-1.8"
    if x <= 2.40 + 1e-12:
        return "1.8-2.4"
    return ">2.4"


def engineer_weight_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    out["interval_cat"] = (
        pd.to_numeric(out["interval_index"], errors="coerce")
        .fillna(-1)
        .astype(int)
        .astype(str)
    )

    for col in ("ne0", "avp0", "current_ne", "current_avp"):
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out["ne0_sq"] = out["ne0"] ** 2
    out["avp0_sq"] = out["avp0"] ** 2
    out["current_ne_sq"] = out["current_ne"] ** 2
    out["current_avp_sq"] = out["current_avp"] ** 2
    out["ne0_x_avp0"] = out["ne0"] * out["avp0"]
    out["current_ne_x_avp"] = out["current_ne"] * out["current_avp"]

    out["ne_fraction_t0"] = np.where(
        out["ne0"].abs() > 1e-12,
        out["current_ne"] / out["ne0"],
        np.nan,
    )
    out["avp_fraction_t0"] = np.where(
        out["avp0"].abs() > 1e-12,
        out["current_avp"] / out["avp0"],
        np.nan,
    )

    out["ne0_bin"] = [ne_dose_bin_g7(x) for x in out["ne0"]]
    out["avp0_bin"] = [avp_dose_bin_g7(x) for x in out["avp0"]]
    out["current_ne_bin"] = [ne_dose_bin_g7(x) for x in out["current_ne"]]
    out["current_avp_bin"] = [avp_dose_bin_g7(x) for x in out["current_avp"]]

    rmap = pd.to_numeric(out["recent_map_min"], errors="coerce")
    out["recent_map_state"] = np.select(
        [
            rmap.isna(),
            rmap < 65,
            (rmap >= 65) & (rmap < 70),
        ],
        ["missing", "<65", "65-69"],
        default=">=70",
    )

    rlac = pd.to_numeric(out["recent_lactate"], errors="coerce")
    out["recent_lactate_state"] = np.select(
        [
            rlac.isna(),
            rlac > 2,
        ],
        ["missing", ">2"],
        default="<=2",
    )

    return out


def fit_weighted_logistic_ridge(
    df: pd.DataFrame,
    y_col: str,
    sample_weight_col: str,
    continuous: Sequence[str],
    categorical: Sequence[str],
    ridge_lambda: float = RIDGE_LAMBDA,
) -> LogisticModel:
    X, design = _prepare_design_fit(
        df,
        continuous,
        categorical,
    )

    y = (
        pd.to_numeric(df[y_col], errors="coerce")
        .fillna(0)
        .to_numpy(dtype=float)
    )

    sw = (
        pd.to_numeric(df[sample_weight_col], errors="coerce")
        .fillna(0)
        .clip(lower=0)
        .to_numpy(dtype=float)
    )

    if not np.any(sw > 0):
        raise RuntimeError("Outcome model has no positive analysis weights.")

    beta = np.zeros(X.shape[1], dtype=float)

    penalty = np.eye(X.shape[1], dtype=float)
    penalty[0, 0] = 0.0

    converged = False

    for iteration in range(1, LOGIT_MAX_ITER + 1):
        eta = X @ beta
        p = sigmoid(eta)

        v = np.clip(
            sw * p * (1.0 - p),
            1e-8,
            None,
        )

        grad = (
            X.T @ (sw * (y - p))
            - ridge_lambda * penalty @ beta
        )

        hess_pos = (
            X.T @ (X * v[:, None])
            + ridge_lambda * penalty
            + 1e-8 * np.eye(X.shape[1])
        )

        try:
            step = np.linalg.solve(hess_pos, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(hess_pos) @ grad

        beta_new = beta + step

        if np.max(np.abs(beta_new - beta)) < 1e-7:
            beta = beta_new
            converged = True
            break

        beta = beta_new

    return LogisticModel(
        beta=beta,
        design=design,
        ridge_lambda=ridge_lambda,
        converged=converged,
        iterations=iteration,
    )


def fit_g7_censor_weights(
    weight_df: pd.DataFrame,
    spec_name: str,
) -> Tuple[pd.DataFrame, List[dict]]:
    """
    Fit one denominator plus two numerator specifications per arm.

    baseline-stabilized:
        P(uncensored | baseline, interval)
        -----------------------------------
        P(uncensored | baseline, interval, time-varying history)

    marginal/time-only stabilized:
        P(uncensored | interval)
        -----------------------------------
        P(uncensored | baseline, interval, time-varying history)
    """
    out = engineer_weight_features(weight_df)

    for c in (
        "p_num_baseline_raw",
        "p_num_marginal_raw",
        "p_den_raw",
    ):
        out[c] = np.nan

    diagnostics: List[dict] = []

    for arm in ("NE_first", "AVP_first"):
        mask = out["arm"] == arm
        d = out.loc[mask].copy()

        y = d["remain_uncensored"].astype(int)

        if y.nunique() < 2:
            raise RuntimeError(
                f"{spec_name}/{arm}: censor-model outcome has one class."
            )

        model_base = fit_logistic_ridge(
            d,
            "remain_uncensored",
            G7_BASELINE_NUM_CONT,
            G7_BASELINE_NUM_CAT,
        )
        model_margin = fit_logistic_ridge(
            d,
            "remain_uncensored",
            G7_MARGINAL_NUM_CONT,
            G7_MARGINAL_NUM_CAT,
        )
        model_den = fit_logistic_ridge(
            d,
            "remain_uncensored",
            G7_DEN_CONT,
            G7_DEN_CAT,
        )

        p_base = predict_logistic(model_base, d)
        p_margin = predict_logistic(model_margin, d)
        p_den = predict_logistic(model_den, d)

        out.loc[mask, "p_num_baseline_raw"] = p_base
        out.loc[mask, "p_num_marginal_raw"] = p_margin
        out.loc[mask, "p_den_raw"] = p_den

        diagnostics.append(
            dict(
                spec_name=spec_name,
                arm=arm,
                n_weight_rows=len(d),
                censor_events=int((y == 0).sum()),
                baseline_num_converged=int(model_base.converged),
                marginal_num_converged=int(model_margin.converged),
                denominator_converged=int(model_den.converged),
                baseline_num_iterations=model_base.iterations,
                marginal_num_iterations=model_margin.iterations,
                denominator_iterations=model_den.iterations,
                min_raw_p_den=float(np.min(p_den)),
                p01_raw_p_den=float(np.quantile(p_den, 0.01)),
                p05_raw_p_den=float(np.quantile(p_den, 0.05)),
                median_raw_p_den=float(np.median(p_den)),
                max_raw_p_den=float(np.max(p_den)),
                pct_den_below_0_01=float(np.mean(p_den < 0.01) * 100.0),
                pct_den_below_0_05=float(np.mean(p_den < 0.05) * 100.0),
            )
        )

    for base in (
        "p_num_baseline_raw",
        "p_num_marginal_raw",
        "p_den_raw",
    ):
        out[base.replace("_raw", "_clip")] = out[base].clip(
            PROB_CLIP_LO,
            PROB_CLIP_HI,
        )

    for c in (
        "sw_baseline_start",
        "sw_baseline_end",
        "sw_marginal_start",
        "sw_marginal_end",
    ):
        out[c] = np.nan

    for clone_id, idx in out.groupby("clone_id", sort=False).groups.items():
        indices = list(idx)
        indices.sort(key=lambda i: out.at[i, "interval_index"])

        sw_base = 1.0
        sw_margin = 1.0

        for i in indices:
            out.at[i, "sw_baseline_start"] = sw_base
            out.at[i, "sw_marginal_start"] = sw_margin

            # Contributions in an interval that ends in artificial censoring
            # or an outcome use the weight accumulated through interval start.
            if int(out.at[i, "artificial_censor"]) == 1:
                out.at[i, "sw_baseline_end"] = sw_base
                out.at[i, "sw_marginal_end"] = sw_margin
                continue

            if int(out.at[i, "outcome_this_interval"]) == 1:
                out.at[i, "sw_baseline_end"] = sw_base
                out.at[i, "sw_marginal_end"] = sw_margin
                continue

            den = float(out.at[i, "p_den_clip"])
            sw_base *= float(out.at[i, "p_num_baseline_clip"]) / den
            sw_margin *= float(out.at[i, "p_num_marginal_clip"]) / den

            out.at[i, "sw_baseline_end"] = sw_base
            out.at[i, "sw_marginal_end"] = sw_margin

    return out, diagnostics


def clone_weight_paths_g7(
    weight_df: pd.DataFrame,
    scheme: str,
) -> Dict[str, dict]:
    if scheme not in {"baseline", "marginal"}:
        raise ValueError(scheme)

    start_col = f"sw_{scheme}_start"
    end_col = f"sw_{scheme}_end"

    paths = {}

    for clone_id, g in weight_df.groupby("clone_id", sort=False):
        g = g.sort_values("interval_index")

        rows = []

        for _, r in g.iterrows():
            rows.append(
                dict(
                    start=r["interval_start"],
                    end=r["interval_end"],
                    sw_start=float(r[start_col]),
                    sw_end=float(r[end_col]),
                    artificial_censor=int(r["artificial_censor"]),
                    strategy_completed=int(r["strategy_completed"]),
                    outcome_this_interval=int(r["outcome_this_interval"]),
                )
            )

        paths[clone_id] = dict(
            rows=rows,
            final_sw=(float(rows[-1]["sw_end"]) if rows else 1.0),
        )

    return paths


def add_baseline_to_analysis_rows(
    analysis_df: pd.DataFrame,
    records: Sequence[dict],
) -> pd.DataFrame:
    rec_df = pd.DataFrame(
        [
            {
                "sid": r["sid"],
                "age": r["age"],
                "sofa_score": r["sofa_score"],
                "ne0": r["ne0"],
                "avp0": r["avp0"],
                "delay_h": r["delay_h"],
                "baseline_lactate": r["baseline_lactate"],
                "baseline_map_mean": r["baseline_map_mean"],
                "baseline_map_min": r["baseline_map_min"],
                "baseline_hr_mean": r["baseline_hr_mean"],
                "steroid_prior24h": r["steroid_prior24h"],
                "inotrope_at_t0": r["inotrope_at_t0"],
                "gender": r["gender"],
                "race_group": r["race_group"],
                "careunit_group": r["careunit_group"],
            }
            for r in records
        ]
    )

    out = analysis_df.merge(
        rec_df,
        on="sid",
        how="left",
        validate="many_to_one",
    )

    out["time_h"] = (
        pd.to_numeric(out["interval_index"], errors="coerce")
        * STEP_MIN
        / 60.0
    )
    out["time_h2"] = out["time_h"] ** 2
    out["time_h3"] = out["time_h"] ** 3
    out["arm_ne"] = (out["arm"] == "NE_first").astype(float)
    out["arm_time_h"] = out["arm_ne"] * out["time_h"]
    out["arm_time_h2"] = out["arm_ne"] * out["time_h2"]
    out["arm_time_h3"] = out["arm_ne"] * out["time_h3"]

    return out


def standardize_outcome_risk(
    analysis_df: pd.DataFrame,
    records: Sequence[dict],
    weight_col: str,
    horizon_h: int,
) -> Tuple[List[dict], dict]:
    """
    Fit a weighted pooled-logistic outcome model, then standardize predicted
    risks over the original (uncloned) CCW baseline distribution.
    """
    d = add_baseline_to_analysis_rows(
        analysis_df,
        records,
    )

    model = fit_weighted_logistic_ridge(
        d,
        "event",
        weight_col,
        G7_OUTCOME_CONT,
        G7_OUTCOME_CAT,
    )

    base = pd.DataFrame(
        [
            {
                "age": r["age"],
                "sofa_score": r["sofa_score"],
                "ne0": r["ne0"],
                "avp0": r["avp0"],
                "delay_h": r["delay_h"],
                "baseline_lactate": r["baseline_lactate"],
                "baseline_map_mean": r["baseline_map_mean"],
                "baseline_map_min": r["baseline_map_min"],
                "baseline_hr_mean": r["baseline_hr_mean"],
                "steroid_prior24h": r["steroid_prior24h"],
                "inotrope_at_t0": r["inotrope_at_t0"],
                "gender": r["gender"],
                "race_group": r["race_group"],
                "careunit_group": r["careunit_group"],
            }
            for r in records
        ]
    )

    n_intervals = int(math.ceil(horizon_h * 60 / STEP_MIN))
    results = []

    for arm in ("NE_first", "AVP_first"):
        surv = np.ones(len(base), dtype=float)
        arm_ne = 1.0 if arm == "NE_first" else 0.0

        for k in range(n_intervals):
            pred = base.copy()
            t = k * STEP_MIN / 60.0

            pred["time_h"] = t
            pred["time_h2"] = t ** 2
            pred["time_h3"] = t ** 3
            pred["arm_ne"] = arm_ne
            pred["arm_time_h"] = arm_ne * t
            pred["arm_time_h2"] = arm_ne * (t ** 2)
            pred["arm_time_h3"] = arm_ne * (t ** 3)

            haz = predict_logistic(model, pred)
            haz = np.clip(haz, 0.0, 1.0)
            surv *= (1.0 - haz)

        indiv_risk = 1.0 - surv

        results.append(
            dict(
                arm=arm,
                risk=float(np.mean(indiv_risk)),
                median_individual_risk=float(np.median(indiv_risk)),
            )
        )

    diagnostics = dict(
        outcome_model_converged=int(model.converged),
        outcome_model_iterations=model.iterations,
        outcome_model_n_rows=len(d),
        outcome_model_events=int(d["event"].sum()),
    )

    return results, diagnostics


def g7_weight_distribution_rows(
    weight_df: pd.DataFrame,
    clone_protocol: Sequence[dict],
    spec_name: str,
) -> List[dict]:
    completion_ids = {
        p["clone_id"]
        for p in clone_protocol
        if p["completion_time"] is not None
    }

    rows = []

    for scheme in ("baseline", "marginal"):
        start_col = f"sw_{scheme}_start"
        end_col = f"sw_{scheme}_end"

        for arm in ("NE_first", "AVP_first"):
            d = weight_df[weight_df["arm"] == arm]

            vals = pd.to_numeric(d[start_col], errors="coerce").dropna()

            if not vals.empty:
                rows.append(
                    dict(
                        spec_name=spec_name,
                        weight_scheme=scheme,
                        arm=arm,
                        scope="person_period_start",
                        n=len(vals),
                        mean_weight=float(vals.mean()),
                        median_weight=float(vals.median()),
                        p95_weight=float(vals.quantile(0.95)),
                        p99_weight=float(vals.quantile(0.99)),
                        max_weight=float(vals.max()),
                        ess=effective_sample_size(vals),
                    )
                )

            comp = d[d["clone_id"].isin(completion_ids)].copy()
            if not comp.empty:
                final = (
                    comp.sort_values(["clone_id", "interval_index"])
                    .groupby("clone_id", as_index=False)
                    .tail(1)
                )
                vals2 = pd.to_numeric(final[end_col], errors="coerce").dropna()

                rows.append(
                    dict(
                        spec_name=spec_name,
                        weight_scheme=scheme,
                        arm=arm,
                        scope="strategy_completers_final",
                        n=len(vals2),
                        mean_weight=float(vals2.mean()),
                        median_weight=float(vals2.median()),
                        p95_weight=float(vals2.quantile(0.95)),
                        p99_weight=float(vals2.quantile(0.99)),
                        max_weight=float(vals2.max()),
                        ess=effective_sample_size(vals2),
                    )
                )

    return rows


def riskset_balance_rows(
    weight_df: pd.DataFrame,
    spec_name: str,
    weight_col: str = "sw_marginal_start",
) -> List[dict]:
    out = []

    for k in sorted(weight_df["interval_index"].unique()):
        d = weight_df[weight_df["interval_index"] == k]

        a = d[d["arm"] == "NE_first"]
        b = d[d["arm"] == "AVP_first"]

        if a.empty or b.empty:
            continue

        for var in G7_BALANCE_VARS:
            s = smd_cont(
                a[var],
                b[var],
                a[weight_col],
                b[weight_col],
            )

            out.append(
                dict(
                    spec_name=spec_name,
                    interval_index=int(k),
                    minutes_from_t0=int(k) * STEP_MIN,
                    variable=var,
                    n_NE=len(a),
                    n_AVP=len(b),
                    weighted_SMD=s,
                    abs_weighted_SMD=(abs(s) if s is not None else None),
                )
            )

    return out


def completion_balance_g7(
    records: Sequence[dict],
    clone_protocol: Sequence[dict],
    paths: Dict[str, dict],
    spec_name: str,
) -> List[dict]:
    rec_lookup = {r["sid"]: r for r in records}
    rows = []

    for p in clone_protocol:
        if p["completion_time"] is None:
            continue

        rec = rec_lookup[p["sid"]]
        w = float(paths[p["clone_id"]]["final_sw"])

        row = dict(
            arm=p["arm"],
            weight=w,
        )
        for var in G7_BALANCE_VARS:
            row[var] = rec[var]
        rows.append(row)

    if not rows:
        return []

    d = pd.DataFrame(rows)
    out = []

    for var in G7_BALANCE_VARS:
        a = d[d["arm"] == "NE_first"]
        b = d[d["arm"] == "AVP_first"]

        s_un = smd_cont(a[var], b[var])
        s_w = smd_cont(a[var], b[var], a["weight"], b["weight"])

        out.append(
            dict(
                spec_name=spec_name,
                variable=var,
                n_NE=len(a),
                n_AVP=len(b),
                unweighted_SMD=s_un,
                marginal_weighted_SMD=s_w,
                abs_marginal_weighted_SMD=(abs(s_w) if s_w is not None else None),
            )
        )

    return out


def natural_censor_audit_rows(
    records: Sequence[dict],
) -> List[dict]:
    rows = []

    for h in OUTCOME_WINDOWS_H:
        cutoff_key = f"{h}h"

        discharged = 0
        died = 0

        for r in records:
            horizon = r["t0"] + timedelta(hours=h)

            if r["outtime"] is not None and r["outtime"] < horizon:
                discharged += 1

            dt = r.get("deathtime")
            if dt is not None and r["t0"] < dt <= horizon:
                died += 1

        rows.append(
            dict(
                horizon=cutoff_key,
                n=len(records),
                ICU_out_before_horizon=discharged,
                pct_ICU_out_before_horizon=(100.0 * discharged / len(records)),
                death_before_horizon=died,
                pct_death_before_horizon=(100.0 * died / len(records)),
            )
        )

    return rows


def run_g7_synthetic_tests() -> None:
    # Weighted ridge must remain finite and converge on a simple toy example.
    toy = pd.DataFrame(
        {
            "y": [0, 0, 1, 0, 1, 1, 0, 1],
            "x": [-2, -1, -0.5, 0, 0.5, 1, 2, 3],
            "cat": ["a", "a", "b", "b", "a", "b", "a", "b"],
            "w": [1, 2, 1, 1, 3, 1, 1, 2],
        }
    )
    model = fit_weighted_logistic_ridge(
        toy,
        "y",
        "w",
        ["x"],
        ["cat"],
    )
    p = predict_logistic(model, toy)
    assert np.all(np.isfinite(p))
    assert np.all((p > 0) & (p < 1))

    # Engineered weight features should respect exact dose-bin boundaries.
    x = pd.DataFrame(
        {
            "interval_index": [0, 1],
            "ne0": [0.10, 0.10001],
            "avp0": [1.80, 1.80001],
            "current_ne": [0.05, 0.20],
            "current_avp": [1.20, 2.40],
            "recent_map_min": [64.9, 70.0],
            "recent_lactate": [2.01, 2.0],
        }
    )
    e = engineer_weight_features(x)
    assert e.loc[0, "ne0_bin"] == "0.05-0.10"
    assert e.loc[1, "ne0_bin"] == "0.10-0.20"
    assert e.loc[0, "avp0_bin"] == "1.2-1.8"
    assert e.loc[1, "avp0_bin"] == "1.8-2.4"
    assert e.loc[0, "recent_map_state"] == "<65"
    assert e.loc[0, "recent_lactate_state"] == ">2"
    assert e.loc[1, "recent_lactate_state"] == "<=2"


# ---------------------------------------------------------------------
# G8 non-anticipating treatment rule and diagnostics
# ---------------------------------------------------------------------
def immediate_event_for_drug(
    tl: Timeline,
    t0: datetime,
    grace_min: int,
    drug: str,
    baseline: float,
    mode: str,
    observation_end: datetime,
) -> Optional[ActionEvent]:
    """Return the first observed qualifying rate transition; no future look-ahead."""
    if baseline <= 0:
        return None

    if mode == "cumulative30_immediate":
        threshold = baseline * 0.70
    elif mode == "dovss_scale_immediate":
        abs_drop = DOVSS_NE_ABS_DROP if drug == "NE" else DOVSS_AVP_ABS_DROP
        threshold = max(0.0, baseline - abs_drop)
    else:
        raise ValueError(f"Unknown immediate mode: {mode}")

    grace_end = min(t0 + timedelta(minutes=grace_min), observation_end)

    for tt, fr, to, _persist in tl.transitions(drug):
        if tt <= t0 or tt > grace_end or fr <= 0:
            continue
        if to <= 0 or to <= threshold + 1e-12:
            return ActionEvent(
                time=tt,
                drug=drug,
                kind=("stop" if to <= 0 else "taper"),
                from_rate=baseline,
                to_rate=to,
            )
    return None


def get_spec_events(
    tl: Timeline,
    rec: dict,
    spec: dict,
) -> Tuple[Optional[ActionEvent], Optional[ActionEvent]]:
    mode = spec["mode"]

    if mode in {"dovss_scale_immediate", "cumulative30_immediate"}:
        ne = immediate_event_for_drug(
            tl, rec["t0"], spec["grace_minutes"], "NE",
            rec["ne0"], mode, rec["outtime"]
        )
        avp = immediate_event_for_drug(
            tl, rec["t0"], spec["grace_minutes"], "AVP",
            rec["avp0"], mode, rec["outtime"]
        )
        return ne, avp

    # Historical fallback only.
    ne = cumulative_event_for_drug(
        tl, rec["t0"], spec["grace_minutes"], "NE",
        rec["ne0"], mode, rec["outtime"]
    )
    avp = cumulative_event_for_drug(
        tl, rec["t0"], spec["grace_minutes"], "AVP",
        rec["avp0"], mode, rec["outtime"]
    )
    return ne, avp


def positivity_relevant_rows(
    weight_df: pd.DataFrame,
    spec_name: str,
) -> List[dict]:
    rows: List[dict] = []

    for arm in ("NE_first", "AVP_first"):
        d0 = weight_df[weight_df["arm"] == arm].copy()

        scopes = {
            "all_rows": d0,
            "artificial_censor_rows": d0[d0["artificial_censor"] == 1],
            "outcome_rows": d0[d0["outcome_this_interval"] == 1],
            "weight_update_rows": d0[
                (d0["artificial_censor"] == 0)
                & (d0["outcome_this_interval"] == 0)
            ],
            "ongoing_uncensored_rows": d0[
                (d0["artificial_censor"] == 0)
                & (d0["outcome_this_interval"] == 0)
                & (d0["strategy_completed"] == 0)
            ],
            "strategy_completion_rows": d0[d0["strategy_completed"] == 1],
        }

        for scope, d in scopes.items():
            if d.empty:
                continue

            p = pd.to_numeric(d["p_den_raw"], errors="coerce").dropna()
            num = pd.to_numeric(d["p_num_marginal_raw"], errors="coerce")
            den = pd.to_numeric(d["p_den_raw"], errors="coerce")
            valid = num.notna() & den.notna() & (den > 0)
            ratio = num[valid] / den[valid]

            rows.append(
                dict(
                    spec_name=spec_name,
                    arm=arm,
                    scope=scope,
                    n=len(d),
                    p_den_min=(float(p.min()) if len(p) else None),
                    p_den_p01=(float(p.quantile(0.01)) if len(p) else None),
                    p_den_p05=(float(p.quantile(0.05)) if len(p) else None),
                    p_den_median=(float(p.median()) if len(p) else None),
                    pct_p_den_lt_0_01=(
                        float((p < 0.01).mean() * 100.0) if len(p) else None
                    ),
                    pct_p_den_lt_0_05=(
                        float((p < 0.05).mean() * 100.0) if len(p) else None
                    ),
                    marginal_ratio_p95=(
                        float(ratio.quantile(0.95)) if len(ratio) else None
                    ),
                    marginal_ratio_p99=(
                        float(ratio.quantile(0.99)) if len(ratio) else None
                    ),
                    marginal_ratio_max=(
                        float(ratio.max()) if len(ratio) else None
                    ),
                )
            )

    return rows


def action_timing_and_outcome_rows(
    records: Sequence[dict],
    timelines: Dict[int, Timeline],
    spec: dict,
    by_stay_map: Dict[int, list],
) -> Tuple[List[dict], List[dict]]:
    timing = defaultdict(list)
    counts = Counter()
    counts["eligible"] = len(records)

    for rec in records:
        ne, avp = get_spec_events(timelines[rec["sid"]], rec, spec)
        strategy = classify_events(ne, avp)

        if ne is not None and avp is not None:
            first_time = min(ne.time, avp.time)
        elif ne is not None:
            first_time = ne.time
        elif avp is not None:
            first_time = avp.time
        else:
            first_time = None

        if strategy in {"NE_first", "AVP_first"} and first_time is not None:
            timing[strategy].append(
                (first_time - rec["t0"]).total_seconds() / 60.0
            )

        event6 = first_pair_lt65_time(
            by_stay_map.get(rec["sid"], []),
            rec["t0"],
            6,
        )

        if event6 is None:
            counts["no_6h_outcome"] += 1
        elif first_time is None:
            counts["6h_outcome_no_qualifying_action"] += 1
        elif event6 <= first_time:
            counts["6h_outcome_before_or_at_first_action"] += 1
        else:
            counts["6h_outcome_after_first_action"] += 1

    timing_rows = []
    for strategy in ("NE_first", "AVP_first"):
        vals = sorted(timing.get(strategy, []))
        timing_rows.append(
            dict(
                spec_name=spec["name"],
                strategy=strategy,
                n=len(vals),
                p25_minutes_to_action=(qtile_or_none(vals, 0.25) if vals else None),
                median_minutes_to_action=(median_or_none(vals) if vals else None),
                p75_minutes_to_action=(qtile_or_none(vals, 0.75) if vals else None),
                pct_within_30min=(
                    100.0 * sum(x <= 30 for x in vals) / len(vals)
                    if vals else None
                ),
                pct_within_60min=(
                    100.0 * sum(x <= 60 for x in vals) / len(vals)
                    if vals else None
                ),
            )
        )

    outcome_rows = [
        dict(
            spec_name=spec["name"],
            eligible=counts["eligible"],
            no_6h_outcome=counts["no_6h_outcome"],
            outcome_no_qualifying_action=counts["6h_outcome_no_qualifying_action"],
            outcome_before_or_at_first_action=counts[
                "6h_outcome_before_or_at_first_action"
            ],
            outcome_after_first_action=counts["6h_outcome_after_first_action"],
            pct_all_with_outcome_before_or_at_action=(
                100.0
                * counts["6h_outcome_before_or_at_first_action"]
                / counts["eligible"]
                if counts["eligible"]
                else None
            ),
        )
    ]

    return timing_rows, outcome_rows


def run_g8_synthetic_tests() -> None:
    base = datetime(2026, 1, 1, 0, 0)

    def row(itemid, rate, unit, st, en):
        return dict(
            itemid=itemid,
            rate=rate,
            rateuom=unit,
            patientweight=70,
            starttime=base + timedelta(minutes=st),
            endtime=base + timedelta(minutes=en),
            statusdescription="ChangeDose/Rate",
            ordercategorydescription="Continuous",
        )

    # Exact DOVSS NE threshold qualifies immediately even if reversed 10 min later.
    rows = [
        row(NE_ITEMID, 0.30, "mcg/kg/min", 0, 60),
        row(NE_ITEMID, 0.20, "mcg/kg/min", 60, 70),
        row(NE_ITEMID, 0.30, "mcg/kg/min", 70, 240),
        row(AVP_ITEMID, 2.40, "units/hour", 0, 240),
    ]
    segs, _ = segments_from_rows(rows)
    tl, _ = build_atomic_timeline(segs)
    ev = immediate_event_for_drug(
        tl, base + timedelta(minutes=30), 120,
        "NE", 0.30, "dovss_scale_immediate",
        base + timedelta(minutes=240),
    )
    assert ev is not None and ev.time == base + timedelta(minutes=60)

    # Just below the exact threshold does not qualify.
    rows = [
        row(NE_ITEMID, 0.30, "mcg/kg/min", 0, 60),
        row(NE_ITEMID, 0.2001, "mcg/kg/min", 60, 240),
        row(AVP_ITEMID, 2.40, "units/hour", 0, 240),
    ]
    segs, _ = segments_from_rows(rows)
    tl, _ = build_atomic_timeline(segs)
    ev = immediate_event_for_drug(
        tl, base + timedelta(minutes=30), 120,
        "NE", 0.30, "dovss_scale_immediate",
        base + timedelta(minutes=240),
    )
    assert ev is None

    # Exact AVP 0.60 U/h drop qualifies immediately.
    rows = [
        row(NE_ITEMID, 0.20, "mcg/kg/min", 0, 240),
        row(AVP_ITEMID, 2.40, "units/hour", 0, 60),
        row(AVP_ITEMID, 1.80, "units/hour", 60, 240),
    ]
    segs, _ = segments_from_rows(rows)
    tl, _ = build_atomic_timeline(segs)
    ev = immediate_event_for_drug(
        tl, base + timedelta(minutes=30), 120,
        "AVP", 2.40, "dovss_scale_immediate",
        base + timedelta(minutes=240),
    )
    assert ev is not None and ev.time == base + timedelta(minutes=60)

    # 15 minutes apart is ordered in G8; exact same timestamp is unresolved.
    ne = ActionEvent(
        base + timedelta(minutes=60), "NE", "taper", 0.30, 0.20
    )
    avp_later = ActionEvent(
        base + timedelta(minutes=75), "AVP", "taper", 2.40, 1.80
    )
    avp_same = ActionEvent(
        base + timedelta(minutes=60), "AVP", "taper", 2.40, 1.80
    )
    assert classify_events(ne, avp_later) == "NE_first"
    assert classify_events(ne, avp_same) == "near_simultaneous"


def t0_boundary_audit(
    records: Sequence[dict],
    timelines: Dict[int, Timeline],
) -> List[dict]:
    """
    Audit whether an infusion rate transition occurs exactly at the grid-defined
    t0. G8 does not silently change the G5 fixed t0; it quantifies this boundary
    issue so the final bootstrap implementation can decide whether a left-limit
    t0 representation is required.
    """
    c = Counter()
    c["eligible_t0"] = len(records)

    eps = timedelta(microseconds=1)

    for rec in records:
        tl = timelines[rec["sid"]]
        t0 = rec["t0"]

        pre_ne = tl.rate_at(t0 - eps, "NE")
        pre_avp = tl.rate_at(t0 - eps, "AVP")
        post_ne = tl.rate_at(t0, "NE")
        post_avp = tl.rate_at(t0, "AVP")

        vals = (pre_ne, pre_avp, post_ne, post_avp)
        if any(v is None for v in vals):
            c["conflict_or_missing_at_boundary"] += 1
            continue

        pre_ne = float(pre_ne)
        pre_avp = float(pre_avp)
        post_ne = float(post_ne)
        post_avp = float(post_avp)

        ne_changed = abs(post_ne - pre_ne) > 1e-12
        avp_changed = abs(post_avp - pre_avp) > 1e-12

        if ne_changed or avp_changed:
            c["any_rate_change_exactly_at_t0"] += 1

        if post_ne < pre_ne - 1e-12 or post_avp < pre_avp - 1e-12:
            c["any_downward_rate_change_exactly_at_t0"] += 1

        dovss_ne = (
            pre_ne > 0
            and (
                post_ne <= 0
                or pre_ne - post_ne >= DOVSS_NE_ABS_DROP - 1e-12
            )
        )
        dovss_avp = (
            pre_avp > 0
            and (
                post_avp <= 0
                or pre_avp - post_avp >= DOVSS_AVP_ABS_DROP - 1e-12
            )
        )

        if dovss_ne or dovss_avp:
            c["DOVSS_scale_downward_action_exactly_at_t0"] += 1

        pre_dose_eligible = (
            pre_ne > 0
            and pre_avp > 0
            and pre_ne <= PRIMARY_NE_MAX + 1e-12
            and pre_avp <= PRIMARY_AVP_MAX + 1e-12
        )

        if not pre_dose_eligible:
            c["post_t0_dose_eligible_but_pre_t0_not_dose_eligible"] += 1

    n = max(c["eligible_t0"], 1)

    order = [
        "eligible_t0",
        "conflict_or_missing_at_boundary",
        "any_rate_change_exactly_at_t0",
        "any_downward_rate_change_exactly_at_t0",
        "DOVSS_scale_downward_action_exactly_at_t0",
        "post_t0_dose_eligible_but_pre_t0_not_dose_eligible",
    ]

    return [
        dict(
            metric=k,
            count=int(c[k]),
            percent_of_eligible=(
                100.0 * c[k] / n
                if k != "eligible_t0"
                else 100.0
            ),
        )
        for k in order
    ]



# ---------------------------------------------------------------------
# G10 frozen decision point from G9
# ---------------------------------------------------------------------
def g10_exact_transition(
    tl: Timeline,
    t: datetime,
) -> bool:
    for drug in ("NE", "AVP"):
        for tt, fr, to, _persist in tl.transitions(drug):
            if tt == t and abs(to - fr) > 1e-12:
                return True
    return False


def g10_recent_dovss_scale_change(
    tl: Timeline,
    t: datetime,
    drug: str,
    minutes: int = 120,
) -> bool:
    start = t - timedelta(minutes=minutes)
    abs_thr = DOVSS_NE_ABS_DROP if drug == "NE" else DOVSS_AVP_ABS_DROP

    for tt, fr, to, _persist in tl.transitions(drug):
        if not (start < tt < t):
            continue
        if fr <= 0:
            continue
        if to <= 0:
            return True
        if abs(to - fr) + 1e-12 >= abs_thr:
            return True

    return False


def g10_other_pressor_overlap(
    rows: Sequence[dict],
    start: datetime,
    end: datetime,
    label_classes: Dict[int, str],
) -> bool:
    for r in rows:
        itemid = int(r["itemid"])
        if label_classes.get(itemid, "other") != "other_pressor":
            continue

        st = r.get("starttime")
        en = r.get("endtime")

        if (
            st is not None
            and en is not None
            and st < end
            and en > start
        ):
            return True

    return False


def g10_find_final_t0(
    tl: Timeline,
    map_rows: Sequence[Tuple[datetime, float]],
    other_rows: Sequence[dict],
    label_classes: Dict[int, str],
    meta: dict,
) -> Tuple[Optional[datetime], dict]:
    """
    Frozen G9 primary decision point:
    dovss_stable120_no_otherpressor2h.
    """
    first_dual = tl.first_dual_time()

    stage = dict(
        two_hour_dual=False,
        sepsis_before_t=False,
        adult_ok=False,
        map_stable=False,
        dose_ok=False,
        no_escalation=False,
        boundary_safe=False,
        dovss_stable120=False,
        other_pressor_clean120=False,
        eligible=False,
    )

    if first_dual is None:
        return None, stage

    start = ceil_dt(
        first_dual + timedelta(minutes=DUAL_LOOKBACK_MIN),
        STEP_MIN,
    )
    stop = min(
        meta["outtime"],
        first_dual + timedelta(days=SEARCH_MAX_DAYS),
    )

    eps = timedelta(microseconds=1)
    t = start

    while t <= stop:
        if not tl.continuous_dual(
            t - timedelta(minutes=DUAL_LOOKBACK_MIN),
            t,
        ):
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["two_hour_dual"] = True

        sus = meta.get("suspected_infection_time")
        sofa = meta.get("sofa_time")
        if sus is None or sofa is None or sus > t or sofa > t:
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["sepsis_before_t"] = True

        age = meta.get("age")
        if age is None or float(age) < 18.0:
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["adult_ok"] = True

        if not map_stable(map_rows, t):
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["map_stable"] = True

        # Eligibility doses are explicitly pre-action / left-limit.
        teval = t - eps
        ne = tl.rate_at(teval, "NE")
        avp = tl.rate_at(teval, "AVP")

        if (
            ne is None
            or avp is None
            or ne <= 0
            or avp <= 0
            or ne > PRIMARY_NE_MAX + 1e-12
            or avp > PRIMARY_AVP_MAX + 1e-12
        ):
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["dose_ok"] = True

        # Preserve the prior G8 recovery criterion. Exact-t transitions are
        # separately forbidden below, so this is genuinely pre-action.
        if (
            has_recent_escalation(tl, t, "NE")
            or has_recent_escalation(tl, t, "AVP")
        ):
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["no_escalation"] = True

        if g10_exact_transition(tl, t):
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["boundary_safe"] = True

        if (
            g10_recent_dovss_scale_change(tl, t, "NE", 120)
            or g10_recent_dovss_scale_change(tl, t, "AVP", 120)
        ):
            t += timedelta(minutes=STEP_MIN)
            continue
        stage["dovss_stable120"] = True

        if g10_other_pressor_overlap(
            other_rows,
            t - timedelta(hours=2),
            t,
            label_classes,
        ):
            t += timedelta(minutes=STEP_MIN)
            continue

        # Explicitly catch a third pressor starting exactly at t0.
        if interval_active(
            other_rows,
            t,
            ["other_pressor"],
            label_classes,
        ):
            t += timedelta(minutes=STEP_MIN)
            continue

        stage["other_pressor_clean120"] = True
        stage["eligible"] = True

        return t, stage

    return None, stage


# ---------------------------------------------------------------------
# G10 bootstrap helpers
# ---------------------------------------------------------------------
def fit_g10_marginal_censor_weights(
    weight_df: pd.DataFrame,
    spec_name: str,
) -> Tuple[pd.DataFrame, List[dict]]:
    """
    Marginal/time-only stabilized IPCW only.

    Numerator: flexible interval time, separately by arm.
    Denominator: G7/G8 rich baseline + time-varying censoring model.
    """
    out = engineer_weight_features(weight_df)

    out["p_num_marginal_raw"] = np.nan
    out["p_den_raw"] = np.nan

    diagnostics = []

    for arm in ("NE_first", "AVP_first"):
        mask = out["arm"] == arm
        d = out.loc[mask].copy()

        y = d["remain_uncensored"].astype(int)
        if y.nunique() < 2:
            raise RuntimeError(
                f"{spec_name}/{arm}: censor outcome has one class."
            )

        num_model = fit_logistic_ridge(
            d,
            "remain_uncensored",
            G7_MARGINAL_NUM_CONT,
            G7_MARGINAL_NUM_CAT,
        )
        den_model = fit_logistic_ridge(
            d,
            "remain_uncensored",
            G7_DEN_CONT,
            G7_DEN_CAT,
        )

        p_num = predict_logistic(num_model, d)
        p_den = predict_logistic(den_model, d)

        out.loc[mask, "p_num_marginal_raw"] = p_num
        out.loc[mask, "p_den_raw"] = p_den

        diagnostics.append(
            dict(
                spec_name=spec_name,
                arm=arm,
                numerator_converged=int(num_model.converged),
                denominator_converged=int(den_model.converged),
                n_rows=len(d),
                censor_events=int((y == 0).sum()),
                p_den_min=float(np.min(p_den)),
                p_den_p01=float(np.quantile(p_den, 0.01)),
                p_den_p05=float(np.quantile(p_den, 0.05)),
            )
        )

    out["p_num_marginal_clip"] = out["p_num_marginal_raw"].clip(
        PROB_CLIP_LO,
        PROB_CLIP_HI,
    )
    out["p_den_clip"] = out["p_den_raw"].clip(
        PROB_CLIP_LO,
        PROB_CLIP_HI,
    )

    out["sw_marginal_start"] = np.nan
    out["sw_marginal_end"] = np.nan

    for clone_id, idxs in out.groupby("clone_id", sort=False).groups.items():
        indices = list(idxs)
        indices.sort(key=lambda i: out.at[i, "interval_index"])

        sw = 1.0

        for i in indices:
            out.at[i, "sw_marginal_start"] = sw

            if int(out.at[i, "artificial_censor"]) == 1:
                out.at[i, "sw_marginal_end"] = sw
                continue

            if int(out.at[i, "outcome_this_interval"]) == 1:
                out.at[i, "sw_marginal_end"] = sw
                continue

            sw *= (
                float(out.at[i, "p_num_marginal_clip"])
                / float(out.at[i, "p_den_clip"])
            )

            out.at[i, "sw_marginal_end"] = sw

    return out, diagnostics


def g10_attach_marginal_weights(
    analysis_df: pd.DataFrame,
    weight_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Attach start-of-interval marginal IPC weights.
    After strategy completion, carry the final adherence weight forward.
    """
    wstart = weight_df[
        [
            "clone_id",
            "interval_index",
            "sw_marginal_start",
        ]
    ].copy()

    final = (
        weight_df.sort_values(["clone_id", "interval_index"])
        .groupby("clone_id", as_index=False)
        .tail(1)[["clone_id", "sw_marginal_end"]]
        .rename(columns={"sw_marginal_end": "_final_sw"})
    )

    out = analysis_df.merge(
        wstart,
        on=["clone_id", "interval_index"],
        how="left",
        validate="many_to_one",
    )
    out = out.merge(
        final,
        on="clone_id",
        how="left",
        validate="many_to_one",
    )

    out["weight_marginal_raw"] = (
        out["sw_marginal_start"]
        .fillna(out["_final_sw"])
        .fillna(1.0)
        .astype(float)
    )

    out = out.drop(
        columns=["sw_marginal_start", "_final_sw"]
    )

    out["weight_marginal_trunc"] = truncate_weights(
        out,
        "weight_marginal_raw",
    )

    return out


def g10_bootstrap_copy(
    df: pd.DataFrame,
    draw_map: pd.DataFrame,
) -> pd.DataFrame:
    out = df.merge(
        draw_map,
        on="subject_id",
        how="inner",
        validate="many_to_many",
    )

    out["clone_id"] = (
        out["clone_id"].astype(str)
        + "|boot"
        + out["boot_copy"].astype(str)
    )

    return out.drop(columns=["boot_copy"])


def g10_bootstrap_one_spec(
    weight_base: pd.DataFrame,
    analysis_base_by_horizon: Dict[int, pd.DataFrame],
    draw_map: pd.DataFrame,
    spec_name: str,
) -> Tuple[List[dict], dict]:
    wb = g10_bootstrap_copy(
        weight_base,
        draw_map,
    )

    wb, diag = fit_g10_marginal_censor_weights(
        wb,
        spec_name,
    )

    converged = all(
        x["numerator_converged"] == 1
        and x["denominator_converged"] == 1
        for x in diag
    )

    estimates = []

    for horizon_h, abase in analysis_base_by_horizon.items():
        ab = g10_bootstrap_copy(
            abase,
            draw_map,
        )

        ab = g10_attach_marginal_weights(
            ab,
            wb,
        )

        for col, version in (
            ("weight_marginal_raw", "untruncated"),
            ("weight_marginal_trunc", "truncated_1_99"),
        ):
            km = weighted_km_risk(
                ab,
                col,
            )
            by_arm = {
                r["arm"]: r
                for r in km
            }

            ne = float(by_arm["NE_first"]["risk"])
            avp = float(by_arm["AVP_first"]["risk"])

            estimates.append(
                dict(
                    spec_name=spec_name,
                    horizon_hours=horizon_h,
                    weight_version=version,
                    risk_NE_first=ne,
                    risk_AVP_first=avp,
                    risk_difference_NE_minus_AVP=ne - avp,
                    risk_ratio_NE_over_AVP=(
                        ne / avp
                        if avp > 0
                        else None
                    ),
                )
            )

    if not converged:
        raise RuntimeError(
            f"{spec_name}: bootstrap censor model did not converge."
        )

    if not converged:
        raise RuntimeError(
            f"{spec_name}: bootstrap censor model did not converge."
        )

    dd = pd.DataFrame(diag)

    return estimates, dict(
        censor_models_converged=int(converged),
        min_p01_denominator_probability=(
            float(dd["p_den_p01"].min())
            if len(dd)
            else None
        ),
    )


def g10_percentile_ci(
    vals: Sequence[float],
    alpha: float = 0.05,
) -> Tuple[Optional[float], Optional[float]]:
    x = np.asarray(
        [
            float(v)
            for v in vals
            if v is not None and math.isfinite(float(v))
        ],
        dtype=float,
    )

    if len(x) == 0:
        return None, None

    return (
        float(np.quantile(x, alpha / 2.0)),
        float(np.quantile(x, 1.0 - alpha / 2.0)),
    )


def g10_bootstrap_summary(
    point_rows: Sequence[dict],
    replicate_rows: Sequence[dict],
    requested_reps: int,
) -> List[dict]:
    point_df = pd.DataFrame(point_rows)
    boot = pd.DataFrame(replicate_rows)

    out = []

    keys = [
        "spec_name",
        "horizon_hours",
        "weight_version",
    ]

    for key, pgrp in point_df.groupby(keys, dropna=False):
        spec_name, horizon_h, weight_version = key

        if horizon_h != 6:
            # Primary inferential bootstrap is deliberately limited to 6 h.
            continue

        p = pgrp.iloc[0]

        b = boot[
            (boot["spec_name"] == spec_name)
            & (boot["horizon_hours"] == horizon_h)
            & (boot["weight_version"] == weight_version)
            & (boot["success"] == 1)
        ]

        rd_lo, rd_hi = g10_percentile_ci(
            b["risk_difference_NE_minus_AVP"]
        )
        rr_lo, rr_hi = g10_percentile_ci(
            b["risk_ratio_NE_over_AVP"]
        )
        ne_lo, ne_hi = g10_percentile_ci(
            b["risk_NE_first"]
        )
        avp_lo, avp_hi = g10_percentile_ci(
            b["risk_AVP_first"]
        )

        success_n = len(b)

        out.append(
            dict(
                spec_name=spec_name,
                horizon_hours=horizon_h,
                weight_version=weight_version,
                bootstrap_requested=requested_reps,
                bootstrap_successful=success_n,
                bootstrap_success_fraction=(
                    success_n / requested_reps
                    if requested_reps
                    else None
                ),
                risk_NE_first=float(p["risk_NE_first"]),
                risk_NE_first_CI_low=ne_lo,
                risk_NE_first_CI_high=ne_hi,
                risk_AVP_first=float(p["risk_AVP_first"]),
                risk_AVP_first_CI_low=avp_lo,
                risk_AVP_first_CI_high=avp_hi,
                risk_difference_NE_minus_AVP=float(
                    p["risk_difference_NE_minus_AVP"]
                ),
                RD_CI_low=rd_lo,
                RD_CI_high=rd_hi,
                risk_ratio_NE_over_AVP=float(
                    p["risk_ratio_NE_over_AVP"]
                ) if pd.notna(p["risk_ratio_NE_over_AVP"]) else None,
                RR_CI_low=rr_lo,
                RR_CI_high=rr_hi,
            )
        )

    return out


def run_g10_synthetic_tests() -> None:
    # Percentile CI sanity.
    lo, hi = g10_percentile_ci(range(100))
    assert abs(lo - 2.475) < 1e-9
    assert abs(hi - 96.525) < 1e-9

    # Cluster-copy semantics: a duplicated subject duplicates all its rows.
    toy = pd.DataFrame(
        {
            "subject_id": [1, 1, 2],
            "clone_id": ["a", "b", "c"],
            "interval_index": [0, 0, 0],
        }
    )
    draws = pd.DataFrame(
        {
            "subject_id": [1, 1, 2],
            "boot_copy": [0, 1, 2],
        }
    )
    z = g10_bootstrap_copy(toy, draws)
    assert len(z) == 5
    assert z["clone_id"].nunique() == 5

# ---------------------------------------------------------------------
# Synthetic tests
# ---------------------------------------------------------------------
def run_synthetic_tests() -> None:
    base = datetime(
        2026,
        1,
        1,
        0,
        0,
    )

    def row(
        itemid,
        rate,
        unit,
        st_min,
        en_min,
        status="ChangeDose/Rate",
    ):
        return dict(
            itemid=itemid,
            rate=rate,
            rateuom=unit,
            patientweight=70,
            starttime=base
            + timedelta(
                minutes=st_min
            ),
            endtime=base
            + timedelta(
                minutes=en_min
            ),
            statusdescription=status,
            ordercategorydescription="Continuous",
        )

    # Unit conversion.
    d, r, why = normalize_rate(
        AVP_ITEMID,
        0.03,
        "units/min",
        70,
    )
    assert d == "AVP"
    assert abs(r - 1.8) < 1e-12
    assert why == "accepted"

    # DOVSS-scale NE exact threshold.
    rows = [
        row(
            NE_ITEMID,
            0.30,
            "mcg/kg/min",
            0,
            120,
        ),
        row(
            NE_ITEMID,
            0.20,
            "mcg/kg/min",
            120,
            300,
        ),
        row(
            AVP_ITEMID,
            2.40,
            "units/hour",
            0,
            300,
        ),
    ]

    segs, _ = segments_from_rows(
        rows
    )
    tl, _ = build_atomic_timeline(
        segs
    )

    ev = cumulative_event_for_drug(
        tl,
        base + timedelta(
            minutes=60
        ),
        120,
        "NE",
        0.30,
        "dovss_scale_30m",
        base + timedelta(
            minutes=300
        ),
    )

    assert ev is not None
    assert ev.time == base + timedelta(
        minutes=120
    )

    # DOVSS-scale AVP exact 0.6 U/h threshold.
    rows = [
        row(
            NE_ITEMID,
            0.20,
            "mcg/kg/min",
            0,
            300,
        ),
        row(
            AVP_ITEMID,
            2.40,
            "units/hour",
            0,
            120,
        ),
        row(
            AVP_ITEMID,
            1.80,
            "units/hour",
            120,
            300,
        ),
    ]

    segs, _ = segments_from_rows(
        rows
    )
    tl, _ = build_atomic_timeline(
        segs
    )

    ev = cumulative_event_for_drug(
        tl,
        base + timedelta(
            minutes=60
        ),
        120,
        "AVP",
        2.40,
        "dovss_scale_30m",
        base + timedelta(
            minutes=300
        ),
    )

    assert ev is not None
    assert ev.time == base + timedelta(
        minutes=120
    )

    # Brief 20-min reduction is not sustained.
    rows = [
        row(
            NE_ITEMID,
            0.30,
            "mcg/kg/min",
            0,
            120,
        ),
        row(
            NE_ITEMID,
            0.20,
            "mcg/kg/min",
            120,
            140,
        ),
        row(
            NE_ITEMID,
            0.30,
            "mcg/kg/min",
            140,
            300,
        ),
        row(
            AVP_ITEMID,
            2.40,
            "units/hour",
            0,
            300,
        ),
    ]

    segs, _ = segments_from_rows(
        rows
    )
    tl, _ = build_atomic_timeline(
        segs
    )

    ev = cumulative_event_for_drug(
        tl,
        base + timedelta(
            minutes=60
        ),
        120,
        "NE",
        0.30,
        "dovss_scale_30m",
        base + timedelta(
            minutes=300
        ),
    )

    assert ev is None

    # Hypotension pair boundary.
    maps = [
        (
            base + timedelta(
                minutes=10
            ),
            64.0,
        ),
        (
            base + timedelta(
                minutes=25
            ),
            64.0,
        ),
    ]

    assert first_pair_lt65_time(
        maps,
        base,
        6,
    ) == base + timedelta(
        minutes=25
    )

    maps = [
        (
            base + timedelta(
                minutes=10
            ),
            64.0,
        ),
        (
            base + timedelta(
                minutes=71
            ),
            64.0,
        ),
    ]

    assert first_pair_lt65_time(
        maps,
        base,
        6,
    ) is None

    # Other vasoactive label classification.
    assert (
        classify_vasoactive_label(
            "Epinephrine"
        )
        == "other_pressor"
    )
    assert (
        classify_vasoactive_label(
            "Phenylephrine 80 mg/250 mL"
        )
        == "other_pressor"
    )
    assert (
        classify_vasoactive_label(
            "Dobutamine"
        )
        == "inotrope"
    )
    assert (
        classify_vasoactive_label(
            "Norepinephrine"
        )
        == "study_drug"
    )

    # Ridge logistic should fit a simple non-degenerate dataset.
    toy = pd.DataFrame(
        {
            "y": [
                1,
                1,
                1,
                0,
                0,
                0,
            ],
            "x": [
                -2,
                -1,
                0,
                1,
                2,
                3,
            ],
            "cat": [
                "a",
                "a",
                "b",
                "b",
                "a",
                "b",
            ],
        }
    )

    model = fit_logistic_ridge(
        toy,
        "y",
        ["x"],
        ["cat"],
    )

    p = predict_logistic(
        model,
        toy,
    )

    assert np.all(
        np.isfinite(p)
    )
    assert np.all(
        (p > 0)
        & (p < 1)
    )

    # Weighted KM sanity:
    km = pd.DataFrame(
        {
            "arm": [
                "NE_first",
                "NE_first",
                "AVP_first",
                "AVP_first",
            ],
            "interval_index": [
                0,
                0,
                0,
                0,
            ],
            "event": [
                1,
                0,
                0,
                0,
            ],
            "w": [
                1.0,
                1.0,
                1.0,
                1.0,
            ],
        }
    )

    risk = weighted_km_risk(
        km,
        "w",
    )

    rr = {
        x["arm"]: x["risk"]
        for x in risk
    }

    assert abs(
        rr["NE_first"] - 0.5
    ) < 1e-12
    assert abs(
        rr["AVP_first"] - 0.0
    ) < 1e-12


# ---------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------
def fetch_all(
    conn,
    name: str,
    sql: str,
    params=(),
) -> List[dict]:
    started = time.monotonic()
    log(f"▶ {name} ...")

    try:
        with conn.cursor() as cur:
            cur.execute(
                sql,
                params,
            )
            rows = cur.fetchall()

        log(
            f"  OK: {len(rows)} row(s), "
            f"{time.monotonic() - started:.1f}s"
        )

        return rows

    except Exception as exc:
        first = str(
            exc
        ).splitlines()[0]

        ERRORS.append(
            dict(
                stage=name,
                error=first,
            )
        )

        log(
            f"  FAIL: "
            f"{type(exc).__name__}: "
            f"{first}"
        )

        return []


def fetch_optional(
    conn,
    name: str,
    sql: str,
    params=(),
) -> List[dict]:
    started = time.monotonic()
    log(f"▶ {name} (optional) ...")

    try:
        with conn.cursor() as cur:
            cur.execute(
                sql,
                params,
            )
            rows = cur.fetchall()

        log(
            f"  OK: {len(rows)} row(s), "
            f"{time.monotonic() - started:.1f}s"
        )

        return rows

    except Exception as exc:
        first = str(
            exc
        ).splitlines()[0]

        warn(
            name,
            "failed:",
            f"{type(exc).__name__}:",
            first,
        )

        return []


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    run_synthetic_tests()
    run_g7_synthetic_tests()
    run_g8_synthetic_tests()
    run_g10_synthetic_tests()

    log(
        "=== G10 FROZEN-PROTOCOL CCW + SUBJECT-CLUSTER BOOTSTRAP ==="
    )
    log(
        "Synthetic tests: PASS (G6/G7/G8 inheritance + G10 bootstrap tests)"
    )

    if os.environ.get(
        "G10_SELFTEST_ONLY"
    ) == "1":
        log(
            "G10_SELFTEST_ONLY=1 -> "
            "database work skipped."
        )
        return

    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    import psycopg
    from psycopg.rows import dict_row

    conn = psycopg.connect(
        dbname=DBNAME,
        user=USER,
        host=SOCKET_DIR,
        port=PORT,
        row_factory=dict_row,
    )

    try:
        conn.execute(
            "SET default_transaction_read_only = on"
        )
        conn.execute(
            f"SET statement_timeout = "
            f"{int(TIMEOUT_MS)}"
        )

        # ----------------------------------------------------------
        # 1. Candidate stays and baseline demographics
        # ----------------------------------------------------------
        candidate_sql = """
        WITH ds AS (
          SELECT
            stay_id,
            BOOL_OR(itemid=%s) AS has_ne,
            BOOL_OR(itemid=%s) AS has_avp
          FROM mimiciv_icu.inputevents
          WHERE itemid IN (%s,%s)
            AND stay_id IS NOT NULL
            AND rate IS NOT NULL
            AND rate > 0
            AND endtime > starttime
            AND COALESCE(
                  lower(statusdescription),
                  ''
                )
                NOT IN (
                  'rewritten',
                  'flushed'
                )
            AND NOT (
              itemid=%s
              AND lower(
                    COALESCE(
                      statusdescription,
                      ''
                    )
                  )='bolus'
            )
          GROUP BY stay_id
        )
        SELECT
          ie.subject_id,
          ie.hadm_id,
          ie.stay_id,
          ie.intime,
          ie.outtime,
          ie.first_careunit,
          pat.gender,
          adm.race,
          adm.deathtime,
          adm.hospital_expire_flag,
          (
            pat.anchor_age
            + EXTRACT(
                YEAR
                FROM adm.admittime
              )
            - pat.anchor_year
          )::double precision AS age,
          s3.suspected_infection_time,
          s3.sofa_time,
          s3.sofa_score
        FROM ds
        JOIN mimiciv_icu.icustays ie
          USING (stay_id)
        JOIN mimiciv_hosp.admissions adm
          ON adm.hadm_id=ie.hadm_id
        JOIN mimiciv_hosp.patients pat
          ON pat.subject_id=ie.subject_id
        JOIN mimiciv_derived.sepsis3 s3
          ON s3.stay_id=ie.stay_id
         AND s3.sepsis3
        WHERE ds.has_ne
          AND ds.has_avp
        ORDER BY ie.stay_id
        """

        candidates = fetch_all(
            conn,
            "01_candidate_stays",
            candidate_sql,
            (
                NE_ITEMID,
                AVP_ITEMID,
                NE_ITEMID,
                AVP_ITEMID,
                NE_ITEMID,
            ),
        )

        if not candidates:
            raise RuntimeError(
                "No candidate stays returned."
            )

        stay_ids = [
            int(r["stay_id"])
            for r in candidates
        ]

        hadm_ids = sorted(
            {
                int(r["hadm_id"])
                for r in candidates
                if r.get(
                    "hadm_id"
                )
                is not None
            }
        )

        meta = {
            int(r["stay_id"]): r
            for r in candidates
        }

        # ----------------------------------------------------------
        # 2. NE / AVP inputevents
        # ----------------------------------------------------------
        infusion_sql = """
        SELECT
          stay_id,
          itemid,
          starttime,
          endtime,
          rate,
          rateuom,
          statusdescription,
          ordercategorydescription,
          patientweight
        FROM mimiciv_icu.inputevents
        WHERE stay_id = ANY(%s::bigint[])
          AND itemid IN (%s,%s)
        ORDER BY
          stay_id,
          starttime,
          endtime,
          itemid
        """

        inf_rows = fetch_all(
            conn,
            "02_NE_AVP_infusions",
            infusion_sql,
            (
                stay_ids,
                NE_ITEMID,
                AVP_ITEMID,
            ),
        )

        # ----------------------------------------------------------
        # 3. MAP + heart rate
        # ----------------------------------------------------------
        vital_sql = """
        SELECT
          stay_id,
          charttime,
          mbp,
          heart_rate
        FROM mimiciv_derived.vitalsign
        WHERE stay_id = ANY(%s::bigint[])
          AND (
            (
              mbp IS NOT NULL
              AND mbp > 0
              AND mbp < 300
            )
            OR (
              heart_rate IS NOT NULL
              AND heart_rate > 0
              AND heart_rate < 300
            )
          )
        ORDER BY stay_id, charttime
        """

        vital_rows = fetch_all(
            conn,
            "03_vitals",
            vital_sql,
            (stay_ids,),
        )

        # ----------------------------------------------------------
        # 4. Lactate
        # ----------------------------------------------------------
        lactate_sql = """
        SELECT
          hadm_id,
          charttime,
          valuenum AS lactate
        FROM mimiciv_hosp.labevents
        WHERE hadm_id = ANY(%s::bigint[])
          AND itemid = 50813
          AND valuenum IS NOT NULL
          AND valuenum > 0
          AND valuenum < 50
        ORDER BY hadm_id, charttime
        """

        lac_rows = fetch_all(
            conn,
            "04_lactate",
            lactate_sql,
            (hadm_ids,),
        )

        # ----------------------------------------------------------
        # 5. Discover other vasoactive itemids from d_items
        # ----------------------------------------------------------
        vaso_item_sql = """
        SELECT
          itemid,
          label,
          category
        FROM mimiciv_icu.d_items
        WHERE lower(COALESCE(label,'')) ~
          '(epinephrine|phenylephrine|dopamine|dobutamine|milrinone|angiotensin|methylene blue|isoproterenol|isoprenaline|norepinephrine|vasopressin)'
        ORDER BY itemid
        """

        vaso_items = fetch_optional(
            conn,
            "05_vasoactive_d_items",
            vaso_item_sql,
        )

        label_classes = {
            int(r["itemid"]):
            classify_vasoactive_label(
                r["label"]
            )
            for r in vaso_items
        }

        other_itemids = sorted(
            [
                itemid
                for itemid, cls
                in label_classes.items()
                if cls in {
                    "other_pressor",
                    "inotrope",
                }
            ]
        )

        if other_itemids:
            other_vaso_sql = """
            SELECT
              stay_id,
              itemid,
              starttime,
              endtime,
              rate,
              rateuom,
              statusdescription,
              ordercategorydescription
            FROM mimiciv_icu.inputevents
            WHERE stay_id = ANY(%s::bigint[])
              AND itemid = ANY(%s::int[])
              AND rate IS NOT NULL
              AND rate > 0
              AND endtime > starttime
              AND COALESCE(
                    lower(statusdescription),
                    ''
                  )
                  NOT IN (
                    'rewritten',
                    'flushed'
                  )
            ORDER BY
              stay_id,
              starttime,
              endtime,
              itemid
            """

            other_vaso_rows = fetch_optional(
                conn,
                "06_other_vasoactives",
                other_vaso_sql,
                (
                    stay_ids,
                    other_itemids,
                ),
            )
        else:
            warn(
                "No non-study vasoactive itemids "
                "were discovered from d_items."
            )
            other_vaso_rows = []

        # ----------------------------------------------------------
        # 6. Optional systemic steroid orders
        # ----------------------------------------------------------
        steroid_sql = """
        SELECT
          hadm_id,
          starttime,
          stoptime,
          drug
        FROM mimiciv_hosp.prescriptions
        WHERE hadm_id = ANY(%s::bigint[])
          AND lower(COALESCE(drug,'')) ~
            '(hydrocortisone|methylprednisolone|prednisone|dexamethasone)'
        ORDER BY hadm_id, starttime
        """

        steroid_rows = fetch_optional(
            conn,
            "07_systemic_steroid_orders",
            steroid_sql,
            (hadm_ids,),
        )

        # ----------------------------------------------------------
        # Organize rows
        # ----------------------------------------------------------
        by_stay_inf = defaultdict(list)
        for r in inf_rows:
            by_stay_inf[
                int(r["stay_id"])
            ].append(r)

        by_stay_map = defaultdict(list)
        by_stay_hr = defaultdict(list)

        for r in vital_rows:
            sid = int(
                r["stay_id"]
            )

            if r.get(
                "mbp"
            ) is not None:
                by_stay_map[
                    sid
                ].append(
                    (
                        r["charttime"],
                        float(
                            r["mbp"]
                        ),
                    )
                )

            if r.get(
                "heart_rate"
            ) is not None:
                by_stay_hr[
                    sid
                ].append(
                    (
                        r["charttime"],
                        float(
                            r[
                                "heart_rate"
                            ]
                        ),
                    )
                )

        by_hadm_lac = defaultdict(list)
        for r in lac_rows:
            by_hadm_lac[
                int(r["hadm_id"])
            ].append(
                (
                    r["charttime"],
                    float(
                        r["lactate"]
                    ),
                )
            )

        by_stay_other = defaultdict(list)
        for r in other_vaso_rows:
            by_stay_other[
                int(r["stay_id"])
            ].append(r)

        by_hadm_steroid = defaultdict(list)
        for r in steroid_rows:
            by_hadm_steroid[
                int(r["hadm_id"])
            ].append(r)

        # ----------------------------------------------------------
        # Reconstruct NE/AVP timelines
        # ----------------------------------------------------------
        timelines = {}
        qc_all = Counter()

        for i, sid in enumerate(
            stay_ids,
            1,
        ):
            segs, qc = segments_from_rows(
                by_stay_inf.get(
                    sid,
                    [],
                )
            )

            tl, qc2 = build_atomic_timeline(
                segs
            )

            timelines[
                sid
            ] = tl

            qc_all.update(
                qc
            )
            qc_all.update(
                qc2
            )

            if i % 500 == 0:
                log(
                    "  reconstructed stays:",
                    i,
                )

        qc_rows = [
            dict(
                metric=k,
                count=int(v),
            )
            for k, v in sorted(
                qc_all.items()
            )
        ]

        qc_rows += [
            dict(
                metric="candidate_stays",
                count=len(
                    stay_ids
                ),
            ),
            dict(
                metric="NE_AVP_raw_rows",
                count=len(
                    inf_rows
                ),
            ),
            dict(
                metric="vital_rows",
                count=len(
                    vital_rows
                ),
            ),
            dict(
                metric="lactate_rows",
                count=len(
                    lac_rows
                ),
            ),
            dict(
                metric="discovered_other_vasoactive_itemids",
                count=len(
                    other_itemids
                ),
            ),
            dict(
                metric="other_vasoactive_rows",
                count=len(
                    other_vaso_rows
                ),
            ),
            dict(
                metric="steroid_order_rows",
                count=len(
                    steroid_rows
                ),
            ),
        ]

        save_csv(
            "00_data_QC.csv",
            qc_rows,
        )

        # ----------------------------------------------------------
        # Frozen G10 primary decision point (selected a priori in G9)
        # ----------------------------------------------------------
        final_records = []
        stage_counter = Counter()

        for sid in stay_ids:
            m = meta[sid]
            tl = timelines[sid]

            t0, stage = g10_find_final_t0(
                tl,
                by_stay_map.get(sid, []),
                by_stay_other.get(sid, []),
                label_classes,
                m,
            )

            for k, v in stage.items():
                if v:
                    stage_counter[k] += 1

            if t0 is None:
                continue

            hadm_id = int(m["hadm_id"])
            eps = timedelta(microseconds=1)

            lac_near_icu = [
                (tt, v)
                for tt, v in by_hadm_lac.get(hadm_id, [])
                if (
                    m["intime"] - timedelta(hours=24)
                    <= tt
                    <= m["outtime"]
                )
            ]

            map_mean, map_minimum, _ = recent_summary(
                by_stay_map.get(sid, []),
                t0 + timedelta(microseconds=1),
                120,
            )

            hr_mean, _hr_min, _hr_latest = recent_summary(
                by_stay_hr.get(sid, []),
                t0 + timedelta(microseconds=1),
                120,
            )

            baseline_lac = latest_lactate(
                lac_near_icu,
                t0 + timedelta(microseconds=1),
                24,
            )

            other_rows_this = by_stay_other.get(sid, [])
            steroid_rows_this = by_hadm_steroid.get(hadm_id, [])

            first_dual = tl.first_dual_time()

            ne0 = tl.rate_at(t0 - eps, "NE")
            avp0 = tl.rate_at(t0 - eps, "AVP")

            final_records.append(
                dict(
                    sid=sid,
                    subject_id=int(m["subject_id"]),
                    hadm_id=hadm_id,
                    t0=t0,
                    outtime=m["outtime"],
                    deathtime=m.get("deathtime"),
                    age=float(m["age"]),
                    gender=str(m.get("gender") or "Missing"),
                    race_group=race_group(m.get("race")),
                    careunit=str(
                        m.get("first_careunit") or "Unknown"
                    ),
                    sofa_score=(
                        float(m["sofa_score"])
                        if m.get("sofa_score") is not None
                        else None
                    ),
                    ne0=float(ne0),
                    avp0=float(avp0),
                    delay_h=(
                        (t0 - first_dual).total_seconds() / 3600.0
                        if first_dual is not None
                        else None
                    ),
                    baseline_lactate=baseline_lac,
                    baseline_map_mean=map_mean,
                    baseline_map_min=map_minimum,
                    baseline_hr_mean=hr_mean,
                    steroid_prior24h=int(
                        steroid_prior24h(
                            steroid_rows_this,
                            t0,
                        )
                    ),
                    other_pressor_at_t0=int(
                        interval_active(
                            other_rows_this,
                            t0,
                            ["other_pressor"],
                            label_classes,
                        )
                    ),
                    other_pressor_prior2h=int(
                        g10_other_pressor_overlap(
                            other_rows_this,
                            t0 - timedelta(hours=2),
                            t0,
                            label_classes,
                        )
                    ),
                    inotrope_at_t0=int(
                        interval_active(
                            other_rows_this,
                            t0,
                            ["inotrope"],
                            label_classes,
                        )
                    ),
                    exact_transition_at_t0=int(
                        g10_exact_transition(
                            tl,
                            t0,
                        )
                    ),
                )
            )

        if len(final_records) != EXPECTED_FINAL_ELIGIBLE_N:
            raise RuntimeError(
                "Frozen G9 cohort reconciliation failed: "
                f"expected {EXPECTED_FINAL_ELIGIBLE_N}, "
                f"observed {len(final_records)}."
            )

        if any(r["exact_transition_at_t0"] for r in final_records):
            raise RuntimeError(
                "Frozen cohort unexpectedly contains a rate transition at t0."
            )

        if any(
            r["other_pressor_at_t0"]
            or r["other_pressor_prior2h"]
            for r in final_records
        ):
            raise RuntimeError(
                "Frozen cohort unexpectedly contains non-study pressor exposure."
            )

        # Collapse sparse care units after the cohort is frozen.
        care_map = collapse_rare(
            [r["careunit"] for r in final_records],
            min_count=30,
        )
        for r in final_records:
            r["careunit_group"] = care_map[r["careunit"]]

        ccw_records = final_records

        # Cohort / clustering audit.
        subject_counts = Counter(
            r["subject_id"]
            for r in ccw_records
        )
        n_subjects = len(subject_counts)
        repeated_subjects = sum(
            n > 1
            for n in subject_counts.values()
        )
        max_stays_per_subject = max(
            subject_counts.values()
        ) if subject_counts else 0

        save_csv(
            "01_frozen_protocol_population.csv",
            [
                dict(
                    eligible_stays=len(ccw_records),
                    unique_subjects=n_subjects,
                    subjects_with_multiple_eligible_stays=repeated_subjects,
                    max_eligible_stays_per_subject=max_stays_per_subject,
                    median_hours_first_overlap_to_t0=median_or_none(
                        r["delay_h"] for r in ccw_records
                    ),
                    median_NE_t0=median_or_none(
                        r["ne0"] for r in ccw_records
                    ),
                    median_AVP_t0=median_or_none(
                        r["avp0"] for r in ccw_records
                    ),
                    exact_transition_at_t0=sum(
                        r["exact_transition_at_t0"]
                        for r in ccw_records
                    ),
                    other_pressor_active_at_t0=sum(
                        r["other_pressor_at_t0"]
                        for r in ccw_records
                    ),
                    other_pressor_prior2h=sum(
                        r["other_pressor_prior2h"]
                        for r in ccw_records
                    ),
                )
            ],
        )

        save_csv(
            "02_frozen_protocol_funnel.csv",
            [
                dict(
                    candidate_stays=len(stay_ids),
                    two_hour_dual=stage_counter["two_hour_dual"],
                    sepsis_before_t=stage_counter["sepsis_before_t"],
                    adult_ok=stage_counter["adult_ok"],
                    map_stable=stage_counter["map_stable"],
                    dose_ok=stage_counter["dose_ok"],
                    no_escalation=stage_counter["no_escalation"],
                    boundary_safe=stage_counter["boundary_safe"],
                    dovss_stable120=stage_counter["dovss_stable120"],
                    other_pressor_clean120=stage_counter[
                        "other_pressor_clean120"
                    ],
                    eligible_t0=stage_counter["eligible"],
                )
            ],
        )

        save_csv(
            "03_natural_censoring_audit.csv",
            natural_censor_audit_rows(ccw_records),
        )

        # ----------------------------------------------------------
        # Run repaired CCW estimators
        # ----------------------------------------------------------
        observed_action_rows = []
        weight_model_rows = []
        weight_dist_rows = []
        riskset_balance_all = []
        completion_balance_all = []
        risk_rows = []
        estimator_agreement_rows = []
        outcome_model_rows = []
        positivity_rows_all = []
        action_timing_rows_all = []
        outcome_timing_rows_all = []
        bootstrap_weight_bases = {}
        bootstrap_analysis_bases = {}
        primary_point_rows_for_bootstrap = []

        for spec in SPECS:
            log(
                "▶ G10 CCW spec:",
                spec["name"],
            )

            # ------------------------------------------------------
            # Observed action counts
            # ------------------------------------------------------
            cc = Counter()

            for rec in ccw_records:
                ne, avp = get_spec_events(
                    timelines[rec["sid"]],
                    rec,
                    spec,
                )
                cc[classify_events(ne, avp)] += 1

            clear = cc["NE_first"] + cc["AVP_first"]

            observed_action_rows.append(
                dict(
                    spec_name=spec["name"],
                    mode=spec["mode"],
                    grace_minutes=spec["grace_minutes"],
                    eligible_t0=len(ccw_records),
                    NE_first=cc["NE_first"],
                    AVP_first=cc["AVP_first"],
                    near_simultaneous=cc["near_simultaneous"],
                    no_clear_action=cc["no_clear_action"],
                    clear_actions=clear,
                    clear_action_fraction=(
                        clear / len(ccw_records)
                        if ccw_records
                        else None
                    ),
                    AVP_share_among_clear=(
                        cc["AVP_first"] / clear
                        if clear
                        else None
                    ),
                )
            )

            timing_rows, outcome_timing_rows = action_timing_and_outcome_rows(
                ccw_records,
                timelines,
                spec,
                by_stay_map,
            )

            if spec["name"] == "primary_dovss120_immediate":
                observed_actions = {
                    "NE_first": int(cc["NE_first"]),
                    "AVP_first": int(cc["AVP_first"]),
                    "near_simultaneous": int(cc["near_simultaneous"]),
                    "no_clear_action": int(cc["no_clear_action"]),
                }
                if observed_actions != EXPECTED_FINAL_PRIMARY_ACTIONS:
                    raise RuntimeError(
                        "Frozen primary-action reconciliation failed: "
                        f"expected {EXPECTED_FINAL_PRIMARY_ACTIONS}, "
                        f"observed {observed_actions}."
                    )
            action_timing_rows_all.extend(timing_rows)
            outcome_timing_rows_all.extend(outcome_timing_rows)

            # G8 intentionally changes the action definition, so the historical
            # G6/G7 action counts are not a hard gate.
            if False:
                observed = {
                    k: int(cc[k])
                    for k in EXPECTED_G6_PRIMARY_ACTIONS
                }
                if observed != EXPECTED_G6_PRIMARY_ACTIONS:
                    raise RuntimeError(
                        "G6 primary-action reconciliation failed: "
                        f"expected {EXPECTED_G6_PRIMARY_ACTIONS}, "
                        f"observed {observed}."
                    )

            # ------------------------------------------------------
            # Censoring models and dual weights
            # ------------------------------------------------------
            weight_df, clone_protocol = build_weight_rows_for_spec(
                ccw_records,
                timelines,
                spec,
                by_stay_map,
                by_stay_hr,
                by_hadm_lac,
                by_stay_other,
                label_classes,
                by_hadm_steroid,
                outcome_h=6,
            )

            if weight_df.empty:
                raise RuntimeError(
                    f"{spec['name']}: empty weight dataset."
                )

            # Cache the deterministic person-period data before fitting;
            # bootstrap replicates will resample subjects and refit the censor models.
            bootstrap_weight_bases[spec["name"]] = weight_df.copy()

            weight_df, diag = fit_g7_censor_weights(
                weight_df,
                spec["name"],
            )

            weight_model_rows.extend(diag)

            positivity_rows_all.extend(
                positivity_relevant_rows(
                    weight_df,
                    spec["name"],
                )
            )

            weight_dist_rows.extend(
                g7_weight_distribution_rows(
                    weight_df,
                    clone_protocol,
                    spec["name"],
                )
            )

            # Truncated marginal weights are used only for sensitivity estimation;
            # raw marginal weights are used for balance diagnostics.
            weight_df["sw_marginal_start_trunc"] = truncate_weights(
                weight_df.rename(
                    columns={
                        "sw_marginal_start": "_g7_tmp_weight"
                    }
                ),
                "_g7_tmp_weight",
            ).to_numpy()

            marginal_paths = clone_weight_paths_g7(
                weight_df,
                "marginal",
            )
            baseline_paths = clone_weight_paths_g7(
                weight_df,
                "baseline",
            )

            riskset_balance_all.extend(
                riskset_balance_rows(
                    weight_df,
                    spec["name"],
                    "sw_marginal_start",
                )
            )

            completion_balance_all.extend(
                completion_balance_g7(
                    ccw_records,
                    clone_protocol,
                    marginal_paths,
                    spec["name"],
                )
            )

            # ------------------------------------------------------
            # Outcome estimation at 6 / 12 / 24 h
            # ------------------------------------------------------
            for horizon_h in OUTCOME_WINDOWS_H:
                # Estimator B: marginal/time-only stabilized IPCW + weighted KM.
                analysis_m = build_analysis_rows(
                    ccw_records,
                    timelines,
                    spec,
                    clone_protocol,
                    marginal_paths,
                    by_stay_map,
                    horizon_h,
                )

                if analysis_m.empty:
                    raise RuntimeError(
                        f"{spec['name']}/{horizon_h}h: empty marginal analysis."
                    )

                bootstrap_analysis_bases.setdefault(
                    spec["name"],
                    {},
                )[horizon_h] = analysis_m.copy()

                analysis_m["weight_marginal_raw"] = analysis_m["weight_raw"]
                analysis_m["weight_marginal_trunc"] = truncate_weights(
                    analysis_m,
                    "weight_marginal_raw",
                )

                for col, version in (
                    ("weight_marginal_raw", "raw"),
                    ("weight_marginal_trunc", "truncated_1_99"),
                ):
                    km = weighted_km_risk(
                        analysis_m,
                        col,
                    )
                    by_arm = {
                        x["arm"]: x
                        for x in km
                    }

                    ne_risk = by_arm["NE_first"]["risk"]
                    avp_risk = by_arm["AVP_first"]["risk"]

                    risk_rows.append(
                        dict(
                            spec_name=spec["name"],
                            grace_minutes=spec["grace_minutes"],
                            outcome="two_MAP_lt65_15_to_60min_pair",
                            horizon_hours=horizon_h,
                            estimator="marginal_time_only_IPCW_weighted_KM",
                            weight_version=version,
                            risk_NE_first=ne_risk,
                            risk_AVP_first=avp_risk,
                            risk_difference_NE_minus_AVP=(ne_risk - avp_risk),
                            risk_ratio_NE_over_AVP=(
                                ne_risk / avp_risk
                                if avp_risk > 0
                                else None
                            ),
                            observed_events_NE=by_arm["NE_first"]["observed_events"],
                            observed_events_AVP=by_arm["AVP_first"]["observed_events"],
                        )
                    )

                    if horizon_h == 6:
                        primary_point_rows_for_bootstrap.append(
                            dict(
                                spec_name=spec["name"],
                                horizon_hours=6,
                                weight_version=(
                                    "untruncated"
                                    if version == "raw"
                                    else "truncated_1_99"
                                ),
                                risk_NE_first=ne_risk,
                                risk_AVP_first=avp_risk,
                                risk_difference_NE_minus_AVP=(
                                    ne_risk - avp_risk
                                ),
                                risk_ratio_NE_over_AVP=(
                                    ne_risk / avp_risk
                                    if avp_risk > 0
                                    else None
                                ),
                            )
                        )

                # Estimator A: baseline-stabilized IPCW + weighted outcome model
                # + standardization over original uncloned baseline cohort.
                analysis_b = build_analysis_rows(
                    ccw_records,
                    timelines,
                    spec,
                    clone_protocol,
                    baseline_paths,
                    by_stay_map,
                    horizon_h,
                )

                analysis_b["weight_baseline_raw"] = analysis_b["weight_raw"]
                analysis_b["weight_baseline_trunc"] = truncate_weights(
                    analysis_b,
                    "weight_baseline_raw",
                )

                for col, version in (
                    ("weight_baseline_raw", "raw"),
                    ("weight_baseline_trunc", "truncated_1_99"),
                ):
                    std, omdiag = standardize_outcome_risk(
                        analysis_b,
                        ccw_records,
                        col,
                        horizon_h,
                    )

                    omdiag_row = dict(
                        spec_name=spec["name"],
                        horizon_hours=horizon_h,
                        weight_version=version,
                    )
                    omdiag_row.update(omdiag)
                    outcome_model_rows.append(omdiag_row)

                    by_arm = {
                        x["arm"]: x
                        for x in std
                    }

                    ne_risk = by_arm["NE_first"]["risk"]
                    avp_risk = by_arm["AVP_first"]["risk"]

                    risk_rows.append(
                        dict(
                            spec_name=spec["name"],
                            grace_minutes=spec["grace_minutes"],
                            outcome="two_MAP_lt65_15_to_60min_pair",
                            horizon_hours=horizon_h,
                            estimator="baseline_stabilized_IPCW_standardized_outcome",
                            weight_version=version,
                            risk_NE_first=ne_risk,
                            risk_AVP_first=avp_risk,
                            risk_difference_NE_minus_AVP=(ne_risk - avp_risk),
                            risk_ratio_NE_over_AVP=(
                                ne_risk / avp_risk
                                if avp_risk > 0
                                else None
                            ),
                            observed_events_NE=None,
                            observed_events_AVP=None,
                        )
                    )

                # Agreement metric uses truncated estimates, which are the
                # pre-specified robust point estimates for protocol freezing.
                current = [
                    r for r in risk_rows
                    if (
                        r["spec_name"] == spec["name"]
                        and r["horizon_hours"] == horizon_h
                        and r["weight_version"] == "truncated_1_99"
                    )
                ]

                est_map = {
                    r["estimator"]: r
                    for r in current
                }

                km_rd = est_map[
                    "marginal_time_only_IPCW_weighted_KM"
                ]["risk_difference_NE_minus_AVP"]

                std_rd = est_map[
                    "baseline_stabilized_IPCW_standardized_outcome"
                ]["risk_difference_NE_minus_AVP"]

                estimator_agreement_rows.append(
                    dict(
                        spec_name=spec["name"],
                        horizon_hours=horizon_h,
                        marginal_KM_RD=km_rd,
                        standardized_outcome_RD=std_rd,
                        absolute_RD_difference=abs(km_rd - std_rd),
                        same_direction=int(
                            (km_rd == 0)
                            or (std_rd == 0)
                            or (math.copysign(1.0, km_rd) == math.copysign(1.0, std_rd))
                        ),
                    )
                )

        save_csv(
            "05_observed_actions_in_CCW_cohort.csv",
            observed_action_rows,
        )
        save_csv(
            "06_censor_weight_model_diagnostics.csv",
            weight_model_rows,
        )
        save_csv(
            "07_weight_distributions_dual_estimators.csv",
            weight_dist_rows,
        )
        save_csv(
            "08_marginal_weight_riskset_balance.csv",
            riskset_balance_all,
        )
        save_csv(
            "09_marginal_weight_completion_balance.csv",
            completion_balance_all,
        )
        save_csv(
            "10_point_estimates_dual_estimators.csv",
            risk_rows,
        )
        save_csv(
            "11_outcome_model_diagnostics.csv",
            outcome_model_rows,
        )
        save_csv(
            "12_estimator_agreement.csv",
            estimator_agreement_rows,
        )

        save_csv(
            "13_action_timing.csv",
            action_timing_rows_all,
        )
        save_csv(
            "14_outcome_timing_relative_to_action.csv",
            outcome_timing_rows_all,
        )
        save_csv(
            "15_positivity_relevant_rows.csv",
            positivity_rows_all,
        )

        # ----------------------------------------------------------
        # Subject-cluster bootstrap for PRIMARY 6-h inference
        # ----------------------------------------------------------
        log(
            "▶ subject-cluster bootstrap:",
            BOOTSTRAP_REPS,
            "replicates",
        )

        subjects = np.asarray(
            sorted(
                {
                    int(r["subject_id"])
                    for r in ccw_records
                }
            ),
            dtype=np.int64,
        )

        if len(subjects) == 0:
            raise RuntimeError("No subjects available for bootstrap.")

        rng = np.random.default_rng(
            BOOTSTRAP_SEED
        )

        bootstrap_replicates = []
        bootstrap_diag = []

        for b in range(1, BOOTSTRAP_REPS + 1):
            sampled = rng.choice(
                subjects,
                size=len(subjects),
                replace=True,
            )

            draw_map = pd.DataFrame(
                {
                    "subject_id": sampled,
                    "boot_copy": np.arange(
                        len(sampled),
                        dtype=np.int64,
                    ),
                }
            )

            rep_success = True
            rep_errors = []

            for spec in SPECS:
                name = spec["name"]

                try:
                    est, diag = g10_bootstrap_one_spec(
                        bootstrap_weight_bases[name],
                        {
                            6: bootstrap_analysis_bases[name][6]
                        },
                        draw_map,
                        name,
                    )

                    for r in est:
                        r["bootstrap_rep"] = b
                        r["success"] = 1
                        r["error"] = None
                        bootstrap_replicates.append(r)

                    bootstrap_diag.append(
                        dict(
                            bootstrap_rep=b,
                            spec_name=name,
                            success=1,
                            censor_models_converged=diag[
                                "censor_models_converged"
                            ],
                            min_p01_denominator_probability=diag[
                                "min_p01_denominator_probability"
                            ],
                            error=None,
                        )
                    )

                    if not diag["censor_models_converged"]:
                        rep_success = False
                        rep_errors.append(
                            f"{name}:model_nonconvergence"
                        )

                except Exception as exc:
                    rep_success = False
                    msg = (
                        f"{type(exc).__name__}: "
                        f"{str(exc).splitlines()[0]}"
                    )
                    rep_errors.append(
                        f"{name}:{msg}"
                    )

                    # Add explicit failed rows for both reported weight versions.
                    for version in (
                        "untruncated",
                        "truncated_1_99",
                    ):
                        bootstrap_replicates.append(
                            dict(
                                bootstrap_rep=b,
                                spec_name=name,
                                horizon_hours=6,
                                weight_version=version,
                                risk_NE_first=None,
                                risk_AVP_first=None,
                                risk_difference_NE_minus_AVP=None,
                                risk_ratio_NE_over_AVP=None,
                                success=0,
                                error=msg,
                            )
                        )

                    bootstrap_diag.append(
                        dict(
                            bootstrap_rep=b,
                            spec_name=name,
                            success=0,
                            censor_models_converged=0,
                            min_p01_denominator_probability=None,
                            error=msg,
                        )
                    )

            if (
                b == 1
                or b % 25 == 0
                or b == BOOTSTRAP_REPS
            ):
                log(
                    "  bootstrap",
                    b,
                    "/",
                    BOOTSTRAP_REPS,
                    "completed",
                    "OK" if rep_success else (
                        "with issue(s): " + "; ".join(rep_errors)
                    ),
                )

        save_csv(
            "20_bootstrap_replicates_6h.csv",
            bootstrap_replicates,
        )
        save_csv(
            "21_bootstrap_model_diagnostics.csv",
            bootstrap_diag,
        )

        bootstrap_summary = g10_bootstrap_summary(
            primary_point_rows_for_bootstrap,
            bootstrap_replicates,
            BOOTSTRAP_REPS,
        )

        for r in bootstrap_summary:
            r["minimum_success_fraction_required"] = (
                BOOTSTRAP_MIN_SUCCESS_FRACTION
            )
            r["bootstrap_quality_flag"] = (
                "PASS"
                if (
                    r["bootstrap_success_fraction"]
                    is not None
                    and r["bootstrap_success_fraction"]
                    >= BOOTSTRAP_MIN_SUCCESS_FRACTION
                )
                else "FAIL"
            )

        save_csv(
            "22_bootstrap_summary_6h.csv",
            bootstrap_summary,
        )

        # ----------------------------------------------------------
        # Legacy diagnostic verdicts remain below as supplementary audit.
        # ----------------------------------------------------------
        # ----------------------------------------------------------
        # Legacy G7-style verdict: retained for audit, not decisive in G8
        # ----------------------------------------------------------
        diag_df = pd.DataFrame(weight_model_rows)
        wd_df = pd.DataFrame(weight_dist_rows)
        rb_df = pd.DataFrame(riskset_balance_all)
        cb_df = pd.DataFrame(completion_balance_all)
        ag_df = pd.DataFrame(estimator_agreement_rows)
        act_df = pd.DataFrame(observed_action_rows)
        om_df = pd.DataFrame(outcome_model_rows)

        verdict_rows = []

        for spec in SPECS:
            name = spec["name"]

            dd = diag_df[diag_df["spec_name"] == name]
            ww = wd_df[
                (wd_df["spec_name"] == name)
                & (wd_df["weight_scheme"] == "marginal")
            ]
            rr = rb_df[rb_df["spec_name"] == name]
            cc = cb_df[cb_df["spec_name"] == name]
            aa = act_df[act_df["spec_name"] == name]
            gg = ag_df[
                (ag_df["spec_name"] == name)
                & (ag_df["horizon_hours"] == 6)
            ]
            oo = om_df[om_df["spec_name"] == name]

            censor_models_converged = bool(
                len(dd) == 2
                and (dd["baseline_num_converged"] == 1).all()
                and (dd["marginal_num_converged"] == 1).all()
                and (dd["denominator_converged"] == 1).all()
            )

            outcome_models_converged = bool(
                len(oo) > 0
                and (oo["outcome_model_converged"] == 1).all()
            )

            min_p01 = (
                float(dd["p01_raw_p_den"].min())
                if len(dd)
                else None
            )

            comp_ww = ww[ww["scope"] == "strategy_completers_final"]
            max_p99 = (
                float(comp_ww["p99_weight"].max())
                if len(comp_ww)
                else None
            )
            max_weight = (
                float(comp_ww["max_weight"].max())
                if len(comp_ww)
                else None
            )

            # Ignore interval 0: cloning guarantees exact baseline balance there.
            rr_after0 = rr[rr["interval_index"] > 0]
            max_riskset_smd = (
                float(rr_after0["abs_weighted_SMD"].max())
                if len(rr_after0)
                else None
            )

            max_completion_smd = (
                float(cc["abs_marginal_weighted_SMD"].max())
                if len(cc)
                else None
            )

            min_arm = (
                int(
                    min(
                        aa.iloc[0]["NE_first"],
                        aa.iloc[0]["AVP_first"],
                    )
                )
                if len(aa)
                else None
            )

            clear_fraction = (
                float(aa.iloc[0]["clear_action_fraction"])
                if len(aa)
                else None
            )

            estimator_rd_gap_6h = (
                float(gg.iloc[0]["absolute_RD_difference"])
                if len(gg)
                else None
            )

            same_direction_6h = (
                int(gg.iloc[0]["same_direction"])
                if len(gg)
                else 0
            )

            flags = []

            if not censor_models_converged:
                flags.append("censor_model_nonconvergence")

            if not outcome_models_converged:
                flags.append("outcome_model_nonconvergence")

            if min_p01 is not None and min_p01 < 0.05:
                flags.append("p01_denominator_probability_lt0.05")

            if max_p99 is not None and max_p99 > 10:
                flags.append("marginal_completion_p99_weight_gt10")

            if max_weight is not None and max_weight > 25:
                flags.append("marginal_completion_max_weight_gt25")

            if max_riskset_smd is not None and max_riskset_smd > 0.10:
                flags.append("marginal_riskset_abs_SMD_gt0.10")

            # Completion balance is informative but less decisive than risk-set balance.
            if max_completion_smd is not None and max_completion_smd > 0.20:
                flags.append("marginal_completion_abs_SMD_gt0.20")

            if min_arm is not None and min_arm < 100:
                flags.append("observed_action_min_arm_lt100")

            if clear_fraction is not None and clear_fraction < 0.10:
                flags.append("clear_action_fraction_lt0.10")

            if estimator_rd_gap_6h is not None and estimator_rd_gap_6h > 0.03:
                flags.append("dual_estimators_6h_RD_gap_gt0.03")

            if same_direction_6h != 1:
                flags.append("dual_estimators_6h_direction_disagrees")

            verdict_rows.append(
                dict(
                    spec_name=name,
                    censor_models_converged=int(censor_models_converged),
                    outcome_models_converged=int(outcome_models_converged),
                    min_p01_denominator_probability=min_p01,
                    marginal_completion_p99_weight=max_p99,
                    marginal_completion_max_weight=max_weight,
                    max_marginal_riskset_abs_SMD=max_riskset_smd,
                    max_marginal_completion_abs_SMD=max_completion_smd,
                    observed_action_min_arm=min_arm,
                    clear_action_fraction=clear_fraction,
                    dual_estimators_6h_RD_gap=estimator_rd_gap_6h,
                    dual_estimators_6h_same_direction=same_direction_6h,
                    diagnostic_flags=(";".join(flags) if flags else ""),
                    protocol_freeze_status=(
                        "READY_FOR_BOOTSTRAP_PROTOCOL_FREEZE"
                        if not flags
                        else "REVIEW_BEFORE_BOOTSTRAP"
                    ),
                )
            )

        save_csv(
            "16_legacy_G7_style_verdict_not_decisive.csv",
            verdict_rows,
        )

        # ----------------------------------------------------------
        # Corrected G8 protocol-freeze verdict
        # ----------------------------------------------------------
        pos_df = pd.DataFrame(positivity_rows_all)
        risk_df = pd.DataFrame(risk_rows)

        corrected_rows = []

        for spec in SPECS:
            name = spec["name"]

            dd = diag_df[diag_df["spec_name"] == name]
            ww = wd_df[
                (wd_df["spec_name"] == name)
                & (wd_df["weight_scheme"] == "marginal")
            ]
            rr = rb_df[rb_df["spec_name"] == name]
            aa = act_df[act_df["spec_name"] == name]
            gg = ag_df[
                (ag_df["spec_name"] == name)
                & (ag_df["horizon_hours"] == 6)
            ]
            oo = om_df[om_df["spec_name"] == name]
            pp = pos_df[
                (pos_df["spec_name"] == name)
                & (pos_df["scope"] == "weight_update_rows")
            ]

            censor_models_converged = bool(
                len(dd) == 2
                and (dd["baseline_num_converged"] == 1).all()
                and (dd["marginal_num_converged"] == 1).all()
                and (dd["denominator_converged"] == 1).all()
            )

            outcome_models_converged = bool(
                len(oo) > 0
                and (oo["outcome_model_converged"] == 1).all()
            )

            min_update_p01 = (
                float(pp["p_den_p01"].min())
                if len(pp)
                else None
            )
            max_pct_update_p_lt001 = (
                float(pp["pct_p_den_lt_0_01"].max())
                if len(pp)
                else None
            )
            max_update_ratio_p99 = (
                float(pp["marginal_ratio_p99"].max())
                if len(pp)
                else None
            )

            ppw = ww[ww["scope"] == "person_period_start"]
            compw = ww[ww["scope"] == "strategy_completers_final"]

            max_pp_p99_weight = (
                float(ppw["p99_weight"].max())
                if len(ppw)
                else None
            )
            max_pp_weight = (
                float(ppw["max_weight"].max())
                if len(ppw)
                else None
            )
            max_completion_p99 = (
                float(compw["p99_weight"].max())
                if len(compw)
                else None
            )

            rr_after0 = rr[rr["interval_index"] > 0]
            max_riskset_smd = (
                float(rr_after0["abs_weighted_SMD"].max())
                if len(rr_after0)
                else None
            )

            min_arm = (
                int(min(aa.iloc[0]["NE_first"], aa.iloc[0]["AVP_first"]))
                if len(aa)
                else None
            )
            clear_fraction = (
                float(aa.iloc[0]["clear_action_fraction"])
                if len(aa)
                else None
            )

            rd_gap = (
                float(gg.iloc[0]["absolute_RD_difference"])
                if len(gg)
                else None
            )
            same_direction = (
                int(gg.iloc[0]["same_direction"])
                if len(gg)
                else 0
            )

            six = risk_df[
                (risk_df["spec_name"] == name)
                & (risk_df["horizon_hours"] == 6)
            ]

            def trunc_gap(estimator_name: str):
                z = six[six["estimator"] == estimator_name]
                raw = z[z["weight_version"] == "raw"]
                tr = z[z["weight_version"] == "truncated_1_99"]
                if len(raw) != 1 or len(tr) != 1:
                    return None
                return abs(
                    float(raw.iloc[0]["risk_difference_NE_minus_AVP"])
                    - float(tr.iloc[0]["risk_difference_NE_minus_AVP"])
                )

            marginal_trunc_gap = trunc_gap(
                "marginal_time_only_IPCW_weighted_KM"
            )
            standardized_trunc_gap = trunc_gap(
                "baseline_stabilized_IPCW_standardized_outcome"
            )

            flags = []

            if not censor_models_converged:
                flags.append("censor_model_nonconvergence")
            if not outcome_models_converged:
                flags.append("outcome_model_nonconvergence")

            # Positivity is evaluated on rows that actually update IPCW.
            if (
                max_pct_update_p_lt001 is not None
                and max_pct_update_p_lt001 > 1.0
            ):
                flags.append("gt1pct_weight_updates_have_pden_lt0.01")

            if (
                max_update_ratio_p99 is not None
                and max_update_ratio_p99 > 5.0
            ):
                flags.append("p99_single_interval_marginal_ratio_gt5")

            if (
                max_pp_p99_weight is not None
                and max_pp_p99_weight > 5.0
            ):
                flags.append("person_period_p99_weight_gt5")

            if (
                max_pp_weight is not None
                and max_pp_weight > 10.0
            ):
                flags.append("person_period_max_weight_gt10")

            if (
                max_completion_p99 is not None
                and max_completion_p99 > 10.0
            ):
                flags.append("completion_p99_weight_gt10")

            if (
                max_riskset_smd is not None
                and max_riskset_smd > 0.10
            ):
                flags.append("riskset_abs_SMD_gt0.10")

            if min_arm is not None and min_arm < 100:
                flags.append("observed_action_min_arm_lt100")

            if clear_fraction is not None and clear_fraction < 0.10:
                flags.append("clear_action_fraction_lt0.10")

            if rd_gap is not None and rd_gap > 0.03:
                flags.append("dual_estimators_6h_RD_gap_gt0.03")

            if same_direction != 1:
                flags.append("dual_estimators_6h_direction_disagrees")

            if (
                marginal_trunc_gap is not None
                and marginal_trunc_gap > 0.02
            ):
                flags.append("marginal_raw_vs_trunc_6h_RD_gap_gt0.02")

            if (
                standardized_trunc_gap is not None
                and standardized_trunc_gap > 0.03
            ):
                flags.append("standardized_raw_vs_trunc_6h_RD_gap_gt0.03")

            corrected_rows.append(
                dict(
                    spec_name=name,
                    censor_models_converged=int(censor_models_converged),
                    outcome_models_converged=int(outcome_models_converged),
                    min_p01_denominator_probability_among_weight_updates=min_update_p01,
                    max_pct_weight_updates_pden_lt0_01=max_pct_update_p_lt001,
                    max_p99_single_interval_marginal_ratio=max_update_ratio_p99,
                    marginal_person_period_p99_weight=max_pp_p99_weight,
                    marginal_person_period_max_weight=max_pp_weight,
                    marginal_completion_p99_weight=max_completion_p99,
                    max_marginal_riskset_abs_SMD=max_riskset_smd,
                    observed_action_min_arm=min_arm,
                    clear_action_fraction=clear_fraction,
                    dual_estimators_6h_RD_gap=rd_gap,
                    dual_estimators_6h_same_direction=same_direction,
                    marginal_raw_vs_trunc_6h_RD_gap=marginal_trunc_gap,
                    standardized_raw_vs_trunc_6h_RD_gap=standardized_trunc_gap,
                    diagnostic_flags=(";".join(flags) if flags else ""),
                    protocol_freeze_status=(
                        "READY_FOR_PATIENT_LEVEL_BOOTSTRAP"
                        if not flags
                        else "REVIEW_BEFORE_BOOTSTRAP"
                    ),
                )
            )

        save_csv(
            "17_corrected_protocol_freeze_verdict.csv",
            corrected_rows,
        )

        # ----------------------------------------------------------
        # Manifest
        # ----------------------------------------------------------
        manifest = {
            "stage": "G10_frozen_protocol_CCW_subject_cluster_bootstrap",
            "script_version": "1.1",
            "run_time": datetime.now().astimezone().isoformat(),
            "error": ERRORS or None,
            "warnings": WARNINGS or None,
            "reconciliation": {
                "G9_frozen_eligible_expected": EXPECTED_FINAL_ELIGIBLE_N,
                "G9_frozen_eligible_observed": len(ccw_records),
                "primary_actions_expected": EXPECTED_FINAL_PRIMARY_ACTIONS,
                "primary_actions_observed": {
                    "NE_first": int(
                        next(
                            r["NE_first"]
                            for r in observed_action_rows
                            if r["spec_name"] == "primary_dovss120_immediate"
                        )
                    ),
                    "AVP_first": int(
                        next(
                            r["AVP_first"]
                            for r in observed_action_rows
                            if r["spec_name"] == "primary_dovss120_immediate"
                        )
                    ),
                    "near_simultaneous": int(
                        next(
                            r["near_simultaneous"]
                            for r in observed_action_rows
                            if r["spec_name"] == "primary_dovss120_immediate"
                        )
                    ),
                    "no_clear_action": int(
                        next(
                            r["no_clear_action"]
                            for r in observed_action_rows
                            if r["spec_name"] == "primary_dovss120_immediate"
                        )
                    ),
                },
            },
            "frozen_decision_point": (
                "G9 dovss_stable120_no_otherpressor2h"
            ),
            "population": {
                "adult": True,
                "sepsis3_onset_before_or_at_t0": True,
                "continuous_NE_AVP_minutes": DUAL_LOOKBACK_MIN,
                "MAP_stability_minutes": MAP_LOOKBACK_MIN,
                "left_limit_NE_max_mcg_kg_min": PRIMARY_NE_MAX,
                "left_limit_AVP_max_units_hour": PRIMARY_AVP_MAX,
                "no_meaningful_escalation_lookback_minutes": (
                    NO_ESCALATION_LOOKBACK_MIN
                ),
                "no_NE_AVP_transition_exactly_at_t0": True,
                "no_DOVSS_scale_NE_AVP_change_prior_minutes": 120,
                "no_nonstudy_pressor_exposure_prior_minutes": 120,
            },
            "specifications": SPECS,
            "primary_action_timing_rule": (
                "first observed threshold-crossing rate transition after the "
                "frozen pre-action t0; no future persistence information"
            ),
            "primary_endpoint": {
                "definition": (
                    "two MAP values <65 mmHg separated by 15-60 minutes"
                ),
                "window_hours": 6,
            },
            "secondary_endpoint_windows_hours": [12, 24],
            "primary_estimator": {
                "name": (
                    "marginal/time-only stabilized IPCW + weighted Kaplan-Meier"
                ),
                "censor_weight_numerator": (
                    "flexible interval time only, fitted separately by arm"
                ),
                "censor_weight_denominator": (
                    "baseline covariates + flexible interval time + "
                    "time-varying hemodynamic/treatment history"
                ),
                "probability_clip": [
                    PROB_CLIP_LO,
                    PROB_CLIP_HI,
                ],
                "primary_weight_version": (
                    "untruncated cumulative stabilized weight"
                ),
                "weight_truncation_sensitivity": [
                    WEIGHT_LO_Q,
                    WEIGHT_HI_Q,
                ],
            },
            "parallel_model_based_sensitivity": {
                "name": (
                    "baseline-stabilized IPCW + weighted pooled logistic "
                    "outcome model + standardization"
                ),
                "role": "point-estimate sensitivity analysis",
            },
            "bootstrap": {
                "cluster_unit": "subject_id",
                "requested_replicates": BOOTSTRAP_REPS,
                "seed": BOOTSTRAP_SEED,
                "CI_method": "percentile 95%",
                "minimum_success_fraction_required": (
                    BOOTSTRAP_MIN_SUCCESS_FRACTION
                ),
                "bootstrap_output": "22_bootstrap_summary_6h.csv",
            },
            "privacy": {
                "patient_identifiers_used_internally_for_clustering": True,
                "patient_identifiers_exported": False,
            },
            "critical_limitations": [
                (
                    "Primary inferential endpoint is 6 hours. The 12/24-hour "
                    "estimates are exploratory because natural ICU discharge/death "
                    "censoring is not additionally inverse-probability weighted."
                ),
                (
                    "The septic-shock phenotype is an operational real-world "
                    "target-trial population and does not prove adequate fluid "
                    "resuscitation as required by the full Sepsis-3 shock definition."
                ),
                (
                    "Percentile bootstrap intervals are used for the frozen "
                    "primary estimator; no multiplicity adjustment is applied "
                    "to sensitivity analyses."
                ),
                (
                    "Only aggregate results and bootstrap-replicate summaries "
                    "are exported; subject/stay identifiers are not saved."
                ),
            ],
        }

        (
            OUT / "23_manifest.json"
        ).write_text(
            json.dumps(
                manifest,
                ensure_ascii=False,
                indent=2,
                default=str,
            ),
            encoding="utf-8",
        )

        log(
            "=== G10 complete ==="
        )
        log(
            "Output:",
            OUT,
        )
        log(
            "Errors:",
            len(
                ERRORS
            ),
            "Warnings:",
            len(
                WARNINGS
            ),
        )
        log(
            "Please upload the entire G10 output folder as a ZIP."
        )

    finally:
        try:
            conn.close()
        finally:
            if OUT.exists():
                (
                    OUT
                    / "24_run_log.txt"
                ).write_text(
                    "\n".join(
                        LOG
                    ),
                    encoding="utf-8",
                )

                manifest_path = (
                    OUT
                    / "23_manifest.json"
                )

                if not manifest_path.exists():
                    emergency = {
                        "stage": (
                            "G10_frozen_protocol_CCW_subject_cluster_bootstrap"
                        ),
                        "script_version": "1.1",
                        "run_time": (
                            datetime.now()
                            .astimezone()
                            .isoformat()
                        ),
                        "completed": False,
                        "errors": (
                            ERRORS or None
                        ),
                        "warnings": (
                            WARNINGS or None
                        ),
                        "note": (
                            "Execution stopped before normal manifest creation. "
                            "Inspect 24_run_log.txt."
                        ),
                    }

                    manifest_path.write_text(
                        json.dumps(
                            emergency,
                            ensure_ascii=False,
                            indent=2,
                            default=str,
                        ),
                        encoding="utf-8",
                    )


if __name__ == "__main__":
    main()
