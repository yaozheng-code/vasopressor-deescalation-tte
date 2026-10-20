"""Validate private aggregate G10/G11 outputs against the frozen report.

Works without the MIMIC database and does not ingest individual patient data.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def _csv_rows(path: Path, errors: list[str]) -> list[dict]:
    if not path.is_file():
        errors.append(f"Missing file: {path.name}")
        return []
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def _manifest(path: Path, errors: list[str]) -> dict:
    if not path.is_file():
        errors.append(f"Missing file: {path.name}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"Invalid {path.name}: {exc}")
        return {}


def _eq(label: str, actual, expected, errors: list[str], passed: list[str]) -> None:
    try:
        actual_n = int(actual)
    except (ValueError, TypeError):
        errors.append(f"{label}: noninteger or missing value {actual!r}")
        return
    if actual_n == int(expected):
        passed.append(f"{label}: {actual_n}")
    else:
        errors.append(f"{label}: observed {actual_n}, expected {expected}")


def _close(
    label: str,
    actual,
    expected: float,
    tolerance: float,
    errors: list[str],
    passed: list[str],
) -> None:
    try:
        v = float(actual)
    except (ValueError, TypeError):
        errors.append(f"{label}: missing/non-numeric value {actual!r}")
        return
    if math.isfinite(v) and abs(v - expected) <= tolerance:
        passed.append(f"{label}: {v:.6f}")
    else:
        errors.append(
            f"{label}: observed {v}, expected {expected} "
            f"(tolerance ±{tolerance})"
        )


def audit(
    g10_dir: Path | None = None,
    g11_dir: Path | None = None,
    *,
    allow_alternate_bootstrap: bool = False,
) -> dict:
    expected = json.loads(
        (ROOT / "expected_results.json").read_text(encoding="utf-8")
    )
    errors: list[str] = []
    warnings: list[str] = []
    passed: list[str] = []
    if g10_dir is None and g11_dir is None:
        errors.append("Provide at least one of --g10-dir or --g11-dir")
        return {"passed": passed, "warnings": warnings, "errors": errors}

    if g10_dir is not None:
        d = Path(g10_dir)
        man = _manifest(d / "23_manifest.json", errors)
        if man:
            if man.get("script_version") == "1.1":
                passed.append("G10 version: 1.1")
            else:
                errors.append("G10 manifest version is not 1.1")
            if man.get("error") not in (None, [], {}):
                errors.append(f"G10 manifest error: {man['error']}")
            if man.get("warnings"):
                warnings.append(f"G10 run recorded warnings: {man['warnings']}")
            b = man.get("bootstrap", {})
            if b.get("cluster_unit") == "subject_id":
                passed.append("Bootstrap cluster unit: subject_id")
            else:
                errors.append("Bootstrap cluster unit is not subject_id")
            _eq("Bootstrap seed", b.get("seed"),
                expected["bootstrap_seed"], errors, passed)
            if allow_alternate_bootstrap:
                try:
                    if int(b.get("requested_replicates", 0)) < 1:
                        errors.append("Invalid alternate bootstrap replicate count")
                except (ValueError, TypeError):
                    errors.append("Bootstrap replicate count missing")
                warnings.append(
                    "Noncanonical bootstrap replicate count: results are "
                    "not a literal 500-replicate manuscript reproduction."
                )
            else:
                _eq("Bootstrap requested", b.get("requested_replicates"),
                    expected["bootstrap_replicates"], errors, passed)

        log = d / "24_run_log.txt"
        if log.is_file() and "=== G10 complete ===" in log.read_text(
            encoding="utf-8", errors="replace"
        ):
            passed.append("G10 run completion log")
        else:
            errors.append("G10 completion log is absent/incomplete")

        population = _csv_rows(d / "01_frozen_protocol_population.csv", errors)
        if len(population) == 1:
            for k, exp in expected["population"].items():
                _eq(f"G10 population.{k}", population[0].get(k),
                    exp, errors, passed)
        else:
            errors.append(f"G10 population summary should have 1 row; found {len(population)}")

        actions = _csv_rows(d / "05_observed_actions_in_CCW_cohort.csv", errors)
        qualifying = [
            r for r in actions
            if r.get("spec_name") == expected["primary_untruncated_6h"]["spec_name"]
        ]
        if len(qualifying) == 1:
            for k, exp in expected["primary_action_counts"].items():
                _eq(f"Primary observed actions.{k}", qualifying[0].get(k),
                    exp, errors, passed)
        else:
            errors.append("Primary action support row absent/nonunique")

        results = _csv_rows(d / "22_bootstrap_summary_6h.csv", errors)
        primary_key = expected["primary_untruncated_6h"]
        rows = [
            r for r in results
            if r.get("spec_name") == primary_key["spec_name"]
            and r.get("weight_version") == primary_key["weight_version"]
            and r.get("horizon_hours") == str(primary_key["horizon_hours"])
        ]
        if len(rows) == 1:
            row = rows[0]
            n_reps = int(row.get("bootstrap_requested", "0"))
            success = int(row.get("bootstrap_successful", "0"))
            if success != n_reps or n_reps < 1:
                errors.append(
                    f"Bootstrap success {success}/{n_reps}, expected full success"
                )
            else:
                passed.append(f"Bootstrap successful: {success}/{n_reps}")
            if not allow_alternate_bootstrap:
                _eq("Bootstrap summary requested", n_reps,
                    expected["bootstrap_replicates"], errors, passed)
            tol = expected["tolerances"]
            for k in ["risk_NE_first", "risk_AVP_first"]:
                _close(f"Primary.{k}", row.get(k), primary_key[k],
                       tol["risk_absolute"], errors, passed)
            _close(
                "Primary.risk_difference_NE_minus_AVP",
                row.get("risk_difference_NE_minus_AVP"),
                primary_key["risk_difference_NE_minus_AVP"],
                tol["risk_difference_absolute"], errors, passed,
            )
            _close(
                "Primary.risk_ratio_NE_over_AVP",
                row.get("risk_ratio_NE_over_AVP"),
                primary_key["risk_ratio_NE_over_AVP"],
                tol["risk_ratio_absolute"], errors, passed,
            )
            if not allow_alternate_bootstrap:
                for k in ("RD_CI_low", "RD_CI_high", "RR_CI_low", "RR_CI_high"):
                    _close(f"Primary.{k}", row.get(k), primary_key[k],
                           tol["bootstrap_CI_bound_absolute"], errors, passed)
        else:
            errors.append("Primary untruncated 6-h bootstrap row absent/nonunique")

        for spec, expected_rd in expected["sensitivity_6h_untruncated_RD"].items():
            found = [
                r for r in results
                if r.get("spec_name") == spec
                and r.get("weight_version") == "untruncated"
                and r.get("horizon_hours") == "6"
            ]
            if len(found) == 1:
                _close(f"Sensitivity RD {spec}",
                       found[0].get("risk_difference_NE_minus_AVP"),
                       expected_rd,
                       expected["tolerances"]["sensitivity_RD_absolute"],
                       errors, passed)
            else:
                errors.append(f"Missing/nonunique sensitivity row {spec}")

    if g11_dir is not None:
        d = Path(g11_dir)
        man = _manifest(d / "23_manifest.json", errors)
        if man:
            if man.get("script_version") == "1.0":
                passed.append("G11 version: 1.0")
            else:
                errors.append("G11 manifest version is not 1.0")
            _eq("G11 eligible stays", man.get("eligible_ICU_stays"),
                expected["population"]["eligible_stays"], errors, passed)
            _eq("G11 unique subjects", man.get("unique_subjects"),
                expected["population"]["unique_subjects"], errors, passed)
        baseline = _csv_rows(d / "04_Table1_baseline_characteristics.csv", errors)
        characteristics = {r.get("Characteristic", "").strip() for r in baseline}
        for important in (
            "Age, years", "Norepinephrine dose at t0, mcg/kg/min",
            "Vasopressin dose at t0, U/h",
        ):
            if important in characteristics:
                passed.append("Table 1: " + important)
            else:
                errors.append("G11 Table 1 missing: " + important)

    return {"passed": passed, "warnings": warnings, "errors": errors}


def display_report(report: dict) -> None:
    for s in report["passed"]:
        print("PASS:", s)
    for s in report["warnings"]:
        print("WARNING:", s)
    for s in report["errors"]:
        print("FAIL:", s)
    if not report["errors"]:
        print("\nREPRODUCIBILITY CHECK: PASS")
    else:
        print(f"\nREPRODUCIBILITY CHECK: FAIL ({len(report['errors'])} issues)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--g10-dir", type=Path)
    ap.add_argument("--g11-dir", type=Path)
    ap.add_argument("--allow-alternate-bootstrap", action="store_true")
    args = ap.parse_args()
    report = audit(
        args.g10_dir, args.g11_dir,
        allow_alternate_bootstrap=args.allow_alternate_bootstrap,
    )
    display_report(report)
    return 1 if report["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
