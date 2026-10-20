"""No-database tests for a public, code-only reproducibility release."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from check_integrity import check
from verify_results import audit


def write_csv(path: Path, data: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(data[0]))
        wr.writeheader()
        wr.writerows(data)


class ReleaseTests(unittest.TestCase):
    def test_source_integrity(self):
        check()

    def test_syntax_all_python(self):
        for file in ROOT.rglob("*.py"):
            with self.subTest(file=file.name):
                compile(file.read_text(encoding="utf-8"), str(file), "exec")

    def test_reference_math_consistent(self):
        r = json.loads((ROOT / "expected_results.json").read_text())
        p = r["primary_untruncated_6h"]
        self.assertAlmostEqual(
            p["risk_NE_first"] - p["risk_AVP_first"],
            p["risk_difference_NE_minus_AVP"],
            places=5,
        )
        self.assertAlmostEqual(
            p["risk_NE_first"] / p["risk_AVP_first"],
            p["risk_ratio_NE_over_AVP"],
            places=4,
        )

    def test_built_in_synthetic_self_tests(self):
        process = subprocess.run(
            [
                sys.executable,
                str(ROOT / "run_pipeline.py"),
                "--self-test-only",
                "--step",
                "all",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=90,
        )
        self.assertEqual(
            process.returncode,
            0,
            f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}",
        )
        self.assertIn("Offline synthetic self-tests: PASS", process.stdout)

    def test_output_verifier_on_synthetic_aggregate_fixture(self):
        # This fixture contains only prepublished headline AGGREGATES,
        # never patient-level MIMIC data.
        r = json.loads((ROOT / "expected_results.json").read_text())
        with tempfile.TemporaryDirectory() as td:
            g10 = Path(td) / "g10"
            g11 = Path(td) / "g11"
            g10.mkdir()
            g11.mkdir()
            write_csv(
                g10 / "01_frozen_protocol_population.csv",
                [r["population"]],
            )
            write_csv(
                g10 / "05_observed_actions_in_CCW_cohort.csv",
                [{
                    "spec_name": "primary_dovss120_immediate",
                    **r["primary_action_counts"],
                }],
            )
            main = r["primary_untruncated_6h"]
            rows = [
                {
                    "spec_name": main["spec_name"],
                    "horizon_hours": main["horizon_hours"],
                    "weight_version": "untruncated",
                    "bootstrap_requested": 500,
                    "bootstrap_successful": 500,
                    **{k:v for k,v in main.items()
                       if k not in ("spec_name","horizon_hours","weight_version")},
                }
            ]
            for name, rd in r["sensitivity_6h_untruncated_RD"].items():
                rows.append({
                    "spec_name": name,
                    "horizon_hours": 6,
                    "weight_version": "untruncated",
                    "bootstrap_requested": 500,
                    "bootstrap_successful": 500,
                    "risk_difference_NE_minus_AVP": rd,
                })
            # Unequal row fields are normalized to a common CSV header.
            keys = list(dict.fromkeys(k for row in rows for k in row))
            write_csv(
                g10 / "22_bootstrap_summary_6h.csv",
                [{k:row.get(k,"") for k in keys} for row in rows],
            )
            (g10 / "23_manifest.json").write_text(json.dumps({
                "script_version": "1.1",
                "error": None,
                "warnings": None,
                "bootstrap": {
                    "cluster_unit": "subject_id",
                    "seed": 20261015,
                    "requested_replicates": 500,
                }
            }))
            (g10 / "24_run_log.txt").write_text("=== G10 complete ===")
            (g11 / "23_manifest.json").write_text(json.dumps({
                "script_version": "1.0",
                "eligible_ICU_stays": 2260,
                "unique_subjects": 2168,
            }))
            write_csv(
                g11 / "04_Table1_baseline_characteristics.csv",
                [
                    {"Characteristic": c, "Overall": "synthetic"}
                    for c in (
                        "Age, years",
                        "Norepinephrine dose at t0, mcg/kg/min",
                        "Vasopressin dose at t0, U/h",
                    )
                ],
            )
            ok = audit(g10,g11)
            self.assertEqual(ok["errors"], [])
            # Deliberately corrupt an aggregate field; verifier must fail.
            population = {**r["population"], "eligible_stays": 2200}
            write_csv(
                g10 / "01_frozen_protocol_population.csv",
                [population],
            )
            bad = audit(g10, g11)
            self.assertTrue(
                any("eligible_stays" in x for x in bad["errors"])
            )


if __name__ == "__main__":
    unittest.main()
