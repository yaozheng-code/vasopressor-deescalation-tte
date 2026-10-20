"""Run the unchanged frozen G10 inferential analysis and/or G11 Table 1.

The external wrapper configures local PostgreSQL and file locations without
altering the frozen analysis functions.
"""
from __future__ import annotations

import argparse
import getpass
import importlib.metadata
import importlib.util
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STAGE_FILES = {
    "g10": ROOT / "analysis" / "G10_v11_MIMIC_NE_AVP_final_CCW_subject_bootstrap.py",
    "g11": ROOT / "analysis" / "G11_MIMIC_NE_AVP_manuscript_Table1.py",
}
SELFTEST_ENV = {"g10": "G10_SELFTEST_ONLY", "g11": "G11_SELFTEST_ONLY"}


def load_frozen_module(stage: str):
    path = STAGE_FILES[stage]
    spec = importlib.util.spec_from_file_location(f"_frozen_neavp_{stage}", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run_stage(stage: str, args, output_dir: Path) -> None:
    mod = load_frozen_module(stage)
    # Only machine-dependent configuration values change. Scientific
    # definitions, SQL and outcome/estimator code stay byte-for-byte frozen.
    mod.DBNAME = args.db_name
    mod.SOCKET_DIR = args.db_host
    mod.PORT = args.db_port
    mod.USER = args.db_user
    mod.OUT = output_dir
    # When reproducing manuscript results, the canonical values are 500,
    # 20261015. Other values are exploratory diagnostics only.
    mod.BOOTSTRAP_REPS = args.n_bootstrap
    mod.BOOTSTRAP_SEED = args.seed

    env_name = SELFTEST_ENV[stage]
    old = os.environ.get(env_name)
    try:
        if args.self_test_only:
            os.environ[env_name] = "1"
        else:
            os.environ.pop(env_name, None)
        print(f"\n>>> Running frozen {stage.upper()}")
        if not args.self_test_only:
            print(f"Output directory: {output_dir}")
        mod.main()
    finally:
        if old is None:
            os.environ.pop(env_name, None)
        else:
            os.environ[env_name] = old


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def metadata(args: argparse.Namespace) -> dict:
    from check_integrity import ROOT as _ROOT
    import hashlib

    return {
        "run_started_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": args.step,
        "bootstrap_replicates": args.n_bootstrap,
        "bootstrap_seed": args.seed,
        "db_name": args.db_name,
        "db_host": args.db_host,
        "db_port": args.db_port,
        "db_user": "[not recorded]",
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "versions": {
            x: package_version(x) for x in ("numpy", "pandas", "psycopg")
        },
        "source_sha256": {
            stage: hashlib.sha256(path.read_bytes()).hexdigest()
            for stage, path in STAGE_FILES.items()
        },
        "note": "Original analysis files unchanged; output files are PRIVATE",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", choices=("all", "g10", "g11"), default="all",
                    help="Default: main G10 followed by G11 baseline Table 1")
    ap.add_argument("--db-name", default="mimiciv")
    ap.add_argument("--db-host", default="/var/run/postgresql",
                    help="Unix socket directory, or network hostname")
    ap.add_argument("--db-port", default=5432, type=int)
    ap.add_argument("--db-user", default=getpass.getuser())
    ap.add_argument("--output-dir", type=Path, default=ROOT / "outputs",
                    help="Private local results directory; ignored by Git")
    ap.add_argument("--n-bootstrap", type=int, default=500)
    ap.add_argument("--seed", type=int, default=20261015)
    ap.add_argument("--self-test-only", action="store_true",
                    help="Run built-in synthetic tests; no DB connection")
    ap.add_argument("--skip-preflight", action="store_true")
    ap.add_argument("--skip-verification", action="store_true")
    args = ap.parse_args()

    if args.n_bootstrap < 1:
        ap.error("--n-bootstrap must be >=1")
    from check_integrity import check
    check()

    stages = ["g10", "g11"] if args.step == "all" else [args.step]
    if args.self_test_only:
        for stage in stages:
            run_stage(stage, args, ROOT / "outputs" / "not_created_for_selftest")
        print("\nOffline synthetic self-tests: PASS (no patient data accessed)")
        return 0

    if not args.skip_preflight:
        from preflight import check_database
        check_database(
            args.db_name, args.db_host, args.db_port, args.db_user
        )
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    runroot = args.output_dir.resolve() / timestamp
    runroot.mkdir(parents=True, exist_ok=False)
    (runroot / "run_environment.json").write_text(
        json.dumps(metadata(args), indent=2) + "\n", encoding="utf-8"
    )

    for stage in stages:
        run_stage(stage, args, runroot / stage)

    if not args.skip_verification:
        from verify_results import audit, display_report
        report = audit(
            g10_dir=(runroot/"g10") if "g10" in stages else None,
            g11_dir=(runroot/"g11") if "g11" in stages else None,
            allow_alternate_bootstrap=args.n_bootstrap != 500,
        )
        display_report(report)
        if report["errors"]:
            raise SystemExit(1)

    print(f"\nComplete. PRIVATE local output (do not publish): {runroot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
