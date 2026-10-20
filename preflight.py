"""Read-only availability check for the MIMIC-IV PostgreSQL schemas.

No patient-level values are displayed or saved.
"""
from __future__ import annotations

import argparse
import getpass

REQUIRED_COLUMNS = {
    "mimiciv_icu.inputevents": {
        "stay_id", "subject_id", "itemid", "starttime", "endtime",
        "rate", "rateuom", "statusdescription",
    },
    "mimiciv_icu.icustays": {
        "stay_id", "subject_id", "hadm_id", "intime", "outtime",
        "first_careunit",
    },
    "mimiciv_icu.d_items": {"itemid", "label"},
    "mimiciv_hosp.admissions": {
        "hadm_id", "subject_id", "admittime", "race", "deathtime",
        "hospital_expire_flag",
    },
    "mimiciv_hosp.patients": {
        "subject_id", "gender", "anchor_age", "anchor_year",
    },
    "mimiciv_hosp.labevents": {
        "subject_id", "hadm_id", "itemid", "charttime", "valuenum",
    },
    "mimiciv_hosp.prescriptions": {
        "subject_id", "hadm_id", "starttime", "stoptime",
        "drug",
    },
    "mimiciv_derived.sepsis3": {
        "stay_id", "sepsis3", "suspected_infection_time",
        "sofa_time", "sofa_score",
    },
    "mimiciv_derived.vitalsign": {
        "stay_id", "charttime", "mbp", "heart_rate",
    },
}


def check_database(db_name: str, db_host: str, db_port: int, db_user: str) -> None:
    try:
        import psycopg
        from psycopg.rows import dict_row
    except ImportError as exc:
        raise RuntimeError(
            "psycopg 3 is required: pip install -r requirements.txt"
        ) from exc

    print("Checking PostgreSQL schemas (metadata only, no patient records)...")
    with psycopg.connect(
        dbname=db_name,
        host=db_host,
        port=db_port,
        user=db_user,
        row_factory=dict_row,
    ) as conn:
        conn.execute("SET default_transaction_read_only = on")
        version = conn.execute("SHOW server_version").fetchone()
        print("PostgreSQL:", version["server_version"])
        missing = []
        for full_name, required in REQUIRED_COLUMNS.items():
            schema, table = full_name.split(".", 1)
            rows = conn.execute(
                "SELECT column_name "
                "FROM information_schema.columns "
                "WHERE table_schema = %s AND table_name = %s",
                (schema, table),
            ).fetchall()
            columns = {r["column_name"] for r in rows}
            missing_columns = sorted(required - columns)
            if not columns:
                print("  MISSING TABLE:", full_name)
                missing.append(full_name + " (missing or inaccessible)")
            elif missing_columns:
                print("  INCOMPATIBLE:", full_name, "missing", missing_columns)
                missing.append(full_name + " missing " + ", ".join(missing_columns))
            else:
                print("  OK:", full_name)

        if missing:
            raise RuntimeError(
                "Required MIMIC-IV tables/columns are unavailable:\n"
                + "\n".join(" - " + x for x in missing)
                + "\nBuild/check mimiciv_derived concepts and dataset version."
            )
    print("Read-only schema preflight: PASS")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-name", default="mimiciv")
    ap.add_argument("--db-host", default="/var/run/postgresql")
    ap.add_argument("--db-port", type=int, default=5432)
    ap.add_argument("--db-user", default=getpass.getuser())
    args = ap.parse_args()
    check_database(args.db_name, args.db_host, args.db_port, args.db_user)


if __name__ == "__main__":
    main()
