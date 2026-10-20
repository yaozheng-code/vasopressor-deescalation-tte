"""Check byte-for-byte SHA-256 integrity of the frozen G10/G11 scripts."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def check() -> None:
    record = json.loads((ROOT / "source_sha256.json").read_text(encoding="utf-8"))
    errors = []
    for basename, expected in record["sha256"].items():
        path = ROOT / "analysis" / basename
        if not path.is_file():
            errors.append(f"{basename}: missing")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        print(f"{basename}: {'PASS' if actual == expected else 'FAIL'}")
        if actual != expected:
            errors.append(f"{basename}: hash mismatch ({actual})")
    if errors:
        raise SystemExit("Frozen-source integrity failure:\n" + "\n".join(errors))
    print("Frozen-source SHA-256 checks: PASS")


if __name__ == "__main__":
    check()
