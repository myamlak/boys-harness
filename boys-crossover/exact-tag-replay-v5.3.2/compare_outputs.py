#!/usr/bin/env python3
"""Compare the v5.3.2 replay's measurement records with the deposited v5.3.1
record by SHA-256.

usage: python -I compare_outputs.py

The deposited record is `boys-crossover/exact-tag-replay`. Every deposited
file that is still present is hashed directly; `run/paths.csv` was removed
after its digest was taken, so that row compares the replay's digest with the
digest in the deposited `digests.txt`. The exit code is 0 only when all 13
measurement records are identical.
"""
import hashlib
import re
import sys
from pathlib import Path

RECORD = Path(__file__).resolve().parent
DEPOSIT = RECORD.parent / "exact-tag-replay"

RECORDS = [
    "run/cells-9353.csv",
    "run/admission-9353.csv",
    "run/replica_gain-9353.csv",
    "run/cells-9352.csv",
    "run/admission-9352.csv",
    "run/replica_gain-9352.csv",
    "run/cells-sweep2.csv",
    "run/admission-sweep2.csv",
    "run/paths.csv",
    "run/certificates/cells_cert.csv",
    "run/certificates/band_cert.csv",
    "run/certificates/model_cells.csv",
    "run/certificates/summary.json",
]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def deposited_digests():
    table = {}
    for line in (DEPOSIT / "digests.txt").read_text(encoding="utf-8-sig").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        match = re.search(r"\b([0-9a-f]{64})\b", line)
        if not match:
            continue
        tail = line[match.end():].strip()
        if not tail:
            continue
        table[tail.split()[0]] = match.group(1)
    return table


def main():
    expected = deposited_digests()
    failures = 0
    print("%-34s %-12s %-64s %s" % ("record", "bytes", "sha256", "result"))
    for path in RECORDS:
        produced = RECORD / path
        if path not in expected:
            print("%-34s %-12s %-64s %s" % (path, "-", "-", "NO DEPOSITED DIGEST"))
            failures += 1
            continue
        if not produced.is_file():
            print("%-34s %-12s %-64s %s" % (path, "-", "-", "PRODUCED FILE MISSING"))
            failures += 1
            continue
        got = sha256(produced)
        deposited_file = DEPOSIT / path
        if deposited_file.is_file():
            dep = sha256(deposited_file)
            if dep != expected[path]:
                print("%-34s %-12s %-64s %s" % (
                    path, produced.stat().st_size, got,
                    "DEPOSIT FILE vs ITS digests.txt MISMATCH"))
                failures += 1
                continue
            note = "identical (deposited file compared byte for byte)"
        else:
            note = "identical (deposited raw removed; digest from its digests.txt)"
        if got != expected[path]:
            note = "DIFFERS: deposited %s" % expected[path]
            failures += 1
        print("%-34s %-12d %-64s %s" % (path, produced.stat().st_size, got, note))
    print()
    print("%d of %d measurement records identical" % (len(RECORDS) - failures, len(RECORDS)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
