#!/usr/bin/env python3
"""Reproduce the v5.0.0 comparison for this v5.1.1 record.

The v5.1.1 records are under ``sweeps/``; the deposited v5.0.0 records are
``../`` (first sweep) and ``../../amplification-form/`` (second sweep).
Line endings are normalised before every comparison.  A shared file is

  SAME    the v5.0.0 rows equal the v5.1.1 rows;
  PREFIX  every v5.0.0 row is the corresponding v5.1.1 row, and the
          v5.1.1 file carries the recovered argument's rows at its end.

The script also checks the record's counts and the recovered argument
x = 5.949924076054243 = nextafter(kRegionA1Edge, +inf).

Usage:
    python check_v500_compare.py

Exit code 0 when every comparison passes.
"""

import csv
import json
import lzma
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
V500_PRODUCT = ROOT.parent
V500_AMPFORM = ROOT.parent.parent / "amplification-form"
V511_PRODUCT = ROOT / "sweeps" / "product_amplification"
V511_AMPFORM = ROOT / "sweeps" / "amplification_form"

RECOVERED_X = 5.949924076054243


def read_bytes(path):
    """The file's bytes with CRLF folded to LF; .xz siblings are opened."""
    if path.suffix == ".xz":
        with lzma.open(path, "rb") as handle:
            data = handle.read()
    else:
        data = path.read_bytes()
    return data.replace(b"\r\n", b"\n")


def resolve(path):
    """The plain file, or its .xz sibling when the plain file is absent."""
    if path.is_file():
        return path
    packed = path.with_name(path.name + ".xz")
    return packed if packed.is_file() else None


def rows_of(data):
    """The file's lines without the trailing newline's empty element."""
    lines = data.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    return lines


def compare(name, v500, v511, mode):
    """Print one SAME/PREFIX/DIFF line; return True when it passes."""
    old_path = resolve(v500)
    new_path = resolve(v511)
    if old_path is None:
        print("  %-30s MISSING in the v5.0.0 deposit" % name)
        return False
    if new_path is None:
        print("  %-30s MISSING in the v5.1.1 record" % name)
        return False

    old = read_bytes(old_path)
    new = read_bytes(new_path)
    old_rows = rows_of(old)
    new_rows = rows_of(new)

    if mode == "same":
        if old == new:
            print("  %-30s SAME (%d lines)" % (name, len(new_rows)))
            return True
        for index, (left, right) in enumerate(zip(old_rows, new_rows), 1):
            if left != right:
                print("  %-30s DIFF at line %d" % (name, index))
                return False
        print("  %-30s DIFF in length (%d vs %d lines)"
              % (name, len(old_rows), len(new_rows)))
        return False

    if new_rows[:len(old_rows)] == old_rows:
        print("  %-30s PREFIX (v5.0.0 %d lines; v5.1.1 %d lines)"
              % (name, len(old_rows), len(new_rows)))
        return True
    for index, (left, right) in enumerate(zip(old_rows, new_rows), 1):
        if left != right:
            print("  %-30s DIFF at line %d" % (name, index))
            return False
    print("  %-30s DIFF in length (%d vs %d lines)"
          % (name, len(old_rows), len(new_rows)))
    return False


def compare_validation():
    """The validation block differs only in the appended row's counts."""
    old_path = V500_PRODUCT / "data" / "reference_validation.txt"
    new_path = V511_PRODUCT / "data" / "reference_validation.txt"
    if not old_path.is_file() or not new_path.is_file():
        print("  %-30s MISSING" % "data/reference_validation.txt")
        return False

    def normalised(path):
        lines = []
        for line in read_bytes(path).split(b"\n"):
            if line.startswith(b"arguments:"):
                continue
            line = re.sub(rb"grid: .*", b"grid: PATH", line)
            line = re.sub(rb"over \d+ cells", b"over N cells", line)
            lines.append(line)
        return b"\n".join(lines)

    if normalised(old_path) == normalised(new_path):
        print("  %-30s SAME (route values; counts grow by the recovered row)"
              % "data/reference_validation.txt")
        return True
    print("  %-30s DIFF" % "data/reference_validation.txt")
    return False


def count_rows(path):
    with open(path, newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def recovered_rows(path):
    with open(path, newline="") as handle:
        return [row for row in csv.DictReader(handle)
                if float(row["x"]) == RECOVERED_X]


def check_counts():
    """The record's counts and the recovered argument's rows."""
    failures = 0
    provenance = json.loads(
        (V511_PRODUCT / "results" / "provenance.json").read_text(encoding="utf-8"))
    edge = float(provenance["constants"]["kRegionA1Edge"])
    successor = math.nextafter(edge, math.inf)

    print("recovered argument")
    if successor == RECOVERED_X:
        print("  %-30s OK (%r = nextafter(%r, +inf))"
              % ("x", RECOVERED_X, edge))
    else:
        print("  %-30s FAIL: nextafter(%r, +inf) = %r"
              % ("x", edge, successor))
        failures += 1

    with open(V511_PRODUCT / "data" / "grid.csv", newline="") as handle:
        grid = list(csv.DictReader(handle))
    bands = {"A1": 0, "A2": 0}
    recovered = []
    for row in grid:
        bands[row["band"]] += 1
        if "recovered" in row["sets"].split("+"):
            recovered.append(row)
    total = len(grid)
    designed = total - len(recovered)
    if total == 9353 and bands == {"A1": 4812, "A2": 4541} and designed == 9352:
        print("  %-30s OK (9353 = 4812 A1 + 4541 A2; 9352 designed + %d recovered)"
              % ("arguments", len(recovered)))
    else:
        print("  %-30s FAIL (total %d, bands %r, designed %d)"
              % ("arguments", total, bands, designed))
        failures += 1
    if (len(recovered) == 1 and recovered[0]["band"] == "A2"
            and float(recovered[0]["x"]) == RECOVERED_X
            and grid[-1] is recovered[0]):
        print("  %-30s OK (last grid row, band A2, tag recovered)"
              % "recovered grid row")
    else:
        print("  %-30s FAIL (%r)" % ("recovered grid row", recovered))
        failures += 1

    cells = count_rows(V511_PRODUCT / "results" / "cells.csv")
    if cells == 308649 and cells == total * 33:
        print("  %-30s OK (308649 = 9353 x 33)" % "cells")
    else:
        print("  %-30s FAIL (%d rows)" % ("cells", cells))
        failures += 1

    reference = recovered_rows(V511_PRODUCT / "data" / "reference.csv")
    if len(reference) == 1 and all(reference[0]["f%d" % n] for n in range(33)):
        print("  %-30s OK (one row, 33 orders)" % "recovered reference row")
    else:
        print("  %-30s FAIL (%r)" % ("recovered reference row", reference))
        failures += 1

    recovered_cells = recovered_rows(V511_PRODUCT / "results" / "cells.csv")
    orders = sorted(int(row["n"]) for row in recovered_cells)
    if len(recovered_cells) == 33 and orders == list(range(33)):
        print("  %-30s OK (33 rows, n = 0..32)" % "recovered cells")
    else:
        print("  %-30s FAIL (%d rows)" % ("recovered cells", len(recovered_cells)))
        failures += 1

    admission = recovered_rows(V511_PRODUCT / "results" / "admission.csv")
    if (len(admission) == 1 and admission[0]["direction"] == "upward"
            and admission[0]["checked"] == admission[0]["matched"] == "17"):
        print("  %-30s OK (upward, 17 of 17 bit-identical)"
              % "recovered admission row")
    else:
        print("  %-30s FAIL (%r)" % ("recovered admission row", admission))
        failures += 1

    gains = recovered_rows(V511_PRODUCT / "results" / "replica_gain.csv")
    gain_orders = sorted(int(row["n"]) for row in gains)
    if (len(gains) == 16 and all(row["case"] == "up" for row in gains)
            and gain_orders == list(range(1, 17))):
        print("  %-30s OK (16 rows, case up, n = 1..16)"
              % "recovered gain rows")
    else:
        print("  %-30s FAIL (%d rows)" % ("recovered gain rows", len(gains)))
        failures += 1

    return failures


def main():
    failures = 0

    print("first sweep: v5.0.0 deposit vs v5.1.1 record")
    for name, mode in (
        ("data/grid.csv", "prefix"),
        ("data/reference.csv", "prefix"),
        ("data/documented.csv", "same"),
        ("results/admission.csv", "prefix"),
        ("results/band_summary.csv", "same"),
        ("results/cells.csv", "prefix"),
        ("results/decomposition.csv", "same"),
        ("results/gain_up.csv", "same"),
        ("results/growth.csv", "same"),
        ("results/per_order.csv", "same"),
        ("results/replica_gain.csv", "prefix"),
    ):
        if not compare(name, V500_PRODUCT / name, V511_PRODUCT / name, mode):
            failures += 1
    if not compare_validation():
        failures += 1

    print("second sweep: v5.0.0 deposit vs v5.1.1 record")
    for name in (
        "data/fits.csv",
        "data/fits_validation.txt",
        "data/grid.csv",
        "data/reference.csv",
        "data/reference_validation.txt",
        "results/admission.csv",
        "results/band_summary.csv",
        "results/cells.csv",
        "results/paths.csv",
        "results/per_order.csv",
        "results/prediction.csv",
        "results/prediction_scan.csv",
        "results/transfer_cells.csv",
    ):
        if not compare(name, V500_AMPFORM / name, V511_AMPFORM / name, "same"):
            failures += 1
    new_record = V511_AMPFORM / "results" / "prediction_record.csv"
    print("  %-30s %s" % ("results/prediction_record.csv",
                          "NEW in the v5.1.1 record"
                          if new_record.is_file() else "MISSING"))
    if not new_record.is_file():
        failures += 1

    failures += check_counts()

    print()
    if failures:
        print("FAIL: %d comparison(s) failed" % failures)
        return 1
    print("PASS: every shared row is unchanged; the v5.1.1 record appends "
          "the recovered argument's rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
