#!/usr/bin/env python3
"""Reconstruct the two Paper A input lists that the harness tag v1.0.34 did
not carry as argument lists.

1. R4: the 73 arguments of the kernel's committed 45-digit grid
   (tests/data/boys_reference.csv at boys v5.3.1) below the region-A upper
   edge, with their 33 orders: 73 x 33 = 2409 cells.  The reconstruction is
   verified by recomputing the R4 agreement against the deposited 60-digit
   reference: the worst absolute difference must be the summary's 4.975493e-46
   at (n=2, x=3.162277660168379e-06) over exactly 2409 cells.

2. The 9,353rd argument of the first sweep: the last row of the deposited
   9,353-argument grid, nextafter(kRegionA1Edge, +inf), absent from the
   deposited 9,352-argument designed grid.  The reconstruction is verified
   against both grids and the C99 hexadecimal spelling.

Outputs (beside this script): r4-arguments.csv, r4-cells.csv,
argument-9353.csv, verification.txt, digests.txt.
"""

import csv
import hashlib
import math
import sys
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 80

A1_EDGE = 5.949924076054242      # kRegionA1Edge at v5.3.1
A2_EDGE = 11.899848152108484     # kRegionAEnd at v5.3.1
RECOVERED = 5.949924076054243    # nextafter(A1_EDGE, +inf)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_committed(path):
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append((int(row["n"]), float(row["x"]), mpf(row["value"])))
    return rows


def read_reference(path):
    ref = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            ref[float(row["x"])] = [mpf(row["f%d" % n]) for n in range(33)]
    return ref


def read_grid(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    tag_tree = Path(sys.argv[1])
    paper_tree = Path(sys.argv[2])
    here = Path(__file__).resolve().parent

    committed_path = tag_tree / "tests" / "data" / "boys_reference.csv"
    form_grid_path = paper_tree / "sweeps" / "amplification_form" / "data" / "grid.csv"
    form_ref_path = paper_tree / "sweeps" / "amplification_form" / "data" / "reference.csv"
    designed_grid_path = paper_tree / "sweeps" / "product_amplification" / "data" / "grid.csv"
    deposit_grid_path = (Path(sys.argv[3]) / "product-amplification" / "first-sweep-9353"
                         / "sweeps" / "product_amplification" / "data" / "grid.csv")
    # the designed 9,352-argument record deposited at the harness tag
    designed_9352_path = (Path(sys.argv[3]) / "product-amplification" / "data" / "grid.csv")

    # ---- R4 ---------------------------------------------------------------
    committed = read_committed(committed_path)
    all_x = sorted({x for _, x, _ in committed})
    args = [x for x in all_x if x < A2_EDGE]
    assert len(all_x) == 88, len(all_x)
    assert len(args) == 73, len(args)
    counts = {}
    for n, x, _ in committed:
        counts[x] = counts.get(x, 0) + 1
    assert all(counts[x] == 33 for x in args), "an R4 argument lacks 33 orders"

    band_of = {}
    for row in read_grid(form_grid_path):
        band_of[float(row["x"])] = row["band"]
    missing_band = [x for x in args if x not in band_of]
    assert not missing_band, missing_band

    with open(here / "r4-arguments.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "x", "x_hex", "orders", "cells"])
        for x in args:
            writer.writerow([band_of[x], repr(x), x.hex(), 33, 33])

    cells = [(band_of[x], n, x, value) for n, x, value in committed if x in args]
    assert len(cells) == 2409, len(cells)
    with open(here / "r4-cells.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "n", "x", "value_45_digit"])
        for band, n, x, value in cells:
            writer.writerow([band, n, repr(x), mp.nstr(value, 45)])

    # verification against the deposited 60-digit reference
    ref = read_reference(form_ref_path)
    worst = None
    where = None
    checked = 0
    for n, x, value in committed:
        if x not in args:
            continue
        assert x in ref, x
        diff = abs(ref[x][n] - value)
        checked += 1
        if worst is None or diff > worst:
            worst = diff
            where = (n, x)
    assert checked == 2409, checked
    expected = mpf("4.975493e-46")
    assert abs(worst - expected) / expected < mpf("2e-6"), (worst, expected)
    assert where == (2, 3.162277660168379e-06), where

    # ---- the 9,353rd first-sweep argument ---------------------------------
    designed = read_grid(designed_grid_path)
    deposited = read_grid(deposit_grid_path)
    designed_9352 = read_grid(designed_9352_path)
    assert len(designed) == 9353, len(designed)
    assert len(deposited) == 9353, len(deposited)
    assert len(designed_9352) == 9352, len(designed_9352)
    assert [r["band"] for r in designed] == [r["band"] for r in deposited]
    assert [r["x"] for r in designed] == [r["x"] for r in deposited]
    # the deposited 9,352-argument designed record is the prefix of the 9,353 grid
    assert [r["band"] for r in designed_9352] == [r["band"] for r in designed[:9352]]
    assert [r["x"] for r in designed_9352] == [r["x"] for r in designed[:9352]]
    assert designed[-1]["sets"] == "recovered"
    last = designed[-1]
    assert last["band"] == "A2"
    assert last["sets"] == "recovered"
    x = float(last["x"])
    assert x == RECOVERED, (x, RECOVERED)
    assert math.nextafter(A1_EDGE, math.inf) == RECOVERED
    assert x.hex() == "0x1.7ccb8e5a1580ep+2", x.hex()
    assert designed[-2]["x"] != last["x"]
    assert sum(1 for r in designed if r["sets"] == "recovered") == 1

    with open(here / "argument-9353.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["position", "band", "x", "x_hex", "sets", "rule"])
        writer.writerow([9353, last["band"], repr(x), x.hex(), last["sets"],
                         "nextafter(kRegionA1Edge, +inf); kRegionA1Edge = %r" % A1_EDGE])

    # its 33 reference values and cells, from the deposited 9,353 record
    ref35 = read_reference(paper_tree / "sweeps" / "amplification_form" / "data" / "reference.csv")
    # the first-sweep reference is the 60-digit file beside the grid
    first_ref_path = deposit_grid_path.parent / "reference.csv"
    first_ref = read_reference(first_ref_path)
    assert x in first_ref, "recovered argument absent from the first-sweep reference"
    first_cells_path = deposit_grid_path.parent.parent / "results" / "cells.csv"
    recovered_cells = [r for r in csv.DictReader(open(first_cells_path, newline=""))
                       if float(r["x"]) == x]
    assert len(recovered_cells) == 33, len(recovered_cells)

    # ---- digests ----------------------------------------------------------
    lines = ["sha256 digests of the inputs and the deposited lists", ""]
    for label, path in [
        ("kernel_committed_45_digit_grid", committed_path),
        ("second_sweep_grid", form_grid_path),
        ("second_sweep_60_digit_reference", form_ref_path),
        ("first_sweep_9353_grid", designed_grid_path),
        ("first_sweep_deposited_9353_grid", deposit_grid_path),
        ("first_sweep_deposited_9352_grid", designed_9352_path),
        ("r4_arguments", here / "r4-arguments.csv"),
        ("r4_cells", here / "r4-cells.csv"),
        ("argument_9353", here / "argument-9353.csv"),
    ]:
        lines.append("%-40s %s" % (label, sha256(path)))
    (here / "digests.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ---- verification log -------------------------------------------------
    report = [
        "R4 committed-grid argument list",
        "--------------------------------",
        "source: tests/data/boys_reference.csv at boys v5.3.1",
        "arguments in the committed file: %d" % len(all_x),
        "arguments below kRegionAEnd = %r: %d" % (A2_EDGE, len(args)),
        "orders per argument: 33; cells: %d" % checked,
        "present in the second sweep's 836-argument grid: %d of %d" %
        (len(args) - len(missing_band), len(args)),
        "worst |60-digit reference - committed value|: %s" % mp.nstr(worst, 10),
        "at (n, x) = (%d, %r); deposited summary states 4.975493e-46" % where,
        "",
        "The 9,353rd first-sweep argument",
        "---------------------------------",
        "deposited designed grid rows: %d" % (len(designed_9352)),
        "deposited 9,353-argument grid rows: %d" % (len(deposited)),
        "paper first-sweep grid rows: %d" % (len(designed)),
        "the 9,352-argument designed record is the 9,353 grid's first 9,352 rows: verified",
        "last row: band %s, x = %r = %s, tag %s" %
        (last["band"], x, x.hex(), last["sets"]),
        "nextafter(kRegionA1Edge, +inf) = %r: verified" % RECOVERED,
        "recovered reference row present: 33 orders; 33 cells in the deposited cells record",
        "",
        "Result: both lists reconstructed and verified.",
    ]
    (here / "verification.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
