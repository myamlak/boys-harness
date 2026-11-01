#!/usr/bin/env python3
"""Print the F1 report's markdown tables from a sweep record.

Reads a sweep record (the delivered one or the corrected copy) and prints:

  * the width table (W at x0 and a2);
  * the Chebyshev stored counts (x0/a2);
  * the rational stored counts (x0/a2);
  * the accepted rational pairs at x0 and at a2, with r/l solver letters and
    a '*' on any row whose recorded below check met its target;
  * the per-target second-edge count comparison.

Usage:
    python _knee-f1/make_report_tables.py --sweep <record>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
WIDTHS = [f"2^{j}" for j in range(-7, 6)] + ["ship", "paper"]
TARGETS = ("rel:1e-13", "abs:5e-14", "rel:1e-14")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep",
                    default=str(PAPER / "results-second-instance-f1"
                                / "sweep_all.json"))
    args = ap.parse_args()
    doc = json.loads(Path(args.sweep).read_text(encoding="utf-8"))
    rows = {(r["edge"], r["target"], r["width_name"]): r
            for r in doc["rows"]}

    print("## widths\n")
    print("| width | W x0 | W a2 |")
    print("| --- | --- | --- |")
    for name in WIDTHS:
        r0 = rows[("x0", "rel:1e-13", name)]
        r2 = rows[("a2", "rel:1e-13", name)]
        print(f"| {name} | {r0['width']} | {r2['width']} |")

    print("\n## chebyshev\n")
    print("| width | rel:1e-13 | abs:5e-14 | rel:1e-14 |")
    print("| --- | --- | --- | --- |")
    for name in WIDTHS:
        cells = []
        for t in TARGETS:
            c0 = rows[("x0", t, name)]["cheb"]["stored"]
            c2 = rows[("a2", t, name)]["cheb"]["stored"]
            cells.append(f"{c0}/{c2}")
        print(f"| {name} | " + " | ".join(cells) + " |")

    print("\n## rational\n")
    print("| width | rel:1e-13 | abs:5e-14 | rel:1e-14 |")
    print("| --- | --- | --- | --- |")
    for name in WIDTHS:
        cells = []
        for t in TARGETS:
            a = rows[("x0", t, name)]["rat"]
            b = rows[("a2", t, name)]["rat"]
            mark = "" if not (a.get("below_meets") or b.get("below_meets")) \
                else "*"
            cells.append(f"{a['stored']}/{b['stored']}{mark}")
        print(f"| {name} | " + " | ".join(cells) + " |")

    print("\n## pairs\n")
    for edge in ("x0", "a2"):
        print(f"\n{edge}:\n")
        print("```")
        for name in WIDTHS:
            cells = []
            for t in TARGETS:
                e = rows[(edge, t, name)]["rat"]
                letter = {"remez": "r", "lawson": "l"}[e["solver"]]
                mark = "*" if e.get("below_meets") else ""
                cells.append(f"[{e['m']}/{e['k']}]{letter}{mark}")
            pad = " " * (7 - len(name))
            print(f"{name}{pad}" + "  ".join(cells))
        print("```")

    print("\n## second edge\n")
    print("| target | family | width | x0 | a2 | diff |")
    print("| --- | --- | --- | --- | --- | --- |")
    for t in TARGETS:
        for fam in ("cheb", "rat"):
            for name in WIDTHS:
                r0 = rows[("x0", t, name)][fam]
                r2 = rows[("a2", t, name)][fam]
                d = r2["stored"] - r0["stored"]
                if d:
                    print(f"| {t} | {fam} | {name} | {r0['stored']} | "
                          f"{r2['stored']} | {d:+d} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
