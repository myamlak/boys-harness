#!/usr/bin/env python3
r"""Re-run the paper's own free-line test on the record and the measured counts.

The paper's pre-set free-line test is implemented by
`../analysis/ladder_analysis.py` (the paper tree's copy, byte-identical,
cited in supplement section S10). It reads a kernel-format `sweep_all.json`:
one row per (width, target) with a `rat` (and `cheb`) block carrying
`stored`.

This script builds that record from the corrected counts: rational counts from
`../results-corrected/corrected-counts.json`, Chebyshev counts from the
supplement's Table S5 (hardcoded below). It then runs the paper's
`ladder_analysis.py` on the corrected record and on the record with the
control's measured 60-digit counts substituted:

* paper \(10^{-15}\) = cell A's count,
* ship \(10^{-15}\) = cell B's count,
* \(2^5\) \(10^{-14}\) and \(10^{-15}\) = cells C2's and C1's counts,
* \(2^0\) \(10^{-15}\) = cell D's count.

The corrected record must reproduce supplement Table S6 before the substituted
runs are read. A case without its cell record is skipped and marked.

Usage: python verify_freeline.py
"""

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
SWEEP = PAPER / "records" / "ladder" / "sweep_all.json"
CORRECTED = PAPER / "results-corrected" / "corrected-counts.json"
LADDER_ANALYSIS = PAPER / "analysis" / "ladder_analysis.py"
SCRATCH = HERE / "_freeline-scratch"

SHIP_WIDTH_F = 19.100151847891517
PAPER_WIDTH_F = 17.089489586712254
P5_WIDTH_F = 32.0
W1_WIDTH_F = 1.0

NAME_MAP = {"m7": "2^-7", "p0": "2^0", "p5": "2^5",
            "ship": "ship", "paper": "paper"}
# Supplement Table S5, corrected counts, Chebyshev column, targets 1e-8..1e-15.
CHEB = {
    "2^-7": [3, 3, 3, 3, 5, 5, 5, 5],
    "2^0": [5, 6, 6, 7, 7, 8, 9, 9],
    "2^5": [16, 17, 19, 21, 23, 25, 26, 28],
    "ship": [13, 14, 16, 17, 18, 20, 21, 23],
    "paper": [12, 13, 15, 16, 18, 19, 20, 22],
}
TARGETS = ["rel:1e-08", "rel:1e-09", "rel:1e-10", "rel:1e-11",
           "rel:1e-12", "rel:1e-13", "rel:1e-14", "rel:1e-15"]


def cell_count(key):
    p = HERE / f"cell_{key}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text()).get("count60")


def build(subs):
    """subs: {(width_name, target): count} in the sweep record's names."""
    sweep = json.loads(SWEEP.read_text())
    corr = json.loads(CORRECTED.read_text())
    rat = {(NAME_MAP[r["width_name"]], r["target"]): r["rational_corrected"]
           for r in corr["ladder"]}
    for row in sweep["rows"]:
        key = (row["width_name"], row["target"])
        row["rat"]["stored"] = rat[key]
        row["rat"]["below_meets"] = False
        row["cheb"]["stored"] = CHEB[row["width_name"]][
            TARGETS.index(row["target"])]
        row["cheb"]["below_meets"] = False
        if key in subs:
            row["rat"]["stored"] = subs[key]
    return sweep


def run(name, subs):
    d = SCRATCH / name
    d.mkdir(parents=True, exist_ok=True)
    (d / "sweep_all.json").write_text(json.dumps(build(subs), indent=1))
    out = d / "analysis.json"
    proc = subprocess.run(
        [sys.executable, str(LADDER_ANALYSIS), "--dir", str(d),
         "--out", str(out)], capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"{name}: ladder_analysis.py exit {proc.returncode}")
        print(proc.stdout[-2000:])
        print(proc.stderr[-2000:])
        raise SystemExit(1)
    doc = json.loads(out.read_text())
    print(f"{name}: ladder_analysis.py exit 0; verdict {doc['verdict']}; "
          f"flips {[(f['width'], f['family']) for f in doc['flips']]}")
    for width in (W1_WIDTH_F, PAPER_WIDTH_F, SHIP_WIDTH_F, P5_WIDTH_F):
        e = doc["widths"][repr(width)]["rat"]
        print(f"  {width!r}: counts {e['counts']}; "
              f"residual {e['max_abs_residual']!r}; "
              f"2nd diff {e['max_abs_second_difference']!r}; "
              f"affine {e['affine_within_tolerance']}")
    return doc


def main():
    c1, c2 = cell_count("C1"), cell_count("C2")
    cases = [("record", {})]
    if cell_count("A") is not None:
        cases.append(("paper15", {("paper", "rel:1e-15"): cell_count("A")}))
    if cell_count("B") is not None:
        cases.append(("ship16", {("ship", "rel:1e-15"): cell_count("B")}))
    if cell_count("D") is not None:
        cases.append(("p0_8", {("2^0", "rel:1e-15"): cell_count("D")}))
    if c1 is not None and c2 is not None:
        cases.append(("p5", {("2^5", "rel:1e-14"): c2,
                             ("2^5", "rel:1e-15"): c1}))
    docs = {}
    for name, subs in cases:
        docs[name] = run(name, subs)
    # the paper's Table S6 values, as a self-check on the constructed record
    rec = docs["record"]
    expect = {
        W1_WIDTH_F: (0.8333333333333321, 2, False),
        PAPER_WIDTH_F: (1.75, 2, False),
        SHIP_WIDTH_F: (1.1666666666666679, 1, False),
        P5_WIDTH_F: (3.8333333333333286, 4, False),
    }
    for width, (res, sec, ok) in expect.items():
        e = rec["widths"][repr(width)]["rat"]
        assert e["max_abs_residual"] == res, (width, e)
        assert e["max_abs_second_difference"] == sec, (width, e)
        assert e["affine_within_tolerance"] is ok, (width, e)
    print("record reproduces supplement Table S6; "
          f"{len(cases) - 1} substituted case(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
