#!/usr/bin/env python3
"""Compare the 60-digit control with the record, and run the free-line scenario.

Two reports, both read from the control record
(results-ladder/rational_unrounded_control.json) and the paper's own files:

  --compare   the decisive cells side by side: the record's stored count,
              pair, solver and binary64 delivered error against the 60-digit
              count, pair and error, with the best miss at the failing level;
              then the three model-based increments and their residuals.

  --scenario  the paper's free-line test re-run with the measured 60-digit
              counts substituted for the measured rows. Every other row keeps
              the record's count, so its verdict is a scenario, not a
              measurement, and it says what the pre-set test would return if
              the unmeasured cells were unchanged.

The control record has two revisions. Cells written before the resume carry
the accepted pair only inside their levels; the helpers below derive the same
fields from either shape.

Usage:
  python rational_unrounded_scenario.py --compare
  python rational_unrounded_scenario.py --scenario
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent
CONTROL = PAPER / "results-ladder" / "rational_unrounded_control.json"
SWEEP = PAPER / "results-ladder" / "sweep_all.json"
SCRATCH = PAPER / "_unrounded_scenario"

# The three model-based increments of the paper's Table tab:ladderpred, with
# the model's predicted increment over the span and the number of decades.
SPANS = [
    ("paper", "rel:1e-08", "rel:1e-14", 6, 4.78),
    ("ship", "rel:1e-08", "rel:1e-15", 7, 5.71),
    ("2^5", "rel:1e-08", "rel:1e-15", 7, 6.41),
]

TWO_COEFFICIENT_STEP = 2.0


def load(path):
    return json.loads(path.read_text())


def derive_from_levels(cell):
    """The accepted pair and the failing level, from either record revision."""
    best = None
    fail = None
    for level in cell.get("levels", []):
        for p in level.get("passes") or [level]:
            for r in p.get("rows", []):
                if not r.get("meets60"):
                    continue
                key = (r["stored"], float(r["rel60"]))
                if best is None or key < best[0]:
                    best = (key, r)
        if level.get("passes") is not None and not level.get("meets60"):
            fail = level
    out = {}
    if best is not None:
        r = best[1]
        out.update({"count60_m": r["m"], "count60_k": r["k"],
                    "count60_solver": r["solver"],
                    "count60_pair": f"[{r['m']}/{r['k']}] {r['solver']}",
                    "count60_rel60": r["rel60"],
                    "count60_rel64": r["rel64"],
                    "count60_meets64": r["meets64"]})
    if fail is not None:
        out["failing_level"] = {
            "stored": fail["stored"],
            "routes": [{"pass": p["pass"], "best_rel60": p["best_rel60"],
                        "solved": p["solved"]} for p in fail["passes"]],
        }
    return out


def cell_fields(cell):
    """The fields the comparison needs, from either record revision."""
    out = {k: cell.get(k) for k in
           ("record", "count60", "direction", "count60_pair",
            "count60_solver", "count60_m", "count60_k", "count60_rel60",
            "count60_rel64", "count60_meets64", "failing_level")}
    derived = derive_from_levels(cell)
    for key, value in derived.items():
        if out.get(key) is None:
            out[key] = value
    return out


def by_cell(control):
    return {(c["width_name"], c["target"]): cell_fields(c)
            for c in control["cells"] if "error" not in c}


def fmt(value, width=14):
    return f"{'-' if value is None else value:>{width}}"


def compare(control):
    cells = by_cell(control)
    print("decisive cells: record (binary64) against control (60 digits)\n")
    header = ("width   target    rec  record pair        rec rel64        "
              "ctl  control pair      ctl rel60        ctl rel64    dir   "
              "fail best rel60")
    print(header)
    print("-" * len(header))
    width_order = {name: i for i, (name, *_rest) in enumerate(SPANS)}
    for (name, target), c in sorted(
            cells.items(),
            key=lambda kv: (width_order.get(kv[0][0], 99), kv[0][1])):
        rec = c.get("record") or {}
        fail = c.get("failing_level") or {}
        routes = "; ".join(f"{r['pass']} {r['best_rel60']}"
                           for r in fail.get("routes", []))
        print(f"{name:<7} {target:<9} "
              f"{fmt(rec.get('stored'), 4)} "
              f"[{rec.get('m')}/{rec.get('k')}] {str(rec.get('solver')):<8} "
              f"{str(rec.get('delivered_rel')):<16} "
              f"{fmt(c.get('count60'), 4)} "
              f"{str(c.get('count60_pair')):<16} "
              f"{str(c.get('count60_rel60')):<16} "
              f"{str(c.get('count60_rel64')):<16} "
              f"{str(c.get('direction')):<5} "
              f"{routes}")

    print("\nmodel-based increments (endpoint counts)\n")
    print("width  span              decades  model  record  ctl60  "
          "res64  res60")
    for name, start, end, decades, predicted in SPANS:
        a = cells.get((name, start))
        b = cells.get((name, end))
        rec_a = None if a is None else a["record"]["stored"]
        rec_b = None if b is None else b["record"]["stored"]
        ctl_a = None if a is None else a.get("count60")
        ctl_b = None if b is None else b.get("count60")
        rec_inc = None if None in (rec_a, rec_b) else rec_b - rec_a
        ctl_inc = None if None in (ctl_a, ctl_b) else ctl_b - ctl_a
        res64 = None if rec_inc is None else rec_inc - predicted
        res60 = None if ctl_inc is None else ctl_inc - predicted
        print(f"{name:<6} {start.split(':')[1]} -> {end.split(':')[1]}  "
              f"{decades:>7}  {predicted:>5.2f}  {fmt(rec_inc, 6)} "
              f"{fmt(ctl_inc, 6)}  {fmt(None if res64 is None else f'{res64:+.2f}', 5)} "
              f"{fmt(None if res60 is None else f'{res60:+.2f}', 5)}")
        if res60 is not None:
            beyond = abs(res60) - TWO_COEFFICIENT_STEP
            print(f"       control residual {res60:+.2f} is {beyond:+.2f} "
                  f"coefficients beyond one two-coefficient step "
                  f"({TWO_COEFFICIENT_STEP:.0f})")


def scenario(control):
    rows = load(SWEEP)["rows"]
    cells = by_cell(control)
    subs = {}
    for (name, target), c in cells.items():
        if c.get("count60") is None:
            continue
        subs[(name, target)] = {
            "count60": c["count60"], "record": c["record"]["stored"],
            "pair": c.get("count60_pair"), "rel60": c.get("count60_rel60"),
            "direction": c.get("direction")}
    substituted = 0
    for r in rows:
        key = (r["width_name"], r["target"])
        if key not in subs:
            continue
        s = subs[key]
        r["rat"]["stored"] = s["count60"]
        # The walk certified the count below as a miss, so the substituted
        # row is minimal under the same one-level check the record uses.
        r["rat"]["below_meets"] = False
        r["rat"]["unrounded_control"] = {
            "count60": s["count60"], "record_stored": s["record"],
            "pair60": s["pair"], "rel60": s["rel60"],
            "direction": s["direction"]}
        substituted += 1
    SCRATCH.mkdir(exist_ok=True)
    (SCRATCH / "sweep_all.json").write_text(json.dumps({"rows": rows}, indent=1))
    out = SCRATCH / "ladder_analysis.json"
    subprocess.run([sys.executable, str(PAPER / "ladder_analysis.py"),
                    "--dir", str(SCRATCH), "--out", str(out)], check=True,
                   capture_output=True, text=True)
    analysis = load(out)
    print(f"\nsubstituted {substituted} of {len(rows)} rational rows; "
          f"the rest keep the record's count")
    print(f"free-line verdict: {analysis['verdict']}")
    for width, entry in analysis["widths"].items():
        rat = entry.get("rat")
        if not rat:
            continue
        counts = rat.get("counts")
        print(f"width {width}: points above floor {rat['points_above_floor']}, "
              f"counts {counts}, max|residual| {rat['max_abs_residual']:.3f}, "
              f"max|second difference| {rat['max_abs_second_difference']}, "
              f"affine {rat['affine_within_tolerance']}")
        for ex in rat.get("excluded_non_minimal", []):
            print(f"  excluded non-minimal: {ex['target']} stored {ex['stored']}")
    totals = {}
    for r in rows:
        totals[r["width_name"]] = totals.get(r["width_name"], 0) + 1
    per_width = {}
    for name, _target in subs:
        per_width[name] = per_width.get(name, 0) + 1
    print("\nrows substituted by width (measured / total):")
    for name in sorted(totals):
        print(f"  {name:<8} {per_width.get(name, 0):>2} / {totals[name]}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--scenario", action="store_true")
    ap.add_argument("--control", default=str(CONTROL))
    args = ap.parse_args(argv)
    if not (args.compare or args.scenario):
        ap.error("name --compare or --scenario")
    control = load(Path(args.control))
    if args.compare:
        compare(control)
    if args.scenario:
        scenario(control)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
