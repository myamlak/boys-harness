#!/usr/bin/env python3
"""Apply the pre-registered clean-edge decision rule to edge62_sweep.json.

Reads the ladder record, fits the stored count against -log10(target) on a free
line for each width and family, and reports the residual and second difference.
Primary range: all eight targets (1e-8 .. 1e-15). Secondary: the paper's contract
range (1e-8 .. 1e-13).

The rule (C-clean-edge-preregistration.md):
  Chebyshev holds at a width if max |residual| <= 1 stored coefficient.
  Rational fails at a width if max |residual| > 2 or max |second difference| >= 2.
  Clean-edge failure confirmed if rational fails at >= 3 of 5 widths while
  Chebyshev holds everywhere; transient-specific if rational holds at >= 4 of 5;
  otherwise mixed.

No third-party imports. Writes edge62_analysis.json and prints a readable table.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RECORD = HERE / "edge62_sweep.json"
OUT = HERE / "edge62_analysis.json"

TARGETS = ["rel:1e-08", "rel:1e-09", "rel:1e-10", "rel:1e-11",
           "rel:1e-12", "rel:1e-13", "rel:1e-14", "rel:1e-15"]
WIDTHS = ["2^-7", "2^0", "2^5", "ship", "paper"]
EXPONENTS = [8, 9, 10, 11, 12, 13, 14, 15]
CONTRACT = 6  # first six targets, 1e-8 .. 1e-13


def fit_line(xs: list[int], ys: list[float]) -> tuple[float, float]:
    n = len(xs)
    xbar = math.fsum(xs) / n
    ybar = math.fsum(ys) / n
    sxx = math.fsum((x - xbar) ** 2 for x in xs)
    sxy = math.fsum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = ybar - slope * xbar
    return intercept, slope


def metrics(ys: list[float], count: int) -> dict:
    xs = EXPONENTS[:count]
    ys = ys[:count]
    intercept, slope = fit_line(xs, ys)
    residual = max(abs(y - (intercept + slope * x)) for x, y in zip(xs, ys))
    second = max((abs(ys[i - 1] - 2 * ys[i] + ys[i + 1])
                  for i in range(1, len(ys) - 1)), default=0.0)
    return {"n": count, "intercept": intercept, "slope": slope,
            "max_abs_residual": residual, "max_abs_second_difference": second}


def family_verdict(stats: dict, family: str) -> str:
    if family == "cheb":
        return "holds" if stats["max_abs_residual"] <= 1.0 else "fails"
    if stats["max_abs_residual"] > 2.0 or stats["max_abs_second_difference"] >= 2.0:
        return "fails"
    return "holds"


def main() -> int:
    if not RECORD.exists():
        print(f"missing {RECORD}", file=sys.stderr)
        return 1
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    rows = record.get("rows", [])
    by_key = {(r["target"], r["width_name"]): r for r in rows}
    missing = [(t, w) for t in TARGETS for w in WIDTHS if (t, w) not in by_key]
    if missing:
        print(f"incomplete: {len(rows)}/40 rows; missing {missing}", file=sys.stderr)
        return 2

    analysis: dict = {"edge": record.get("meta", {}).get("edge"),
                      "scale": record.get("meta", {}).get("scale"),
                      "widths": {}, "missing": []}
    for width in WIDTHS:
        entry: dict = {}
        for family in ("cheb", "rat"):
            ys = [float(by_key[(t, width)][family]["stored"]) for t in TARGETS]
            primary = metrics(ys, len(TARGETS))
            secondary = metrics(ys, CONTRACT)
            entry[family] = {
                "stored": ys,
                "primary_all_eight": primary,
                "secondary_contract": secondary,
                "verdict_primary": family_verdict(primary, family),
                "verdict_secondary": family_verdict(secondary, family),
            }
        analysis["widths"][width] = entry

    for key in ("primary_all_eight", "secondary_contract"):
        suffix = "primary" if key == "primary_all_eight" else "secondary"
        cheb_ok = all(analysis["widths"][w]["cheb"][f"verdict_{suffix}"] == "holds"
                      for w in WIDTHS)
        rat_fail = [w for w in WIDTHS
                    if analysis["widths"][w]["rat"][f"verdict_{suffix}"] == "fails"]
        rat_hold = [w for w in WIDTHS
                    if analysis["widths"][w]["rat"][f"verdict_{suffix}"] == "holds"]
        if not cheb_ok:
            verdict = "uninterpretable-chebyshev"
        elif len(rat_fail) >= 3:
            verdict = "clean-edge-failure-confirmed"
        elif len(rat_hold) >= 4:
            verdict = "transient-specific"
        else:
            verdict = "mixed"
        analysis[f"verdict_{suffix}"] = {
            "chebyshev_holds_everywhere": cheb_ok,
            "rational_fail_widths": rat_fail,
            "rational_hold_widths": rat_hold,
            "verdict": verdict,
        }

    OUT.write_text(json.dumps(analysis, indent=2) + "\n", encoding="utf-8")
    for width in WIDTHS:
        e = analysis["widths"][width]
        print(f"{width:6s} cheb r={e['cheb']['primary_all_eight']['max_abs_residual']:8.3f} "
              f"d2={e['cheb']['primary_all_eight']['max_abs_second_difference']:5.1f} "
              f"{e['cheb']['verdict_primary']:5s} | "
              f"rat r={e['rat']['primary_all_eight']['max_abs_residual']:8.3f} "
              f"d2={e['rat']['primary_all_eight']['max_abs_second_difference']:5.1f} "
              f"{e['rat']['verdict_primary']:5s}")
    print()
    print("primary verdict:", analysis["verdict_primary"]["verdict"],
          "| fail widths:", analysis["verdict_primary"]["rational_fail_widths"])
    print("secondary verdict:", analysis["verdict_secondary"]["verdict"],
          "| fail widths:", analysis["verdict_secondary"]["rational_fail_widths"])
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
