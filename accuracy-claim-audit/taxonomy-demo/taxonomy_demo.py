#!/usr/bin/env python3
"""Reproduce the two synthetic demonstrations of the row taxonomy.

The demonstrations are constructed, not found: they show that a row's
property type decides which defects it can see.  A fast path is given a
planted drift that stays inside the lane's accuracy bound while exceeding
the per-band bound and the two-path consistency allowance.

The model is deliberately small and deterministic:

* ``F0`` is the Boys function of order zero, evaluated from ``math.erf``.
* ``scalar`` and ``fast`` are two evaluation paths.  Both carry the same
  baseline error of 9e-14 * F0(x) on [9, 10] (a pre-existing, accepted
  error inside the 1e-13 lane bound), so they agree everywhere before the
  plant.  The plant adds 5e-15 * F0(x) on [2, 3] to ``fast`` only.
* The lane accuracy row bounds the worst error of either path against the
  reference by 1e-13 * F0(x).
* The band row bounds the fast path on [2, 3] by 1e-15 * F0(x).
* The consistency row bounds |fast - scalar| by the sum of the two paths'
  half-lane bounds, 2e-15 * F0(x).
* The exactness row checks F0(0) == 1 on both paths.
* The ordinary test checks only the lane accuracy bound.

Panel A reports the accuracy row, the ordinary test and the consistency
row.  Panel B reports the four-row testbed: the band row, the lane
accuracy row, the consistency row and the exactness row.

Run: python taxonomy_demo.py [--json out.json]
Exit status 0 when the printed verdicts match the expected pattern.
"""

from __future__ import annotations

import argparse
import json
import math
import sys

RIGHT = 11.899848152108484
LANE_BOUND = 1e-13
BAND_BOUND = 1e-15
HALF_LANE_BOUND = 1e-15
ALLOWANCE = 2.0 * HALF_LANE_BOUND
BASE_ERROR = 9e-14
BASE_REGION = (9.0, 10.0)
DEFECT_REGION = (2.0, 3.0)
DRIFT = 5e-15

# A designed grid: the domain ends, the two regions' ends, and 4093 interior
# points.  The regions are covered densely enough that their maxima are
# evaluated, and every value is a binary64 double.
GRID = [0.0, RIGHT]
GRID += [RIGHT * i / 4096.0 for i in range(1, 4096)]
for r in (BASE_REGION, DEFECT_REGION):
    GRID += [r[0], r[1], 0.5 * (r[0] + r[1])]
GRID = sorted(set(GRID))


def f0(x: float) -> float:
    """The Boys function F_0 at x >= 0, with F_0(0) = 1."""
    if x == 0.0:
        return 1.0
    return math.sqrt(math.pi) / (2.0 * math.sqrt(x)) * math.erf(math.sqrt(x))


def baseline(x: float) -> float:
    """The shared baseline error, inside the lane bound but above the band bound."""
    return BASE_ERROR * f0(x) if BASE_REGION[0] <= x <= BASE_REGION[1] else 0.0


def _in_domain(x: float) -> None:
    if not (0.0 <= x <= RIGHT):
        raise ValueError(f"x outside the fitted domain: {x!r}")


def scalar(x: float) -> float:
    _in_domain(x)
    return f0(x) + baseline(x)


def fast(x: float, planted: bool) -> float:
    _in_domain(x)
    drift = DRIFT * f0(x) if planted and DEFECT_REGION[0] <= x <= DEFECT_REGION[1] else 0.0
    return scalar(x) + drift


def max_abs(values) -> float:
    return max(values, default=0.0)


def lane_row(planted: bool) -> dict:
    """The lane accuracy row: worst |path - reference| / F0 over both paths."""
    ratios = []
    for x in GRID:
        ref = f0(x)
        if ref == 0.0:
            continue
        ratios.append(abs(scalar(x) - ref) / ref)
        ratios.append(abs(fast(x, planted) - ref) / ref)
    worst = max_abs(ratios)
    return {"row": "lane.accuracy", "bound": LANE_BOUND, "worst": worst,
            "verdict": "pass" if worst <= LANE_BOUND else "fail"}


def band_row(planted: bool) -> dict:
    """The band row: worst |fast - reference| / F0 on [2, 3], bound 1e-15."""
    ratios = []
    for x in GRID:
        if DEFECT_REGION[0] <= x <= DEFECT_REGION[1]:
            ref = f0(x)
            if ref:
                ratios.append(abs(fast(x, planted) - ref) / ref)
    worst = max_abs(ratios)
    return {"row": "band.bound", "bound": BAND_BOUND, "worst": worst,
            "verdict": "pass" if worst <= BAND_BOUND else "fail"}


def consistency_row(planted: bool) -> dict:
    """The consistency row: worst |fast - scalar| / F0, bound 2e-15."""
    ratios = []
    for x in GRID:
        ref = f0(x)
        if ref:
            ratios.append(abs(fast(x, planted) - scalar(x)) / ref)
    worst = max_abs(ratios)
    return {"row": "path.consistency", "bound": ALLOWANCE, "worst": worst,
            "verdict": "pass" if worst <= ALLOWANCE else "fail"}


def exactness_row(planted: bool) -> dict:
    """F0(0) == 1 on both paths, exactly."""
    ok = scalar(0.0) == 1.0 and fast(0.0, planted) == 1.0
    return {"row": "exact.zero", "bound": 0.0, "worst": 0.0 if ok else 1.0,
            "verdict": "pass" if ok else "fail"}


def ordinary_test(planted: bool) -> dict:
    """An ordinary test checks only the lane's accuracy bound."""
    row = lane_row(planted)
    return {"row": "ordinary.test", "bound": LANE_BOUND, "worst": row["worst"],
            "verdict": row["verdict"]}


def panel(planted: bool) -> dict:
    return {
        "planted": planted,
        "lane": lane_row(planted),
        "band": band_row(planted),
        "consistency": consistency_row(planted),
        "exactness": exactness_row(planted),
        "ordinary": ordinary_test(planted),
    }


def expected() -> dict:
    return {
        "before": {"lane": "pass", "band": "pass", "consistency": "pass",
                   "exactness": "pass", "ordinary": "pass"},
        "after": {"lane": "pass", "band": "fail", "consistency": "fail",
                  "exactness": "pass", "ordinary": "pass"},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)
    before = panel(False)
    after = panel(True)
    result = {
        "constants": {
            "right": RIGHT, "lane_bound": LANE_BOUND, "band_bound": BAND_BOUND,
            "half_lane_bound": HALF_LANE_BOUND, "allowance": ALLOWANCE,
            "base_error": BASE_ERROR, "base_region": BASE_REGION,
            "defect_region": DEFECT_REGION, "drift": DRIFT,
            "grid_points": len(GRID),
        },
        "before": before, "after": after,
    }
    payload = {"schema": "taxonomy-demo-1", "result": result}
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1)
            handle.write("\n")

    exp = expected()

    def line(phase, key, row):
        print(f"{phase:6s} {row['row']:18s} bound={row['bound']:.3e} "
              f"worst={row['worst']:.6e} {row['verdict'].upper()}")

    print(f"grid: {len(GRID)} points on [0, {RIGHT!r}]")
    print(f"drift: {DRIFT:.3e} * F0 on {DEFECT_REGION}; lane bound {LANE_BOUND:.3e}; "
          f"band bound {BAND_BOUND:.3e}; consistency allowance {ALLOWANCE:.3e}")
    for phase, panel_data in (("before", before), ("after", after)):
        for key in ("lane", "band", "consistency", "exactness", "ordinary"):
            line(phase, key, panel_data[key])
        print()

    failures = []
    for phase, panel_data in (("before", before), ("after", after)):
        for key, want in exp[phase].items():
            got = panel_data[key]["verdict"]
            if got != want:
                failures.append(f"{phase} {key}: {got} != {want}")
    if failures:
        print("MISMATCH: " + "; ".join(failures))
        return 1
    print("VERDICTS MATCH THE TWO DEMONSTRATIONS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
