#!/usr/bin/env python3
"""Apply the pre-set C-LAW test to the dense target-ladder records.

For each width and family, above the solver floor: fit the stored count against
L = log(1/epsilon), compare the measured slope with the model's 1/ell or 1/mu,
and report the largest residual and the largest second difference. The verdict
rests on the residual and the second difference; the slope comparison is
reported beside it. The portfolio threshold: the counts are affine in L after
lattice rounding when the largest residual is at most one stored coefficient and
no second difference exceeds that tolerance. One width/family outside it is
C-NOLAW. A row the solver flags below_meets is known not to carry the minimal
count and is excluded from the fit and listed beside it. A record with any
width/family short of three above-floor targets is C-INCOMPLETE and decides
nothing.

Usage:
    python ladder_analysis.py [--dir results-ladder] [--out ladder_analysis.json]
"""

import argparse
import json
import math
from pathlib import Path

from mpmath import mp, mpf, sqrt, log, ellipk, pi

mp.dps = 40

A = mpf("11.899848152108484")
FLOOR_TOLERANCE = 1.0


def kappa(width):
    return (sqrt(A + width) - sqrt(A)) / (sqrt(A + width) + sqrt(A))


def ell_mu(width):
    k = kappa(width)
    m = k * k
    return log(1 / k), pi / 2 * ellipk(1 - m) / ellipk(m)


def parse_L(spec):
    # rel:1e-13 -> L = 13 ln 10
    assert spec.startswith("rel:"), spec
    exponent = -int(spec.rsplit("e", 1)[1])
    return exponent * math.log(10.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results-ladder")
    ap.add_argument("--out", default="ladder_analysis.json")
    args = ap.parse_args()

    rows = json.load(open(Path(args.dir) / "sweep_all.json"))["rows"]
    by_width = {}
    for r in rows:
        by_width.setdefault(round(float(r["width_f"]), 17), []).append(r)

    # Each family's floor is the smallest count that family ever stores; a row
    # where the other family failed is still a point of this one.
    floors = {}
    for fam in ("cheb", "rat"):
        values = [r[fam]["stored"] for r in rows if r.get(fam, {}).get("stored") is not None]
        floors[fam] = min(values) if values else None

    report = {"floor": floors, "widths": {}}
    flips = []
    if not rows:
        report.update({"flips": [], "verdict": "C-INCOMPLETE", "complete": False,
                       "note": "the record has no rows"})
        Path(args.out).write_text(json.dumps(report, indent=1) + "\n")
        print(json.dumps(report, indent=1))
        return
    for width_f in sorted(by_width):
        entries = sorted(by_width[width_f], key=lambda r: parse_L(r["target"]))
        ell, mu = ell_mu(mpf(repr(width_f)))
        family_report = {}
        for fam, rate in (("cheb", ell), ("rat", mu)):
            fam_entries = [r for r in entries if r.get(fam, {}).get("stored") is not None]
            # A row the solver flagged below_meets returns a count known not to be
            # minimal; it is reported and kept out of the fit.
            excluded = [
                {"target": r["target"], "stored": r[fam]["stored"]}
                for r in fam_entries if r[fam].get("below_meets")
            ]
            minimal = [r for r in fam_entries if not r[fam].get("below_meets")]
            pts = [(parse_L(r["target"]), r[fam]["stored"]) for r in minimal]
            used = [(L, n) for L, n in pts if floors[fam] is not None and n > floors[fam]]
            entry = {
                "points_total": len(pts),
                "points_above_floor": len(used),
                "counts": [n for _, n in used],
                "Ls": [round(L, 6) for L, _ in used],
                "excluded_non_minimal": excluded,
            }
            if len(used) >= 3:
                n = len(used)
                mx = sum(L for L, _ in used) / n
                my = sum(c for _, c in used) / n
                sxx = sum((L - mx) ** 2 for L, _ in used)
                sxy = sum((L - mx) * (c - my) for L, c in used)
                slope = sxy / sxx
                intercept = my - slope * mx
                residuals = [c - (intercept + slope * L) for L, c in used]
                second = [
                    used[i + 1][1] - 2 * used[i][1] + used[i - 1][1]
                    for i in range(1, n - 1)
                ]
                slope_se = math.sqrt(
                    sum(r * r for r in residuals) / (n - 2) / sxx
                ) if n > 2 else float("nan")
                predicted = float(1 / rate)
                entry.update(
                    {
                        "slope_measured": slope,
                        "slope_predicted": predicted,
                        "slope_se": slope_se,
                        "slope_sigma": abs(slope - predicted) / slope_se
                        if slope_se else None,
                        "intercept": intercept,
                        "max_abs_residual": max(abs(r) for r in residuals),
                        "max_abs_second_difference": max(abs(s) for s in second)
                        if second
                        else 0.0,
                    }
                )
                ok = (
                    entry["max_abs_residual"] <= FLOOR_TOLERANCE
                    and entry["max_abs_second_difference"] <= FLOOR_TOLERANCE
                )
                entry["affine_within_tolerance"] = ok
                if not ok:
                    flips.append({"width": width_f, "family": fam})
            else:
                entry["affine_within_tolerance"] = None
                entry["note"] = "fewer than three above-floor targets"
            family_report[fam] = entry
        report["widths"][repr(width_f)] = family_report

    report["verdict"] = "C-INCOMPLETE"
    report["flips"] = flips
    report["complete"] = all(
        f["affine_within_tolerance"] is not None
        for w in report["widths"].values()
        for f in w.values()
    )
    if report["complete"]:
        report["verdict"] = "C-NOLAW" if flips else "C-LAW"
    Path(args.out).write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
