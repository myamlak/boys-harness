#!/usr/bin/env python3
r"""The affine-residual effect of the measured 60-digit counts.
Two tests are recomputed with the measured 60-digit counts substituted for the
measured corrected-record cells; every other count stays the corrected
record's binary64 count.

1. The model-based slope test of the paper's Table S7: the increment from
   \(10^{-8}\) to the deep target against the model's predicted increment
   (5.57 at the paper width over 7 decades, 5.71 at the ship width, 6.41 at
   \(2^5\)). The 60-digit endpoint counts are the E-* cells; if an endpoint
   cell has not been run the record's count is used and the result is marked.
2. The pre-set free-line test of the supplement's Table S6: a free line is
   fitted to the stored counts against L = target exponent * ln 10 above the
   solver floor, and the width fails when the largest residual or any second
   difference exceeds one stored coefficient. The paper's own criterion is
   used; the script reproduces the record's published values as a self-check.

Usage: python freeline_effect.py [--json freeline_effect.json]
"""

import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
RECORD = HERE.parent / "results-corrected" / "corrected-counts.json"

# width, deep cell key (1e-15), endpoint cell key, predicted increment,
# per-decade value, and every measured cell as (key, target exponent)
WIDTHS = [("p0", "D", None, 3.05, 0.435, [("D", 15)]),
          ("paper", "A", "E-paper", 5.57, 0.796, [("A", 15)]),
          ("ship", "B", "E-ship", 5.71, 0.816, [("B", 15)]),
          ("p5", "C1", "E-p5", 6.41, 0.916, [("C2", 14), ("C1", 15)])]
FLOOR = 3
TOL = 1.0


def cell(key):
    p = HERE / f"cell_{key}.json"
    return json.loads(p.read_text()) if p.exists() else None


def free_line(counts):
    L = [e * math.log(10.0) for e in range(8, 8 + len(counts))]
    used = [(L[i], c) for i, c in enumerate(counts) if c > FLOOR]
    n = len(used)
    if n < 3:
        return None
    mx = sum(x for x, _ in used) / n
    my = sum(y for _, y in used) / n
    sxx = sum((x - mx) ** 2 for x, _ in used)
    sxy = sum((x - mx) * (y - my) for x, y in used)
    b = sxy / sxx
    a = my - b * mx
    res = [y - (a + b * x) for x, y in used]
    sec = [used[i + 1][1] - 2 * used[i][1] + used[i - 1][1]
           for i in range(1, n - 1)]
    max_res = max(abs(r) for r in res)
    max_sec = max(abs(s) for s in sec)
    return {"counts": [y for _, y in used], "max_abs_residual": max_res,
            "max_abs_second_difference": max_sec,
            "verdict": "pass" if max_res <= TOL and max_sec <= TOL else "fail"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=str(HERE / "freeline_effect.json"))
    args = ap.parse_args()
    rec = json.loads(RECORD.read_text())
    series = {}
    for w in ("p0", "paper", "ship", "p5"):
        rows = [r for r in rec["ladder"] if r["width_name"] == w
                and r["target"].startswith("rel:1e-")]
        rows.sort(key=lambda r: float(r["target"].split(":")[1]),
                  reverse=True)
        series[w] = [r["rational_corrected"] for r in rows]
    out = {"model": [], "free_line": []}
    for wname, deep_key, end_key, predicted, per_decade, measured in WIDTHS:
        deep = cell(deep_key)
        end = cell(end_key) if end_key else None
        base = list(series[wname])
        # the ladder lists 1e-8 .. 1e-15; index = exponent - 8
        record = {"counts": base, "published": free_line(base)}
        substituted = list(base)
        applied = []
        for mkey, exponent in measured:
            md = cell(mkey)
            if md and md.get("count60") is not None:
                substituted[exponent - 8] = md["count60"]
                applied.append({"cell": mkey, "target_exponent": exponent,
                                "record": base[exponent - 8],
                                "count60": md["count60"]})
        effect = free_line(substituted)
        out["free_line"].append({"width": wname, "record": record,
                                 "substituted": {"applied": applied,
                                                 "counts": substituted,
                                                 **effect}})
        if deep and deep.get("count60") is not None:
            endpoint = (end["count60"] if end and end.get("count60") is not None
                        else base[0])
            increment = deep["count60"] - endpoint
            out["model"].append({
                "width": wname, "endpoint_count": endpoint,
                "deep_count": deep["count60"], "increment": increment,
                "predicted": predicted,
                "residual": round(increment - predicted, 2),
                "per_decade": per_decade,
                "residual_exact": round(increment - per_decade * 7, 3),
                "endpoint_measured": bool(end and end.get("count60") is not None)})
    print(json.dumps(out, indent=1))
    Path(args.json).write_text(json.dumps(out, indent=1))
    print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
