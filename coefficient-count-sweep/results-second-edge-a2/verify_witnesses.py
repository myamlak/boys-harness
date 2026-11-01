#!/usr/bin/env python3
"""Re-evaluate every stored a2 witness from its coefficients.

This is a separate implementation of the paper's acceptance rule, not a call
into the driver: the grid is rebuilt from each record's own a and W, and the
witness coefficients are read from the record.  For both the 2001-node
certify grid and the 4001-node check grid it recomputes, at 70 decimal
digits, the worst relative error, the worst absolute error and the
denominator ratio min|den|/max|den|, by Horner in binary64 by both the fused
(math.fma) and the separate route.  A cell passes when the error in its
target's kind is within the target and the guard ratio is at least 1e-3, on
both grids.  The recomputed numbers must also agree with the record's own
delivered numbers to within four binary64 ulps.

Exit 0 when every cell passes, 1 otherwise.
"""
from __future__ import annotations

import glob
import json
import math
import sys
import time

import mpmath as mp

mp.mp.dps = 70


def f0(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def grid(a, W, n):
    ts = [2.0 * i / n - 1.0 for i in range(n + 1)]
    xs = [a + W * (t + 1.0) / 2.0 for t in ts]
    return ts, [f0(mp.mpf(float(x))) for x in xs]


def horner(cs, t, fused):
    acc = 0.0
    for c in cs[::-1]:
        acc = math.fma(acc, t, c) if fused else acc * t + c
    return acc


def evaluate(P, Q, ts, fs):
    pd = [float(c) for c in P]
    qd = [float(c) for c in Q[1:]]
    wr = wa = 0.0
    dmin = math.inf
    dmax = 0.0
    for fused in (True, False):
        for t, fm in zip(ts, fs):
            t = float(t)
            num = horner(pd, t, fused)
            hq = horner(qd, t, fused)
            den = math.fma(hq, t, 1.0) if fused else hq * t + 1.0
            v = num / den
            f = float(fm)
            e = abs(v - f)
            if e > wa:
                wa = e
            r = e / f
            if r > wr:
                wr = r
            d = abs(den)
            if d < dmin:
                dmin = d
            if d > dmax:
                dmax = d
    return wr, wa, dmin / dmax


def near(x, y):
    return abs(x - y) <= 4 * max(abs(x), abs(y), 1e-300) * 2.0 ** -52


def main() -> int:
    t0 = time.time()
    bad = []
    paths = sorted(glob.glob("cells/A2_*.final.json"))
    if len(paths) != 45:
        print(f"expected 45 final records, found {len(paths)}")
        return 1
    for path in paths:
        r = json.loads(open(path, encoding="utf-8").read())
        kind, target = r["kind"], float(r["target"])
        w = r["witness"]
        for n, label in ((2000, "2001"), (4000, "4001")):
            ts, fs = grid(r["a"], r["W"], n)
            rel, ab, guard = evaluate(w["P"], w["Q"], ts, fs)
            delivered = rel if kind == "rel" else ab
            recorded = w[f"rel_{label}"] if kind == "rel" else w[f"abs_{label}"]
            rec_guard = w[f"den_{label}"]
            if delivered > target:
                bad.append(f"{r['cell_id']} {label}: {delivered:.3e} > {target:.0e}")
            if guard < 1e-3:
                bad.append(f"{r['cell_id']} {label}: guard {guard:.3f} < 1e-3")
            if not (near(delivered, recorded) and near(guard, rec_guard)):
                bad.append(f"{r['cell_id']} {label}: record mismatch "
                           f"({delivered!r} vs {recorded!r}, "
                           f"{guard!r} vs {rec_guard!r})")
        print(f"OK {r['cell_id']} count={r['corrected_count']} "
              f"pair={r['pair']} kind={kind} target={target:.0e}")
    print(f"\n{len(paths)} cells, {len(bad)} failures, "
          f"{time.time() - t0:.1f}s")
    for b in bad:
        print("FAIL", b)
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
