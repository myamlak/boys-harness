#!/usr/bin/env python3
"""Step (a): alternation/de la Vallee Poussin bound on the stored [5/6] fit.

The stored 12-coefficient nominal truncation is, by hypothesis, a near-minimax
type-(5,6) fit.  This script:
  * locates the 13 alternating extrema of F_0 - r numerically;
  * proves, with exact rational arithmetic at those 13 nodes, a lower bound
    delta on the best type-(5,6) uniform error (13 = m+k+2, the dVP count);
  * reports the candidate's exact sup at those nodes.

The gap between delta and the candidate's sup measures how far the stored fit
is from the exact type-(5,6) minimax.  This step bounds the fit's own type; it
does not by itself bound any stored count.

Usage (from the deposit root):
  python -I alternation_stored_candidate.py \
      --corrected inputs/paper_cell_abs_corrected.json \
      --out alternation_stored_candidate.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f0exact as fx  # noqa: E402

A = Fraction("11.899848152108484")
W = Fraction("17.089489586712256")
CORRECTED = (r"C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts"
             r"\results-corrected\corrected-counts.json")
TARGET = fx.TARGET
X0_MP = None
W_MP = None


def f0_mp(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def f1_mp(x):
    return (f0_mp(x) - mp.e ** (-x)) / (2 * x)


def horner(cs, t):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * t + c
    return acc


def deriv(cs):
    return [i * cs[i] for i in range(1, len(cs))]


def g_value(P, Q, t):
    x = X0_MP + W_MP * (t + 1) / 2
    return f0_mp(x) - horner(P, t) / horner(Q, t)


def g_prime(P, Q, dP, dQ, t):
    x = X0_MP + W_MP * (t + 1) / 2
    p, q = horner(P, t), horner(Q, t)
    return -(W_MP / 2) * f1_mp(x) - (horner(dP, t) * q -
                                    p * horner(dQ, t)) / (q * q)


def find_extrema(P, Q, n=8001):
    dP, dQ = deriv(P), deriv(Q)
    ts = [mp.mpf(2 * i) / (n - 1) - 1 for i in range(n)]
    gs = [g_prime(P, Q, dP, dQ, t) for t in ts]
    pts = [(ts[0], g_value(P, Q, ts[0])),
           (ts[-1], g_value(P, Q, ts[-1]))]
    for i in range(n - 1):
        if gs[i] == 0:
            pts.append((ts[i], g_value(P, Q, ts[i])))
        elif (gs[i] > 0) != (gs[i + 1] > 0):
            lo, hi = ts[i], ts[i + 1]
            for _ in range(90):
                mid = (lo + hi) / 2
                if (g_prime(P, Q, dP, dQ, mid) > 0) == (gs[i] > 0):
                    lo = mid
                else:
                    hi = mid
            t = (lo + hi) / 2
            pts.append((t, g_value(P, Q, t)))
    pts.sort(key=lambda z: z[0])
    dedup = []
    for t, g in pts:
        if dedup and abs(t - dedup[-1][0]) < mp.mpf("1e-12"):
            continue
        dedup.append((t, g))
    alt = []
    for t, g in dedup:
        if not alt or (g > 0) != (alt[-1][1] > 0):
            alt.append((t, g))
        elif abs(g) > abs(alt[-1][1]):
            alt[-1] = (t, g)
    return dedup, alt


def main(argv=None):
    global X0_MP, W_MP
    started = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--corrected", default=CORRECTED)
    ap.add_argument("--out", default="alternation_stored_candidate.json")
    args = ap.parse_args(argv)
    mp.mp.dps = 60
    X0_MP = mp.mpf(A.numerator) / A.denominator
    W_MP = mp.mpf(W.numerator) / W.denominator

    doc = json.load(open(args.corrected, encoding="utf-8"))
    cell = next(c for c in doc["cells"]
                if c.get("width_name") == "paper"
                and c.get("target") == "abs:5e-14")
    m, k = 5, 6
    P_str = cell["witness_P"][:m + 1]
    Q_str = cell["witness_Q"][:k + 1]
    P = [Fraction(float(s)) for s in P_str]
    Q = [Fraction(float(s)) for s in Q_str]
    P_mp = [mp.mpf(s) for s in P_str]
    Q_mp = [mp.mpf(s) for s in Q_str]

    ext, alt = find_extrema(P_mp, Q_mp)
    print(f"located {len(ext)} extrema, {len(alt)} alternating "
          f"(need {m+k+2})")
    if len(alt) < m + k + 2:
        print("FAIL: too few alternating extrema")
        return 1
    nodes = []
    for t, _g in alt:
        x = X0_MP + W_MP * (t + 1) / 2
        xf = Fraction(mp.nstr(x, 20))
        tf = 2 * (xf - A) / W - 1
        nodes.append((xf, tf))
    ok, delta, rows = fx.dvp_check_mapped(P, Q, nodes,
                                          target=Fraction(1, 10**30),
                                          n_terms=120)
    if not ok:
        print("FAIL: exact dVP check did not produce alternating signs")
        return 1
    signs = [r["sign"] for r in rows]
    alternating = all(signs[i] != signs[i + 1] for i in range(len(signs) - 1))
    upper = Fraction(0)
    for x, t in nodes:
        lo, hi = fx.f0_bounds(x, n_terms=120)
        r, _ = fx.rat_eval(P, Q, t)
        upper = max(upper, abs(lo - r), abs(hi - r))
    print(f"type-(5,6) exact dVP lower bound: {float(delta):.6e}")
    print(f"alternating: {alternating}, nodes: {len(nodes)}")
    print(f"stored candidate's exact sup at the 13 nodes: "
          f"{float(upper):.6e}")
    print(f"ratio candidate sup / lower bound: "
          f"{float(upper/delta):.6f}")
    out = dict(
        schema="paper-c-alternation-stored-candidate-1",
        pair=[m, k], needed_alternations=m + k + 2,
        located_extrema=len(ext), alternating_nodes=len(nodes),
        signs=signs,
        dvp_lower_bound="%.6e" % float(delta),
        dvp_lower_num=str(delta.numerator),
        dvp_lower_den=str(delta.denominator),
        candidate_sup_at_nodes="%.6e" % float(upper),
        sup_over_lower="%.6f" % float(upper / delta),
        target="5e-14",
        note=("The 13-point dVP bound is exact for type (5,6). The gap "
              "between it and the stored candidate's sup shows the stored "
              "BRASIL fit is near-minimax but not exactly equioscillating; "
              "the alternation theorem's equality form does not certify it "
              "as the exact type-(5,6) minimax.")
    )
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {args.out} ({time.time()-started:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
