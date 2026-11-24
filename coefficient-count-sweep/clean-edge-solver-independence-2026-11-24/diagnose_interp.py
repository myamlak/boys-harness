"""Rational-interpolation initialization diagnostic for p5 [6/8].

Build the Chebyshev rational interpolant of type [6/8] (P - f Q = 0 at 15
Chebyshev nodes), then run the Newton-Remez exchange from it.  This is an
independent, generic initialization that should land on the smooth
positive-denominator branch.
"""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

mp.mp.dps = R.DPS
A = 62.0
W = 166.72481653881133
M, K = 6, 8
N = M + K + 2
prob = R.problem_f0(A, W)
ts, tsm, fs = R.build_grid(A, W, 2000)
ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]


def stats(P, Q):
    errs, qvals = R.profile(P, Q, tsm, fs, ws)
    gm = max(abs(e) for e in errs)
    g = R.guard_ratio(qvals)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    return gm, g, merged


def interp_fit(nodes):
    A_ = []
    b_ = []
    for t in nodes:
        f = prob.f(t)
        row = []
        tp = mp.mpf(1)
        for _ in range(M + 1):
            row.append(tp)
            tp *= t
        tp = mp.mpf(1)
        for _ in range(K):
            tp *= t
            row.append(-f * tp)
        A_.append(row)
        b_.append(f)
    u = mp.lu_solve(mp.matrix(A_), mp.matrix(b_))
    P = [u[j] for j in range(M + 1)]
    Q = [mp.mpf(1)] + [u[M + 1 + j] for j in range(K)]
    return P, Q


def newton_remez(P, Q, E, iters=30, tag=""):
    best = None
    hist = []
    for it in range(iters):
        gm, g, merged = stats(P, Q)
        hist.append((gm, g, len(merged)))
        if g >= R.DEN_FLOOR and (best is None or gm < best[0]):
            best = (gm, g, len(merged), P, Q, E, it)
        ref = R.choose_window(merged, N)
        if ref is None:
            ref = R.chebyshev_reference(N)
        Pn, Qn, En, ok, note, res = R.newton_solve(M, K, ref, prob, P, Q, E)
        if not ok:
            print(f"    {tag}: newton stopped at it={it}: {note} res={res}")
            break
        P, Q, E = Pn, Qn, En
    gm, g, merged = stats(P, Q)
    hist.append((gm, g, len(merged)))
    if g >= R.DEN_FLOOR and (best is None or gm < best[0]):
        best = (gm, g, len(merged), P, Q, E, iters)
    gm, g, ne, Pb, Qb, Eb, itb = best
    print(f"  {tag:>28s}: best grid_max={mp.nstr(gm, 12)} guard="
          f"{mp.nstr(g, 6)} extrema={ne} at it={itb}")
    print("     history: " + ", ".join(mp.nstr(h[0], 4) for h in hist[::5]))
    return best


print("rational interpolation initializations for [6/8]:")
for label, nodes in (
    ("T_15 roots", [mp.cos(mp.pi * (2 * i + 1) / (2 * (N - 1)))
                    for i in range(N - 1)]),
    ("T_15 extrema-minus-one",
     [mp.cos(mp.pi * i / (N - 1)) for i in range(N - 1)]),
    ("equispaced 15", [mp.mpf(-1) + 2 * mp.mpf(i) / (N - 2)
                       for i in range(N - 1)]),
):
    P, Q = interp_fit(nodes)
    gm, g, merged = stats(P, Q)
    print(f"\n{label}: interpolant grid_max={mp.nstr(gm, 10)} "
          f"guard={mp.nstr(g, 6)} extrema={len(merged)}")
    newton_remez(P, Q, mp.mpf(0), tag=label)
