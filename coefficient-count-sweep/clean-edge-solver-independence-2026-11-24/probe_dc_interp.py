"""Probe: differential-correction exchange started from the [6/8] rational
interpolant's own extrema (16 points).  Pure DC steps, no Newton.
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

# rational interpolant at the T_15 roots
nodes = [mp.cos(mp.pi * (2 * i + 1) / (2 * (N - 1))) for i in range(N - 1)]
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


def stats(P, Q):
    errs, qv = R.profile(P, Q, tsm, fs, ws)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    return max(abs(e) for e in errs), R.guard_ratio(qv), merged, errs


gm, g, merged, errs = stats(P, Q)
print(f"interpolant: grid_max={mp.nstr(gm, 8)} guard={mp.nstr(g, 6)} "
      f"extrema={len(merged)}")
ref = R.choose_window(merged, N)
print(f"reference from its extrema: {len(ref)}")
for mode in ("global", "maehly"):
    rr = list(ref)
    pp, qq = list(P), list(Q)
    print(f"\n--- DC exchange mode={mode}")
    for it in range(40):
        pp, qq, ee = R.dc_solve(M, K, rr, prob)
        gm, g, merged, errs = stats(pp, qq)
        if mode == "maehly":
            new = R.maehly_exchange(merged, rr)
        else:
            new = R.choose_window(merged, N)
        if it % 2 == 0 or it < 10:
            print(f"  it={it:2d} grid_max={mp.nstr(gm, 8)} "
                  f"guard={mp.nstr(g, 5)} extrema={len(merged)}", flush=True)
        if new is None:
            print("  no window; stop")
            break
        if new == rr:
            print("  reference fixed; stop")
            break
        rr = new
