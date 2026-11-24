"""Strategy trials for p5 [6/8]: Maehly vs global exchange, cold vs
rational-interpolation warm start."""

from __future__ import annotations

import sys
import time
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


def interp_fit():
    nodes = [mp.cos(mp.pi * (2 * i + 1) / (2 * (N - 1)))
             for i in range(N - 1)]
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


def report(tag, r, t0):
    print(f"{tag}: conv={r['converged']} status={r['status']} "
          f"grid_max={r.get('grid_max')} guard={r.get('guard')} "
          f"ratio={r.get('ratio')} newton={r.get('newton_note')} "
          f"({time.perf_counter()-t0:.0f}s)", flush=True)


# A: cold DC, Maehly exchange
t0 = time.perf_counter()
r = R.solve_pair(M, K, tsm, fs, prob, mp.mpf("1e-15"),
                 dc_iters=150, max_outer=150, exchange_mode="maehly")
report("A cold-DC maehly", r, t0)

# B: interpolation init, Newton, Maehly
P0, Q0 = interp_fit()
t0 = time.perf_counter()
r = R.solve_pair(M, K, tsm, fs, prob, mp.mpf("1e-15"),
                 dc_iters=0, max_outer=80, exchange_mode="maehly",
                 init=(P0, Q0, mp.mpf(0)))
report("B interp Newton maehly", r, t0)

# C: interpolation init, Newton, global exchange
t0 = time.perf_counter()
r = R.solve_pair(M, K, tsm, fs, prob, mp.mpf("1e-15"),
                 dc_iters=0, max_outer=80, exchange_mode="global",
                 init=(P0, Q0, mp.mpf(0)))
report("C interp Newton global", r, t0)

# D: cold DC, global exchange (baseline)
t0 = time.perf_counter()
r = R.solve_pair(M, K, tsm, fs, prob, mp.mpf("1e-15"),
                 dc_iters=30, max_outer=60, exchange_mode="global")
report("D cold-DC global", r, t0)
