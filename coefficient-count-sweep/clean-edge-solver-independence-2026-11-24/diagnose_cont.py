"""Degree continuation diagnostic for p5 [6/k], k = 0..8.

Start from the exact polynomial minimax [6/0] (known-answer-tested solver),
then increase the denominator degree one at a time, warm-starting each step
from the previous fit with the exact Newton solve on a Chebyshev reference.
The continuation stays on the positive-denominator branch instead of jumping
to the near-pole spurious solution found by a cold differential correction.
"""

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
M, KMAX = 6, 8
prob = R.problem_f0(A, W)
ts, tsm, fs = R.build_grid(A, W, 2000)
ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]

started = time.perf_counter()
# ---- [6/0]: cold solve, the polynomial branch
res = R.solve_pair(M, 0, tsm, fs, prob, mp.mpf("1e-15"), verbose=False)
print(f"[6/0] conv={res['converged']} grid_max={res['grid_max']} "
      f"guard={res['guard']} extrema_ref={len(res['ref'])} "
      f"({time.perf_counter()-started:.1f}s)", flush=True)
P = [mp.mpf(c) for c in res["P"]]
Q = [mp.mpf(c) for c in res["Q"]]
E = mp.mpf(res["E"])
best = (mp.mpf(res["grid_max"]), res["converged"], 0, res)

for k in range(1, KMAX + 1):
    t0 = time.perf_counter()
    Qpad = list(Q) + [mp.mpf(0)] * (k + 1 - len(Q))
    ref = R.chebyshev_reference(M + k + 2)
    r = R.solve_pair(M, k, tsm, fs, prob, mp.mpf("1e-15"),
                     dc_iters=0, max_outer=60, verbose=False,
                     ref0=ref, init=(P, Qpad, E))
    print(f"[{M}/{k}] conv={r['converged']} status={r['status']} "
          f"grid_max={r['grid_max']} guard={r['guard']} ratio={r['ratio']} "
          f"newton={r['newton_note']} ({time.perf_counter()-t0:.1f}s)",
          flush=True)
    P = [mp.mpf(c) for c in r["P"]]
    Q = [mp.mpf(c) for c in r["Q"]]
    E = mp.mpf(r["E"])
    if mp.mpf(r["grid_max"]) < best[0]:
        best = (mp.mpf(r["grid_max"]), r["converged"], k, r)

gm, conv, k, r = best
print(f"\nbest along continuation: [6/{k}] grid_max={mp.nstr(gm, 17)} "
      f"converged={conv}")
Pb = [mp.mpf(c) for c in r["P"]]
Qb = [mp.mpf(c) for c in r["Q"]]
s = R.score60(Pb, Qb, A, W, 1e-15, grid=(ts, tsm, fs))
print("60-digit score of the best fit: e60=", s["e60"],
      "guard60=", s["guard60"], "meets=", s["meets60"],
      "at t=", s["e60_at_t"])
print(f"total {time.perf_counter()-started:.1f}s")
