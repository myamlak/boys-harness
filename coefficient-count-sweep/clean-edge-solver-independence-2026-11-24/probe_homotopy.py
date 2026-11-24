"""Function homotopy for p5 [6/8].

Start from f_0(t) = A/(B+t), a rational with a positive linear denominator
matching the good-basin structure (B = 1.8; the BRASIL fit has Q ~ 1+0.556 t).
The [6/8] family contains f_0 exactly.  Continue f_lambda = (1-lambda) f_0 +
lambda f in lambda, warm-starting each step from the previous minimax fit, and
run the Remez/Newton solver.  If the path stays on the normal branch the
lambda=1 fit is the minimax of f.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

mp.mp.dps = R.DPS
A_INTERVAL = 62.0
W = 166.72481653881133
M, K = 6, 8
N = M + K + 2
prob_f = R.problem_f0(A_INTERVAL, W)
ts, tsm, fs = R.build_grid(A_INTERVAL, W, 2000)
f_mid = prob_f.f(mp.mpf(0))
B = mp.mpf("1.8")
A0 = f_mid * B
print(f"f_0(t) = {mp.nstr(A0, 10)} / ({B} + t)")


def f0_eval(t):
    return A0 / (B + t)


def interp_init(prob, fs_lam):
    """Rational interpolant of the current f at the T_15 roots."""
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
    return P, Q, mp.mpf(0)


def make_problem(lam):
    def eval_f(t):
        return (1 - lam) * f0_eval(t) + lam * prob_f.f(t)

    return R.Problem(eval_f, lambda t, f: 1 / f)


def stats(P, Q, prob, fs_lam):
    ws = [prob.weight(t, f) for t, f in zip(tsm, fs_lam)]
    errs, qv = R.profile(P, Q, tsm, fs_lam, ws)
    gm = max(abs(e) for e in errs)
    g = R.guard_ratio(qv)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    return gm, g, len(merged), ws


lams = [mp.mpf(s) for s in
        ("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9",
         "0.95", "1.0")]
P = Q = E = None
ref = None
t_start = time.perf_counter()
best = None
for lam in lams:
    prob = make_problem(lam)
    fs_lam = [prob.f(t) for t in tsm]
    if P is None:
        P0, Q0, E0 = interp_init(prob, fs_lam)
        gm0, g0, ne0, _ = stats(P0, Q0, prob, fs_lam)
        print(f"lambda={mp.nstr(lam, 4)} interpolant grid_max="
              f"{mp.nstr(gm0, 6)} guard={mp.nstr(g0, 5)} extrema={ne0}",
              flush=True)
        init = (P0, Q0, E0)
        ref0 = None
    else:
        init = (P, Q, E)
        ref0 = ref
    t0 = time.perf_counter()
    r = R.solve_pair(M, K, tsm, fs_lam, prob, mp.mpf("1e-15"),
                     dc_iters=0, max_outer=80, exchange_mode="global",
                     ref0=ref0, init=init)
    print(f"lambda={mp.nstr(lam, 4)} conv={r['converged']} "
          f"grid_max={r['grid_max']} guard={r['guard']} ratio={r['ratio']} "
          f"it={r['iterations']} newton={r['newton_note']} "
          f"({time.perf_counter()-t0:.0f}s)", flush=True)
    P = [mp.mpf(c) for c in r["P"]]
    Q = [mp.mpf(c) for c in r["Q"]]
    E = mp.mpf(r["E"])
    ref = [mp.mpf(t) for t in r["ref"]]
    if r["converged"] or best is None or mp.mpf(r["grid_max"]) < best[0]:
        best = (mp.mpf(r["grid_max"]), r["converged"], lam, P, Q)

gm, conv, lam, Pb, Qb = best
print(f"\nbest: lambda={mp.nstr(lam, 4)} grid_max={mp.nstr(gm, 17)} "
      f"converged={conv}")
s = R.score60(Pb, Qb, A_INTERVAL, W, 1e-15, grid=(ts, tsm, fs))
print("score60:", {k: s[k] for k in ("e60", "guard60", "meets60",
                                     "e60_at_t")})
print("Q=", [mp.nstr(c, 20) for c in Qb])
print(f"total {time.perf_counter()-t_start:.0f}s")
