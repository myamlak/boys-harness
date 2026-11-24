"""Instrument the [6/1] continuation step: Newton residual/step history and
Jacobian condition, with and without line search."""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

mp.mp.dps = R.DPS
A = 62.0
W = 166.72481653881133
prob = R.problem_f0(A, W)
ts, tsm, fs = R.build_grid(A, W, 2000)
ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]

# [6/0] first
r0 = R.solve_pair(6, 0, tsm, fs, prob, mp.mpf("1e-15"))
P0 = [mp.mpf(c) for c in r0["P"]]
E0 = mp.mpf(r0["E"])
print("polynomial [6/0]:", r0["grid_max"], r0["converged"])
print("P0 =", [mp.nstr(c, 8) for c in P0], "E0 =", mp.nstr(E0, 6))

K = 1
M = 6
N = M + K + 2
ref = R.chebyshev_reference(N)
P = P0
Q = [mp.mpf(1), mp.mpf(0)]
E = E0


def jac_res(P, Q, E, ref):
    J = []
    Rv = []
    for i, t in enumerate(ref):
        f = prob.f(t)
        w = prob.weight(t, f)
        s = 1 if i % 2 == 0 else -1
        Pi = R.poly_eval(P, t)
        Qi = R.poly_eval(Q, t)
        Rv.append(w * (Pi / Qi - f) - s * E)
        row = []
        tp = mp.mpf(1)
        for _ in range(M + 1):
            row.append(w * tp)
            tp *= t
        tp = mp.mpf(1)
        for _ in range(K):
            tp *= t
            row.append(-w * f * tp - s * E * tp)
        row.append(mp.mpf(-s) * Qi)
        J.append(row)
    return J, Rv


def grid_stats(P, Q):
    errs, qv = R.profile(P, Q, tsm, fs, ws)
    return max(abs(e) for e in errs), R.guard_ratio(qv)


print(f"\nreference N={N}:", [mp.nstr(t, 8) for t in ref])
errs, qv = R.profile(P, Q, tsm, fs, ws)
print(f"init grid_max={mp.nstr(max(abs(e) for e in errs), 8)} "
      f"guard={mp.nstr(R.guard_ratio(qv), 6)}")
print("init ref errors:", [mp.nstr(e, 5) for e in
                           [R.error_at(P, Q, t, prob)[0] for t in ref]])

J, Rv = jac_res(P, Q, E, ref)
U, S, V = mp.svd(mp.matrix(J))
print("cond=", mp.nstr(max(S) / min(S), 6),
      "sv range", mp.nstr(min(S), 4), mp.nstr(max(S), 4))

for it in range(15):
    J, Rv = jac_res(P, Q, E, ref)
    res = max(abs(r) for r in Rv)
    d = mp.lu_solve(mp.matrix(J), mp.matrix([-r for r in Rv]))
    step = max(abs(x) for x in d)
    print(f"it={it:2d} res={mp.nstr(res, 6)} step={mp.nstr(step, 6)} "
          f"dE={mp.nstr(d[M + 1 + K], 6)}")
    # plain Newton
    for j in range(M + 1):
        P[j] += d[j]
    for j in range(K):
        Q[1 + j] += d[M + 1 + j]
    E += d[M + 1 + K]
    gm, g = grid_stats(P, Q)
    print(f"      grid_max={mp.nstr(gm, 6)} guard={mp.nstr(g, 6)} "
          f"Q1={mp.nstr(Q[1], 6)}")
