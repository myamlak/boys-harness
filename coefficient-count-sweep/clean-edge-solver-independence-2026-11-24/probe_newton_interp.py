"""Probe: backtracking Newton from the [6/8] rational interpolant on its own
16 extrema, printing the residual and the scored fit at every step."""

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
E = mp.mpf(0)


def stats(P, Q):
    errs, qv = R.profile(P, Q, tsm, fs, ws)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    return max(abs(e) for e in errs), R.guard_ratio(qv), merged


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


gm, g, merged = stats(P, Q)
ref = R.choose_window(merged, N)
print(f"start: grid_max={mp.nstr(gm, 8)} guard={mp.nstr(g, 6)} "
      f"extrema={len(merged)} ref={len(ref)}")

for it in range(30):
    J, Rv = jac_res(P, Q, E, ref)
    res = max(abs(r) for r in Rv)
    try:
        d = mp.lu_solve(mp.matrix(J), mp.matrix([-r for r in Rv]))
    except Exception as ex:
        print("solve failed", repr(ex))
        break
    step = max(abs(x) for x in d)
    # backtracking on the exact reference residual
    alpha = mp.mpf(1)
    accepted = False
    for _ in range(14):
        Pt = [P[j] + alpha * d[j] for j in range(M + 1)]
        Qt = [mp.mpf(1)] + [Q[1 + j] + alpha * d[M + 1 + j]
                            for j in range(K)]
        Et = E + alpha * d[M + 1 + K]
        nr = R._newton_residual(M, K, ref, prob, Pt, Qt, Et)
        if nr < res:
            accepted = True
            break
        alpha /= 2
    gm, g, merged = stats(P, Q)
    print(f"it={it:2d} res={mp.nstr(res, 6)} step={mp.nstr(step, 6)} "
          f"alpha={mp.nstr(alpha, 4)} accepted={accepted} "
          f"grid_max={mp.nstr(gm, 6)} guard={mp.nstr(g, 6)}")
    if not accepted:
        print("line search failed")
        break
    P, Q, E = Pt, Qt, Et
    gm, g, merged = stats(P, Q)
    print(f"      -> grid_max={mp.nstr(gm, 6)} guard={mp.nstr(g, 6)} "
          f"extrema={len(merged)}")
    if res < mp.mpf("1e-30"):
        break

gm, g, merged = stats(P, Q)
print(f"\nfinal: grid_max={mp.nstr(gm, 17)} guard={mp.nstr(g, 8)} "
      f"extrema={len(merged)}")
s = R.score60(P, Q, A, W, 1e-15, grid=(ts, tsm, fs))
print("score60:", {k: s[k] for k in ("e60", "guard60", "meets60",
                                     "e60_at_t")})
print("P=", [mp.nstr(c, 25) for c in P])
print("Q=", [mp.nstr(c, 25) for c in Q])
