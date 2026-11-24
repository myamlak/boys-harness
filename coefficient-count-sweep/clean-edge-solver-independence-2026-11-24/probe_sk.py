"""Sanathanan-Koerner (weighted least-squares) rational fit for p5 [6/8].

A different algorithm family from minimax: linearised weighted least squares
iterated on the denominator.  It is used only as an initializer; the minimax
certificate still comes from the Remez/Newton solver.
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


def sk_fit(Q0, iters=20, tag=""):
    Q = list(Q0) + [mp.mpf(0)] * (K + 1 - len(Q0))
    for it in range(iters):
        nb = M + 1 + K
        Mtx = [[mp.mpf(0)] * nb for _ in range(nb)]
        rhs = [mp.mpf(0)] * nb
        for t, f, w in zip(tsm, fs, ws):
            q = R.poly_eval(Q, t)
            rho = 1 / (f * q) ** 2
            phi = []
            tp = mp.mpf(1)
            for _ in range(M + 1):
                phi.append(tp)
                tp *= t
            tp = mp.mpf(1)
            for _ in range(K):
                tp *= t
                phi.append(-f * tp)
            for a in range(nb):
                rhs[a] += rho * phi[a] * f
                for b in range(a, nb):
                    Mtx[a][b] += rho * phi[a] * phi[b]
        for a in range(nb):
            for b in range(a):
                Mtx[a][b] = Mtx[b][a]
        u = mp.lu_solve(mp.matrix(Mtx), mp.matrix(rhs))
        P = [u[j] for j in range(M + 1)]
        Q = [mp.mpf(1)] + [u[M + 1 + j] for j in range(K)]
        errs, qv = R.profile(P, Q, tsm, fs, ws)
        gm = max(abs(e) for e in errs)
        g = R.guard_ratio(qv)
        if it % 4 == 0 or it == iters - 1:
            print(f"  {tag} it={it:2d} grid_max={mp.nstr(gm, 8)} "
                  f"guard={mp.nstr(g, 6)} Q1={mp.nstr(Q[1], 6)} "
                  f"maxQ={mp.nstr(max(abs(q) for q in Q), 6)}", flush=True)
    return P, Q


def newton_try(P, Q, tag, iters=40):
    E = mp.mpf(0)
    best = None
    for it in range(iters):
        errs, qv = R.profile(P, Q, tsm, fs, ws)
        gm = max(abs(e) for e in errs)
        g = R.guard_ratio(qv)
        merged = R.refine_extrema(P, Q, ts, errs, prob)
        if g >= R.DEN_FLOOR and (best is None or gm < best[0]):
            best = (gm, g, len(merged), list(P), list(Q))
        ref = R.choose_window(merged, N)
        if ref is None:
            ref = R.chebyshev_reference(N)
        Pn, Qn, En, ok, note, res = R.newton_solve(M, K, ref, prob, P, Q, E)
        # even a "failed" Newton with a tiny residual is useful: accept a fit
        # whose reference residual is far below the fit error
        if not ok and res is not None and res < gm * mp.mpf("1e-3"):
            ok = True
            note = f"accepted residual {mp.nstr(res, 5)}"
        if not ok:
            print(f"  {tag}: newton stop it={it}: {note}")
            break
        P, Q, E = Pn, Qn, En
    errs, qv = R.profile(P, Q, tsm, fs, ws)
    gm = max(abs(e) for e in errs)
    g = R.guard_ratio(qv)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    if g >= R.DEN_FLOOR and (best is None or gm < best[0]):
        best = (gm, g, len(merged), list(P), list(Q))
    gm, g, ne, Pb, Qb = best
    print(f"  {tag}: best grid_max={mp.nstr(gm, 10)} guard={mp.nstr(g, 6)} "
          f"extrema={ne}")
    return Pb, Qb


for c0 in ("0", "0.5", "0.55", "0.6"):
    print(f"\nSK start Q=1+{c0}t")
    P, Q = sk_fit([mp.mpf(1), mp.mpf(c0)], tag=f"c={c0}")
    Pb, Qb = newton_try(P, Q, tag=f"newton c={c0}")
    s = R.score60(Pb, Qb, A, W, 1e-15, grid=(ts, tsm, fs))
    print(f"  score60 e60={s['e60']} guard60={s['guard60']} "
          f"meets={s['meets60']}")
