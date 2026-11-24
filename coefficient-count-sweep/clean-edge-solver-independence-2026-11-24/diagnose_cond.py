"""Conditioning probe for the [6/8] rational reference equations.

Run the same rational-interpolation + Newton step at 70 and 90 digits and
print the Jacobian condition estimate and the Newton residual history.
"""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

A = 62.0
W = 166.72481653881133
M, K = 6, 8
N = M + K + 2


def run(dps):
    mp.mp.dps = dps
    prob = R.problem_f0(A, W)
    ts, tsm, fs = R.build_grid(A, W, 2000)
    ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]
    # rational interpolant at T_{N-1} roots
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
    errs, qvals = R.profile(P, Q, tsm, fs, ws)
    gm = max(abs(e) for e in errs)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    ref = R.choose_window(merged, N)
    print(f"dps={dps}: interpolant grid_max={mp.nstr(gm, 8)} "
          f"guard={mp.nstr(R.guard_ratio(qvals), 8)} extrema={len(merged)} "
          f"ref={len(ref) if ref else None}")

    # build the Jacobian at the reference and estimate its condition
    J = []
    Rv = []
    for i, t in enumerate(ref):
        f = prob.f(t)
        w = prob.weight(t, f)
        s = 1 if i % 2 == 0 else -1
        Pi = R.poly_eval(P, t)
        Qi = R.poly_eval(Q, t)
        Rv.append(w * (Pi / Qi - f) - s * mp.mpf(0))
        row = []
        tp = mp.mpf(1)
        for _ in range(M + 1):
            row.append(w * tp)
            tp *= t
        tp = mp.mpf(1)
        for _ in range(K):
            tp *= t
            row.append(-w * f * tp - s * mp.mpf(0) * tp)
        row.append(mp.mpf(-s) * Qi)
        J.append(row)
    try:
        U, S, V = mp.svd(mp.matrix(J))
        cond = max(S) / min(S)
        print(f"  Jacobian singular values: max={mp.nstr(max(S), 6)} "
              f"min={mp.nstr(min(S), 6)} cond={mp.nstr(cond, 6)}")
    except Exception as ex:
        print("  svd failed:", repr(ex))

    # Newton residual history
    Pn, Qn, En = P, Q, mp.mpf(0)
    for it in range(12):
        J2 = []
        R2 = []
        for i, t in enumerate(ref):
            f = prob.f(t)
            w = prob.weight(t, f)
            s = 1 if i % 2 == 0 else -1
            Pi = R.poly_eval(Pn, t)
            Qi = R.poly_eval(Qn, t)
            R2.append(w * (Pi / Qi - f) - s * En)
            row = []
            tp = mp.mpf(1)
            for _ in range(M + 1):
                row.append(w * tp)
                tp *= t
            tp = mp.mpf(1)
            for _ in range(K):
                tp *= t
                row.append(-w * f * tp - s * En * tp)
            row.append(mp.mpf(-s) * Qi)
            J2.append(row)
        res = max(abs(r) for r in R2)
        try:
            d = mp.lu_solve(mp.matrix(J2), mp.matrix([-r for r in R2]))
        except Exception as ex:
            print("  newton solve failed", repr(ex))
            break
        step = max(abs(x) for x in d)
        print(f"  it={it:2d} res={mp.nstr(res, 5)} step={mp.nstr(step, 5)}")
        if res < mp.mpf(10) ** (-(dps - 8)):
            break
        for j in range(M + 1):
            Pn[j] += d[j]
        for j in range(K):
            Qn[1 + j] += d[M + 1 + j]
        En += d[M + 1 + K]


for dps in (70, 90, 110):
    run(dps)
    print()
