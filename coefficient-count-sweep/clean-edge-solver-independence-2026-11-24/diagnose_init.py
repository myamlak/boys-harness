"""Multi-start diagnostic for p5 [6/8]: which initialization reaches a
normal (16-extremum) minimax?

Initializations:
  * the BRASIL accepted fit (diagnostic reference, not used for the answer);
  * Q(t) = 1 + c t for a sweep of c, with P the weighted least-squares fit of
    f*Q at the grid nodes (a generic positive-denominator start).

For each start: Newton-Remez loop (refine extrema -> exchange -> exact Newton),
report the best grid max, the guard, the number of alternating extrema and the
equioscillation ratio.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

mp.mp.dps = R.DPS
HERE = Path(__file__).resolve().parent
CELL = Path(r"C:\Users\myaml\dev\_c60-clean-edge-full-2026-11-24\cells"
            r"\p5_rel_1e-15.json")

doc = json.loads(CELL.read_text(encoding="utf-8"))
acc = doc["accepted"]
a = float(doc["cell"]["a"])
W = float(doc["cell"]["width"])
prob = R.problem_f0(a, W)
ts, tsm, fs = R.build_grid(a, W, 2000)
ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]
n = 2000
M, K = 6, 8
N = M + K + 2


def weighted_ls_P(Q):
    """Weighted least squares for P with Q fixed: min sum w^2 (P - f Q)^2."""
    rows = []
    rhs = []
    for t, f, w in zip(tsm, fs, ws):
        q = R.poly_eval(Q, t)
        row = []
        tp = mp.mpf(1)
        for _ in range(M + 1):
            row.append(w * tp)
            tp *= t
        rows.append(row)
        rhs.append(w * f * q)
    A = mp.matrix(rows)
    b = mp.matrix(rhs)
    AtA = A.T * A
    Atb = A.T * b
    P = mp.lu_solve(AtA, Atb)
    return [P[j] for j in range(M + 1)]


def stats(P, Q):
    errs, qvals = R.profile(P, Q, tsm, fs, ws)
    gm = max(abs(e) for e in errs)
    g = R.guard_ratio(qvals)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    return gm, g, merged


def newton_remez(P, Q, E, iters=25, tag=""):
    best = None
    for it in range(iters):
        gm, g, merged = stats(P, Q)
        if g >= R.DEN_FLOOR and (best is None or gm < best[0]):
            best = (gm, g, len(merged), P, Q, E, it)
        ref = R.choose_window(merged, N)
        if ref is None:
            ref = R.chebyshev_reference(N)
        Pn, Qn, En, ok, note, res = R.newton_solve(M, K, ref, prob, P, Q, E)
        if not ok:
            break
        P, Q, E = Pn, Qn, En
    gm, g, merged = stats(P, Q)
    if g >= R.DEN_FLOOR and (best is None or gm < best[0]):
        best = (gm, g, len(merged), P, Q, E, iters)
    gm, g, ne, Pb, Qb, Eb, itb = best
    print(f"  {tag:>26s}: best grid_max={mp.nstr(gm, 10)} guard="
          f"{mp.nstr(g, 6)} extrema={ne} at it={itb}")
    return best


print("BRASIL [6/8] as a start:")
P = [mp.mpf(c) for c in acc["P70"]]
Q = [mp.mpf(c) for c in acc["Q70"]]
gm, g, merged = stats(P, Q)
print(f"  stored: grid_max={mp.nstr(gm, 10)} guard={mp.nstr(g, 6)} "
      f"extrema={len(merged)}")
newton_remez(P, Q, mp.mpf("0.1"), tag="newton from BRASIL")

print("\npositive-denominator least-squares starts, Q = 1 + c t:")
for c in ("0", "0.2", "0.4", "0.5", "0.55", "0.6", "0.7", "0.8", "0.9"):
    Q0 = [mp.mpf(1), mp.mpf(c)]
    P0 = weighted_ls_P(Q0)
    errs, _ = R.profile(P0, Q0, tsm, fs, ws)
    E0 = max(abs(e) for e in errs)
    newton_remez(P0, Q0, E0, tag=f"Q=1+{c}t")
