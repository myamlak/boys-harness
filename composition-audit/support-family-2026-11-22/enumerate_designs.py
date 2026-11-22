"""Enumerate every 0/1 design with Omega_{T,R0} = 0 exactly.

Homogeneous system A n = 0 (98 alias equations), n in {0,1}^36.
Iteratively solve MILP and exclude found supports.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from instance import N_CELLS, PHI, T_IDX, R0_IDX, CELLS, AXES, SUPPORTS, CHOICES

PAIRS = [(a, b) for a in T_IDX for b in R0_IDX]
A = np.array([PHI[:, a] * PHI[:, b] for a, b in PAIRS])
scale = np.abs(A).max(axis=1)
scale[scale == 0] = 1.0
A = A / scale[:, None]
# add the homogeneous count row? no: keep A n = 0 only
CONS = [LinearConstraint(A, np.full(len(A), -1e-9), np.full(len(A), 1e-9))]


def solve(N=None, exclude=(), objective=None, seed=0):
    cons = list(CONS)
    if N is not None:
        ones = np.ones((1, N_CELLS))
        cons.append(LinearConstraint(ones, [N], [N]))
    # exclude every previously found set
    Aex = []
    bex = []
    for S in exclude:
        row = np.zeros(N_CELLS)
        row[list(S)] = 1.0
        # exclude exactly S, allow supersets: sum_S n - sum_notS n <= |S|-1
        row[list(set(range(N_CELLS)) - set(S))] = -1.0
        Aex.append(row)
        bex.append(len(S) - 1.0)
    if Aex:
        cons.append(LinearConstraint(np.array(Aex), np.full(len(Aex), -np.inf), np.array(bex)))
    if objective is None:
        c = np.zeros(N_CELLS)
    else:
        c = objective
    res = milp(
        c,
        constraints=cons,
        integrality=np.ones(N_CELLS),
        bounds=Bounds(np.zeros(N_CELLS), np.ones(N_CELLS)),
        options={"time_limit": 60, "disp": False},
    )
    if res.status != 0 or res.x is None:
        return None, res
    n = np.round(res.x).astype(int)
    return frozenset(np.flatnonzero(n)), res


def rank_of(S):
    idx = sorted(S)
    M = PHI[np.ix_(idx, T_IDX)]
    # Gram with weights 1/N
    G = M.T @ M / len(idx)
    return np.linalg.matrix_rank(G, tol=1e-9), np.linalg.eigvalsh(G)


if __name__ == "__main__":
    found = []
    for N in range(1, 37):
        print(f"--- N={N}")
        while True:
            S, res = solve(N=N, exclude=found)
            if S is None:
                print("   no (more) solutions")
                break
            found.append(S)
            r, ev = rank_of(S)
            print(f"   design size {len(S)}, rank {r}, evmin {ev.min():.3e}, "
                  f"cells {sorted(S)}")
            if len(found) > 300:
                raise SystemExit("too many")
    print(f"\ntotal {len(found)} designs")
    for S in sorted(found, key=len):
        r, ev = rank_of(S)
        print(f"N={len(S):2d} rank={r:2d} evmin={ev.min():.3e} cells={sorted(S)}")
