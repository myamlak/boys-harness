"""Exact rank and certified bound for the 75-run pooled design.

D75 = all 72 cells with pack=0 (the full five-axis factorial on
route x scheme x partition x form x exp), plus the three cells with pack=1
at partition levels 0,1,2 and route=scheme=form=exp=0.

Rank proof (analytic):
  * the 72 pack-mean columns restricted to the 72 pack=0 rows are the
    five-axis orthonormal basis, invertible;
  * subtracting the matching pack-mean column from each of the three
    pack-contrast columns leaves a 75x75 matrix supported on the three
    pack=1 rows in the last three columns, with the 3x3 matrix
        [1, l(p_i), q(p_i)],  p_i = 0,1,2,
    which is invertible.
We also verify the rank modulo primes p = 1 mod 24.
"""

from __future__ import annotations

import json
import math

import numpy as np

import pooled_allowance as pa
from pooled_d75 import evaluate, five_axis_cells, extra_cells

P_IDX = np.flatnonzero(pa.PERM)
Q_IDX = np.flatnonzero(pa.Q)
T_IDX = pa.T_IDX
G = pa.G_ALPHA
N = pa.N


def modular_rank(p):
    # sqrt(2), sqrt(3) mod p
    r2 = next(x for x in range(2, p) if (x * x - 2) % p == 0)
    r3 = next(x for x in range(2, p) if (x * x - 3) % p == 0)
    half = pow(2, -1, p)
    r6 = (r2 * r3) % p

    def val(level_count, row, level):
        if level_count == 2:
            return 1 if row == 0 else (1 if level == 0 else -1)
        if row == 0:
            return 1
        if row == 1:
            return (r6 * half * (1 if level == 0 else (0 if level == 1 else -1))) % p
        return (r2 * half * (1 if level == 0 else (-2 if level == 1 else 1))) % p

    cells = five_axis_cells(0) + extra_cells([0, 1, 2], base=(0, 0, 0, 0))
    levels = [
        tuple((c // s) % cnt for cnt, s in zip(pa.LEVELS, pa.STR))
        for c in cells
    ]
    A = []
    for lv in levels:
        row = []
        for idx in P_IDX:
            v = 1
            for axis, ch in enumerate(pa.CHOICES[idx]):
                v = (v * val(pa.LEVELS[axis], ch, lv[axis])) % p
            row.append(v)
        A.append(row)
    # Gaussian elimination mod p
    A = np.array(A, dtype=np.int64) % p
    m = len(A)
    rank = 0
    for col in range(m):
        piv = None
        for i in range(rank, m):
            if A[i, col] % p:
                piv = i
                break
        if piv is None:
            continue
        A[[rank, piv]] = A[[piv, rank]]
        inv = pow(int(A[rank, col]), -1, p)
        A[rank] = (A[rank] * inv) % p
        for i in range(m):
            if i != rank and A[i, col] % p:
                A[i] = (A[i] - A[i, col] * A[rank]) % p
        rank += 1
        if rank == m:
            break
    return rank


def main():
    cells = five_axis_cells(0) + extra_cells([0, 1, 2], base=(0, 0, 0, 0))
    e = evaluate(cells)
    print(f"D75: N={e['N']}, numeric rank {e['rank']}/75, "
          f"cond {e['cond']:.6f}, kappa {e['kappa_max']:.6f}, "
          f"amp {e['amp_max']:.3f}")
    for p in [1009, 2017, 3001]:
        r = modular_rank(p)
        print(f"  rank mod {p}: {r}/75")
    # alias spectrum
    B = e["B"]
    ratio = np.abs(B / G[Q_IDX][None, :])
    idx = np.dstack(np.unravel_index(np.argsort(ratio, axis=None)[::-1][:12],
                                     ratio.shape))[0]
    print("\nlargest |B_alpha,beta| / G_beta entries:")
    for i, j in idx:
        a = P_IDX[i]
        b = Q_IDX[j]
        print(f"  alpha {pa.SUPPORTS[a]} {pa.CHOICES[a]}  "
              f"beta {pa.SUPPORTS[b]} {pa.CHOICES[b]}  "
              f"|B|/G = {ratio[i, j]:.6f}  B={B[i, j]:.6f}")
    out = {
        "N": e["N"], "cells": e["cells"], "rank": e["rank"],
        "cond": e["cond"], "kappa": e["kappa_max"], "amp": e["amp_max"],
        "mod_ranks": {str(p): modular_rank(p) for p in [1009, 2017, 3001]},
    }
    with open("pooled_d75_exact.json", "w") as h:
        json.dump(out, h, indent=1)
    # certified bound for the fused fp64 allowance
    etaQ64 = 3.9541e-16
    etaQ32 = 2.3758e-7
    print(f"\nfused fp64: kappa*eta_Q = {e['kappa_max']*etaQ64:.3e}; "
          f"target signal 1.329e-14; published figure 2.558e-14")
    print(f"fused fp32: kappa*eta_Q = {e['kappa_max']*etaQ32:.3e}; "
          f"target signal 5.220e-08; published figure 1.5e-07")


if __name__ == "__main__":
    main()
