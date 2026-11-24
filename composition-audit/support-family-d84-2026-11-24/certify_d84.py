"""Certificate for the 84-run pooled design D84.

Axes (pooled six-axis cube): route, scheme, partition, pack, form, exp.
The corrected packed-sum cone is {route, scheme, partition, pack}; the upward
cone is {route, scheme, partition, form, exp}.  Their intersection is the
three partition rows times route and scheme, so the permitted set is

    (2*2*3*3*2) + (2*2*3*2) - (2*2*3) = 72 + 24 - 12 = 84 channels.

D84 = every cell with pack = arguments (the full 72-cell five-axis factorial on
route x scheme x partition x form x exp), plus the 12 cells with pack = orders
at route in {0,1}, scheme in {0,1}, partition in {0,1,2}, form = exp = 0.

Rank proof:
  * the 72 pack-mean columns on the 72 pack=arguments rows are the five-axis
    orthonormal basis, invertible;
  * subtracting the matching pack-mean column from each of the twelve
    pack-contrast columns leaves a matrix supported on the twelve
    pack=orders rows, whose 12x12 block is the Kronecker product of the
    route (2x2), scheme (2x2) and partition (3x3) orthonormal bases evaluated
    at the selected levels.  Each factor is invertible, so the block is.
The rank is also checked modulo three primes.

The printed kappa is the target-restricted alias constant, written as
`kappa_target`.  It is not the joint contamination constant; that constant is
four and is verified over the full cube by verify_cert.py.
"""

from __future__ import annotations

import json

import numpy as np

import pooled_allowance as pa
import pooled_d75 as d75

# Correct the packed cone: route and scheme reach the packed sum.
pa.PACKED = {0, 1, 2, 3}
pa.PERM = np.array([pa.in_pooled(s) for s in pa.SUPPORTS])
pa.Q = ~pa.PERM
pa.TARGET = np.array([pa.PERM[i] and 2 <= pa.ORDERS[i] <= pa.W_P
                      for i in range(pa.N)])
pa.Q_IDX = np.flatnonzero(pa.Q)
pa.T_IDX = np.flatnonzero(pa.TARGET)
d75.P_IDX = np.flatnonzero(pa.PERM)
d75.Q_IDX = pa.Q_IDX
d75.T_IDX = pa.T_IDX

P_IDX = d75.P_IDX
Q_IDX = d75.Q_IDX
T_IDX = d75.T_IDX
G = pa.G_ALPHA

CELLS = sorted(d75.five_axis_cells(0) + [
    d75.cell_index(r, s, p, 1, 0, 0)
    for r in range(2) for s in range(2) for p in range(3)])


def modular_rank(p):
    r2 = next(x for x in range(2, p) if (x * x - 2) % p == 0)
    r3 = next(x for x in range(2, p) if (x * x - 3) % p == 0)
    half = pow(2, -1, p)
    r6 = (r2 * r3) % p

    def val(count, row, level):
        if count == 2:
            return 1 if row == 0 else (1 if level == 0 else -1)
        if row == 0:
            return 1
        if row == 1:
            return (r6 * half * (1 if level == 0 else
                                 (0 if level == 1 else -1))) % p
        return (r2 * half * (1 if level == 0 else
                             (-2 if level == 1 else 1))) % p

    levels = [tuple((c // s) % cnt for cnt, s in zip(pa.LEVELS, pa.STR))
              for c in CELLS]
    a = []
    for lv in levels:
        row = []
        for idx in P_IDX:
            v = 1
            for axis, ch in enumerate(pa.CHOICES[idx]):
                v = (v * val(pa.LEVELS[axis], ch, lv[axis])) % p
            row.append(v)
        a.append(row)
    a = np.array(a, dtype=np.int64) % p
    m = len(a)
    rank = 0
    for col in range(m):
        piv = next((i for i in range(rank, m) if a[i, col] % p), None)
        if piv is None:
            continue
        a[[rank, piv]] = a[[piv, rank]]
        inv = pow(int(a[rank, col]), -1, p)
        a[rank] = (a[rank] * inv) % p
        for i in range(m):
            if i != rank and a[i, col] % p:
                a[i] = (a[i] - a[i, col] * a[rank]) % p
        rank += 1
        if rank == m:
            break
    return rank


def main():
    print(f"permitted {len(P_IDX)}, complement {len(Q_IDX)}, "
          f"targets {len(T_IDX)}, cells {len(CELLS)}")
    e = d75.evaluate(CELLS)
    print(f"rank {e['rank']}/{len(CELLS)}, cond {e['cond']:.12f}, "
          f"kappa {e['kappa_max']:.15f}, amp {e['amp_max']:.12f}")
    mods = {}
    for p in (1009, 2017, 3001):
        mods[str(p)] = modular_rank(p)
        print(f"rank mod {p}: {mods[str(p)]}/{len(CELLS)}")
    pos = {int(p): i for i, p in enumerate(P_IDX)}
    b = e["B"]
    gq = G[Q_IDX]
    gp = G[P_IDX]
    ratio = np.abs(b / gq[None, :])
    g_kappa = []
    for t in T_IDX:
        i = pos[int(t)]
        g_kappa.append(float(gp[i]) * float(ratio[i].max()))
    print(f"sum G kappa {sum(g_kappa):.15f}, "
          f"min {min(g_kappa):.15f}, max {max(g_kappa):.15f}")
    out = {
        "convention": "packed cone {route, scheme, partition, pack}",
        "permitted": int(len(P_IDX)),
        "complement": int(len(Q_IDX)),
        "targets": int(len(T_IDX)),
        "cells": [int(c) for c in CELLS],
        "rank": int(e["rank"]),
        "cond": float(e["cond"]),
        "kappa_target": float(e["kappa_max"]),
        "amplification": float(e["amp_max"]),
        "mod_ranks": mods,
        "G_kappa": g_kappa,
        "sum_G_kappa": float(sum(g_kappa)),
    }
    with open("certify_d84.json", "w") as h:
        json.dump(out, h, indent=1)


if __name__ == "__main__":
    main()
