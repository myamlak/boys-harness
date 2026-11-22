"""File-backed certificate for the 75-run pooled design D75.

Axes (pooled six-axis cube): route, scheme, partition, pack, form, exp.
D75 = every cell with pack = arguments (the full five-axis factorial on
route x scheme x partition x form x exp, 72 cells), plus the three cells
with pack = orders at partition = 0,1,2 and route = scheme = form = exp = 0.

Rank proof:
  * the 72 pack-mean columns on the 72 pack=arguments cells are the
    five-axis orthonormal basis, invertible;
  * subtracting the matching pack-mean column from each of the three
    pack-contrast columns leaves a matrix supported on the three pack=orders
    rows with the 3x3 matrix [1, linear(p), quadratic(p)] at p = 0,1,2,
    whose determinant is -3*sqrt(3).
The rank is also checked modulo three primes p = 1 mod 24.
"""

from __future__ import annotations

import json

import numpy as np

import pooled_allowance as pa
import pooled_d75 as d75
from pooled_d75_exact import modular_rank

P_IDX = np.flatnonzero(pa.PERM)
Q_IDX = np.flatnonzero(pa.Q)
T_IDX = pa.T_IDX
G = pa.G_ALPHA
AXES = pa.AXES


def main():
    cells = sorted(d75.five_axis_cells(0)
                   + d75.extra_cells([0, 1, 2], base=(0, 0, 0, 0)))
    print("D75 construction: all 72 cells with pack = arguments, plus "
          "pack = orders at partition 0,1,2 with route=scheme=form=exp=0")
    print("pack=orders cells:", [c for c in cells if c >= 0 and pa.CHOICES[c][3] == 1])
    print("cell count:", len(cells))
    print("cells:", cells)
    e = d75.evaluate(cells)
    print(f"\nnumeric rank of the 75x75 permitted block: {e['rank']} / 75")
    print(f"condition number: {e['cond']:.9f}")
    print(f"kappa_max = {e['kappa_max']:.15f}")
    print(f"threshold amplification max_alpha sum_gamma |invPP| G = "
          f"{e['amp_max']:.6f}")
    for p in [1009, 2017, 3001]:
        print(f"rank mod {p}: {modular_rank(p)} / 75")
    B = e["B"]
    GQ = G[Q_IDX]
    GP = G[P_IDX]
    ratio = np.abs(B / GQ[None, :])
    pos = {int(p): i for i, p in enumerate(P_IDX)}
    vals = []
    print("\nthe sixty-six G_alpha kappa_alpha values (target: value):")
    for t in T_IDX:
        i = pos[int(t)]
        k = float(ratio[i].max())
        vals.append(float(GP[i]) * k)
        supp = [AXES[a] for a in sorted(pa.SUPPORTS[t])]
        print(f"  {supp} {pa.CHOICES[t]}: {float(GP[i]) * k:.15f}")
    print(f"\nsum G_alpha kappa_alpha = {sum(vals):.15f}")
    print(f"min = {min(vals):.15f}, max = {max(vals):.15f}")
    print("\nlargest |B_alpha,beta|/G_beta entries:")
    flat = np.argsort(ratio, axis=None)[::-1][:8]
    for fi in flat:
        i, j = np.unravel_index(fi, ratio.shape)
        print(f"  alpha {pa.SUPPORTS[P_IDX[i]]} {pa.CHOICES[P_IDX[i]]}, "
              f"beta {pa.SUPPORTS[Q_IDX[j]]} {pa.CHOICES[Q_IDX[j]]}, "
              f"ratio {ratio[i, j]:.12f}, B {B[i, j]:.12f}")
    with open("certify_d75.json", "w") as h:
        json.dump({
            "cells": cells,
            "rank": e["rank"],
            "mod_ranks": {str(p): modular_rank(p) for p in [1009, 2017, 3001]},
            "kappa_max": e["kappa_max"],
            "G_kappa": vals,
            "sum_G_kappa": sum(vals),
            "amp": e["amp_max"],
            "cond": e["cond"],
        }, h, indent=1)


if __name__ == "__main__":
    main()
