"""File-backed certificate for the structured 21-run pathwise design D21.

D21 construction (axes: 0 partition/granularity, 1 pack, 2 form, 3 exp):
  * the full 3x3 partition x form grid at (pack, exp) = (0, 0);
  * the full 3x3 partition x form grid at (pack, exp) = (1, 1);
  * the three cells at (pack, exp) = (0, 1) with partition = 0,1,2 and
    form = 0.

The script prints the exact cells and levels, the exact rank of the 21x21
permitted block over Q(sqrt2, sqrt3), kappa_max, and the fourteen
G_alpha kappa_alpha values.  It writes certify_d21.json.
"""

from __future__ import annotations

import json

import numpy as np

import exact_verify as ev
from instance import AXES, CELLS, G_ALPHA, SUPPORTS, CHOICES, T_IDX


def cell(partition, pack, form, exp):
    return partition * 12 + pack * 6 + form * 2 + exp


D21 = (
    [cell(p, 0, f, 0) for p in range(3) for f in range(3)]      # (0,0) grid
    + [cell(p, 1, f, 1) for p in range(3) for f in range(3)]    # (1,1) grid
    + [cell(p, 0, 0, 1) for p in range(3)]                      # (0,1) extras
)
D21 = sorted(D21)


def main():
    print("D21 cells:", D21)
    print("D21 cell count:", len(D21))
    for c in D21:
        print(f"  cell {c:2d} levels (partition,pack,form,exp) = {CELLS[c]}")
    out, OmPP, B, kappa = ev.exact_design(D21)
    print(f"\nexact rank of the 21x21 permitted block over Q(sqrt2,sqrt3): "
          f"{out['rank_P']} / 21")
    print(f"kappa_max = {out['kappa_max']:.15f}")
    M = ev.EXACT_PHI.submatrix(D21, range(36)).to_float()
    F = M[:, ev.P_IDX]
    invPP = OmPP.inverse().to_float()
    amp = float((np.abs(invPP) @ ev.G_ALPHA[ev.P_IDX]).max())
    print(f"cond(F_P) = {np.linalg.cond(F):.9f}")
    print(f"threshold amplification max_alpha sum_gamma |invPP| G = {amp:.6f}")
    print(f"eigenvalues of Omega_PP: {np.round(out['eigPP'], 9).tolist()}")
    print("\nthe fourteen G_alpha kappa_alpha values:")
    vals = []
    for t in T_IDX:
        g = float(G_ALPHA[t])
        k = float(kappa[int(t)])
        vals.append(g * k)
        print(f"  target support {[AXES[a] for a in sorted(SUPPORTS[t])]} "
              f"choices {CHOICES[t]}  G={g:.12f}  kappa={k:.12f}  "
              f"G*kappa={g*k:.15f}")
    print(f"\nsum G_alpha kappa_alpha = {sum(vals):.15f}")
    print(f"min = {min(vals):.15f}, max = {max(vals):.15f}")
    with open("certify_d21.json", "w") as h:
        json.dump({
            "cells": D21,
            "rank": out["rank_P"],
            "kappa_max": out["kappa_max"],
            "G_kappa": vals,
            "sum_G_kappa": sum(vals),
            "amp": amp,
            "eigPP": out["eigPP"].tolist(),
        }, h, indent=1)


if __name__ == "__main__":
    main()
