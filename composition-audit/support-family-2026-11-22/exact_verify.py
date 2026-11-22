"""Exact verification of family-aware designs on the Boys path.

Builds the exact Q(sqrt2,sqrt3) basis, screens random designs numerically,
then verifies rank, Omega_PP invertibility, exact zero pattern and kappa in
exact arithmetic.
"""

from __future__ import annotations

import json
import math

import numpy as np

from field import F4, Mat, basis_rows, ONE, ZERO
from instance import (
    N_CELLS, PHI, T_IDX, R0_IDX, Q_IDX, G_ALPHA, AXES, SUPPORTS, ORDERS,
    CHOICES, CELLS, LEVELS, STRIDES,
)

P_IDX = np.sort(np.concatenate([T_IDX, R0_IDX]))


def build_exact_phi():
    rows = [basis_rows(c) for c in LEVELS]
    phi = [[ZERO] * N_CELLS for _ in range(N_CELLS)]
    for cell_i, cell in enumerate(CELLS):
        for idx in range(N_CELLS):
            v = ONE
            for axis, row in enumerate(CHOICES[idx]):
                v = v * rows[axis][row][cell[axis]]
            phi[cell_i][idx] = v
    return Mat(phi)


EXACT_PHI = build_exact_phi()


def check_float_agreement():
    M = EXACT_PHI.to_float()
    print("exact vs float Phi max deviation:", np.abs(M - PHI).max())


def exact_design(cells):
    cells = sorted(cells)
    N = len(cells)
    F = EXACT_PHI.submatrix(cells, range(N_CELLS))
    FP = F.submatrix(range(N), P_IDX)
    FQ = F.submatrix(range(N), Q_IDX)
    rank = FP.rank()
    out = {"N": N, "cells": cells, "rank_P": rank}
    if rank < len(P_IDX):
        return out, None, None, None
    OmPP = FP.transpose().matmul(FP).scale(F4(1, 0, 0, 0) / N) if False else None
    # division by N: F4(1/N)
    invN = F4(1) / F4(N)
    OmPP = FP.transpose().matmul(FP).scale(invN)
    OmPQ = FP.transpose().matmul(FQ).scale(invN)
    invPP = OmPP.inverse()
    B = invPP.matmul(OmPQ)
    Bf = B.to_float()
    omf = OmPP.to_float()
    # kappa per target row (P ordering: T then R0? P_IDX is sorted mix)
    pos = {int(p): i for i, p in enumerate(P_IDX)}
    kappa = {}
    for t in T_IDX:
        row = np.abs(Bf[pos[int(t)], :] / G_ALPHA[Q_IDX])
        kappa[int(t)] = float(row.max())
    zero_pattern = B.zero_entries()
    out["kappa_max"] = max(kappa.values())
    out["kappa"] = kappa
    out["exact_zero_entries"] = int(sum(v for row in zero_pattern for v in row))
    out["total_entries"] = len(B.rows) * len(B.rows[0])
    out["eigPP"] = np.linalg.eigvalsh(omf)
    return out, OmPP, B, kappa


def screen_random(N, need=3, tries=200000, seed=0, min_smin=1e-3):
    rng = np.random.default_rng(seed)
    found = []
    for t in range(tries):
        cells = rng.choice(N_CELLS, size=N, replace=False)
        F = PHI[np.ix_(cells, P_IDX)]
        s = np.linalg.svd(F, compute_uv=False)
        if s.min() > min_smin * s.max():
            found.append((int(s.min() / s.max()), tuple(sorted(int(c) for c in cells))))
            if len(found) >= need:
                break
    found.sort(reverse=True)
    return found, t + 1


def main():
    check_float_agreement()
    results = {}
    for N, need, tries in [(21, 5, 200000), (22, 3, 50000), (25, 3, 20000)]:
        found, used = screen_random(N, need=need, tries=tries, seed=100 + N)
        print(f"N={N}: screened {used} random designs, found {len(found)}")
        best = None
        for cond, cells in found:
            out, OmPP, B, kappa = exact_design(cells)
            if out["rank_P"] != 21:
                print(f"  candidate {cells}: exact rank {out['rank_P']} (reject)")
                continue
            print(f"  candidate {cells}: exact rank 21, kappa_max {out['kappa_max']:.6f}, "
                  f"cond {out['eigPP'].max()/out['eigPP'].min():.3f}, "
                  f"exact-zero B entries {out['exact_zero_entries']}/{out['total_entries']}")
            if best is None or out["kappa_max"] < best[0]["kappa_max"]:
                best = (out, cells)
        results[N] = best
        if best:
            print(f"  BEST N={N}: kappa {best[0]['kappa_max']:.6f} cells {best[1]}")
    with open("exact_designs.json", "w") as h:
        json.dump({str(k): (v[1] if v else None) for k, v in results.items()}, h)

    # verify the two 18-cell exact-orthogonality designs
    d18 = [
        [0, 2, 4, 7, 9, 11, 12, 14, 16, 19, 21, 23, 24, 26, 28, 31, 33, 35],
        [1, 3, 5, 6, 8, 10, 13, 15, 17, 18, 20, 22, 25, 27, 29, 30, 32, 34],
    ]
    print("\n18-cell orthogonal-to-R0 designs:")
    for cells in d18:
        out, _, _, _ = exact_design_exact_only(cells)
        print(f"  cells {cells}: rank_P {out['rank_P']}, rank_TT {out['rank_TT']}")

    # full cube
    out, _, _, _ = exact_design(list(range(N_CELLS)))
    print(f"\nfull cube: rank_P {out['rank_P']}, kappa_max {out['kappa_max']:.3e}")


def exact_design_exact_only(cells):
    cells = sorted(cells)
    N = len(cells)
    F = EXACT_PHI.submatrix(cells, range(N_CELLS))
    FP = F.submatrix(range(N), P_IDX)
    FT = F.submatrix(range(N), T_IDX)
    rankP = FP.rank()
    rankT = FT.rank()
    return {"N": N, "rank_P": rankP, "rank_TT": rankT}, None, None, None


if __name__ == "__main__":
    main()
