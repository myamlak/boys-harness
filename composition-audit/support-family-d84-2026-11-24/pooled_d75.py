"""The 75-run saturated design on the pooled six-axis cube.

The 72 permitted components that do not involve the pack contrast are the
complete basis on the five axes of the pooled upward-step cone
(route, scheme, partition, form, exp).  Their restrictions are independent
only on the full 72-cell five-axis factorial.  The remaining three permitted
components carry the pack contrast (pack alone, and partition linear/quadratic
times pack).  Three cells at the opposite packing level, at three partition
levels, identify them.

So D75 = {(u, pack=0) : u in the 72-cell five-axis cube}
         + {(u_i, pack=1) : i = 1,2,3}.

We verify rank, kappa and the amplification, and compare with the full 144.
"""

from __future__ import annotations

import json
import math

import numpy as np

import pooled_allowance as pa

PHI = pa.PHI
P_IDX = np.flatnonzero(pa.PERM)
Q_IDX = np.flatnonzero(pa.Q)
T_IDX = pa.T_IDX
G = pa.G_ALPHA
LEVELS = pa.LEVELS
STR = pa.STR


def cell_index(route, scheme, partition, pack, form, exp):
    return (route * 72 + scheme * 36 + partition * 12 + pack * 6
            + form * 2 + exp)


def five_axis_cells(pack_level):
    return [cell_index(r, s, p, pack_level, f, e)
            for r in range(2) for s in range(2) for p in range(3)
            for f in range(3) for e in range(2)]


def evaluate(cells):
    cells = np.sort(np.asarray(cells, dtype=int))
    N = len(cells)
    F = PHI[cells, :]
    PP = F[:, P_IDX]
    QQ = F[:, Q_IDX]
    s = np.linalg.svd(PP, compute_uv=False)
    rank = int((s > 1e-9 * max(1.0, s.max())).sum())
    OmPP = PP.T @ PP / N
    OmPQ = PP.T @ QQ / N
    inv = np.linalg.inv(OmPP)
    B = inv @ OmPQ
    kappa = np.abs(B / G[Q_IDX][None, :])
    amp = np.abs(inv) @ G[P_IDX]
    pos = {int(p): i for i, p in enumerate(P_IDX)}
    return {
        "N": N, "cells": [int(c) for c in cells], "rank": rank,
        "cond": float(s.max() / s.min()) if s.min() > 0 else float("inf"),
        "kappa_max": float(kappa.max()),
        "kappa_targets_max": float(kappa[[pos[int(t)] for t in T_IDX]].max()),
        "amp_max": float(amp.max()),
        "amp_targets_max": float(amp[[pos[int(t)] for t in T_IDX]].max()),
        "OmPP": OmPP, "B": B,
    }


def extra_cells(partition_levels, base=(0, 0, 0, 0)):
    """Cells at pack=1 for the given partition levels, fixed other axes."""
    r, s, f, e = base
    return [cell_index(r, s, p, 1, f, e) for p in partition_levels]


def main():
    base72 = five_axis_cells(0)
    print("five-axis factorial cells:", len(base72))
    full = evaluate(list(range(144)))
    print(f"full 144: rank {full['rank']}, kappa {full['kappa_max']:.3f}, "
          f"amp {full['amp_max']:.3f}")
    results = {}
    for extra in [
        extra_cells([0, 1, 2], base=(0, 0, 0, 0)),
        extra_cells([0, 1, 2], base=(0, 0, 0, 1)),
        extra_cells([0, 1, 2], base=(0, 0, 2, 0)),
        extra_cells([0, 2, 1], base=(0, 0, 0, 0)),
    ]:
        e = evaluate(base72 + extra)
        tag = tuple(extra)
        results[str(extra)] = {k: v for k, v in e.items()
                               if k not in ("OmPP", "B")}
        print(f"\nD75 extra {extra}:")
        print(f"  rank {e['rank']}/75, cond {e['cond']:.3f}, "
              f"kappa {e['kappa_max']:.4f} (targets {e['kappa_targets_max']:.4f}), "
              f"amp {e['amp_max']:.3f} (targets {e['amp_targets_max']:.3f})")
    # add more pack=1 cells to see the trade-off
    for k in [1, 2, 3, 4, 6, 9, 12]:
        extras = five_axis_cells(1)[:k]
        e = evaluate(base72 + extras)
        print(f"\nN={72 + k}: extra {k} pack-1 cells: "
              f"rank {e['rank']}, kappa {e['kappa_max']:.4f}, "
              f"amp {e['amp_max']:.3f}, cond {e['cond']:.3f}")
    with open("pooled_d75.json", "w") as h:
        json.dump(results, h, indent=1)


if __name__ == "__main__":
    main()
