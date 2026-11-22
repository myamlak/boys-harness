"""Family-aware star design on the four-axis second instance.

The record's four-axis instance has singleton cones {A_i}, so w_p = 1 and
the permitted family is the mean plus the four main effects.  The 11
order-two..four interactions are structural zeros.  A saturated star design
(baseline + one cell per axis) estimates the five permitted coefficients
with 5 runs, against 8 for the strength-2 half fraction and 16 for the full
factorial.

This script reads the recorded full-factorial cells (read-only), builds the
star design matrix, and exhibits the alias operator B = Omega_PP^{-1} Omega_PQ,
the contamination kappa, and the estimates against the recorded full-factorial
coefficients.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import os

import numpy as np

def _harness_root():
    """Locate the harness repository holding the second-instance record."""
    override = os.environ.get("BOYS_HARNESS_REPO")
    if override:
        return os.path.abspath(override)
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    probe = os.path.join(root, "composition-audit", "second-instance",
                         "results", "four-axis", "full-cells.csv")
    if os.path.isfile(probe):
        return root
    raise SystemExit("set BOYS_HARNESS_REPO to the harness repository root")


BASE = os.path.join(_harness_root(), "composition-audit", "second-instance",
                    "results", "four-axis")
N_AXES = 4
LEVELS = [2] * N_AXES


def read_csv(name):
    with open(os.path.join(BASE, name), newline="") as h:
        return list(csv.DictReader(h))


def phi(mask, levels):
    """Record convention: bit i of mask is axis A(i+1); level 0 -> +1, 1 -> -1."""
    v = 1.0
    for axis in range(N_AXES):
        if (mask >> axis) & 1:
            v *= (1.0 if levels[axis] == 0 else -1.0)
    return v


def cell_levels(cell):
    """full_cell numbering A1 A2 A3 A4, A1 most significant."""
    return [(cell >> (N_AXES - 1 - axis)) & 1 for axis in range(N_AXES)]


def full_phi():
    M = np.zeros((16, 16))
    for cell in range(16):
        levels = cell_levels(cell)
        for mask in range(16):
            M[cell, mask] = phi(mask, levels)
    return M


PHI = full_phi()
SUPPORTS = [frozenset(i for i in range(N_AXES) if (mask >> i) & 1)
            for mask in range(16)]
P = np.array([len(s) <= 1 for s in SUPPORTS])
Q = ~P
P_IDX = np.flatnonzero(P)
Q_IDX = np.flatnonzero(Q)


def design_matrix(cells):
    F = PHI[np.ix_(cells, P_IDX)]
    N = len(cells)
    OmPP = F.T @ F / N
    OmPQ = F.T @ PHI[np.ix_(cells, Q_IDX)] / N
    return F, OmPP, OmPQ


def stars_and_subsets():
    return {
        "star": [0, 8, 4, 2, 1],          # baseline + A1..A4 at level 1
        "half-fraction": [0, 3, 5, 6, 9, 10, 12, 15],  # even parity
        "full": list(range(16)),
    }


def main():
    cells_rows = read_csv("full-cells.csv")
    coeff_rows = read_csv("coefficients.csv")
    # index recorded full-factorial coefficients by mask
    full_c = {int(r["mask"]): float(r["full_factorial_coefficient"])
              for r in coeff_rows}
    error = {}
    delivered = {}
    exact = {}
    for r in cells_rows:
        fc = int(r["full_cell"])
        error[fc] = float(r["error"])
        delivered[fc] = float(r["joint_delivered"])
        exact[fc] = float(r["joint_exact"])
    print("recorded full coefficients (order <= 1):",
          {k: full_c[k] for k in range(16) if len(SUPPORTS[k]) <= 1})
    print("recorded interactions, max |c|:",
          max(abs(full_c[k]) for k in range(16) if len(SUPPORTS[k]) >= 2))

    out = {}
    for name, cells in stars_and_subsets().items():
        F, OmPP, OmPQ = design_matrix(cells)
        N = len(cells)
        rank = np.linalg.matrix_rank(F)
        inv = np.linalg.inv(OmPP)
        B = inv @ OmPQ
        # contamination kappa per permitted component; G=1 for all
        kappa = np.abs(B).max(axis=1)
        # raw contrasts and estimate from the realised joint error
        y = np.array([error[c] for c in cells])
        raw = F.T @ y / N
        est = inv @ raw
        # compare with full-factorial coefficients on P components
        truth = np.array([full_c[int(i)] for i in P_IDX])
        out[name] = {
            "cells": cells,
            "runs": N,
            "rank_P": int(rank),
            "cond": float(np.linalg.cond(F)),
            "kappa_max": float(kappa.max()),
            "kappa_P": [float(v) for v in kappa],
            "max_abs_est_minus_full": float(np.abs(est - truth).max()),
            "estimates": [float(v) for v in est],
            "truth": [float(v) for v in truth],
            "B_nonzero": int(np.sum(np.abs(B) > 1e-12)),
            "B_entries": int(B.size),
        }
        print(f"\n{name}: {N} runs, rank(P) {rank}, cond {out[name]['cond']:.3f}, "
              f"kappa_max {kappa.max():.3f}, max|est-full| "
              f"{out[name]['max_abs_est_minus_full']:.3e}")
        print("  kappa per permitted component:",
              {("mean" if i == 0 else f"A{i}"): round(float(v), 3)
               for i, v in zip(P_IDX, kappa)})
    with open("second_instance_star.json", "w") as h:
        json.dump(out, h, indent=1)


if __name__ == "__main__":
    main()
