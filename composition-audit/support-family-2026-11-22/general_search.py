"""General support-family-aware design search.

Model P = all 21 components with support in a cone (mean + 6 mains +
14 interactions); Q = the 15 structural zeros.  For a 0/1 design D with
N cells:

    Omega_PP = (1/N) Phi_D,P^T Phi_D,P       (must be invertible)
    Omega_PQ = (1/N) Phi_D,P^T Phi_D,Q
    B = Omega_PP^{-1} Omega_PQ
    kappa_alpha = max_{beta in Q} |B_{alpha,beta}| / G_beta

The certified bound is delta_alpha + kappa_alpha * eta, where delta_alpha
is the measurement threshold propagated through Omega_PP^{-1}.
"""

from __future__ import annotations

import itertools
import time

import numpy as np

from instance import (
    N_CELLS,
    PHI,
    T_IDX,
    R0_IDX,
    Q_IDX,
    G_ALPHA,
    AXES,
    SUPPORTS,
    CHOICES,
    ORDERS,
    CELLS,
)

P_IDX = np.sort(np.concatenate([T_IDX, R0_IDX]))
assert len(P_IDX) == 21


def evaluate(cells):
    """Return rank, target kappa vector, and diagnostics for a cell subset."""
    cells = np.sort(np.asarray(cells, dtype=int))
    N = len(cells)
    F = PHI[cells, :]
    PP = F[:, P_IDX]
    QQ = F[:, Q_IDX]
    OmPP = PP.T @ PP / N
    OmPQ = PP.T @ QQ / N
    # condition / rank
    s = np.linalg.svd(PP, compute_uv=False)
    rank = int((s > 1e-9 * max(1.0, s.max())).sum())
    if rank < len(P_IDX):
        return {
            "N": N,
            "cells": cells,
            "rank": rank,
            "kappa": None,
            "kappa_max": None,
            "kappa_mean": None,
            "smin": s.min(),
            "smax": s.max(),
            "cond": None,
            "amp": None,
            "amp_target_max": None,
            "omPQ_max": None,
        }
    B = np.linalg.solve(OmPP, OmPQ)
    kappa_all = np.abs(B / G_ALPHA[Q_IDX][None, :]).max(axis=1)
    # mapping from P_IDX position to target position
    pos = {p: i for i, p in enumerate(P_IDX)}
    kappa_t = {int(p): float(kappa_all[pos[p]]) for p in T_IDX}
    # noise amplification: max row l1 of OmPP^{-1} * mean |Phi_P|
    inv = np.linalg.inv(OmPP)
    amp = np.abs(inv).sum(axis=1)
    return {
        "N": N,
        "cells": cells,
        "rank": rank,
        "kappa": kappa_t,
        "kappa_max": max(kappa_t.values()),
        "kappa_mean": float(np.mean(list(kappa_t.values()))),
        "smin": s.min(),
        "smax": s.max(),
        "cond": s.max() / s.min(),
        "amp": amp,
        "amp_target_max": float(amp[[pos[p] for p in T_IDX]].max()),
        "omPQ_max": np.abs(OmPQ).max(),
    }


def full_cube_minus(remove):
    remove = set(remove)
    return [a for a in range(N_CELLS) if a not in remove]


def search_single_removal():
    out = []
    for a in range(N_CELLS):
        e = evaluate(full_cube_minus([a]))
        out.append((e["kappa_max"], a, e))
    out.sort()
    return out


def search_pairs():
    out = []
    bad = 0
    for a, b in itertools.combinations(range(N_CELLS), 2):
        e = evaluate(full_cube_minus([a, b]))
        if e["kappa_max"] is None:
            bad += 1
            continue
        out.append((e["kappa_max"], (a, b), e))
    out.sort()
    return out, bad


def local_search(N, restarts=30, iters=1500, seed=0, kappa_weight=1.0,
                 amp_weight=0.1):
    """Simulated annealing over subsets of size N; cost = kappa + amp weight."""
    rng = np.random.default_rng(seed)
    best = None
    for r in range(restarts):
        cells = set(rng.choice(N_CELLS, size=N, replace=False).tolist())
        cur = evaluate(sorted(cells))
        if cur["rank"] < 21:
            cost = np.inf
        else:
            cost = kappa_weight * cur["kappa_max"] + amp_weight * cur["amp_target_max"]
        T = 0.5
        for it in range(iters):
            inn = int(rng.integers(0, N_CELLS))
            if inn in cells:
                continue
            trial = set(cells)
            trial.remove(rng.choice(list(cells)))
            trial.add(inn)
            ev = evaluate(sorted(trial))
            if ev["rank"] < 21:
                cst = np.inf
            else:
                cst = kappa_weight * ev["kappa_max"] + amp_weight * ev["amp_target_max"]
            if cst < cost or (np.isfinite(cst) and rng.random() < np.exp((cost - cst) / max(T, 1e-9))):
                cells, cur, cost = trial, ev, cst
            T *= 0.998
            if best is None or (np.isfinite(cost) and cost < best[0]):
                if cur["rank"] == 21:
                    best = (cost, cur)
    return best


if __name__ == "__main__":
    print("=== single removals (N=35) ===")
    for k, a, e in search_single_removal()[:10]:
        print(f"remove {a:2d} {CELLS[a]}: kappa_max {k:.4f} "
              f"amp {e['amp_target_max']:.3f} cond {e['cond']:.2f} "
              f"|PQ| {e['omPQ_max']:.3f}")
    print("\n=== pairs (N=34), best 15 ===")
    t0 = time.time()
    pr, bad = search_pairs()
    print(f"pairs evaluated in {time.time()-t0:.1f}s; "
          f"{bad} rank-deficient")
    for k, ab, e in pr[:15]:
        print(f"remove {ab}: kappa_max {k:.4f} amp {e['amp_target_max']:.3f} "
              f"cond {e['cond']:.2f}")
    print("\n=== local search by N ===")
    for N in range(21, 36):
        best = local_search(N, restarts=12, iters=600, seed=1000 + N)
        if best is None:
            print(f"N={N}: no full-rank design found")
        else:
            c, e = best
            print(f"N={N}: cost {c:.3f} kappa_max {e['kappa_max']:.4f} "
                  f"amp {e['amp_target_max']:.3f} cond {e['cond']:.2f}")
