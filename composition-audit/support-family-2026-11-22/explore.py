"""Explore the polytope of designs orthogonal to mean/main effects.

Condition B: Omega_{alpha,beta} = 0 for every target alpha (permitted
interaction, order 2..3) and every beta in R0 (mean + main effects).
The normalized weight vector m lives in a 3-dimensional affine polytope.
We enumerate its vertices and test the rank of the target Gram matrix.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.optimize import linprog

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


def build_system():
    pairs = [(a, b) for a in T_IDX for b in R0_IDX]
    A = np.array([PHI[:, a] * PHI[:, b] for a, b in pairs])
    scale = np.abs(A).max(axis=1)
    scale[scale == 0] = 1.0
    A = A / scale[:, None]
    Aeq = np.vstack([A, np.ones(N_CELLS)])
    beq = np.concatenate([np.zeros(len(A)), [1.0]])
    return Aeq, beq


AEQ, BEQ = build_system()


def gram(m):
    return (PHI * np.asarray(m)[:, None]).T @ PHI


def rank_and_alias(m):
    Om = gram(m)
    TT = Om[np.ix_(T_IDX, T_IDX)]
    ev = np.linalg.eigvalsh(TT)
    rank = np.linalg.matrix_rank(TT, tol=1e-9)
    # null vector if deficient
    nulls = []
    if rank < len(T_IDX):
        _, s, vh = np.linalg.svd(TT)
        for k in range(len(T_IDX)):
            if s[k] < 1e-9:
                nulls.append(vh[k])
    return Om, TT, ev, rank, nulls


def support_of(m, tol=1e-9):
    return [i for i in range(N_CELLS) if m[i] > tol]


def describe_design(m, label):
    nz = support_of(m)
    mpos = m[nz]
    eq = np.allclose(mpos, mpos[0], rtol=1e-6, atol=1e-9)
    N = 1.0 / mpos[0] if eq else None
    Om, TT, ev, rank, nulls = rank_and_alias(m)
    print(f"\n=== {label}: support {len(nz)}, equal weights {eq}"
          + (f", simple N={N:.6f}" if eq else "")
          + f", rank(TT)={rank}/14 ===")
    print(f"  cells: {[(i, CELLS[i], round(float(m[i]), 6)) for i in nz]}")
    print(f"  eig(TT): {np.round(ev, 6)}")
    print(f"  |TR0|max {np.abs(Om[np.ix_(T_IDX, R0_IDX)]).max():.2e}, "
          f"|TQ|max {np.abs(Om[np.ix_(T_IDX, Q_IDX)]).max():.2e}")
    for k, v in enumerate(nulls):
        print(f"  null vector {k}: " + ", ".join(
            f"{np.round(v[j], 4)}*{tuple(AXES[a] for a in sorted(SUPPORTS[T_IDX[j]])) or ('mean',)}"
            + f"{CHOICES[T_IDX[j]]}" for j in range(len(T_IDX)) if abs(v[j]) > 1e-6))
    return rank


def polytope_vertices(n_random=400, seed=1):
    rng = np.random.default_rng(seed)
    bounds = [(0.0, None)] * N_CELLS
    seen = {}
    objectives = [np.eye(N_CELLS)[i] for i in range(N_CELLS)]
    objectives += [-np.eye(N_CELLS)[i] for i in range(N_CELLS)]
    objectives += [rng.standard_normal(N_CELLS) for _ in range(n_random)]
    for c in objectives:
        res = linprog(c, A_eq=AEQ, b_eq=BEQ, bounds=bounds, method="highs")
        if not res.success:
            continue
        m = np.round(res.x, 8)
        m[m < 1e-8] = 0.0
        key = tuple(np.flatnonzero(m > 0))
        # also key on weights sign pattern
        if key not in seen:
            seen[key] = m
    return list(seen.values())


def vertex_properties(m):
    nz = np.array(support_of(m))
    mpos = m[nz]
    equal = np.allclose(mpos, mpos[0], rtol=1e-6, atol=1e-9)
    Om, TT, ev, rank, nulls = rank_and_alias(m)
    return {
        "support": len(nz),
        "equal": equal,
        "N": 1.0 / mpos[0] if equal else None,
        "rank": rank,
        "evmin": ev.min(),
        "tr0": np.abs(Om[np.ix_(T_IDX, R0_IDX)]).max(),
        "tq": np.abs(Om[np.ix_(T_IDX, Q_IDX)]).max(),
        "nulls": nulls,
    }


if __name__ == "__main__":
    print("rank AEQ:", np.linalg.matrix_rank(AEQ, tol=1e-10), "of", AEQ.shape)
    vs = polytope_vertices(n_random=600)
    print(f"{len(vs)} distinct vertices\n")
    for j, m in enumerate(vs):
        p = vertex_properties(m)
        print(f"V{j:02d}: sup={p['support']:2d} equal={p['equal']} "
              f"N={p['N'] if p['N'] else '-'} rank={p['rank']} "
              f"evmin={p['evmin']:.2e} |TR0|={p['tr0']:.1e} |TQ|={p['tq']:.2e}")
    # detailed description of every vertex with rank 14
    good = [(j, m) for j, m in enumerate(vs) if vertex_properties(m)["rank"] == 14]
    print(f"\nvertices with full target rank: {[j for j, _ in good]}")
    for j, m in good:
        describe_design(m, f"V{j:02d}")
    # whole-support uniform vertex (full cube)
    describe_design(np.ones(N_CELLS) / N_CELLS, "full cube")
