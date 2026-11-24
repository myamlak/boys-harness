"""Exact mixed-level instance for the Boys resolving fp64/fp32 path.

Axes and levels, in the order (partition, pack, form, exp):
    partition / granularity : 3
    pack                    : 2
    form / division form    : 3
    exp / region-B expon.   : 2

Cone table on the path:
    upward step  C1 = {partition, form, exp} = {0,2,3}
    packed sum   C2 = {partition, pack}       = {0,1}
    region-A coefficient sum {partition} = {0} (subset of C1)
    far ladder, sort: empty cones (irrelevant)

Basis: orthonormal under the uniform measure, per axis_basis().
"""

from __future__ import annotations

import itertools
import math

import numpy as np

LEVELS = [3, 2, 3, 2]
AXES = ["partition", "pack", "form", "exp"]
CONES = [
    {"partition", "form", "exp"},   # upward step
    {"partition", "pack"},          # packed single-precision sum
]
CONE_IDX = [frozenset({0, 2, 3}), frozenset({0, 1})]
W_P = 3


def axis_basis(level_count):
    if level_count == 2:
        rows = np.array([[1.0, 1.0], [1.0, -1.0]])
    elif level_count == 3:
        rows = np.array(
            [
                [1.0, 1.0, 1.0],
                math.sqrt(1.5) * np.array([1.0, 0.0, -1.0]),
                (1.0 / math.sqrt(2.0)) * np.array([1.0, -2.0, 1.0]),
            ]
        )
    else:
        raise ValueError(level_count)
    return rows


def strides(level_counts):
    out, s = [], 1
    for c in reversed(level_counts):
        out.append(s)
        s *= c
    return list(reversed(out))


N_CELLS = math.prod(LEVELS)
STRIDES = strides(LEVELS)
ROWS = [axis_basis(c) for c in LEVELS]
# convenience: sups of |row| per axis
ROW_SUP = [np.abs(r).max(axis=1) for r in ROWS]


def cell_levels(cell):
    return tuple((cell // s) % c for c, s in zip(LEVELS, STRIDES))


CELLS = [cell_levels(c) for c in range(N_CELLS)]
CELL_LEVEL_ARR = np.array(CELLS, dtype=np.int64)


def component_choices(index):
    """Multi-index (one basis row per axis) of basis component `index`."""
    return tuple((index // s) % c for c, s in zip(LEVELS, STRIDES))


def support_axes(choices):
    return frozenset(i for i, ch in enumerate(choices) if ch != 0)


def build_phi():
    phi = np.ones((N_CELLS, N_CELLS))
    for k, cell in enumerate(CELLS):
        for idx in range(N_CELLS):
            ch = component_choices(idx)
            v = 1.0
            for axis, row in enumerate(ch):
                v *= ROWS[axis][row, cell[axis]]
            phi[k, idx] = v
    return phi


PHI = build_phi()
CHOICES = [component_choices(i) for i in range(N_CELLS)]
SUPPORTS = [support_axes(c) for c in CHOICES]
ORDERS = [len(s) for s in SUPPORTS]
G_ALPHA = np.array(
    [
        math.prod(ROW_SUP[axis][row] for axis, row in enumerate(ch))
        for ch in CHOICES
    ]
)


def in_permitted(supp):
    return any(supp <= cone for cone in CONE_IDX)


PERMITTED = np.array([in_permitted(s) for s in SUPPORTS])
TARGET = np.array(
    [PERMITTED[i] and 2 <= ORDERS[i] <= W_P for i in range(N_CELLS)]
)
R0 = np.array([PERMITTED[i] and ORDERS[i] <= 1 for i in range(N_CELLS)])
Q = np.array([not PERMITTED[i] for i in range(N_CELLS)])

T_IDX = np.flatnonzero(TARGET)
R0_IDX = np.flatnonzero(R0)
Q_IDX = np.flatnonzero(Q)

assert len(T_IDX) == 14, len(T_IDX)
assert len(R0_IDX) == 7, len(R0_IDX)
assert len(Q_IDX) == 15, len(Q_IDX)
assert len(T_IDX) + len(R0_IDX) + len(Q_IDX) == N_CELLS


def describe():
    print(f"cells={N_CELLS} permitted components={PERMITTED.sum()} "
          f"(targets={len(T_IDX)}, low-order permitted={len(R0_IDX)}) "
          f"structural zeros={len(Q_IDX)}")
    for name, idxs in (("T", T_IDX), ("R0", R0_IDX), ("Q", Q_IDX)):
        print(f"{name}:")
        for i in idxs:
            print(
                f"  idx {i:2d} order {ORDERS[i]} support "
                f"{[AXES[a] for a in sorted(SUPPORTS[i])] or ['mean']} "
                f"choices {CHOICES[i]} G={G_ALPHA[i]:.6f}"
            )
    # channel counts per support
    from collections import Counter

    counts = Counter(SUPPORTS[i] for i in T_IDX)
    print("target channel counts per support:")
    for supp, n in sorted(counts.items(), key=lambda kv: (len(kv[0]), sorted(kv[0]))):
        print(f"  {[AXES[a] for a in sorted(supp)]}: {n}")


if __name__ == "__main__":
    describe()
    # full-cube sanity
    print("Phi^T Phi / N = I:", np.allclose(PHI.T @ PHI / N_CELLS, np.eye(N_CELLS)))
    print("G range:", G_ALPHA.min(), G_ALPHA.max())
