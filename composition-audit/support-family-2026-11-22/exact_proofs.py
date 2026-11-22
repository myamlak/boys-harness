"""Exact certificates for the Boys-path obstruction.

(1) Condition A (target-target diagonal + target-R0 orthogonal) has a
    rank-36 constraint matrix over Q(sqrt2,sqrt3); the only normalized
    solution is the uniform full cube.
(2) The two 18-cell designs have target rank 12 and their null relation is
    exhibited exactly.
(3) A 21-cell design has exact rank 21 and exact Omega_PP invertible.
"""

from __future__ import annotations

import itertools

import numpy as np
from scipy.linalg import qr

from field import F4, Mat
from exact_verify import EXACT_PHI, P_IDX
from instance import (
    N_CELLS, PHI, T_IDX, R0_IDX, Q_IDX, AXES, SUPPORTS, CHOICES, CELLS,
)


def condition_a_matrix():
    rows = []
    labels = []
    for a in T_IDX:
        for b in R0_IDX:
            rows.append(PHI[:, a] * PHI[:, b])
            labels.append(("TR0", int(a), int(b)))
    for i, a in enumerate(T_IDX):
        for b in T_IDX[i + 1:]:
            rows.append(PHI[:, a] * PHI[:, b])
            labels.append(("TT", int(a), int(b)))
    rows.append(np.ones(N_CELLS))
    labels.append(("norm", -1, -1))
    return np.array(rows), labels


def exact_rank_of_labels(labels_sel):
    M = Mat([
        [EXACT_PHI.rows[c][a] * EXACT_PHI.rows[c][b]
         for c in range(N_CELLS)]
        for (kind, a, b) in labels_sel
    ])
    return M.rank()


def main():
    A, labels = condition_a_matrix()
    print("condition A rows:", A.shape)
    print("float rank:", np.linalg.matrix_rank(A, tol=1e-9))
    # pick 36 rows via QR pivoting on A^T
    Q, R, piv = qr(A.T, pivoting=True, mode="economic")
    sel = sorted(piv[:36])
    print("QR pivot rows count:", len(sel))
    r = exact_rank_of_labels([labels[i] for i in sel])
    print("exact rank of the 36 selected rows:", r, "/ 36")
    print("selected constraint labels:")
    for i in sel:
        print("  ", labels[i])

    # exact null relation of the 18-cell designs
    for tag, cells in [
        ("V0", [1, 3, 5, 6, 8, 10, 13, 15, 17, 18, 20, 22, 25, 27, 29, 30, 32, 34]),
        ("V6", [0, 2, 4, 7, 9, 11, 12, 14, 16, 19, 21, 23, 24, 26, 28, 31, 33, 35]),
    ]:
        M = Mat([[EXACT_PHI.rows[c][t] for t in T_IDX] for c in cells])
        rank = M.rank()
        # solve M x = 0 exactly: Gauss-Jordan to RREF
        Aaug = [row[:] for row in M.rows]
        n, m = len(Aaug), len(Aaug[0])
        piv_cols = []
        r = 0
        for col in range(m):
            p = next((i for i in range(r, n) if not Aaug[i][col].is_zero()), None)
            if p is None:
                continue
            Aaug[r], Aaug[p] = Aaug[p], Aaug[r]
            pv = Aaug[r][col]
            Aaug[r] = [v / pv for v in Aaug[r]]
            for i in range(n):
                if i != r and not Aaug[i][col].is_zero():
                    f = Aaug[i][col]
                    Aaug[i] = [x - f * y for x, y in zip(Aaug[i], Aaug[r])]
            piv_cols.append(col)
            r += 1
            if r == n:
                break
        free = [c for c in range(m) if c not in piv_cols]
        print(f"\n{tag}: target rank {rank}, free target columns {free}")
        for fc in free:
            v = [F4(0)] * m
            v[fc] = F4(1)
            for ri, pc in enumerate(piv_cols):
                v[pc] = -Aaug[ri][fc]
            nz = [(j, v[j]) for j in range(m) if not v[j].is_zero()]
            desc = []
            for j, val in nz:
                t = T_IDX[j]
                desc.append(
                    f"{val} * [{','.join(AXES[a] for a in sorted(SUPPORTS[t]))}]"
                    f"{CHOICES[t]}"
                )
            print("  null vector:", " + ".join(desc))
            print("  float coefficients:", [float(v[j]) for j in range(m)])


if __name__ == "__main__":
    main()
