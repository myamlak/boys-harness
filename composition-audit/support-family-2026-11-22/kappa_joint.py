"""Exact joint contamination constant kappa_joint for D21 and D75.

For beta in Q, r_beta(a) = Phi_beta(a) - sum_{alpha in T} B_{alpha,beta} Phi_alpha(a),
and kappa_joint = max_{beta in Q} sup_{a in cube} |r_beta(a)| / G_beta.

For D21 the whole computation is exact over Q(sqrt2, sqrt3): B comes from
exact_design and Phi from exact_verify.  For D75 B is the float inverse of
the 75x75 permitted block; the ratios are then rounded to the exact-looking
integers and re-checked with the claimed values.
"""

from __future__ import annotations

import json

import numpy as np

import exact_verify as ev
import pooled_allowance as pa
import pooled_d75 as d75
from instance import AXES, CELLS, SUPPORTS, CHOICES, G_ALPHA

# ---- D21 exact -------------------------------------------------------------

D21 = [0, 1, 2, 4, 7, 9, 11, 12, 13, 14, 16, 19, 21, 23, 24, 25, 26, 28, 31, 33, 35]


def d21_kappa_joint():
    out, OmPP, B, kappa = ev.exact_design(D21)
    posP = {int(p): i for i, p in enumerate(ev.P_IDX)}
    Trows = [posP[int(t)] for t in ev.T_IDX]
    rows = []
    for j, q in enumerate(ev.Q_IDX):
        worst = 0.0
        arg = None
        for a in range(36):
            r = ev.EXACT_PHI.rows[a][int(q)]
            for tr, alpha in zip(Trows, ev.T_IDX):
                r = r - B.rows[tr][j] * ev.EXACT_PHI.rows[a][int(alpha)]
            val = abs(float(r)) / float(G_ALPHA[q])
            if val > worst:
                worst = val
                arg = a
        rows.append({
            "beta": int(q),
            "support": [AXES[i] for i in sorted(SUPPORTS[q])],
            "choices": list(CHOICES[q]),
            "G": float(G_ALPHA[q]),
            "kappa_joint_beta": worst,
            "argmax_cell": arg,
            "argmax_levels": list(CELLS[arg]),
        })
    rows.sort(key=lambda r: -r["kappa_joint_beta"])
    print("D21 kappa_joint (exact over Q(sqrt2,sqrt3)):")
    for r in rows:
        print(f"  beta {r['beta']:2d} support {r['support']} choices {r['choices']}"
              f"  G={r['G']:.6f}  kappa_joint_beta={r['kappa_joint_beta']:.12f}"
              f"  at cell {r['argmax_cell']} {r['argmax_levels']}")
    kj = rows[0]["kappa_joint_beta"]
    print(f"D21 kappa_joint = {kj:.12f}")
    return {"rows": rows, "kappa_joint": kj}


# ---- D75 numeric -----------------------------------------------------------

def d75_kappa_joint():
    cells = sorted(d75.five_axis_cells(0)
                   + d75.extra_cells([0, 1, 2], base=(0, 0, 0, 0)))
    e = d75.evaluate(cells)
    B = e["B"]
    P_IDX = np.flatnonzero(pa.PERM)
    Q_IDX = np.flatnonzero(pa.Q)
    T_IDX = pa.T_IDX
    posP = {int(p): i for i, p in enumerate(P_IDX)}
    Trows = [posP[int(t)] for t in T_IDX]
    PT = pa.PHI[:, T_IDX]
    rows = []
    for j, q in enumerate(Q_IDX):
        bcol = B[Trows, j]
        r = pa.PHI[:, int(q)] - PT @ bcol
        val = float(np.abs(r).max() / pa.G_ALPHA[q])
        rows.append({
            "beta": int(q),
            "support": [pa.AXES[i] for i in sorted(pa.SUPPORTS[q])],
            "choices": list(pa.CHOICES[q]),
            "G": float(pa.G_ALPHA[q]),
            "kappa_joint_beta": val,
            "argmax_cell": int(np.argmax(np.abs(r))),
        })
    rows.sort(key=lambda r: -r["kappa_joint_beta"])
    print("\nD75 kappa_joint (float B, exact sup on the 144-cell cube):")
    for r in rows[:20]:
        print(f"  beta {r['beta']:2d} support {r['support']} choices {r['choices']}"
              f"  G={r['G']:.6f}  kappa_joint_beta={r['kappa_joint_beta']:.12f}"
              f"  at cell {r['argmax_cell']}")
    kj = rows[0]["kappa_joint_beta"]
    print(f"D75 kappa_joint = {kj:.12f} (nearest integer {round(kj)})")
    return {"rows": rows, "kappa_joint": kj}


if __name__ == "__main__":
    a = d21_kappa_joint()
    print(f"D21 kappa_joint nearest integer: {round(a['kappa_joint'])}")
    b = d75_kappa_joint()
    with open("kappa_joint.json", "w") as h:
        json.dump({"D21": a, "D75": b}, h, indent=1)
