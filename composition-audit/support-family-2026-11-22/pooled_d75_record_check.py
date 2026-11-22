"""Run the 75-cell pooled design on the deposited 576-cell record.

For every point of the fused (or separate) fp64/fp32 record:
  * compute the full-cube coefficient vector c over the 144 cells;
  * apply the D75 design matrix: raw = (1/75) sum_D y Phi_P, chat = inv(Omega_PP) raw;
  * compare chat with c_P and with the alias identity chat - c_P = B c_Q;
  * report the certified bound chat + kappa * eta_Q(point).
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pooled_allowance as pa
import pooled_d75 as d75

def _harness_root():
    """Locate the harness repository holding the full-cube record."""
    override = os.environ.get("BOYS_HARNESS_REPO")
    if override:
        return os.path.abspath(override)
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    probe = os.path.join(root, "composition-audit", "fullcube",
                         "fullcube_decomposition.py")
    if os.path.isfile(probe):
        return root
    raise SystemExit("set BOYS_HARNESS_REPO to the harness repository root")


REPO = _harness_root()
SCRIPT = os.path.join(REPO, "composition-audit", "fullcube",
                      "fullcube_decomposition.py")
spec = importlib.util.spec_from_file_location("fc", SCRIPT)
fc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fc)

P_IDX = np.flatnonzero(pa.PERM)
Q_IDX = np.flatnonzero(pa.Q)
T_IDX = pa.T_IDX
G = pa.G_ALPHA
PHI = pa.PHI
N = pa.N
CELLS = d75.five_axis_cells(0) + d75.extra_cells([0, 1, 2], base=(0, 0, 0, 0))
CELLS = np.array(sorted(CELLS))
ND = len(CELLS)
F = PHI[CELLS, :]
PP = F[:, P_IDX]
QQ = F[:, Q_IDX]
OmPP = PP.T @ PP / ND
OmPQ = PP.T @ QQ / ND
invPP = np.linalg.inv(OmPP)
B = invPP @ OmPQ
kappa_all = np.abs(B / G[Q_IDX][None, :])
pos = {int(p): i for i, p in enumerate(P_IDX)}
kappa_t = kappa_all[[pos[int(t)] for t in T_IDX]].max(axis=1)


def run(path, label):
    data = fc.load_route(path, label)
    R = data["R"]
    out = {}
    for lane_name in ["fp64", "fp32"]:
        lane = fc.LANES.index(lane_name)
        worst_alias = 0.0
        worst_identity = 0.0
        worst_target_est = 0.0
        ratio_max = 0.0
        n_above = 0
        examples = []
        for point in range(R.shape[1]):
            vec = R[lane, point, :]
            if np.isnan(vec).any():
                raise SystemExit("missing")
            c = PHI.T @ vec / N
            raw = PP.T @ vec[CELLS] / ND
            chat = invPP @ raw
            # alias identity on targets
            pred = c[P_IDX] + B @ c[Q_IDX]
            worst_identity = max(worst_identity, float(np.abs(chat - pred).max()))
            # target comparison
            d = np.abs(chat - c[P_IDX])
            worst_alias = max(worst_alias, float(d.max()))
            etaQ = float(np.sum(G[Q_IDX] * np.abs(c[Q_IDX])))
            bound = kappa_t * etaQ
            dt = d[[pos[int(t)] for t in T_IDX]]
            bt = bound
            if np.any(bt > 1e-30):
                ratio_max = max(ratio_max, float((dt[bt > 1e-30] / bt[bt > 1e-30]).max()))
            est_t = np.abs(chat[[pos[int(t)] for t in T_IDX]]).max()
            worst_target_est = max(worst_target_est, float(est_t))
            if bt.size and est_t > bt.max():
                n_above += 1
            if len(examples) < 3:
                i = int(np.argmax(np.abs(chat[[pos[int(t)] for t in T_IDX]])))
                x = data["xs"][point % len(data["xs"])]
                n = data["ns"][point // len(data["xs"])]
                examples.append({
                    "n": int(n), "x": float(x),
                    "target_est_max": float(est_t),
                    "eta_Q": etaQ, "bound_max": float(bound.max()),
                })
        out[lane_name] = {
            "max_chat_minus_c": worst_alias,
            "max_alias_identity_residual": worst_identity,
            "max_target_est": worst_target_est,
            "max_chat_minus_c_over_bound": ratio_max,
            "points_est_above_bound": n_above,
            "examples": examples,
        }
        print(f"{label} {lane_name}: max|chat-c| {worst_alias:.3e}, "
              f"identity residual {worst_identity:.2e}, "
              f"max target est {worst_target_est:.3e}, "
              f"max ratio to kappa*eta_Q {ratio_max:.3f}, "
              f"points above bound {n_above}/{R.shape[1]}")
        print("   examples:", examples)
    return out


def main():
    arms = {
        "fused": REPO + "/composition-audit/records/cells/cells-msvc-fused.csv.xz",
        "separate": REPO + "/composition-audit/records/cells/cells-msvc-separate.csv.xz",
    }
    out = {}
    for label, path in arms.items():
        print("loading", label, file=sys.stderr)
        out[label] = run(path, label)
    with open("pooled_d75_record_check.json", "w") as h:
        json.dump(out, h, indent=1)


if __name__ == "__main__":
    main()
