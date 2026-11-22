"""Compact pathwise summary for both arms and the resolving lanes.

For every non-far point and every fixed (route, scheme) path, decompose the
four free axes (partition, pack, form, exp), check the reconstruction, and
report the structural-zero allowance eta_Q, the four-axis high-order allowance
eta_old (order 4), and the permitted interaction signal.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from instance import (N_CELLS, PHI, T_IDX, Q_IDX, G_ALPHA, ORDERS, AXES,
                      SUPPORTS, CHOICES)

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

ORDER4 = np.array([i for i in range(N_CELLS) if ORDERS[i] == 4])


def analyse_arm(path, label):
    data = fc.load_route(path, label)
    R = data["R"]
    xs, ns = data["xs"], data["ns"]
    out = {"label": label, "sha256": data["checks"]["sha256"],
           "lanes": {}}
    for lane_name in ["fp64", "fp32"]:
        lane = fc.LANES.index(lane_name)
        summ = {
            "eta_Q_max": 0.0, "eta_Q_at": None,
            "eta_old_max": 0.0, "eta_old_at": None,
            "q_coef_max": 0.0, "q_coef_at": None,
            "target_coef_max": 0.0, "target_coef_at": None,
            "recon_max": 0.0, "ratio_max": 0.0,
            "n_points_paths": 0,
        }
        for route in (0, 1):
            for scheme in (0, 1):
                base = route * 72 + scheme * 36
                for point in range(R.shape[1]):
                    x = xs[point % len(xs)]
                    if fc.region_of(x) == "far":
                        continue
                    n = ns[point // len(xs)]
                    vec = R[lane, point, base:base + N_CELLS]
                    if np.isnan(vec).any():
                        raise SystemExit("missing cell")
                    c = PHI.T @ vec / N_CELLS
                    recon = PHI @ c
                    summ["n_points_paths"] += 1
                    summ["recon_max"] = max(
                        summ["recon_max"], float(np.abs(recon - vec).max()))
                    eta_Q = float(np.sum(G_ALPHA[Q_IDX] * np.abs(c[Q_IDX])))
                    eta_old = float(np.sum(G_ALPHA[ORDER4] * np.abs(c[ORDER4])))
                    ratio = eta_Q / eta_old if eta_old > 0 else float("inf")
                    summ["ratio_max"] = max(summ["ratio_max"], ratio)
                    qmax_i = Q_IDX[np.argmax(np.abs(c[Q_IDX]) * G_ALPHA[Q_IDX])]
                    tmax_i = T_IDX[np.argmax(np.abs(c[T_IDX]))]
                    at = {"n": int(n), "x": float(x), "region": fc.region_of(x),
                          "route": fc.LEVELS["route"][route],
                          "scheme": fc.LEVELS["scheme"][scheme]}
                    if eta_Q > summ["eta_Q_max"]:
                        summ["eta_Q_max"] = eta_Q
                        summ["eta_Q_at"] = {**at, "component_index": int(qmax_i),
                            "support": [AXES[a] for a in sorted(SUPPORTS[qmax_i])],
                            "choices": list(CHOICES[qmax_i])}
                    if eta_old > summ["eta_old_max"]:
                        summ["eta_old_max"] = eta_old
                        summ["eta_old_at"] = at
                    if abs(c[qmax_i]) > summ["q_coef_max"]:
                        summ["q_coef_max"] = float(abs(c[qmax_i]))
                        summ["q_coef_at"] = {**at, "component_index": int(qmax_i),
                            "support": [AXES[a] for a in sorted(SUPPORTS[qmax_i])],
                            "choices": list(CHOICES[qmax_i])}
                    if abs(c[tmax_i]) > summ["target_coef_max"]:
                        summ["target_coef_max"] = float(abs(c[tmax_i]))
                        summ["target_coef_at"] = {**at, "component_index": int(tmax_i),
                            "support": [AXES[a] for a in sorted(SUPPORTS[tmax_i])],
                            "choices": list(CHOICES[tmax_i])}
        summ["ratio_max"] = (None if summ["ratio_max"] == float("inf")
                             else summ["ratio_max"])
        # order-4 exact zeros?
        out["lanes"][lane_name] = summ
    return out


def main():
    arms = {
        "fused": REPO + "/composition-audit/records/cells/cells-msvc-fused.csv.xz",
        "separate": REPO + "/composition-audit/records/cells/cells-msvc-separate.csv.xz",
    }
    out = {}
    for label, path in arms.items():
        print(f"loading {label} ...", file=sys.stderr)
        out[label] = analyse_arm(path, label)
    with open("pathwise_summary.json", "w") as h:
        json.dump(out, h, indent=1)
    for label, arm in out.items():
        print(f"\n=== {label} ===")
        for lane, s in arm["lanes"].items():
            print(f"{lane}: eta_Q {s['eta_Q_max']:.6e}  eta_old {s['eta_old_max']:.3e}  "
                  f"ratio {s['ratio_max']}  q_coef {s['q_coef_max']:.3e}  "
                  f"target_coef {s['target_coef_max']:.3e}  recon {s['recon_max']:.2e}")
            print(f"   eta_Q at {s['eta_Q_at']}")
            print(f"   q_coef at {s['q_coef_at']}")
            print(f"   target at {s['target_coef_at']}")


if __name__ == "__main__":
    main()
