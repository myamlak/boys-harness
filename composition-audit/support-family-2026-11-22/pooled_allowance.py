"""Pooled six-axis allowance from the deposited 576-cell record.

Axes (fullcube order): route(2), scheme(2), partition(3), pack(2), form(3),
exp(2).  Pooled cones from the source audit:

  upward step union: {route, scheme, partition, form, exp}   (5 axes)
  packed sum:        {partition, pack}                       (2 axes)

so the permitted supports are every support without pack, plus {pack} and
{partition,pack}; the structural-zero supports Q contain pack together with
at least one of route, scheme, form, exp.

For every lane/point we decompose the 144-cell response and report
eta_Q, eta_old (order > 5, i.e. order 6), and the permitted target signal.
"""

from __future__ import annotations

import importlib.util
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from instance import axis_basis

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

LEVELS = [2, 2, 3, 2, 3, 2]
AXES = ["route", "scheme", "partition", "pack", "form", "exp"]
W_P = 5
UPWARD = {0, 1, 2, 4, 5}
PACKED = {2, 3}


def strides(counts):
    out, s = [], 1
    for c in reversed(counts):
        out.append(s)
        s *= c
    return list(reversed(out))


STR = strides(LEVELS)
N = math.prod(LEVELS)  # 144
ROWS = [axis_basis(c) for c in LEVELS]
G_ROW = [np.abs(r).max(axis=1) for r in ROWS]


def choices(idx):
    return tuple((idx // s) % c for c, s in zip(LEVELS, STR))


def support(idx):
    return frozenset(i for i, ch in enumerate(choices(idx)) if ch != 0)


def in_pooled(supp):
    return supp <= UPWARD or supp <= PACKED


CHOICES = [choices(i) for i in range(N)]
SUPPORTS = [support(i) for i in range(N)]
ORDERS = [len(s) for s in SUPPORTS]
PERM = np.array([in_pooled(s) for s in SUPPORTS])
Q = ~PERM
TARGET = np.array([PERM[i] and 2 <= ORDERS[i] <= W_P for i in range(N)])
Q_IDX = np.flatnonzero(Q)
T_IDX = np.flatnonzero(TARGET)
ORDER6 = np.array([i for i in range(N) if ORDERS[i] == 6])
G_ALPHA = np.array([
    math.prod(G_ROW[i][ch] for i, ch in enumerate(chs)) for chs in CHOICES
])
print(f"components {N}, permitted {PERM.sum()}, targets {len(T_IDX)}, "
      f"Q {len(Q_IDX)}, main effects {sum(ORDERS[i]==1 for i in range(N))}")


def build_phi():
    phi = np.ones((N, N))
    levels = [tuple((c // s) % cnt for cnt, s in zip(LEVELS, STR))
              for c in range(N)]
    for idx, chs in enumerate(CHOICES):
        v = np.ones(N)
        for axis, ch in enumerate(chs):
            v *= ROWS[axis][ch][np.array([lv[axis] for lv in levels])]
        phi[:, idx] = v
    return phi


PHI = build_phi()
assert np.allclose(PHI.T @ PHI / N, np.eye(N))


def arm_report(path, label):
    data = fc.load_route(path, label)
    R = data["R"]
    out = {"sha256": data["checks"]["sha256"], "lanes": {}}
    for lane_name in ["fp64", "fp32"]:
        lane = fc.LANES.index(lane_name)
        stats = {"eta_Q_max": 0.0, "eta_Q_at": None, "eta_old_max": 0.0,
                 "target_max": 0.0, "target_at": None, "q_coef_max": 0.0,
                 "eta_Q_max_all_points": []}
        for point in range(R.shape[1]):
            vec = R[lane, point, :]
            if np.isnan(vec).any():
                raise SystemExit("missing")
            c = PHI.T @ vec / N
            eta_Q = float(np.sum(G_ALPHA[Q_IDX] * np.abs(c[Q_IDX])))
            eta_old = float(np.sum(G_ALPHA[ORDER6] * np.abs(c[ORDER6])))
            tmax_i = T_IDX[np.argmax(np.abs(c[T_IDX]))]
            qmax_i = Q_IDX[np.argmax(np.abs(c[Q_IDX]) * G_ALPHA[Q_IDX])]
            x = data["xs"][point % len(data["xs"])]
            n = data["ns"][point // len(data["xs"])]
            at = {"n": int(n), "x": float(x), "region": fc.region_of(x)}
            if eta_Q > stats["eta_Q_max"]:
                stats["eta_Q_max"] = eta_Q
                stats["eta_Q_at"] = {**at,
                    "support": [AXES[a] for a in sorted(SUPPORTS[qmax_i])],
                    "component_index": int(qmax_i)}
            stats["eta_old_max"] = max(stats["eta_old_max"], eta_old)
            if abs(c[tmax_i]) > stats["target_max"]:
                stats["target_max"] = float(abs(c[tmax_i]))
                stats["target_at"] = {**at, "component_index": int(tmax_i),
                    "support": [AXES[a] for a in sorted(SUPPORTS[tmax_i])]}
            stats["q_coef_max"] = max(stats["q_coef_max"],
                                      float(abs(c[qmax_i])))
        stats.pop("eta_Q_max_all_points")
        out["lanes"][lane_name] = stats
    return out


def main():
    arms = {
        "fused": REPO + "/composition-audit/records/cells/cells-msvc-fused.csv.xz",
        "separate": REPO + "/composition-audit/records/cells/cells-msvc-separate.csv.xz",
    }
    out = {}
    for label, path in arms.items():
        print("loading", label, file=sys.stderr)
        out[label] = arm_report(path, label)
    with open("pooled_allowance.json", "w") as h:
        json.dump(out, h, indent=1)
    for label, arm in out.items():
        for lane, s in arm["lanes"].items():
            print(f"\n{label} {lane}: eta_Q {s['eta_Q_max']:.4e}  "
                  f"eta_old {s['eta_old_max']:.3e}  target {s['target_max']:.3e}  "
                  f"q_coef {s['q_coef_max']:.3e}")
            print("  eta_Q at", s["eta_Q_at"])
            print("  target at", s["target_at"])


if __name__ == "__main__":
    main()
