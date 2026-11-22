"""Corrected frontier with the joint contamination constant.

For every point and design, compute

  M_hat   = max_a |sum_{alpha in T} chat_alpha Phi_alpha(a)|   (measured)
  eta_Q   = sum_{beta in Q} G_beta |c_beta|                    (measured full cube)
  N_full  = max_a |sum_{alpha in T} c_alpha Phi_alpha(a)
                    + sum_{beta in Q} c_beta Phi_beta(a)|      (full cube)

The corrected non-additive bound is M_hat + kappa_joint * eta_Q
(+ measurement, negligible here).  The full-cube bound is figure + N_full.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact_verify as ev
import pooled_allowance as pa
import pooled_d75
from instance import (PHI as PHI36, T_IDX as T36, Q_IDX as Q36,
                      G_ALPHA as G36, N_CELLS as N36)

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

D21 = [0, 1, 2, 4, 7, 9, 11, 12, 13, 14, 16, 19, 21, 23, 24, 25, 26, 28,
       31, 33, 35]
D75 = sorted(pooled_d75.five_axis_cells(0)
             + pooled_d75.extra_cells([0, 1, 2], base=(0, 0, 0, 0)))

P36 = ev.P_IDX
posP36 = {int(p): i for i, p in enumerate(P36)}
posT36 = [posP36[int(t)] for t in T36]

P144 = np.flatnonzero(pa.PERM)
posP144 = {int(p): i for i, p in enumerate(P144)}
posT144 = [posP144[int(t)] for t in pa.T_IDX]

# D21 matrices
PP21 = PHI36[np.ix_(D21, P36)]
Om21 = PP21.T @ PP21 / len(D21)
inv21 = np.linalg.inv(Om21)
PT21 = PHI36[:, T36]
PQ21 = PHI36[:, Q36]

# D75 matrices
PP75 = pa.PHI[np.ix_(D75, P144)]
Om75 = PP75.T @ PP75 / len(D75)
inv75 = np.linalg.inv(Om75)
PT144 = pa.PHI[:, pa.T_IDX]
PQ144 = pa.PHI[:, pa.Q_IDX]

KJ21 = 4.0
KJ75 = 3.0


def measurement_inflation(design_cells, P_IDX, T_IDX, PHI, inv, G):
    """max_alpha sum_gamma |inv[alpha,gamma]| * (1/N) sum_a |Phi_gamma(a)|."""
    cells = np.asarray(design_cells)
    avg_abs = np.abs(PHI[cells, :]).mean(axis=0)
    amp = np.abs(inv) @ avg_abs[P_IDX]
    pos = {int(p): i for i, p in enumerate(P_IDX)}
    return {
        "max_all_permitted": float(amp.max()),
        "max_targets": float(amp[[pos[int(t)] for t in T_IDX]].max()),
        "sum_G_amp_targets": float(
            (G[T_IDX] * amp[[pos[int(t)] for t in T_IDX]]).sum()),
    }


def run_arm(path, label):
    data = fc.load_route(path, label)
    R = data["R"]
    out = {}
    for lane_name in ["fp64", "fp32"]:
        lane = fc.LANES.index(lane_name)
        acc = {
            "path": {"M_hat": 0.0, "etaQ": 0.0, "Nfull": 0.0,
                     "post": 0.0, "at_post": None},
            "pooled": {"M_hat": 0.0, "etaQ": 0.0, "Nfull": 0.0,
                       "post": 0.0, "at_post": None},
        }
        for point in range(R.shape[1]):
            vec = R[lane, point, :]
            if np.isnan(vec).any():
                raise SystemExit("missing")
            x = data["xs"][point % len(data["xs"])]
            n = data["ns"][point // len(data["xs"])]
            at = (int(n), float(x), fc.region_of(x))

            # pooled D75
            c = pa.PHI.T @ vec / 144.0
            chat = inv75 @ (PP75.T @ vec[D75] / 75.0)
            mhat = float(np.abs(PT144 @ chat[posT144]).max())
            etaq = float(np.sum(pa.G_ALPHA[pa.Q_IDX] * np.abs(c[pa.Q_IDX])))
            nfull = float(np.abs(PT144 @ c[pa.T_IDX]
                                 + PQ144 @ c[pa.Q_IDX]).max())
            post = mhat + KJ75 * etaq
            P = acc["pooled"]
            P["M_hat"] = max(P["M_hat"], mhat)
            P["etaQ"] = max(P["etaQ"], etaq)
            P["Nfull"] = max(P["Nfull"], nfull)
            if post > P["post"]:
                P["post"] = post
                P["at_post"] = {"n": at[0], "x": at[1], "region": at[2],
                                "M_hat": mhat, "etaQ": etaq, "Nfull": nfull}

            # pathwise D21, four fixed (route, scheme) paths
            for base in (0, 36, 72, 108):
                v36 = vec[base:base + N36]
                c36 = PHI36.T @ v36 / N36
                chat36 = inv21 @ (PP21.T @ v36[D21] / len(D21))
                mhat = float(np.abs(PT21 @ chat36[posT36]).max())
                etaq = float(np.sum(G36[Q36] * np.abs(c36[Q36])))
                nfull = float(np.abs(PT21 @ c36[T36]
                                     + PQ21 @ c36[Q36]).max())
                post = mhat + KJ21 * etaq
                A = acc["path"]
                A["M_hat"] = max(A["M_hat"], mhat)
                A["etaQ"] = max(A["etaQ"], etaq)
                A["Nfull"] = max(A["Nfull"], nfull)
                if post > A["post"]:
                    A["post"] = post
                    A["at_post"] = {"n": at[0], "x": at[1], "region": at[2],
                                    "path_base": base, "M_hat": mhat,
                                    "etaQ": etaq, "Nfull": nfull}
        out[lane_name] = acc
    return out


def main():
    fig = {"fp64": 2.558314305293831e-14, "fp32": 1.5e-7}
    old_sumGk = {"path": 14.0, "pooled": 66.0}
    kj = {"path": KJ21, "pooled": KJ75}
    arms = {
        "fused": REPO + "/composition-audit/records/cells/cells-msvc-fused.csv.xz",
        "separate": REPO + "/composition-audit/records/cells/cells-msvc-separate.csv.xz",
    }
    out = {}
    for label, path in arms.items():
        print("loading", label, file=sys.stderr)
        out[label] = run_arm(path, label)

    # measurement inflation (formula of the corrected theorem)
    mi21 = measurement_inflation(D21, P36, T36, PHI36, inv21, G36)
    mi75 = measurement_inflation(D75, P144, pa.T_IDX, pa.PHI, inv75,
                                 pa.G_ALPHA)
    print("\nmeasurement inflation factors (max_alpha sum_gamma |inv| * mean|Phi|):")
    print("  D21:", mi21)
    print("  D75:", mi75)

    print("\ncorrected frontier, per-term breakdown")
    print("design/lane  full_nonadd  full_bound  eta_Q_max  sumGk(old)  "
          "kappa_joint  charge  M_hat  posterior  family_bound  ratio")
    table = {}
    for label, arm in out.items():
        for lane in ["fp64", "fp32"]:
            for kind in ["path", "pooled"]:
                a = arm[lane][kind]
                full_nonadd = a["Nfull"]
                full_bound = fig[lane] + full_nonadd
                eta = a["etaQ"]
                charge = kj[kind] * eta
                posterior = a["post"]
                family_bound = fig[lane] + posterior
                ratio = family_bound / full_bound
                key = f"{label}/{lane}/{kind}"
                table[key] = {
                    "full_nonadd": full_nonadd,
                    "full_bound": full_bound,
                    "eta_Q_max": eta,
                    "sumGk_old": old_sumGk[kind],
                    "old_charge": (1.0 + old_sumGk[kind]) * eta,
                    "kappa_joint": kj[kind],
                    "charge": charge,
                    "M_hat_max": a["M_hat"],
                    "posterior": posterior,
                    "family_bound": family_bound,
                    "ratio": ratio,
                    "at_post": a["at_post"],
                }
                print(f"{key:28s} {full_nonadd:.4e} {full_bound:.4e} "
                      f"{eta:.4e} {old_sumGk[kind]:.0f} {kj[kind]:.0f} "
                      f"{charge:.4e} {a['M_hat']:.4e} {posterior:.4e} "
                      f"{family_bound:.4e} {ratio:.4f}")
    with open("corrected_frontier.json", "w") as h:
        json.dump({"measurement_inflation": {"D21": mi21, "D75": mi75},
                   "table": table}, h, indent=1)


if __name__ == "__main__":
    main()
