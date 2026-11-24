"""The D84 frontier ratio and its distribution from the measured record.

The certified frontier at kappa_joint = 4 is the maximum over the
record's 858 points of the pointwise family/full ratio

    (tilde_figure(t) + M_hat(t) + 4 eta_Q(t))
    ------------------------------------------
    (tilde_figure(t) + N_full(t))

with tilde_figure(t) = max_cells figure(cell,t) / F0(x).  The maxima of
the family and full bounds need not occur at the same point.  The
propagated measurement term is reported separately; it is below
2.5e-27 (fp64) and 2.5e-20 (fp32), far below the printed precision.

The script also recomputes the per-term maxima and the one-allowance
ratio that the deposited `d84_frontier_recompute.json` reports, so the
two quantities are not confused.  It writes `frontier_d84.json` beside
itself.

Usage:

    python -I frontier_d84.py
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import d84_fixture as fx             # noqa: E402
import record_figures as rf          # noqa: E402

fc = fx.fc
pa = fx.pa

KAPPA = 4.0
EPS = {"fp64": 1e-29, "fp32": 1e-22}
LANES = ["fp64", "fp32"]
EXPECTED = {
    ("fused", "fp64"): 1.009916,
    ("separate", "fp64"): 1.009916,
    ("fused", "fp32"): 1.635887,
    ("separate", "fp32"): 1.705403,
}


def arm_paths():
    return {
        "fused": os.path.join(fx.REPO, "composition-audit", "records",
                              "cells", "cells-msvc-fused.csv.xz"),
        "separate": os.path.join(fx.REPO, "composition-audit", "records",
                                 "cells", "cells-msvc-separate.csv.xz"),
    }


def percentiles(values):
    v = np.asarray(values, dtype=float)
    return {
        "min": float(v.min()),
        "median": float(np.median(v)),
        "p95": float(np.percentile(v, 95)),
        "max": float(v.max()),
        "frac_le_1.01": float((v <= 1.01).mean()),
        "frac_le_1.05": float((v <= 1.05).mean()),
        "frac_le_1.10": float((v <= 1.10).mean()),
    }


def lane_terms(R_lane, d):
    c = R_lane @ pa.PHI / 144.0
    Nd = len(d["cells"])
    chat = (R_lane[:, d["cells"]] @ d["PP"] / Nd) @ d["inv"].T
    M = chat[:, d["posT"]] @ d["PT"].T
    Mhat = np.abs(M).max(axis=1)
    eta = (d["GQ"][None, :] * np.abs(c[:, d["Q_IDX"]])).sum(axis=1)
    N = c[:, d["T_IDX"]] @ d["PT"].T + c[:, d["Q_IDX"]] @ d["PQ"].T
    Nfull = np.abs(N).max(axis=1)
    return {"Mhat": Mhat, "eta": eta, "Nfull": Nfull, "M": M, "N": N}


def main():
    d = fx.build("d84")
    fx.verify_against_certificate(d)
    measure = {lane: EPS[lane]
               * d["measurement_inflation"]["sum_G_amp_targets"]
               for lane in LANES}
    arms = arm_paths()
    out = {
        "meta": {
            "design": "d84",
            "cells": int(len(d["cells"])),
            "permitted": int(len(d["P_IDX"])),
            "complement": int(len(d["Q_IDX"])),
            "targets": int(len(d["T_IDX"])),
            "kappa_joint": KAPPA,
            "convention": (
                "maximum over the 858 points of "
                "(tilde_figure + M_hat + kappa_joint*eta_Q) / "
                "(tilde_figure + N_full), kappa_joint = 4"),
            "figure": "tilde_figure(t) = max_cells figure(cell,t)/F0(x)",
            "measurement_term": (
                "reported separately; at most 2.5e-27 (fp64) and "
                "2.5e-20 (fp32), below the printed precision"),
            "record": {arm: rf.sha256_file(path)
                       for arm, path in arms.items()},
        },
        "per_arm": {},
        "certificate_check": {},
    }
    print(f"design D84: {len(d['cells'])} cells, {len(d['P_IDX'])} "
          f"permitted, {len(d['Q_IDX'])} complement, "
          f"{len(d['T_IDX'])} targets, kappa_joint {KAPPA}")
    print(f"{'arm':9s} {'lane':5s} {'ratio':>10s} {'median':>10s} "
          f"{'p95':>10s} {'within1%':>9s} {'within5%':>9s} "
          f"{'block':>10s} {'ok':>3s}")
    for arm, path in arms.items():
        run = fc.load_route(path, arm)
        cfig = rf.load_figures(path)
        for lane_name in LANES:
            lane = fc.LANES.index(lane_name)
            scaled = rf.scaled_figure(cfig["FIG"], cfig["FN"], lane_name,
                                      len(run["xs"]))
            tilde = scaled["tilde"]
            terms = lane_terms(run["R"][lane], d)
            full = tilde + terms["Nfull"]
            family = tilde + terms["Mhat"] + KAPPA * terms["eta"]
            ratio = family / full
            ratio_measure = (family + measure[lane_name]) / full
            one = (tilde + terms["Mhat"] + terms["eta"]) / full
            maxima = float(family.max() / full.max())
            p = int(np.argmax(ratio))
            nx = len(run["xs"])
            row = {
                "tilde_min": float(tilde.min()),
                "tilde_max": float(tilde.max()),
                "Nfull_max": float(terms["Nfull"].max()),
                "Mhat_max": float(terms["Mhat"].max()),
                "eta_Q_max": float(terms["eta"].max()),
                "charge_max": float(KAPPA * terms["eta"].max()),
                "full_max": float(full.max()),
                "family_max": float(family.max()),
                "ratio_pointwise": float(ratio.max()),
                "ratio_pointwise_with_measurement": float(
                    ratio_measure.max()),
                "ratio_maxima": maxima,
                "ratio_one_allowance": {
                    "max": float(one.max()),
                    "at": int(np.argmax(one)),
                    "note": (
                        "reported by d84_frontier_recompute.json as "
                        "ratio_pointwise_new; it is not the certified "
                        "frontier and does not use kappa_joint = 4"),
                },
                "percentiles": percentiles(ratio),
                "worst_point": {
                    "point_index": p,
                    "n": int(run["ns"][p // nx]),
                    "x": float(run["xs"][p % nx]),
                    "region": fc.region_of(run["xs"][p % nx]),
                    "tilde_figure": float(tilde[p]),
                    "Nfull": float(terms["Nfull"][p]),
                    "M_hat": float(terms["Mhat"][p]),
                    "eta_Q": float(terms["eta"][p]),
                    "charge": float(KAPPA * terms["eta"][p]),
                    "ratio": float(ratio[p]),
                },
                "signal_maximum": {
                    "Nfull": float(terms["Nfull"].max()),
                    "at": int(np.argmax(terms["Nfull"])),
                },
            }
            out["per_arm"][f"{arm}/{lane_name}"] = row
            block = EXPECTED[(arm, lane_name)]
            ok = round(row["ratio_pointwise"], 6) == block
            out["certificate_check"][f"{arm}/{lane_name}"] = {
                "computed": row["ratio_pointwise"],
                "computed_6dp": round(row["ratio_pointwise"], 6),
                "certificate_block": block,
                "agrees": bool(ok),
            }
            print(f"{arm:9s} {lane_name:5s} "
                  f"{row['ratio_pointwise']:10.6f} "
                  f"{row['percentiles']['median']:10.6f} "
                  f"{row['percentiles']['p95']:10.6f} "
                  f"{100 * row['percentiles']['frac_le_1.01']:8.2f}% "
                  f"{100 * row['percentiles']['frac_le_1.05']:8.2f}% "
                  f"{block:10.6f} {'yes' if ok else 'NO':>3s}")
    out["all_agree"] = all(v["agrees"]
                           for v in out["certificate_check"].values())
    with open(os.path.join(HERE, "frontier_d84.json"), "w") as handle:
        json.dump(out, handle, indent=1, default=float)
    print("\ncertificate frontier block "
          f"{'confirmed' if out['all_agree'] else 'CORRECTED'}")
    return 0 if out["all_agree"] else 1


if __name__ == "__main__":
    sys.exit(main())
