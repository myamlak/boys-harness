"""The three prospective checks on the D84 pooled design.

A. Split-sample transport of the measured allowance: 100 random half
   splits of the 858 points; a scalar allowance is estimated on one half
   and the certified inequality is checked on the other.
B. A planted outside-family order-two channel supported on pack and
   form: the null allowance is re-measured, the certificate ratio is
   reported at delta0 and the frozen null allowance's failure threshold
   is bracketed per precision.
C. The distribution of the pointwise family/full ratio over the 858
   points: median, 95th percentile, maximum and the fractions within one
   and five per cent.

The protocol is the one in the article's reproduction package
(`worked-check/prospective/additions.py`, SHA-256 546e0304...), run with
the corrected 84-run pooled design instead of the earlier 75-run one.
The D21 path rows are unchanged; only the pooled rows differ.  With
`--design d75` the same script reproduces the package's D75 numbers as a
self-check.

The pool fixture, the six-axis basis and the measured record are read
from the harness repository named by BOYS_HARNESS_REPO (or the one two
directories above this script).  No harness or article file is written.

Usage, from this directory:

    python -I prospective_d84.py
    python -I prospective_d84.py --design d75 --out-dir repro_d75
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import d84_fixture as fx             # noqa: E402
import record_figures as rf          # noqa: E402

cf = fx.cf
fc = fx.fc
pa = fx.pa
d75 = fx.d75

SEED = 20261123
N_SPLITS = 100
N_CELL_SPLITS = 100
EPS = {"fp64": 1e-29, "fp32": 1e-22}
PLANT_IDX = 8
CHUNK = 25
LANES = ["fp64", "fp32"]


def arm_paths():
    return {
        "fused": os.path.join(fx.REPO, "composition-audit", "records",
                              "cells", "cells-msvc-fused.csv.xz"),
        "separate": os.path.join(fx.REPO, "composition-audit", "records",
                                 "cells", "cells-msvc-separate.csv.xz"),
    }


def measure_factors(d):
    mi21 = cf.measurement_inflation(cf.D21, cf.P36, cf.T36, cf.PHI36,
                                    cf.inv21, cf.G36)
    mi_pooled = d["measurement_inflation"]
    return {
        "path": {lane: EPS[lane] * mi21["sum_G_amp_targets"]
                 for lane in LANES},
        "pooled": {lane: EPS[lane] * mi_pooled["sum_G_amp_targets"]
                   for lane in LANES},
        "mi21": mi21,
        "mi_pooled": mi_pooled,
    }


def precompute(R_lane, d):
    path = []
    for base in (0, 36, 72, 108):
        v = R_lane[:, base:base + 36]
        c = v @ cf.PHI36 / 36.0
        chat = (v[:, cf.D21] @ cf.PP21 / len(cf.D21)) @ cf.inv21.T
        S = (cf.G36[cf.T36][None, :]
             * np.abs(chat[:, cf.posT36])).sum(axis=1)
        eta = (cf.G36[cf.Q36][None, :]
               * np.abs(c[:, cf.Q36])).sum(axis=1)
        N = c[:, cf.T36] @ cf.PT21.T + c[:, cf.Q36] @ cf.PQ21.T
        M = chat[:, cf.posT36] @ cf.PT21.T
        path.append({"c": c, "chat": chat, "S": S, "eta": eta,
                     "N": N, "M": M, "Nfull": np.abs(N).max(axis=1)})
    c = R_lane @ pa.PHI / 144.0
    Nd = len(d["cells"])
    chat = (R_lane[:, d["cells"]] @ d["PP"] / Nd) @ d["inv"].T
    S = (d["GT"][None, :] * np.abs(chat[:, d["posT"]])).sum(axis=1)
    eta = (d["GQ"][None, :] * np.abs(c[:, d["Q_IDX"]])).sum(axis=1)
    N = c[:, d["T_IDX"]] @ d["PT"].T + c[:, d["Q_IDX"]] @ d["PQ"].T
    M = chat[:, d["posT"]] @ d["PT"].T
    pooled = {"c": c, "chat": chat, "S": S, "eta": eta, "N": N, "M": M,
              "Nfull": np.abs(N).max(axis=1)}
    return path, pooled


def region_array(xs, ns):
    nx = len(xs)
    return np.array([fc.region_of(xs[p % nx])
                     for p in range(len(ns) * nx)])


def stats(values):
    v = np.asarray(values, dtype=float)
    return {"min": float(v.min()), "median": float(np.median(v)),
            "p95": float(np.percentile(v, 95)), "max": float(v.max()),
            "mean": float(v.mean()), "std": float(v.std()),
            "n": int(v.size)}


# --------------------------------------------------------------------------
# A. Prospective transport
# --------------------------------------------------------------------------


def transport_test(path, pooled, d, lane_name, meas, xs, ns, splits,
                   regions):
    mp = meas["path"][lane_name]
    md = meas["pooled"][lane_name]
    eta_path = np.max([b["eta"] for b in path], axis=0)
    eta_pooled = pooled["eta"]
    n = len(eta_path)
    n_splits = len(splits)
    out = {}
    for kind in ("path", "pooled"):
        ratios = np.empty((n_splits, n // 2))
        ratios_ct = np.empty((n_splits, n // 2))
        split_cov = np.zeros(n_splits)
        split_worst = np.zeros(n_splits)
        split_cov_ct = np.zeros(n_splits)
        worst = {"ratio": -1.0}
        worst_ct = {"ratio": -1.0}
        cell_ok = cell_total = 0
        ct_point_ok = 0
        for si, perm in enumerate(splits):
            train, test = perm[:n // 2], perm[n // 2:]
            if kind == "path":
                eta_train = float(eta_path[train].max())
                bound = np.array([b["S"][test] + cf.KJ21 * eta_train + mp
                                  for b in path])
                bound_cell = np.array([np.abs(b["M"][test])
                                       + cf.KJ21 * eta_train + mp
                                       for b in path])
                true = np.array([b["Nfull"][test] for b in path])
                true_cell = np.array([np.abs(b["N"][test]) for b in path])
                r = (true / bound).max(axis=0)
                rct = (true_cell / bound_cell).max(axis=(0, 2))
                cell_ok += int((true_cell <= bound_cell).sum())
                cell_total += int(true_cell.size)
            else:
                eta_train = float(eta_pooled[train].max())
                bound = pooled["S"][test] + d["kappa_joint"] * eta_train + md
                bound_cell = np.abs(pooled["M"][test]) \
                    + d["kappa_joint"] * eta_train + md
                true = pooled["Nfull"][test]
                true_cell = np.abs(pooled["N"][test])
                r = true / bound
                rct = (true_cell / bound_cell).max(axis=1)
                cell_ok += int((true_cell <= bound_cell).sum())
                cell_total += int(true_cell.size)
            ratios[si] = r
            ratios_ct[si] = rct
            ct_point_ok += int((rct <= 1.0).sum())
            split_cov[si] = float((r <= 1.0).mean())
            split_worst[si] = float(r.max())
            split_cov_ct[si] = float((rct <= 1.0).mean())
            if r.max() > worst["ratio"]:
                j = int(np.argmax(r))
                p = int(test[j])
                worst = {"ratio": float(r[j]), "split": si,
                         "point_index": p,
                         "n": int(ns[p // len(xs)]),
                         "x": float(xs[p % len(xs)]),
                         "region": str(regions[p])}
            if rct.max() > worst_ct["ratio"]:
                j = int(np.argmax(rct))
                p = int(test[j])
                worst_ct = {"ratio": float(rct[j]), "split": si,
                            "point_index": p,
                            "n": int(ns[p // len(xs)]),
                            "x": float(xs[p % len(xs)]),
                            "region": str(regions[p])}
        flat = ratios.ravel()
        flat_ct = ratios_ct.ravel()
        test_points = np.array([p for perm in splits
                                for p in perm[n // 2:]])
        test_regions = regions[test_points]
        per_region = {}
        for region in ("A", "B", "far"):
            mask = test_regions == region
            if not mask.any():
                per_region[region] = None
                continue
            per_region[region] = {
                "n_trials": int(mask.sum()),
                "coverage": float((flat[mask] <= 1.0).mean()),
                "worst_ratio": float(flat[mask].max()),
                "median_ratio": float(np.median(flat[mask])),
                "p95_ratio": float(np.percentile(flat[mask], 95)),
            }
        out[kind] = {
            "allowance": (
                "scalar, max over training points of the row's pointwise "
                "eta_Q = max_b eta_b (path) / eta_pooled; kappa applied "
                "separately"),
            "inequality": (
                "max_a |N(a)| <= S(t) + kappa_joint*eta_train + "
                "measurement, S(t)=sum_T G_a|chat_a(t)|, per design block; "
                "cell-tight variant replaces S(t) by |sum_T chat_a Phi_a(a)| "
                "at each cell"),
            "coverage": float((flat <= 1.0).mean()),
            "cell_tight_point_coverage":
                ct_point_ok / (n_splits * (n // 2)),
            "cell_tight_cell_coverage": cell_ok / cell_total,
            "cellwise_coverage": cell_ok / cell_total,
            "worst_ratio": worst,
            "ratio_stats": stats(flat),
            "ratio_ct_stats": stats(flat_ct),
            "split_coverage_min": float(split_cov.min()),
            "split_coverage_max": float(split_cov.max()),
            "split_coverage_mean": float(split_cov.mean()),
            "split_coverage_std": float(split_cov.std()),
            "split_worst_ratio_min": float(split_worst.min()),
            "split_worst_ratio_max": float(split_worst.max()),
            "split_coverage_ct_min": float(split_cov_ct.min()),
            "split_coverage_ct_max": float(split_cov_ct.max()),
            "cell_tight_worst_ratio": worst_ct,
            "per_region": per_region,
        }
    return out


def cell_transport_test(R_lane, path, pooled, d, lane_name, meas, xs, ns,
                        masks36, masks144):
    mp = meas["path"][lane_name]
    md = meas["pooled"][lane_name]
    out = {}
    for kind in ("path", "pooled"):
        point_ratios = []
        cell_ok = cell_total = 0
        worst = {"ratio": -1.0}
        if kind == "path":
            for mi, mask in enumerate(masks36):
                use = np.zeros(36, dtype=bool)
                use[mask] = True
                held = ~use
                block_ratio = []
                for b in path:
                    cq = (b["c"][:, use]
                          @ cf.PHI36[np.ix_(np.flatnonzero(use),
                                            cf.Q36)]) / use.sum()
                    eta_s = (cf.G36[cf.Q36][None, :]
                             * np.abs(cq)).sum(axis=1)
                    bound = b["S"] + cf.KJ21 * eta_s + mp
                    r = np.abs(b["N"][:, held]) / bound[:, None]
                    block_ratio.append(r.max(axis=1))
                    cell_ok += int((np.abs(b["N"][:, held])
                                    <= bound[:, None]).sum())
                    cell_total += int(r.size)
                r = np.max(block_ratio, axis=0)
                point_ratios.append(r)
                if r.max() > worst["ratio"]:
                    j = int(np.argmax(r))
                    worst = {"ratio": float(r[j]), "mask": mi,
                             "point_index": j,
                             "n": int(ns[j // len(xs)]),
                             "x": float(xs[j % len(xs)]),
                             "region": fc.region_of(xs[j % len(xs)])}
        else:
            for mi, mask in enumerate(masks144):
                use = np.zeros(144, dtype=bool)
                use[mask] = True
                held = ~use
                cq = (R_lane[:, use]
                      @ pa.PHI[np.ix_(np.flatnonzero(use),
                                      d["Q_IDX"])]) / use.sum()
                eta_s = (d["GQ"][None, :] * np.abs(cq)).sum(axis=1)
                bound = pooled["S"] + d["kappa_joint"] * eta_s + md
                r = np.abs(pooled["N"][:, held]) / bound[:, None]
                cell_ok += int((np.abs(pooled["N"][:, held])
                                <= bound[:, None]).sum())
                cell_total += int(r.size)
                rr = r.max(axis=1)
                point_ratios.append(rr)
                if rr.max() > worst["ratio"]:
                    j = int(np.argmax(rr))
                    worst = {"ratio": float(rr[j]), "mask": mi,
                             "point_index": j,
                             "n": int(ns[j // len(xs)]),
                             "x": float(xs[j % len(xs)]),
                             "region": fc.region_of(xs[j % len(xs)])}
        flat = np.concatenate(point_ratios)
        out[kind] = {
            "cell_splits": N_CELL_SPLITS,
            "inequality": (
                "|N(a)| <= S(t) + kappa_joint*eta_train_cells + measurement "
                "for every held-out cell a; eta_train_cells = "
                "sum_Q G_b|c_b^S| with c_b^S the 1/|S| subset projection"),
            "coverage": float((flat <= 1.0).mean()),
            "cellwise_coverage": cell_ok / cell_total,
            "worst_ratio": worst,
            "ratio_stats": stats(flat),
        }
    return out


# --------------------------------------------------------------------------
# B. Planted outside-family coefficient
# --------------------------------------------------------------------------


def alias_table(b, idx_map, supports, axis_names, G, orders, choices):
    rows = []
    for pos, value in enumerate(b):
        comp = int(idx_map[pos])
        if abs(value) > 1e-12:
            rows.append({
                "component_index": comp,
                "support": [axis_names[i] for i in sorted(supports[comp])],
                "channel_choices": [int(c) for c in choices[comp]],
                "order": int(orders[comp]),
                "G": float(G[comp]),
                "alias_coefficient": float(value),
                "G_times_alias": float(G[comp] * abs(value)),
            })
    rows.sort(key=lambda r: -abs(r["alias_coefficient"]))
    return rows


def planted_test(path, pooled, d, lane_name, meas, xs, ns, arm):
    G21 = cf.G36[cf.T36]
    G_idx = float(pa.G_ALPHA[PLANT_IDX])
    phi21 = cf.PHI36[:, PLANT_IDX]
    phi_pooled = pa.PHI[:, PLANT_IDX]
    b21 = cf.inv21 @ (cf.PP21.T @ phi21[cf.D21] / len(cf.D21))
    b_pooled = d["inv"] @ (d["PP"].T @ phi_pooled[d["cells"]]
                           / len(d["cells"]))
    out = {}

    def run_design(kind):
        if kind == "path":
            blocks = path
            eta_b = [b["eta"] for b in blocks]
            eta = np.max(eta_b, axis=0)
            kappa = cf.KJ21
            bvec, posT, G_T, phi = b21, cf.posT36, G21, phi21
            chat_T = [b["chat"][:, cf.posT36] for b in blocks]
            N_base = [b["N"] for b in blocks]
            S_base = [b["S"] for b in blocks]
            c8 = [b["c"][:, PLANT_IDX] for b in blocks]
            measure = meas["path"][lane_name]
            axis_names, Gvec = fx.inst.AXES, fx.inst.G_ALPHA
            orders, choices = fx.inst.ORDERS, fx.inst.CHOICES
            supports = fx.inst.SUPPORTS
            idx_map = cf.P36
            PT = cf.PT21
        else:
            blocks = [pooled]
            eta_b = [pooled["eta"]]
            eta = pooled["eta"]
            kappa = d["kappa_joint"]
            bvec, posT, G_T, phi = b_pooled, d["posT"], d["GT"], phi_pooled
            chat_T = [pooled["chat"][:, d["posT"]]]
            N_base = [pooled["N"]]
            S_base = [pooled["S"]]
            c8 = [pooled["c"][:, PLANT_IDX]]
            measure = meas["pooled"][lane_name]
            axis_names, Gvec = pa.AXES, pa.G_ALPHA
            orders, choices = pa.ORDERS, pa.CHOICES
            supports = pa.SUPPORTS
            idx_map = d["P_IDX"]
            PT = d["PT"]
        bT = bvec[posT]
        comp_of_pos = np.asarray(idx_map)
        tmask = np.isin(comp_of_pos, d["T_IDX"] if kind == "pooled"
                        else cf.T36)
        r0mask = np.array([orders[int(c)] <= 1 for c in comp_of_pos])
        r0_alias = bvec[r0mask]
        r0_G = np.array([Gvec[int(c)] for c in comp_of_pos[r0mask]])
        r0_components = [
            {"component_index": int(c),
             "support": [axis_names[i] for i in sorted(supports[int(c)])],
             "order": int(orders[int(c)]),
             "alias_coefficient": float(v)}
            for c, v in zip(comp_of_pos[r0mask], r0_alias)
            if abs(v) > 1e-12]
        t_components = [
            {"component_index": int(c),
             "support": [axis_names[i] for i in sorted(supports[int(c)])],
             "order": int(orders[int(c)]),
             "alias_coefficient": float(v)}
            for c, v in zip(comp_of_pos[tmask], bvec[tmask])
            if abs(v) > 1e-12]
        d0 = kappa * float(eta.max())
        max_cT = max(float(np.abs(b["c"][:, d["T_IDX"]
                                        if kind == "pooled" else cf.T36]
                                  ).max()) for b in blocks)
        max_S = max(float(S.max()) for S in S_base)
        max_M = max(float(np.abs(b["M"]).max()) for b in blocks)
        support = sorted(supports[PLANT_IDX])
        residual = (phi - PT @ bT)
        kappa_beta = float(np.abs(residual).max() / G_idx)

        def one_delta(delta, variant):
            best = {"ratio": -1.0}
            nfail = 0
            max_true = 0.0
            for bi in range(len(blocks)):
                tr = np.abs(N_base[bi] + delta * phi).max(axis=1)
                max_true = max(max_true, float(tr.max()))
                sp = (G_T[None, :]
                      * np.abs(chat_T[bi] + delta * bT)).sum(axis=1)
                if variant == "global_frozen":
                    eta_used = np.full_like(sp, float(eta.max()))
                elif variant == "pointwise_frozen":
                    eta_used = np.asarray(eta_b[bi])
                elif variant == "self_consistent_exact":
                    eta_used = eta_b[bi] + G_idx * (
                        np.abs(c8[bi] + delta) - np.abs(c8[bi]))
                else:
                    eta_used = eta_b[bi] + G_idx * abs(delta)
                bound = sp + kappa * eta_used + measure
                r = tr / bound
                nfail += int((r > 1.0).sum())
                j = int(np.argmax(r))
                if r[j] > best["ratio"]:
                    best = {"ratio": float(r[j]), "true": float(tr[j]),
                            "bound": float(bound[j]), "point_index": j,
                            "block": bi, "sp": float(sp[j]),
                            "eta_used": float(eta_used[j]),
                            "n": int(ns[j // len(xs)]),
                            "x": float(xs[j % len(xs)]),
                            "region": fc.region_of(xs[j % len(xs)])}
            best["failing_points"] = nfail
            best["max_true_record"] = max_true
            return best

        def sweep(deltas, variant):
            nd = len(deltas)
            ratio = np.zeros(nd)
            for s in range(0, nd, CHUNK):
                dv = deltas[s:s + CHUNK]
                for bi in range(len(blocks)):
                    tr = np.abs(N_base[bi][None, :, :]
                                + dv[:, None, None]
                                * phi[None, None, :]).max(axis=2)
                    sp = (G_T[None, None, :]
                          * np.abs(chat_T[bi][None, :, :]
                                   + dv[:, None, None]
                                   * bT[None, None, :])).sum(axis=2)
                    if variant == "global_frozen":
                        eta_used = np.full_like(sp, float(eta.max()))
                    elif variant == "pointwise_frozen":
                        eta_used = np.broadcast_to(
                            eta_b[bi][None, :], sp.shape)
                    elif variant == "self_consistent_exact":
                        eta_used = eta_b[bi][None, :] + G_idx * (
                            np.abs(c8[bi][None, :] + dv[:, None])
                            - np.abs(c8[bi][None, :]))
                    else:
                        eta_used = eta_b[bi][None, :] \
                            + G_idx * np.abs(dv)[:, None]
                    bound = sp + kappa * eta_used + measure
                    rr = (tr / bound).max(axis=1)
                    ratio[s:s + CHUNK] = np.maximum(
                        ratio[s:s + CHUNK], rr)
            return ratio

        variants = ("global_frozen", "pointwise_frozen",
                    "self_consistent_exact", "self_consistent_upper")
        results = {}
        for variant in variants:
            for sign in (+1, -1):
                absd = np.concatenate([[0.0],
                                       np.geomspace(1e-3 * d0, 1e3 * d0,
                                                    400)])
                deltas = sign * absd
                ratio = sweep(deltas, variant)
                fail = ratio > 1.0
                threshold = None
                first_fail = last_hold = None
                if fail.any():
                    i0 = int(np.argmax(fail))
                    f0 = one_delta(float(deltas[i0]), variant)
                    first_fail = {"delta": float(deltas[i0]),
                                  "delta_over_delta0": float(absd[i0] / d0),
                                  "ratio": float(ratio[i0]), **{
                                      k: f0[k] for k in
                                      ("true", "bound", "failing_points",
                                       "max_true_record", "point_index",
                                       "block", "n", "x", "region",
                                       "sp", "eta_used")}}
                    if i0 > 0:
                        h0 = one_delta(float(deltas[i0 - 1]), variant)
                        last_hold = {"delta": float(deltas[i0 - 1]),
                                     "delta_over_delta0":
                                         float(absd[i0 - 1] / d0),
                                     "ratio": float(ratio[i0 - 1]), **{
                                         k: h0[k] for k in
                                         ("true", "bound", "failing_points",
                                          "max_true_record", "point_index",
                                          "block", "n", "x", "region",
                                          "sp", "eta_used")}}
                        lo, hi = absd[i0 - 1], absd[i0]
                        for _ in range(60):
                            mid = 0.5 * (lo + hi)
                            rmid = float(sweep(np.array([sign * mid]),
                                               variant)[0])
                            if rmid > 1.0:
                                hi = mid
                            else:
                                lo = mid
                        threshold = {"abs_delta": float(0.5 * (lo + hi)),
                                     "delta_over_delta0":
                                         float(0.5 * (lo + hi) / d0)}
                    last = int(np.max(np.flatnonzero(fail)))
                    recovered = last < len(deltas) - 1
                else:
                    recovered = None
                table = []
                for frac in (0.0, 0.1, 1.0, 3.0, 10.0):
                    dv = sign * frac * d0
                    hit = one_delta(dv, variant)
                    table.append({
                        "delta": float(dv),
                        "delta_over_delta0": float(frac),
                        "max_true": hit["max_true_record"],
                        "true_at_binding": hit["true"],
                        "bound_at_binding": hit["bound"],
                        "S_at_binding": hit["sp"],
                        "eta_used": hit["eta_used"],
                        "kappa_eta_used": kappa * hit["eta_used"],
                        "ratio": hit["ratio"],
                        "failing_points": hit["failing_points"],
                        "pass": bool(hit["ratio"] <= 1.0),
                        "n": hit["n"], "x": hit["x"],
                        "region": hit["region"]})
                for tag, item in (("threshold_prev", last_hold),
                                  ("threshold_fail", first_fail)):
                    if item:
                        table.append({
                            "delta": item["delta"],
                            "delta_over_delta0": item["delta_over_delta0"],
                            "max_true": item["max_true_record"],
                            "true_at_binding": item["true"],
                            "bound_at_binding": item["bound"],
                            "S_at_binding": item["sp"],
                            "eta_used": item["eta_used"],
                            "kappa_eta_used": kappa * item["eta_used"],
                            "ratio": item["ratio"],
                            "failing_points": item["failing_points"],
                            "pass": tag == "threshold_prev",
                            "n": item["n"], "x": item["x"],
                            "region": item["region"], "tag": tag})
                end = one_delta(sign * 1000.0 * d0, variant)
                results[f"{variant}/{'+' if sign > 0 else '-'}"] = {
                    "delta0": d0,
                    "delta_range_over_delta0": [1e-3, 1e3],
                    "threshold": threshold,
                    "first_failure": first_fail,
                    "largest_holding_before_failure": last_hold,
                    "no_failure_up_to_1000_delta0": not fail.any(),
                    "max_ratio": float(ratio.max()),
                    "max_ratio_delta_over_delta0":
                        float(absd[int(np.argmax(ratio))] / d0),
                    "ratio_at_1000_delta0": float(ratio[-1]),
                    "at_1000_delta0": {
                        "max_true": end["max_true_record"],
                        "true_at_binding": end["true"],
                        "S_at_binding": end["sp"],
                        "eta_used_at_binding": end["eta_used"],
                        "kappa_eta_used": kappa * end["eta_used"],
                        "true_over_kappa_eta": (
                            end["true"] / (kappa * end["eta_used"])
                            if end["eta_used"] > 0 else None),
                        "bound": end["bound"], "ratio": end["ratio"]},
                    "recovered_after_failure": recovered,
                    "table": table,
                }
        return {
            "arm": arm, "lane": lane_name, "design": kind,
            "planted_component_index": PLANT_IDX,
            "support_4axis": [fx.inst.AXES[i]
                              for i in sorted(fx.inst.SUPPORTS[PLANT_IDX])],
            "support_6axis": [pa.AXES[i]
                              for i in sorted(pa.SUPPORTS[PLANT_IDX])],
            "channels_4axis": [int(fx.inst.CHOICES[PLANT_IDX][i])
                               for i in sorted(fx.inst.SUPPORTS[PLANT_IDX])],
            "channels_6axis": [int(pa.CHOICES[PLANT_IDX][i])
                               for i in sorted(pa.SUPPORTS[PLANT_IDX])],
            "order": int(pa.ORDERS[PLANT_IDX]),
            "G": G_idx,
            "why_outside": (
                "support {pack, form} is a subset neither of the upward-step "
                "cone nor of the corrected packed-sum cone "
                "{route,scheme,partition,pack}; it is a structural-zero "
                "(Q) channel"),
            "kappa_joint": kappa,
            "kappa_joint_beta_residual_sup": kappa_beta,
            "delta0_definition": "kappa_joint * max_t eta_Q(t)",
            "delta0": d0,
            "eta_scalar": float(eta.max()),
            "max_abs_c_target": max_cT,
            "max_S": max_S,
            "max_M": max_M,
            "alias_into_permitted": {
                "target_charge_per_delta": float((G_T * np.abs(bT)).sum()),
                "G_planted": G_idx,
                "max_abs_alias_T": float(np.abs(bT).max()),
                "max_abs_alias_R0": (
                    float(np.abs(r0_alias).max())
                    if r0_alias.size else 0.0),
                "sum_G_abs_alias_R0": (
                    float((r0_G * np.abs(r0_alias)).sum())
                    if r0_alias.size else 0.0),
                "target_components_aliased": t_components,
                "r0_components_aliased": r0_components,
                "alias_components": alias_table(
                    bvec, idx_map, supports, axis_names, Gvec, orders,
                    choices),
            },
            "variants": results,
        }

    out["path"] = run_design("path")
    out["pooled"] = run_design("pooled")
    return out


# --------------------------------------------------------------------------
# C. Frontier ratio distribution
# --------------------------------------------------------------------------


def distribution(path, pooled, d, lane_name, fig_scaled, meas):
    out = {}
    for kind in ("path", "pooled"):
        tol = fig_scaled
        if kind == "path":
            Nfull = np.max([b["Nfull"] for b in path], axis=0)
            post = np.max([np.abs(b["M"]).max(axis=1)
                           + cf.KJ21 * b["eta"] for b in path], axis=0)
            measure = meas["path"][lane_name]
            full = tol + Nfull
            family = tol + post + measure
            cells = []
            for b in path:
                cells.append((tol[:, None] + np.abs(b["M"])
                              + cf.KJ21 * b["eta"][:, None] + measure)
                             / (tol[:, None] + np.abs(b["N"])))
            cell_ratio = np.max(cells, axis=0).max(axis=1)
            Mmax = np.max([np.abs(b["M"]).max(axis=1) for b in path],
                          axis=0)
            etamax = np.max([b["eta"] for b in path], axis=0)
        else:
            Nfull = pooled["Nfull"]
            Mmax = np.abs(pooled["M"]).max(axis=1)
            etamax = pooled["eta"]
            measure = meas["pooled"][lane_name]
            post = Mmax + d["kappa_joint"] * etamax
            full = tol + Nfull
            family = tol + post + measure
            cell_ratio = ((tol[:, None] + np.abs(pooled["M"])
                           + d["kappa_joint"] * pooled["eta"][:, None]
                           + measure)
                          / (tol[:, None] + np.abs(pooled["N"]))).max(axis=1)
        ratio = family / full
        p = int(np.argmax(ratio))
        nfull_over_tol = float(Nfull[p] / tol[p]) if tol[p] > 0 else np.inf
        out[kind] = {
            "pointwise": {
                "min": float(ratio.min()),
                "median": float(np.median(ratio)),
                "p95": float(np.percentile(ratio, 95)),
                "max": float(ratio.max()),
                "frac_le_1.01": float((ratio <= 1.01).mean()),
                "frac_le_1.05": float((ratio <= 1.05).mean()),
                "frac_le_1.10": float((ratio <= 1.10).mean()),
            },
            "cellwise_max": {
                "min": float(cell_ratio.min()),
                "median": float(np.median(cell_ratio)),
                "p95": float(np.percentile(cell_ratio, 95)),
                "max": float(cell_ratio.max()),
                "frac_le_1.01": float((cell_ratio <= 1.01).mean()),
                "frac_le_1.05": float((cell_ratio <= 1.05).mean()),
                "frac_le_1.10": float((cell_ratio <= 1.10).mean()),
            },
            "worst_point": {
                "n": None, "x": None, "region": None,
                "ratio": float(ratio[p]),
                "tilde_figure": float(tol[p]),
                "Nfull": float(Nfull[p]),
                "M_max": float(Mmax[p]),
                "eta_Q": float(etamax[p]),
                "Nfull_over_tilde": nfull_over_tol,
                "signal_bearing": bool(Nfull[p] > 0.01 * tol[p]),
                "label": ("signal-bearing" if Nfull[p] > 0.01 * tol[p]
                          else "figure-only (<1% of the figure)"),
            },
            "family_bound_max": float(family.max()),
            "full_bound_max": float(full.max()),
            "family_bound_min": float(family.min()),
            "full_bound_min": float(full.min()),
            "point_ratios": ratio,
            "cell_ratios": cell_ratio,
            "Nfull_over_tilde": Nfull / tol,
        }
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--design", choices=("d84", "d75"), default="d84")
    parser.add_argument("--out-dir", default=HERE)
    parser.add_argument("--prefix", default=None)
    args = parser.parse_args()
    prefix = args.prefix or f"prospective_{args.design}"
    os.makedirs(args.out_dir, exist_ok=True)

    started = time.perf_counter()
    d = fx.build(args.design)
    if args.design == "d84":
        fx.verify_against_certificate(d)
    meas = measure_factors(d)
    out = {
        "meta": {
            "design": args.design,
            "cells": int(len(d["cells"])),
            "permitted": int(len(d["P_IDX"])),
            "complement": int(len(d["Q_IDX"])),
            "targets": int(len(d["T_IDX"])),
            "kappa_joint": float(d["kappa_joint"]),
            "harness": (
                "harness repository resolved from BOYS_HARNESS_REPO or the "
                "script's parent checkout"),
            "records": {arm: rf.sha256_file(path)
                        for arm, path in arm_paths().items()},
            "seed": SEED,
            "n_point_splits": N_SPLITS,
            "n_cell_splits": N_CELL_SPLITS,
            "cell_split_seed": SEED + 1,
            "figure": "tilde_figure(t) = max_cells figure(cell,t)/F0(x)",
            "theorem_17": (
                "|N(a)| <= sum_T G_a|chat_a| + kappa_joint eta_Q^+ + "
                "sum_T G_a delta_a; N(a)=sum_T c_a Phi_a(a)+"
                "sum_Q c_b Phi_b(a)"),
            "measurement_epsilon": EPS,
            "measurement_terms": {"path": meas["path"],
                                  "pooled": meas["pooled"]},
            "protocol": (
                "worked-check/prospective/additions.py at SHA-256 "
                "546e0304b438dd1d029ddc45a3d5b0ab016ab340f317d842a22b913"
                "84690c87a, run with the corrected pooled design"),
            "runtime_seconds": None,
        },
        "A_transport": {}, "B_planted": {}, "C_distribution": {},
    }
    rng = np.random.default_rng(SEED)
    splits = [rng.permutation(858) for _ in range(N_SPLITS)]
    rng2 = np.random.default_rng(SEED + 1)
    masks36 = [rng2.permutation(36)[:18] for _ in range(N_CELL_SPLITS)]
    masks144 = [rng2.permutation(144)[:72] for _ in range(N_CELL_SPLITS)]
    dist_rows = []
    fig_checks = {}
    for arm, path in arm_paths().items():
        print(f"loading {arm}", flush=True)
        run = fc.load_route(path, arm)
        cfig = rf.load_figures(path)
        regions = region_array(run["xs"], run["ns"])
        for lane_name in LANES:
            lane = fc.LANES.index(lane_name)
            scaled = rf.scaled_figure(cfig["FIG"], cfig["FN"], lane_name,
                                      len(run["xs"]))
            tilde = scaled["tilde"]
            fig_checks[f"{arm}/{lane_name}"] = {
                "max_rel_dev": scaled["max_rel_dev"],
                "max_abs_dev": scaled["max_abs_dev"],
            }
            path_blocks, pooled = precompute(run["R"][lane], d)
            print(f"  {arm}/{lane_name}: A", flush=True)
            a = transport_test(path_blocks, pooled, d, lane_name, meas,
                               run["xs"], run["ns"], splits, regions)
            a["cell_subsample"] = cell_transport_test(
                run["R"][lane], path_blocks, pooled, d, lane_name, meas,
                run["xs"], run["ns"], masks36, masks144)
            out["A_transport"][f"{arm}/{lane_name}"] = a
            print(f"  {arm}/{lane_name}: B", flush=True)
            b = planted_test(path_blocks, pooled, d, lane_name, meas,
                             run["xs"], run["ns"], arm)
            out["B_planted"][f"{arm}/{lane_name}"] = b
            print(f"  {arm}/{lane_name}: C", flush=True)
            dist = distribution(path_blocks, pooled, d, lane_name, tilde,
                                meas)
            nx = len(run["xs"])
            for kind in ("path", "pooled"):
                wp = dist[kind]["worst_point"]
                p = int(np.argmax(dist[kind]["point_ratios"]))
                wp["n"] = int(run["ns"][p // nx])
                wp["x"] = float(run["xs"][p % nx])
                wp["region"] = fc.region_of(run["xs"][p % nx])
                for q in range(len(tilde)):
                    dist_rows.append([
                        arm, lane_name, kind,
                        int(run["ns"][q // nx]),
                        repr(float(run["xs"][q % nx])),
                        fc.region_of(run["xs"][q % nx]),
                        repr(float(tilde[q])),
                        repr(float(dist[kind]["point_ratios"][q])),
                        repr(float(dist[kind]["cell_ratios"][q])),
                        repr(float(dist[kind]["Nfull_over_tilde"][q])),
                    ])
            out["C_distribution"][f"{arm}/{lane_name}"] = dist
        print(f"  done {arm} ({time.perf_counter() - started:.1f} s)",
              flush=True)

    out["meta"]["runtime_seconds"] = time.perf_counter() - started
    out["meta"]["figure_checks"] = fig_checks
    for v in out["C_distribution"].values():
        for kind in ("path", "pooled"):
            v[kind].pop("point_ratios", None)
            v[kind].pop("cell_ratios", None)
            v[kind].pop("Nfull_over_tilde", None)
    with open(os.path.join(args.out_dir, f"{prefix}_results.json"), "w",
              encoding="utf-8") as handle:
        json.dump(out, handle, indent=1, default=float)

    with open(os.path.join(args.out_dir, f"{prefix}_transport.csv"), "w",
              encoding="utf-8", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(["arm", "lane", "design", "variant", "coverage",
                    "cellwise_coverage", "worst_ratio", "worst_n",
                    "worst_x", "worst_region", "median_ratio", "p95_ratio",
                    "split_coverage_min", "split_coverage_max",
                    "split_coverage_std", "cell_tight_point_coverage",
                    "cell_tight_cell_coverage",
                    "cell_tight_worst_ratio"])
        for key, v in out["A_transport"].items():
            arm, lane = key.split("/")
            for kind in ("path", "pooled"):
                dd = v[kind]
                w.writerow([arm, lane, kind, "scalar_cross_point",
                            repr(dd["coverage"]),
                            repr(dd["cellwise_coverage"]),
                            repr(dd["worst_ratio"]["ratio"]),
                            dd["worst_ratio"]["n"],
                            repr(dd["worst_ratio"]["x"]),
                            dd["worst_ratio"]["region"],
                            repr(dd["ratio_stats"]["median"]),
                            repr(dd["ratio_stats"]["p95"]),
                            repr(dd["split_coverage_min"]),
                            repr(dd["split_coverage_max"]),
                            repr(dd["split_coverage_std"]),
                            repr(dd["cell_tight_point_coverage"]),
                            repr(dd["cell_tight_cell_coverage"]),
                            repr(dd["cell_tight_worst_ratio"]["ratio"])])
                c = v["cell_subsample"][kind]
                w.writerow([arm, lane, kind, "cell_subsample",
                            repr(c["coverage"]),
                            repr(c["cellwise_coverage"]),
                            repr(c["worst_ratio"]["ratio"]),
                            c["worst_ratio"]["n"],
                            repr(c["worst_ratio"]["x"]),
                            c["worst_ratio"]["region"],
                            repr(c["ratio_stats"]["median"]),
                            repr(c["ratio_stats"]["p95"]), "", "", "",
                            "", "", ""])

    header = ["arm", "lane", "design", "variant", "sign",
              "delta_over_delta0", "delta", "max_true",
              "true_at_binding", "S_at_binding", "eta_used",
              "kappa_eta_used", "bound", "ratio", "failing_points",
              "pass", "n", "x", "region", "tag"]
    frozen = ("global_frozen", "pointwise_frozen")
    consistent = ("self_consistent_exact", "self_consistent_upper")
    with open(os.path.join(args.out_dir, f"{prefix}_planted.csv"), "w",
              encoding="utf-8", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(header)
        for key, v in out["B_planted"].items():
            arm, lane = key.split("/")
            for kind in ("path", "pooled"):
                for vk, r in v[kind]["variants"].items():
                    variant, sign = vk.rsplit("/", 1)
                    if not variant.startswith(frozen):
                        continue
                    for row in r["table"]:
                        w.writerow([arm, lane, kind, variant, sign,
                                    repr(row["delta_over_delta0"]),
                                    repr(row["delta"]),
                                    repr(row["max_true"]),
                                    repr(row["true_at_binding"]),
                                    repr(row["S_at_binding"]),
                                    repr(row["eta_used"]),
                                    repr(row["kappa_eta_used"]),
                                    repr(row["bound_at_binding"]),
                                    repr(row["ratio"]),
                                    row["failing_points"],
                                    int(row["pass"]), row["n"],
                                    repr(row["x"]), row["region"],
                                    row.get("tag", "")])
    with open(os.path.join(args.out_dir,
                           f"{prefix}_planted_selfconsistent.csv"), "w",
              encoding="utf-8", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(header)
        for key, v in out["B_planted"].items():
            arm, lane = key.split("/")
            for kind in ("path", "pooled"):
                dd = v[kind]
                for vk, r in dd["variants"].items():
                    variant, sign = vk.rsplit("/", 1)
                    if not variant.startswith(consistent):
                        continue
                    end = r["at_1000_delta0"]
                    w.writerow([
                        arm, lane, kind, variant, sign,
                        "1000.0", repr(1000.0 * dd["delta0"]),
                        repr(end["max_true"]),
                        repr(end["true_at_binding"]),
                        repr(end["S_at_binding"]),
                        repr(end["eta_used_at_binding"]),
                        repr(end["kappa_eta_used"]), repr(end["bound"]),
                        repr(end["ratio"]), "",
                        int(end["ratio"] <= 1.0),
                        "", "", "", "at_1000_delta0"])
                    for row in r["table"]:
                        w.writerow([arm, lane, kind, variant, sign,
                                    repr(row["delta_over_delta0"]),
                                    repr(row["delta"]),
                                    repr(row["max_true"]),
                                    repr(row["true_at_binding"]),
                                    repr(row["S_at_binding"]),
                                    repr(row["eta_used"]),
                                    repr(row["kappa_eta_used"]),
                                    repr(row["bound_at_binding"]),
                                    repr(row["ratio"]),
                                    row["failing_points"],
                                    int(row["pass"]), row["n"],
                                    repr(row["x"]), row["region"],
                                    row.get("tag", "")])

    with open(os.path.join(args.out_dir, f"{prefix}_distribution.csv"), "w",
              encoding="utf-8", newline="") as handle:
        w = csv.writer(handle)
        w.writerow(["arm", "lane", "design", "n", "x", "region",
                    "tilde_figure", "ratio_pointwise",
                    "ratio_cellwise_max", "Nfull_over_tilde"])
        w.writerows(dist_rows)

    print("\nA. scalar transport coverage (100 splits)")
    for key, v in out["A_transport"].items():
        for kind in ("path", "pooled"):
            dd = v[kind]
            print(f"  {key:16s} {kind:7s} cov {dd['coverage']:.6f} "
                  f"ctpt {dd['cell_tight_point_coverage']:.6f} "
                  f"ctcell {dd['cell_tight_cell_coverage']:.6f} "
                  f"worst {dd['worst_ratio']['ratio']:.6f} "
                  f"ctworst {dd['cell_tight_worst_ratio']['ratio']:.6f} "
                  f"med {dd['ratio_stats']['median']:.6f} "
                  f"p95 {dd['ratio_stats']['p95']:.6f}")
    print("\nA. cell-subsample coverage")
    for key, v in out["A_transport"].items():
        for kind in ("path", "pooled"):
            dd = v["cell_subsample"][kind]
            print(f"  {key:16s} {kind:7s} cov {dd['coverage']:.6f} "
                  f"cell {dd['cellwise_coverage']:.6f} "
                  f"worst {dd['worst_ratio']['ratio']:.6f} "
                  f"med {dd['ratio_stats']['median']:.6f}")
    print("\nB. planted: frozen allowance at delta0 and thresholds")
    for key, v in out["B_planted"].items():
        for kind in ("path", "pooled"):
            dd = v[kind]
            for variant in ("global_frozen", "pointwise_frozen"):
                for sign in ("+", "-"):
                    r = dd["variants"][f"{variant}/{sign}"]
                    row0 = [t for t in r["table"]
                            if t["delta_over_delta0"] == 1.0][0]
                    th = r["threshold"]
                    print(f"  {key:16s} {kind:7s} {variant:17s} {sign} "
                          f"delta0 {dd['delta0']:.6e} "
                          f"ratio(1d0) {row0['ratio']:.6f} "
                          f"maxratio {r['max_ratio']:.6f} "
                          f"frozen*={'none' if th is None else format(th['delta_over_delta0'], '.6f')}")
            print(f"  {key:16s} {kind:7s} aliasT "
                  f"{dd['alias_into_permitted']['max_abs_alias_T']:.6e} "
                  f"charge/delta "
                  f"{dd['alias_into_permitted']['target_charge_per_delta']:.6f}")
    print("\nB. planted: self-consistent allowance")
    for key, v in out["B_planted"].items():
        for kind in ("path", "pooled"):
            dd = v[kind]
            for variant in ("self_consistent_exact", "self_consistent_upper"):
                r = dd["variants"][f"{variant}/+"]
                row0 = [t for t in r["table"]
                        if t["delta_over_delta0"] == 1.0][0]
                th = r["threshold"]
                end = r["at_1000_delta0"]
                print(f"  {key:16s} {kind:7s} {variant:22s} "
                      f"ratio(1d0) {row0['ratio']:.6f} "
                      f"th={'none' if th is None else format(th['delta_over_delta0'], '.6f')} "
                      f"maxratio {r['max_ratio']:.6f} "
                      f"at1000 {end['ratio']:.6f}")
    print("\nC. pointwise ratio distribution")
    for key, v in out["C_distribution"].items():
        for kind in ("path", "pooled"):
            dd = v[kind]["pointwise"]
            print(f"  {key:16s} {kind:7s} min {dd['min']:.6f} "
                  f"med {dd['median']:.6f} p95 {dd['p95']:.6f} "
                  f"max {dd['max']:.6f} <=1.01 {dd['frac_le_1.01']:.4f} "
                  f"<=1.05 {dd['frac_le_1.05']:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
