#!/usr/bin/env python3
"""Definitive resolution of the tensor-vs-direct order-3 materiality disagreement.

Reads the deposited 576-cell record with the same loader, basis and transforms
as C:\\Users\\myaml\\dev\\_d-fullcube\\fullcube_decomposition.py (imported, not
retyped), reproduces both the per-axis tensor path and the direct 144-term
path, and then computes the materiality of every order-3 component at every
A/B point in 60-digit arithmetic on the exact recorded doubles.

Outputs verification.json beside this script and a human-readable summary on
stdout.  Read-only with respect to the record, the decomposition and the
papers.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from mpmath import mp

LANDING = Path(__file__).resolve().parent
DEFAULT_SCRIPT = Path(__file__).resolve().parent / "fullcube_decomposition.py"
DEFAULT_FUSED = (Path(__file__).resolve().parents[1]
                / "records" / "cells" / "cells-msvc-fused.csv.xz")
DEFAULT_SEPARATE = (Path(__file__).resolve().parents[1]
                   / "records" / "cells" / "cells-msvc-separate.csv.xz")
DEFAULT_SUMMARY = Path(__file__).resolve().parent / "fullcube_summary.json"


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("fcd", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mp_axis_matrices(fcd):
    matrices = []
    for count in fcd.LEVEL_COUNTS:
        if count == 2:
            rows = [
                [mp.mpf(1), mp.mpf(1)],
                [mp.mpf(1), mp.mpf(-1)],
            ]
        elif count == 3:
            s = mp.sqrt(mp.mpf(3) / 2)
            q = mp.mpf(1) / mp.sqrt(mp.mpf(2))
            rows = [
                [mp.mpf(1), mp.mpf(1), mp.mpf(1)],
                [s, mp.mpf(0), -s],
                [q, -2 * q, q],
            ]
        else:
            raise ValueError(count)
        matrices.append(rows)
    return matrices


def axis_maps(fcd):
    maps = []
    for axis, count in enumerate(fcd.LEVEL_COUNTS):
        stride = fcd.STRIDES[axis]
        entries = []
        for index in range(fcd.N_CELLS):
            outer = index // (stride * count)
            j = (index // stride) % count
            inner = index % stride
            base = outer * stride * count + inner
            entries.append((index, j, base))
        maps.append(entries)
    return maps


def mp_forward(values, matrices, maps, fcd):
    """Exact (60-digit) per-axis transform, identical algebra to tensor_forward."""
    array = [mp.mpf(v) for v in values]
    for axis, (matrix, entries) in enumerate(zip(matrices, maps)):
        count = fcd.LEVEL_COUNTS[axis]
        stride = fcd.STRIDES[axis]
        new = [None] * fcd.N_CELLS
        for index, j, base in entries:
            row = matrix[j]
            acc = mp.mpf(0)
            for level in range(count):
                coefficient = row[level]
                if coefficient == 0:
                    continue
                acc += coefficient * array[base + level * stride]
            new[index] = acc
        array = new
    return array


def exact_point(values, matrices, maps, fcd, one_over):
    c = mp_forward(values, matrices, maps, fcd)
    c = [v * one_over for v in c]
    d0 = mp.fsum([v * v for v in c[1:]])
    return c, d0


def max_phi_weights(phi):
    return np.abs(phi).max(axis=1)


def component_description(fcd, index):
    comp = fcd.COMPONENTS[index]
    return {
        "component_index": int(index),
        "support": list(comp["support_axes"]),
        "channels": list(comp["channels"]),
        "channels_1based": list(comp["channels_1based"]),
        "order": int(comp["order"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--script", default=str(DEFAULT_SCRIPT))
    parser.add_argument("--fused", default=str(DEFAULT_FUSED))
    parser.add_argument("--separate", default=str(DEFAULT_SEPARATE))
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY))
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--out", default=str(LANDING / "verification.json"))
    args = parser.parse_args()

    started = time.perf_counter()
    mp.dps = args.dps
    fcd = load_module(Path(args.script))
    summary = json.loads(Path(args.summary).read_text(encoding="utf-8"))

    phi = fcd.build_phi()
    components = fcd.COMPONENTS
    mask3 = fcd.component_selectors(components, 3)
    order3_indices = [int(v) for v in mask3]
    weights = max_phi_weights(phi)
    matrices = mp_axis_matrices(fcd)
    maps = axis_maps(fcd)
    one_over = mp.mpf(1) / fcd.N_CELLS

    report = {
        "generated_local": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "dps": args.dps,
        "script": str(Path(args.script)),
        "inputs": {
            "fused": str(Path(args.fused)),
            "separate": str(Path(args.separate)),
            "summary": str(Path(args.summary)),
        },
        "definitions": {
            "denominator": "D0 = sum_{k != 0} c_k^2, the point's centred (non-mean) total sum of squares",
            "threshold": "ALPHA = 0.1, strict > in both code paths",
            "components": "the 44 order-3 components; a component is material if its own share exceeds ALPHA at any A/B point",
            "tensor_path": "full_cube[route][lane].orders['3'].material_components_10pct (858 points; A/B mask equals all here)",
            "direct_path": "checks.coefficient_transform_agreement[route].order3_share_stability[lane].direct_matmul (825 A/B points)",
        },
        "reproduction": {},
        "exact_counts": {},
        "per_route_lane": {},
    }

    for route, path in (("fused", args.fused), ("separate", args.separate)):
        load_started = time.perf_counter()
        run = fcd.load_route(str(path), route)
        points = run["R"].shape[1]
        tensor_raw = fcd.tensor_forward(run["R"])
        direct_raw = run["R"].reshape(-1, fcd.N_CELLS) @ phi.T / fcd.N_CELLS
        tensor = tensor_raw.reshape(len(fcd.LANES), points, fcd.N_CELLS)
        direct = direct_raw.reshape(len(fcd.LANES), points, fcd.N_CELLS)
        regions = np.array(
            [fcd.region_of(run["xs"][p % len(run["xs"])]) for p in range(points)]
        )
        nonfar = regions != "far"
        nonfar_indices = np.where(nonfar)[0]
        report["reproduction"][route] = {"load_seconds": time.perf_counter() - load_started}
        report["per_route_lane"][route] = {}

        for lane in fcd.RESOLVING_LANES:
            lane_index = fcd.LANES.index(lane)
            lane_started = time.perf_counter()
            ct = tensor[lane_index]
            cd = direct[lane_index]
            dt = np.square(ct[:, 1:]).sum(axis=1)
            dd = np.square(cd[:, 1:]).sum(axis=1)
            st = np.zeros((points, len(mask3)))
            sd = np.zeros((points, len(mask3)))
            ok_t = dt > 0
            ok_d = dd > 0
            st[ok_t] = np.square(ct[ok_t][:, mask3]) / dt[ok_t, None]
            sd[ok_d] = np.square(cd[ok_d][:, mask3]) / dd[ok_d, None]

            order_share_t = st.sum(axis=1)
            order_share_d = sd.sum(axis=1)

            component_max_t = st[nonfar].max(axis=0)
            component_max_d = sd[nonfar].max(axis=0)
            material_t = np.where(component_max_t > fcd.ALPHA)[0]
            material_d = np.where(component_max_d > fcd.ALPHA)[0]

            # Summary cross-check.
            summary_tensor_count = summary["full_cube"][route][lane]["orders"]["3"][
                "material_components_10pct"
            ]
            summary_tensor_nonfar = summary["full_cube"][route][lane]["orders"]["3"][
                "material_components_10pct_nonfar"
            ]
            stability = summary["checks"]["coefficient_transform_agreement"][route][
                "order3_share_stability"
            ][lane]
            summary_direct_count = stability["direct_matmul"][
                "order3_material_components_AB"
            ]
            summary_tensor_max = stability["tensor"]["max_order3_share_AB"]
            summary_direct_max = stability["direct_matmul"]["max_order3_share_AB"]

            entry = {
                "tensor": {
                    "max_order_share_AB": float(order_share_t[nonfar].max()),
                    "max_order_share_AB_at": fcd.point_location(
                        run, int(nonfar_indices[int(np.argmax(order_share_t[nonfar]))])
                    ),
                    "component_max_shares": {
                        int(mask3[c]): float(component_max_t[c]) for c in material_t
                    },
                    "material_component_count": int(len(material_t)),
                    "material_component_count_all_points": int(
                        (st.max(axis=0) > fcd.ALPHA).sum()
                    ),
                    "summary_count": summary_tensor_count,
                    "summary_count_nonfar": summary_tensor_nonfar,
                    "summary_max_order_share": summary_tensor_max,
                },
                "direct": {
                    "max_order_share_AB": float(order_share_d[nonfar].max()),
                    "max_order_share_AB_at": fcd.point_location(
                        run, int(nonfar_indices[int(np.argmax(order_share_d[nonfar]))])
                    ),
                    "component_max_shares": {
                        int(mask3[c]): float(component_max_d[c]) for c in material_d
                    },
                    "material_component_count": int(len(material_d)),
                    "summary_count": summary_direct_count,
                    "summary_max_order_share": summary_direct_max,
                },
                "candidate_readings": [],
            }

            # Every (component, point) reading above ALPHA in either path.
            seen = set()
            for path_name, shares in (("tensor", st), ("direct", sd)):
                for local_component in sorted(set(material_t) | set(material_d)):
                    points_above = np.where(shares[nonfar, local_component] > fcd.ALPHA)[0]
                    for point_local in points_above:
                        key = (int(local_component), int(point_local))
                        if key in seen:
                            continue
                        seen.add(key)
                        point = int(nonfar_indices[point_local])
                        component_index = int(mask3[local_component])
                        reading = {
                            "found_by": path_name,
                            **component_description(fcd, component_index),
                            "point": fcd.point_location(run, point),
                            "tensor": {
                                "c": float(ct[point, component_index]),
                                "D0": float(dt[point]),
                                "share": float(st[point, local_component]),
                            },
                            "direct": {
                                "c": float(cd[point, component_index]),
                                "D0": float(dd[point]),
                                "share": float(sd[point, local_component]),
                            },
                        }
                        entry["candidate_readings"].append(reading)

            # Exact pass over every A/B point for all 44 order-3 components.
            exact_component_max = {index: mp.mpf(0) for index in order3_indices}
            exact_component_at = {index: None for index in order3_indices}
            exact_order_max = mp.mpf(0)
            exact_order_at = None
            exact_zero_denominators = 0
            exact_nonzero_min_d0 = None
            constant_points = []
            for point in nonfar_indices:
                point = int(point)
                c, d0 = exact_point(
                    run["R"][lane_index, point], matrices, maps, fcd, one_over
                )
                unique_ratios = np.unique(run["R"][lane_index, point]).size
                if unique_ratios == 1:
                    constant_points.append(fcd.point_location(run, point))
                if d0 == 0:
                    exact_zero_denominators += 1
                    continue
                if exact_nonzero_min_d0 is None or d0 < exact_nonzero_min_d0:
                    exact_nonzero_min_d0 = d0
                order_total = mp.mpf(0)
                for local, index in enumerate(order3_indices):
                    share = c[index] * c[index] / d0
                    order_total += share
                    if share > exact_component_max[index]:
                        exact_component_max[index] = share
                        exact_component_at[index] = fcd.point_location(run, point)
                if order_total > exact_order_max:
                    exact_order_max = order_total
                    exact_order_at = fcd.point_location(run, point)

            exact_material = sorted(
                index
                for index in order3_indices
                if exact_component_max[index] > mp.mpf("0.1")
            )
            entry["exact"] = {
                "dps": args.dps,
                "material_component_count": len(exact_material),
                "material_components": {
                    str(index): {
                        **component_description(fcd, index),
                        "max_share": mp.nstr(exact_component_max[index], 20),
                        "max_share_float": float(exact_component_max[index]),
                        "at": exact_component_at[index],
                    }
                    for index in exact_material
                },
                "all_component_maxima": {
                    str(index): {
                        "max_share": mp.nstr(exact_component_max[index], 20),
                        "at": exact_component_at[index],
                    }
                    for index in order3_indices
                },
                "max_order3_share": mp.nstr(exact_order_max, 20),
                "max_order3_share_float": float(exact_order_max),
                "max_order3_share_at": exact_order_at,
                "zero_denominator_points": exact_zero_denominators,
                "smallest_nonzero_D0": mp.nstr(exact_nonzero_min_d0, 10)
                if exact_nonzero_min_d0 is not None
                else None,
                "constant_points": constant_points,
                "seconds": time.perf_counter() - lane_started,
            }

            # Detailed exact adjudication of every candidate reading.
            for reading in entry["candidate_readings"]:
                point = None
                loc = reading["point"]
                for candidate in range(points):
                    candidate_loc = fcd.point_location(run, candidate)
                    if (
                        candidate_loc["n"] == loc["n"]
                        and candidate_loc["x"] == loc["x"]
                    ):
                        point = candidate
                        break
                row = run["R"][lane_index, point]
                c, d0 = exact_point(row, matrices, maps, fcd, one_over)
                unique_ratios = np.unique(row)
                exact_ss = (
                    "0 (all 144 recorded ratios are the same double)"
                    if unique_ratios.size == 1
                    else mp.nstr(d0 * fcd.N_CELLS, 20)
                )
                index = reading["component_index"]
                norm_sq = mp.fsum([v * v for v in c])
                max_r = mp.mpf(float(np.abs(row).max()))
                reconstructed = float(np.abs(ct[point]).max())
                share_hp = c[index] * c[index] / d0 if d0 != 0 else mp.mpf(0)
                reading["exact"] = {
                    "unique_recorded_ratios": int(unique_ratios.size),
                    "exact_centred_ss": exact_ss,
                    "D0": mp.nstr(d0, 20),
                    "D0_float": float(d0),
                    "c": mp.nstr(c[index], 20),
                    "c_float": float(c[index]),
                    "share_non_mean_total": mp.nstr(share_hp, 20),
                    "share_non_mean_total_float": float(share_hp),
                }
                reading["alternative_denominators"] = {
                    "coefficient_norm_sq_incl_mean": mp.nstr(norm_sq, 20),
                    "share_of_coefficient_norm_sq": mp.nstr(
                        c[index] * c[index] / norm_sq if norm_sq != 0 else mp.mpf(0),
                        20,
                    ),
                    "max_abs_Phi_component": float(weights[index]),
                    "max_abs_r_point": float(np.abs(row).max()),
                    "max_abs_r_hat_point": float(reconstructed),
                    "response_ratio_over_max_r": mp.nstr(
                        abs(c[index]) * mp.mpf(float(weights[index])) / max_r
                        if max_r != 0
                        else mp.mpf(0),
                        20,
                    ),
                    "response_ratio_over_max_r_hat": mp.nstr(
                        abs(c[index]) * mp.mpf(float(weights[index])) / reconstructed
                        if reconstructed != 0
                        else mp.mpf(0),
                        20,
                    ),
                }
                reading["tensor"]["share_at_point"] = float(st[point, list(mask3).index(index)])
                reading["direct"]["share_at_point"] = float(sd[point, list(mask3).index(index)])

            # Full-record counts under the alternative denominators (numpy,
            # tensor coefficients; exact counts confirm the D0 entry).
            alt = {}
            # coefficient norm squared, mean included
            norm_all = np.square(ct).sum(axis=1)
            safe_norm = np.where(norm_all > 0, norm_all, np.inf)
            norm_share = np.square(ct[:, mask3]) / safe_norm[:, None]
            alt["coefficient_norm_sq_incl_mean"] = {
                "material_component_count": int(
                    (norm_share[nonfar].max(axis=0) > fcd.ALPHA).sum()
                ),
                "max_component_share": float(norm_share[nonfar].max()),
            }
            # response ratio: |c| max|Phi| / max|r| at the point
            max_r = np.abs(run["R"][lane_index]).max(axis=1)
            safe_r = np.where(max_r > 0, max_r, np.inf)
            response_ratio = (
                np.abs(ct[:, mask3]) * weights[mask3][None, :] / safe_r[:, None]
            )
            alt["max_phi_weight_over_max_r_point"] = {
                "material_component_count": int(
                    (response_ratio[nonfar].max(axis=0) > fcd.ALPHA).sum()
                ),
                "max_component_ratio": float(response_ratio[nonfar].max()),
            }
            entry["alternative_denominators_full_record"] = alt

            report["reproduction"][route][lane] = {
                "tensor_reproduced": int(len(material_t)),
                "tensor_summary": summary_tensor_count,
                "tensor_reproduced_nonfar": int(
                    (component_max_t > fcd.ALPHA).sum()
                ),
                "tensor_summary_nonfar": summary_tensor_nonfar,
                "direct_reproduced": int(len(material_d)),
                "direct_summary": summary_direct_count,
                "tensor_max_order_share_reproduced": float(order_share_t[nonfar].max()),
                "tensor_max_order_share_summary": summary_tensor_max,
                "direct_max_order_share_reproduced": float(order_share_d[nonfar].max()),
                "direct_max_order_share_summary": summary_direct_max,
            }
            report["exact_counts"].setdefault(route, {})[lane] = len(exact_material)
            report["per_route_lane"][route][lane] = entry
            print(
                f"{route} {lane}: tensor {len(material_t)}/44 direct {len(material_d)}/44 "
                f"exact {len(exact_material)}/44 "
                f"(summary {summary_tensor_count}/{summary_direct_count}), "
                f"exact max order-3 share {mp.nstr(exact_order_max, 8)} "
                f"at n={exact_order_at['n']}, x={exact_order_at['x']:g}; "
                f"{time.perf_counter() - lane_started:.1f} s",
                flush=True,
            )

    report["runtime_seconds"] = time.perf_counter() - started
    Path(args.out).write_text(json.dumps(report, indent=1, default=str), encoding="utf-8")
    print(f"wrote {args.out} in {report['runtime_seconds']:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
