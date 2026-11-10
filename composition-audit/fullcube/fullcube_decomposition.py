#!/usr/bin/env python3
"""Full-cube mixed-level factorial decomposition of the composition record.

Paper D, full-cube decomposition lane.  This script reads the deposited
576-cell composition-audit record and computes, for every lane and every
(n, x) point, the exact decomposition of the signed F0-scaled error

    r(c) = ((value - ref_hi) - ref_lo) / F0

on the six call axes route, scheme, partition, pack, form, exp with level
counts (2, 2, 3, 2, 3, 2).  The product basis Phi_alpha uses the orthonormal
rows of the task:

    two levels      mean 1, contrast (1, -1);
    three levels    mean 1, linear sqrt(3/2)(1, 0, -1),
                    quadratic (1/sqrt(2))(1, -2, 1).

Every row has mean square 1 under the uniform measure over levels, zero mean
(the mean row excepted), and the rows of one axis are mutually orthogonal.
The coefficients and the reconstruction are

    c_alpha = (1/144) sum_cells r(c) Phi_alpha(c),
    r_hat(c) = sum_alpha c_alpha Phi_alpha(c).

The d79 corner-design estimator is also computed because the verification
note D-cube-contrasts.md reports its maxima.  That note's estimator uses
Euclidean-normalised rows psi (mean row 1/sqrt(L), contrast rows of unit
Euclidean length) and the support-only product:

    P_alpha = sum_{x in block(T_alpha)} r(x) psi_alpha(x),

where block(T_alpha) is the complete factorial on the support axes with every
other axis at its baseline level.  With B = prod_{i in T} L_i,

    P_alpha = sqrt(B) * (1/B) sum_block r(x) Phi_alpha(x),

so the note's P and this script's c are different scalings of the same block
projection.  Both are reported, and the note's values are compared against P,
which is the quantity its table calls the contrast.

Nothing here re-runs the driver or edits a record.

Usage:
    python -I fullcube_decomposition.py \
        --fused records/cells/cells-msvc-fused.csv.xz \
        --separate records/cells/cells-msvc-separate.csv.xz \
        --out-dir results
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# --- design -----------------------------------------------------------------

AXES = ["route", "scheme", "partition", "pack", "form", "exp"]
LEVELS = {
    "route": ["chebyshev", "rational-minimax"],
    "scheme": ["split-clenshaw", "horner"],
    "partition": ["shipped", "narrow", "uniform"],
    "pack": ["arguments", "orders"],
    "form": ["exact-division", "plain-reciprocal", "refined-reciprocal"],
    "exp": ["accurate", "fast"],
}
LEVEL_COUNTS = [len(LEVELS[a]) for a in AXES]          # [2, 2, 3, 2, 3, 2]
N_CELLS = math.prod(LEVEL_COUNTS)                       # 144
LANES = ["fp64", "fp32", "fp16", "bf16"]
RESOLVING_LANES = ["fp64", "fp32"]
STORE_LANES = ["fp16", "bf16"]

# Region selectors from the note (boys kX0 / kFlatHi).
KX0 = 11.899848152108484
KFLATHI = 31.0

MEASUREMENT_FLOOR = 1e-32
ALPHA = 0.10

COLUMNS = [
    "precision", "route", "scheme", "partition", "pack", "form", "exp",
    "n", "x_nominal", "x", "value", "ref_hi", "ref_lo", "abs_err", "f0",
    "fn", "figure", "half_digit", "limit", "accessor_ratio", "figure_ratio",
    "pack_live", "partition_in_domain", "argument_kind", "status",
]

# The note's reported d79 maxima, used only for comparison, never adjusted.
NOTE = {
    "order3_max": {
        ("fused", "fp64"): 2.26e-16,
        ("separate", "fp64"): 2.76e-16,
        ("fused", "fp32"): 1.40e-7,
        ("separate", "fp32"): 1.40e-7,
    },
    "order2_max": {
        ("fused", "fp64"): 3.26e-14,
        ("separate", "fp64"): 3.26e-14,
        ("fused", "fp32"): 1.98e-7,
        ("separate", "fp32"): 1.98e-7,
    },
    "order1_max": {
        ("fused", "fp64"): 4.45e-14,
        ("separate", "fp64"): 4.43e-14,
        ("fused", "fp32"): 2.79e-7,
        ("separate", "fp32"): 2.79e-7,
    },
    "order3_exact_zero_components": {
        ("fused", "fp64"): 12,
        ("separate", "fp64"): 11,
        ("fused", "fp32"): 28,
        ("separate", "fp32"): 28,
    },
    "order3_material_components": {
        ("fused", "fp64"): 17,
        ("separate", "fp64"): 20,
        ("fused", "fp32"): 13,
        ("separate", "fp32"): 14,
    },
    "order3_material_readings": {
        ("fused", "fp64"): 172,
        ("separate", "fp64"): 267,
        ("fused", "fp32"): 2271,
        ("separate", "fp32"): 2485,
    },
    "order3_max_relative": {
        ("fused", "fp64"): 0.631,
        ("separate", "fp64"): 0.703,
        ("fused", "fp32"): 1.382,
        ("separate", "fp32"): 1.382,
    },
}


# --- basis ------------------------------------------------------------------


def mixed_strides(level_counts):
    strides = []
    stride = 1
    for count in reversed(level_counts):
        strides.append(stride)
        stride *= count
    return list(reversed(strides))


STRIDES = mixed_strides(LEVEL_COUNTS)


def axis_basis(level_count):
    """Task rows: orthonormal under the uniform measure over levels."""
    if level_count == 2:
        rows = np.array([[1.0, 1.0], [1.0, -1.0]])
        names = ["mean", "contrast"]
    elif level_count == 3:
        rows = np.array(
            [
                [1.0, 1.0, 1.0],
                math.sqrt(1.5) * np.array([1.0, 0.0, -1.0]),
                (1.0 / math.sqrt(2.0)) * np.array([1.0, -2.0, 1.0]),
            ]
        )
        names = ["mean", "linear", "quadratic"]
    else:
        raise ValueError(f"unsupported level count {level_count}")
    return rows, names


def euclidean_rows(level_count):
    """Note rows: every row of unit Euclidean length (sum of squares 1)."""
    rows, names = axis_basis(level_count)
    return rows / math.sqrt(level_count), names


def build_phi():
    phi = None
    for count in LEVEL_COUNTS:
        rows, _ = axis_basis(count)
        phi = rows if phi is None else np.kron(phi, rows)
    return phi


def tensor_forward(R):
    """Coefficient transform applied one axis at a time.

    Each axis is a sum of at most three terms, so the roundoff is a few
    machine epsilons rather than the 144-term condensation of one dot product.
    """
    array = R.reshape(-1, *LEVEL_COUNTS)
    for axis, count in enumerate(LEVEL_COUNTS):
        matrix, _ = axis_basis(count)
        array = np.tensordot(matrix, array, axes=([1], [axis + 1]))
        array = np.moveaxis(array, 0, axis + 1)
    return array.reshape(-1, N_CELLS) / N_CELLS


def tensor_inverse(C):
    """Reconstruct r from its coefficients, one axis at a time."""
    array = C.reshape(-1, *LEVEL_COUNTS)
    for axis, count in enumerate(LEVEL_COUNTS):
        matrix, _ = axis_basis(count)
        array = np.tensordot(matrix.T, array, axes=([1], [axis + 1]))
        array = np.moveaxis(array, 0, axis + 1)
    return array.reshape(-1, N_CELLS)


def build_euclidean_axes():
    return [euclidean_rows(count)[0] for count in LEVEL_COUNTS]


def build_components():
    """Every component of the product basis with its support and channels."""
    rows_per_axis = [axis_basis(count)[1] for count in LEVEL_COUNTS]
    components = []
    for index in range(N_CELLS):
        choices = tuple(
            (index // stride) % count
            for count, stride in zip(LEVEL_COUNTS, STRIDES)
        )
        support = tuple(i for i, choice in enumerate(choices) if choice != 0)
        components.append(
            {
                "index": index,
                "choices": choices,
                "order": len(support),
                "support": support,
                "support_axes": [AXES[i] for i in support],
                "channels": [rows_per_axis[i][choices[i]] for i in support],
                "channels_1based": [choices[i] for i in support],
            }
        )
    return components


def cell_levels_table():
    table = np.zeros((N_CELLS, len(AXES)), dtype=np.int64)
    for cell in range(N_CELLS):
        for axis, (count, stride) in enumerate(zip(LEVEL_COUNTS, STRIDES)):
            table[cell, axis] = (cell // stride) % count
    return table


CELL_LEVELS = cell_levels_table()
COMPONENTS = build_components()


def component_selectors(components, order):
    return np.array([c["index"] for c in components if c["order"] == order])


# --- data -------------------------------------------------------------------


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def region_of(x):
    if x >= KFLATHI:
        return "far"
    if x >= KX0:
        return "B"
    return "A"


def load_route(path, route_name):
    """Read one cell file into R[lane, point, cell] of signed F0 ratios."""
    started = time.perf_counter()
    lane_index = {lane: i for i, lane in enumerate(LANES)}
    level_index = {
        axis: {name: i for i, name in enumerate(LEVELS[axis])} for axis in AXES
    }
    x_index = {}
    kind_by_x = {}
    lanes = []
    cells = []
    n_values = []
    x_values = []
    ratios = []
    row_count = 0
    skipped = 0
    zero_f0 = 0
    worst_abs_err_identity = 0.0
    worst_f0_identity = 0.0
    worst_ratio_identity = 0.0
    kind_conflicts = 0
    level_conflicts = 0

    with lzma.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header != COLUMNS:
            raise SystemExit(
                f"{path}: unexpected header\n  got {header}\n  want {COLUMNS}"
            )
        for row in reader:
            row_count += 1
            status = row[24]
            if status != "ok":
                skipped += 1
                continue
            lane = lane_index[row[0]]
            cell = 0
            for position, (axis, stride) in enumerate(zip(AXES, STRIDES)):
                choice = level_index[axis].get(row[1 + position])
                if choice is None:
                    level_conflicts += 1
                    break
                cell += choice * stride
            n = int(row[7])
            x = float(row[8])
            kind = row[23]
            previous_kind = kind_by_x.get(x)
            if previous_kind is None:
                kind_by_x[x] = kind
            elif previous_kind != kind:
                kind_conflicts += 1

            f0 = float(row[14])
            value = float(row[10])
            ref_hi = float(row[11])
            ref_lo = float(row[12])
            abs_err = float(row[13])
            signed = (value - ref_hi) - ref_lo

            worst_abs_err_identity = max(
                worst_abs_err_identity, abs(abs_err - abs(signed))
            )
            # The C++ driver sets f0 = ref.hi[0] + ref.lo[0] once per group.
            if n == 0:
                worst_f0_identity = max(
                    worst_f0_identity, abs(f0 - (ref_hi + ref_lo))
                )
            ratio = signed / f0 if f0 != 0.0 else 0.0
            if f0 == 0.0:
                zero_f0 += 1
            worst_ratio_identity = max(
                worst_ratio_identity, abs(abs(ratio) - abs_err / f0)
                if f0 != 0.0 else 0.0
            )

            previous_x = x_index.get(x)
            if previous_x is None:
                previous_x = len(x_index)
                x_index[x] = previous_x

            lanes.append(lane)
            cells.append(cell)
            n_values.append(n)
            x_values.append(x)
            ratios.append(ratio)

    xs = sorted(x_index)
    x_remap = {value: index for index, value in enumerate(xs)}
    n_sorted = sorted(set(n_values))
    n_remap = {old: new for new, old in enumerate(n_sorted)}
    point_count = len(n_sorted) * len(xs)

    point_ids = np.empty(len(lanes), dtype=np.int64)
    lane_ids = np.empty(len(lanes), dtype=np.int64)
    cell_ids = np.empty(len(lanes), dtype=np.int64)
    ratio_values = np.empty(len(lanes), dtype=np.float64)
    for i in range(len(lanes)):
        point_ids[i] = n_remap[n_values[i]] * len(xs) + x_remap[x_values[i]]
        lane_ids[i] = lanes[i]
        cell_ids[i] = cells[i]
        ratio_values[i] = ratios[i]

    R = np.full((len(LANES), point_count, N_CELLS), np.nan, dtype=np.float64)
    R[lane_ids, point_ids, cell_ids] = ratio_values

    far_points = [
        p
        for p in range(point_count)
        if region_of(xs[p % len(xs)]) == "far"
    ]
    far_max_spread = 0.0
    for lane in range(len(LANES)):
        for point in far_points:
            far_max_spread = max(
                far_max_spread, float(np.abs(R[lane, point] - R[lane, point, 0]).max())
            )

    keys = (lane_ids * point_count + point_ids) * N_CELLS + cell_ids
    unique_keys = np.unique(keys).size
    missing = int(np.isnan(R).sum())
    expected_rows = len(LANES) * point_count * N_CELLS

    checks = {
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": Path(path).stat().st_size,
        "rows": row_count,
        "skipped_rows": skipped,
        "zero_f0_rows": zero_f0,
        "lanes_present": sorted({LANES[i] for i in lanes}),
        "points_per_lane": point_count,
        "n_count": len(n_sorted),
        "x_count": len(xs),
        "cells_per_lane": N_CELLS,
        "unique_row_keys": int(unique_keys),
        "expected_rows": expected_rows,
        "missing_values": missing,
        "kind_conflicts": kind_conflicts,
        "level_conflicts": level_conflicts,
        "far_points": len(far_points),
        "far_max_spread_across_cells": far_max_spread,
        "max_abs_err_identity_residual": worst_abs_err_identity,
        "max_f0_identity_residual": worst_f0_identity,
        "max_ratio_identity_residual": worst_ratio_identity,
    }
    if (
        missing
        or unique_keys != row_count
        or row_count != expected_rows
        or skipped
        or kind_conflicts
        or level_conflicts
    ):
        raise SystemExit(f"{path}: incomplete or inconsistent record: {checks}")

    kinds = [kind_by_x[x] for x in xs]
    return {
        "route": route_name,
        "R": R,
        "xs": xs,
        "ns": n_sorted,
        "kinds": kinds,
        "checks": checks,
        "load_seconds": time.perf_counter() - started,
    }


# --- locations --------------------------------------------------------------


def point_location(run, point):
    nx = len(run["xs"])
    n_index, x_index = divmod(int(point), nx)
    n = run["ns"][n_index]
    x = run["xs"][x_index]
    return {
        "n": int(n),
        "x": float(x),
        "x_nominal": float(x),
        "region": region_of(x),
        "kind": run["kinds"][x_index],
    }


def component_location(component):
    return {
        "component_index": int(component["index"]),
        "order": int(component["order"]),
        "support": list(component["support_axes"]),
        "channels": list(component["channels"]),
        "channels_1based": [int(c) for c in component["channels_1based"]],
    }


def describe(run, component, point):
    out = component_location(component)
    out.update(point_location(run, point))
    return out


# --- full cube --------------------------------------------------------------


def full_cube_stats(run, C, Phi, components, R_lane):
    """Exact full-cube decomposition statistics for one lane."""
    point_count = C.shape[0]
    regions = np.array(
        [
            region_of(run["xs"][point % len(run["xs"])])
            for point in range(point_count)
        ]
    )
    region_masks = {region: regions == region for region in ("A", "B", "far")}
    nonfar = regions != "far"

    centered_ss = np.zeros(point_count)
    order_ss = {k: np.zeros(point_count) for k in range(len(AXES) + 1)}
    masks = {}
    for order in range(len(AXES) + 1):
        mask = component_selectors(components, order)
        masks[order] = mask
        if order >= 1:
            block_ss = np.square(C[:, mask]).sum(axis=1)
            order_ss[order] += block_ss
            centered_ss += block_ss

    safe_centered = np.where(centered_ss > 0.0, centered_ss, np.inf)
    shares = np.zeros_like(C)
    np.divide(np.square(C), safe_centered[:, None], out=shares)

    orders = {}
    for order in range(1, len(AXES) + 1):
        block = C[:, masks[order]]
        abs_block = np.abs(block)
        flat = int(np.argmax(abs_block))
        point, local = np.unravel_index(flat, abs_block.shape)
        component = components[int(masks[order][local])]
        max_component = float(abs_block[point, local])
        component_max = abs_block.max(axis=0)
        component_share_max = shares[:, masks[order]].max(axis=0)
        order_share = np.zeros(point_count)
        np.divide(
            order_ss[order], safe_centered, out=order_share, where=centered_ss > 0.0
        )
        max_order_point = int(np.argmax(order_share))

        component_max_nonfar = (
            np.abs(C[np.ix_(nonfar, masks[order])]).max(axis=0)
            if nonfar.any()
            else component_max
        )
        component_share_nonfar = shares[np.ix_(nonfar, masks[order])]
        order_share_nonfar = order_share[nonfar]
        max_order_point_nonfar = (
            int(np.where(nonfar)[0][int(np.argmax(order_share_nonfar))])
            if nonfar.any()
            else None
        )
        region_share = {}
        for region, selection in region_masks.items():
            if not selection.any():
                region_share[region] = None
                continue
            if region == "far":
                region_share[region] = {
                    "max_order_share": None,
                    "at": None,
                    "material_points_10pct": 0,
                    "note": (
                        "mathematical zero: at x = 32 the response is "
                        "bit-identical across all 144 cells, so every contrast "
                        "is zero and a share is not defined.  The residual "
                        "coefficients are summation roundoff."
                    ),
                }
                continue
            values = order_share[selection]
            region_point = int(np.where(selection)[0][int(np.argmax(values))])
            region_share[region] = {
                "max_order_share": float(values.max()),
                "at": point_location(run, region_point),
                "material_points_10pct": int((values > ALPHA).sum()),
            }

        material_components = component_share_max > ALPHA
        readings = int((shares[:, masks[order]] > ALPHA).sum())
        orders[str(order)] = {
            "components": int(len(masks[order])),
            "max_abs_c": max_component,
            "max_at": describe(run, component, point),
            "local_max_abs_r": float(np.abs(R_lane[point, :]).max()),
            "component_share_at_max": float(shares[point, component["index"]]),
            "nonzero_above_floor": int((component_max > MEASUREMENT_FLOOR).sum()),
            "nonzero_above_floor_nonfar": int(
                (component_max_nonfar > MEASUREMENT_FLOOR).sum()
            ),
            "exactly_zero": int((component_max == 0.0).sum()),
            "exactly_zero_nonfar": int((component_max_nonfar == 0.0).sum()),
            "subfloor_nonzero": int(
                ((component_max > 0.0) & (component_max <= MEASUREMENT_FLOOR)).sum()
            ),
            "material_components_10pct": int(material_components.sum()),
            "material_components_10pct_nonfar": int(
                (component_share_nonfar.max(axis=0) > ALPHA).sum()
            ),
            "material_readings_10pct": readings,
            "material_readings_10pct_nonfar": int(
                (component_share_nonfar > ALPHA).sum()
            ),
            "readings_total": int(point_count * len(masks[order])),
            "readings_total_nonfar": int(
                int(nonfar.sum()) * len(masks[order])
            ),
            "max_order_share": float(order_share[max_order_point]),
            "max_order_share_at": point_location(run, max_order_point),
            "max_order_share_nonfar": (
                float(order_share_nonfar.max()) if nonfar.any() else None
            ),
            "max_order_share_nonfar_at": (
                point_location(run, max_order_point_nonfar)
                if max_order_point_nonfar is not None
                else None
            ),
            "max_order_share_by_region": region_share,
            "material_points_10pct": int((order_share > ALPHA).sum()),
            "material_points_nonfar_10pct": int((order_share_nonfar > ALPHA).sum()),
            "global_ss_share": (
                float(order_ss[order].sum() / centered_ss.sum())
                if centered_ss.sum() > 0.0
                else None
            ),
            "global_ss_share_nonfar": (
                float(order_ss[order][nonfar].sum() / centered_ss[nonfar].sum())
                if centered_ss[nonfar].sum() > 0.0
                else None
            ),
        }

    # Order-0 (mean) coefficient is reported for completeness.
    mean_c = C[:, 0]
    mean_flat = int(np.argmax(np.abs(mean_c)))
    orders["0"] = {
        "components": 1,
        "max_abs_c": float(abs(mean_c[mean_flat])),
        "max_at": point_location(run, mean_flat),
        "global_mean": float(mean_c.mean()),
    }

    # Non-additive part: every component of order two and above.
    ge2 = np.concatenate([masks[k] for k in range(2, len(AXES) + 1)])
    residual = C[:, ge2] @ Phi[ge2, :]
    abs_residual = np.abs(residual)
    flat = int(np.argmax(abs_residual))
    point, cell = np.unravel_index(flat, abs_residual.shape)
    order_ge2_share = np.zeros(point_count)
    order_ge2_ss = sum(order_ss[k] for k in range(2, len(AXES) + 1))
    np.divide(
        order_ge2_ss, safe_centered, out=order_ge2_share, where=centered_ss > 0.0
    )
    max_share_point = int(np.argmax(order_ge2_share))
    max_share_point_nonfar = (
        int(np.where(nonfar)[0][int(np.argmax(order_ge2_share[nonfar]))])
        if nonfar.any()
        else None
    )

    non_additive = {
        "max_abs_residual": float(abs_residual[point, cell]),
        "at": {
            **point_location(run, point),
            "cell_index": int(cell),
            "cell_levels": [int(v) for v in CELL_LEVELS[cell]],
        },
        "r_at": float(R_lane[point, cell]),
        "order_ge2_share_at": float(order_ge2_share[point]),
        "max_order_ge2_share": float(order_ge2_share[max_share_point]),
        "max_order_ge2_share_at": point_location(run, max_share_point),
        "max_order_ge2_share_nonfar": (
            float(order_ge2_share[nonfar].max()) if nonfar.any() else None
        ),
        "max_order_ge2_share_nonfar_at": (
            point_location(run, max_share_point_nonfar)
            if max_share_point_nonfar is not None
            else None
        ),
        "material_points_10pct": int((order_ge2_share > ALPHA).sum()),
        "material_points_nonfar_10pct": int(
            (order_ge2_share[nonfar] > ALPHA).sum()
        ),
        "global_ss_share_ge2": (
            float(order_ge2_ss.sum() / centered_ss.sum())
            if centered_ss.sum() > 0.0
            else None
        ),
        "global_ss_share_ge2_nonfar": (
            float(order_ge2_ss[nonfar].sum() / centered_ss[nonfar].sum())
            if centered_ss[nonfar].sum() > 0.0
            else None
        ),
    }

    # Per-component table for the CSV output.
    component_rows = []
    for component in components:
        index = component["index"]
        column = C[:, index]
        abs_column = np.abs(column)
        point = int(np.argmax(abs_column))
        column_nonfar = np.abs(C[np.ix_(np.where(nonfar)[0], [index])])[:, 0]
        share_nonfar = shares[np.ix_(np.where(nonfar)[0], [index])][:, 0]
        component_rows.append(
            {
                "component_index": int(index),
                "order": int(component["order"]),
                "support": "+".join(component["support_axes"]) or "(mean)",
                "channels": "+".join(component["channels"]) or "(mean)",
                "channels_1based": ":".join(
                    str(c) for c in component["channels_1based"]
                )
                or "-",
                "max_abs_c": float(abs_column[point]),
                "max_abs_c_nonfar": float(column_nonfar.max()),
                "at_n": int(point_location(run, point)["n"]),
                "at_x": float(point_location(run, point)["x"]),
                "region": point_location(run, point)["region"],
                "kind": point_location(run, point)["kind"],
                "component_share_max": float(shares[:, index].max()),
                "component_share_max_nonfar": float(share_nonfar.max()),
                "material_10pct": bool(shares[:, index].max() > ALPHA),
                "material_10pct_nonfar": bool(share_nonfar.max() > ALPHA),
                "nonzero_above_floor": bool(
                    abs_column.max() > MEASUREMENT_FLOOR
                ),
                "nonzero_above_floor_nonfar": bool(
                    column_nonfar.max() > MEASUREMENT_FLOOR
                ),
                "exactly_zero": bool(abs_column.max() == 0.0),
                "exactly_zero_nonfar": bool(column_nonfar.max() == 0.0),
            }
        )

    return {
        "points": int(point_count),
        "orders": orders,
        "non_additive": non_additive,
        "component_rows": component_rows,
    }


def region_maxima(run, C, components, order):
    mask = component_selectors(components, order)
    out = {}
    for region in ("A", "B", "far"):
        points = [
            p
            for p in range(C.shape[0])
            if point_location(run, p)["region"] == region
        ]
        if not points:
            out[region] = None
            continue
        block = np.abs(C[np.ix_(points, mask)])
        flat = int(np.argmax(block))
        point_local, comp_local = np.unravel_index(flat, block.shape)
        point = points[point_local]
        component = components[int(mask[comp_local])]
        out[region] = {
            "max_abs_c": float(block[point_local, comp_local]),
            "at": {
                **point_location(run, point),
                "support": list(component["support_axes"]),
                "channels": list(component["channels"]),
            },
        }
    return out


def order3_stability(C_lane, run, components):
    """Order-3 share summary used to compare transform implementations.

    Near-zero points make an order share ill-conditioned; the tensor transform
    (three-term sums per axis) is more accurate than one 144-term dot product,
    so comparing the two exposes points whose share is roundoff, not signal.
    """
    point_count = C_lane.shape[0]
    regions = np.array(
        [
            region_of(run["xs"][point % len(run["xs"])])
            for point in range(point_count)
        ]
    )
    nonfar = regions != "far"
    centered = np.square(C_lane[:, 1:]).sum(axis=1)
    safe = np.where(centered > 0.0, centered, np.inf)
    mask = component_selectors(components, 3)
    order_ss = np.square(C_lane[:, mask]).sum(axis=1)
    share = np.divide(
        order_ss, safe, out=np.zeros(point_count), where=centered > 0.0
    )
    share_ab = share[nonfar]
    component_share_ab = (np.square(C_lane[:, mask]) / safe[:, None])[nonfar]
    point = int(np.where(nonfar)[0][int(np.argmax(share_ab))])
    return {
        "max_order3_share_AB": float(share_ab.max()),
        "max_order3_share_AB_at": point_location(run, point),
        "order3_material_components_AB": int(
            (component_share_ab.max(axis=0) > ALPHA).sum()
        ),
        "order3_material_readings_AB": int((component_share_ab > ALPHA).sum()),
        "order3_material_points_AB": int((share_ab > ALPHA).sum()),
    }


# --- d79 corner blocks ------------------------------------------------------


def d79_tables(R_lane, Phi, euclidean_axes, components):
    """Note-convention P and task-convention c for every d79 block."""
    point_count = R_lane.shape[0]
    P = np.zeros((point_count, N_CELLS), dtype=np.float64)
    c_block = np.zeros((point_count, N_CELLS), dtype=np.float64)

    for component in components:
        block = component_block_cells(component)
        size = len(block)
        phi_block = Phi[component["index"], block]
        c_values = (R_lane[:, block] @ phi_block) / size
        c_block[:, component["index"]] = c_values
        P[:, component["index"]] = c_values * math.sqrt(size)
    return P, c_block


def component_block_cells(component):
    support = set(component["support"])
    return [
        cell
        for cell in range(N_CELLS)
        if all(
            CELL_LEVELS[cell, axis] == 0
            for axis in range(len(AXES))
            if axis not in support
        )
    ]


def d79_stats(run, C, R_lane, P, c_block, components, Phi, euclidean_axes):
    """Statistics of the d79 corner blocks, with the note's P convention."""
    point_count = C.shape[0]
    masks = {
        order: component_selectors(components, order)
        for order in range(1, len(AXES) + 1)
    }
    d79_cells = sorted(
        {
            cell
            for component in components
            if component["order"] <= 3
            for cell in component_block_cells(component)
        }
    )

    # The note's materiality: a contrast relative to the largest response on
    # its own block at the same point.  This reproduces the note's max relative
    # and its 10% material counts exactly.
    block_scale = np.zeros((point_count, N_CELLS), dtype=np.float64)
    for component in components:
        block = component_block_cells(component)
        block_scale[:, component["index"]] = np.abs(R_lane[:, block]).max(axis=1)
    safe_scale = np.where(block_scale > 0.0, block_scale, np.inf)
    relative = np.abs(P) / safe_scale

    # The same measure applied to the pure full-cube component: its self term
    # in the note's P units is c_alpha sqrt(B).
    sqrt_block = np.array(
        [math.sqrt(len(component_block_cells(c))) for c in components]
    )
    relative_pure = (np.abs(C) * sqrt_block[None, :]) / safe_scale

    def alias_split(point, target):
        """Decompose a d79 block functional into its pure and alias orders.

        The full-cube coefficients are in the task (mean-square) basis.  The
        note's Euclidean basis coefficient is 12 times larger, hence the 12.
        The self term is c_alpha sqrt(B); higher orders are the triangular
        aliases whose support strictly contains the target's.
        """
        target_support = set(target["support"])
        block_size = len(component_block_cells(target))
        self_term = float(C[point, target["index"]] * math.sqrt(block_size))
        by_order = {}
        for beta in components:
            if beta["order"] < 1:
                continue
            support = set(beta["support"])
            if not target_support.issubset(support):
                continue
            if not all(
                beta["choices"][axis] == target["choices"][axis]
                for axis in target_support
            ):
                continue
            weight = 12.0
            for axis in support - target_support:
                weight *= euclidean_axes[axis][beta["choices"][axis]][0]
            for axis in range(len(AXES)):
                if axis not in support:
                    weight /= math.sqrt(LEVEL_COUNTS[axis])
            by_order[beta["order"]] = by_order.get(beta["order"], 0.0) + float(
                C[point, beta["index"]] * weight
            )
        return self_term, by_order

    orders = {}
    for order in range(1, 4):
        mask = masks[order]
        block = P[:, mask]
        abs_block = np.abs(block)
        flat = int(np.argmax(abs_block))
        point, local = np.unravel_index(flat, abs_block.shape)
        component = components[int(mask[local])]
        component_max = abs_block.max(axis=0)

        rel_block = relative[:, mask]
        rel_flat = int(np.argmax(rel_block))
        rel_point, rel_local = np.unravel_index(rel_flat, rel_block.shape)
        rel_component = components[int(mask[rel_local])]

        rel_pure = relative_pure[:, mask]
        pure_flat = int(np.argmax(rel_pure))
        pure_point, pure_local = np.unravel_index(pure_flat, rel_pure.shape)
        pure_component = components[int(mask[pure_local])]

        self_term, alias_by_order = alias_split(point, component)
        direct = float(block[point, local])

        orders[str(order)] = {
            "components": int(len(mask)),
            "max_abs_P_note": abs(direct),
            "P_signed_at_max": direct,
            "max_abs_P_at": describe(run, component, point),
            "local_max_abs_r": float(np.abs(R_lane[point, :]).max()),
            "block_scale_at_max": float(block_scale[point, component["index"]]),
            "relative_at_max": float(relative[point, component["index"]]),
            "max_abs_c_block": float(np.abs(c_block[:, mask]).max()),
            "max_abs_c_full": float(np.abs(C[:, mask]).max()),
            "exactly_zero_P": int((component_max == 0.0).sum()),
            "at_or_below_floor_P": int(
                (component_max <= MEASUREMENT_FLOOR).sum()
            ),
            "nonzero_above_floor_P": int(
                (component_max > MEASUREMENT_FLOOR).sum()
            ),
            "note_relative": {
                "max_relative": float(rel_block[rel_point, rel_local]),
                "max_relative_at": describe(run, rel_component, rel_point),
                "material_components_10pct": int(
                    (rel_block.max(axis=0) > ALPHA).sum()
                ),
                "material_readings_10pct": int((rel_block > ALPHA).sum()),
                "material_points_10pct": int(
                    (rel_block.max(axis=1) > ALPHA).sum()
                ),
            },
            "pure_full_cube_relative": {
                "max_relative": float(rel_pure[pure_point, pure_local]),
                "max_relative_at": describe(run, pure_component, pure_point),
                "material_components_10pct": int(
                    (rel_pure.max(axis=0) > ALPHA).sum()
                ),
                "material_readings_10pct": int((rel_pure > ALPHA).sum()),
                "material_points_10pct": int(
                    (rel_pure.max(axis=1) > ALPHA).sum()
                ),
            },
            "alias_split": {
                "P_direct": direct,
                "self_term": self_term,
                "alias_by_order": {
                    str(k): float(v)
                    for k, v in alias_by_order.items()
                    if k != order
                },
                "expansion_sum": float(sum(alias_by_order.values())),
                "self_fraction": (
                    self_term / direct if direct != 0.0 else None
                ),
            },
        }

        # Alternative share definitions, kept for transparency.
        def share_stats(values, total_mask):
            total = np.square(values[:, total_mask]).sum(axis=1)
            safe = np.where(total > 0.0, total, np.inf)
            sub = np.square(values[:, mask])
            per_point = sub.sum(axis=1) / safe
            per_component_max = (sub / safe[:, None]).max(axis=0)
            return {
                "material_components_10pct": int(
                    (per_component_max > ALPHA).sum()
                ),
                "material_readings_10pct": int(
                    ((sub / safe[:, None]) > ALPHA).sum()
                ),
                "max_order_share": float(per_point.max()),
                "material_points_10pct": int((per_point > ALPHA).sum()),
            }

        orders[str(order)]["shares"] = {
            "P": share_stats(P, np.arange(N_CELLS)[1:]),
            "c_block": share_stats(c_block, np.arange(N_CELLS)[1:]),
            "c_full": share_stats(C, np.arange(N_CELLS)[1:]),
        }

    # Region maxima of the note's P for orders 1-3.
    regions = {}
    for order in range(1, 4):
        mask = masks[order]
        regions[str(order)] = {}
        for region in ("A", "B", "far"):
            points = [
                p
                for p in range(point_count)
                if point_location(run, p)["region"] == region
            ]
            if not points:
                regions[str(order)][region] = None
                continue
            block = np.abs(P[np.ix_(points, mask)])
            flat = int(np.argmax(block))
            point_local, comp_local = np.unravel_index(flat, block.shape)
            point = points[point_local]
            component = components[int(mask[comp_local])]
            regions[str(order)][region] = {
                "max_abs_P_note": float(block[point_local, comp_local]),
                "at": {
                    **point_location(run, point),
                    "support": list(component["support_axes"]),
                    "channels": list(component["channels"]),
                },
            }

    # The d79 corner model: mean + orders 1-3, each component evaluated on the
    # d79 cells that contain its block.  Orders above 3 are not estimated.
    estimated = set()
    for order in range(0, 4):
        estimated.update(int(i) for i in component_selectors(components, order))
    model_matrix = np.zeros((len(d79_cells), N_CELLS), dtype=np.float64)
    nonadd_matrix = np.zeros_like(model_matrix)
    for row, cell in enumerate(d79_cells):
        nonbase = {
            axis
            for axis in range(len(AXES))
            if CELL_LEVELS[cell, axis] != 0
        }
        for component in components:
            if component["order"] > 3:
                continue
            if set(component["support"]) <= nonbase:
                value = Phi[component["index"], cell]
                model_matrix[row, component["index"]] = value
                if component["order"] >= 2:
                    nonadd_matrix[row, component["index"]] = value

    model_values = c_block[:, :] @ model_matrix.T
    nonadd_values = c_block[:, :] @ nonadd_matrix.T
    reconstruction = model_values - R_lane[:, d79_cells]
    flat = int(np.argmax(np.abs(nonadd_values)))
    point, row = np.unravel_index(flat, nonadd_values.shape)

    corner_model = {
        "d79_cells": len(d79_cells),
        "estimated_components": len(estimated),
        "max_abs_nonadd": float(np.abs(nonadd_values)[point, row]),
        "max_abs_nonadd_at": {
            **point_location(run, point),
            "cell_index": int(d79_cells[row]),
            "cell_levels": [int(v) for v in CELL_LEVELS[d79_cells[row]]],
        },
        "max_abs_reconstruction_error": float(np.abs(reconstruction).max()),
        "max_abs_r_on_d79": float(np.abs(R_lane[:, d79_cells]).max()),
    }

    return {
        "blocks": [component_block_cells(c) for c in components],
        "d79_cells": d79_cells,
        "orders": orders,
        "regions": regions,
        "corner_model": corner_model,
        "P": P,
        "c_block": c_block,
    }


# --- precision --------------------------------------------------------------


def precision_decomposition(C, components):
    """Decompose every call-axis component across the four lanes."""
    psi = np.array(
        [
            [1.0, 1.0, 1.0, 1.0],
            [-3.0, -1.0, 1.0, 3.0],
            [1.0, -1.0, -1.0, 1.0],
            [-1.0, 3.0, -3.0, 1.0],
        ],
        dtype=np.float64,
    )
    psi[1] /= math.sqrt(5.0)
    psi[2] /= 2.0
    psi[3] /= math.sqrt(5.0)
    gram = (psi @ psi.T) / 4.0
    check = float(np.abs(gram - np.eye(4)).max())
    contrasts = np.einsum("jl,lpc->jpc", psi, C) / 4.0
    return {
        "rows": psi,
        "row_names": ["mean", "linear", "quadratic", "cubic"],
        "orthonormality_residual": check,
        "Q": contrasts,
    }


def precision_stats(run, C, decomposition, components):
    Q = decomposition["Q"]
    point_count = C.shape[0]
    names = decomposition["row_names"]

    def best(component_mask, label):
        block = np.abs(Q[1:, :, :][:, :, component_mask])
        flat = int(np.argmax(block))
        contrast, point, local = np.unravel_index(flat, block.shape)
        component = components[int(np.array(component_mask)[local])]
        return {
            "kind": label,
            "coefficient": float(block[contrast, point, local]),
            "precision_contrast": names[contrast + 1],
            "precision_row": int(contrast + 1),
            "at": describe(run, component, point),
            "lane_values": [float(v) for v in C[:, point, component["index"]]],
        }

    order_masks = {
        order: component_selectors(components, order)
        for order in range(1, len(AXES) + 1)
    }
    marginal = np.array([0])
    result = {
        "precision_main_effect": best(marginal, "precision main effect"),
        "precision_by_call": best(
            np.concatenate(list(order_masks.values())), "precision-by-call"
        ),
        "precision_by_call_order": {},
    }
    for order, mask in order_masks.items():
        result["precision_by_call_order"][str(order)] = best(
            mask, f"precision-by-call order {order}"
        )

    # Per-component summary: the largest precision contrast of each call axis
    # component, which is the precision interaction of that component.
    rows = []
    for component in components:
        column = Q[1:, :, component["index"]]
        flat = int(np.argmax(np.abs(column)))
        contrast, point = np.unravel_index(flat, column.shape)
        rows.append(
            {
                "component_index": int(component["index"]),
                "order": int(component["order"]),
                "support": "+".join(component["support_axes"]) or "(mean)",
                "channels": "+".join(component["channels"]) or "(mean)",
                "max_abs_precision_contrast": float(
                    abs(column[contrast, point])
                ),
                "precision_contrast": names[contrast + 1],
                "at_n": int(point_location(run, point)["n"]),
                "at_x": float(point_location(run, point)["x"]),
                "region": point_location(run, point)["region"],
                "lane_values": [float(v) for v in C[:, point, component["index"]]],
            }
        )
    result["component_rows"] = rows
    return result


# --- output -----------------------------------------------------------------


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def write_csv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def fmt(value):
    if isinstance(value, float):
        return f"{value:.17g}"
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fused", required=True, help="msvc-fused cell file")
    parser.add_argument("--separate", required=True, help="msvc-separate cell file")
    parser.add_argument("--out-dir", default="results")
    args = parser.parse_args()

    started = time.perf_counter()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    Phi = build_phi()
    euclidean_axes = build_euclidean_axes()
    components = COMPONENTS

    # Orthonormality checks.
    per_axis = {}
    for axis, count in zip(AXES, LEVEL_COUNTS):
        rows, names = axis_basis(count)
        gram = (rows @ rows.T) / count
        per_axis[axis] = {
            "levels": count,
            "rows": names,
            "max_abs_gram_minus_identity": float(
                np.abs(gram - np.eye(count)).max()
            ),
        }
    phi_gram = (Phi @ Phi.T) / N_CELLS
    orthonormality = {
        "per_axis": per_axis,
        "phi_max_abs_gram_minus_identity": float(
            np.abs(phi_gram - np.eye(N_CELLS)).max()
        ),
        "measure": "uniform over levels; mean square 1, zero mean, "
        "mutually orthogonal",
    }

    routes = {}
    for route_name, path in (("fused", args.fused), ("separate", args.separate)):
        print(f"loading {route_name}: {path}", flush=True)
        run = load_route(path, route_name)
        print(
            f"  {run['checks']['rows']} rows, {run['checks']['points_per_lane']} "
            f"points, {run['checks']['cells_per_lane']} cells, "
            f"{run['load_seconds']:.1f} s",
            flush=True,
        )

        C_flat = tensor_forward(run["R"])
        C = C_flat.reshape(len(LANES), run["R"].shape[1], N_CELLS)
        R_hat_flat = tensor_inverse(C_flat)
        R_hat = R_hat_flat.reshape(len(LANES), run["R"].shape[1], N_CELLS)
        direct_matmul = (
            run["R"].reshape(-1, N_CELLS) @ Phi.T
        ) / N_CELLS
        transform_agreement = {
            "max_abs_difference_tensor_vs_matmul": float(
                np.abs(C_flat - direct_matmul).max()
            ),
            "max_abs_coefficient": float(np.abs(C_flat).max()),
            "max_relative_difference": float(
                np.abs(C_flat - direct_matmul).max() / np.abs(C_flat).max()
            ),
            "order3_share_stability": {},
        }
        direct_C = direct_matmul.reshape(
            len(LANES), run["R"].shape[1], N_CELLS
        )
        for lane_name in RESOLVING_LANES:
            lane_index = LANES.index(lane_name)
            transform_agreement["order3_share_stability"][lane_name] = {
                "tensor": order3_stability(C[lane_index], run, components),
                "direct_matmul": order3_stability(
                    direct_C[lane_index], run, components
                ),
            }
        reconstruction = np.abs(R_hat - run["R"])

        # Independent check of two coefficients: group means do not use the
        # transform code at all.  At n = 0, x = 12:
        #   partition linear  c = sqrt(3/2) (m0 - m2) / 3
        #   route contrast    c = (m0 - m1) / 2
        point_12 = run["ns"].index(0) * len(run["xs"]) + run["xs"].index(12.0)
        partition_linear = next(
            component
            for component in components
            if component["support_axes"] == ["partition"]
            and component["channels"] == ["linear"]
        )
        route_contrast = next(
            component
            for component in components
            if component["support_axes"] == ["route"]
        )
        independent = {}
        for lane, lane_name in enumerate(LANES):
            values = run["R"][lane, point_12]
            partition_means = [
                values[CELL_LEVELS[:, 2] == level].mean() for level in range(3)
            ]
            route_means = [
                values[CELL_LEVELS[:, 0] == level].mean() for level in range(2)
            ]
            partition_direct = (
                math.sqrt(1.5) * (partition_means[0] - partition_means[2]) / 3.0
            )
            route_direct = (route_means[0] - route_means[1]) / 2.0
            partition_transform = float(
                C[lane, point_12, partition_linear["index"]]
            )
            route_transform = float(C[lane, point_12, route_contrast["index"]])
            independent[lane_name] = {
                "partition_linear_group_mean": partition_direct,
                "partition_linear_transform": partition_transform,
                "route_contrast_group_mean": route_direct,
                "route_contrast_transform": route_transform,
                "max_abs_difference": max(
                    abs(partition_direct - partition_transform),
                    abs(route_direct - route_transform),
                ),
            }
        recon = {}
        for lane, lane_name in enumerate(LANES):
            recon[lane_name] = {
                "max_abs_reconstruction_error": float(
                    reconstruction[lane].max()
                ),
                "max_abs_r": float(np.abs(run["R"][lane]).max()),
                "max_abs_reconstruction_error_over_max_abs_r": float(
                    reconstruction[lane].max() / np.abs(run["R"][lane]).max()
                ),
                "theoretical_scale_n_eps_max_r": float(
                    N_CELLS
                    * np.finfo(np.float64).eps
                    * np.abs(run["R"][lane]).max()
                ),
            }

        full = {}
        d79 = {}
        precision = {}
        for lane, lane_name in enumerate(LANES):
            Cl = C[lane]
            R_lane = run["R"][lane]
            full[lane_name] = full_cube_stats(
                run, Cl, Phi, components, R_lane
            )
            full[lane_name]["region_maxima"] = {
                str(order): region_maxima(run, Cl, components, order)
                for order in range(1, len(AXES) + 1)
            }
            P, c_block = d79_tables(
                R_lane, Phi, euclidean_axes, components
            )
            d79[lane_name] = d79_stats(
                run, Cl, R_lane, P, c_block, components, Phi, euclidean_axes
            )
            # Keep only the summary in memory; the tables are large.
            d79[lane_name].pop("P")
            d79[lane_name].pop("c_block")
            d79[lane_name].pop("blocks")
            d79[lane_name].pop("d79_cells")

        decomposition = precision_decomposition(C, components)
        precision = precision_stats(run, C, decomposition, components)
        precision["orthonormality_residual"] = decomposition[
            "orthonormality_residual"
        ]
        precision["rows"] = decomposition["row_names"]
        precision["fully_crossed"] = {
            "lanes": list(LANES),
            "lanes_present": run["checks"]["lanes_present"],
            "cells_per_lane": run["checks"]["cells_per_lane"],
            "points_per_lane": run["checks"]["points_per_lane"],
            "all_rows_ok": run["checks"]["skipped_rows"] == 0,
            "statement": "every lane is measured at every cell and every point",
        }

        routes[route_name] = {
            "checks": run["checks"],
            "xs": run["xs"],
            "ns": run["ns"],
            "kinds": run["kinds"],
            "transform_agreement": transform_agreement,
            "independent_coefficient_check": independent,
            "reconstruction": recon,
            "full_cube": full,
            "d79": d79,
            "precision": precision,
            "C": C,
        }

    # --- comparisons against the note ---------------------------------------
    # The note prints two or three significant figures; a difference within
    # 0.5 per cent is agreement at the note's own precision.
    note_tolerance = 5e-3
    comparisons = {}
    for route_name, route in routes.items():
        comparisons[route_name] = {}
        for lane_name in RESOLVING_LANES:
            key = (route_name, lane_name)
            d79_lane = route["d79"][lane_name]
            entry = {}
            for name, field, note_key in (
                ("order1", "1", "order1_max"),
                ("order2", "2", "order2_max"),
                ("order3", "3", "order3_max"),
            ):
                note_value = NOTE[note_key].get(key)
                ours = d79_lane["orders"][field]["max_abs_P_note"]
                entry[name] = {
                    "ours_P": ours,
                    "note": note_value,
                    "ratio": (ours / note_value) if note_value else None,
                    "agrees_at_note_precision": (
                        abs(ours - note_value) <= note_tolerance * abs(note_value)
                        if note_value
                        else None
                    ),
                }
            note_relative = NOTE["order3_max_relative"].get(key)
            ours_relative = d79_lane["orders"]["3"]["note_relative"]["max_relative"]
            entry["order3_max_relative"] = {
                "ours": ours_relative,
                "note": note_relative,
                "agrees_at_note_precision": (
                    abs(ours_relative - note_relative)
                    <= note_tolerance * abs(note_relative)
                    if note_relative
                    else None
                ),
            }
            note_zero = NOTE["order3_exact_zero_components"].get(key)
            entry["order3_exact_zero_components"] = {
                "ours_bit_exact": d79_lane["orders"]["3"]["exactly_zero_P"],
                "ours_at_or_below_floor": d79_lane["orders"]["3"][
                    "at_or_below_floor_P"
                ],
                "note": note_zero,
                "note_agrees_with_bit_exact": (
                    d79_lane["orders"]["3"]["exactly_zero_P"] == note_zero
                ),
                "note_agrees_with_at_or_below_floor": (
                    d79_lane["orders"]["3"]["at_or_below_floor_P"] == note_zero
                ),
                "status": (
                    "reproduced"
                    if d79_lane["orders"]["3"]["exactly_zero_P"] == note_zero
                    or d79_lane["orders"]["3"]["at_or_below_floor_P"] == note_zero
                    else "not reproduced; the note's zero convention was not "
                    "identified, so both counts are reported"
                ),
            }
            note_material = NOTE["order3_material_components"].get(key)
            note_readings = NOTE["order3_material_readings"].get(key)
            material = d79_lane["orders"]["3"]["note_relative"]
            entry["order3_material_components"] = {
                "ours_note_relative": material["material_components_10pct"],
                "note": note_material,
                "agrees": material["material_components_10pct"] == note_material,
            }
            entry["order3_material_readings"] = {
                "ours_note_relative": material["material_readings_10pct"],
                "note": note_readings,
                "agrees": material["material_readings_10pct"] == note_readings,
            }
            pure = d79_lane["orders"]["3"]["pure_full_cube_relative"]
            entry["order3_pure_full_cube_under_note_measure"] = {
                "max_relative": pure["max_relative"],
                "material_components_10pct": pure[
                    "material_components_10pct"
                ],
                "material_readings_10pct": pure["material_readings_10pct"],
                "note_block_contrast_components": note_material,
                "note_block_contrast_readings": note_readings,
            }
            entry["order3_alias_split"] = d79_lane["orders"]["3"]["alias_split"]
            comparisons[route_name][lane_name] = entry

    # --- write outputs -------------------------------------------------------
    full_rows = []
    for route_name, route in routes.items():
        for lane_name, lane in route["full_cube"].items():
            for row in lane["component_rows"]:
                full_rows.append(
                    [
                        route_name,
                        lane_name,
                        row["component_index"],
                        row["order"],
                        row["support"],
                        row["channels"],
                        row["channels_1based"],
                        fmt(row["max_abs_c"]),
                        fmt(row["max_abs_c_nonfar"]),
                        row["at_n"],
                        fmt(row["at_x"]),
                        row["region"],
                        row["kind"],
                        fmt(row["component_share_max"]),
                        fmt(row["component_share_max_nonfar"]),
                        int(row["material_10pct"]),
                        int(row["material_10pct_nonfar"]),
                        int(row["nonzero_above_floor"]),
                        int(row["nonzero_above_floor_nonfar"]),
                        int(row["exactly_zero"]),
                        int(row["exactly_zero_nonfar"]),
                    ]
                )
    write_csv(
        out_dir / "full_cube_components.csv",
        [
            "route", "lane", "component_index", "order", "support", "channels",
            "channels_1based", "max_abs_c", "max_abs_c_nonfar", "at_n",
            "at_x", "region", "kind", "component_share_max",
            "component_share_max_nonfar", "material_10pct",
            "material_10pct_nonfar", "nonzero_above_floor",
            "nonzero_above_floor_nonfar", "exactly_zero",
            "exactly_zero_nonfar",
        ],
        full_rows,
    )

    order_rows = []
    for route_name, route in routes.items():
        for lane_name, lane in route["full_cube"].items():
            for order in ("1", "2", "3", "4", "5", "6"):
                stat = lane["orders"][order]
                order_rows.append(
                    [
                        route_name, lane_name, int(order),
                        fmt(stat["max_abs_c"]),
                        "+".join(stat["max_at"]["support"]),
                        "+".join(stat["max_at"]["channels"]),
                        stat["max_at"]["n"], fmt(stat["max_at"]["x"]),
                        stat["max_at"]["region"], stat["max_at"]["kind"],
                        stat["components"],
                        stat["nonzero_above_floor"],
                        stat["exactly_zero"],
                        stat["subfloor_nonzero"],
                        stat["material_components_10pct"],
                        stat["material_readings_10pct"],
                        stat["readings_total"],
                        fmt(stat["max_order_share"]),
                        stat["max_order_share_at"]["n"],
                        fmt(stat["max_order_share_at"]["x"]),
                        fmt(stat["global_ss_share"])
                        if stat["global_ss_share"] is not None else "",
                        stat["nonzero_above_floor_nonfar"],
                        stat["exactly_zero_nonfar"],
                        stat["material_components_10pct_nonfar"],
                        stat["material_readings_10pct_nonfar"],
                        stat["readings_total_nonfar"],
                        fmt(stat["max_order_share_nonfar"])
                        if stat["max_order_share_nonfar"] is not None else "",
                        stat["max_order_share_nonfar_at"]["n"]
                        if stat["max_order_share_nonfar_at"] is not None else "",
                        fmt(stat["max_order_share_nonfar_at"]["x"])
                        if stat["max_order_share_nonfar_at"] is not None else "",
                        stat["material_points_nonfar_10pct"],
                        fmt(stat["global_ss_share_nonfar"])
                        if stat["global_ss_share_nonfar"] is not None else "",
                    ]
                )
    write_csv(
        out_dir / "full_cube_orders.csv",
        [
            "route", "lane", "order", "max_abs_c", "support", "channels",
            "at_n", "at_x", "region", "kind", "components",
            "nonzero_above_floor", "exactly_zero", "subfloor_nonzero",
            "material_components_10pct", "material_readings_10pct",
            "readings_total", "max_order_share", "max_order_share_at_n",
            "max_order_share_at_x", "global_ss_share",
            "nonzero_above_floor_nonfar", "exactly_zero_nonfar",
            "material_components_10pct_nonfar",
            "material_readings_10pct_nonfar", "readings_total_nonfar",
            "max_order_share_nonfar", "max_order_share_nonfar_at_n",
            "max_order_share_nonfar_at_x", "material_points_nonfar_10pct",
            "global_ss_share_nonfar",
        ],
        order_rows,
    )

    d79_rows = []
    for route_name, route in routes.items():
        for lane_name, lane in route["d79"].items():
            for order in ("1", "2", "3"):
                stat = lane["orders"][order]
                note_rel = stat["note_relative"]
                pure_rel = stat["pure_full_cube_relative"]
                alias = stat["alias_split"]
                d79_rows.append(
                    [
                        route_name, lane_name, int(order),
                        fmt(stat["max_abs_P_note"]),
                        "+".join(stat["max_abs_P_at"]["support"]),
                        "+".join(stat["max_abs_P_at"]["channels"]),
                        stat["max_abs_P_at"]["n"],
                        fmt(stat["max_abs_P_at"]["x"]),
                        stat["max_abs_P_at"]["region"],
                        stat["max_abs_P_at"]["kind"],
                        fmt(stat["block_scale_at_max"]),
                        fmt(stat["relative_at_max"]),
                        fmt(stat["max_abs_c_block"]),
                        fmt(stat["max_abs_c_full"]),
                        stat["components"],
                        stat["exactly_zero_P"],
                        stat["at_or_below_floor_P"],
                        stat["nonzero_above_floor_P"],
                        fmt(note_rel["max_relative"]),
                        note_rel["material_components_10pct"],
                        note_rel["material_readings_10pct"],
                        note_rel["material_points_10pct"],
                        fmt(pure_rel["max_relative"]),
                        pure_rel["material_components_10pct"],
                        pure_rel["material_readings_10pct"],
                        fmt(alias["self_term"]),
                        fmt(alias["alias_by_order"].get("4", 0.0)),
                        fmt(alias["alias_by_order"].get("5", 0.0)),
                        fmt(alias["alias_by_order"].get("6", 0.0)),
                        fmt(alias["self_fraction"])
                        if alias["self_fraction"] is not None else "",
                    ]
                )
    write_csv(
        out_dir / "d79_orders.csv",
        [
            "route", "lane", "order", "max_abs_P_note", "support", "channels",
            "at_n", "at_x", "region", "kind", "block_scale_at_max",
            "relative_at_max", "max_abs_c_block", "max_abs_c_full",
            "components", "exactly_zero_P", "at_or_below_floor_P",
            "nonzero_above_floor_P", "note_max_relative",
            "note_material_components_10pct",
            "note_material_readings_10pct",
            "note_material_points_10pct", "pure_max_relative",
            "pure_material_components_10pct", "pure_material_readings_10pct",
            "self_term_P",
            "alias_order4_P", "alias_order5_P", "alias_order6_P",
            "self_fraction",
        ],
        d79_rows,
    )

    precision_rows = []
    for route_name, route in routes.items():
        precision = route["precision"]
        for kind in ("precision_main_effect", "precision_by_call"):
            stat = precision[kind]
            precision_rows.append(
                [
                    route_name, kind, stat["precision_contrast"],
                    stat["coefficient"],
                    "+".join(stat["at"]["support"]) or "(call mean)",
                    "+".join(stat["at"]["channels"]) or "(call mean)",
                    stat["at"]["order"], stat["at"]["n"],
                    fmt(stat["at"]["x"]), stat["at"]["region"],
                    "|".join(fmt(v) for v in stat["lane_values"]),
                ]
            )
        for order, stat in precision["precision_by_call_order"].items():
            precision_rows.append(
                [
                    route_name, f"precision_by_call_order_{order}",
                    stat["precision_contrast"], stat["coefficient"],
                    "+".join(stat["at"]["support"]) or "(call mean)",
                    "+".join(stat["at"]["channels"]) or "(call mean)",
                    stat["at"]["order"], stat["at"]["n"],
                    fmt(stat["at"]["x"]), stat["at"]["region"],
                    "|".join(fmt(v) for v in stat["lane_values"]),
                ]
            )
    write_csv(
        out_dir / "precision_maxima.csv",
        [
            "route", "kind", "precision_contrast", "coefficient",
            "call_support", "call_channels", "call_order", "at_n", "at_x",
            "region", "lane_values_fp64|fp32|fp16|bf16",
        ],
        precision_rows,
    )

    precision_component_rows = []
    for route_name, route in routes.items():
        for row in route["precision"]["component_rows"]:
            precision_component_rows.append(
                [
                    route_name, row["component_index"], row["order"],
                    row["support"], row["channels"],
                    fmt(row["max_abs_precision_contrast"]),
                    row["precision_contrast"], row["at_n"], fmt(row["at_x"]),
                    row["region"],
                    "|".join(fmt(v) for v in row["lane_values"]),
                ]
            )
    write_csv(
        out_dir / "precision_by_component.csv",
        [
            "route", "component_index", "order", "support", "channels",
            "max_abs_precision_contrast", "precision_contrast", "at_n",
            "at_x", "region", "lane_values_fp64|fp32|fp16|bf16",
        ],
        precision_component_rows,
    )

    runtime = time.perf_counter() - started
    summary = {
        "schema": 1,
        "title": "Full-cube mixed-level factorial decomposition of the "
        "576-cell composition record",
        "design": {
            "axes": AXES,
            "levels": {a: LEVELS[a] for a in AXES},
            "level_counts": LEVEL_COUNTS,
            "cells_per_lane": N_CELLS,
            "lanes": LANES,
            "points": "33 orders x 26 arguments",
            "regions": {
                "A": f"x < {KX0}",
                "B": f"{KX0} <= x < {KFLATHI}",
                "far": f"x >= {KFLATHI}",
            },
            "measurement_floor": MEASUREMENT_FLOOR,
            "materiality_alpha": ALPHA,
        },
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
            "runtime_seconds": runtime,
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "inputs": {
                route: routes[route]["checks"] for route in routes
            },
        },
        "checks": {
            "orthonormality": orthonormality,
            "reconstruction": {
                route: routes[route]["reconstruction"] for route in routes
            },
            "coefficient_transform_agreement": {
                route: routes[route]["transform_agreement"] for route in routes
            },
            "independent_coefficient_check": {
                route: routes[route]["independent_coefficient_check"]
                for route in routes
            },
            "data": {
                route: routes[route]["checks"] for route in routes
            },
        },
        "full_cube": {
            route: {
                lane: routes[route]["full_cube"][lane]
                for lane in LANES
            }
            for route in routes
        },
        "d79": {
            route: {
                lane: routes[route]["d79"][lane] for lane in LANES
            }
            for route in routes
        },
        "precision": {
            route: routes[route]["precision"] for route in routes
        },
        "note_comparison": comparisons,
        "notes": {
            "signed_error": "r = ((value - ref_hi) - ref_lo) / F0",
            "full_cube_coefficient": "c_alpha = (1/144) sum_cells r Phi_alpha",
            "d79_note_estimator": (
                "P_alpha = sum_block r psi_alpha with Euclidean-normalised "
                "rows; P_alpha = sqrt(B) * (1/B) sum_block r Phi_alpha, "
                "where B is the block size"
            ),
            "d79_block": (
                "complete factorial on the support axes with every other axis "
                "at its baseline level; the union over orders <= 3 is d79"
            ),
            "materiality": (
                "component share of the point's non-mean (centred) total sum "
                "of squares; alpha = 0.10"
            ),
            "precision_axis": (
                "four-level factor with rows mean, linear (-3,-1,1,3)/sqrt(5), "
                "quadratic (1,-1,-1,1)/2, cubic (-1,3,-3,1)/sqrt(5), "
                "orthonormal under the uniform measure over the four lanes"
            ),
            "far_region": (
                "x >= 31 is the single argument x = 32.  There the measured "
                "response is bit-identical across all 144 cells (checked and "
                "recorded as far_max_spread_across_cells = 0), so every "
                "contrast is mathematically zero; the ~1e-31 fp64 residues are "
                "summation roundoff.  Nonzero counts and order shares are "
                "reported both over all points and over regions A/B, and the "
                "A/B figures are the physical ones."
            ),
            "d79_zero_counts": (
                "the verification note's 'exactly zero' counts (12/44, 11/44, "
                "28/44, 28/44) are not reproduced by bit-exact zeros or by "
                "at-or-below-floor zeros of the note's P; both our counts and "
                "the note's are reported.  Its maxima, max relative and 10% "
                "material counts are reproduced exactly under the per-block "
                "maximum scale."
            ),
            "d79_alias": (
                "the note's P_alpha is a block functional, not the pure "
                "component: it sums c_beta over every beta whose support "
                "contains alpha, weighted by baseline contrast values on the "
                "extra axes.  The alias split at each order-3 maximum is "
                "reported as self_term, alias_by_order, and self_fraction."
            ),
        },
    }
    with open(out_dir / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(json_safe(summary), handle, indent=2)

    # Console headline.
    print("\northonormality:")
    for axis, stat in per_axis.items():
        print(
            f"  {axis}: levels {stat['levels']}, max |G-I| = "
            f"{stat['max_abs_gram_minus_identity']:.3e}"
        )
    print(
        "  full Phi: max |Phi Phi^T/144 - I| = "
        f"{orthonormality['phi_max_abs_gram_minus_identity']:.3e}"
    )
    for route_name, route in routes.items():
        print(f"\n{route_name} reconstruction:")
        for lane_name in LANES:
            stat = route["reconstruction"][lane_name]
            print(
                f"  {lane_name}: max |r - r_hat| = "
                f"{stat['max_abs_reconstruction_error']:.3e}, "
                f"max |r| = {stat['max_abs_r']:.3e}, "
                f"relative = "
                f"{stat['max_abs_reconstruction_error_over_max_abs_r']:.3e}"
            )
        agreement = route["transform_agreement"]
        print(
            f"  transform check: max |C_tensor - C_matmul| = "
            f"{agreement['max_abs_difference_tensor_vs_matmul']:.3e} "
            f"(relative {agreement['max_relative_difference']:.3e})"
        )
        independent = route["independent_coefficient_check"]
        print(
            "  independent group-mean coefficient check: max |diff| = "
            f"{max(v['max_abs_difference'] for v in independent.values()):.3e}"
        )
        for lane_name in RESOLVING_LANES:
            stability = agreement["order3_share_stability"][lane_name]
            tens = stability["tensor"]
            matm = stability["direct_matmul"]
            print(
                f"    {lane_name} order-3 share stability: tensor "
                f"{tens['max_order3_share_AB']:.6f} "
                f"({tens['order3_material_components_AB']}/44, "
                f"{tens['order3_material_readings_AB']} readings) vs "
                f"direct matmul {matm['max_order3_share_AB']:.6f} "
                f"({matm['order3_material_components_AB']}/44, "
                f"{matm['order3_material_readings_AB']} readings)"
            )
        print(f"\n{route_name} full cube:")
        for lane_name in RESOLVING_LANES:
            lane = route["full_cube"][lane_name]
            print(f"  {lane_name}:")
            for order in ("1", "2", "3"):
                stat = lane["orders"][order]
                at = stat["max_at"]
                share_max_ab = 0.0
                for region_name in ("A", "B"):
                    region_stat = stat["max_order_share_by_region"][region_name]
                    if region_stat:
                        share_max_ab = max(
                            share_max_ab, region_stat["max_order_share"]
                        )
                print(
                    f"    order {order}: max |c| = "
                    f"{stat['max_abs_c']:.6e} at "
                    f"{'+'.join(at['support'])} "
                    f"({'+'.join(at['channels'])}), n={at['n']}, "
                    f"x={at['x']:g}, region {at['region']}, "
                    f"{stat['material_components_10pct_nonfar']}/"
                    f"{stat['components']} material components (A/B), "
                    f"{stat['material_readings_10pct_nonfar']}/"
                    f"{stat['readings_total_nonfar']} material readings (A/B); "
                    f"nonzero A/B {stat['nonzero_above_floor_nonfar']}/"
                    f"{stat['components']}, max share A/B "
                    f"{share_max_ab:.4f}"
                )
            for order in ("4", "5", "6"):
                stat = lane["orders"][order]
                print(
                    f"    order {order}: max |c| = {stat['max_abs_c']:.6e} "
                    f"(global SS share {stat['global_ss_share']:.3e})"
                )
            nonadd = lane["non_additive"]
            print(
                f"    max non-additive residual = "
                f"{nonadd['max_abs_residual']:.6e} at n="
                f"{nonadd['at']['n']}, x={nonadd['at']['x']:g}, "
                f"{nonadd['at']['region']}; max order>=2 share A/B "
                f"{nonadd['max_order_ge2_share_nonfar']:.4f}; "
                f"material points A/B {nonadd['material_points_nonfar_10pct']}"
            )
            print(
                f"    far check: response spread across cells at x>=31 = "
                f"{route['checks']['far_max_spread_across_cells']:.3e}"
            )
        print(f"\n{route_name} d79 (note P):")
        for lane_name in RESOLVING_LANES:
            lane = route["d79"][lane_name]
            print(f"  {lane_name}:")
            for order in ("1", "2", "3"):
                stat = lane["orders"][order]
                at = stat["max_abs_P_at"]
                note_rel = stat["note_relative"]
                pure_rel = stat["pure_full_cube_relative"]
                print(
                    f"    order {order}: max |P| = "
                    f"{stat['max_abs_P_note']:.6e} at "
                    f"{'+'.join(at['support'])} "
                    f"({'+'.join(at['channels'])}), n={at['n']}, "
                    f"x={at['x']:g}, region {at['region']}, "
                    f"max rel {note_rel['max_relative']:.4f}, "
                    f"material {note_rel['material_components_10pct']}/"
                    f"{stat['components']} components, "
                    f"{note_rel['material_readings_10pct']} readings; "
                    f"pure component under the same measure "
                    f"{pure_rel['material_components_10pct']}/"
                    f"{stat['components']}, "
                    f"{pure_rel['material_readings_10pct']} readings; "
                    f"bit-exact zeros {stat['exactly_zero_P']}/"
                    f"{stat['components']}"
                )
            alias = lane["orders"]["3"]["alias_split"]
            print(
                f"    order-3 max alias split: self {alias['self_term']:.6e} "
                f"({alias['self_fraction']:.3f} of P), "
                f"4FI {alias['alias_by_order'].get('4', 0.0):.3e}, "
                f"5FI {alias['alias_by_order'].get('5', 0.0):.3e}, "
                f"6FI {alias['alias_by_order'].get('6', 0.0):.3e}"
            )
        print(f"\n{route_name} precision:")
        for kind in ("precision_main_effect", "precision_by_call"):
            stat = route["precision"][kind]
            print(
                f"  {kind}: {stat['coefficient']:.6e} "
                f"({stat['precision_contrast']}) at "
                f"{'+'.join(stat['at']['support']) or '(call mean)'}, "
                f"n={stat['at']['n']}, x={stat['at']['x']:g}, "
                f"lanes={['%.3e' % v for v in stat['lane_values']]}"
            )

    print(f"\nwrote {out_dir}/summary.json and CSVs")
    print(f"runtime {runtime:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
