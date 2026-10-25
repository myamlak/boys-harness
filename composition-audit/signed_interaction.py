#!/usr/bin/env python3
"""Recompute the interaction verdict on signed pointwise errors.

The audit's verdict statistic is built from the non-negative response
R = |e|/F_0 and from cell maxima over (n, x). A referee showed that the
absolute value can manufacture an interaction: an exactly additive signed
error, e = a on the first axis and e = b on the second, has

    e00 = 0, e10 = a, e01 = b, e11 = a + b,

so |e| is 0, |a|, |b|, |a+b| and the contrast of |e| is generally non-zero
even though the signed errors add exactly. The triangle inequality the rule
rests on applies to signed errors, so the interaction that matters for C2 is
the one the signed errors show.

This script reads the four deposited cells files and recomputes the pairwise
and multi-way statistics on

    e = (value - ref_hi) - ref_lo,     r = e / F_0,

using the same reference pair the driver subtracted. The ref_hi/ref_lo
columns of each cells file are the double-double split of the single 60-digit
mpmath reference; the reference.csv file in each run directory carries the
same pair and is cross-checked here row by row.

The pairwise contrast is pointwise: at one (n, x), with the same completion
and the same two level indices as the paper's cell-max contrast,

    pi = r11 - r10 - r01 + r00.

It is informative where |r10 - r00| + |r01 - r00| > 0, and material at margin
alpha where |pi| > alpha * max(|r00|, |r10|, |r01|, |r11|), the pointwise
analogue of the paper's |I| > alpha * max(E). The multi-way statistic
decomposes the signed response at each point with the same orthonormal
contrasts and reports each order's share of the point's centred total sum of
squares.

Two sensitivity statistics are reported beside the pointwise one:
  * a signed cell-max contrast, the paper's construction with the signed
    response in place of |e|;
  * a pointwise numerator over the paper's cell-max denominator, which asks
    whether a signed pointwise interaction clears the scale the paper used.
The script also reproduces the paper's absolute-value counts from the same
data, so a difference between the two is a difference of statistic and not of
input.

Usage, from the paper directory:

    python scripts/signed_interaction.py

It writes records/signed-interaction.json and records/signed-interaction.tsv.
Nothing here re-runs the driver.
"""

import csv
import hashlib
import itertools
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import summarise_interaction as si

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "records"

AXES = si.AXES
WITHIN = si.WITHIN
LEVELS = si.LEVELS
ALPHAS = si.ALPHAS
PRIMARY_ALPHA = si.PRIMARY_ALPHA
RESOLVING = si.RESOLVING_LANES

INPUTS = {
    "msvc-fused": ROOT / "results" / "cells.csv",
    "msvc-separate": ROOT / "results-msvc-separate" / "cells.csv",
    "gcc-fused": ROOT / "results-gcc-fused" / "cells.csv",
    "gcc-separate": ROOT / "results-gcc-separate" / "cells.csv",
}

REFERENCE = ROOT / "results" / "reference.csv"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_reference(path):
    """The single reference pair, keyed by (lane, n, x)."""
    reference = {}
    with open(path, newline="", encoding="utf-8") as f:
        for raw in csv.DictReader(f):
            reference[(raw["lane"], int(raw["n"]), float(raw["x"]))] = (
                float(raw["hi"]),
                float(raw["lo"]),
            )
    return reference


def read_cells(path, reference):
    """Resolving-lane rows as (signed, absolute, flags) per cell and point."""
    cells = defaultdict(dict)
    rows = 0
    refusals = 0
    reference_checked = 0
    reference_mismatch = []
    kind_counts = defaultdict(int)

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        column = {name: index for index, name in enumerate(header)}
        get = column.__getitem__
        for raw in reader:
            rows += 1
            lane = raw[get("precision")]
            if lane not in RESOLVING:
                continue
            if raw[get("status")] != "ok":
                refusals += 1
                continue
            key = tuple(raw[get(axis)] for axis in WITHIN)
            n = int(raw[get("n")])
            x_nominal = float(raw[get("x_nominal")])
            x = float(raw[get("x")])
            value = float(raw[get("value")])
            ref_hi = float(raw[get("ref_hi")])
            ref_lo = float(raw[get("ref_lo")])
            reference_checked += 1
            expected = reference.get((lane, n, x))
            if expected != (ref_hi, ref_lo):
                if len(reference_mismatch) < 12:
                    reference_mismatch.append(
                        {"lane": lane, "n": n, "x": x,
                         "cells": [ref_hi, ref_lo],
                         "reference.csv": list(expected) if expected else None}
                    )
            f0 = float(raw[get("f0")])
            # The driver forms (value - hi) - lo; do the same two operations.
            signed = (value - ref_hi) - ref_lo
            absolute = abs(signed)
            recorded = float(raw[get("abs_err")])
            if absolute != recorded:
                raise SystemExit(
                    f"{path}: abs_err column disagrees with (value - hi) - lo "
                    f"at {lane} n={n} x={x}: {recorded} vs {absolute}"
                )
            signed_r = signed / f0 if f0 > 0.0 else 0.0
            absolute_r = absolute / f0 if f0 > 0.0 else 0.0
            cells[(lane, key)][(n, x_nominal)] = (
                signed_r,
                absolute_r,
                int(raw[get("pack_live")]),
                int(raw[get("partition_in_domain")]),
                raw[get("argument_kind")],
            )
            kind_counts[raw[get("argument_kind")]] += 1

    if reference_mismatch:
        raise SystemExit(
            "reference mismatch between cells.csv and reference.csv: "
            + json.dumps(reference_mismatch)
        )
    return {
        "cells": cells,
        "rows": rows,
        "refusals": refusals,
        "reference_checked": reference_checked,
        "kind_counts": dict(kind_counts),
    }


def domain_ok(rec, pair):
    if "pack" in pair and not rec[2]:
        return False
    if "partition" in pair and not rec[3]:
        return False
    return True


def lane_view(cells, lane):
    """Cells of one lane, keyed by the six-axis key."""
    return {
        key: points
        for (row_lane, key), points in cells.items()
        if row_lane == lane
    }


def combined_view(cells):
    """Both resolving lanes, keyed by (lane,) + six-axis key."""
    view = {}
    for (lane, key), points in cells.items():
        view[(lane,) + key] = points
    return view


def pair_contrasts(view, axes):
    """Signed pointwise, signed cell-max and paper-statistic contrasts.

    `view` maps a cell key to its points. When `lane` is given the key is the
    six-axis tuple and the lane is fixed; otherwise precision is one of
    `axes` and the key is (lane,) + six-axis tuple.
    """
    contrasts = []
    for a, b in itertools.combinations(axes, 2):
        other = [axis for axis in axes if axis not in (a, b)]
        for completion in itertools.product(*(LEVELS[axis] for axis in other)):
            base = dict(zip(other, completion))
            level_a0 = LEVELS[a][0]
            level_b0 = LEVELS[b][0]

            def key_for(level_a, level_b):
                setting = dict(base)
                setting[a] = level_a
                setting[b] = level_b
                if "precision" in axes:
                    cell_lane = setting["precision"]
                    return (cell_lane,) + tuple(setting[axis] for axis in WITHIN)
                return tuple(setting[axis] for axis in WITHIN)

            def get(level_a, level_b):
                return view.get(key_for(level_a, level_b))

            for level_a in LEVELS[a][1:]:
                for level_b in LEVELS[b][1:]:
                    blocks = {
                        "00": get(level_a0, level_b0),
                        "10": get(level_a, level_b0),
                        "01": get(level_a0, level_b),
                        "11": get(level_a, level_b),
                    }
                    if any(points is None for points in blocks.values()):
                        continue
                    pair = (a, b)
                    material_points = {str(alpha): 0 for alpha in ALPHAS}
                    informative_points = 0
                    point_count = 0
                    max_pointwise = None
                    cell_max_abs = {name: 0.0 for name in blocks}
                    cell_max_signed = {name: None for name in blocks}

                    # First pass: the four cell maxima, needed by both the
                    # paper's statistic and the guard on the local scale.
                    valid = []
                    for point, r00 in blocks["00"].items():
                        recs = {"00": r00}
                        missing = False
                        for name in ("10", "01", "11"):
                            rec = blocks[name].get(point)
                            if rec is None:
                                missing = True
                                break
                            recs[name] = rec
                        if missing:
                            continue
                        if not all(domain_ok(recs[name], pair) for name in recs):
                            continue
                        valid.append((point, recs))
                        for name, rec in recs.items():
                            cell_max_abs[name] = max(cell_max_abs[name], rec[1])
                            if (cell_max_signed[name] is None
                                    or rec[0] > cell_max_signed[name]):
                                cell_max_signed[name] = rec[0]

                    scale_cell_abs = max(cell_max_abs.values())
                    guard_scale = 0.1 * scale_cell_abs

                    # Second pass: the signed pointwise statistic, with a
                    # guarded count that ignores points whose four responses
                    # are all below one tenth of the largest cell response,
                    # the floor the supplement uses for its held-out margins.
                    for point, recs in valid:
                        point_count += 1
                        r = {name: recs[name][0] for name in recs}
                        pi = r["11"] - r["10"] - r["01"] + r["00"]
                        m_a = r["10"] - r["00"]
                        m_b = r["01"] - r["00"]
                        scale = max(abs(r[name]) for name in r)
                        if abs(m_a) + abs(m_b) > 0.0:
                            informative_points += 1
                        if scale > 0.0 and abs(m_a) + abs(m_b) > 0.0:
                            relative = abs(pi) / scale
                            for alpha in ALPHAS:
                                if abs(pi) > alpha * scale:
                                    material_points[str(alpha)] += 1
                                if scale >= guard_scale and abs(pi) > alpha * scale:
                                    material_points.setdefault("guarded", {})
                                    material_points["guarded"].setdefault(
                                        str(alpha), 0)
                                    material_points["guarded"][str(alpha)] += 1
                            if max_pointwise is None or relative > max_pointwise["relative"]:
                                max_pointwise = {
                                    "relative": relative,
                                    "interaction": pi,
                                    "at": {"n": point[0], "x_nominal": point[1]},
                                    "scale": scale,
                                    "signed": r,
                                }

                    paper = None
                    if scale_cell_abs > 0.0:
                        e = cell_max_abs
                        interaction = e["11"] - e["10"] - e["01"] + e["00"]
                        paper = {
                            "interaction": interaction,
                            "relative": abs(interaction) / scale_cell_abs,
                            "material_at": {
                                str(alpha): bool(
                                    abs(e["10"] - e["00"]) + abs(e["01"] - e["00"]) > 0.0
                                    and abs(interaction) > alpha * scale_cell_abs
                                )
                                for alpha in ALPHAS
                            },
                        }

                    scale_cell_signed = (
                        max(abs(value) for value in cell_max_signed.values())
                        if all(value is not None for value in cell_max_signed.values())
                        else 0.0
                    )
                    signed_max = None
                    if scale_cell_signed > 0.0:
                        s = cell_max_signed
                        interaction = s["11"] - s["10"] - s["01"] + s["00"]
                        signed_max = {
                            "interaction": interaction,
                            "relative": abs(interaction) / scale_cell_signed,
                            "material_at": {
                                str(alpha): bool(
                                    abs(s["10"] - s["00"]) + abs(s["01"] - s["00"]) > 0.0
                                    and abs(interaction) > alpha * scale_cell_signed
                                )
                                for alpha in ALPHAS
                            },
                        }

                    cell_scale_pointwise = None
                    if max_pointwise is not None and scale_cell_abs > 0.0:
                        cell_scale_pointwise = {
                            "relative": abs(max_pointwise["interaction"]) / scale_cell_abs,
                            "material_at": {
                                str(alpha): bool(
                                    abs(max_pointwise["interaction"]) > alpha * scale_cell_abs
                                )
                                for alpha in ALPHAS
                            },
                        }

                    contrasts.append(
                        {
                            "pair": list(pair),
                            "levels": {a: [level_a0, level_a],
                                       b: [level_b0, level_b]},
                            "completion": base,
                            "points": point_count,
                            "informative_points": informative_points,
                            "material_points": material_points,
                            "max_pointwise": max_pointwise,
                            "paper_cell_max": paper,
                            "signed_cell_max": signed_max,
                            "pointwise_over_cell_scale": cell_scale_pointwise,
                        }
                    )
    return contrasts


def multiway(view, axes, points_filter=None, levels_override=None):
    """Signed multi-way decomposition at every common point."""
    levels_override = levels_override or {}
    level_lists = [levels_override.get(axis, LEVELS[axis]) for axis in axes]
    levels = [len(level_list) for level_list in level_lists]
    metadata = si.component_metadata(levels)
    combinations = list(itertools.product(*level_lists))

    if "precision" in axes:
        keys = [
            (dict(zip(axes, combo))["precision"],)
            + tuple(dict(zip(axes, combo))[axis] for axis in WITHIN)
            for combo in combinations
        ]
    else:
        keys = [
            tuple(dict(zip(axes, combo))[axis] for axis in WITHIN)
            for combo in combinations
        ]

    if any(key not in view for key in keys):
        missing = sum(1 for key in keys if key not in view)
        return {"points": 0, "note": "incomplete factorial",
                "missing_cells": missing}

    common = set(view[keys[0]])
    for key in keys[1:]:
        common &= set(view[key])
    points = sorted(common)
    if points_filter is not None:
        points = [point for point in points if points_filter(point)]

    material_points = {str(alpha): 0 for alpha in ALPHAS}
    max_order_share = {str(order): 0.0 for order in range(2, len(axes) + 1)}
    max_order_share_at = {}
    max_interaction_coefficient = 0.0
    max_interaction_at = None

    for point in points:
        raw = [view[key][point][0] for key in keys]
        effects = si.decompose_values(raw, levels)
        coefficients = [
            effects[index] / meta["norm"] for index, meta in enumerate(metadata)
        ]
        centered_ss = sum(
            coefficient * coefficient
            for coefficient, meta in zip(coefficients, metadata)
            if meta["order"] >= 1
        )
        if centered_ss <= 0.0:
            continue
        order_ss = defaultdict(float)
        for coefficient, meta in zip(coefficients, metadata):
            if meta["order"] >= 2:
                order_ss[meta["order"]] += coefficient * coefficient
                if abs(coefficient) > max_interaction_coefficient:
                    max_interaction_coefficient = abs(coefficient)
                    max_interaction_at = {
                        "n": point[0], "x_nominal": point[1],
                        "order": meta["order"],
                        "support": [axes[i] for i in meta["support"]],
                    }
        for order, value in order_ss.items():
            share = value / centered_ss
            if share > max_order_share[str(order)]:
                max_order_share[str(order)] = share
                max_order_share_at[str(order)] = {
                    "n": point[0], "x_nominal": point[1]
                }
            for alpha in ALPHAS:
                if share > alpha:
                    material_points[str(alpha)] += 1

    return {
        "factors": list(axes),
        "points": len(points),
        "material_points": material_points,
        "max_order_share": max_order_share,
        "max_order_share_at": max_order_share_at,
        "max_interaction_coefficient": max_interaction_coefficient,
        "max_interaction_at": max_interaction_at,
    }


def verdict(pairwise, multiway_blocks, alphas=ALPHAS):
    pairwise_material = {str(alpha): 0 for alpha in alphas}
    informative = 0
    for contrast in pairwise:
        if contrast["informative_points"] > 0:
            informative += 1
        for alpha in alphas:
            if contrast["material_points"][str(alpha)] > 0:
                pairwise_material[str(alpha)] += 1
    multiway_material = {str(alpha): 0 for alpha in alphas}
    for block in multiway_blocks:
        for alpha in alphas:
            multiway_material[str(alpha)] += block["material_points"][str(alpha)]
    outcome = {}
    for alpha in alphas:
        key = str(alpha)
        if pairwise_material[key] > 0 or multiway_material[key] > 0:
            outcome[key] = "interaction-material"
        elif informative == 0:
            outcome[key] = "inconclusive"
        else:
            outcome[key] = f"additive-within-{int(round(alpha * 100))}pct"
    return {
        "pairwise_material_contrasts": pairwise_material,
        "multiway_material_points": multiway_material,
        "informative_contrasts": informative,
        "verdict_by_alpha": outcome,
    }


def summarise_pair_counts(pairwise):
    """Per pair: material contrasts and the largest signed pointwise contrast."""
    out = defaultdict(
        lambda: {
            "contrasts": 0,
            "material": {str(alpha): 0 for alpha in ALPHAS},
            "max_relative": 0.0,
            "max_at": None,
        }
    )
    for contrast in pairwise:
        name = "+".join(contrast["pair"])
        entry = out[name]
        entry["contrasts"] += 1
        for alpha in ALPHAS:
            if contrast["material_points"][str(alpha)] > 0:
                entry["material"][str(alpha)] += 1
        largest = contrast["max_pointwise"]
        if largest is not None and largest["relative"] > entry["max_relative"]:
            entry["max_relative"] = largest["relative"]
            entry["max_at"] = {
                "levels": contrast["levels"],
                "completion": contrast["completion"],
                "at": largest["at"],
                "interaction": largest["interaction"],
            }
    return dict(out)


def statistic_block(pairwise, multiway_blocks):
    return {
        "pairwise": {
            "material_contrasts": {
                str(alpha): sum(
                    1 for c in pairwise if c["material_points"][str(alpha)] > 0
                )
                for alpha in ALPHAS
            },
            "material_point_readings": {
                str(alpha): sum(c["material_points"][str(alpha)] for c in pairwise)
                for alpha in ALPHAS
            },
            "informative_contrasts": sum(
                1 for c in pairwise if c["informative_points"] > 0
            ),
            "per_pair": summarise_pair_counts(pairwise),
            "absolute_value_reproduction": {
                str(alpha): sum(
                    1 for c in pairwise
                    if c["paper_cell_max"] is not None
                    and c["paper_cell_max"]["material_at"][str(alpha)]
                )
                for alpha in ALPHAS
            },
        },
        "multiway": multiway_blocks,
        "verdict": verdict(pairwise, multiway_blocks),
        "sensitivity": {
            "signed_cell_max_material_contrasts": {
                str(alpha): sum(
                    1 for c in pairwise
                    if c["signed_cell_max"] is not None
                    and c["signed_cell_max"]["material_at"][str(alpha)]
                )
                for alpha in ALPHAS
            },
            "pointwise_over_cell_max_scale_material_contrasts": {
                str(alpha): sum(
                    1 for c in pairwise
                    if c["pointwise_over_cell_scale"] is not None
                    and c["pointwise_over_cell_scale"]["material_at"][str(alpha)]
                )
                for alpha in ALPHAS
            },
            "guarded_pointwise_material_contrasts": {
                str(alpha): sum(
                    1 for c in pairwise
                    if c["material_points"].get("guarded", {}).get(str(alpha), 0) > 0
                )
                for alpha in ALPHAS
            },
        },
    }


def self_test():
    """The referee's exactly additive example, computed both ways."""
    a, b = 1.0, -0.6
    signed = {"00": 0.0, "10": a, "01": b, "11": a + b}
    absolute = {name: abs(value) for name, value in signed.items()}
    signed_contrast = signed["11"] - signed["10"] - signed["01"] + signed["00"]
    absolute_contrast = (
        absolute["11"] - absolute["10"] - absolute["01"] + absolute["00"]
    )
    scale = max(absolute.values())
    return {
        "signed_error": signed,
        "absolute_error": absolute,
        "signed_contrast": signed_contrast,
        "absolute_contrast": absolute_contrast,
        "absolute_relative": abs(absolute_contrast) / scale,
        "material_at_10pct_by_absolute_value": abs(absolute_contrast) > 0.10 * scale,
        "material_at_10pct_by_signed_value": abs(signed_contrast) > 0.10 * scale,
    }


def analyse(path, reference):
    loaded = read_cells(path, reference)
    cells = loaded["cells"]
    lanes = {lane: lane_view(cells, lane) for lane in RESOLVING}
    combined = combined_view(cells)

    result = {
        "path": path.relative_to(ROOT).as_posix(),
        "sha256": sha256(path),
        "rows": loaded["rows"],
        "refusals": loaded["refusals"],
        "reference_checked": loaded["reference_checked"],
        "rows_by_argument_kind": loaded["kind_counts"],
        "lanes": {},
    }

    lane_blocks = {}
    for lane, view in lanes.items():
        pairwise = pair_contrasts(view, WITHIN)
        multi = multiway(view, WITHIN)
        lane_blocks[lane] = multi
        result["lanes"][lane] = {
            "cells": len(view),
            "points": len(next(iter(view.values()))),
            "statistic": statistic_block(pairwise, [multi]),
        }

    pairwise_all = pair_contrasts(combined, AXES)
    cross = multiway(
        combined, AXES,
        points_filter=lambda point: combined[next(iter(combined))][point][4] == "exact",
        levels_override={"precision": RESOLVING},
    )
    result["paper_wide"] = statistic_block(
        pairwise_all, [lane_blocks["fp64"], lane_blocks["fp32"], cross]
    )
    result["paper_wide"]["multiway_blocks"] = {
        "fp64": lane_blocks["fp64"],
        "fp32": lane_blocks["fp32"],
        "cross_lane_exact": cross,
    }
    return result


def tsv_rows(results):
    header = [
        "scope", "configuration", "lane", "alpha",
        "pairwise_material_contrasts", "pairwise_material_point_readings",
        "multiway_material_points", "signed_verdict",
        "paper_statistic_material_contrasts",
        "signed_cell_max_material_contrasts",
        "pointwise_over_cell_scale_material_contrasts",
        "guarded_pointwise_material_contrasts",
    ]
    rows = []
    for configuration, result in results.items():
        for lane in RESOLVING:
            block = result["lanes"][lane]["statistic"]
            for alpha in ALPHAS:
                key = str(alpha)
                rows.append([
                    "route-lane", configuration, lane, key,
                    block["pairwise"]["material_contrasts"][key],
                    block["pairwise"]["material_point_readings"][key],
                    block["multiway"][0]["material_points"][key],
                    block["verdict"]["verdict_by_alpha"][key],
                    block["pairwise"]["absolute_value_reproduction"][key],
                    block["sensitivity"]["signed_cell_max_material_contrasts"][key],
                    block["sensitivity"][
                        "pointwise_over_cell_max_scale_material_contrasts"][key],
                    block["sensitivity"]["guarded_pointwise_material_contrasts"][key],
                ])
        block = result["paper_wide"]
        for alpha in ALPHAS:
            key = str(alpha)
            rows.append([
                "paper-wide", configuration, "fp64+fp32", key,
                block["pairwise"]["material_contrasts"][key],
                block["pairwise"]["material_point_readings"][key],
                block["verdict"]["multiway_material_points"][key],
                block["verdict"]["verdict_by_alpha"][key],
                block["pairwise"]["absolute_value_reproduction"][key],
                block["sensitivity"]["signed_cell_max_material_contrasts"][key],
                block["sensitivity"][
                    "pointwise_over_cell_max_scale_material_contrasts"][key],
                block["sensitivity"]["guarded_pointwise_material_contrasts"][key],
            ])
    return header, rows


def main():
    reference = read_reference(REFERENCE)
    reference_digest = sha256(REFERENCE)
    results = {}
    for configuration, path in INPUTS.items():
        print(f"reading {path.relative_to(ROOT)}", flush=True)
        results[configuration] = analyse(path, reference)
        wide = results[configuration]["paper_wide"]
        print(
            "  paper-wide signed verdict: "
            + json.dumps(wide["verdict"]["verdict_by_alpha"]),
            flush=True,
        )

    route_checks = {}
    for route in ("fused", "separate"):
        digests = {
            results[name]["sha256"]
            for name in (f"msvc-{route}", f"gcc-{route}")
        }
        route_checks[route] = {
            "byte_identical": len(digests) == 1,
            "sha256": sorted(digests),
        }

    record = {
        "schema": "signed-interaction/1",
        "paths_relative_to": "papers/composed-bounds",
        "command": "python scripts/signed_interaction.py",
        "python": platform.python_version(),
        "script": {"path": "scripts/signed_interaction.py"},
        "statistic": {
            "signed_error": "e = (value - ref_hi) - ref_lo",
            "response": "r = e / F_0",
            "pairwise": (
                "at each (n, x): pi = r11 - r10 - r01 + r00; informative "
                "where |r10-r00|+|r01-r00| > 0; material at alpha where "
                "|pi| > alpha * max(|r00|,|r10|,|r01|,|r11|)"
            ),
            "multiway": (
                "orthonormal decomposition of r at each point; an interaction "
                "order's share of the point's centred total sum of squares > "
                "alpha"
            ),
        },
        "reference": {
            "source": "ref_hi and ref_lo columns of each cells.csv",
            "why": (
                "these are the two halves the driver subtracts, copied from the "
                "run's reference.csv; the reference is the single 60-digit "
                "mpmath evaluation at each lane's narrowed argument"
            ),
            "reference_csv": REFERENCE.relative_to(ROOT).as_posix(),
            "reference_csv_sha256": reference_digest,
            "rows_checked": sum(
                result["reference_checked"] for result in results.values()
            ),
            "mismatches": 0,
        },
        "inputs": {
            name: {
                "path": result["path"],
                "sha256": result["sha256"],
                "rows": result["rows"],
                "refusals": result["refusals"],
            }
            for name, result in results.items()
        },
        "route_checks": route_checks,
        "self_test": self_test(),
        "configuration_results": results,
    }

    RECORDS.mkdir(exist_ok=True)
    (RECORDS / "signed-interaction.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8")
    header, rows = tsv_rows(results)
    (RECORDS / "signed-interaction.tsv").write_text(
        "\n".join(["\t".join(header)]
                  + ["\t".join(str(value) for value in row) for row in rows])
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {RECORDS / 'signed-interaction.json'}")
    print(f"wrote {RECORDS / 'signed-interaction.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
