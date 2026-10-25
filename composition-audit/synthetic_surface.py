"""Synthetic surfaces with known interaction order.

The paper's interaction verdict is a statement about signed errors, and its
design envelope is a statement about absolute errors. This script separates the
two on surfaces whose structure is known: an exactly additive signed-error
surface, the same surface plus one planted pairwise interaction, and the same
surface plus one planted three-factor interaction. It runs the paper's
statistics and the paper's 1, 9, 35, 79 and 144-cell designs over a fixed set of
seeds and writes the aggregate rates, coverage and margins.

The surfaces live on the paper's six axes with levels (2, 2, 3, 2, 3, 2). Cell
values are linear combinations of the orthonormal factorial basis, so an
additive surface has exactly zero signed pairwise contrast at every point and
every completion. The absolute-value statistic is applied to the cell maxima of
|e|, and the signed statistic to the signed pointwise errors, both as in
scripts/summarise_interaction.py and scripts/signed_interaction.py.

Two regimes bracket the sign structure of the paper's own errors. In the
base-dominated regime the common term keeps every cell positive, so |e| = e
pointwise; in the sign-changing regime there is no common term and the cells
cross zero, so |e| is not e. The record tests three predictions fixed before
the run: that the absolute-value statistic reports interactions on an additive
signed surface while the signed statistic abstains, and that the 35-cell design
covers a pairwise surface but not a three-factor one.

Run from the paper directory:
  python -I scripts/synthetic_surface.py

Outputs:
  records/synthetic-surface.json   seeds, parameters and full-precision results

No external input is read. The script is deterministic for a fixed seed list.
"""

import argparse
import hashlib
import itertools
import json
import platform
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "records" / "synthetic-surface.json"

FACTORS = ["route", "scheme", "partition", "pack", "form", "exp"]
LEVELS = [2, 2, 3, 2, 3, 2]
ALPHAS = (0.05, 0.10, 0.20)
PRIMARY = 0.10
SEED_BASE = 20261023
DEFAULT_SEEDS = 200
DEFAULT_POINTS = 8

# A base-dominated regime has a common term much larger than the axis terms,
# so the signed surface keeps one sign and |e| is e pointwise. A sign-changing
# regime has no common term, so the signed errors cross zero and |e| is not e
# pointwise. The pair brackets the paper's own errors, which are small signed
# differences of both signs.
REGIMES = {
    "base-dominated": {"base": 1.0, "main_amplitude": 0.2},
    "sign-changing": {"base": 0.0, "main_amplitude": 0.5},
}
MAGNITUDES = (0.02, 0.05, 0.10, 0.20, 0.40, 0.80)

CONTRASTS = {
    2: ((1.0, 1.0), (1.0, -1.0)),
    3: ((1.0, 1.0, 1.0), (-1.0, 0.0, 1.0), (1.0, -2.0, 1.0)),
}

STRIDES = []
for a in range(6):
    stride = 1
    for b in range(a + 1, 6):
        stride *= LEVELS[b]
    STRIDES.append(stride)

CELLS = list(itertools.product(*[range(level) for level in LEVELS]))
PAIRS = list(itertools.combinations(range(6), 2))
TRIPLES = list(itertools.combinations(range(6), 3))
PAIR_INDEX = {pair: i for i, pair in enumerate(PAIRS)}


def axis_matrix(level_count):
    matrix = np.array(CONTRASTS[level_count], dtype=float)
    matrix /= np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix


def basis_matrix():
    """Orthonormal factorial basis, one row per cell and one column per component.

    A Kronecker product would transpose the factor matrices on the row side;
    the explicit product keeps the entry as M[channel, level] for every axis.
    """
    matrices = [axis_matrix(level) for level in LEVELS]
    matrix = np.zeros((len(CELLS), len(CELLS)))
    for cell_index, cell in enumerate(CELLS):
        for component in range(len(CELLS)):
            channels = tuple((component // stride) % level
                             for stride, level in zip(STRIDES, LEVELS))
            value = 1.0
            for axis in range(6):
                value *= matrices[axis][channels[axis], cell[axis]]
            matrix[cell_index, component] = value
    return matrix


PHI = basis_matrix()


def flat(cell):
    return sum(level * stride for level, stride in zip(cell, STRIDES))


def component_index(channels):
    return sum(channel * stride for channel, stride in zip(channels, STRIDES))


def cell_with(assign):
    cell = [0] * 6
    for axis, level in assign.items():
        cell[axis] = level
    return flat(cell)


def nonbase(index):
    return sum(1 for level in CELLS[index] if level != 0)


def design_keys(max_nonbase):
    return [index for index in range(len(CELLS)) if nonbase(index) <= max_nonbase]


def functionals(design, order):
    """Coefficient functionals and model basis for one hierarchical design.

    This mirrors fit_point and model_at in scripts/certification_curve.py: the
    mean is the design average, a main coefficient is the raw projection on a
    one-factor block, and pair and triple coefficients are raw projections on
    complete blocks. A block absent from the design contributes no column.
    """
    design_set = set(design)
    rows = []
    cols = []

    mean_row = np.zeros(len(CELLS))
    mean_row[design] = 1.0 / len(design)
    rows.append(mean_row)
    cols.append(np.ones(len(CELLS)))

    if order >= 1:
        for axis in range(6):
            block = [cell_with({axis: level}) for level in range(LEVELS[axis])]
            if not all(cell in design_set for cell in block):
                continue
            for channel in range(1, LEVELS[axis]):
                row = CONTRASTS[LEVELS[axis]][channel]
                norm = float(np.linalg.norm(row))
                functional = np.zeros(len(CELLS))
                for level, cell in enumerate(block):
                    functional[cell] = row[level] / norm
                rows.append(functional)
                cols.append(np.array([row[level] / norm
                                      for level in (cell[axis] for cell in CELLS)]))

    if order >= 2:
        for i, j in PAIRS:
            block = {(a, b): cell_with({i: a, j: b})
                     for a in range(LEVELS[i]) for b in range(LEVELS[j])}
            if not all(cell in design_set for cell in block.values()):
                continue
            for ci in range(1, LEVELS[i]):
                ri = CONTRASTS[LEVELS[i]][ci]
                ni = float(np.linalg.norm(ri))
                for cj in range(1, LEVELS[j]):
                    rj = CONTRASTS[LEVELS[j]][cj]
                    nj = float(np.linalg.norm(rj))
                    functional = np.zeros(len(CELLS))
                    for (a, b), cell in block.items():
                        functional[cell] = ri[a] * rj[b] / (ni * nj)
                    rows.append(functional)
                    cols.append(np.array([
                        ri[cell[i]] * rj[cell[j]] / (ni * nj) for cell in CELLS]))

    if order >= 3:
        for i, j, k in TRIPLES:
            block = {(a, b, c): cell_with({i: a, j: b, k: c})
                     for a in range(LEVELS[i])
                     for b in range(LEVELS[j])
                     for c in range(LEVELS[k])}
            if not all(cell in design_set for cell in block.values()):
                continue
            for ci in range(1, LEVELS[i]):
                ri = CONTRASTS[LEVELS[i]][ci]
                ni = float(np.linalg.norm(ri))
                for cj in range(1, LEVELS[j]):
                    rj = CONTRASTS[LEVELS[j]][cj]
                    nj = float(np.linalg.norm(rj))
                    for ck in range(1, LEVELS[k]):
                        rk = CONTRASTS[LEVELS[k]][ck]
                        nk = float(np.linalg.norm(rk))
                        functional = np.zeros(len(CELLS))
                        for (a, b, c), cell in block.items():
                            functional[cell] = ri[a] * rj[b] * rk[c] / (ni * nj * nk)
                        rows.append(functional)
                        cols.append(np.array([
                            ri[cell[i]] * rj[cell[j]] * rk[cell[k]] / (ni * nj * nk)
                            for cell in CELLS]))

    return np.array(rows), np.array(cols).T


def design_maps():
    designs = {
        "d1": (design_keys(0), 0),
        "d9": (design_keys(1), 1),
        "d35": (design_keys(2), 2),
        "d79": (design_keys(3), 3),
    }
    maps = {}
    for name, (design, order) in designs.items():
        functional, model_basis = functionals(design, order)
        maps[name] = {
            "design": np.array(design, dtype=np.int64),
            "order": order,
            "B": functional,
            "A": model_basis,
            "cells": len(design),
            "components": len(functional),
        }
    maps["d144"] = {
        "design": np.arange(len(CELLS), dtype=np.int64),
        "order": 6,
        "B": PHI.T,
        "A": PHI,
        "cells": len(CELLS),
        "components": len(CELLS),
    }
    return maps


DESIGN_MAPS = design_maps()


def envelope(values, entry):
    """Uniform bound, design allowance and fitted model for one design.

    values has shape (points, cells). The bound charges every non-mean
    coefficient by the largest value of any orthonormal basis function, 1/6.
    """
    coefficients = values @ entry["B"].T
    mean_bound = float(np.max(np.abs(coefficients[:, 0])))
    if coefficients.shape[1] > 1:
        component_bound = float(
            np.sum(np.max(np.abs(coefficients[:, 1:]), axis=0)))
    else:
        component_bound = 0.0
    bound = mean_bound + component_bound / 6.0
    model = coefficients @ entry["A"].T
    design = entry["design"]
    if design.size:
        allowance = float(np.max(np.abs(values[:, design] - model[:, design])))
    else:
        allowance = 0.0
    return bound, allowance, model


def score_design(values_abs, entry, bound, allowance):
    """Coverage and margin of one bound at one design.

    values_abs is |e|. A held reading is covered when it does not exceed the
    bound; the margin is the smallest bound-to-response ratio over held
    readings at or above one tenth of the largest held response, as in
    scripts/certification_curve.py.
    """
    cells = entry["cells"]
    design = set(entry["design"].tolist())
    held = np.array([index for index in range(len(CELLS)) if index not in design],
                    dtype=np.int64)
    total = bound + allowance
    record = {
        "cells": cells,
        "held_cells": int(held.size),
        "readings": int(held.size * values_abs.shape[0]),
        "bound": bound,
        "allowance": allowance,
        "violations": 0,
        "violations_no_allowance": 0,
        "worst": 0.0,
        "margin": None,
        "covered_cells": int(held.size),
    }
    if held.size == 0:
        return record
    held_values = values_abs[:, held]
    record["violations"] = int(np.sum(held_values > total))
    record["violations_no_allowance"] = int(np.sum(held_values > bound))
    record["worst"] = float(np.max(held_values / total)) if total > 0 else float("inf")
    largest = float(np.max(held_values))
    mask = (held_values > 0.0) & (held_values >= 0.1 * largest)
    if np.any(mask):
        record["margin"] = float(np.min(total / held_values[mask]))
    record["covered_cells"] = int(np.sum(np.all(held_values <= total, axis=0)))
    return record


def pair_blocks():
    blocks = []
    for i, j in PAIRS:
        others = [axis for axis in range(6) if axis not in (i, j)]
        completions = list(itertools.product(*[range(LEVELS[axis]) for axis in others]))
        index = np.zeros((len(completions), LEVELS[i], LEVELS[j]), dtype=np.int64)
        for position, completion in enumerate(completions):
            cell = [0] * 6
            for axis, level in zip(others, completion):
                cell[axis] = level
            for level_i in range(LEVELS[i]):
                cell[i] = level_i
                for level_j in range(LEVELS[j]):
                    cell[j] = level_j
                    index[position, level_i, level_j] = flat(cell)
        blocks.append((i, j, index))
    return blocks


PAIR_BLOCKS = pair_blocks()


def absolute_statistic(values, alphas):
    """Cell-maximum difference in differences on |e|.

    values is the signed error array, shape (points, cells). One contrast is
    one completion and one pair of non-baseline levels; it is material when
    |I| > alpha * max(E_ab).
    """
    response = np.abs(values)
    material = {f"{alpha:.2f}": 0 for alpha in alphas}
    pair_mask = {f"{alpha:.2f}": 0 for alpha in alphas}
    max_relative = 0.0
    informative = 0
    for position, (i, j, index) in enumerate(PAIR_BLOCKS):
        block = response[:, index]
        cell_max = block.max(axis=0)
        for level_i in range(1, LEVELS[i]):
            for level_j in range(1, LEVELS[j]):
                e00 = cell_max[:, 0, 0]
                e10 = cell_max[:, level_i, 0]
                e01 = cell_max[:, 0, level_j]
                e11 = cell_max[:, level_i, level_j]
                marginal_i = e10 - e00
                marginal_j = e01 - e00
                interaction = e11 - e10 - e01 + e00
                scale = np.maximum.reduce([np.abs(e00), np.abs(e10),
                                           np.abs(e01), np.abs(e11)])
                good = ((np.abs(marginal_i) + np.abs(marginal_j)) > 0.0) & (scale > 0.0)
                informative += int(np.sum(good))
                if np.any(good):
                    relative = np.abs(interaction[good]) / scale[good]
                    max_relative = max(max_relative, float(np.max(relative)))
                for alpha in alphas:
                    key = f"{alpha:.2f}"
                    hit = good & (np.abs(interaction) > alpha * scale)
                    count = int(np.sum(hit))
                    if count:
                        material[key] += count
                        pair_mask[key] |= 1 << position
    return {
        "material_contrasts": material,
        "pair_mask": pair_mask,
        "max_relative": max_relative,
        "informative_readings": informative,
    }


def signed_statistic(values, alphas):
    """Signed pointwise difference in differences on e.

    values is the signed error array, shape (points, cells). One point of one
    completion and one pair of non-baseline levels is material when
    |pi| > alpha * max|r_ab| at that point.
    """
    material = {f"{alpha:.2f}": 0 for alpha in alphas}
    readings = {f"{alpha:.2f}": 0 for alpha in alphas}
    pair_mask = {f"{alpha:.2f}": 0 for alpha in alphas}
    max_relative = 0.0
    informative = 0
    for position, (i, j, index) in enumerate(PAIR_BLOCKS):
        block = values[:, index]
        for level_i in range(1, LEVELS[i]):
            for level_j in range(1, LEVELS[j]):
                r00 = block[:, :, 0, 0]
                r10 = block[:, :, level_i, 0]
                r01 = block[:, :, 0, level_j]
                r11 = block[:, :, level_i, level_j]
                marginal_i = r10 - r00
                marginal_j = r01 - r00
                interaction = r11 - r10 - r01 + r00
                scale = np.maximum.reduce([np.abs(r00), np.abs(r10),
                                           np.abs(r01), np.abs(r11)])
                good = ((np.abs(marginal_i) + np.abs(marginal_j)) > 0.0) & (scale > 0.0)
                informative += int(np.sum(good))
                if np.any(good):
                    relative = np.abs(interaction[good]) / scale[good]
                    max_relative = max(max_relative, float(np.max(relative)))
                for alpha in alphas:
                    key = f"{alpha:.2f}"
                    hit = good & (np.abs(interaction) > alpha * scale)
                    count = int(np.sum(hit))
                    if count:
                        readings[key] += count
                        material[key] += int(np.sum(np.any(hit, axis=0)))
                        pair_mask[key] |= 1 << position
    return {
        "material_contrasts": material,
        "material_readings": readings,
        "pair_mask": pair_mask,
        "max_relative": max_relative,
        "informative_readings": informative,
    }


def base_coefficients(rng, points, base, amplitude):
    """Coefficients of one additive draw.

    The mean basis function is 1/sqrt(144), so the mean coefficient is scaled
    to make base the constant level of the surface.
    """
    coefficients = np.zeros((points, len(CELLS)))
    coefficients[:, 0] = base * np.sqrt(len(CELLS))
    for axis in range(6):
        for channel in range(1, LEVELS[axis]):
            column = component_index([channel if a == axis else 0 for a in range(6)])
            coefficients[:, column] = rng.uniform(-amplitude, amplitude, size=points)
    return coefficients


def support_channels(rng, support):
    channels = []
    for axis in support:
        channels.append(1 if LEVELS[axis] == 2 else int(rng.integers(1, LEVELS[axis])))
    return channels


def expected_pairs(surface, support):
    if surface == "pair":
        return 1 << PAIR_INDEX[tuple(support)]
    if surface == "triple":
        mask = 0
        for pair in itertools.combinations(support, 2):
            mask |= 1 << PAIR_INDEX[tuple(pair)]
        return mask
    return 0


def quantiles(values):
    if not values:
        return {"n": 0, "min": None, "q05": None, "median": None,
                "q95": None, "max": None}
    array = np.array(values, dtype=float)
    return {
        "n": int(array.size),
        "min": float(np.min(array)),
        "q05": float(np.quantile(array, 0.05)),
        "median": float(np.median(array)),
        "q95": float(np.quantile(array, 0.95)),
        "max": float(np.max(array)),
    }


def aggregate_statistics(seeds, abs_result, signed_result, expected_masks):
    """Rates and support checks over seeds for one surface."""
    out = {}
    for alpha in ALPHAS:
        key = f"{alpha:.2f}"
        abs_flags = [int(result["material_contrasts"][key] > 0) for result in abs_result]
        signed_flags = [int(result["material_contrasts"][key] > 0)
                        for result in signed_result]
        abs_masks = [result["pair_mask"][key] for result in abs_result]
        signed_masks = [result["pair_mask"][key] for result in signed_result]
        abs_outside = sum(1 for mask, expected in zip(abs_masks, expected_masks)
                          if mask & ~expected)
        signed_outside = sum(1 for mask, expected in zip(signed_masks, expected_masks)
                             if mask & ~expected)
        abs_inside = sum(1 for mask, expected in zip(abs_masks, expected_masks)
                         if mask & expected)
        signed_inside = sum(1 for mask, expected in zip(signed_masks, expected_masks)
                            if mask & expected)
        oracle_mismatch = sum(
            1 for result, flag in zip(signed_result, signed_flags)
            if (result["max_relative"] > alpha) != bool(flag)
        )
        abs_false_negative = sum(
            1 for abs_flag, signed_flag in zip(abs_flags, signed_flags)
            if signed_flag and not abs_flag
        )
        out[key] = {
            "abs_material_seeds": int(sum(abs_flags)),
            "abs_material_rate": float(np.mean(abs_flags)) if seeds else None,
            "signed_material_seeds": int(sum(signed_flags)),
            "signed_material_rate": float(np.mean(signed_flags)) if seeds else None,
            "abs_material_contrasts_mean": float(np.mean(
                [result["material_contrasts"][key] for result in abs_result])),
            "abs_material_contrasts_total": int(sum(
                result["material_contrasts"][key] for result in abs_result)),
            "signed_material_contrasts_mean": float(np.mean(
                [result["material_contrasts"][key] for result in signed_result])),
            "signed_material_contrasts_total": int(sum(
                result["material_contrasts"][key] for result in signed_result)),
            "signed_material_readings_mean": float(np.mean(
                [result["material_readings"][key] for result in signed_result])),
            "abs_pair_outside_seeds": int(abs_outside),
            "abs_pair_outside_rate": float(abs_outside / seeds) if seeds else None,
            "signed_pair_outside_seeds": int(signed_outside),
            "signed_pair_outside_rate": float(signed_outside / seeds) if seeds else None,
            "abs_expected_pair_seeds": int(abs_inside),
            "abs_expected_pair_rate": float(abs_inside / seeds) if seeds else None,
            "signed_expected_pair_seeds": int(signed_inside),
            "signed_expected_pair_rate": float(signed_inside / seeds) if seeds else None,
            "signed_oracle_mismatches": int(oracle_mismatch),
            "abs_false_negative_given_signed": int(abs_false_negative),
            "abs_false_negative_rate_given_signed": (
                float(abs_false_negative / sum(signed_flags)) if sum(signed_flags) else None),
        }
    out["max_relative"] = {
        "abs": quantiles([result["max_relative"] for result in abs_result]),
        "signed": quantiles([result["max_relative"] for result in signed_result]),
    }
    return out


def aggregate_design(records, seed_ids):
    """Coverage and margin over seeds for one design and one fit domain."""
    violation = [record["violations"] for record in records]
    no_allowance = [record["violations_no_allowance"] for record in records]
    margins = [record["margin"] for record in records if record["margin"] is not None]
    covered = [record["covered_cells"] / record["held_cells"]
               for record in records if record["held_cells"]]
    bounds = [record["bound"] for record in records]
    allowances = [record["allowance"] for record in records]
    worst_index = max((index for index, record in enumerate(records)
                       if record["held_cells"]),
                      key=lambda index: records[index]["worst"], default=None)
    margin_index = None
    for index, record in enumerate(records):
        if record["margin"] is None:
            continue
        if margin_index is None or record["margin"] < records[margin_index]["margin"]:
            margin_index = index
    first = records[0]
    return {
        "cells": first["cells"],
        "held_cells": first["held_cells"],
        "readings": first["readings"],
        "seeds": len(records),
        "violation_seeds": sum(1 for count in violation if count > 0),
        "violation_rate": (sum(1 for count in violation if count > 0) / len(records)
                           if records else None),
        "violations_total": int(sum(violation)),
        "no_allowance_seeds": sum(1 for count in no_allowance if count > 0),
        "no_allowance_rate": (sum(1 for count in no_allowance if count > 0) / len(records)
                              if records else None),
        "no_allowance_total": int(sum(no_allowance)),
        "violations_mean": float(np.mean(violation)) if records else None,
        "covered_fraction_mean": float(np.mean(covered)) if covered else None,
        "worst_max": (records[worst_index]["worst"]
                      if worst_index is not None else None),
        "worst_case": (None if worst_index is None else {
            "seed": seed_ids[worst_index],
            "ratio": records[worst_index]["worst"],
        }),
        "closest_margin": (None if margin_index is None else {
            "seed": seed_ids[margin_index],
            "margin": records[margin_index]["margin"],
        }),
        "margin": quantiles(margins),
        "bound": quantiles(bounds),
        "allowance": quantiles(allowances),
    }


def run_surface(seed_values, seed_ids, surface, magnitude, expected_masks, alphas):
    results = []
    for values in seed_values:
        results.append({
            "abs": absolute_statistic(values, alphas),
            "signed": signed_statistic(values, alphas),
        })
    abs_result = [result["abs"] for result in results]
    signed_result = [result["signed"] for result in results]
    statistics = aggregate_statistics(len(seed_values), abs_result, signed_result,
                                      expected_masks)
    designs = {}
    for name, entry in DESIGN_MAPS.items():
        signed_records = []
        absolute_records = []
        for values in seed_values:
            values_abs = np.abs(values)
            bound, allowance, _ = envelope(values, entry)
            signed_records.append(score_design(values_abs, entry, bound, allowance))
            bound, allowance, _ = envelope(values_abs, entry)
            absolute_records.append(score_design(values_abs, entry, bound, allowance))
        designs[name] = {
            "signed": dict(aggregate_design(signed_records, seed_ids),
                           components=entry["components"], order=entry["order"]),
            "absolute": dict(aggregate_design(absolute_records, seed_ids),
                             components=entry["components"], order=entry["order"]),
        }
    return {
        "surface": surface,
        "magnitude": magnitude,
        "expected_pair_masks": expected_masks,
        "statistics": statistics,
        "pair_masks_at_primary": {
            "alpha": PRIMARY,
            "abs": [result["pair_mask"][f"{PRIMARY:.2f}"] for result in abs_result],
            "signed": [result["pair_mask"][f"{PRIMARY:.2f}"] for result in signed_result],
        },
        "max_relative_per_seed": {
            "abs": [result["max_relative"] for result in abs_result],
            "signed": [result["max_relative"] for result in signed_result],
        },
        "designs": designs,
    }


def self_test(points=4):
    """Checks of the basis, the statistics and the design ladder."""
    out = {}
    out["basis_orthonormality_error"] = float(np.max(np.abs(PHI.T @ PHI - np.eye(len(CELLS)))))
    sizes = [DESIGN_MAPS[name]["cells"] for name in ("d1", "d9", "d35", "d79", "d144")]
    components = [DESIGN_MAPS[name]["components"] for name in ("d1", "d9", "d35", "d79", "d144")]
    out["design_cells"] = sizes
    out["design_components"] = components
    assert sizes == [1, 9, 35, 79, 144], sizes
    assert components == sizes, components

    # The referee's 2x2 example: signed errors 0, a, b, a+b with a = 1,
    # b = -0.6 give a zero signed contrast and an absolute cell-max contrast
    # of -1.2, relative 1.2. The surface depends only on axes 0 and 1, so the
    # other axes contribute no contrast.
    example_values = {(0, 0): 0.0, (1, 0): 1.0, (0, 1): -0.6, (1, 1): 0.4}
    example = np.zeros((1, len(CELLS)))
    for index, cell in enumerate(CELLS):
        example[0, index] = example_values[(cell[0], cell[1])]
    signed_example = signed_statistic(example, ALPHAS)
    abs_example = absolute_statistic(example, ALPHAS)
    out["referee_example"] = {
        "signed_max_relative": signed_example["max_relative"],
        "absolute_max_relative": abs_example["max_relative"],
        "signed_material_at_10": signed_example["material_contrasts"]["0.10"],
        "absolute_material_at_10": abs_example["material_contrasts"]["0.10"],
    }
    assert abs(signed_example["max_relative"]) < 1e-14
    assert abs(abs_example["max_relative"] - 1.2) < 1e-12
    assert signed_example["material_contrasts"]["0.10"] == 0
    assert abs_example["material_contrasts"]["0.10"] > 0

    rng = np.random.default_rng(731)
    base = base_coefficients(rng, points, 1.0, 0.1)
    additive = base @ PHI.T
    additive_signed = signed_statistic(additive, ALPHAS)
    additive_abs = absolute_statistic(additive, ALPHAS)
    out["additive"] = {
        "signed_max_relative": additive_signed["max_relative"],
        "absolute_max_relative": additive_abs["max_relative"],
        "signed_material_at_10": additive_signed["material_contrasts"]["0.10"],
        "absolute_material_at_10": additive_abs["material_contrasts"]["0.10"],
    }
    assert additive_signed["max_relative"] < 1e-13
    assert additive_signed["material_contrasts"]["0.10"] == 0

    for name in ("d9", "d35", "d79"):
        bound, allowance, _ = envelope(additive, DESIGN_MAPS[name])
        record = score_design(np.abs(additive), DESIGN_MAPS[name], bound, allowance)
        out[f"additive_{name}"] = {"allowance": allowance, "violations": record["violations"]}

    # The same additive construction with no common term: the signed errors
    # cross zero, so |e| is not e pointwise and the two fit domains diverge.
    rng = np.random.default_rng(733)
    centred = base_coefficients(rng, points, 0.0, 0.5) @ PHI.T
    centred_signed = signed_statistic(centred, ALPHAS)
    centred_absolute = absolute_statistic(centred, ALPHAS)
    signed_bound, signed_allowance, _ = envelope(centred, DESIGN_MAPS["d35"])
    absolute_bound, absolute_allowance, _ = envelope(np.abs(centred), DESIGN_MAPS["d35"])
    out["sign_changing_additive"] = {
        "negative_fraction": float(np.mean(centred < 0)),
        "signed_max_relative": centred_signed["max_relative"],
        "absolute_max_relative": centred_absolute["max_relative"],
        "signed_material_at_10": centred_signed["material_contrasts"]["0.10"],
        "absolute_material_at_10": centred_absolute["material_contrasts"]["0.10"],
        "d35_signed_bound": signed_bound,
        "d35_signed_allowance": signed_allowance,
        "d35_absolute_bound": absolute_bound,
        "d35_absolute_allowance": absolute_allowance,
    }
    assert out["sign_changing_additive"]["negative_fraction"] > 0.2
    assert centred_signed["max_relative"] < 1e-12
    assert centred_signed["material_contrasts"]["0.10"] == 0
    assert abs(signed_allowance - absolute_allowance) > 1e-6

    pair = base.copy()
    pair[:, component_index([1, 1, 0, 0, 0, 0])] += 0.8
    pair_values = pair @ PHI.T
    pair_signed = signed_statistic(pair_values, ALPHAS)
    out["pair"] = {
        "signed_max_relative": pair_signed["max_relative"],
        "signed_material_at_10": pair_signed["material_contrasts"]["0.10"],
        "pair_mask": pair_signed["pair_mask"]["0.10"],
    }
    assert pair_signed["max_relative"] > 0.1
    assert pair_signed["material_contrasts"]["0.10"] > 0
    assert pair_signed["pair_mask"]["0.10"] == 1 << PAIR_INDEX[(0, 1)]
    for name in ("d35", "d79"):
        bound, allowance, _ = envelope(pair_values, DESIGN_MAPS[name])
        record = score_design(np.abs(pair_values), DESIGN_MAPS[name], bound, allowance)
        out[f"pair_{name}"] = {"allowance": allowance, "violations": record["violations"]}
    bound, allowance, _ = envelope(pair_values, DESIGN_MAPS["d9"])
    record = score_design(np.abs(pair_values), DESIGN_MAPS["d9"], bound, allowance)
    out["pair_d9"] = {"allowance": allowance, "violations": record["violations"]}
    assert allowance > 0.0

    triple = base.copy()
    triple[:, component_index([1, 1, 1, 0, 0, 0])] += 0.8
    triple_values = triple @ PHI.T
    triple_signed = signed_statistic(triple_values, ALPHAS)
    out["triple"] = {
        "signed_max_relative": triple_signed["max_relative"],
        "signed_material_at_10": triple_signed["material_contrasts"]["0.10"],
        "pair_mask": triple_signed["pair_mask"]["0.10"],
    }
    assert triple_signed["max_relative"] > 0.1
    bound, allowance, _ = envelope(triple_values, DESIGN_MAPS["d35"])
    record = score_design(np.abs(triple_values), DESIGN_MAPS["d35"], bound, allowance)
    out["triple_d35"] = {"allowance": allowance, "violations": record["violations"]}
    assert allowance > 0.0
    bound, allowance, _ = envelope(triple_values, DESIGN_MAPS["d79"])
    record = score_design(np.abs(triple_values), DESIGN_MAPS["d79"], bound, allowance)
    out["triple_d79"] = {"allowance": allowance, "violations": record["violations"]}
    return out


def script_digest():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def prediction_tests(results, magnitudes):
    """Test the three stated predictions against the measured rates.

    The criteria are fixed here and stated in the record, so the verdicts are
    not chosen after seeing the numbers. False positive means a material
    verdict on an additive signed surface; coverage means the held-out
    absolute errors do not exceed the design bound.
    """
    alpha_keys = [f"{alpha:.2f}" for alpha in ALPHAS]
    primary = f"{PRIMARY:.2f}"
    p1_observed = {}
    p2_observed = {}
    for regime, data in results.items():
        statistics = data["additive"]["statistics"]
        p1_observed[regime] = {key: statistics[key]["abs_material_rate"]
                               for key in alpha_keys}
        p2_observed[regime] = {key: statistics[key]["signed_material_rate"]
                               for key in alpha_keys}

    p1_per_regime = {}
    for regime, rates in p1_observed.items():
        if rates[primary] >= 0.5:
            p1_per_regime[regime] = "supported"
        elif any(rate > 0 for rate in rates.values()):
            p1_per_regime[regime] = "partial"
        else:
            p1_per_regime[regime] = "refuted"
    if "supported" in p1_per_regime.values():
        p1_overall = "supported"
    elif "partial" in p1_per_regime.values():
        p1_overall = "partial"
    else:
        p1_overall = "refuted"

    p2_per_regime = {}
    for regime, rates in p2_observed.items():
        p2_per_regime[regime] = (
            "supported" if all(rate == 0 for rate in rates.values()) else "refuted")
    p2_overall = "supported" if all(value == "supported"
                                    for value in p2_per_regime.values()) else "refuted"

    p3_observed = {}
    for regime, data in results.items():
        p3_observed[regime] = {}
        for magnitude in magnitudes:
            key = f"{magnitude:.4f}"
            entry = {}
            for fit in ("signed", "absolute"):
                pair = data["pair"][key]["designs"]["d35"][fit]
                triple = data["triple"][key]["designs"]["d35"][fit]
                entry[fit] = {
                    "pair_violation_rate": pair["violation_rate"],
                    "triple_violation_rate": triple["violation_rate"],
                    "pair_margin_median": pair["margin"]["median"],
                    "triple_margin_median": triple["margin"]["median"],
                }
            p3_observed[regime][key] = entry

    largest = f"{max(magnitudes):.4f}"
    pair_rates = {regime: p3_observed[regime][largest]["signed"]["pair_violation_rate"]
                  for regime in p3_observed}
    triple_rates = {regime: p3_observed[regime][largest]["signed"]["triple_violation_rate"]
                    for regime in p3_observed}
    if all(rate == 0.0 for rate in pair_rates.values()) and any(
            rate > 0.0 for rate in triple_rates.values()):
        p3_overall = "supported"
    elif not any(rate > 0.0 for rate in triple_rates.values()):
        p3_overall = "refuted"
    else:
        p3_overall = "mixed"
    p3_per_regime = {
        regime: {
            "pair_violation_rate": pair_rates[regime],
            "triple_violation_rate": triple_rates[regime],
            "pair_covered": pair_rates[regime] == 0.0,
            "triple_covered": triple_rates[regime] == 0.0,
        }
        for regime in p3_observed
    }

    return {
        "p1_absolute_manufactures": {
            "statement": "the absolute-value cell-maximum statistic reports material "
                         "interactions on a signed-additive surface",
            "observed_abs_false_positive_rate": p1_observed,
            "observed_signed_false_positive_rate": p2_observed,
            "criterion": "per regime: supported if the primary-margin absolute rate is "
                         "at least 0.5, partial if any margin is positive, refuted if "
                         "every rate is zero; overall takes the strongest regime",
            "per_regime": p1_per_regime,
            "overall": p1_overall,
        },
        "p2_signed_abstains": {
            "statement": "the signed pointwise statistic abstains on a signed-additive "
                         "surface",
            "observed_signed_false_positive_rate": p2_observed,
            "criterion": "supported if the signed rate is zero at every margin in every "
                         "regime",
            "per_regime": p2_per_regime,
            "overall": p2_overall,
        },
        "p3_design_35_covers_pair_not_triple": {
            "statement": "the 35-cell design covers a one-pairwise surface but not a "
                         "three-factor one",
            "observed": p3_observed,
            "criterion": "at the largest magnitude and the signed fit: supported if the "
                         "pair violation rate is zero in every regime and the triple rate "
                         "is positive in at least one regime; refuted if the triple rate is "
                         "zero in every regime; mixed otherwise",
            "per_regime": p3_per_regime,
            "overall": p3_overall,
        },
    }


def mask_pairs(mask):
    """The axis pairs whose bit is set in a pair mask."""
    return [list(PAIRS[index]) for index in range(len(PAIRS))
            if mask & (1 << index)]


def coverage_curve(results, fit):
    """Violation rate and closest margin by design, class and fit domain."""
    curve = {}
    for design in ("d1", "d9", "d35", "d79", "d144"):
        entry = {
            "cells": DESIGN_MAPS[design]["cells"],
            "components": DESIGN_MAPS[design]["components"],
            "classes": {},
        }
        for surf in ("additive", "pair", "triple"):
            blocks = []
            for regime, data in results.items():
                if surf == "additive":
                    blocks.append(data["additive"]["designs"][design][fit])
                else:
                    blocks.extend(data[surf][magnitude]["designs"][design][fit]
                                  for magnitude in data[surf])
            surfaces = sum(block["seeds"] for block in blocks)
            violating = sum(block["violation_seeds"] for block in blocks)
            readings = sum(block["readings"] * block["seeds"] for block in blocks)
            margins = [block["margin"]["min"] for block in blocks
                       if block["margin"]["min"] is not None]
            entry["classes"][surf] = {
                "surfaces": surfaces,
                "violating_surfaces": violating,
                "violation_rate": violating / surfaces if surfaces else None,
                "readings": readings,
                "closest_margin": min(margins) if margins else None,
            }
        curve[design] = entry
    return curve


def decision_summary(results, seeds):
    """Score the two portfolio thresholds on this run's own numbers.

    Threshold A: no false material absolute-value interaction at the primary
    margin on the known-additive surfaces. Threshold B: no held-out violation
    of the 35-cell design with its measured allowance, and every margin at
    least one, on every synthetic surface. Each threshold names its worst
    case, so one failure cannot hide inside an average.
    """
    primary = f"{PRIMARY:.2f}"
    per_regime = {}
    for regime, data in results.items():
        statistics = data["additive"]["statistics"][primary]
        abs_relative = data["additive"]["max_relative_per_seed"]["abs"]
        signed_relative = data["additive"]["max_relative_per_seed"]["signed"]
        abs_mask = data["additive"]["pair_masks_at_primary"]["abs"]
        seed_ids = [support["seed"] for support in data["supports"]]
        abs_index = max(range(len(abs_relative)), key=lambda i: abs_relative[i])
        signed_index = max(range(len(signed_relative)),
                           key=lambda i: signed_relative[i])
        per_regime[regime] = {
            "surfaces": len(seed_ids),
            "abs_material_seeds": statistics["abs_material_seeds"],
            "abs_material_contrasts": statistics["abs_material_contrasts_total"],
            "abs_material_rate": statistics["abs_material_rate"],
            "signed_material_seeds": statistics["signed_material_seeds"],
            "signed_material_contrasts": statistics["signed_material_contrasts_total"],
            "signed_material_rate": statistics["signed_material_rate"],
            "worst_abs": {
                "seed": seed_ids[abs_index],
                "max_relative": abs_relative[abs_index],
                "pairs": mask_pairs(abs_mask[abs_index]),
            },
            "worst_signed": {
                "seed": seed_ids[signed_index],
                "max_relative": signed_relative[signed_index],
            },
        }
    surfaces = seeds * len(results)
    abs_seeds = sum(entry["abs_material_seeds"] for entry in per_regime.values())
    abs_contrasts = sum(entry["abs_material_contrasts"] for entry in per_regime.values())
    signed_seeds = sum(entry["signed_material_seeds"] for entry in per_regime.values())
    signed_contrasts = sum(entry["signed_material_contrasts"]
                           for entry in per_regime.values())
    worst_regime = max(per_regime,
                       key=lambda name: per_regime[name]["worst_abs"]["max_relative"])
    threshold_a = {
        "statement": "zero false material absolute-value interactions at the "
                     "primary margin on the known-additive surfaces",
        "surfaces": surfaces,
        "abs_material_seeds": abs_seeds,
        "abs_material_contrasts": abs_contrasts,
        "abs_material_rate": abs_seeds / surfaces if surfaces else None,
        "signed_material_seeds": signed_seeds,
        "signed_material_contrasts": signed_contrasts,
        "signed_material_rate": signed_seeds / surfaces if surfaces else None,
        "met": abs_seeds == 0 and abs_contrasts == 0,
        "per_regime": per_regime,
        "worst_case": {
            "regime": worst_regime,
            **per_regime[worst_regime]["worst_abs"],
        },
    }

    def d35_blocks(fit):
        for regime, data in results.items():
            yield regime, "additive", "", data["additive"]["designs"]["d35"][fit]
            for surf in ("pair", "triple"):
                for magnitude in MAGNITUDES:
                    key = f"{magnitude:.4f}"
                    yield regime, surf, key, data[surf][key]["designs"]["d35"][fit]

    threshold_b = {
        "statement": "zero held-out violations of the 35-cell design with its "
                     "measured allowance, and every margin at least one, on "
                     "every synthetic surface",
        "per_fit": {},
    }
    for fit in ("absolute", "signed"):
        classes = {
            name: {
                "surfaces": 0,
                "violating_surfaces": 0,
                "violations": 0,
                "readings": 0,
                "worst": {"ratio": None},
                "closest": {"margin": None},
            }
            for name in ("additive", "pair", "triple")
        }
        for regime, surf, magnitude, block in d35_blocks(fit):
            entry = classes[surf]
            entry["surfaces"] += block["seeds"]
            entry["violating_surfaces"] += block["violation_seeds"]
            entry["violations"] += block["violations_total"]
            entry["readings"] += block["readings"] * block["seeds"]
            worst = block["worst_max"]
            if worst is not None and (entry["worst"]["ratio"] is None
                                      or worst > entry["worst"]["ratio"]):
                entry["worst"] = {
                    "ratio": worst,
                    "regime": regime,
                    "surface": surf,
                    "magnitude": magnitude,
                    "seed": block["worst_case"]["seed"],
                }
            margin = block["margin"]["min"]
            if margin is not None and (entry["closest"]["margin"] is None
                                       or margin < entry["closest"]["margin"]):
                entry["closest"] = {
                    "margin": margin,
                    "regime": regime,
                    "surface": surf,
                    "magnitude": magnitude,
                    "seed": block["closest_margin"]["seed"],
                }
        overall = {
            "surfaces": sum(entry["surfaces"] for entry in classes.values()),
            "violating_surfaces": sum(entry["violating_surfaces"]
                                      for entry in classes.values()),
            "violations": sum(entry["violations"] for entry in classes.values()),
            "readings": sum(entry["readings"] for entry in classes.values()),
            "worst": max((entry["worst"] for entry in classes.values()
                          if entry["worst"]["ratio"] is not None),
                         key=lambda item: item["ratio"], default=None),
            "closest": min((entry["closest"] for entry in classes.values()
                            if entry["closest"]["margin"] is not None),
                           key=lambda item: item["margin"], default=None),
        }
        overall["met"] = (
            overall["violations"] == 0
            and overall["closest"] is not None
            and overall["closest"]["margin"] >= 1.0
        )
        threshold_b["per_fit"][fit] = {"classes": classes, "overall": overall}
    threshold_b["met"] = all(entry["overall"]["met"]
                             for entry in threshold_b["per_fit"].values())

    return {
        "threshold_a_additive_abs_at_primary": threshold_a,
        "threshold_b_d35_coverage": threshold_b,
        "coverage_curve": {
            fit: coverage_curve(results, fit) for fit in ("signed", "absolute")
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=DEFAULT_SEEDS)
    parser.add_argument("--points", type=int, default=DEFAULT_POINTS)
    args = parser.parse_args()
    points = args.points
    seeds = list(range(args.seeds))

    checks = self_test()

    contrast_count = sum(
        len(index) * (LEVELS[i] - 1) * (LEVELS[j] - 1)
        for i, j, index in PAIR_BLOCKS)

    results = {}
    for regime_index, (regime_name, regime) in enumerate(REGIMES.items()):
        regime_result = {}
        additive = []
        pair_surfaces = {f"{magnitude:.4f}": [] for magnitude in MAGNITUDES}
        triple_surfaces = {f"{magnitude:.4f}": [] for magnitude in MAGNITUDES}
        supports = []
        seed_ids = [SEED_BASE + 10000 * regime_index + seed for seed in seeds]
        for seed, seed_id in zip(seeds, seed_ids):
            rng = np.random.default_rng(seed_id)
            base = base_coefficients(rng, points, regime["base"],
                                     regime["main_amplitude"])
            pair_support = PAIRS[seed % len(PAIRS)]
            triple_support = TRIPLES[seed % len(TRIPLES)]
            pair_channels = support_channels(rng, pair_support)
            triple_channels = support_channels(rng, triple_support)
            pair_column = component_index(
                [pair_channels[pair_support.index(axis)] if axis in pair_support else 0
                 for axis in range(6)])
            triple_column = component_index(
                [triple_channels[triple_support.index(axis)] if axis in triple_support else 0
                 for axis in range(6)])
            supports.append({
                "seed": seed_id,
                "pair": {"axes": list(pair_support), "channels": pair_channels},
                "triple": {"axes": list(triple_support), "channels": triple_channels},
            })
            additive.append(base @ PHI.T)
            for magnitude in MAGNITUDES:
                values = base.copy()
                values[:, pair_column] += magnitude
                pair_surfaces[f"{magnitude:.4f}"].append(values @ PHI.T)
                values = base.copy()
                values[:, triple_column] += magnitude
                triple_surfaces[f"{magnitude:.4f}"].append(values @ PHI.T)

        regime_result["additive"] = run_surface(
            additive, seed_ids, "additive", 0.0, [0] * len(seeds), ALPHAS)
        regime_result["pair"] = {}
        for magnitude in MAGNITUDES:
            key = f"{magnitude:.4f}"
            expected = [expected_pairs("pair", tuple(support["pair"]["axes"]))
                        for support in supports]
            regime_result["pair"][key] = run_surface(
                pair_surfaces[key], seed_ids, "pair", magnitude, expected, ALPHAS)
        regime_result["triple"] = {}
        for magnitude in MAGNITUDES:
            key = f"{magnitude:.4f}"
            expected = [expected_pairs("triple", tuple(support["triple"]["axes"]))
                        for support in supports]
            regime_result["triple"][key] = run_surface(
                triple_surfaces[key], seed_ids, "triple", magnitude, expected, ALPHAS)
        regime_result["supports"] = supports
        results[regime_name] = regime_result

    record = {
        "schema": "synthetic-surface/1",
        "command": f"python -I scripts/synthetic_surface.py --seeds {args.seeds} "
                   f"--points {points}",
        "script": {"path": "scripts/synthetic_surface.py", "sha256": script_digest()},
        "python": platform.python_version(),
        "numpy": np.__version__,
        "paths_relative_to": "papers/composed-bounds",
        "parameters": {
            "factors": FACTORS,
            "levels": LEVELS,
            "cells": len(CELLS),
            "points_per_surface": points,
            "seed_base": SEED_BASE,
            "seeds": args.seeds,
            "seed_rule": "seed = seed_base + 10000*regime_index + i, i in range(seeds)",
            "unit": "synthetic signed error in the paper's F0-scaled unit, F0 = 1",
            "alphas": list(ALPHAS),
            "primary_alpha": PRIMARY,
            "regimes": REGIMES,
            "magnitudes": list(MAGNITUDES),
            "contrast_rows": {str(level): [list(row) for row in rows]
                              for level, rows in CONTRASTS.items()},
            "designs": {name: {"cells": entry["cells"],
                               "order": entry["order"],
                               "components": entry["components"]}
                        for name, entry in DESIGN_MAPS.items()},
            "pair_contrasts_per_surface": contrast_count,
            "pairs": [list(pair) for pair in PAIRS],
            "triples": [list(triple) for triple in TRIPLES],
            "support_rule": "pair is PAIRS[seed % 15]; triple is TRIPLES[seed % 20]; "
                            "channels are uniform over the non-baseline rows",
        },
        "checks": checks,
        "results": results,
        "predictions": prediction_tests(results, MAGNITUDES),
        "decision": decision_summary(results, args.seeds),
    }
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
