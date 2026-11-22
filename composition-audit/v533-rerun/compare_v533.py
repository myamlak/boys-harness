#!/usr/bin/env python3
"""Compare the deposited v5.0.0 cube with the v5.3.3 re-measurement.

Reads the two deposited cell archives (../../records/cells) and the two
compressed re-measurement cell files (results-fused, results-separate), and
recomputes every comparison number summary.json states.  It exits 0 only when
all of them reproduce; a mismatch is printed and returns 1.

The signed-interaction reduction is the harness's own signed_interaction.py,
imported from the parent directory, so the contrast definition is the one the
deposited record uses.
"""
import csv
import itertools
import json
import lzma
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import signed_interaction as sig
import summarise_interaction as si

WITHIN = sig.WITHIN
AXES = sig.AXES
RESOLVING = list(sig.RESOLVING)
OLD = {
    "fused": HERE.parent / "records" / "cells" / "cells-msvc-fused.csv.xz",
    "separate": HERE.parent / "records" / "cells" / "cells-msvc-separate.csv.xz",
}
NEW = {
    "fused": HERE / "results-fused" / "cells.csv.xz",
    "separate": HERE / "results-separate" / "cells.csv.xz",
}
EXPECTED = {
    "fused": {
        "changed_cells": 7128, "identical_rows": 487080,
        "max_abs_value_difference": 1.1157741397482823e-14,
        "affected_combinations": 24,
        "material_old": {"0.05": 912, "0.1": 901, "0.2": 899},
        "material_new": {"0.05": 897, "0.1": 885, "0.2": 883},
        "fp64_old": {"0.05": 336, "0.1": 325, "0.2": 323},
        "fp64_new": {"0.05": 321, "0.1": 309, "0.2": 307},
        "fp32_old": {"0.05": 276, "0.1": 276, "0.2": 276},
        "fp32_new": {"0.05": 276, "0.1": 276, "0.2": 276},
        "max_statement_ratio": {"fp64": 0.47303043999577227,
                                "fp32": 0.8191245254432976,
                                "fp16": 0.9934027438291815,
                                "bf16": 0.9985908700338076},
    },
    "separate": {
        "changed_cells": 7128, "identical_rows": 487080,
        "max_abs_value_difference": 1.1157741397482823e-14,
        "affected_combinations": 24,
        "material_old": {"0.05": 914, "0.1": 903, "0.2": 901},
        "material_new": {"0.05": 895, "0.1": 891, "0.2": 885},
        "fp64_old": {"0.05": 338, "0.1": 327, "0.2": 325},
        "fp64_new": {"0.05": 319, "0.1": 315, "0.2": 309},
        "fp32_old": {"0.05": 276, "0.1": 276, "0.2": 276},
        "fp32_new": {"0.05": 276, "0.1": 276, "0.2": 276},
        "max_statement_ratio": {"fp64": 0.47303043999577227,
                                "fp32": 0.8191245254432976,
                                "fp16": 0.9934027438291815,
                                "bf16": 0.9985908700338076},
    },
    "contrast_counts": {
        "total": 1808, "informative": 1648, "per_lane": 664,
        "cross_lane": 480, "order_three_components": 44,
    },
}


def load(path):
    rows = 0
    values = {}
    with lzma.open(path, "rt", newline="") as f:
        for r in csv.DictReader(f):
            rows += 1
            if r["status"] != "ok":
                continue
            key = (r["precision"],) + tuple(r[a] for a in WITHIN) + (
                int(r["n"]), float(r["x_nominal"]))
            values[key] = (float(r["value"]), float(r["abs_err"]),
                           float(r["figure"]), float(r["half_digit"]),
                           float(r["limit"]), float(r["figure_ratio"]))
    return values, rows


def signed_cells(path):
    cells = defaultdict(dict)
    with lzma.open(path, "rt", newline="") as f:
        for r in csv.DictReader(f):
            lane = r["precision"]
            if lane not in sig.RESOLVING or r["status"] != "ok":
                continue
            key = tuple(r[a] for a in WITHIN)
            n = int(r["n"])
            xnom = float(r["x_nominal"])
            e = (float(r["value"]) - float(r["ref_hi"])) - float(r["ref_lo"])
            f0 = float(r["f0"])
            cells[(lane, key)][(n, xnom)] = (
                e / f0 if f0 > 0 else 0.0,
                abs(e) / f0 if f0 > 0 else 0.0,
                int(r["pack_live"]), int(r["partition_in_domain"]),
                r["argument_kind"])
    return cells


def material(view, axes):
    contrasts = sig.pair_contrasts(view, axes)
    informative = 0
    mats = {str(a): 0 for a in sig.ALPHAS}
    for c in contrasts:
        if c["informative_points"] > 0 and c["max_pointwise"] is not None:
            informative += 1
            rel = c["max_pointwise"]["relative"]
            for a in sig.ALPHAS:
                if rel > a:
                    mats[str(a)] += 1
    return contrasts, informative, mats


def combined(cells):
    return {(lane,) + key: pts for (lane, key), pts in cells.items()}


def lane(cells, name):
    return {key: pts for (row_lane, key), pts in cells.items() if row_lane == name}


def main():
    problems = []
    out = {"contrast_counts": {}}
    # Design counts: 664 per-lane contrasts, 480 crossing the precision axis,
    # 44 order-three components in the six-axis decomposition.
    meta = si.component_metadata([2, 2, 3, 2, 3, 2])
    order_three = sum(1 for m in meta if m["order"] == 3)
    out["contrast_counts"]["order_three_components"] = order_three
    if order_three != EXPECTED["contrast_counts"]["order_three_components"]:
        problems.append(f"order-three components {order_three}")

    old_all = {r: signed_cells(OLD[r]) for r in OLD}
    new_all = {r: signed_cells(NEW[r]) for r in NEW}
    for route in ("fused", "separate"):
        old, old_rows = load(OLD[route])
        new, new_rows = load(NEW[route])
        common = set(old) & set(new)
        changed = [k for k in common if old[k][0] != new[k][0]]
        changed_err = sum(1 for k in common if old[k][1] != new[k][1])
        changed_figure = sum(1 for k in common
                             if old[k][2:5] != new[k][2:5])
        maxdiff = max((abs(old[k][0] - new[k][0]) for k in changed), default=0.0)
        combos = Counter((k[0],) + k[1:7] for k in changed)
        points = sorted({k[8] for k in changed})
        orders = sorted({k[7] for k in changed})
        ratios = {name: max(v[5] for k, v in new.items() if k[0] == name)
                  for name in ("fp64", "fp32", "fp16", "bf16")}
        old_ratios = {name: max(v[5] for k, v in old.items() if k[0] == name)
                      for name in ("fp64", "fp32", "fp16", "bf16")}

        contrasts_all, informative_all, mats_old = material(
            combined(old_all[route]), AXES)
        _, _, mats_new = material(combined(new_all[route]), AXES)
        lane_mats = {}
        for name in RESOLVING:
            lane_contrasts, lane_inf, lo = material(
                lane(old_all[route], name), WITHIN)
            _, _, ln = material(lane(new_all[route], name), WITHIN)
            lane_mats[name] = {"old": lo, "new": ln, "informative": lane_inf,
                               "contrasts": len(lane_contrasts)}
        cross = sum(1 for c in contrasts_all if "precision" in c["pair"])
        per_lane = lane_mats["fp64"]["contrasts"]
        out["contrast_counts"].update({
            "total": len(contrasts_all),
            "informative": informative_all,
            "per_lane": per_lane,
            "cross_lane": cross,
        })

        block = {
            "rows_old": old_rows,
            "rows_new": new_rows,
            "cells_compared": len(common),
            "changed_cells": len(changed),
            "changed_abs_err": changed_err,
            "changed_figure_half_digit_limit": changed_figure,
            "identical_rows": len(common) - len(changed),
            "max_abs_value_difference": maxdiff,
            "affected_combinations": len(combos),
            "affected_lane": sorted({c[0] for c in combos}),
            "affected_route": sorted({c[1] for c in combos}),
            "affected_partition": sorted({c[3] for c in combos}),
            "affected_schemes": sorted({c[2] for c in combos}),
            "affected_packs": sorted({c[4] for c in combos}),
            "affected_forms": sorted({c[5] for c in combos}),
            "affected_exps": sorted({c[6] for c in combos}),
            "changed_arguments": points,
            "changed_orders": {"count": len(orders),
                               "first": orders[0] if orders else None,
                               "last": orders[-1] if orders else None},
            "material_old": mats_old,
            "material_new": mats_new,
            "per_lane_material": lane_mats,
            "max_statement_ratio_new": ratios,
            "max_statement_ratio_old": old_ratios,
        }
        out[route] = block

        exp = EXPECTED[route]
        if len(changed) != exp["changed_cells"]:
            problems.append(f"{route}: changed cells {len(changed)}")
        if block["identical_rows"] != exp["identical_rows"]:
            problems.append(f"{route}: identical rows {block['identical_rows']}")
        if maxdiff != exp["max_abs_value_difference"]:
            problems.append(f"{route}: max diff {maxdiff!r}")
        if len(combos) != exp["affected_combinations"]:
            problems.append(f"{route}: combinations {len(combos)}")
        if mats_old != exp["material_old"]:
            problems.append(f"{route}: material old {mats_old}")
        if mats_new != exp["material_new"]:
            problems.append(f"{route}: material new {mats_new}")
        for name in ("fp64", "fp32"):
            if lane_mats[name]["old"] != exp[f"{name}_old"]:
                problems.append(f"{route}: {name} old {lane_mats[name]['old']}")
            if lane_mats[name]["new"] != exp[f"{name}_new"]:
                problems.append(f"{route}: {name} new {lane_mats[name]['new']}")
        if ratios != exp["max_statement_ratio"]:
            problems.append(f"{route}: maxima {ratios}")

    cc = out["contrast_counts"]
    if cc["total"] != EXPECTED["contrast_counts"]["total"]:
        problems.append(f"contrasts {cc['total']}")
    if cc["informative"] != EXPECTED["contrast_counts"]["informative"]:
        problems.append(f"informative {cc['informative']}")
    if cc["per_lane"] != EXPECTED["contrast_counts"]["per_lane"]:
        problems.append(f"per-lane contrasts {cc['per_lane']}")
    if cc["cross_lane"] != EXPECTED["contrast_counts"]["cross_lane"]:
        problems.append(f"cross-lane contrasts {cc['cross_lane']}")

    out["problems"] = problems
    print(json.dumps(out, indent=1))
    if problems:
        print("MISMATCH:", "; ".join(problems), file=sys.stderr)
        return 1
    print("all comparison numbers reproduce", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
