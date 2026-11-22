#!/usr/bin/env python3
"""Distribution of the signed interaction-to-response ratio.

Recomputes the 1,808 four-cell contrasts from the deposited v5.0.0 cube
(../../records/cells) with the harness's own pair_contrasts(), and writes
b7_ratio_distribution.json beside this script.  It reports how saturated the
materiality counts are: how many informative contrasts sit above one, the
quantiles, and how many move between the 5, 10 and 20 per cent margins.
"""
import csv
import json
import lzma
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import signed_interaction as sig

FILES = {
    "msvc-fused": HERE.parent / "records" / "cells" / "cells-msvc-fused.csv.xz",
    "msvc-separate": HERE.parent / "records" / "cells" / "cells-msvc-separate.csv.xz",
}
RESOLVING = set(sig.RESOLVING)
WITHIN = sig.WITHIN


def load(path):
    cells = defaultdict(dict)
    rows = 0
    with lzma.open(path, "rt", newline="") as f:
        for raw in csv.DictReader(f):
            rows += 1
            lane = raw["precision"]
            if lane not in RESOLVING or raw["status"] != "ok":
                continue
            key = tuple(raw[a] for a in WITHIN)
            n = int(raw["n"])
            xnom = float(raw["x_nominal"])
            value = float(raw["value"])
            hi = float(raw["ref_hi"])
            lo = float(raw["ref_lo"])
            f0 = float(raw["f0"])
            signed = (value - hi) - lo
            sr = signed / f0 if f0 > 0.0 else 0.0
            ar = abs(signed) / f0 if f0 > 0.0 else 0.0
            cells[(lane, key)][(n, xnom)] = (
                sr, ar, int(raw["pack_live"]), int(raw["partition_in_domain"]),
                raw["argument_kind"])
    return cells, rows


def quantiles(values, ps):
    v = sorted(values)
    out = {}
    for p in ps:
        if not v:
            out[str(p)] = None
            continue
        idx = (len(v) - 1) * p
        lo = int(idx)
        hi = min(lo + 1, len(v) - 1)
        frac = idx - lo
        out[str(p)] = v[lo] * (1 - frac) + v[hi] * frac
    return out


result = {}
for cfg, path in FILES.items():
    cells, rows = load(path)
    view = {(lane,) + key: pts for (lane, key), pts in cells.items()}
    contrasts = sig.pair_contrasts(view, sig.AXES)
    rels = []
    paper_rels = []
    informative = 0
    for c in contrasts:
        if c["informative_points"] > 0 and c["max_pointwise"] is not None:
            informative += 1
            rels.append(c["max_pointwise"]["relative"])
            if c["paper_cell_max"] is not None:
                paper_rels.append(c["paper_cell_max"]["relative"])
    mats = {a: sum(1 for r in rels if r > a) for a in sig.ALPHAS}
    bands = {
        "<=0.05": sum(1 for r in rels if r <= 0.05),
        "(0.05,0.10]": sum(1 for r in rels if 0.05 < r <= 0.10),
        "(0.10,0.20]": sum(1 for r in rels if 0.10 < r <= 0.20),
        "(0.20,0.50]": sum(1 for r in rels if 0.20 < r <= 0.50),
        "(0.50,1.00]": sum(1 for r in rels if 0.50 < r <= 1.00),
        ">1": sum(1 for r in rels if r > 1.0),
    }
    result[cfg] = {
        "rows_read": rows,
        "contrasts": len(contrasts),
        "informative": informative,
        "material": {str(a): mats[a] for a in sig.ALPHAS},
        "bands": bands,
        "quantiles": quantiles(rels, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0]),
        "share_gt_1": bands[">1"] / len(rels) if rels else None,
        "paper_cell_max_quantiles": quantiles(paper_rels, [0, 0.5, 0.9, 0.99, 1.0]),
        "paper_cell_max_gt_1": sum(1 for r in paper_rels if r > 1.0),
    }
    print(cfg, json.dumps(result[cfg], indent=1))

with open(HERE / "b7_ratio_distribution.json", "w", newline="") as f:
    json.dump(result, f, indent=1)
    f.write("\n")
