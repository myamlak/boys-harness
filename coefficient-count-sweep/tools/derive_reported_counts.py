#!/usr/bin/env python3
"""Derive the paper's reported coefficient counts from the 60-digit audit.

The audit record (``results-corrected/corrected-counts-60dig.json``) stores,
for every one of the 45 Group A cells, three numbers: the nominal count
``(m + 1) + k`` of the recorded pair, the full stored-witness count
``len(P) + len(Q) - 1`` of that same pair, and a re-solve at the control's
count ``count60`` that sweeps every admissible ``(m, k)`` pair.

This script applies one rule to those measurements and writes the derived
record ``results-corrected/reported-counts.json``.  The rule is stated once,
here and in the record's ``definition`` field:

* If the recorded pair's nominal truncation meets the cell target on both
  grids (2,001 certifying and 4,001 checking nodes) under both evaluations
  (60 decimal digits and binary64), the reported count is that nominal count
  and the mode is ``nominal``.
* Otherwise, if the re-solve at ``count60`` exhibits a pair whose nominal
  truncation meets the target on both grids under at least one of the two
  evaluations, the reported count is ``count60`` and the mode is ``genuine``;
  the pairs are recorded, binary64 first because binary64 is the paper's
  acceptance rule.
* Otherwise no genuine fit exists at that count, and the reported count is
  the full stored witness of the recorded pair (the solver zero-pads both
  coefficient arrays, so the stored representation is longer than its
  nominal truncation); the mode is ``full``.

The script verifies the audit record's digest, checks the resulting
aggregate against the constants the revision was verified with, and prints a
one-line manifest.  Stdlib only; no network, no git.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

SOURCE = "results-corrected/corrected-counts-60dig.json"
LEGACY = "results-corrected/corrected-counts.json"
OUTPUT = "results-corrected/reported-counts.json"

# The audit record this derivation was verified against.
SOURCE_SHA256 = (
    "18a41a42798f0cd1fead4a143309ded3112d1bba94727b093bf27df02e3256e7")

POW = {
    "m7": "2^-7", "m6": "2^-6", "m5": "2^-5", "m4": "2^-4",
    "m3": "2^-3", "m2": "2^-2", "m1": "2^-1", "p0": "2^0",
    "p1": "2^1", "p2": "2^2", "p3": "2^3", "p4": "2^4", "p5": "2^5",
}
TARGETS = ("rel:1e-13", "abs:5e-14", "rel:1e-14")

DEFINITION = (
    "Reported stored coefficient count: the nominal count (m+1)+k of the "
    "recorded pair when that truncation meets the cell target on both the "
    "2,001-node certifying and 4,001-node checking grids under both the "
    "60-digit and binary64 evaluations (mode 'nominal', 40 cells); "
    "otherwise the nominal count of a genuine re-solve at the control's "
    "count count60 whose truncation meets the target on both grids under "
    "binary64 or 60-digit evaluation (mode 'genuine', 1 cell); otherwise "
    "the full stored witness len(P)+len(Q)-1 of the recorded pair, whose "
    "nominal truncation does not meet the target and for which no genuine "
    "re-solve exists at count60 (mode 'full', 4 cells). Each stored "
    "coefficient counts once; a cell's pair is the pair that meets."
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _meets(block, which: str, grid: int, base: int) -> bool:
    return bool(block[which]["grid_{}".format(grid)]["meets{}".format(base)])


def classify(cell):
    """The rule, applied to one audit cell.  Returns (mode, count, pair)."""
    fit = cell["measurements"]["recorded_fit"]
    nominal = all(
        _meets(fit, "nominal", grid, base)
        for grid in (2001, 4001) for base in (60, 64))
    if nominal:
        return "nominal", int(cell["nominal_count"]), list(
            cell["recorded_pair"])
    level = cell["measurements"].get("count60_level") or {}
    checks = level.get("checks") or []
    good64 = [list(check["pair"]) for check in checks
              if check.get("nominal_meets64_2001")
              and check.get("nominal_meets64_4001")]
    good60 = [list(check["pair"]) for check in checks
              if check.get("nominal_meets60_2001")
              and check.get("nominal_meets60_4001")]
    if good64 or good60:
        pair = good64[0] if good64 else good60[0]
        return "genuine", int(level.get("count60")
                              or cell["nominal_count"]), pair
    return "full", int(fit["full_count"]), list(cell["recorded_pair"])


def main(argv):
    source_path = ROOT / SOURCE
    digest = _sha256(source_path)
    if digest != SOURCE_SHA256:
        raise SystemExit(
            "audit record digest {} != expected {}".format(
                digest, SOURCE_SHA256))
    audit = json.loads(source_path.read_text(encoding="utf-8"))

    legacy = {}
    legacy_path = ROOT / LEGACY
    if legacy_path.exists():
        old = json.loads(legacy_path.read_text(encoding="utf-8"))
        legacy = {(cell["target"], cell["width_name"]): cell
                  for cell in old["cells"]}

    cells = []
    for cell in audit["cells"]:
        measurement = cell["measurements"]
        fit = measurement["recorded_fit"]
        mode, count, pair = classify(cell)
        key = (cell["target"], cell["width_name"])
        item = {
            "index": int(cell["index"]),
            "target": cell["target"],
            "target_kind": cell["target_kind"],
            "target_value": cell["target_value"],
            "width": cell["width"],
            "width_name": cell["width_name"],
            "nominal_count": int(cell["nominal_count"]),
            "full_count": int(fit["full_count"]),
            "chebyshev_count": int(
                measurement["chebyshev"]["count_new"]),
            "recorded_pair": list(cell["recorded_pair"]),
            "count_mode": mode,
            "reported_count": count,
            "reported_pair": pair,
            "delivered_2001": cell["delivered_2001"],
            "delivered_4001": cell["delivered_4001"],
        }
        level = measurement.get("count60_level") or {}
        if level:
            item["count60"] = level.get("count60")
            item["genuine_pairs_64"] = [
                list(check["pair"]) for check in level.get("checks", [])
                if check.get("nominal_meets64_2001")
                and check.get("nominal_meets64_4001")]
            item["genuine_pairs_60"] = [
                list(check["pair"]) for check in level.get("checks", [])
                if check.get("nominal_meets60_2001")
                and check.get("nominal_meets60_4001")]
        item["recorded_fit_meets"] = {
            "nominal": {
                "meets60_2001": _meets(fit, "nominal", 2001, 60),
                "meets64_2001": _meets(fit, "nominal", 2001, 64),
                "meets60_4001": _meets(fit, "nominal", 4001, 60),
                "meets64_4001": _meets(fit, "nominal", 4001, 64),
            },
            "full": {
                "meets60_2001": _meets(fit, "full", 2001, 60),
                "meets64_2001": _meets(fit, "full", 2001, 64),
                "meets60_4001": _meets(fit, "full", 4001, 60),
                "meets64_4001": _meets(fit, "full", 4001, 64),
            },
        }
        # Fields the ledger already reads from the old cell record, kept
        # under their old names and meanings: ``corrected_count`` is the
        # record's nominal count, ``pair`` the pair that meets, and the
        # below-floor diagnostics describe the search.
        item["corrected_count"] = int(cell["nominal_count"])
        item["pair"] = pair
        old = legacy.get(key)
        for field in ("below_ratio", "best_below", "best_below_pair"):
            if old is not None and field in old:
                item[field] = old[field]
        cells.append(item)

    by_key = {(item["target"], item["width_name"]): item for item in cells}
    if len(cells) != 45 or len(by_key) != 45:
        raise SystemExit("expected 45 distinct Group A cells, got {}".format(
            len(by_key)))

    for target in TARGETS:
        for name in POW:
            if (target, name) not in by_key:
                raise SystemExit("missing cell {} {}".format(target, name))

    modes = {}
    for item in cells:
        modes[item["count_mode"]] = modes.get(item["count_mode"], 0) + 1
    if modes != {"nominal": 40, "genuine": 1, "full": 4}:
        raise SystemExit("unexpected mode counts: {}".format(modes))

    total_nominal = sum(item["nominal_count"] for item in cells)
    total_full = sum(item["full_count"] for item in cells)
    total_reported = sum(item["reported_count"] for item in cells)
    total_cheb = sum(item["chebyshev_count"] for item in cells)
    if (total_nominal, total_full, total_reported, total_cheb) != (
            361, 454, 392, 494):
        raise SystemExit("unexpected count sums: {}".format(
            (total_nominal, total_full, total_reported, total_cheb)))

    fewer = []
    equal = []
    more = []
    for item in cells:
        label = "{} {}".format(item["target"], item["width_name"])
        if item["reported_count"] < item["chebyshev_count"]:
            fewer.append(label)
        elif item["reported_count"] == item["chebyshev_count"]:
            equal.append(label)
        else:
            more.append(label)
    if (len(fewer), len(equal), len(more)) != (39, 3, 3):
        raise SystemExit("unexpected comparison: {}/{}/{}".format(
            len(fewer), len(equal), len(more)))
    if equal != ["rel:1e-13 m3", "abs:5e-14 m3", "rel:1e-14 m4"]:
        raise SystemExit("unexpected equal cells: {}".format(equal))
    expected_more = ["rel:1e-14 p4", "rel:1e-14 p5", "rel:1e-14 paper"]
    if more != expected_more:
        raise SystemExit("unexpected more cells: {}".format(more))

    charged = [(item["target"], item["width_name"], item["reported_count"],
                item["recorded_pair"])
               for item in cells if item["count_mode"] == "full"]
    expected_charged = [
        ("rel:1e-13", "p5", 21, [10, 4]),
        ("rel:1e-14", "p4", 21, [3, 10]),
        ("rel:1e-14", "p5", 29, [14, 3]),
        ("rel:1e-14", "paper", 21, [3, 10]),
    ]
    if charged != expected_charged:
        raise SystemExit("unexpected charged cells: {}".format(charged))

    genuine = [(item["target"], item["width_name"], item["reported_count"],
                item["reported_pair"])
               for item in cells if item["count_mode"] == "genuine"]
    if genuine != [("abs:5e-14", "p5", 14, [7, 6])]:
        raise SystemExit("unexpected genuine cells: {}".format(genuine))

    series = {}
    for target in TARGETS:
        series[target] = {
            "widths": list(POW),
            "reported": [by_key[(target, name)]["reported_count"]
                         for name in POW],
            "nominal": [by_key[(target, name)]["nominal_count"]
                        for name in POW],
            "chebyshev": [by_key[(target, name)]["chebyshev_count"]
                          for name in POW],
            "d": [by_key[(target, name)]["chebyshev_count"]
                  - by_key[(target, name)]["reported_count"]
                  for name in POW],
        }
    expected_series = {
        "rel:1e-13": {
            "reported": [4, 4, 4, 4, 5, 5, 6, 6, 7, 9, 11, 12, 21],
            "chebyshev": [5, 5, 5, 5, 5, 6, 7, 8, 9, 11, 14, 18, 25],
            "d": [1, 1, 1, 1, 0, 1, 1, 2, 2, 2, 3, 6, 4],
        },
        "abs:5e-14": {
            "reported": [4, 4, 4, 4, 5, 5, 6, 6, 7, 9, 10, 12, 14],
            "chebyshev": [5, 5, 5, 5, 5, 6, 7, 8, 9, 11, 14, 18, 24],
            "d": [1, 1, 1, 1, 0, 1, 1, 2, 2, 2, 4, 6, 10],
        },
        "rel:1e-14": {
            "reported": [4, 4, 4, 5, 5, 5, 6, 6, 9, 10, 11, 21, 29],
            "chebyshev": [5, 5, 5, 5, 6, 6, 7, 9, 10, 12, 15, 20, 26],
            "d": [1, 1, 1, 0, 1, 1, 1, 3, 1, 2, 4, -1, -3],
        },
    }
    for target in TARGETS:
        got = series[target]
        want = expected_series[target]
        if (got["reported"] != want["reported"]
                or got["chebyshev"] != want["chebyshev"]
                or got["d"] != want["d"]):
            raise SystemExit("series mismatch for {}: {}".format(target, got))

    full_vs_nominal = [
        {
            "target": item["target"],
            "width_name": item["width_name"],
            "nominal_count": item["nominal_count"],
            "full_count": item["full_count"],
            "chebyshev_count": item["chebyshev_count"],
            "count_mode": item["count_mode"],
        }
        for item in cells
        if not item["recorded_fit_meets"]["nominal"]["meets64_2001"]
    ]
    if len(full_vs_nominal) != 5:
        raise SystemExit("expected five full-vs-nominal cells, got {}".format(
            len(full_vs_nominal)))

    record = {
        "schema": "qcx-paper-c-reported-counts/v1",
        "generated_utc": audit.get("generated_utc"),
        "definition": DEFINITION,
        "source": {
            "record": SOURCE,
            "sha256": SOURCE_SHA256,
            "audit_generated_utc": audit.get("generated_utc"),
            "harness_tag": audit.get("method", {}).get("harness_tag"),
            "harness_commit": audit.get("method", {}).get("harness_commit"),
        },
        "cells": cells,
        "summary": {
            "cells": len(cells),
            "counts": {
                "nominal": total_nominal,
                "reported": total_reported,
                "full": total_full,
                "chebyshev": total_cheb,
            },
            "modes": modes,
            "comparison": {
                "fewer": len(fewer),
                "equal": len(equal),
                "more": len(more),
                "equal_cells": equal,
                "more_cells": more,
            },
            "full_vs_nominal": full_vs_nominal,
            "charged_cells": [
                {"target": target, "width_name": name, "reported_count":
                 count, "pair": pair}
                for target, name, count, pair in charged],
            "genuine_cells": [
                {"target": target, "width_name": name, "reported_count":
                 count, "pair": pair}
                for target, name, count, pair in genuine],
            "series": series,
            "worst_charge": {
                "cell": ["rel:1e-14", "p5"],
                "nominal": by_key[("rel:1e-14", "p5")]["nominal_count"],
                "reported": by_key[("rel:1e-14", "p5")]["reported_count"],
                "chebyshev": by_key[("rel:1e-14", "p5")][
                    "chebyshev_count"],
            },
        },
    }

    output_path = ROOT / OUTPUT
    output_path.write_text(
        json.dumps(record, indent=1, ensure_ascii=False) + "\n",
        encoding="utf-8")
    print("wrote {} ({} cells, modes {})".format(
        OUTPUT, len(cells),
        ", ".join("{}={}".format(key, modes[key]) for key in sorted(modes))))
    print("reported counts: {}; sums nominal/reported/full/chebyshev = "
          "{}/{}/{}/{}".format(
              OUTPUT, total_nominal, total_reported, total_full, total_cheb))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
