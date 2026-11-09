#!/usr/bin/env python3
"""Re-count the coefficient-count paper's Chebyshev cells from this deposit.

The Chebyshev fit itself is unchanged by the split-Clenshaw repair; only the
evaluator that scores it is.  This driver rebuilds every cell's grid, target
and degree ladder from the deposited source records and re-solves the
Chebyshev fit alone.  It never touches the rational solvers.

The grid for a cell is ``a = X0``, ``b = a + width``, the record's certifying
grid and the record's 60-digit Boys reference.  A cell is
``(width, target, degree ladder)``; the ladder is taken from each source
record's own settings (the ladder record admits degrees 0-2 as well).

Node convention
---------------
The source records store ``settings.nodes = 2001``: the number of points on
the certifying grid.  The sweep API is interval-based::

    make_nodes(g, a, b, n)  ->  n + 1 points

so this driver passes ``n = settings.nodes - 1 = 2000``.  The public record's
per-cell ``nodes`` field stores that API parameter (2000), not the point
count.  Both numbers are asserted at every grid construction.

Inputs
------
The v5.1.1 generator (the corrected split-Clenshaw evaluator) and the five
source records are deposited beside this file and in ``records/``.  The
recorded sweep harness is ``count_sweep_recorded.py`` at the campaign root;
the section S1 record inventory names it by digest.  The v5.0.0 generator is
the control: it re-solves the same cells with the pre-repair evaluator and
must reproduce every recorded pre-repair count.

Examples
--------
Corrected run, verified against the deposited record::

    python recount_cheb.py \\
        --out recount-run-v5.1.1.json --verify ../chebyshev_recount.json \\
        --report verify-new.json --metadata run-metadata.json

Pre-repair control, verified against the record's old counts::

    python recount_cheb.py --generator generators/gen_boys_coefficients_v5.0.0.py \\
        --out recount-control-v5.0.0.json --verify ../chebyshev_recount.json \\
        --expect old --report verify-old.json

Do not run with ``python -I``: the generator's mpmath comes from user-site.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import shutil
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parents[1]
REPO = CAMPAIGN.parent

DEFAULT_GENERATOR = HERE / "generators" / "gen_boys_coefficients_v5.1.1.py"
CONTROL_GENERATOR = HERE / "generators" / "gen_boys_coefficients_v5.0.0.py"
DEFAULT_COUNT_SWEEP = CAMPAIGN / "count_sweep_recorded.py"

SOURCE_RECORDS = (
    ("core", CAMPAIGN / "records" / "core" / "sweep_all.json"),
    ("abs", CAMPAIGN / "records" / "core" / "sweep_abs_5e-14.json"),
    ("rel13", CAMPAIGN / "records" / "core" / "sweep_rel_1e-13.json"),
    ("rel14", CAMPAIGN / "records" / "core" / "sweep_rel_1e-14.json"),
    ("ladder", CAMPAIGN / "records" / "ladder" / "sweep_all.json"),
)
GROUP_B = {"2.8567491415463806", "3.405194117991412", "4.091271749574425",
           "4.947466945561672", "3.7994698932176263"}
FALSIFIER_WIDTHS = ("2^1", "2^2", "2^3", "2^5")

# Provenance printed in the record; the generator digests are recomputed from
# the deposited files, the other three are the revision the record names.
BOYS_VERSION = "5.1.1"
BOYS_DESCRIBE = "v5.1.0-1-g3b187d5"
COEFFICIENT_HEADER_SHA256 = (
    "8ce5e17157a4be5f621d695f377094d212cb839732e278d47f23aa225b4d9204")

NODE_CONVENTION = {
    "source_record_field": "settings.nodes = 2001 is the number of grid points",
    "sweep_api": "make_nodes(g, a, b, n) takes n intervals and returns n + 1 points",
    "driver_parameter": "n = settings.nodes - 1 = 2000",
    "public_record_field": "cells[].nodes stores the API parameter 2000, not the point count",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()


SCRIPT_SHA256 = sha256(Path(__file__))


def rel(path: Path) -> str:
    """A repository-relative path when possible, else the absolute path."""
    p = Path(path).resolve()
    try:
        return p.relative_to(REPO).as_posix()
    except ValueError:
        return p.as_posix()


def load_count_sweep(path: Path):
    spec = importlib.util.spec_from_file_location("count_sweep_recorded", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["count_sweep_recorded"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_generator(cs, generator: Path):
    """Load a deposited generator file under the sweep's expected path.

    ``count_sweep.load_generator`` imports ``tools/gen_boys_coefficients.py``
    from a boys root.  The deposited generator is a single tagged file, so it
    is copied into a temporary root at that path and the digest is checked
    before and after.
    """
    expected = sha256(generator)
    tmp = tempfile.TemporaryDirectory(prefix="recount-boys-")
    root = Path(tmp.name)
    (root / "tools").mkdir()
    shutil.copyfile(generator, root / "tools" / "gen_boys_coefficients.py")
    g = cs.load_generator(root)
    if sha256(generator) != expected:
        raise RuntimeError(f"generator changed while loading: {generator}")
    return g, tmp


def collect_cells(source_records):
    """Unique (width, target, ladder, intervals) cells from the source rows."""
    cells = {}
    order = []
    for name, path in source_records:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
        for row in doc["rows"]:
            settings = row.get("settings", {})
            points = int(settings.get("nodes", 2001))
            intervals = points - 1
            mode = settings.get("cheb_degree_mode", "records")
            key = (row["width"], row["target"], mode, intervals)
            if key not in cells:
                cells[key] = {
                    "width": row["width"],
                    "target": row["target"],
                    "target_kind": row["target_kind"],
                    "target_value": row["target_value"],
                    "mode": mode,
                    "n_nodes": intervals,
                    "nodes_actual": points,
                    "rows": [],
                }
                order.append(key)
            cells[key]["rows"].append({
                "record": name,
                "width_name": row["width_name"],
                "source_stored": row["cheb"].get("stored"),
                "source_m": row["cheb"].get("m"),
            })
    return cells, order


def run(cs, generator, cells, order, progress=False):
    g, tmp = load_generator(cs, generator)
    try:
        a = g.X0
        grids = {}
        results = []
        started = time.perf_counter()
        for index, key in enumerate(order, 1):
            cell = cells[key]
            points = cell["nodes_actual"]
            intervals = cell["n_nodes"]
            if points != cs.N_NODES + 1:
                raise SystemExit(
                    f"{cell['width']} {cell['target']}: source record has "
                    f"{points} nodes; this deposit's convention is "
                    f"{cs.N_NODES + 1}")
            grid_key = (cell["width"], intervals)
            if grid_key not in grids:
                width = g.mpf(cell["width"])
                b = a + width
                _ts, xs = cs.make_nodes(g, a, b, intervals)
                if len(xs) != points:
                    raise SystemExit(
                        f"grid convention: {intervals} intervals gave "
                        f"{len(xs)} points, expected {points}")
                reference_started = time.perf_counter()
                ref = cs.reference(g, xs)
                grids[grid_key] = (width, b, xs, ref,
                                   time.perf_counter() - reference_started)
            width, b, xs, ref, ref_seconds = grids[grid_key]
            target = {"kind": cell["target_kind"],
                      "value": cell["target_value"],
                      "spec": cell["target"]}
            cheb = cs.cheb_count(g, a, b, xs, ref, target,
                                 degrees=cs.cheb_degrees(cell["mode"]))
            entry = dict(cell)
            entry["new"] = cheb
            entry["reference_seconds"] = round(ref_seconds, 3)
            results.append(entry)
            if progress:
                print(f"[{index}/{len(order)}] {cell['width']:>18} "
                      f"{cell['target']:>12} {cell['mode']:>7} "
                      f"new {cheb.get('stored')} ({cheb.get('seconds')} s)",
                      flush=True)
        return g, results, time.perf_counter() - started
    finally:
        tmp.cleanup()


def cell_key(cell):
    intervals = cell["n_nodes"] if "n_nodes" in cell else cell["nodes"]
    return (cell["width"], cell["target"], cell["mode"], intervals)


def annotate_with_record(results, record):
    """Add the deposited record's old count and changed flag to each cell."""
    expected = {cell_key(c): c for c in record["cells"]}
    for cell in results:
        rc = expected.get(cell_key(cell))
        if rc is None:
            continue
        cell["old_stored"] = rc["old_stored"]
        cell["changed"] = cell["new"].get("stored") != rc["old_stored"]
        for row in cell["rows"]:
            row["old_stored"] = rc["old_stored"]


def compare(record, results, expect):
    """Field-by-field verification of the run against the deposited record."""
    expected = {cell_key(c): c for c in record["cells"]}
    computed = {cell_key(c): c for c in results}
    report = {
        "record": {"cells": len(expected), "computed": len(computed)},
        "expect": expect,
        "missing": sorted(list(k) for k in set(expected) - set(computed)),
        "extra": sorted(list(k) for k in set(computed) - set(expected)),
        "fields": {},
        "changed": {},
        "summary": {},
        "provenance": {},
        "mismatches": [],
    }
    if expect == "new":
        values = {
            "new_stored": lambda c: c["new"].get("stored"),
            "new_m": lambda c: c["new"].get("m"),
            "delivered_rel": lambda c: c["new"].get("delivered_rel"),
            "delivered_abs": lambda c: c["new"].get("delivered_abs"),
            "at_x_rel": lambda c: c["new"].get("at_x_rel"),
            "at_x_abs": lambda c: c["new"].get("at_x_abs"),
        }
    else:
        values = {"old_stored": lambda c: c["new"].get("stored")}
    for name, get in values.items():
        checked = matched = 0
        for key, rc in expected.items():
            got_cell = computed.get(key)
            if got_cell is None:
                continue
            checked += 1
            got = get(got_cell)
            if expect == "old":
                want = rc["old_stored"]
            elif name == "new_stored":
                want = rc["new_stored"]
            elif name == "new_m":
                want = rc["new_m"]
            else:
                want = rc[name]
            if got == want:
                matched += 1
            else:
                report["mismatches"].append(
                    {"field": name, "cell": list(key), "record": want,
                     "computed": got})
        report["fields"][name] = {"checked": checked, "matched": matched}

    if expect == "new":
        changed = changed_matched = changed_mismatches = 0
        for key, rc in expected.items():
            got_cell = computed.get(key)
            if got_cell is None:
                continue
            if rc["changed"]:
                changed += 1
            is_changed = got_cell["new"].get("stored") != rc["old_stored"]
            if is_changed == rc["changed"]:
                changed_matched += 1
            else:
                changed_mismatches += 1
                report["mismatches"].append(
                    {"field": "changed", "cell": list(key),
                     "record": rc["changed"], "computed": is_changed})
        report["changed"] = {
            "record": sum(1 for c in record["cells"] if c["changed"]),
            "computed": sum(1 for key, rc in expected.items()
                            if computed.get(key) is not None
                            and computed[key]["new"].get("stored") != rc["old_stored"]),
            "checked": len(expected), "matched": changed_matched,
            "mismatches": changed_mismatches,
        }
        report["summary"] = compare_summary(record, computed)
    report["result"] = (
        "PASS" if not report["mismatches"] and not report["missing"]
        and not report["extra"]
        and all(v["checked"] == v["matched"]
                for v in report["fields"].values())
        and report["summary"].get("result", "PASS") in ("PASS", "SKIPPED")
        else "FAIL")
    return report


def compare_summary(record, computed):
    """Recompute the record's summary from the measured cells and source rows."""
    summary = record["summary"]
    if len(computed) != summary["unique_cells"]:
        return {"result": "SKIPPED", "checked": [],
                "reason": "partial run; the summary describes the full grid"}
    cells = list(computed.values())
    changed = sum(1 for c in cells if c["new"].get("stored") != c["old_stored"])
    deltas = sorted({c["new"].get("stored") - c["old_stored"] for c in cells})
    checks = {
        "unique_cells": (len(cells), summary["unique_cells"]),
        "changed_cells": (changed, summary["changed_cells"]),
        "unchanged_cells": (len(cells) - changed, summary["unchanged_cells"]),
        "deltas": (deltas, summary["deltas"]),
    }

    by_key = {cell_key(c): c for c in cells}

    def source_key(row):
        mode = row.get("settings", {}).get("cheb_degree_mode", "records")
        return (row["width"], row["target"], mode,
                int(row.get("settings", {}).get("nodes", 2001)) - 1)

    arm_parity = {}
    for name in ("rel13", "abs", "rel14"):
        rows = json.loads(dict(SOURCE_RECORDS)[name].read_text(
            encoding="utf-8"))["rows"]
        odd = even = 0
        for row in rows:
            if row["width_name"] in GROUP_B:
                continue
            stored = by_key[source_key(row)]["new"].get("stored")
            if stored % 2:
                odd += 1
            else:
                even += 1
        arm_parity[name] = {"odd": odd, "even": even}
    checks["group_a_parity"] = ({
        "odd": sum(a["odd"] for a in arm_parity.values()),
        "even": sum(a["even"] for a in arm_parity.values()),
        "of": sum(a["odd"] + a["even"] for a in arm_parity.values()),
        "by_arm": arm_parity,
    }, summary["group_a_parity"])

    ladder_odd = ladder_even = 0
    ladder_rows = json.loads(
        dict(SOURCE_RECORDS)["ladder"].read_text(encoding="utf-8"))["rows"]
    for row in ladder_rows:
        stored = by_key[source_key(row)]["new"].get("stored")
        if stored % 2:
            ladder_odd += 1
        else:
            ladder_even += 1
    checks["ladder_parity"] = (
        {"odd": ladder_odd, "even": ladder_even,
         "of": ladder_odd + ladder_even}, summary["ladder_parity"])

    group_b = {}
    for name in ("rel13", "abs"):
        rows = json.loads(dict(SOURCE_RECORDS)[name].read_text(
            encoding="utf-8"))["rows"]
        old = new = rational = 0
        for row in rows:
            if row["width_name"] not in GROUP_B:
                continue
            cell = by_key[source_key(row)]
            old += cell["old_stored"]
            new += cell["new"].get("stored")
            rational += row["rat"]["stored"]
        group_b[name] = {"chebyshev_old": old, "chebyshev_new": new,
                         "rational": rational}
    checks["group_b_sums"] = (group_b, summary["group_b_sums"])

    falsifier = {}
    for width_name in FALSIFIER_WIDTHS:
        values = {}
        for name in ("rel13", "rel14"):
            rows = json.loads(dict(SOURCE_RECORDS)[name].read_text(
                encoding="utf-8"))["rows"]
            row = next(r for r in rows if r["width_name"] == width_name)
            cell = by_key[source_key(row)]
            values[name + "_new"] = cell["new"].get("stored")
            values[name + "_old"] = cell["old_stored"]
        values["increment_old"] = values["rel14_old"] - values["rel13_old"]
        values["increment_new"] = values["rel14_new"] - values["rel13_new"]
        falsifier[width_name] = values
    checks["falsifier4_increments"] = (
        falsifier, summary["falsifier4_increments"])

    static_path = CAMPAIGN / "records" / "core" / "static_counts.json"
    if static_path.exists():
        static = json.loads(static_path.read_text(encoding="utf-8"))
        narrow = static["derived_counts"]["region_b_narrow"]
        shipped = {"chebyshev": narrow["chebyshev"]["stored"],
                   "rational": narrow["rational"]["stored"]}
        checks["shipped_totals"] = (shipped, summary["shipped_totals"])

    mismatches = []
    for name, (got, want) in checks.items():
        if got != want:
            mismatches.append({"field": name, "record": want, "computed": got})
    return {"checked": sorted(checks),
            "matched": len(checks) - len(mismatches),
            "mismatches": mismatches,
            "result": "PASS" if not mismatches else "FAIL"}


def compare_provenance(record, corrected_generator, control):
    """Check the record's paper-record and generator provenance against the deposit."""
    report = {"paper_records": {}, "boys": {}, "mismatches": []}
    digests = HERE / "paper-record-digests.tsv"
    if digests.exists():
        rows = digests.read_text(encoding="utf-8").splitlines()
        index = {name: i for i, name in enumerate(rows[0].split("\t"))}
        for line in rows[1:]:
            if not line.strip():
                continue
            parts = line.split("\t")
            name = parts[index["paper_records_key"]]
            original = parts[index["original_sha256"]]
            want = record.get("paper_records", {}).get(name, {}).get("sha256")
            ok = want == original
            report["paper_records"][name] = {
                "record": want, "original_sha256": original, "match": ok}
            if not ok:
                report["mismatches"].append(
                    {"field": f"paper_records.{name}", "record": want,
                     "computed": original})
    boys = record.get("boys", {})
    expected = {
        "generator_sha256": sha256(corrected_generator),
        "previous_generator_sha256": sha256(control),
        "coefficient_header_sha256": COEFFICIENT_HEADER_SHA256,
        "version": BOYS_VERSION,
        "describe": BOYS_DESCRIBE,
    }
    for name, value in expected.items():
        ok = boys.get(name) == value
        report["boys"][name] = {"record": boys.get(name), "deposit": value,
                                "match": ok}
        if not ok:
            report["mismatches"].append(
                {"field": f"boys.{name}", "record": boys.get(name),
                 "computed": value})
    report["result"] = "PASS" if not report["mismatches"] else "FAIL"
    return report


def metadata_document(args, generator, g, results, seconds, record_path,
                      source_records, count_sweep, report):
    import mpmath
    doc = {
        "what": "Per-run metadata for the Paper C Chebyshev recount driver",
        "driver": {
            "path": rel(Path(__file__)),
            "sha256": SCRIPT_SHA256,
        },
        "invocation": {
            "argv": [Path(sys.executable).name, *sys.argv],
            "cwd": rel(Path.cwd()),
        },
        "node_convention": NODE_CONVENTION,
        "generator": {
            "path": rel(generator),
            "sha256": sha256(generator),
        },
        "corrected_generator": {
            "path": rel(args.corrected_generator),
            "sha256": sha256(args.corrected_generator),
        },
        "count_sweep": {
            "path": rel(count_sweep),
            "sha256": sha256(count_sweep),
        },
        "source_records": [
            {"name": name, "path": rel(path), "sha256": sha256(path)}
            for name, path in source_records
        ],
        "output": {
            "path": rel(Path(args.out)),
            "sha256": sha256(Path(args.out)),
            "bytes": Path(args.out).stat().st_size,
        },
        "environment": {
            "python": platform.python_version(),
            "mpmath": mpmath.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
        },
        "seconds": round(seconds, 3),
        "cells": len(results),
        "x0": g.mp.nstr(g.X0, 17),
    }
    if record_path is not None:
        doc["record"] = {
            "path": rel(record_path),
            "sha256": sha256(record_path),
        }
    if args.report is not None and Path(args.report).exists():
        doc["report"] = {
            "path": rel(Path(args.report)),
            "sha256": sha256(Path(args.report)),
        }
    if report is not None:
        doc["verification"] = {
            "expect": report["expect"],
            "result": report["result"],
            "record_cells": report["record"]["cells"],
            "computed_cells": report["record"]["computed"],
            "fields": report["fields"],
            "changed": report.get("changed", {}),
            "summary": {
                "result": report.get("summary", {}).get("result", "PASS"),
                "checked": len(report.get("summary", {}).get("checked", [])),
            },
            "provenance": report.get("provenance", {}).get("result", "PASS"),
            "mismatches": (report["mismatches"]
                           + report.get("missing", [])
                           + report.get("extra", [])
                           + report.get("summary", {}).get("mismatches", [])
                           + report.get("provenance", {}).get("mismatches", [])),
        }
    return doc


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generator", type=Path, default=DEFAULT_GENERATOR)
    parser.add_argument("--corrected-generator", type=Path,
                        default=DEFAULT_GENERATOR,
                        help="the pre-repair record's corrected generator")
    parser.add_argument("--control-generator", type=Path,
                        default=CONTROL_GENERATOR,
                        help="the pre-repair generator, for the provenance check")
    parser.add_argument("--count-sweep", type=Path, default=DEFAULT_COUNT_SWEEP)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--verify", type=Path, default=None,
                        help="the deposited records/chebyshev_recount.json")
    parser.add_argument("--expect", choices=("new", "old"), default="new",
                        help="compare against the record's new or old counts")
    parser.add_argument("--report", type=Path, default=None,
                        help="write the verification report here")
    parser.add_argument("--metadata", type=Path, default=None,
                        help="write the per-run metadata here")
    parser.add_argument("--records", default="all",
                        help="comma list of record names, or all")
    parser.add_argument("--limit", type=int, default=0,
                        help="stop after this many unique cells (0 = all)")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args(argv)

    source_records = tuple(SOURCE_RECORDS)
    if args.records != "all":
        wanted = set(args.records.split(","))
        source_records = tuple(r for r in source_records if r[0] in wanted)
    missing = [str(p) for _n, p in source_records if not p.exists()]
    if missing:
        raise SystemExit(f"missing source records: {missing}")

    cs = load_count_sweep(args.count_sweep)
    cells, order = collect_cells(source_records)
    if args.limit:
        order = order[:args.limit]
    g, results, seconds = run(cs, args.generator, cells, order,
                              progress=args.progress)

    record = None
    if args.verify is not None:
        record = json.loads(args.verify.read_text(encoding="utf-8"))
        annotate_with_record(results, record)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    run_doc = {
        "what": "Paper C Chebyshev re-count from the deposited records",
        "driver": {"path": rel(Path(__file__)),
                   "sha256": SCRIPT_SHA256},
        "generator": {"path": rel(args.generator),
                      "sha256": sha256(args.generator)},
        "count_sweep": {"path": rel(args.count_sweep),
                        "sha256": sha256(args.count_sweep)},
        "x0": g.mp.nstr(g.X0, 17),
        "node_convention": NODE_CONVENTION,
        "source_records": [
            {"name": name, "path": rel(path), "sha256": sha256(path)}
            for name, path in source_records
        ],
        "cells": results,
        "seconds": round(seconds, 3),
    }
    args.out.write_text(json.dumps(run_doc, indent=1) + "\n",
                        encoding="utf-8")
    print(f"wrote {args.out}: {len(results)} cells in {seconds:.1f} s")

    report = None
    if record is not None:
        report = compare(record, results, args.expect)
        report["provenance"] = compare_provenance(
            record, args.corrected_generator, args.control_generator)
        if report["provenance"]["result"] != "PASS":
            report["result"] = "FAIL"
        if args.report is not None:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(json.dumps(report, indent=1) + "\n",
                                   encoding="utf-8")
        counts = report["fields"]
        print(f"verify {args.expect}: {report['result']}; "
              f"cells {report['record']['computed']}/{report['record']['cells']}; "
              + ", ".join(f"{k} {v['matched']}/{v['checked']}"
                          for k, v in counts.items()))

    if args.metadata is not None:
        doc = metadata_document(args, args.generator, g, results, seconds,
                                args.verify, source_records, args.count_sweep,
                                report)
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(doc, indent=1) + "\n",
                                 encoding="utf-8")
        print(f"wrote {args.metadata}")

    if report is not None and report["result"] != "PASS":
        for mismatch in report["mismatches"][:20]:
            print(f"  MISMATCH {mismatch}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
