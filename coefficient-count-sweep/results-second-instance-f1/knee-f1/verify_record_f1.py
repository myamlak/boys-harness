#!/usr/bin/env python3
"""Verify the F1 second-instance record before analysing it.

Checks, against the record itself and the paper's own records:

  1. coverage: 90 unique (edge, target, width) cells, 45 per edge, 30 per
     target, all fifteen Group A widths at both edges and all three targets;
  2. the recorded settings are the recorded ones, uniform in every row;
  3. the a2 edge is exactly twice x0: b - a = 2 * (b_x0 - a_x0), a = 2 x0;
  4. every cell's delivered error meets its own target;
  5. every row is internally consistent (stored = m + 1 [+ k], parity);
  6. the file digests in the meta equal the files on disk;
  7. the F1 reference validation and the fit-path check records;
  8. the lane's own three-cell F0 control reproduces the paper's recount
     (Chebyshev) and core (rational) records on six fields;
  9. the two knee rows are the only rows with below_meets true.

Exit 0 only when every check passes.  Writes _knee-f1/verify_record.json.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
RECORD = PAPER / "results-second-instance-f1" / "sweep_all.json"
CONTROL = PAPER / "results-second-instance-f1-control" / "sweep_all.json"
OUT = PAPER / "_knee-f1" / "verify_record.json"
RECOUNT = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
               r"\records\chebyshev_recount.json")
CORE = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
            r"\records\core")
TARGETS = ("rel:1e-13", "abs:5e-14", "rel:1e-14")
WIDTHS = [f"2^{j}" for j in range(-7, 6)] + ["ship", "paper"]
FIELDS = ("stored", "m", "delivered_rel", "delivered_abs",
          "at_x_rel", "at_x_abs")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def field(entry: dict, name: str):
    if name == "stored":
        return entry.get("stored")
    if name == "m":
        return entry.get("m")
    return entry.get(name)


def main() -> int:
    doc = json.loads(RECORD.read_text(encoding="utf-8"))
    meta = doc["meta"]
    rows = doc["rows"]
    fails = []
    out = {"record": str(RECORD), "checks": {}, "fails": fails}

    def check(name, ok, detail):
        out["checks"][name] = {"ok": bool(ok), "detail": detail}
        if not ok:
            fails.append(name)

    # 1. coverage
    keys = [(r["edge"], r["target"], r["width_name"]) for r in rows]
    check("rows_90", len(rows) == 90, len(rows))
    check("unique_cells", len(set(keys)) == len(keys), len(set(keys)))
    for edge, n in (("x0", 45), ("a2", 45)):
        check(f"edge_{edge}_45", sum(r["edge"] == edge for r in rows) == n,
              sum(r["edge"] == edge for r in rows))
    for target, n in zip(TARGETS, (30, 30, 30)):
        check(f"target_{target}_30",
              sum(r["target"] == target for r in rows) == n,
              sum(r["target"] == target for r in rows))
    for edge in ("x0", "a2"):
        for target in TARGETS:
            names = {r["width_name"] for r in rows
                     if r["edge"] == edge and r["target"] == target}
            check(f"complete_{edge}_{target}", names == set(WIDTHS),
                  sorted(set(WIDTHS) - names))

    # 2. settings uniform and equal to the meta
    want = {"nodes": meta["grid"]["nodes"],
            "coarse_stride": meta["grid"]["coarse_stride"],
            "scan_outer": meta["solver"]["scan_outer"],
            "refine_outer": meta["solver"]["refine_outer"],
            "remez_iters": meta["solver"]["remez_iters"],
            "rat_k_max": meta["solver"]["rat_k_max"],
            "cheb_degree_mode": meta["solver"]["cheb_degree_mode"]}
    bad = []
    for r in rows:
        s = r["settings"]
        for k, v in want.items():
            if s.get(k) != v:
                bad.append([r["edge"], r["target"], r["width_name"], k,
                            s.get(k), v])
    check("settings_uniform", not bad, bad[:5])
    check("grid_2001_5", meta["grid"]["nodes"] == 2001
          and meta["grid"]["coarse_stride"] == 5, meta["grid"])

    # 3. a2 is twice x0
    by = {(r["edge"], r["target"], r["width_name"]): r for r in rows}
    dist = []
    for r in rows:
        if r["edge"] != "a2":
            continue
        r0 = by[("x0", r["target"], r["width_name"])]
        a_ratio = float(r["a"]) / float(r0["a"])
        w_ratio = ((float(r["b"]) - float(r["a"]))
                   / (float(r0["b"]) - float(r0["a"])))
        dist.append([r["width_name"], a_ratio, w_ratio])
    check("a2_double", all(abs(d[1] - 2) < 1e-15 for d in dist)
          and all(abs(d[2] - 2) < 1e-14 for d in dist), dist[:3])

    # 4. every cell meets its target
    misses = []
    for r in rows:
        eps = float(r["target_value"])
        for fam in ("cheb", "rat"):
            entry = r[fam]
            if entry.get("stored") is None:
                misses.append([fam, r["edge"], r["target"], r["width_name"],
                               "no stored count"])
                continue
            if r["target_kind"] == "rel":
                got = float(entry["delivered_rel"])
            else:
                got = float(entry["delivered_abs"])
            if not got <= eps:
                misses.append([fam, r["edge"], r["target"],
                               r["width_name"], got, eps])
    check("all_cells_meet", not misses, misses[:5])

    # 5. internal consistency
    bad = []
    for r in rows:
        c, t = r["cheb"], r["rat"]
        if c["stored"] != c["m"] + 1:
            bad.append(["cheb", r["edge"], r["target"], r["width_name"]])
        if t["stored"] != t["m"] + 1 + t["k"]:
            bad.append(["rat", r["edge"], r["target"], r["width_name"]])
        if r.get("parity") is not None:
            # The parity field is a per-integer label record, not a type name.
            par = r["parity"]
            for key, value in (("cheb_degree", c["m"]),
                               ("rat_m", t["m"]), ("rat_k", t["k"])):
                want = "even" if value % 2 == 0 else "odd"
                if par.get(key + "_parity") != want:
                    bad.append(["parity", r["edge"], r["target"],
                                r["width_name"], key, par.get(key + "_parity"),
                                want])
    check("rows_consistent", not bad, bad[:5])

    # 6. digests
    for name, entry in (("header", meta["header"]),
                        ("generator", meta["generator"]),
                        ("harness", meta["harness"]),
                        ("f1_reference", meta["f1_reference"]),
                        ("driver", meta["driver"])):
        path = Path(entry["path"])
        got = sha(path) if path.exists() else None
        check(f"digest_{name}", got == entry["sha256"],
              {"path": str(path), "record": entry["sha256"], "disk": got})

    # 7. reference validation and fit-path check
    val = json.loads((PAPER / "results-second-instance-f1"
                      / "f1_validation.json").read_text(encoding="utf-8"))
    fit = json.loads((PAPER / "results-second-instance-f1"
                      / "fitpath_check.json").read_text(encoding="utf-8"))
    check("f1_reference_validation",
          "57" in str(val.get("verdict")) or "at least" in str(val.get("verdict")),
          {"verdict": val.get("verdict"), "worst_digits": val.get("worst_digits")})
    fit_rows = [c for bands in fit.values() if isinstance(bands, list)
                for c in bands]
    fit_max = max((c["max_abs_diff"] for c in fit_rows), default=None)
    fit_rel = max((c["max_rel_to_max_coeff"] for c in fit_rows), default=None)
    check("fitpath_check", fit_max is not None and fit_max < 1e-30,
          {"max_abs_diff": fit_max, "max_rel_to_max_coeff": fit_rel,
           "n_bands": len(fit_rows)})

    # 8. the lane's three-cell F0 control
    control = json.loads(CONTROL.read_text(encoding="utf-8"))
    recount = json.loads(RECOUNT.read_text(encoding="utf-8"))
    rec_cells = {}
    for cell in recount["cells"]:
        if (cell.get("mode") == "records" and cell.get("nodes") == 2000
                and len(cell.get("width_names") or []) == 1):
            rec_cells[(cell["target"], cell["width_names"][0])] = cell
    core = {t: json.loads((CORE / f"sweep_{t.replace(':', '_')}.json")
                          .read_text(encoding="utf-8"))
            for t in TARGETS}
    core_rows = {(t, row["target"], row.get("width_name"), row["a"]): row
                 for t, doc_ in core.items()
                 for row in doc_["rows"]}
    diffs = []
    for row in control["rows"]:
        target, name = row["target"], row["width_name"]
        want = rec_cells.get((target, name))
        for f in FIELDS:
            if want is not None and field(row["cheb"], f) != \
                    (want.get("new_stored") if f == "stored" else
                     want.get("new_m") if f == "m" else want.get(f)):
                diffs.append(["cheb", target, name, f])
        key = [t for t in TARGETS if any(
            k[0] == t and k[1] == target and k[2] == name
            and k[3].startswith(row["a"][:12]) for k in core_rows)]
        core_row = next((core_rows[k] for k in core_rows
                         if k[1] == target and k[2] == name
                         and k[3].startswith(row["a"][:12])), None)
        for f in FIELDS:
            if core_row is None or field(row["rat"], f) != \
                    field(core_row["rat"], f):
                diffs.append(["rat", target, name, f])
    check("lane_control_3", len(control["rows"]) == 3 and not diffs, diffs[:6])

    # 9. below_meets rows
    below = [[r["edge"], r["target"], r["width_name"], r["rat"]["stored"],
              r["rat"]["m"], r["rat"]["k"], r["rat"]["below_meets"]]
             for r in rows if r["rat"].get("below_meets")]
    check("below_rows_two", len(below) == 2, below)
    out["below_rows"] = below
    out["meta"] = {"generated_utc": meta.get("generated_utc"),
                   "kernel_version": meta.get("kernel_version"),
                   "python": meta.get("python"),
                   "platform": meta.get("platform")}

    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
    for name, c in out["checks"].items():
        print(f"{'PASS' if c['ok'] else 'FAIL'} {name}: {c['detail']}")
    print(f"\nchecks {len(out['checks'])}; fails {len(fails)}")
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
