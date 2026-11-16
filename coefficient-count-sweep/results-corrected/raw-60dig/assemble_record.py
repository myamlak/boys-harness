#!/usr/bin/env python3
"""Assemble the 45 raw cell records into corrected-counts-60dig.json."""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SWEEP = ROOT / "sweep"
OUT = ROOT / "corrected-counts-60dig.json"
CORRECTED = (ROOT / "coefficient-count-sweep" / "results-corrected"
             / "corrected-counts.json")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    corrected = json.loads(CORRECTED.read_text())
    n = len(corrected["cells"])
    cells = []
    missing = []
    for i in range(n):
        p = SWEEP / f"cell_g{i:02d}.json"
        if not p.exists():
            missing.append(i)
            continue
        cells.append(json.loads(p.read_text()))
    if missing:
        print(f"missing raw cells: {missing}", file=sys.stderr)
        return 1

    # ---- supplementary count-level sweeps for reproduction-failed cells.
    # The tagged control aborts before its count ladder when delivered_2001
    # does not reproduce, so `count60_level` is absent.  A sidecar file (see
    # analysis/supplement_count60.py) swept the level at the recorded nominal
    # count with the same scoring.  It is merged here and marked as such.
    supplements = []
    for c in cells:
        meas = c["measurements"]
        if meas.get("count60_level") is not None:
            continue
        sp = SWEEP / f"cell_g{c['index']:02d}_count60-supplement.json"
        if not sp.exists():
            continue
        sup = json.loads(sp.read_text())
        if sup["index"] != c["index"]:
            raise RuntimeError(f"supplement index {sup['index']} != "
                               f"{c['index']}")
        if sup["count60"] != c["nominal_count"]:
            raise RuntimeError(
                f"supplement count {sup['count60']} != nominal "
                f"{c['nominal_count']} for g{c['index']:02d}")
        lv = {k: sup[k] for k in (
            "count60", "swept_by_control", "extra_level_run", "rows_solved",
            "checks", "any_nominal_meets60_2001",
            "any_nominal_meets64_2001", "any_nominal_meets60_4001",
            "any_nominal_meets64_4001", "genuine_pairs_60",
            "genuine_pairs_64")}
        lv["supplement"] = True
        lv["supplement_reason"] = sup["reason"]
        lv["supplement_script_sha256"] = sup["meta"][
            "supplement_script_sha256"]
        lv["supplement_generated_utc"] = sup["meta"]["generated_utc"]
        meas["count60_level"] = lv
        supplements.append({
            "index": c["index"],
            "cell": f"{c['width_name']} {c['target']}",
            "count60": sup["count60"],
            "file": sp.name,
            "sha256": sha256(sp),
            "reason": sup["reason"],
        })

    # ---- summary counts
    failed = []
    repro_failed = []
    count60_none = []
    for c in cells:
        st = c["measurements"].get("control_status")
        name = f"{c['width_name']} {c['target']}"
        if st != "ok" or c["measurements"].get("recorded_fit") is None:
            failed.append({
                "index": c["index"], "width_name": c["width_name"],
                "target": c["target"], "status": st,
                "control": c["control"].get("status")
                if isinstance(c["control"], dict) else None,
                "recorded_fit_present":
                    c["measurements"].get("recorded_fit") is not None,
            })
        if st == "reproduction-failed":
            rep = c["measurements"].get("reproduction") or {}
            repro_failed.append({
                "index": c["index"], "cell": name,
                "recorded_delivered_2001": rep.get(
                    "record_delivered_2001"),
                "rerun_e64": rep.get("e64"),
                "settings": rep.get("settings"),
            })
        if c["measurements"].get("count60") is None:
            count60_none.append({
                "index": c["index"], "cell": name,
                "status": st,
                "control_error": c["control"].get("error")
                if isinstance(c["control"], dict) else None,
                "supplemented":
                    c["measurements"].get("count60_level") is not None,
            })

    def cell_flag(c, block, key):
        rf = c["measurements"].get("recorded_fit")
        if not rf:
            return None
        return rf[block][key]

    disagree = {
        "full_vs_nominal_60_2001": [],
        "full_vs_nominal_64_2001": [],
        "full_vs_nominal_60_4001": [],
        "full_vs_nominal_64_4001": [],
    }
    disagree_src = {
        "full_vs_nominal_60_2001": "disagree_60_2001",
        "full_vs_nominal_64_2001": "disagree_64_2001",
        "full_vs_nominal_60_4001": "disagree_60_4001",
        "full_vs_nominal_64_4001": "disagree_64_4001",
    }
    for c in cells:
        rf = c["measurements"].get("recorded_fit")
        if not rf:
            continue
        for key, src in disagree_src.items():
            if rf.get(src):
                disagree[key].append(f"{c['width_name']} {c['target']}")

    genuine = {}
    for basis in ("any_nominal_meets60_2001", "any_nominal_meets64_2001",
                  "any_nominal_meets60_4001", "any_nominal_meets64_4001"):
        hits = []
        for c in cells:
            lv = c["measurements"].get("count60_level")
            if lv and lv.get(basis):
                hits.append({
                    "width_name": c["width_name"], "target": c["target"],
                    "pairs": lv["genuine_pairs_60"] if "60" in basis
                    else lv["genuine_pairs_64"],
                    "count60": lv["count60"],
                })
        genuine[basis] = hits

    # recorded-fit decisive test (the Q-429 question), both bases
    decisive = []
    for c in cells:
        rf = c["measurements"].get("recorded_fit")
        if not rf:
            continue
        decisive.append({
            "width_name": c["width_name"],
            "target": c["target"],
            "pair": rf["pair"],
            "nominal_count": rf["nominal_count"],
            "full_count": rf["full_count"],
            "stored_witness_full_count":
                c["measurements"]["stored_witness"]["full_evaluated_count"],
            "chebyshev_count": c["measurements"]["chebyshev"]["count_new"],
            "full_meets60_2001": rf["full"]["grid_2001"]["meets60"],
            "nominal_meets60_2001": rf["nominal"]["grid_2001"]["meets60"],
            "full_meets64_2001": rf["full"]["grid_2001"]["meets64"],
            "nominal_meets64_2001": rf["nominal"]["grid_2001"]["meets64"],
            "full_meets60_4001": rf["full"]["grid_4001"]["meets60"],
            "nominal_meets60_4001": rf["nominal"]["grid_4001"]["meets60"],
            "full_meets64_4001": rf["full"]["grid_4001"]["meets64"],
            "nominal_meets64_4001": rf["nominal"]["grid_4001"]["meets64"],
            "full_e60_2001": rf["full"]["grid_2001"]["e60"],
            "nominal_e60_2001": rf["nominal"]["grid_2001"]["e60"],
            "full_guard60_2001": rf["full"]["grid_2001"]["guard60"],
            "nominal_guard60_2001": rf["nominal"]["grid_2001"]["guard60"],
        })

    # ---- count ordering: nominal, full witness, Chebyshev
    order = {
        "cells_nominal_lt_full": [],
        "cells_nominal_lt_cheb": [],
        "cells_cheb_lt_full": [],
        "reversals_full_gt_cheb": [],
    }
    order_stored = {
        "cells_cheb_lt_stored_witness": [],
        "reversals_stored_witness_gt_cheb": [],
    }
    sums = {"nominal": 0, "full": 0, "chebyshev": 0}
    for c in decisive:
        nm, fu, ch = (c["nominal_count"], c["full_count"],
                      c["chebyshev_count"])
        sw = c["stored_witness_full_count"]
        sums["nominal"] += nm
        sums["full"] += fu
        sums["chebyshev"] += ch
        name = f"{c['width_name']} {c['target']}"
        if nm < fu:
            order["cells_nominal_lt_full"].append(name)
        if nm < ch:
            order["cells_nominal_lt_cheb"].append(name)
        if ch < fu:
            order["cells_cheb_lt_full"].append(name)
            order["reversals_full_gt_cheb"].append(
                {"cell": name, "nominal": nm, "full": fu, "chebyshev": ch})
        if ch < sw:
            order_stored["cells_cheb_lt_stored_witness"].append(name)
            order_stored["reversals_stored_witness_gt_cheb"].append(
                {"cell": name, "nominal": nm, "stored_witness": sw,
                 "chebyshev": ch})

    # ---- direction of every disagreement, both bases
    directions = {}
    for basis, fkey, nkey in (
            ("60_2001", "full_meets60_2001", "nominal_meets60_2001"),
            ("64_2001", "full_meets64_2001", "nominal_meets64_2001"),
            ("60_4001", "full_meets60_4001", "nominal_meets60_4001"),
            ("64_4001", "full_meets64_4001", "nominal_meets64_4001")):
        directions[f"full_meets_nominal_fails_{basis}"] = [
            f"{d['width_name']} {d['target']}" for d in decisive
            if d[fkey] and not d[nkey]]
        directions[f"full_fails_nominal_meets_{basis}"] = [
            f"{d['width_name']} {d['target']}" for d in decisive
            if d[nkey] and not d[fkey]]

    # ---- self-check: the recorded fit's full 60-digit score must equal the
    # control's own confirm-row score (same computation, two paths)
    inconsistent = []
    for c in cells:
        rf = c["measurements"].get("recorded_fit")
        if not rf:
            continue
        if (rf["full"]["grid_2001"]["e60"]
                != rf["control_row"]["e60"]):
            inconsistent.append(
                f"{c['width_name']} {c['target']}: "
                f"{rf['full']['grid_2001']['e60']} != "
                f"{rf['control_row']['e60']}")

    # ---- genuine nominal fit at the recorded pair (the Q-429 decisive test)
    recorded_pair_genuine = {
        "nominal_meets60_2001": [
            f"{d['width_name']} {d['target']} [{d['pair'][0]}/{d['pair'][1]}]"
            for d in decisive if d["nominal_meets60_2001"]],
        "nominal_meets64_2001": [
            f"{d['width_name']} {d['target']} [{d['pair'][0]}/{d['pair'][1]}]"
            for d in decisive if d["nominal_meets64_2001"]],
        "nominal_meets60_4001": [
            f"{d['width_name']} {d['target']} [{d['pair'][0]}/{d['pair'][1]}]"
            for d in decisive if d["nominal_meets60_4001"]],
        "nominal_meets64_4001": [
            f"{d['width_name']} {d['target']} [{d['pair'][0]}/{d['pair'][1]}]"
            for d in decisive if d["nominal_meets64_4001"]],
    }

    doc = {
        "schema": "corrected-counts-60dig-1",
        "generated_utc": datetime.now(timezone.utc)
        .strftime("%Y-%m-%dT%H:%M:%SZ"),
        "method": {
            "purpose": (
                "Q-429 remedy: every Group A cell re-solved with the v1.0.39 "
                "60-digit control; the stored witness and its nominal "
                "truncation are scored separately against the cell target"),
            "harness_tag": "v1.0.39",
            "harness_commit": "309cf873d04e3466bf072b2b19bb8ee9718d4546",
            "harness_tag_commit": "85856bbe49cee5a915c86ea0fa07bcd7f1e4a399",
            "control_script": (
                "coefficient-count-sweep/control-60dig/"
                "control_60dig_brasil.py at v1.0.39"),
            "control_script_sha256": sha256(
                ROOT / "coefficient-count-sweep" / "control-60dig"
                / "control_60dig_brasil.py"),
            "driver_script": "sweep_cells.py",
            "driver_script_sha256": (
                sorted({c["meta"]["script_sha256"] for c in cells})[0]
                if len({c["meta"]["script_sha256"] for c in cells}) == 1
                else sorted({c["meta"]["script_sha256"] for c in cells})),
            "corrected_record_sha256": sha256(CORRECTED),
            "chebyshev_record_sha256": sha256(
                ROOT / "coefficient-count-sweep" / "records"
                / "chebyshev_recount.json"),
            "root_isolation_sha256": sha256(
                ROOT / "coefficient-count-sweep" / "strong-optional"
                / "root-isolation" / "records" / "root_isolation.json"),
            "versions": cells[0]["meta"]["versions"],
            "convert_dps": 70,
            "evaluation_dps": 60,
            "rounding": "binary64 after the 70-digit conversion",
            "horner": "both routes, fused (math.fma) and separate; "
                      "the worse error is reported",
            "grids": {"certifying_nodes": 2001, "check_nodes": 4001,
                      "t_distribution": "uniform on [-1, 1]"},
            "reference": "F_0 at 70 decimal digits (mpmath)",
            "denominator_guard": "min|Q|/max|Q| >= 1e-3",
            "target_kinds": ["rel", "abs"],
            "score_note": (
                "the tagged control scores the relative error; for the "
                "abs:5e-14 cells the wrapper scores the absolute error, "
                "matching the corrected record's own acceptance"),
            "control_walk": (
                "reproduce the record pair; walk downward with a full sweep "
                "of every admissible pair at each count (default route, "
                "tightened only if the default route found no meeting fit), "
                "stopping at the first count where neither route meets; "
                "accepted = meeting fit with smallest e60; 4001 check on "
                "the accepted fit"),
            "supplement_note": (
                "cells whose record delivered_2001 did not reproduce have no "
                "control count60 ladder; for those the level at the recorded "
                "nominal count was swept by analysis/supplement_count60.py "
                "with the same scoring and routes, recorded in "
                "summary.supplementary_count60_sweeps and marked "
                "count60_level.supplement=true; the raw cell files are "
                "unchanged"),
        },
        "summary": {
            "cells_in_corrected_record": n,
            "cells_measured": sum(
                1 for c in cells
                if c["measurements"].get("recorded_fit") is not None),
            "cells_swept": len(cells) - len(failed),
            "cells_failed": failed,
            "cells_reproduction_failed": repro_failed,
            "cells_without_count60": count60_none,
            "supplementary_count60_sweeps": supplements,
            "disagreements_recorded_fit": {
                k: {"count": len(v), "cells": v} for k, v in disagree.items()},
            "disagreement_directions": {
                k: {"count": len(v), "cells": v}
                for k, v in directions.items()},
            "genuine_nominal_at_recorded_pair": {
                k: {"count": len(v), "cells": v}
                for k, v in recorded_pair_genuine.items()},
            "genuine_nominal_at_count60": {
                k: {"count": len(v), "cells": v} for k, v in genuine.items()},
            "recorded_fit_score_consistency": {
                "mismatches": inconsistent,
                "count": len(inconsistent)},
            "count_ordering": {
                k: {"count": len(v), "cells": v} for k, v in order.items()},
            "count_ordering_stored_witness": {
                k: {"count": len(v), "cells": v}
                for k, v in order_stored.items()},
            "count_sums": sums,
            "recorded_fit_decisive": decisive,
        },
        "cells": cells,
    }
    OUT.write_text(json.dumps(doc, indent=1))
    print(f"wrote {OUT}")
    print(f"sha256 {sha256(OUT)}")
    print(json.dumps({
        "cells_swept": len(cells) - len(failed),
        "cells_measured": sum(
            1 for c in cells
            if c["measurements"].get("recorded_fit") is not None),
        "cells_failed": len(failed),
        "cells_without_count60": len(count60_none),
        "supplementary_count60_sweeps": len(supplements),
        "disagree_60_2001": len(disagree["full_vs_nominal_60_2001"]),
        "disagree_64_2001": len(disagree["full_vs_nominal_64_2001"]),
        "recorded_pair_genuine_60_2001":
            len(recorded_pair_genuine["nominal_meets60_2001"]),
        "recorded_pair_genuine_64_2001":
            len(recorded_pair_genuine["nominal_meets64_2001"]),
        "genuine_60_2001":
            len(genuine["any_nominal_meets60_2001"]),
        "genuine_64_2001":
            len(genuine["any_nominal_meets64_2001"]),
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
