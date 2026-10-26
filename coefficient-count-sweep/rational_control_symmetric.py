#!/usr/bin/env python3
"""The symmetric control: coefficient rounding and the Chebyshev endpoint cells.

WHAT THIS ADDS.  rational_unrounded_control.py re-solved the seven decisive
rational cells at 60 digits and scored each accepted fit twice: at 60 digits
(e60) and through the full binary64 route (e64).  That separates the solver's
precision from the stored arithmetic as a whole, but it cannot say how much of
the binary64 loss is the coefficient storage and how much is the evaluation.
This script adds three measurements:

  rational    every accepted 60-digit fit of the seven cells, re-solved from
              its own pair and solver, scored a third way: the coefficients
              rounded to binary64 and the series evaluated at 60 digits on the
              exact t grid (e64r).  The chain e60 -> e64r -> e64 attributes the
              loss: e64r - e60 is the coefficient storage, e64 - e64r is the
              evaluation arithmetic and the x -> t map.

  chebyshev   the same seven (width, target) cells on the Chebyshev side.  A
              Chebyshev fit is a coefficient vector with no solver search, so
              the control is direct: the record's corrected count, the count
              the same interpolation reaches when it is evaluated at 60 digits,
              and the full binary64 score of the 60-digit fit.  The record's
              own degree is re-solved at 60 digits beside it, and its binary64
              score must reproduce the recount record digit for digit.

  flip        the one accepted rational pair whose delivered margin moves
              between the stability pass and its follow-up: [7/8] Lawson at
              2^5, rel 1e-13.  The pass's two cached fits are re-scored on
              every grid, and the pair is solved again serially, one process,
              at the default and the doubled settings.  Each solve runs on two
              grids: the one the pass rebuilt from the record's 17-digit
              endpoint strings, and the one the sweep and the re-certification
              build from the generator's kX0.  The two print the same 17 digits
              and differ in the 32nd decimal place; whether that difference
              moves the delivered margin is the question the record leaves
              open.

  table       the six wide-end rows, assembled from the committed records: no
              fit is solved for it.

Every number is written to results-ladder/rational_control_symmetric.json.
Parts run independently (--part), and a completed part is reused on a rerun
unless --force is given.

Run WITHOUT -I: isolated mode hides the user-site mpmath the generator needs.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import re
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

PAPER = Path(__file__).resolve().parent
OUT_DEFAULT = PAPER / "results-ladder" / "rational_control_symmetric.json"
SWEEP = PAPER / "results-ladder" / "sweep_all.json"
RECERT = PAPER / "results-ladder" / "recertification.json"
CONTROL = PAPER / "results-ladder" / "rational_unrounded_control.json"
CORE_1E13 = PAPER / "results" / "sweep_rel_1e-13.json"
PASS_LATEST = PAPER / "results-stability" / "pass_2026-10-22.json"
CHEB_RECOUNT = Path(
    r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
    r"\records\chebyshev_recount.json")
RAT_BOYS = Path(r"C:\Users\myaml\dev\bladawiec\external\boys")
CHEB_BOYS = Path(r"C:\Users\myaml\dev\_mainonly")

EVAL_DPS = 60
K_MAX = 12
SCHEMA = "rational-control-symmetric-1"

# The seven cells the existing control uses, in its order.
CELLS = [
    ("paper", "rel:1e-08"),
    ("ship", "rel:1e-08"),
    ("2^5", "rel:1e-08"),
    ("paper", "rel:1e-14"),
    ("2^5", "rel:1e-14"),
    ("ship", "rel:1e-15"),
    ("2^5", "rel:1e-15"),
]

# The six wide-end rows of the paper's rational ladder.
WIDE_END = [
    ("2^5", "rel:1e-13"),
    ("2^5", "rel:1e-14"),
    ("2^5", "rel:1e-15"),
    ("ship", "rel:1e-15"),
    ("paper", "rel:1e-14"),
    ("paper", "rel:1e-15"),
]

FLIP_ROW = "rel:1e-13|2^5"
FLIP_M, FLIP_K = 7, 8

# The five scoring grids of the flip: two constructions at the default node
# count, two at the doubled count, and the pass's own interleaved doubled grid
# rebuilt from the record-string construction.
GRID_KEYS = [
    ("record_strings", "default"),
    ("record_strings", "doubled"),
    ("x0_based", "default"),
    ("x0_based", "doubled"),
    ("pass_interleaved", "doubled"),
]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def nstr(g, value, digits=17):
    return g.mp.nstr(value, digits)


def ratio(g, a, b):
    if b == 0:
        return None
    return nstr(g, a / b)


def first_differing_digit(g, x, y, hi=60):
    """First decimal digit at which the two printed values disagree."""
    for digits in range(1, hi + 1):
        if nstr(g, x, digits) != nstr(g, y, digits):
            return digits
    return None


def boys_version(root: Path):
    """The CMake project version of a boys checkout, read from its top line."""
    try:
        text = (root / "CMakeLists.txt").read_text(errors="replace")
    except OSError:
        return None
    m = re.search(r"^project\(\s*boys\s+VERSION\s+([0-9.]+)", text, re.M)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# the Chebyshev 60-digit evaluation


def cheb_value_60(g, coeffs, t):
    """Clenshaw's recurrence for sum c_k T_k(t), at the working precision.

    The kernel's split form is an algebraically identical rearrangement of this
    sum; evaluating the sum itself is the control's point, so the control does
    not reuse the split recursion it is checking.
    """
    b1 = g.mpf(0)
    b2 = g.mpf(0)
    for c in reversed(coeffs[1:]):
        b0 = 2 * t * b1 - b2 + c
        b2 = b1
        b1 = b0
    return t * b1 - b2 + coeffs[0]


def cheb_direct_60(g, coeffs, t):
    """The same sum written out through mpmath's own T_k."""
    total = g.mpf(0)
    for k, c in enumerate(coeffs):
        total += c * g.mp.chebyt(k, t)
    return total


def cheb_clenshaw_check(g, coeffs, ts, probes=7):
    """The Clenshaw recurrence against the direct sum at a few nodes."""
    mp = g.mp
    worst = g.mpf(0)
    where = None
    step = max(1, (len(ts) - 1) // (probes - 1))
    with mp.workdps(EVAL_DPS):
        for t in ts[::step][:probes]:
            a = cheb_value_60(g, coeffs, t)
            b = cheb_direct_60(g, coeffs, t)
            den = abs(b) if b != 0 else g.mpf(1)
            rel = abs(a - b) / den
            if rel > worst:
                worst, where = rel, t
    return worst, where


def cheb_errors_60(g, coeffs, ts, ref):
    """Worst relative and absolute error of a Chebyshev series at 60 digits."""
    mp = g.mp
    worst_rel = g.mpf(0)
    worst_abs = g.mpf(0)
    at_rel = None
    at_abs = None
    with mp.workdps(EVAL_DPS):
        for t, r in zip(ts, ref):
            v = cheb_value_60(g, coeffs, t)
            e = abs(v - r)
            if e > worst_abs:
                worst_abs, at_abs = e, t
            rel = e / r
            if rel > worst_rel:
                worst_rel, at_rel = rel, t
    return worst_rel, worst_abs, at_rel, at_abs


def cheb_coeffs_60(g, a, b, deg):
    """The interpolant's coefficients at 60 digits, not rounded."""
    with g.mp.workdps(g.RAT_DPS):
        return g.cheb_coeffs(0, a, b, deg)


def rounded(g, values):
    """The same values rounded to binary64 and carried exactly at 60 digits."""
    return [g.mpf(float(v)) for v in values]


def score_60(g, cs, target, wr, wa, at_r, at_a):
    return {
        "rel": nstr(g, wr),
        "abs": nstr(g, wa),
        "at_t": nstr(g, at_r) if at_r is not None else None,
        "at_t_abs": nstr(g, at_a) if at_a is not None else None,
        "meets": bool(cs.meets(target, wr, wa)),
    }


def score_64(g, cs, target, err):
    if err is None:
        return {"rel": None, "abs": None, "meets": None}
    wr, wa, at_rel, at_abs = err
    return {
        "rel": nstr(g, wr),
        "abs": nstr(g, wa),
        "at_x_rel": nstr(g, at_rel),
        "at_x_abs": nstr(g, at_abs),
        "meets": bool(cs.meets(target, wr, wa)),
    }


# ---------------------------------------------------------------------------
# reading the records


def control_cells():
    data = json.loads(CONTROL.read_text())
    out = {}
    for c in data["cells"]:
        out[(c["width_name"], c["target"])] = c
    return data, out


def control_accepted(cell):
    """The accepted row of the existing control, at its 60-digit count.

    Cells 2 to 7 carry the count60_* scalars.  Cell 1 was written before those
    fields existed and keeps the accepted row inside its top level.
    """
    if cell.get("count60_m") is not None:
        return {
            "m": cell["count60_m"], "k": cell["count60_k"],
            "solver": cell["count60_solver"],
            "pair": cell.get("count60_pair"),
            "rel60": cell.get("count60_rel60"),
            "abs60": cell.get("count60_abs60"),
            "rel64": cell.get("count60_rel64"),
        }
    for level in cell.get("levels", []):
        if level.get("stored") != cell.get("count60"):
            continue
        for p in (level.get("passes") or [level]):
            for r in p.get("rows", []):
                if r.get("meets60"):
                    return {
                        "m": r["m"], "k": r["k"], "solver": r["solver"],
                        "pair": f"[{r['m']}/{r['k']}] {r['solver']}",
                        "rel60": r["rel60"], "abs60": r["abs60"],
                        "rel64": r["rel64"],
                    }
    return None


def sweep_rows():
    data = json.loads(SWEEP.read_text())
    return {(r["width_name"], r["target"]): r for r in data["rows"]}


def recert_cells():
    data = json.loads(RECERT.read_text())
    return data, {(c["width_name"], c["target"]): c for c in data["cells"]}


def cheb_records():
    """The seven recount cells, every mode of each kept separately.

    The recount writes one cell per (width, target, mode).  The mode named
    "records" is the paper's degree ladder, 4 to 64; "all" adds degrees 0 to 2.
    A cell present in both carries the same count in both.
    """
    data = json.loads(CHEB_RECOUNT.read_text())
    wanted = set(CELLS)
    out = {}
    for c in data["cells"]:
        for name in c.get("width_names", []):
            key = (name, c["target"])
            if key not in wanted:
                continue
            slot = out.setdefault(key, {})
            if c["mode"] in slot:
                raise RuntimeError(f"duplicate recount cell {key} {c['mode']}")
            slot[c["mode"]] = c
    missing = sorted(wanted - set(out))
    if missing:
        raise RuntimeError(f"recount record lacks cells: {missing}")
    return data, out


def check(name, ok, detail=None):
    return {"check": name, "ok": bool(ok), "detail": detail}


# ---------------------------------------------------------------------------
# part 1: the rational cells with the rounded-coefficient column


def rational_part(cs, g, rc, cells, checkpoint):
    widths = dict(cs.width_set(g, "all"))
    targets = {cs.target_tag(t): t for t in cs.parse_targets(
        "rel:1e-08,rel:1e-14,rel:1e-15")}
    _, control_by_key = control_cells()
    record_cells = {c["id"]: c for c in rc.load_record_cells()}
    for index, (name, tag) in enumerate(CELLS, 1):
        cell = record_cells[index]
        if (cell["width_name"], cell["target"]) != (name, tag):
            raise RuntimeError(f"cell {index} is not {name} {tag}")
        control = control_by_key[(name, tag)]
        accepted = control_accepted(control)
        if accepted is None:
            raise RuntimeError(f"no accepted fit in the control for cell {index}")
        m, k = accepted["m"], accepted["k"]
        solver = accepted["solver"]
        target = targets[tag]
        width = widths[name]
        a = g.X0
        b = a + width
        ts, xs = cs.make_nodes(g, a, b, cs.N_NODES)
        started = time.perf_counter()
        ref = cs.reference(g, xs)
        reference_seconds = time.perf_counter() - started
        full_idx = list(range(len(ts)))
        fits = []
        if solver == "remez":
            # The control confirmed this cell on the minimax route; a Lawson
            # solve would add nothing the control does not already record.
            try:
                r = g.rat_remez(m, k, full_idx, ts, ref, None,
                                iters=cs.RAT_REMEZ_ITERS)
            except Exception:
                r = None
            if r is not None:
                fits.append(("remez", r[1], r[2]))
        else:
            # The control ran this cell's Lawson route through the worker job,
            # which ships the grid as 62-digit strings and parses it there.  At
            # four of the 2,001 nodes that round trip returns the neighbouring
            # 60-digit value, and the Lawson map amplifies the difference, so an
            # in-process solve is not the control's accepted fit.  Re-run the
            # control's own job; the accepted candidate is the Lawson fit.
            spec = (cs.BOYS_ROOT_STR, m, k, full_idx, cs.ship_text(g, ts),
                    cs.ship_text(g, ref), cs.RAT_REFINE_OUTER,
                    cs.RAT_REMEZ_ITERS)
            fits = cs.rat_candidate_job(spec)
        solve_seconds = time.perf_counter() - started - reference_seconds
        chosen = [f for f in fits if f[0] == solver]
        if not chosen:
            cells.append({
                "id": index, "width_name": name, "target": tag,
                "error": f"no {solver} fit returned for [{m}/{k}]",
                "candidates": [f[0] for f in fits],
                "seconds": round(time.perf_counter() - started, 1),
            })
            checkpoint()
            print(f"[rational {index}/7] {name} {tag}: no fit; recorded",
                  flush=True)
            continue
        _, p, q = chosen[0]
        e60 = rc.rat_errors_unrounded(g, p, q, a, b, ts, ref)
        e64r = rc.rat_errors_unrounded(g, rounded(g, p), rounded(g, q),
                                       a, b, ts, ref)
        e64 = cs.rat_errors(g, p, q, a, b, xs, ref)
        s60 = score_60(g, cs, target, *e60)
        s64r = score_60(g, cs, target, *e64r)
        s64 = score_64(g, cs, target, e64)
        checks = [
            check("e60 reproduces the control's rel60",
                  s60["rel"] == accepted["rel60"],
                  {"control": accepted["rel60"], "now": s60["rel"]}),
            check("e64 reproduces the control's rel64",
                  s64["rel"] == accepted["rel64"],
                  {"control": accepted["rel64"], "now": s64["rel"]}),
        ]
        # Where the control's count confirms the record's stored count, the
        # accepted pair is the record's pair; at the excluded 2^5 1e-14 row the
        # control found a lower count with another pair, so equality is not
        # required there and the fact is recorded instead.
        same_pair = (m, k) == tuple(cell["record_pair"])
        accepted_vs_record = {
            "record_pair": cell["record_pair"],
            "record_solver": cell.get("record_solver"),
            "start_stored": cell.get("start_stored"),
            "control_count60": control["count60"],
            "control_confirms_the_record":
                control["count60"] == cell.get("start_stored"),
            "accepted_pair_is_the_record_pair": same_pair,
        }
        if accepted_vs_record["control_confirms_the_record"] and not same_pair:
            checks.append(check("the accepted pair is the record's pair",
                                False, accepted_vs_record))
        others = []
        for fname, fp, fq in fits:
            if fname == solver:
                continue
            other64 = cs.rat_errors(g, fp, fq, a, b, xs, ref)
            others.append({
                "solver": fname,
                "e60": score_60(g, cs, target,
                                *rc.rat_errors_unrounded(g, fp, fq, a, b,
                                                         ts, ref)),
                "e64": score_64(g, cs, target, other64),
            })
        entry = {
            "id": index,
            "width_name": name,
            "width": nstr(g, width),
            "target": tag,
            "a": nstr(g, a),
            "b": nstr(g, b),
            "record": cell["record"],
            "control": {"count60": control["count60"],
                        "accepted": accepted},
            "accepted": {"m": m, "k": k, "solver": solver,
                         "stored": (m + 1) + k,
                         "pair": f"[{m}/{k}] {solver}"},
            "accepted_vs_record": accepted_vs_record,
            "e60": s60,
            "e64r": s64r,
            "e64": s64,
            "ratios": {
                "e64r_over_e60": ratio(g, e64r[0], e60[0]),
                "e64_over_e64r": (None if e64 is None
                                  else ratio(g, e64[0], e64r[0])),
                "e64_over_e60": (None if e64 is None
                                 else ratio(g, e64[0], e60[0])),
            },
            "other_candidates": others,
            "coefficients": {
                "p": [nstr(g, c, 62) for c in p],
                "q": [nstr(g, c, 62) for c in q],
            },
            "checks": checks,
            "reference_seconds": round(reference_seconds, 1),
            "solve_seconds": round(solve_seconds, 1),
            "seconds": round(time.perf_counter() - started, 1),
        }
        cells.append(entry)
        checkpoint()
        print(f"[rational {index}/7] {name} {tag}: [{m}/{k}] {solver} "
              f"e60={s60['rel']} e64r={s64r['rel']} e64={s64['rel']} "
              f"meets60/64r/64={s60['meets']}/{s64r['meets']}/{s64['meets']} "
              f"e64r/e60={entry['ratios']['e64r_over_e60']} "
              f"e64/e64r={entry['ratios']['e64_over_e64r']} "
              f"[{entry['seconds']} s]", flush=True)
        for c in checks:
            if not c["ok"]:
                print(f"    CHECK FAILED: {c['check']}: {c['detail']}",
                      flush=True)


# ---------------------------------------------------------------------------
# part 2: the Chebyshev endpoint cells


def chebyshev_part(cs, g, cells, checkpoint):
    _, cheb_by_key = cheb_records()
    targets = {cs.target_tag(t): t for t in cs.parse_targets(
        "rel:1e-08,rel:1e-14,rel:1e-15")}
    for index, (name, tag) in enumerate(CELLS, 1):
        modes = cheb_by_key[(name, tag)]
        mode = "records" if "records" in modes else "all"
        rec = modes[mode]
        target = targets[tag]
        # The recount's own construction: the width is the record's decimal
        # string; the left edge is the generator's kX0.  Measured against the
        # record, this reproduces the corrected binary64 score digit for digit.
        width = g.mpf(rec["width"])
        a = g.X0
        b = a + width
        ts, xs = cs.make_nodes(g, a, b, cs.N_NODES)
        started = time.perf_counter()
        ref = cs.reference(g, xs)
        reference_seconds = time.perf_counter() - started
        deg_record = rec["new_m"]
        # The 60-digit walk over the same degree ladder the recount used: the
        # first degree whose interpolant, evaluated at 60 digits on the exact t
        # grid, meets the target.
        walk = []
        count60 = None
        fit60 = None
        e60_fit = None
        for deg in cs.cheb_degrees(mode):
            coeffs = cheb_coeffs_60(g, a, b, deg)
            err = cheb_errors_60(g, coeffs, ts, ref)
            row = {"degree": deg, "stored": deg + 1,
                   "rel60": nstr(g, err[0]), "abs60": nstr(g, err[1]),
                   "meets60": bool(cs.meets(target, err[0], err[1]))}
            walk.append(row)
            if row["meets60"]:
                count60 = deg + 1
                fit60 = coeffs
                e60_fit = err
                break
        if fit60 is None:
            cells.append({
                "id": index, "width_name": name, "target": tag,
                "error": "no degree met at 60 digits on the mode's ladder",
                "walk": walk,
                "seconds": round(time.perf_counter() - started, 1),
            })
            checkpoint()
            print(f"[chebyshev {index}/7] {name} {tag}: no count60; recorded",
                  flush=True)
            continue
        clenshaw_rel, clenshaw_at = cheb_clenshaw_check(g, fit60, ts)
        s60 = score_60(g, cs, target, *e60_fit)
        e64r_fit = cheb_errors_60(g, rounded(g, fit60), ts, ref)
        s64r = score_60(g, cs, target, *e64r_fit)
        e64_fit = cs.cheb_errors(g, [float(c) for c in fit60], a, b, xs, ref)
        s64 = score_64(g, cs, target, e64_fit)
        # The record's own degree, re-solved at 60 digits and measured every
        # way; its binary64 score must reproduce the recount record.
        coeffs_rec = cheb_coeffs_60(g, a, b, deg_record)
        e60_rec = cheb_errors_60(g, coeffs_rec, ts, ref)
        e64r_rec = cheb_errors_60(g, rounded(g, coeffs_rec), ts, ref)
        e64_rec = cs.cheb_errors(g, [float(c) for c in coeffs_rec],
                                 a, b, xs, ref)
        s60_rec = score_60(g, cs, target, *e60_rec)
        s64r_rec = score_60(g, cs, target, *e64r_rec)
        s64_rec = score_64(g, cs, target, e64_rec)
        checks = [
            check("the record-degree binary64 score reproduces the recount",
                  s64_rec["rel"] == rec["delivered_rel"],
                  {"record": rec["delivered_rel"], "now": s64_rec["rel"],
                   "abs_record": rec["delivered_abs"],
                   "abs_now": s64_rec["abs"]}),
            check("the record degree is the record's stored count minus one",
                  deg_record + 1 == rec["new_stored"],
                  {"new_m": deg_record, "new_stored": rec["new_stored"]}),
            check("the ladder modes agree where both exist",
                  all(c["new_stored"] == rec["new_stored"]
                      for c in modes.values()),
                  {m2: c["new_stored"] for m2, c in modes.items()}),
            check("the Clenshaw recurrence equals the direct sum",
                  clenshaw_rel < g.mpf("1e-45"),
                  {"worst_rel": nstr(g, clenshaw_rel, 8)}),
        ]
        entry = {
            "id": index,
            "width_name": name,
            "width": nstr(g, width),
            "target": tag,
            "a": nstr(g, a),
            "b": nstr(g, b),
            "mode": mode,
            "record": {
                "old_stored": rec["old_stored"],
                "new_stored": rec["new_stored"],
                "degree": deg_record,
                "delivered_rel": rec["delivered_rel"],
                "delivered_abs": rec["delivered_abs"],
                "at_x_rel": rec["at_x_rel"],
                "at_x_abs": rec["at_x_abs"],
                "changed": rec.get("changed"),
                "records": rec.get("records"),
                "modes": {m2: {"old_stored": c["old_stored"],
                               "new_stored": c["new_stored"],
                               "delivered_rel": c["delivered_rel"]}
                          for m2, c in modes.items()},
                "source": str(CHEB_RECOUNT),
            },
            "count60": {"stored": count60, "degree": count60 - 1},
            "delta_vs_record": count60 - rec["new_stored"],
            "walk": walk,
            "clenshaw_check": {"worst_rel": nstr(g, clenshaw_rel, 8),
                               "at_t": nstr(g, clenshaw_at)},
            "fit60": {
                "degree": count60 - 1,
                "stored": count60,
                "e60": s60,
                "e64r": s64r,
                "e64": s64,
                "ratios": {
                    "e64r_over_e60": ratio(g, e64r_fit[0], e60_fit[0]),
                    "e64_over_e64r": (None if e64_fit is None
                                      else ratio(g, e64_fit[0], e64r_fit[0])),
                },
                "coefficients": [nstr(g, c, 62) for c in fit60],
            },
            "record_degree_fit": {
                "degree": deg_record,
                "stored": deg_record + 1,
                "e60": s60_rec,
                "e64r": s64r_rec,
                "e64": s64_rec,
                "ratios": {
                    "e64r_over_e60": ratio(g, e64r_rec[0], e60_rec[0]),
                    "e64_over_e64r": (None if e64_rec is None
                                      else ratio(g, e64_rec[0], e64r_rec[0])),
                },
                "coefficients": [nstr(g, c, 62) for c in coeffs_rec],
            },
            "checks": checks,
            "reference_seconds": round(reference_seconds, 1),
            "seconds": round(time.perf_counter() - started, 1),
        }
        cells.append(entry)
        checkpoint()
        print(f"[chebyshev {index}/7] {name} {tag} ({mode}): record "
              f"{rec['new_stored']} (deg {deg_record}) count60 {count60} "
              f"(deg {count60 - 1}) delta {count60 - rec['new_stored']:+d} "
              f"fit60 e60={s60['rel']} e64r={s64r['rel']} e64={s64['rel']} "
              f"record-fit e64={s64_rec['rel']} "
              f"clenshaw={nstr(g, clenshaw_rel, 4)} "
              f"[{entry['seconds']} s]", flush=True)
        for c in checks:
            if not c["ok"]:
                print(f"    CHECK FAILED: {c['check']}: {c['detail']}",
                      flush=True)


# ---------------------------------------------------------------------------
# part 3: the [7/8] flip


def flip_part(cs, g, ssv, out, checkpoint, solve=True):
    target = cs.parse_targets("rel:1e-13")[0]
    sweep = sweep_rows()
    row = sweep[("2^5", "rel:1e-13")]
    width = g.mpf(2) ** 5

    constructions = {
        "record_strings": {
            "a": g.mpf(row["a"]), "b": g.mpf(row["b"]),
            "source": ("the row's a and b decimal strings, parsed at 60 digits: "
                       "the stability pass's construction (stability_save_verify "
                       "pair_job and evaluate_fits)")},
        "x0_based": {
            "a": g.X0, "b": g.X0 + width,
            "source": ("a = g.X0 and b = a + width, the width being 2^5 exactly: "
                       "the sweep's and the re-certification's construction "
                       "(count_sweep.sweep_one and recertify_ladder.cell_result)")},
    }
    a_r, b_r = constructions["record_strings"]["a"], constructions["record_strings"]["b"]
    a_x, b_x = constructions["x0_based"]["a"], constructions["x0_based"]["b"]
    out["row"] = FLIP_ROW
    out["m"], out["k"] = FLIP_M, FLIP_K
    out["target"] = "rel:1e-13"
    out["width"] = nstr(g, width)
    out["constructions"] = {
        name: {"a": nstr(g, c["a"], 62), "b": nstr(g, c["b"], 62),
               "a17": nstr(g, c["a"]), "b17": nstr(g, c["b"]),
               "source": c["source"]}
        for name, c in constructions.items()}
    out["endpoint_difference"] = {
        "a_abs": nstr(g, a_x - a_r, 8),
        "a_rel": nstr(g, (a_x - a_r) / a_r, 8),
        "b_abs": nstr(g, b_x - b_r, 8),
        "b_rel": nstr(g, (b_x - b_r) / b_r, 8),
        "a_first_differing_decimal_digit": first_differing_digit(g, a_x, a_r),
        "b_first_differing_decimal_digit": first_differing_digit(g, b_x, b_r),
        "a_prints_the_same_at_17": nstr(g, a_x) == nstr(g, a_r),
        "b_prints_the_same_at_17": nstr(g, b_x) == nstr(g, b_r),
        "float_a_equal": float(a_x) == float(a_r),
        "float_b_equal": float(b_x) == float(b_r),
    }
    checkpoint()
    print(f"[flip] endpoint difference a {out['endpoint_difference']['a_abs']} "
          f"b {out['endpoint_difference']['b_abs']} "
          f"(first differing digit a "
          f"{out['endpoint_difference']['a_first_differing_decimal_digit']}, "
          f"b {out['endpoint_difference']['b_first_differing_decimal_digit']})",
          flush=True)

    grids = {}
    grid_specs = {"record_strings": (a_r, b_r), "x0_based": (a_x, b_x)}
    for cname, (ca, cb) in grid_specs.items():
        for nodes, label in ((cs.N_NODES, "default"),
                             (2 * cs.N_NODES, "doubled")):
            t0 = time.perf_counter()
            ts, xs = cs.make_nodes(g, ca, cb, nodes)
            ref = cs.reference(g, xs)
            grids[(cname, label)] = {
                "a": ca, "b": cb, "nodes": nodes, "points": len(ts),
                "ts": ts, "xs": xs, "ref": ref,
                "seconds": round(time.perf_counter() - t0, 1)}
    # The pass's doubled grid was built by interleaving the row-string default
    # grid; rebuild it, because that is the grid its 1.223e-13 was measured on.
    t0 = time.perf_counter()
    pass_row = {"a": nstr(g, a_r, 62), "b": nstr(g, b_r, 62)}
    base = grids[("record_strings", "default")]
    ts_p, xs_p, ref_p = ssv.doubled_grid(
        g, pass_row, base["ts"], base["xs"], base["ref"])
    grids[("pass_interleaved", "doubled")] = {
        "a": a_r, "b": b_r, "nodes": 2 * cs.N_NODES, "points": len(ts_p),
        "ts": ts_p, "xs": xs_p, "ref": ref_p,
        "seconds": round(time.perf_counter() - t0, 1)}
    out["grids"] = {
        f"{cname}/{label}": {"nodes": G["nodes"], "points": G["points"],
                             "reference_seconds": G["seconds"]}
        for (cname, label), G in grids.items()}
    # Do the two constructions' grids differ exactly where their endpoints do?
    diffs = {}
    for label in ("default", "doubled"):
        g1 = grids[("record_strings", label)]
        g2 = grids[("x0_based", label)]
        dx = max(abs(u - v) for u, v in zip(g1["xs"], g2["xs"]))
        dr = max(abs(u - v) for u, v in zip(g1["ref"], g2["ref"]))
        rr = max(abs(u - v) / abs(v) for u, v in zip(g1["ref"], g2["ref"]))
        diffs[label] = {
            "max_abs_dx": nstr(g, dx, 8),
            "max_abs_dref": nstr(g, dr, 8),
            "max_rel_dref": nstr(g, rr, 8),
            "xs_equal_exact": all(u == v for u, v in zip(g1["xs"], g2["xs"])),
            "xs_equal_float": all(float(u) == float(v)
                                  for u, v in zip(g1["xs"], g2["xs"])),
        }
    out["reference_difference"] = diffs
    # The interleaved doubled grid, rebuilt from each construction's own
    # default grid, must equal that construction's direct 4,001-point grid.
    interleave = {}
    for cname in grid_specs:
        row_obj = {"a": nstr(g, grid_specs[cname][0], 62),
                   "b": nstr(g, grid_specs[cname][1], 62)}
        b0 = grids[(cname, "default")]
        ts_i, xs_i, ref_i = ssv.doubled_grid(
            g, row_obj, b0["ts"], b0["xs"], b0["ref"])
        direct = grids[(cname, "doubled")]
        interleave[cname] = {
            "points": len(ts_i),
            "ts_equal": all(nstr(g, x, 40) == nstr(g, y, 40)
                            for x, y in zip(ts_i, direct["ts"])),
            "xs_equal": all(nstr(g, x, 40) == nstr(g, y, 40)
                            for x, y in zip(xs_i, direct["xs"])),
            "ref_equal": all(nstr(g, x, 40) == nstr(g, y, 40)
                             for x, y in zip(ref_i, direct["ref"])),
            "ts_equal_exact": all(x == y for x, y in zip(ts_i, direct["ts"])),
            "xs_equal_exact": all(x == y for x, y in zip(xs_i, direct["xs"])),
            "ref_equal_exact": all(x == y for x, y in zip(ref_i, direct["ref"])),
        }
    out["interleave_check"] = interleave
    out["pass_interleaved_check"] = {
        "points": len(ts_p),
        "same_as_i__record_strings_doubled_at_40_digits": all(
            nstr(g, x, 40) == nstr(g, y, 40)
            for x, y in zip(ts_p, grids[("record_strings", "doubled")]["ts"])),
    }
    checkpoint()
    print(f"[flip] interleave checks: {interleave}; pass grid 40-digit equal: "
          f"{out['pass_interleaved_check']['same_as_i__record_strings_doubled_at_40_digits']}",
          flush=True)

    def score_on_all(p, q):
        scores = {}
        for cname, label in GRID_KEYS:
            G = grids[(cname, label)]
            scores[f"{cname}/{label}"] = score_64(
                g, cs, target,
                cs.rat_errors(g, p, q, G["a"], G["b"], G["xs"], G["ref"]))
        return scores

    default_spec = {"row": FLIP_ROW, "role": "accepted", "m": FLIP_M,
                    "k": FLIP_K, "nodes": cs.N_NODES,
                    "outer": cs.RAT_REFINE_OUTER,
                    "remez_iters": cs.RAT_REMEZ_ITERS, "doubled": False}
    doubled_spec = {"row": FLIP_ROW, "role": "accepted", "m": FLIP_M,
                    "k": FLIP_K, "nodes": 2 * cs.N_NODES,
                    "outer": ssv.DOUBLED_OUTER,
                    "remez_iters": ssv.DOUBLED_REMEZ, "doubled": True}

    cached = {}
    for label, spec in (("default_settings", default_spec),
                        ("doubled_settings", doubled_spec)):
        t0 = time.perf_counter()
        path = ssv.cache_path(spec)
        entry = {"settings": spec, "cache_path": str(path),
                 "cached_utc": (datetime.fromtimestamp(
                     path.stat().st_mtime, timezone.utc).strftime(
                         "%Y-%m-%dT%H:%M:%SZ") if path.exists() else None)}
        if not path.exists():
            entry["error"] = "cache entry missing"
            cached[label] = entry
            continue
        data = json.loads(path.read_text())
        entry["fits"] = []
        for fit in data["fits"]:
            p = [g.mpf(s) for s in fit["p"]]
            q = [g.mpf(s) for s in fit["q"]]
            entry["fits"].append({
                "solver": fit["solver"],
                "p": fit["p"], "q": fit["q"],
                "scores": score_on_all(p, q),
            })
        entry["seconds"] = round(time.perf_counter() - t0, 1)
        cached[label] = entry
        for fit in entry["fits"]:
            s = fit["scores"]
            print(f"[flip] cached {label} {fit['solver']}: "
                  f"record/default {s['record_strings/default']['rel']} "
                  f"record/doubled {s['record_strings/doubled']['rel']} "
                  f"pass-interleaved {s['pass_interleaved/doubled']['rel']} "
                  f"x0/default {s['x0_based/default']['rel']} "
                  f"x0/doubled {s['x0_based/doubled']['rel']}", flush=True)
    out["cached"] = cached
    checkpoint()

    # The recorded values this experiment is asked to explain.
    pass_data = json.loads(PASS_LATEST.read_text())
    pass_default_rel = None
    for pair in pass_data["saved"].get(FLIP_ROW, []):
        if (pair.get("role") == "accepted" and pair.get("m") == FLIP_M
                and pair.get("k") == FLIP_K):
            for e in pair.get("default", []):
                if e.get("solver") == "lawson" and e.get("residual"):
                    pass_default_rel = e["residual"]["rel"]
    pass_doubled_rel = None
    for t in pass_data.get("thin_evaluated", []):
        if (t.get("row") == FLIP_ROW and t.get("role") == "accepted"
                and t.get("m") == FLIP_M and t.get("k") == FLIP_K):
            for e in t.get("fits", []):
                if e.get("solver") == "lawson" and e.get("residual"):
                    pass_doubled_rel = e["residual"]["rel"]
    recert = json.loads(RECERT.read_text())
    recert_cell = next(c for c in recert["cells"]
                       if (c["width_name"], c["target"]) == ("2^5", "rel:1e-13"))
    core = json.loads(CORE_1E13.read_text())
    core_rel = None
    for r in core["rows"]:
        if ((r["width_name"], r["target"]) == ("2^5", "rel:1e-13")
                and r.get("rat", {}).get("m") == FLIP_M
                and r["rat"].get("k") == FLIP_K):
            core_rel = r["rat"]["delivered_rel"]
    recorded = {
        "ladder_sweep_rel": row["rat"]["delivered_rel"],
        "ladder_sweep_source": "results-ladder/sweep_all.json, the row's solve",
        "core_1e13_rel": core_rel,
        "core_1e13_source": "results/sweep_rel_1e-13.json, the core run",
        "pass_default_rel": pass_default_rel,
        "pass_default_source": (f"{PASS_LATEST.name}, the saved fit solved at "
                                "the row's own settings, evaluated on the row "
                                "grid"),
        "pass_doubled_rel": pass_doubled_rel,
        "pass_doubled_source": (f"{PASS_LATEST.name}, the thin fit solved at "
                                "the doubled settings, evaluated on the pass's "
                                "interleaved doubled grid"),
        "recert_doubled_rel": recert_cell["certifying_pair"]["delivered_rel"],
        "recert_doubled_source": ("results-ladder/recertification.json cell 6, "
                                  "the cold doubled-settings solve"),
    }
    out["recorded_values"] = recorded

    # The pass's two cached fits, scored on each construction's own grid and on
    # the other construction's grid.  This is the decisive test of whether the
    # recorded flip is the grid or the solve, and it costs one evaluation per
    # grid.
    def cached_fit(label, solver="lawson"):
        for fit in cached.get(label, {}).get("fits", []):
            if fit["solver"] == solver:
                return fit
        return None

    cached_resolution = {}
    for label, own_key, rec_key in (
            ("default_settings", "record_strings/default", "pass_default_rel"),
            ("default_settings", "x0_based/default", "ladder_sweep_rel"),
            ("doubled_settings", "pass_interleaved/doubled",
             "pass_doubled_rel"),
            ("doubled_settings", "x0_based/doubled", "recert_doubled_rel")):
        fit = cached_fit(label)
        if fit is None:
            cached_resolution[f"{label}:{own_key}"] = {
                "cached_rel": None, "recorded_key": rec_key,
                "recorded_rel": recorded[rec_key], "matches": None}
            continue
        cached_resolution[f"{label}:{own_key}"] = {
            "cached_rel": fit["scores"][own_key]["rel"],
            "cached_abs": fit["scores"][own_key]["abs"],
            "cached_meets": fit["scores"][own_key]["meets"],
            "own_grid": own_key,
            "recorded_key": rec_key,
            "recorded_rel": recorded[rec_key],
            "matches": fit["scores"][own_key]["rel"] == recorded[rec_key],
        }
    out["cached_resolution"] = cached_resolution
    checkpoint()
    print(f"[flip] cached resolution: {json.dumps(cached_resolution, indent=1)}",
          flush=True)

    if not solve:
        out["serial_skipped"] = ("cached-only run: the four serial solves were "
                                 "not run")
        checkpoint()
        return

    serial = {}
    for cname in ("record_strings", "x0_based"):
        for nodes, label, outer, remez in (
                (cs.N_NODES, "default_settings", cs.RAT_REFINE_OUTER,
                 cs.RAT_REMEZ_ITERS),
                (2 * cs.N_NODES, "doubled_settings", ssv.DOUBLED_OUTER,
                 ssv.DOUBLED_REMEZ)):
            variant = "default" if nodes == cs.N_NODES else "doubled"
            key = f"{cname}/{label}"
            G = grids[(cname, variant)]
            t0 = time.perf_counter()
            idx = list(range(len(G["ts"])))
            if cname == "record_strings":
                # The pass's pair_job built the grid inside the worker from the
                # row's a/b strings; rebuilding it in this process is the same
                # solve on bitwise the same nodes.
                cands = cs.rational_candidates(g, FLIP_M, FLIP_K, idx,
                                               G["ts"], G["ref"], outer,
                                               remez)
                solve_path = ("grid from the record a/b strings, in-process "
                              "solve (the pass's pair_job shape)")
            else:
                # The sweep and the re-certification reach the solvers through
                # refine_pairs, which ships the grid as 62-digit strings and
                # parses it in the worker.  Four of the 2,001 nodes do not
                # round-trip at 60 digits, so that trip is part of the
                # recorded construction and is repeated here.
                spec = (cs.BOYS_ROOT_STR, FLIP_M, FLIP_K, idx,
                        cs.ship_text(g, G["ts"]), cs.ship_text(g, G["ref"]),
                        outer, remez)
                cands = cs.rat_candidate_job(spec)
                solve_path = ("62-digit ts/ref strings through the worker job "
                              "(refine_pairs' path)")
            entry = {
                "construction": cname,
                "solve_path": solve_path,
                "settings": {"nodes": nodes, "points": nodes + 1,
                             "outer": outer, "remez_iters": remez},
                "solve_seconds": round(time.perf_counter() - t0, 1),
                "fits": [],
            }
            for fname, p, q in cands:
                entry["fits"].append({
                    "solver": fname,
                    "p": [nstr(g, c, 62) for c in p],
                    "q": [nstr(g, c, 62) for c in q],
                    "scores": score_on_all(p, q),
                })
            serial[key] = entry
            for fit in entry["fits"]:
                s = fit["scores"]
                print(f"[flip] serial {key} {fit['solver']}: "
                      f"own-grid {s[f'{cname}/{variant}']['rel']} "
                      f"record/default {s['record_strings/default']['rel']} "
                      f"record/doubled {s['record_strings/doubled']['rel']} "
                      f"pass-interleaved {s['pass_interleaved/doubled']['rel']} "
                      f"x0/default {s['x0_based/default']['rel']} "
                      f"x0/doubled {s['x0_based/doubled']['rel']} "
                      f"[{entry['solve_seconds']} s]", flush=True)
            out["serial"] = serial
            checkpoint()

    # Do the fresh coefficients reproduce the cached fits exactly?
    reproduce = {}
    for label in ("default_settings", "doubled_settings"):
        old = {f["solver"]: f for f in cached.get(label, {}).get("fits", [])}
        for cname in ("record_strings", "x0_based"):
            fresh = {f["solver"]: f
                     for f in serial.get(f"{cname}/{label}", {}).get("fits", [])}
            for solver in sorted(set(fresh) & set(old)):
                reproduce[f"{cname}/{label}:{solver}"] = {
                    "coefficients_identical": (
                        fresh[solver]["p"] == old[solver]["p"]
                        and fresh[solver]["q"] == old[solver]["q"]),
                    "fresh_p": fresh[solver]["p"],
                    "fresh_q": fresh[solver]["q"],
                    "cached_p": old[solver]["p"],
                    "cached_q": old[solver]["q"],
                }
    out["reproduce"] = reproduce

    def fresh_fit(cname, label, solver="lawson"):
        for fit in serial.get(f"{cname}/{label}", {}).get("fits", []):
            if fit["solver"] == solver:
                return fit
        return None

    resolution = {}
    for cname, label, own_key, rec_key in (
            ("record_strings", "default_settings", "record_strings/default",
             "pass_default_rel"),
            ("record_strings", "doubled_settings", "pass_interleaved/doubled",
             "pass_doubled_rel"),
            ("x0_based", "default_settings", "x0_based/default",
             "ladder_sweep_rel"),
            ("x0_based", "doubled_settings", "x0_based/doubled",
             "recert_doubled_rel")):
        fit = fresh_fit(cname, label)
        if fit is None:
            resolution[f"{cname}/{label}"] = {
                "fresh_rel": None, "recorded_key": rec_key,
                "recorded_rel": recorded[rec_key], "matches": None}
            continue
        resolution[f"{cname}/{label}"] = {
            "fresh_rel": fit["scores"][own_key]["rel"],
            "fresh_abs": fit["scores"][own_key]["abs"],
            "fresh_meets": fit["scores"][own_key]["meets"],
            "own_grid": own_key,
            "recorded_key": rec_key,
            "recorded_rel": recorded[rec_key],
            "matches": fit["scores"][own_key]["rel"] == recorded[rec_key],
        }
    out["resolution"] = resolution
    checkpoint()
    print(f"[flip] fresh resolution: {json.dumps(resolution, indent=1)}",
          flush=True)


# ---------------------------------------------------------------------------
# part 4: the wide-end table, from the records only


def table_part(cs, g, out, checkpoint):
    sweep = sweep_rows()
    _, recert = recert_cells()
    _, control = control_cells()
    rows = []
    for name, tag in WIDE_END:
        s = sweep[(name, tag)]
        rat = s["rat"]
        c = recert.get((name, tag))
        ctl = control.get((name, tag))
        recertified = None
        if c is not None:
            pair = c["certifying_pair"]
            recertified = {
                "certified_stored": c["certified_stored"],
                "pair": f"[{pair['m']}/{pair['k']}] {pair['solver']}",
                "delivered_rel": pair["delivered_rel"],
                "delivered_abs": pair["delivered_abs"],
                "verdict": c["verdict"],
                "walk": c.get("walk"),
                "record_stored": c["record_stored"],
                "old_count_status": c.get("old_count_status"),
            }
        reason = None
        if c is None:
            if (name, tag) == ("paper", "rel:1e-15"):
                reason = ("not re-certified: the paper excludes this row as "
                          "non-minimal at stored 18 and does not report its "
                          "true minimum; the doubling run stopped after cell 3 "
                          "(recertify_ladder.py, scope.note)")
            else:
                reason = "not one of the six re-certification cells"
        rows.append({
            "width_name": name,
            "width": s.get("width"),
            "target": tag,
            "record": {
                "stored": rat.get("stored"),
                "m": rat.get("m"), "k": rat.get("k"),
                "solver": rat.get("solver"),
                "delivered_rel": rat.get("delivered_rel"),
                "delivered_abs": rat.get("delivered_abs"),
                "below_meets": rat.get("below_meets"),
                "at_x_rel": rat.get("at_x_rel"),
                "at_x_abs": rat.get("at_x_abs"),
            },
            "recertified": recertified,
            "recertified_absent_because": reason,
            "control": (None if ctl is None else {
                "count60": ctl.get("count60"),
                "pair": ctl.get("count60_pair"),
                "rel60": ctl.get("count60_rel60"),
                "rel64": ctl.get("count60_rel64"),
            }),
            "control_absent_because": (None if ctl is not None else
                                       "not one of the control's seven cells"),
        })
        print(f"[table] {name} {tag}: record {rat.get('stored')} "
              f"[{rat.get('m')}/{rat.get('k')}] {rat.get('solver')} "
              f"rel={rat.get('delivered_rel')} | recertified "
              f"{recertified['certified_stored'] if recertified else '-'} "
              f"{recertified['pair'] if recertified else '-'} "
              f"{recertified['delivered_rel'] if recertified else '-'} "
              f"verdict={recertified['verdict'] if recertified else '-'} | "
              f"control {ctl.get('count60') if ctl else '-'}", flush=True)
    out["rows"] = rows
    checkpoint()


# ---------------------------------------------------------------------------
# driver


def meta_block(cs):
    return {
        "schema": SCHEMA,
        "generated_utc": now(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "harness": str(PAPER / "count_sweep.py"),
        "harness_sha256": sha256(PAPER / "count_sweep.py"),
        "control_script": str(PAPER / "rational_unrounded_control.py"),
        "control_script_sha256": sha256(
            PAPER / "rational_unrounded_control.py"),
        "stability_script": str(PAPER / "stability_save_verify.py"),
        "stability_script_sha256": sha256(
            PAPER / "stability_save_verify.py"),
        "rational_boys": {
            "root": str(RAT_BOYS),
            "version": boys_version(RAT_BOYS),
            "generator_sha256": sha256(
                RAT_BOYS / "tools" / "gen_boys_coefficients.py"),
            "header_sha256": sha256(
                RAT_BOYS / "include" / "boys" / "boys_coefficients.hpp"),
        },
        "chebyshev_boys": {
            "root": str(CHEB_BOYS),
            "version": boys_version(CHEB_BOYS),
            "generator_sha256": sha256(
                CHEB_BOYS / "tools" / "gen_boys_coefficients.py"),
            "header_sha256": sha256(
                CHEB_BOYS / "include" / "boys" / "boys_coefficients.hpp"),
        },
        "records": {
            "sweep_all": {"path": str(SWEEP), "sha256": sha256(SWEEP)},
            "recertification": {"path": str(RECERT),
                                "sha256": sha256(RECERT)},
            "control": {"path": str(CONTROL), "sha256": sha256(CONTROL)},
            "core_1e13": {"path": str(CORE_1E13),
                          "sha256": sha256(CORE_1E13)},
            "pass_latest": {"path": str(PASS_LATEST),
                            "sha256": sha256(PASS_LATEST)},
            "chebyshev_recount": {"path": str(CHEB_RECOUNT),
                                  "sha256": sha256(CHEB_RECOUNT)},
        },
        "settings": {
            "nodes": cs.N_NODES,
            "grid_points": cs.N_NODES + 1,
            "remez_iters": cs.RAT_REMEZ_ITERS,
            "refine_outer": cs.RAT_REFINE_OUTER,
            "rat_k_max": K_MAX,
            "reference_dps": 60,
            "evaluation_dps": EVAL_DPS,
            "grid": "uniform in t = 2*(x-a)/(b-a) - 1, endpoints included",
            "acceptance": ("relative to F_0(x) at each node; the binary64 "
                           "route is the worse of the fused and separate "
                           "multiply-add evaluations"),
            "rational_chain": ("e60: 60-digit coefficients, 60-digit "
                               "evaluation on the exact t grid; e64r: the same "
                               "coefficients rounded to binary64, evaluated at "
                               "60 digits on the exact t grid; e64: binary64 "
                               "coefficients and evaluation, x -> t mapped in "
                               "binary64"),
            "chebyshev_chain": ("record count: the corrected binary64 recount; "
                                "count60: the first degree of the cell's own "
                                "mode ladder whose interpolant, evaluated at 60 "
                                "digits on the exact t grid, meets the target"),
            "doubled": {"nodes": 2 * cs.N_NODES,
                        "grid_points": 2 * cs.N_NODES + 1,
                        "remez_iters": 120, "refine_outer": 600},
        },
        "python": sys.version,
        "platform": platform.platform(),
    }


def load_previous(path: Path):
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text())
    except Exception as exc:
        print(f"warning: cannot read {path} ({exc}); starting fresh",
              flush=True)
        return {}
    if data.get("meta", {}).get("schema") != SCHEMA:
        print(f"warning: {path.name} has another schema; starting fresh",
              flush=True)
        return {}
    return data


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--part", action="append",
                    choices=["rational", "chebyshev", "flip", "table"],
                    help="run only this part; repeatable; default all")
    ap.add_argument("--out", default=str(OUT_DEFAULT))
    ap.add_argument("--force", action="store_true",
                    help="recompute parts already marked done")
    ap.add_argument("--flip-cached-only", action="store_true",
                    help="flip part: re-score the two cached fits on every "
                         "grid, but skip the four fresh serial solves")
    args = ap.parse_args(argv)
    parts = args.part or ["chebyshev", "table", "rational", "flip"]

    if str(PAPER) not in sys.path:
        sys.path.insert(0, str(PAPER))
    cs = load_module("count_sweep", PAPER / "count_sweep.py")
    rc = load_module("rational_unrounded_control",
                     PAPER / "rational_unrounded_control.py")
    ssv = load_module("stability_save_verify",
                      PAPER / "stability_save_verify.py")

    out_path = Path(args.out)
    record = load_previous(out_path)
    record["meta"] = meta_block(cs)
    record.setdefault("commands", [])
    record["commands"].append({
        "command": " ".join([sys.executable, str(Path(__file__).resolve())]
                            + list(argv if argv is not None else sys.argv[1:])),
        "utc": now(),
        "parts": parts,
    })

    def dump():
        out_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = out_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(record, indent=1), encoding="utf-8")
        tmp.replace(out_path)

    started = time.perf_counter()
    errors = []
    for part in parts:
        done_key = f"{part}_done"
        if not args.force and record.get(done_key):
            print(f"[{part}] already done in the record; reused "
                  f"(use --force to recompute)", flush=True)
            continue
        try:
            if part == "rational":
                g = cs.load_generator(RAT_BOYS)
                record["rational_cells"] = []
                rational_part(cs, g, rc, record["rational_cells"], dump)
                record["rational_done"] = True
            elif part == "chebyshev":
                g = cs.load_generator(CHEB_BOYS)
                record["chebyshev_cells"] = []
                chebyshev_part(cs, g, record["chebyshev_cells"], dump)
                record["chebyshev_done"] = True
            elif part == "flip":
                g = cs.load_generator(RAT_BOYS)
                record["flip"] = {}
                flip_part(cs, g, ssv, record["flip"], dump,
                          solve=not args.flip_cached_only)
                record["flip_done"] = not args.flip_cached_only
            elif part == "table":
                g = cs.load_generator(RAT_BOYS)
                record["wide_end_table"] = {}
                table_part(cs, g, record["wide_end_table"], dump)
                record["table_done"] = True
        except Exception:
            errors.append(part)
            record.setdefault("errors", []).append({
                "part": part,
                "traceback": traceback.format_exc(),
                "utc": now(),
            })
            print(f"[{part}] ERROR\n{traceback.format_exc()}", flush=True)
        dump()
        print(f"[{part}] written to {out_path.name}", flush=True)

    record["runtime_seconds"] = round(time.perf_counter() - started, 1)
    dump()
    print(f"wrote {out_path}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
