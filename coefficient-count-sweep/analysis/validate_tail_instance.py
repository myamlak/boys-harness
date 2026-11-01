#!/usr/bin/env python3
"""Independent Chebyshev recount for the pure-tail instance.

This is a witness, not the instrument.  Nothing here imports the kernel's
generator, its Clenshaw evaluation or the harness's fitting code: the
coefficients are formed from the textbook discrete Chebyshev transform of
x^-1/2 at the zeros of T_{d+1}, and the polynomial is summed with a plain
three-term recurrence at 70-digit precision, on the harness's 2,001-node
uniform grid.  The degree ladder and the acceptance rule are the paper's
(degrees 4..64 at both parities; a row meets a relative target when the worst
relative error on the grid does, and an absolute target when the worst
absolute error does).

Every one of the 90 Group A cells -- 15 widths x 2 edges x 3 targets -- is
recounted and compared with the recorded stored count.

Usage:
    python validate_tail_instance.py [--out results-second-instance]
Exit 0 when every recounted cell equals its recorded count.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

PAPER = Path(__file__).resolve().parent
NAMES = [f"2^{j}" for j in range(-7, 6)] + ["ship", "paper"]
TARGETS = ["rel:1e-13", "abs:5e-14", "rel:1e-14"]
X0 = mp.mpf("11.899848152108484")
A2 = 2 * X0
DPS = 70
mp.mp.dps = DPS
X0 = mp.mpf("11.899848152108484")
A2 = 2 * X0
N_NODES = 2000
# The paper's record ladder: degrees 4..64 at both parities.  Degree 3 is
# never read by the kernel's split Clenshaw; degrees 0..2 belong to the
# unrestricted path, which the record mode does not use.
LADDER = [d for d in range(4, 65)]


def width_of(name: str):
    if name.startswith("2^"):
        return mp.mpf(2) ** int(name.split("^")[1])
    if name == "ship":
        return mp.mpf("31") - X0
    if name == "paper":
        return mp.mpf("28.989337738820740") - X0
    raise ValueError(name)


def coefficients(a, b, deg):
    mid = (a + b) / 2
    half = (b - a) / 2
    n = deg + 1
    values = []
    for k in range(n):
        x = mid + half * mp.cos(mp.pi * (k + mp.mpf("0.5")) / n)
        values.append(1 / mp.sqrt(x))
    coeffs = []
    for j in range(n):
        s = mp.mpf(0)
        for k in range(n):
            s += values[k] * mp.cos(mp.pi * j * (k + mp.mpf("0.5")) / n)
        coeffs.append((1 if j == 0 else 2) * s / n)
    return coeffs


def errors(a, b, coeffs, xs, ref):
    mid = (a + b) / 2
    half = (b - a) / 2
    worst_rel = mp.mpf(0)
    worst_abs = mp.mpf(0)
    for x, f in zip(xs, ref):
        t = (x - mid) / half
        t_prev = mp.mpf(1)
        p = coeffs[0]
        if len(coeffs) > 1:
            t_cur = t
            p += coeffs[1] * t_cur
            for c in coeffs[2:]:
                t_next = 2 * t * t_cur - t_prev
                p += c * t_next
                t_prev, t_cur = t_cur, t_next
        diff = abs(p - f)
        worst_abs = max(worst_abs, diff)
        worst_rel = max(worst_rel, diff / abs(f))
    return worst_rel, worst_abs


def meets(target, wr, wa):
    kind, value = target.split(":", 1)
    return wr <= mp.mpf(value) if kind == "rel" else wa <= mp.mpf(value)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(PAPER / "results-second-instance"))
    args = parser.parse_args(argv)
    mp.mp.dps = DPS

    recorded = {}
    cell_path = Path(args.out) / "chebyshev_only.json"
    doc = json.loads(cell_path.read_text(encoding="utf-8"))
    for cell in doc["cells"]:
        recorded[(cell["a"], cell["target"], cell["width_name"])] = \
            cell["cheb"]["stored"]

    edges = [(X0, "x0"), (A2, "a2")]
    mismatches = []
    rows = []
    for a, label in edges:
        for name in NAMES:
            width = width_of(name) * (a / X0)
            b = a + width
            xs = [a + (b - a) * (mp.mpf(i) / N_NODES) for i in range(N_NODES + 1)]
            ref = [1 / mp.sqrt(x) for x in xs]
            accepted = {}
            for deg in LADDER:
                coeffs = coefficients(a, b, deg)
                wr, wa = errors(a, b, coeffs, xs, ref)
                for target in TARGETS:
                    if target in accepted:
                        continue
                    if meets(target, wr, wa):
                        accepted[target] = {
                            "stored": deg + 1, "m": deg,
                            "delivered_rel": mp.nstr(wr, 17),
                            "delivered_abs": mp.nstr(wa, 17)}
                if len(accepted) == len(TARGETS):
                    break
            for target in TARGETS:
                got = accepted.get(target, {}).get("stored")
                want = recorded.get((mp.nstr(a, 17), target, name))
                ok = got == want
                if not ok:
                    mismatches.append({"edge": label, "width": name,
                                       "target": target,
                                       "recorded": want, "recount": got})
                rows.append({"edge": label, "width": name, "target": target,
                             "recorded": want, "recount": got, "match": ok,
                             **{k: accepted.get(target, {}).get(k)
                                for k in ("m", "delivered_rel",
                                          "delivered_abs")}})
            print(f"{label} {name:>6}: " + " ".join(
                f"{t}={accepted.get(t, {}).get('stored')}"
                f"/{recorded.get((mp.nstr(a, 17), t, name))}"
                for t in TARGETS), flush=True)

    out = Path(args.out) / "validation_recount.json"
    out.write_text(json.dumps({
        "witness": "independent 70-digit Chebyshev transform and recurrence",
        "grid": N_NODES + 1,
        "ladder": [LADDER[0], LADDER[-1]],
        "cells": len(rows),
        "mismatches": mismatches,
        "rows": rows,
    }, indent=1) + "\n", encoding="utf-8")
    print(f"\n{len(rows)} cells recounted; {len(mismatches)} mismatch(es)")
    for item in mismatches:
        print("  MISMATCH", item)
    print(f"wrote {out}")
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
