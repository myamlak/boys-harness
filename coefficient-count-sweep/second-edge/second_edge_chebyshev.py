#!/usr/bin/env python3
"""The Chebyshev half of the second-left-edge sweep, for all 45 Group A cells.

This is the paper's own Chebyshev re-count instrument (the one that produced
the corrected counts at x0) run at a = 2 * X0 with every Group A width
doubled.  It writes results-second-edge/chebyshev_only.json: the same cells
and fields the full sweep will carry, without the rational half.  The full
sweep's sweep_all.json remains the record of record; this file exists so the
Chebyshev comparison is complete even while the rational fits are running.

Nothing is fitted here beyond the library's cheb_coeffs; the acceptance test
is cheb_errors, as in count_sweep.cheb_count.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

PAPER = Path(__file__).resolve().parent
sys.path.insert(0, str(PAPER))

import count_sweep as cs  # noqa: E402

DEFAULT_BOYS = Path(r"C:\Users\myaml\dev\_mainonly")
DEFAULT_OUT = PAPER / "results-second-edge" / "chebyshev_only.json"
TARGETS = ("rel:1e-13", "abs:5e-14", "rel:1e-14")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boys-root", default=str(DEFAULT_BOYS))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--scale", type=int, default=2)
    parser.add_argument("--nodes", type=int, default=cs.N_NODES)
    parser.add_argument("--targets", default=",".join(TARGETS))
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args(argv)

    boys = Path(args.boys_root)
    g = cs.load_generator(boys)
    cs.BOYS_ROOT_STR = str(boys)
    a = g.X0 * args.scale
    widths = [(name, width * args.scale)
              for name, width in cs.width_set(g, "all")]
    targets = cs.parse_targets(args.targets)
    cells = []
    started = time.perf_counter()
    for name, width in widths:
        b = a + width
        ts, xs = cs.make_nodes(g, a, b, args.nodes)
        reference_started = time.perf_counter()
        ref = cs.reference(g, xs)
        ref_seconds = time.perf_counter() - reference_started
        for target in targets:
            cheb = cs.cheb_count(g, a, b, xs, ref, target,
                                 degrees=cs.cheb_degrees("records"),
                                 verbose=args.progress)
            cell = {
                "target": cs.target_tag(target),
                "target_kind": target["kind"],
                "target_value": target["value"],
                "width_name": name,
                "width": g.mp.nstr(width, 17),
                "width_x0": g.mp.nstr(width / args.scale, 17),
                "width_f": float(width),
                "a": g.mp.nstr(a, 17),
                "b": g.mp.nstr(b, 17),
                "cheb": cheb,
                "reference_seconds": round(ref_seconds, 3),
                "settings": {"nodes": args.nodes + 1,
                             "cheb_degree_mode": "records"},
            }
            cells.append(cell)
            print(f"{cs.target_tag(target):>10} {name:<7} "
                  f"W={float(width):<12.6g} cheb={cheb.get('stored')} "
                  f"(deg={cheb.get('m')}, "
                  f"{'odd' if (cheb.get('stored') or 0) % 2 else 'even'}) "
                  f"rel={cheb.get('delivered_rel')} "
                  f"({cheb.get('seconds')} s, ref {ref_seconds:.1f} s)",
                  flush=True)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    header = boys / cs.HEADER_REL
    generator = boys / cs.GENERATOR_REL
    doc = {
        "meta": {
            "schema": 2,
            "generated_utc": __import__("datetime").datetime.now(
                __import__("datetime").timezone.utc).isoformat(),
            "boys_root": str(boys),
            "git": cs.git_rev(boys),
            "header": {"path": str(header), "sha256": cs.sha256(header)},
            "generator": {"path": str(generator),
                          "sha256": cs.sha256(generator)},
            "python": sys.version.split()[0],
            "platform": __import__("platform").platform(),
            "left_edge": g.mp.nstr(a, 17),
            "base_left_edge": g.mp.nstr(g.X0, 17),
            "width_scale": args.scale,
            "grid": {"nodes": args.nodes + 1},
            "solver": {"cheb_degree_mode": "records"},
        },
        "cells": cells,
    }
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(cells)} cells) in "
          f"{time.perf_counter() - started:.1f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
