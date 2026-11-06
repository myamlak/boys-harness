"""Bit-for-bit replay of the sweep's cells from the model, and the contraction
questions the paper left to the delivered cells.

usage: python replay.py [<sweep-root>] [--both-downward] [--limit N]

Compares model.ship_batch / model.product_batch against results/cells.csv at
every (band, x, order).  With --both-downward, tries the separate and the
fused downward step and reports both mismatch counts.  --limit N replays the
first N arguments of each band, for a clean-checkout smoke test.  The sweep
root defaults to --sweep / CERT_SWEEP (record_paths.py).
"""
import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model
import record_paths


def load_cells(root):
    cells = defaultdict(dict)
    with open(Path(root) / "results" / "cells.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            x = float(row["x"])
            cells[row["band"]][(x, int(row["n"]))] = (
                float(row["recursed"]), float(row["direct"]))
    return cells


def bits(v):
    import struct
    return struct.unpack(">Q", struct.pack(">d", v))[0]


def ulps(a, b):
    ia, ib = bits(a), bits(b)
    if (ia >> 63) == 1:
        ia = (1 << 64) - ia
    if (ib >> 63) == 1:
        ib = (1 << 64) - ib
    return abs(ia - ib)


def replay(root, exp, limit=None):
    cells = load_cells(root)
    mism = defaultdict(int)
    total = defaultdict(int)
    first = {}
    for band, per_band in cells.items():
        bi = 0 if band == "A1" else 1
        xs = sorted({x for x, _ in per_band})
        if limit is not None:
            xs = xs[:limit]
        for x in xs:
            half = exp[x][1]
            ship = model.ship_batch(x, half)
            prod = model.product_batch(x, bi)
            for n in range(33):
                got_s, got_p = per_band[(x, n)]
                for kind, got, mod in (("ship", got_s, ship[n]["value"]),
                                       ("product", got_p, prod[n]["value"])):
                    total[kind] += 1
                    if bits(got) != bits(mod):
                        mism[kind] += 1
                        if kind not in first:
                            first[kind] = (band, x, n, got, mod, ulps(got, mod))
    return dict(total=dict(total), mism=dict(mism),
                first={k: (v[0], v[1], v[2], v[3], v[4], v[5])
                       for k, v in first.items()})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sweep_positional", nargs="?", metavar="SWEEP",
                        help="sweep root (default: --sweep / CERT_SWEEP)")
    record_paths.add_kernel(parser)
    record_paths.add_inputs(parser, exp_bits_arg=True)
    parser.add_argument("--both-downward", action="store_true",
                        help="try the separate and the fused downward step")
    parser.add_argument("--limit", type=int, metavar="N",
                        help="replay the first N arguments of each band")
    args = parser.parse_args(argv)
    model.configure(kernel=args.kernel, exp_bits=args.exp_bits)
    root = record_paths.sweep(args.sweep or args.sweep_positional)
    exp = model.load_exp()
    if args.both_downward:
        for fused in (False, True):
            model.DOWNWARD_FUSED = fused
            out = replay(root, exp, args.limit)
            print("downward fused =", fused)
            print("  total", out["total"], "mismatches", out["mism"])
            for k, v in out["first"].items():
                print("  first", k, v)
    else:
        model.DOWNWARD_FUSED = False
        out = replay(root, exp, args.limit)
        print("total", out["total"], "mismatches", out["mism"])
        for k, v in out["first"].items():
            print("first", k, v)
    return 0


if __name__ == "__main__":
    sys.exit(main())


if __name__ == "__main__":
    main()
