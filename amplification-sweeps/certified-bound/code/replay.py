"""Bit-for-bit replay of the sweep's cells from the model, and the contraction
questions the paper left to the delivered cells.

usage: python replay.py <sweep-root> [--both-downward]

Compares model.ship_batch / model.product_batch against results/cells.csv at
every (band, x, order).  With --both-downward, tries the separate and the
fused downward step and reports both mismatch counts.
"""
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model


def load_cells(root):
    cells = defaultdict(dict)
    with open(root + "/results/cells.csv", newline="") as handle:
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


def replay(root, exp):
    cells = load_cells(root)
    mism = defaultdict(int)
    total = defaultdict(int)
    first = {}
    for band, per_band in cells.items():
        bi = 0 if band == "A1" else 1
        xs = sorted({x for x, _ in per_band})
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


def main():
    root = sys.argv[1]
    both = "--both-downward" in sys.argv
    exp = model.load_exp()
    if both:
        for fused in (False, True):
            model.DOWNWARD_FUSED = fused
            out = replay(root, exp)
            print("downward fused =", fused)
            print("  total", out["total"], "mismatches", out["mism"])
            for k, v in out["first"].items():
                print("  first", k, v)
    else:
        model.DOWNWARD_FUSED = False
        out = replay(root, exp)
        print("total", out["total"], "mismatches", out["mism"])
        for k, v in out["first"].items():
            print("first", k, v)


if __name__ == "__main__":
    main()
