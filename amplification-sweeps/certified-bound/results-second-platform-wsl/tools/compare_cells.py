"""Cell-by-cell comparison of the second-platform sweep cells with the deposit.

Usage: compare_cells.py DEPOSITED SECOND_PLATFORM
Compares (band, x, n) row identity and the bit patterns of the recursed and
direct delivered doubles.  Prints differing entry counts, per-column counts,
the largest ULP difference and the differing rows (first 40).
"""
import csv
import struct
import sys


def bits(v):
    return struct.unpack(">Q", struct.pack(">d", v))[0]


def ulp(a, b):
    ia, ib = bits(a), bits(b)
    if (ia >> 63) != (ib >> 63):
        return abs(ia - ib)
    return abs(ia - ib)


def load(path):
    rows = {}
    order = []
    with open(path, newline="") as h:
        for r in csv.DictReader(h):
            key = (r["band"], float(r["x"]), int(r["n"]))
            rows[key] = (r["x"], float(r["recursed"]), float(r["direct"]))
            order.append(key)
    return rows, order


def main():
    a, oa = load(sys.argv[1])
    b, ob = load(sys.argv[2])
    print("rows", len(a), len(b), "keys equal", oa == ob)
    rec = dirn = both = 0
    max_rec = max_dir = 0
    examples = []
    args = set()
    for key in oa:
        if key not in b:
            continue
        xa, ra, da = a[key]
        xb, rb, db = b[key]
        if xa != xb:
            print("x text differs", key, xa, xb)
        dr = ra != rb
        dd = da != db
        if dr:
            rec += 1
            max_rec = max(max_rec, ulp(ra, rb))
        if dd:
            dirn += 1
            max_dir = max(max_dir, ulp(da, db))
        if dr and dd:
            both += 1
        if (dr or dd) and len(examples) < 40:
            examples.append((key, dr, ulp(ra, rb) if dr else 0,
                             dd, ulp(da, db) if dd else 0))
        if dr or dd:
            args.add((key[0], key[1]))
    print("recursed entries differ", rec, "of", len(a),
          "max ulp", max_rec)
    print("direct entries differ", dirn, "of", len(a),
          "max ulp", max_dir)
    print("cells with both entries differing", both)
    print("arguments with any differing entry", len(args))
    for key, dr, ur, dd, ud in examples:
        print("  %s x=%.17g n=%d recursed_diff=%s(%d ulp) direct_diff=%s(%d ulp)"
              % (key[0], key[1], key[2], dr, ur, dd, ud))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
