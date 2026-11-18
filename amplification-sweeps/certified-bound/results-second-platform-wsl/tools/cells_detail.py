"""Per-argument detail of the second-platform delivered-value differences."""
import csv
import struct
import sys
import math


def bits(v):
    return struct.unpack(">Q", struct.pack(">d", v))[0]


def load(path):
    rows = {}
    with open(path, newline="") as h:
        for r in csv.DictReader(h):
            rows[(r["band"], float(r["x"]), int(r["n"]))] = (
                float(r["recursed"]), float(r["direct"]))
    return rows


def main():
    a, b = load(sys.argv[1]), load(sys.argv[2])
    by_arg = {}
    max_ulp = (0, None)
    max_rel = (0.0, None)
    for key, (ra, da) in a.items():
        rb, db = b[key]
        for name, va, vb in (("recursed", ra, rb), ("direct", da, db)):
            if va != vb:
                u = abs(bits(va) - bits(vb))
                rel = abs(vb - va) / abs(va) if va else math.inf
                by_arg.setdefault((key[0], key[1]), []).append((name, key[2], u, rel, va, vb))
                if u > max_ulp[0]:
                    max_ulp = (u, (key, name, va, vb, rel))
                if rel > max_rel[0]:
                    max_rel = (rel, (key, name, va, vb, u))
    print("arguments with differences", len(by_arg))
    for (band, x), items in sorted(by_arg.items()):
        print("%s x=%.17g entries=%d methods=%s max_ulp=%d max_rel=%.3e"
              % (band, x, len(items),
                 sorted({i[0] for i in items}),
                 max(i[2] for i in items), max(i[3] for i in items)))
    print()
    (key, name, va, vb, rel) = max_ulp[1]
    print("largest ULP difference: %s x=%.17g n=%d %s" % (key[0], key[1], key[2], name))
    print("  deposited %.17g  second %.17g  ulp %d  relative %.3e"
          % (va, vb, max_ulp[0], rel))
    rel = max_rel[0]
    (key, name, va, vb, u) = max_rel[1]
    print("largest relative difference: %s x=%.17g n=%d %s" % (key[0], key[1], key[2], name))
    print("  deposited %.17g  second %.17g  relative %.3e  ulp %d"
          % (va, vb, rel, u))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
