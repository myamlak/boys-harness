"""Compare the second-platform exponential probe bits with the deposited bits.

Usage: compare_exp.py DEPOSITED SECOND_PLATFORM
Prints the argument count, how many exp(-x) and half-scaling bit patterns
differ, the largest absolute integer-bit difference of the exponentials, and
lists the differing arguments.
"""
import struct
import sys


def load(path):
    out = {}
    with open(path) as h:
        for line in h:
            xb, eb, hb = (int(t, 16) for t in line.split())
            out[xb] = (eb, hb)
    return out


def main():
    a, b = load(sys.argv[1]), load(sys.argv[2])
    only_a = set(a) - set(b)
    only_b = set(b) - set(a)
    diff_e = diff_h = 0
    max_bitdiff = 0
    rows = []
    for xb in sorted(set(a) & set(b)):
        ea, ha = a[xb]
        eb, hb = b[xb]
        if ea != eb:
            diff_e += 1
            d = abs(ea - eb)
            max_bitdiff = max(max_bitdiff, d)
            x = struct.unpack(">d", struct.pack(">Q", xb))[0]
            rows.append((x, ea, eb, d))
        if ha != hb:
            diff_h += 1
    print("arguments", len(a), len(b), "only_deposited", len(only_a),
          "only_second", len(only_b))
    print("exp bits differ", diff_e, "half bits differ", diff_h,
          "max bit difference", max_bitdiff)
    for x, ea, eb, d in rows[:40]:
        print("  x=%.17g deposited=%016x second=%016x diff=%d" % (x, ea, eb, d))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
