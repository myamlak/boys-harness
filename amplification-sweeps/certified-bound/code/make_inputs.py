"""Write the unique tested arguments as 64-bit patterns for the exp probe.

Sources: every reference argument of both sweeps, plus every per-order worst
cell's x for either metric.  The reference and per-order values are parsed by
Python's correctly rounded decimal reader, the same double strtod yields.
"""
import csv
import struct

ROOTS = [
    r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
]
OUT = r"C:/Users/myaml/dev/_a-cert-run/inputs/x_bits.txt"


def bits(x):
    return struct.unpack(">Q", struct.pack(">d", x))[0]


seen = {}
for root in ROOTS:
    with open(root + "/data/reference.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            x = float(row["x"])
            seen[bits(x)] = x
    with open(root + "/results/per_order.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            for col in ("at_x_abs", "at_x_rel"):
                x = float(row[col])
                seen[bits(x)] = x

with open(OUT, "w", newline="") as handle:
    for b in sorted(seen):
        handle.write("%016x\n" % b)

print("unique arguments", len(seen))
print("min", min(seen.values()), "max", max(seen.values()))
