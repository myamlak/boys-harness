"""Write the unique tested arguments as 64-bit patterns for the exp probe.

Sources: every reference argument of the deposited sweeps, plus every
per-order worst cell's x for either metric.  The reference and per-order
values are parsed by Python's correctly rounded decimal reader, the same
double strtod yields.  A reference file may be deposited as reference.csv or
reference.csv.xz.

usage: python make_inputs.py --out inputs/x_bits.txt [--sweep DIR ...]

The sweep roots default to the two deposited siblings,
../amplification-form and ../product-amplification (record_paths.py).
"""
import argparse
import csv
import lzma
import os
import struct
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import record_paths

DEFAULT_SWEEPS = [record_paths.RECORD_DIR.parent / "amplification-form",
                  record_paths.RECORD_DIR.parent / "product-amplification"]


def bits(x):
    return struct.unpack(">Q", struct.pack(">d", x))[0]


def open_reference(root):
    plain = root / "data" / "reference.csv"
    if plain.is_file():
        return open(plain, newline="", encoding="utf-8")
    packed = root / "data" / "reference.csv.xz"
    if packed.is_file():
        return lzma.open(packed, "rt", newline="", encoding="utf-8")
    raise SystemExit(f"no reference.csv or reference.csv.xz under {root}/data")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sweep", action="append", metavar="DIR",
                        help="sweep root; repeat for several (default: the "
                             "two deposited siblings)")
    parser.add_argument("--out", required=True, metavar="FILE",
                        help="output bit file")
    args = parser.parse_args(argv)
    roots = [Path(v).expanduser().resolve() for v in args.sweep] \
        if args.sweep else DEFAULT_SWEEPS

    seen = {}
    for root in roots:
        with open_reference(root) as handle:
            for row in csv.DictReader(handle):
                x = float(row["x"])
                seen[bits(x)] = x
        with open(root / "results" / "per_order.csv", newline="") as handle:
            for row in csv.DictReader(handle):
                for col in ("at_x_abs", "at_x_rel"):
                    x = float(row[col])
                    seen[bits(x)] = x

    out = Path(args.out).expanduser().resolve()
    with open(out, "w", newline="") as handle:
        for b in sorted(seen):
            handle.write("%016x\n" % b)

    print("unique arguments", len(seen))
    print("min", min(seen.values()), "max", max(seen.values()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
