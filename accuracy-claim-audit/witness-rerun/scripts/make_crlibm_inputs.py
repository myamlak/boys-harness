#!/usr/bin/env python3
"""Generate the input list for the crlibm sinh witness.

Writes:
  raw/crlibm_sinh_inputs.txt       one 64-bit double pattern per line (hex)
  raw/crlibm_sinh_inputs_meta.csv  bits, x, label

The named witness is the crlibm constant max_input_csh = 0x408633ce8fb9f87e
(csh_fast.h), the double nearest asinh(DBL_MAX).  The list covers its 64
neighbours on each side, the negative mirror, and a coarse scan through the
whole positive and negative overflow tail up to DBL_MAX.
"""
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw"

BASE = 0x408633CE8FB9F87E
SIGN = 1 << 63
DBL_MAX_BITS = 0x7FEFFFFFFFFFFFFF


def bits_to_float(b: int) -> float:
    return struct.unpack(">d", struct.pack(">Q", b))[0]


def main() -> None:
    rows = []  # (bits, label)

    for d in range(-64, 65):
        rows.append((BASE + d, f"positive neighbour {d:+d}"))
    for d in range(-64, 65):
        rows.append((SIGN | (BASE + d), f"negative neighbour {d:+d}"))

    # Coarse scans through the overflow tail, in the form used by the
    # audit's own m3_defects.py.
    for k in range(1, 41):
        b = BASE + k * 0x1000000000
        rows.append((b, f"positive tail k={k}"))
        rows.append((SIGN | b, f"negative tail k={k}"))
    for i in range(1, 41):
        b = DBL_MAX_BITS - i * 0x100000000000
        rows.append((b, f"positive near DBL_MAX i={i}"))
        rows.append((SIGN | b, f"negative near DBL_MAX i={i}"))

    # A few further points below the boundary and the two format endpoints.
    for b, label in (
        (DBL_MAX_BITS, "DBL_MAX"),
        (BASE - 0x100000000, "base - 2^32"),
        (BASE - 0x10000000000, "base - 2^40"),
        (0x3FF0000000000000, "1.0 sanity point"),
    ):
        rows.append((b, label))

    seen = {}
    for b, label in rows:
        try:
            bits_to_float(b)
        except Exception:
            continue
        seen.setdefault(b, label)

    RAW.mkdir(parents=True, exist_ok=True)
    with (RAW / "crlibm_sinh_inputs.txt").open("w", newline="\n") as f:
        f.write("# one 64-bit IEEE-754 double pattern per line, hex\n")
        for b in sorted(seen):
            f.write(f"{b:016x}\n")
    with (RAW / "crlibm_sinh_inputs_meta.csv").open("w", newline="\n") as f:
        f.write("input_bits,x,label\n")
        for b in sorted(seen):
            f.write(f"{b:016x},{bits_to_float(b)!r},{seen[b]}\n")

    print(f"wrote {len(seen)} distinct inputs")
    print(f"named witness {BASE:016x} = {bits_to_float(BASE)!r}")
    print(f"negative mirror {SIGN | BASE:016x} = "
          f"{bits_to_float(SIGN | BASE)!r}")


if __name__ == "__main__":
    main()
