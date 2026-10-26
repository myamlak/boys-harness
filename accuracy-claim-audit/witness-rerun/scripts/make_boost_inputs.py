#!/usr/bin/env python3
"""Generate the argument list for the Boost.Math Airy witness.

Writes:
  raw/boost_airy_inputs.txt       one 64-bit double pattern per line (hex),
                                  the format boost_airy_witness.cpp reads
  raw/boost_airy_inputs_meta.csv  bits, exact decimal, neighbourhood label

The named arguments are the decimal literals the paper reports (35.9 and
78.1252), the paper's widest sampled argument 78.13, and the external gate's
bi maximum 88.258434442270058.  The neighbourhood grids are dense enough
(step 0.001) that the reported error is resolved as a smooth function over the
window, and wide enough (plus or minus 2 for the coarse grid) to show whether
the named point is a local worst case.
"""
from pathlib import Path
import struct

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw"


def bits_of_decimal(text: str) -> int:
    x = float(text)
    return struct.unpack(">Q", struct.pack(">d", x))[0]


def main() -> None:
    rows = []  # (bits, label)
    named = [
        ("35.9", "named airy_ai witness (paper)"),
        ("78.1252", "named airy_bi witness (paper)"),
        ("78.13", "paper's widest sampled argument"),
        ("88.258434442270058", "external gate's airy_bi maximum (context)"),
    ]
    for text, label in named:
        b = bits_of_decimal(text)
        rows.append((b, label))
        # exact neighbours of the named double
        for delta in (-1, 1):
            rows.append((b + delta, f"{label} neighbour {delta:+d} ulp"))

    dense_windows = [
        (35.5, 36.3, 0.001, "airy_ai dense window"),
        (77.7, 78.5, 0.001, "airy_bi dense window"),
    ]
    for lo, hi, step, label in dense_windows:
        n = round((hi - lo) / step)
        for i in range(n + 1):
            text = f"{lo + i * step:.6f}"
            rows.append((bits_of_decimal(text), label))

    coarse_windows = [
        (34.0, 38.0, 0.01, "airy_ai coarse window"),
        (76.0, 80.0, 0.01, "airy_bi coarse window"),
    ]
    for lo, hi, step, label in coarse_windows:
        n = round((hi - lo) / step)
        for i in range(n + 1):
            text = f"{lo + i * step:.6f}"
            rows.append((bits_of_decimal(text), label))

    # Audit grids reconstructed from the row table in the deposit:
    #   A2  linear:65 over [0.001, 100]
    #   D1  linear:33 over [30, 35.9]
    # The A2 point i=50 is 78.12521875 and is the argument behind the paper's
    # displayed "x = 78.1252"; D1's last point is exactly 35.9.  Computing the
    # grid point in binary64 here reproduces the audit's argument, which the
    # printed four-decimal label does not.
    def grid(lo, hi, n, label):
        for i in range(n + 1):
            x = lo + i * (hi - lo) / n
            b = struct.unpack(">Q", struct.pack(">d", x))[0]
            rows.append((b, f"{label} i={i}"))

    grid(0.001, 100.0, 64, "audit A2 grid")
    grid(30.0, 35.9, 32, "audit D1 grid")
    # Pin the reconstructed A2 witness by its exact decimal as well, so the
    # argument is checkable without re-running the grid arithmetic.
    rows.append((bits_of_decimal("78.12521875"),
                 "audit A2 i=50 exact decimal (paper prints 78.1252)"))

    # de-duplicate by bit pattern, keeping the first label
    seen = {}
    for b, label in rows:
        seen.setdefault(b, label)

    RAW.mkdir(parents=True, exist_ok=True)
    inp = RAW / "boost_airy_inputs.txt"
    meta = RAW / "boost_airy_inputs_meta.csv"
    with inp.open("w", newline="\n") as f:
        f.write("# one 64-bit IEEE-754 double pattern per line, hex\n")
        for b in sorted(seen):
            f.write(f"{b:016x}\n")
    with meta.open("w", newline="\n") as f:
        f.write("input_bits,x,label\n")
        for b in sorted(seen):
            x = struct.unpack(">d", struct.pack(">Q", b))[0]
            f.write(f"{b:016x},{x!r},{seen[b]}\n")
    print(f"wrote {len(seen)} distinct arguments to {inp}")
    print(f"wrote {meta}")
    for text, label in named:
        b = bits_of_decimal(text)
        x = struct.unpack(">d", struct.pack(">Q", b))[0]
        print(f"  {label}: {text} -> {b:016x} = {x!r}")


if __name__ == "__main__":
    main()
