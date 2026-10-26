#!/usr/bin/env python3
"""Compare crlibm's four sinh routines with correctly rounded mpmath values.

Reads  raw/crlibm_sinh_raw.csv   (written by crlibm_sinh_witness)
Writes raw/crlibm_sinh_reference.csv
       raw/crlibm_sinh_independent_check.txt

The reference is mpmath 1.4.1 at 60 significant decimal digits.  The exact
mpf value is converted to a Fraction and rounded to binary64 under each of the
four IEEE 754 rounding directions, with the overflow rule of clause 7.4:
a magnitude above DBL_MAX overflows; round-to-nearest carries it to infinity
only at or above the midpoint 2^1024 - 2^970, while toward-zero and toward
minus-infinity carry it to DBL_MAX and toward-plus-infinity carries a positive
overflow to +inf.
"""
from fractions import Fraction
import csv
import math
from pathlib import Path
import struct

import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw"
DPS = 60

DBL_MAX_F = Fraction(2**1024 - 2**971)
T_RN_F = Fraction(2**1024 - 2**970)
DBL_MAX_MPF = mp.mpf(2) ** 1024 - mp.mpf(2) ** 971
T_RN_MPF = mp.mpf(2) ** 1024 - mp.mpf(2) ** 970
INF_BITS = 0x7FF0000000000000
DBL_MAX_BITS = 0x7FEFFFFFFFFFFFFF
MODES = ("rn", "ru", "rd", "rz")


def bits_to_float(b: int) -> float:
    return struct.unpack(">d", struct.pack(">Q", b))[0]


def float_to_bits(x: float) -> int:
    return struct.unpack(">Q", struct.pack(">d", x))[0]


def mpf_to_fraction(v: mp.mpf) -> Fraction:
    sign, man, exp, bc = v._mpf_
    if man == 0:
        return Fraction(0)
    if exp >= 0:
        fr = Fraction(man) * (1 << exp)
    else:
        fr = Fraction(man, 1 << (-exp))
    return -fr if sign else fr


def round_to_bits(v: mp.mpf, mode: str) -> int:
    """Correctly round the exact value v to binary64 under `mode`."""
    sign = -1 if v < 0 else 1
    a = abs(v)

    if a > DBL_MAX_MPF:
        if mode == "rn":
            overflow = a >= T_RN_MPF
            if overflow:
                return INF_BITS if sign > 0 else (INF_BITS | (1 << 63))
            return float_to_bits(sign * float(DBL_MAX_F))
        if mode == "rz":
            return float_to_bits(sign * float(DBL_MAX_F))
        if mode == "ru":
            return (INF_BITS if sign > 0
                    else float_to_bits(-float(DBL_MAX_F)))
        if mode == "rd":
            return (INF_BITS | (1 << 63) if sign < 0
                    else float_to_bits(float(DBL_MAX_F)))
        raise ValueError(mode)

    # No overflow, so the exact value fits a Fraction without a huge exponent.
    fr = mpf_to_fraction(v)
    a = abs(fr)
    nearest = float(a)  # CPython rounds Fraction -> float correctly
    nf = Fraction(nearest)
    lower = nearest if nf <= a else math.nextafter(nearest, 0.0)
    upper = nearest if nf >= a else math.nextafter(nearest, math.inf)
    if mode == "rn":
        mag = nearest
    elif sign > 0:
        mag = upper if mode == "ru" else lower
    else:
        mag = lower if mode == "ru" else (upper if mode == "rd" else lower)
    return float_to_bits(sign * mag)


def read_raw(path: Path):
    rows = []
    with path.open(newline="") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split(",")
            if len(parts) != 10:
                continue
            try:
                int(parts[0], 16)
            except ValueError:
                continue
            rows.append(parts)
    return rows


def read_meta(path: Path) -> dict[str, str]:
    labels = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            labels[row["input_bits"]] = row["label"]
    return labels


def main() -> None:
    mp.mp.dps = DPS
    labels = read_meta(RAW / "crlibm_sinh_inputs_meta.csv")
    raw = read_raw(RAW / "crlibm_sinh_raw.csv")
    out_rows = []
    mismatches = {m: [] for m in MODES}
    inf_rows = {m: [] for m in MODES}
    exact_by_bits = {}

    for parts in raw:
        btext = parts[0]
        x = bits_to_float(int(btext, 16))
        xm = mp.mpf(x)
        exact = mp.sinh(xm)
        exact_by_bits[btext] = exact
        expect = {m: round_to_bits(exact, m) for m in MODES}
        got = {
            "rn": int(parts[2], 16),
            "ru": int(parts[4], 16),
            "rd": int(parts[6], 16),
            "rz": int(parts[8], 16),
        }
        row = {
            "input_bits": btext,
            "x": repr(x),
            "label": labels.get(btext, ""),
            "exact_sinh_60d": mp.nstr(exact, 30),
            "exact_gt_DBL_MAX": str(abs(exact) > DBL_MAX_MPF),
        }
        for m in MODES:
            row[f"got_{m}_bits"] = f"{got[m]:016x}"
            row[f"expected_{m}_bits"] = f"{expect[m]:016x}"
            row[f"match_{m}"] = str(got[m] == expect[m])
            if got[m] != expect[m]:
                mismatches[m].append((btext, x, exact, expect[m], got[m]))
            if got[m] == INF_BITS:
                inf_rows[m].append(btext)
        out_rows.append(row)

    out = RAW / "crlibm_sinh_reference.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {len(out_rows)} rows to {out}")

    # 60-vs-100 digit stability of every rounding decision
    unstable = 0
    for btext, exact60 in exact_by_bits.items():
        mp.mp.dps = 100
        exact100 = mp.sinh(mp.mpf(bits_to_float(int(btext, 16))))
        for m in MODES:
            if round_to_bits(exact60, m) != round_to_bits(exact100, m):
                unstable += 1
        mp.mp.dps = DPS
    print(f"60-vs-100 digit rounding decisions that differ: {unstable}")

    base = "408633ce8fb9f87e"
    negbase = "c08633ce8fb9f87e"
    print("\n== named witness and its two neighbours ==")
    print(f"{'input_bits':16s} {'x':>22s} {'exact_sinh':>18s} "
          f"{'rn_expected':16s} {'rn_got':16s} {'ru_got':16s} "
          f"{'rd_expected':16s} {'rd_got':16s} {'rz_got':16s}")
    for parts in raw:
        btext = parts[0]
        if btext not in (f"{int(base,16)-1:016x}", base,
                         f"{int(base,16)+1:016x}", negbase):
            continue
        x = bits_to_float(int(btext, 16))
        exact = exact_by_bits[btext]
        exp = {m: round_to_bits(exact, m) for m in MODES}
        got = {"rn": int(parts[2], 16), "ru": int(parts[4], 16),
               "rd": int(parts[6], 16), "rz": int(parts[8], 16)}
        print(f"{btext:16s} {x!r:>22} {mp.nstr(exact, 12):>18} "
              f"{exp['rn']:016x} {got['rn']:016x} {got['ru']:016x} "
              f"{exp['rd']:016x} {got['rd']:016x} {got['rz']:016x}")

    print("\n== mismatches per routine ==")
    for m in MODES:
        print(f"sinh_{m}: {len(mismatches[m])} of {len(raw)} inputs")
        for btext, x, exact, want, got_bits in mismatches[m][:8]:
            print(f"  x={x!r} ({btext}) exact={mp.nstr(exact, 20)} "
                  f"expected={want:016x} got={got_bits:016x}")

    print("\n== sampled inputs where the routine returned +inf ==")
    for m in MODES:
        bits_list = sorted(int(b, 16) for b in inf_rows[m])
        if not bits_list:
            print(f"sinh_{m}: none")
            continue
        # contiguity in the sampled set
        contiguous = all(bits_list[i + 1] - bits_list[i] == 1
                         for i in range(len(bits_list) - 1))
        print(f"sinh_{m}: {len(bits_list)} sampled inputs, bits "
              f"{bits_list[0]:016x}..{bits_list[-1]:016x}, "
              f"sampled-contiguous={contiguous}")
        # is the returned infinity ever the correct result?
        wrong = [b for b in bits_list if any(
            bb == f"{b:016x}" for bb, *_ in mismatches[m])]
        print(f"  of these, {len(wrong)} are wrong")

    lines = []
    lines.append("crlibm sinh reference: independent checks")
    lines.append(f"mpmath {mp.__version__}, mp.dps = {DPS}")
    x_named = bits_to_float(int(base, 16))
    xm = mp.mpf(x_named)
    y = mp.sinh(xm)
    y1 = mp.diff(mp.sinh, xm, 1)
    y2 = mp.diff(mp.sinh, xm, 2)
    residual = abs(y2 - y) / max(abs(y2), abs(y))
    lines.append(f"x = 0x{base} = {x_named!r}")
    lines.append(f"exact sinh = {mp.nstr(y, 40)}")
    lines.append(f"ODE residual sinh'' - sinh = {mp.nstr(residual, 6)} "
                 f"(relative)")
    # independent identity: sinh(x) = (exp(x) - exp(-x))/2, and
    # cosh(x)^2 - sinh(x)^2 = 1 at a moderate argument where the identity is
    # resolvable at 60 digits.
    x_id = mp.mpf(1) / 3
    lhs = mp.cosh(x_id) ** 2 - mp.sinh(x_id) ** 2
    lines.append(f"cosh(1/3)^2 - sinh(1/3)^2 - 1 = {mp.nstr(lhs - 1, 6)}")
    log_sinh = mp.log(y)
    lines.append(f"log(sinh(x)) = {mp.nstr(log_sinh, 30)}")
    lines.append(f"asinh(exact) - x = {mp.nstr(mp.asinh(y) - xm, 6)}")
    text = "\n".join(lines) + "\n"
    (RAW / "crlibm_sinh_independent_check.txt").write_text(text)
    print("\n" + text)


if __name__ == "__main__":
    main()
