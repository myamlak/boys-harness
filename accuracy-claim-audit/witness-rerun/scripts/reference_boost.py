#!/usr/bin/env python3
"""Compare the Boost.Math Airy witness output with a 60-digit mpmath reference.

Reads  raw/boost_airy_raw.csv   (written by boost_airy_witness.exe)
Writes raw/boost_airy_reference.csv  (one row per function and argument)
       raw/boost_independent_check.txt

The relative error is reported in units of eps = 2**-52, i.e.
    eps_units = |library - reference| / |reference| / 2**-52.

The reference is mpmath 1.4.1 at 60 significant decimal digits.  Its argument
is the exact binary64 value of the input, never the printed decimal.  Two
independent checks are run at x = 35.9: the Airy ODE residual y'' - x y and the
Wronskian Ai Bi' - Ai' Bi = 1/pi.  A 60-vs-100 digit stability check confirms
that the reference is far inside half an ulp at the reported points.
"""
from pathlib import Path
import csv
import struct

import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw"
DPS = 60
EPS = mp.mpf(2) ** -52

NAMED = {
    "35.9": 121.395,
    "78.1252": None,
    "78.12521875": 414.111,
    "78.13": None,
    "88.258434442270058": None,
}


def bits_to_float(b: str) -> float:
    return struct.unpack(">d", struct.pack(">Q", int(b, 16)))[0]


def read_meta(path: Path) -> dict[str, str]:
    labels = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            labels[row["input_bits"]] = row["label"]
    return labels


def read_raw(path: Path):
    rows = []
    with path.open(newline="") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split(",")
            if len(parts) != 6:
                continue
            try:
                int(parts[0], 16)
            except ValueError:
                continue
            rows.append(parts)
    return rows


def error_eps(result_bits: str, reference: mp.mpf) -> tuple[float, float]:
    got = mp.mpf(bits_to_float(result_bits))
    err = abs(got - reference) / abs(reference) / EPS
    rel = abs(got - reference) / abs(reference)
    return float(err), float(rel)


def main() -> None:
    mp.mp.dps = DPS
    raw = read_raw(RAW / "boost_airy_raw.csv")
    labels = read_meta(RAW / "boost_airy_inputs_meta.csv")
    out_rows = []
    per_function = {}
    for input_bits, x_text, ai_bits, ai_text, bi_bits, bi_text in raw:
        x = bits_to_float(input_bits)
        xm = mp.mpf(x)
        refs = {
            "airy_ai": mp.airyai(xm),
            "airy_bi": mp.airybi(xm),
        }
        for name, ref in refs.items():
            res_bits = ai_bits if name == "airy_ai" else bi_bits
            err, rel = error_eps(res_bits, ref)
            out_rows.append(
                {
                    "input_bits": input_bits,
                    "x": repr(x),
                    "label": labels.get(input_bits, ""),
                    "function": name,
                    "result_bits": res_bits,
                    "result": repr(bits_to_float(res_bits)),
                    "reference_60d": mp.nstr(ref, 30),
                    "rel_error": f"{rel:.17g}",
                    "error_eps": f"{err:.12g}",
                }
            )
            per_function.setdefault(name, []).append((err, x, input_bits))

    out = RAW / "boost_airy_reference.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {len(out_rows)} rows to {out}")

    print("\n== named arguments: error in eps = 2^-52 ==")
    print(f"{'function':8s} {'x':>22s} {'input_bits':16s} {'error_eps':>12s}"
          f" {'paper_eps':>10s}")
    named_hits = {}
    for name in ("airy_ai", "airy_bi"):
        for err, x, bits in per_function[name]:
            key = None
            for text in NAMED:
                if bits == f"{struct.unpack('>Q', struct.pack('>d', float(text)))[0]:016x}":
                    key = text
            if key is None:
                continue
            named_hits[(name, key)] = (err, x, bits)
            paper = NAMED[key]
            paper_s = "n/a" if paper is None else f"{paper:.3f}"
            print(f"{name:8s} {x!r:>22} {bits:16s} {err:12.6f} {paper_s:>10s}")

    print("\n== worst in dense windows (step 0.001) ==")
    windows = {
        "airy_ai": (35.4, 36.4),
        "airy_bi": (77.6, 78.6),
    }
    for name, (lo, hi) in windows.items():
        inside = [(e, x, b) for e, x, b in per_function[name] if lo <= x <= hi]
        e, x, b = max(inside)
        print(f"{name:8s} window [{lo}, {hi}]: max {e:.6f} eps at x={x!r}"
              f" ({b})")
        # named-point value for comparison
        for text in NAMED:
            nb = f"{struct.unpack('>Q', struct.pack('>d', float(text)))[0]:016x}"
            for e2, x2, b2 in inside:
                if b2 == nb:
                    print(f"         named {text}: {e2:.6f} eps at x={x2!r}")
        # worst within 0.5 of the named point
        named_text = "35.9" if name == "airy_ai" else "78.1252"
        c = float(named_text)
        near = [(e3, x3, b3) for e3, x3, b3 in inside if abs(x3 - c) <= 0.5]
        e3, x3, b3 = max(near)
        print(f"         max within +/-0.5 of {named_text}: {e3:.6f} eps at"
              f" x={x3!r} ({b3})")

    print("\n== reconstructed audit grids (row table in the deposit) ==")
    for grid_label in ("audit A2 grid", "audit D1 grid"):
        for name in ("airy_ai", "airy_bi"):
            grid_rows = [
                (e, x, b) for e, x, b in per_function[name]
                if labels.get(b, "").startswith(grid_label)
            ]
            e, x, b = max(grid_rows)
            print(f"{name:8s} {grid_label}: max {e:.6f} eps at x={x!r}"
                  f" ({b}) [{labels.get(b, '')}]")

    print("\n== overall maxima in the sampled set ==")
    for name in ("airy_ai", "airy_bi"):
        e, x, b = max(per_function[name])
        print(f"{name:8s}: {e:.6f} eps at x={x!r} ({b})")

    # ---------------------------------------------------------------- checks
    lines = []
    lines.append("Boost.Math Airy reference: independent checks")
    lines.append(f"mpmath {mp.__version__}, mp.dps = {DPS}")
    x_named = bits_to_float(
        f"{struct.unpack('>Q', struct.pack('>d', 35.9))[0]:016x}")
    lines.append(f"x = 35.9 exact double = {x_named!r}, bits "
                 f"{struct.unpack('>Q', struct.pack('>d', x_named))[0]:016x}")
    for label, fn in (("Ai", mp.airyai), ("Bi", mp.airybi)):
        y = fn(mp.mpf(x_named))
        y1 = mp.diff(fn, mp.mpf(x_named), 1)
        y2 = mp.diff(fn, mp.mpf(x_named), 2)
        residual = abs(y2 - mp.mpf(x_named) * y) / max(
            abs(y2), abs(mp.mpf(x_named) * y))
        lines.append(f"ODE residual {label}'' - x {label}: |res| = "
                     f"{mp.nstr(residual, 6)} (relative)")
    ai = mp.airyai(mp.mpf(x_named))
    bi = mp.airybi(mp.mpf(x_named))
    ai1 = mp.diff(mp.airyai, mp.mpf(x_named), 1)
    bi1 = mp.diff(mp.airybi, mp.mpf(x_named), 1)
    wronskian = ai * bi1 - ai1 * bi
    wres = abs(wronskian - 1 / mp.pi) / (1 / mp.pi)
    lines.append(f"Wronskian Ai Bi' - Ai' Bi = {mp.nstr(wronskian, 30)}")
    lines.append(f"Wronskian residual vs 1/pi = {mp.nstr(wres, 6)} (relative)")

    # 60-vs-100 digit stability at the named points and window maxima
    check_points = [float(t) for t in NAMED]
    for name, (lo, hi) in windows.items():
        inside = [(e, x, b) for e, x, b in per_function[name] if lo <= x <= hi]
        check_points.append(max(inside)[1])
    lines.append("")
    lines.append("60-vs-100 digit stability (does the reference round to the "
                 "same binary64?)")
    stable = True
    for x in sorted(set(check_points)):
        xm = mp.mpf(x)
        mp.mp.dps = 60
        a60 = mp.airyai(xm)
        b60 = mp.airybi(xm)
        mp.mp.dps = 100
        a100 = mp.airyai(xm)
        b100 = mp.airybi(xm)
        same_a = float(a60) == float(a100)
        same_b = float(b60) == float(b100)
        stable = stable and same_a and same_b
        lines.append(f"  x={x!r}: Ai same={same_a} Bi same={same_b}")
    mp.mp.dps = DPS
    lines.append(f"all stable: {stable}")

    text = "\n".join(lines) + "\n"
    (RAW / "boost_independent_check.txt").write_text(text)
    print("\n" + text)


if __name__ == "__main__":
    main()
