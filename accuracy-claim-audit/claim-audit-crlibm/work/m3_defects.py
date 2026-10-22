"""M3: map the exact extent of the two measured defects.

sinh:  directed-rounding overflow near the first input whose sinh overflows.
atanpi: directed rounding where the exact value approaches 1/2 from below.
"""
import sys
import os
import struct
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref
import run_crlibm


def d2hex(x):
    return struct.pack(">d", x).hex()


def bits_to_float(b):
    return struct.unpack(">d", struct.pack(">Q", b))[0]


def scan(func, bitlist, modes=("rn", "ru", "rd", "rz")):
    queries, labels = [], []
    for b in bitlist:
        x = bits_to_float(b)
        h = d2hex(x)
        for m in modes:
            queries.append((func, m, h))
            labels.append((b, m))
    got = run_crlibm.evaluate(queries)
    bad = []
    for (b, m), g in zip(labels, got):
        x = bits_to_float(b)
        rn, ru, rd, rz, st = ref.reference(func, d2hex(x))
        want = {"rn": rn, "ru": ru, "rd": rd, "rz": rz}[m]
        if st != "decided":
            continue
        if g != want:
            bad.append((b, m, want, g))
    return bad


def main():
    print("== sinh: dense scan around the overflow boundary ==")
    base = 0x408633ce8fb9f870
    bl = list(range(base - 4, base + 16))
    bad = scan("sinh", bl)
    print("  inputs scanned:", len(bl), " wrong results:", len(bad))
    for b, m, want, g in bad[:12]:
        print("    x=%016x sinh_%s want=%s got=%s" % (b, m, want, g))
    print("== sinh: negative mirror ==")
    nbase = base | 0x8000000000000000
    bl2 = [nbase - 4 + i for i in range(20)]
    bad2 = scan("sinh", bl2)
    print("  wrong results:", len(bad2))
    for b, m, want, g in bad2[:12]:
        print("    x=%016x sinh_%s want=%s got=%s" % (b, m, want, g))
    print("== sinh: coarse scan across the whole overflow region ==")
    bl3 = [0x408633ce8fb9f87e + i * 0x1000000000 for i in range(1, 40)]
    bl3 += [0x7fefffffffffffff - i * 0x100000000000 for i in range(1, 40)]
    bad3 = scan("sinh", bl3)
    print("  inputs:", len(bl3), " wrong results:", len(bad3))
    for b, m, want, g in bad3[:8]:
        print("    x=%016x sinh_%s want=%s got=%s" % (b, m, want, g))

    print("== atanpi: threshold scan ==")
    # smallest x with atan(x)/pi > 0.5 - 2**-54
    import mpmath
    mpmath.mp.dps = 60
    xb = mpmath.tan(mpmath.pi * (mpmath.mpf(1) / 2 - mpmath.mpf(2) ** -54))
    xf = float(xb)
    print("  exact boundary x* =", mpmath.nstr(xb, 20), "=", xf.hex())
    bl4 = [struct.unpack(">Q", struct.pack(">d", xf))[0] + i for i in range(-3, 4)]
    bad4 = scan("atanpi", bl4)
    print("  inputs scanned:", len(bl4), " wrong results:", len(bad4))
    for b, m, want, g in bad4:
        print("    x=%016x atanpi_%s want=%s got=%s" % (b, m, want, g))
    print("  documented proven interval for atanpi: [tan(2^-25 pi), tan(2^-5 pi)] =")
    lo = mpmath.tan(mpmath.mpf(2) ** -25 * mpmath.pi)
    hi = mpmath.tan(mpmath.mpf(2) ** -5 * mpmath.pi)
    print("     [%s, %s]" % (mpmath.nstr(lo, 12), mpmath.nstr(hi, 12)))


if __name__ == "__main__":
    main()
