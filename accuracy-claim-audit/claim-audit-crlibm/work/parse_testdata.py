"""Parser for crlibm's tests/*.testdata vector files.

Line shapes seen in the tree:
    <TAG> <in_hi> <in_lo> <out_hi> <out_lo>      # comment
    <TAG> 0x<in_hi> 0x<in_lo> 0x<out_hi> 0x<out_lo>
    <TAG> <x_hi> <x_lo> <y_hi> <y_lo> <out_hi> <out_lo>     # pow only
Words are high-then-low, big-endian hex.
"""
import os
import re

TD = r"C:\Users\myaml\dev\_claim-audit-crlibm\src\crlibm\tests"
MODES = {"N": "rn", "P": "ru", "M": "rd", "Z": "rz"}
WORD = re.compile(r"^(?:0x)?([0-9a-fA-F]{1,8})$")


def _words(tokens):
    """Words are 32-bit values printed without zero padding."""
    out = []
    for t in tokens:
        m = WORD.match(t)
        if not m:
            return None
        out.append(m.group(1).lower().rjust(8, "0"))
    return out


def load(fn, path=TD):
    """Yield (mode, inhex, outhex) or (mode, xhex, yhex, outhex); None-safe."""
    unparsed = 0
    with open(os.path.join(path, fn)) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if not line:
                continue
            p = line.split()
            if p[0] not in MODES:
                continue
            w = _words(p[1:])
            if w is None:
                unparsed += 1
                continue
            if len(w) == 4:
                yield MODES[p[0]], w[0] + w[1], w[2] + w[3]
            elif len(w) == 6:
                yield MODES[p[0]], w[0] + w[1], w[2] + w[3], w[4] + w[5]
            else:
                unparsed += 1
    if unparsed:
        raise SystemExit("%s: %d unparsed lines" % (fn, unparsed))


def counts(path=TD):
    """Per file: (vectors, by_mode)."""
    out = {}
    for fn in sorted(f for f in os.listdir(path) if f.endswith(".testdata")):
        rows = list(load(fn, path))
        by = {}
        for r in rows:
            by[r[0]] = by.get(r[0], 0) + 1
        out[fn[:-len(".testdata")]] = (len(rows), by)
    return out


if __name__ == "__main__":
    tot = 0
    for name, (n, by) in counts().items():
        print("%-10s %6d   %s" % (name, n, by))
        tot += n
    print("TOTAL machine-checked vector assertions:", tot)
