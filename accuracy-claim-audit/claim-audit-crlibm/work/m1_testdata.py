"""M1: run crlibm against its own shipped test vectors.

tests/*.testdata hold input/output pairs for the four rounding modes; the
expected outputs were produced by the crlibm authors with MPFR.  Each line is
a pointwise assertion that crlibm returns the correctly rounded double, so the
files are the library's machine-checked bounds -- published in the source
tree, not in the manual.

This measures whether crlibm satisfies them.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parse_testdata
import run_crlibm


def main():
    files = sorted(f for f in os.listdir(parse_testdata.TD)
                   if f.endswith(".testdata"))
    total = bad_total = 0
    rows = []
    for fn in files:
        name = fn[:-len(".testdata")]
        cases = list(parse_testdata.load(fn))
        queries, expect = [], []
        for c in cases:
            if len(c) == 3:
                m, i, o = c
                queries.append((name, m, i))
            else:
                m, i, j, o = c
                queries.append((name, m, i, j))
            expect.append(o)
        got = run_crlibm.evaluate(queries)
        bad = 0
        for q, e, g in zip(queries, expect, got):
            total += 1
            if g != e:
                bad += 1
                bad_total += 1
                rows.append((name, q[1], " ".join(q[2:]), e, g))
        print("  %-8s vectors=%-6d disagree=%-4d" % (name, len(cases), bad))
    print("TOTAL vectors: %d   disagreements: %d" % (total, bad_total))
    with open("m1_results.tsv", "w") as fh:
        fh.write("func\tmode\tinput\texpected\tgot\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("wrote m1_results.tsv (%d rows)" % len(rows))


if __name__ == "__main__":
    main()
