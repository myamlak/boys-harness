"""Which evaluator is inaccurate at the oracle outlier?

Cell: A1, x = 0.00040666245819265257 (the float key), n = 24, product.
Compares, at increasing working precision:
  - the fast route's truncation t = fit(x) - reference F_n(x),
  - the same fit against an independent hyp1f1 evaluation of F_n,
  - the deposit's identity oracle (exact.entry_truncation_exact).
"""
import argparse
import os
import sys
from pathlib import Path

from mpmath import mp, mpf, hyp1f1

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
import exactroute as er
import model
import record_paths

X = float("0.00040666245819265257")
N = 24


def fn_hyp(x, n):
    return hyp1f1(n + mpf(1) / 2, n + mpf(3) / 2, -x) / (2 * n + 1)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    record_paths.add_kernel(parser)
    record_paths.add_inputs(parser, deposit_arg=True)
    args = parser.parse_args(argv)
    model.configure(kernel=args.kernel, exp_bits=args.exp_bits)
    root = record_paths.sweep(args.sweep)
    dep = record_paths.deposit(args.deposit)
    with record_paths.deposit_on_path(dep):
        import exact
    mp.dps = 70

    ref = er.load_reference(root)
    for dps in (70, 100, 140):
        mp.dps = dps
        x = mpf(X)
        refv = ref[("A1", X)][N]
        fit = er.fit_value(N, 0, x)
        hyp = fn_hyp(x, N)
        print("dps", dps)
        print("  fit-ref   ", mp.nstr(fit - refv, 12))
        print("  fit-hyp   ", mp.nstr(fit - hyp, 12))
        print("  ref-hyp   ", mp.nstr(refv - hyp, 12))
        t = exact.entry_truncation_exact("A1", N, "product", x)
        print("  oracle    ", mp.nstr(t, 20))
    return 0


if __name__ == "__main__":
    sys.exit(main())
