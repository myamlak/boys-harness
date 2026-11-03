"""Which evaluator is inaccurate at the oracle outlier?

Cell: A1, x = 0.00040666245819265257 (the float key), n = 24, product.
Compares, at increasing working precision:
  - the fast route's truncation t = fit(x) - reference F_n(x),
  - the same fit against an independent hyp1f1 evaluation of F_n,
  - the deposit's identity oracle (exact.entry_truncation_exact).
"""
import sys

from mpmath import mp, mpf, hyp1f1

sys.path.insert(0, r"C:/Users/myaml/dev/_a-cert-run/code")
sys.path.insert(0, r"C:/Users/myaml/dev/_harness-v1023/amplification-sweeps/apriori-bound")
import exactroute as er
import exact

ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
X = float("0.00040666245819265257")
N = 24


def fn_hyp(x, n):
    return hyp1f1(n + mpf(1) / 2, n + mpf(3) / 2, -x) / (2 * n + 1)


def main():
    for dps in (70, 100, 140):
        mp.dps = dps
        x = mpf(X)
        ref = er.load_reference(ROOT)[("A1", X)][N]
        fit = er.fit_value(N, 0, x)
        hyp = fn_hyp(x, N)
        print("dps", dps)
        print("  fit-ref   ", mp.nstr(fit - ref, 12))
        print("  fit-hyp   ", mp.nstr(fit - hyp, 12))
        print("  ref-hyp   ", mp.nstr(ref - hyp, 12))
        t = exact.entry_truncation_exact("A1", N, "product", x)
        print("  oracle    ", mp.nstr(t, 20))


if __name__ == "__main__":
    main()
