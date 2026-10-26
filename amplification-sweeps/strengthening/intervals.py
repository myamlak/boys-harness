"""Intervals for the manuscript's headline counts, with the dependence stated.

    python intervals.py > intervals.log 2>&1

Imports the deposited surrogate control and reproduces its margin-1 and
margin-0.94 numbers before computing anything.  The control names the paper
tree it was first run from in its own source; the import block below binds it
to the deposit's analyze.py and data files, which are byte-identical to that
tree's.  The 132 band-metric comparisons are four groups over the same 33
orders and the same argument sets, so the comparisons are not independent; a
binomial interval treats them as independent trials.  Both are reported, plus
an order-block jackknife and bootstrap that resample whole orders.
"""

import csv
import importlib.util
import math
import random
import sys
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 70

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE.parent

# surrogate_control.py names the paper tree it was first run from in a
# module-level TREE.  A clone of this deposit has no paper tree, and every
# record the control reads is here, so register the deposit's analyze.py under
# the name it imports, then point its TREE at the deposit.  The data files are
# the same bytes the paper tree held, so no computed number changes.
_spec = importlib.util.spec_from_file_location(
    "analyze", DEPOSIT / "amplification-form" / "analyze.py"
)
_analyze = importlib.util.module_from_spec(_spec)
sys.path.insert(0, str(DEPOSIT / "common"))
sys.path.insert(0, str(DEPOSIT / "amplification-form"))
sys.modules["analyze"] = _analyze
_spec.loader.exec_module(_analyze)

sys.path.insert(0, str(DEPOSIT / "surrogate-control"))
import surrogate_control as S  # noqa: E402

S.TREE = DEPOSIT / "amplification-form"
if S.A is not _analyze:
    raise SystemExit("surrogate_control did not import the deposit's analyze.py")

NMAX = S.NMAX
ORDERS = list(range(NMAX + 1))
BANDS = ("A1", "A2")
METRICS = ("absolute", "relative")


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    denominator = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return (max(0.0, centre - half), min(1.0, centre + half))


def beta_regularized(a, b, x):
    return mp.betainc(a, b, 0, x, regularized=True)


def beta_quantile(a, b, p):
    """Inverse regularized incomplete beta by bisection."""
    lo, hi = mpf(0), mpf(1)
    for _ in range(200):
        mid = (lo + hi) / 2
        value = beta_regularized(a, b, mid)
        if value < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clopper_pearson(k, n, alpha=0.05):
    if k == 0:
        lower = mpf(0)
    else:
        lower = beta_quantile(mpf(k), mpf(n - k + 1), mpf(alpha / 2))
    if k == n:
        upper = mpf(1)
    else:
        upper = beta_quantile(mpf(k + 1), mpf(n - k), mpf(1 - alpha / 2))
    return (float(lower), float(upper))


def frozen_min_factory():
    cache = {}

    def nmin_of(x):
        return min(NMAX, max(0, math.ceil(x - 0.5)))

    def frozen_min(n, x):
        m = nmin_of(x)
        if n <= m:
            return S.A.amplification(n, x)
        key = (m, x)
        if key not in cache:
            cache[key] = S.A.amplification(m, x)
        return cache[key]

    return frozen_min


def table_rows(rec_abs, prod_abs, arith_abs, rec_rel, prod_rel, arith_rel, winner, margin):
    return S.all_rows(rec_abs, prod_abs, arith_abs, rec_rel, prod_rel, arith_rel, winner, margin)


def flatten(table):
    """(band, metric, order) -> row."""
    out = {}
    for (band, metric), rows in table.items():
        for row in rows:
            out[(band, metric, row["n"])] = row
    return out


def jackknife(values):
    """Delete-one estimate list -> (mean, standard error)."""
    m = len(values)
    mean = sum(values) / m
    variance = (m - 1) / m * sum((v - mean) ** 2 for v in values)
    return mean, math.sqrt(variance)


def block_ratio_bootstrap(by_order, numerator, denominator, reps=20000, seed=20261026):
    """Resample whole orders; numerator and denominator are sums of the blocks."""
    rng = random.Random(seed)
    num = {o: numerator(o) for o in ORDERS}
    den = {o: denominator(o) for o in ORDERS}
    ratios = []
    m = len(ORDERS)
    for _ in range(reps):
        total_num = total_den = 0
        for _ in range(m):
            o = rng.choice(ORDERS)
            total_num += num[o]
            total_den += den[o]
        ratios.append(total_num / total_den if total_den else float("nan"))
    ratios.sort()
    return ratios


def percentile(sorted_values, q):
    index = q * (len(sorted_values) - 1)
    lo = int(math.floor(index))
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = index - lo
    return sorted_values[lo] * (1 - frac) + sorted_values[hi] * frac


def main():
    reference, cells, error0, winner, *_ = S.build_measured()
    entries, population, model_from = S.build_cells(reference, cells, error0)
    print("population:", population, "model_from:", model_from)

    true_abs = S.predicted(entries, S.A.amplification, "absolute")
    true_rel = S.predicted(entries, S.A.amplification, "relative")
    frozen = frozen_min_factory()
    sur_abs = S.predicted(entries, frozen, "absolute")
    sur_rel = S.predicted(entries, frozen, "relative")

    # ---- validation against the deposited control ------------------------
    true1 = flatten(table_rows(*true_abs, *true_rel, winner, mpf(1)))
    sur1 = flatten(table_rows(*sur_abs, *sur_rel, winner, mpf(1)))
    true094 = flatten(table_rows(*true_abs, *true_rel, winner, mpf("0.94")))
    sur094 = flatten(table_rows(*sur_abs, *sur_rel, winner, mpf("0.94")))
    true094b = flatten(table_rows(*true_abs, *true_rel, winner, mpf("0.946")))
    sur094b = flatten(table_rows(*sur_abs, *sur_rel, winner, mpf("0.946")))
    true0 = flatten(table_rows(*true_abs, *true_rel, winner, mpf(0)))
    sur0 = flatten(table_rows(*sur_abs, *sur_rel, winner, mpf(0)))

    def score_map(rows):
        committed = [k for k, r in rows.items() if r["state"] != "abstain"]
        right = sum(1 for k in committed if rows[k]["right"])
        return len(committed), right, len(committed) - right

    print("validation (deposited control numbers in brackets):")
    print("  true margin 1:", score_map(true1), "[67, 67, 0]")
    print("  surrogate margin 1:", score_map(sur1), "[21, 21, 0]")
    print("  true margin 0.94:", score_map(true094), "[74, 74, 0]")
    print("  surrogate margin 0.94:", score_map(sur094), "[68, 68, 0]")
    print("  true margin 0.946:", score_map(true094b), "[73, 73, 0]")
    print("  surrogate margin 0.946:", score_map(sur094b), "[67, 67, 0]")
    t0 = score_map(true0)
    print("  true margin 0:", t0, "[100, 83, 17]")

    true_set1 = {k for k, r in true1.items() if r["state"] != "abstain"}
    sur_set1 = {k for k, r in sur1.items() if r["state"] != "abstain"}
    print("  surrogate 21 inside true 67:", len(sur_set1 - true_set1) == 0,
          "; lost =", len(true_set1 - sur_set1))
    true_set094 = {k for k, r in true094.items() if r["state"] != "abstain"}
    sur_set094 = {k for k, r in sur094.items() if r["state"] != "abstain"}
    print("  74 set, 68 set, overlap:", len(true_set094), len(sur_set094),
          len(true_set094 & sur_set094), "; true only", len(true_set094 - sur_set094),
          "; surrogate only", len(sur_set094 - true_set094))
    print()

    # ---- intervals -------------------------------------------------------
    first_rows = list(csv.DictReader(open(HERE / "out" / "first" / "prediction.csv")))
    first_commits = [r for r in first_rows if r["state"] != "abstain"]
    first_right = sum(1 for r in first_commits if r["right_where_committed"] == "1")
    first_wrong = len(first_commits) - first_right
    print("first experiment, read from out/first/prediction.csv: commits %d, right %d, wrong %d"
          % (len(first_commits), first_right, first_wrong))
    counts = (
        ("model commitments right, second experiment, margin 1", 67, 67),
        ("wrong commitments, second experiment, margin 0", 17, 100),
        ("  of which not a measured tie", 14, 100),
        ("true-model commitments lost by the surrogate, margin 1", 46, 67),
        ("model commitments right, first experiment, margin 1",
         first_right, len(first_commits)),
    )
    print("binomial intervals (the independence assumption is false; the block rows below")
    print("resample whole orders):")
    print("count                                   k/n      Wilson 95%          Clopper-Pearson 95%")
    for label, k, n in counts:
        w = wilson(k, n)
        cp = clopper_pearson(k, n)
        print("  %-38s %3d/%-4d [%.4f, %.4f]   [%.4f, %.4f]"
              % (label, k, n, w[0], w[1], cp[0], cp[1]))
    print()

    # ---- per-order structure --------------------------------------------
    def by_order(rows, predicate):
        out = {o: 0 for o in ORDERS}
        for (band, metric, n), row in rows.items():
            if predicate(row):
                out[n] += 1
        return out

    wrong0 = by_order(true0, lambda r: r["state"] != "abstain" and not r["right"])
    commits0 = by_order(true0, lambda r: r["state"] != "abstain")
    print("margin 0: wrong comparisons by order (committed comparisons in brackets):")
    print("  " + "  ".join("%d:%d(%d)" % (o, wrong0[o], commits0[o]) for o in ORDERS
                           if wrong0[o] or commits0[o]))
    print("  distinct orders carrying at least one wrong commitment:",
          sum(1 for o in ORDERS if wrong0[o] > 0), "of", sum(1 for o in ORDERS if commits0[o] > 0),
          "orders with commitments")
    print()

    # order-block jackknife for the margin-0 wrong rate
    rates = []
    for o in ORDERS:
        num = sum(wrong0[x] for x in ORDERS if x != o)
        den = sum(commits0[x] for x in ORDERS if x != o)
        if den:
            rates.append(num / den)
    mean, se = jackknife(rates)
    print("margin 0 wrong rate, order-block jackknife: mean %.4f  SE %.4f  (naive binomial SE %.4f)"
          % (mean, se, math.sqrt(0.17 * 0.83 / 100)))
    boot = block_ratio_bootstrap(
        None,
        lambda o: wrong0[o],
        lambda o: commits0[o],
        reps=5000,
    )
    print("  order-block bootstrap 95%% interval for the wrong rate: [%.4f, %.4f]"
          % (percentile(boot, 0.025), percentile(boot, 0.975)))
    print()

    # the 46 lost commitments, by order and group
    lost = true_set1 - sur_set1
    lost_by_order = {o: 0 for o in ORDERS}
    for (band, metric, n) in lost:
        lost_by_order[n] += 1
    print("the 46 lost commitments by order (all right for the true model):")
    print("  " + "  ".join("%d:%d" % (o, lost_by_order[o]) for o in ORDERS if lost_by_order[o]))
    print("  orders with at least one lost commitment:", sum(1 for o in ORDERS if lost_by_order[o]))
    true_commit1 = by_order(true1, lambda r: r["state"] != "abstain")
    lost_boot = block_ratio_bootstrap(
        None,
        lambda o: lost_by_order[o],
        lambda o: true_commit1[o],
        reps=5000,
    )
    print("  order-block bootstrap 95%% interval for the lost fraction: [%.4f, %.4f]"
          % (percentile(lost_boot, 0.025), percentile(lost_boot, 0.975)))
    right1 = by_order(true1, lambda r: r["state"] != "abstain" and r["right"])
    right_boot = block_ratio_bootstrap(
        None,
        lambda o: right1[o],
        lambda o: true_commit1[o],
        reps=5000,
    )
    print("  the 67/67 right rate under the same block resample: [%.4f, %.4f]"
          " (every committed comparison is right, so resampling orders cannot"
          " produce a failure; the lower bound is a statement about the sampling"
          " of comparisons, not about this sample)"
          % (percentile(right_boot, 0.025), percentile(right_boot, 0.975)))
    print()

    # the six-decision gap: true 74 vs surrogate 68 at margin 0.94
    true_commit094 = by_order(true094, lambda r: r["state"] != "abstain")
    sur_commit094 = by_order(sur094, lambda r: r["state"] != "abstain")
    diff_blocks = {o: true_commit094[o] - sur_commit094[o] for o in ORDERS}
    print("margin 0.94 per-order commitment difference (true minus surrogate):")
    print("  " + "  ".join("%d:%+d" % (o, diff_blocks[o]) for o in ORDERS if diff_blocks[o]))
    print("  total difference:", sum(diff_blocks.values()))
    # order-block jackknife of the difference
    diffs = [sum(diff_blocks[x] for x in ORDERS if x != o) * len(ORDERS) / (len(ORDERS) - 1)
             for o in ORDERS]
    dmean, dse = jackknife(diffs)
    print("  order-block jackknife of the difference: mean %.3f SE %.3f" % (dmean, dse))
    rng = random.Random(20261026)
    boots = []
    for _ in range(5000):
        total = 0
        for _ in range(len(ORDERS)):
            total += diff_blocks[rng.choice(ORDERS)]
        boots.append(total)
    boots.sort()
    print("  order-block bootstrap 95%% interval for the difference: [%d, %d]"
          % (percentile(boots, 0.025), percentile(boots, 0.975)))

    # the null that the six is a threshold coincidence: at the same margin the
    # true model commits at 74 and the surrogate at 68; both are error-free.


if __name__ == "__main__":
    main()
