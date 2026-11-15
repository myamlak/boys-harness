"""Q-417 interval analysis for Paper A, one experiment per invocation.

Reads only deposited records under the repository root that contains this
record (read-only).  Writes only under this record's work/ directory.  The
repository root is derived from this file's location; set BOYS_HARNESS to
override it.

Model (the paper's revision, deposited at
amplification-sweeps/product-amplification/analysis/sweeps/amplification_form/
analyze.py, the same functions copied in held-out/analyze.py):

    S(n,x)   = eps0(x) * A_n(x),   eps0(x) = |recursed_0(x) - F_0(x)|
    tau(n,x) = 3*2^-53 * sum_{l=1..n} F_l(x) A_n(x)/A_l(x)
             + 2^-52 * exp(-x) * sum_{l=0..n-1} A_n(x)/((2l+1) A_l(x))
    W = max rec,  D = max direct error,  T = max tau    (population x >= theta0)

The model is served at order n when n == 0 or x >= kTierThresholds[n]; at an
unserved order inside the population both entries read the same stored fit,
and the deposited model sets the predicted recursed error to the measured
direct error (arithmetic term zero).  The enclosure reproduces that as the
degenerate interval [d(x), d(x)]; see README.md beside this record.

Referee Q-417 enclosure, same population and same metric units:

    L = max_x max(0, S - tau),  U = max_x (S + tau),  Dpop = max_x d
    direct   if Dpop < L
    recursed if U < Dpop
    abstain  otherwise

Since the deposited W is S+tau at served cells and d at unserved cells, W is
exactly U; the invalid step is W - T as a lower bound.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
from mpmath import mp, mpf

mp.dps = 70
NMAX = 32
U53 = mpf(2) ** -53
ETA_E = mpf(2) ** -52
HARNESS = Path(os.environ.get(
    "BOYS_HARNESS", Path(__file__).resolve().parents[3]))
WORK = Path(__file__).resolve().parent.parent

CASES = {
    "first9353": dict(
        cells=HARNESS / "boys-crossover/exact-tag-replay/run/cells-9353.csv",
        reference=HARNESS / "amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/data/reference.csv",
        provenance=HARNESS / "amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/results/provenance.json",
        per_order=HARNESS / "amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/results/per_order.csv",
        per_order_kind="method",
        record=None,
    ),
    "first9352": dict(
        cells=HARNESS / "boys-crossover/exact-tag-replay/run/cells-9352.csv",
        reference=HARNESS / "amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/data/reference.csv",
        provenance=HARNESS / "amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/results/provenance.json",
        per_order=HARNESS / "amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/results/per_order.csv",
        per_order_kind="method",
        record=None,
    ),
    "second": dict(
        cells=HARNESS / "amplification-sweeps/amplification-form/results/cells.csv",
        reference=HARNESS / "amplification-sweeps/amplification-form/data/reference.csv",
        provenance=HARNESS / "amplification-sweeps/amplification-form/results/provenance.json",
        per_order=HARNESS / "amplification-sweeps/amplification-form/results/per_order.csv",
        per_order_kind="method",
        record=HARNESS / "amplification-sweeps/certified-bound/inputs/model-record/prediction_record.csv",
        record_labels=("predicted_ship", "predicted_product", "arithmetic", "ship", "product"),
    ),
    "heldout": dict(
        cells=HARNESS / "amplification-sweeps/held-out/results/cells.csv",
        reference=HARNESS / "amplification-sweeps/held-out/data/reference.csv",
        provenance=HARNESS / "amplification-sweeps/held-out/results/provenance.json",
        per_order=HARNESS / "amplification-sweeps/held-out/results/per_order.csv",
        per_order_kind="band_metric",
        record=HARNESS / "amplification-sweeps/held-out/results/prediction_heldout.csv",
        record_labels=("predicted_recursed", "predicted_direct", "arithmetic", "recursed", "direct"),
    ),
}

PAPER_ANALYZE = (HARNESS / "amplification-sweeps/product-amplification/analysis/"
                 "sweeps/amplification_form/analyze.py")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_float(text):
    return float(text.strip().replace("nan(ind)", "nan"))


def read_reference(path):
    reference = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header[:2] != ["band", "x"]:
            raise SystemExit("unexpected reference header: %r" % (header[:2],))
        for row in reader:
            values = [mpf(value) for value in row[2:]]
            if len(values) != NMAX + 1:
                raise SystemExit("reference row has %d values" % len(values))
            reference[(row[0], float(row[1]))] = values
    return reference


def read_cells(path):
    data = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            band, x, n = row["band"], float(row["x"]), int(row["n"])
            entry = data.get((band, x))
            if entry is None:
                entry = [None] * (NMAX + 1)
                data[(band, x)] = entry
            entry[n] = (csv_float(row["recursed"]), csv_float(row["direct"]))
    return data


def update(store, key, value, x):
    current = store.get(key)
    if current is None or value > current[0]:
        store[key] = (value, x)


def winner_of(ship, prod):
    if ship < prod:
        return "recursed"
    if prod < ship:
        return "direct"
    return "tie"


def nstr(value):
    if value is None:
        return ""
    return mp.nstr(value, 17)


def per_order_check(path, kind, ship_abs, prod_abs, ship_rel, prod_rel):
    """Compare recomputed whole-band order maxima with the deposited file."""
    if not path.is_file():
        return {"file": str(path), "present": False}
    checked = mismatched = 0
    worst = mpf(0)
    examples = []
    if kind == "method":
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                band, n = row["band"], int(row["n"])
                method = row["method"]
                store_abs = ship_abs if method == "ship" else prod_abs
                store_rel = ship_rel if method == "ship" else prod_rel
                for metric, store, col, at_col in (
                    ("absolute", store_abs, "worst_abs", "at_x_abs"),
                    ("relative", store_rel, "worst_rel", "at_x_rel"),
                ):
                    value, x = store[(band, metric, n)]
                    dep_value = mpf(row[col])
                    dep_x = float(row[at_col])
                    checked += 1
                    rel = abs(value - dep_value) / abs(dep_value) if dep_value != 0 else abs(value)
                    worst = max(worst, rel)
                    if rel > mpf("1e-15") or x != dep_x:
                        mismatched += 1
                        if len(examples) < 5:
                            examples.append((band, metric, n, method,
                                             nstr(value), row[col],
                                             repr(x), row[at_col]))
    else:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                band, n, metric = row["band"], int(row["order"]), row["metric"]
                for store, col, at_col in (
                    (ship_abs if metric == "absolute" else ship_rel,
                     "measured_recursed_worst", "measured_recursed_at_x"),
                    (prod_abs if metric == "absolute" else prod_rel,
                     "measured_direct_worst", "measured_direct_at_x"),
                ):
                    rec_key = (band, metric, n)
                    value, x = store[rec_key]
                    dep_value = mpf(row[col])
                    dep_x = float(row[at_col])
                    checked += 1
                    rel = abs(value - dep_value) / abs(dep_value) if dep_value != 0 else abs(value)
                    worst = max(worst, rel)
                    if rel > mpf("1e-15") or x != dep_x:
                        mismatched += 1
                        if len(examples) < 5:
                            examples.append((band, metric, n, col,
                                             nstr(value), row[col],
                                             repr(x), row[at_col]))
    return {"file": str(path), "present": True, "checked": checked,
            "mismatched": mismatched, "max_rel": nstr(worst),
            "examples": examples}


def formula_self_check(cells, reference, thresholds):
    """Compare the cumulative-sum tau and A_n with the paper's own functions."""
    import importlib.util
    common = HARNESS / "amplification-sweeps/common"
    if str(common) not in sys.path:
        sys.path.insert(0, str(common))
    spec = importlib.util.spec_from_file_location("paper_a_analyze_q417", PAPER_ANALYZE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    args = {}
    for (band, x), rows in cells.items():
        if x >= thresholds[0]:
            args.setdefault(band, []).append(x)
    sample = []
    for band, xs in args.items():
        xs.sort()
        if not xs:
            continue
        for index in (0, len(xs) // 2, len(xs) - 1):
            x = xs[index]
            for n in range(NMAX + 1):
                if n == 0 or x >= thresholds[n]:
                    sample.append((band, x, n))
    worst_tau = mpf(0)
    worst_amp = mpf(0)
    count = 0
    for band, x, n in sample:
        values = reference[(band, x)]
        A = nstr(module.amplification(n, mpf(x)))
        # the recurrence form used by the main loop
        amp = mpf(1)
        twox = 2 * mpf(x)
        for k in range(1, n + 1):
            amp = amp * (2 * k - 1) / twox
        dep_amp = module.amplification(n, mpf(x))
        rel = abs(amp - dep_amp) / abs(dep_amp) if dep_amp != 0 else abs(amp)
        worst_amp = max(worst_amp, rel)
        tau_paper = module.recurrence_arithmetic(n, mpf(x), values)
        # cumulative-sum evaluation, the main loop's arithmetic
        invA = [mpf(1)] * (NMAX + 1)
        for k in range(1, NMAX + 1):
            invA[k] = invA[k - 1] * twox / (2 * k - 1)
        P = mpf(0)
        Q = mpf(0)
        for k in range(1, n + 1):
            P += values[k] * invA[k]
            Q += invA[k - 1] / (2 * (k - 1) + 1)
        tau = 3 * U53 * amp * P + ETA_E * mp.exp(-mpf(x)) * amp * Q
        rel = abs(tau - tau_paper) / abs(tau_paper) if tau_paper != 0 else abs(tau)
        worst_tau = max(worst_tau, rel)
        count += 1
    return {"samples": count, "max_rel_tau": nstr(worst_tau),
            "max_rel_amplification": nstr(worst_amp),
            "paper_analyze": str(PAPER_ANALYZE),
            "paper_analyze_sha256": sha256(PAPER_ANALYZE)}


def record_check(path, labels, rows, winner_whole):
    """Compare recomputed W, D, T, state and winners with a deposited record."""
    if path is None or not Path(path).is_file():
        return {"file": str(path), "present": False}
    col_w, col_d, col_t, lab_rec, lab_dir = labels
    checked = 0
    state_mismatch = winner_mismatch = choice_mismatch = right_mismatch = 0
    max_rel = {col_w: mpf(0), col_d: mpf(0), col_t: mpf(0)}
    examples = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["band"], row["metric"], int(row["order"]))
            mine = rows[key]
            checked += 1
            for col, value in ((col_w, mine["W"]), (col_d, mine["D"]), (col_t, mine["T"])):
                dep = mpf(row[col])
                rel = abs(value - dep) / abs(dep) if dep != 0 else abs(value)
                max_rel[col] = max(max_rel[col], rel)
                if rel > mpf("1e-15"):
                    if len(examples) < 5:
                        examples.append((key, col, nstr(value), row[col]))
            dep_state = row["state"]
            if dep_state in (lab_rec, lab_dir):
                mapped = {"recursed": "recursed", "direct": "direct",
                          "ship": "recursed", "product": "direct"}[dep_state]
            else:
                mapped = "abstain"
            if mine["dep_state"] != mapped:
                state_mismatch += 1
            dep_winner = row["measured_winner"]
            mapped_winner = {"recursed": "recursed", "direct": "direct",
                             "ship": "recursed", "product": "direct",
                             "tie": "tie"}.get(dep_winner, dep_winner)
            if mine["winner_whole"] != mapped_winner:
                winner_mismatch += 1
            dep_choice = row["choice"]
            mapped_choice = {"recursed": "recursed", "direct": "direct",
                             "ship": "recursed", "product": "direct",
                             "": "abstain", "abstain": "abstain"}.get(dep_choice, dep_choice)
            if mine["dep_state"] != mapped_choice:
                choice_mismatch += 1
            dep_right = row["right_where_committed"].strip()
            mine_right = mine["dep_right_whole"]
            expected_right = "1" if mine_right == "1" else "0"
            if dep_right != expected_right:
                right_mismatch += 1
    return {"file": str(path), "present": True, "checked": checked,
            "state_mismatches": state_mismatch,
            "winner_mismatches": winner_mismatch,
            "choice_mismatches": choice_mismatch,
            "right_mismatches": right_mismatch,
            "max_rel": {col: nstr(value) for col, value in max_rel.items()},
            "examples": examples}


def run_case(name):
    case = CASES[name]
    reference = read_reference(case["reference"])
    cells = read_cells(case["cells"])
    constants = json.load(open(case["provenance"], encoding="utf-8"))["constants"]
    thresholds = [float(v) for v in constants["kTierThresholds"]]
    if len(thresholds) <= NMAX:
        raise SystemExit("threshold table too short")
    model_from = thresholds[0]

    bands = ("A1", "A2")
    metrics = ("absolute", "relative")
    orders = range(NMAX + 1)

    # measured maxima, whole band and population, per (band, metric, n)
    ship_abs, ship_rel, prod_abs, prod_rel = {}, {}, {}, {}
    pop_ship_abs, pop_ship_rel, pop_prod_abs, pop_prod_rel = {}, {}, {}, {}

    # model and enclosure accumulators
    W = {(b, m, n): mpf(0) for b in bands for m in metrics for n in orders}
    D = {(b, m, n): mpf(0) for b in bands for m in metrics for n in orders}
    T = {(b, m, n): mpf(0) for b in bands for m in metrics for n in orders}
    L = {(b, m, n): mpf(0) for b in bands for m in metrics for n in orders}
    U = {(b, m, n): mpf(0) for b in bands for m in metrics for n in orders}
    Dpop = {(b, m, n): mpf(0) for b in bands for m in metrics for n in orders}
    pop_args = {b: 0 for b in bands}
    pop_cells = {(b, n): 0 for b in bands for n in orders}
    served_cells = {(b, n): 0 for b in bands for n in orders}

    missing = []

    for band in bands:
        xs = sorted(x for (b, x) in cells if b == band)
        for x in xs:
            entry = cells[(band, x)]
            if any(item is None for item in entry):
                missing.append((band, x))
                continue
            F = reference[(band, x)]
            xm = mpf(x)

            # measured error maxima, whole band
            for n in orders:
                recursed, direct = entry[n]
                expected = F[n]
                rec_err = abs(mpf(recursed) - expected)
                dir_err = abs(mpf(direct) - expected)
                rec_rel = rec_err / expected if expected != 0 else mpf("inf")
                dir_rel = dir_err / expected if expected != 0 else mpf("inf")
                update(ship_abs, (band, "absolute", n), rec_err, x)
                update(prod_abs, (band, "absolute", n), dir_err, x)
                update(ship_rel, (band, "relative", n), rec_rel, x)
                update(prod_rel, (band, "relative", n), dir_rel, x)

            if x < model_from:
                continue
            pop_args[band] += 1
            eps0 = abs(mpf(entry[0][0]) - F[0])

            # A_n and 1/A_l by recurrence; cumulative sums for tau
            A = [mpf(1)] * (NMAX + 1)
            invA = [mpf(1)] * (NMAX + 1)
            twox = 2 * xm
            for k in range(1, NMAX + 1):
                A[k] = A[k - 1] * (2 * k - 1) / twox
                invA[k] = invA[k - 1] * twox / (2 * k - 1)
            ex = mp.exp(-xm)
            P = mpf(0)
            Q = mpf(0)
            for n in orders:
                if n > 0:
                    P += F[n] * invA[n]
                    Q += invA[n - 1] / (2 * (n - 1) + 1)
                expected = F[n]
                recursed, direct = entry[n]
                dir_err = abs(mpf(direct) - expected)
                dir_rel = dir_err / expected if expected != 0 else mpf("inf")
                pop_cells[(band, n)] += 1
                update(pop_ship_abs, (band, "absolute", n),
                       abs(mpf(recursed) - expected), x)
                update(pop_prod_abs, (band, "absolute", n), dir_err, x)
                update(pop_ship_rel, (band, "relative", n),
                       abs(mpf(recursed) - expected) / expected if expected != 0 else mpf("inf"), x)
                update(pop_prod_rel, (band, "relative", n), dir_rel, x)

                served = (n == 0) or (x >= thresholds[n])
                if served:
                    served_cells[(band, n)] += 1
                    tau_abs = 3 * U53 * A[n] * P + ETA_E * ex * A[n] * Q
                    s_abs = eps0 * A[n]
                    lo_abs = s_abs - tau_abs
                    if lo_abs < 0:
                        lo_abs = mpf(0)
                    hi_abs = s_abs + tau_abs
                    if expected != 0:
                        lo_rel = lo_abs / expected
                        hi_rel = hi_abs / expected
                        tau_rel = tau_abs / expected
                    else:
                        lo_rel = hi_rel = tau_rel = mpf("inf")
                else:
                    tau_abs = mpf(0)
                    tau_rel = mpf(0)
                    lo_abs = hi_abs = dir_err
                    lo_rel = hi_rel = dir_rel

                for metric, lo, hi, dval, tauval in (
                    ("absolute", lo_abs, hi_abs, dir_err, tau_abs),
                    ("relative", lo_rel, hi_rel, dir_rel, tau_rel),
                ):
                    key = (band, metric, n)
                    if hi > U[key]:
                        U[key] = hi
                    if lo > L[key]:
                        L[key] = lo
                    if dval > Dpop[key]:
                        Dpop[key] = dval
                    if hi > W[key]:
                        W[key] = hi
                    if dval > D[key]:
                        D[key] = dval
                    if tauval > T[key]:
                        T[key] = tauval

    # winners
    winner_whole, winner_pop = {}, {}
    for band in bands:
        for metric in metrics:
            for n in orders:
                key = (band, metric, n)
                store_ship = ship_abs if metric == "absolute" else ship_rel
                store_prod = prod_abs if metric == "absolute" else prod_rel
                winner_whole[key] = winner_of(store_ship[key][0], store_prod[key][0])
                p_ship = (pop_ship_abs if metric == "absolute" else pop_ship_rel)[key]
                p_prod = (pop_prod_abs if metric == "absolute" else pop_prod_rel)[key]
                winner_pop[key] = winner_of(p_ship[0], p_prod[0])

    rows = {}
    for band in bands:
        for metric in metrics:
            for n in orders:
                key = (band, metric, n)
                gap = abs(W[key] - D[key])
                if gap <= T[key]:
                    dep_state = "abstain"
                else:
                    dep_state = "recursed" if W[key] < D[key] else "direct"
                if Dpop[key] < L[key]:
                    enc_state = "direct"
                elif U[key] < Dpop[key]:
                    enc_state = "recursed"
                else:
                    enc_state = "abstain"

                def right(state, winner):
                    if state == "abstain":
                        return "n/a"
                    if winner == "tie":
                        return "tie"
                    return "1" if state == winner else "0"

                rows[key] = dict(
                    experiment=name, band=band, metric=metric, order=n,
                    population_arguments=pop_args[band],
                    population_cells=pop_cells[(band, n)],
                    served_cells=served_cells[(band, n)],
                    unserved_cells=pop_cells[(band, n)] - served_cells[(band, n)],
                    W=W[key], D=D[key], T=T[key],
                    W_minus_T=W[key] - T[key],
                    L=L[key], U=U[key], Dpop=Dpop[key],
                    dep_state=dep_state,
                    dep_choice=dep_state if dep_state != "abstain" else "",
                    enc_state=enc_state,
                    winner_whole=winner_whole[key],
                    winner_pop=winner_pop[key],
                    dep_right_whole=right(dep_state, winner_whole[key]),
                    dep_right_pop=right(dep_state, winner_pop[key]),
                    enc_right_whole=right(enc_state, winner_whole[key]),
                    enc_right_pop=right(enc_state, winner_pop[key]),
                    dep_justified_by_enc=(
                        "n/a" if dep_state == "abstain"
                        else ("1" if enc_state == dep_state else "0")),
                    enc_commits_dep_abstains=(
                        "1" if dep_state == "abstain" and enc_state != "abstain" else "0"),
                )

    def count(rows, state_field, right_field):
        out = {"commits": 0, "abstains": 0, "recursed": 0, "direct": 0,
               "right": 0, "wrong": 0, "tie": 0}
        for row in rows:
            state = row[state_field]
            if state == "abstain":
                out["abstains"] += 1
                continue
            out["commits"] += 1
            out[state] += 1
            r = row[right_field]
            if r == "1":
                out["right"] += 1
            elif r == "0":
                out["wrong"] += 1
            else:
                out["tie"] += 1
        return out

    rowlist = [rows[(band, metric, n)] for band in bands
               for metric in metrics for n in orders]

    dep_direct_unjust = [(r["band"], r["metric"], r["order"])
                         for r in rowlist
                         if r["dep_state"] == "direct" and r["dep_justified_by_enc"] == "0"]
    dep_rec_unjust = [(r["band"], r["metric"], r["order"])
                      for r in rowlist
                      if r["dep_state"] == "recursed" and r["dep_justified_by_enc"] == "0"]
    enc_extra = [(r["band"], r["metric"], r["order"], r["enc_state"])
                 for r in rowlist if r["enc_commits_dep_abstains"] == "1"]
    w_below_l = [(r["band"], r["metric"], r["order"],
                  nstr(r["W_minus_T"]), nstr(r["L"]))
                 for r in rowlist if r["W_minus_T"] < r["L"]]
    u_ne_w = max(abs(r["U"] - r["W"]) for r in rowlist)

    summary = dict(
        case=name,
        inputs={key: {"path": str(case[key]), "sha256": sha256(case[key])}
                for key in ("cells", "reference", "provenance", "per_order")
                if case.get(key)},
        thresholds=[repr(v) for v in thresholds],
        model_from=repr(model_from),
        missing_cells=len(missing),
        rows=len(rowlist),
        counts=dict(
            deposited_whole=count(rowlist, "dep_state", "dep_right_whole"),
            deposited_pop=count(rowlist, "dep_state", "dep_right_pop"),
            enclosure_whole=count(rowlist, "enc_state", "enc_right_whole"),
            enclosure_pop=count(rowlist, "enc_state", "enc_right_pop"),
        ),
        state_agreement=sum(1 for r in rowlist if r["dep_state"] == r["enc_state"]),
        dep_direct_unjustified=dep_direct_unjust,
        dep_recursed_unjustified=dep_rec_unjust,
        enc_commits_where_dep_abstains=enc_extra,
        W_minus_T_below_L=w_below_l,
        max_abs_U_minus_W=nstr(u_ne_w),
        per_order=per_order_check(case["per_order"], case["per_order_kind"],
                                  ship_abs, prod_abs, ship_rel, prod_rel),
        formula_self_check=formula_self_check(cells, reference, thresholds),
        record=record_check(case.get("record"), case.get("record_labels"),
                            rows, winner_whole),
    )
    return summary, rowlist


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, choices=sorted(CASES))
    args = parser.parse_args()
    summary, rowlist = run_case(args.case)
    out = WORK / "work" / (args.case + ".json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    csv_path = WORK / "work" / (args.case + ".csv")
    columns = ["experiment", "band", "metric", "order",
               "population_arguments", "population_cells", "served_cells",
               "unserved_cells", "W", "D", "T", "W_minus_T", "L", "U", "Dpop",
               "dep_state", "dep_choice", "enc_state", "winner_whole",
               "winner_pop", "dep_right_whole", "dep_right_pop",
               "enc_right_whole", "enc_right_pop", "dep_justified_by_enc",
               "enc_commits_dep_abstains"]
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row in rowlist:
            writer.writerow([
                row[column] if not isinstance(row[column], mpf) else nstr(row[column])
                for column in columns])
    print(json.dumps({key: summary[key] for key in
                      ("case", "rows", "counts", "state_agreement",
                       "dep_direct_unjustified", "dep_recursed_unjustified",
                       "enc_commits_where_dep_abstains", "W_minus_T_below_L",
                       "max_abs_U_minus_W", "per_order", "formula_self_check",
                       "record")}, indent=1, default=str))
    print("wrote", out)
    print("wrote", csv_path)


if __name__ == "__main__":
    main()
