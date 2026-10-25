"""All measurements behind APRIORI-BOUND.md.

Primary: at the paper's per-order worst cells, can the a priori signed
intervals order the two entries?  Both sweeps, the 67 committed comparisons
of the model's sweep, and every variant of the coefficient bound.
Secondary: band-maximum certification over the paper's population, envelope
tightness, residual signedness, and the rational route's lower bound.
"""

import csv
import json
import math
import sys
import time
from collections import Counter, defaultdict

from mpmath import mp, mpf, inf

import bounds

mp.dps = 60

THETA0 = bounds.THETA0

SWEEPS = {
    "amplification_form": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    "product_amplification": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
}


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------

def read_reference_args(root):
    args = defaultdict(list)
    with open(f"{root}/data/reference.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            args[row["band"]].append(float(row["x"]))
    for band in args:
        args[band].sort()
    return args


def read_per_order(root):
    rows = {}
    with open(f"{root}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            rows[(r["band"], int(r["n"]), r["method"])] = dict(
                worst_abs=mpf(r["worst_abs"]), at_x_abs=float(r["at_x_abs"]),
                worst_rel=mpf(r["worst_rel"]), at_x_rel=float(r["at_x_rel"]))
    return rows


def read_prediction(root):
    path = f"{root}/results/prediction.csv"
    try:
        handle = open(path, newline="")
    except FileNotFoundError:
        return None
    out = {}
    with handle:
        for r in csv.DictReader(handle):
            out[(r["band"], r["metric"], int(r["order"]))] = dict(
                state=r["state"], winner=r["measured_winner"])
    return out


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def stats(values):
    vals = sorted(float(v) for v in values if v != mp.inf)
    n = len(vals)
    if n == 0:
        return dict(n=0)
    def q(p):
        return vals[min(n - 1, max(0, int(round(p * (n - 1)))))]
    return dict(n=n, min=vals[0], median=q(0.5), q10=q(0.1), q90=q(0.9),
                max=vals[-1])


def certify_magnitude(lo_s, hi_s, lo_p, hi_p):
    return bounds.certify(lo_s, hi_s, lo_p, hi_p)


def certify_relative(band, n, x_s, x_p, eb_s, eb_p):
    ls, hs = bounds.relative_interval(band, n, x_s, eb_s, signed_method="coef")
    lp, hp = bounds.relative_interval(band, n, x_p, eb_p, signed_method="coef")
    if hs < lp:
        return "ship", (ls, hs), (lp, hp)
    if hp < ls:
        return "product", (ls, hs), (lp, hp)
    return "", (ls, hs), (lp, hp)


def trunc_envelope(band, n, method, x, kind="sharp"):
    """Symmetric truncation envelope of one entry at one argument (cheap)."""
    x = mpf(x)
    bi = 0 if band == "A1" else 1
    if method == "product":
        return bounds.fit_sym(n, bi, x, kind)
    if x < THETA0:
        transfer = bounds.amplification(n, x) / bounds.amplification(bounds.NMAX, x)
        return transfer * bounds.fit_sym(bounds.NMAX, 0, x, kind)
    if n == 0 or x >= mpf(bounds.TABLE["thresholds"][n]):
        return bounds.amplification(n, x) * bounds.seed_sym(x, kind)
    return bounds.fit_sym(n, bi, x, kind)


# ---------------------------------------------------------------------------
# the cell comparison, one variant
# ---------------------------------------------------------------------------

VARIANTS = {
    "primary_sampled400_paper": dict(signed_method="sampled", samples=400,
                                     size_kind="apriori", charge_kind="paper"),
    "sampled100_paper": dict(signed_method="sampled", samples=100,
                             size_kind="apriori", charge_kind="paper"),
    "sampled1600_paper": dict(signed_method="sampled", samples=1600,
                              size_kind="apriori", charge_kind="paper"),
    "sampled400_paper_strictsize": dict(signed_method="sampled", samples=400,
                                        size_kind="bound", charge_kind="paper"),
    "sampled400_worstcharge": dict(signed_method="sampled", samples=400,
                                   size_kind="apriori", charge_kind="worst"),
    "coef_paper": dict(signed_method="coef", samples=0,
                       size_kind="apriori", charge_kind="paper"),
    "crude_paper": dict(signed_method="crude", samples=0,
                        size_kind="apriori", charge_kind="paper"),
}


def cell_pair(band, metric, n, order, variant):
    ship = order[(band, n, "ship")]
    prod = order[(band, n, "product")]
    if metric == "absolute":
        x_s, x_p = ship["at_x_abs"], prod["at_x_abs"]
        del_s, del_p = ship["worst_abs"], prod["worst_abs"]
    else:
        x_s, x_p = ship["at_x_rel"], prod["at_x_rel"]
        del_s, del_p = ship["worst_rel"], prod["worst_rel"]
    kw = VARIANTS[variant]
    bs = bounds.entry_bounds(band, n, "ship", x_s, **kw)
    bp = bounds.entry_bounds(band, n, "product", x_p, **kw)
    if metric == "absolute":
        cert = certify_magnitude(bs["lo"], bs["hi"], bp["lo"], bp["hi"])
        rels = relp = None
    else:
        cert, rels, relp = certify_relative(band, n, x_s, x_p, bs, bp)
    return dict(band=band, metric=metric, n=n, x_ship=x_s, x_prod=x_p,
                delivered_ship=del_s, delivered_prod=del_p,
                bs=bs, bp=bp, cert=cert, rel_ship=rels, rel_prod=relp)


def run_sweep(name, root, do_variants=True):
    args = read_reference_args(root)
    order = read_per_order(root)
    pred = read_prediction(root)
    out = dict(name=name, population={b: len([x for x in xs if mpf(x) >= THETA0])
                                      for b, xs in args.items()},
               total_args={b: len(xs) for b, xs in args.items()})

    comparisons = [(band, metric, n)
                   for band in ("A1", "A2")
                   for metric in ("absolute", "relative")
                   for n in range(33)]
    primary = {}
    for band, metric, n in comparisons:
        primary[(band, metric, n)] = cell_pair(band, metric, n, order,
                                               "primary_sampled400_paper")

    # --- counts by variant over the 132 comparisons -------------------------
    variant_counts = {}
    for variant in VARIANTS:
        counts = Counter()
        directions = Counter()
        for band, metric, n in comparisons:
            row = primary[(band, metric, n)] if variant == "primary_sampled400_paper" \
                else cell_pair(band, metric, n, order, variant)
            cert = row["cert"]
            counts["certified" if cert else "abstain"] += 1
            if cert:
                directions[cert] += 1
        variant_counts[variant] = dict(counts=dict(counts),
                                       directions=dict(directions))

    # --- committed comparisons and agreement --------------------------------
    committed = [k for k in primary if pred and pred[k]["state"] in ("ship", "product")]
    cert_committed = [k for k in committed if primary[k]["cert"]]
    agree = [k for k in cert_committed if primary[k]["cert"] == pred[k]["state"]]
    cert_all = [k for k in primary if primary[k]["cert"]]
    cert_all_vs_delivered = [
        k for k in cert_all
        if primary[k]["cert"] == (
            "ship" if primary[k]["delivered_ship"] < primary[k]["delivered_prod"]
            else "product")]
    out["counts"] = dict(
        comparisons=len(comparisons),
        committed=len(committed),
        certified_among_committed=len(cert_committed),
        certified_among_committed_agree=len(agree),
        certified_at_all=len(cert_all),
        certified_at_all_agree_with_delivered=len(cert_all_vs_delivered),
        by_metric=dict(Counter(k[1] for k in committed)),
    )
    out["variant_counts"] = variant_counts
    out["committed_detail"] = [
        dict(band=k[0], metric=k[1], n=k[2], state=pred[k]["state"],
             measured_winner=pred[k]["winner"], cert=primary[k]["cert"],
             x_ship=primary[k]["x_ship"], x_prod=primary[k]["x_prod"])
        for k in committed]

    # --- entry modes at the cells -------------------------------------------
    modes = Counter()
    same_fit = 0
    both_signed = 0
    for k in primary:
        row = primary[k]
        modes[(row["bs"]["mode"], row["bp"]["mode"])] += 1
        if row["bs"]["mode"] == "same-fit":
            same_fit += 1
        if (bounds.magnitude_interval(row["bs"]["lo"], row["bs"]["hi"])[0] > 0
                and bounds.magnitude_interval(row["bp"]["lo"], row["bp"]["hi"])[0] > 0):
            both_signed += 1
    out["modes"] = {"%s/%s" % k: v for k, v in modes.items()}
    out["same_fit_comparisons"] = same_fit
    out["both_intervals_one_signed"] = both_signed
    out["same_fit_is_identical"] = all(
        primary[k]["bs"]["trunc"] == primary[k]["bp"]["trunc"]
        for k in primary if primary[k]["bs"]["mode"] == "same-fit")

    # --- envelope / delivered ratios at every cell --------------------------
    ratios = defaultdict(list)
    one_signed = Counter()
    one_signed_total = Counter()
    for k in primary:
        row = primary[k]
        metric = k[1]
        for meth in ("ship", "product"):
            b = row["bs"] if meth == "ship" else row["bp"]
            delivered = row["delivered_ship"] if meth == "ship" else row["delivered_prod"]
            if delivered == 0:
                continue
            total = max(abs(b["lo"]), abs(b["hi"]))
            ratios[("total_" + metric)].append(total / delivered)
            ratios[("trunc_" + metric)].append(b["trunc"] / delivered)
            charge = b["charge"]
            ratios[("symsharp_" + metric)].append(
                (b["sym_sharp"] + charge) / delivered)
            ratios[("symcrude_" + metric)].append(
                (b["sym_crude"] + charge) / delivered)
            key = (meth, metric, b["mode"])
            one_signed_total[key] += 1
            if bounds.magnitude_interval(b["lo"], b["hi"])[0] > 0:
                one_signed[key] += 1
    out["ratios"] = {k: stats(v) for k, v in ratios.items()}
    out["one_signed_fraction"] = {
        "%s/%s/%s" % k: [one_signed[k], one_signed_total[k]]
        for k in one_signed_total}

    # --- band-maximum certification, absolute metric ------------------------
    band_max = {}
    for band in ("A1", "A2"):
        pop = [x for x in args[band] if mpf(x) >= THETA0]
        for n in range(33):
            for meth in ("ship", "product"):
                peak = mpf(0)
                at = None
                for x in pop:
                    e = trunc_envelope(band, n, meth, x, "sharp")
                    if e > peak:
                        peak, at = e, x
                charge = bounds.charge_uniform(n, band, meth, "paper")
                band_max[(band, n, meth)] = dict(
                    trunc=peak, at=at, charge=charge, total=peak + charge)
    cert_bm = Counter()
    bm_rows = []
    for band, metric, n in comparisons:
        if metric != "absolute":
            continue
        row = primary[(band, metric, n)]
        Bs = band_max[(band, n, "ship")]["total"]
        Bp = band_max[(band, n, "product")]["total"]
        ls = bounds.magnitude_interval(row["bs"]["lo"], row["bs"]["hi"])[0]
        lp = bounds.magnitude_interval(row["bp"]["lo"], row["bp"]["hi"])[0]
        cert = ""
        if Bs < lp:
            cert = "ship"
        if Bp < ls:
            cert = "product"
        cert_bm[cert or "abstain"] += 1
        if cert:
            bm_rows.append(dict(band=band, n=n, cert=cert,
                                B_ship=mp.nstr(Bs, 4), B_prod=mp.nstr(Bp, 4),
                                lower_ship=mp.nstr(ls, 4), lower_prod=mp.nstr(lp, 4)))
    out["band_max"] = dict(counts=dict(cert_bm), rows=bm_rows,
                           peaks={("%s/%d/%s" % k): dict(total=mp.nstr(v["total"], 4),
                                                         at=v["at"])
                                  for k, v in band_max.items()})
    return out, primary


# ---------------------------------------------------------------------------
# residual signedness over every stored fit
# ---------------------------------------------------------------------------

def signedness():
    rows = []
    polys = [bounds.fit_poly(n, b) for n in range(33) for b in (0, 1)]
    polys.append(bounds.ext_poly())
    for p in polys:
        r = p.R
        spread = sum(abs(v) for v in r[1:]) + p.tail
        ratio = abs(r[0]) / spread if spread > 0 else mp.inf
        rows.append(dict(name=p.name, order=p.order, r0=mp.nstr(r[0], 6),
                         spread=mp.nstr(spread, 6), ratio=float(ratio)))
    ratios = sorted(x["ratio"] for x in rows)
    return dict(
        count=len(rows),
        one_signed=sum(1 for x in rows if x["ratio"] >= 1),
        max=rows[max(range(len(rows)), key=lambda i: rows[i]["ratio"])],
        median=ratios[len(ratios) // 2],
        detail=[x for x in rows if x["ratio"] > 0.1],
    )


# ---------------------------------------------------------------------------
# the rational route's |q| lower bound
# ---------------------------------------------------------------------------

def rational_bounds(samples=2000):
    table = bounds.TABLE
    rows = []
    for index in range(len(table["pieces"])):
        piece = table["pieces"][index]
        off = table["rat_offset"][index]
        m = table["rat_num_deg"][index]
        k = table["rat_den_deg"][index]
        c = [mpf(v) for v in table["rat_coeffs"][off:off + m + 1 + k]]
        p = c[:m + 1]
        q = [mpf(1)] + c[m + 1:m + 1 + k]  # q_0 = 1, q_1..q_k stored
        crude = 1 - sum(abs(v) for v in q[1:])
        dq = sum(j * abs(q[j]) for j in range(1, len(q)))
        # sampled min of |q| on [-1,1] with a rigorous slope correction
        qv = None
        if k > 0:
            step = 2.0 / samples
            vals = [abs(bounds.core.poly_value(q, mpf(-1 + step * i)))
                    for i in range(samples + 1)]
            qv = min(vals) - dq * step / 2
        else:
            qv = mpf(1)
        p_l1 = sum(abs(v) for v in p)
        rows.append(dict(index=index, num_deg=m, den_deg=k,
                         crude=float(crude), sampled=float(qv),
                         value_bound=(float(p_l1 / crude) if crude > 0 else None),
                         sampled_value_bound=(float(p_l1 / qv) if qv > 0 else None)))
    return dict(
        pieces=len(rows),
        crude_positive=sum(1 for r in rows if r["crude"] > 0),
        sampled_positive=sum(1 for r in rows if r["sampled"] > 0),
        crude_min=min(r["crude"] for r in rows),
        sampled_min=min(r["sampled"] for r in rows),
        worst=sorted(rows, key=lambda r: r["sampled"])[:5],
    )


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    report = {}
    primaries = {}
    for name, root in SWEEPS.items():
        print(f"[{name}] starting at {time.time() - t0:.1f}s", flush=True)
        out, primary = run_sweep(name, root)
        report[name] = out
        primaries[name] = primary
        print(f"[{name}] done at {time.time() - t0:.1f}s", flush=True)
        print(f"  population {out['population']} of {out['total_args']}")
        print("  counts:", out["counts"])
        for variant, vc in out["variant_counts"].items():
            print("   ", variant, vc)
        print("  modes:", out["modes"], "same-fit:", out["same_fit_comparisons"],
              "both one-signed:", out["both_intervals_one_signed"],
              "identical trunc at same-fit:", out["same_fit_is_identical"])
        for key in sorted(out["ratios"]):
            print("   ratio", key, {k: (round(v, 3) if isinstance(v, float) else v)
                                    for k, v in out["ratios"][key].items()})
        print("  band-max counts:", out["band_max"]["counts"])
    report["signedness"] = signedness()
    print("signedness:", {k: v for k, v in report["signedness"].items()
                          if k != "detail"})
    report["rational"] = rational_bounds()
    print("rational:", {k: v for k, v in report["rational"].items() if k != "worst"})
    with open(r"C:\Users\myaml\AppData\Local\Temp\apriori-bound\results.json", "w") as handle:
        json.dump(report, handle, indent=1, default=str)
    print("total %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
