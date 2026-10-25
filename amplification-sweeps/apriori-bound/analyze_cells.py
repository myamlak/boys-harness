"""The 132 band-metric-order comparisons, certified or not."""

import csv
from mpmath import mp, mpf

import core
import measure

mp.dps = 60

ROOT = measure.ROOT
PER_ORDER = f"{ROOT}/results/per_order.csv"


def read_prediction():
    pred = {}
    with open(f"{ROOT}/results/prediction.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            pred[(row["band"], row["metric"], int(row["order"]))] = dict(
                state=row["state"], measured=row["measured_winner"],
                ps=mpf(row["predicted_ship"]), pd=mpf(row["predicted_product"]),
                arith=mpf(row["arithmetic"]))
    return pred


def key_rows():
    """(band, metric, n) -> dict(ship=row, product=row)."""
    out = {}
    for row in measure.read_per_order(PER_ORDER):
        for metric, xkey, wkey in (("absolute", "at_x_abs", "worst_abs"),
                                   ("relative", "at_x_rel", "worst_rel")):
            k = (row["band"], metric, row["n"])
            out.setdefault(k, {})[row["method"]] = dict(
                x=row[xkey], worst=row[wkey])
    return out


def certify(lo_s, hi_s, lo_p, hi_p):
    ms = measure.magnitude_interval(lo_s, hi_s)
    mp_ = measure.magnitude_interval(lo_p, hi_p)
    if ms[1] < mp_[0]:
        return "ship"
    if mp_[1] < ms[0]:
        return "product"
    return ""


def main(samples=300, bound_kind="sharp", signed_method="sampled"):
    rows = key_rows()
    pred = read_prediction()
    out = []
    for (band, metric, n), pair in sorted(rows.items()):
        ship, prod = pair["ship"], pair["product"]
        bs = measure.entry_bounds(band, n, "ship", ship["x"], bound_kind, signed_method,
                                  samples=samples)
        bp = measure.entry_bounds(band, n, "product", prod["x"], bound_kind, signed_method,
                                  samples=samples)
        if bs is None:
            cert = "out-of-population"
        else:
            cert = certify(bs["lo"], bs["hi"], bp["lo"], bp["hi"])
        delivered = "ship" if ship["worst"] < prod["worst"] else (
            "product" if prod["worst"] < ship["worst"] else "tie")
        p = pred[(band, metric, n)]
        out.append(dict(
            band=band, metric=metric, n=n,
            x_ship=ship["x"], deliv_ship=ship["worst"],
            half_ship=(bs["half"] if bs else mpf(0)),
            lo_ship=(bs["lo"] if bs else mpf(0)), hi_ship=(bs["hi"] if bs else mpf(0)),
            mode_ship=(bs["mode"] if bs else "downward"),
            arith_ship=(bs["arithmetic"] if bs else mpf(0)),
            x_prod=prod["x"], deliv_prod=prod["worst"],
            half_prod=bp["half"], lo_prod=bp["lo"], hi_prod=bp["hi"],
            arith_prod=bp["arithmetic"],
            certified=cert, delivered=delivered, paper_state=p["state"],
        ))
    return out


def summarize(out):
    n = len(out)
    committed = [r for r in out if r["paper_state"] in ("ship", "product")]
    print("comparisons: %d, paper commits: %d" % (n, len(committed)))
    for r in committed:
        print("%s %-8s n=%2d paper=%-7s delivered=%-7s certified=%-16s "
              "x_s=%-9.4g x_p=%-9.4g deliv_s=%.3e deliv_p=%.3e half_s=%.3e half_p=%.3e"
              % (r["band"], r["metric"], r["n"], r["paper_state"], r["delivered"],
                 r["certified"] or "-", r["x_ship"], r["x_prod"],
                 float(r["deliv_ship"]), float(r["deliv_prod"]),
                 float(r["half_ship"]), float(r["half_prod"])))
    from collections import Counter
    print(Counter(r["certified"] for r in out))
    both = [r for r in out if r["certified"] == r["paper_state"] and r["certified"]]
    print("certified in the paper's committed direction:", len(both))


if __name__ == "__main__":
    import sys
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    out = main(samples)
    summarize(out)
