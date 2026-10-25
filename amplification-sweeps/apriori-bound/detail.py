"""Per-comparison detail of the primary certificate and the ordering gap."""

import csv
import json
from collections import Counter

from mpmath import mp, mpf

import bounds
import run_all

mp.dps = 60


def magnitude(lo, hi):
    return bounds.magnitude_interval(lo, hi)


def pack(entry):
    return {k: (mp.nstr(v, 8) if isinstance(v, mpf) else v)
            for k, v in entry.items()}


def main():
    out = {}
    for name, root in run_all.SWEEPS.items():
        args = run_all.read_reference_args(root)
        order = run_all.read_per_order(root)
        pred = run_all.read_prediction(root)
        rows = {}
        for band in ("A1", "A2"):
            for metric in ("absolute", "relative"):
                for n in range(33):
                    row = run_all.cell_pair(band, metric, n, order,
                                            "primary_sampled400_paper")
                    bs, bp = row["bs"], row["bp"]
                    us = magnitude(bs["lo"], bs["hi"])
                    up = magnitude(bp["lo"], bp["hi"])
                    # width of the signed interval under the other two routes
                    widths = {}
                    for route in ("coef", "crude"):
                        wa = bounds.entry_bounds(band, n, "ship", row["x_ship"],
                                                 signed_method=route, samples=0)
                        wb = bounds.entry_bounds(band, n, "product", row["x_prod"],
                                                 signed_method=route, samples=0)
                        widths[route] = dict(
                            ship=str(wa["hi"] - wa["lo"]),
                            product=str(wb["hi"] - wb["lo"]))
                    widths["sampled"] = dict(
                        ship=str(bs["hi"] - bs["lo"]),
                        product=str(bp["hi"] - bp["lo"]))
                    winner = ("ship" if row["delivered_ship"] < row["delivered_prod"]
                              else "product" if row["delivered_prod"] < row["delivered_ship"]
                              else "tie")
                    if winner == "ship":
                        U, L = us[1], up[0]
                        delivered_winner = row["delivered_ship"]
                        delivered_loser = row["delivered_prod"]
                    else:
                        U, L = up[1], us[0]
                        delivered_winner = row["delivered_prod"]
                        delivered_loser = row["delivered_ship"]
                    p = pred[(band, metric, n)] if pred else dict(state="", winner="")
                    d = dict(band=band, metric=metric, n=n, x_ship=row["x_ship"],
                             x_prod=row["x_prod"],
                             delivered_ship=str(row["delivered_ship"]),
                             delivered_prod=str(row["delivered_prod"]),
                             state=p["state"], paper_winner=p["winner"],
                             delivered_winner=winner, cert=row["cert"],
                             ship=pack(bs), product=pack(bp), widths=widths,
                             U_winner=mp.nstr(U, 8), L_loser=mp.nstr(L, 8),
                             envelope_gap_ratio=(float(U / L) if L > 0 else None),
                             envelope_over_winner_error=float(U / delivered_winner)
                             if delivered_winner else None,
                             lower_over_loser_error=float(L / delivered_loser)
                             if delivered_loser else None)
                    rows[f"{band}|{metric}|{n}"] = d
        out[name] = rows
        committed = [k for k, r in rows.items() if r["state"] in ("ship", "product")]
        print(f"== {name}: committed {len(committed)}")
        ratios = [rows[k]["envelope_gap_ratio"] for k in committed
                  if rows[k]["envelope_gap_ratio"] is not None]
        zero = sum(1 for k in committed if rows[k]["envelope_gap_ratio"] is None)
        print("   committed with L_loser = 0 (sign obstruction):", zero)
        if ratios:
            print("   gap ratio among the rest: min %.3g median %.3g max %.3g"
                  % (min(ratios), sorted(ratios)[len(ratios) // 2], max(ratios)))
        env = [rows[k]["envelope_over_winner_error"] for k in committed]
        low = [rows[k]["lower_over_loser_error"] for k in committed]
        print("   U_winner / delivered_winner: min %.3g median %.3g max %.3g"
              % (min(env), sorted(env)[len(env) // 2], max(env)))
        print("   L_loser / delivered_loser: min %.3g median %.3g max %.3g"
              % (min(low), sorted(low)[len(low) // 2], max(low)))
        print("   certified among committed:",
              sum(1 for k in committed if rows[k]["cert"]))
        for k in committed:
            r = rows[k]
            print("    %s %-8s n=%2d %-7s winner=%-7s cert=%-7s x_s=%-10.4g "
                  "x_p=%-10.4g Uw=%s Ll=%s ratio=%s"
                  % (r["band"], r["metric"], r["n"], r["state"], r["delivered_winner"],
                     r["cert"] or "-", r["x_ship"], r["x_prod"], r["U_winner"],
                     r["L_loser"],
                     "%.3g" % r["envelope_gap_ratio"]
                     if r["envelope_gap_ratio"] is not None else "inf"))
    with open(r"C:\Users\myaml\AppData\Local\Temp\apriori-bound\detail.json", "w") as handle:
        json.dump(out, handle, indent=1)


if __name__ == "__main__":
    main()
