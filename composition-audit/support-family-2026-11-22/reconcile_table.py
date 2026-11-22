"""Per-term breakdown of the old and corrected frontier rows."""

from __future__ import annotations

import json

FIG = {"fp64": 2.558314305293831e-14, "fp32": 1.5e-7}
OLD_SUM = {"path": 14.0, "pooled": 66.0}
KJ = {"path": 4.0, "pooled": 3.0}


def main():
    d = json.load(open("corrected_frontier.json"))
    t = d["table"]
    rows = []
    for key in sorted(t):
        v = t[key]
        lane = key.split("/")[1]
        kind = key.split("/")[2]
        fig = FIG[lane]
        full_nonadd = v["full_nonadd"]
        full = fig + full_nonadd
        eta = v["eta_Q_max"]
        old_charge = (1.0 + OLD_SUM[kind]) * eta
        old_family = fig + old_charge
        old_inc = old_family - full
        new_charge = KJ[kind] * eta
        new_family = fig + v["posterior"]
        new_inc = new_family - full
        rows.append({
            "row": key,
            "figure": fig,
            "full_nonadd": full_nonadd,
            "full_bound": full,
            "eta_Q_max": eta,
            "old_sum_Gk": OLD_SUM[kind],
            "old_charge": old_charge,
            "old_family_bound": old_family,
            "old_increment": old_inc,
            "old_increment_in_eta": old_inc / eta,
            "kappa_joint": KJ[kind],
            "new_charge": new_charge,
            "M_hat_max": v["M_hat_max"],
            "posterior": v["posterior"],
            "new_family_bound": new_family,
            "new_increment": new_inc,
            "new_increment_in_eta": new_inc / eta,
            "new_ratio": new_family / full,
            "at_post": v["at_post"],
        })
    print("row | full_nonadd | full_bound | eta_Q | old_sumGk | old_charge | "
          "old_family | old_inc | old_inc/eta | kj | new_charge | M_hat | "
          "posterior | new_family | new_inc | new_inc/eta | new_ratio")
    for r in rows:
        print(f"{r['row']} | {r['full_nonadd']:.6e} | {r['full_bound']:.6e} | "
              f"{r['eta_Q_max']:.6e} | {r['old_sum_Gk']:.0f} | "
              f"{r['old_charge']:.6e} | {r['old_family_bound']:.6e} | "
              f"{r['old_increment']:.6e} | {r['old_increment_in_eta']:.3f} | "
              f"{r['kappa_joint']:.0f} | {r['new_charge']:.6e} | "
              f"{r['M_hat_max']:.6e} | {r['posterior']:.6e} | "
              f"{r['new_family_bound']:.6e} | {r['new_increment']:.6e} | "
              f"{r['new_increment_in_eta']:.3f} | {r['new_ratio']:.6f}")
    with open("reconcile_table.json", "w") as h:
        json.dump(rows, h, indent=1)


if __name__ == "__main__":
    main()
