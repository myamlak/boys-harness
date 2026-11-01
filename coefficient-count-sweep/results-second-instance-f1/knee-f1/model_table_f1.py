#!/usr/bin/env python3
"""Model secants and the class comparison for the F1 instance.

Computes, with the paper's own formulas:

  * the fixed-constant model secants over the thirteen scaling widths at both
    edges (Chebyshev at Cc = 0.02, 0.03 and 1; the rational leading law at
    Cr = 4; the diagonal-odd and off-diagonal-even stored models), using the
    same functions the tail report used, imported from analyze_tail_instance;
  * the measured slope minus each model secant, in standard errors;
  * the F1-minus-F0 slope difference at x0 with an unpaired standard error,
    the class comparison (the rate should not depend on alpha);
  * the rational delivered/target margin per row, from the record, to show
    which cells meet their target only just.

Writes _knee-f1/model_table.json and prints the tables.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PAPER))

import mpmath as mp  # noqa: E402

import analyze_second_instance_f1 as af  # noqa: E402
import analyze_tail_instance as at  # noqa: E402

RECORD = PAPER / "results-second-instance-f1" / "sweep_all.json"
ANALYSIS = PAPER / "results-second-instance-f1" / "analysis.json"
OUT = PAPER / "_knee-f1" / "model_table.json"


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--analysis", default=str(ANALYSIS))
    ap.add_argument("--record", default=str(RECORD))
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    doc = json.loads(Path(args.analysis).read_text(encoding="utf-8"))
    slopes = doc["slopes"]
    f0 = doc["f0_bands"]["slopes"]

    mp.mp.dps = 50
    a_x0 = af.X0
    a_a2 = af.X0 * 2
    f_x0 = mp.mpf(doc["F1_at_x0"])
    f_a2 = mp.mpf(doc["F1_at_a2"])

    out = {"source_analysis": str(Path(args.analysis).resolve()),
           "model": {}, "class": {},
           "f0_slopes": f0, "margins": {}}

    for target in af.TARGETS:
        out["model"][target] = {}
        for fam in ("cheb", "rat"):
            out["model"][target][fam] = {}
            for edge, a, fval in (("x0", a_x0, f_x0), ("a2", a_a2, f_a2)):
                L = af.target_L(target, a, fval)
                entry = {"L": float(L)}
                if fam == "cheb":
                    for label, c in (("Cc=0.02", "0.02"),
                                     ("Cc=0.03", "0.03"),
                                     ("Cc=1", "1")):
                        entry[label] = at.model_secant("cheb", a, L, mp.mpf(c))
                else:
                    entry["Cr=4"] = at.model_secant("rat", a, L, mp.mpf(4))
                    entry["diag_odd"] = at.model_secant_rat_stored(
                        a, L, "diag_odd")
                    entry["offdiag_even"] = at.model_secant_rat_stored(
                        a, L, "offdiag_even")
                measured = slopes[target][fam][edge]
                entry["measured"] = measured["beta"]
                entry["se"] = measured["se"]
                entry["n"] = measured["n"]
                entry["delta"] = {
                    k: {"delta": measured["beta"] - v,
                        "se": (measured["beta"] - v) / measured["se"]}
                    for k, v in entry.items()
                    if k.startswith(("Cc=", "Cr=")) or k in
                    ("diag_odd", "offdiag_even")
                }
                out["model"][target][fam][edge] = entry

    for target in af.TARGETS:
        out["class"][target] = {}
        for fam in ("cheb", "rat"):
            f1 = slopes[target][fam]["x0"]
            f0x = f0[target][fam]["x0"]
            diff = f1["beta"] - f0x["beta"]
            se = (f1["se"] ** 2 + f0x["se"] ** 2) ** 0.5
            out["class"][target][fam] = {
                "f1_x0": f1["beta"], "f1_se": f1["se"], "f1_n": f1["n"],
                "f0_x0": f0x["beta"], "f0_se": f0x["se"], "f0_n": f0x["n"],
                "diff": diff, "se_unpaired": se, "z": diff / se,
            }

    record = json.loads(Path(args.record).read_text(encoding="utf-8"))
    rows = []
    for row in record["rows"]:
        rat = row["rat"]
        kind, value = row["target"].split(":", 1)
        eps = float(value)
        if kind == "rel":
            delivered = float(rat["delivered_rel"])
        else:
            delivered = float(rat["delivered_abs"])
        rows.append({
            "edge": row["edge"], "width": row["width_name"],
            "target": row["target"], "stored": rat["stored"],
            "delivered": delivered, "eps": eps,
            "margin": delivered / eps,
            "below_meets": rat.get("below_meets"),
        })
    out["margins"] = rows

    Path(args.out).write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"wrote {args.out}")

    print("\nfixed-constant model secants (measured +/- s.e.; delta/se)")
    for target in af.TARGETS:
        for fam in ("cheb", "rat"):
            for edge in ("x0", "a2"):
                e = out["model"][target][fam][edge]
                keys = [k for k in e if k.startswith(("Cc=", "Cr=")) or
                        k in ("diag_odd", "offdiag_even")]
                cells = "; ".join(
                    f"{k} {e[k]:.5f} ({e['delta'][k]['delta']:+.4f}, "
                    f"{e['delta'][k]['se']:+.2f} s.e.)" for k in keys)
                print(f"{target:10s} {fam:4s} {edge}: "
                      f"{e['measured']:.5f}+/-{e['se']:.5f} n={e['n']} | "
                      f"{cells}")
    print("\nF1 - F0 at x0 (class comparison)")
    for target in af.TARGETS:
        for fam in ("cheb", "rat"):
            c = out["class"][target][fam]
            print(f"{target:10s} {fam:4s}: F1 {c['f1_x0']:.5f}+/-"
                  f"{c['f1_se']:.5f} (n={c['f1_n']}) minus F0 "
                  f"{c['f0_x0']:.5f}+/-{c['f0_se']:.5f} (n={c['f0_n']}) = "
                  f"{c['diff']:+.5f} +/- {c['se_unpaired']:.5f} "
                  f"({c['z']:+.2f} s.e.)")
    print("\nrational margins delivered/target")
    for edge in ("x0", "a2"):
        for target in af.TARGETS:
            vals = [r["margin"] for r in rows
                    if r["edge"] == edge and r["target"] == target]
            vals_sorted = sorted(vals)
            print(f"{edge} {target:10s}: n={len(vals)} min "
                  f"{vals_sorted[0]:.3f} median {vals_sorted[len(vals)//2]:.3f} "
                  f"max {vals_sorted[-1]:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
