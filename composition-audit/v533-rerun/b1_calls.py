#!/usr/bin/env python3
"""Per joint call, the largest measured error over candidate per-axis sums.

Reduces the deposited v5.0.0 cube (../../records/cells) to one row per call and
writes b1_calls.tsv beside this script.  The four denominators are the accessor
figure (err/figure), the row's whole published statement (err/limit), six
copies of that statement (err/(6*limit)), and the same six over the lane base.
The exact sum over the six separate axis rows is not defined: the published
axis rows are heterogeneous, so this record quotes the ratios it can form and
does not claim a per-axis sum.
"""
import csv
import lzma
from pathlib import Path

HERE = Path(__file__).resolve().parent
FILES = {
    "fused": HERE.parent / "records" / "cells" / "cells-msvc-fused.csv.xz",
    "separate": HERE.parent / "records" / "cells" / "cells-msvc-separate.csv.xz",
}
FLOATISH = {"fp32", "fp16", "bf16"}
HALF = {"fp16", "bf16"}

rows_out = []
for cfg, path in FILES.items():
    best = {}
    with lzma.open(path, "rt", newline="") as f:
        for r in csv.DictReader(f):
            if r["status"] != "ok":
                continue
            lane = r["precision"]
            key = (lane, r["route"], r["scheme"], r["partition"], r["pack"],
                   r["form"], r["exp"])
            n = int(r["n"])
            x = float(r["x"])
            xnom = float(r["x_nominal"])
            err = float(r["abs_err"])
            fig = float(r["figure"])
            lim = float(r["limit"])
            formterm = 1e-7 if (lane in FLOATISH and r["form"] == "plain-reciprocal") else 0.0
            base = fig - formterm
            statement = lim if lane in HALF else fig
            dens = {
                "fig_ratio": fig,
                "statement_ratio": statement,
                "six_statement_ratio": 6.0 * statement,
                "six_base_ratio": 6.0 * base,
                "base_ratio": base,
            }
            b = best.get(key)
            if b is None:
                b = best[key] = {"max_abs_err": -1.0, "cell": None}
                for d in dens:
                    b[d] = -1.0
            if err > b["max_abs_err"]:
                b["max_abs_err"] = err
                b["cell"] = (n, x)
            for d, den in dens.items():
                if den > 0:
                    q = err / den
                    if q > b[d]:
                        b[d] = q
    for key, b in sorted(best.items()):
        n, x = b["cell"]
        rows_out.append([cfg] + list(key) + [repr(b["max_abs_err"]), n, repr(x)]
                        + [repr(b[d]) for d in ("fig_ratio", "statement_ratio",
                                                "six_statement_ratio", "six_base_ratio",
                                                "base_ratio")])

header = ["config", "lane", "route", "scheme", "partition", "pack", "form", "exp",
          "max_abs_err", "cell_n", "cell_x", "fig_ratio", "statement_ratio",
          "six_statement_ratio", "six_base_ratio", "base_ratio"]
with open(HERE / "b1_calls.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(header)
    w.writerows(rows_out)

print("rows", len(rows_out))
for cfg in ("fused", "separate"):
    print("===", cfg)
    for lane in ("fp64", "fp32", "fp16", "bf16"):
        sub = [r for r in rows_out if r[0] == cfg and r[1] == lane]
        print(" ", lane,
              "fig", max(float(r[11]) for r in sub),
              "stmt", max(float(r[12]) for r in sub),
              "sixstmt", max(float(r[13]) for r in sub),
              "sixbase", max(float(r[14]) for r in sub),
              "base", max(float(r[15]) for r in sub))
