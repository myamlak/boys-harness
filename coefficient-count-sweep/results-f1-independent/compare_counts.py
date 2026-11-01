#!/usr/bin/env python3
"""Compare independent rational counts with the kernel's delivered record.

Usage:
  python compare_counts.py --kind f1 --independent out_f1_x0/run_f1_x0.json \
      out_f1_a2/run_f1_a2.json
  python compare_counts.py --kind tail --independent out_tail_x0/run_tail_x0.json \
      out_tail_a2/run_tail_a2.json
  python compare_counts.py --kind f0 --independent out_f0_a62/run_f0_a62.json

The kernel record is the corrected F1 record for F1, the pure-tail record for
the tail, and the corrected F0 record for F0. Prints one row per cell with
both counts, the pair, and the difference (independent - kernel).
"""
import argparse
import json
from pathlib import Path

SWEEP = Path(__file__).resolve().parent.parent
RECORDS = {
    "f1": SWEEP / "results-second-instance-f1" / "knee-f1" / "corrected" / "sweep_all.json",
    "tail": SWEEP / "results-second-instance" / "sweep_all.json",
    "f0": SWEEP / "results-corrected" / "corrected-counts.json",
}
TAG = {"rel13": "rel:1e-13", "abs514": "abs:5e-14", "rel14": "rel:1e-14"}
EDGE_TAG = {"x0": "x0", "a2": "a2"}
NAME = {"2^-7": "m7", "2^-6": "m6", "2^-5": "m5", "2^-4": "m4",
        "2^-3": "m3", "2^-2": "m2", "2^-1": "m1", "2^0": "p0",
        "2^1": "p1", "2^2": "p2", "2^3": "p3", "2^4": "p4",
        "2^5": "p5", "ship": "ship", "paper": "paper"}


def load_kernel(kind):
    doc = json.loads(RECORDS[kind].read_text())
    rows = doc.get("rows", doc.get("cells", []))
    out = {}
    for r in rows:
        if kind == "f0":
            edge = r.get("edge", "x0")
            target = r.get("target")
            name = r.get("width_name")
            rat = r.get("corrected_count")
            pair = r.get("pair")
        else:
            edge = r.get("edge")
            target = r.get("target")
            name = NAME.get(r.get("width_name"), r.get("width_name"))
            rat = (r.get("rat") or {}).get("stored") if isinstance(r.get("rat"), dict) else None
            pair = None
            if isinstance(r.get("rat"), dict):
                m, k = r["rat"].get("m"), r["rat"].get("k")
                pair = [m, k] if m is not None else None
        out[(edge, target, name)] = {"kernel": rat, "kernel_pair": pair}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=("f0", "f1", "tail"), required=True)
    ap.add_argument("--independent", nargs="+", required=True)
    args = ap.parse_args()
    kernel = load_kernel(args.kind)
    independent = {}
    for path in args.independent:
        doc = json.loads(Path(path).read_text())
        edge = doc.get("left_tag")
        for r in doc.get("cells", []):
            target = TAG.get(r.get("target_tag"), r.get("target_tag"))
            independent[(edge, target, r.get("width_name"))] = r
    print(f"{'edge':>4} {'target':>10} {'width':>6} {'kernel':>7} {'indep':>6} "
          f"{'diff':>5} {'pair':>10}")
    diffs = []
    for key in sorted(independent):
        ind = independent[key]
        ker = kernel.get(key, {})
        k = ker.get("kernel")
        i = ind.get("corrected_count")
        d = None if k is None or i is None else i - k
        if d is not None:
            diffs.append(d)
        print(f"{key[0]:>4} {key[1]:>10} {key[2]:>6} {str(k):>7} {str(i):>6} "
              f"{str(d):>5} {str(ind.get('pair')):>10}")
    if diffs:
        print(f"\ncells {len(diffs)}, differing {sum(1 for d in diffs if d)}, "
              f"sum {sum(diffs)}, max drop {min(diffs)}, max rise {max(diffs)}")


if __name__ == "__main__":
    main()
