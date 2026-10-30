import sys
from pathlib import Path

paper = Path(r"C:\Users\myaml\dev\papers-a-v5\papers\boys-crossover")
run = paper / "sweeps" / "v5.2.0-check" / "run"

pairs = [
    ("product_amplification", "results/cells.csv"),
    ("product_amplification", "results/admission.csv"),
    ("product_amplification", "results/replica_gain.csv"),
    ("product_amplification", "results/per_order.csv"),
    ("product_amplification", "results/band_summary.csv"),
    ("product_amplification", "results/gain_up.csv"),
    ("product_amplification", "results/growth.csv"),
    ("product_amplification", "results/decomposition.csv"),
    ("product_amplification", "results/summary.txt"),
    ("product_amplification", "results/provenance.json"),
    ("product_amplification", "data/documented.csv"),
    ("amplification_form", "results/paths.csv"),
    ("amplification_form", "results/cells.csv"),
    ("amplification_form", "results/admission.csv"),
    ("amplification_form", "results/per_order.csv"),
    ("amplification_form", "results/band_summary.csv"),
    ("amplification_form", "results/prediction.csv"),
    ("amplification_form", "results/prediction_scan.csv"),
    ("amplification_form", "results/prediction_record.csv"),
    ("amplification_form", "results/transfer_cells.csv"),
    ("amplification_form", "results/summary.txt"),
    ("amplification_form", "results/provenance.json"),
]

def lines(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n").split(b"\n")

print("sweep tree / file: byte comparison, CRLF folded to LF")
for tree, relative in pairs:
    a = paper / "sweeps" / tree / relative
    b = run / "sweeps" / tree / relative
    if not b.is_file():
        print("MISSING NEW  %-22s %s" % (tree, relative))
        continue
    la, lb = lines(a), lines(b)
    if la == lb:
        print("SAME        %-22s %s  (%d lines, %d bytes deposited)" % (tree, relative, len(la) - 1, a.stat().st_size))
        continue
    print("DIFF        %-22s %s  (%d vs %d lines)" % (tree, relative, len(la) - 1, len(lb) - 1))
    shown = 0
    for i, (x, y) in enumerate(zip(la, lb), 1):
        if x != y:
            print("    line %d" % i)
            print("      deposited: %s" % x.decode(errors="replace")[:220])
            print("      new      : %s" % y.decode(errors="replace")[:220])
            shown += 1
            if shown >= 6:
                break
    if len(la) != len(lb):
        print("    length differs beyond the common prefix")
