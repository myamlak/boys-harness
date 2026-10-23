#!/usr/bin/env python3
"""Print the six headline numbers of the paper from the deposited records.

Run from this directory, or point it at a copy of the audit material:

    python headline_numbers.py
    python headline_numbers.py --root /path/to/boys-harness/accuracy-claim-audit

The first five numbers come from the two claim registries. Both registries
are in the boys-harness repository at tag v1.0.0, so those five print from
that repository alone. The sixth needs the captured cell summary,
records/summarize_cells.txt, which is supplied with the paper and is not in
the repository; when it is absent the script says so rather than guessing.

Only the Python standard library is used.
"""

import argparse
import csv
import sys
from pathlib import Path


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def boost_counts(path):
    """Total documented statements, and the strictly row-shaped count."""
    rows = read_tsv(path)
    total = sum(int(row["count"]) for row in rows)
    row_shaped = sum(
        int(row["count"])
        for row in rows
        if row["row_shaped"].strip().lower() == "yes"
    )
    return total, row_shaped


def crlibm_counts(path):
    """Documented claims, and the strictly row-shaped count.

    A documented claim is a tier-A or tier-B entry for a named function. The
    tier-B entry for "all functions" is a global bound over the library, not
    a claim about one entry point, and is excluded by the paper's counting
    rule. A strictly row-shaped entry is class ROW.
    """
    rows = read_tsv(path)
    claims = [
        row
        for row in rows
        if row["tier"] in ("A", "B")
        and row["function"].strip() != "all functions"
    ]
    strict = sum(1 for row in rows if row["class"].strip() == "ROW")
    return len(claims), strict


def boost_false_claims(cells_path, rows_path):
    """Boost claims whose registered row has a judged worst value over bound."""
    bounds = {}
    with rows_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            bounds[row["id"]] = (row["function"], float(row["bound"]))

    failed_functions = set()
    with cells_path.open(encoding="utf-8") as handle:
        next(handle)  # header
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            row_id = fields[0]
            worst = float(fields[4])
            if row_id not in bounds:
                continue
            function, bound = bounds[row_id]
            if worst > bound:
                failed_functions.add(function)
    return len(failed_functions)


def crlibm_false_claims(path):
    """Strictly row-shaped entries whose note records a measured falsehood."""
    rows = read_tsv(path)
    return sum(
        1
        for row in rows
        if row["class"].strip() == "ROW"
        and "MEASURED FALSE" in (row["note"] or "").upper()
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory holding the audit material (default: this directory)",
    )
    args = parser.parse_args()
    root = args.root

    boost_tsv = root / "boost_accuracy_rows" / "claims" / "claims_audit.tsv"
    rows_csv = root / "boost_accuracy_rows" / "claims" / "rows.csv"
    crlibm_tsv = root / "claim-audit-crlibm" / "work" / "claims_registry.tsv"
    cells_txt = root / "records" / "summarize_cells.txt"

    missing = [path for path in (boost_tsv, rows_csv, crlibm_tsv) if not path.is_file()]
    if missing:
        for path in missing:
            print(f"missing: {path}", file=sys.stderr)
        return 2

    boost_total, boost_strict = boost_counts(boost_tsv)
    crlibm_claims, crlibm_strict = crlibm_counts(crlibm_tsv)
    together = boost_strict + crlibm_strict

    if cells_txt.is_file():
        false_total = boost_false_claims(cells_txt, rows_csv) + crlibm_false_claims(
            crlibm_tsv
        )
        false_note = ""
    else:
        false_total = None
        false_note = (
            "  [unavailable: records/summarize_cells.txt is not in the "
            "boys-harness repository]"
        )

    print(f"{boost_total:>3}  documented accuracy statements in Boost.Math")
    print(f"{boost_strict:>3}  strictly row-shaped among them")
    print(f"{crlibm_claims:>3}  documented accuracy claims in crlibm")
    print(f"{crlibm_strict:>3}  strictly row-shaped among them")
    print(f"{together:>3}  strictly row-shaped across the two")
    if false_total is None:
        print(f"  ?  row-shaped claims measured false{false_note}")
    else:
        print(f"{false_total:>3}  row-shaped claims measured false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
