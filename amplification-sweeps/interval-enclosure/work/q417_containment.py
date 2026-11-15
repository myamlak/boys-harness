"""Diagnostic: does the referee interval contain the measured recursed error?

For every served population cell, compare |measured recursed error| with the
model's per-cell enclosure [max(0, S-tau), S+tau].  Also compare each row's
measured recursed band maximum with the band enclosure [L, U].

Read-only against the harness; writes only under this record's work/
directory.  The repository root is derived from this file's location; set
BOYS_HARNESS to override it.
"""
import argparse
import csv
import json
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
WORK = Path(__file__).resolve().parent

sys.path.insert(0, str(WORK))
from q417_interval import CASES, read_cells, read_reference  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, choices=sorted(CASES))
    args = parser.parse_args()
    case = CASES[args.case]
    reference = read_reference(case["reference"])
    cells = read_cells(case["cells"])
    thresholds = [float(v) for v in json.load(
        open(case["provenance"], encoding="utf-8"))["constants"]["kTierThresholds"]]
    model_from = thresholds[0]

    cells_checked = 0
    mag_below_lo = 0
    mag_above_hi = 0
    signed_outside = 0
    worst_below = mpf(0)
    worst_above = mpf(0)
    # per row measured recursed max and band enclosure from the case CSV
    measured_max = {}
    for band in ("A1", "A2"):
        for (b, x), entry in cells.items():
            if b != band or x < model_from:
                continue
            F = reference[(band, x)]
            xm = mpf(x)
            eps0 = abs(mpf(entry[0][0]) - F[0])
            A = [mpf(1)] * (NMAX + 1)
            invA = [mpf(1)] * (NMAX + 1)
            twox = 2 * xm
            for k in range(1, NMAX + 1):
                A[k] = A[k - 1] * (2 * k - 1) / twox
                invA[k] = invA[k - 1] * twox / (2 * k - 1)
            ex = mp.exp(-xm)
            P = mpf(0)
            Q = mpf(0)
            for n in range(NMAX + 1):
                if n > 0:
                    P += F[n] * invA[n]
                    Q += invA[n - 1] / (2 * (n - 1) + 1)
                recursed = mpf(entry[n][0])
                e = recursed - F[n]
                mag = abs(e)
                for metric, scale in (("absolute", mpf(1)),
                                      ("relative", 1 / F[n])):
                    key = (band, metric, n)
                    current = measured_max.get(key)
                    if current is None or mag * scale > current:
                        measured_max[key] = mag * scale
                if (n == 0) or (x >= thresholds[n]):
                    tau = 3 * U53 * A[n] * P + ETA_E * ex * A[n] * Q
                    s = eps0 * A[n]
                    lo = s - tau
                    if lo < 0:
                        lo = mpf(0)
                    hi = s + tau
                    cells_checked += 1
                    if mag < lo:
                        mag_below_lo += 1
                        if lo != 0:
                            worst_below = max(worst_below, (lo - mag) / lo)
                    if mag > hi:
                        mag_above_hi += 1
                        if hi != 0:
                            worst_above = max(worst_above, (mag - hi) / hi)
                    if e < lo or e > hi:
                        signed_outside += 1

    rows = list(csv.DictReader(open(WORK / (args.case + ".csv"), newline="")))
    # L and U travel through the record at 17 significant digits; compare the
    # measured maximum at the same precision so an exact-equality row (the
    # order-0 interval is a point) is not split by representation round-off.
    measured17 = {key: mpf(mp.nstr(value, 17))
                  for key, value in measured_max.items()}
    row_below = row_above = 0
    worst_row_below = mpf(0)
    worst_row_above = mpf(0)
    below_list = []
    above_list = []
    for row in rows:
        key = (row["band"], row["metric"], int(row["order"]))
        mag = measured17[key]
        lo, hi = mpf(row["L"]), mpf(row["U"])
        if mag < lo:
            row_below += 1
            if lo != 0:
                worst_row_below = max(worst_row_below, (lo - mag) / lo)
            below_list.append(dict(band=key[0], metric=key[1], order=key[2],
                                   measured_max=mp.nstr(mag, 17),
                                   L=row["L"], U=row["U"],
                                   unserved_cells=int(row["unserved_cells"])))
        if mag > hi:
            row_above += 1
            if hi != 0:
                worst_row_above = max(worst_row_above, (mag - hi) / hi)
            above_list.append(dict(band=key[0], metric=key[1], order=key[2],
                                   measured_max=mp.nstr(mag, 17),
                                   L=row["L"], U=row["U"],
                                   unserved_cells=int(row["unserved_cells"])))

    # Do the enclosure rule's own commitments rest on a valid bound?  A
    # recursed commitment (U < Dpop) needs the measured recursed band maximum
    # at or below U; a direct commitment (Dpop < L) needs it at or above L.
    recursed_rows = []
    direct_rows = []
    for row in rows:
        key = (row["band"], row["metric"], int(row["order"]))
        mag = measured17[key]
        lo, hi = mpf(row["L"]), mpf(row["U"])
        if row["enc_state"] == "recursed":
            recursed_rows.append((key, mag <= hi, mag, lo, hi))
        elif row["enc_state"] == "direct":
            direct_rows.append((key, mag >= lo, mag, lo, hi))

    def summarise(entries, endpoint):
        failed = [e for e in entries if not e[1]]
        return dict(
            commits=len(entries),
            endpoint_holds=len(entries) - len(failed),
            endpoint_fails=len(failed),
            failing=[dict(band=k[0], metric=k[1], order=k[2],
                          measured_max=mp.nstr(e[2], 17),
                          L=mp.nstr(e[3], 17), U=mp.nstr(e[4], 17))
                     for e in failed
                     for k in [e[0]]],
            endpoint=endpoint,
        )

    out = dict(
        case=args.case,
        served_cells=cells_checked,
        magnitude_below_lower=mag_below_lo,
        magnitude_above_upper=mag_above_hi,
        signed_outside=signed_outside,
        worst_below_relative=mp.nstr(worst_below, 17),
        worst_above_relative=mp.nstr(worst_above, 17),
        rows=len(rows),
        measured_recursed_max_below_L=row_below,
        measured_recursed_max_above_U=row_above,
        rows_below=below_list,
        rows_above=above_list,
        worst_row_below_relative=mp.nstr(worst_row_below, 17),
        worst_row_above_relative=mp.nstr(worst_row_above, 17),
        enclosure_recursed_commitments=summarise(
            recursed_rows, "measured_max <= U"),
        enclosure_direct_commitments=summarise(
            direct_rows, "measured_max >= L"),
    )
    path = WORK / (args.case + ".containment.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(out, handle, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
