"""Reproduce and diagnose the three secondary checks in
amplification-form/results/summary.txt that the manuscript does not explain.

    python check_transfers.py > check_transfers.log 2>&1

Reads only the decompressed copies of the deposited records under inputs/ and
the deposited amplification-form reference and transfer_cells records.  The
sampling, the gates and the arithmetic are those of the deposited
amplification-form/analyze.py:

  (1) exact-perturbation transfer: rows at row_index % 97 == 0 with n != start,
      measured transfer_ulp against A_n(x)/A_start(x);
  (2) start-independence: starts not in (0, 32) at every 16th argument,
      (value_n - ref_n) / (value_start - ref_start) against A_n(x)/A_start(x);
  (3) the two-parameter n^p q^n fit to log(transfer) over the measurable
      canonical cells with n >= 1.
"""

import csv
import json
import math
from pathlib import Path

from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
DEPOSIT = HERE.parent
NMAX = 32
EPS = 2.220446049250313e-16


def amplification(n, x):
    if n == 0:
        return mpf(1)
    odd = math.factorial(2 * n) // (2**n * math.factorial(n))
    return mpf(odd) / (2 * mpf(x)) ** n


def amplification_float(n, x):
    if n == 0:
        return 1.0
    odd = math.factorial(2 * n) // (2**n * math.factorial(n))
    return float(mpf(odd) / (2 * mpf(x)) ** n)


def read_reference(path):
    reference = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            reference[(row[0], float(row[1]))] = [mpf(v) for v in row[2:]]
    return reference


def least_squares_two(data):
    s11 = s12 = s22 = t1 = t2 = 0.0
    for n, y in data:
        ln = math.log(n)
        s11 += ln * ln
        s12 += ln * n
        s22 += float(n) * n
        t1 += ln * y
        t2 += n * y
    determinant = s11 * s22 - s12 * s12
    if determinant == 0.0:
        return None
    p = (t1 * s22 - t2 * s12) / determinant
    q = (s11 * t2 - s12 * t1) / determinant
    return p, q


def main():
    mp.dps = 70
    reference = read_reference(DEPOSIT / "amplification-form" / "data" / "reference.csv")

    rows_total = 0
    ulp_cells = 0
    ulp_worst = None
    ulp_records = []
    cross_cells = 0
    cross_worst = None
    cross_records = []
    cross_unmeasurable = 0

    current_key = None
    argument_index = -1
    seed_cross = {}

    with open(INPUTS / "paths.csv", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            rows_total += 1
            route, band, x_text, start_s, degree_s, n_s, value_s, pert_s, seed_s, transfer_s = row
            x = float(x_text)
            n = int(n_s)
            start = int(start_s)
            degree = int(degree_s)
            key = (band, x)

            if key != current_key:
                current_key = key
                argument_index += 1
                seed_cross = {}

            if route == "seed" or x == 0.0:
                continue

            value = float(value_s.replace("nan(ind)", "nan"))
            pert_delta = float(pert_s.replace("nan(ind)", "nan"))
            transfer_ulp = float(transfer_s.replace("nan(ind)", "nan"))

            # ---- check 1: exact-perturbation transfer ----------------------
            if rows_total % 97 == 0 and n != start:
                expected = amplification(n, x) / amplification(start, x)
                resolved = abs(pert_delta) > 16.0 * EPS * abs(value)
                if expected > 0 and resolved and math.isfinite(transfer_ulp):
                    deviation = float(abs(mpf(transfer_ulp) - expected) / expected)
                    ulp_cells += 1
                    ulp_records.append((deviation, band, x_text, start, degree, n,
                                        float(mpf(transfer_ulp)), float(expected),
                                        abs(pert_delta) / abs(value) if value else float("inf"),
                                        abs(pert_delta) / (EPS * abs(value)) if value else float("inf")))
                    if ulp_worst is None or deviation > ulp_worst:
                        ulp_worst = deviation

            # ---- check 2: start-independence ------------------------------
            if start not in (0, 32) and argument_index % 16 == 0:
                refs = reference[(band, x)]
                if n == start:
                    seed_cross[(start, degree)] = mpf(value) - refs[start]
                else:
                    seed_error = seed_cross.get((start, degree))
                    if seed_error is not None and seed_error != 0:
                        numerator = mpf(value) - refs[n]
                        transfer = numerator / seed_error
                        expected = amplification(n, x) / amplification(start, x)
                        if expected > 0 and mp.isfinite(transfer):
                            cross_cells += 1
                            deviation = abs(abs(transfer) - expected) / expected
                            expected_numerator = abs(mpf(seed_error)) * expected
                            floor = 16.0 * mpf(EPS) * abs(mpf(value))
                            unmeasurable = expected_numerator <= floor
                            cross_unmeasurable += 1 if unmeasurable else 0
                            cross_records.append((deviation, band, x_text, start, degree, n,
                                                  float(transfer), float(expected),
                                                  float(seed_error), value,
                                                  float(expected_numerator / floor) if floor else float("inf"),
                                                  bool(unmeasurable)))
                            if cross_worst is None or deviation > cross_worst:
                                cross_worst = deviation

    print("paths rows scanned:", rows_total)
    print()
    print("(1) exact-perturbation transfer, every 97th row, n != start")
    print("    cells:", ulp_cells)
    print("    worst relative deviation: %.16g" % ulp_worst)
    ulp_records.sort(reverse=True, key=lambda item: item[0])
    print("    records kept: %d; max by key %.16g; above 2: %d"
          % (len(ulp_records), ulp_records[0][0] if ulp_records else float("nan"),
             sum(1 for item in ulp_records if item[0] > 2.0)))
    print("    ten worst rows (deviation, band, x, start, degree, n, measured, closed form,")
    print("    |pert_delta|/|value|, that over eps.  The last ratio is the gate's own")
    print("    margin: a value near 16 is a perturbation one rounding away from noise.)")
    for rec in ulp_records[:12]:
        print("      %.6g  %s x=%s start=%d deg=%d n=%d  meas=%.6g closed=%.6g  pert/value=%.3g  /eps=%.3g" % rec)
    print("    the summary's worst cell, found whether or not it is in the sample:")
    for rec in ulp_records:
        if rec[1] == "A2" and rec[4] == 2 and rec[5] == 25:
            print("      %.6g  %s x=%s start=%d deg=%d n=%d  meas=%.6g closed=%.6g  pert/value=%.3g  /eps=%.3g" % rec)
    deviations = sorted(item[0] for item in ulp_records)
    percentile = lambda q: deviations[int(q * (len(deviations) - 1))]
    far = [item for item in ulp_records if item[9] > 4500.0]
    near = [item for item in ulp_records if item[9] <= 32.0]
    print("    deviation distribution: median %.3g  p90 %.3g  p99 %.3g" %
          (percentile(0.5), percentile(0.9), percentile(0.99)))
    print("    cells with the perturbation above 4500 eps of the value: %d, worst %.3g"
          % (len(far), max(item[0] for item in far) if far else float("nan")))
    print("    cells within a factor two of the 16-eps gate (<= 32 eps): %d, worst %.3g"
          % (len(near), max(item[0] for item in near) if near else float("nan")))


    print()
    print("(2) start-independence, non-canonical starts, every 16th argument")
    print("    cells:", cross_cells)
    print("    worst relative deviation: %.16g" % cross_worst)
    print("    unmeasurable cells by the 16-eps gate on the expected numerator:", cross_unmeasurable,
          "of", cross_cells)
    cross_records.sort(reverse=True, key=lambda r: r[0])
    print("    ten worst rows (deviation, band, x, start, degree, n, measured, closed form,")
    print("    seed error, |value_n|, expected numerator / 16-eps floor)")
    for rec in cross_records[:10]:
        dev, band, x_text, start, degree, n, transfer, expected, seed, value, ratio, unmeas = rec
        print("      %.6g  %s x=%s start=%d deg=%d n=%d  meas=%.6g closed=%.6g  seed=%.3g  value=%.3g  ratio=%.3g %s"
              % (dev, band, x_text, start, degree, n, transfer, expected, seed, value, ratio,
                 "UNMEASURABLE" if unmeas else ""))
    measurable = [rec for rec in cross_records if not rec[11]]
    measurable.sort(reverse=True, key=lambda r: r[0])
    print("    measurable cells: %d of %d" % (len(measurable), len(cross_records)))
    if measurable:
        print("    worst measurable deviation: %.6g at %s start=%d deg=%d n=%d x=%s"
              % (measurable[0][0], measurable[0][1], measurable[0][3], measurable[0][4],
                 measurable[0][5], measurable[0][2]))
        print("    measurable cells with deviation above 1: %d" %
              sum(1 for rec in measurable if rec[0] > 1))
    top100 = cross_records[:100]
    print("    unmeasurable among the 100 worst cells: %d" %
          sum(1 for rec in top100 if rec[11]))

    # ---- check 3: the n^p q^n fit -----------------------------------------
    fit_rows = []
    with open(DEPOSIT / "amplification-form" / "results" / "transfer_cells.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["measurable"] != "1":
                continue
            n = int(row["n"])
            transfer = float(row["transfer"])
            if n >= 1 and transfer > 0:
                fit_rows.append((row["band"], float(row["x"]), n, transfer))
    data = [(n, math.log(transfer)) for band, x, n, transfer in fit_rows]
    p, q = least_squares_two(data)
    residuals = [y - (p * math.log(n) + q * n) for n, y in data]
    rms = math.sqrt(sum(r * r for r in residuals) / len(residuals))
    print()
    print("(3) n^p q^n with two free parameters")
    print("    cells:", len(data), " p=%.16g q=%.16g" % (p, q))
    print("    maximum |residual| in log units: %.16g" % max(abs(r) for r in residuals))
    print("    rms: %.16g" % rms)
    print("    range of log(transfer): %.6g .. %.6g" % (min(y for n, y in data), max(y for n, y in data)))
    print("    range of n*log(x): %.6g .. %.6g"
          % (min(n * math.log(x) for band, x, n, t in fit_rows),
             max(n * math.log(x) for band, x, n, t in fit_rows)))

    # per-argument refit of the same two-parameter family
    by_argument = {}
    for band, x, n, transfer in fit_rows:
        by_argument.setdefault((band, x), []).append((n, math.log(transfer)))
    per_arg_worst = []
    for key, rows in by_argument.items():
        if len(rows) < 4:
            continue
        solution = least_squares_two(rows)
        if solution is None:
            continue
        pp, qq = solution
        worst = max(abs(y - (pp * math.log(n) + qq * n)) for n, y in rows)
        per_arg_worst.append((worst, key, pp, qq))
    per_arg_worst.sort(reverse=True)
    print("    per-argument refit of the same family: %d arguments, worst residual %.6g"
          % (len(per_arg_worst), per_arg_worst[0][0] if per_arg_worst else float("nan")))
    print("    worst five arguments:", [(round(w, 6), k) for w, k, pp, qq in per_arg_worst[:5]])

    # three-parameter fit with the x slope free
    s11 = s12 = s13 = s22 = s23 = s33 = t1 = t2 = t3 = 0.0
    for band, x, n, transfer in fit_rows:
        y = math.log(transfer)
        a1 = math.log(n)
        a2 = float(n)
        a3 = n * math.log(x)
        s11 += a1 * a1; s12 += a1 * a2; s13 += a1 * a3
        s22 += a2 * a2; s23 += a2 * a3; s33 += a3 * a3
        t1 += a1 * y; t2 += a2 * y; t3 += a3 * y
    # solve by Cramer's rule
    det = (s11 * (s22 * s33 - s23 * s23)
           - s12 * (s12 * s33 - s23 * s13)
           + s13 * (s12 * s23 - s22 * s13))
    det1 = (t1 * (s22 * s33 - s23 * s23)
            - s12 * (t2 * s33 - s23 * t3)
            + s13 * (t2 * s23 - s22 * t3))
    det2 = (s11 * (t2 * s33 - s23 * t3)
            - t1 * (s12 * s33 - s23 * s13)
            + s13 * (s12 * t3 - t2 * s13))
    det3 = (s11 * (s22 * t3 - t2 * s23)
            - s12 * (s12 * t3 - t2 * s13)
            + t1 * (s12 * s23 - s22 * s13))
    c1, c2, c3 = det1 / det, det2 / det, det3 / det
    residuals3 = [math.log(transfer) - (c1 * math.log(n) + c2 * n + c3 * n * math.log(x))
                  for band, x, n, transfer in fit_rows]
    print("    three-parameter refit log T = c1 log n + c2 n + c3 n log x:")
    print("      c1=%.16g c2=%.16g c3=%.16g (c3 = -1 would be the derived x^-n)" % (c1, c2, c3))
    print("      maximum |residual|: %.6g  rms: %.6g"
          % (max(abs(r) for r in residuals3),
             math.sqrt(sum(r * r for r in residuals3) / len(residuals3))))


if __name__ == "__main__":
    main()
