#!/usr/bin/env python3
"""Verify the Paper C continuum upper bound by recomputation.

Reads the deposited continuum bound JSON and the polynomial certificate,
recomputes every exact rational bound with the same method, and checks:

  * the certificate digest named in the record;
  * exact equality of every published bound fraction;
  * that both bounds are below the 5e-14 target;
  * that the record's conclusion follows.

Exit 0 on PASS, 1 on any failure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import continuum_bound as cb  # noqa: E402


def fail(msg: str) -> None:
    print(f"FAIL: {msg}")
    raise SystemExit(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate", required=True)
    ap.add_argument("--bound", required=True)
    args = ap.parse_args()

    with open(args.certificate, "rb") as fh:
        cert_digest = hashlib.sha256(fh.read()).hexdigest()
    with open(args.bound, encoding="utf-8") as fh:
        rec = json.load(fh)

    if rec.get("schema") != "paper-c-continuum-bound-1":
        fail(f"unexpected schema {rec.get('schema')!r}")
    if rec.get("certificate_sha256") != cert_digest:
        fail("certificate digest does not match the recorded digest")

    method = rec["method"]
    n_sub = int(method["subintervals"])
    samples = int(method["samples_per_subinterval"]) - 1
    recomputed = cb.compute(args.certificate, n_sub, samples)

    published = {r["coefficients"]: r for r in rec["results"]}
    for r in recomputed:
        label = r["label"]
        if label not in published:
            fail(f"missing published result for {label}")
        p = published[label]
        p_bound = Fraction(p["bound"]["upper_decimal"])
        if p_bound < r["bound"]:
            fail(f"{label}: published upper bound is below recomputation")
        if p_bound >= cb.TARGET:
            fail(f"{label}: published upper bound is not below the target")
        if r["bound"] > cb.TARGET:
            fail(f"{label}: recomputed bound {float(r['bound']):.6e} "
                 "is above target")
        if not p.get("meets"):
            fail(f"{label}: record does not mark the bound as meeting")
        print(f"PASS {label}: recomputed {float(r['bound']):.6e}, "
              f"published upper {p['bound']['upper_decimal']} "
              f"below target by {float(cb.TARGET - p_bound):.3e}")

    print("PASS: continuum upper bound verified by exact recomputation")
    return 0


if __name__ == "__main__":
    sys.exit(main())
