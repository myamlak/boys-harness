"""Write _a-cert-run/provenance.json: every input this record reads, its
digest, and the evidence behind each assumption.
"""
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

ROOT = r"C:/Users/myaml/dev"
CERT = ROOT + "/_a-cert-run"
BOYS = ROOT + "/_boys-tags"
SWEEP = ROOT + "/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
DEPOSIT = ROOT + "/_harness-v1023/amplification-sweeps/apriori-bound"


def sha(path):
    with open(path, "rb") as h:
        return hashlib.sha256(h.read()).hexdigest()


def files(base, names):
    out = {}
    for name in names:
        p = base + "/" + name
        out[name] = {"bytes": os.path.getsize(p), "sha256": sha(p)}
    return out


def main():
    tag = subprocess.run(["git", "-C", BOYS, "rev-parse", "HEAD"],
                         capture_output=True, text=True, check=True).stdout.strip()
    describe = subprocess.run(["git", "-C", BOYS, "describe", "--tags"],
                              capture_output=True, text=True, check=True).stdout.strip()
    sweep_prov = json.load(open(SWEEP + "/results/provenance.json"))
    record = {
        "record": "certified per-call rounding bound for the two region-A Boys "
                  "entries, Paper A (Boys route comparison)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel": {
            "tag": describe,
            "commit": tag,
            "tree": BOYS,
            "headers": files(BOYS + "/include/boys", [
                "boys_coefficients.hpp",
                "boys_coefficients_d_coarse_cheb.hpp",
                "boys_impl.hpp"]),
            "generator": files(BOYS + "/tools", ["gen_boys_coefficients.py"]),
            "coefficient_identity_with_sweep": (
                "the sweep's provenance.json records boys_coefficients.hpp "
                "sha256 8ce5e17157a4be5f621d695f377094d212cb839732e278d47f23aa225b4d9204 "
                "in its v5.1.1 tree; the v5.2.4 header has the same digest, so "
                "the stored fits and thresholds the sweep evaluated are the "
                "ones this record bounds"),
        },
        "sweep": {
            "path": SWEEP,
            "provenance": files(SWEEP + "/results", ["provenance.json"]),
            "boys_commit": sweep_prov["boys"]["commit"],
            "boys_describe": sweep_prov["boys"]["describe"],
            "compiler": sweep_prov["compiler"]["path"],
            "flags": sweep_prov["flags"],
            "machine": sweep_prov["machine"],
            "files": files(SWEEP + "/results", [
                "cells.csv", "per_order.csv", "prediction.csv",
                "prediction_record.csv", "band_summary.csv", "summary.txt"]),
            "reference": files(SWEEP + "/data", ["reference.csv"]),
            "model_code": files(SWEEP, ["analyze.py"]),
            "model_functions_used": [
                "prediction.csv / prediction_record.csv: the 132 comparisons "
                "and the model's 65 commitments",
                "per_order.csv: each entry's worst absolute and relative cells",
                "cells.csv: the delivered recursed and direct doubles",
                "reference.csv: 60-digit F_n at 836 arguments",
            ],
        },
        "deposit": {
            "path": DEPOSIT,
            "files": files(DEPOSIT, ["bounds.py", "exact.py", "run_all.py",
                                     "exact_trunc_cert.py",
                                     "exact_amplification_form.csv",
                                     "README.md"]),
            "role": "the existing a priori-bound method's certified truncation "
                    "intervals and its seven variants; this record keeps them "
                    "and replaces the heuristic charge with the per-call "
                    "running-error radius",
            "not_modified": True,
        },
        "manuscript": {
            "path": SWEEP.rsplit("/sweeps/", 1)[0] + "/main.tex",
            "sha256": sha(SWEEP.rsplit("/sweeps/", 1)[0] + "/main.tex"),
            "read_only": True,
            "section_read": "5.6 (lines 786-814), the a priori failure "
                            "paragraph whose 16/49, 9.6 and 93 are checked",
        },
        "exponential": {
            "probe_source": files(CERT + "/exp", ["exp_probe.cpp", "exp_check.py"]),
            "probe_binary": files(CERT + "/exp", ["exp_probe.exe"]),
            "bit_file": files(CERT + "/exp", ["exp_bits.txt"]),
            "bit_input": files(CERT + "/exp", ["xbits_in.txt"]),
            "reproduced_output": files(CERT + "/exp",
                                       ["exp_bits_reproduced.txt"]),
            "analysis_log": files(CERT + "/exp", ["exp_check.log"]),
            "reproduction": "the deposited probe binary rerun on the input bits "
                            "recovered from exp_bits.txt reproduces the file "
                            "byte for byte (cmp exit 0, 2026-11-03); the "
                            "build flags above are the measuring session's",
            "arguments": 9541,
            "compiler": "MSVC 14.51.36231, host x64",
            "flags": "/nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2",
            "worst_relative_over_2^-53": 0.999662996081,
            "worst_argument": 11.084208097312864,
            "arguments_above_half_unit": 2585,
            "arguments_above_one_unit": 0,
            "python_math_exp_bit_mismatches": 0,
            "half_scaling_mismatches": 0,
            "conditional_eta_E": "2^-52",
            "evidence": "the probe prints std::exp(-x) bits at every argument; "
                        "the worst relative departure is 0.999662996081 * 2^-53 "
                        "(0.4998 ulp), so eta_E <= 2^-52 holds with a 2x margin; "
                        "Python math.exp matches bit for bit at all 9541 "
                        "arguments and 0.5*exp(-x) is exact",
            "assumption": "the bound is conditional on eta_E <= 2^-52, which "
                          "the probe supports on this machine and toolchain; a "
                          "different libm is outside the evidence",
        },
        "oracle_cross_check": {
            "files": files(CERT + "/results",
                           ["oracle_check.csv", "oracle_check.log",
                            "dps_check.log"]),
            "scripts": files(CERT + "/code", ["oracle_check.py"]) | files(
                CERT + "/explore", ["dps_check.py"]),
            "cells": 60,
            "findings": "59 of 60 cells agree with the deposit's independent "
                        "identity oracle to 8.2e-50 or better; one small-x "
                        "product cell (A1, x=4.0666245819265257e-4, n=24) "
                        "differs by 7.99e-24 because the oracle's dps-70 "
                        "quadrature is inaccurate there.  Our route's "
                        "truncation is stable at dps 70/100/140 and agrees "
                        "with an independent hyp1f1 evaluation of F_n to "
                        "4.4e-65; the oracle converges to our value at dps "
                        "140.  The record's truncation values are unaffected.",
            "relevance": "evidence that the exact-real route evaluates the "
                         "same truncation the deposit certifies; the outlier "
                         "is an oracle quadrature limit, not a route error",
        },
        "results": {
            "per_cell": files(CERT + "/results",
                              ["cells_cert.csv", "band_cert.csv",
                               "model_cells.csv", "harness_cells.csv",
                               "tiers_cells.csv"]),
            "summaries": files(CERT + "/results",
                               ["summary.json", "tiers_summary.json",
                                "commands.md"]),
            "logs": files(CERT + "/results",
                          ["certify.log", "harness.log", "replay.log",
                           "tiers.log", "oracle_check.log", "dps_check.log",
                           "paper93.log"]),
        },
        "assumptions": {
            "rounding": "round-to-nearest, IEEE-754 binary64, no flush-to-zero",
            "evaluation": "the kernel's exact operation sequence at "
                          "EvalScheme::kSplitClenshaw, DivisionForm::kRefinedReciprocal, "
                          "with the hoisted exp(-x)/2 shared by the recurrence",
            "overflow": "no overflow at the tested cells; the largest delivered "
                        "magnitude is recorded in results/summary.json",
            "underflow": "the reference grid includes x=0 and the subnormal "
                         "x=5e-324; the x=0 branch is one correctly rounded "
                         "division per order, and at the subnormal argument no "
                         "term carrying the result underflows (the delivered "
                         "value is inside the interval, as validation shows)",
            "radius": "first-order running-error sums with SAFETY=1+2^-40; "
                      "every rounding at u=2^-53, accumulation in double",
            "truncation": "the exact-real route value at 70-digit mpmath minus "
                          "the 60-digit reference, with EPS_BASE=1e-45 and "
                          "2^-51|t| for the delivered comparison; the deposit's "
                          "independent identity oracle agrees on the checked cells",
        },
        "record_files": {
            "code": files(CERT + "/code", [
                "model.py", "exactroute.py", "certify.py", "tiers.py",
                "replay.py", "oracle_check.py", "make_provenance.py"]),
            "explore": files(CERT + "/explore", [
                "harnesscells.py", "dps_check.py", "paper93.py"]),
            "exp": files(CERT + "/exp", ["exp_probe.cpp", "exp_check.py"]),
        },
    }
    with open(CERT + "/provenance.json", "w") as h:
        json.dump(record, h, indent=1)
    print("wrote", CERT + "/provenance.json")


if __name__ == "__main__":
    main()
