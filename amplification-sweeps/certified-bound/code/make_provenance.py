"""Write <run-dir>/provenance.json: every input this record reads, its
digest, and the evidence behind each assumption.

All paths are resolved by record_paths.py from the command line, the
environment, or defaults relative to the record.  The results section lists
what the selected run directory actually contains, so a provenance written
after the certificate, the replay, the oracle check and the diagnostics
covers those files.
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model
import record_paths

CORE_SWEEP = ["results/cells.csv", "results/per_order.csv",
              "results/provenance.json", "data/reference.csv"]
OPTIONAL_SWEEP = ["results/prediction.csv", "results/band_summary.csv",
                  "results/summary.txt", "analyze.py"]
CORE_DEPOSIT = ["bounds.py", "exact.py", "coeffs.py", "core.py"]
OPTIONAL_DEPOSIT = ["run_all.py", "exact_trunc_cert.py",
                    "exact_amplification_form.csv", "README.md"]
EXP_NAMES = ["exp_probe.cpp", "exp_check.py", "exp_bits.txt",
             "exp_bits_reproduced.txt", "exp_check.log"]


def sha(path):
    with open(path, "rb") as h:
        return hashlib.sha256(h.read()).hexdigest()


def files(base, names):
    out = {}
    for name in names:
        p = Path(base) / name
        out[name] = {"bytes": p.stat().st_size, "sha256": sha(p)}
    return out


def files_present(base, names):
    out = {}
    for name in names:
        p = Path(base) / name
        if p.is_file():
            out[name] = {"bytes": p.stat().st_size, "sha256": sha(p)}
        else:
            out[name] = None
    return out


def git(root, *args):
    return subprocess.run(["git", "-C", str(root), *args],
                          capture_output=True, text=True,
                          check=True).stdout.strip()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    record_paths.add_kernel(parser)
    record_paths.add_inputs(parser, model_record_arg=True, deposit_arg=True)
    record_paths.add_out(parser)
    parser.add_argument("--manuscript", metavar="FILE",
                        help="manuscript file read for the derivation (optional)")
    args = parser.parse_args(argv)

    kernel = record_paths.kernel(args.kernel)
    sweep = record_paths.sweep(args.sweep)
    model_record = record_paths.model_record(args.model_record)
    deposit = record_paths.deposit(args.deposit)
    exp_bits = record_paths.exp_bits(args.exp_bits)
    out = record_paths.out_dir(args.out)

    model.configure(kernel=args.kernel, exp_bits=args.exp_bits)
    tables = model.table()

    sweep_prov = {}
    sweep_prov_path = Path(sweep) / "results" / "provenance.json"
    if sweep_prov_path.is_file():
        sweep_prov = json.load(open(sweep_prov_path))

    record = {
        "record": "certified per-call rounding bound for the two region-A Boys "
                  "entries, Paper A (Boys route comparison), repaired "
                  "underflow radius",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": {
            "python": sys.version.split()[0],
            "mpmath": mpmath.__version__,
            "platform": platform.platform(),
        },
        "model_constants": {
            "u": 2.0 ** -53,
            "u_add": 2.0 ** -53 / (1.0 - 2.0 ** -53),
            "min_normal": model.MIN_NORMAL,
            "underflow_allowance": model.UNDERFLOW,
            "safety": model.SAFETY,
        },
        "paths": {
            "resolution": "record_paths.py: command line, then environment, "
                          "then a default relative to the record",
            "kernel": str(kernel),
            "sweep": str(sweep),
            "model_record": str(model_record),
            "deposit": str(deposit),
            "exp_bits": str(exp_bits),
            "run_directory": str(out),
        },
        "kernel": {
            "tag": git(kernel, "describe", "--tags"),
            "commit": git(kernel, "rev-parse", "HEAD"),
            "release": tables["release"],
            "parsed_table_digest": model.FIT_CONTENT_DIGEST,
            "headers": files(kernel / "include" / "boys", [
                "boys_coefficients.hpp",
                "boys_coefficients_d_coarse_cheb.hpp",
                "boys_impl.hpp"]),
            "generator": files(kernel / "tools", ["gen_boys_coefficients.py"]),
            "coefficient_identity": (
                "the model parses kPieces, kPieceStart, kTierThresholds, "
                "kExtendedBcoeffs, kExtendedBX0 and kX0 from "
                "boys_coefficients.hpp and the Chebyshev table from "
                "boys_coefficients_d_coarse_cheb.hpp; the parsed content "
                "digest is the same in the v5.1.1 sweep tree and the v5.2.5 "
                "and v5.2.6 releases, and the coarse table's whole-file "
                "digest is unchanged, so the fits, thresholds and extended "
                "seed the sweep evaluated are the ones this record bounds.  "
                "The v5.2.5 header differs from v5.1.1 only in the rational "
                "region-B seed, which the named policy does not read; the "
                "v5.2.6 header is byte-identical to the v5.2.5 header."),
        },
        "sweep": {
            "path": str(sweep),
            "files": files(sweep, CORE_SWEEP),
            "optional_files": files_present(sweep, OPTIONAL_SWEEP),
            "boys_commit": sweep_prov.get("boys", {}).get("commit"),
            "boys_describe": sweep_prov.get("boys", {}).get("describe"),
            "compiler": sweep_prov.get("compiler", {}).get("path"),
            "flags": sweep_prov.get("flags"),
            "machine": sweep_prov.get("machine"),
            "model_functions_used": [
                "cells.csv: the delivered recursed and direct doubles",
                "per_order.csv: each entry's worst absolute and relative cells",
                "reference.csv: 60-digit F_n at 836 arguments",
            ],
        },
        "model_record": {
            "path": str(model_record),
            "files": files(model_record, ["prediction.csv",
                                          "prediction_record.csv"]),
            "role": "the paper model's 132 comparisons, its 65 commitments "
                    "and the band-max arguments of each entry; the deposited "
                    "sweep's prediction.csv is a different model revision "
                    "(67 commitments), so the certificate reads the paper's "
                    "record deposited here",
        },
        "deposit": {
            "path": str(deposit),
            "files": files(deposit, CORE_DEPOSIT),
            "optional_files": files_present(deposit, OPTIONAL_DEPOSIT),
            "role": "the existing a priori-bound method's certified truncation "
                    "intervals and its seven variants; this record keeps them "
                    "and replaces the heuristic charge with the per-call "
                    "running-error radius",
            "not_modified": True,
        },
        "exponential": {
            "bit_file": str(exp_bits),
            "files": files(Path(exp_bits).parent, [Path(exp_bits).name]),
            "probe_files": files_present(Path(exp_bits).parent, EXP_NAMES),
            "bit_input": files(record_paths.RECORD_DIR / "inputs",
                               ["x_bits.txt"]),
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
            "evidence": "the deposited probe output prints std::exp(-x) bits "
                        "at every argument; the worst relative departure is "
                        "0.999662996081 * 2^-53 (0.4998 ulp), so eta_E <= "
                        "2^-52 holds with a 2x margin; Python math.exp matches "
                        "bit for bit at all 9541 arguments and 0.5*exp(-x) is "
                        "exact",
            "assumption": "the bound is conditional on eta_E <= 2^-52, which "
                          "the probe supports on this machine and toolchain; a "
                          "different libm is outside the evidence",
        },
        "assumptions": {
            "rounding": "round-to-nearest, IEEE-754 binary64, no flush-to-zero",
            "evaluation": "the kernel's exact operation sequence at "
                          "EvalScheme::kSplitClenshaw, "
                          "DivisionForm::kRefinedReciprocal, with the hoisted "
                          "exp(-x)/2 shared by the recurrence",
            "overflow": "no overflow at the tested cells; the largest delivered "
                        "magnitude is recorded in the run's summary.json",
            "underflow": (
                "the grid includes the least positive subnormal x=2**-1074. "
                "In the downward chain the product x*F_32 has an exact value "
                "about 7.60e-326, which rounds to zero, so the multiplication "
                "loses a nonzero exact result below the subnormal range.  The "
                "repaired primitive rule carries the absolute allowance "
                "2**-1071 = 16*2**-1075 whenever the computed result or the "
                "radius bracket is below 2**-1022: one half-ulp for the "
                "IEEE-754 round-to-nearest underflow bound 2**-1075, which is "
                "not itself a binary64 number, and fifteen for the absolute "
                "roundings of the radius arithmetic in the subnormal range.  "
                "The allowance propagates through the chain and the regression "
                "test records the exact product, the computed zero and the "
                "radius at the affected step."),
            "radius": "first-order running-error sums with the absolute "
                      "underflow allowance above and SAFETY=1+2^-40; every "
                      "rounding at u=2^-53, accumulation in double",
            "truncation": "the exact-real route value at 70-digit mpmath minus "
                          "the 60-digit reference, with EPS_BASE=1e-45 and "
                          "2^-51|t| for the delivered comparison; the deposit's "
                          "independent identity oracle is cross-checked on 60 "
                          "cells",
        },
    }

    if args.manuscript:
        manuscript = Path(args.manuscript).expanduser().resolve()
        record["manuscript"] = {
            "path": str(manuscript),
            "sha256": sha(manuscript),
            "read_only": True,
            "section_read": "S18, the certified enclosure derived",
        }

    record["record_files"] = {
        "code": files(record_paths.CODE_DIR,
                      sorted(p.name for p in record_paths.CODE_DIR.glob("*.py"))),
        "explore": files(record_paths.EXPLORE_DIR,
                         sorted(p.name for p in record_paths.EXPLORE_DIR.glob("*.py"))),
        "exp": files(record_paths.RECORD_DIR / "exp", EXP_NAMES),
    }
    results = {}
    for p in sorted(out.iterdir()):
        if p.is_file() and p.name != "provenance.json":
            results[p.name] = {"bytes": p.stat().st_size, "sha256": sha(p)}
    record["results"] = results

    with open(out / "provenance.json", "w") as h:
        json.dump(record, h, indent=1)
    print("wrote", out / "provenance.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
