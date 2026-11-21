#!/usr/bin/env python3
"""The class test with one solver on both functions.

The alpha-independence test compares F1's fitted width-law slope with F0's at
the same edge and target.  The comparison is clean only when both functions'
counts come from the same solver.  This script builds that matched comparison
from the independent solver's records at both edges:

  F1 at x0 and a2   results-f1-independent/run_f1_{x0,a2}.json
  F0 at x0          results-corrected/corrected-counts.json
  F0 at a2          results-second-edge-a2/analysis.json

It also carries the earlier mixed variant, in which F0's rational rows are the
kernel's records (results-f1-independent/f1_independent_summary.json, the
class_test block).  That variant compares solvers and is reported as a
sensitivity, not as the alpha test.

The statistic is the paired slope difference: the same thirteen scaling widths
are used for both fits, the two residuals correlate, and the paired standard
error is the right scale.  The script writes
results-f1-independent/class_test_matched.json and prints the table.

    python class_test_matched.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGETS = (
    ("rel:1e-13", "rel", 1e-13),
    ("abs:5e-14", "abs", 5e-14),
    ("rel:1e-14", "rel", 1e-14),
)
DYAD = range(-7, 6)


def ols(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    beta = sxy / sxx
    alpha = my - beta * mx
    residuals = [y - (alpha + beta * x) for x, y in zip(xs, ys)]
    s2 = sum(r * r for r in residuals) / (n - 2)
    return {
        "alpha": alpha,
        "beta": beta,
        "se": math.sqrt(s2 / sxx),
        "residuals": residuals,
        "n": n,
        "sxx": sxx,
    }


def paired_slope_test(f1, f0):
    diff = f1["beta"] - f0["beta"]
    se_unpaired = math.sqrt(f1["se"] ** 2 + f0["se"] ** 2)
    n = f1["n"]
    cov = sum(a * b for a, b in zip(f1["residuals"], f0["residuals"])) / (n - 2)
    se_paired = math.sqrt(max(
        0.0,
        f1["se"] ** 2 + f0["se"] ** 2 - 2 * cov / f1["sxx"],
    ))
    return {
        "diff": diff,
        "se_unpaired": se_unpaired,
        "se_paired": se_paired,
        "z_unpaired": diff / se_unpaired if se_unpaired else None,
        "z_paired": diff / se_paired if se_paired else None,
        "residual_cov": cov,
    }


def load_f1(path, target):
    """One target's F1 independent counts over the thirteen scaling widths."""
    document = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for cell in document["cells"]:
        if cell.get("group") != "A":
            continue
        if abs(cell["target"] - target[2]) > 1e-30:
            continue
        if cell.get("kind") != target[1]:
            continue
        width_name = cell["width_name"]
        if width_name in ("ship", "paper"):
            continue
        exponent = -int(width_name[1:]) if width_name.startswith("m") else int(width_name[1:])
        out[str(exponent)] = {
            "stored": cell["corrected_count"],
            "below_meets": bool(cell.get("below_meets")),
        }
    return out


def load_f0_x0(path):
    """The independent solver's F0 x0 counts, by target and width."""
    document = json.loads(path.read_text(encoding="utf-8"))
    out = {target[0]: {} for target in TARGETS}
    for cell in document["cells"]:
        width_name = cell["width_name"]
        if width_name in ("ship", "paper"):
            continue
        exponent = -int(width_name[1:]) if width_name.startswith("m") else int(width_name[1:])
        for target in TARGETS:
            if cell["target"] != target[0]:
                continue
            best_below = cell.get("best_below")
            out[target[0]][str(exponent)] = {
                "stored": cell["corrected_count"],
                "below_meets": (
                    best_below is not None and best_below <= target[2]
                ),
            }
    return out


def load_f0_a2(path):
    """The independent solver's F0 a2 counts, by target and width."""
    document = json.loads(path.read_text(encoding="utf-8"))
    out = {target[0]: {} for target in TARGETS}
    for target in TARGETS:
        for exponent in DYAD:
            name = "p{}".format(exponent) if exponent >= 0 else "m{}".format(-exponent)
            cell = document["cells"]["{}/{}".format(target[0], name)]
            below_ratio = cell.get("below_ratio_a2")
            out[target[0]][str(exponent)] = {
                "stored": cell["count_a2"],
                "below_meets": below_ratio is not None and below_ratio <= 1.0,
            }
    return out


def matched_blocks(f1, f0):
    """The twelve blocks, with both functions from the independent solver."""
    blocks = {}
    for edge in ("x0", "a2"):
        blocks[edge] = {}
        for target in TARGETS:
            blocks[edge][target[0]] = {}
            for family in ("cheb", "rat"):
                if family == "rat":
                    members = []
                    for exponent in DYAD:
                        key = str(exponent)
                        e1 = f1[edge][target[0]][key]
                        e0 = f0[edge][target[0]][key]
                        members.append((exponent, e1, e0))
                    xs = [float(exponent) for exponent, _, _ in members]
                    y1 = [math.log2(e1["stored"]) for _, e1, _ in members]
                    y0 = [math.log2(e0["stored"]) for _, _, e0 in members]
                    f1fit = ols(xs, y1)
                    f0fit = ols(xs, y0)
                    paired = paired_slope_test(f1fit, f0fit)
                    blocks[edge][target[0]][family] = {
                        "n": len(members),
                        "f1_slope": f1fit["beta"],
                        "f1_se": f1fit["se"],
                        "f0_slope": f0fit["beta"],
                        "f0_se": f0fit["se"],
                        **paired,
                    }
    return blocks


def main() -> int:
    f1 = {
        "x0": {
            target[0]: load_f1(
                HERE.parent / "results-f1-independent" / "run_f1_x0.json", target)
            for target in TARGETS
        },
        "a2": {
            target[0]: load_f1(
                HERE.parent / "results-f1-independent" / "run_f1_a2.json", target)
            for target in TARGETS
        },
    }
    f0 = {
        "x0": load_f0_x0(HERE.parent / "results-corrected" / "corrected-counts.json"),
        "a2": load_f0_a2(HERE.parent / "results-second-edge-a2" / "analysis.json"),
    }
    matched = matched_blocks(f1, f0)
    summary = json.loads(
        (HERE.parent / "results-f1-independent" / "f1_independent_summary.json")
        .read_text(encoding="utf-8"))
    for edge in ("x0", "a2"):
        for target in TARGETS:
            block = summary["class_test"][edge][target[0]]["cheb"]
            matched[edge][target[0]]["cheb"] = {
                "n": block["n"],
                "f1_slope": block["f1_slope"],
                "f1_se": block["f1_se"],
                "f0_slope": block["f0_slope"],
                "f0_se": block["f0_se"],
                "diff": block["diff"],
                "z_paired": block["z_paired"],
                "z_unpaired": block["z_unpaired"],
                "dropped": block["dropped"],
                "source": "kernel records for both functions",
            }
    mixed = {}
    for edge in ("x0", "a2"):
        mixed[edge] = {}
        for target in TARGETS:
            block = summary["class_test"][edge][target[0]]["rat"]
            mixed[edge][target[0]] = {
                "n": block["n"],
                "f1_slope": block["f1_slope"],
                "f1_se": block["f1_se"],
                "f0_slope": block["f0_slope"],
                "f0_se": block["f0_se"],
                "diff": block["diff"],
                "z_paired": block["z_paired"],
                "z_unpaired": block["z_unpaired"],
                "dropped": block["dropped"],
                "source": "F1 independent; F0 kernel record",
            }

    relative = [
        (edge, target[0], matched[edge][target[0]][family])
        for edge in ("x0", "a2")
        for target in TARGETS
        if target[1] == "rel"
        for family in ("cheb", "rat")
    ]
    above_two = [
        (edge, target, block["diff"], block["z_paired"])
        for edge, target, block in relative
        if block["z_paired"] is not None and block["z_paired"] > 2.0
    ]
    out = {
        "schema": "coefficient-counts/class-test-matched/1",
        "description": (
            "The alpha-independence class test with the independent solver's "
            "records for F1 and F0 at both edges. The mixed variant, in which "
            "F0's rational rows are the kernel's, is carried as a sensitivity."
        ),
        "method": (
            "OLS of log2(stored count) on log2(width) over the thirteen "
            "scaling widths, paired slope difference and paired standard "
            "error."
        ),
        "targets": [target[0] for target in TARGETS],
        "matched": matched,
        "mixed_sensitivity": mixed,
        "relative_tally": {
            "blocks": len(relative),
            "above_two_paired_se": len(above_two),
            "which": [
                {"edge": edge, "target": target,
                 "diff": diff, "z_paired": z}
                for edge, target, diff, z in above_two
            ],
            "min_diff": min(block["diff"] for _, _, block in relative),
            "max_diff": max(block["diff"] for _, _, block in relative),
        },
    }
    out_path = HERE.parent / "results-f1-independent" / "class_test_matched.json"
    out_path.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")

    for edge in ("x0", "a2"):
        for target in TARGETS:
            for family in ("cheb", "rat"):
                block = matched[edge][target[0]][family]
                print(
                    "{:2s} {:10s} {:4s}: n={:2d} F1 {:.5f} F0 {:.5f} "
                    "diff {:+.5f} z {:+.2f}".format(
                        edge, target[0], family, block["n"],
                        block["f1_slope"], block["f0_slope"],
                        block["diff"], block["z_paired"]))
    print("wrote", out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
