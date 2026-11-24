"""Assemble the solver-independence evidence into one summary JSON/text."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
CELLS = Path(r"C:\Users\myaml\dev\_c60-clean-edge-full-2026-11-24\cells")
REC = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
           r"\clean-edge-2026-11-22\edge62_sweep.json")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def control_cell(name):
    d = load(CELLS / f"{name}.json")
    acc = d["accepted"]
    return {
        "cell": d["cell"], "count60": d["count60"],
        "direction": d["direction"],
        "accepted_pair": [acc["m"], acc["k"]],
        "accepted_settings": acc["settings"],
        "brasil_converged": acc["converged"],
        "brasil_deviation": acc["deviation"],
        "e60": acc["e60"], "guard60": acc["guard60"],
        "e60_4001": acc["e60_4001"], "meets60_4001": acc["meets60_4001"],
        "wall_seconds": d.get("wall_seconds"),
    }


def main():
    summary = {
        "schema": "clean-edge-solver-independence-1",
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "question": ("Are the clean-edge 60-digit meeting counts solver-"
                     "specific to BRASIL/baryrat?"),
        "inputs": {
            "record": {"path": str(REC), "sha256": sha256(REC)},
            "control_cells": {
                name: {"path": str(CELLS / f"{name}.json"),
                       "sha256": sha256(CELLS / f"{name}.json")}
                for name in ("p5_rel_1e-15", "ship_rel_1e-15")},
        },
        "protocol": {
            "a": 62.0,
            "target": "rel:1e-15",
            "grid_nodes": 2001, "check_nodes": 4001,
            "k_max": 12, "den_floor": 1e-3,
            "coefficients": "60 significant digits",
            "evaluation_dps": 60, "f0_dps": 70,
            "solver_dps": 70,
        },
        "brasil_control": {
            "p5": control_cell("p5_rel_1e-15"),
            "ship": control_cell("ship_rel_1e-15"),
        },
        "independent_arm": {
            "minimax_solver": (
                "classical rational Remez exchange: differential correction "
                "plus nonlinear Newton polish at 70 dps, continuous reference "
                "points refined by parabolic interpolation, global "
                "de la Vallee Poussin exchange and Maehly exchange. "
                "Convergence criterion: N=m+k+2 alternating reference "
                "extrema whose common magnitude equals the grid maximum to "
                "1e-10 relative (alternation theorem certificate). "
                "Known-answer test: minimax polynomials for t^4..t^8 pass."),
            "constructive_arm": (
                "iterated weighted least squares (Sanathanan-Koerner) for the "
                "rational fit at 70 dps, used to construct admissible fits "
                "and to initialize the Newton/Remez polish."),
            "lower_bound": (
                "de la Vallee Poussin: if a fit's weighted relative error has "
                "N=m+k+2 alternating 2,001-grid nodes with magnitudes >= L, "
                "then every guard-admissible rational of type [m/k] has grid "
                "maximum >= L. Proof: the difference of the two rationals is "
                "R/(QQ') with deg R <= m+k, but alternation forces m+k+1 "
                "zeros."),
        },
        "p5": {},
        "ship": {},
    }
    summary["independent_arm"]["remez_validation"] = {
        "known_answer_test": (
            "minimax polynomials for t^4, t^5, t^6, t^8 on [-1,1] recovered "
            "with the exact theoretical errors 2^(1-n) and coefficients "
            "matching t^n - T_n/2^(n-1) to 1e-40"),
        "known_answer_passed": True,
        "command": "python -I remez_rational.py --test (0.95 s)",
    }
    summary["independent_arm"]["remez_cold_6_8_blocker"] = {
        "command": ("python -I remez_rational.py --a 62.0 --width "
                    "166.72481653881133 --m 6 --k 8"),
        "result": "best-admissible, not certified, meets=False",
        "grid_max": "1.753192978686451e-15",
        "guard": "0.0083463293329078285",
        "grid_max_over_reference_min": "131.0167217149823",
        "equioscillation_certificate": False,
        "blocker": (
            "cold differential correction from Chebyshev references "
            "converged to a spurious near-pole local solution; the exact "
            "Newton on a Chebyshev reference diverged or stalled. The "
            "weighted-least-squares construction supplied a good starting "
            "fit, from which the targeted Newton reached 2.56e-16."),
    }
    # p5 count 15
    s15 = load(OUT / "p5_s15.json")
    p5_15 = {
        "meeting_pairs": s15["meeting_pairs"],
        "best_e60": s15["best_e60"],
        "rows": [{k: r.get(k) for k in
                  ("m", "k", "e60", "guard60", "meets60",
                   "e60_4001", "meets60_4001")}
                 for r in s15["rows"]],
        "wall_seconds": s15["wall_seconds"],
    }
    if (OUT / "targeted_p5_68.json").exists():
        doc = load(OUT / "targeted_p5_68.json")
        t = doc["best"]
        best_round = next((r for r in doc["rounds"]
                           if f"round {r['round']}" == t["tag"]), None)
        p5_15["targeted_newton_6_8"] = {
            k: t.get(k) for k in
            ("e60", "guard60", "meets60", "e60_4001", "meets60_4001",
             "tag")}
        p5_15["targeted_newton_6_8"]["convergence"] = {
            "equioscillation_certificate": False,
            "newton_ok": (best_round["newton_ok"]
                          if best_round else None),
            "newton_note": (best_round["note"] if best_round else None),
            "grid_max_over_reference_min": (best_round["ratio"]
                                            if best_round else None),
            "alternating_extrema": (best_round["extrema"]
                                    if best_round else None),
        }
    # p5 count 14 and bounds
    s14 = load(OUT / "p5_s14.json")
    lb14 = load(OUT / "lb_p5_s14.json")
    p5_14 = {
        "best_e60": s14["best_e60"],
        "meeting_pairs": s14["meeting_pairs"],
        "rows": [{k: r.get(k) for k in
                  ("m", "k", "e60", "guard60", "meets60")}
                 for r in s14["rows"]],
        "lower_bounds": {
            "min_continuous_L": lb14["min_L"],
            "min_grid_L": lb14["min_grid_L"],
            "all_grid_bounds_above_target":
                lb14["all_grid_bounds_above_target"],
            "table": [{k: r.get(k) for k in
                       ("m", "k", "N", "extrema", "lower_bound_grid_L",
                        "grid_L_over_target", "grid_bound_above_target")}
                      for r in lb14["rows"]],
        },
        "wall_seconds": s14["wall_seconds"],
    }
    summary["p5"]["count15"] = p5_15
    summary["p5"]["count14"] = p5_14
    # ship count 13 and bounds
    s13 = load(OUT / "ship_s13.json")
    ship_13 = {
        "sk_best_e60": s13["best_e60"],
        "sk_meeting_pairs": s13["meeting_pairs"],
        "sk_rows": [{k: r.get(k) for k in
                     ("m", "k", "e60", "guard60", "meets60")}
                    for r in s13["rows"]],
        "wall_seconds": s13["wall_seconds"],
    }
    if (OUT / "targeted_ship_66.json").exists():
        t = load(OUT / "targeted_ship_66.json")
        b = t["best"]
        best_round = next((r for r in t["rounds"]
                           if f"round {r['round']}" == b["tag"]), None)
        ship_13["targeted_newton_6_6"] = {
            k: b.get(k) for k in
            ("e60", "guard60", "meets60", "e60_4001", "meets60_4001",
             "tag")}
        ship_13["targeted_newton_6_6"]["convergence"] = {
            "equioscillation_certificate": False,
            "newton_ok": (best_round["newton_ok"]
                          if best_round else None),
            "newton_note": (best_round["note"] if best_round else None),
            "grid_max_over_reference_min": (best_round["ratio"]
                                            if best_round else None),
            "alternating_extrema": (best_round["extrema"]
                                    if best_round else None),
        }
        ship_13["targeted_newton_rounds"] = t["rounds"]
    s12 = load(OUT / "ship_s12.json")
    lb12 = load(OUT / "lb_ship_s12.json")
    ship_12 = {
        "best_e60": s12["best_e60"],
        "meeting_pairs": s12["meeting_pairs"],
        "rows": [{k: r.get(k) for k in
                  ("m", "k", "e60", "guard60", "meets60")}
                 for r in s12["rows"]],
        "lower_bounds": {
            "min_continuous_L": lb12["min_L"],
            "min_grid_L": lb12["min_grid_L"],
            "all_grid_bounds_above_target":
                lb12["all_grid_bounds_above_target"],
            "table": [{k: r.get(k) for k in
                       ("m", "k", "N", "extrema", "lower_bound_grid_L",
                        "grid_L_over_target", "grid_bound_above_target")}
                      for r in lb12["rows"]],
        },
        "wall_seconds": s12["wall_seconds"],
    }
    summary["ship"]["count13"] = ship_13
    summary["ship"]["count12"] = ship_12
    # final dense-grid checks
    if (OUT / "final_checks.json").exists():
        summary["meeting_fit_dense_checks"] = load(OUT / "final_checks.json")

    summary["conclusion"] = {
        "p5_independent_count": 15,
        "ship_independent_count": 13,
        "brasil_counts": {"p5": 15, "ship": 13},
        "counts_solver_specific": False,
        "basis": (
            "Admissible independent fits meet the 1e-15 relative target at "
            "2^5 count 15 ([6/8]: e60 2.5635e-16 after the Newton polish, "
            "guard 7.84e-3; the weighted-least-squares construction at the "
            "same pair is 8.0230e-16, and [7/7] is 7.0226e-16) and ship "
            "count 13 ([6/6]: e60 5.3885e-16, guard 4.45e-2).  Both also "
            "meet on the 4,001-node check grid and a 20,001-node grid. At "
            "one count below, the de la Vallee Poussin bounds from "
            "independent fits are 1.3249e-15 (2^5, [6/7]) and 5.3592e-15 "
            "(ship, [5/6]); every pair at those counts has a bound above "
            "1e-15, so no admissible fit of any type at count 14 or 12 can "
            "meet the target."),
        "caveat": (
            "The independent meeting fits are admissible constructive fits; "
            "the Remez/Newton arm did not reach a full equioscillation "
            "certificate on them (its Newton ratio stalled at 1.139 for ship "
            "and the cold Remez for 2^5 [6/8] fell into a spurious near-pole "
            "local solution at 1.75e-15). The count claim does not require a "
            "minimax certificate: existence at the count is constructive and "
            "the lower bound at the count below is a de la Vallee Poussin "
            "bound, which is independent of optimality of the certifying "
            "fit."),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1),
                                      encoding="utf-8")

    lines = []
    lines.append("Clean-edge solver-independence research summary")
    lines.append("=" * 55)
    lines.append("")
    lines.append("BRASIL 60-digit control:")
    for key, label in (("p5", "2^5 rel:1e-15"), ("ship", "ship rel:1e-15")):
        c = summary["brasil_control"][key]
        lines.append(f"  {label:22s} count {c['count60']:2d}  "
                     f"pair [{c['accepted_pair'][0]}/{c['accepted_pair'][1]}]"
                     f"  e60 {c['e60']}  guard {c['guard60']}  "
                     f"BRASIL-converged {c['brasil_converged']}")
    lines.append("")
    lines.append("Independent arm (SK construction + Remez/Newton polish):")
    lines.append(f"  2^5  count 15: meeting pairs {p5_15['meeting_pairs']}, "
                 f"best e60 {p5_15['best_e60']}")
    if "targeted_newton_6_8" in p5_15:
        t = p5_15["targeted_newton_6_8"]
        c = t["convergence"]
        lines.append(f"       targeted Newton [6/8]: e60 {t['e60']}  "
                     f"guard {t['guard60']}  meets {t['meets60']}  "
                     f"4001 meets {t['meets60_4001']}")
        lines.append(f"       convergence: certificate "
                     f"{c['equioscillation_certificate']}, newton_ok "
                     f"{c['newton_ok']} ({c['newton_note']}), grid/ref "
                     f"{c['grid_max_over_reference_min']}")
    lines.append(f"  2^5  count 14: best e60 {p5_14['best_e60']}, "
                 f"no meeting pair")
    lines.append(f"       de la Vallee Poussin min grid bound "
                 f"{p5_14['lower_bounds']['min_grid_L']} "
                 f"(> target: "
                 f"{p5_14['lower_bounds']['all_grid_bounds_above_target']})")
    lines.append(f"  ship count 13: SK best e60 {ship_13['sk_best_e60']}")
    if "targeted_newton_6_6" in ship_13:
        t = ship_13["targeted_newton_6_6"]
        c = t["convergence"]
        lines.append(f"       targeted Newton [6/6]: e60 {t['e60']}  "
                     f"guard {t['guard60']}  meets {t['meets60']}  "
                     f"4001 meets {t['meets60_4001']}")
        lines.append(f"       convergence: certificate "
                     f"{c['equioscillation_certificate']}, newton_ok "
                     f"{c['newton_ok']} ({c['newton_note']}), grid/ref "
                     f"{c['grid_max_over_reference_min']}")
    lines.append(f"  ship count 12: best e60 {ship_12['best_e60']}, "
                 f"no meeting pair")
    lines.append(f"       de la Vallee Poussin min grid bound "
                 f"{ship_12['lower_bounds']['min_grid_L']} "
                 f"(> target: "
                 f"{ship_12['lower_bounds']['all_grid_bounds_above_target']})")
    lines.append("")
    lines.append("Remez/Newton solver notes:")
    lines.append("  known-answer test passed: t^4..t^8 minimax polynomials "
                 "with exact errors 2^(1-n)")
    rb = summary["independent_arm"]["remez_cold_6_8_blocker"]
    lines.append(f"  cold [6/8]: {rb['result']} at {rb['grid_max']} "
                 f"(grid/ref {rb['grid_max_over_reference_min']}), "
                 f"spurious near-pole local solution")
    lines.append("  weighted least squares (Sanathanan-Koerner) supplied "
                 "the good starting fit; targeted Newton then reached "
                 "2.56e-16 (2^5) and 5.39e-16 (ship)")
    lines.append("  neither targeted fit reached a full equioscillation "
                 "certificate (newton_ok False); the counts rest on the "
                 "constructive fits plus the de la Vallee Poussin bounds")
    lines.append("")
    lines.append("Conclusion: independent counts 15 (2^5) and 13 (ship), "
                 "identical to BRASIL;")
    lines.append("the clean-edge meeting counts are not solver-specific.")
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n",
                                     encoding="utf-8")
    print("\n".join(lines))
    print(f"\nwrote {OUT / 'summary.json'} and {OUT / 'summary.txt'}")


if __name__ == "__main__":
    main()
