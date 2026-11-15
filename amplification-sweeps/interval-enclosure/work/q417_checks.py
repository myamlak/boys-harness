"""Assemble the Q-417 checks and the record files from the four case runs.

Reads work/<case>.{json,csv} written by q417_interval.py and the deposited
charge record under
amplification-sweeps/certified-bound/results-underflow-fix/.  Writes
results.csv, results.json and outcome_summary.json beside this record's
README.  The repository root is derived from this file's location; set
BOYS_HARNESS to override it.
"""
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True

WORK = Path(__file__).resolve().parent
RECORD = WORK.parent
HARNESS = Path(os.environ.get(
    "BOYS_HARNESS", Path(__file__).resolve().parents[3]))
CASES = ("first9353", "first9352", "second", "heldout")
CERT = HARNESS / "amplification-sweeps/certified-bound/results-underflow-fix"


def load_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def case_checks(case):
    rows = load_csv(WORK / (case + ".csv"))
    summary = json.load(open(WORK / (case + ".json"), encoding="utf-8"))
    dep = [r for r in rows if r["dep_state"] != "abstain"]
    enc = [r for r in rows if r["enc_state"] != "abstain"]
    dep_by_state = {s: [r for r in dep if r["dep_state"] == s]
                    for s in ("recursed", "direct")}
    enc_by_state = {s: [r for r in enc if r["enc_state"] == s]
                    for s in ("recursed", "direct")}
    order0 = [r for r in rows if r["order"] == "0"]
    model = [r for r in rows if r["order"] != "0"]

    def split(part):
        return {"rows": len(part),
                "commits": sum(1 for r in part if r["dep_state"] != "abstain"),
                "recursed": sum(1 for r in part if r["dep_state"] == "recursed"),
                "direct": sum(1 for r in part if r["dep_state"] == "direct")}

    justified = [r for r in dep if r["dep_justified_by_enc"] == "1"]
    unjustified = [r for r in dep if r["dep_justified_by_enc"] == "0"]
    wrong_whole = [r for r in dep if r["dep_right_whole"] == "0"]
    wrong_pop = [r for r in dep if r["dep_right_pop"] == "0"]
    dep_winner_whole = {
        "committed_rows": len(dep),
        "winner_direct": sum(1 for r in dep if r["winner_whole"] == "direct"),
        "winner_recursed": sum(1 for r in dep if r["winner_whole"] == "recursed"),
        "winner_tie": sum(1 for r in dep if r["winner_whole"] == "tie"),
    }
    enc_wrong_whole = [r for r in enc if r["enc_right_whole"] == "0"]
    enc_extra = [r for r in rows if r["enc_commits_dep_abstains"] == "1"]

    always_direct_whole = sum(1 for r in rows if r["winner_whole"] == "direct")
    always_direct_pop = sum(1 for r in rows if r["winner_pop"] == "direct")
    always_recursed_whole = sum(1 for r in rows if r["winner_whole"] == "recursed")
    always_recursed_pop = sum(1 for r in rows if r["winner_pop"] == "recursed")
    ties_whole = sum(1 for r in rows if r["winner_whole"] == "tie")
    ties_pop = sum(1 for r in rows if r["winner_pop"] == "tie")

    served_zero = [r for r in rows if int(r["served_cells"]) == 0]
    structural = {
        "rows": len(served_zero),
        "orders": sorted({(r["band"], int(r["order"])) for r in served_zero}),
        "all_abstain": all(r["dep_state"] == "abstain" for r in served_zero),
        "W_equals_D": all(r["W"] == r["D"] for r in served_zero),
    }
    if case == "heldout":
        structural["abstentions_total"] = sum(1 for r in rows
                                              if r["dep_state"] == "abstain")
        structural["abstentions_structural"] = sum(
            1 for r in served_zero if r["dep_state"] == "abstain")

    shortfall = []
    for r in rows:
        wmt = float(r["W_minus_T"])
        low = float(r["L"])
        if wmt < low:
            shortfall.append(((low - wmt) / low if low else None,
                              r["band"], r["metric"], int(r["order"]),
                              r["W_minus_T"], r["L"]))
    shortfall.sort(reverse=True)
    shortfall_max = shortfall[0][0] if shortfall else None

    return {
        "case": case,
        "counts": summary["counts"],
        "state_agreement": summary["state_agreement"],
        "max_abs_U_minus_W": summary["max_abs_U_minus_W"],
        "deposited_split": {
            "order0": split(order0), "model": split(model),
        },
        "deposited_by_state": {s: len(v) for s, v in dep_by_state.items()},
        "enclosure_by_state": {s: len(v) for s, v in enc_by_state.items()},
        "deposited_winner_whole": dep_winner_whole,
        "justified": [(r["band"], r["metric"], int(r["order"]), r["dep_state"])
                      for r in justified],
        "unjustified_count": len(unjustified),
        "unjustified_by_state": {
            s: sum(1 for r in unjustified if r["dep_state"] == s)
            for s in ("recursed", "direct")},
        "wrong_whole": [(r["band"], r["metric"], int(r["order"]),
                         r["dep_state"], r["winner_whole"])
                        for r in wrong_whole],
        "wrong_pop": [(r["band"], r["metric"], int(r["order"]),
                       r["dep_state"], r["winner_pop"])
                      for r in wrong_pop],
        "enclosure_wrong_whole": [(r["band"], r["metric"], int(r["order"]),
                                   r["enc_state"], r["winner_whole"])
                                  for r in enc_wrong_whole],
        "enclosure_extra": [(r["band"], r["metric"], int(r["order"]),
                             r["enc_state"]) for r in enc_extra],
        "measured_whole": {"direct": always_direct_whole,
                           "recursed": always_recursed_whole,
                           "tie": ties_whole},
        "measured_pop": {"direct": always_direct_pop,
                         "recursed": always_recursed_pop,
                         "tie": ties_pop},
        "structural": structural,
        "W_minus_T_below_L": {
            "rows": len(shortfall),
            "max_relative_shortfall": shortfall_max,
            "worst": shortfall[:5],
        },
        "population_score_pair": {
            "whole": summary["counts"]["deposited_whole"],
            "restricted": summary["counts"]["deposited_pop"],
        },
    }


def charge_check():
    log = CERT / "harness.log"
    cells = CERT / "harness_cells.csv"
    rows = load_csv(cells)
    from mpmath import mp, mpf
    mp.dps = 70
    violations = []
    for row in rows:
        if mpf(row["r_ship"]) > mpf(row["paper_ship"]):
            violations.append((row["band"], row["metric"], int(row["n"]), "ship"))
        if mpf(row["r_prod"]) > mpf(row["paper_prod"]):
            violations.append((row["band"], row["metric"], int(row["n"]), "product"))
    text = log.read_text(encoding="utf-8")
    line = [ln for ln in text.splitlines() if ln.startswith("charge violations")]
    return {
        "record": str(cells),
        "sha256": sha256(cells),
        "log": str(log),
        "log_sha256": sha256(log),
        "log_line": line,
        "rows_commitments": len(rows),
        "entry_cells": 2 * len(rows),
        "violations": len(violations),
        "violating_cells": violations,
    }


def inputs_block():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from q417_interval import CASES, PAPER_ANALYZE
    cache = {}

    def digest(path):
        path = Path(path)
        if path not in cache:
            cache[path] = {"path": str(path), "sha256": sha256(path)}
        return cache[path]

    block = {"paper_model_analyze_py": digest(PAPER_ANALYZE)}
    for case, spec in CASES.items():
        block[case] = {key: digest(spec[key])
                       for key in ("cells", "reference", "provenance",
                                   "per_order", "record")
                       if spec.get(key) is not None}
    block["model_record_prediction_csv"] = digest(
        HARNESS / "amplification-sweeps/certified-bound/inputs/model-record/prediction.csv")
    block["certified_bound_readme"] = digest(
        HARNESS / "amplification-sweeps/certified-bound/README.md")
    block["heldout_readme"] = digest(
        HARNESS / "amplification-sweeps/held-out/README.md")
    block["exact_tag_replay_readme"] = digest(
        HARNESS / "boys-crossover/exact-tag-replay/README.md")
    block["older_revision_analyze_py"] = digest(
        HARNESS / "amplification-sweeps/amplification-form/analyze.py")
    return block


def main():
    RECORD.mkdir(parents=True, exist_ok=True)
    checks = {case: case_checks(case) for case in CASES}
    charge = charge_check()
    containment = {}
    for case in ("first9353", "second", "heldout"):
        path = WORK / (case + ".containment.json")
        if path.is_file():
            containment[case] = json.load(open(path, encoding="utf-8"))

    second = checks["second"]
    held = checks["heldout"]
    first = checks["first9353"]
    five_checks = {
        "first_control_pair": {
            "expected": "68 commitments, 68 right restricted, 67 of 68 whole-band",
            "first9353": first["population_score_pair"],
            "first9352": checks["first9352"]["population_score_pair"],
            "confirmed": (
                first["population_score_pair"]["restricted"]["commits"] == 68
                and first["population_score_pair"]["whole"]["right"] == 67
                and first["population_score_pair"]["whole"]["wrong"] == 1
                and first["population_score_pair"]["restricted"]["right"] == 68
                and first["population_score_pair"]["restricted"]["wrong"] == 0),
        },
        "second_decomposition": {
            "expected": "65 = 4 order-0 measured + 61 model (57 direct, 4 recursed)",
            "order0": second["deposited_split"]["order0"],
            "model": second["deposited_split"]["model"],
            "confirmed": (
                second["deposited_split"]["order0"] == {
                    "rows": 4, "commits": 4, "recursed": 2, "direct": 2}
                and second["deposited_split"]["model"] == {
                    "rows": 128, "commits": 61, "recursed": 4, "direct": 57}),
        },
        "heldout_decomposition": {
            "expected": "58 = 53 direct + 5 recursed; always-direct right on 53 of 58",
            "deposited_by_state": held["deposited_by_state"],
            "commit_winners_whole": held["deposited_winner_whole"],
            "recursed_commitments": [
                r for r in held["justified"] if r[3] == "recursed"],
            "confirmed": (
                held["deposited_by_state"] == {"recursed": 5, "direct": 53}
                and held["deposited_winner_whole"]["winner_direct"] == 53
                and held["deposited_winner_whole"]["winner_recursed"] == 5
                and held["wrong_whole"] == []
                and sum(1 for r in held["justified"] if r[3] == "recursed") == 5),
        },
        "heldout_structural_abstentions": {
            "expected": "32 of 74 abstentions structural, orders 17-32, same stored fit",
            "structural": held["structural"],
            "confirmed": (
                held["structural"]["rows"] == 32
                and held["structural"]["all_abstain"]
                and held["structural"]["W_equals_D"]
                and held["structural"]["abstentions_structural"] == 32
                and held["structural"]["abstentions_total"] == 74),
        },
        "heuristic_charge": {
            "expected": "paper charge exceeded at 5 of 130 harness entry-cells",
            "confirmed": (charge["violations"] == 5
                          and charge["entry_cells"] == 130),
            "detail": charge,
        },
    }

    all_rows = []
    columns = ["experiment", "band", "metric", "order",
               "population_arguments", "population_cells", "served_cells",
               "unserved_cells", "W", "D", "T", "W_minus_T", "L", "U", "Dpop",
               "dep_state", "dep_choice", "enc_state", "winner_whole",
               "winner_pop", "dep_right_whole", "dep_right_pop",
               "enc_right_whole", "enc_right_pop", "dep_justified_by_enc",
               "enc_commits_dep_abstains"]
    for case in CASES:
        all_rows.extend(load_csv(WORK / (case + ".csv")))
    with open(RECORD / "results.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in all_rows:
            writer.writerow({column: row[column] for column in columns})

    payload = {
        "checks": checks,
        "five_checks": five_checks,
        "charge_check": charge,
        "containment": containment,
        "inputs": inputs_block(),
        "cases": {case: json.load(open(WORK / (case + ".json"), encoding="utf-8"))
                  for case in CASES},
    }
    with open(RECORD / "results.json", "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1, default=str)

    compact = {case: {
        "population_score_pair": checks[case]["population_score_pair"],
        "counts": checks[case]["counts"],
        "deposited_split": checks[case]["deposited_split"],
        "deposited_winner_whole": checks[case]["deposited_winner_whole"],
        "justified": checks[case]["justified"],
        "unjustified_count": checks[case]["unjustified_count"],
        "unjustified_by_state": checks[case]["unjustified_by_state"],
        "worst_whole": checks[case]["wrong_whole"],
        "enclosure_wrong_whole": checks[case]["enclosure_wrong_whole"],
        "enclosure_extra": checks[case]["enclosure_extra"],
        "measured_whole": checks[case]["measured_whole"],
        "structural": checks[case]["structural"],
        "W_minus_T_below_L": checks[case]["W_minus_T_below_L"],
    } for case in CASES}
    with open(RECORD / "outcome_summary.json", "w", encoding="utf-8") as handle:
        json.dump({"cases": compact,
                   "five_checks_confirmed": {
                       name: entry["confirmed"]
                       for name, entry in five_checks.items()},
                   "charge_check": charge}, handle, indent=1, default=str)
    print(json.dumps({"cases": compact, "charge": charge}, indent=1, default=str))


if __name__ == "__main__":
    main()
