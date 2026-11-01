#!/usr/bin/env python3
"""Build the corrected F1 record from the two knee walks.

The F1 sweep recorded two rational rows whose below-count check also reached
the target (below_meets true).  The paper's rule makes the stored count the
smallest stored total that reaches the target on the certifying grid, so a row
whose below count also reaches is not minimal.  knee_walk_f1.py solved every
level below each of those rows with count_sweep.refine_pairs at the recorded
settings and stopped at the first level that misses.

This script replaces each flagged row's rational entry with the corrected
level's best meeting fit, clears below_meets, updates saving and parity, and
records the original values under a top-level corrections block.  Everything
else in the record is copied verbatim.  The output is a derived record, not a
replacement for the sweep record.

Usage:
    python make_corrected_record.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
RECORD = PAPER / "results-second-instance-f1" / "sweep_all.json"
OUT_DIR = PAPER / "_knee-f1" / "corrected"
OUT = OUT_DIR / "sweep_all.json"

WALKS = {
    ("a2", "rel:1e-13", "paper"): PAPER / "_knee-f1" / "knee_walk_a2-paper-rel13.json",
    ("a2", "rel:1e-14", "2^5"): PAPER / "_knee-f1" / "knee_walk_a2-2^5-rel14.json",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    doc = json.loads(RECORD.read_text(encoding="utf-8"))
    corrections = []
    for key, walk_path in WALKS.items():
        row = next(r for r in doc["rows"]
                   if (r["edge"], r["target"], r["width_name"]) == key)
        walk = json.loads(walk_path.read_text(encoding="utf-8"))
        corr = walk.get("corrected") or {}
        if corr.get("stored") is None:
            raise SystemExit(f"{walk_path}: walk found no meeting level")
        if corr.get("first_missing_level") != corr["stored"] - 1:
            raise SystemExit(
                f"{walk_path}: unexpected walk shape {corr}")
        if not walk.get("generator_digest_matches_record"):
            raise SystemExit(f"{walk_path}: generator digest mismatch")
        best = corr["best"]
        entry = row["rat"]
        original = {k: entry.get(k) for k in
                    ("stored", "m", "k", "solver", "delivered_rel",
                     "delivered_abs", "at_x_rel", "at_x_abs", "seconds",
                     "below_meets")}
        level = next(e for e in walk["levels"]
                     if e["stored"] == corr["stored"])
        entry.update({
            "stored": corr["stored"],
            "m": best["m"],
            "k": best["k"],
            "solver": best["solver"],
            "delivered_rel": best["rel"],
            "delivered_abs": best["abs"],
            "at_x_rel": best["at_x_rel"],
            "at_x_abs": best["at_x_abs"],
            "refined": True,
            "below_meets": False,
            "seconds": level["seconds"],
            "seconds_note": ("level total for eight pairs on four workers; "
                             "not a single-pair solver time"),
        })
        row["saving"] = row["cheb"]["stored"] - entry["stored"]
        if "parity" in row:
            row["parity"]["rat_m"] = entry["m"]
            row["parity"]["rat_k"] = entry["k"]
            row["parity"]["rat_m_parity"] = ("even" if entry["m"] % 2 == 0
                                             else "odd")
            row["parity"]["rat_k_parity"] = ("even" if entry["k"] % 2 == 0
                                             else "odd")
        corrections.append({
            "cell": walk["cell"],
            "edge": key[0], "target": key[1], "width_name": key[2],
            "walk": str(walk_path), "walk_digest": sha(walk_path),
            "recorded": original,
            "corrected": {
                "stored": entry["stored"], "m": entry["m"], "k": entry["k"],
                "solver": entry["solver"],
                "delivered_rel": entry["delivered_rel"],
                "delivered_abs": entry["delivered_abs"],
                "at_x_rel": entry["at_x_rel"], "at_x_abs": entry["at_x_abs"],
                "below_meets": False,
            },
            "first_missing_level": corr["first_missing_level"],
            "moved": corr["moved"],
            "walk_seconds": walk.get("walk_seconds"),
        })
    doc["corrections"] = {
        "note": ("Derived record: the two rational rows with below_meets true "
                 "replaced by the corrected minimal stored count from "
                 "knee_walk_f1.py at the recorded settings.  All other rows, "
                 "meta and settings are verbatim from the sweep record."),
        "source_record": str(RECORD),
        "source_record_digest": sha(RECORD),
        "entries": corrections,
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print(f"wrote {OUT} ({sha(OUT)})")
    for c in corrections:
        print(f"{c['cell']}: recorded {c['recorded']['stored']} "
              f"[{c['recorded']['m']}/{c['recorded']['k']}] "
              f"{c['recorded']['solver']} -> corrected {c['corrected']['stored']} "
              f"[{c['corrected']['m']}/{c['corrected']['k']}] "
              f"{c['corrected']['solver']} (moved {c['moved']}); "
              f"level {c['first_missing_level']} misses")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
