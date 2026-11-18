#!/usr/bin/env python3
"""Extract the two paper-cell records that the certificate scripts consume.

Inputs (read-only, Paper C's results-corrected records):
  results-corrected/corrected-counts-60dig.json   (the 45-cell 60-digit audit)
  results-corrected/corrected-counts.json         (the corrected BRASIL record)

The source digests are the ones the Paper C supplement names.  The script
writes one-cell documents so the deposit is self-contained; the deposited
certificates were generated from the full records, and only the recorded
input path/digest fields differ when regenerated from these copies.

Usage (from the deposit root):
  python -I inputs/extract_cells.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = Path(r"C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts")
AUDIT = PAPER / "results-corrected" / "corrected-counts-60dig.json"
CORRECTED = PAPER / "results-corrected" / "corrected-counts.json"
AUDIT_SHA256 = "18a41a42798f0cd1fead4a143309ded3112d1bba94727b093bf27df02e3256e7"
CORRECTED_SHA256 = (
    "9bb8ac0633b10e8a7922684ca8a419135103f2c7825b4d06b9c71603399b7195")
TARGET = "abs:5e-14"
WIDTH_NAME = "paper"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_cell(doc, width_name, target):
    for cell in doc["cells"]:
        if cell.get("width_name") == width_name and cell.get("target") == target:
            return cell
    raise SystemExit(f"cell ({width_name}, {target}) not found")


def main():
    digest = sha256(AUDIT)
    if digest != AUDIT_SHA256:
        raise SystemExit(f"audit digest {digest} != {AUDIT_SHA256}")
    digest = sha256(CORRECTED)
    if digest != CORRECTED_SHA256:
        raise SystemExit(f"corrected digest {digest} != {CORRECTED_SHA256}")
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    corrected = json.loads(CORRECTED.read_text(encoding="utf-8"))
    cells = {
        "paper_cell_abs_audit.json": dict(
            schema="paper-c-optimality-input-cell-1",
            meta=dict(source=str(AUDIT), sha256=AUDIT_SHA256,
                      width_name=WIDTH_NAME, target=TARGET),
            cells=[find_cell(audit, WIDTH_NAME, TARGET)],
        ),
        "paper_cell_abs_corrected.json": dict(
            schema="paper-c-optimality-input-cell-1",
            meta=dict(source=str(CORRECTED), sha256=CORRECTED_SHA256,
                      width_name=WIDTH_NAME, target=TARGET),
            cells=[find_cell(corrected, WIDTH_NAME, TARGET)],
        ),
    }
    for name, doc in cells.items():
        path = HERE / name
        path.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
        print(f"wrote {path} sha256={sha256(path)}")


if __name__ == "__main__":
    main()
