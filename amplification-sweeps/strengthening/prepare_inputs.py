"""Decompress the input streams the strengthening scripts read.

    python prepare_inputs.py

The deposit ships its three sizeable inputs compressed with xz.  The analysis
scripts read plain CSV from inputs/, so run this once first.  It decompresses,
in a streaming loop,

    product-amplification/results/cells.csv.xz   -> inputs/cells.csv
    product-amplification/data/reference.csv.xz  -> inputs/reference.csv
    amplification-form/results/paths.csv.xz      -> inputs/paths.csv

and checks each output against the sha256 recorded for the original
decompression in verify_provenance.log:1-3.  The .xz streams are opened
read-only and are never modified.  inputs/ is generated and is not part of the
deposit; re-running this script rewrites it from the deposit's streams.
"""

import hashlib
import lzma
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE.parent
INPUTS = HERE / "inputs"

# (deposit source relative to DEPOSIT, target name, decompressed sha256)
STREAMS = (
    (
        "product-amplification/results/cells.csv.xz",
        "cells.csv",
        "dd26d0b2a4352113d9f15433ec991adfacc211607ae9586c0bb51cd12b3f3369",
    ),
    (
        "product-amplification/data/reference.csv.xz",
        "reference.csv",
        "1e554a3540c6a15535b80f7f2597360bfc3a75919fa71d1496b284f348deb5b9",
    ),
    (
        "amplification-form/results/paths.csv.xz",
        "paths.csv",
        "81001a06d80c88e5b3ddac1f3dd011cae6f8ce2e806c8bf1f7f26e85dd8c8219",
    ),
)


def prepare(source, name, expected):
    packed = DEPOSIT / source
    target = INPUTS / name
    temporary = INPUTS / (name + ".tmp")
    digest = hashlib.sha256()
    size = 0
    with lzma.open(packed, "rb") as reader, open(temporary, "wb") as writer:
        while True:
            block = reader.read(1 << 22)
            if not block:
                break
            writer.write(block)
            digest.update(block)
            size += len(block)
    actual = digest.hexdigest()
    if actual != expected:
        temporary.unlink()
        print("FAIL inputs/%s: sha256 %s, expected %s" % (name, actual, expected))
        return False
    os.replace(temporary, target)
    print("prepared inputs/%-13s %10d bytes  sha256 %s" % (name, size, actual))
    return True


def main():
    INPUTS.mkdir(exist_ok=True)
    ok = True
    for source, name, expected in STREAMS:
        packed = DEPOSIT / source
        if not packed.is_file():
            print("FAIL missing deposit stream %s" % source)
            ok = False
            continue
        ok = prepare(source, name, expected) and ok
    if not ok:
        print("FAILED")
        return 1
    print("all three inputs prepared under inputs/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
