# The v5.3.1 release replay

This directory is the check that carries Paper A's measurements to the
submission release. It rebuilds both sweep programs against the release's
headers, runs them on the deposited grids and fits, and compares every
produced record with the deposited one. It also re-runs the certified-bound
record's model replay and its full bound validation at the release.

The harness tag of the submission is `v1.0.34`.

## What it proves

- Both sweeps, compiled against the release's headers, reproduce the deposited
  cells bit for bit: 308,616 first-sweep rows (617,232 entry values) and
  27,588 second-sweep rows (55,176 entry values), with zero mismatches. The
  first sweep's admission and replica-gain records and the second sweep's
  admission and path records are byte-identical as well.
- The certificate's model replay re-derives all 27,588 recursed and 27,588
  direct values bit for bit with the separate downward step. The fused
  downward variant mismatches the same 208 recursed values as the recorded
  run, the first at (A1, 0.029916890803630927, 0), one ulp.
- The bound validation reports 0 violations at all 55,176 entry-cells, with
  the recorded worst measured-to-certified ratio 0.71875 at (A1, 0, 24, ship).
  Every certificate output file is byte-identical to the deposited run.

The source-identity files beside this record, `digests.txt` and
`identity-v531.txt`, state what remains: 45 of the 62 `include/boys` headers,
10 of the 14 `src` files and 55 of the 76 evaluated files are byte-identical
to the measured tree. The applicability of the measurements to the release
therefore rests on execution, not on source identity: the release executes the
measured values at every deposited cell.

## Reproducing

The tag was not cut when this replay ran. The build used the release tree
whose `include/` and `src/` content differs from the tagged tree only in
`include/boys/version.hpp`; `digests.txt` records the version header at both
identities and every other evaluated file once, and `provenance.json` records
the tag objects. Re-run against the tag after it is cut: a difference is a
failed check, not a number to adjust.

The two sweep programs are the paper's `sweeps/product_amplification/sweep1.cpp`
and `sweeps/amplification_form/sweep2.cpp` with `sweeps/common/replica.hpp`.
The compile line is the runbook's own:

    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 \
       /I <release>/include /Fobuild\sweep1.obj /Febuild\sweep1.exe sweep1.cpp
    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 \
       /I <release>/include /Fobuild\sweep2.obj /Febuild\sweep2.exe sweep2.cpp

run on the deposited inputs:

    build\sweep1.exe --grid data/grid.csv --cells results/cells.csv \
        --admission results/admission.csv --gain results/replica_gain.csv \
        --batch full
    build\sweep2.exe --grid data/grid.csv --fits data/fits.csv \
        --paths results/paths.csv --cells results/cells.csv \
        --admission results/admission.csv

The certificate record's own commands, from the certificate directory:

    python -I code/replay.py --kernel <release> --both-downward
    python -I code/certify.py --kernel <release> --sweep ../amplification-form \
        --model-record inputs/model-record --exp-bits exp/exp_bits.txt --out <run>

All four programs exited 0. `sweep-replay-v531.log` and `replay-v531.log` are
the captured runs, with the byte comparison of every produced record.

## Files

| File | What it records |
| --- | --- |
| `README.md` | this file |
| `digests.txt` | sha256 of the 76 evaluated files at the three identities, and the counts |
| `identity-v531.txt` | the same counts as prose, one block per identity |
| `sweep-replay-v531.log` | the two builds, the two runs, and the byte comparison of every sweep record |
| `replay-v531.log` | the certificate's model replay, the full bound validation, and the output comparison |
| `provenance.json` | machine-readable identity: tag objects, input digests, run exits and results |
