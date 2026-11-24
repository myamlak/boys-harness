# The qcx integration bridge

The adapter through which the qcx integral engine calls the Boys library, and
the two compile-time switches the integral-engine record was run with. The
probes and their build stanzas are deposited beside this directory:
`../scripts/` holds the water probe and `../benzene-pvdz/` the benzene probe.
This directory carries the engine-side half those probes are built against.

The files are taken from the qcx framework at commit
`ac8540b7c23d3818c9e67f8d59bf6cfbbef7fa0d` (v0.1.0-36, project version 0.2.0);
the adapter and the module build file carry the bridge's two switches. They
mirror their paths in that tree and are applied to a checkout of that revision.
The engine tree itself is not deposited.

* `integrals/src/internal/md_boys.hpp` is the engine's Boys adapter, with the
  route switch and the argument capture used by the water and cc-pVTZ runs.
* `benzene/integrals/src/internal/md_boys.hpp` is the capture variant used by
  the benzene run: it keeps one record per distinct argument instead of one per
  call.
* `integrals/include/qcx/integrals/boys.hpp` and `f16.hpp` are the front door
  the adapter includes; `integrals/src/boys_build_record.cpp` is the record
  translation unit the probe writes.
* `integrals/CMakeLists.txt` compiles the adapter and defines the two switches
  on `qcx-boys-route`.

The switches are `QCX_BRIDGE_BOYS_POLICY` (`default` or `paper-a`), which
selects the route the engine's all-orders Boys call evaluates - the library's
host fp64 class row, or the route-comparison study's named recursed policy -
and `QCX_BRIDGE_BOYS_CAPTURE`, which records the argument and ladder of every
engine call so the probe can evaluate both routes on the run's own stream.

The engine's root build file adds the probe target:

    option(QCX_BUILD_BRIDGE_PROBE "Build the chemistry bridge probe" OFF)
    ...
    if (QCX_BUILD_BRIDGE_PROBE)
        add_subdirectory(bridge)
    endif()

The dense tensors are not deposited; the probes rebuild them.

## Build

Place the probe and its build stanza from `../scripts/` (or
`../benzene-pvdz/` for the benzene record) at `bridge/probe.cpp` and
`bridge/CMakeLists.txt`, then configure and build one arm:

    cmake -S <qcx-checkout> -B build-bridge -DCMAKE_BUILD_TYPE=Release \
        -DQCX_BUILD_TESTS=OFF -DQCX_BUILD_BENCHMARKS=OFF -DQCX_ENABLE_IO=OFF \
        -DQCX_INTEGRALS_LMAX=3 -DQCX_BUILD_BRIDGE_PROBE=ON \
        -DQCX_BRIDGE_BOYS_CAPTURE=ON -DQCX_BRIDGE_BOYS_POLICY=default
    cmake --build build-bridge --config Release --target qcx-bridge-probe

The second arm is the same configure in a second build directory with
`-DQCX_BRIDGE_BOYS_POLICY=paper-a`.

## Run

    qcx-bridge-probe.exe <out-prefix> data/basis/cc-pvtz

The benzene record runs the same target with the molecule named:

    qcx-bridge-probe.exe <out-prefix> data/basis/cc-pvdz benzene

`../scripts/commands.md` gives the exact settings and exit codes for the water
runs; `../benzene-pvdz/commands.md` gives the benzene run and its comparison.
The water records were measured against Boys 5.0.0 and the benzene record
against Boys 5.2.0; `../provenance.json` is the record.
