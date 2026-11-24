"""D84 fixture: the corrected 84-run pooled design on the six-axis cube.

The pooled cube has axes route(2), scheme(2), partition(3), pack(2),
form(3), exp(2).  The corrected packed-sum cone is
{route, scheme, partition, pack}, so the permitted set is the 72
pack-free channels of the upward cone plus the 24 packed-cone channels
less the 12 channels in both cones, 84 in all.  The complement has 60
structural zeros and the targets are the 75 permitted channels of order
2 to 5.

The design D84 is the 72-cell five-axis factorial at pack = arguments
plus the 12 cells at pack = orders with route, scheme in {0,1},
partition in {0,1,2} and form = exp = 0.

This module returns the matrices and constants the prospective checks
need and verifies them against the deposited D84 certificate.  It reads
the basis and cell conventions of the D75 fixture and the full-cube
loader; it does not regenerate the measured record.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def harness_root():
    """Locate the harness repository holding the composition audit."""
    override = os.environ.get("BOYS_HARNESS_REPO")
    if override:
        return os.path.abspath(override)
    root = os.path.abspath(os.path.join(HERE, "..", ".."))
    probe = os.path.join(root, "composition-audit", "fullcube",
                         "fullcube_decomposition.py")
    if os.path.isfile(probe):
        return root
    raise SystemExit("set BOYS_HARNESS_REPO to the harness repository root")


REPO = harness_root()
D84_DIR = os.path.join(REPO, "composition-audit",
                       "support-family-d84-2026-11-24")
D75_DIR = os.path.join(REPO, "composition-audit",
                       "support-family-2026-11-22")
sys.path.insert(0, D75_DIR)

import corrected_frontier as cf          # noqa: E402
import instance as inst                  # noqa: E402
import pooled_allowance as pa            # noqa: E402
import pooled_d75 as d75                 # noqa: E402

fc = cf.fc

# The D84 design cells, in the certificate's order.
D84 = np.array(sorted(
    d75.five_axis_cells(0)
    + [d75.cell_index(r, s, p, 1, 0, 0)
       for r in range(2) for s in range(2) for p in range(3)]), dtype=np.int64)

AGAINST_D84 = "certify_d84.json"


def build(design="d84"):
    """Return the pooled design's matrices and constants.

    design is "d84" (the corrected cone) or "d75" (the earlier record's
    cone, used only to reproduce the paper's D75 numbers).
    """
    if design == "d84":
        packed = {0, 1, 2, 3}          # route, scheme, partition, pack
        cells = D84
        kappa_joint = 4.0
    elif design == "d75":
        packed = {2, 3}                # partition, pack
        cells = np.array(sorted(
            d75.five_axis_cells(0)
            + d75.extra_cells([0, 1, 2], base=(0, 0, 0, 0))),
            dtype=np.int64)
        kappa_joint = 3.0
    else:
        raise ValueError(design)

    # Correct the packed cone in the shared allowance module.  Nothing
    # below reads the module globals again except through this dict.
    pa.PACKED = set(packed)
    perm = np.array([pa.in_pooled(s) for s in pa.SUPPORTS])
    q = ~perm
    target = np.array([perm[i] and 2 <= pa.ORDERS[i] <= pa.W_P
                       for i in range(pa.N)])
    P_IDX = np.flatnonzero(perm)
    Q_IDX = np.flatnonzero(q)
    T_IDX = np.flatnonzero(target)

    phi = pa.PHI
    F = phi[cells, :]
    PP = F[:, P_IDX]
    FQ = F[:, Q_IDX]
    Om = PP.T @ PP / len(cells)
    inv = np.linalg.inv(Om)
    PT = phi[:, T_IDX]
    PQ = phi[:, Q_IDX]
    posP = {int(p): i for i, p in enumerate(P_IDX)}
    posT = np.array([posP[int(t)] for t in T_IDX])
    G = pa.G_ALPHA[P_IDX]
    GQ = pa.G_ALPHA[Q_IDX]
    GT = pa.G_ALPHA[T_IDX]

    B = inv @ (PP.T @ FQ / len(cells))
    kappa = np.abs(B / GQ[None, :])
    kappa_target = float(kappa[posT].max())
    g_kappa = G[posT] * kappa[posT].max(axis=1)
    amplification = float((np.abs(inv) @ G).max())

    mi = cf.measurement_inflation(cells, P_IDX, T_IDX, phi, inv,
                                  pa.G_ALPHA)

    # Numeric joint constant: sup over the full 144-cell cube.
    residual = PQ - PT @ B[posT, :]
    ratio = np.abs(residual) / GQ[None, :]
    kappa_joint_values = ratio.max(axis=0)
    kappa_joint_numeric = float(kappa_joint_values.max())

    return {
        "design": design,
        "cells": cells,
        "P_IDX": P_IDX, "Q_IDX": Q_IDX, "T_IDX": T_IDX,
        "G": G, "GQ": GQ, "GT": GT,
        "PP": PP, "inv": inv, "B": B,
        "PT": PT, "PQ": PQ, "posT": posT,
        "kappa_joint": kappa_joint,
        "kappa_target": kappa_target,
        "g_kappa": g_kappa,
        "amplification": amplification,
        "measurement_inflation": mi,
        "kappa_joint_numeric": kappa_joint_numeric,
        "kappa_joint_values": kappa_joint_values,
    }


def verify_against_certificate(d, path=None):
    """Compare a d84 build with certify_d84.json; raise on mismatch."""
    path = path or os.path.join(D84_DIR, AGAINST_D84)
    with open(path) as handle:
        cert = json.load(handle)
    checks = {}
    checks["cells"] = bool(np.array_equal(d["cells"], cert["cells"]))
    checks["permitted"] = len(d["P_IDX"]) == cert["permitted"]
    checks["complement"] = len(d["Q_IDX"]) == cert["complement"]
    checks["targets"] = len(d["T_IDX"]) == cert["targets"]
    checks["rank"] = int(np.linalg.matrix_rank(d["PP"])) == cert["rank"]
    checks["kappa_target"] = abs(d["kappa_target"] - cert["kappa_target"]) < 1e-12
    checks["sum_G_kappa"] = abs(float(d["g_kappa"].sum())
                                - cert["sum_G_kappa"]) < 1e-9
    checks["G_kappa"] = bool(np.allclose(d["g_kappa"], cert["G_kappa"],
                                        rtol=0, atol=1e-12))
    checks["amplification"] = abs(d["amplification"]
                                  - cert["amplification"]) < 1e-9
    checks["kappa_joint"] = abs(d["kappa_joint_numeric"]
                                - cert["kappa_joint"]) < 1e-9
    values = sorted({round(float(v), 6)
                     for v in d["kappa_joint_values"]})
    checks["kappa_joint_values"] = values == [
        float(v) for v in cert["kappa_joint_values"]]
    checks["measurement_max_target"] = abs(
        d["measurement_inflation"]["max_targets"]
        - cert["measurement_max_target"]) < 1e-9
    checks["measurement_g_sum"] = abs(
        d["measurement_inflation"]["sum_G_amp_targets"]
        - cert["measurement_g_sum"]) < 1e-9
    failed = [k for k, ok in checks.items() if not ok]
    if failed:
        raise SystemExit(f"certificate mismatch: {failed}")
    return checks


def main():
    print("harness repository resolved from BOYS_HARNESS_REPO or the "
          "script's parent checkout")

    def plain(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, dict):
            return {k: plain(v) for k, v in value.items()}
        return value

    out = {}
    for design in ("d75", "d84"):
        d = build(design)
        out[design] = {
            k: plain(v) for k, v in d.items()
            if k not in ("PP", "inv", "B", "PT", "PQ", "cells",
                         "P_IDX", "Q_IDX", "T_IDX")
        }
        out[design]["cells"] = [int(c) for c in d["cells"]]
        out[design]["counts"] = {
            "cells": int(len(d["cells"])),
            "permitted": int(len(d["P_IDX"])),
            "complement": int(len(d["Q_IDX"])),
            "targets": int(len(d["T_IDX"])),
        }
        print(f"{design}: cells {len(d['cells'])}, permitted "
              f"{len(d['P_IDX'])}, complement {len(d['Q_IDX'])}, targets "
              f"{len(d['T_IDX'])}, kappa_joint {d['kappa_joint']} "
              f"(numeric {d['kappa_joint_numeric']:.15f}), "
              f"kappa_target {d['kappa_target']:.15f}, "
              f"amplification {d['amplification']:.12f}")
        print(f"  kappa_joint values "
              f"{sorted({round(float(v), 6) for v in d['kappa_joint_values']})}")
        print(f"  measurement {d['measurement_inflation']}")
        if design == "d84":
            checks = verify_against_certificate(d)
            print("  certificate checks:", "all pass"
                  if all(checks.values()) else checks)
            # The D75 and D84 builds are different designs; compare only
            # the common matrices that must be identical.
            d75 = build("d75")
            agree = {
                "D84 vs certificate cells": bool(np.array_equal(
                    d["cells"], D84)),
                "D84 cells distinct": len(set(d["cells"].tolist())) == 84,
                "D75 cells match corrected_frontier": bool(np.array_equal(
                    d75["cells"], np.array(cf.D75))),
                "D75 inv matches corrected_frontier": bool(np.allclose(
                    d75["inv"], cf.inv75, rtol=1e-12, atol=1e-15)),
                "D75 PP matches corrected_frontier": bool(np.allclose(
                    d75["PP"], cf.PP75, rtol=1e-12, atol=1e-15)),
            }
            for key, ok in agree.items():
                print(f"  {key}: {'yes' if ok else 'NO'}")
            out["agreement"] = agree
    with open(os.path.join(HERE, "d84_fixture.json"), "w") as handle:
        json.dump(out, handle, indent=1, default=float)
    return 0


if __name__ == "__main__":
    sys.exit(main())
