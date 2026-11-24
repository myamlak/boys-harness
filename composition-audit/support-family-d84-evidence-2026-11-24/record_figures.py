"""Read the deposited cell records' figure column and rebuild the scaled
figure used by the prospective checks.

The record's `figure` column is the per-cell accuracy figure.  The scaled
figure is the maximum over cells divided by the reference F0 at the point,
as in the article's reproduction package:

    tilde_figure(t) = max_cells figure(cell,t) / F0(x).

The accessor formula of boys v5.3.2 is checked against the recorded
column as a consistency test, not as an input.
"""

from __future__ import annotations

import csv
import hashlib
import lzma
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def harness_root():
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
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "fc", os.path.join(REPO, "composition-audit", "fullcube",
                       "fullcube_decomposition.py"))
fc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fc)

COLUMNS = [
    "precision", "route", "scheme", "partition", "pack", "form", "exp",
    "n", "x_nominal", "x", "value", "ref_hi", "ref_lo", "abs_err", "f0",
    "fn", "figure", "half_digit", "limit", "accessor_ratio", "figure_ratio",
    "pack_live", "partition_in_domain", "argument_kind", "status",
]
F64_REL = 1e-13
KFLAT_FIG = 2.0 ** -22
F32_FLOOR = 1.5e-7
PLAIN_TERM = 1e-7


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_figures(path):
    """Read figure, fn and the signed ratio per (lane, point, cell)."""
    lane_index = {lane: i for i, lane in enumerate(fc.LANES)}
    level_index = {
        axis: {name: i for i, name in enumerate(fc.LEVELS[axis])}
        for axis in fc.AXES
    }
    xs, ns = [], set()
    recs = []
    with lzma.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        if next(reader) != COLUMNS:
            raise SystemExit(f"{path}: unexpected header")
        for row in reader:
            if row[24] != "ok":
                continue
            cell = 0
            for pos, (axis, stride) in enumerate(zip(fc.AXES, fc.STRIDES)):
                cell += level_index[axis][row[1 + pos]] * stride
            n = int(row[7])
            x = float(row[8])
            xs.append(x)
            ns.add(n)
            value = float(row[10])
            ref_hi = float(row[11])
            ref_lo = float(row[12])
            f0 = float(row[14])
            signed = (value - ref_hi) - ref_lo
            recs.append((lane_index[row[0]], n, x, cell, float(row[16]),
                         float(row[15]),
                         signed / f0 if f0 != 0.0 else 0.0))
    ns = sorted(ns)
    xs = sorted(set(xs))
    x_index = {v: i for i, v in enumerate(xs)}
    n_index = {v: i for i, v in enumerate(ns)}
    npoints = len(ns) * len(xs)
    FIG = np.full((len(fc.LANES), npoints, fc.N_CELLS), np.nan)
    FN = np.full_like(FIG, np.nan)
    RATIO = np.full_like(FIG, np.nan)
    for lane, n, x, cell, fig, fn, ratio in recs:
        p = n_index[n] * len(xs) + x_index[x]
        FIG[lane, p, cell] = fig
        FN[lane, p, cell] = fn
        RATIO[lane, p, cell] = ratio
    return {"FIG": FIG, "FN": FN, "RATIO": RATIO, "xs": xs, "ns": ns,
            "sha256": sha256_file(path)}


def scaled_figure(FIG, FN, lane_name, nx):
    """Return fig_max / F0 at each point and the accessor comparison."""
    lane = fc.LANES.index(lane_name)
    fn = FN[lane, :, 0]
    f0 = FN[lane, np.arange(FIG.shape[1]) % nx, 0]
    fig_max = FIG[lane].max(axis=1)
    if lane_name == "fp64":
        predict_max = F64_REL * f0
    else:
        predict_max = np.maximum(KFLAT_FIG * fn, F32_FLOOR) + PLAIN_TERM
    return {
        "tilde": fig_max / f0,
        "fig_max": fig_max,
        "f0": f0,
        "fn": fn,
        "max_rel_dev": float((np.abs(fig_max - predict_max)
                              / np.abs(fig_max)).max()),
        "max_abs_dev": float(np.abs(fig_max - predict_max).max()),
    }
