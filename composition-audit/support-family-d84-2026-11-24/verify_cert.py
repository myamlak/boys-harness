"""Independent D84 certificate verification with the checker's primitives.

Set DESIGN_RULE_CHECK_DIR to the directory holding design_rule_check.py; that
module is external to this deposit.
"""

import importlib.util
import itertools
import math
import os
import sys

HERE = os.environ.get("DESIGN_RULE_CHECK_DIR")
if not HERE:
    raise SystemExit(
        "set DESIGN_RULE_CHECK_DIR to the directory holding "
        "design_rule_check.py")
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location(
    "drc", os.path.join(HERE, "design_rule_check.py"))
drc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(drc)

levels = [2, 2, 3, 2, 3, 2]
cones = [frozenset([0, 1, 2, 4, 5]), frozenset([0, 1, 2, 3])]


def d84_cells():
    cells = []
    for r in range(2):
        for s in range(2):
            for p in range(3):
                for f in range(3):
                    for e in range(2):
                        cells.append((r, s, p, 0, f, e))
    for r in range(2):
        for s in range(2):
            for p in range(3):
                cells.append((r, s, p, 1, 0, 0))
    return sorted(set(cells))


cells = d84_cells()
print("cells", len(cells), "distinct", len(set(cells)))
ranks = []
for p in (1009, 2017, 3001):
    s2 = drc.sqrt_mod(2, p)
    s3 = drc.sqrt_mod(3, p)
    rows_m, comp_m, perm_m = drc.build_family(
        levels, cones, lambda n: drc.axis_rows_mod(n, s2, s3, p))
    matrix = [[drc.phi_mod(rows_m, comp_m, cell, index, p)
               for index in perm_m] for cell in cells]
    ranks.append(drc.rank_mod(matrix, p))
print("mod ranks", ranks)

rows = [drc.axis_rows_float(n) for n in levels]
components = []
for choices in itertools.product(*[range(len(r)) for r in rows]):
    support = frozenset(i for i, c in enumerate(choices) if c != 0)
    components.append((choices, support))
permitted = [i for i, (_, support) in enumerate(components)
             if any(support <= cone for cone in cones)]
permitted, outside, targets = drc.family_split(components, permitted)
print("permitted", len(permitted), "outside", len(outside),
      "targets", len(targets))
omega_ii, omega_iq = drc.design_matrices(
    cells, cells, rows, components, permitted, outside)
inverse = drc.invert(omega_ii)
alias = drc.matmul(inverse, omega_iq)
position = {index: i for i, index in enumerate(permitted)}
kappa_target = max(
    max(abs(alias[position[alpha]][j]) / drc.g_float(rows, components, beta)
        for j, beta in enumerate(outside))
    for alpha in targets)
g_kappa = [drc.g_float(rows, components, alpha)
           * max(abs(alias[position[alpha]][j])
                 / drc.g_float(rows, components, beta)
                 for j, beta in enumerate(outside))
           for alpha in targets]
print("kappa_target", kappa_target)
print("g_kappa range", min(g_kappa), max(g_kappa),
      "sum", math.fsum(g_kappa))
cube = [(r, s, p, k, f, e) for r in range(2) for s in range(2)
        for p in range(3) for k in range(2) for f in range(3)
        for e in range(2)]
kappa_joint_values = []
for j, beta in enumerate(outside):
    ratios = []
    for cell in cube:
        residual = drc.phi_float(rows, components, cell, beta) - math.fsum(
            alias[position[alpha]][j]
            * drc.phi_float(rows, components, cell, alpha)
            for alpha in targets)
        ratios.append(abs(residual) / drc.g_float(rows, components, beta))
    kappa_joint_values.append(max(ratios))
print("kappa_joint max", max(kappa_joint_values),
      "values", sorted(set(round(v, 6) for v in kappa_joint_values)))
amp = max(
    math.fsum(abs(inverse[position[alpha]][i])
              * drc.g_float(rows, components, permitted[i])
              for i in range(len(permitted)))
    for alpha in targets)
print("amplification", amp)
print("gram condition", drc.gram_condition(omega_ii),
      "model condition", math.sqrt(drc.gram_condition(omega_ii)))
average_abs = [math.fsum(
    abs(drc.phi_float(rows, components, cell, index))
    for cell in cells) / len(cells) for index in permitted]
measurement = [math.fsum(abs(inverse[position[alpha]][i])
                         * average_abs[i] for i in range(len(permitted)))
               for alpha in targets]
print("measurement max target", max(measurement))
print("measurement G-sum",
      math.fsum(drc.g_float(rows, components, alpha)
                * measurement[targets.index(alpha)] for alpha in targets))

# 12x12 Kronecker block determinant
route = [[1.0, 1.0], [1.0, -1.0]]
scheme = [[1.0, 1.0], [1.0, -1.0]]
part = [[1.0, 1.0, 1.0],
        [drc.SQRT_3_2, 0.0, -drc.SQRT_3_2],
        [-drc.INV_SQRT_2, 2.0 * drc.INV_SQRT_2, -drc.INV_SQRT_2]]


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


det_r = route[0][0] * route[1][1] - route[0][1] * route[1][0]
det_s = scheme[0][0] * scheme[1][1] - scheme[0][1] * scheme[1][0]
det_p = det3(part)
print("det route", det_r, "det scheme", det_s, "det partition", det_p,
      "=3sqrt3", 3 * math.sqrt(3))
block_det = ((-2.0) ** 12) * (det_r ** 6) * (det_s ** 6) * (det_p ** 4)
print("12x12 block determinant", block_det)
