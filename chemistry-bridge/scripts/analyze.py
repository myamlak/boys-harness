#!/usr/bin/env python3
"""Paper-A chemistry bridge: offline analysis of the two probe arms.

Reads the files written by qcx-bridge-probe for the default class-row arm and
the Paper-A policy arm, all built from the same source tree and the same input:

    <prefix>.summary.txt     key=value run summary
    <prefix>.eri.txt         n, then "i j k l value" over the full dense tensor
    <prefix>.matrices.txt    "matrix name rows cols" then row-major values
    <prefix>.args.txt        one line per unique Boys argument:
                             nmax count x engine[0..nmax] class[0..nmax] papera[0..nmax]

It reports:
  * the captured argument stream is the same in both arms (nmax, x bits, count);
  * each arm's engine ladder is reproduced bit-for-bit by the policy it was
    built with, and the other policy column agrees across arms;
  * per-order |F_n - reference| for both routes against an 80-digit mpmath
    reference (validated against mp.quad at selected arguments), each error
    also as a multiple of the kernel's 1e-13 * F_0 contract;
  * per-order |F_n(class) - F_n(paper-a)|, the difference that drove the ERIs;
  * max |delta ERI| and delta E_2e at fixed densities (the GWH seed is
    route-independent by construction; converged densities are cross-evaluated);
  * delta E_total for the fixed-density comparison and for the two SCF runs;
  * every energy number against chemical accuracy (1.594e-3 Ha) and against the
    kernel's 1e-13 relative contract scaled by |E_2e|.

Usage: python -I analyze.py <default-prefix> <paper-a-prefix> <out.txt>

The script is analysis only: no qcx code and no boys code is imported or run.
"""

import struct
import sys

import mpmath as mp

mp.mp.dps = 80

CHEMICAL_ACCURACY = 1.594e-3   # Ha
KERNEL_RELATIVE = 1e-13        # per-value F_n contract is 1e-13 * F_0(x)


def kahan_sum(values):
    total = 0.0
    compensation = 0.0
    for value in values:
        adjusted = value - compensation
        updated = total + adjusted
        compensation = (updated - total) - adjusted
        total = updated
    return total


def read_summary(path):
    summary = {}
    with open(path, "r", encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            key, value = line.split("=", 1)
            summary[key] = value
    return summary


def read_matrices(path):
    matrices = {}
    with open(path, "r", encoding="ascii") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    index = 0
    while index < len(lines):
        header = lines[index].split()
        assert header[0] == "matrix", lines[index]
        name, rows, cols = header[1], int(header[2]), int(header[3])
        index += 1
        matrix = []
        for _ in range(rows):
            matrix.append([float(token) for token in lines[index].split()])
            index += 1
        assert len(matrix[0]) == cols
        matrices[name] = matrix
    return matrices


def read_eri(path):
    """The dense tensor as {(i,j,k,l): value}, in the file's row-major order."""
    values = {}
    with open(path, "r", encoding="ascii") as handle:
        n = int(handle.readline())
        for line in handle:
            i, j, k, l, value = line.split()
            values[(int(i), int(j), int(k), int(l))] = float(value)
    assert len(values) == n ** 4, (len(values), n)
    return n, values


def read_eri_dimension(path):
    with open(path, "r", encoding="ascii") as handle:
        return int(handle.readline())


# Above this many entries the two tensors are compared in a lockstep stream
# instead of two Python dictionaries: cc-pVTZ (58 functions) is 11.3M entries
# and two dicts of 4-tuples would cost several GB. The streamed path computes
# the same numbers from the same file order.
STREAM_THRESHOLD = 2_000_000


class Kahan:
    """Compensated accumulator, the same variant the probe uses."""

    __slots__ = ("total", "compensation")

    def __init__(self):
        self.total = 0.0
        self.compensation = 0.0

    def add(self, value):
        adjusted = value - self.compensation
        updated = self.total + adjusted
        self.compensation = (updated - self.total) - adjusted
        self.total = updated


def stream_eri(path_default, path_paper, density):
    """One lockstep pass over both tensor dumps at a fixed density.

    Returns the comparison statistics and the density-weighted difference of
    the Coulomb and exchange sums. The exchange contraction uses the identity

        sum_{uvls} D_uv D_ls (ul|sv) = sum_{abcd} D_ad D_bc (ab|cd),

    a relabelling that lets the exchange term read the file's own order. The
    difference is taken per entry before weighting, so it is the density-
    weighted integral error, not the difference of two large sums.
    """
    max_eri = 0.0
    max_delta = 0.0
    max_delta_index = None
    sum_square = 0.0
    differing = 0
    count = 0
    coulomb = Kahan()
    exchange = Kahan()
    with open(path_default, "r", encoding="ascii") as left, \
            open(path_paper, "r", encoding="ascii") as right:
        n = int(left.readline())
        assert int(right.readline()) == n
        for left_line, right_line in zip(left, right):
            i, j, k, l, left_value = left_line.split()
            i2, j2, k2, l2, right_value = right_line.split()
            i, j, k, l = int(i), int(j), int(k), int(l)
            i2, j2, k2, l2 = int(i2), int(j2), int(k2), int(l2)
            assert (i2, j2, k2, l2) == (i, j, k, l), (left_line, right_line)
            a = float(left_value)
            b = float(right_value)
            count += 1
            max_eri = max(max_eri, abs(a), abs(b))
            if bits(a) != bits(b):
                differing += 1
            difference = a - b
            sum_square += difference * difference
            if abs(difference) > max_delta:
                max_delta = abs(difference)
                max_delta_index = (i, j, k, l, a, b)
            dij = density[i][j]
            dkl = density[k][l]
            coulomb.add(0.5 * dij * dkl * difference)
            exchange.add(-0.25 * density[i][l] * density[j][k] * difference)
    assert count == n ** 4, (count, n)
    if max_delta_index is None:  # identical tensors
        max_delta_index = (0, 0, 0, 0, 0.0, 0.0)
    return {
        "n": n,
        "entries": count,
        "differing": differing,
        "max_eri": max_eri,
        "max_delta": max_delta,
        "max_delta_index": max_delta_index,
        "rms": (sum_square / count) ** 0.5,
        "delta_coulomb": coulomb.total,
        "delta_exchange": exchange.total,
        "delta_2e": coulomb.total + exchange.total,
    }


def read_args(path):
    arguments = []
    with open(path, "r", encoding="ascii") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            nmax = int(fields[0])
            count = int(fields[1])
            x = float(fields[2])
            width = nmax + 1
            body = [float(token) for token in fields[3:]]
            assert len(body) == 3 * width, (nmax, len(body))
            arguments.append({
                "nmax": nmax,
                "count": count,
                "x": x,
                "x_bits": struct.pack("<d", x),
                "engine": body[0:width],
                "class": body[width:2 * width],
                "papera": body[2 * width:3 * width],
            })
    return arguments


def bits(value):
    return struct.pack("<d", value)


def series_ladder(nmax, x_float):
    """F_0(x)..F_nmax(x) from the Taylor series about x = 0.

    F_n(x) = sum_{k>=0} (-x)^k / (k! (2n + 2k + 1)),

    which converges for every x and is the accurate route at small x: the
    upward recurrence amplifies each step's error by about (2n+1)/(2x), so at
    x = 1e-6 the error reaches about 1e-2 by order 12 even at 80 digits. The
    alternating terms shrink by about x/(k+1); the loop stops once a term is
    below 1e-110, which needs about 45 terms at the x = 1e-2 branch boundary.
    """
    x = mp.mpf(x_float)
    threshold = mp.mpf("1e-110")
    ladders = []
    for n in range(nmax + 1):
        total = mp.mpf(0)
        power = mp.mpf(1)
        for k in range(250):
            term = power / (mp.factorial(k) * (2 * n + 2 * k + 1))
            total += -term if k % 2 else term
            if abs(term) < threshold:
                break
            power *= x
        else:
            raise RuntimeError("small-x series did not converge")
        ladders.append(total)
    return ladders


def reference_ladder(nmax, x_float):
    """F_0(x)..F_nmax(x) at mp.dps: series below 1e-2, recurrence above it.

    F_0(x) = sqrt(pi) / (2 sqrt(x)) * erf(sqrt(x)),
    F_{n+1}(x) = ((2n+1) F_n(x) - exp(-x)) / (2x).

    At x = 1e-2 the recurrence still carries about 40 digits at order 12
    (each step amplifies by (2n+1)/(2x); the product over 12 steps is about
    1e39 against 80 digits), while below it the series is exact. The former
    switch at 1e-6 left orders above 8 with a broken reference between the
    two branches, found by the validation below on cc-pVTZ.
    """
    x = mp.mpf(x_float)
    if x == 0:
        return [mp.mpf(1) / (2 * n + 1) for n in range(nmax + 1)]
    if x < mp.mpf("1e-2"):
        return series_ladder(nmax, x_float)
    exp_minus_x = mp.exp(-x)
    ladder = [mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))]
    for n in range(nmax):
        ladder.append(((2 * n + 1) * ladder[-1] - exp_minus_x) / (2 * x))
    return ladder


def reference_by_quad(nmax, x_float):
    """The same ladder by direct quadrature, for validating the recurrence."""
    x = mp.mpf(x_float)
    ladder = []
    for n in range(nmax + 1):
        if x == 0:
            ladder.append(mp.mpf(1) / (2 * n + 1))
        else:
            ladder.append(mp.quad(lambda t, n=n: t ** (2 * n) * mp.exp(-x * t * t), [0, 1]))
    return ladder


def rel(value, scale):
    return value / scale if scale else float("inf")


def fixed_density_energy(eri_values, density, core, n):
    """The closed-shell E_2e and Tr[DH] at one density, in a fixed order.

    E_2e = 1/2 sum_uvls D_uv D_ls (uv|ls) - 1/4 sum_uvls D_uv D_ls (ul|sv).

    The one-particle labels are spatial and D is the spin-summed density, so
    this is the same expression the SCF's Coulomb and exchange supermatrices
    build; the trace identity against HfResult.electronicEnergy is checked
    below. eri_values is the dictionary in the file's row-major order, so the
    summation order is the probe's own.
    """
    coulomb_terms = []
    exchange_terms = []
    core_terms = []
    for (i, j, k, l), value in eri_values.items():
        pair = density[i][j] * density[k][l]
        coulomb_terms.append(0.5 * pair * value)
        exchange_terms.append(-0.25 * pair * eri_values[(i, k, l, j)])
    for i in range(n):
        for j in range(n):
            core_terms.append(density[i][j] * core[i][j])
    return kahan_sum(coulomb_terms), kahan_sum(exchange_terms), kahan_sum(core_terms)


def fixed_density_delta(eri_default, eri_papera, density):
    """Delta of the Coulomb and exchange sums at a fixed density.

    Differencing the two Kahan energies directly would difference two large
    numbers; this form keeps only the integral's difference, so the number is
    the density-weighted integral error itself. Returns
    (delta_coulomb, delta_exchange); their sum is delta E_2e.
    """
    coulomb_terms = []
    exchange_terms = []
    for key, a in eri_default.items():
        i, j, k, l = key
        pair = density[i][j] * density[k][l]
        coulomb_terms.append(0.5 * pair * (a - eri_papera[key]))
        exchange_terms.append(-0.25 * pair
                              * (eri_default[(i, k, l, j)] - eri_papera[(i, k, l, j)]))
    return kahan_sum(coulomb_terms), kahan_sum(exchange_terms)


def matrix_delta(a, b):
    rows = len(a)
    cols = len(a[0])
    largest = 0.0
    identical = True
    for i in range(rows):
        for j in range(cols):
            if bits(a[i][j]) != bits(b[i][j]):
                identical = False
            largest = max(largest, abs(a[i][j] - b[i][j]))
    return largest, identical


def main():
    if len(sys.argv) != 4:
        print("usage: python -I analyze.py <default-prefix> <paper-a-prefix> <out.txt>",
              file=sys.stderr)
        return 2

    default_prefix, paper_prefix, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    out = open(out_path, "w", encoding="ascii")

    def emit(text=""):
        out.write(text + "\n")
        print(text)

    default_summary = read_summary(default_prefix + ".summary.txt")
    paper_summary = read_summary(paper_prefix + ".summary.txt")
    default_matrices = read_matrices(default_prefix + ".matrices.txt")
    paper_matrices = read_matrices(paper_prefix + ".matrices.txt")
    n = read_eri_dimension(default_prefix + ".eri.txt")
    n_paper = read_eri_dimension(paper_prefix + ".eri.txt")
    assert n == n_paper
    streaming = n ** 4 > STREAM_THRESHOLD
    default_eri = None
    paper_eri = None
    if not streaming:
        n, default_eri = read_eri(default_prefix + ".eri.txt")
        _, paper_eri = read_eri(paper_prefix + ".eri.txt")

    emit("Paper-A chemistry bridge analysis")
    emit("=" * 72)
    emit("")
    emit("Arms")
    emit("-" * 72)
    for label, summary in (("default", default_summary), ("paper-a", paper_summary)):
        emit(f"  {label:8s} variant={summary['variant']:20s} "
             f"captured={summary['capture_calls']:>7s} unique={summary['capture_unique']:>6s} "
             f"ladder_mismatches={summary['capture_ladder_mismatches']}")
    basis = default_summary.get("basis_source", "sto-3g-fixture")
    basis_dir = default_summary.get("basis_directory", "")
    emit(f"  n = {n} ({basis}{', ' + basis_dir if basis_dir else ''}; "
         f"spin-summed RHF)")
    # 'shells' counts the engine's contracted shells, which are not the
    # function count: one NWChem general-contraction block with several rows
    # is one engine shell. cc-pVTZ's 58 functions live in 10 such shells, the
    # same number as 6-31G*'s 18 - the coincidence is real, not stale data.
    emit(f"  engine shells={default_summary.get('shells', '?')} "
         f"pairs={default_summary.get('shell_pairs', '?')} "
         f"canonical quartets={default_summary.get('canonical_quartets', '?')} "
         f"screened survivors={default_summary.get('screened_quartets', '?')}")
    emit(f"  tensor comparison path: {'lockstep stream' if streaming else 'full dictionary'}")
    emit("")

    # ------------------------------------------------------------------ #
    # The argument stream: the same arguments in both arms, and each
    # arm's engine ladder reproduced by the policy it was built with.
    # ------------------------------------------------------------------ #
    emit("Argument stream and policy reproduction")
    emit("-" * 72)
    default_args = read_args(default_prefix + ".args.txt")
    paper_args = read_args(paper_prefix + ".args.txt")

    same_stream = len(default_args) == len(paper_args)
    if same_stream:
        for a, b in zip(default_args, paper_args):
            if (a["nmax"] != b["nmax"] or a["x_bits"] != b["x_bits"]
                    or a["count"] != b["count"]):
                same_stream = False
                break

    class_agree = True
    papera_agree = True
    engine_matches_own = True
    engine_other_matches = 0
    for a, b in zip(default_args, paper_args):
        if a["class"] != b["class"]:
            class_agree = False
        if a["papera"] != b["papera"]:
            papera_agree = False
        # The default arm ran the class row; the paper-a arm ran Paper A.
        if a["engine"] != a["class"] or b["engine"] != b["papera"]:
            engine_matches_own = False
        if b["engine"] == b["class"]:
            engine_other_matches += 1

    calls = sum(argument["count"] for argument in default_args)
    nmax_values = sorted({argument["nmax"] for argument in default_args})
    xs = [argument["x"] for argument in default_args]
    positive_xs = [x for x in xs if x > 0.0]

    emit(f"  unique arguments: default={len(default_args)} paper-a={len(paper_args)}; "
         f"same (nmax, x bits, multiplicity): {same_stream}")
    emit(f"  total engine calls in the stream: {calls}")
    emit(f"  nmax values present: {nmax_values}")
    emit(f"  x range: min={min(xs):.17g} max={max(xs):.17g}; "
         f"smallest positive x={min(positive_xs):.17g}")
    emit(f"  class-row column identical across arms: {class_agree}")
    emit(f"  paper-a column identical across arms: {papera_agree}")
    emit(f"  engine ladder reproduced by its own arm's policy, every argument: "
         f"{engine_matches_own}")
    emit(f"  (for contrast, paper-a arm engine ladder equals the class row on "
         f"{engine_other_matches}/{len(paper_args)} arguments)")
    emit(f"  repeated-call ladder inconsistencies within an arm: "
         f"default={default_summary['capture_ladder_mismatches']} "
         f"paper-a={paper_summary['capture_ladder_mismatches']}")
    emit("")

    # ------------------------------------------------------------------ #
    # The 80-digit reference. Validate whichever branch reference_ladder
    # takes against quadrature at selected arguments before using it.
    # ------------------------------------------------------------------ #
    emit("Reference validation (against direct mp.quad and the closed form, 80 dps)")
    emit("-" * 72)
    largest_nmax = max(nmax_values)

    # Quadrature resolves the exponential only while x is large enough to
    # matter at 80 digits; below about 1e-6 the integrand is 1 to machine
    # precision at quadrature's own tolerance and mp.quad stops being a
    # reference. The series is validated down there in its own right below.
    quad_candidates = []
    quad_safe = [x for x in sorted(set(positive_xs)) if x > 1e-6]
    if quad_safe:
        picks = [quad_safe[0], quad_safe[len(quad_safe) // 4],
                 quad_safe[len(quad_safe) // 2], quad_safe[(3 * len(quad_safe)) // 4],
                 quad_safe[-1]]
        for x in picks:
            if x not in quad_candidates:
                quad_candidates.append(x)

    worst_quad = 0.0
    for x in quad_candidates:
        ladder = reference_ladder(largest_nmax, x)
        quadrature = reference_by_quad(largest_nmax, x)
        difference = max(abs(a - b) for a, b in zip(ladder, quadrature))
        worst_quad = max(worst_quad, difference)
        branch = "series" if mp.mpf(x) < mp.mpf("1e-2") else "recurrence"
        emit(f"  x={x:.17g} ({branch}): max |reference - quad| = "
             f"{mp.nstr(difference, 8)} (orders 0..{largest_nmax})")
    emit(f"  worst reference-vs-quad over the checked arguments: "
         f"{mp.nstr(worst_quad, 8)}")

    for x in (1e-3, 1e-5, 1e-8):
        ladder = series_ladder(largest_nmax, x)
        quadrature = reference_by_quad(largest_nmax, x)
        difference = max(abs(a - b) for a, b in zip(ladder, quadrature))
        emit(f"  x={x:.0e} (series branch): max |series - quad| = "
             f"{mp.nstr(difference, 8)}")
    if positive_xs:
        x_small = min(positive_xs)
        ladder = reference_ladder(largest_nmax, x_small)
        closed = (mp.sqrt(mp.pi) / (2 * mp.sqrt(mp.mpf(x_small)))
                  * mp.erf(mp.sqrt(mp.mpf(x_small))))
        emit(f"  smallest x={x_small:.17g} (series branch): "
             f"|F_0 - closed form| = {mp.nstr(abs(ladder[0] - closed), 8)}")
    emit("")

    # ------------------------------------------------------------------ #
    # Per-order errors: each route against the reference, and the two
    # routes against each other, normalised by the 1e-13 * F_0 contract.
    # ------------------------------------------------------------------ #
    emit("Per-order F_n errors against the 80-digit reference")
    emit("-" * 72)
    emit("  orders are the Boys order n; errors are max over the captured arguments")
    emit("  that request that order (nmax >= n). 'bound' = 1e-13 * F_0(x) per value.")
    emit("")
    header = ("  n    max|class-ref|   x(class max)     class/bound"
              "      max|papera-ref|   x(papera max)   papera/bound"
              "      max|class-papera|   x(dF max)      dF/bound")
    emit(header)

    reference_rows = []
    for argument in default_args:
        ladder = reference_ladder(argument["nmax"], argument["x"])
        reference_rows.append((argument, ladder,
                               mp.mpf(KERNEL_RELATIVE) * ladder[0]))

    per_order = []
    for order in range(largest_nmax + 1):
        class_err = 0.0
        class_at = 0.0
        class_ratio = 0.0
        paper_err = 0.0
        paper_at = 0.0
        paper_ratio = 0.0
        delta = 0.0
        delta_at = 0.0
        delta_ratio = 0.0
        for argument, ladder, bound in reference_rows:
            if argument["nmax"] < order:
                continue
            reference = ladder[order]
            x = argument["x"]

            err_class = abs(mp.mpf(argument["class"][order]) - reference)
            err_paper = abs(mp.mpf(argument["papera"][order]) - reference)
            err_delta = abs(mp.mpf(argument["class"][order])
                            - mp.mpf(argument["papera"][order]))

            if err_class > class_err:
                class_err, class_at = err_class, x
            class_ratio = max(class_ratio, float(err_class / bound))
            if err_paper > paper_err:
                paper_err, paper_at = err_paper, x
            paper_ratio = max(paper_ratio, float(err_paper / bound))
            if err_delta > delta:
                delta, delta_at = err_delta, x
            delta_ratio = max(delta_ratio, float(err_delta / bound))

        per_order.append({
            "order": order,
            "class_err": float(class_err), "class_at": class_at,
            "class_ratio": class_ratio, "paper_err": float(paper_err),
            "paper_at": paper_at, "paper_ratio": paper_ratio,
            "delta": float(delta), "delta_at": delta_at, "delta_ratio": delta_ratio,
        })
        emit(f"  {order:<3d} {float(class_err):>16.6e} {class_at:>16.6e} "
             f"{class_ratio:>14.3e} {float(paper_err):>18.6e} {paper_at:>16.6e} "
             f"{paper_ratio:>13.3e} {float(delta):>18.6e} {delta_at:>16.6e} "
             f"{delta_ratio:>13.3e}")

    worst_class = max(per_order, key=lambda row: row["class_ratio"])
    worst_paper = max(per_order, key=lambda row: row["paper_ratio"])
    worst_delta = max(per_order, key=lambda row: row["delta_ratio"])
    emit("")
    emit(f"  worst class-row error / bound: order {worst_class['order']}, "
         f"{worst_class['class_ratio']:.3e} x bound, at x={worst_class['class_at']:.17g}")
    emit(f"  worst paper-a error / bound:   order {worst_paper['order']}, "
         f"{worst_paper['paper_ratio']:.3e} x bound, at x={worst_paper['paper_at']:.17g}")
    emit(f"  worst |class - paper-a| / bound: order {worst_delta['order']}, "
         f"{worst_delta['delta_ratio']:.3e} x bound, at x={worst_delta['delta_at']:.17g}")
    emit("")

    # ------------------------------------------------------------------ #
    # The ERI tensors.
    # ------------------------------------------------------------------ #
    emit("ERI tensor comparison (same input, same engine, two Boys routes)")
    emit("-" * 72)
    streamed = {}
    if streaming:
        # The comparison and the shared-seed density contraction in one pass.
        stats = stream_eri(default_prefix + ".eri.txt", paper_prefix + ".eri.txt",
                           default_matrices["seed_density"])
        streamed["seed (GWH, shared)"] = stats
        entries = stats["entries"]
        differing = stats["differing"]
        max_eri = stats["max_eri"]
        max_delta = stats["max_delta"]
        rms = stats["rms"]
        max_delta_index = stats["max_delta_index"]
    else:
        max_eri = 0.0
        max_delta = 0.0
        max_delta_index = None
        sum_square = 0.0
        differing = 0
        assert list(default_eri.keys()) == list(paper_eri.keys())
        for key, a in default_eri.items():
            b = paper_eri[key]
            max_eri = max(max_eri, abs(a))
            if bits(a) != bits(b):
                differing += 1
            difference = a - b
            sum_square += difference * difference
            if abs(difference) > max_delta:
                max_delta = abs(difference)
                max_delta_index = (*key, a, b)

        entries = len(default_eri)
        rms = (sum_square / entries) ** 0.5
        if max_delta_index is None:  # identical tensors
            key = next(iter(default_eri))
            max_delta_index = (*key, default_eri[key], paper_eri[key])
    i, j, k, l, a, b = max_delta_index
    emit(f"  entries: {entries}; bitwise differing: {differing}")
    emit(f"  max |ERI| = {max_eri:.6e}")
    emit(f"  max |delta ERI| = {max_delta:.6e} at ({i},{j},{k},{l}): "
         f"{a:.17g} vs {b:.17g}")
    emit(f"  RMS delta ERI = {rms:.6e}")
    emit(f"  max |delta ERI| / max|ERI| = {max_delta / max_eri:.3e}")
    emit("")

    # ------------------------------------------------------------------ #
    # Fixed-density energies. S and H are Boys-free: prove they are
    # identical bit-for-bit before treating the seed density as shared.
    # ------------------------------------------------------------------ #
    emit("Fixed-density energy comparison")
    emit("-" * 72)
    for name in ("overlap", "core"):
        largest, identical = matrix_delta(default_matrices[name], paper_matrices[name])
        emit(f"  {name}: bitwise identical = {identical}, max |delta| = {largest:.3e}")
    for name in ("seed_density", "default_density", "tight_density"):
        largest, identical = matrix_delta(default_matrices[name], paper_matrices[name])
        emit(f"  {name}: bitwise identical = {identical}, max |delta| = {largest:.3e}")
    emit("")

    core = default_matrices["core"]
    e_nuc = float(default_summary["e_nuc"])
    densities = {
        "seed (GWH, shared)": default_matrices["seed_density"],
        "default-arm SCF density": default_matrices["default_density"],
        "paper-a-arm SCF density": paper_matrices["default_density"],
        "default-arm tight density": default_matrices["tight_density"],
        "paper-a-arm tight density": paper_matrices["tight_density"],
    }

    emit("  delta E_2e at fixed densities, full expression (Coulomb - exchange, Kahan):")
    fixed_density = {}
    for label, density in densities.items():
        if streaming:
            stats = streamed.get(label)
            if stats is None:
                stats = stream_eri(default_prefix + ".eri.txt",
                                   paper_prefix + ".eri.txt", density)
                streamed[label] = stats
            delta_coulomb = stats["delta_coulomb"]
            delta_exchange = stats["delta_exchange"]
        else:
            delta_coulomb, delta_exchange = fixed_density_delta(
                default_eri, paper_eri, density)
        delta_2e = delta_coulomb + delta_exchange
        fixed_density[label] = delta_2e
        emit(f"    at {label:28s}: delta E_2e = {delta_2e:+.6e} Ha"
             f"   ({abs(delta_2e) / CHEMICAL_ACCURACY:.3e} x chemical accuracy;"
             f" Coulomb {delta_coulomb:+.3e}, exchange {delta_exchange:+.3e})")
    emit("")

    # Cross-check: recompute each arm's full E_2e at the shared seed density
    # and compare with the probe's own Kahan value. This proves the offline
    # convention matches the probe's and that the probe's numbers are sound.
    probe_seed_default = float(default_summary["seed_e2e"])
    probe_seed_paper = float(paper_summary["seed_e2e"])
    if streaming:
        emit("  probe vs streamed recomputation at the shared seed density:")
        emit(f"    default arm E_2e (probe) = {probe_seed_default:.17g}")
        emit(f"    paper-a arm E_2e (probe) = {probe_seed_paper:.17g}")
        seed_delta_probe = probe_seed_default - probe_seed_paper
        seed_delta_stream = fixed_density["seed (GWH, shared)"]
        emit(f"    seed delta E_2e (probe difference) = {seed_delta_probe:+.6e} Ha")
        emit(f"    seed delta E_2e (streamed)         = {seed_delta_stream:+.6e} Ha")
        emit(f"    |difference| = "
             f"{abs(seed_delta_probe - seed_delta_stream):.3e}")
        seed_2e_default = probe_seed_default
    else:
        emit("  probe vs offline recomputation at the shared seed density:")
        seed_density = default_matrices["seed_density"]
        seed_c_default, seed_x_default, seed_core_default = fixed_density_energy(
            default_eri, seed_density, core, n)
        seed_c_paper, seed_x_paper, seed_core_paper = fixed_density_energy(
            paper_eri, seed_density, core, n)
        seed_2e_default = seed_c_default + seed_x_default
        seed_2e_paper = seed_c_paper + seed_x_paper
        emit(f"    default arm: probe {probe_seed_default:.17g}, offline "
             f"{seed_2e_default:.17g}, |difference| = "
             f"{abs(probe_seed_default - seed_2e_default):.3e}"
             f" (Coulomb {seed_c_default:.6e}, exchange {seed_x_default:.6e})")
        emit(f"    paper-a arm: probe {probe_seed_paper:.17g}, offline "
             f"{seed_2e_paper:.17g}, |difference| = "
             f"{abs(probe_seed_paper - seed_2e_paper):.3e}"
             f" (Coulomb {seed_c_paper:.6e}, exchange {seed_x_paper:.6e})")
        emit(f"    seed delta E_2e from the two recomputations: "
             f"{seed_2e_default - seed_2e_paper:+.6e} Ha")
        emit(f"    seed total (e_nuc + Tr[DH] + E_2e): default "
             f"{e_nuc + seed_core_default + seed_2e_default:.17g}, paper-a "
             f"{e_nuc + seed_core_paper + seed_2e_paper:.17g}")
    emit("")

    # The trace identity: at each arm's own returned density, Tr[DH] + E_2e
    # must equal HfResult.electronicEnergy. It closes only if the exchange
    # term is present, so this is the convention's own check. In the streamed
    # path the probe has already evaluated it with its own tensor, and its
    # residuals are quoted directly.
    emit("  trace identity Tr[DH] + E_2e = E_electronic at each returned density:")
    if streaming:
        for gate in ("default", "tight"):
            emit(f"    {gate} gate: probe residual default="
                 f"{float(default_summary[f'{gate}_trace_identity_residual']):.3e}, "
                 f"paper-a={float(paper_summary[f'{gate}_trace_identity_residual']):.3e}")
    else:
        identity_checks = (
            ("default gate, default arm", default_eri, "default_density",
             default_summary["scf_default_electronic_energy"]),
            ("default gate, paper-a arm", paper_eri, "default_density",
             paper_summary["scf_default_electronic_energy"]),
            ("tight gate, default arm", default_eri, "tight_density",
             default_summary["scf_tight_electronic_energy"]),
            ("tight gate, paper-a arm", paper_eri, "tight_density",
             paper_summary["scf_tight_electronic_energy"]),
        )
        for label, eri, density_name, electronic in identity_checks:
            arm_density = (default_matrices if "default arm" in label
                           else paper_matrices)[density_name]
            c, x, tr = fixed_density_energy(eri, arm_density, core, n)
            offline = tr + c + x
            emit(f"    {label}: offline {offline:.17g}, result {float(electronic):.17g}, "
                 f"|residual| = {abs(offline - float(electronic)):.3e}")
    emit("")

    # The density itself changes between routes; report how much, so the
    # relaxed part of any SCF energy difference is not mistaken for the
    # integral error.
    emit("  density relaxation between the arms:")
    for label in ("default_density", "tight_density"):
        largest, identical = matrix_delta(default_matrices[label],
                                         paper_matrices[label])
        emit(f"    {label}: max |delta D| = {largest:.3e}, bitwise identical = {identical}")
    emit("")

    # ------------------------------------------------------------------ #
    # SCF energies as the probes returned them, and the contract ratios.
    # ------------------------------------------------------------------ #
    emit("SCF comparison (as returned by each arm)")
    emit("-" * 72)
    for gate in ("default", "tight"):
        d_total = float(default_summary[f"scf_{gate}_total_energy"])
        p_total = float(paper_summary[f"scf_{gate}_total_energy"])
        d_electronic = float(default_summary[f"scf_{gate}_electronic_energy"])
        p_electronic = float(paper_summary[f"scf_{gate}_electronic_energy"])
        delta_total = p_total - d_total
        delta_electronic = p_electronic - d_electronic
        emit(f"  {gate} gate:")
        emit(f"    default: E_total={d_total:.17g} Ha, converged="
             f"{default_summary[f'scf_{gate}_converged']}, "
             f"iterations={default_summary[f'scf_{gate}_iterations']}, "
             f"achieved dE={float(default_summary[f'scf_{gate}_achieved_energy_delta']):.3e}")
        emit(f"    paper-a: E_total={p_total:.17g} Ha, converged="
             f"{paper_summary[f'scf_{gate}_converged']}, "
             f"iterations={paper_summary[f'scf_{gate}_iterations']}, "
             f"achieved dE={float(paper_summary[f'scf_{gate}_achieved_energy_delta']):.3e}")
        emit(f"    delta E_total (paper-a - default) = {delta_total:+.6e} Ha")
        emit(f"      = {abs(delta_total) / CHEMICAL_ACCURACY:.3e} x chemical accuracy")
        emit(f"      = {abs(delta_total) / (KERNEL_RELATIVE * abs(d_total)):.3e} "
             f"x the 1e-13 relative contract scaled by |E_total|")
        emit(f"    delta E_electronic = {delta_electronic:+.6e} Ha")
        emit("")

    # The headline fixed-density number: at the route-independent GWH seed,
    # the only difference between the arms is the integral tensor.
    emit("Headline numbers")
    emit("-" * 72)
    headline = fixed_density["seed (GWH, shared)"]
    emit(f"  max |delta F_n| between the routes:        "
         f"{max(row['delta'] for row in per_order):.6e} (see per-order table)")
    emit(f"  max |delta F_n| / contract:                "
         f"{max(row['delta_ratio'] for row in per_order):.3e} x 1e-13 F_0")
    emit(f"  max |delta ERI|:                           {max_delta:.6e}")
    emit(f"  delta E_2e at the shared seed density:     {headline:+.6e} Ha")
    emit(f"  delta E_total at the shared seed density:  {headline:+.6e} Ha "
         f"(S and H are Boys-free and bit-identical)")
    emit(f"  delta E_2e / chemical accuracy:            "
         f"{abs(headline) / CHEMICAL_ACCURACY:.3e}")
    emit(f"  delta E_2e / (1e-13 |E_2e|):               "
         f"{abs(headline) / (KERNEL_RELATIVE * abs(seed_2e_default)):.3e} "
         f"(|E_2e(seed)| = {abs(seed_2e_default):.6e} Ha)")
    tight_key = "default-arm tight density"
    emit(f"  delta E_2e at the default arm's tight density: "
         f"{fixed_density[tight_key]:+.6e} Ha"
         f" ({abs(fixed_density[tight_key]) / CHEMICAL_ACCURACY:.3e} x chemical accuracy)")
    emit("")

    out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
