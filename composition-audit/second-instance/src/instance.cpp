// d-second-instance: real, small configuration cubes on which a fractional
// orthogonal array of strength 2 certifies the additive composition of
// per-axis accuracy bounds, recorded against claimgate 0.1.0.
//
// ---------------------------------------------------------------------------
// The instances
// ---------------------------------------------------------------------------
// Both instances evaluate a sum of independent series approximations:
//
//   Y(a) = T1(a1) + T2(a2) + ... + Tn(an),
//
// where axis Ai selects one of two truncation levels for term Ti.  Each term
// is computed in long double from its alternating series, then rounded to a
// common fixed-point quantum Q = 2^-40.  The final addition is performed on
// the integer numerators, so it is exact: N = n1 + ... + nn and Y = N*Q.
//
//   three-axis: A1 exp(-1), A2 ln 2, A3 pi/4      full 8 cells,  fraction 4
//   four-axis:  the same three plus A4 cos(1)     full 16 cells, fraction 8
//
//   A1: exp(-1) = sum_{k>=0} (-1)^k / k!          levels: K=5, K=8
//   A2: ln 2    = sum_{k>=1} (-1)^(k+1) / k       levels: K=8, K=12
//   A3: pi/4    = sum_{k>=0} (-1)^k / (2k+1)      levels: K=6, K=10
//   A4: cos(1)  = sum_{k>=0} (-1)^k / (2k)!       levels: K=4, K=6
//
// The exact reference is the platform's long-double elementary function for
// each term and the sum of them for the joint call.  The per-axis bound b_i is
// the alternating-series remainder bound at the coarser level plus Q/2 for the
// fixed-point rounding and a 1e-15 evaluation/reference margin.  The joint
// additive bound is B = sum b_i; the reported allowance is eta = 1e-12 (the
// measurement floor; the structural high-order remainder is zero), so the
// certified joint bound is B + eta.
//
// Dependency cones (the error-carrying operations):
//   op_i: the series evaluation and its fixed-point rounding -> cone {Ai}.
//   op_add: the integer addition and the power-of-two scaling have a local
//           error that is identically zero, so they add no axis to the cone;
//           the paper applies the same reading to an exact final store, whose
//           operand may depend on every live axis but which adds no axis.
//           The operation therefore contributes no term to the first-order
//           error decomposition and is not a cone in the table below.
// Hence w_p = 1 and U_p = {A1,...,An} for both instances.
//
// The designs:
//   full factorial: 2^n cells;
//   fractional: the even-parity half fraction (I = A1*...*An), 2^(n-1) cells:
//   OA(4,3,2,2) for three axes and OA(8,4,2,2) for four, strength 2.
//   Main effects alias with complementary (n-1)-factor interactions; in the
//   four-axis fraction the two-factor interactions also alias with each other.
//   Every interaction coefficient is zero by the cone structure, so the
//   fractional mean and main effects estimate the full-factorial ones within
//   the allowance, and half the cells are saved.
//
// claimgate 0.1.0 is used for the parts of the claim model it has: each
// per-axis bound and the joint bound are registered as claims and every cell
// of each design is measured against the same independent reference.  The
// library has no axes, cones, orthogonal arrays, contrasts, aliases, or
// allowance term; the design computation below is therefore each instance's
// own, and its result is recorded in the CSV files beside the claimgate
// reports.  A zero-target row (an interaction coefficient against 0) is
// always "vacuous" under claimgate's bound-vs-reference rule, so those rows
// are registered as records, not judged claims; the certificate checks them
// numerically below and says so.
//
// Output layout for `instance all --out DIR`:
//   DIR/three-axis/{axis-cone-table,design-matrices,full-cells,
//                   fractional-cells,coefficients,fractional-contrasts,
//                   alias-matrix,bounds}.csv, certificate.txt
//   DIR/four-axis/{the same files}
//   DIR/instances.csv              one parameterised row per instance
//
// Build (against an installed claimgate 0.1.0 package):
//   cmake -S record -B build/instance -DCMAKE_PREFIX_PATH=<install>/claimgate-0.1.0
//   cmake --build build/instance --config Release
// Run:
//   instance all --out record/results
//   instance three-axis --out record/results
//   instance four-axis --out record/results

#include "claimgate/claimgate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dsi {
namespace {

// ---------------------------------------------------------------------------
// The shared arithmetic quantum and the reported allowance.
// ---------------------------------------------------------------------------

// Q = 2^-40.  Every delivered term is n*Q with |n| < 2^41, so the integer
// numerator and the term are exactly representable in double, and the joint
// integer sum is exact.
constexpr long double kQ = 1.0L / 1099511627776.0L; // 2^-40
constexpr double kQd = 1.0 / 1099511627776.0;       // 2^-40, the double value
constexpr long double kQuantizationMargin = kQ / 2.0L;

// The long-double series evaluation and the platform reference carry a
// rounding error far below the fixed-point quantum; 1e-15 is a conservative
// cover for both at the magnitudes used here.
constexpr long double kEvaluationMargin = 1e-15L;

// The reported high-order allowance.  The cone structure makes the
// structural remainder zero; eta is the measurement floor carried so the
// certificate is robust to the floating-point evaluation of the contrasts.
constexpr long double kEta = 1e-12L;

struct AxisSpec {
    std::string name;
    std::string operation;
    std::string function;
    int terms[2];
    long double remainder[2]; // alternating-series remainder bound per level
    long double exact;        // the exact term the axis is measured against
    long double bound;        // documented per-axis bound, worst over levels
};

std::vector<AxisSpec> MakeAxes(int axisCount) {
    if (axisCount != 3 && axisCount != 4) {
        throw std::invalid_argument("axis count must be 3 or 4");
    }

    std::vector<AxisSpec> a;
    a.push_back({"A1", "exp(-1) Taylor series + fixed-point rounding",
                 "sum_{k>=0} (-1)^k / k!", {5, 8},
                 {1.0L / 120.0L, 1.0L / 40320.0L}, std::exp(-1.0L), 0.0L});
    a.push_back({"A2", "ln 2 alternating harmonic series + fixed-point rounding",
                 "sum_{k>=1} (-1)^(k+1) / k", {8, 12},
                 {1.0L / 9.0L, 1.0L / 13.0L}, std::log(2.0L), 0.0L});
    a.push_back({"A3", "pi/4 arctangent series + fixed-point rounding",
                 "sum_{k>=0} (-1)^k / (2k+1)", {6, 10},
                 {1.0L / 13.0L, 1.0L / 21.0L}, std::atan(1.0L), 0.0L});
    if (axisCount == 4) {
        a.push_back({"A4", "cos(1) Taylor series + fixed-point rounding",
                     "sum_{k>=0} (-1)^k / (2k)!", {4, 6},
                     {1.0L / 40320.0L, 1.0L / 479001600.0L}, std::cos(1.0L), 0.0L});
    }

    for (AxisSpec& axis : a) {
        axis.bound = std::max(axis.remainder[0], axis.remainder[1]) +
                     kQuantizationMargin + kEvaluationMargin;
    }

    return a;
}

long double PartialSum(int axis, int level, const std::vector<AxisSpec>& axes) {
    if (axis < 0 || axis >= static_cast<int>(axes.size())) {
        throw std::out_of_range("partial sum: axis out of range");
    }
    if (level < 0 || level > 1) {
        throw std::out_of_range("partial sum: level out of range");
    }

    const int terms = axes[static_cast<std::size_t>(axis)].terms[level];
    long double sum = 0.0L;

    switch (axis) {
    case 0: {
        long double factorial = 1.0L;
        for (int k = 0; k < terms; ++k) {
            if (k > 0) {
                factorial *= static_cast<long double>(k);
            }
            const long double term = 1.0L / factorial;
            sum += (k % 2 == 0) ? term : -term;
        }
        break;
    }
    case 1: {
        for (int k = 1; k <= terms; ++k) {
            const long double term = 1.0L / static_cast<long double>(k);
            sum += (k % 2 == 1) ? term : -term;
        }
        break;
    }
    case 2: {
        for (int k = 0; k < terms; ++k) {
            const long double term = 1.0L / static_cast<long double>(2 * k + 1);
            sum += (k % 2 == 0) ? term : -term;
        }
        break;
    }
    case 3: {
        long double factorial = 1.0L; // (2k)!
        for (int k = 0; k < terms; ++k) {
            if (k > 0) {
                factorial *= static_cast<long double>(2 * k - 1) *
                             static_cast<long double>(2 * k);
            }
            const long double term = 1.0L / factorial;
            sum += (k % 2 == 0) ? term : -term;
        }
        break;
    }
    default:
        throw std::out_of_range("partial sum: unknown axis");
    }

    return sum;
}

long double JointAdditiveBound(const std::vector<AxisSpec>& axes) {
    long double b = 0.0L;
    for (const AxisSpec& axis : axes) {
        b += axis.bound;
    }
    return b;
}

// ---------------------------------------------------------------------------
// One cell and one design.
// ---------------------------------------------------------------------------

struct CellResult {
    int cell = 0; // design-local cell number
    std::vector<int> level;
    std::vector<long long> n; // fixed-point numerators
    std::vector<long double> delivered;
    std::vector<long double> exactTerm;
    std::vector<long double> m; // delivered - exact, per axis
    double jointDelivered = 0.0;
    long double jointExact = 0.0L;
    long double error = 0.0L;            // joint delivered - joint exact
    long double additiveResidual = 0.0L; // error - sum_i m_i
};

CellResult EvaluateCell(int cell, const std::vector<int>& level,
                        const std::vector<AxisSpec>& axes) {
    if (level.size() != axes.size()) {
        throw std::invalid_argument("cell level count does not match axis count");
    }

    CellResult r;
    r.cell = cell;
    r.level = level;

    const std::size_t nAxes = axes.size();
    r.n.assign(nAxes, 0);
    r.delivered.assign(nAxes, 0.0L);
    r.exactTerm.assign(nAxes, 0.0L);
    r.m.assign(nAxes, 0.0L);

    long long numeratorSum = 0;
    long double jointExact = 0.0L;
    long double sumOfParts = 0.0L;

    for (std::size_t i = 0; i < nAxes; ++i) {
        const long double partial = PartialSum(static_cast<int>(i), level[i], axes);
        const long long n = std::llround(partial / kQ);
        const double delivered = static_cast<double>(n) * kQd;

        r.n[i] = n;
        r.delivered[i] = static_cast<long double>(delivered);
        r.exactTerm[i] = axes[i].exact;
        r.m[i] = r.delivered[i] - r.exactTerm[i];

        numeratorSum += n;
        jointExact += r.exactTerm[i];
        sumOfParts += r.m[i];
    }

    r.jointDelivered = static_cast<double>(numeratorSum) * kQd;
    r.jointExact = jointExact;
    r.error = static_cast<long double>(r.jointDelivered) - jointExact;
    r.additiveResidual = r.error - sumOfParts;

    return r;
}

struct DesignRow {
    int cell = 0;
    std::vector<int> level;
};

std::vector<DesignRow> FullFactorialDesign(int axisCount) {
    const int cells = 1 << axisCount;
    std::vector<DesignRow> rows;
    rows.reserve(static_cast<std::size_t>(cells));

    for (int ci = 0; ci < cells; ++ci) {
        DesignRow row;
        row.cell = ci;
        row.level.resize(static_cast<std::size_t>(axisCount));
        for (int i = 0; i < axisCount; ++i) {
            row.level[static_cast<std::size_t>(i)] =
                (ci >> (axisCount - 1 - i)) & 1;
        }
        rows.push_back(row);
    }

    return rows;
}

// The even-parity half fraction with defining relation I = A1*...*An: the
// cells whose n level bits sum to an even number.  It is an OA(2^(n-1),n,2,2)
// of strength 2: every column is balanced and every pair of columns contains
// all four level combinations equally often.
std::vector<DesignRow> FractionalDesign(int axisCount) {
    const int fullCells = 1 << axisCount;
    std::vector<DesignRow> rows;
    rows.reserve(static_cast<std::size_t>(fullCells / 2));

    int local = 0;
    for (int ci = 0; ci < fullCells; ++ci) {
        std::vector<int> level(static_cast<std::size_t>(axisCount), 0);
        int parity = 0;
        for (int i = 0; i < axisCount; ++i) {
            level[static_cast<std::size_t>(i)] =
                (ci >> (axisCount - 1 - i)) & 1;
            parity ^= level[static_cast<std::size_t>(i)];
        }

        if (parity == 0) {
            DesignRow row;
            row.cell = local++;
            row.level = level;
            rows.push_back(row);
        }
    }

    return rows;
}

std::vector<CellResult> EvaluateDesign(const std::vector<DesignRow>& design,
                                       const std::vector<AxisSpec>& axes) {
    std::vector<CellResult> cells;
    cells.reserve(design.size());

    for (const DesignRow& row : design) {
        cells.push_back(EvaluateCell(row.cell, row.level, axes));
    }

    return cells;
}

// ---------------------------------------------------------------------------
// The two-level contrast basis and the factorial algebra.
// ---------------------------------------------------------------------------

// phi_i,1(level 0) = +1, phi_i,1(level 1) = -1; the mean row is +1.
int Phi(int mask, const std::vector<int>& level) {
    int value = 1;

    for (std::size_t i = 0; i < level.size(); ++i) {
        if ((mask & (1 << static_cast<int>(i))) != 0) {
            value *= 1 - 2 * level[i];
        }
    }

    return value;
}

int Order(int mask) {
    int order = 0;
    while (mask != 0) {
        order += mask & 1;
        mask >>= 1;
    }
    return order;
}

int Components(int axisCount) { return 1 << axisCount; }

int FullCellNumber(const std::vector<int>& level) {
    int full = 0;
    for (const int bit : level) {
        full = (full << 1) | bit;
    }
    return full;
}

std::string SupportName(int mask, const std::vector<AxisSpec>& axes) {
    if (mask == 0) {
        return "{}";
    }

    std::string name = "{";
    bool first = true;

    for (std::size_t i = 0; i < axes.size(); ++i) {
        if ((mask & (1 << static_cast<int>(i))) != 0) {
            if (!first) {
                name += ",";
            }
            name += axes[i].name;
            first = false;
        }
    }

    name += "}";
    return name;
}

std::string UnionName(const std::vector<AxisSpec>& axes) {
    std::string name = "{";
    for (std::size_t i = 0; i < axes.size(); ++i) {
        if (i != 0) {
            name += ",";
        }
        name += axes[i].name;
    }
    name += "}";
    return name;
}

std::vector<long double> Coefficients(const std::vector<CellResult>& cells,
                                      int components) {
    std::vector<long double> c(static_cast<std::size_t>(components), 0.0L);

    for (const CellResult& cell : cells) {
        for (int mask = 0; mask < components; ++mask) {
            c[static_cast<std::size_t>(mask)] +=
                cell.error * static_cast<long double>(Phi(mask, cell.level));
        }
    }

    const long double count = static_cast<long double>(cells.size());

    for (long double& value : c) {
        value /= count;
    }

    return c;
}

using AliasMatrix = std::vector<std::vector<long double>>;

AliasMatrix ComputeAliasMatrix(const std::vector<CellResult>& cells,
                               int components) {
    AliasMatrix omega(static_cast<std::size_t>(components),
                      std::vector<long double>(static_cast<std::size_t>(components),
                                               0.0L));
    const long double count = static_cast<long double>(cells.size());

    for (int a = 0; a < components; ++a) {
        for (int b = 0; b < components; ++b) {
            long double sum = 0.0L;
            for (const CellResult& cell : cells) {
                sum += static_cast<long double>(Phi(a, cell.level) *
                                                Phi(b, cell.level));
            }
            omega[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] =
                sum / count;
        }
    }

    return omega;
}

// The design is strength 2 when every column is balanced and every pair of
// columns contains each of the four level combinations the same number of
// times.
bool IsStrengthTwo(const std::vector<DesignRow>& design, int axisCount) {
    for (int i = 0; i < axisCount; ++i) {
        int counts[2] = {0, 0};
        for (const DesignRow& row : design) {
            ++counts[row.level[static_cast<std::size_t>(i)]];
        }
        if (counts[0] != counts[1]) {
            return false;
        }
    }

    for (int i = 0; i < axisCount; ++i) {
        for (int j = i + 1; j < axisCount; ++j) {
            int counts[2][2] = {{0, 0}, {0, 0}};
            for (const DesignRow& row : design) {
                ++counts[row.level[static_cast<std::size_t>(i)]]
                        [row.level[static_cast<std::size_t>(j)]];
            }
            const int expected = static_cast<int>(design.size()) / 4;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    if (counts[a][b] != expected) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// Formatting and file helpers.
// ---------------------------------------------------------------------------

std::string Num(long double value) {
    std::ostringstream os;
    os << std::setprecision(21) << value;
    return os.str();
}

std::string Num(double value) {
    std::ostringstream os;
    os << std::setprecision(17) << value;
    return os.str();
}

std::string Quote(const std::string& text) {
    std::string out = "\"";
    for (const char c : text) {
        if (c == '"') {
            out += "\"\"";
        } else {
            out += c;
        }
    }
    out += "\"";
    return out;
}

void WriteTextFile(const std::filesystem::path& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("cannot write " + path.string());
    }
    out << text;
}

void WriteCsv(const std::filesystem::path& path,
              const std::vector<std::string>& header,
              const std::vector<std::vector<std::string>>& rows) {
    std::ostringstream os;

    for (std::size_t i = 0; i < header.size(); ++i) {
        if (i != 0) {
            os << ',';
        }
        os << header[i];
    }
    os << '\n';

    for (const std::vector<std::string>& row : rows) {
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i != 0) {
                os << ',';
            }
            os << row[i];
        }
        os << '\n';
    }

    WriteTextFile(path, os.str());
}

std::string FullDesignName(int axisCount) {
    return "full " + std::to_string(1 << axisCount) + "-cell factorial";
}

std::string FractionalDesignName(int axisCount) {
    return "OA(" + std::to_string(1 << (axisCount - 1)) + "," +
           std::to_string(axisCount) + ",2,2) half fraction";
}

std::string JointCallName(int axisCount) {
    std::string name = "joint call (";
    for (int i = 0; i < axisCount; ++i) {
        if (i != 0) {
            name += "+";
        }
        name += "T" + std::to_string(i + 1);
    }
    name += ")";
    return name;
}

std::string CommonCellCsv(const std::string& designName,
                          const std::vector<CellResult>& cells,
                          const std::vector<AxisSpec>& axes,
                          long double jointBound) {
    std::ostringstream os;
    os << "design,cell,full_cell";
    for (std::size_t i = 0; i < axes.size(); ++i) {
        os << ",A" << (i + 1) << "_level";
    }
    for (std::size_t i = 0; i < axes.size(); ++i) {
        os << ",n" << (i + 1);
    }
    for (std::size_t i = 0; i < axes.size(); ++i) {
        os << ",delivered_A" << (i + 1);
    }
    for (std::size_t i = 0; i < axes.size(); ++i) {
        os << ",exact_A" << (i + 1);
    }
    for (std::size_t i = 0; i < axes.size(); ++i) {
        os << ",m_A" << (i + 1);
    }
    os << ",joint_delivered,joint_exact,error,additive_residual,"
          "per_axis_ok,joint_ok\n";

    for (const CellResult& cell : cells) {
        const int fullCell = FullCellNumber(cell.level);

        bool perAxisOk = true;
        for (std::size_t i = 0; i < axes.size(); ++i) {
            perAxisOk = perAxisOk && std::abs(cell.m[i]) <= axes[i].bound;
        }
        const bool jointOk = std::abs(cell.error) <= jointBound;

        os << Quote(designName) << ',' << cell.cell << ',' << fullCell;
        for (const int level : cell.level) {
            os << ',' << level;
        }
        for (const long long n : cell.n) {
            os << ',' << n;
        }
        for (const long double delivered : cell.delivered) {
            os << ',' << Num(delivered);
        }
        for (const long double exact : cell.exactTerm) {
            os << ',' << Num(exact);
        }
        for (const long double m : cell.m) {
            os << ',' << Num(m);
        }
        os << ',' << Num(cell.jointDelivered) << ',' << Num(cell.jointExact)
           << ',' << Num(cell.error) << ',' << Num(cell.additiveResidual)
           << ',' << (perAxisOk ? 1 : 0) << ',' << (jointOk ? 1 : 0) << '\n';
    }

    return os.str();
}

// ---------------------------------------------------------------------------
// The claimgate books.
// ---------------------------------------------------------------------------

int RunAccuracyBook(const std::string& designName,
                    const std::vector<CellResult>& cells,
                    const std::vector<AxisSpec>& axes,
                    std::ostream& out) {
    claimgate::Book book(static_cast<int>(cells.size()) - 1);

    std::vector<claimgate::Claim> axisClaims;
    axisClaims.reserve(axes.size());
    for (const AxisSpec& axis : axes) {
        const std::string lane = axis.name + " series term";
        const std::string region = designName + ", both levels";
        axisClaims.push_back(book.AddClaim(claimgate::PropertyType::Accuracy,
                                           lane, region,
                                           static_cast<double>(axis.bound)));
    }

    const long double jointBound = JointAdditiveBound(axes) + kEta;
    const claimgate::Claim joint = book.AddClaim(
        claimgate::PropertyType::Accuracy, JointCallName(static_cast<int>(axes.size())),
        designName + ", joint additive bound B+eta",
        static_cast<double>(jointBound));

    for (const CellResult& cell : cells) {
        for (std::size_t i = 0; i < axes.size(); ++i) {
            const double got = static_cast<double>(cell.delivered[i]);
            const double ref = static_cast<double>(cell.exactTerm[i]);
            book.Measure(axisClaims[i],
                         claimgate::Cell{cell.cell,
                                         static_cast<double>(cell.level[i]), got,
                                         ref, false},
                         static_cast<double>(axes[i].bound));
        }

        book.Measure(joint,
                     claimgate::Cell{cell.cell, static_cast<double>(cell.cell),
                                     cell.jointDelivered,
                                     static_cast<double>(cell.jointExact), false},
                     static_cast<double>(jointBound));
    }

    out << "\n=== claimgate book: " << designName << " (" << cells.size()
        << " cells) ===\n";
    const int status = book.PrintReport(out);
    out << "book exit status: " << status << "\n";
    return status;
}

// The fractional mean and main-effect contrasts against the full-factorial
// coefficients.  The bound is the allowance eta: the theorem's alias
// contamination is at most G*eta, and G=1 for the two-level +-1 basis.
int RunContrastBook(const std::vector<long double>& full,
                    const std::vector<long double>& fractional,
                    const std::vector<AxisSpec>& axes,
                    std::ostream& out) {
    const int axisCount = static_cast<int>(axes.size());
    claimgate::Book book(axisCount);

    std::vector<claimgate::Claim> claims;
    claims.reserve(static_cast<std::size_t>(axisCount + 1));

    const std::string region =
        FractionalDesignName(axisCount) + " vs " + FullDesignName(axisCount);

    for (int k = 0; k <= axisCount; ++k) {
        const int mask = (k == 0) ? 0 : (1 << (k - 1));
        const std::string lane =
            (k == 0) ? "fractional mean contrast"
                     : "fractional main contrast " + axes[static_cast<std::size_t>(k - 1)].name;

        claims.push_back(book.AddClaim(claimgate::PropertyType::Accuracy, lane,
                                       region, static_cast<double>(kEta)));
        book.Measure(claims.back(),
                     claimgate::Cell{k, 0.0,
                                     static_cast<double>(
                                         fractional[static_cast<std::size_t>(mask)]),
                                     static_cast<double>(
                                         full[static_cast<std::size_t>(mask)]),
                                     false},
                     static_cast<double>(kEta));
    }

    out << "\n=== claimgate book: fractional contrasts against the full factorial "
           "(mean and " << axisCount << " main effects, bound eta) ===\n";
    const int status = book.PrintReport(out);
    out << "book exit status: " << status << "\n";
    return status;
}

// The unsupported-order coefficients are measured against a zero target.
// claimgate's rule counts every such cell as non-discriminating (the bound is
// at least the reference magnitude, 0), so these rows are records, not judged
// claims: a clean row would otherwise be reported "vacuous only" and would not
// count as met.  The record still fails the row if the coefficient is over the
// bound, and the certificate checks the same inequality numerically.
int RunInteractionRecordBook(const std::vector<long double>& full,
                             const std::vector<AxisSpec>& axes,
                             std::ostream& out) {
    const int components = Components(static_cast<int>(axes.size()));
    int recordCount = 0;
    for (int mask = 0; mask < components; ++mask) {
        if (Order(mask) >= 2) {
            ++recordCount;
        }
    }

    claimgate::Book book(recordCount - 1);
    std::vector<claimgate::Claim> claims;
    claims.reserve(static_cast<std::size_t>(recordCount));
    int slot = 0;

    for (int mask = 0; mask < components; ++mask) {
        if (Order(mask) < 2) {
            continue;
        }

        const std::string lane =
            "full-factorial coefficient " + SupportName(mask, axes);
        claims.push_back(book.AddClaim(
            claimgate::PropertyType::Accuracy, lane,
            FullDesignName(static_cast<int>(axes.size())) + ", zero target",
            static_cast<double>(kEta), false));
        book.Measure(claims.back(),
                     claimgate::Cell{slot, 0.0,
                                     static_cast<double>(
                                         full[static_cast<std::size_t>(mask)]),
                                     0.0, false},
                     static_cast<double>(kEta));
        ++slot;
    }

    out << "\n=== claimgate book: unsupported-order coefficients as records "
           "(" << recordCount << " rows, zero target, not judged) ===\n";
    const int status = book.PrintReport(out);
    out << "book exit status: " << status << "\n";
    return status;
}

// ---------------------------------------------------------------------------
// Tables shared by the modes.
// ---------------------------------------------------------------------------

void WriteAxisAndConeTable(const std::filesystem::path& outDir,
                           const std::vector<AxisSpec>& axes) {
    std::vector<std::vector<std::string>> rows;

    for (const AxisSpec& axis : axes) {
        std::ostringstream levels;
        levels << "0=" << axis.terms[0] << " terms, 1=" << axis.terms[1]
               << " terms";
        std::ostringstream remainder;
        remainder << "level0 " << Num(axis.remainder[0]) << "; level1 "
                  << Num(axis.remainder[1]);

        rows.push_back({Quote(axis.name), Quote(axis.operation),
                        Quote(levels.str()), Quote("{" + axis.name + "}"),
                        Quote(remainder.str()), Num(axis.bound),
                        Num(axis.exact)});
    }

    rows.push_back({Quote("op_add"),
                    Quote("exact fixed-point integer addition and 2^-40 scaling"),
                    Quote("fixed, no level"), Quote("{}"),
                    Quote("no local error; exact in double"), "0", Quote("-")});

    WriteCsv(outDir / "axis-cone-table.csv",
             {"axis", "operation", "levels", "cone", "series_remainder_bound",
              "documented_bound", "exact_reference"},
             rows);
}

void WriteDesignMatrices(const std::filesystem::path& outDir,
                         const std::vector<DesignRow>& full,
                         const std::vector<DesignRow>& fractional,
                         int axisCount) {
    std::vector<std::string> header{"design", "cell", "full_cell"};
    for (int i = 0; i < axisCount; ++i) {
        header.push_back("A" + std::to_string(i + 1) + "_level");
    }
    for (int i = 0; i < axisCount; ++i) {
        header.push_back("A" + std::to_string(i + 1) + "_coded");
    }

    std::vector<std::vector<std::string>> rows;

    const auto addRows = [&](const std::string& designName,
                             const std::vector<DesignRow>& design) {
        for (const DesignRow& row : design) {
            std::vector<std::string> csvRow{
                Quote(designName), std::to_string(row.cell),
                std::to_string(FullCellNumber(row.level))};
            for (const int level : row.level) {
                csvRow.push_back(std::to_string(level));
            }
            for (const int level : row.level) {
                csvRow.push_back(std::to_string(1 - 2 * level));
            }
            rows.push_back(csvRow);
        }
    };

    addRows("full", full);
    addRows("fractional", fractional);

    WriteCsv(outDir / "design-matrices.csv", header, rows);
}

void WriteCoefficientTable(const std::filesystem::path& outDir,
                           const std::vector<long double>& full,
                           const std::vector<long double>& fractional,
                           const std::vector<AxisSpec>& axes) {
    const int components = static_cast<int>(full.size());
    const int complement = components - 1;

    std::vector<std::vector<std::string>> rows;
    std::vector<std::vector<std::string>> contrastRows;

    for (int mask = 0; mask < components; ++mask) {
        const int order = Order(mask);
        rows.push_back({std::to_string(mask), std::to_string(order),
                        Quote(SupportName(mask, axes)),
                        Num(full[static_cast<std::size_t>(mask)]),
                        Num(fractional[static_cast<std::size_t>(mask)]),
                        Num(fractional[static_cast<std::size_t>(mask)] -
                            full[static_cast<std::size_t>(mask)])});

        const int partner = mask ^ complement;
        const bool supported = (order <= 1);
        const std::string status =
            supported ? (mask == 0 ? "supported (mean)"
                                   : "supported (main effect)")
                      : "aliased estimate, not a coefficient";
        std::ostringstream identity;
        identity << Num(full[static_cast<std::size_t>(mask)]) << " + "
                 << Num(full[static_cast<std::size_t>(partner)]);

        contrastRows.push_back(
            {std::to_string(mask), std::to_string(order),
             Quote(SupportName(mask, axes)),
             Num(fractional[static_cast<std::size_t>(mask)]),
             Num(full[static_cast<std::size_t>(mask)]),
             Num(fractional[static_cast<std::size_t>(mask)] -
                 full[static_cast<std::size_t>(mask)]),
             std::to_string(partner), Quote(SupportName(partner, axes)),
             identity.str(),
             Num(fractional[static_cast<std::size_t>(mask)] -
                 (full[static_cast<std::size_t>(mask)] +
                  full[static_cast<std::size_t>(partner)])),
             Quote(status)});
    }

    WriteCsv(outDir / "coefficients.csv",
             {"mask", "order", "support", "full_factorial_coefficient",
              "fractional_estimate", "difference"},
             rows);

    WriteCsv(outDir / "fractional-contrasts.csv",
             {"mask", "order", "support", "fractional_estimate",
              "full_coefficient", "difference", "alias_partner",
              "alias_partner_support", "alias_identity",
              "fractional_minus_alias_sum", "status"},
             contrastRows);
}

void WriteAliasMatrix(const std::filesystem::path& outDir,
                      const AliasMatrix& omega,
                      const std::vector<AxisSpec>& axes) {
    const int components = static_cast<int>(axes.size()) == 0
                               ? 0
                               : Components(static_cast<int>(axes.size()));
    std::vector<std::string> header;
    header.push_back("alpha_beta");

    for (int b = 0; b < components; ++b) {
        header.push_back(std::to_string(b) + " " + SupportName(b, axes));
    }

    std::vector<std::vector<std::string>> rows;

    for (int a = 0; a < components; ++a) {
        std::vector<std::string> row;
        row.push_back(std::to_string(a) + " " + SupportName(a, axes));

        for (int b = 0; b < components; ++b) {
            row.push_back(
                Num(omega[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)]));
        }

        rows.push_back(row);
    }

    WriteCsv(outDir / "alias-matrix.csv", header, rows);
}

// ---------------------------------------------------------------------------
// The certificate.
// ---------------------------------------------------------------------------

struct CertificateValues {
    int axisCount = 0;
    int components = 0;
    long double B = 0.0L;
    long double eta = kEta;
    long double certifiedJointBound = 0.0L;
    long double delta = 0.0L;
    long double nonAdditiveBound = 0.0L;
    long double aliasContaminationBound = 0.0L;
    long double maxFullInteraction = 0.0L;
    long double maxMainDifference = 0.0L;
    long double maxMeanDifference = 0.0L;
    long double maxPerAxisRatio = 0.0L;
    long double maxJointRatio = 0.0L;
    long double maxAdditiveResidual = 0.0L;
    int supportedComponents = 0;
    int unsupportedComponents = 0;
    bool strengthTwo = false;
    bool aliasStructureOk = false;
    bool perAxisBoundsOk = false;
    bool jointBoundOk = false;
    bool unsupportedOrdersOk = false;
    bool mainEffectsMatchOk = false;
    bool livenessOk = false;
};

CertificateValues ComputeCertificateValues(
    const std::vector<AxisSpec>& axes, const std::vector<CellResult>& fullCells,
    const std::vector<CellResult>& fractionalCells,
    const std::vector<long double>& full,
    const std::vector<long double>& fractional, const AliasMatrix& omega,
    const std::vector<DesignRow>& fractionalDesign) {
    CertificateValues v;

    v.axisCount = static_cast<int>(axes.size());
    v.components = Components(v.axisCount);
    const int complement = v.components - 1;

    v.B = JointAdditiveBound(axes);
    v.eta = kEta;
    v.certifiedJointBound = v.B + v.eta;

    // delta is the largest supported non-mean contrast; the mean is reported
    // separately.  For w_p=1 the theorem's coefficient bound over 2<=|S|<=w_p
    // is empty and the non-additive bound is eta, but delta is reported as the
    // audit protocol asks.
    for (int mask = 1; mask < v.components; ++mask) {
        if (Order(mask) == 1) {
            v.delta = std::max(
                v.delta, std::abs(fractional[static_cast<std::size_t>(mask)]));
        }
    }

    // w_p=1: N_w=0, G_w=1 (the two-level +-1 basis has sup 1), so the
    // non-additive bound N_w*G_w*(delta+G_w*eta)+eta is eta.
    v.nonAdditiveBound = v.eta;
    v.aliasContaminationBound = v.eta;

    for (int mask = 0; mask < v.components; ++mask) {
        if (Order(mask) >= 2) {
            ++v.unsupportedComponents;
            v.maxFullInteraction = std::max(
                v.maxFullInteraction,
                std::abs(full[static_cast<std::size_t>(mask)]));
        } else {
            ++v.supportedComponents;
        }
    }

    for (int i = 0; i < v.axisCount; ++i) {
        const int mask = 1 << i;
        v.maxMainDifference = std::max(
            v.maxMainDifference,
            std::abs(fractional[static_cast<std::size_t>(mask)] -
                     full[static_cast<std::size_t>(mask)]));
    }
    v.maxMeanDifference = std::abs(fractional[0] - full[0]);

    v.strengthTwo = IsStrengthTwo(fractionalDesign, v.axisCount);

    // Alias structure: for this half fraction the only nonzero entries are
    // omega[a][a]=1 and omega[a][a^complement]=1; every other entry is zero.
    v.aliasStructureOk = true;
    for (int a = 0; a < v.components; ++a) {
        for (int b = 0; b < v.components; ++b) {
            const long double expected =
                (b == a || b == (a ^ complement)) ? 1.0L : 0.0L;
            if (std::abs(
                    omega[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] -
                    expected) > 1e-12L) {
                v.aliasStructureOk = false;
            }
        }
    }

    v.perAxisBoundsOk = true;
    v.jointBoundOk = true;
    v.maxAdditiveResidual = 0.0L;

    for (const std::vector<CellResult>* design :
         {&fullCells, &fractionalCells}) {
        for (const CellResult& cell : *design) {
            for (int i = 0; i < v.axisCount; ++i) {
                const long double ratio =
                    std::abs(cell.m[static_cast<std::size_t>(i)]) /
                    axes[static_cast<std::size_t>(i)].bound;
                v.maxPerAxisRatio = std::max(v.maxPerAxisRatio, ratio);
                if (ratio > 1.0L) {
                    v.perAxisBoundsOk = false;
                }
            }

            const long double jointRatio =
                std::abs(cell.error) / v.certifiedJointBound;
            v.maxJointRatio = std::max(v.maxJointRatio, jointRatio);
            if (jointRatio > 1.0L) {
                v.jointBoundOk = false;
            }

            v.maxAdditiveResidual = std::max(
                v.maxAdditiveResidual, std::abs(cell.additiveResidual));
        }
    }

    v.unsupportedOrdersOk = (v.maxFullInteraction <= v.eta);
    v.mainEffectsMatchOk =
        (v.maxMainDifference <= v.eta) && (v.maxMeanDifference <= v.eta);

    // C3 liveness: both levels of each axis change the delivered term.
    v.livenessOk = true;
    for (int i = 0; i < v.axisCount; ++i) {
        long double low = 0.0L;
        long double high = 0.0L;
        bool haveLow = false;
        bool haveHigh = false;

        for (const CellResult& cell : fullCells) {
            if (cell.level[static_cast<std::size_t>(i)] == 0) {
                low = cell.delivered[static_cast<std::size_t>(i)];
                haveLow = true;
            } else {
                high = cell.delivered[static_cast<std::size_t>(i)];
                haveHigh = true;
            }
        }

        if (!haveLow || !haveHigh || std::abs(high - low) <= 1e-6L) {
            v.livenessOk = false;
        }
    }

    return v;
}

void WriteBoundsTable(const std::filesystem::path& outDir,
                      const std::vector<AxisSpec>& axes,
                      const CertificateValues& v) {
    std::vector<std::vector<std::string>> rows;

    for (const AxisSpec& axis : axes) {
        rows.push_back({"b_" + axis.name, Num(axis.bound),
                        Quote("series remainder at the coarser level + Q/2 + 1e-15")});
    }

    rows.push_back({"axes", std::to_string(v.axisCount), Quote("number of two-level axes")});
    rows.push_back({"components", std::to_string(v.components), Quote("2^axes contrast components")});
    rows.push_back({"w_p", "1", Quote("widest dependency cone on the path")});
    rows.push_back({"U_p", Quote(UnionName(axes)), Quote("union of the path's cones")});
    rows.push_back({"required_strength", "2", Quote("2*w_p")});
    rows.push_back({"B", Num(v.B), Quote("sum of the per-axis bounds")});
    rows.push_back({"eta", Num(v.eta), Quote("reported allowance (measurement floor; structural high-order remainder zero)")});
    rows.push_back({"certified_joint_bound", Num(v.certifiedJointBound), Quote("B + eta")});
    rows.push_back({"delta", Num(v.delta), Quote("largest supported non-mean fractional contrast")});
    rows.push_back({"N_w", "0", Quote("w_p = 1: no supported interaction order")});
    rows.push_back({"G_w", "1", Quote("sup of the two-level +-1 basis row")});
    rows.push_back({"non_additive_bound", Num(v.nonAdditiveBound), Quote("N_w*G_w*(delta+G_w*eta)+eta = eta for w_p=1")});
    rows.push_back({"alias_contamination_bound", Num(v.aliasContaminationBound), Quote("G_w*eta")});
    rows.push_back({"max_full_interaction", Num(v.maxFullInteraction), Quote("largest |full-factorial coefficient| of order >= 2")});
    rows.push_back({"max_main_difference", Num(v.maxMainDifference), Quote("largest |fractional main - full main|")});
    rows.push_back({"max_mean_difference", Num(v.maxMeanDifference), Quote("|fractional mean - full mean|")});
    rows.push_back({"max_per_axis_ratio", Num(v.maxPerAxisRatio), Quote("worst |m_i|/b_i over both designs")});
    rows.push_back({"max_joint_ratio", Num(v.maxJointRatio), Quote("worst |error|/(B+eta) over both designs")});
    rows.push_back({"max_additive_residual", Num(v.maxAdditiveResidual), Quote("largest |error - sum m_i| (floating-point floor)")});
    rows.push_back({"full_cells", std::to_string(1 << v.axisCount), Quote("full factorial")});
    rows.push_back({"fractional_cells", std::to_string(1 << (v.axisCount - 1)), Quote("even-parity half fraction")});
    rows.push_back({"saved_cells", std::to_string(1 << (v.axisCount - 1)), Quote("cells not run")});
    rows.push_back({"saved_percent", "50", Quote("half of the full factorial")});

    WriteCsv(outDir / "bounds.csv", {"quantity", "value", "note"}, rows);
}

std::string CertificateText(const std::string& instanceName,
                            const std::vector<AxisSpec>& axes,
                            const std::vector<long double>& full,
                            const std::vector<long double>& fractional,
                            const CertificateValues& v) {
    std::ostringstream os;

    const std::string header = "d-second-instance certificate: " + instanceName;
    os << header << "\n" << std::string(header.size(), '=') << "\n\n";
    os << "Instance: " << instanceName << " (" << v.axisCount
       << " axes), " << v.axisCount
       << "-term series bank, fixed-point quantum Q = 2^-40.\n"
          "Axes and cones:\n";
    for (const AxisSpec& axis : axes) {
        os << "  " << axis.name << ": " << axis.function
           << "; levels 0=" << axis.terms[0] << " terms, 1=" << axis.terms[1]
           << " terms; cone {" << axis.name << "}; bound b_" << axis.name
           << " = " << Num(axis.bound) << "\n";
    }
    os << "  op_add: exact fixed-point integer addition; no local error, so no\n"
          "          first-order error term and no cone in the error decomposition.\n";
    os << "  w_p = 1; U_p = " << UnionName(axes) << "; |U_p| = " << v.axisCount
       << "; required strength 2*w_p = 2.\n\n";

    const int fullCells = 1 << v.axisCount;
    const int fractionalCells = 1 << (v.axisCount - 1);
    std::string relation = "I = ";
    for (std::size_t i = 0; i < axes.size(); ++i) {
        if (i != 0) {
            relation += "*";
        }
        relation += axes[i].name;
    }

    os << "Design:\n"
       << "  full factorial: " << fullCells << " cells\n"
       << "  fractional: the even-parity half fraction " << relation << ", "
       << fractionalCells << " cells, an " << FractionalDesignName(v.axisCount)
       << "\n"
       << "  strength-2 check (all columns balanced, all pairs balanced): "
       << (v.strengthTwo ? "PASS" : "FAIL") << "\n"
       << "  cells saved: " << fractionalCells << " of " << fullCells
       << " (50.0%)\n\n";

    os << "Bounds:\n";
    for (const AxisSpec& axis : axes) {
        os << "  b_" << axis.name << " = " << Num(axis.bound) << "\n";
    }
    os << "  B = sum b_i = " << Num(v.B) << "\n"
       << "  eta (allowance) = " << Num(v.eta) << "\n"
       << "  certified joint bound B + eta = " << Num(v.certifiedJointBound) << "\n"
       << "  delta (largest supported non-mean contrast) = " << Num(v.delta) << "\n"
       << "  N_w = 0, G_w = 1; non-additive bound = " << Num(v.nonAdditiveBound)
       << " (for w_p=1 the formula reads eta)\n\n";

    os << "Full-factorial coefficients:\n"
          "  mask order support                  coefficient\n";
    for (int mask = 0; mask < v.components; ++mask) {
        os << "  " << mask << "    " << Order(mask) << "     " << std::left
           << std::setw(22) << SupportName(mask, axes) << std::right
           << Num(full[static_cast<std::size_t>(mask)]) << "\n";
    }
    os << "\nFractional estimates:\n"
          "  mask order support                  estimate                    full                        difference\n";
    for (int mask = 0; mask < v.components; ++mask) {
        os << "  " << mask << "    " << Order(mask) << "     " << std::left
           << std::setw(22) << SupportName(mask, axes) << std::right
           << std::setw(28) << Num(fractional[static_cast<std::size_t>(mask)])
           << std::setw(28) << Num(full[static_cast<std::size_t>(mask)])
           << std::setw(28) << Num(fractional[static_cast<std::size_t>(mask)] -
                                   full[static_cast<std::size_t>(mask)])
           << (Order(mask) <= 1 ? ""
                                : "   (aliased estimate, not a coefficient)")
           << "\n";
    }

    const int complement = v.components - 1;
    os << "\nAlias structure:\n"
          "  omega[alpha][alpha] = 1 for every alpha; all other nonzero entries are\n"
          "  the complement pairs:\n";
    for (int mask = 0; mask < v.components / 2; ++mask) {
        os << "    " << SupportName(mask, axes) << " <-> "
           << SupportName(mask ^ complement, axes) << "\n";
    }
    os << "  (the full " << v.components << "x" << v.components
       << " matrix is in alias-matrix.csv)\n";

    os << "\nCertificate checks:\n";
    const auto check = [&os](bool ok, const std::string& text) {
        os << "  " << (ok ? "PASS" : "FAIL") << ": " << text << "\n";
    };

    check(v.perAxisBoundsOk,
          "per-axis bounds hold on all full and fractional cells (worst |m_i|/b_i = " +
              Num(v.maxPerAxisRatio) + ")");
    check(v.jointBoundOk,
          "joint bound B+eta holds on all full and fractional cells (worst |error|/(B+eta) = " +
              Num(v.maxJointRatio) + ")");
    check(v.unsupportedOrdersOk,
          "every unsupported order is zero by the cone structure within eta (largest |full interaction| = " +
              Num(v.maxFullInteraction) + " over " +
              std::to_string(v.unsupportedComponents) + " unsupported components; eta = " +
              Num(v.eta) + ")");
    check(v.mainEffectsMatchOk,
          "fractional main effects (and mean) match the full factorial within eta (largest main difference = " +
              Num(v.maxMainDifference) + ", mean difference = " +
              Num(v.maxMeanDifference) + ")");
    check(v.aliasStructureOk,
          "alias matrix is the complement-pair structure of the half fraction");
    check(v.livenessOk, "both levels of every axis are live (C3)");
    check(v.strengthTwo, "the fractional design has orthogonal-array strength 2");
    check(v.maxAdditiveResidual <= v.eta,
          "measured additive residual is within the allowance (largest |error - sum m_i| = " +
              Num(v.maxAdditiveResidual) + ")");

    const bool allOk = v.perAxisBoundsOk && v.jointBoundOk &&
                       v.unsupportedOrdersOk && v.mainEffectsMatchOk &&
                       v.aliasStructureOk && v.livenessOk && v.strengthTwo &&
                       (v.maxAdditiveResidual <= v.eta);

    os << "\nRESULT: " << (allOk ? "CERTIFIED" : "NOT CERTIFIED")
       << " - full factorial " << fullCells << " cells, fractional "
       << FractionalDesignName(v.axisCount) << " " << fractionalCells
       << " cells, " << fractionalCells << " cells saved (50.0%).\n";

    return os.str();
}

// ---------------------------------------------------------------------------
// Modes.
// ---------------------------------------------------------------------------

struct InstanceSummary {
    std::string name;
    int axisCount = 0;
    int fullCells = 0;
    int fractionalCells = 0;
    int savedCells = 0;
    long double B = 0.0L;
    long double eta = 0.0L;
    long double certifiedJointBound = 0.0L;
    long double delta = 0.0L;
    long double maxFullInteraction = 0.0L;
    long double maxMainDifference = 0.0L;
    long double maxMeanDifference = 0.0L;
    long double maxPerAxisRatio = 0.0L;
    long double maxJointRatio = 0.0L;
    long double maxAdditiveResidual = 0.0L;
    int unsupportedComponents = 0;
    bool strengthTwo = false;
    bool certified = false;
};

int RunFull(const std::filesystem::path& outDir,
            const std::vector<AxisSpec>& axes,
            const std::vector<CellResult>& fullCells) {
    const long double jointBound = JointAdditiveBound(axes) + kEta;
    WriteTextFile(outDir / "full-cells.csv",
                  CommonCellCsv("full", fullCells, axes, jointBound));
    return RunAccuracyBook(FullDesignName(static_cast<int>(axes.size())),
                           fullCells, axes, std::cout);
}

int RunFractional(const std::filesystem::path& outDir,
                  const std::vector<AxisSpec>& axes,
                  const std::vector<CellResult>& fractionalCells) {
    const long double jointBound = JointAdditiveBound(axes) + kEta;
    WriteTextFile(outDir / "fractional-cells.csv",
                  CommonCellCsv("fractional", fractionalCells, axes, jointBound));
    return RunAccuracyBook(FractionalDesignName(static_cast<int>(axes.size())),
                           fractionalCells, axes, std::cout);
}

int RunCertify(const std::filesystem::path& outDir,
               const std::string& instanceName,
               const std::vector<AxisSpec>& axes,
               const std::vector<DesignRow>& fractionalDesign,
               const std::vector<CellResult>& fullCells,
               const std::vector<CellResult>& fractionalCells,
               InstanceSummary* summary) {
    const int components = Components(static_cast<int>(axes.size()));
    const std::vector<long double> full = Coefficients(fullCells, components);
    const std::vector<long double> fractional =
        Coefficients(fractionalCells, components);
    const AliasMatrix omega = ComputeAliasMatrix(fractionalCells, components);

    WriteCoefficientTable(outDir, full, fractional, axes);
    WriteAliasMatrix(outDir, omega, axes);

    const CertificateValues v = ComputeCertificateValues(
        axes, fullCells, fractionalCells, full, fractional, omega,
        fractionalDesign);

    WriteBoundsTable(outDir, axes, v);
    WriteTextFile(outDir / "certificate.txt",
                  CertificateText(instanceName, axes, full, fractional, v));

    std::cout << "\n"
              << CertificateText(instanceName, axes, full, fractional, v);

    const int contrastStatus = RunContrastBook(full, fractional, axes, std::cout);
    const int recordStatus =
        RunInteractionRecordBook(full, axes, std::cout);

    const bool certificateOk =
        v.perAxisBoundsOk && v.jointBoundOk && v.unsupportedOrdersOk &&
        v.mainEffectsMatchOk && v.aliasStructureOk && v.livenessOk &&
        v.strengthTwo && (v.maxAdditiveResidual <= v.eta);

    if (summary != nullptr) {
        summary->name = instanceName;
        summary->axisCount = v.axisCount;
        summary->fullCells = static_cast<int>(fullCells.size());
        summary->fractionalCells = static_cast<int>(fractionalCells.size());
        summary->savedCells = summary->fullCells - summary->fractionalCells;
        summary->B = v.B;
        summary->eta = v.eta;
        summary->certifiedJointBound = v.certifiedJointBound;
        summary->delta = v.delta;
        summary->maxFullInteraction = v.maxFullInteraction;
        summary->maxMainDifference = v.maxMainDifference;
        summary->maxMeanDifference = v.maxMeanDifference;
        summary->maxPerAxisRatio = v.maxPerAxisRatio;
        summary->maxJointRatio = v.maxJointRatio;
        summary->maxAdditiveResidual = v.maxAdditiveResidual;
        summary->unsupportedComponents = v.unsupportedComponents;
        summary->strengthTwo = v.strengthTwo;
        summary->certified = certificateOk;
    }

    std::cout << "\ncertificate checks: " << (certificateOk ? "PASS" : "FAIL")
              << "\n"
              << "contrast book exit status: " << contrastStatus << "\n"
              << "interaction record book exit status: " << recordStatus
              << "\n";

    return (certificateOk && contrastStatus == 0 && recordStatus == 0) ? 0 : 1;
}

void WriteInstanceSummary(const std::filesystem::path& outDir,
                          const std::vector<InstanceSummary>& summaries) {
    std::vector<std::vector<std::string>> rows;

    for (const InstanceSummary& s : summaries) {
        rows.push_back(
            {Quote(s.name), std::to_string(s.axisCount),
             std::to_string(s.fullCells), std::to_string(s.fractionalCells),
             std::to_string(s.savedCells),
             Num(100.0L * static_cast<long double>(s.savedCells) /
                 static_cast<long double>(s.fullCells)),
             "1", "2", Num(s.B), Num(s.eta), Num(s.certifiedJointBound),
             Num(s.delta), Num(s.maxFullInteraction), Num(s.maxMainDifference),
             Num(s.maxMeanDifference), Num(s.maxPerAxisRatio),
             Num(s.maxJointRatio), Num(s.maxAdditiveResidual),
             std::to_string(s.unsupportedComponents),
             s.strengthTwo ? "1" : "0", s.certified ? "1" : "0"});
    }

    WriteCsv(outDir / "instances.csv",
             {"instance", "axes", "full_cells", "fractional_cells",
              "saved_cells", "saved_percent", "w_p", "strength", "B", "eta",
              "certified_joint_bound", "delta", "max_full_interaction",
              "max_main_difference", "max_mean_difference",
              "max_per_axis_ratio", "max_joint_ratio",
              "max_additive_residual", "unsupported_components",
              "strength_two_check", "certified"},
             rows);
}

int RunOneInstance(const std::filesystem::path& outDir, const std::string& name,
                   int axisCount, const std::string& mode,
                   InstanceSummary* summary) {
    std::filesystem::create_directories(outDir);

    const std::vector<AxisSpec> axes = MakeAxes(axisCount);
    const std::vector<DesignRow> fullDesign = FullFactorialDesign(axisCount);
    const std::vector<DesignRow> fractionalDesign = FractionalDesign(axisCount);
    const std::vector<CellResult> fullCells = EvaluateDesign(fullDesign, axes);
    const std::vector<CellResult> fractionalCells =
        EvaluateDesign(fractionalDesign, axes);

    WriteAxisAndConeTable(outDir, axes);
    WriteDesignMatrices(outDir, fullDesign, fractionalDesign, axisCount);

    if (mode == "full") {
        return RunFull(outDir, axes, fullCells);
    }
    if (mode == "fractional") {
        return RunFractional(outDir, axes, fractionalCells);
    }
    if (mode == "certify") {
        return RunCertify(outDir, name, axes, fractionalDesign, fullCells,
                          fractionalCells, summary);
    }
    if (mode == "all") {
        const int fullStatus = RunFull(outDir, axes, fullCells);
        const int fractionalStatus =
            RunFractional(outDir, axes, fractionalCells);
        const int certifyStatus =
            RunCertify(outDir, name, axes, fractionalDesign, fullCells,
                       fractionalCells, summary);
        return (fullStatus == 0 && fractionalStatus == 0 &&
                certifyStatus == 0)
                   ? 0
                   : 1;
    }

    std::cerr << "usage: instance {all|three-axis|four-axis|full|fractional|"
                 "certify} [--out DIR]\n";
    return 2;
}

int Run(const std::string& mode, const std::filesystem::path& outDir) {
    std::filesystem::create_directories(outDir);

    if (mode == "all") {
        std::vector<InstanceSummary> summaries;
        int status = 0;

        for (const int axisCount : {3, 4}) {
            const std::string name =
                (axisCount == 3) ? "three-axis" : "four-axis";
            std::cout
                << "\n================================================================================"
                   "\ninstance: "
                << name << " (" << axisCount << " axes)"
                << "\n================================================================================"
                   "\n";
            InstanceSummary s;
            const int st =
                RunOneInstance(outDir / name, name, axisCount, "all", &s);
            summaries.push_back(s);
            if (st != 0) {
                status = 1;
            }
        }

        WriteInstanceSummary(outDir, summaries);
        std::cout << "\noverall exit status: " << status << "\n";
        return status;
    }

    if (mode == "three-axis" || mode == "four-axis") {
        const int axisCount = (mode == "three-axis") ? 3 : 4;
        std::cout
            << "\n================================================================================"
               "\ninstance: "
            << mode << " (" << axisCount << " axes)"
            << "\n================================================================================"
               "\n";
        return RunOneInstance(outDir / mode, mode, axisCount, "all", nullptr);
    }

    if (mode == "full" || mode == "fractional" || mode == "certify") {
        int status = 0;
        for (const int axisCount : {3, 4}) {
            const std::string name =
                (axisCount == 3) ? "three-axis" : "four-axis";
            std::cout
                << "\n================================================================================"
                   "\ninstance: "
                << name << " (" << axisCount << " axes), mode " << mode
                << "\n================================================================================"
                   "\n";
            const int st =
                RunOneInstance(outDir / name, name, axisCount, mode, nullptr);
            if (st != 0) {
                status = 1;
            }
        }
        return status;
    }

    std::cerr << "usage: instance {all|three-axis|four-axis|full|fractional|"
                 "certify} [--out DIR]\n";
    return 2;
}

} // namespace
} // namespace dsi

int main(int argc, char** argv) {
    std::string mode = (argc > 1) ? argv[1] : "all";
    std::filesystem::path outDir = ".";

    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            outDir = argv[++i];
        }
    }

    try {
        return dsi::Run(mode, outDir);
    } catch (const std::exception& e) {
        std::cerr << "instance error: " << e.what() << "\n";
        return 2;
    }
}
