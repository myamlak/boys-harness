// Accuracy comparison of the Boys function: boys v5.3.2 double default lanes
// against libint2 2.13.1's double evaluators, on the Boys repository's
// committed reference grids.
//
// Reads a committed grid CSV with columns n,x,value (the gate grid's extra
// columns are ignored).  For every argument x the probe evaluates the full
// order ladder 0..32 once per implementation and compares each (n,x) cell with
// the high-precision reference value.  The reference is parsed into long
// double, never into double, so the reference's own rounding is not charged to
// the implementation under test.
//
// The Boys library's documented double-lane contract is
//     |F_hat - F_n(x)| <= max(1e-13 * F_n(x), 1e-13 * F_0(x))
//                        = 1e-13 * F_0(x)          (F_n <= F_0),
// and the bound is evaluated with the grid's reference F_0(x), not with the
// implementation's own F_0.  Per the specification, the contract binds only
// where the reference value is at or above the double lane's smallest normal
// (2.2250738585072014e-308); cells below that floor are counted separately and
// excluded from pass/fail.
//
// libint2 publishes no guaranteed accuracy bound for its production evaluator,
// FmEval_Chebyshev7<double>.  Its header documents the constructor's precision
// parameter as a desired relative precision defaulting to
// std::numeric_limits<double>::epsilon(), which is a target, not a bound.  Its
// reference evaluator, FmEval_Reference2<double>, has the header statement
// "When Real is double, can maintain absolute precision of epsilon for up to
// m=40" -- an approximate statement, reported here as such.
//
// Build: see record/build-command.txt.  No timing is measured.

#include <boys/boys.hpp>
#include <libint2/boys.h>

#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace {

constexpr int kMaxOrder = 32;
constexpr long double kBoysBoundScale = 1e-13L;
// Smallest normal of binary64 and binary64 epsilon, written exactly.
constexpr long double kDoubleMin = 2.2250738585072014e-308L;
constexpr long double kDoubleEps = 2.2204460492503131e-16L;

struct Row {
    int n = 0;
    long double ref = 0;
    std::string ref_text;
    bool ref_parse_underflow = false;
};

struct Group {
    double x = 0;
    std::vector<Row> rows;
    long double f0 = 0;
    std::string f0_text;
    bool has_f0 = false;
    bool f0_parse_underflow = false;
};

struct CellInfo {
    long double err = 0;
    long double bound = 0;
    long double ratio = 0;
    long double rel = 0;
    int n = 0;
    double x = 0;
    long double ref = 0;
    double value = 0;
    std::string ref_text;
};

struct ImplStats {
    std::string key;
    std::string documented_bound;
    bool has_guaranteed_bound = false;
    bool has_absolute_figure = false;  // libint2 Reference2's epsilon statement
    long double absolute_figure = 0;

    long long rows = 0;
    long long contract_cells = 0;  // reference >= binary64 smallest normal
    long long below_floor = 0;
    long long ref_underflow_long_double = 0;
    long long meet_contract = 0;
    long long over_contract = 0;
    long long meet_all_rows = 0;
    long long meet_absolute_figure = 0;
    long long non_finite = 0;
    long long zero_forbidden = 0;     // returned zero while the true value exceeds its bound
    long long zero_within_bound = 0;  // returned zero while the true value is inside its bound

    CellInfo worst_abs;
    CellInfo worst_ratio;
    CellInfo worst_rel;
    CellInfo worst_absolute_figure_ratio;
    bool have_abs = false;
    bool have_ratio = false;
    bool have_rel = false;
    bool have_absolute_figure_ratio = false;
};

bool HasNonzeroDigit(const std::string& s) {
    for (char c : s) {
        if (c == 'e' || c == 'E') break;
        if (c >= '1' && c <= '9') return true;
    }
    return false;
}

// Parses a decimal literal into long double and reports whether a nonzero
// literal underflowed to zero (the value is below the long double range).
bool ParseLongDouble(const std::string& text, long double* out) {
    errno = 0;
    char* end = nullptr;
    const long double value = std::strtold(text.c_str(), &end);
    if (end == text.c_str()) {
        std::fprintf(stderr, "cannot parse reference value '%s'\n", text.c_str());
        std::exit(2);
    }
    *out = value;
    return std::isinf(value) || (value == 0 && HasNonzeroDigit(text) && errno == ERANGE);
}

std::string Trim(std::string s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' ')) s.pop_back();
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(s.begin());
    return s;
}

std::vector<Group> LoadGrid(const char* path, std::size_t* row_count) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "cannot open grid '%s'\n", path);
        std::exit(2);
    }

    std::vector<Group> groups;
    std::map<double, std::size_t> index;

    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        line = Trim(line);
        if (line.empty()) continue;
        if (header && line.rfind("n,", 0) == 0) {
            header = false;
            continue;
        }
        header = false;

        const std::size_t c1 = line.find(',');
        const std::size_t c2 = c1 == std::string::npos ? std::string::npos : line.find(',', c1 + 1);
        if (c1 == std::string::npos || c2 == std::string::npos) {
            std::fprintf(stderr, "malformed grid line: '%s'\n", line.c_str());
            std::exit(2);
        }

        const std::string n_text = line.substr(0, c1);
        const std::string x_text = line.substr(c1 + 1, c2 - c1 - 1);
        const std::string ref_text = line.substr(c2 + 1, line.find(',', c2 + 1) - c2 - 1);

        char* end = nullptr;
        errno = 0;
        const long n_long = std::strtol(n_text.c_str(), &end, 10);
        if (end == n_text.c_str() || n_long < 0 || n_long > kMaxOrder) {
            std::fprintf(stderr, "bad order '%s'\n", n_text.c_str());
            std::exit(2);
        }
        errno = 0;
        const double x = std::strtod(x_text.c_str(), &end);
        if (end == x_text.c_str() || !std::isfinite(x) || x < 0.0) {
            std::fprintf(stderr, "bad argument '%s'\n", x_text.c_str());
            std::exit(2);
        }

        Row row;
        row.n = static_cast<int>(n_long);
        row.ref_text = ref_text;
        row.ref_parse_underflow = ParseLongDouble(ref_text, &row.ref);

        auto it = index.find(x);
        std::size_t gi;
        if (it == index.end()) {
            gi = groups.size();
            index.emplace(x, gi);
            groups.push_back(Group{});
            groups.back().x = x;
        } else {
            gi = it->second;
        }
        Group& g = groups[gi];
        for (const Row& existing : g.rows) {
            if (existing.n == row.n) {
                std::fprintf(stderr, "duplicate cell n=%d x=%.17g\n", row.n, x);
                std::exit(3);
            }
        }
        g.rows.push_back(row);
        if (row.n == 0) {
            g.f0 = row.ref;
            g.f0_text = row.ref_text;
            g.has_f0 = true;
            g.f0_parse_underflow = row.ref_parse_underflow;
        }
    }

    for (const Group& g : groups) {
        if (!g.has_f0) {
            std::fprintf(stderr, "argument x=%.17g has no n=0 reference\n", g.x);
            std::exit(3);
        }
        if (g.f0_parse_underflow || !(g.f0 > 0)) {
            std::fprintf(stderr, "argument x=%.17g has no representable positive F_0 reference\n", g.x);
            std::exit(3);
        }
        if (g.rows.size() != static_cast<std::size_t>(kMaxOrder + 1)) {
            std::fprintf(stderr, "argument x=%.17g has %zu orders, expected %d\n", g.x, g.rows.size(),
                         kMaxOrder + 1);
            std::exit(3);
        }
    }

    *row_count = groups.size() * static_cast<std::size_t>(kMaxOrder + 1);
    return groups;
}

void NoteCell(const CellInfo& cell, CellInfo& worst, bool& have, bool use_ratio) {
    const long double key = use_ratio ? cell.ratio : cell.err;
    const long double old = use_ratio ? worst.ratio : worst.err;
    if (!have || key > old) {
        worst = cell;
        have = true;
    }
}

void Accumulate(ImplStats& s, const CellInfo& cell, bool in_contract, bool ref_underflow) {
    ++s.rows;
    if (ref_underflow) ++s.ref_underflow_long_double;

    if (!std::isfinite(cell.value)) ++s.non_finite;

    if (cell.err <= cell.bound) ++s.meet_all_rows;

    if (!in_contract) {
        ++s.below_floor;
        return;
    }

    ++s.contract_cells;
    if (cell.value == 0.0) {
        if (cell.ref > cell.bound) {
            ++s.zero_forbidden;
        } else {
            ++s.zero_within_bound;
        }
    }
    if (cell.err <= cell.bound) {
        ++s.meet_contract;
    } else {
        ++s.over_contract;
    }
    if (s.has_absolute_figure && cell.err <= s.absolute_figure) ++s.meet_absolute_figure;

    // A non-finite value is a failure, not a candidate for the worst cell:
    // keep it out of the numeric maxima and let the counter report it.
    if (!std::isfinite(cell.value)) return;

    NoteCell(cell, s.worst_abs, s.have_abs, false);
    NoteCell(cell, s.worst_ratio, s.have_ratio, true);
    if (cell.ref > 0 && (!s.have_rel || cell.rel > s.worst_rel.rel)) {
        s.worst_rel = cell;
        s.have_rel = true;
    }
    if (s.has_absolute_figure) {
        const long double figure_ratio = cell.err / s.absolute_figure;
        if (!s.have_absolute_figure_ratio ||
            figure_ratio > s.worst_absolute_figure_ratio.err) {
            s.worst_absolute_figure_ratio = cell;
            s.worst_absolute_figure_ratio.err = figure_ratio;
            s.have_absolute_figure_ratio = true;
        }
    }
}

void PrintCell(std::FILE* out, const char* label, const CellInfo& c) {
    std::fprintf(out, "%s n=%d x=%.17g ref=%s value=%.17g abs_err=%.17Le rel_err=%.17Le bound=%.17Le "
                      "err_over_bound=%.17Le\n",
                 label, c.n, c.x, c.ref_text.c_str(), c.value, c.err, c.rel, c.bound, c.ratio);
}

void PrintRow(std::FILE* out, const std::string& grid, const ImplStats& s) {
    std::fprintf(out,
                 "%s\t%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%.17Le\t%d\t%.17g\t%s\t%.17g\t"
                 "%.17Le\t%d\t%.17g\t%s\t%.17g\t%.17Le\t%d\t%.17g\t%s\t%.17g\t%s\t",
                 grid.c_str(), s.key.c_str(), s.rows, s.contract_cells, s.below_floor,
                 s.ref_underflow_long_double, s.meet_contract, s.over_contract, s.meet_all_rows,
                 s.non_finite, s.zero_forbidden, s.zero_within_bound,
                 s.worst_abs.err, s.worst_abs.n, s.worst_abs.x, s.worst_abs.ref_text.c_str(),
                 s.worst_abs.value, s.worst_ratio.ratio, s.worst_ratio.n, s.worst_ratio.x,
                 s.worst_ratio.ref_text.c_str(), s.worst_ratio.value, s.worst_rel.rel,
                 s.worst_rel.n, s.worst_rel.x, s.worst_rel.ref_text.c_str(), s.worst_rel.value,
                 s.documented_bound.c_str());
    if (s.has_guaranteed_bound || s.has_absolute_figure) {
        const long double documented_ratio =
            s.has_absolute_figure ? (s.worst_abs.err / s.absolute_figure) : s.worst_ratio.ratio;
        std::fprintf(out, "%.17Le\t", documented_ratio);
    } else {
        std::fprintf(out, "\t");
    }
    if (s.has_absolute_figure) {
        std::fprintf(out, "%lld\n", s.meet_absolute_figure);
    } else if (s.has_guaranteed_bound) {
        std::fprintf(out, "%lld\n", s.meet_contract);
    } else {
        std::fprintf(out, "\n");
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <grid.csv> <grid-label>\n", argv[0]);
        return 2;
    }
    const char* grid_path = argv[1];
    const std::string grid_label = argv[2];

    std::size_t row_count = 0;
    const std::vector<Group> groups = LoadGrid(grid_path, &row_count);

    ImplStats boys_stats;
    boys_stats.key = "boys_v5.3.2_double_single_BoysAllOrders";
    boys_stats.documented_bound = "max(1e-13*F_n(x),1e-13*F_0(x)) = 1e-13*F_0(x)";
    boys_stats.has_guaranteed_bound = true;

    ImplStats boys_batch_stats;
    boys_batch_stats.key = "boys_v5.3.2_double_batch_BoysAllN";
    boys_batch_stats.documented_bound = "max(1e-13*F_n(x),1e-13*F_0(x)) = 1e-13*F_0(x)";
    boys_batch_stats.has_guaranteed_bound = true;

    ImplStats cheb_stats;
    cheb_stats.key = "libint2_2.13.1_FmEval_Chebyshev7_double";
    cheb_stats.documented_bound = "none published";

    ImplStats ref2_stats;
    ref2_stats.key = "libint2_2.13.1_FmEval_Reference2_double";
    ref2_stats.documented_bound = "header statement: absolute precision ~ 1*DBL_EPSILON for m<=40";
    ref2_stats.has_absolute_figure = true;
    ref2_stats.absolute_figure = kDoubleEps;

    libint2::FmEval_Chebyshev7<double> chebyshev(kMaxOrder);

    // The batch lane is called once over every grid argument, in the grid's own
    // (unsorted) order, exactly as its documented contract allows.
    const std::size_t count = groups.size();
    std::vector<double> xs(count);
    for (std::size_t i = 0; i < count; ++i) xs[i] = groups[i].x;
    std::vector<double> batch_values(count * static_cast<std::size_t>(kMaxOrder + 1));
    std::vector<std::size_t> workspace(boys::BoysAllNWorkspaceSize(count));
    boys::BoysAllN(kMaxOrder, xs.data(), batch_values.data(), count, workspace.data());

    std::array<double, kMaxOrder + 1> boys_values{};
    std::array<double, kMaxOrder + 1> cheb_values{};
    std::array<double, kMaxOrder + 1> ref2_values{};
    long long batch_differences = 0;

    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
        const Group& g = groups[gi];
        boys::BoysAllOrders(kMaxOrder, g.x, boys_values.data());

        chebyshev.eval(cheb_values.data(), g.x, kMaxOrder);
        libint2::FmEval_Reference2<double>::eval(ref2_values.data(), g.x, kMaxOrder);

        const long double f0 = g.f0;
        for (const Row& row : g.rows) {
            const std::size_t n = static_cast<std::size_t>(row.n);

            CellInfo cell;
            cell.n = row.n;
            cell.x = g.x;
            cell.ref = row.ref;
            cell.ref_text = row.ref_text;
            cell.bound = kBoysBoundScale * f0;
            const bool in_contract = !row.ref_parse_underflow && row.ref >= kDoubleMin;

            if (boys_values[n] != batch_values[n * count + gi]) ++batch_differences;

            cell.value = boys_values[n];
            cell.err = std::fabs(static_cast<long double>(cell.value) - row.ref);
            cell.ratio = cell.bound > 0 ? cell.err / cell.bound : 0;
            cell.rel = row.ref > 0 ? cell.err / row.ref : 0;
            Accumulate(boys_stats, cell, in_contract, row.ref_parse_underflow);

            cell.value = batch_values[n * count + gi];
            cell.err = std::fabs(static_cast<long double>(cell.value) - row.ref);
            cell.ratio = cell.bound > 0 ? cell.err / cell.bound : 0;
            cell.rel = row.ref > 0 ? cell.err / row.ref : 0;
            Accumulate(boys_batch_stats, cell, in_contract, row.ref_parse_underflow);

            cell.value = cheb_values[n];
            cell.err = std::fabs(static_cast<long double>(cell.value) - row.ref);
            cell.ratio = cell.bound > 0 ? cell.err / cell.bound : 0;
            cell.rel = row.ref > 0 ? cell.err / row.ref : 0;
            Accumulate(cheb_stats, cell, in_contract, row.ref_parse_underflow);

            cell.value = ref2_values[n];
            cell.err = std::fabs(static_cast<long double>(cell.value) - row.ref);
            cell.ratio = cell.bound > 0 ? cell.err / cell.bound : 0;
            cell.rel = row.ref > 0 ? cell.err / row.ref : 0;
            Accumulate(ref2_stats, cell, in_contract, row.ref_parse_underflow);
        }
    }

    std::printf("# grid: %s\n", grid_path);
    std::printf("# grid label: %s\n", grid_label.c_str());
    std::printf("# rows: %zu | unique x: %zu | orders: 0..%d\n", row_count, groups.size(), kMaxOrder);
    std::printf("# boys single vs batch differing cells: %lld\n", batch_differences);
    std::printf("# contract domain: reference value >= DBL_MIN (2.2250738585072014e-308)\n");
    std::printf("# columns: grid, impl, rows, contract_cells, below_double_min, ref_underflow, "
                "meet_contract, over_contract, meet_all_rows, non_finite, zero_forbidden, "
                "zero_within_bound, worst_abs_err, worst_abs_n, worst_abs_x, "
                "worst_abs_ref, worst_abs_value, worst_bound_ratio, worst_ratio_n, worst_ratio_x, "
                "worst_ratio_ref, worst_ratio_value, worst_rel_err, worst_rel_n, worst_rel_x, "
                "worst_rel_ref, worst_rel_value, documented_bound, worst_over_documented_bound, "
                "documented_bound_cells_met\n");
    PrintRow(stdout, grid_label, boys_stats);
    PrintRow(stdout, grid_label, boys_batch_stats);
    PrintRow(stdout, grid_label, cheb_stats);
    PrintRow(stdout, grid_label, ref2_stats);

    std::printf("# worst cells\n");
    PrintCell(stdout, "# boys single worst_abs", boys_stats.worst_abs);
    PrintCell(stdout, "# boys single worst_ratio", boys_stats.worst_ratio);
    PrintCell(stdout, "# boys single worst_rel", boys_stats.worst_rel);
    PrintCell(stdout, "# boys batch worst_abs", boys_batch_stats.worst_abs);
    PrintCell(stdout, "# boys batch worst_ratio", boys_batch_stats.worst_ratio);
    PrintCell(stdout, "# boys batch worst_rel", boys_batch_stats.worst_rel);
    PrintCell(stdout, "# cheb7 worst_abs", cheb_stats.worst_abs);
    PrintCell(stdout, "# cheb7 worst_ratio_to_boys_contract", cheb_stats.worst_ratio);
    PrintCell(stdout, "# cheb7 worst_rel", cheb_stats.worst_rel);
    PrintCell(stdout, "# ref2 worst_abs", ref2_stats.worst_abs);
    PrintCell(stdout, "# ref2 worst_ratio_to_boys_contract", ref2_stats.worst_ratio);
    PrintCell(stdout, "# ref2 worst_rel", ref2_stats.worst_rel);
    std::printf("# reference cells below the double floor: %lld (of which long-double-underflow: %lld)\n",
                boys_stats.below_floor, boys_stats.ref_underflow_long_double);

    return 0;
}
