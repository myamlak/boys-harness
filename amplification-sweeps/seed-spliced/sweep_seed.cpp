// Seed-spliced sweep: the upward recurrence seeded by the direct entry's own
// order-0 value.
//
// The paper's per-order comparison confounds two things: the recurrence's
// amplification of a seed error, and the different stored fits the two entries
// read.  This driver removes the second: it keeps the direct entry's own F_0
// as the seed and steps it upward through the kernel's own recurrence body, so
// the seed error and the direct error come from the same fit family.
//
// For every argument of the second sweep's grid this driver emits
//
//   * the shipped recursed entry (the paper's policy, from its band seed),
//   * the direct matrix-product entry, read once per band, and
//   * the spliced ladder: F_0 = the direct entry's order-0 value, then
//     F_l = boys::detail::UpwardStep<PaperPolicy::kDivision>(l, F_{l-1}, ...)
//     for l = 1, ..., 32.
//
// Every step is also run through the replica body of ../common/replica.hpp and
// the two are compared bit for bit, so "the recurrence this sweep ran" is the
// kernel's own body and the replica is shown to agree with it everywhere.
// The replica admission of each shipped ladder is the same check the first two
// sweeps record.
//
// x = 0 is present in the grid and has no upward recurrence (the step divides
// by x); its spliced column is written `nan` and the analysis excludes it.
//
// Outputs:
//   cells.csv          band,x,n,shipped,direct,spliced
//   admission.csv      band,x,direction,checked,matched
//   replica_steps.csv  band,x,steps,matched

#include "../common/replica.hpp"

#include "boys/boys_fit_core.hpp"
#include "boys/boys_transform.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace {

using paper_a::kNmax;
using paper_a::PaperPolicy;
using paper_a::SameBits;

struct Row {
    std::string band;
    double x;
};

std::string Argument(int argc, char** argv, const char* name, const char* fallback) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::string(argv[i]) == name) {
            return argv[i + 1];
        }
    }
    return fallback;
}

int IntegerArgument(int argc, char** argv, const char* name, int fallback) {
    const std::string text = Argument(argc, argv, name, "");
    return text.empty() ? fallback : std::atoi(text.c_str());
}

std::vector<Row> ReadGrid(const std::string& path) {
    std::FILE* file = std::fopen(path.c_str(), "r");

    if (!file) {
        std::fprintf(stderr, "cannot open %s\n", path.c_str());
        std::exit(2);
    }

    std::vector<Row> rows;
    char line[512];
    std::fgets(line, sizeof line, file); // header

    while (std::fgets(line, sizeof line, file)) {
        const std::string text(line);
        const std::size_t first = text.find(',');
        const std::size_t second = text.find(',', first + 1);

        if (first == std::string::npos || second == std::string::npos) {
            continue;
        }

        rows.push_back({text.substr(0, first), std::stod(text.substr(first + 1, second - first - 1))});
    }

    std::fclose(file);
    return rows;
}

std::FILE* Open(const std::string& path, const char* header) {
    std::FILE* file = std::fopen(path.c_str(), "w");

    if (!file) {
        std::fprintf(stderr, "cannot write %s\n", path.c_str());
        std::exit(2);
    }

    std::fputs(header, file);
    return file;
}

} // namespace

int main(int argc, char** argv) {
    const std::string grid_path = Argument(argc, argv, "--grid", "data/grid.csv");
    const std::string cells_path = Argument(argc, argv, "--cells", "results/cells.csv");
    const std::string admission_path = Argument(argc, argv, "--admission", "results/admission.csv");
    const std::string steps_path = Argument(argc, argv, "--steps", "results/replica_steps.csv");
    const int limit = IntegerArgument(argc, argv, "--max-arguments", 0);

    std::vector<Row> rows = ReadGrid(grid_path);

    if (limit > 0 && static_cast<int>(rows.size()) > limit) {
        rows.resize(static_cast<std::size_t>(limit));
    }

    struct Band {
        const char* name;
        boys::RegionABand region;
        std::vector<Row> rows;
    };

    std::array<Band, 2> bands = {{
        {"A1", boys::RegionABand::kA1, {}},
        {"A2", boys::RegionABand::kA2, {}},
    }};

    for (const Row& row : rows) {
        for (Band& band : bands) {
            if (row.band == band.name) {
                band.rows.push_back(row);
            }
        }
    }

    std::FILE* cells = Open(cells_path, "band,x,n,shipped,direct,spliced\n");
    std::FILE* admission = Open(admission_path, "band,x,direction,checked,matched\n");
    std::FILE* steps = Open(steps_path, "band,x,steps,matched\n");

    for (Band& band : bands) {
        const std::size_t count = band.rows.size();

        if (count == 0) {
            continue;
        }

        std::vector<double> xs(count);
        for (std::size_t i = 0; i < count; ++i) {
            xs[i] = band.rows[i].x;
        }

        std::vector<double> direct(((kNmax + 1) * count), 0.0);
        boys::BoysRegionAProduct<boys::ProductMode::kFp64>(
            band.region, kNmax, xs.data(), direct.data(), count);

        for (std::size_t i = 0; i < count; ++i) {
            const double x = xs[i];
            std::array<double, kNmax + 1> shipped{};
            boys::BoysAllOrders<PaperPolicy>(kNmax, x, shipped.data());

            const paper_a::Admission admitted = paper_a::AdmitReplica(x, shipped.data());
            std::fprintf(admission, "%s,%.17g,%s,%d,%d\n", band.name, x, admitted.direction,
                         admitted.checked, admitted.matched);

            std::array<double, kNmax + 1> spliced{};
            int step_count = 0;
            int step_matched = 0;

            if (x != 0.0) {
                const double seed = direct[i]; // order 0, same layout as the first sweep
                const double half_exp = 0.5 * std::exp(-x);
                const double reciprocal = 1.0 / x;
                double f = seed;
                spliced[0] = f;

                for (int l = 1; l <= kNmax; ++l) {
                    const double kernel =
                        boys::detail::UpwardStep<PaperPolicy::kDivision>(l, f, x, reciprocal, half_exp);
                    const double replica = paper_a::UpwardStep(l, f, x, reciprocal, half_exp);
                    ++step_count;
                    step_matched += SameBits(kernel, replica) ? 1 : 0;
                    f = kernel;
                    spliced[static_cast<std::size_t>(l)] = f;
                }
            } else {
                const double nan = std::numeric_limits<double>::quiet_NaN();
                for (double& value : spliced) {
                    value = nan;
                }
            }

            std::fprintf(steps, "%s,%.17g,%d,%d\n", band.name, x, step_count, step_matched);

            for (int n = 0; n <= kNmax; ++n) {
                const double direct_value =
                    direct[static_cast<std::size_t>(n) * count + i];
                std::fprintf(cells, "%s,%.17g,%d,%.17g,%.17g,", band.name, x, n, shipped[n],
                             direct_value);

                if (std::isfinite(spliced[static_cast<std::size_t>(n)])) {
                    std::fprintf(cells, "%.17g\n", spliced[static_cast<std::size_t>(n)]);
                } else {
                    std::fputs("nan\n", cells);
                }
            }
        }
    }

    std::fclose(cells);
    std::fclose(admission);
    std::fclose(steps);
    return 0;
}
