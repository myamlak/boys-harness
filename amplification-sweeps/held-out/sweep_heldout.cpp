// Held-out sweep: the shipped recursed entry and the direct matrix-product
// entry on an argument set disjoint from both earlier sweeps.
//
// The two evaluations are the second sweep's, unchanged: one matrix-product
// read per band, one all-orders call per argument, and the bitwise replica
// admission.  Nothing here knows the prediction rule; the analysis applies
// the rule frozen from tables/prediction-record.tex to these cells.
//
// Outputs:
//   cells.csv      band,x,n,recursed,direct
//   admission.csv  band,x,direction,checked,matched

#include "../common/replica.hpp"

#include "boys/boys_transform.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

using paper_a::kNmax;
using paper_a::PaperPolicy;

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

    const std::vector<Row> rows = ReadGrid(grid_path);

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

    std::FILE* cells = Open(cells_path, "band,x,n,recursed,direct\n");
    std::FILE* admission = Open(admission_path, "band,x,direction,checked,matched\n");

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

            for (int n = 0; n <= kNmax; ++n) {
                std::fprintf(cells, "%s,%.17g,%d,%.17g,%.17g\n", band.name, x, n, shipped[n],
                             direct[static_cast<std::size_t>(n) * count + i]);
            }
        }
    }

    std::fclose(cells);
    std::fclose(admission);
    return 0;
}
