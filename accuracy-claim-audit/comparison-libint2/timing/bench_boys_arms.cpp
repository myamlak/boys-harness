// Wall-clock comparison of the Boys kernel's double-precision entries with
// libint2 2.13.1's double evaluators on the committed gate grid.
//
// Four arms, all evaluated at orders 0..32 for every argument:
//   boys ladder   boys::BoysAllOrders
//   boys batch    boys::BoysAllN
//   cheb7         libint2::FmEval_Chebyshev7<double>
//   reference2    libint2::FmEval_Reference2<double>
//
// Each arm is warmed up twice, then timed five times.  The median is reported.
// No accuracy claim is made here.  Build with -O3.

#include <boys/boys.hpp>
#include <libint2/boys.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr int kMaxOrder = 32;
constexpr int kRepeats = 5;

double Now() {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

std::vector<double> LoadUniqueX(const char* path) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", path);
        std::exit(2);
    }
    std::vector<double> xs;
    std::set<unsigned long long> seen;
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (header) { header = false; continue; }
        const char* p = line.c_str();
        char* end = nullptr;
        std::strtod(p, &end);
        if (end == p) continue;
        p = end;
        while (*p == ',' || *p == ' ') ++p;
        const double x = std::strtod(p, &end);
        if (end == p) continue;
        unsigned long long bits = 0;
        std::memcpy(&bits, &x, sizeof bits);
        if (seen.insert(bits).second) xs.push_back(x);
    }
    return xs;
}

double Median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <grid.csv>\n", argv[0]);
        return 2;
    }
    const std::vector<double> xs = LoadUniqueX(argv[1]);
    const std::size_t count = xs.size();
    std::printf("# arguments: %zu, orders: 0..%d\n", count, kMaxOrder);

    std::array<double, kMaxOrder + 1> out{};
    std::vector<double> batch(count * (kMaxOrder + 1));
    std::vector<std::size_t> workspace(boys::BoysAllNWorkspaceSize(count));
    libint2::FmEval_Chebyshev7<double> chebyshev(kMaxOrder);

    double sink = 0.0;
    auto ladder = [&]() {
        for (double x : xs) {
            boys::BoysAllOrders(kMaxOrder, x, out.data());
            sink += out[kMaxOrder];
        }
    };
    auto batch_call = [&]() {
        boys::BoysAllN(kMaxOrder, xs.data(), batch.data(), count,
                       workspace.data());
        for (std::size_t i = 0; i < count; ++i) sink += batch[i];
    };
    auto cheb7 = [&]() {
        for (double x : xs) {
            chebyshev.eval(out.data(), x, kMaxOrder);
            sink += out[kMaxOrder];
        }
    };
    auto reference2 = [&]() {
        for (double x : xs) {
            libint2::FmEval_Reference2<double>::eval(out.data(), x, kMaxOrder);
            sink += out[kMaxOrder];
        }
    };

    // Two warm-up passes per arm.
    ladder(); ladder();
    batch_call(); batch_call();
    cheb7(); cheb7();
    reference2(); reference2();

    auto time_arm = [&](auto&& fn, const char* name) {
        std::vector<double> times;
        for (int r = 0; r < kRepeats; ++r) {
            const double t0 = Now();
            fn();
            times.push_back(Now() - t0);
        }
        const double med = Median(times);
        const double ns_arg = med * 1e9 / count;
        const double ns_cell = ns_arg / (kMaxOrder + 1);
        std::printf("%-12s total %.6f s  ns/argument %10.2f  ns/order %7.2f\n",
                    name, med, ns_arg, ns_cell);
    };

    time_arm(ladder, "boys_ladder");
    time_arm(batch_call, "boys_batch");
    time_arm(cheb7, "cheb7");
    time_arm(reference2, "reference2");

    std::printf("# checksum %.17g\n", sink);
    return 0;
}
