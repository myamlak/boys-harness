// The platform exponential at every tested argument, with the sweep's own
// flags and compiler.  Reads "xbits" lines (16 hex digits), writes
// "xbits expbits halfbits" where exp = std::exp(-x) and half = 0.5 * exp.
#include <cmath>
#include <cstdio>
#include <cstring>

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::fprintf(stderr, "usage: exp_probe <in> <out>\n");
        return 2;
    }

    std::FILE* in = std::fopen(argv[1], "r");
    std::FILE* out = std::fopen(argv[2], "w");
    if (in == nullptr || out == nullptr)
    {
        std::fprintf(stderr, "cannot open input or output\n");
        return 2;
    }

    unsigned long long xb = 0;
    while (std::fscanf(in, "%llx", &xb) == 1)
    {
        double x = 0.0;
        std::memcpy(&x, &xb, sizeof x);
        const double e = std::exp(-x);
        const double h = 0.5 * std::exp(-x);
        unsigned long long eb = 0;
        unsigned long long hb = 0;
        std::memcpy(&eb, &e, sizeof e);
        std::memcpy(&hb, &h, sizeof h);
        std::fprintf(out, "%016llx %016llx %016llx\n", xb, eb, hb);
    }

    std::fclose(in);
    std::fclose(out);
    return 0;
}
