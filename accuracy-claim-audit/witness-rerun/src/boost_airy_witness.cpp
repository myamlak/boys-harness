/* Boost.Math Airy witness.
 *
 * Reads one 64-bit double pattern per line from the file named on the command
 * line, evaluates boost::math::airy_ai and boost::math::airy_bi in double, and
 * prints the input and both results as bit patterns (authoritative) and as
 * 17-significant-digit decimals.
 *
 * Built header-only against the Boost 1.92.0 release tree with MSVC Release,
 * /fp:precise.  See ../README.md for the build command.
 */
#include <cstdint>
#include <cstdio>
#include <cstring>

#include <boost/math/special_functions/airy.hpp>
#include <boost/version.hpp>

static double bits_to_double(std::uint64_t b) {
  double d;
  std::memcpy(&d, &b, sizeof d);
  return d;
}

static std::uint64_t double_to_bits(double d) {
  std::uint64_t b;
  std::memcpy(&b, &d, sizeof b);
  return b;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s inputs.txt\n", argv[0]);
    return 2;
  }
  std::FILE* f = std::fopen(argv[1], "r");
  if (f == nullptr) {
    std::perror(argv[1]);
    return 2;
  }

  std::printf("# BOOST_VERSION=%d BOOST_LIB_VERSION=%s\n", BOOST_VERSION,
              BOOST_LIB_VERSION);
#ifdef _MSC_VER
  std::printf("# _MSC_VER=%d _MSC_FULL_VER=%d\n", _MSC_VER, _MSC_FULL_VER);
#endif
#if defined(_M_FP_PRECISE)
  std::printf("# _M_FP_PRECISE=1\n");
#elif defined(_M_FP_FAST)
  std::printf("# _M_FP_FAST=1\n");
#elif defined(_M_FP_STRICT)
  std::printf("# _M_FP_STRICT=1\n");
#else
  std::printf("# _M_FP_* not predefined\n");
#endif
  std::printf("# sizeof(double)=%zu FLT_EVAL_METHOD=%d\n", sizeof(double),
              (int)FLT_EVAL_METHOD);
  std::printf("input_bits,x,airy_ai_bits,airy_ai,airy_bi_bits,airy_bi\n");

  char line[256];
  while (std::fgets(line, sizeof line, f) != nullptr) {
    std::uint64_t b = 0;
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
      continue;
    }
    if (std::sscanf(line, "%llx", (unsigned long long*)&b) != 1) {
      continue;
    }
    const double x = bits_to_double(b);
    const double ai = boost::math::airy_ai(x);
    const double bi = boost::math::airy_bi(x);
    std::printf("%016llx,%.17g,%016llx,%.17g,%016llx,%.17g\n",
                (unsigned long long)b, x,
                (unsigned long long)double_to_bits(ai), ai,
                (unsigned long long)double_to_bits(bi), bi);
  }
  std::fclose(f);
  return 0;
}
