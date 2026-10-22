// Prints the identifiers the library's own test suite keys its expected-error bounds on, so a
// bound transcribed out of that suite can be the one selected for this configuration.

#include <cstdio>
#include <limits>

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

int main() {
  std::printf("BOOST_LIB_VERSION  = %s\n", BOOST_LIB_VERSION);
  std::printf("BOOST_COMPILER     = %s\n", BOOST_COMPILER);
  std::printf("BOOST_STDLIB       = %s\n", BOOST_STDLIB);
  std::printf("BOOST_PLATFORM     = %s\n", BOOST_PLATFORM);
  std::printf("double digits      = %d\n", std::numeric_limits<double>::digits);
  std::printf("long double digits = %d\n", std::numeric_limits<long double>::digits);
  std::printf("double epsilon     = %.17g\n", std::numeric_limits<double>::epsilon());
  std::printf("boost::math::tools::epsilon<double>() = %.17g\n",
              boost::math::tools::epsilon<double>());
  std::printf("cpp_bin_float_100 digits10 = %d\n",
              std::numeric_limits<boost::multiprecision::cpp_bin_float_100>::digits10);
  return 0;
}
