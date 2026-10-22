// A claim gate over a third-party library's documented accuracy claims.
//
// The gate reads a row table: each row is one documented accuracy claim written down in
// machine-readable form, carrying a function, a domain, a bound, and how the domain is
// sampled.  It measures every sampled point of every row against a high-precision reference,
// gives each row exactly one verdict, and fails while any row is unmet.  Its headline reports
// how many rows were registered and how many of its comparison cells can discriminate.
//
// The rule for a cell is the one the method uses: a cell can discriminate when the row's bound,
// expressed in the units the comparison is made at, is at least eight ulps of the value the
// comparison is made at -- above the rounding floor of double precision, where a breach of the
// bound and the arithmetic that computed the value are distinguishable.  Cells below that floor
// are counted, reported, and not judged either way.
//
// The metric is the one the library's own test suite uses, so that a bound transcribed out of
// that suite means the same thing here: relative_difference(computed, reference) / epsilon,
// which is the "errors in units of epsilon" the library's tables and bounds are stated in.
//
// Exit status: 0 if no row is unmet, 1 if any row is unmet, 2 if the table could not be read.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/version.hpp>

namespace {

using ref_t = boost::multiprecision::cpp_bin_float_100;
const double kEps = std::numeric_limits<double>::epsilon();

double ulp_of(double v) {
  v = std::fabs(v);
  return std::nextafter(v, std::numeric_limits<double>::infinity()) - v;
}

// One documented claim, machine-readable.  The first seven fields are the method's row; the
// rest carry the provenance of the claim and are printed so that no verdict can be read without
// the words it came from.
struct Row {
  std::string id;
  std::string function;
  double param = 0.0;      // second argument, where the function takes one (Bessel order)
  std::string quantity;    // "eps_relative" or "absolute"
  double lo = 0.0;
  double hi = 0.0;
  double bound = 0.0;      // in units of epsilon for eps_relative rows
  std::string sampling;    // "linear:n" | "geometric:n:first" | "data:<file>"
  std::string source;      // "manual" | "test-bound"
  std::string interpreted; // "no" | a short note on what had to be decided
  std::string quote;
};

std::vector<std::string> split(const std::string& line, char sep) {
  std::vector<std::string> out;
  std::string item;
  std::istringstream in(line);
  while (std::getline(in, item, sep)) out.push_back(item);
  return out;
}

std::vector<Row> read_rows(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open the row table '" + path + "'");
  std::vector<Row> rows;
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (header) {
      header = false;
      continue;
    }
    const std::vector<std::string> f = split(line, '\t');
    if (f.size() < 12) throw std::runtime_error("malformed row: " + line);
    Row r;
    r.id = f[0];
    r.function = f[1];
    r.param = std::stod(f[2]);
    r.quantity = f[3];
    r.lo = std::stod(f[4]);
    r.hi = std::stod(f[5]);
    r.bound = std::stod(f[6]);
    r.sampling = f[7];
    r.source = f[8];
    r.interpreted = f[9];
    r.quote = f[11];
    rows.push_back(r);
  }
  if (rows.empty()) throw std::runtime_error("the row table '" + path + "' has no rows");
  return rows;
}

// ---- the library under verification, in double and at 100 decimal digits -------------------

double evaluate_double(const std::string& f, double nu, double x) {
  if (f == "erf") return boost::math::erf(x);
  if (f == "erfc") return boost::math::erfc(x);
  if (f == "tgamma") return boost::math::tgamma(x);
  if (f == "lgamma") return boost::math::lgamma(x);
  if (f == "digamma") return boost::math::digamma(x);
  if (f == "cyl_bessel_j") return boost::math::cyl_bessel_j(nu, x);
  if (f == "cyl_bessel_i") return boost::math::cyl_bessel_i(nu, x);
  if (f == "expint") return boost::math::expint(x);
  if (f == "airy_ai") return boost::math::airy_ai(x);
  if (f == "airy_bi") return boost::math::airy_bi(x);
  if (f == "zeta") return boost::math::zeta(x);
  throw std::runtime_error("row table names an unknown function: " + f);
}

ref_t evaluate_ref(const std::string& f, double nu, double x) {
  const ref_t rx(x);
  if (f == "erf") return boost::math::erf(rx);
  if (f == "erfc") return boost::math::erfc(rx);
  if (f == "tgamma") return boost::math::tgamma(rx);
  if (f == "lgamma") return boost::math::lgamma(rx);
  if (f == "digamma") return boost::math::digamma(rx);
  if (f == "cyl_bessel_j") return boost::math::cyl_bessel_j(ref_t(nu), rx);
  if (f == "cyl_bessel_i") return boost::math::cyl_bessel_i(ref_t(nu), rx);
  if (f == "expint") return boost::math::expint(rx);
  if (f == "airy_ai") return boost::math::airy_ai(rx);
  if (f == "airy_bi") return boost::math::airy_bi(rx);
  if (f == "zeta") return boost::math::zeta(rx);
  throw std::runtime_error("row table names an unknown function: " + f);
}

// A sampled point: the argument, and the order where the function takes one.  A row whose
// domain is a group of published test points in two variables carries no single order, and
// takes it from the point itself.
struct Sample {
  double x = 0.0;
  double nu = 0.0;
};

std::vector<Sample> sample_points(const Row& row) {
  std::vector<Sample> out;
  const std::vector<std::string> s = split(row.sampling, ':');
  if (s[0] == "linear" || s[0] == "geometric") {
    std::vector<double> xs;
    if (s[0] == "linear") {
      const int n = std::stoi(s[1]);
      for (int i = 0; i < n; ++i)
        xs.push_back(row.lo + (row.hi - row.lo) * i / static_cast<double>(n - 1));
    } else {
      const int n = std::stoi(s[1]);
      const double first = std::stod(s[2]);
      xs.push_back(row.lo);
      const double ratio = row.hi / first;
      for (int i = 0; i < n - 1; ++i)
        xs.push_back(first * std::pow(ratio, i / static_cast<double>(n - 2)));
    }
    for (double x : xs) out.push_back(Sample{x, row.param});
  } else if (s[0] == "data") {
    std::ifstream in(s[1]);
    if (!in) throw std::runtime_error("cannot open the sampled data file '" + s[1] + "'");
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      const std::vector<std::string> c = split(line, ',');
      const double x = std::stod(c[0]);
      if (x < row.lo || x > row.hi) continue;
      double nu = row.param;
      if (std::isnan(row.param)) {
        if (c.size() < 2) throw std::runtime_error("row " + row.id + " wants an order per point");
        nu = std::stod(c[1]);
      }
      out.push_back(Sample{x, nu});
    }
  } else {
    throw std::runtime_error("unknown sampling plan: " + row.sampling);
  }
  return out;
}

struct RefDump {
  std::string function;
  double nu = 0.0;
  double x = 0.0;
  double value = 0.0;
  std::string reference;
};

std::vector<RefDump> g_dump;
bool g_dump_enabled = false;

struct Cell {
  double x = 0.0;
  double nu = 0.0;
  double observed_eps = 0.0;  // relative_difference(computed, reference) / epsilon
  double scale = 0.0;         // the magnitude the relative difference is taken against
  double floor_eps = 0.0;     // the method's rounding floor, in the same units
  bool judged = false;
};

struct RowResult {
  const Row* row = nullptr;
  std::vector<Cell> cells;
  int judged = 0;
  double worst_eps = 0.0;
  double worst_x = 0.0;
  double worst_exceeding = 0.0;
  double worst_exceeding_x = 0.0;
  bool unmet = false;
  bool unverifiable = false;
};

RowResult evaluate(const Row& row) {
  RowResult result;
  result.row = &row;
  for (const Sample& s : sample_points(row)) {
    Cell c;
    c.x = s.x;
    c.nu = s.nu;
    const double got = evaluate_double(row.function, s.nu, s.x);
    const ref_t ref = evaluate_ref(row.function, s.nu, s.x);
    const double ref_d = static_cast<double>(ref);
    c.scale = std::fmin(std::fabs(got), std::fabs(ref_d));
    c.observed_eps = static_cast<double>(boost::math::relative_difference(got, ref) / kEps);
    // The method's floor: eight ulps of the value the comparison is made at, carried into the
    // row's units.  For a value of magnitude v, a last-bit change of the result moves the
    // relative difference by about ulp(v)/v, so the floor in epsilon units is that ratio times
    // eight.
    c.floor_eps = (c.scale > 0.0) ? 8.0 * ulp_of(ref_d) / (c.scale * kEps)
                                  : std::numeric_limits<double>::infinity();
    c.judged = row.bound >= c.floor_eps;
    if (c.judged) ++result.judged;
    if (c.observed_eps > result.worst_eps) {
      result.worst_eps = c.observed_eps;
      result.worst_x = c.x;
    }
    if (c.observed_eps > row.bound && (c.judged || c.observed_eps > c.floor_eps) &&
        c.observed_eps > result.worst_exceeding) {
      result.worst_exceeding = c.observed_eps;
      result.worst_exceeding_x = c.x;
    }
    result.cells.push_back(c);
    if (g_dump_enabled) {
      // The reference exactly as the gate used it, so an outside implementation can be asked
      // whether it agrees with these digits.
      std::ostringstream os;
      os << std::setprecision(50) << ref;
      g_dump.push_back(RefDump{row.function, s.nu, s.x, got, os.str()});
    }
  }
  result.unmet = result.worst_exceeding > 0.0;
  result.unverifiable = (result.judged == 0);
  return result;
}

const char* verdict_of(const RowResult& r) {
  if (r.unmet) return "UNMET";
  if (r.unverifiable) return "unverifiable";
  return "met";
}

}  // namespace

int main(int argc, char** argv) {
  std::string rows_path = "claims/rows.csv";
  std::string refout;
  bool verbose = false;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg.rfind("--rows=", 0) == 0) rows_path = arg.substr(7);
    if (arg.rfind("--refout=", 0) == 0) refout = arg.substr(9);
    if (arg == "--cells") verbose = true;
  }

  std::printf("claim gate over a third-party library's documented accuracy claims\n");
  std::printf("  library     : Boost %s (Boost.Math)\n", BOOST_LIB_VERSION);
  std::printf("  build       : %s | %s | %s | double, %d mantissa bits\n", BOOST_COMPILER,
              BOOST_STDLIB, BOOST_PLATFORM, std::numeric_limits<double>::digits);
  std::printf("  reference   : boost::multiprecision::cpp_bin_float_100, %d decimal digits\n",
              std::numeric_limits<ref_t>::max_digits10);
  std::printf("  metric      : relative_difference(value, reference) / epsilon, epsilon = %.3e\n",
              kEps);
  std::printf("  table       : %s\n\n", rows_path.c_str());

  std::vector<Row> rows;
  try {
    rows = read_rows(rows_path);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "gate: %s\n", e.what());
    return 2;
  }

  g_dump_enabled = !refout.empty();
  std::vector<RowResult> results;
  for (const Row& row : rows) results.push_back(evaluate(row));
  if (g_dump_enabled) {
    std::ofstream out(refout);
    if (!out) {
      std::fprintf(stderr, "gate: cannot write the reference dump '%s'\n", refout.c_str());
      return 2;
    }
    out << "function,order,x,value_double,reference\n";
    for (const RefDump& d : g_dump)
      out << d.function << "," << d.nu << "," << std::setprecision(17) << d.x << ","
          << std::setprecision(17) << d.value << "," << d.reference << "\n";
  }

  int met = 0, unmet = 0, unverifiable = 0, cells = 0, judged = 0;
  for (const RowResult& r : results) {
    if (r.unmet) ++unmet;
    else if (r.unverifiable) ++unverifiable;
    else ++met;
    cells += static_cast<int>(r.cells.size());
    judged += r.judged;
  }
  std::printf("  claims      : %d registered   (met %d, unmet %d, unverifiable %d)\n",
              static_cast<int>(results.size()), met, unmet, unverifiable);
  std::printf("  cells       : %d compared, %d can discriminate (%.0f%%)\n\n", cells, judged,
              cells ? 100.0 * judged / cells : 0.0);

  std::printf(" %-4s %-14s %-13s %-20s %9s %6s %7s %10s %8s  %s\n", "id", "function",
              "quantity", "domain", "bound", "cells", "judged", "worst", "bound at", "verdict");
  bool any_unmet = false;
  for (const RowResult& r : results) {
    const Row& row = *r.row;
    char domain[64];
    std::snprintf(domain, sizeof(domain), "[%g, %g]", row.lo, row.hi);
    std::printf(" %-4s %-14s %-13s %-20s %9.4g %6d %7d %10.4g %8s  %s\n", row.id.c_str(),
                row.function.c_str(), row.quantity.c_str(), domain, row.bound,
                static_cast<int>(r.cells.size()), r.judged, r.worst_eps,
                row.source == "test-bound" ? "eps" : "eps", verdict_of(r));
    if (r.unmet) {
      any_unmet = true;
      std::printf("      breach: x = %g, observed %.4g eps against the bound %.4g eps\n",
                  r.worst_exceeding_x, r.worst_exceeding, row.bound);
      std::printf("      claim : %s\n", row.quote.c_str());
    }
    if (verbose) {
      for (const Cell& c : r.cells)
        std::printf("      cell x=%-14g observed=%-12.4g floor=%-10.4g scale=%-12.4g %s\n", c.x,
                    c.observed_eps, c.floor_eps, c.scale, c.judged ? "judged" : "below floor");
    }
  }

  std::printf("\n  input: the gate reads the row table and nothing else. A documented accuracy\n");
  std::printf("         claim that was never registered as a row is not checked, and is not\n");
  std::printf("         reported as unchecked: the gate has no way to know that it exists.\n\n");
  std::printf("GATE: %s\n", any_unmet ? "FAIL (at least one row is unmet)"
                                        : "PASS (no row is unmet)");
  return any_unmet ? 1 : 0;
}
