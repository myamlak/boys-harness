// The Paper-A chemistry bridge probe (scratch copy only).
//
// One executable, built twice: once with the MD engine's fp64 Boys ladder on
// the library's host fp64 all-orders class row (the default arm) and once on
// Paper A's measured recursed policy (QCX_BRIDGE_BOYS_POLICY=paper-a). Both
// arms run the same H2O/STO-3G calculation here:
//
//   1. the one-electron matrices S, T, V (no Boys call in this tree reaches
//      the bridge switch, so they are route-independent by construction);
//   2. the dense ERI tensor through the MD engine, with every Boys argument
//      and the ladder it returned captured in MdBoysBatch;
//   3. the captured argument stream re-evaluated under BOTH policies, so the
//      offline reference sees the two ladders on exactly the run's arguments;
//   4. RHF twice - the library's default tolerances and a tight gate - with
//      the default run's density dump giving the route-independent starting
//      density (the GWH guess is built from S and H only).
//
// Everything is written as text at 17 significant digits, which round-trips a
// double exactly, so the offline analysis reads the same bits the engine
// computed. The one binary file is the SCF density dump, parsed here.
//
// Usage: qcx-bridge-probe <out-prefix> [basis-directory]
//
// With no basis directory the probe runs the in-tree H2O/STO-3G fixture. With
// one, it parses that NWChem basis directory (filtered to the water molecule's
// elements, Z = 1 and 8) and runs the same water geometry in that basis. The
// molecule is always the fixture's published geometry, so a basis series
// differs in the basis alone.

// Paper A's policy reads the coarse split-Clenshaw Chebyshev table, which the
// library's front door does not open. It must precede boys/boys.hpp in the
// translation unit that instantiates the policy (the fit's static_assert names
// the missing include otherwise).
#include "boys/boys_coefficients_d_coarse_cheb.hpp"

#include "internal/md_boys.hpp"

#include "benzene_sto3g.hpp"
#include "h2o_sto3g.hpp"

#include "qcx/integrals/boys.hpp"
#include "qcx/integrals/eri_dense.hpp"
#include "qcx/integrals/one_electron.hpp"
#include "qcx/integrals/screening.hpp"
#include "qcx/integrals/shell_pairs.hpp"
#include "qcx/molecule/mass_properties.hpp"
#include "qcx/scf/rhf.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

using Tensor2 = qcx::memory::Tensor<double, 2, qcx::backend::CpuTag>;
using Tensor4 = qcx::memory::Tensor<double, 4, qcx::backend::CpuTag>;

/// The comparison policy: the host fp64 all-orders class row, the combination
/// an unnamed BoysAllOrders call resolves to. The static asserts below pin the
/// seven cells to the table's row, so a library revision that moves the row
/// fails this probe's build rather than silently comparing other policies.
using BoysClassRowPolicy = boys::DefaultPolicy<boys::Precision::kFp64, boys::Shape::kAllOrders>;

/// Paper A's measured recursed policy, named as the paper's option row states
/// it: Chebyshev fits, split Clenshaw summation, the coarsest partition, the
/// refined reciprocal, the fast region-B exponential. The budget cell is inert
/// on the double lane and the arguments packing axis does not change the
/// scalar all-orders body this shape runs; the other five cells are the
/// measured ones.
using BoysPaperAPolicy =
    boys::StatedEvalPolicy<boys::FitRoute::kChebyshev,
                           boys::EvalScheme::kSplitClenshaw,
                           boys::BoysBudget::kFloat,
                           boys::PackAxis::kArguments,
                           boys::FitGranularity::kCoarsest,
                           boys::DivisionForm::kRefinedReciprocal,
                           boys::RegionBExp::kFast>;

static_assert(boys::EvalPolicyLike<BoysClassRowPolicy>);
static_assert(boys::EvalPolicyLike<BoysPaperAPolicy>);

// The class row, cell by cell (boys_build_defaults.hpp, the kHost/kFp64/
// kAllOrders row).
static_assert(BoysClassRowPolicy::kRoute == boys::FitRoute::kChebyshev);
static_assert(BoysClassRowPolicy::kScheme == boys::EvalScheme::kHorner);
static_assert(BoysClassRowPolicy::kBudget == boys::BoysBudget::kFloat);
static_assert(BoysClassRowPolicy::kPack == boys::PackAxis::kOrders);
static_assert(BoysClassRowPolicy::kGranularity == boys::FitGranularity::kUniform);
static_assert(BoysClassRowPolicy::kDivision == boys::DivisionForm::kRefinedReciprocal);
static_assert(BoysClassRowPolicy::kRegionBExp == boys::RegionBExp::kFast);

// Paper A's row, cell by cell.
static_assert(BoysPaperAPolicy::kRoute == boys::FitRoute::kChebyshev);
static_assert(BoysPaperAPolicy::kScheme == boys::EvalScheme::kSplitClenshaw);
static_assert(BoysPaperAPolicy::kPack == boys::PackAxis::kArguments);
static_assert(BoysPaperAPolicy::kGranularity == boys::FitGranularity::kCoarsest);
static_assert(BoysPaperAPolicy::kDivision == boys::DivisionForm::kRefinedReciprocal);
static_assert(BoysPaperAPolicy::kRegionBExp == boys::RegionBExp::kFast);

constexpr std::size_t kLadderSize = static_cast<std::size_t>(boys::kMaxBoysOrder) + 1;

#if defined(QCX_BRIDGE_BOYS_POLICY_PAPER_A)
constexpr const char* kEnginePolicy = "paper-a";
#else
constexpr const char* kEnginePolicy = "default-class-row";
#endif

#if defined(QCX_BRIDGE_BOYS_CAPTURE)
constexpr bool kCaptureEnabled = true;
#else
constexpr bool kCaptureEnabled = false;
#endif

[[noreturn]] void Fail(const std::string& what) {
    std::cerr << "qcx-bridge-probe: " << what << '\n';
    std::exit(1);
}

template <typename T> T& Require(qcx::Result<T>& result, const char* what) {
    if (!result.has_value())
    {
        Fail(std::string(what) + ": " + result.error().message);
    }

    return *result;
}

Eigen::MatrixXd ToEigen(const Tensor2& tensor) {
    const auto shape = tensor.Shape();
    Eigen::MatrixXd matrix(static_cast<Eigen::Index>(shape[0]),
                           static_cast<Eigen::Index>(shape[1]));

    for (std::size_t i = 0; i < shape[0]; ++i)
    {
        for (std::size_t j = 0; j < shape[1]; ++j)
        {
            matrix(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = tensor(i, j);
        }
    }

    return matrix;
}

std::ofstream Open(const std::string& path) {
    std::ofstream stream(path, std::ios::trunc);

    if (!stream)
    {
        Fail("cannot open " + path);
    }

    stream << std::setprecision(17);
    return stream;
}

void WriteMatrix(std::ostream& out, const char* name, const Eigen::MatrixXd& matrix) {
    out << "matrix " << name << ' ' << matrix.rows() << ' ' << matrix.cols() << '\n';

    for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    {
        for (Eigen::Index col = 0; col < matrix.cols(); ++col)
        {
            out << matrix(row, col) << (col + 1 == matrix.cols() ? '\n' : ' ');
        }
    }
}

/// The two-electron energy and the core-Hamiltonian trace at one fixed
/// density. The closed-shell two-electron energy has a Coulomb and an
/// exchange part, in the index order of the SCF's own supermatrices:
///
///   E_2e = 1/2 sum D_uv D_ls (uv|ls) - 1/4 sum D_uv D_ls (ul|sv).
///
/// Both sums are compensated, so the quantity is not limited by the
/// accumulation of 2,401 terms, and both arms run this same loop over the same
/// density, so their difference is the ERI response and nothing else. The
/// trace identity against HfResult.electronicEnergy is checked in the summary;
/// it closes to the last bit because the exchange term is included.
struct FixedDensityEnergy {
    double coulomb = 0.0;
    double exchange = 0.0;
    double twoElectron = 0.0; // coulomb + exchange
    double coreTrace = 0.0;
    double total = 0.0;
};

FixedDensityEnergy EnergyAtDensity(const Tensor4& eri,
                                   const Eigen::MatrixXd& density,
                                   const Eigen::MatrixXd& core,
                                   double nuclearRepulsion,
                                   std::size_t n) {
    double coulomb = 0.0;
    double coulombCompensation = 0.0;
    double exchange = 0.0;
    double exchangeCompensation = 0.0;

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            const double dij = density(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));

            for (std::size_t k = 0; k < n; ++k)
            {
                for (std::size_t l = 0; l < n; ++l)
                {
                    const double dkl =
                        density(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(l));

                    const double coulombTerm = 0.5 * dij * dkl * eri(i, j, k, l);
                    const double coulombY = coulombTerm - coulombCompensation;
                    const double coulombT = coulomb + coulombY;
                    coulombCompensation = (coulombT - coulomb) - coulombY;
                    coulomb = coulombT;

                    const double exchangeTerm = -0.25 * dij * dkl * eri(i, k, l, j);
                    const double exchangeY = exchangeTerm - exchangeCompensation;
                    const double exchangeT = exchange + exchangeY;
                    exchangeCompensation = (exchangeT - exchange) - exchangeY;
                    exchange = exchangeT;
                }
            }
        }
    }

    double coreTrace = 0.0;

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            coreTrace += density(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) *
                         core(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
        }
    }

    FixedDensityEnergy energy;
    energy.coulomb = coulomb;
    energy.exchange = exchange;
    energy.twoElectron = coulomb + exchange;
    energy.coreTrace = coreTrace;
    energy.total = nuclearRepulsion + coreTrace + energy.twoElectron;
    return energy;
}

[[maybe_unused]] void EvaluateClassRow(int nmax, double x, double* out) {
    boys::BoysAllOrders<BoysClassRowPolicy>(nmax, x, out);
}

[[maybe_unused]] void EvaluatePaperA(int nmax, double x, double* out) {
    boys::BoysAllOrders<BoysPaperAPolicy>(nmax, x, out);
}

std::uint64_t ReadLeU64(std::istream& in) {
    unsigned char bytes[8] = {};
    in.read(reinterpret_cast<char*>(bytes), 8);
    std::uint64_t value = 0;

    for (int i = 0; i < 8; ++i)
    {
        value |= static_cast<std::uint64_t>(bytes[i]) << (8 * i);
    }

    return value;
}

/// The first density record of an SCF density dump: iteration 0's density is
/// the starting guess, built from the overlap and core Hamiltonian alone, so
/// it is the same matrix in both arms bit for bit.
Eigen::MatrixXd ReadSeedDensity(const std::string& path, std::size_t expectedN) {
    std::ifstream in(path, std::ios::binary);

    if (!in)
    {
        Fail("cannot open the density dump " + path);
    }

    char magic[8] = {};
    in.read(magic, 8);

    if (!in || std::string(magic, 8) != "QXCDFDMP")
    {
        Fail("the density dump magic is missing: " + path);
    }

    const std::uint64_t n = ReadLeU64(in);

    if (!in || n != expectedN)
    {
        Fail("the density dump dimension does not match the basis");
    }

    for (;;)
    {
        unsigned char tag = 0;
        in.read(reinterpret_cast<char*>(&tag), 1);

        if (!in)
        {
            break;
        }

        const std::uint64_t iteration = ReadLeU64(in);
        const std::uint64_t dimension = ReadLeU64(in);
        static_cast<void>(iteration);

        if (!in || dimension != expectedN)
        {
            Fail("a density dump record is malformed");
        }

        Eigen::MatrixXd matrix(static_cast<Eigen::Index>(dimension),
                               static_cast<Eigen::Index>(dimension));

        for (std::uint64_t row = 0; row < dimension; ++row)
        {
            for (std::uint64_t col = 0; col < dimension; ++col)
            {
                const std::uint64_t bits = ReadLeU64(in);
                double value = 0.0;
                std::memcpy(&value, &bits, sizeof(value));
                matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) = value;
            }
        }

        if (!in)
        {
            Fail("the density dump is truncated");
        }

        // Tag 2 is a density record; the first one is iteration 0, the seed.
        if (tag == 2)
        {
            return matrix;
        }
    }

    Fail("the density dump holds no density record");

    return Eigen::MatrixXd{};
}

/// Sorts the capture's one-record-per-argument output into (nmax, x, xBits)
/// order: the order the argument file is written in, so a rerun of the same
/// sweep reproduces it.
void SortCapturedArguments(std::vector<qcx::integrals::internal::BoysCaptureRecord>& records) {
    std::sort(records.begin(), records.end(), [](const auto& a, const auto& b) {
        if (a.nmax != b.nmax)
        {
            return a.nmax < b.nmax;
        }

        if (a.x != b.x)
        {
            return a.x < b.x;
        }

        return a.xBits < b.xBits;
    });
}

struct ScfArm {
    bool converged = false;
    int iterations = 0;
    double totalEnergy = 0.0;
    double electronicEnergy = 0.0;
    double achievedEnergyDelta = 0.0;
    double achievedRmsDensityChange = 0.0;
    Eigen::MatrixXd density;
};

ScfArm RunScf(const qcx::molecule::Molecule& molecule,
              const Eigen::MatrixXd& overlap,
              const Eigen::MatrixXd& core,
              const Tensor4& eri,
              qcx::scf::RhfOptions options) {
    auto result = qcx::scf::RunRhfScf(molecule, overlap, core, eri, options);

    if (!result.has_value())
    {
        Fail("RunRhfScf: " + result.error().message);
    }

    ScfArm arm;
    arm.converged = result->converged;
    arm.iterations = result->iterations;
    arm.totalEnergy = result->totalEnergy;
    arm.electronicEnergy = result->electronicEnergy;
    arm.achievedEnergyDelta = result->achievedEnergyDelta;
    arm.achievedRmsDensityChange = result->achievedRmsDensityChange;
    arm.density = result->density;
    return arm;
}

void WriteScfSummary(std::ostream& out, const char* name, const ScfArm& arm) {
    out << name << "_converged=" << (arm.converged ? 1 : 0) << '\n';
    out << name << "_iterations=" << arm.iterations << '\n';
    out << name << "_total_energy=" << arm.totalEnergy << '\n';
    out << name << "_electronic_energy=" << arm.electronicEnergy << '\n';
    out << name << "_achieved_energy_delta=" << arm.achievedEnergyDelta << '\n';
    out << name << "_achieved_rms_density_change=" << arm.achievedRmsDensityChange << '\n';
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 4)
    {
        std::cerr << "usage: qcx-bridge-probe <out-prefix> [basis-directory] [water|benzene]\n";
        return 2;
    }

    const std::string prefix = argv[1];
    const std::string basisDirectory = argc >= 3 ? argv[2] : std::string{};
    const bool externalBasis = argc >= 3 && !basisDirectory.empty();
    const std::string moleculeName = argc == 4 ? argv[3] : "water";
    const bool benzene = moleculeName == "benzene";

    if (moleculeName != "water" && !benzene)
    {
        std::cerr << "unknown molecule '" << moleculeName << "': water or benzene\n";
        return 2;
    }

    auto moleculeResult = benzene ? qcx::testing::MakeBenzeneSto3g() : qcx::testing::MakeH2oSto3g();
    const qcx::molecule::Molecule& molecule =
        Require(moleculeResult, benzene ? "MakeBenzeneSto3g" : "MakeH2oSto3g");

    qcx::Result<qcx::basisset::BasisSet> basisResult =
        externalBasis
            ? qcx::basisset::ParseNwchemDirectoryFiltered(
                  basisDirectory, benzene ? std::array<int, 2>{1, 6} : std::array<int, 2>{1, 8})
            : (benzene ? qcx::testing::MakeBenzeneSto3gBasis()
                       : qcx::testing::MakeH2oSto3gBasis());
    const qcx::basisset::BasisSet& basis =
        Require(basisResult,
                externalBasis ? "ParseNwchemDirectoryFiltered"
                              : (benzene ? "MakeBenzeneSto3gBasis" : "MakeH2oSto3gBasis"));

    auto overlapResult = qcx::integrals::BuildOverlapMatrix(molecule, basis);
    auto kineticResult = qcx::integrals::BuildKineticMatrix(molecule, basis);
    auto nuclearResult = qcx::integrals::BuildNuclearAttractionMatrix(molecule, basis);

    const Eigen::MatrixXd overlap = ToEigen(Require(overlapResult, "BuildOverlapMatrix"));
    const Eigen::MatrixXd kinetic = ToEigen(Require(kineticResult, "BuildKineticMatrix"));
    const Eigen::MatrixXd nuclear = ToEigen(Require(nuclearResult, "BuildNuclearAttractionMatrix"));
    const Eigen::MatrixXd core = kinetic + nuclear;
    const std::size_t n = static_cast<std::size_t>(overlap.rows());
    const double nuclearRepulsion = qcx::molecule::NuclearRepulsionEnergy(molecule);

    qcx::integrals::EriDenseOptions eriOptions;

    if (benzene)
    {
        // The dense build admits by its own size estimate; benzene/cc-pVDZ is
        // 2.52 GiB against the 2 GiB default, so the molecule the probe was
        // extended for raises the cap it needs.
        eriOptions.maxTensorBytes = 4ull * 1024 * 1024 * 1024;
    }

    auto eriResult = qcx::integrals::BuildEriTensorGeneral(molecule, basis, eriOptions);
    const Tensor4& eri = Require(eriResult, "BuildEriTensorGeneral");

#if defined(QCX_BRIDGE_BOYS_CAPTURE)
    std::size_t capturedCalls = 0;
    std::size_t ladderMismatches = 0;
    std::vector<qcx::integrals::internal::BoysCaptureRecord> unique =
        qcx::integrals::internal::BoysCaptureTake(&capturedCalls, &ladderMismatches);
    SortCapturedArguments(unique);
#else
    const std::size_t ladderMismatches = 0;
#endif

    // The shell/pair/quartet census, and the screened survivor count under the
    // dense builder's own defaults (kNormal, screen on). It runs after the
    // capture is taken, so the sweep's own diagonal quartets - the same
    // arguments the ERI build already formed - never double-count in the
    // recorded stream.
    std::size_t shells = 0;
    std::size_t pairs = 0;
    std::size_t canonicalQuartets = 0;
    std::size_t screenedQuartets = 0;
    std::size_t pairFunctionCount = 0;
    {
        auto pairList = qcx::integrals::BuildShellPairs(molecule, basis);

        if (!pairList.has_value())
        {
            Fail("BuildShellPairs: " + pairList.error().message);
        }

        pairFunctionCount = pairList->functionCount;

        auto bounds = qcx::integrals::ComputeSchwarzBounds(molecule, basis);

        if (!bounds.has_value())
        {
            Fail("ComputeSchwarzBounds: " + bounds.error().message);
        }

        shells = pairList->shells.size();
        pairs = pairList->pairs.size();
        const double threshold = qcx::integrals::SchwarzThreshold(eriOptions.accuracy);

        for (std::size_t bra = 0; bra < pairs; ++bra)
        {
            for (std::size_t ket = 0; ket <= bra; ++ket)
            {
                ++canonicalQuartets;

                if ((*bounds)[bra] * (*bounds)[ket] >= threshold)
                {
                    ++screenedQuartets;
                }
            }
        }
    }

    std::size_t captureCalls = 0;
    std::size_t captureUnique = 0;
    std::array<std::size_t, kLadderSize> callsByOrder{};
    std::array<std::size_t, kLadderSize> uniqueByOrder{};

#if defined(QCX_BRIDGE_BOYS_CAPTURE)
    captureCalls = capturedCalls;
    captureUnique = unique.size();

    for (const auto& argument : unique)
    {
        callsByOrder[static_cast<std::size_t>(argument.nmax)] += argument.count;
    }

    for (const auto& argument : unique)
    {
        ++uniqueByOrder[static_cast<std::size_t>(argument.nmax)];
    }

    // The captured stream under both policies. The engine's own ladder is
    // reproduced by the arm it was built with; the other column is the policy
    // the other arm would have run. Writing both here means one variant's
    // output already carries the pair, and the two variants cross-check each
    // other.
    {
        std::ofstream args = Open(prefix + ".args.txt");
        args << "# nmax count x engine[0..nmax] class-row[0..nmax] paper-a[0..nmax]\n";

        std::array<double, kLadderSize> classRow{};
        std::array<double, kLadderSize> paperA{};

        for (const auto& argument : unique)
        {
            classRow.fill(0.0);
            paperA.fill(0.0);
            EvaluateClassRow(argument.nmax, argument.x, classRow.data());
            EvaluatePaperA(argument.nmax, argument.x, paperA.data());

            args << argument.nmax << ' ' << argument.count << ' ' << argument.x;

            for (int l = 0; l <= argument.nmax; ++l)
            {
                args << ' ' << argument.ladder[static_cast<std::size_t>(l)];
            }

            for (int l = 0; l <= argument.nmax; ++l)
            {
                args << ' ' << classRow[static_cast<std::size_t>(l)];
            }

            for (int l = 0; l <= argument.nmax; ++l)
            {
                args << ' ' << paperA[static_cast<std::size_t>(l)];
            }

            args << '\n';
        }
    }

    unique.clear();
    unique.shrink_to_fit();
#endif

    // The ERI dump: every value, in the tensor's own index order. 17 digits is
    // a round-trip, so the offline difference is the engine's difference.
    {
        std::ofstream eriOut = Open(prefix + ".eri.txt");
        eriOut << n << '\n';
        double maxAbs = 0.0;

        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                for (std::size_t k = 0; k < n; ++k)
                {
                    for (std::size_t l = 0; l < n; ++l)
                    {
                        const double value = eri(i, j, k, l);
                        maxAbs = std::max(maxAbs, std::abs(value));
                        eriOut << i << ' ' << j << ' ' << k << ' ' << l << ' ' << value << '\n';
                    }
                }
            }
        }

        eriOut.close();
        std::ofstream maxOut = Open(prefix + ".eri_max.txt");
        maxOut << maxAbs << '\n';
    }

    // The two RHF runs. The default-tolerance run carries the density dump;
    // its iteration-0 density is the GWH seed, built from S and H alone and so
    // identical in both arms.
    qcx::scf::RhfOptions defaultOptions;
    defaultOptions.densityDumpFile = prefix + ".default.density.bin";
    defaultOptions.traceFile = prefix + ".default.trace.txt";
    const ScfArm defaultArm = RunScf(molecule, overlap, core, eri, defaultOptions);

    qcx::scf::RhfOptions tightOptions;
    tightOptions.maxIterations = 500;
    tightOptions.energyTolerance = 1e-12;
    tightOptions.densityTolerance = 1e-10;
    tightOptions.traceFile = prefix + ".tight.trace.txt";
    const ScfArm tightArm = RunScf(molecule, overlap, core, eri, tightOptions);

    const Eigen::MatrixXd seedDensity = ReadSeedDensity(defaultOptions.densityDumpFile, n);

    const FixedDensityEnergy seedEnergy =
        EnergyAtDensity(eri, seedDensity, core, nuclearRepulsion, n);
    const FixedDensityEnergy defaultEnergy =
        EnergyAtDensity(eri, defaultArm.density, core, nuclearRepulsion, n);
    const FixedDensityEnergy tightEnergy =
        EnergyAtDensity(eri, tightArm.density, core, nuclearRepulsion, n);

    // Matrices: S, H, the seed density and both converged densities. The
    // offline analysis cross-evaluates each density against both arms' ERIs.
    {
        std::ofstream matrices = Open(prefix + ".matrices.txt");
        WriteMatrix(matrices, "overlap", overlap);
        WriteMatrix(matrices, "core", core);
        WriteMatrix(matrices, "seed_density", seedDensity);
        WriteMatrix(matrices, "default_density", defaultArm.density);
        WriteMatrix(matrices, "tight_density", tightArm.density);
    }


    {
        std::ofstream summary = Open(prefix + ".summary.txt");
        summary << "variant=" << kEnginePolicy << '\n';
        summary << "capture_enabled=" << (kCaptureEnabled ? 1 : 0) << '\n';
        summary << "basis_source=" << (externalBasis ? "nwchem-directory" : "sto-3g-fixture")
                << '\n';

        if (externalBasis)
        {
            summary << "basis_directory=" << basisDirectory << '\n';
        }

        summary << "n=" << n << '\n';
        summary << "molecule_atoms=" << molecule.AtomCount() << '\n';
        summary << "pair_function_count=" << pairFunctionCount << '\n';

        for (const qcx::basisset::ElementBasis& element : basis.Elements())
        {
            std::size_t elementFunctions = 0;

            for (const qcx::basisset::Shell& shell : element.shells)
            {
                const std::size_t angular = shell.isSpherical
                                                ? 2 * static_cast<std::size_t>(shell.angularMomentum) + 1
                                                : static_cast<std::size_t>(shell.angularMomentum + 1) *
                                                      static_cast<std::size_t>(shell.angularMomentum + 2) / 2;
                elementFunctions += angular * shell.coefficients.size();
            }

            summary << "element_Z_" << element.atomicNumber << "_shells=" << element.shells.size()
                    << " functions=" << elementFunctions << '\n';
        }

        summary << "shells=" << shells << '\n';
        summary << "shell_pairs=" << pairs << '\n';
        summary << "canonical_quartets=" << canonicalQuartets << '\n';
        summary << "screened_quartets=" << screenedQuartets << '\n';
        summary << "e_nuc=" << nuclearRepulsion << '\n';
        summary << "capture_calls=" << captureCalls << '\n';
        summary << "capture_unique=" << captureUnique << '\n';
        summary << "capture_ladder_mismatches=" << ladderMismatches << '\n';

        for (std::size_t order = 0; order < callsByOrder.size(); ++order)
        {
            if (callsByOrder[order] != 0 || uniqueByOrder[order] != 0)
            {
                summary << "calls_nmax_" << order << '=' << callsByOrder[order] << '\n';
                summary << "unique_nmax_" << order << '=' << uniqueByOrder[order] << '\n';
            }
        }
        summary << "seed_coulomb=" << seedEnergy.coulomb << '\n';
        summary << "seed_exchange=" << seedEnergy.exchange << '\n';
        summary << "seed_e2e=" << seedEnergy.twoElectron << '\n';
        summary << "seed_trdh=" << seedEnergy.coreTrace << '\n';
        summary << "seed_total=" << seedEnergy.total << '\n';
        summary << "default_coulomb=" << defaultEnergy.coulomb << '\n';
        summary << "default_exchange=" << defaultEnergy.exchange << '\n';
        summary << "default_e2e=" << defaultEnergy.twoElectron << '\n';
        summary << "default_trdh=" << defaultEnergy.coreTrace << '\n';
        summary << "default_total=" << defaultEnergy.total << '\n';
        summary << "tight_coulomb=" << tightEnergy.coulomb << '\n';
        summary << "tight_exchange=" << tightEnergy.exchange << '\n';
        summary << "tight_e2e=" << tightEnergy.twoElectron << '\n';
        summary << "tight_trdh=" << tightEnergy.coreTrace << '\n';
        summary << "tight_total=" << tightEnergy.total << '\n';
        summary << "default_trace_identity_residual="
                << std::abs(defaultEnergy.coreTrace + defaultEnergy.twoElectron -
                            defaultArm.electronicEnergy)
                << '\n';
        summary << "tight_trace_identity_residual="
                << std::abs(tightEnergy.coreTrace + tightEnergy.twoElectron -
                            tightArm.electronicEnergy)
                << '\n';
        WriteScfSummary(summary, "scf_default", defaultArm);
        WriteScfSummary(summary, "scf_tight", tightArm);
    }

    {
        std::ofstream record(prefix + ".boys_record.txt");
        record << qcx::integrals::BoysBuildRecord();
    }

    std::cout << "variant=" << kEnginePolicy << " n=" << n << " captured=" << captureCalls
              << " unique=" << captureUnique << '\n';
    std::cout << std::setprecision(17) << "seed_e2e=" << seedEnergy.twoElectron
              << " seed_total=" << seedEnergy.total << '\n';
    std::cout << "scf_default total=" << defaultArm.totalEnergy
              << " converged=" << (defaultArm.converged ? 1 : 0)
              << " iterations=" << defaultArm.iterations << '\n';
    std::cout << "scf_tight total=" << tightArm.totalEnergy
              << " converged=" << (tightArm.converged ? 1 : 0)
              << " iterations=" << tightArm.iterations << '\n';
    return 0;
}
