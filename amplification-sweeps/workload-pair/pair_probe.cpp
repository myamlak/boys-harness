// The Paper-A workload pair probe (scratch copy only).
//
// One build, six ERI tensors: the named recursed policy at every order, the
// matrix-product direct entry at every order, and the paper's per-order model
// under two metrics and two readings of its population. The engine reaches the
// routes through the runtime selection in integrals/src/internal/md_boys.hpp
// (QCX_BRIDGE_BOYS_PAIR_MODE), which is this scratch tree's only seam and does
// not exist in the host tree.
//
// Each arm builds the dense ERI tensor through the MD engine, captures every
// Boys argument and the ladder the engine actually consumed, runs RHF at the
// default and tight tolerances, and dumps the tensor, densities, one-electron
// matrices, the per-argument ladder columns and the fixed-density energies.
// The offline analysis compares the emitted ladders with a 60-digit reference
// and the tensors and energies with each other.
//
// Usage: qcx-pair-probe <out-prefix> [basis-directory]
//
// With no basis directory the probe runs the in-tree H2O/STO-3G fixture.

// Paper A's policy reads the coarse split-Clenshaw Chebyshev table, which the
// library's front door does not open. It must precede boys/boys.hpp in the
// translation unit that instantiates the policy.
#include "boys/boys_coefficients_d_coarse_cheb.hpp"

#include "internal/md_boys.hpp"

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

using qcx::integrals::internal::BoysPairArm;
using qcx::integrals::internal::BoysPairDirectAt;
using qcx::integrals::internal::BoysPairSetArm;
using qcx::integrals::internal::PaperABoysPolicy;
using qcx::integrals::internal::BoysCaptureRecord;
using qcx::integrals::internal::BoysCaptureTake;

constexpr std::size_t kLadderSize = static_cast<std::size_t>(boys::kMaxBoysOrder) + 1;

struct ArmSpec {
    const char* name;
    BoysPairArm arm;
    bool populationLimited;
};

constexpr std::array<ArmSpec, 6> kArms = {{
    {"paper-a", BoysPairArm::kPaperA, false},
    {"direct", BoysPairArm::kDirect, false},
    {"model-abs-band", BoysPairArm::kModelAbsolute, false},
    {"model-abs-pop", BoysPairArm::kModelAbsolute, true},
    {"model-rel-band", BoysPairArm::kModelRelative, false},
    {"model-rel-pop", BoysPairArm::kModelRelative, true},
}};

[[noreturn]] void Fail(const std::string& what) {
    std::cerr << "qcx-pair-probe: " << what << '\n';
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

/// The closed-shell two-electron energy and the core-Hamiltonian trace at one
/// fixed density, in the SCF's own index order and compensated:
///
///   E_2e = 1/2 sum D_uv D_ls (uv|ls) - 1/4 sum D_uv D_ls (ul|sv).
struct FixedDensityEnergy {
    double coulomb = 0.0;
    double exchange = 0.0;
    double twoElectron = 0.0;
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

/// The named recursed policy on one argument, every order.
void EvaluatePaperA(int nmax, double x, double* out) {
    boys::BoysAllOrders<PaperABoysPolicy>(nmax, x, out);
}

/// The matrix-product direct entry on one argument, every order. Above region
/// A the entry has no lane; the value there is the named policy's, which is
/// what every arm emits for such an argument.
void EvaluateDirect(int nmax, double x, double* out) {
    if (x < boys::kRegionAEnd)
    {
        const boys::RegionABand band =
            x < boys::kRegionA1Edge ? boys::RegionABand::kA1 : boys::RegionABand::kA2;
        boys::BoysRegionAProduct<boys::ProductMode::kFp64>(band, nmax, &x, out, 1);
    } else
    {
        boys::BoysAllOrders<PaperABoysPolicy>(nmax, x, out);
    }
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

        static_cast<void>(ReadLeU64(in)); // iteration
        const std::uint64_t dimension = ReadLeU64(in);

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

        if (tag == 2)
        {
            return matrix;
        }
    }

    Fail("the density dump holds no density record");

    return Eigen::MatrixXd{};
}

struct UniqueArgument {
    int nmax = 0;
    double x = 0.0;
    std::uint64_t xBits = 0;
    std::size_t count = 0;
    std::array<double, kLadderSize> engineLadder{};
};

std::vector<UniqueArgument>
UniqueArguments(const std::vector<BoysCaptureRecord>& records, std::size_t& ladderMismatches) {
    std::map<std::pair<int, std::uint64_t>, std::size_t> index;
    std::vector<UniqueArgument> unique;
    ladderMismatches = 0;

    for (const auto& record : records)
    {
        const auto key = std::make_pair(record.nmax, std::bit_cast<std::uint64_t>(record.x));
        const auto found = index.find(key);

        if (found == index.end())
        {
            UniqueArgument argument;
            argument.nmax = record.nmax;
            argument.x = record.x;
            argument.xBits = key.second;
            argument.count = 1;
            argument.engineLadder = record.ladder;
            index.emplace(key, unique.size());
            unique.push_back(argument);
            continue;
        }

        UniqueArgument& argument = unique[found->second];
        ++argument.count;

        for (int l = 0; l <= record.nmax; ++l)
        {
            if (std::bit_cast<std::uint64_t>(argument.engineLadder[static_cast<std::size_t>(l)]) !=
                std::bit_cast<std::uint64_t>(record.ladder[static_cast<std::size_t>(l)]))
            {
                ++ladderMismatches;
                break;
            }
        }
    }

    std::sort(unique.begin(), unique.end(), [](const UniqueArgument& a, const UniqueArgument& b) {
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

    return unique;
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
    if (argc < 2 || argc > 3)
    {
        std::cerr << "usage: qcx-pair-probe <out-prefix> [basis-directory]\n";
        return 2;
    }

    const std::string prefix = argv[1];
    const std::string basisDirectory = argc == 3 ? argv[2] : std::string{};
    const bool externalBasis = argc == 3;

    auto moleculeResult = qcx::testing::MakeH2oSto3g();
    const qcx::molecule::Molecule& molecule = Require(moleculeResult, "MakeH2oSto3g");

    qcx::Result<qcx::basisset::BasisSet> basisResult =
        externalBasis
            ? qcx::basisset::ParseNwchemDirectoryFiltered(basisDirectory,
                                                          std::array<int, 2>{1, 8})
            : qcx::testing::MakeH2oSto3gBasis();
    const qcx::basisset::BasisSet& basis =
        Require(basisResult, externalBasis ? "ParseNwchemDirectoryFiltered" : "MakeH2oSto3gBasis");

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

    for (const ArmSpec& spec : kArms)
    {
        BoysPairSetArm(spec.arm, spec.populationLimited);

        auto eriResult = qcx::integrals::BuildEriTensorGeneral(molecule, basis, eriOptions);
        const Tensor4& eri = Require(eriResult, "BuildEriTensorGeneral");

        const std::vector<BoysCaptureRecord> records = BoysCaptureTake();
        std::size_t ladderMismatches = 0;
        const std::vector<UniqueArgument> unique = UniqueArguments(records, ladderMismatches);

        std::array<std::size_t, kLadderSize> callsByOrder{};
        std::array<std::size_t, kLadderSize> uniqueByOrder{};
        std::array<std::size_t, kLadderSize> directOrders{};
        std::array<std::size_t, kLadderSize> recursedOrders{};
        std::size_t regionAArguments = 0;
        std::size_t aboveRegionAArguments = 0;

        for (const auto& record : records)
        {
            ++callsByOrder[static_cast<std::size_t>(record.nmax)];
        }

        std::array<double, kLadderSize> paperA{};
        std::array<double, kLadderSize> direct{};

        std::ofstream args = Open(prefix + "." + spec.name + ".args.txt");
        args << "# arm=" << spec.name << " nmax count x paper-a[0..nmax] direct[0..nmax]"
             << " engine[0..nmax]\n";

        for (const UniqueArgument& argument : unique)
        {
            ++uniqueByOrder[static_cast<std::size_t>(argument.nmax)];

            if (argument.x < boys::kRegionAEnd)
            {
                ++regionAArguments;
            } else
            {
                ++aboveRegionAArguments;
            }

            paperA.fill(0.0);
            direct.fill(0.0);
            EvaluatePaperA(argument.nmax, argument.x, paperA.data());
            EvaluateDirect(argument.nmax, argument.x, direct.data());

            args << argument.nmax << ' ' << argument.count << ' ' << argument.x;

            for (int l = 0; l <= argument.nmax; ++l)
            {
                args << ' ' << paperA[static_cast<std::size_t>(l)];
            }

            for (int l = 0; l <= argument.nmax; ++l)
            {
                args << ' ' << direct[static_cast<std::size_t>(l)];
            }

            for (int l = 0; l <= argument.nmax; ++l)
            {
                args << ' ' << argument.engineLadder[static_cast<std::size_t>(l)];
            }

            args << '\n';

            for (int order = 0; order <= argument.nmax; ++order)
            {
                if (BoysPairDirectAt(argument.x, order))
                {
                    ++directOrders[static_cast<std::size_t>(order)];
                } else
                {
                    ++recursedOrders[static_cast<std::size_t>(order)];
                }
            }
        }

        args.close();

        {
            std::ofstream eriOut = Open(prefix + "." + spec.name + ".eri.txt");
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
                            eriOut << i << ' ' << j << ' ' << k << ' ' << l << ' ' << value
                                   << '\n';
                        }
                    }
                }
            }

            eriOut.close();
            std::ofstream maxOut = Open(prefix + "." + spec.name + ".eri_max.txt");
            maxOut << maxAbs << '\n';
        }

        qcx::scf::RhfOptions defaultOptions;
        defaultOptions.densityDumpFile = prefix + "." + spec.name + ".default.density.bin";
        defaultOptions.traceFile = prefix + "." + spec.name + ".default.trace.txt";
        const ScfArm defaultArm = RunScf(molecule, overlap, core, eri, defaultOptions);

        qcx::scf::RhfOptions tightOptions;
        tightOptions.maxIterations = 500;
        tightOptions.energyTolerance = 1e-12;
        tightOptions.densityTolerance = 1e-10;
        tightOptions.traceFile = prefix + "." + spec.name + ".tight.trace.txt";
        const ScfArm tightArm = RunScf(molecule, overlap, core, eri, tightOptions);

        const Eigen::MatrixXd seedDensity =
            ReadSeedDensity(defaultOptions.densityDumpFile, n);

        const FixedDensityEnergy seedEnergy =
            EnergyAtDensity(eri, seedDensity, core, nuclearRepulsion, n);
        const FixedDensityEnergy defaultEnergy =
            EnergyAtDensity(eri, defaultArm.density, core, nuclearRepulsion, n);
        const FixedDensityEnergy tightEnergy =
            EnergyAtDensity(eri, tightArm.density, core, nuclearRepulsion, n);

        {
            std::ofstream matrices = Open(prefix + "." + spec.name + ".matrices.txt");
            WriteMatrix(matrices, "overlap", overlap);
            WriteMatrix(matrices, "core", core);
            WriteMatrix(matrices, "seed_density", seedDensity);
            WriteMatrix(matrices, "default_density", defaultArm.density);
            WriteMatrix(matrices, "tight_density", tightArm.density);
        }

        {
            std::ofstream summary = Open(prefix + "." + spec.name + ".summary.txt");
            summary << "variant=" << spec.name << '\n';
            summary << "arm=" << static_cast<int>(spec.arm) << '\n';
            summary << "population_limited=" << (spec.populationLimited ? 1 : 0) << '\n';
            summary << "n=" << n << '\n';
            summary << "e_nuc=" << nuclearRepulsion << '\n';
            summary << "capture_calls=" << records.size() << '\n';
            summary << "capture_unique=" << unique.size() << '\n';
            summary << "capture_ladder_mismatches=" << ladderMismatches << '\n';
            summary << "region_a_arguments=" << regionAArguments << '\n';
            summary << "above_region_a_arguments=" << aboveRegionAArguments << '\n';

            for (std::size_t order = 0; order < callsByOrder.size(); ++order)
            {
                if (callsByOrder[order] != 0 || uniqueByOrder[order] != 0)
                {
                    summary << "calls_nmax_" << order << '=' << callsByOrder[order] << '\n';
                    summary << "unique_nmax_" << order << '=' << uniqueByOrder[order] << '\n';
                    summary << "direct_orders_nmax_" << order << '=' << directOrders[order] << '\n';
                    summary << "recursed_orders_nmax_" << order << '=' << recursedOrders[order]
                            << '\n';
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

        std::cout << std::setprecision(17) << spec.name << " n=" << n
                  << " calls=" << records.size() << " unique=" << unique.size()
                  << " seed_e2e=" << seedEnergy.twoElectron
                  << " scf_default=" << defaultArm.totalEnergy
                  << " scf_tight=" << tightArm.totalEnergy << '\n';
    }

    {
        std::ofstream record = Open(prefix + ".boys_record.txt");
        record << qcx::integrals::BoysBuildRecord();
    }

    return 0;
}
