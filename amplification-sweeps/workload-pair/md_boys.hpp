#pragma once

// The Boys seed adapter of the MD engine: the fp64 pipeline consumes
// BoysAllOrders directly, the certified fp32 pipeline consumes
// BoysAllOrdersF32 - whose region-A seeds are computed in double internally,
// so the amplification-sensitive seeds are safe. The fp32-vs-fp64
// choice is a per-call precision policy, never a per-quartet branch inside
// the kernel: the caller's gate decides which pipeline each quartet goes
// through.

// Paper-A chemistry bridge (scratch copy only): the measured recursed policy
// reads the coarse split-Clenshaw Chebyshev table, and the library's front
// door opens only the default class row's coefficient groups. The table's
// header must therefore precede boys/boys.hpp in every unit that instantiates
// this adapter, or the fit's static_assert names the missing include.
#if defined(QCX_BRIDGE_BOYS_POLICY_PAPER_A) || defined(QCX_BRIDGE_BOYS_PAIR_MODE)
#include "boys/boys_coefficients_d_coarse_cheb.hpp"
#endif

#include "qcx/integrals/boys.hpp"

#include <cstddef>
#include <type_traits>

// Paper-A workload pair (scratch copy only): a runtime, per-order selection
// between the paper's named recursed policy and its matrix-product direct
// entry. The engine has no such switch; this is the smallest seam that lets a
// single ERI tensor mix the two routes order by order. The route is set before
// a tensor build and read at each ladder call, so the arithmetic of every
// other engine step is the one build's own.
#if defined(QCX_BRIDGE_BOYS_PAIR_MODE)
#include "boys/boys_transform.hpp"

#include <array>
#include <atomic>
#endif

#if defined(QCX_BRIDGE_BOYS_CAPTURE)
#include <array>
#include <mutex>
#include <vector>
#endif

// The direct-timing instrument of the Boys-share probe
// (benchmarks/boys_share_probe.cpp). It is defined only where
// QCX_BOYS_SHARE_INSTRUMENT is, which the build sets on qcx-integrals and on
// that one benchmark and on nothing else, so every other build compiles this
// header to the text it always had.
//
// Two instruments run at once, because they fail differently. The stamp
// counter reads the calls' own cost in place, and the elision flag beside it
// takes the library call out altogether so the same spans can be differenced
// with no per-call model at all. The first is what the second checks.
#if defined(QCX_BOYS_SHARE_INSTRUMENT)

#include <cstdint>

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

namespace qcx::integrals::internal {

/// The time-stamp counter, read at the instrument's two ends.
inline std::uint64_t BoysShareStamp() noexcept {
    return static_cast<std::uint64_t>(__rdtsc());
}

/// Cycles spent inside the instrumented calls, and how many of them, on the
/// calling thread. Valid only while \c gBoysShareSkip is false: the elided
/// arm times its own zero-fill instead.
inline thread_local std::uint64_t gBoysShareCycles = 0;
inline thread_local std::uint64_t gBoysShareCalls = 0;

/// True while the probe runs the Boys-elided arm. The call is then skipped and
/// the ladder left zeroed, so the arm executes the same instruction stream
/// with the library's arithmetic taken out. It is one predictable branch per
/// call and it is present in both arms, so it cancels in the difference.
inline bool gBoysShareSkip = false;

/// The ladder an elided call leaves: zeros, which no later step can turn into
/// a denormal, so the rest of the kernel runs at its own speed.
template <typename T> inline void BoysShareZero(int nmax, T* out) noexcept {
    for (int l = 0; l <= nmax; ++l)
    {
        out[l] = T{};
    }
}

} // namespace qcx::integrals::internal

#endif

namespace qcx::integrals::internal {

#if defined(QCX_BRIDGE_BOYS_POLICY_PAPER_A) || defined(QCX_BRIDGE_BOYS_PAIR_MODE)

/// Paper A's measured recursed policy, named exactly as its option row states
/// it: the Chebyshev route summed by the split Clenshaw recurrence, the
/// coarsest partition, the refined reciprocal, the fast region-B exponential.
/// The budget cell is inert on the double lane; the packing axis names the
/// arguments axis the paper's measurement used, which this all-orders entry
/// does not vectorise over (one argument, every order) - it selects the same
/// scalar body either way.
using PaperABoysPolicy =
    boys::StatedEvalPolicy<boys::FitRoute::kChebyshev,
                           boys::EvalScheme::kSplitClenshaw,
                           boys::BoysBudget::kFloat,
                           boys::PackAxis::kArguments,
                           boys::FitGranularity::kCoarsest,
                           boys::DivisionForm::kRefinedReciprocal,
                           boys::RegionBExp::kFast>;

#endif

#if defined(QCX_BRIDGE_BOYS_PAIR_MODE)

/// The four arms the workload probe runs, and the two readings of the model's
/// population it builds tensors for.
///
/// The model's state is the paper's own: at margin 1 the second experiment's
/// record names, per band, metric and order, the evaluation with the smaller
/// truncation term, and abstains where the two modelled errors differ by no
/// more than the recursed bound's own accumulated-arithmetic term. The tables
/// below are that record, transcribed; a state of 0 is an abstention, and the
/// caller's named fallback is the direct entry, so a 0 and a -1 both take the
/// direct route.
///
/// Two readings of the model's population are built because the paper's
/// derivation restricts it to x >= the first dispatch threshold (below it the
/// entry seeds downward from order 32, so A_n is not the transfer), while the
/// paper's own end-to-end table applies a band-order commitment to the whole
/// band. Both tensors are formed; neither is presented as the other.
enum class BoysPairArm : int {
    kPaperA = 0,
    kDirect = 1,
    kModelAbsolute = 2,
    kModelRelative = 3,
};

/// +1: the model names the recursed evaluation. -1: it names the direct entry.
/// 0: it abstains (fallback: direct).
inline constexpr std::array<std::array<signed char, boys::kMaxBoysOrder + 1>, 2> kPairAbsolute = {{
    {1, 0, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0},
    {-1, -1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1, -1, -1},
}};

inline constexpr std::array<std::array<signed char, boys::kMaxBoysOrder + 1>, 2> kPairRelative = {{
    {1, 1, 0, 0, -1, 0, 0, -1, -1, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0},
    {-1, -1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1},
}};

/// The first dispatch threshold, which is also the paper's model population
/// starts. Pinned against the kernel's own table: a revision that moves it
/// fails this probe's build rather than silently moving the population.
inline constexpr double kPairPopulationStart = 1.08552523453493333;
static_assert(kPairPopulationStart == boys::detail::kTierThresholds[0]);

/// The arm one tensor build runs, and whether the model's population is
/// restricted to the paper's. Both are read at every ladder call; they are set
/// before a build and never during one.
inline std::atomic<int> gPairArm{static_cast<int>(BoysPairArm::kPaperA)};
inline std::atomic<int> gPairPopulationLimited{1};

inline void BoysPairSetArm(BoysPairArm arm, bool populationLimited) noexcept {
    gPairArm.store(static_cast<int>(arm), std::memory_order_relaxed);
    gPairPopulationLimited.store(populationLimited ? 1 : 0, std::memory_order_relaxed);
}

/// True when this (argument, order) cell takes the direct entry.
inline bool BoysPairDirectAt(double x, int order) noexcept {
    // Outside the two fitted bands the direct entry has no lane:
    // BoysRegionAProduct takes region-A arguments only. Such an argument is
    // forced to the named policy in every arm and cannot carry a difference.
    if (!(x < boys::kRegionAEnd))
    {
        return false;
    }

    const int arm = gPairArm.load(std::memory_order_relaxed);

    if (arm == static_cast<int>(BoysPairArm::kPaperA))
    {
        return false;
    }

    if (arm == static_cast<int>(BoysPairArm::kDirect))
    {
        return true;
    }

    // The paper's model has no state below its population start; the named
    // fallback there is the direct entry.
    if (gPairPopulationLimited.load(std::memory_order_relaxed) != 0 &&
        x < kPairPopulationStart)
    {
        return true;
    }

    const int band = x < boys::kRegionA1Edge ? 0 : 1;
    const bool absolute = arm == static_cast<int>(BoysPairArm::kModelAbsolute);
    const int state =
        static_cast<int>((absolute ? kPairAbsolute : kPairRelative)[static_cast<std::size_t>(band)]
                                       [static_cast<std::size_t>(order)]);
    return state <= 0;
}

/// The mixed ladder: the named policy's ladder and the direct entry's, copied
/// per order. Both are formed for every argument so that the two route columns
/// exist beside the emitted one; the selection alone decides what the engine
/// sees.
inline void BoysPairLadder(int nmax, double x, double* out) noexcept {
    std::array<double, boys::kMaxBoysOrder + 1> recursed{};
    std::array<double, boys::kMaxBoysOrder + 1> direct{};

    boys::BoysAllOrders<PaperABoysPolicy>(nmax, x, recursed.data());

    if (x < boys::kRegionAEnd)
    {
        const boys::RegionABand band =
            x < boys::kRegionA1Edge ? boys::RegionABand::kA1 : boys::RegionABand::kA2;
        boys::BoysRegionAProduct<boys::ProductMode::kFp64>(band, nmax, &x, direct.data(), 1);
    }

    for (int order = 0; order <= nmax; ++order)
    {
        out[order] = BoysPairDirectAt(x, order)
                         ? direct[static_cast<std::size_t>(order)]
                         : recursed[static_cast<std::size_t>(order)];
    }
}

#endif

#if defined(QCX_BRIDGE_BOYS_CAPTURE)

/// One argument the engine formed and the ladder it emitted for it: the
/// ladder's top order, the reduced argument x = p q / (p + q) * |P - Q|^2,
/// and F_0(x)..F_nmax(x) as the engine's own call returned them.
struct BoysCaptureRecord {
    int nmax;
    double x;
    std::array<double, kMaxBoysOrder + 1> ladder;
};

/// The capture sink every translation unit of the binary shares. The static
/// lives in an inline function, so the program holds one instance however
/// many class translation units call the adapter.
struct BoysCaptureSink {
    std::mutex mutex;
    std::vector<BoysCaptureRecord> records;
};

inline BoysCaptureSink& BoysCapture() {
    static BoysCaptureSink sink;
    return sink;
}

/// Records one argument and the ladder it returned. It never changes a
/// computed value and is compiled only where the bridge asks for it.
inline void BoysCaptureAppend(int nmax, double x, const double* ladder) {
    BoysCaptureRecord record;
    record.nmax = nmax;
    record.x = x;

    for (int l = 0; l <= nmax; ++l)
    {
        record.ladder[static_cast<std::size_t>(l)] = ladder[l];
    }

    BoysCaptureSink& sink = BoysCapture();
    const std::lock_guard<std::mutex> lock(sink.mutex);
    sink.records.push_back(record);
}

/// Takes every record made so far and leaves the sink empty.
inline std::vector<BoysCaptureRecord> BoysCaptureTake() {
    BoysCaptureSink& sink = BoysCapture();
    const std::lock_guard<std::mutex> lock(sink.mutex);
    std::vector<BoysCaptureRecord> taken = std::move(sink.records);
    sink.records.clear();
    return taken;
}

#endif

/// F_0(x)..F_nmax(x) into out, in the working scalar type T.
///
/// The shape is per-argument: one argument, the whole ladder up to \p nmax.
/// The engine forms its argument inside the primitive-quadruple loop and
/// consumes the ladder on the next line, and the library's many-argument
/// entries buy their grouping only from a gather the caller has to perform
/// first - see the seed-entry note in qcx/integrals/boys.hpp for the
/// measurement behind the choice and for what a batched shape was measured to
/// be worth on this engine's own argument stream.
/// \tparam T double (the library's published bound is the larger of
/// 1e-13 * F_n(x) and 1e-13 * F_0(x) per value) or float (its float lane's is
/// the larger of 2^-22 * F_n(x) and 1.5e-7 per value, plus 1e-7 under the
/// plain reciprocal).
/// \tparam kFp64Seeds false, the default, leaves the float arm on the fp32
/// entry's own seeds; true forms the ladder in double and rounds each value
/// once, so that route's seed error is the fp64 bound plus one float rounding
/// per value. Which of the two the engine's argument stream favours is a
/// measurement question, not an assumption.
template <typename T, bool kFp64Seeds = false>
inline void MdBoysBatch(int nmax, double x, T* out) noexcept {
    static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>);

#if defined(QCX_BOYS_SHARE_INSTRUMENT)
    const std::uint64_t boysShareStart = BoysShareStamp();

    if (gBoysShareSkip)
    {
        BoysShareZero(nmax, out);
        gBoysShareCycles += BoysShareStamp() - boysShareStart;
        ++gBoysShareCalls;

        return;
    }
#endif

    if constexpr (std::is_same_v<T, double>)
    {
#if defined(QCX_BRIDGE_BOYS_PAIR_MODE)
        BoysPairLadder(nmax, x, out);
#elif defined(QCX_BRIDGE_BOYS_POLICY_PAPER_A)
        boys::BoysAllOrders<PaperABoysPolicy>(nmax, x, out);
#else
        BoysAllOrders(nmax, x, out);
#endif

#if defined(QCX_BRIDGE_BOYS_CAPTURE)
        BoysCaptureAppend(nmax, x, out);
#endif
    } else if constexpr (kFp64Seeds)
    {
        // The fp64 entry writes double, so its ladder lands in a stack buffer
        // and is rounded once per value on the way into the caller's array -
        // the whole of what this arm adds over the one below.
        double ladder[kMaxBoysOrder + 1];
        BoysAllOrders(nmax, x, ladder);

        for (int l = 0; l <= nmax; ++l)
        {
            out[l] = static_cast<float>(ladder[l]);
        }
    } else
    {
        // The ladder is written straight into the caller's array: the fp32
        // entry takes and returns float, so an intermediate buffer would be a
        // float-to-float copy and nothing else. The argument is the one
        // conversion, and it is the call's own (the region-A seeds are
        // computed in double inside the entry, so the
        // amplification-sensitive seeds are safe).
        BoysAllOrdersF32(nmax, static_cast<float>(x), out);
    }

#if defined(QCX_BOYS_SHARE_INSTRUMENT)
    gBoysShareCycles += BoysShareStamp() - boysShareStart;
    ++gBoysShareCalls;
#endif
}

} // namespace qcx::integrals::internal
