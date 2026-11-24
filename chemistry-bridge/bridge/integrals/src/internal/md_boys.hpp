#pragma once

// The Boys seed adapter of the MD engine: the fp64 pipeline consumes
// BoysAllOrders directly, the certified fp32 pipeline consumes
// BoysAllOrdersF32 - whose region-A seeds are computed in double internally,
// so the amplification-sensitive seeds are safe. The fp32-vs-fp64
// choice is a per-call precision policy, never a per-quartet branch inside
// the kernel: the caller's gate decides which pipeline each quartet goes
// through.

// The route-comparison chemistry bridge: the measured recursed policy
// reads the coarse split-Clenshaw Chebyshev table, and the library's front
// door opens only the default class row's coefficient groups. The table's
// header must therefore precede boys/boys.hpp in every unit that instantiates
// this adapter, or the fit's static_assert names the missing include.
#if defined(QCX_BRIDGE_BOYS_POLICY_PAPER_A)
#include "boys/boys_coefficients_d_coarse_cheb.hpp"
#endif

#include "qcx/integrals/boys.hpp"

#include <cstddef>
#include <type_traits>

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

#if defined(QCX_BRIDGE_BOYS_POLICY_PAPER_A)

/// The route-comparison study's measured recursed policy, named exactly as
/// its option row states it: the Chebyshev route summed by the split Clenshaw
/// recurrence, the coarsest partition, the refined reciprocal, the fast
/// region-B exponential.
/// The budget cell is inert on the double lane; the packing axis names the
/// arguments axis the route-comparison study's measurement used, which this all-orders entry
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
#if defined(QCX_BRIDGE_BOYS_POLICY_PAPER_A)
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
