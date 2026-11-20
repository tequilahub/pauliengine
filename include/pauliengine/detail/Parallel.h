#pragma once
/// Internal: when and how PauliEngine starts OpenMP threads.

#include <cstddef>

#ifdef PAULIENGINE_HAS_OPENMP
#include <omp.h>
#endif

// Work below which no thread team is started: records for sorting and merging,
// pairs times words for forming products. Starting a team costs a few microseconds
// with libgomp (g++), much less with the MSVC runtime.
#ifndef PAULIENGINE_OMP_THRESHOLD
#if defined(_MSC_VER) && !defined(__clang__)
#define PAULIENGINE_OMP_THRESHOLD 4096
#else
#define PAULIENGINE_OMP_THRESHOLD 32768
#endif
#endif

namespace pauliengine_detail {

/// True if `work` is large enough to start a thread team.
inline bool use_threads(size_t work) {
#ifdef PAULIENGINE_HAS_OPENMP
        return work >= PAULIENGINE_OMP_THRESHOLD && omp_get_max_threads() > 1;
#else
        (void)work;
        return false;
#endif
}

/// Calls f(i) for i in [0, n), on all threads if `parallel`. A serial call opens no
/// OpenMP region at all; even an inactive one costs microseconds with libgomp.
template<typename F>
inline void for_each_index(size_t n, bool parallel, F&& f) {
        const std::ptrdiff_t count = static_cast<std::ptrdiff_t>(n);  // signed for OpenMP 2.0
#ifdef PAULIENGINE_HAS_OPENMP
        if (parallel) {
                #pragma omp parallel for schedule(static)
                for (std::ptrdiff_t i = 0; i < count; ++i) f(static_cast<size_t>(i));
                return;
        }
#else
        (void)parallel;
#endif
        for (std::ptrdiff_t i = 0; i < count; ++i) f(static_cast<size_t>(i));
}

}  // namespace pauliengine_detail
