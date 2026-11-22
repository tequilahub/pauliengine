#pragma once
#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

#include "pauliengine/detail/FlatTerms.h"
#include "pauliengine/detail/Parallel.h"

namespace pauliengine_detail {

/// The pairs (i, j) of two sums whose strings share a qubit, row by row with
/// ascending j. Only these pairs can anticommute, so for strings with small
/// supports (local Hamiltonians) a commutator needs to look at them alone instead
/// of at all |a| |b| pairs. Found through an index of b by qubit: row i collects the
/// terms of b on each qubit of its own support, without repeats.
class OverlappingPairs {
    public:
        std::vector<size_t> start;    ///< row i owns cols[start[i] .. start[i + 1])
        std::vector<uint32_t> cols;   ///< indices into b

        size_t size() const { return cols.size(); }

        /// The overlapping pairs, or nothing when testing all pairs costs less:
        /// collecting a pair (index lookup, repeat check, sort) costs about
        /// kCollectCost tests, nearly independent of the string width.
        static std::optional<OverlappingPairs> find(const FlatTerms& a, const FlatTerms& b) {
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                if (n2 >= std::numeric_limits<uint32_t>::max()) {
                        return std::nullopt;
                }
                // qubits beyond the shorter width carry nothing on one side
                const size_t w = std::min(a.width(), b.width());
                const size_t nq = 64 * w;
                std::vector<size_t> on_a(nq, 0);
                std::vector<size_t> at(nq + 1, 0);  // terms of b per qubit, then offsets
                for (size_t i = 0; i < n1; ++i) {
                        for_each_site(a, i, w, [&](size_t q) { ++on_a[q]; });
                }
                for (size_t j = 0; j < n2; ++j) {
                        for_each_site(b, j, w, [&](size_t q) { ++at[q + 1]; });
                }
                // pairs sharing a qubit, counted once per shared qubit
                double collected = 0.0;
                for (size_t q = 0; q < nq; ++q) {
                        collected += static_cast<double>(on_a[q]) * static_cast<double>(at[q + 1]);
                }
                if (kCollectCost * collected >= static_cast<double>(n1) * static_cast<double>(n2)) {
                        return std::nullopt;
                }
                for (size_t q = 0; q < nq; ++q) {
                        at[q + 1] += at[q];
                }
                std::vector<uint32_t> by_qubit(at[nq]);
                {
                        std::vector<size_t> fill(at.begin(), at.end() - 1);
                        for (size_t j = 0; j < n2; ++j) {
                                for_each_site(b, j, w, [&](size_t q) { by_qubit[fill[q]++] = static_cast<uint32_t>(j); });
                        }
                }
                // Row i: the terms of b on its qubits, each once (marked with a per-row
                // stamp); counted first, then written into the reserved range and sorted.
                const bool parallel = use_threads(static_cast<size_t>(collected) + n1);
                const auto collect = [&](size_t i, uint32_t* out) {
                        thread_local std::vector<uint64_t> seen;
                        thread_local uint64_t stamp = 0;
                        if (seen.size() < n2) {
                                seen.assign(n2, 0);
                                stamp = 0;
                        }
                        ++stamp;
                        size_t found = 0;
                        for_each_site(a, i, w, [&](size_t q) {
                                for (size_t k = at[q]; k < at[q + 1]; ++k) {
                                        const uint32_t j = by_qubit[k];
                                        if (seen[j] != stamp) {
                                                seen[j] = stamp;
                                                if (out) out[found] = j;
                                                ++found;
                                        }
                                }
                        });
                        return found;
                };
                OverlappingPairs pairs;
                pairs.start.assign(n1 + 1, 0);
                for_each_index(n1, parallel, [&](size_t i) { pairs.start[i + 1] = collect(i, nullptr); });
                for (size_t i = 0; i < n1; ++i) {
                        pairs.start[i + 1] += pairs.start[i];
                }
                pairs.cols.resize(pairs.start[n1]);
                for_each_index(n1, parallel, [&](size_t i) {
                        uint32_t* row = pairs.cols.data() + pairs.start[i];
                        std::sort(row, row + collect(i, row));
                });
                return pairs;
        }

    private:
        /// Cost of collecting one pair (counted once per shared qubit), in pair tests;
        /// measured break-even between 8 and 16.
        static constexpr double kCollectCost = 16.0;

        template <typename F>
        static void for_each_site(const FlatTerms& t, size_t i, size_t w, F&& f) {
                const uint64_t* x = t.x(i);
                const uint64_t* y = t.y(i);
                for (size_t k = 0; k < w; ++k) {
                        uint64_t m = x[k] | y[k];
                        while (m) {
                                f(64 * k + static_cast<size_t>(std::countr_zero(m)));
                                m &= m - 1;
                        }
                }
        }
};

}  // namespace pauliengine_detail
