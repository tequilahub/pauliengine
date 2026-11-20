#pragma once

#include <algorithm>
#include <bit>
#include <complex>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/detail/PrefixPartitionedProduct.h"
#include "pauliengine/UnweightedPauliString.h"

/// An ordered set of Pauli strings without coefficients, each with its own phase i^k,
/// to be used only as the basis of a Lie algebra. It is not a linear combination: every operator
/// occurs once and the first phase seen for an operator wins, so A + A = A.
///
/// An element keeps its index for good, and an open-addressing table finds an operator
/// without touching most other elements.
class unweighted_QubitHamiltonian {
        public:
                /// The empty set.
                unweighted_QubitHamiltonian() = default;

                /// From a list of strings; later duplicates of an operator are ignored.
                explicit unweighted_QubitHamiltonian(const std::vector<unweighted_pauliString>& terms) {
                        for (const auto& term : terms) {
                                insert(term);
                        }
                }

                /// The operators of a Hamiltonian with phase 0; its coefficients are dropped.
                template<typename Coeff>
                explicit unweighted_QubitHamiltonian(const QubitHamiltonian<Coeff>& hamiltonian) {
                        for (const auto& term : hamiltonian.terms()) {
                                insert(unweighted_pauliString(term));
                        }
                }

                /// Number of elements.
                size_t size() const noexcept { return data_.size(); }

                /// Adds the string unless its operator is present; returns whether it was added.
                bool insert(const unweighted_pauliString& term) {
                        return insert_words(term.x.data(), term.y.data(), term.x.size(), term.phase());
                }

                /// True if the operator of the string is present, whatever the phases are.
                bool contains(const unweighted_pauliString& term) const {
                        return find(term).has_value();
                }

                /// Index of the operator of the string, if it is present.
                std::optional<size_t> find(const unweighted_pauliString& term) const {
                        if (slots_.empty()) {
                                return std::nullopt;
                        }
                        const size_t n = term.x.size();
                        const uint64_t h = hash_words(term.x.data(), term.y.data(), n);
                        const uint32_t tag = static_cast<uint32_t>(h >> 32);
                        const size_t mask = slots_.size() - 1;
                        for (size_t pos = h & mask; slots_[pos].index != 0; pos = (pos + 1) & mask) {
                                const size_t i = slots_[pos].index - 1;
                                if (slots_[pos].tag == tag && element_equals(i, term.x.data(), term.y.data(), n)) {
                                        return i;
                                }
                        }
                        return std::nullopt;
                }

                /// All elements in insertion order.
                const std::vector<unweighted_pauliString>& terms() const noexcept { return data_; }

                /// Union of two sets: the elements of this one, then the new ones of other.
                unweighted_QubitHamiltonian operator+(const unweighted_QubitHamiltonian& other) const {
                        unweighted_QubitHamiltonian result(*this);
                        result += other;
                        return result;
                }

                unweighted_QubitHamiltonian& operator+=(const unweighted_QubitHamiltonian& other) {
                        if (this == &other) {
                                return *this;
                        }
                        for (const auto& term : other.data_) {
                                insert(term);
                        }
                        return *this;
                }

                /// The set of all products a * b in the order of a, then b.
                unweighted_QubitHamiltonian operator*(const unweighted_QubitHamiltonian& other) const {
                        return pairwise(other, false);
                }

                /// The set of commutators without their factor 2: a * b for every
                /// anti-commuting pair, in the order of a, then b. See
                /// unweighted_pauliString::commutator().
                unweighted_QubitHamiltonian commutator(const unweighted_QubitHamiltonian& other) const {
                        return pairwise(other, true);
                }

                /// Equal sets with equal phases, whatever the order of the elements is.
                bool operator==(const unweighted_QubitHamiltonian& other) const {
                        if (size() != other.size()) {
                                return false;
                        }
                        for (const auto& term : data_) {
                                const std::optional<size_t> j = other.find(term);
                                if (!j || other.data_[*j].phase() != term.phase()) {
                                        return false;
                                }
                        }
                        return true;
                }

                /// The ordinary Hamiltonian with coefficient coeff * i^k on every element.
                template<typename Coeff>
                QubitHamiltonian<Coeff> to_qubit_hamiltonian(Coeff coeff) const {
                        std::vector<PauliString<Coeff>> result;
                        result.reserve(size());
                        for (const auto& term : data_) {
                                result.push_back(term.to_pauli_string(coeff));
                        }
                        return QubitHamiltonian<Coeff>(result);
                }

        private:
                /// term index + 1 (0 marks a free slot) and the upper hash bits, which
                /// reject most mismatches before the bits have to be compared.
                struct Slot {
                        uint32_t index{0};
                        uint32_t tag{0};
                };

                /// The table refers to these elements, so they change through insert() only.
                std::vector<unweighted_pauliString> data_;
                std::vector<Slot> slots_;

                static uint64_t hash_words(const uint64_t* x, const uint64_t* y, size_t n) noexcept {
                        uint64_t h = 0x9E3779B97F4A7C15ULL;
                        const auto mix = [&h](uint64_t value) {
                                h = (h ^ value) * 0xFF51AFD7ED558CCDULL;
                                h ^= h >> 32;
                        };
                        for (size_t i = 0; i < n; ++i) {
                                mix(x[i]);
                                mix(y[i]);
                        }
                        return h;
                }

                /// Compares with n trimmed words; the elements are trimmed as well.
                bool element_equals(size_t i, const uint64_t* x, const uint64_t* y, size_t n) const noexcept {
                        const unweighted_pauliString& element = data_[i];
                        if (element.x.size() != n) {
                                return false;
                        }
                        for (size_t k = 0; k < n; ++k) {
                                if (element.x[k] != x[k] || element.y[k] != y[k]) {
                                        return false;
                                }
                        }
                        return true;
                }

                void place(uint64_t h, uint32_t index) noexcept {
                        const size_t mask = slots_.size() - 1;
                        size_t pos = h & mask;
                        while (slots_[pos].index != 0) {
                                pos = (pos + 1) & mask;
                        }
                        slots_[pos] = Slot{index + 1, static_cast<uint32_t>(h >> 32)};
                }

                void grow_table() {
                        slots_.assign(std::max<size_t>(16, slots_.size() * 2), Slot{});
                        for (size_t i = 0; i < size(); ++i) {
                                const unweighted_pauliString& element = data_[i];
                                place(hash_words(element.x.data(), element.y.data(), element.x.size()), static_cast<uint32_t>(i));
                        }
                }

                /// x and y must not point into this object.
                bool insert_words(const uint64_t* x, const uint64_t* y, size_t n, int phase) {
                        while (n > 0 && x[n - 1] == 0 && y[n - 1] == 0) {
                                --n;
                        }
                        if ((size() + 1) * 2 > slots_.size()) {
                                grow_table();
                        }
                        const uint64_t h = hash_words(x, y, n);
                        const uint32_t tag = static_cast<uint32_t>(h >> 32);
                        const size_t mask = slots_.size() - 1;
                        for (size_t pos = h & mask; slots_[pos].index != 0; pos = (pos + 1) & mask) {
                                if (slots_[pos].tag == tag && element_equals(slots_[pos].index - 1, x, y, n)) {
                                        return false;
                                }
                        }
                        if (size() >= std::numeric_limits<uint32_t>::max() - 1) {
                                throw std::length_error("unweighted_QubitHamiltonian: too many elements");
                        }
                        place(h, static_cast<uint32_t>(size()));
                        data_.emplace_back(x, y, n, phase);
                        return true;
                }

                static size_t longest(const std::vector<unweighted_pauliString>& terms) noexcept {
                        size_t n = 0;
                        for (const auto& term : terms) {
                                n = std::max(n, term.x.size());
                        }
                        return n;
                }

                /// From this many pairs on, and with several threads, products go through
                /// the prefix-partitioned product of QubitHamiltonian
                /// (pauliengine_detail::PrefixPartitionedProduct).
                /// With one thread the loop below is as fast: building the set, which
                /// both need, takes most of the time.
                static constexpr size_t kBlockPairs = 4096;

                /// The terms as flat rows of n_words words; the coefficient of row i is
                /// (i, phase exponent), as PrefixPartitionedProduct::first_pairs() expects.
                static pauliengine_detail::FlatTerms flat_rows(const std::vector<unweighted_pauliString>& terms, size_t n_words) {
                        pauliengine_detail::FlatTerms rows(terms.size(), n_words);
                        for (size_t i = 0; i < terms.size(); ++i) {
                                const unweighted_pauliString& term = terms[i];
                                uint64_t* x = rows.x(i);
                                uint64_t* y = rows.y(i);
                                for (size_t k = 0; k < n_words; ++k) {
                                        x[k] = k < term.x.size() ? term.x[k] : 0;
                                        y[k] = k < term.y.size() ? term.y[k] : 0;
                                }
                                rows.set_coeff(i, {static_cast<double>(i), static_cast<double>(term.phase())});
                        }
                        return rows;
                }

                /// The same set as the loop over all pairs gives, in the same order and
                /// with the same phases: the prefix-partitioned product finds, for every operator, its
                /// pair of lowest rank i * other.size() + j, which is the pair the loop
                /// would insert first. The operators are then put in the order of these
                /// ranks, by row i first (counting) and by j within a row.
                std::optional<unweighted_QubitHamiltonian> pairwise_by_blocks(const unweighted_QubitHamiltonian& other,
                                                                             bool anti_commuting_only) const {
                        using pauliengine_detail::FlatTerms;
                        using pauliengine_detail::PrefixPartitionedProduct;
                        const size_t n1 = size();
                        const size_t n2 = other.size();
                        // The kernel forms a commutator only where a fair share of the pairs
                        // anti-commute, judged as for QubitHamiltonian by a grid of pairs.
                        if (anti_commuting_only) {
                                const size_t side = n1 * n2 >= 65536 ? 32 : 16;
                                const size_t rows_a = std::min(n1, side);
                                const size_t rows_b = std::min(n2, side);
                                size_t found = 0;
                                for (size_t s = 0; s < rows_a; ++s) {
                                        const unweighted_pauliString& left = data_[s * n1 / rows_a];
                                        for (size_t t = 0; t < rows_b; ++t) {
                                                found += !left.commutes(other.data_[t * n2 / rows_b]);
                                        }
                                }
                                if (8 * found < rows_a * rows_b) {
                                        return std::nullopt;
                                }
                        }
                        const size_t n_words = std::max(longest(data_), longest(other.data_));
                        const FlatTerms a = flat_rows(data_, n_words);
                        const FlatTerms b = flat_rows(other.data_, n_words);
                        const FlatTerms firsts = PrefixPartitionedProduct::first_pairs(a, b, anti_commuting_only);
                        const size_t m = firsts.size();
                        if (m >= std::numeric_limits<uint32_t>::max() - 1) {
                                throw std::length_error("unweighted_QubitHamiltonian: too many elements");
                        }
                        std::vector<uint64_t> rank(m);
                        std::vector<uint32_t> row_start(n1 + 1, 0);
                        for (size_t r = 0; r < m; ++r) {
                                rank[r] = static_cast<uint64_t>(firsts.coeff(r).real());
                                ++row_start[rank[r] / n2 + 1];
                        }
                        for (size_t i = 0; i < n1; ++i) {
                                row_start[i + 1] += row_start[i];
                        }
                        std::vector<uint32_t> order(m);
                        {
                                std::vector<uint32_t> next(row_start.begin(), row_start.end() - 1);
                                for (size_t r = 0; r < m; ++r) {
                                        order[next[rank[r] / n2]++] = static_cast<uint32_t>(r);
                                }
                        }
                        for (size_t i = 0; i < n1; ++i) {
                                std::sort(order.begin() + row_start[i], order.begin() + row_start[i + 1],
                                          [&rank](uint32_t p, uint32_t q) { return rank[p] < rank[q]; });
                        }

                        // The strings and their hashes, on all threads.
                        unweighted_QubitHamiltonian result;
                        result.data_.resize(m);
                        std::vector<uint64_t> hashes(m);
                        pauliengine_detail::for_each_index(m, pauliengine_detail::use_threads(m * n_words), [&](size_t p) {
                                const uint32_t r = order[p];
                                unweighted_pauliString& added = result.data_[p];
                                added = unweighted_pauliString(firsts.x(r), firsts.y(r), n_words, static_cast<int>(firsts.coeff(r).imag()));
                                hashes[p] = hash_words(added.x.data(), added.y.data(), added.x.size());
                        });
                        // The table is filled in the order of the home slots, a stretch of
                        // 256 slots at a time, instead of all over it.
                        const size_t n_slots = std::bit_ceil(std::max<size_t>(16, 2 * (m + 1)));
                        result.slots_.assign(n_slots, Slot{});
                        const unsigned shift = n_slots > 256 ? 8 : 0;
                        std::vector<uint32_t> stretch_start((n_slots >> shift) + 1, 0);
                        for (size_t p = 0; p < m; ++p) {
                                ++stretch_start[((hashes[p] & (n_slots - 1)) >> shift) + 1];
                        }
                        for (size_t g = 1; g < stretch_start.size(); ++g) {
                                stretch_start[g] += stretch_start[g - 1];
                        }
                        std::vector<uint32_t> by_slot(m);
                        for (size_t p = 0; p < m; ++p) {
                                by_slot[stretch_start[(hashes[p] & (n_slots - 1)) >> shift]++] = static_cast<uint32_t>(p);
                        }
                        for (const uint32_t p : by_slot) {
                                result.place(hashes[p], p);
                        }
                        return result;
                }

                unweighted_QubitHamiltonian pairwise(const unweighted_QubitHamiltonian& other, bool anti_commuting_only) const {
                        const size_t pairs = size() * other.size();
                        if (pairs >= kBlockPairs && pauliengine_detail::use_threads(pairs)) {
                                if (auto result = pairwise_by_blocks(other, anti_commuting_only)) {
                                        return std::move(*result);
                                }
                        }
                        unweighted_QubitHamiltonian result;
                        const size_t n_buffer = std::max(longest(data_), longest(other.data_));
                        std::vector<uint64_t> out_x(n_buffer);
                        std::vector<uint64_t> out_y(n_buffer);
                        for (const auto& a : data_) {
                                for (const auto& b : other.data_) {
                                        if (anti_commuting_only && a.commutes(b)) {
                                                continue;
                                        }
                                        const int power_of_i = multiply_symplectic_words(
                                                a.x.data(), a.y.data(), a.x.size(),
                                                b.x.data(), b.y.data(), b.x.size(),
                                                out_x.data(), out_y.data());
                                        result.insert_words(out_x.data(), out_y.data(), std::max(a.x.size(), b.x.size()),
                                                            a.phase() + b.phase() + power_of_i);
                                }
                        }
                        return result;
                }
};
