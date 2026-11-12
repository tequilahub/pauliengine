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
                        for (const auto& term : hamiltonian.data) {
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

                unweighted_QubitHamiltonian pairwise(const unweighted_QubitHamiltonian& other, bool anti_commuting_only) const {
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
