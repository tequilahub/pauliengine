#pragma once
#include <unordered_map>
#include <string>
#include <complex>
#include <vector>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <memory>
#include <numeric>
#include <sstream>
#include <utility>

#include "pauliengine/PauliString.h"

#include <symengine/expression.h>
#include <symengine/functions.h>
using SymEngine::Expression;

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


#define BITS_IN_INTEGER (sizeof(uint64_t) * 8)

template<typename Coeff>
using Pauli_structure = std::pair<Coeff, std::unordered_map<int, std::string>>;

template<typename Coeff>
using Hamiltonian_structure = std::vector<Pauli_structure<Coeff>>;

using Matrix2D = std::vector<std::vector<std::complex<double>>>;

const std::complex<double> neg_I(0.0 , -1.0);

/// Tag type used by the fast-path constructor below: pass this when the caller
/// already guarantees the input vector is compacted (deduplicated + zero-filtered).
/// Lets operations like operator* skip a redundant compact() at the end.
struct already_compacted_t {};
inline constexpr already_compacted_t already_compacted{};

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

/// One term with a complex coefficient, read from a string or from flat storage.
/// x and y hold nw words each.
struct TermView {
        const uint64_t* x;
        const uint64_t* y;
        size_t nw;
        std::complex<double> coeff;
};

inline TermView view_of(const PauliString<std::complex<double>>& ps) {
        return {ps.x.data(), ps.y.data(), ps.x.size(), ps.coeff};
}

/// scale * i^power_of_i * a * b.
inline std::complex<double> product_coeff(std::complex<double> a, std::complex<double> b,
                                          int power_of_i, double scale) {
        static constexpr double kRealSign[4] = {1.0, -1.0, -1.0, 1.0};
        static constexpr double kImagSign[4] = {1.0, 1.0, -1.0, -1.0};
        const double parts[2] = {a.real() * b.real() - a.imag() * b.imag(),
                                 a.real() * b.imag() + a.imag() * b.real()};
        const int k = power_of_i & 3;
        return {scale * kRealSign[k] * parts[k & 1], scale * kImagSign[k] * parts[(k & 1) ^ 1]};
}

/// Whether two strings anti-commute: the symplectic product has odd parity. Much
/// cheaper than forming the product, which a commutator needs only for these pairs.
inline bool anticommute(const TermView& a, const TermView& b) {
        const size_t n = std::min(a.nw, b.nw);
        uint64_t anti = 0;
        for (size_t k = 0; k < n; ++k) anti ^= (a.x[k] & b.y[k]) ^ (b.x[k] & a.y[k]);
        return std::popcount(anti) & 1;
}

/// The terms of a Hamiltonian with complex coefficients, stored flat. Row i holds
/// the x words and then the y words of term i, both padded with zeros to the common
/// width; the coefficients are a second array of (real, imaginary) pairs. A sum of
/// n terms takes two allocations, instead of one string object per term plus heap
/// words for every string beyond 128 qubits.
class FlatTerms {
    public:
        FlatTerms() = default;

        /// n terms of nw words each, contents uninitialised.
        FlatTerms(size_t n, size_t nw)
                : n_(n), nw_(nw),
                  words_(n * nw > 0 ? new uint64_t[2 * n * nw] : nullptr),
                  coeffs_(n > 0 ? new double[2 * n] : nullptr) {}

        FlatTerms(const FlatTerms& o) : FlatTerms(o.n_, o.nw_) {
                std::copy(o.words_.get(), o.words_.get() + 2 * n_ * nw_, words_.get());
                std::copy(o.coeffs_.get(), o.coeffs_.get() + 2 * n_, coeffs_.get());
        }

        FlatTerms(FlatTerms&& o) noexcept
                : n_(std::exchange(o.n_, 0)), nw_(std::exchange(o.nw_, 0)),
                  words_(std::move(o.words_)), coeffs_(std::move(o.coeffs_)) {}

        FlatTerms& operator=(const FlatTerms& o) {
                if (this != &o) *this = FlatTerms(o);
                return *this;
        }

        FlatTerms& operator=(FlatTerms&& o) noexcept {
                n_ = std::exchange(o.n_, 0);
                nw_ = std::exchange(o.nw_, 0);
                words_ = std::move(o.words_);
                coeffs_ = std::move(o.coeffs_);
                return *this;
        }

        size_t size() const noexcept { return n_; }
        /// Words per term, the same for x and y.
        size_t width() const noexcept { return nw_; }

        uint64_t* x(size_t i) noexcept { return words_.get() + 2 * nw_ * i; }
        uint64_t* y(size_t i) noexcept { return x(i) + nw_; }
        const uint64_t* x(size_t i) const noexcept { return words_.get() + 2 * nw_ * i; }
        const uint64_t* y(size_t i) const noexcept { return x(i) + nw_; }

        std::complex<double> coeff(size_t i) const noexcept { return {coeffs_[2 * i], coeffs_[2 * i + 1]}; }
        void set_coeff(size_t i, std::complex<double> c) noexcept {
                coeffs_[2 * i] = c.real();
                coeffs_[2 * i + 1] = c.imag();
        }

        TermView view(size_t i) const noexcept { return {x(i), y(i), nw_, coeff(i)}; }

        /// Words of term i without its trailing zero words, as a PauliString keeps them.
        size_t used_words(size_t i) const noexcept {
                size_t n = nw_;
                while (n > 0 && x(i)[n - 1] == 0 && y(i)[n - 1] == 0) --n;
                return n;
        }

        /// Sets row i to the term `v`, which has at most width() words.
        void set_row(size_t i, const TermView& v) noexcept {
                uint64_t* rx = x(i);
                uint64_t* ry = y(i);
                for (size_t k = 0; k < nw_; ++k) {
                        rx[k] = k < v.nw ? v.x[k] : 0;
                        ry[k] = k < v.nw ? v.y[k] : 0;
                }
                set_coeff(i, v.coeff);
        }

        /// Sets row i to the operator part of a * b and returns the exponent k of the
        /// phase i^k; the coefficient is left to the caller.
        int set_row_product(size_t i, const TermView& a, const TermView& b) noexcept {
                uint64_t* rx = x(i);
                uint64_t* ry = y(i);
                const int power_of_i = multiply_symplectic_words(a.x, a.y, a.nw, b.x, b.y, b.nw, rx, ry);
                for (size_t k = std::max(a.nw, b.nw); k < nw_; ++k) {
                        rx[k] = 0;
                        ry[k] = 0;
                }
                return power_of_i;
        }

        PauliString<std::complex<double>> term(size_t i) const {
                PauliString<std::complex<double>> ps;
                const size_t nw = used_words(i);
                ps.x.assign(x(i), x(i) + nw);
                ps.y.assign(y(i), y(i) + nw);
                ps.coeff = coeff(i);
                return ps;
        }

        std::vector<PauliString<std::complex<double>>> terms() const {
                std::vector<PauliString<std::complex<double>>> out(n_);
                for_each_index(n_, use_threads(n_ * nw_), [&](size_t i) { out[i] = term(i); });
                return out;
        }

        /// Multiplies every coefficient by s; a real s scales both parts alone.
        template<typename S>
        void scale(S s) noexcept {
                for (size_t i = 0; i < n_; ++i) set_coeff(i, coeff(i) * s);
        }

        /// Removes the terms with coefficient zero, keeping the order of the others.
        void drop_zero_terms() noexcept {
                size_t write = 0;
                for (size_t read = 0; read < n_; ++read) {
                        if (coeffs_[2 * read] == 0.0 && coeffs_[2 * read + 1] == 0.0) continue;
                        if (write != read) {
                                std::copy(x(read), x(read) + 2 * nw_, x(write));
                                coeffs_[2 * write] = coeffs_[2 * read];
                                coeffs_[2 * write + 1] = coeffs_[2 * read + 1];
                        }
                        ++write;
                }
                n_ = write;
        }

        bool has_zero_coeff() const noexcept {
                for (size_t i = 0; i < n_; ++i) {
                        if (coeffs_[2 * i] == 0.0 && coeffs_[2 * i + 1] == 0.0) return true;
                }
                return false;
        }

    private:
        size_t n_ = 0;
        size_t nw_ = 0;
        std::unique_ptr<uint64_t[]> words_;
        std::unique_ptr<double[]> coeffs_;
};

}  // namespace pauliengine_detail

/// A sum of Pauli strings. Kept compacted: no two terms share an operator, and
/// zero terms are dropped (only scaling by zero leaves them in place).
///
/// With complex coefficients the terms are stored flat (pauliengine_detail::FlatTerms),
/// and products, commutators and sums work on that storage directly; terms() builds
/// Pauli strings from it on request. Symbolic coefficients are kept as a vector of
/// strings, since SymEngine expressions cannot be stored that way.
template<typename Coeff>
class QubitHamiltonian{
        static constexpr bool kFlat = std::is_same_v<Coeff, std::complex<double>>;
        using FlatTerms = pauliengine_detail::FlatTerms;
        using TermView = pauliengine_detail::TermView;
        using Storage = std::conditional_t<kFlat, FlatTerms, std::vector<PauliString<Coeff>>>;

        Storage store_;

        struct from_storage_t {};
        QubitHamiltonian(Storage&& store, from_storage_t) : store_(std::move(store)) {}

    public:
        QubitHamiltonian(const std::vector<PauliString<Coeff>>& terms) {
                if constexpr (kFlat) {
                        store_ = merge_strings(terms);
                } else {
                        store_ = terms;
                        compact_in_place();
                }
        }

        QubitHamiltonian(std::vector<PauliString<Coeff>>&& terms) {
                if constexpr (kFlat) {
                        store_ = merge_strings(terms);
                } else {
                        store_ = std::move(terms);
                        compact_in_place();
                }
        }

        /// Fast-path constructor: skips compact(). Only use when you can prove
        /// the vector has no duplicates and no zero-coefficient terms.
        QubitHamiltonian(std::vector<PauliString<Coeff>>&& terms, already_compacted_t) {
                if constexpr (kFlat) {
                        store_ = copy_strings(terms);
                } else {
                        store_ = std::move(terms);
                }
        }

        QubitHamiltonian(const Hamiltonian_structure<Coeff>& data)
                : QubitHamiltonian(strings_of(data)) {}

        /// Number of terms.
        size_t size() const { return store_.size(); }

        /// Term i as a Pauli string.
        PauliString<Coeff> term(size_t i) const {
                if constexpr (kFlat) {
                        return store_.term(i);
                } else {
                        return store_[i];
                }
        }

        /// All terms as Pauli strings. Complex Hamiltonians build them from the flat
        /// storage and return them by value; symbolic ones return their own vector.
        decltype(auto) terms() const {
                if constexpr (kFlat) {
                        return store_.terms();
                } else {
                        return (store_);
                }
        }

        QubitHamiltonian trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state) const {
                // The constructor drops zero-coefficient terms, including the sentinel
                // default-constructed PauliString returned for tracing X/Y against a
                // computational-basis state.
                return map_terms([&](const PauliString<Coeff>& ps) { return ps.trace_out_qubits(qubits, state); });
        }

        /// Overload accepting arbitrary single-qubit amplitudes (a, b) for each
        /// traced qubit (general superposition states).
        QubitHamiltonian trace_out_qubits(
                const std::vector<int>& qubits,
                const std::vector<std::pair<std::complex<double>, std::complex<double>>>& states
        ) const {
                return map_terms([&](const PauliString<Coeff>& ps) { return ps.trace_out_qubits(qubits, states); });
        }

        std::string to_string() const {
                std::string result = "";
                for (const auto& ps : terms()) {
                        result.append(ps.to_string() + "\n");
                }
                return result;
        }


        QubitHamiltonian operator*(std::complex<double> scale) const { return scaled(scale); }

        QubitHamiltonian operator*(int scale) const {
                if constexpr (kFlat) {
                        return scaled(static_cast<double>(scale));
                } else {
                        return scaled(scale);
                }
        }

        QubitHamiltonian operator*(double scale) const { return scaled(scale); }

        QubitHamiltonian operator*(const QubitHamiltonian& other) const {
                if constexpr (kFlat) {
                        return QubitHamiltonian(product_flat(store_, other.store_), from_storage_t{});
                } else {
                        std::vector<PauliString<Coeff>> product;
                        product.reserve(store_.size() * other.store_.size());
                        for (const auto& a : store_) {
                                for (const auto& b : other.store_) {
                                        product.push_back(a * b);
                                }
                        }
                        return QubitHamiltonian(std::move(product));
                }
        }

        QubitHamiltonian operator+(const QubitHamiltonian& other) const {
                if constexpr (kFlat) {
                        return QubitHamiltonian(add_flat(store_, other.store_, false), from_storage_t{});
                } else {
                        std::vector<PauliString<Coeff>> combined;
                        combined.reserve(store_.size() + other.store_.size());
                        combined.insert(combined.end(), store_.begin(), store_.end());
                        combined.insert(combined.end(), other.store_.begin(), other.store_.end());
                        return QubitHamiltonian(std::move(combined));  // constructor calls compact()
                }
        }

        QubitHamiltonian operator-(const QubitHamiltonian& other) const {
                if constexpr (kFlat) {
                        return QubitHamiltonian(add_flat(store_, other.store_, true), from_storage_t{});
                } else {
                        return *this + (other * -1.0);
                }
        }

        QubitHamiltonian operator-() const {
                return *this * -1.0;
        }

        bool operator==(const QubitHamiltonian& other) const {
                using Map = std::unordered_map<PauliString<Coeff>, Coeff,
                        PauliStringHash<Coeff>, PauliStringOperatorEqual<Coeff>>;
                Map first_data, second_data;
                for (const auto& ps : terms()) {
                        first_data[ps] = first_data[ps] + ps.coeff;
                }
                for (const auto& ps : other.terms()) {
                        second_data[ps] = second_data[ps] + ps.coeff;
                }
                if (first_data.size() != second_data.size()) {
                        return false;
                }
                for (const auto& [ps, coeff] : first_data) {
                        auto it = second_data.find(ps);
                        if (it == second_data.end() || !coeff_equal(coeff, it->second)) {
                                return false;
                        }
                }
                return true;
        }

        /// Coefficient of a single target Pauli operator in the commutator
        /// [this, other], computed as a hash-join instead of the full O(|A||B|)
        /// pair scan (PauliEngine Features, "targeted commutators"). For a fixed
        /// term P and target R there is exactly one operator Q that can produce R,
        /// namely v(Q) = v(P) XOR v(R) (equivalently Q ∝ P·R). We hash one side by
        /// its operator signature and, for each term of the other side, look up the
        /// unique matching partner, accumulating the commutator coefficient of the
        /// target. Cost O(min(|A|,|B|)) lookups after building the map.
        Coeff targeted_commutator(const QubitHamiltonian& other,
                                  const PauliString<Coeff>& target) const {
                using Map = std::unordered_map<PauliString<Coeff>, const PauliString<Coeff>*,
                        PauliStringHash<Coeff>, PauliStringOperatorEqual<Coeff>>;
                Coeff total(0.0);

                // Iterate the smaller side; hash the larger one (keyed by operator).
                const bool iter_this = size() <= other.size();
                const auto& these = terms();
                const auto& those = other.terms();
                const std::vector<PauliString<Coeff>>& outer = iter_this ? these : those;
                const std::vector<PauliString<Coeff>>& inner = iter_this ? those : these;

                Map lut;
                lut.reserve(inner.size());
                for (const auto& term : inner) {
                        lut.emplace(term, &term);
                }

                for (const auto& t : outer) {
                        // Operator of the unique partner that can produce `target`.
                        const PauliString<Coeff> partner_op = t * target;
                        const auto it = lut.find(partner_op);
                        if (it == lut.end()) {
                                continue;
                        }
                        // Preserve the [this, other] ordering (the commutator is
                        // antisymmetric): the term from `this` is always the left factor.
                        const PauliString<Coeff> c = iter_this
                                ? t.commutator(*it->second)
                                : it->second->commutator(t);
                        total = total + c.coeff;
                }
                return total;
        }

        QubitHamiltonian set_all_coeff(std::complex<double> value) const {
                Coeff coeff;
                if constexpr (kFlat) {
                        coeff = value;
                } else {
                        coeff = Expression(value);
                }
                return map_terms([&](const PauliString<Coeff>& ps) { return PauliString<Coeff>(ps.x, ps.y, coeff); });
        }

        QubitHamiltonian diff(std::string symbol) const {
                return map_terms([&](const PauliString<Coeff>& ps) { return ps.diff(symbol); });
        }

        QubitHamiltonian substitute(const std::unordered_map<std::string, std::complex<double>>& substitution_map) const {
                return map_terms([&](const PauliString<Coeff>& ps) { return ps.substitute(substitution_map); });
        }

        // Coefficient of the term whose operator part equals `ps`, or zero if the
        // Hamiltonian has no such term. The coefficient of `ps` itself is ignored,
        // so a query string of either coefficient type works. compact() keeps the
        // operator parts unique, hence at most one term can match.
        template<typename C>
        Coeff contains(const PauliString<C>& ps) const {
                if constexpr (kFlat) {
                        const size_t nw = ps.x.size();
                        for (size_t i = 0; i < store_.size(); ++i) {
                                if (store_.used_words(i) == nw
                                    && std::equal(ps.x.begin(), ps.x.end(), store_.x(i))
                                    && std::equal(ps.y.begin(), ps.y.end(), store_.y(i))) {
                                        return store_.coeff(i);
                                }
                        }
                } else {
                        for (const auto& term : store_) {
                                if (term.x == ps.x && term.y == ps.y) {
                                        return term.coeff;
                                }
                        }
                }
                return make_zero();
        }

        // {qubit: "X"} as in PauliString(coeff, map).
        Coeff contains(const std::unordered_map<int, std::string>& ops) const {
                std::unordered_map<int, std::string> normalized;
                for (const auto& [qubit, pauli] : ops) {
                        if (pauli.size() != 1) {
                                throw std::invalid_argument("contains: invalid Pauli operator '" + pauli + "'");
                        }
                        normalized[qubit] = std::string(1, normalize_pauli(pauli[0]));
                }
                return contains(PauliString<Coeff>(make_zero(), normalized));
        }

        // OpenFermion style [('X', 0), ('Z', 3)] as in PauliString((coeff, ops)).
        Coeff contains(const std::vector<std::pair<char, int>>& ops) const {
                std::vector<std::pair<char, int>> normalized;
                normalized.reserve(ops.size());
                for (const auto& [pauli, qubit] : ops) {
                        normalized.emplace_back(normalize_pauli(pauli), qubit);
                }
                return contains(PauliString<Coeff>(std::make_pair(make_zero(), normalized)));
        }

        // Either a dense string "XIZY" (character i acts on qubit i, as in
        // PauliString(coeff, string)) or a sparse one "X0 Z3" (as in the Python
        // PauliString factory). Any digit in the string selects the sparse form.
        Coeff contains(const std::string& ops) const {
                std::vector<std::pair<char, int>> parsed;
                const bool sparse = std::any_of(ops.begin(), ops.end(),
                        [](unsigned char c) { return std::isdigit(c); });
                if (sparse) {
                        std::istringstream tokens(ops);
                        std::string token;
                        while (tokens >> token) {
                                if (token.size() < 2 || !std::all_of(token.begin() + 1, token.end(),
                                                [](unsigned char c) { return std::isdigit(c); })) {
                                        throw std::invalid_argument("contains: invalid Pauli term '" + token + "', expected e.g. 'X0 Z3'");
                                }
                                parsed.emplace_back(token[0], std::stoi(token.substr(1)));
                        }
                } else {
                        for (size_t i = 0; i < ops.size(); ++i) {
                                if (!std::isspace(static_cast<unsigned char>(ops[i]))) {
                                        parsed.emplace_back(ops[i], static_cast<int>(i));
                                }
                        }
                }
                return contains(parsed);
        }


        QubitHamiltonian compact(){
                compact_in_place();
                return *this;
        }

        Hamiltonian_structure<Coeff> to_dictionary() const{
                const auto& ts = terms();
                Hamiltonian_structure<Coeff> output;
                output.reserve(ts.size());
                for (const auto& entry : ts) {
                        std::unordered_map<int, std::string> temp;
                        temp.reserve(entry.x.size() * BITS_IN_INTEGER);
                        for (size_t i = 0; i < entry.x.size(); ++i) {
                                uint64_t x_word = entry.x[i];
                                uint64_t y_word = entry.y[i];
                                for (int j = 0; j < BITS_IN_INTEGER; ++j) {
                                        uint64_t bit = (uint64_t)1 << j;
                                        bool x_bit = x_word & bit;
                                        bool y_bit = y_word & bit;
                                        if (!x_bit && !y_bit) continue;
                                        int index = i * BITS_IN_INTEGER + j;
                                        if (x_bit && y_bit) {
                                                temp[(int)index] = "Z";
                                        } else if (x_bit) {
                                                temp[(int)index] = "X";
                                        } else {
                                                temp[(int)index] = "Y";
                                        }
                                }
                        }
                        output.emplace_back(entry.coeff, std::move(temp));
                }
                return output;
        }

        std::vector<PauliString<Coeff>> to_list() const {
                return terms();
        }

        QubitHamiltonian commutator(const QubitHamiltonian& other) const {
                if constexpr (kFlat) {
                        return QubitHamiltonian(commutator_flat(store_, other.store_), from_storage_t{});
                } else {
                        std::vector<PauliString<Coeff>> terms;
                        terms.reserve(store_.size() * other.store_.size() / 8 + 1);
                        for (const auto& a : store_) {
                                for (const auto& b : other.store_) {
                                        PauliString<Coeff> t = a.commutator(b);
                                        if (!t.x.empty()) {  // commuting pairs give the empty string
                                                terms.push_back(std::move(t));
                                        }
                                }
                        }
                        return QubitHamiltonian(std::move(terms));  // constructor calls compact()
                }
        }

        std::vector<int> qubits() const {
                std::set<int> qs;
                const auto add = [&qs](const uint64_t* x, const uint64_t* y, size_t nw) {
                        for (size_t i = 0; i < nw; ++i) {
                                uint64_t xy = x[i] | y[i];
                                while (xy) {
                                        int j = std::countr_zero(xy);
                                        qs.insert((int)(i * BITS_IN_INTEGER + j));
                                        xy &= xy - 1;
                                }
                        }
                };
                if constexpr (kFlat) {
                        for (size_t i = 0; i < store_.size(); ++i) add(store_.x(i), store_.y(i), store_.width());
                } else {
                        for (const auto& ps : store_) add(ps.x.data(), ps.y.data(), ps.x.size());
                }
                return std::vector<int>(qs.begin(), qs.end());
        }

        int n_qubits() const { return (int)qubits().size(); }

        bool is_all_z() const {
                if constexpr (kFlat) {
                        for (size_t i = 0; i < store_.size(); ++i) {
                                if (!std::equal(store_.x(i), store_.x(i) + store_.width(), store_.y(i))) return false;
                        }
                } else {
                        for (const auto& ps : store_) {
                                if (!ps.is_all_z()) return false;
                        }
                }
                return true;
        }

        int count_measurements() const {
                return is_all_z() ? 1 : (int)size();
        }

        bool is_hermitian() const {
                if constexpr (kFlat) {
                        for (size_t i = 0; i < store_.size(); ++i) {
                                if (std::abs(store_.coeff(i).imag()) > 1e-12) return false;
                        }
                        return true;
                } else {
                        try {
                                for (const auto& ps : store_) {
                                        auto c = PauliString<Coeff>::to_complex(ps.coeff);
                                        if (std::abs(c.imag()) > 1e-12) return false;
                                }
                                return true;
                        } catch (...) {
                                return false;
                        }
                }
        }

        bool is_antihermitian() const {
                if constexpr (kFlat) {
                        for (size_t i = 0; i < store_.size(); ++i) {
                                if (std::abs(store_.coeff(i).real()) > 1e-12) return false;
                        }
                        return true;
                } else {
                        try {
                                for (const auto& ps : store_) {
                                        auto c = PauliString<Coeff>::to_complex(ps.coeff);
                                        if (std::abs(c.real()) > 1e-12) return false;
                                }
                                return true;
                        } catch (...) {
                                return false;
                        }
                }
        }

        std::vector<PauliString<Coeff>> paulistrings() const {
                return terms();
        }

        static QubitHamiltonian zero() {
                return QubitHamiltonian(std::vector<PauliString<Coeff>>{});
        }

        static QubitHamiltonian unit() {
                Coeff one = make_one();
                PauliString<Coeff> identity_ps(std::vector<uint64_t>{}, std::vector<uint64_t>{}, one);
                return QubitHamiltonian(std::vector<PauliString<Coeff>>{identity_ps});
        }

        QubitHamiltonian dagger() const {
                return map_terms([](const PauliString<Coeff>& ps) {
                        return PauliString<Coeff>(ps.x, ps.y, conjugate_coeff(ps.coeff));
                });
        }

        QubitHamiltonian conjugate() const {
                return map_terms([](const PauliString<Coeff>& ps) {
                        Coeff new_coeff = conjugate_coeff(ps.coeff);
                        if (ps.count_y() % 2 == 1) {
                                new_coeff = new_coeff * -1.0;
                        }
                        return PauliString<Coeff>(ps.x, ps.y, new_coeff);
                });
        }

        QubitHamiltonian transpose() const {
                return map_terms([](const PauliString<Coeff>& ps) {
                        Coeff new_coeff = ps.coeff;
                        if (ps.count_y() % 2 == 1) {
                                new_coeff = new_coeff * -1.0;
                        }
                        return PauliString<Coeff>(ps.x, ps.y, new_coeff);
                });
        }

        QubitHamiltonian simplify(double threshold = 0.0) const {
                std::vector<PauliString<Coeff>> result;
                for (const auto& ps : terms()) {
                        bool keep;
                        if constexpr (kFlat) {
                                keep = std::abs(ps.coeff) > threshold;
                        } else {
                                try {
                                        auto c = PauliString<Coeff>::to_complex(ps.coeff);
                                        keep = std::abs(c) > threshold;
                                } catch (...) {
                                        keep = true;
                                }
                        }
                        if (keep) result.push_back(ps);
                }
                return QubitHamiltonian(std::move(result));
        }

        std::pair<QubitHamiltonian, QubitHamiltonian> split() const {
                std::vector<PauliString<Coeff>> herm_terms;
                std::vector<PauliString<Coeff>> anti_terms;
                for (const auto& ps : terms()) {
                        std::complex<double> c;
                        if constexpr (kFlat) {
                                c = ps.coeff;
                        } else {
                                try {
                                        c = PauliString<Coeff>::to_complex(ps.coeff);
                                } catch (...) {
                                        throw std::runtime_error("split() requires fully evaluated numeric coefficients");
                                }
                        }
                        if (c.real() != 0.0) {
                                Coeff real_coeff;
                                if constexpr (kFlat) {
                                        real_coeff = std::complex<double>(c.real(), 0.0);
                                } else {
                                        real_coeff = Expression(c.real());
                                }
                                herm_terms.push_back(PauliString<Coeff>(ps.x, ps.y, real_coeff));
                        }
                        if (c.imag() != 0.0) {
                                Coeff imag_coeff;
                                if constexpr (kFlat) {
                                        imag_coeff = std::complex<double>(0.0, c.imag());
                                } else {
                                        imag_coeff = Expression(std::complex<double>(0.0, c.imag()));
                                }
                                anti_terms.push_back(PauliString<Coeff>(ps.x, ps.y, imag_coeff));
                        }
                }
                return {QubitHamiltonian(std::move(herm_terms)), QubitHamiltonian(std::move(anti_terms))};
        }

        QubitHamiltonian map_qubits(std::unordered_map<int, int> qubit_map) const {
                return map_terms([&](const PauliString<Coeff>& ps) { return ps.map_qubits(qubit_map); });
        }

        QubitHamiltonian power(int n) const {
                if (n < 0) {
                        throw std::runtime_error("Negative power not supported for QubitHamiltonian");
                }
                if (n == 0) return unit();
                QubitHamiltonian result = *this;
                for (int i = 1; i < n; ++i) {
                        result = result * *this;
                }
                return result;
        }

        Matrix2D to_matrix(bool ignore_unused_qubits = true) const {
                std::vector<int> q_list = qubits();
                int nq;
                if (ignore_unused_qubits) {
                        nq = (int)q_list.size();
                } else if (q_list.empty()) {
                        nq = 0;
                } else {
                        nq = q_list.back() + 1;
                }

                size_t dim = (nq == 0) ? 1 : ((size_t)1 << nq);
                Matrix2D total(dim, std::vector<std::complex<double>>(dim, 0.0));

                std::unordered_map<int, int> q_to_pos;
                if (ignore_unused_qubits) {
                        for (int i = 0; i < (int)q_list.size(); ++i) {
                                q_to_pos[q_list[i]] = i;
                        }
                }

                for (const auto& ps : terms()) {
                        std::complex<double> c;
                        if constexpr (kFlat) {
                                c = ps.coeff;
                        } else {
                                c = PauliString<Coeff>::to_complex(ps.coeff);
                        }

                        auto dict = ps.to_dictionary().second;
                        std::unordered_map<int, std::string> remapped;
                        for (const auto& [q, op] : dict) {
                                int pos = ignore_unused_qubits ? q_to_pos.at(q) : q;
                                remapped[pos] = op;
                        }

                        Matrix2D mat = {{std::complex<double>(1.0, 0.0)}};
                        for (int q = 0; q < nq; ++q) {
                                auto it = remapped.find(q);
                                std::string op = (it != remapped.end()) ? it->second : "I";
                                mat = tensor_product(mat, get_pauli_matrix(op));
                        }

                        for (size_t i = 0; i < dim; ++i) {
                                for (size_t j = 0; j < dim; ++j) {
                                        total[i][j] += c * mat[i][j];
                                }
                        }
                }
                return total;
        }

    private:
        static std::vector<PauliString<Coeff>> strings_of(const Hamiltonian_structure<Coeff>& data) {
                std::vector<PauliString<Coeff>> strings;
                strings.reserve(data.size());
                for (const Pauli_structure<Coeff>& entry : data) {
                        strings.push_back(PauliString<Coeff>(entry.first, entry.second));
                }
                return strings;
        }

        /// The Hamiltonian of f(term) over all terms, compacted.
        template<typename F>
        QubitHamiltonian map_terms(F&& f) const {
                const auto& ts = terms();
                std::vector<PauliString<Coeff>> mapped;
                mapped.reserve(ts.size());
                for (const auto& ps : ts) {
                        mapped.push_back(f(ps));
                }
                return QubitHamiltonian(std::move(mapped));
        }

        template<typename S>
        QubitHamiltonian scaled(S scale) const {
                if constexpr (kFlat) {
                        // Stored terms are compacted, apart from the zeros an earlier scaling
                        // by zero leaves; those are dropped as compact() would.
                        FlatTerms out = store_;
                        if (out.has_zero_coeff()) out.drop_zero_terms();
                        out.scale(scale);
                        return QubitHamiltonian(std::move(out), from_storage_t{});
                } else {
                        QubitHamiltonian out(store_);
                        for (auto& entry : out.store_) {
                                entry.coeff *= scale;
                        }
                        return out;
                }
        }

        static Coeff conjugate_coeff(const Coeff& c) {
                if constexpr (kFlat) {
                        return std::conj(c);
                } else {
                        return Expression(SymEngine::conjugate(c.get_basic()));
                }
        }

        void compact_in_place() {
                if constexpr (kFlat) {
                        store_ = merge_views(store_.size(), store_.width(), [this](size_t i) { return store_.view(i); });
                } else {
                        compact_symbolic();
                }
        }

        static bool coeff_is_zero(const Coeff& c) {
                if constexpr (kFlat) {
                        return c == std::complex<double>(0.0, 0.0);
                } else {
                    const auto& basic = *c.get_basic();
                        return SymEngine::is_a_Number(basic)
                                && SymEngine::down_cast<const SymEngine::Number&>(basic).is_zero();
                }
        }

        static bool coeff_equal(const Coeff& a, const Coeff& b) {
                if constexpr (kFlat) {
                        return a == b;
                } else {
                        return coeff_is_zero(SymEngine::expand(a - b));
                }
        }

        // ------------------------------------------------------------------
        // Complex coefficients: sorting and merging on flat storage.
        //
        // Every operation forms records (operator words plus coefficient), sorts
        // them by operator, adds the coefficients of equal operators and writes the
        // non-zero sums as flat storage. Up to kMaxSmallKeyWords words the records
        // are fixed-size structs sorted by radix sort; longer strings are sorted
        // through an index, since moving whole records would cost too much.
        // ------------------------------------------------------------------

        /// Largest word count handled with fixed-size records (8 words = 512 qubits).
        static constexpr size_t kMaxSmallKeyWords = 8;

        /// Compile-time dispatch: calls f.operator()<W>() for the smallest W in
        /// {1, 2, 4, ..., kMaxSmallKeyWords} that fits max_words.
        /// Returns false when max_words exceeds the fast-path limit.
        template<size_t W = 1, typename F>
        static bool with_small_key(size_t max_words, F&& f) {
                if constexpr (W > kMaxSmallKeyWords) {
                        return false;
                } else {
                        if (max_words <= W) {
                                f.template operator()<W>();
                                return true;
                        }
                        return with_small_key<W * 2>(max_words, std::forward<F>(f));
                }
        }

        /// Fixed-size, trivially copyable record of one term.
        template<size_t W>
        struct KeyedTerm {
                uint64_t key[2 * W];
                double re;
                double im;
        };

        template<size_t W>
        static void set_keyed_term(KeyedTerm<W>& term, const TermView& v) {
                for (size_t k = 0; k < W; ++k) {
                        term.key[k]     = (k < v.nw) ? v.x[k] : 0;
                        term.key[W + k] = (k < v.nw) ? v.y[k] : 0;
                }
                term.re = v.coeff.real();
                term.im = v.coeff.imag();
        }

        /// Sets `term` to scale * a * b.
        template<size_t W>
        static void set_keyed_product(KeyedTerm<W>& term, const TermView& a, const TermView& b, double scale) {
                const int power_of_i = multiply_symplectic_words(a.x, a.y, a.nw, b.x, b.y, b.nw, term.key, term.key + W);
                for (size_t k = std::max(a.nw, b.nw); k < W; ++k) {
                        term.key[k] = 0;
                        term.key[W + k] = 0;
                }
                const std::complex<double> c = pauliengine_detail::product_coeff(a.coeff, b.coeff, power_of_i, scale);
                term.re = c.real();
                term.im = c.imag();
        }

        /// Below this many terms a range is sorted by insertion, which needs no spare memory.
        static constexpr size_t kInsertionSortLimit = 32;

        /// Stable MSD radix sort by key, starting at key byte `byte`. The n terms are in
        /// `src`, `other` is scratch space of the same size, and the result ends up in
        /// `src` if final_in_src, else in `other`. Levels where all terms share the digit
        /// move nothing, and the two arrays swap roles per level instead of copying back.
        /// Stability keeps equal operators in input order, so that adding them rounds the
        /// same way for any number of threads.
        template<size_t W>
        static void radix_sort_keyed(KeyedTerm<W>* src, KeyedTerm<W>* other, size_t n, size_t byte, bool final_in_src) {
                using Term = KeyedTerm<W>;
                if (n <= kInsertionSortLimit) {
                        const auto key_less = [](const Term& a, const Term& b) {
                                for (size_t k = 0; k < 2 * W; ++k) {
                                        if (a.key[k] != b.key[k]) return a.key[k] < b.key[k];
                                }
                                return false;
                        };
                        for (size_t i = 1; i < n; ++i) {
                                const Term moved = src[i];
                                size_t j = i;
                                for (; j > 0 && key_less(moved, src[j - 1]); --j) {
                                        src[j] = src[j - 1];
                                }
                                src[j] = moved;
                        }
                        if (!final_in_src) std::copy(src, src + n, other);
                        return;
                }

                const auto digit = [](const Term& t, size_t b) {
                        return static_cast<size_t>((t.key[b / 8] >> (56 - 8 * (b % 8))) & 0xFF);
                };
                size_t counts[256];
                for (;; ++byte) {
                        if (byte == sizeof(Term::key)) {  // all keys are equal
                                if (!final_in_src) std::copy(src, src + n, other);
                                return;
                        }
                        std::fill(counts, counts + 256, size_t{0});
                        for (size_t i = 0; i < n; ++i) ++counts[digit(src[i], byte)];
                        if (counts[digit(src[0], byte)] != n) break;
                }

                size_t next[256];
                size_t position = 0;
                for (size_t d = 0; d < 256; ++d) {
                        next[d] = position;
                        position += counts[d];
                }
                for (size_t i = 0; i < n; ++i) {
                        other[next[digit(src[i], byte)]++] = src[i];
                }
                position = 0;
                for (size_t d = 0; d < 256; ++d) {
                        if (counts[d] > 0) {
                                radix_sort_keyed<W>(other + position, src + position, counts[d], byte + 1, !final_in_src);
                                position += counts[d];
                        }
                }
        }

        /// Sorts the n terms at `first` and adds the coefficients of equal operators.
        /// The non-zero results are left at `first`; returns how many there are.
        template<size_t W>
        static size_t sort_and_merge_keyed(KeyedTerm<W>* first, KeyedTerm<W>* spare, size_t n) {
                using Term = KeyedTerm<W>;
                radix_sort_keyed<W>(first, spare, n, 0, true);
                const auto key_equal = [](const Term& a, const Term& b) {
                        for (size_t k = 0; k < 2 * W; ++k) {
                                if (a.key[k] != b.key[k]) return false;
                        }
                        return true;
                };
                Term* const last = first + n;
                Term* write = first;
                for (Term* read = first; read != last;) {
                        double re = read->re;
                        double im = read->im;
                        Term* next = read + 1;
                        for (; next != last && key_equal(*read, *next); ++next) {
                                re += next->re;
                                im += next->im;
                        }
                        if (re != 0.0 || im != 0.0) {
                                if (write != read) *write = *read;
                                write->re = re;
                                write->im = im;
                                ++write;
                        }
                        read = next;
                }
                return static_cast<size_t>(write - first);
        }

        /// Writes `count` records to rows first, first + 1, ... of `out`.
        template<size_t W>
        static void write_keyed_terms(const KeyedTerm<W>* terms, size_t count, FlatTerms& out, size_t first) {
                const size_t nw = out.width();
                for (size_t i = 0; i < count; ++i) {
                        const KeyedTerm<W>& term = terms[i];
                        uint64_t* const x = out.x(first + i);
                        uint64_t* const y = out.y(first + i);
                        for (size_t k = 0; k < nw; ++k) {
                                x[k] = term.key[k];
                                y[k] = term.key[W + k];
                        }
                        out.set_coeff(first + i, {term.re, term.im});
                }
        }

        /// Sorts `terms` by operator, adds the coefficients of equal operators and
        /// returns the non-zero results, ordered by operator, as flat storage of width
        /// nw (no record has a non-zero word beyond nw). Large inputs use a sample
        /// sort: regular samples give splitters, every term moves to the bucket of its
        /// operator, and buckets are then sorted and merged independently. Equal
        /// operators always share a bucket, so no serial merge step is left. The
        /// result does not depend on the number of threads. `terms` is overwritten.
        template<size_t W>
        static FlatTerms merge_keyed_terms(KeyedTerm<W>* terms, size_t n, size_t nw) {
                using Term = KeyedTerm<W>;
#ifdef PAULIENGINE_HAS_OPENMP
                if (pauliengine_detail::use_threads(n)) {
                        const int nthreads = omp_get_max_threads();
                        const auto key_less = [](const Term& a, const Term& b) {
                                for (size_t k = 0; k < 2 * W; ++k) {
                                        if (a.key[k] != b.key[k]) return a.key[k] < b.key[k];
                                }
                                return false;
                        };
                        constexpr int kBucketsPerThread = 4;  // balances uneven buckets
                        constexpr size_t kOversampling = 32;
                        const int n_buckets = kBucketsPerThread * nthreads;
                        const size_t n_samples = static_cast<size_t>(n_buckets) * kOversampling;

                        std::vector<Term> samples(n_samples);
                        for (size_t i = 0; i < n_samples; ++i) {
                                samples[i] = terms[i * n / n_samples];
                        }
                        std::sort(samples.begin(), samples.end(), key_less);
                        std::vector<Term> splitters(n_buckets - 1);
                        for (int b = 0; b + 1 < n_buckets; ++b) {
                                splitters[b] = samples[(b + 1) * kOversampling];
                        }

                        std::unique_ptr<uint16_t[]> bucket_of(new uint16_t[n]);
                        std::vector<size_t> cursor(static_cast<size_t>(nthreads) * n_buckets, 0);
                        #pragma omp parallel for schedule(static)
                        for (int t = 0; t < nthreads; ++t) {
                                size_t* counts = cursor.data() + static_cast<size_t>(t) * n_buckets;
                                const size_t last = n / nthreads * (t + 1) + (t + 1 == nthreads ? n % nthreads : 0);
                                for (size_t i = n / nthreads * t; i < last; ++i) {
                                        const auto b = std::upper_bound(splitters.begin(), splitters.end(), terms[i], key_less)
                                                - splitters.begin();
                                        bucket_of[i] = static_cast<uint16_t>(b);
                                        ++counts[b];
                                }
                        }

                        std::vector<size_t> bucket_begin(static_cast<size_t>(n_buckets) + 1, 0);
                        size_t position = 0;
                        for (int b = 0; b < n_buckets; ++b) {
                                bucket_begin[b] = position;
                                for (int t = 0; t < nthreads; ++t) {
                                        size_t& entry = cursor[static_cast<size_t>(t) * n_buckets + b];
                                        const size_t count = entry;
                                        entry = position;
                                        position += count;
                                }
                        }
                        bucket_begin[n_buckets] = position;

                        // Chunks are walked in order, so every bucket keeps the input order.
                        std::unique_ptr<Term[]> partitioned(new Term[n]);
                        Term* const sorted = partitioned.get();
                        #pragma omp parallel for schedule(static)
                        for (int t = 0; t < nthreads; ++t) {
                                size_t* next = cursor.data() + static_cast<size_t>(t) * n_buckets;
                                const size_t last = n / nthreads * (t + 1) + (t + 1 == nthreads ? n % nthreads : 0);
                                for (size_t i = n / nthreads * t; i < last; ++i) {
                                        sorted[next[bucket_of[i]]++] = terms[i];
                                }
                        }

                        std::vector<size_t> kept(n_buckets, 0);
                        #pragma omp parallel for schedule(dynamic)
                        for (int b = 0; b < n_buckets; ++b) {
                                kept[b] = sort_and_merge_keyed<W>(sorted + bucket_begin[b], terms + bucket_begin[b],
                                                                  bucket_begin[b + 1] - bucket_begin[b]);
                        }

                        std::vector<size_t> out_begin(static_cast<size_t>(n_buckets) + 1, 0);
                        for (int b = 0; b < n_buckets; ++b) out_begin[b + 1] = out_begin[b] + kept[b];
                        // The storage is left untouched until here, so the threads fault its pages in.
                        FlatTerms out(out_begin[n_buckets], nw);
                        #pragma omp parallel for schedule(dynamic)
                        for (int b = 0; b < n_buckets; ++b) {
                                write_keyed_terms<W>(sorted + bucket_begin[b], kept[b], out, out_begin[b]);
                        }
                        return out;
                }
#endif
                std::unique_ptr<Term[]> spare(n > kInsertionSortLimit ? new Term[n] : nullptr);
                const size_t kept = sort_and_merge_keyed<W>(terms, spare.get(), n);
                FlatTerms out(kept, nw);
                write_keyed_terms<W>(terms, kept, out, 0);
                return out;
        }

        /// Sorts the rows of `records` listed in `order` by operator, adds the
        /// coefficients of equal operators and returns the non-zero results. For
        /// strings beyond kMaxSmallKeyWords words: only the index is sorted.
        static FlatTerms merge_long(const FlatTerms& records, std::vector<size_t>& order) {
                using pauliengine_detail::for_each_index;
                using pauliengine_detail::use_threads;
                const size_t row_words = 2 * records.width();
                const auto same_operator = [&](size_t i, size_t j) {
                        return std::equal(records.x(i), records.x(i) + row_words, records.x(j));
                };
                // Equal operators stay in input order, so their sum rounds the same way
                // for any number of threads.
                const auto less = [&](size_t i, size_t j) {
                        const uint64_t* a = records.x(i);
                        const uint64_t* b = records.x(j);
                        for (size_t k = 0; k < row_words; ++k) {
                                if (a[k] != b[k]) return a[k] < b[k];
                        }
                        return i < j;
                };
                parallel_sort_vec(order, less, use_threads(order.size()));

                std::vector<size_t> heads;
                std::vector<std::complex<double>> sums;
                for (size_t k = 0; k < order.size();) {
                        std::complex<double> sum = records.coeff(order[k]);
                        size_t next = k + 1;
                        for (; next < order.size() && same_operator(order[k], order[next]); ++next) {
                                sum += records.coeff(order[next]);
                        }
                        if (sum != 0.0) {
                                heads.push_back(order[k]);
                                sums.push_back(sum);
                        }
                        k = next;
                }

                FlatTerms out(heads.size(), records.width());
                for_each_index(heads.size(), use_threads(heads.size() * records.width()), [&](size_t i) {
                        std::copy(records.x(heads[i]), records.x(heads[i]) + row_words, out.x(i));
                        out.set_coeff(i, sums[i]);
                });
                return out;
        }

        /// Merges the n terms get(0), ..., get(n - 1), each of at most nw words.
        template<typename Get>
        static FlatTerms merge_views(size_t n, size_t nw, Get&& get) {
                using pauliengine_detail::for_each_index;
                using pauliengine_detail::use_threads;
                FlatTerms out;
                if (with_small_key(nw, [&]<size_t W>() {
                        std::unique_ptr<KeyedTerm<W>[]> records(new KeyedTerm<W>[n]);
                        KeyedTerm<W>* const raw = records.get();
                        for_each_index(n, use_threads(n * W), [&](size_t i) { set_keyed_term<W>(raw[i], get(i)); });
                        out = merge_keyed_terms<W>(raw, n, nw);
                })) {
                        return out;
                }
                FlatTerms records(n, nw);
                for_each_index(n, use_threads(n * nw), [&](size_t i) { records.set_row(i, get(i)); });
                std::vector<size_t> order(n);
                std::iota(order.begin(), order.end(), size_t{0});
                return merge_long(records, order);
        }

        /// The largest word count among the strings.
        static size_t widest(const std::vector<PauliString<Coeff>>& strings) {
                size_t nw = 0;
                for (const auto& ps : strings) nw = std::max(nw, ps.x.size());
                return nw;
        }

        /// The strings, compacted.
        static FlatTerms merge_strings(const std::vector<PauliString<Coeff>>& strings) {
                return merge_views(strings.size(), widest(strings),
                                   [&](size_t i) { return pauliengine_detail::view_of(strings[i]); });
        }

        /// The strings in their order, without compacting.
        static FlatTerms copy_strings(const std::vector<PauliString<Coeff>>& strings) {
                const size_t nw = widest(strings);
                FlatTerms out(strings.size(), nw);
                pauliengine_detail::for_each_index(strings.size(), pauliengine_detail::use_threads(strings.size() * nw),
                        [&](size_t i) { out.set_row(i, pauliengine_detail::view_of(strings[i])); });
                return out;
        }

        /// a + b, or a - b if `subtract`.
        static FlatTerms add_flat(const FlatTerms& a, const FlatTerms& b, bool subtract) {
                const size_t na = a.size();
                return merge_views(na + b.size(), std::max(a.width(), b.width()), [&](size_t i) {
                        if (i < na) return a.view(i);
                        TermView v = b.view(i - na);
                        if (subtract) v.coeff = -v.coeff;
                        return v;
                });
        }

        static FlatTerms product_flat(const FlatTerms& a, const FlatTerms& b) {
                using pauliengine_detail::for_each_index;
                using pauliengine_detail::use_threads;
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                const size_t pairs = n1 * n2;
                const size_t nw = std::max(a.width(), b.width());
                FlatTerms out;
                if (with_small_key(nw, [&]<size_t W>() {
                        // new[] leaves the records untouched, so the threads fault their own pages in.
                        std::unique_ptr<KeyedTerm<W>[]> records(new KeyedTerm<W>[pairs]);
                        KeyedTerm<W>* const raw = records.get();
                        for_each_index(n1, use_threads(pairs * W), [&](size_t i) {
                                const TermView left = a.view(i);
                                KeyedTerm<W>* const row = raw + i * n2;
                                for (size_t j = 0; j < n2; ++j) {
                                        set_keyed_product<W>(row[j], left, b.view(j), 1.0);
                                }
                        });
                        out = merge_keyed_terms<W>(raw, pairs, nw);
                })) {
                        return out;
                }
                FlatTerms records(pairs, nw);
                for_each_index(n1, use_threads(pairs * nw), [&](size_t i) {
                        const TermView left = a.view(i);
                        for (size_t j = 0; j < n2; ++j) {
                                const TermView right = b.view(j);
                                const int power_of_i = records.set_row_product(i * n2 + j, left, right);
                                records.set_coeff(i * n2 + j, pauliengine_detail::product_coeff(left.coeff, right.coeff, power_of_i, 1.0));
                        }
                });
                std::vector<size_t> order(pairs);
                std::iota(order.begin(), order.end(), size_t{0});
                return merge_long(records, order);
        }

        /// [a, b] for strings of at most W words, written with width nw.
        template<size_t W>
        static FlatTerms commutator_small_key(const FlatTerms& a, const FlatTerms& b, size_t nw) {
                using Term = KeyedTerm<W>;
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                const size_t pairs = n1 * n2;
                // Appends the commutators of row i with all of b that do not vanish.
                const auto add_row = [&](size_t i, std::vector<Term>& mine) {
                        const TermView left = a.view(i);
                        for (size_t j = 0; j < n2; ++j) {
                                const TermView right = b.view(j);
                                if (pauliengine_detail::anticommute(left, right)) {
                                        set_keyed_product<W>(mine.emplace_back(), left, right, 2.0);
                                }
                        }
                };
#ifdef PAULIENGINE_HAS_OPENMP
                if (pauliengine_detail::use_threads(pairs * W)) {
                        const int nthreads = omp_get_max_threads();
                        const std::ptrdiff_t rows = static_cast<std::ptrdiff_t>(n1);  // signed for OpenMP 2.0
                        std::vector<std::vector<Term>> local(nthreads);
                        #pragma omp parallel
                        {
                                // Filled as a local object: the vectors in `local` share cache lines.
                                std::vector<Term> mine;
                                mine.reserve(pairs / (8 * static_cast<size_t>(nthreads)) + 1);
                                #pragma omp for schedule(static)
                                for (std::ptrdiff_t i = 0; i < rows; ++i) add_row(static_cast<size_t>(i), mine);
                                local[omp_get_thread_num()] = std::move(mine);
                        }
                        // A static schedule hands out ascending rows, so this is the serial order.
                        std::vector<size_t> offset(static_cast<size_t>(nthreads) + 1, 0);
                        for (int t = 0; t < nthreads; ++t) offset[t + 1] = offset[t] + local[t].size();
                        std::unique_ptr<Term[]> records(new Term[offset[nthreads]]);
                        Term* const raw = records.get();
                        pauliengine_detail::for_each_index(static_cast<size_t>(nthreads), true, [&](size_t t) {
                                std::copy(local[t].begin(), local[t].end(), raw + offset[t]);
                                std::vector<Term>().swap(local[t]);
                        });
                        return merge_keyed_terms<W>(raw, offset[nthreads], nw);
                }
#endif
                std::vector<Term> mine;
                mine.reserve(pairs / 8 + 1);
                for (size_t i = 0; i < n1; ++i) add_row(i, mine);
                return merge_keyed_terms<W>(mine.data(), mine.size(), nw);
        }

        /// [a, b] = ab - ba: twice the product of every anti-commuting pair.
        static FlatTerms commutator_flat(const FlatTerms& a, const FlatTerms& b) {
                using pauliengine_detail::anticommute;
                using pauliengine_detail::for_each_index;
                using pauliengine_detail::use_threads;
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                const size_t pairs = n1 * n2;
                const size_t nw = std::max(a.width(), b.width());
                FlatTerms out;
                if (with_small_key(nw, [&]<size_t W>() { out = commutator_small_key<W>(a, b, nw); })) {
                        return out;
                }
                // Long strings: row i of a writes its results to the start of its own
                // block of n2 records; the rest of the block is never touched.
                FlatTerms records(pairs, nw);
                std::vector<size_t> found(n1, 0);
                for_each_index(n1, use_threads(pairs * nw), [&](size_t i) {
                        const TermView left = a.view(i);
                        size_t r = i * n2;
                        for (size_t j = 0; j < n2; ++j) {
                                const TermView right = b.view(j);
                                if (!anticommute(left, right)) continue;
                                const int power_of_i = records.set_row_product(r, left, right);
                                records.set_coeff(r, pauliengine_detail::product_coeff(left.coeff, right.coeff, power_of_i, 2.0));
                                ++r;
                        }
                        found[i] = r - i * n2;
                });
                std::vector<size_t> order;
                order.reserve(std::accumulate(found.begin(), found.end(), size_t{0}));
                for (size_t i = 0; i < n1; ++i) {
                        for (size_t r = i * n2; r < i * n2 + found[i]; ++r) order.push_back(r);
                }
                return merge_long(records, order);
        }

        // ------------------------------------------------------------------
        // Symbolic coefficients
        // ------------------------------------------------------------------

        void compact_symbolic() {
                if (store_.empty()) return;
                if (with_small_key(widest(store_), [this]<size_t W>() { compact_small_key_symbolic<W>(); })) {
                        return;
                }

                // sort-based compact. Single-pass ordering: by word
                // count first (constructors trim trailing zero words, so equal
                // operators always have equal sizes), then by content.
                auto op_less = [](const PauliString<Coeff>& a, const PauliString<Coeff>& b) {
                        const size_t na = a.x.size();
                        const size_t nb = b.x.size();
                        if (na != nb) return na < nb;
                        for (size_t i = 0; i < na; ++i) {
                                if (a.x[i] != b.x[i]) return a.x[i] < b.x[i];
                        }
                        for (size_t i = 0; i < na; ++i) {
                                if (a.y[i] != b.y[i]) return a.y[i] < b.y[i];
                        }
                        return false;
                };
                auto op_equal = [](const PauliString<Coeff>& a, const PauliString<Coeff>& b) {
                        return a.x == b.x && a.y == b.y;
                };

                // symengine coefficients are not safe to move across threads
                std::sort(store_.begin(), store_.end(), op_less);

                size_t write = 0;
                for (size_t read = 1; read < store_.size(); ++read) {
                        if (op_equal(store_[write], store_[read])) {
                                store_[write].coeff = store_[write].coeff + store_[read].coeff;
                        } else {
                                ++write;
                                if (write != read) {
                                        store_[write] = std::move(store_[read]);
                                }
                        }
                }
                store_.resize(write + 1);

                store_.erase(std::remove_if(store_.begin(), store_.end(),
                        [](const PauliString<Coeff>& ps) { return coeff_is_zero(ps.coeff); }),
                        store_.end());
        }

        template<size_t W>
        void compact_small_key_symbolic() {
                struct Entry {
                        uint64_t key[2 * W];
                        uint32_t idx;
                };
                const size_t n = store_.size();
                std::vector<Entry> entries(n);
                for (size_t i = 0; i < n; ++i) {
                        Entry& e = entries[i];
                        const PauliString<Coeff>& ps = store_[i];
                        const size_t nw = ps.x.size();
                        for (size_t k = 0; k < W; ++k) {
                                e.key[k]     = (k < nw) ? ps.x[k] : 0;
                                e.key[W + k] = (k < nw) ? ps.y[k] : 0;
                        }
                        e.idx = static_cast<uint32_t>(i);
                }

                auto key_less = [](const Entry& a, const Entry& b) {
                        for (size_t k = 0; k < 2 * W; ++k) {
                                if (a.key[k] != b.key[k]) return a.key[k] < b.key[k];
                        }
                        return false;
                };
                auto key_equal = [](const Entry& a, const Entry& b) {
                        for (size_t k = 0; k < 2 * W; ++k) {
                                if (a.key[k] != b.key[k]) return false;
                        }
                        return true;
                };

                parallel_sort_vec(entries, key_less, pauliengine_detail::use_threads(n));

                std::vector<PauliString<Coeff>> out;
                out.reserve(n);
                size_t i = 0;
                while (i < n) {
                        Coeff acc = std::move(store_[entries[i].idx].coeff);
                        size_t j = i + 1;
                        while (j < n && key_equal(entries[i], entries[j])) {
                                acc = acc + store_[entries[j].idx].coeff;
                                ++j;
                        }
                        if (!coeff_is_zero(acc)) {
                                PauliString<Coeff> ps = std::move(store_[entries[i].idx]);
                                ps.coeff = std::move(acc);
                                out.push_back(std::move(ps));
                        }
                        i = j;
                }
                store_ = std::move(out);
        }

        /// Chunked parallel merge sort if `parallel`, else std::sort. Elements must be
        /// safe to move across threads.
        template<typename T, typename Cmp>
        static void parallel_sort_vec(std::vector<T>& v, Cmp cmp, bool parallel) {
#ifdef PAULIENGINE_HAS_OPENMP
                if (parallel) {
                        const size_t n = v.size();
                        const int nthreads = omp_get_max_threads();
                        const size_t chunk = (n + nthreads - 1) / nthreads;

                        #pragma omp parallel for schedule(static)
                        for (int t = 0; t < nthreads; ++t) {
                                const size_t begin = static_cast<size_t>(t) * chunk;
                                if (begin >= n) continue;
                                const size_t end = std::min(begin + chunk, n);
                                std::sort(v.begin() + begin, v.begin() + end, cmp);
                        }

                        for (int step = 1; step < nthreads; step *= 2) {
                                #pragma omp parallel for schedule(dynamic)
                                for (int t = 0; t < nthreads; t += 2 * step) {
                                        const size_t begin = static_cast<size_t>(t) * chunk;
                                        const size_t mid   = std::min(static_cast<size_t>(t + step) * chunk, n);
                                        const size_t end   = std::min(static_cast<size_t>(t + 2 * step) * chunk, n);
                                        if (mid >= end || begin >= mid) continue;
                                        std::inplace_merge(v.begin() + begin,
                                                           v.begin() + mid,
                                                           v.begin() + end,
                                                           cmp);
                                }
                        }
                        return;
                }
#else
                (void)parallel;
#endif
                std::sort(v.begin(), v.end(), cmp);
        }

        static Coeff make_one() {
                if constexpr (kFlat) {
                        return std::complex<double>(1.0, 0.0);
                } else {
                        return Coeff(1);
                }
        }

        static Coeff make_zero() {
                if constexpr (kFlat) {
                        return std::complex<double>(0.0, 0.0);
                } else {
                        return Coeff(0);
                }
        }

        // Upper-cases a single Pauli label and rejects anything but I, X, Y, Z.
        static char normalize_pauli(char c) {
                const char upper = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (upper != 'I' && upper != 'X' && upper != 'Y' && upper != 'Z') {
                        throw std::invalid_argument(std::string("contains: invalid Pauli operator '") + c + "'");
                }
                return upper;
        }

    static Matrix2D get_pauli_matrix(const std::string& p){static const std::complex<double> I(0,1);
        if (p == "I") return {{1,0},{0,1}};
        else if (p == "X") return {{0,1},{1,0}};
        else if (p == "Y") return {{0,-I},{I,0}};
        else if (p == "Z") return {{1,0},{0,-1}};
        else throw std::runtime_error("Unknown Pauli operator");
}
    static Matrix2D tensor_product(const Matrix2D& A, const Matrix2D& B){
        int rows = A.size() * B.size();
        int cols = A[0].size() * B[0].size();
        Matrix2D result(rows, std::vector<std::complex<double>>(cols, 0));

        for (int i=0; i<(int)A.size(); ++i) {
                for (int j=0; j<(int)A[0].size(); ++j) {
                        for (int k=0; k<(int)B.size(); ++k) {
                                for (int l=0; l<(int)B[0].size(); ++l) {
                                        result[i*B.size()+k][j*B[0].size()+l] = A[i][j] * B[k][l];
                                }
                        }
                }
        }
        return result;
}
    static Matrix2D pauli_string_to_matrix(const std::unordered_map<int, std::string>& pauli_map, int num_qubits){
        Matrix2D result = get_pauli_matrix("I");
        for (int q=0; q < num_qubits; ++q) {
                std::string op = "I";
                auto it = pauli_map.find(q);
                if (it != pauli_map.end()) op = it->second;
                Matrix2D pmat = get_pauli_matrix(op);
                result = tensor_product(result, pmat);
        }
        return result;
}

};
