#pragma once
#include <array>
#include <unordered_map>
#include <string>
#include <complex>
#include <vector>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <algorithm>
#include <cctype>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <new>
#include <numeric>
#include <optional>
#include <sstream>
#include <utility>

#include "pauliengine/detail/Coefficients.h"
#include "pauliengine/PauliString.h"
#include "pauliengine/detail/PrefixPartitionedProduct.h"
#include "pauliengine/detail/OverlappingPairs.h"


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

/// A sum of Pauli strings. Kept compacted: no two terms share an operator, and
/// zero terms are dropped (only scaling by zero leaves them in place).
///
/// With complex coefficients the terms are stored flat (pauliengine_detail::FlatTerms),
/// and products, commutators and sums work on that storage directly; terms() builds
/// Pauli strings from it on request. Symbolic (polynomial) coefficients have a
/// variable size and are kept as a vector of strings; their products, commutators
/// and compaction run in parallel as well.
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

        /// Product this * other. For h * h (the same object on both sides) the
        /// square formula below needs only about a quarter of the pairs.
        QubitHamiltonian operator*(const QubitHamiltonian& other) const {
                if (this == &other) {
                        return square();
                }
                if constexpr (kFlat) {
                        return QubitHamiltonian(product_flat(store_, other.store_), from_storage_t{});
                } else {
                        Storage packed;
                        if (packed_pairs(store_, other.store_, false, packed)) {
                                return QubitHamiltonian(std::move(packed), from_storage_t{});
                        }
                        if (auto terms = symbolic_pairs(store_, other.store_, PairKind::Product)) {
                                return QubitHamiltonian(std::move(*terms), already_compacted);
                        }
                        // Fallback (monomials do not fit 64 bits, or very wide strings):
                        // polynomial coefficients are plain values, so every row of the
                        // product table can be formed on its own thread.
                        const size_t na = store_.size();
                        const size_t nb = other.store_.size();
                        std::vector<PauliString<Coeff>> product(na * nb);
                        pauliengine_detail::for_each_index(na, pauliengine_detail::use_threads(na * nb), [&](size_t i) {
                                const PauliString<Coeff>& a = store_[i];
                                for (size_t j = 0; j < nb; ++j) {
                                        product[i * nb + j] = a * other.store_[j];
                                }
                        });
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
        /// target. Cost O(min(|A|,|B|) L) expected word operations for the lookups
        /// plus O(max(|A|,|B|) L) for building the index, where L is the number of
        /// words per string. Works on the stored terms directly: the candidate
        /// operator is a plain XOR of words, and one phase pass per found partner
        /// gives both the commutation parity and i^k, without forming the product.
        Coeff targeted_commutator(const QubitHamiltonian& other,
                                  const PauliString<Coeff>& target) const {
                Coeff total(0.0);

                // Iterate the smaller side; index the larger one by operator.
                const bool iter_this = size() <= other.size();
                const QubitHamiltonian& outer = iter_this ? *this : other;
                const QubitHamiltonian& inner = iter_this ? other : *this;
                const size_t n_inner = inner.size();
                if (n_inner == 0 || outer.size() == 0) return total;

                // Hash index over the operators of `inner`: buckets keyed by the hash
                // of the trimmed words, chained through `next`.
                constexpr uint32_t kNone = static_cast<uint32_t>(-1);
                std::unordered_map<size_t, uint32_t> head;
                head.reserve(n_inner);
                std::vector<uint32_t> next(n_inner, kNone);
                for (size_t j = 0; j < n_inner; ++j) {
                        const OperatorWords q = inner.operator_words(j);
                        const auto [it, inserted] = head.try_emplace(hash_words(q.x, q.y, q.nw), static_cast<uint32_t>(j));
                        if (!inserted) {
                                next[j] = it->second;
                                it->second = static_cast<uint32_t>(j);
                        }
                }

                const size_t nr = target.x.size();
                std::vector<uint64_t> cand;
                std::vector<Coeff> parts;  // symbolic contributions, summed at the end
                for (size_t i = 0; i < outer.size(); ++i) {
                        const OperatorWords p = outer.operator_words(i);
                        // Operator of the unique partner that can produce `target`:
                        // v(Q) = v(P) XOR v(R), trimmed of trailing zero words.
                        size_t nc = std::max(p.nw, nr);
                        cand.assign(2 * nc, 0);
                        for (size_t k = 0; k < nc; ++k) {
                                cand[k]      = (k < p.nw ? p.x[k] : 0) ^ (k < nr ? target.x[k] : 0);
                                cand[nc + k] = (k < p.nw ? p.y[k] : 0) ^ (k < nr ? target.y[k] : 0);
                        }
                        const uint64_t* cx = cand.data();
                        const uint64_t* cy = cand.data() + nc;
                        while (nc > 0 && cx[nc - 1] == 0 && cy[nc - 1] == 0) --nc;

                        const auto bucket = head.find(hash_words(cx, cy, nc));
                        if (bucket == head.end()) continue;
                        for (uint32_t j = bucket->second; j != kNone; j = next[j]) {
                                const OperatorWords q = inner.operator_words(j);
                                if (q.nw != nc || !std::equal(cx, cx + nc, q.x) || !std::equal(cy, cy + nc, q.y)) {
                                        continue;
                                }
                                // Preserve the [this, other] ordering (the commutator is
                                // antisymmetric): the term from `this` is the left factor.
                                const OperatorWords& left = iter_this ? p : q;
                                const OperatorWords& right = iter_this ? q : p;
                                const int power_of_i = symplectic_product_phase(
                                        left.x, left.y, left.nw, right.x, right.y, right.nw);
                                if (power_of_i & 1) {  // anti-commuting pair
                                        const size_t li = iter_this ? i : j;
                                        const size_t ri = iter_this ? j : i;
                                        Coeff part = PauliString<Coeff>::commutator_coeff(
                                                term_coeff(li), other.term_coeff(ri), power_of_i);
                                        if constexpr (kFlat) {
                                                total = total + part;
                                        } else {
                                                parts.push_back(std::move(part));
                                        }
                                }
                                break;  // operators are unique within a compacted sum
                        }
                }
                if constexpr (!kFlat) {
                        // One merge of all contributions; adding them one by one is quadratic.
                        total = Coeff::sum_of(parts.size(), [&](size_t k) -> const Coeff& { return parts[k]; });
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
                if constexpr (kFlat) {
                        return map_terms([&](const PauliString<Coeff>& ps) { return ps.diff(symbol); });
                } else {
                        // Resolve the name once, not once per term under the table lock.
                        typename Coeff::DiffTarget target;
                        if (!Coeff::prepare_diff(symbol, target)) return zero();
                        return map_terms([&](const PauliString<Coeff>& ps) {
                                return PauliString<Coeff>(ps.x, ps.y, ps.coeff.diff(target));
                        });
                }
        }

        /// Values ({name: complex}) or polynomials ({name: Expression}) substituted for
        /// symbols of every coefficient.
        template <typename Values>
        QubitHamiltonian substitute(const Values& substitution_map) const {
                if constexpr (kFlat) {
                        return map_terms([&](const PauliString<Coeff>& ps) { return ps.substitute(substitution_map); });
                } else {
                        const auto values = Coeff::prepare_subs(substitution_map);
                        return map_terms([&](const PauliString<Coeff>& ps) {
                                return PauliString<Coeff>(ps.x, ps.y, ps.coeff.subs(values));
                        });
                }
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
                                for (size_t j = 0; j < BITS_IN_INTEGER; ++j) {
                                        uint64_t bit = (uint64_t)1 << j;
                                        bool x_bit = x_word & bit;
                                        bool y_bit = y_word & bit;
                                        if (!x_bit && !y_bit) continue;
                                        const int index = static_cast<int>(i * BITS_IN_INTEGER + j);
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
                        Storage packed;
                        if (packed_pairs(store_, other.store_, true, packed)) {
                                return QubitHamiltonian(std::move(packed), from_storage_t{});
                        }
                        if (auto terms = symbolic_pairs(store_, other.store_, PairKind::Commutator)) {
                                return QubitHamiltonian(std::move(*terms), already_compacted);
                        }
                        // Fallback: one list per row of the pair table, filled in parallel and
                        // joined in row order, so the result does not depend on the thread count.
                        const size_t na = store_.size();
                        const size_t nb = other.store_.size();
                        std::vector<std::vector<PauliString<Coeff>>> rows(na);
                        pauliengine_detail::for_each_index(na, pauliengine_detail::use_threads(na * nb), [&](size_t i) {
                                const PauliString<Coeff>& a = store_[i];
                                for (const auto& b : other.store_) {
                                        PauliString<Coeff> t = a.commutator(b);
                                        if (!t.x.empty()) {  // commuting pairs give the empty string
                                                rows[i].push_back(std::move(t));
                                        }
                                }
                        });
                        size_t total = 0;
                        for (const auto& r : rows) total += r.size();
                        std::vector<PauliString<Coeff>> terms;
                        terms.reserve(total);
                        for (auto& r : rows) {
                                std::move(r.begin(), r.end(), std::back_inserter(terms));
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
                const ConjMap conj = conjugation_map();
                return map_terms([&](const PauliString<Coeff>& ps) {
                        return PauliString<Coeff>(ps.x, ps.y, conjugate_coeff(ps.coeff, conj));
                });
        }

        QubitHamiltonian conjugate() const {
                const ConjMap conj = conjugation_map();
                return map_terms([&](const PauliString<Coeff>& ps) {
                        Coeff new_coeff = conjugate_coeff(ps.coeff, conj);
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
                if (n == 1) return *this;
                QubitHamiltonian result = square();
                for (int i = 2; i < n; ++i) {
                        result = result * *this;
                }
                return result;
        }

        /// this * this, from the pairs i <= j only. P_i P_j and P_j P_i are equal if
        /// the strings commute and cancel if they anti-commute, and P_i P_i = I, so
        ///
        ///     H^2 = (sum_i a_i^2) I + 2 sum_{i < j, P_i P_j = P_j P_i} a_i a_j P_i P_j .
        ///
        /// About a quarter of the pairs of the general product remain, and the
        /// anti-commuting pairs leave no rounding residue.
        QubitHamiltonian square() const {
                if constexpr (kFlat) {
                        return QubitHamiltonian(square_flat(store_), from_storage_t{});
                } else {
                        if (auto terms = symbolic_pairs(store_, store_, PairKind::Square)) {
                                return QubitHamiltonian(std::move(*terms), already_compacted);
                        }
                        const size_t n = store_.size();
                        Coeff diagonal = make_zero();
                        for (const auto& ps : store_) {
                                diagonal = diagonal + ps.coeff * ps.coeff;
                        }
                        std::vector<PauliString<Coeff>> terms;
                        terms.emplace_back(std::vector<uint64_t>{}, std::vector<uint64_t>{}, diagonal);
                        for (size_t i = 0; i < n; ++i) {
                                const PauliString<Coeff>& a = store_[i];
                                for (size_t j = i + 1; j < n; ++j) {
                                        const PauliString<Coeff>& b = store_[j];
                                        if (pauliengine_detail::anticommute_words(a.x.data(), a.y.data(), a.x.size(),
                                                                                  b.x.data(), b.y.data(), b.x.size())) {
                                                continue;
                                        }
                                        PauliString<Coeff> p = a * b;
                                        p.coeff = p.coeff * Coeff(2);
                                        terms.push_back(std::move(p));
                                }
                        }
                        return QubitHamiltonian(std::move(terms));
                }
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
        /// Operator part of a stored term, trimmed of trailing zero words.
        struct OperatorWords {
                const uint64_t* x;
                const uint64_t* y;
                size_t nw;
        };

        OperatorWords operator_words(size_t i) const {
                if constexpr (kFlat) {
                        return {store_.x(i), store_.y(i), store_.used_words(i)};
                } else {
                        const PauliString<Coeff>& ps = store_[i];
                        return {ps.x.data(), ps.y.data(), ps.x.size()};
                }
        }

        Coeff term_coeff(size_t i) const {
                if constexpr (kFlat) {
                        return store_.coeff(i);
                } else {
                        return store_[i].coeff;
                }
        }

        /// Hash of an operator given by its trimmed words.
        static size_t hash_words(const uint64_t* x, const uint64_t* y, size_t nw) {
                size_t h = 0xcbf29ce484222325ULL ^ nw;
                const auto mix = [&h](uint64_t w) {
                        h ^= std::hash<uint64_t>{}(w) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                };
                for (size_t k = 0; k < nw; ++k) mix(x[k]);
                for (size_t k = 0; k < nw; ++k) mix(y[k]);
                return h;
        }

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
                std::vector<PauliString<Coeff>> mapped(ts.size());
                // Symbolic terms are expensive enough to share out even in small numbers.
                const size_t work = kFlat ? ts.size() : ts.size() * 64;
                pauliengine_detail::for_each_index(ts.size(), pauliengine_detail::use_threads(work),
                        [&](size_t i) { mapped[i] = f(ts[i]); });
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
                        // Stored terms are compacted apart from zeros of an earlier scaling by
                        // zero; dropping those is all compact() would do.
                        Storage copy;
                        copy.reserve(store_.size());
                        for (const auto& ps : store_) {
                                if (!ps.coeff.is_zero()) copy.push_back(ps);
                        }
                        QubitHamiltonian out(std::move(copy), from_storage_t{});
                        pauliengine_detail::for_each_index(out.store_.size(), pauliengine_detail::use_threads(out.store_.size() * 16),
                                [&](size_t i) { out.store_[i].coeff *= scale; });
                        return out;
                }
        }

        /// Symbol conjugation map, fetched once per operation (empty: all symbols real).
        using ConjMap = std::vector<uint32_t>;
        static ConjMap conjugation_map() {
                if constexpr (kFlat) {
                        return {};
                } else {
                        return Coeff::conjugation_map();
                }
        }

        static Coeff conjugate_coeff(const Coeff& c, const ConjMap& conj) {
                if constexpr (kFlat) {
                        (void)conj;
                        return std::conj(c);
                } else {
                        return c.conjugate(conj);
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
                        return c.is_zero();
                }
        }

        /// Whether a sum of several coefficients is zero (polynomials are canonical,
        /// so this is the same test as for a single coefficient).
        static bool sum_is_zero(const Coeff& c) { return coeff_is_zero(c); }

        static bool coeff_equal(const Coeff& a, const Coeff& b) {
                return a == b;  // complex numbers, or canonical polynomials
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

        /// Up to this many pairs all records are formed at once and sorted; setting up
        /// the blocks of PrefixPartitionedProduct would take longer than that.
        static constexpr size_t kDirectPairs = 256;

        /// a * b, or a * a from the pairs i <= j that commute if `square` (b is a then),
        /// for few pairs of strings of at most W words.
        template<size_t W>
        static FlatTerms direct_small_key(const FlatTerms& a, const FlatTerms& b, bool square, size_t nw) {
                using Term = KeyedTerm<W>;
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                std::vector<Term> records;
                records.reserve(square ? n1 * (n1 + 1) / 2 : n1 * n2);
                for (size_t i = 0; i < n1; ++i) {
                        const TermView left = a.view(i);
                        if (square) set_keyed_product<W>(records.emplace_back(), left, left, 1.0);
                        for (size_t j = square ? i + 1 : 0; j < n2; ++j) {
                                const TermView right = b.view(j);
                                if (square && pauliengine_detail::anticommute(left, right)) continue;
                                set_keyed_product<W>(records.emplace_back(), left, right, square ? 2.0 : 1.0);
                        }
                }
                return merge_keyed_terms<W>(records.data(), records.size(), nw);
        }

        static FlatTerms product_flat(const FlatTerms& a, const FlatTerms& b) {
                if (a.size() * b.size() <= kDirectPairs) {
                        FlatTerms out;
                        const size_t nw = std::max(a.width(), b.width());
                        if (with_small_key(nw, [&]<size_t W>() { out = direct_small_key<W>(a, b, false, nw); })) return out;
                }
                return pauliengine_detail::PrefixPartitionedProduct::run(a, b, false);
        }

        /// a * a (see square()).
        static FlatTerms square_flat(const FlatTerms& a) {
                if (a.size() * (a.size() + 1) / 2 <= kDirectPairs) {
                        FlatTerms out;
                        if (with_small_key(a.width(), [&]<size_t W>() { out = direct_small_key<W>(a, a, true, a.width()); })) return out;
                }
                return pauliengine_detail::PrefixPartitionedProduct::square(a);
        }

        /// [a, b] for strings of at most W words, written with width nw. With
        /// `overlap`, only the pairs listed there are tested.
        template<size_t W>
        static FlatTerms commutator_small_key(const FlatTerms& a, const FlatTerms& b, size_t nw,
                                              const pauliengine_detail::OverlappingPairs* overlap) {
                using Term = KeyedTerm<W>;
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                const size_t pairs = overlap ? overlap->size() : n1 * n2;
                // Appends the commutators of row i with all of b that do not vanish.
                // Most pairs are only tested here, so the test reads nothing but the words.
                const size_t wa = a.width();
                const size_t wb = b.width();
                // The same for the listed pairs only (kept apart from add_row, whose
                // loops over all of b are the hot path otherwise).
                const auto add_listed = [&](size_t i, std::vector<Term>& mine) {
                        const TermView left = a.view(i);
                        const size_t w = std::min(wa, wb);
                        const uint32_t* col = overlap->cols.data();
                        for (size_t k = overlap->start[i]; k < overlap->start[i + 1]; ++k) {
                                const size_t j = col[k];
                                const uint64_t* row = b.x(j);  // x words, then y words
                                uint64_t anti = 0;
                                for (size_t t = 0; t < w; ++t) anti ^= (left.x[t] & row[wb + t]) ^ (row[t] & left.y[t]);
                                if (std::popcount(anti) & 1) {
                                        set_keyed_product<W>(mine.emplace_back(), left, b.view(j), 2.0);
                                }
                        }
                };
                const auto add_row = [&](size_t i, std::vector<Term>& mine) {
                        const TermView left = a.view(i);
                        if (wa == 1 && wb == 1) {
                                const uint64_t ax = left.x[0];
                                const uint64_t ay = left.y[0];
                                const uint64_t* row = b.x(0);  // x word, then y word
                                for (size_t j = 0; j < n2; ++j, row += 2) {
                                        if (std::popcount((ax & row[1]) ^ (row[0] & ay)) & 1) {
                                                set_keyed_product<W>(mine.emplace_back(), left, b.view(j), 2.0);
                                        }
                                }
                                return;
                        }
                        const size_t w = std::min(wa, wb);
                        const uint64_t* row = b.x(0);  // x words, then y words
                        for (size_t j = 0; j < n2; ++j, row += 2 * wb) {
                                uint64_t anti = 0;
                                for (size_t k = 0; k < w; ++k) anti ^= (left.x[k] & row[wb + k]) ^ (row[k] & left.y[k]);
                                if (std::popcount(anti) & 1) {
                                        set_keyed_product<W>(mine.emplace_back(), left, b.view(j), 2.0);
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
                                for (std::ptrdiff_t i = 0; i < rows; ++i) {
                                        if (overlap) add_listed(static_cast<size_t>(i), mine); else add_row(static_cast<size_t>(i), mine);
                                }
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
                for (size_t i = 0; i < n1; ++i) {
                        if (overlap) add_listed(i, mine); else add_row(i, mine);
                }
                return merge_keyed_terms<W>(mine.data(), mine.size(), nw);
        }

        /// [a, b] = ab - ba: twice the product of every anti-commuting pair.
        static FlatTerms commutator_flat(const FlatTerms& a, const FlatTerms& b) {
                // Below a few thousand pairs the two ways differ by less than the sample costs.
                using pauliengine_detail::PrefixPartitionedProduct;
                if (a.size() * b.size() >= 4096 && PrefixPartitionedProduct::many_anticommute(a, b)) {
                        return PrefixPartitionedProduct::run(a, b, true);
                }
                // Few pairs contribute: pick them out first, then sort and merge them.
                // If few pairs even share a qubit, only those are tested.
                std::optional<pauliengine_detail::OverlappingPairs> overlap;
                if (a.size() * b.size() >= 4096) {
                        overlap = pauliengine_detail::OverlappingPairs::find(a, b);
                }
                const pauliengine_detail::OverlappingPairs* listed = overlap ? &*overlap : nullptr;
                using pauliengine_detail::anticommute;
                using pauliengine_detail::for_each_index;
                using pauliengine_detail::use_threads;
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                const size_t pairs = listed ? listed->size() : n1 * n2;
                const size_t nw = std::max(a.width(), b.width());
                FlatTerms out;
                if (with_small_key(nw, [&]<size_t W>() { out = commutator_small_key<W>(a, b, nw, listed); })) {
                        return out;
                }
                // Long strings, listed pairs: row i writes its results to the start of its
                // own block of listed pairs; the rest of the block is never touched.
                if (listed) {
                        FlatTerms records(pairs, nw);
                        std::vector<size_t> found(n1, 0);
                        for_each_index(n1, use_threads(pairs * nw), [&](size_t i) {
                                const TermView left = a.view(i);
                                size_t r = listed->start[i];
                                for (size_t k = listed->start[i]; k < listed->start[i + 1]; ++k) {
                                        const TermView right = b.view(listed->cols[k]);
                                        if (!anticommute(left, right)) continue;
                                        const int power_of_i = records.set_row_product(r, left, right);
                                        records.set_coeff(r, pauliengine_detail::product_coeff(left.coeff, right.coeff, power_of_i, 2.0));
                                        ++r;
                                }
                                found[i] = r - listed->start[i];
                        });
                        std::vector<size_t> order;
                        order.reserve(std::accumulate(found.begin(), found.end(), size_t{0}));
                        for (size_t i = 0; i < n1; ++i) {
                                for (size_t r = listed->start[i]; r < listed->start[i] + found[i]; ++r) order.push_back(r);
                        }
                        return merge_long(records, order);
                }
                // Long strings, all pairs: a block per row would take memory for all
                // n1 * n2 pairs, so each thread appends its results to a growing buffer
                // instead, and memory follows the pairs found.
                struct Found {
                        FlatTerms terms;  // the first `used` rows are written
                        size_t used = 0;
                };
                const auto scan_row = [&](size_t i, Found& mine) {
                        const TermView left = a.view(i);
                        for (size_t j = 0; j < n2; ++j) {
                                const TermView right = b.view(j);
                                if (!anticommute(left, right)) continue;
                                if (mine.used == mine.terms.size()) {
                                        mine.terms.resize(std::max<size_t>(64, 2 * mine.terms.size()));
                                }
                                const size_t r = mine.used++;
                                const int power_of_i = mine.terms.set_row_product(r, left, right);
                                mine.terms.set_coeff(r, pauliengine_detail::product_coeff(left.coeff, right.coeff, power_of_i, 2.0));
                        }
                };
                std::vector<Found> local(1);
                local[0].terms = FlatTerms(0, nw);
#ifdef PAULIENGINE_HAS_OPENMP
                if (use_threads(pairs * nw)) {
                        local.resize(static_cast<size_t>(omp_get_max_threads()));
                        const std::ptrdiff_t rows = static_cast<std::ptrdiff_t>(n1);  // signed for OpenMP 2.0
                        #pragma omp parallel
                        {
                                // Filled as a local object: the buffers in `local` share cache lines.
                                Found mine;
                                mine.terms = FlatTerms(0, nw);
                                #pragma omp for schedule(static)
                                for (std::ptrdiff_t i = 0; i < rows; ++i) scan_row(static_cast<size_t>(i), mine);
                                local[static_cast<size_t>(omp_get_thread_num())] = std::move(mine);
                        }
                } else
#endif
                {
                        for (size_t i = 0; i < n1; ++i) scan_row(i, local[0]);
                }
                FlatTerms records;
                if (local.size() == 1) {
                        records = std::move(local[0].terms);
                        records.shrink_to(local[0].used);
                } else {
                        // A static schedule hands out ascending rows, so this is the serial order.
                        std::vector<size_t> offset(local.size() + 1, 0);
                        for (size_t t = 0; t < local.size(); ++t) offset[t + 1] = offset[t] + local[t].used;
                        records = FlatTerms(offset.back(), nw);
                        for_each_index(local.size(), true, [&](size_t t) {
                                if (local[t].used) {
                                        records.assign_rows(offset[t], local[t].terms.x(0), local[t].terms.coeff_data(0), local[t].used);
                                }
                                local[t].terms = FlatTerms();
                        });
                }
                std::vector<size_t> order(records.size());
                std::iota(order.begin(), order.end(), size_t{0});
                return merge_long(records, order);
        }

        // ------------------------------------------------------------------
        // Symbolic coefficients
        // ------------------------------------------------------------------

        enum class PairKind { Product, Commutator, Square };

        /// a * b, [a, b] or a * a (b is a then) with symbolic coefficients, or nothing
        /// if the strings are too wide for fixed-size keys.
        static std::optional<std::vector<PauliString<Coeff>>> symbolic_pairs(const Storage& a, const Storage& b, PairKind kind) {
                std::optional<std::vector<PauliString<Coeff>>> out;
                with_small_key(std::max(widest(a), widest(b)), [&]<size_t W>() { out = symbolic_pairs_small_key<W>(a, b, kind); });
                return out;
        }

        /// The pairs are first grouped by the operator they give, with flat keys
        /// only, as in the products of complex sums. The phase i^k and the factor 2 of
        /// a pair are taken into the coefficient of its left factor, which has only
        /// eight such variants: every pair then costs one polynomial product, and every
        /// operator one merged sum, instead of up to three operations per pair and a
        /// chain of additions. Used for squares, and for products and commutators
        /// whose monomials do not fit packed_pairs. Phase and factor are exact (sign changes, swapped parts, a
        /// power of two), so the order in which they are applied does not matter.
        ///
        /// The terms come out ordered by operator and compacted, the pairs of an
        /// operator added in the order (i, j): the result does not depend on the
        /// number of threads.
        template<size_t W>
        static std::vector<PauliString<Coeff>> symbolic_pairs_small_key(const Storage& a, const Storage& b, PairKind kind) {
                struct Pair {
                        uint64_t key[2 * W];
                        uint32_t i, j;
                        uint32_t variant;  // phase exponent k, plus 4 if the pair counts twice
                };
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                if (n1 > UINT32_MAX || n2 > UINT32_MAX) {
                        throw std::length_error("QubitHamiltonian: more than 2^32 - 1 terms in a product");
                }
                const bool square = kind == PairKind::Square;
                std::vector<Pair> pairs;
                pairs.reserve(square ? n1 * (n1 + 1) / 2 : (kind == PairKind::Commutator ? n1 * n2 / 2 : n1 * n2));
                for (size_t i = 0; i < n1; ++i) {
                        const PauliString<Coeff>& left = a[i];
                        for (size_t j = square ? i : 0; j < n2; ++j) {
                                const PauliString<Coeff>& right = b[j];
                                const int k = symplectic_product_phase(left.x.data(), left.y.data(), left.x.size(),
                                                                       right.x.data(), right.y.data(), right.x.size()) & 3;
                                // odd k: the strings anti-commute
                                if (kind == PairKind::Commutator && !(k & 1)) continue;
                                if (square && (k & 1)) continue;
                                Pair& p = pairs.emplace_back();
                                for (size_t w = 0; w < W; ++w) {
                                        p.key[w] = (w < left.x.size() ? left.x[w] : 0) ^ (w < right.x.size() ? right.x[w] : 0);
                                        p.key[W + w] = (w < left.y.size() ? left.y[w] : 0) ^ (w < right.y.size() ? right.y[w] : 0);
                                }
                                p.i = static_cast<uint32_t>(i);
                                p.j = static_cast<uint32_t>(j);
                                const bool twice = kind == PairKind::Commutator || (square && i != j);
                                p.variant = static_cast<uint32_t>(k) | (twice ? 4u : 0u);
                        }
                }
                const auto pair_less = [](const Pair& x, const Pair& y) {
                        for (size_t w = 0; w < 2 * W; ++w) {
                                if (x.key[w] != y.key[w]) return x.key[w] < y.key[w];
                        }
                        return x.i != y.i ? x.i < y.i : x.j < y.j;
                };
                parallel_sort_vec(pairs, pair_less, pauliengine_detail::use_threads(pairs.size() * W));

                // coefficient of a[i] times phase and factor, formed when first needed
                std::vector<Coeff> variants(8 * n1);
                std::vector<uint8_t> formed(8 * n1, 0);
                const auto left_coeff = [&](uint32_t i, uint32_t variant) -> const Coeff& {
                        const size_t slot = 8 * static_cast<size_t>(i) + variant;
                        if (!formed[slot]) {
                                Coeff c = a[i].coeff;
                                if (variant & 1) c *= std::complex<double>(0.0, 1.0);
                                if (variant & 2) c = -c;
                                if (variant & 4) c = square ? c * Coeff(2) : c * 2.0;
                                variants[slot] = std::move(c);
                                formed[slot] = 1;
                        }
                        return variants[slot];
                };
                const auto same_key = [](const Pair& x, const Pair& y) {
                        for (size_t w = 0; w < 2 * W; ++w) {
                                if (x.key[w] != y.key[w]) return false;
                        }
                        return true;
                };
                size_t operators = 0;
                for (size_t r = 0; r < pairs.size(); ++r) operators += r == 0 || !same_key(pairs[r - 1], pairs[r]);
                std::vector<PauliString<Coeff>> terms;
                terms.reserve(operators);
                std::vector<Coeff> parts;
                WordVec x(W), y(W);
                for (size_t r = 0; r < pairs.size();) {
                        size_t next = r + 1;
                        while (next < pairs.size() && same_key(pairs[r], pairs[next])) ++next;
                        Coeff coeff;
                        if (next - r == 1) {
                                coeff = left_coeff(pairs[r].i, pairs[r].variant) * b[pairs[r].j].coeff;
                        } else {
                                parts.clear();
                                for (size_t q = r; q < next; ++q) {
                                        parts.push_back(left_coeff(pairs[q].i, pairs[q].variant) * b[pairs[q].j].coeff);
                                }
                                coeff = Coeff::sum_of(parts.size(), [&](size_t k) -> const Coeff& { return parts[k]; });
                        }
                        if (!(next - r == 1 ? coeff_is_zero(coeff) : sum_is_zero(coeff))) {
                                std::copy(pairs[r].key, pairs[r].key + W, x.data());
                                std::copy(pairs[r].key + W, pairs[r].key + 2 * W, y.data());
                                terms.emplace_back(x, y, std::move(coeff));
                        }
                        r = next;
                }
                return terms;
        }

        /// a * b, or [a, b] if `commutator`, for polynomial coefficients, on flat
        /// records. Every term of a coefficient becomes one record (operator words,
        /// packed monomial) -> complex number; the product of two records multiplies
        /// the operators, adds the monomial words and multiplies the numbers. The
        /// records then go through the same parallel sort and merge as complex
        /// Hamiltonians (the monomial sits in the spare last key word, so equal
        /// operators end up adjacent with ascending monomials), and each group of
        /// equal operators is turned back into one string with a polynomial.
        /// Returns false, leaving `result` alone, if the monomials of the product do
        /// not fit 64 bits or the strings are too wide for fixed-size records.
        static bool packed_pairs(const std::vector<PauliString<Coeff>>& a, const std::vector<PauliString<Coeff>>& b,
                                 bool commutator, std::vector<PauliString<Coeff>>& result) {
                using pauliengine::symbolic::MonomialPacking;
                using pauliengine_detail::for_each_index;
                using pauliengine_detail::use_threads;
                const size_t na = a.size();
                const size_t nb = b.size();
                if (na * nb <= 4) return false;  // not worth the setup

                std::vector<uint32_t> max_a, max_b;
                for (const auto& ps : a) MonomialPacking::note_exponents(ps.coeff, max_a);
                for (const auto& ps : b) MonomialPacking::note_exponents(ps.coeff, max_b);
                MonomialPacking packing;
                if (!packing.init(max_a, max_b)) return false;

                // The terms of all coefficients of one operand: those of string i are
                // [start[i], start[i + 1]).
                struct Monomials {
                        std::vector<size_t> start;
                        std::vector<uint64_t> word;
                        std::vector<std::complex<double>> coeff;
                };
                const auto flatten = [&packing](const std::vector<PauliString<Coeff>>& strings) {
                        Monomials m;
                        m.start.reserve(strings.size() + 1);
                        m.word.reserve(strings.size() * 2);
                        m.coeff.reserve(strings.size() * 2);
                        for (const auto& ps : strings) {
                                m.start.push_back(m.word.size());
                                ps.coeff.for_each_term([&](std::complex<double> c, const uint64_t* f, size_t nf) {
                                        m.word.push_back(packing.pack(f, nf));
                                        m.coeff.push_back(c);
                                });
                        }
                        m.start.push_back(m.word.size());
                        return m;
                };
                const Monomials ma = flatten(a);
                const Monomials mb = flatten(b);
                const size_t nw = std::max(widest(a), widest(b));
                const auto view = [](const PauliString<Coeff>& ps) {
                        return TermView{ps.x.data(), ps.y.data(), ps.x.size(), std::complex<double>(1.0, 0.0)};
                };

                // One word more than the strings need, for the monomial.
                return with_small_key(nw + 1, [&]<size_t W>() {
                        using Term = KeyedTerm<W>;
                        // Row i of `a` writes its records to [offset[i], offset[i + 1]).
                        std::vector<size_t> offset(na + 1, 0);
                        if (commutator) {
                                for_each_index(na, use_threads(na * nb), [&](size_t i) {
                                        const TermView left = view(a[i]);
                                        size_t partners = 0;
                                        for (size_t j = 0; j < nb; ++j) {
                                                if (pauliengine_detail::anticommute(left, view(b[j]))) {
                                                        partners += mb.start[j + 1] - mb.start[j];
                                                }
                                        }
                                        offset[i + 1] = (ma.start[i + 1] - ma.start[i]) * partners;
                                });
                        } else {
                                for (size_t i = 0; i < na; ++i) {
                                        offset[i + 1] = (ma.start[i + 1] - ma.start[i]) * mb.word.size();
                                }
                        }
                        for (size_t i = 0; i < na; ++i) offset[i + 1] += offset[i];
                        const size_t total = offset[na];

                        // new[] leaves the records untouched, so the threads fault their own pages in.
                        std::unique_ptr<Term[]> records(new Term[total]);
                        Term* const raw = records.get();
                        const double scale = commutator ? 2.0 : 1.0;
                        for_each_index(na, use_threads(total * W), [&](size_t i) {
                                const TermView left = view(a[i]);
                                Term* write = raw + offset[i];
                                Term pair;
                                for (size_t j = 0; j < nb; ++j) {
                                        const TermView right = view(b[j]);
                                        if (commutator && !pauliengine_detail::anticommute(left, right)) continue;
                                        set_keyed_product<W>(pair, left, right, scale);  // operator and scale * i^k
                                        const std::complex<double> phase(pair.re, pair.im);
                                        for (size_t s = ma.start[i]; s < ma.start[i + 1]; ++s) {
                                                const std::complex<double> cs = phase * ma.coeff[s];
                                                for (size_t t = mb.start[j]; t < mb.start[j + 1]; ++t) {
                                                        *write = pair;
                                                        write->key[2 * W - 1] = ma.word[s] + mb.word[t];
                                                        const std::complex<double> c = cs * mb.coeff[t];
                                                        write->re = c.real();
                                                        write->im = c.imag();
                                                        ++write;
                                                }
                                        }
                                }
                        });

                        // Sorted by operator, then monomial; equal ones added, zeros dropped.
                        const FlatTerms merged = merge_keyed_terms<W>(raw, total, W);
                        records.reset();
                        const size_t n = merged.size();
                        const auto same_operator = [&](size_t i, size_t j) {
                                return std::equal(merged.x(i), merged.x(i) + W, merged.x(j))
                                        && std::equal(merged.y(i), merged.y(i) + W - 1, merged.y(j));
                        };
                        std::vector<size_t> heads;
                        heads.reserve(n + 1);
                        for (size_t i = 0; i < n; ++i) {
                                if (i == 0 || !same_operator(i - 1, i)) heads.push_back(i);
                        }
                        heads.push_back(n);
                        const size_t groups = heads.size() - 1;

                        std::vector<PauliString<Coeff>> out(groups);
                        for_each_index(groups, use_threads(n * 16), [&](size_t g) {
                                const size_t begin = heads[g];
                                const uint64_t* const x = merged.x(begin);
                                const uint64_t* const y = merged.y(begin);
                                size_t used = nw;
                                while (used > 0 && x[used - 1] == 0 && y[used - 1] == 0) --used;
                                PauliString<Coeff>& ps = out[g];
                                ps.x.assign(x, x + used);
                                ps.y.assign(y, y + used);
                                Coeff c;
                                for (size_t r = begin; r < heads[g + 1]; ++r) {
                                        packing.append(merged.y(r)[W - 1], merged.coeff(r), c);
                                }
                                ps.coeff = std::move(c);
                        });
                        result = std::move(out);
                });
        }

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

                // Only for strings wider than kMaxSmallKeyWords words; serial, since the
                // strings themselves are sorted and have no index to break ties with.
                std::sort(store_.begin(), store_.end(), op_less);

                size_t write = 0;
                std::vector<char> merged(store_.size(), 0);
                for (size_t read = 1; read < store_.size(); ++read) {
                        if (op_equal(store_[write], store_[read])) {
                                store_[write].coeff = store_[write].coeff + store_[read].coeff;
                                merged[write] = 1;
                        } else {
                                ++write;
                                if (write != read) {
                                        store_[write] = std::move(store_[read]);
                                }
                        }
                }
                size_t kept = 0;
                for (size_t r = 0; r <= write; ++r) {
                        if (merged[r] ? sum_is_zero(store_[r].coeff) : coeff_is_zero(store_[r].coeff)) continue;
                        if (kept != r) store_[kept] = std::move(store_[r]);
                        ++kept;
                }
                store_.resize(kept);
        }

        template<size_t W>
        void compact_small_key_symbolic() {
                struct Entry {
                        uint64_t key[2 * W];
                        uint32_t idx;
                };
                using pauliengine_detail::for_each_index;
                using pauliengine_detail::use_threads;
                const size_t n = store_.size();
                std::vector<Entry> entries(n);
                for_each_index(n, use_threads(n * W), [&](size_t i) {
                        Entry& e = entries[i];
                        const PauliString<Coeff>& ps = store_[i];
                        const size_t nw = ps.x.size();
                        for (size_t k = 0; k < W; ++k) {
                                e.key[k]     = (k < nw) ? ps.x[k] : 0;
                                e.key[W + k] = (k < nw) ? ps.y[k] : 0;
                        }
                        e.idx = static_cast<uint32_t>(i);
                });

                // Ties in input order, so sums add up the same way for any thread count.
                auto key_less = [](const Entry& a, const Entry& b) {
                        for (size_t k = 0; k < 2 * W; ++k) {
                                if (a.key[k] != b.key[k]) return a.key[k] < b.key[k];
                        }
                        return a.idx < b.idx;
                };
                auto key_equal = [](const Entry& a, const Entry& b) {
                        for (size_t k = 0; k < 2 * W; ++k) {
                                if (a.key[k] != b.key[k]) return false;
                        }
                        return true;
                };

                parallel_sort_vec(entries, key_less, use_threads(n));

                // Groups of equal operators: [heads[g], heads[g + 1]).
                std::vector<size_t> heads;
                heads.reserve(n + 1);
                for (size_t i = 0; i < n; ++i) {
                        if (i == 0 || !key_equal(entries[i - 1], entries[i])) heads.push_back(i);
                }
                heads.push_back(n);
                const size_t n_groups = heads.size() - 1;

                std::vector<PauliString<Coeff>> out(n_groups);
                std::vector<char> keep(n_groups, 0);
                for_each_index(n_groups, use_threads(n * 16), [&](size_t g) {
                        const size_t begin = heads[g];
                        const size_t count = heads[g + 1] - begin;
                        Coeff acc = count == 1
                                ? std::move(store_[entries[begin].idx].coeff)
                                : Coeff::sum_of(count, [&](size_t t) -> const Coeff& { return store_[entries[begin + t].idx].coeff; });
                        if (!coeff_is_zero(acc)) {
                                PauliString<Coeff> ps = std::move(store_[entries[begin].idx]);
                                ps.coeff = std::move(acc);
                                out[g] = std::move(ps);
                                keep[g] = 1;
                        }
                });
                size_t write = 0;
                for (size_t g = 0; g < n_groups; ++g) {
                        if (!keep[g]) continue;
                        if (write != g) out[write] = std::move(out[g]);
                        ++write;
                }
                out.resize(write);
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
        const size_t rows = A.size() * B.size();
        const size_t cols = A[0].size() * B[0].size();
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
