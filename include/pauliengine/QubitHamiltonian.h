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
#include <memory>
#include <sstream>

#include "pauliengine/PauliString.h"

#include <symengine/expression.h>
#include <symengine/functions.h>
using SymEngine::Expression;

#ifdef PAULIENGINE_HAS_OPENMP
#include <omp.h>
#endif

// Threshold for switching on OpenMP. Below this many pair-multiplications the
// thread-spawn overhead dominates.
#ifndef PAULIENGINE_OMP_THRESHOLD
#define PAULIENGINE_OMP_THRESHOLD 4096
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

template<typename Coeff>
class QubitHamiltonian{
    public:
    std::vector<PauliString<Coeff>> data;

        QubitHamiltonian(const std::vector<PauliString<Coeff>>& data){
                this->data = data;
                this->compact_in_place();
        }


        QubitHamiltonian(std::vector<PauliString<Coeff>>&& data){
                this->data = std::move(data);
                this->compact_in_place();
        }

        /// Fast-path constructor: skips compact(). Only use when you can prove
        /// the vector has no duplicates and no zero-coefficient terms.
        QubitHamiltonian(std::vector<PauliString<Coeff>>&& data, already_compacted_t){
                this->data = std::move(data);
        }

        QubitHamiltonian(const Hamiltonian_structure<Coeff>& data){
                std::vector<PauliString<Coeff>> converted_data;
                converted_data.reserve(data.size());
                for (Pauli_structure<Coeff> entry : data) {
                        converted_data.push_back(PauliString<Coeff>(entry.first, entry.second));
                }
                this->data = converted_data;
                this->compact_in_place();
        }

        QubitHamiltonian trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state) const {
                std::vector<PauliString<Coeff>> reduced;
                reduced.reserve(this->data.size());
                for (const PauliString<Coeff>& current : this->data) {
                        reduced.push_back(current.trace_out_qubits(qubits, state));
                }
                // The constructor runs compact() which drops zero-coefficient terms,
                // including the sentinel default-constructed PauliString returned for
                // tracing X/Y against a computational-basis state.
                return QubitHamiltonian(std::move(reduced));
        }

        /// Overload accepting arbitrary single-qubit amplitudes (a, b) for each
        /// traced qubit (general superposition states).
        QubitHamiltonian trace_out_qubits(
                const std::vector<int>& qubits,
                const std::vector<std::pair<std::complex<double>, std::complex<double>>>& states
        ) const {
                std::vector<PauliString<Coeff>> reduced;
                reduced.reserve(this->data.size());
                for (const PauliString<Coeff>& current : this->data) {
                        reduced.push_back(current.trace_out_qubits(qubits, states));
                }
                return QubitHamiltonian(std::move(reduced));
        }
        std::string to_string(){
                std::string result = "";
                for (int i = 0; i < this->data.size(); i++) {
                        result.append(data[i].to_string() + "\n");
                }
                return result;
        }


        QubitHamiltonian operator*(std::complex<double> scale) const{
            QubitHamiltonian to_scale = QubitHamiltonian(this->data);
            for (auto &entry : to_scale.data) {
                    entry.coeff *= scale;
                }
                return to_scale;
        }

        
        QubitHamiltonian operator*(int scale) const{
            QubitHamiltonian to_scale = QubitHamiltonian(this->data);
            for (auto &entry : to_scale.data) {
                    entry.coeff *= scale;
                }
                return to_scale;
        }

        QubitHamiltonian operator*(double scale) const{
            QubitHamiltonian to_scale = QubitHamiltonian(this->data);
            for (auto &entry : to_scale.data) {
                    entry.coeff *= scale;
                }
                return to_scale;
        }

        QubitHamiltonian operator*(const QubitHamiltonian& other) const{ 
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        std::vector<PauliString<Coeff>> merged;
                        if (with_small_key(max_word_count(data, other.data),
                                           [&]<size_t W>() { product_small_key<W>(other, merged); })) {
                                return QubitHamiltonian(std::move(merged), already_compacted);
                        }
                }

                const int n1 = static_cast<int>(this->data.size());  // signed for OpenMP 2.0
                const size_t n2 = other.data.size();
                std::vector<PauliString<Coeff>> product(static_cast<size_t>(n1) * n2);

                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
#ifdef PAULIENGINE_HAS_OPENMP
                        #pragma omp parallel for schedule(static) if(static_cast<size_t>(n1) * n2 >= PAULIENGINE_OMP_THRESHOLD)
#endif
                        for (int i = 0; i < n1; ++i) {
                                for (size_t j = 0; j < n2; ++j) {
                                        product[static_cast<size_t>(i) * n2 + j] = this->data[i] * other.data[j];
                                }
                        }
                } else {
                        for (int i = 0; i < n1; ++i) {
                                for (size_t j = 0; j < n2; ++j) {
                                        product[static_cast<size_t>(i) * n2 + j] = this->data[i] * other.data[j];
                                }
                        }
                }
                return QubitHamiltonian(std::move(product));
        }

        QubitHamiltonian operator+(const QubitHamiltonian& other) const {
                std::vector<PauliString<Coeff>> combined;
                combined.reserve(this->data.size() + other.data.size());
                combined.insert(combined.end(), this->data.begin(), this->data.end());
                combined.insert(combined.end(), other.data.begin(), other.data.end());
                return QubitHamiltonian(std::move(combined));  // constructor calls compact()
        }

        QubitHamiltonian operator-(const QubitHamiltonian& other) const {
                return *this + (other * -1.0);
        }

        QubitHamiltonian operator-() const {
                return *this * -1.0;
        }

        bool operator==(const QubitHamiltonian& other) const {
                using Map = std::unordered_map<PauliString<Coeff>, Coeff,
                        PauliStringHash<Coeff>, PauliStringOperatorEqual<Coeff>>;
                Map first_data, second_data;
                for (const auto& ps : this->data) {
                        first_data[ps] = first_data[ps] + ps.coeff;
                }
                for (const auto& ps : other.data) {
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
                const bool iter_this = this->data.size() <= other.data.size();
                const std::vector<PauliString<Coeff>>& outer = iter_this ? this->data : other.data;
                const std::vector<PauliString<Coeff>>& inner = iter_this ? other.data : this->data;

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
                std::vector<PauliString<Coeff>> temp_data;
                temp_data.reserve(this->data.size());
                Coeff coeff;
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        coeff = value;
                } else {
                        coeff = Expression(value);
                }
                for (const auto& ps : this->data) {
                        temp_data.push_back(PauliString<Coeff>(ps.x, ps.y, coeff));
                }
                return QubitHamiltonian(std::move(temp_data));
        }

        QubitHamiltonian diff(std::string symbol) const{
                std::vector<PauliString<Coeff>> data;
                for (const PauliString<Coeff>& ps : this->data) {
                        PauliString<Coeff> temp = ps.diff(symbol);
                        if (!coeff_is_zero(temp.coeff)) {
                                data.push_back(temp);
                        }
                }
                return QubitHamiltonian(std::move(data));
        }

        QubitHamiltonian substitute(const std::unordered_map<std::string, std::complex<double>>& substitution_map) const {
                std::vector<PauliString<Coeff>> temp_data;
                temp_data.reserve(this->data.size());
                for (const auto &ps : data) {
                        temp_data.push_back(ps.substitute(substitution_map));
                }
                return QubitHamiltonian(std::move(temp_data));
        }

        // Coefficient of the term whose operator part equals `ps`, or zero if the
        // Hamiltonian has no such term. The coefficient of `ps` itself is ignored,
        // so a query string of either coefficient type works. compact() keeps the
        // operator parts unique, hence at most one term can match.
        template<typename C>
        Coeff contains(const PauliString<C>& ps) const {
                for (const auto& term : this->data) {
                        if (term.x == ps.x && term.y == ps.y) {
                                return term.coeff;
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

    private:
        void compact_in_place() {
                if (data.empty()) return;

                size_t max_words = 0;
                for (const auto& ps : data) max_words = std::max(max_words, ps.x.size());
                if (with_small_key(max_words, [this]<size_t W>() { compact_small_key<W>(); })) {
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

                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        parallel_sort_vec(data, op_less);
                } else {
                        // symengine coefficients are not safe to move across threads
                        std::sort(data.begin(), data.end(), op_less);
                }

                size_t write = 0;
                for (size_t read = 1; read < data.size(); ++read) {
                        if (op_equal(data[write], data[read])) {
                                data[write].coeff = data[write].coeff + data[read].coeff;
                        } else {
                                ++write;
                                if (write != read) {
                                        data[write] = std::move(data[read]);
                                }
                        }
                }
                data.resize(write + 1);

                data.erase(std::remove_if(data.begin(), data.end(),
                        [](const PauliString<Coeff>& ps) { return coeff_is_zero(ps.coeff); }),
                        data.end());
        }

        static bool coeff_is_zero(const Coeff& c) {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        return c == std::complex<double>(0.0, 0.0);
                } else {
                    const auto& basic = *c.get_basic();
                        return SymEngine::is_a_Number(basic)
                                && SymEngine::down_cast<const SymEngine::Number&>(basic).is_zero();
                }
        }

        static bool coeff_equal(const Coeff& a, const Coeff& b) {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        return a == b;
                } else {
                        return coeff_is_zero(SymEngine::expand(a - b));
                }
        }

        /// Largest word count handled by the small-key fast path (8 words =
        /// 512 qubits). Beyond this, compact() falls back to the generic sort.
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

        static size_t max_word_count(const std::vector<PauliString<Coeff>>& a,
                                     const std::vector<PauliString<Coeff>>& b) {
                size_t max_words = 0;
                for (const auto& ps : a) max_words = std::max(max_words, ps.x.size());
                for (const auto& ps : b) max_words = std::max(max_words, ps.x.size());
                return max_words;
        }

        /// Flat, trivially copyable form of a complex-coefficient term.
        template<size_t W>
        struct KeyedTerm {
                uint64_t key[2 * W];
                double re;
                double im;
        };

        template<size_t W>
        static void set_keyed_term(KeyedTerm<W>& term, const PauliString<Coeff>& ps) {
                const size_t nw = ps.x.size();
                for (size_t k = 0; k < W; ++k) {
                        term.key[k]     = (k < nw) ? ps.x[k] : 0;
                        term.key[W + k] = (k < nw) ? ps.y[k] : 0;
                }
                term.re = ps.coeff.real();
                term.im = ps.coeff.imag();
        }

        /// Writes the operator part of a * b to `term` without building a PauliString
        /// and returns the exponent k of the phase i^k (odd: a and b anti-commute).
        template<size_t W>
        static int set_keyed_product_key(KeyedTerm<W>& term, const PauliString<Coeff>& a, const PauliString<Coeff>& b) {
                const int power_of_i = multiply_symplectic_words(
                        a.x.data(), a.y.data(), a.x.size(),
                        b.x.data(), b.y.data(), b.x.size(),
                        term.key, term.key + W);
                for (size_t k = std::max(a.x.size(), b.x.size()); k < W; ++k) {
                        term.key[k] = 0;
                        term.key[W + k] = 0;
                }
                return power_of_i;
        }

        /// Sets the coefficient of `term` to scale * i^power_of_i * a.coeff * b.coeff.
        template<size_t W>
        static void set_keyed_product_coeff(KeyedTerm<W>& term, const PauliString<Coeff>& a, const PauliString<Coeff>& b,
                                            int power_of_i, double scale) {
                static constexpr double kRealSign[4] = {1.0, -1.0, -1.0, 1.0};
                static constexpr double kImagSign[4] = {1.0, 1.0, -1.0, -1.0};
                const double ar = a.coeff.real(), ai = a.coeff.imag();
                const double br = b.coeff.real(), bi = b.coeff.imag();
                const double parts[2] = {ar * br - ai * bi, ar * bi + ai * br};
                const int k = power_of_i & 3;
                term.re = scale * kRealSign[k] * parts[k & 1];
                term.im = scale * kImagSign[k] * parts[(k & 1) ^ 1];
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

        /// Writes terms into existing strings, reusing their storage.
        template<size_t W>
        static void write_keyed_terms(const KeyedTerm<W>* terms, size_t count, PauliString<Coeff>* out) {
                for (size_t i = 0; i < count; ++i) {
                        const KeyedTerm<W>& term = terms[i];
                        size_t nw = W;
                        while (nw > 0 && term.key[nw - 1] == 0 && term.key[W + nw - 1] == 0) --nw;
                        out[i].x.assign(term.key, term.key + nw);
                        out[i].y.assign(term.key + W, term.key + W + nw);
                        out[i].coeff = Coeff(term.re, term.im);
                }
        }

        /// Sorts `terms` by operator, adds the coefficients of equal operators and
        /// writes the non-zero results to `out`, ordered by operator. Large inputs
        /// use a sample sort: regular samples give splitters, every term moves to the
        /// bucket of its operator, and buckets are then sorted and merged independently.
        /// Equal operators always share a bucket, so no serial merge step is left.
        /// The result does not depend on the number of threads. `terms` is overwritten.
        template<size_t W>
        static void merge_keyed_terms(KeyedTerm<W>* terms, size_t n, std::vector<PauliString<Coeff>>& out) {
                using Term = KeyedTerm<W>;
#ifdef PAULIENGINE_HAS_OPENMP
                const int nthreads = omp_get_max_threads();
                if (n >= PAULIENGINE_OMP_THRESHOLD && nthreads > 1) {
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
                                samples[i] = terms[i * (n / n_samples)];
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
                        reserve_output(out, out_begin[n_buckets]);
                        out.resize(out_begin[n_buckets]);
                        #pragma omp parallel for schedule(dynamic)
                        for (int b = 0; b < n_buckets; ++b) {
                                write_keyed_terms<W>(sorted + bucket_begin[b], kept[b], out.data() + out_begin[b]);
                        }
                        return;
                }
#endif
                std::unique_ptr<Term[]> spare(n > kInsertionSortLimit ? new Term[n] : nullptr);
                const size_t kept = sort_and_merge_keyed<W>(terms, spare.get(), n);
                out.resize(kept);
                write_keyed_terms<W>(terms, kept, out.data());
        }

#ifdef PAULIENGINE_HAS_OPENMP
        /// Grows the capacity of `out` and lets all threads fault its pages in. The
        /// serial resize() afterwards only has to write to memory that is already mapped.
        static void reserve_output(std::vector<PauliString<Coeff>>& out, size_t total) {
                if (total <= out.capacity()) return;
                out.reserve(total);
                constexpr std::ptrdiff_t kPageSize = 4096;
                char* const storage = reinterpret_cast<char*>(out.data());
                const std::ptrdiff_t pages = static_cast<std::ptrdiff_t>(total * sizeof(PauliString<Coeff>)) / kPageSize;
                #pragma omp parallel for schedule(static)
                for (std::ptrdiff_t page = 0; page < pages; ++page) {
                        *static_cast<volatile char*>(storage + page * kPageSize) = 0;
                }
        }
#endif

        template<size_t W>
        void product_small_key(const QubitHamiltonian& other, std::vector<PauliString<Coeff>>& out) const {
                const int n1 = static_cast<int>(data.size());  // signed for OpenMP 2.0
                const size_t n2 = other.data.size();
                const size_t pairs = static_cast<size_t>(n1) * n2;
                // new[] leaves the terms untouched, so the threads fault their own pages in.
                std::unique_ptr<KeyedTerm<W>[]> terms(new KeyedTerm<W>[pairs]);
                KeyedTerm<W>* const raw = terms.get();
#ifdef PAULIENGINE_HAS_OPENMP
                #pragma omp parallel for schedule(static) if(pairs >= PAULIENGINE_OMP_THRESHOLD)
#endif
                for (int i = 0; i < n1; ++i) {
                        const PauliString<Coeff>& a = data[i];
                        KeyedTerm<W>* const row = raw + static_cast<size_t>(i) * n2;
                        for (size_t j = 0; j < n2; ++j) {
                                const PauliString<Coeff>& b = other.data[j];
                                const int power_of_i = set_keyed_product_key<W>(row[j], a, b);
                                set_keyed_product_coeff<W>(row[j], a, b, power_of_i, 1.0);
                        }
                }
                merge_keyed_terms<W>(raw, pairs, out);
        }

        template<size_t W>
        void commutator_small_key(const QubitHamiltonian& other, std::vector<PauliString<Coeff>>& out) const {
                const int n1 = static_cast<int>(data.size());  // signed for OpenMP 2.0
                const size_t n2 = other.data.size();
                const size_t pairs = static_cast<size_t>(n1) * n2;
#ifdef PAULIENGINE_HAS_OPENMP
                const int nthreads = pairs >= PAULIENGINE_OMP_THRESHOLD ? omp_get_max_threads() : 1;
#else
                const int nthreads = 1;
#endif
                std::vector<std::vector<KeyedTerm<W>>> local(nthreads);
#ifdef PAULIENGINE_HAS_OPENMP
                #pragma omp parallel if(nthreads > 1)
#endif
                {
                        // Filled as a local object: the vectors in `local` share cache lines.
                        std::vector<KeyedTerm<W>> mine;
                        mine.reserve(pairs / (8 * static_cast<size_t>(nthreads)) + 1);
#ifdef PAULIENGINE_HAS_OPENMP
                        #pragma omp for schedule(static)
#endif
                        for (int i = 0; i < n1; ++i) {
                                const PauliString<Coeff>& a = data[i];
                                for (size_t j = 0; j < n2; ++j) {
                                        const PauliString<Coeff>& b = other.data[j];
                                        KeyedTerm<W> term;
                                        const int power_of_i = set_keyed_product_key<W>(term, a, b);
                                        if (power_of_i & 1) {  // commuting pairs contribute nothing
                                                set_keyed_product_coeff<W>(term, a, b, power_of_i, 2.0);
                                                mine.push_back(term);
                                        }
                                }
                        }
#ifdef PAULIENGINE_HAS_OPENMP
                        local[omp_get_thread_num()] = std::move(mine);
#else
                        local[0] = std::move(mine);
#endif
                }

                if (nthreads == 1) {
                        merge_keyed_terms<W>(local[0].data(), local[0].size(), out);
                        return;
                }
                std::vector<size_t> offset(nthreads + 1, 0);
                for (int t = 0; t < nthreads; ++t) offset[t + 1] = offset[t] + local[t].size();
                std::unique_ptr<KeyedTerm<W>[]> terms(new KeyedTerm<W>[offset[nthreads]]);
                KeyedTerm<W>* const raw = terms.get();
#ifdef PAULIENGINE_HAS_OPENMP
                #pragma omp parallel for schedule(static)
#endif
                for (int t = 0; t < nthreads; ++t) {
                        // A static schedule hands out ascending i, so this is the serial order.
                        std::copy(local[t].begin(), local[t].end(), raw + offset[t]);
                        std::vector<KeyedTerm<W>>().swap(local[t]);
                }
                merge_keyed_terms<W>(raw, offset[nthreads], out);
        }

        template<size_t W>
        void compact_small_key() {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(data.size());
                        std::unique_ptr<KeyedTerm<W>[]> terms(new KeyedTerm<W>[data.size()]);
                        KeyedTerm<W>* const raw = terms.get();
#ifdef PAULIENGINE_HAS_OPENMP
                        #pragma omp parallel for schedule(static) if(data.size() >= PAULIENGINE_OMP_THRESHOLD)
#endif
                        for (std::ptrdiff_t i = 0; i < n; ++i) {
                                set_keyed_term<W>(raw[i], data[i]);
                        }
                        merge_keyed_terms<W>(raw, data.size(), data);
                } else {
                        compact_small_key_symbolic<W>();
                }
        }

        template<size_t W>
        void compact_small_key_symbolic() {
                struct Entry {
                        uint64_t key[2 * W];
                        uint32_t idx;
                };
                const size_t n = data.size();
                std::vector<Entry> entries(n);
                for (size_t i = 0; i < n; ++i) {
                        Entry& e = entries[i];
                        const PauliString<Coeff>& ps = data[i];
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

                parallel_sort_vec(entries, key_less);

                std::vector<PauliString<Coeff>> out;
                out.reserve(n);
                size_t i = 0;
                while (i < n) {
                        Coeff acc = std::move(data[entries[i].idx].coeff);
                        size_t j = i + 1;
                        while (j < n && key_equal(entries[i], entries[j])) {
                                acc = acc + data[entries[j].idx].coeff;
                                ++j;
                        }
                        if (!coeff_is_zero(acc)) {
                                PauliString<Coeff> ps = std::move(data[entries[i].idx]);
                                ps.coeff = std::move(acc);
                                out.push_back(std::move(ps));
                        }
                        i = j;
                }
                data = std::move(out);
        }

        /// Chunked parallel merge sort; falls back to std::sort without OpenMP
        /// or for small inputs. Elements must be safe to move across threads.
        template<typename T, typename Cmp>
        static void parallel_sort_vec(std::vector<T>& v, Cmp cmp) {
#ifdef PAULIENGINE_HAS_OPENMP
                const size_t n = v.size();
                const int nthreads = omp_get_max_threads();
                if (n >= PAULIENGINE_OMP_THRESHOLD && nthreads > 1) {
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
#endif
                std::sort(v.begin(), v.end(), cmp);
        }
    public:

        Hamiltonian_structure<Coeff> to_dictionary() const{
                Hamiltonian_structure<Coeff> output;
                output.reserve(this->data.size());
                for (const auto& entry : this->data) {
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

        std::vector<PauliString<Coeff>> to_list(){
                return this->data;
        }

        QubitHamiltonian commutator(const QubitHamiltonian& other) const {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        std::vector<PauliString<Coeff>> merged;
                        if (with_small_key(max_word_count(data, other.data),
                                           [&]<size_t W>() { commutator_small_key<W>(other, merged); })) {
                                return QubitHamiltonian(std::move(merged), already_compacted);
                        }
                }

                const int n1 = static_cast<int>(this->data.size()); 
                const size_t n2 = other.data.size();
                const size_t pairs = static_cast<size_t>(n1) * n2;


                std::vector<PauliString<Coeff>> terms;

                const bool parallel_ok =
#ifdef PAULIENGINE_HAS_OPENMP
                        std::is_same_v<Coeff, std::complex<double>>
                        && pairs >= PAULIENGINE_OMP_THRESHOLD;
#else
                        false;
#endif

                if (!parallel_ok) {
                        terms.reserve(pairs / 8 + 1);
                        for (int i = 0; i < n1; ++i) {
                                for (size_t j = 0; j < n2; ++j) {
                                        PauliString<Coeff> t = this->data[i].commutator(other.data[j]);
                                        if (!t.x.empty()) {
                                                terms.push_back(std::move(t));
                                        }
                                }
                        }
                } else {
#ifdef PAULIENGINE_HAS_OPENMP
                        const int nthreads = omp_get_max_threads();
                        std::vector<std::vector<PauliString<Coeff>>> local(nthreads);
                        for (auto& v : local) v.reserve(pairs / (8 * nthreads) + 1);

                        #pragma omp parallel for schedule(static)
                        for (int i = 0; i < n1; ++i) {
                                auto& mine = local[omp_get_thread_num()];
                                for (size_t j = 0; j < n2; ++j) {
                                        PauliString<Coeff> t = this->data[i].commutator(other.data[j]);
                                        if (!t.x.empty()) {
                                                mine.push_back(std::move(t));
                                        }
                                }
                        }

                        size_t total = 0;
                        for (const auto& v : local) total += v.size();
                        terms.reserve(total);
                        for (auto& v : local) {
                                for (auto& t : v) terms.push_back(std::move(t));
                                v.clear();
                        }
#endif
                }
                return QubitHamiltonian(std::move(terms));  // constructor calls compact()
        }

        size_t size() const { return data.size(); }

        std::vector<int> qubits() const {
                std::set<int> qs;
                for (const auto& ps : data) {
                        for (size_t i = 0; i < ps.x.size(); ++i) {
                                uint64_t xy = ps.x[i] | ps.y[i];
                                while (xy) {
                                        int j = std::countr_zero(xy);
                                        qs.insert((int)(i * BITS_IN_INTEGER + j));
                                        xy &= xy - 1;
                                }
                        }
                }
                return std::vector<int>(qs.begin(), qs.end());
        }

        int n_qubits() const { return (int)qubits().size(); }

        bool is_all_z() const {
                for (const auto& ps : data) {
                        if (!ps.is_all_z()) return false;
                }
                return true;
        }

        int count_measurements() const {
                return is_all_z() ? 1 : (int)data.size();
        }

        bool is_hermitian() const {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        for (const auto& ps : data) {
                                if (std::abs(ps.coeff.imag()) > 1e-12) return false;
                        }
                        return true;
                } else {
                        try {
                                for (const auto& ps : data) {
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
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        for (const auto& ps : data) {
                                if (std::abs(ps.coeff.real()) > 1e-12) return false;
                        }
                        return true;
                } else {
                        try {
                                for (const auto& ps : data) {
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
                return this->data;
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
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (const auto& ps : data) {
                        Coeff conj_coeff;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                conj_coeff = std::conj(ps.coeff);
                        } else {
                                conj_coeff = Expression(SymEngine::conjugate(ps.coeff.get_basic()));
                        }
                        result.push_back(PauliString<Coeff>(ps.x, ps.y, conj_coeff));
                }
                return QubitHamiltonian(std::move(result));
        }

        QubitHamiltonian conjugate() const {
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (const auto& ps : data) {
                        int y_count = ps.count_y();
                        Coeff new_coeff;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                new_coeff = std::conj(ps.coeff);
                        } else {
                                new_coeff = Expression(SymEngine::conjugate(ps.coeff.get_basic()));
                        }
                        if (y_count % 2 == 1) {
                                new_coeff = new_coeff * -1.0;
                        }
                        result.push_back(PauliString<Coeff>(ps.x, ps.y, new_coeff));
                }
                return QubitHamiltonian(std::move(result));
        }

        QubitHamiltonian transpose() const {
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (const auto& ps : data) {
                        int y_count = ps.count_y();
                        Coeff new_coeff = ps.coeff;
                        if (y_count % 2 == 1) {
                                new_coeff = new_coeff * -1.0;
                        }
                        result.push_back(PauliString<Coeff>(ps.x, ps.y, new_coeff));
                }
                return QubitHamiltonian(std::move(result));
        }

        QubitHamiltonian simplify(double threshold = 0.0) const {
                std::vector<PauliString<Coeff>> result;
                for (const auto& ps : data) {
                        bool keep;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
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
                for (const auto& ps : data) {
                        std::complex<double> c;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
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
                                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                        real_coeff = std::complex<double>(c.real(), 0.0);
                                } else {
                                        real_coeff = Expression(c.real());
                                }
                                herm_terms.push_back(PauliString<Coeff>(ps.x, ps.y, real_coeff));
                        }
                        if (c.imag() != 0.0) {
                                Coeff imag_coeff;
                                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
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
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (auto ps : data) {
                        result.push_back(ps.map_qubits(qubit_map));
                }
                return QubitHamiltonian(std::move(result));
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

                for (const auto& ps : data) {
                        std::complex<double> c;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
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
        static Coeff make_one() {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        return std::complex<double>(1.0, 0.0);
                } else {
                        return Coeff(1);
                }
        }

        static Coeff make_zero() {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
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
