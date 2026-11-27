#pragma once
/// @file
/// A sum of Pauli strings stored as lists of sites (see SparsePauliString.h).
///
/// QubitHamiltonian stores every term with the bit vectors of the widest one, so
/// a local Hamiltonian on n qubits takes O(n / 64) words per term and every pair
/// in a product or commutator costs as much. Here a term takes only its w sites,
/// and pairs cost O(w_a + w_b):
///
/// - products and squares form all pairs, merge the site lists and add up equal
///   operators in hash buckets (detail/SparseSites.h), in parallel and with a
///   result that does not depend on the number of threads;
/// - commutators visit only pairs that share a qubit when those are few, found
///   through an index of the other sum by qubit;
/// - targeted commutators, coefficient lookup and inner products join on the hash
///   of the operator.
///
/// The sum is kept compacted: no two terms share an operator and zero terms are
/// dropped. Terms are ordered by the hash of their operator, deterministically.

#include <algorithm>
#include <complex>
#include <cstdint>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/SparsePauliString.h"
#include "pauliengine/detail/SparseSites.h"

template <typename Coeff>
class SparsePauliStringSum {
    public:
        using Site = pauliengine_detail::Site;
        using Terms = pauliengine_detail::SparseTerms<Coeff>;
        using String = SparsePauliString<Coeff>;

        /// The zero operator.
        SparsePauliStringSum() = default;

        /// From terms that may repeat operators: equal operators are merged, zeros dropped.
        explicit SparsePauliStringSum(const std::vector<String>& terms) {
                size_t heaviest = 0;
                for (const auto& t : terms) heaviest = std::max(heaviest, t.sites.size());
                terms_ = pauliengine_detail::accumulate<Coeff>(heaviest, terms.size(), 1, 4, [&](size_t i, auto& emit) {
                        emit(terms[i].sites.data(), terms[i].sites.size(), terms[i].coeff);
                });
        }

        /// From Pauli strings with bit vectors.
        explicit SparsePauliStringSum(const std::vector<PauliString<Coeff>>& terms) {
                std::vector<String> sparse;
                sparse.reserve(terms.size());
                for (const auto& ps : terms) sparse.emplace_back(ps);
                *this = SparsePauliStringSum(sparse);
        }

        /// From a list of (coefficient, {qubit: "X"/"Y"/"Z"}) pairs.
        explicit SparsePauliStringSum(const Hamiltonian_structure<Coeff>& data) {
                std::vector<String> sparse;
                sparse.reserve(data.size());
                for (const auto& [c, ops] : data) sparse.emplace_back(c, ops);
                *this = SparsePauliStringSum(sparse);
        }

        /// The same operator as a QubitHamiltonian.
        explicit SparsePauliStringSum(const QubitHamiltonian<Coeff>& dense)
                : SparsePauliStringSum(dense.paulistrings()) {}

        /// Already compacted terms.
        static SparsePauliStringSum from_terms(Terms terms) {
                SparsePauliStringSum h;
                h.terms_ = std::move(terms);
                return h;
        }

        static SparsePauliStringSum zero() { return SparsePauliStringSum(); }
        static SparsePauliStringSum unit() {
                return SparsePauliStringSum(std::vector<String>{String(pauliengine_detail::coeff_one<Coeff>())});
        }

        const Terms& terms() const { return terms_; }
        size_t size() const { return terms_.size(); }

        String term(size_t i) const {
                String s(terms_.coeffs[i]);
                s.sites.assign(terms_.sites_of(i), terms_.sites_of(i) + terms_.weight_of(i));
                return s;
        }

        std::vector<String> paulistrings() const {
                std::vector<String> out;
                out.reserve(size());
                for (size_t i = 0; i < size(); ++i) out.push_back(term(i));
                return out;
        }

        /// The operator as a QubitHamiltonian (bit vectors up to the highest qubit).
        QubitHamiltonian<Coeff> to_dense() const {
                std::vector<PauliString<Coeff>> dense;
                dense.reserve(size());
                std::vector<uint64_t> x, y;
                for (size_t i = 0; i < size(); ++i) {
                        pauliengine_detail::words_of_sites(terms_.sites_of(i), terms_.weight_of(i), x, y);
                        dense.emplace_back(x, y, terms_.coeffs[i]);
                }
                return QubitHamiltonian<Coeff>(std::move(dense));
        }

        Hamiltonian_structure<Coeff> to_dictionary() const {
                Hamiltonian_structure<Coeff> out;
                out.reserve(size());
                for (size_t i = 0; i < size(); ++i) out.push_back(term(i).to_dictionary());
                return out;
        }

        std::string to_string() const {
                std::ostringstream oss;
                for (size_t i = 0; i < size(); ++i) oss << term(i).to_string() << "\n";
                return oss.str();
        }

        // ------------------------------------------------------------------ arithmetic

        SparsePauliStringSum operator+(const SparsePauliStringSum& other) const { return combine(other, false); }
        SparsePauliStringSum operator-(const SparsePauliStringSum& other) const { return combine(other, true); }
        SparsePauliStringSum operator-() const { return scaled(Coeff(-1.0)); }

        SparsePauliStringSum operator*(std::complex<double> s) const { return scaled(Coeff(s)); }
        SparsePauliStringSum operator*(double s) const { return scaled(Coeff(s)); }
        SparsePauliStringSum operator*(int s) const { return scaled(Coeff(static_cast<double>(s))); }

        /// Product this * other; h * h (the same object) uses the square formula.
        SparsePauliStringSum operator*(const SparsePauliStringSum& other) const {
                if (this == &other) return square();
                const Terms& a = terms_;
                const Terms& b = other.terms_;
                const size_t n2 = b.size();
                const size_t work = n2 * (average_weight() + other.average_weight() + 1);
                return from_terms(pauliengine_detail::accumulate<Coeff>(max_weight() + other.max_weight(), a.size(), n2, work, [&](size_t i, auto& emit) {
                        const Site* sa = a.sites_of(i);
                        const size_t wa = a.weight_of(i);
                        auto& scratch = emit.scratch();
                        for (size_t j = 0; j < n2; ++j) {
                                const size_t wb = b.weight_of(j);
                                if (scratch.size() < wa + wb) scratch.resize(wa + wb);
                                size_t n = 0;
                                const int k = pauliengine_detail::multiply_sites(sa, wa, b.sites_of(j), wb, scratch.data(), n);
                                emit(scratch.data(), n, pauliengine_detail::phase_product(a.coeffs[i], b.coeffs[j], k, 1.0));
                        }
                }));
        }

        /// this * this = sum_i a_i^2 I + 2 sum_{i < j, commuting} a_i a_j P_i P_j: the
        /// anti-commuting pairs cancel, so each unordered pair is visited once.
        SparsePauliStringSum square() const {
                const Terms& a = terms_;
                const size_t n = a.size();
                const size_t work = n * (2 * average_weight() + 1) / 2;
                return from_terms(pauliengine_detail::accumulate<Coeff>(2 * max_weight(), n, n / 2 + 1, work, [&](size_t i, auto& emit) {
                        const Site* si = a.sites_of(i);
                        const size_t wi = a.weight_of(i);
                        emit(nullptr, 0, pauliengine_detail::phase_product(a.coeffs[i], a.coeffs[i], 0, 1.0));
                        auto& scratch = emit.scratch();
                        for (size_t j = i + 1; j < n; ++j) {
                                const size_t wj = a.weight_of(j);
                                if (scratch.size() < wi + wj) scratch.resize(wi + wj);
                                size_t m = 0;
                                const int k = pauliengine_detail::multiply_sites(si, wi, a.sites_of(j), wj, scratch.data(), m);
                                if (k & 1) continue;  // anti-commuting pairs cancel
                                emit(scratch.data(), m, pauliengine_detail::phase_product(a.coeffs[i], a.coeffs[j], k, 2.0));
                        }
                }));
        }

        SparsePauliStringSum power(int n) const {
                if (n < 0) throw std::runtime_error("Negative power not supported for SparsePauliStringSum");
                if (n == 0) return unit();
                SparsePauliStringSum result;
                bool have = false;
                SparsePauliStringSum base = *this;
                while (n > 0) {
                        if (n & 1) {
                                result = have ? result * base : base;
                                have = true;
                        }
                        n >>= 1;
                        if (n > 0) base = base.square();
                }
                return result;
        }

        /// [this, other] = 2 sum over anti-commuting pairs of a_i b_j P_i P_j. When
        /// few pairs share a qubit, only those are visited.
        SparsePauliStringSum commutator(const SparsePauliStringSum& other) const {
                const Terms& a = terms_;
                const Terms& b = other.terms_;
                const size_t n1 = a.size();
                const size_t n2 = b.size();
                if (n1 == 0 || n2 == 0) return SparsePauliStringSum();
                const QubitIndex index(b);
                const double shared = index.shared_pairs(a);
                const bool listed = 8.0 * shared < static_cast<double>(n1) * static_cast<double>(n2);
                const size_t per_row = listed ? static_cast<size_t>(shared / static_cast<double>(n1)) + 1 : n2;
                const size_t work = per_row * (average_weight() + other.average_weight() + 1);
                return from_terms(pauliengine_detail::accumulate<Coeff>(max_weight() + other.max_weight(), n1, per_row, work, [&](size_t i, auto& emit) {
                        const Site* sa = a.sites_of(i);
                        const size_t wa = a.weight_of(i);
                        auto& scratch = emit.scratch();
                        const auto visit = [&](size_t j) {
                                const Site* sb = b.sites_of(j);
                                const size_t wb = b.weight_of(j);
                                if (!pauliengine_detail::anticommute_sites(sa, wa, sb, wb)) return;
                                if (scratch.size() < wa + wb) scratch.resize(wa + wb);
                                size_t n = 0;
                                const int k = pauliengine_detail::multiply_sites(sa, wa, sb, wb, scratch.data(), n);
                                emit(scratch.data(), n, pauliengine_detail::phase_product(a.coeffs[i], b.coeffs[j], k, 2.0));
                        };
                        if (listed) {
                                index.for_each_partner(sa, wa, visit);
                        } else {
                                for (size_t j = 0; j < n2; ++j) visit(j);
                        }
                }));
        }

        /// Coefficient of target's operator in [this, other], via a hash join:
        /// for each term P of the smaller sum the only partner that can give target
        /// is Q = P * target (up to phase). O((|A| + |B|) w) expected.
        Coeff targeted_commutator(const SparsePauliStringSum& other, const String& target) const {
                Coeff total = pauliengine_detail::coeff_zero<Coeff>();
                const bool iter_this = size() <= other.size();
                const SparsePauliStringSum& outer = iter_this ? *this : other;
                const SparsePauliStringSum& inner = iter_this ? other : *this;
                if (outer.size() == 0 || inner.size() == 0) return total;
                const OperatorIndex index(inner.terms_);
                std::vector<Site> partner;
                const Site* st = target.sites.data();
                const size_t wt = target.sites.size();
                for (size_t i = 0; i < outer.size(); ++i) {
                        const Site* sp = outer.terms_.sites_of(i);
                        const size_t wp = outer.terms_.weight_of(i);
                        partner.resize(wp + wt);
                        size_t n = 0;
                        pauliengine_detail::multiply_sites(sp, wp, st, wt, partner.data(), n);
                        const size_t j = index.find(partner.data(), n);
                        if (j == OperatorIndex::npos) continue;
                        const Site* sq = inner.terms_.sites_of(j);
                        const size_t wq = inner.terms_.weight_of(j);
                        // keep the order [this, other]
                        const int k = iter_this ? pauliengine_detail::product_phase_sites(sp, wp, sq, wq)
                                                : pauliengine_detail::product_phase_sites(sq, wq, sp, wp);
                        if (!(k & 1)) continue;
                        const Coeff& left = iter_this ? outer.terms_.coeffs[i] : inner.terms_.coeffs[j];
                        const Coeff& right = iter_this ? inner.terms_.coeffs[j] : outer.terms_.coeffs[i];
                        total = total + pauliengine_detail::phase_product(left, right, k, 2.0);
                }
                return total;
        }

        // ------------------------------------------------------------------ queries

        /// Coefficient of the operator of `ps` (0 if absent); its coefficient is ignored.
        Coeff contains(const String& ps) const {
                const OperatorIndex index(terms_);
                const size_t i = index.find(ps.sites.data(), ps.sites.size());
                return i == OperatorIndex::npos ? pauliengine_detail::coeff_zero<Coeff>() : terms_.coeffs[i];
        }

        /// Normalized Hilbert-Schmidt inner product sum_P conj(a_P) b_P.
        Coeff inner_product(const SparsePauliStringSum& other) const {
                Coeff total = pauliengine_detail::coeff_zero<Coeff>();
                const bool index_other = other.size() >= size();
                const Terms& probe = index_other ? terms_ : other.terms_;
                const Terms& indexed = index_other ? other.terms_ : terms_;
                const OperatorIndex index(indexed);
                for (size_t i = 0; i < probe.size(); ++i) {
                        const size_t j = index.find(probe.sites_of(i), probe.weight_of(i));
                        if (j == OperatorIndex::npos) continue;
                        const Coeff& mine = index_other ? probe.coeffs[i] : indexed.coeffs[j];
                        const Coeff& theirs = index_other ? indexed.coeffs[j] : probe.coeffs[i];
                        total = total + pauliengine_detail::conj_of(mine) * theirs;
                }
                return total;
        }

        bool operator==(const SparsePauliStringSum& other) const {
                if (size() != other.size()) return false;
                const OperatorIndex index(other.terms_);
                for (size_t i = 0; i < size(); ++i) {
                        const size_t j = index.find(terms_.sites_of(i), terms_.weight_of(i));
                        if (j == OperatorIndex::npos
                            || !pauliengine_detail::coeff_equal_to(terms_.coeffs[i], other.terms_.coeffs[j])) {
                                return false;
                        }
                }
                return true;
        }

        std::vector<int> qubits() const {
                std::vector<int> q;
                q.reserve(terms_.sites.size());
                for (Site s : terms_.sites) q.push_back(static_cast<int>(pauliengine_detail::site_qubit(s)));
                std::sort(q.begin(), q.end());
                q.erase(std::unique(q.begin(), q.end()), q.end());
                return q;
        }

        int n_qubits() const { return static_cast<int>(qubits().size()); }

        /// Highest qubit index + 1 (0 for the zero operator or the identity).
        size_t span() const {
                pauliengine_detail::Site top = 0;
                bool any = false;
                for (size_t i = 0; i < size(); ++i) {
                        if (terms_.weight_of(i) == 0) continue;
                        top = std::max(top, terms_.sites_of(i)[terms_.weight_of(i) - 1]);
                        any = true;
                }
                return any ? pauliengine_detail::site_qubit(top) + size_t{1} : 0;
        }

        /// Whether the operator still fits the sparse representation (see
        /// pauliengine_detail::sparse_stays); otherwise to_dense() is the better form.
        bool fits_sparse() const {
                return pauliengine_detail::sparse_stays(static_cast<double>(span()), static_cast<double>(max_weight()));
        }

        bool is_all_z() const {
                return std::all_of(terms_.sites.begin(), terms_.sites.end(),
                                   [](Site s) { return pauliengine_detail::site_pauli(s) == 3; });
        }

        int count_measurements() const { return is_all_z() ? 1 : static_cast<int>(size()); }

        // ------------------------------------------------------------------ coefficients

        /// Terms with |coefficient| <= threshold dropped (symbolic ones that do not
        /// evaluate to a number are kept).
        SparsePauliStringSum simplify(double threshold = 0.0) const {
                return filtered([&](const Coeff& c) {
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                return std::abs(c) > threshold;
                        } else {
                                try {
                                        return std::abs(PauliString<Coeff>::to_complex(c)) > threshold;
                                } catch (...) {
                                        return true;
                                }
                        }
                });
        }

        SparsePauliStringSum set_all_coeff(const Coeff& value) const {
                SparsePauliStringSum out = *this;
                if (pauliengine_detail::coeff_is_zero(value)) return SparsePauliStringSum();
                for (auto& c : out.terms_.coeffs) c = value;
                return out;
        }

        /// Hermitian conjugate: every coefficient conjugated.
        SparsePauliStringSum dagger() const {
                return mapped_coeffs([](size_t, const Coeff& c) { return pauliengine_detail::conj_of(c); });
        }

        /// Complex conjugate: conjugated coefficients, sign flipped for an odd number of Y.
        SparsePauliStringSum conjugate() const {
                return mapped_coeffs([&](size_t i, const Coeff& c) {
                        const Coeff cc = pauliengine_detail::conj_of(c);
                        return odd_y(i) ? pauliengine_detail::coeff_zero<Coeff>() - cc : cc;
                });
        }

        /// Transpose: sign flipped for an odd number of Y.
        SparsePauliStringSum transpose() const {
                return mapped_coeffs([&](size_t i, const Coeff& c) {
                        return odd_y(i) ? pauliengine_detail::coeff_zero<Coeff>() - c : c;
                });
        }

        SparsePauliStringSum diff(const std::string& name) const {
                SparsePauliStringSum out = mapped_coeffs([&](size_t, const Coeff& c) {
                        return pauliengine_detail::coeff_diff(c, name);
                });
                return out.filtered([](const Coeff& c) { return !pauliengine_detail::coeff_is_zero(c); });
        }

        /// Values or polynomials substituted for symbols of every coefficient.
        template <typename Values>
        SparsePauliStringSum substitute(const Values& values) const {
                const auto prepared = Coeff::prepare_subs(values);
                SparsePauliStringSum out = mapped_coeffs([&](size_t, const Coeff& c) {
                        return pauliengine_detail::coeff_subs(c, prepared);
                });
                return out.filtered([](const Coeff& c) { return !pauliengine_detail::coeff_is_zero(c); });
        }

        /// Copy with qubit indices remapped (see PauliString::map_qubits); terms that
        /// become equal are merged.
        SparsePauliStringSum map_qubits(const std::unordered_map<int, int>& qubit_map) const {
                std::vector<String> mapped;
                mapped.reserve(size());
                for (size_t i = 0; i < size(); ++i) mapped.push_back(term(i).map_qubits(qubit_map));
                return SparsePauliStringSum(mapped);
        }

    private:
        Terms terms_;

        size_t average_weight() const { return size() ? terms_.sites.size() / size() : 0; }

        size_t max_weight() const {
                size_t w = 0;
                for (size_t i = 0; i < size(); ++i) w = std::max(w, terms_.weight_of(i));
                return w;
        }

        bool odd_y(size_t i) const {
                const Site* s = terms_.sites_of(i);
                size_t count = 0;
                for (size_t k = 0; k < terms_.weight_of(i); ++k) count += pauliengine_detail::site_pauli(s[k]) == 2;
                return count & 1;
        }

        SparsePauliStringSum scaled(const Coeff& s) const {
                if (pauliengine_detail::coeff_is_zero(s)) return SparsePauliStringSum();
                return mapped_coeffs([&](size_t, const Coeff& c) { return c * s; });
        }

        template <typename F>
        SparsePauliStringSum mapped_coeffs(F f) const {
                SparsePauliStringSum out = *this;
                for (size_t i = 0; i < size(); ++i) out.terms_.coeffs[i] = f(i, terms_.coeffs[i]);
                return out;
        }

        template <typename Keep>
        SparsePauliStringSum filtered(Keep keep) const {
                Terms out;
                for (size_t i = 0; i < size(); ++i) {
                        if (keep(terms_.coeffs[i])) out.push_back(terms_.sites_of(i), terms_.weight_of(i), terms_.coeffs[i]);
                }
                return from_terms(std::move(out));
        }

        SparsePauliStringSum combine(const SparsePauliStringSum& other, bool subtract) const {
                const size_t n1 = size();
                return from_terms(pauliengine_detail::accumulate<Coeff>(std::max(max_weight(), other.max_weight()), n1 + other.size(), 1, 4, [&](size_t r, auto& emit) {
                        if (r < n1) {
                                emit(terms_.sites_of(r), terms_.weight_of(r), terms_.coeffs[r]);
                        } else {
                                const size_t j = r - n1;
                                const Coeff& c = other.terms_.coeffs[j];
                                emit(other.terms_.sites_of(j), other.terms_.weight_of(j),
                                     subtract ? pauliengine_detail::coeff_zero<Coeff>() - c : c);
                        }
                }));
        }

        /// Open-addressing index of the operators of a compacted sum.
        class OperatorIndex {
            public:
                static constexpr size_t npos = SIZE_MAX;

                explicit OperatorIndex(const Terms& t) : terms_(t) {
                        size_t cap = 16;
                        while (cap < 2 * t.size()) cap *= 2;
                        slots_.assign(cap, kEmpty);
                        hashes_.resize(t.size());
                        for (size_t i = 0; i < t.size(); ++i) {
                                const uint64_t h = pauliengine_detail::hash_sites(t.sites_of(i), t.weight_of(i));
                                hashes_[i] = h;
                                size_t pos = static_cast<size_t>(h) & (cap - 1);
                                while (slots_[pos] != kEmpty) pos = (pos + 1) & (cap - 1);
                                slots_[pos] = static_cast<uint32_t>(i);
                        }
                }

                size_t find(const Site* s, size_t n) const {
                        const uint64_t h = pauliengine_detail::hash_sites(s, n);
                        const size_t mask = slots_.size() - 1;
                        for (size_t pos = static_cast<size_t>(h) & mask; slots_[pos] != kEmpty; pos = (pos + 1) & mask) {
                                const uint32_t i = slots_[pos];
                                if (hashes_[i] == h && terms_.weight_of(i) == n && std::equal(s, s + n, terms_.sites_of(i))) return i;
                        }
                        return npos;
                }

            private:
                static constexpr uint32_t kEmpty = UINT32_MAX;
                const Terms& terms_;
                std::vector<uint32_t> slots_;
                std::vector<uint64_t> hashes_;
        };

        /// The terms of a sum by qubit: for a term of another sum, its partners are
        /// the terms sharing a qubit with it, each listed once, by ascending index.
        class QubitIndex {
            public:
                explicit QubitIndex(const Terms& b) : n_terms_(b.size()) {
                        std::vector<std::pair<uint32_t, uint32_t>> pairs;  // (qubit, term)
                        pairs.reserve(b.sites.size());
                        for (size_t j = 0; j < b.size(); ++j) {
                                const Site* s = b.sites_of(j);
                                for (size_t k = 0; k < b.weight_of(j); ++k) {
                                        pairs.emplace_back(pauliengine_detail::site_qubit(s[k]), static_cast<uint32_t>(j));
                                }
                        }
                        std::sort(pairs.begin(), pairs.end());
                        for (size_t p = 0; p < pairs.size(); ++p) {
                                if (p == 0 || pairs[p].first != pairs[p - 1].first) {
                                        qubits_.push_back(pairs[p].first);
                                        start_.push_back(static_cast<uint32_t>(p));
                                }
                                terms_.push_back(pairs[p].second);
                        }
                        start_.push_back(static_cast<uint32_t>(pairs.size()));
                }

                /// Number of (term of a, term of b) pairs sharing a qubit, counted once
                /// per shared qubit; an upper bound for the pairs visited.
                double shared_pairs(const Terms& a) const {
                        double total = 0.0;
                        for (Site s : a.sites) {
                                const size_t q = locate(pauliengine_detail::site_qubit(s));
                                if (q != SIZE_MAX) total += static_cast<double>(start_[q + 1] - start_[q]);
                        }
                        return total;
                }

                /// Calls f(j) for every term j sharing a qubit with the string (sites s,
                /// n), once each and by ascending j.
                template <typename F>
                void for_each_partner(const Site* s, size_t n, F&& f) const {
                        thread_local std::vector<uint32_t> seen;
                        thread_local std::vector<uint32_t> found;
                        thread_local uint32_t stamp = 0;
                        if (seen.size() < n_terms_) {
                                seen.assign(n_terms_, 0);
                                stamp = 0;
                        }
                        if (++stamp == 0) {
                                std::fill(seen.begin(), seen.end(), 0);
                                stamp = 1;
                        }
                        found.clear();
                        for (size_t k = 0; k < n; ++k) {
                                const size_t q = locate(pauliengine_detail::site_qubit(s[k]));
                                if (q == SIZE_MAX) continue;
                                for (uint32_t p = start_[q]; p < start_[q + 1]; ++p) {
                                        const uint32_t j = terms_[p];
                                        if (seen[j] != stamp) {
                                                seen[j] = stamp;
                                                found.push_back(j);
                                        }
                                }
                        }
                        std::sort(found.begin(), found.end());
                        for (uint32_t j : found) f(j);
                }

            private:
                size_t n_terms_;
                std::vector<uint32_t> qubits_;
                std::vector<uint32_t> start_;
                std::vector<uint32_t> terms_;

                size_t locate(uint32_t qubit) const {
                        const auto it = std::lower_bound(qubits_.begin(), qubits_.end(), qubit);
                        if (it == qubits_.end() || *it != qubit) return SIZE_MAX;
                        return static_cast<size_t>(it - qubits_.begin());
                }
        };
};
