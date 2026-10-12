#pragma once

#include <algorithm>
#include <bit>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/PauliCycleSum.h"


/// Represents a Pauli string (with a complex coefficient) together with all of
/// its cyclic rotations on `n_qubits` qubits, i.e. the operator
///
///     C(P) = coeff * sum_{k=0}^{n-1} R^k P        (R = cyclic shift by one qubit)
///
/// For example, n = 3 and the base string XIY stands for XIY + YXI + IYX.
/// The sum always runs over all n rotations, so a base with a shorter cyclic
/// period d (n % d == 0) counts each distinct rotation n/d times: the cycle of
/// ZZZ on 3 qubits is 3*ZZZ, and that of XIXI on 4 qubits is 2*(XIXI + IXIX).
/// This is what makes C(P) C(Q) = sum_k C(P R^k Q) hold for every base.
/// https://arxiv.org/pdf/2604.16701
class PauliCycle {
        public:
                using Coeff = std::complex<double>;

                /// Fixed at construction: the cached unified form below depends on
                /// it, so it must not change afterwards (exposed read-only to Python).
                int n_qubits;
                /// Base holds one representative of the rotation orbit: exactly the
                /// Pauli string the cycle was constructed from. Construction does NOT
                /// unify. operator== and hash() compare unified forms on demand (via a
                /// cached representative), so two cycles from the same rotation orbit
                /// compare equal without an explicit unify() call; the latter
                /// just rewrites Base itself to the unified representative.
                PauliString<Coeff> Base;

                /// Cycle of the base {qubit: "X"/"Y"/"Z"} on `qubits` qubits with coefficient `coeff`.
                /// Throws std::invalid_argument if qubits < 1 or the base acts on a qubit >= qubits.
                PauliCycle(int qubits, const std::unordered_map<int, std::string>& data, Coeff coeff = Coeff(1.0, 0.0))
                        : n_qubits(qubits), Base(coeff, data) {
                        validate();
                }

                /// Cycle of a dense base string such as "XIY" (character i acts on qubit i).
                PauliCycle(int qubits, const std::string& pauli_string, Coeff coeff = Coeff(1.0, 0.0))
                        : n_qubits(qubits) {
                        Base = PauliString<Coeff>(coeff, pauli_string);
                        validate();
                }

                /// Cycle of an existing PauliString; its coefficient becomes the cycle's.
                PauliCycle(int qubits, const PauliString<Coeff>& base)
                        : n_qubits(qubits), Base(base) {
                        validate();
                }

                /// Rewrite Base to the unified representative of the orbit (see
                /// ensure_unified()). Optional and idempotent. operator==/hash are
                /// rotation-invariant on their own, so this is only needed when you
                /// want the stored Base (and hence rotate()/repr) to be unified.
                void unify() {
                        ensure_unified();
                        Base = PauliString<Coeff>(unified_x_, unified_y_, Base.coeff);
                }


                /// Cyclic shift by k qubits: new[(i + k) mod n] = old[i]. Word-parallel,
                /// O(n / 64).
                PauliString<Coeff> rotate(int k) const {
                        const int n = n_qubits;
                        int shift = k % n;
                        if (shift < 0) {
                                shift += n;
                        }
                        return PauliString<Coeff>(rotate_words(Base.x, n, shift),
                                rotate_words(Base.y, n, shift), Base.coeff);
                }

                /// All n rotations rotate(0), ..., rotate(n - 1) of the stored base.
                std::vector<PauliString<Coeff>> rotations() const {
                        std::vector<PauliString<Coeff>> out;
                        out.reserve(n_qubits > 0 ? n_qubits : 0);
                        for (int k = 0; k < n_qubits; ++k) {
                                out.push_back(rotate(k));
                        }
                        return out;
                }

                /// The cycle expanded into an explicit QubitHamiltonian (sum of all rotations).
                QubitHamiltonian<Coeff> to_qubit_hamiltonian() const {
                        return QubitHamiltonian<Coeff>(rotations());
                }


                /// Orbit-level equality: two cycles are equal iff they describe the
                /// same rotation orbit. Both sides are unified on demand (cached),
                /// so no explicit unify() call is required. The coefficient is
                /// ignored (structural equality of the cycle).
                bool operator==(const PauliCycle& other) const {
                        if (n_qubits != other.n_qubits) {
                                return false;
                        }
                        ensure_unified();
                        other.ensure_unified();
                        const size_t nwords = std::max({n_words(), unified_x_.size(),
                                unified_y_.size(), other.unified_x_.size(), other.unified_y_.size()});
                        return cmp_words(unified_x_, other.unified_x_, nwords) == 0
                                && cmp_words(unified_y_, other.unified_y_, nwords) == 0;
                }

                /// Negation of operator==.
                bool operator!=(const PauliCycle& other) const {
                        return !(*this == other);
                }

                /// Hash of the unified representative, consistent with operator==.
                size_t hash() const {
                        ensure_unified();
                        size_t h = std::hash<int>{}(n_qubits);
                        mix_words(h, unified_x_);
                        mix_words(h, unified_y_);
                        return h;
                }

                /// Invertible, coefficient-independent key of the cycle, consistent
                /// with operator==: key = string_key * 2^32 + n_qubits, where
                /// string_key is PauliString::get_hash() of the unified representative
                /// (so all rotations share it). Little-endian 64-bit limbs.
                std::vector<uint64_t> get_hash() const {
                        ensure_unified();
                        const std::vector<uint64_t> s =
                                PauliString<Coeff>(unified_x_, unified_y_, Coeff(1.0, 0.0)).get_hash();
                        std::vector<uint64_t> out(s.size() + 1);
                        out[0] = static_cast<uint64_t>(static_cast<uint32_t>(n_qubits));
                        for (size_t i = 0; i < s.size(); ++i) {
                                out[i] |= s[i] << 32;
                                out[i + 1] = s[i] >> 32;
                        }
                        while (out.size() > 1 && out.back() == 0) {
                                out.pop_back();
                        }
                        return out;
                }

                /// Inverse of get_hash(): the cycle with that key and the given
                /// coefficient (its base is the unified representative).
                static PauliCycle from_hash(const std::vector<uint64_t>& limbs, Coeff coeff = Coeff(1.0, 0.0)) {
                        const uint64_t n = limbs.empty() ? 0 : (limbs[0] & 0xFFFFFFFFULL);
                        if (n < 1 || n > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
                                throw std::invalid_argument("PauliCycle.from_hash: invalid hash (n_qubits = "
                                        + std::to_string(n) + ")");
                        }
                        std::vector<uint64_t> s(limbs.size());
                        for (size_t i = 0; i < limbs.size(); ++i) {
                                s[i] = (limbs[i] >> 32) | (i + 1 < limbs.size() ? limbs[i + 1] << 32 : 0ULL);
                        }
                        // The constructor rejects a string acting beyond n_qubits.
                        return PauliCycle(static_cast<int>(n), PauliString<Coeff>::from_hash(s, coeff));
                }

                /// Product C(P) C(Q) = sum_{b=0}^{n-1} C(P R^b Q), as a PauliCycleSum.
                /// Products that are rotations of each other are merged.
                PauliCycleSum multiply(const PauliCycle& other) const {
                        require_same_n_qubits(other, "multiply");
                        const int n = this->n_qubits;
                        std::vector<PauliCycle> result;
                        result.reserve(n);
                        for (int b = 0; b < n; ++b) {
                                result.emplace_back(n, this->Base * other.rotate(b));
                        }
                        return PauliCycleSum(std::move(result));
                }

                /// [C(P), C(Q)] = sum_k C([P, R^k Q]) (Eq. 40). Only shifts k at which
                /// the supports overlap can contribute (Corollary 3); different shifts
                /// landing in the same rotation orbit are merged by PauliCycleSum.
                PauliCycleSum commutator(const PauliCycle& other) const {
                        require_same_n_qubits(other, "commutator");
                        std::vector<PauliCycle> result;
                        const int n = this->n_qubits;

                        const std::vector<int> shifts = overlap_shifts(other);
                        result.reserve(shifts.size());
                        for (int k : shifts) {
                                PauliString<Coeff> term = this->Base.commutator(other.rotate(k));

                                if (term.coeff != Coeff(0.0, 0.0)) {
                                        result.emplace_back(n, term);
                                }
                        }
                        return PauliCycleSum(std::move(result));
                }

                /// Coefficient of a single target cycle in [this, other], i.e. the sum
                /// of the coefficients of the cycles in the full commutator whose
                /// rotation orbit equals `target`. Uses the same support-difference
                /// shift pruning as commutator() but, instead of constructing and
                /// merging every output cycle, only accumulates the shifts whose
                /// product belongs to the target's cyclic class (PauliEngine Features).
                Coeff targeted_commutator(const PauliCycle& other, const PauliCycle& target) const {
                        require_same_n_qubits(other, "targeted_commutator");
                        require_same_n_qubits(target, "targeted_commutator");
                        const int n = this->n_qubits;
                        Coeff total(0.0, 0.0);
                        // Rotations preserve the weight, so a product of a different
                        // weight can never lie in the target's orbit; checking that
                        // first skips building and unifying most candidates.
                        const size_t target_weight = target.Base.size();

                        for (int k : overlap_shifts(other)) {
                                PauliString<Coeff> term = this->Base.commutator(other.rotate(k));
                                if (term.coeff != Coeff(0.0, 0.0) && term.size() == target_weight) {
                                        // Belongs to the target's cyclic class? (operator== unifies.)
                                        if (PauliCycle(n, term) == target) {
                                                total += term.coeff;
                                        }
                                }
                        }
                        return total;
                }


        private:
                /// A cycle needs at least one qubit, and its base may only act on
                /// qubits 0..n_qubits-1: anything beyond would be dropped by rotate().
                void validate() const {
                        if (n_qubits < 1) {
                                throw std::invalid_argument("PauliCycle: n_qubits must be >= 1, got "
                                        + std::to_string(n_qubits));
                        }
                        const std::vector<int> supp = support(Base);
                        if (!supp.empty() && supp.back() >= n_qubits) {
                                throw std::invalid_argument("PauliCycle: base acts on qubit "
                                        + std::to_string(supp.back()) + ", but n_qubits is "
                                        + std::to_string(n_qubits));
                        }
                }

                void require_same_n_qubits(const PauliCycle& other, const char* op) const {
                        if (other.n_qubits != n_qubits) {
                                throw std::invalid_argument(std::string("PauliCycle.") + op
                                        + ": n_qubits mismatch (" + std::to_string(n_qubits)
                                        + " vs " + std::to_string(other.n_qubits) + ")");
                        }
                }

                /// Shifts k at which the supports of Base and R^k other.Base overlap
                /// (Corollary 3), deduplicated via a bitmap: O(|supp| * |supp'| + n).
                std::vector<int> overlap_shifts(const PauliCycle& other) const {
                        const int n = n_qubits;
                        const std::vector<int> supp_p = support(this->Base);
                        const std::vector<int> supp_pp = support(other.Base);
                        std::vector<char> seen(static_cast<size_t>(n), 0);
                        std::vector<int> shifts;
                        shifts.reserve(std::min(supp_p.size() * supp_pp.size(), static_cast<size_t>(n)));
                        for (int qp : supp_p) {
                                for (int qpp : supp_pp) {
                                        int k = qp - qpp;
                                        if (k < 0) {
                                                k += n;
                                        }
                                        if (!seen[k]) {
                                                seen[k] = 1;
                                                shifts.push_back(k);
                                        }
                                }
                        }
                        return shifts;
                }

                /// Lazily-computed unified representative, cached so repeated ==/hash
                /// on the same cycle don't recompute it. Depends only on Base's operator
                /// part (x, y), which is set at construction and never mutated
                /// afterwards (multiply only touches the coefficient); unify() keeps it
                /// consistent.
                mutable bool unified_valid_ = false;
                mutable WordVec unified_x_;
                mutable WordVec unified_y_;

                /// A rotation orbit is determined by the cyclic sequence of
                /// (gap to the next non-identity site, Pauli) pairs along the support:
                /// two strings are rotations of each other iff these sequences are
                /// cyclic shifts of each other. The unified representative places
                /// the lexicographically least rotation of that sequence at qubit 0.
                /// Cost O(n / 64 + w) for weight w, instead of comparing all n
                /// rotations (O(n^2 / 64)).
                void ensure_unified() const {
                        if (unified_valid_) {
                                return;
                        }
                        const int n = n_qubits;
                        unified_x_ = WordVec(n_words());
                        unified_y_ = WordVec(n_words());

                        const std::vector<int> pos = support(Base);
                        const size_t w = pos.size();
                        if (w > 0) {
                                // item = gap << 2 | pauli, with pauli = x_bit | y_bit << 1
                                std::vector<uint64_t> items(w);
                                for (size_t i = 0; i < w; ++i) {
                                        const int gap = (i + 1 < w ? pos[i + 1] : pos[0] + n) - pos[i];
                                        const uint64_t pauli = (bit_at(Base.x, pos[i]) ? 1ULL : 0ULL)
                                                | (bit_at(Base.y, pos[i]) ? 2ULL : 0ULL);
                                        items[i] = (static_cast<uint64_t>(gap) << 2) | pauli;
                                }
                                const size_t start = least_rotation(items);
                                int q = 0;
                                for (size_t t = 0; t < w; ++t) {
                                        const uint64_t item = items[(start + t) % w];
                                        if (item & 1ULL) {
                                                set_bit(unified_x_, q);
                                        }
                                        if (item & 2ULL) {
                                                set_bit(unified_y_, q);
                                        }
                                        q += static_cast<int>(item >> 2);
                                }
                        }
                        unified_valid_ = true;
                }

                /// Start index of the lexicographically least rotation of s, O(|s|)
                /// (two-pointer minimum-expression algorithm).
                static size_t least_rotation(const std::vector<uint64_t>& s) {
                        const size_t m = s.size();
                        size_t i = 0, j = 1, k = 0;
                        while (i < m && j < m && k < m) {
                                const uint64_t a = s[(i + k) % m];
                                const uint64_t b = s[(j + k) % m];
                                if (a == b) {
                                        ++k;
                                        continue;
                                }
                                if (a > b) {
                                        i += k + 1;
                                } else {
                                        j += k + 1;
                                }
                                if (i == j) {
                                        ++j;
                                }
                                k = 0;
                        }
                        return std::min(i, j);
                }

                /// Rotate the n-bit vector `in` by `shift` (0 <= shift < n) bits:
                /// out = ((in << shift) | (in >> (n - shift))) restricted to n bits.
                static WordVec rotate_words(const WordVec& in, int n, int shift) {
                        const size_t nw = (static_cast<size_t>(n) + BITS_IN_INTEGER - 1) / BITS_IN_INTEGER;
                        WordVec out(nw);
                        const auto get = [&in](long long i) -> uint64_t {
                                return (i >= 0 && static_cast<size_t>(i) < in.size()) ? in[static_cast<size_t>(i)] : 0ULL;
                        };
                        // Left shift by `shift` and right shift by `n - shift`, OR-ed.
                        const long long ls_w = shift / static_cast<int>(BITS_IN_INTEGER);
                        const int ls_b = shift % static_cast<int>(BITS_IN_INTEGER);
                        const long long rs_w = (n - shift) / static_cast<int>(BITS_IN_INTEGER);
                        const int rs_b = (n - shift) % static_cast<int>(BITS_IN_INTEGER);
                        for (size_t i = 0; i < nw; ++i) {
                                const long long ii = static_cast<long long>(i);
                                uint64_t v = get(ii - ls_w) << ls_b;
                                if (ls_b != 0) {
                                        v |= get(ii - ls_w - 1) >> (BITS_IN_INTEGER - ls_b);
                                }
                                if (shift != 0) {
                                        v |= get(ii + rs_w) >> rs_b;
                                        if (rs_b != 0) {
                                                v |= get(ii + rs_w + 1) << (BITS_IN_INTEGER - rs_b);
                                        }
                                }
                                out[i] = v;
                        }
                        const int tail = n % static_cast<int>(BITS_IN_INTEGER);
                        if (tail != 0) {
                                out[nw - 1] &= (1ULL << tail) - 1ULL;
                        }
                        return out;
                }

                size_t n_words() const {
                        return (static_cast<size_t>(n_qubits) + BITS_IN_INTEGER - 1) / BITS_IN_INTEGER;
                }

                static int cmp_words(const WordVec& a, const WordVec& b, size_t nwords) {
                        for (size_t i = nwords; i-- > 0; ) {
                                const uint64_t av = (i < a.size()) ? a[i] : 0ULL;
                                const uint64_t bv = (i < b.size()) ? b[i] : 0ULL;
                                if (av != bv) {
                                        return (av < bv) ? -1 : 1;
                                }
                        }
                        return 0;
                }


                static void mix_words(size_t& h, const WordVec& v) {
                        for (size_t i = 0; i < v.size(); ++i) {
                                h ^= std::hash<uint64_t>{}(v[i]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                        }
                }

                static std::vector<int> support(const PauliString<Coeff>& ps) {
                        std::vector<int> out;
                        for (size_t w = 0; w < ps.x.size(); ++w) {
                                uint64_t bits = ps.x[w] | ps.y[w];
                                while (bits) {
                                        const int bit = std::countr_zero(bits);
                                        out.push_back(static_cast<int>(w * BITS_IN_INTEGER + bit));
                                        bits &= bits - 1;
                                }
                        }
                        return out;
                }

                static bool bit_at(const WordVec& v, int pos) {
                        const size_t word = static_cast<size_t>(pos) / BITS_IN_INTEGER;
                        if (word >= v.size()) {
                                return false;
                        }
                        return (v[word] >> (static_cast<size_t>(pos) % BITS_IN_INTEGER)) & 1ULL;
                }

                static void set_bit(WordVec& v, int pos) {
                        const size_t word = static_cast<size_t>(pos) / BITS_IN_INTEGER;
                        v[word] |= (1ULL << (static_cast<size_t>(pos) % BITS_IN_INTEGER));
                }
};


inline QubitHamiltonian<PauliCycleSum::Coeff> PauliCycleSum::to_qubit_hamiltonian() const {
        std::vector<PauliString<Coeff>> terms;
        for (const auto& cycle : data) {
                std::vector<PauliString<Coeff>> rots = cycle.rotations();
                terms.insert(terms.end(),
                        std::make_move_iterator(rots.begin()),
                        std::make_move_iterator(rots.end()));
        }
        return QubitHamiltonian<Coeff>(std::move(terms));
}

inline PauliCycleSum::PauliCycleSum(std::vector<PauliCycle> cycles) {
        // Hash-join on the rotation orbit: O(1) per summand instead of a linear
        // scan. The map keys point into `data`; reserving up front guarantees
        // they stay valid (data never grows beyond cycles.size()).
        struct OrbitHash {
                size_t operator()(const PauliCycle* c) const { return c->hash(); }
        };
        struct OrbitEqual {
                bool operator()(const PauliCycle* a, const PauliCycle* b) const { return *a == *b; }
        };
        std::unordered_map<const PauliCycle*, size_t, OrbitHash, OrbitEqual> index;
        index.reserve(cycles.size());
        data.reserve(cycles.size());

        for (auto& c : cycles) {
                if (!data.empty() && c.n_qubits != data.front().n_qubits) {
                        throw std::invalid_argument("PauliCycleSum: all cycles must have the same n_qubits ("
                                + std::to_string(data.front().n_qubits) + " vs "
                                + std::to_string(c.n_qubits) + ")");
                }
                const auto it = index.find(&c);
                if (it != index.end()) {
                        data[it->second].Base.coeff += c.Base.coeff;
                        continue;
                }
                data.push_back(std::move(c));
                index.emplace(&data.back(), data.size() - 1);
        }

        // Drop summands whose merged coefficient cancelled to zero.
        data.erase(
                std::remove_if(data.begin(), data.end(),
                        [](const PauliCycle& c) { return c.Base.coeff == Coeff(0.0, 0.0); }),
                data.end());
}

inline PauliCycleSum::Coeff PauliCycleSum::contains(const PauliCycle& cycle) const {
        // The sum is merged, so at most one summand lies in the cycle's orbit.
        for (const auto& c : data) {
                if (c == cycle) {
                        return c.Base.coeff;
                }
        }
        return Coeff(0.0, 0.0);
}

inline PauliCycleSum PauliCycleSum::commutator(const PauliCycleSum& other) const {
        // [ sum_i c_i , sum_j d_j ] = sum_{i,j} [c_i, d_j]. Each cycle carries its
        // own coefficient, so the pairwise results are collected and merged.
        std::vector<PauliCycle> combined;
        for (const auto& ci : data) {
                for (const auto& dj : other.data) {
                        PauliCycleSum part = ci.commutator(dj);
                        combined.insert(combined.end(),
                                std::make_move_iterator(part.data.begin()),
                                std::make_move_iterator(part.data.end()));
                }
        }
        return PauliCycleSum(std::move(combined));
}

inline PauliCycleSum PauliCycleSum::multiply(const PauliCycleSum& other) const {
        // ( sum_i c_i )( sum_j d_j ) = sum_{i,j} c_i d_j, merged across all pairs.
        std::vector<PauliCycle> combined;
        for (const auto& ci : data) {
                for (const auto& dj : other.data) {
                        PauliCycleSum part = ci.multiply(dj);
                        combined.insert(combined.end(),
                                std::make_move_iterator(part.data.begin()),
                                std::make_move_iterator(part.data.end()));
                }
        }
        return PauliCycleSum(std::move(combined));
}

inline PauliCycleSum::Coeff PauliCycleSum::targeted_commutator(const PauliCycleSum& other,
                                                               const PauliCycle& target) const {
        // Cycles carry their own coefficient, so the pairwise targeted coefficients
        // are simply summed (no extra scaling).
        Coeff total(0.0, 0.0);
        for (const auto& ci : data) {
                for (const auto& dj : other.data) {
                        total += ci.targeted_commutator(dj, target);
                }
        }
        return total;
}
