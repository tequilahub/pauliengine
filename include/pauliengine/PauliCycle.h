#pragma once

#include <algorithm>
#include <bit>
#include <complex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/PauliCycleSum.h"


// Represents a Pauli string (with a complex coefficient) together with all of
// its cyclic rotations on `n_qubits` qubits. For example, n = 3 and the base
// string XIY represents { XIY, YXI, IYX }
// https://arxiv.org/pdf/2604.16701 
class PauliCycle {
        public:
                using Coeff = std::complex<double>;

                int n_qubits;
                // Base holds one representative of the rotation orbit: exactly the
                // Pauli string the cycle was constructed from. Construction does NOT
                // canonicalize. operator== and hash() canonicalize on demand (via a
                // cached representative), so two cycles from the same rotation orbit
                // compare equal without an explicit canonicalize() call; the latter
                // just rewrites Base itself to the canonical representative.
                PauliString<Coeff> Base;

                PauliCycle(int qubits, const std::unordered_map<int, std::string>& data, Coeff coeff = Coeff(1.0, 0.0))
                        : n_qubits(qubits), Base(coeff, data) {}

                PauliCycle(int qubits, const std::string& pauli_string, Coeff coeff = Coeff(1.0, 0.0))
                        : n_qubits(qubits) {

                        Base = PauliString<Coeff>(coeff, pauli_string);
                }

                PauliCycle(int qubits, const PauliString<Coeff>& base)
                        : n_qubits(qubits), Base(base) {}

                // Rewrite Base to the canonical representative of the orbit (the
                // rotation with the largest x value, ties broken by the largest y
                // value). Optional and idempotent. operator==/hash canonicalize on
                // their own, so this is only needed when you want the stored Base
                // (and hence rotate()/repr) to be the canonical representative.
                void canonicalize() {
                        if (n_qubits <= 0) {
                                return;
                        }
                        ensure_canonical();
                        Base = PauliString<Coeff>(canon_x_, canon_y_, Base.coeff);
                }


                PauliString<Coeff> rotate(int k) const {
                        const int n = n_qubits;
                        if (n <= 0) {
                                return Base;
                        }
                        int shift = k % n;
                        if (shift < 0) {
                                shift += n;
                        }

                        const size_t n_words = (static_cast<size_t>(n) + BITS_IN_INTEGER - 1) / BITS_IN_INTEGER;
                        WordVec new_x(n_words);
                        WordVec new_y(n_words);

                        for (int i = 0; i < n; ++i) {
                                const int j = (i + shift) % n;
                                if (bit_at(Base.x, i)) {
                                        set_bit(new_x, j);
                                }
                                if (bit_at(Base.y, i)) {
                                        set_bit(new_y, j);
                                }
                        }
                        return PauliString<Coeff>(std::move(new_x), std::move(new_y), Base.coeff);
                }

                std::vector<PauliString<Coeff>> rotations() const {
                        std::vector<PauliString<Coeff>> out;
                        out.reserve(n_qubits > 0 ? n_qubits : 0);
                        for (int k = 0; k < n_qubits; ++k) {
                                out.push_back(rotate(k));
                        }
                        return out;
                }

                QubitHamiltonian<Coeff> to_qubit_hamiltonian() const {
                        return QubitHamiltonian<Coeff>(rotations());
                }


                // Orbit-level equality: two cycles are equal iff they describe the
                // same rotation orbit. Both sides are canonicalized on demand (cached),
                // so no explicit canonicalize() call is required. The coefficient is
                // ignored (structural equality of the cycle).
                bool operator==(const PauliCycle& other) const {
                        if (n_qubits != other.n_qubits) {
                                return false;
                        }
                        ensure_canonical();
                        other.ensure_canonical();
                        const size_t nwords = std::max({n_words(), canon_x_.size(),
                                canon_y_.size(), other.canon_x_.size(), other.canon_y_.size()});
                        return cmp_words(canon_x_, other.canon_x_, nwords) == 0
                                && cmp_words(canon_y_, other.canon_y_, nwords) == 0;
                }

                bool operator!=(const PauliCycle& other) const {
                        return !(*this == other);
                }

                // Hash of the canonical representative, consistent with operator==.
                size_t hash() const {
                        ensure_canonical();
                        size_t h = std::hash<int>{}(n_qubits);
                        mix_words(h, canon_x_);
                        mix_words(h, canon_y_);
                        return h;
                }

                PauliCycleSum multiply(const PauliCycle& other) const {
                        std::vector<PauliCycle> result;
                        const int n = this->n_qubits;
                        if (n <= 0) {
                                return PauliCycleSum(1.0, result);
                        }

                        result.reserve(n);
                        for (int b = 0; b < n; ++b) {
                                PauliString<Coeff> term = this->Base * other.rotate(b);
                                PauliCycle cycle(n, term);

                                // Merge cycles from the same rotation orbit; operator==
                                // canonicalizes on demand, so products that are rotations
                                // of each other are recognised and their coefficients
                                // summed.
                                bool merged = false;
                                for (auto& existing : result) {
                                        if (existing == cycle) {
                                                existing.Base.coeff += cycle.Base.coeff;
                                                merged = true;
                                                break;
                                        }
                                }
                                if (!merged) {
                                        result.push_back(std::move(cycle));
                                }
                        }

                        // Drop terms whose merged coefficient cancelled to zero.
                        result.erase(
                                std::remove_if(result.begin(), result.end(),
                                        [](const PauliCycle& c) { return c.Base.coeff == Coeff(0.0, 0.0); }),
                                result.end());

                        return PauliCycleSum(1.0, result);
                }

                PauliCycleSum commutator(const PauliCycle& other) const {
                        std::vector<PauliCycle> result;
                        const int n = this->n_qubits;
                        if (n <= 0) {
                                return PauliCycleSum(1.0 / (n != 0 ? n : 1), result);
                        }

                        const std::vector<int> supp_p = support(this->Base);
                        const std::vector<int> supp_pp = support(other.Base);

                        // Candidate shifts where the supports can overlap 
                        std::unordered_set<int> shifts;
                        shifts.reserve(supp_p.size() * supp_pp.size());
                        for (int qp : supp_p) {
                                for (int qpp : supp_pp) {
                                        shifts.insert(((qp - qpp) % n + n) % n);
                                }
                        }

                        result.reserve(shifts.size());
                        for (int k : shifts) {
                                PauliString<Coeff> term = this->Base.commutator(other.rotate(k));

                                if (term.coeff != Coeff(0.0, 0.0)) {
                                        result.push_back(PauliCycle(n, term));
                                }
                        }
                        return PauliCycleSum(1.0 / n, result);
                }

                // Coefficient of a single target cycle in [this, other], i.e. the sum
                // of the coefficients of the cycles in the full commutator whose
                // rotation orbit equals `target`. Uses the same support-difference
                // shift pruning as commutator() but, instead of constructing and
                // merging every output cycle, only accumulates the shifts whose
                // product belongs to the target's cyclic class (PauliEngine Features).
                Coeff targeted_commutator(const PauliCycle& other, const PauliCycle& target) const {
                        const int n = this->n_qubits;
                        Coeff total(0.0, 0.0);
                        if (other.n_qubits != n || target.n_qubits != n || n <= 0) {
                                return total;
                        }

                        const std::vector<int> supp_p = support(this->Base);
                        const std::vector<int> supp_pp = support(other.Base);
                        std::unordered_set<int> shifts;
                        shifts.reserve(supp_p.size() * supp_pp.size());
                        for (int qp : supp_p) {
                                for (int qpp : supp_pp) {
                                        shifts.insert(((qp - qpp) % n + n) % n);
                                }
                        }

                        for (int k : shifts) {
                                PauliString<Coeff> term = this->Base.commutator(other.rotate(k));
                                if (term.coeff != Coeff(0.0, 0.0)) {
                                        // Belongs to the target's cyclic class? (operator== canonicalizes.)
                                        if (PauliCycle(n, term) == target) {
                                                total += term.coeff;
                                        }
                                }
                        }
                        return total;
                }


        private:
                // Lazily-computed canonical representative (max-value rotation), cached
                // so repeated ==/hash on the same cycle don't re-scan all n rotations.
                // Depends only on Base's operator part (x, y), which is set at
                // construction and never mutated afterwards (multiply only touches the
                // coefficient); canonicalize() keeps it consistent.
                mutable bool canon_valid_ = false;
                mutable WordVec canon_x_;
                mutable WordVec canon_y_;

                void ensure_canonical() const {
                        if (canon_valid_) {
                                return;
                        }
                        if (n_qubits <= 0) {
                                canon_x_ = Base.x;
                                canon_y_ = Base.y;
                                canon_valid_ = true;
                                return;
                        }
                        PauliString<Coeff> best = rotate(0);
                        for (int k = 1; k < n_qubits; ++k) {
                                PauliString<Coeff> cand = rotate(k);
                                if (rotation_greater(cand, best)) {
                                        best = std::move(cand);
                                }
                        }
                        canon_x_ = std::move(best.x);
                        canon_y_ = std::move(best.y);
                        canon_valid_ = true;
                }

                size_t n_words() const {
                        return n_qubits > 0
                                ? (static_cast<size_t>(n_qubits) + BITS_IN_INTEGER - 1) / BITS_IN_INTEGER
                                : 0;
                }


                bool rotation_greater(const PauliString<Coeff>& a, const PauliString<Coeff>& b) const {
                        const size_t nwords = n_words();
                        const int cx = cmp_words(a.x, b.x, nwords);
                        if (cx != 0) {
                                return cx > 0;
                        }
                        return cmp_words(a.y, b.y, nwords) > 0;
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

inline bool PauliCycleSum::contains(const PauliCycle& cycle) const {
        for (const auto& c : data) {
                if (c == cycle) {
                        return true;
                }
        }
        return false;
}

inline PauliCycleSum PauliCycleSum::commutator(const PauliCycleSum& other) const {
        // [ sum_i c_i , sum_j d_j ] = sum_{i,j} [c_i, d_j]. Each cycle carries its
        // own coefficient, so the pairwise results are simply concatenated.
        std::vector<PauliCycle> combined;
        for (const auto& ci : data) {
                for (const auto& dj : other.data) {
                        PauliCycleSum part = ci.commutator(dj);
                        for (auto& c : part.data) {
                                combined.push_back(std::move(c));
                        }
                }
        }
        return PauliCycleSum(1.0, combined);
}

inline PauliCycleSum PauliCycleSum::multiply(const PauliCycleSum& other) const {
        // ( sum_i c_i )( sum_j d_j ) = sum_{i,j} c_i d_j. PauliCycle::multiply
        // canonicalizes its output cycles, so equal ones can be merged afterwards.
        std::vector<PauliCycle> combined;
        for (const auto& ci : data) {
                for (const auto& dj : other.data) {
                        PauliCycleSum part = ci.multiply(dj);
                        for (auto& c : part.data) {
                                combined.push_back(std::move(c));
                        }
                }
        }

        std::vector<PauliCycle> merged;
        for (auto& c : combined) {
                bool done = false;
                for (auto& e : merged) {
                        if (e == c) {
                                e.Base.coeff += c.Base.coeff;
                                done = true;
                                break;
                        }
                }
                if (!done) {
                        merged.push_back(std::move(c));
                }
        }
        merged.erase(
                std::remove_if(merged.begin(), merged.end(),
                        [](const PauliCycle& c) { return c.Base.coeff == PauliCycle::Coeff(0.0, 0.0); }),
                merged.end());

        return PauliCycleSum(1.0, merged);
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
