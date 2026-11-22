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

#include "pauliengine/detail/Coefficients.h"
#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/PauliCycleSum.h"


/// Represents a Pauli string with a coefficient of type Coeff (std::complex<double>
/// or SymEngine::Expression) together with all of its cyclic rotations on
/// `n_qubits` qubits, i.e. the operator
///
///     C(P) = coeff * sum_{k=0}^{n-1} R^k P        (R = cyclic shift by one qubit)
///
/// For example, n = 3 and the base string XIY stands for XIY + YXI + IYX.
/// The sum always runs over all n rotations, so a base with a shorter cyclic
/// period d (n % d == 0) counts each distinct rotation n/d times: the cycle of
/// ZZZ on 3 qubits is 3*ZZZ, and that of XIXI on 4 qubits is 2*(XIXI + IXIX).
/// This is what makes C(P) C(Q) = sum_k C(P R^k Q) hold for every base.
/// https://arxiv.org/pdf/2604.16701
template <typename Coeff>
class PauliCycle {
        template <typename> friend class PauliCycle;
        template <typename> friend class PauliCycleSum;

        public:
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
                PauliCycle(int qubits, const std::unordered_map<int, std::string>& data,
                           Coeff coeff = pauliengine_detail::coeff_one<Coeff>())
                        : n_qubits(qubits), Base(coeff, data) {
                        validate();
                }

                /// Cycle of a dense base string such as "XIY" (character i acts on qubit i).
                PauliCycle(int qubits, const std::string& pauli_string, Coeff coeff = pauliengine_detail::coeff_one<Coeff>())
                        : n_qubits(qubits) {
                        Base = PauliString<Coeff>(coeff, pauli_string);
                        validate();
                }

                /// Cycle of an existing PauliString; its coefficient becomes the cycle's.
                PauliCycle(int qubits, const PauliString<Coeff>& base)
                        : n_qubits(qubits), Base(base) {
                        validate();
                }

        private:
                /// A base known to act on qubits 0..qubits-1 only, such as a product of
                /// two rotations of valid bases: no validation.
                struct trusted_t {};
                PauliCycle(int qubits, PauliString<Coeff>&& base, trusted_t)
                        : n_qubits(qubits), Base(std::move(base)) {}

        public:

                /// Rewrite Base to the unified representative of the orbit (see
                /// ensure_unified()). Optional and idempotent. operator==/hash are
                /// rotation-invariant on their own, so this is only needed when you
                /// want the stored Base (and hence rotate()/repr) to be unified.
                void unify() {
                        if (base_unified_) {
                                return;
                        }
                        ensure_unified();
                        Base = PauliString<Coeff>(unified_x_, unified_y_, Base.coeff);
                        // the cache is the base now
                        base_unified_ = true;
                        unified_valid_ = false;
                        unified_x_ = WordVec();
                        unified_y_ = WordVec();
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
                /// so no explicit unify() call is required. The coefficient (and its
                /// type) is ignored (structural equality of the cycle).
                template <typename C>
                bool operator==(const PauliCycle<C>& other) const {
                        if (n_qubits != other.n_qubits) {
                                return false;
                        }
                        ensure_unified();
                        other.ensure_unified();
                        const size_t nwords = std::max({n_words(), ux().size(),
                                uy().size(), other.ux().size(), other.uy().size()});
                        return cmp_words(ux(), other.ux(), nwords) == 0
                                && cmp_words(uy(), other.uy(), nwords) == 0;
                }

                /// Period d of the base: the smallest d > 0 with R^d P = P, i.e. the
                /// number of distinct rotations (d divides n_qubits). The sum over all
                /// n rotations counts each distinct one n/d times. Computed in
                /// O(n/64 + w) from the cyclic (gap, Pauli) sequence: if that sequence
                /// consists of r repetitions of a block, then d = n / r.
                int period() const {
                        const std::vector<int> pos = support(Base);
                        const size_t w = pos.size();
                        if (w == 0) {
                                return 1;                       // the identity is invariant
                        }
                        std::vector<uint64_t> items(w);
                        for (size_t i = 0; i < w; ++i) {
                                const int gap = (i + 1 < w ? pos[i + 1] : pos[0] + n_qubits) - pos[i];
                                const uint64_t pauli = (bit_at(Base.x, pos[i]) ? 1ULL : 0ULL)
                                        | (bit_at(Base.y, pos[i]) ? 2ULL : 0ULL);
                                items[i] = (static_cast<uint64_t>(gap) << 2) | pauli;
                        }
                        // smallest block length via the prefix function (KMP)
                        std::vector<size_t> pi(w, 0);
                        for (size_t i = 1; i < w; ++i) {
                                size_t k = pi[i - 1];
                                while (k > 0 && items[i] != items[k]) {
                                        k = pi[k - 1];
                                }
                                if (items[i] == items[k]) {
                                        ++k;
                                }
                                pi[i] = k;
                        }
                        const size_t block = w - pi[w - 1];
                        const size_t repetitions = (w % block == 0) ? w / block : 1;
                        return n_qubits / static_cast<int>(repetitions);
                }

                /// Negation of operator==.
                template <typename C>
                bool operator!=(const PauliCycle<C>& other) const {
                        return !(*this == other);
                }

                /// Hash of the unified representative, consistent with operator==.
                size_t hash() const {
                        ensure_unified();
                        size_t h = std::hash<int>{}(n_qubits);
                        mix_words(h, ux());
                        mix_words(h, uy());
                        return h;
                }

                /// Whether rotation_signature() is cheaper than the unified form: the base
                /// has many sites for its number of words. The same for all rotations of
                /// the base.
                bool prefers_signature() const {
                        const size_t used = std::min(Base.x.size(), Base.y.size());
                        size_t weight = 0;
                        for (size_t w = 0; w < used; ++w) weight += static_cast<size_t>(std::popcount(Base.x[w] | Base.y[w]));
                        return weight > 4 * n_words() + 16;
                }

                /// Hash of quantities that do not change under rotation: n_qubits, the
                /// numbers of X, Y and Z sites and of every pair of cyclically adjacent
                /// Paulis. Equal cycles (bases that are rotations of each other) always
                /// share it; cycles with different signatures differ. O(n / 64), against
                /// O(n / 64 + w) for the unified form, so merging can skip unifying cycles
                /// whose signatures differ.
                size_t rotation_signature() const {
                        const int n = n_qubits;
                        const size_t nw = n_words();
                        thread_local WordVec masks[3];
                        thread_local WordVec shifted[3];
                        for (WordVec& m : masks) m.resize(nw);
                        const size_t used = std::min(Base.x.size(), Base.y.size());
                        uint64_t counts[3] = {0, 0, 0};
                        for (size_t w = 0; w < nw; ++w) {
                                const uint64_t x = w < used ? Base.x[w] : 0;
                                const uint64_t y = w < used ? Base.y[w] : 0;
                                masks[0][w] = x & ~y;  // X
                                masks[1][w] = y & ~x;  // Y
                                masks[2][w] = x & y;   // Z
                                for (int p = 0; p < 3; ++p) counts[p] += static_cast<uint64_t>(std::popcount(masks[p][w]));
                        }
                        uint64_t h = 0x9E3779B97F4A7C15ULL ^ static_cast<uint64_t>(n);
                        const auto mix = [&h](uint64_t v) {
                                h ^= v + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
                                h *= 0xBF58476D1CE4E5B9ULL;
                                h ^= h >> 31;
                        };
                        for (uint64_t c : counts) mix(c);
                        if (n > 1) {
                                for (int p = 0; p < 3; ++p) rotate_words_into(masks[p], n, 1, shifted[p]);
                                for (int p = 0; p < 3; ++p) {
                                        for (int q = 0; q < 3; ++q) {
                                                uint64_t pairs = 0;
                                                for (size_t w = 0; w < nw; ++w) {
                                                        pairs += static_cast<uint64_t>(std::popcount(masks[p][w] & shifted[q][w]));
                                                }
                                                mix(pairs);
                                        }
                                }
                        }
                        return static_cast<size_t>(h);
                }

                /// Invertible, coefficient-independent key of the cycle, consistent
                /// with operator==: key = string_key * 2^32 + n_qubits, where
                /// string_key is PauliString::get_hash() of the unified representative
                /// (so all rotations share it). Little-endian 64-bit limbs.
                std::vector<uint64_t> get_hash() const {
                        ensure_unified();
                        const std::vector<uint64_t> s =
                                PauliString<Coeff>(ux(), uy(), pauliengine_detail::coeff_one<Coeff>()).get_hash();
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
                static PauliCycle from_hash(const std::vector<uint64_t>& limbs,
                                            Coeff coeff = pauliengine_detail::coeff_one<Coeff>()) {
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
                PauliCycleSum<Coeff> operator*(const PauliCycle& other) const {
                        require_same_n_qubits(other, "operator*");
                        const int n = this->n_qubits;
                        std::vector<PauliCycle> result;
                        result.reserve(n);
                        const size_t nw = n_words();
                        const size_t na = std::min(this->Base.x.size(), this->Base.y.size());
                        const Coeff coeff = this->Base.coeff * other.Base.coeff;
                        WordVec rx, ry;
                        std::vector<uint64_t> px(nw), py(nw);
                        for (int b = 0; b < n; ++b) {
                                rotate_words_into(other.Base.x, n, b, rx);
                                rotate_words_into(other.Base.y, n, b, ry);
                                const int power_of_i = multiply_symplectic_words(
                                        this->Base.x.data(), this->Base.y.data(), na, rx.data(), ry.data(), nw,
                                        px.data(), py.data());
                                // as in PauliString::operator*
                                Coeff c = coeff;
                                if (power_of_i & 1) {
                                        c *= std::complex<double>(0.0, 1.0);
                                }
                                if (power_of_i & 2) {
                                        c = -c;
                                }
                                result.push_back(PauliCycle(n, PauliString<Coeff>(WordVec(px), WordVec(py), c), trusted_t{}));
                        }
                        return PauliCycleSum<Coeff>(std::move(result));
                }

                /// The cycle with its coefficient multiplied by `scale`.
                PauliCycle operator*(double scale) const {
                        PauliCycle out = *this;
                        out.Base.coeff *= scale;
                        return out;
                }

                /// The cycle with its coefficient multiplied by `scale`.
                PauliCycle operator*(int scale) const { return *this * static_cast<double>(scale); }

                /// [C(P), C(Q)] = sum_k C([P, R^k Q]) (Eq. 40). Only shifts k at which
                /// the supports overlap can contribute (Corollary 3); different shifts
                /// landing in the same rotation orbit are merged by PauliCycleSum.
                PauliCycleSum<Coeff> commutator(const PauliCycle& other) const {
                        require_same_n_qubits(other, "commutator");
                        std::vector<PauliCycle> result;
                        const int n = this->n_qubits;

                        const std::vector<int> shifts = overlap_shifts(other);
                        result.reserve(shifts.size());
                        const size_t nw = n_words();
                        WordVec rx, ry;
                        std::vector<uint64_t> px(nw), py(nw);
                        const size_t na = std::min(this->Base.x.size(), this->Base.y.size());
                        for (int k : shifts) {
                                rotate_words_into(other.Base.x, n, k, rx);
                                rotate_words_into(other.Base.y, n, k, ry);
                                uint64_t anti = 0;
                                for (size_t w = 0; w < na; ++w) {
                                        anti ^= (this->Base.x[w] & ry[w]) ^ (rx[w] & this->Base.y[w]);
                                }
                                if ((std::popcount(anti) & 1) == 0) {
                                        continue;  // commuting: no commutator term
                                }
                                const int power_of_i = multiply_symplectic_words(
                                        this->Base.x.data(), this->Base.y.data(), na, rx.data(), ry.data(), nw,
                                        px.data(), py.data());
                                PauliString<Coeff> term(WordVec(px), WordVec(py),
                                        PauliString<Coeff>::commutator_coeff(this->Base.coeff, other.Base.coeff, power_of_i));
                                if (!pauliengine_detail::coeff_is_zero(term.coeff)) {
                                        result.push_back(PauliCycle(n, std::move(term), trusted_t{}));
                                }
                        }
                        return PauliCycleSum<Coeff>(std::move(result));
                }

                /// Coefficient of a single target cycle in [this, other], i.e. the sum
                /// of the coefficients of the cycles in the full commutator whose
                /// rotation orbit equals `target`. Uses the same support-difference
                /// shift pruning as commutator() but, instead of constructing and
                /// merging every output cycle, only accumulates the shifts whose
                /// product belongs to the target's cyclic class (PauliEngine Features).
                /// Per shift, the operator of P R^k Q is the XOR of the words; its weight
                /// is checked first (rotations preserve it), then one phase pass gives
                /// the commutation parity and i^k, and only the remaining candidates
                /// are unified and compared with the target. No product string, rotated
                /// cycle or coefficient is formed for a shift that cannot contribute.
                template <typename C>
                Coeff targeted_commutator(const PauliCycle& other, const PauliCycle<C>& target) const {
                        require_same_n_qubits(other, "targeted_commutator");
                        require_same_n_qubits(target, "targeted_commutator");
                        const int n = this->n_qubits;
                        const size_t nw = n_words();
                        Coeff total = pauliengine_detail::coeff_zero<Coeff>();
                        const size_t target_weight = target.Base.size();
                        target.ensure_unified();

                        WordVec cx(nw), cy(nw), ux, uy;
                        for (int k : overlap_shifts(other)) {
                                const WordVec rx = rotate_words(other.Base.x, n, k);
                                const WordVec ry = rotate_words(other.Base.y, n, k);
                                // Operator of the product and its weight in one pass.
                                size_t weight = 0;
                                for (size_t w = 0; w < nw; ++w) {
                                        const uint64_t bx = w < Base.x.size() ? Base.x[w] : 0ULL;
                                        const uint64_t by = w < Base.y.size() ? Base.y[w] : 0ULL;
                                        cx[w] = bx ^ rx[w];
                                        cy[w] = by ^ ry[w];
                                        weight += static_cast<size_t>(std::popcount(cx[w] | cy[w]));
                                }
                                if (weight != target_weight) {
                                        continue;
                                }
                                const int power_of_i = symplectic_product_phase(
                                        Base.x.data(), Base.y.data(), Base.x.size(),
                                        rx.data(), ry.data(), nw);
                                if ((power_of_i & 1) == 0) {
                                        continue;  // commuting: no commutator term
                                }
                                unify_words(cx, cy, n, ux, uy);
                                if (cmp_words(ux, target.ux(), nw) != 0
                                    || cmp_words(uy, target.uy(), nw) != 0) {
                                        continue;
                                }
                                total = total + PauliString<Coeff>::commutator_coeff(
                                        Base.coeff, other.Base.coeff, power_of_i);
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

                template <typename C>
                void require_same_n_qubits(const PauliCycle<C>& other, const char* op) const {
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
                                                // Every shift found (typical for dense bases):
                                                // the remaining pairs can add nothing.
                                                if (shifts.size() == static_cast<size_t>(n)) {
                                                        return shifts;
                                                }
                                        }
                                }
                        }
                        return shifts;
                }

                /// Lazily-computed unified representative, cached so repeated ==/hash
                /// on the same cycle don't recompute it. Depends only on Base's operator
                /// part (x, y), which is set at construction and never mutated
                /// afterwards (operator* only touches the coefficient); unify() keeps it
                /// consistent.
                mutable bool unified_valid_ = false;
                mutable WordVec unified_x_;
                mutable WordVec unified_y_;
                /// Base itself is the unified representative (after unify(), and for the
                /// summands of sparse sum products), so no copy of it is cached.
                bool base_unified_ = false;
                /// The unified words, from Base or the cache; call ensure_unified() first.
                const WordVec& ux() const { return base_unified_ ? Base.x : unified_x_; }
                const WordVec& uy() const { return base_unified_ ? Base.y : unified_y_; }

                /// A rotation orbit is determined by the cyclic sequence of
                /// (gap to the next non-identity site, Pauli) pairs along the support:
                /// two strings are rotations of each other iff these sequences are
                /// cyclic shifts of each other. The unified representative places
                /// the lexicographically least rotation of that sequence at qubit 0.
                /// Cost O(n / 64 + w) for weight w, instead of comparing all n
                /// rotations (O(n^2 / 64)).
                void ensure_unified() const {
                        if (unified_valid_ || base_unified_) {
                                return;
                        }
                        unify_words(Base.x, Base.y, n_qubits, unified_x_, unified_y_);
                        unified_valid_ = true;
                }

                /// Unified representative (ux, uy) of the rotation orbit of the operator
                /// (x, y) on n qubits, with exactly ceil(n / 64) words each; see
                /// ensure_unified(). Two operators lie in the same orbit iff their
                /// unified words are equal.
                static void unify_words(const WordVec& x, const WordVec& y, int n, WordVec& ux, WordVec& uy) {
                        const size_t nw = (static_cast<size_t>(n) + BITS_IN_INTEGER - 1) / BITS_IN_INTEGER;
                        ux = WordVec(nw);
                        uy = WordVec(nw);

                        // item = gap << 2 | pauli, with pauli = x_bit | y_bit << 1, in one pass
                        // over the words; the last gap wraps around to the first site.
                        thread_local std::vector<uint64_t> items;
                        items.clear();
                        const size_t words = std::min(x.size(), y.size());
                        long long first = -1;
                        long long prev = -1;
                        uint64_t prev_pauli = 0;
                        for (size_t w = 0; w < words; ++w) {
                                uint64_t bits = x[w] | y[w];
                                while (bits) {
                                        const int b = std::countr_zero(bits);
                                        const long long pos = static_cast<long long>(w * BITS_IN_INTEGER) + b;
                                        const uint64_t pauli = ((x[w] >> b) & 1ULL) | (((y[w] >> b) & 1ULL) << 1);
                                        if (prev >= 0) {
                                                items.push_back((static_cast<uint64_t>(pos - prev) << 2) | prev_pauli);
                                        } else {
                                                first = pos;
                                        }
                                        prev = pos;
                                        prev_pauli = pauli;
                                        bits &= bits - 1;
                                }
                        }
                        if (prev < 0) {
                                return;
                        }
                        items.push_back((static_cast<uint64_t>(first + n - prev) << 2) | prev_pauli);
                        const size_t w = items.size();
                        const size_t start = least_rotation(items);
                        int q = 0;
                        for (size_t t = 0; t < w; ++t) {
                                const size_t at = start + t < w ? start + t : start + t - w;
                                const uint64_t item = items[at];
                                if (item & 1ULL) {
                                        set_bit(ux, q);
                                }
                                if (item & 2ULL) {
                                        set_bit(uy, q);
                                }
                                q += static_cast<int>(item >> 2);
                        }
                }

                /// Hash of unified words, so that operators of the same orbit collide.
                static size_t hash_unified(const WordVec& ux, const WordVec& uy) {
                        size_t h = 0xcbf29ce484222325ULL;
                        mix_words(h, ux);
                        mix_words(h, uy);
                        return h;
                }

                /// Start index of the lexicographically least rotation of s, O(|s|)
                /// (two-pointer minimum-expression algorithm).
                static size_t least_rotation(const std::vector<uint64_t>& s) {
                        const size_t m = s.size();
                        size_t i = 0, j = 1, k = 0;
                        while (i < m && j < m && k < m) {
                                // i, j, k < m, so one subtraction replaces the modulo
                                const size_t ia = i + k < m ? i + k : i + k - m;
                                const size_t jb = j + k < m ? j + k : j + k - m;
                                const uint64_t a = s[ia];
                                const uint64_t b = s[jb];
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
                        WordVec out;
                        rotate_words_into(in, n, shift, out);
                        return out;
                }

                /// rotate_words() into `out`, which is resized to ceil(n / 64) words and
                /// keeps its memory between calls.
                static void rotate_words_into(const WordVec& in, int n, int shift, WordVec& out) {
                        const size_t nw = (static_cast<size_t>(n) + BITS_IN_INTEGER - 1) / BITS_IN_INTEGER;
                        out.resize(nw);
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


                /// Mixes the words of v up to its last non-zero one, so a trimmed base and
                /// full-length unified words of the same operator hash alike.
                static void mix_words(size_t& h, const WordVec& v) {
                        size_t used = v.size();
                        while (used > 0 && v[used - 1] == 0) {
                                --used;
                        }
                        for (size_t i = 0; i < used; ++i) {
                                h ^= std::hash<uint64_t>{}(v[i]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                        }
                }

                static std::vector<int> support(const PauliString<Coeff>& ps) {
                        return support_words(ps.x, ps.y);
                }

                static std::vector<int> support_words(const WordVec& x, const WordVec& y) {
                        std::vector<int> out;
                        const size_t nw = std::min(x.size(), y.size());
                        for (size_t w = 0; w < nw; ++w) {
                                uint64_t bits = x[w] | y[w];
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


template <typename Coeff>
PauliCycle<Coeff> operator*(double scale, const PauliCycle<Coeff>& c) { return c * scale; }
template <typename Coeff>
PauliCycle<Coeff> operator*(int scale, const PauliCycle<Coeff>& c) { return c * scale; }
