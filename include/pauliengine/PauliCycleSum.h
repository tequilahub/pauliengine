#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <iterator>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "pauliengine/detail/Coefficients.h"
#include "pauliengine/detail/SparseCycles.h"
#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"


// Every use of PauliCycle below depends on Coeff, so the member functions are
// compiled only where they are used, after PauliCycle.h (included at the end).
template <typename Coeff>
class PauliCycle;


/// A sum of PauliCycles, sum_i data[i]. Each cycle carries its own coefficient
/// (in its base string); there is no separate global prefactor. The sum is kept
/// merged: no two summands describe the same rotation orbit, and summands whose
/// coefficient cancelled to zero are dropped. Coeff is std::complex<double> or
/// symbolic::Expression.
template <typename Coeff>
class PauliCycleSum {
        public:
                /// The summands; no two lie in the same rotation orbit.
                std::vector<PauliCycle<Coeff>> data;

                /// The empty sum (zero operator).
                PauliCycleSum() = default;

                /// Merges summands from the same rotation orbit and drops zeros.
                explicit PauliCycleSum(std::vector<PauliCycle<Coeff>> cycles);

                /// The sum expanded into a QubitHamiltonian (all rotations of all summands).
                QubitHamiltonian<Coeff> to_qubit_hamiltonian() const;

                /// Coefficient of the summand in the same rotation orbit as `cycle`
                /// (its own coefficient is ignored), or 0 if there is none.
                template <typename C>
                Coeff contains(const PauliCycle<C>& cycle) const;

                /// Commutator [this, other] = sum_{i,j} [c_i, d_j], merged.
                PauliCycleSum commutator(const PauliCycleSum& other) const;

                /// Product this * other = sum_{i,j} c_i d_j, merged.
                PauliCycleSum operator*(const PauliCycleSum& other) const;

                /// Every coefficient multiplied by `scale`; scaling by 0 gives the empty sum.
                PauliCycleSum operator*(double scale) const;
                PauliCycleSum operator*(int scale) const { return *this * static_cast<double>(scale); }
                /// Every coefficient multiplied by a complex `scale`; scaling by 0 gives the empty sum.
                PauliCycleSum operator*(std::complex<double> scale) const {
                        if (scale == std::complex<double>(0.0)) return PauliCycleSum();
                        PauliCycleSum out = *this;
                        for (auto& c : out.data) c.Base.coeff = c.Base.coeff * Coeff(scale);
                        return out;
                }

                /// Derivative of every coefficient by the symbol `name`; summands whose
                /// derivative is zero are dropped (symbolic coefficients).
                PauliCycleSum diff(const std::string& name) const {
                        std::vector<PauliCycle<Coeff>> terms;
                        terms.reserve(data.size());
                        for (const auto& c : data) terms.push_back(c.diff(name));
                        return PauliCycleSum(std::move(terms));
                }

                /// The sum with values or polynomials substituted for symbols; summands
                /// that become zero are dropped.
                template <typename Values>
                PauliCycleSum substitute(const Values& values) const {
                        const auto prepared = Coeff::prepare_subs(values);
                        std::vector<PauliCycle<Coeff>> terms;
                        terms.reserve(data.size());
                        for (const auto& c : data) terms.push_back(c.substitute(prepared));
                        return PauliCycleSum(std::move(terms));
                }

                /// Coefficient of `target`'s orbit in [this, other], without building
                /// the full commutator.
                template <typename C>
                Coeff targeted_commutator(const PauliCycleSum& other, const PauliCycle<C>& target) const;

                /// Rewrite every summand's base to the unified representative of its
                /// orbit (in place). Cheap: the unified forms are already cached by the
                /// merging in the constructor.
                void unify();

                /// Number of summands.
                size_t size() const { return data.size(); }

        private:
                /// sum_ij c_i d_j, or sum_ij [c_i, d_j], formed on site lists (see
                /// detail/SparseCycles.h); nothing if a base acts on too many sites for
                /// that, or the summands differ in n_qubits.
                std::optional<PauliCycleSum> sparse_pairs(const PauliCycleSum& other, bool commutator) const;
};


template <typename Coeff>
QubitHamiltonian<Coeff> PauliCycleSum<Coeff>::to_qubit_hamiltonian() const {
        std::vector<PauliString<Coeff>> terms;
        for (const auto& cycle : data) {
                std::vector<PauliString<Coeff>> rots = cycle.rotations();
                terms.insert(terms.end(),
                        std::make_move_iterator(rots.begin()),
                        std::make_move_iterator(rots.end()));
        }
        return QubitHamiltonian<Coeff>(std::move(terms));
}

template <typename Coeff>
PauliCycleSum<Coeff>::PauliCycleSum(std::vector<PauliCycle<Coeff>> cycles) {
        // Hash-join on the cycle (the base up to rotation). For cycles with many sites a cheap
        // rotation-invariant signature comes first: a cycle whose signature no other
        // summand has starts a new summand without being unified. Summands sharing a
        // signature, and sparse cycles, whose unified form is cheap anyway, are joined
        // on the unified form, so the merging decisions are those of a join on the
        // unified form alone. Whether a cycle takes the signature is the same for all
        // rotations of its base, so equal cycles always meet. The map keys point into `data`; reserving up
        // front guarantees they stay valid (data never grows beyond cycles.size()).
        struct OrbitHash {
                size_t operator()(const PauliCycle<Coeff>* c) const { return c->hash(); }
        };
        struct OrbitEqual {
                bool operator()(const PauliCycle<Coeff>* a, const PauliCycle<Coeff>* b) const { return *a == *b; }
        };
        constexpr size_t kInIndex = static_cast<size_t>(-1);
        // signature -> the one summand with it that is not in `index` yet, or kInIndex
        std::unordered_map<size_t, size_t> by_signature;
        std::unordered_map<const PauliCycle<Coeff>*, size_t, OrbitHash, OrbitEqual> index;
        index.reserve(cycles.size());
        data.reserve(cycles.size());

        for (auto& c : cycles) {
                if (!data.empty() && c.n_qubits != data.front().n_qubits) {
                        throw std::invalid_argument("PauliCycleSum: all cycles must have the same n_qubits ("
                                + std::to_string(data.front().n_qubits) + " vs "
                                + std::to_string(c.n_qubits) + ")");
                }
                if (c.prefers_signature()) {
                        const auto [slot, fresh] = by_signature.try_emplace(c.rotation_signature(), data.size());
                        if (fresh) {
                                data.push_back(std::move(c));
                                continue;
                        }
                        if (slot->second != kInIndex) {
                                index.emplace(&data[slot->second], slot->second);
                                slot->second = kInIndex;
                        }
                }
                const auto it = index.find(&c);
                if (it != index.end()) {
                        Coeff& merged = data[it->second].Base.coeff;
                        merged = merged + c.Base.coeff;
                        continue;
                }
                data.push_back(std::move(c));
                index.emplace(&data.back(), data.size() - 1);
        }

        // Drop summands whose merged coefficient cancelled to zero.
        data.erase(
                std::remove_if(data.begin(), data.end(),
                        [](const PauliCycle<Coeff>& c) { return pauliengine_detail::coeff_is_zero(c.Base.coeff); }),
                data.end());
}

template <typename Coeff>
void PauliCycleSum<Coeff>::unify() {
        for (auto& c : data) {
                c.unify();
        }
}

template <typename Coeff>
template <typename C>
Coeff PauliCycleSum<Coeff>::contains(const PauliCycle<C>& cycle) const {
        // The sum is merged, so at most one summand lies in the cycle's orbit.
        for (const auto& c : data) {
                if (c == cycle) {
                        return c.Base.coeff;
                }
        }
        return pauliengine_detail::coeff_zero<Coeff>();
}

template <typename Coeff>
PauliCycleSum<Coeff> PauliCycleSum<Coeff>::commutator(const PauliCycleSum& other) const {
        // [ sum_i c_i , sum_j d_j ] = sum_{i,j} [c_i, d_j]. Each cycle carries its
        // own coefficient, so the pairwise results are collected and merged.
        if (auto sparse = sparse_pairs(other, true)) {
                return std::move(*sparse);
        }
        std::vector<PauliCycle<Coeff>> combined;
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

template <typename Coeff>
PauliCycleSum<Coeff> PauliCycleSum<Coeff>::operator*(const PauliCycleSum& other) const {
        // ( sum_i c_i )( sum_j d_j ) = sum_{i,j} c_i d_j, merged across all pairs.
        if (auto sparse = sparse_pairs(other, false)) {
                return std::move(*sparse);
        }
        std::vector<PauliCycle<Coeff>> combined;
        for (const auto& ci : data) {
                for (const auto& dj : other.data) {
                        PauliCycleSum part = ci * dj;
                        combined.insert(combined.end(),
                                std::make_move_iterator(part.data.begin()),
                                std::make_move_iterator(part.data.end()));
                }
        }
        return PauliCycleSum(std::move(combined));
}

template <typename Coeff>
std::optional<PauliCycleSum<Coeff>> PauliCycleSum<Coeff>::sparse_pairs(const PauliCycleSum& other, bool commutator) const {
        using Kernel = pauliengine_detail::SparseCycleProducts<Coeff>;
        using pauliengine_detail::SparseString;
        if (data.empty() || other.data.empty()) {
                return std::nullopt;
        }
        const int n = data.front().n_qubits;
        const size_t nw = (static_cast<size_t>(n) + 63) / 64;
        std::vector<SparseString> sa, sb;
        std::vector<Coeff> ca, cb;
        const auto load = [&](const PauliCycleSum& sum, std::vector<SparseString>& sites, std::vector<Coeff>& coeffs) {
                sites.reserve(sum.data.size());
                coeffs.reserve(sum.data.size());
                for (const auto& c : sum.data) {
                        if (c.n_qubits != n) {
                                return false;  // the pairwise route reports the mismatch
                        }
                        sites.push_back(SparseString::of(c.Base.x, c.Base.y));
                        if (!Kernel::suitable(sites.back().pos.size(), nw)) {
                                return false;
                        }
                        coeffs.push_back(c.Base.coeff);
                }
                return true;
        };
        if (!load(*this, sa, ca) || !load(other, sb, cb)) {
                return std::nullopt;
        }
        Kernel kernel;
        kernel.run(sa, ca, sb, cb, static_cast<uint32_t>(n), commutator ? Kernel::Op::Commutator : Kernel::Op::Product);
        // Every summand's base is the unified representative of its class, which the
        // merging has found: no second copy of it, and unify() has nothing to do.
        const auto build = [&](const typename Kernel::Summand& s) {
                WordVec ux(nw), uy(nw);
                Kernel::unified_words(kernel.items.data() + s.off, s.len, ux.data(), uy.data(), nw);
                PauliCycle<Coeff> c(n, PauliString<Coeff>(std::move(ux), std::move(uy), s.coeff),
                                    typename PauliCycle<Coeff>::trusted_t{});
                c.base_unified_ = true;
                return c;
        };
        const size_t m = kernel.summands.size();
        PauliCycleSum out;
        out.data.reserve(m);
#ifdef PAULIENGINE_HAS_OPENMP
        // a summand costs about 100 ns to build, so a team pays off from tens of thousands
        if (std::is_same_v<Coeff, std::complex<double>> && pauliengine_detail::use_threads(m / 8)) {
                // contiguous ranges per thread, appended in order
                const size_t threads = static_cast<size_t>(omp_get_max_threads());
                std::vector<std::vector<PauliCycle<Coeff>>> parts(threads);
                pauliengine_detail::for_each_index(threads, true, [&](size_t t) {
                        const size_t lo = m * t / threads, hi = m * (t + 1) / threads;
                        parts[t].reserve(hi - lo);
                        for (size_t k = lo; k < hi; ++k) parts[t].push_back(build(kernel.summands[k]));
                });
                for (auto& part : parts) {
                        std::move(part.begin(), part.end(), std::back_inserter(out.data));
                }
                return out;
        }
#endif
        for (const auto& s : kernel.summands) {
                out.data.push_back(build(s));
        }
        return out;
}

template <typename Coeff>
PauliCycleSum<Coeff> PauliCycleSum<Coeff>::operator*(double scale) const {
        if (scale == 0.0) {
                return PauliCycleSum();
        }
        PauliCycleSum out = *this;
        for (auto& c : out.data) {
                c.Base.coeff *= scale;
        }
        return out;
}

template <typename Coeff>
PauliCycleSum<Coeff> operator*(double scale, const PauliCycleSum<Coeff>& s) { return s * scale; }
template <typename Coeff>
PauliCycleSum<Coeff> operator*(int scale, const PauliCycleSum<Coeff>& s) { return s * scale; }

template <typename Coeff>
template <typename C>
Coeff PauliCycleSum<Coeff>::targeted_commutator(const PauliCycleSum& other, const PauliCycle<C>& target) const {
        // Cycles carry their own coefficient, so the pairwise targeted coefficients
        // are simply summed (no extra scaling).
        Coeff total = pauliengine_detail::coeff_zero<Coeff>();
        if (data.empty() || other.data.empty()) {
                return total;
        }
        // All summands of a sum share n_qubits, so checking one pair suffices.
        data.front().require_same_n_qubits(other.data.front(), "targeted_commutator");
        data.front().require_same_n_qubits(target, "targeted_commutator");

        // XOR rule on rotation orbits: the pair (P, Q) can contribute only if
        // P R^k Q = R^m T for some shifts k, m, i.e. R^k Q = P XOR R^m T. Hence the
        // orbit of Q is the orbit of P XOR R^m T for one of the d_T distinct
        // rotations of the target (and symmetrically with the roles swapped).
        // Traverse the smaller sum, index the larger one by unified orbit, and
        // evaluate only the pairs found. If the target has more distinct rotations
        // than the larger sum has summands, the lookups would cost more than the
        // plain pair scan, which is used instead.
        const bool iter_this = data.size() <= other.data.size();
        const std::vector<PauliCycle<Coeff>>& outer = iter_this ? data : other.data;
        const std::vector<PauliCycle<Coeff>>& inner = iter_this ? other.data : data;
        const int n = target.n_qubits;
        const int d_target = target.period();
        const auto pair_value = [&](size_t i, size_t j) {
                return iter_this ? outer[i].targeted_commutator(inner[j], target)
                                 : inner[j].targeted_commutator(outer[i], target);
        };

        if (static_cast<size_t>(d_target) > inner.size()) {
                for (size_t i = 0; i < outer.size(); ++i) {
                        for (size_t j = 0; j < inner.size(); ++j) {
                                total = total + pair_value(i, j);
                        }
                }
                return total;
        }

        // Index of `inner` by unified orbit; the sum is merged, so each orbit
        // occurs at most once, but hash collisions are chained.
        constexpr uint32_t kNone = static_cast<uint32_t>(-1);
        std::unordered_map<size_t, uint32_t> head;
        head.reserve(inner.size());
        std::vector<uint32_t> next(inner.size(), kNone);
        for (size_t j = 0; j < inner.size(); ++j) {
                inner[j].ensure_unified();
                const size_t h = PauliCycle<Coeff>::hash_unified(inner[j].ux(), inner[j].uy());
                const auto [it, inserted] = head.try_emplace(h, static_cast<uint32_t>(j));
                if (!inserted) {
                        next[j] = it->second;
                        it->second = static_cast<uint32_t>(j);
                }
        }

        // The d_T distinct rotations of the target.
        const size_t nw = (static_cast<size_t>(n) + BITS_IN_INTEGER - 1) / BITS_IN_INTEGER;
        std::vector<WordVec> rot_x, rot_y;
        rot_x.reserve(d_target);
        rot_y.reserve(d_target);
        for (int m = 0; m < d_target; ++m) {
                rot_x.push_back(PauliCycle<C>::rotate_words(target.Base.x, n, m));
                rot_y.push_back(PauliCycle<C>::rotate_words(target.Base.y, n, m));
        }

        // Partners of the current outer summand; `stamp` deduplicates them.
        std::vector<size_t> stamp(inner.size(), static_cast<size_t>(-1));
        std::vector<uint32_t> partners;
        WordVec cx(nw), cy(nw), ux, uy;
        for (size_t i = 0; i < outer.size(); ++i) {
                const PauliString<Coeff>& p = outer[i].Base;
                partners.clear();
                for (int m = 0; m < d_target; ++m) {
                        for (size_t w = 0; w < nw; ++w) {
                                cx[w] = (w < p.x.size() ? p.x[w] : 0ULL) ^ rot_x[m][w];
                                cy[w] = (w < p.y.size() ? p.y[w] : 0ULL) ^ rot_y[m][w];
                        }
                        PauliCycle<Coeff>::unify_words(cx, cy, n, ux, uy);
                        const auto bucket = head.find(PauliCycle<Coeff>::hash_unified(ux, uy));
                        if (bucket == head.end()) {
                                continue;
                        }
                        for (uint32_t j = bucket->second; j != kNone; j = next[j]) {
                                if (stamp[j] == i) {
                                        continue;
                                }
                                if (PauliCycle<Coeff>::cmp_words(ux, inner[j].ux(), nw) == 0
                                    && PauliCycle<Coeff>::cmp_words(uy, inner[j].uy(), nw) == 0) {
                                        stamp[j] = i;
                                        partners.push_back(j);
                                }
                        }
                }
                for (uint32_t j : partners) {
                        total = total + pair_value(i, j);
                }
        }
        return total;
}


// The cycles themselves, so that this header can be included on its own.
#include "pauliengine/PauliCycle.h"
