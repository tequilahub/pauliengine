#pragma once

/// @file
/// Arithmetic between all operator representations.
///
/// PauliString, QubitHamiltonian, PauliCycle, PauliCycleSum, PauliOrbit and
/// PauliOrbitSum all describe linear combinations of Pauli strings. This header
/// overloads `+`, `-`, `*` and provides a free `commutator(a, b)` for every pair:
///
/// - Within the cycle family (PauliCycle, PauliCycleSum) the result stays a
///   PauliCycleSum, within the orbit family (PauliOrbit, PauliOrbitSum) a
///   PauliOrbitSum, so nothing is expanded unnecessarily.
/// - PauliString * PauliString and its commutator stay a PauliString;
///   QubitHamiltonian with QubitHamiltonian of the same coefficient type uses the
///   class's own operators.
/// - Every other combination expands both operands with as_hamiltonian() and
///   returns a QubitHamiltonian. Its coefficients are symbolic as soon as one
///   operand is symbolic, otherwise complex.

#include <algorithm>
#include <bit>
#include <complex>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/PauliCycle.h"
#include "pauliengine/PauliOrbit.h"

// ---------------------------------------------------------------------------
// Conversion of every representation to a QubitHamiltonian
// ---------------------------------------------------------------------------

/// A single Pauli string as a one-term Hamiltonian.
template <typename C>
QubitHamiltonian<C> as_hamiltonian(const PauliString<C>& ps) {
        return QubitHamiltonian<C>(std::vector<PauliString<C>>{ps});
}

/// A Hamiltonian is already one.
template <typename C>
const QubitHamiltonian<C>& as_hamiltonian(const QubitHamiltonian<C>& h) {
        return h;
}

/// A cycle expanded into the sum of all its rotations.
inline QubitHamiltonian<std::complex<double>> as_hamiltonian(const PauliCycle& c) {
        return c.to_qubit_hamiltonian();
}

/// A cycle sum expanded into the sum of all rotations of all summands.
inline QubitHamiltonian<std::complex<double>> as_hamiltonian(const PauliCycleSum& s) {
        return s.to_qubit_hamiltonian();
}

/// An orbit expanded into all of its Pauli strings (see PauliOrbit::to_qubit_hamiltonian).
inline QubitHamiltonian<std::complex<double>> as_hamiltonian(const PauliOrbit& o) {
        return o.to_qubit_hamiltonian();
}

/// An orbit sum expanded into all Pauli strings of all its orbits.
inline QubitHamiltonian<std::complex<double>> as_hamiltonian(const PauliOrbitSum& s) {
        return s.to_qubit_hamiltonian();
}

namespace pauliengine_detail {

using Complex = std::complex<double>;
using Symbolic = SymEngine::Expression;

template <typename C>
C coeff_of(const QubitHamiltonian<C>&);

/// Coefficient type of as_hamiltonian(T).
template <typename T>
using hamiltonian_coeff_t = decltype(coeff_of(as_hamiltonian(std::declval<const T&>())));

/// Common coefficient type: symbolic if either side is symbolic.
template <typename A, typename B>
using promoted_coeff_t = std::conditional_t<
        std::is_same_v<hamiltonian_coeff_t<A>, Symbolic> || std::is_same_v<hamiltonian_coeff_t<B>, Symbolic>,
        Symbolic, Complex>;

/// A Hamiltonian with its coefficients converted to `To` (only complex -> symbolic is needed).
template <typename To, typename From>
QubitHamiltonian<To> convert_coeffs(const QubitHamiltonian<From>& h) {
        if constexpr (std::is_same_v<To, From>) {
                return h;
        } else {
                std::vector<PauliString<To>> terms;
                terms.reserve(h.data.size());
                for (const auto& ps : h.data) {
                        terms.emplace_back(ps.x, ps.y, To(ps.coeff));
                }
                return QubitHamiltonian<To>(std::move(terms), already_compacted);
        }
}

/// `a` as a Hamiltonian with the coefficient type common to A and B.
template <typename A, typename B>
QubitHamiltonian<promoted_coeff_t<A, B>> lift(const A& a) {
        return convert_coeffs<promoted_coeff_t<A, B>>(as_hamiltonian(a));
}

template <typename T> inline constexpr bool is_pauli_string = false;
template <typename C> inline constexpr bool is_pauli_string<PauliString<C>> = true;
template <typename T> inline constexpr bool is_hamiltonian = false;
template <typename C> inline constexpr bool is_hamiltonian<QubitHamiltonian<C>> = true;
template <typename T> inline constexpr bool is_cycle_family =
        std::is_same_v<T, PauliCycle> || std::is_same_v<T, PauliCycleSum>;
template <typename T> inline constexpr bool is_orbit_family =
        std::is_same_v<T, PauliOrbit> || std::is_same_v<T, PauliOrbitSum>;

template <typename A, typename B> inline constexpr bool cycle_pair = is_cycle_family<A> && is_cycle_family<B>;
template <typename A, typename B> inline constexpr bool orbit_pair = is_orbit_family<A> && is_orbit_family<B>;
template <typename A, typename B> inline constexpr bool same_hamiltonian =
        is_hamiltonian<A> && std::is_same_v<A, B>;

/// Pairs whose `+` / `-` is defined natively (class members or the family operators below).
template <typename A, typename B> inline constexpr bool native_add =
        cycle_pair<A, B> || orbit_pair<A, B> || same_hamiltonian<A, B>;
/// Pairs whose `*` is defined natively.
template <typename A, typename B> inline constexpr bool native_mul =
        native_add<A, B> || (is_pauli_string<A> && is_pauli_string<B>);

inline PauliCycleSum as_cycle_sum(const PauliCycle& c) { return PauliCycleSum(std::vector<PauliCycle>{c}); }
inline const PauliCycleSum& as_cycle_sum(const PauliCycleSum& s) { return s; }

inline PauliOrbitSum as_orbit_sum(const PauliOrbit& o) {
        PauliOrbitSum s(o.n_qubits);
        s.add(o, Complex(1.0, 0.0));
        return s;
}
inline const PauliOrbitSum& as_orbit_sum(const PauliOrbitSum& s) { return s; }

/// a + sign * b within the cycle family; summands of the same orbit are merged.
template <typename A, typename B>
PauliCycleSum cycle_combine(const A& a, const B& b, double sign) {
        const PauliCycleSum& sa = as_cycle_sum(a);   // binds a temporary for a single cycle
        const PauliCycleSum& sb = as_cycle_sum(b);
        std::vector<PauliCycle> terms = sa.data;
        for (PauliCycle c : sb.data) {
                c.Base.coeff *= sign;
                terms.push_back(std::move(c));
        }
        return PauliCycleSum(std::move(terms));
}

/// a + sign * b within the orbit family.
template <typename A, typename B>
PauliOrbitSum orbit_combine(const A& a, const B& b, double sign) {
        PauliOrbitSum out = as_orbit_sum(a);
        const auto& sb = as_orbit_sum(b);
        if (sb.n_qubits != out.n_qubits) {
                throw std::invalid_argument("PauliOrbitSum: n_qubits mismatch (" + std::to_string(out.n_qubits)
                        + " vs " + std::to_string(sb.n_qubits) + ")");
        }
        for (size_t k = 0; k < sb.data.size(); ++k) {
                out.add(sb.data[k], sign * sb.coeffs[k]);
        }
        out.prune_zeros();
        return out;
}

}  // namespace pauliengine_detail

/// Every type that as_hamiltonian() can expand.
template <typename T>
concept OperatorLike = requires(const T& t) { as_hamiltonian(t); };

// ---------------------------------------------------------------------------
// Same-family operators: stay in the compressed representation
// ---------------------------------------------------------------------------

/// Sum of two cycle-family operators, as a PauliCycleSum.
template <typename A, typename B>
        requires pauliengine_detail::cycle_pair<A, B>
PauliCycleSum operator+(const A& a, const B& b) {
        return pauliengine_detail::cycle_combine(a, b, 1.0);
}

/// Difference of two cycle-family operators, as a PauliCycleSum.
template <typename A, typename B>
        requires pauliengine_detail::cycle_pair<A, B>
PauliCycleSum operator-(const A& a, const B& b) {
        return pauliengine_detail::cycle_combine(a, b, -1.0);
}

/// Product of two cycle-family operators, as a PauliCycleSum.
template <typename A, typename B>
        requires pauliengine_detail::cycle_pair<A, B>
PauliCycleSum operator*(const A& a, const B& b) {
        using namespace pauliengine_detail;
        return as_cycle_sum(a) * as_cycle_sum(b);
}

/// Sum of two orbit-family operators, as a PauliOrbitSum.
template <typename A, typename B>
        requires pauliengine_detail::orbit_pair<A, B>
PauliOrbitSum operator+(const A& a, const B& b) {
        return pauliengine_detail::orbit_combine(a, b, 1.0);
}

/// Difference of two orbit-family operators, as a PauliOrbitSum.
template <typename A, typename B>
        requires pauliengine_detail::orbit_pair<A, B>
PauliOrbitSum operator-(const A& a, const B& b) {
        return pauliengine_detail::orbit_combine(a, b, -1.0);
}

/// Product of two orbit-family operators, as a PauliOrbitSum.
template <typename A, typename B>
        requires pauliengine_detail::orbit_pair<A, B>
PauliOrbitSum operator*(const A& a, const B& b) {
        using namespace pauliengine_detail;
        return as_orbit_sum(a) * as_orbit_sum(b);
}

// ---------------------------------------------------------------------------
// Mixed operators: expand both sides and return a QubitHamiltonian
// ---------------------------------------------------------------------------

/// a + b for operators of different kinds, as a QubitHamiltonian.
template <typename A, typename B>
        requires(OperatorLike<A> && OperatorLike<B> && !pauliengine_detail::native_add<A, B>)
auto operator+(const A& a, const B& b) {
        using namespace pauliengine_detail;
        return lift<A, B>(a) + lift<B, A>(b);
}

/// a - b for operators of different kinds, as a QubitHamiltonian.
template <typename A, typename B>
        requires(OperatorLike<A> && OperatorLike<B> && !pauliengine_detail::native_add<A, B>)
auto operator-(const A& a, const B& b) {
        using namespace pauliengine_detail;
        return lift<A, B>(a) - lift<B, A>(b);
}

/// a * b for operators of different kinds, as a QubitHamiltonian.
template <typename A, typename B>
        requires(OperatorLike<A> && OperatorLike<B> && !pauliengine_detail::native_mul<A, B>)
auto operator*(const A& a, const B& b) {
        using namespace pauliengine_detail;
        return lift<A, B>(a) * lift<B, A>(b);
}

/// Commutator [a, b] = a b - b a for any two operators.
///
/// Same type with a commutator member: that member (PauliString stays a
/// PauliString, cycles and orbits stay sums). Same family: the family's sum
/// commutator. Otherwise both sides are expanded and the result is a
/// QubitHamiltonian.
template <typename A, typename B>
        requires(OperatorLike<A> && OperatorLike<B>)
auto commutator(const A& a, const B& b) {
        using namespace pauliengine_detail;
        if constexpr (std::is_same_v<A, B> && requires { a.commutator(b); }) {
                return a.commutator(b);
        } else if constexpr (cycle_pair<A, B>) {
                return as_cycle_sum(a).commutator(as_cycle_sum(b));
        } else if constexpr (orbit_pair<A, B>) {
                return as_orbit_sum(a).commutator(as_orbit_sum(b));
        } else {
                return lift<A, B>(a).commutator(lift<B, A>(b));
        }
}

// ---------------------------------------------------------------------------
// Hilbert-Schmidt inner product
// ---------------------------------------------------------------------------
//
// inner_product(a, b) = Tr(a^dagger b) / 2^n, the Hilbert-Schmidt inner product
// normalized so that Pauli strings are orthonormal (arXiv:2604.16701, Cor. 1 and
// Lemma 1). It is conjugate-linear in a and does not depend on the number of
// qubits the operators are embedded in; the unnormalized value Tr(a^dagger b) on
// n qubits is 2^n times it.

namespace pauliengine_detail {

inline Complex conj_coeff(const Complex& c) { return std::conj(c); }
inline Symbolic conj_coeff(const Symbolic& c) { return Symbolic(SymEngine::conjugate(c.get_basic())); }

/// Hash of a Pauli string's operator part, independent of the coefficient type.
struct OperatorKeyHash {
        size_t operator()(const std::vector<uint64_t>& v) const {
                size_t h = v.size();
                for (uint64_t w : v) {
                        h ^= std::hash<uint64_t>{}(w) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                }
                return h;
        }
};

/// sum over common operators of conj(a_P) b_P (hash join on the smaller side).
template <typename R, typename CA, typename CB>
R hamiltonian_inner(const QubitHamiltonian<CA>& a, const QubitHamiltonian<CB>& b) {
        R acc(0);
        if (a.data.size() <= b.data.size()) {
                std::unordered_map<std::vector<uint64_t>, const PauliString<CA>*, OperatorKeyHash> index;
                index.reserve(a.data.size());
                for (const auto& p : a.data) index.emplace(p.get_hash(), &p);
                for (const auto& q : b.data) {
                        const auto it = index.find(q.get_hash());
                        if (it != index.end()) acc = acc + R(conj_coeff(it->second->coeff)) * R(q.coeff);
                }
        } else {
                std::unordered_map<std::vector<uint64_t>, const PauliString<CB>*, OperatorKeyHash> index;
                index.reserve(b.data.size());
                for (const auto& q : b.data) index.emplace(q.get_hash(), &q);
                for (const auto& p : a.data) {
                        const auto it = index.find(p.get_hash());
                        if (it != index.end()) acc = acc + R(conj_coeff(p.coeff)) * R(it->second->coeff);
                }
        }
        return acc;
}

/// (number of X, Y, Z sites) of a Pauli string and the highest qubit it acts on (-1 for I).
struct StringLabel {
        int p, q, r, top;
};
template <typename C>
StringLabel string_label(const PauliString<C>& s) {
        StringLabel l{0, 0, 0, -1};
        for (size_t w = 0; w < s.x.size(); ++w) {
                const uint64_t x = s.x[w], y = s.y[w];
                l.p += std::popcount(x & ~y);
                l.q += std::popcount(~x & y);
                l.r += std::popcount(x & y);
                if (x | y) l.top = static_cast<int>(w * BITS_IN_INTEGER) + 63 - std::countl_zero(x | y);
        }
        return l;
}

/// <orbit sum, H>: the orbit B_{pqr} on n qubits is (i/N) times the sum of all
/// strings with label (p,q,r) on qubits 0..n-1, so only H's strings with that
/// label contribute. The orbit is never expanded.
template <typename R, typename C>
R orbit_sum_inner(const PauliOrbitSum& a, const QubitHamiltonian<C>& h) {
        // label -> conj(a_k) * conj(i / N_k) = conj(a_k) * (-i / N_k)
        std::unordered_map<int, Complex> weight;
        for (size_t k = 0; k < a.data.size(); ++k) {
                const PauliOrbit& o = a.data[k];
                const int key = (o.p * (PauliOrbit::kMaxQubits + 1) + o.q) * (PauliOrbit::kMaxQubits + 1) + o.r;
                weight[key] += std::conj(a.coeffs[k]) * Complex(0.0, -1.0 / o.n_terms());
        }
        R acc(0);
        for (const auto& s : h.data) {
                const StringLabel l = string_label(s);
                if (l.top >= a.n_qubits) continue;
                const int key = (l.p * (PauliOrbit::kMaxQubits + 1) + l.q) * (PauliOrbit::kMaxQubits + 1) + l.r;
                const auto it = weight.find(key);
                if (it != weight.end()) acc = acc + R(it->second) * R(s.coeff);
        }
        return acc;
}

/// <A, B> for two orbit sums (Lemma 1): orbits with equal labels on n_a and n_b
/// qubits share the N(label; min(n_a, n_b)) strings on the smaller register.
inline Complex orbit_family_inner(const PauliOrbitSum& a, const PauliOrbitSum& b) {
        const int n_min = std::min(a.n_qubits, b.n_qubits);
        Complex acc(0.0, 0.0);
        for (size_t k = 0; k < a.data.size(); ++k) {
                const PauliOrbit& oa = a.data[k];
                if (oa.p + oa.q + oa.r > n_min) continue;      // no common strings
                const PauliOrbit ob(b.n_qubits, oa.p, oa.q, oa.r);
                const Complex cb = b.contains(ob);
                if (cb == Complex(0.0, 0.0)) continue;
                const double shared = PauliOrbit(n_min, oa.p, oa.q, oa.r).n_terms();
                // conj(i / N_a) * (i / N_b) = 1 / (N_a N_b), times the shared strings
                acc += std::conj(a.coeffs[k]) * cb * (shared / (oa.n_terms() * ob.n_terms()));
        }
        return acc;
}

/// <A, B> for two cycle sums on the same n: equal orbits give conj(c_P) c_Q n^2 / d.
inline Complex cycle_family_inner(const PauliCycleSum& a, const PauliCycleSum& b) {
        std::unordered_map<std::vector<uint64_t>, Complex, OperatorKeyHash> index;
        index.reserve(b.data.size());
        for (const auto& c : b.data) index.emplace(c.get_hash(), c.Base.coeff);
        Complex acc(0.0, 0.0);
        for (const auto& c : a.data) {
                const auto it = index.find(c.get_hash());
                if (it == index.end()) continue;
                const double n = c.n_qubits;
                acc += std::conj(c.Base.coeff) * it->second * (n * n / c.period());
        }
        return acc;
}

template <typename A, typename B>
using inner_t = promoted_coeff_t<A, B>;

}  // namespace pauliengine_detail

/// Normalized Hilbert-Schmidt inner product Tr(a^dagger b) / 2^n of any two
/// operators (arXiv:2604.16701, Corollary 1 and Lemma 1). Conjugate-linear in
/// `a`; symbolic if either operand is symbolic, otherwise complex. Orbits are
/// never expanded, so this also works for orbits too large to write out.
template <typename A, typename B>
        requires(OperatorLike<A> && OperatorLike<B>)
pauliengine_detail::inner_t<A, B> inner_product(const A& a, const B& b) {
        using namespace pauliengine_detail;
        using R = inner_t<A, B>;
        if constexpr (orbit_pair<A, B>) {
                return orbit_family_inner(as_orbit_sum(a), as_orbit_sum(b));
        } else if constexpr (is_orbit_family<A>) {
                return orbit_sum_inner<R>(as_orbit_sum(a), lift<B, A>(b));
        } else if constexpr (is_orbit_family<B>) {
                return conj_coeff(orbit_sum_inner<R>(as_orbit_sum(b), lift<A, B>(a)));
        } else if constexpr (cycle_pair<A, B>) {
                const PauliCycleSum& sa = as_cycle_sum(a);
                const PauliCycleSum& sb = as_cycle_sum(b);
                if (sa.data.empty() || sb.data.empty()) return R(0);
                if (sa.data.front().n_qubits == sb.data.front().n_qubits) return cycle_family_inner(sa, sb);
                return hamiltonian_inner<R>(as_hamiltonian(sa), as_hamiltonian(sb));
        } else {
                return hamiltonian_inner<R>(lift<A, B>(a), lift<B, A>(b));
        }
}
