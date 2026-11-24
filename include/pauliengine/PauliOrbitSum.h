#pragma once

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pauliengine/detail/Coefficients.h"
#include "pauliengine/QubitHamiltonian.h"

class PauliOrbit;

namespace pauliengine_detail {

/// PauliOrbit, named through a template parameter. PauliOrbit is not a template, so
/// without this every use of it in PauliOrbitSum's member functions would be compiled
/// at once and would need PauliOrbit to be complete; this way they are compiled only
/// where they are used, after PauliOrbit.h (included at the end of this header).
template <typename>
struct orbit_type {
        using type = PauliOrbit;
};

}  // namespace pauliengine_detail

/// A scalar-weighted sum of Pauli orbits, sum_i coeffs[i] * data[i]
/// (arXiv:2604.16701, Section V). See PauliOrbit.h for the orbit data structure.
///
/// Coeff is std::complex<double> or SymEngine::Expression. Commutators of the
/// (anti-Hermitian) orbits have real structure constants, but ordinary products
/// carry a full complex Pauli phase (see PauliOrbit::operator*); a symbolic sum
/// multiplies these numeric constants into its symbolic coefficients.
template <typename Coeff>
class PauliOrbitSum {
        public:
                /// PauliOrbit (see pauliengine_detail::orbit_type).
                using Orbit = typename pauliengine_detail::orbit_type<Coeff>::type;

                /// Fixed at construction (exposed read-only to Python).
                int n_qubits;
                std::vector<Orbit> data;        ///< the orbits, each label at most once
                std::vector<Coeff> coeffs;      ///< coefficients, parallel to data

                /// Empty sum on n qubits. Throws std::invalid_argument unless 0 <= n <= PauliOrbit::kMaxQubits.
                explicit PauliOrbitSum(int n);

                /// Adds coeff * orbit; merges into an existing term with the same label in O(1).
                void add(const Orbit& orbit, Coeff coeff);

                /// Coefficient of `orbit` in the sum, or 0 if it does not occur.
                Coeff contains(const Orbit& orbit) const;

                /// Drop terms whose coefficient is zero.
                void prune_zeros();

                /// Commutator [this, other] = sum_{i,j} a_i b_j [B_i, B_j].
                PauliOrbitSum commutator(const PauliOrbitSum& other) const;

                /// Product this * other = sum_{i,j} a_i b_j B_i B_j.
                PauliOrbitSum operator*(const PauliOrbitSum& other) const;

                /// Every coefficient multiplied by `scale`; scaling by 0 gives the empty sum.
                PauliOrbitSum operator*(double scale) const;
                PauliOrbitSum operator*(int scale) const { return *this * static_cast<double>(scale); }
                /// Every coefficient multiplied by a complex `scale`; scaling by 0 gives the empty sum.
                PauliOrbitSum operator*(std::complex<double> scale) const {
                        if (scale == std::complex<double>(0.0)) return PauliOrbitSum(n_qubits);
                        PauliOrbitSum out = *this;
                        for (auto& c : out.coeffs) c = c * Coeff(scale);
                        return out;
                }

                /// Coefficient of `target` in [this, other], without building the full commutator.
                Coeff targeted_commutator(const PauliOrbitSum& other, const Orbit& target) const;

                /// The sum expanded into an explicit QubitHamiltonian.
                QubitHamiltonian<Coeff> to_qubit_hamiltonian() const;

                /// Number of terms.
                size_t size() const { return data.size(); }

        private:
                /// label key (p, q, r) -> position in data/coeffs; n_qubits is fixed
                /// per sum, so the label alone identifies an orbit.
                std::unordered_map<int, size_t> index_;

                static int label_key(int p, int q, int r);
};


template <typename Coeff>
PauliOrbitSum<Coeff>::PauliOrbitSum(int n) : n_qubits(n) {
        if (n < 0 || n > Orbit::kMaxQubits) {
                throw std::invalid_argument("PauliOrbitSum: n_qubits must be in [0, "
                        + std::to_string(Orbit::kMaxQubits) + "], got " + std::to_string(n));
        }
}

template <typename Coeff>
int PauliOrbitSum<Coeff>::label_key(int p, int q, int r) {
        // p, q, r <= kMaxQubits, so this packing is collision-free and fits an int.
        constexpr int base = Orbit::kMaxQubits + 1;
        return (p * base + q) * base + r;
}

template <typename Coeff>
void PauliOrbitSum<Coeff>::add(const Orbit& orbit, Coeff coeff) {
        if (orbit.n_qubits != n_qubits) {
                throw std::invalid_argument("PauliOrbitSum.add: n_qubits mismatch ("
                        + std::to_string(n_qubits) + " vs " + std::to_string(orbit.n_qubits) + ")");
        }
        const auto [it, inserted] = index_.try_emplace(label_key(orbit.p, orbit.q, orbit.r), data.size());
        if (!inserted) {
                coeffs[it->second] = coeffs[it->second] + coeff;
                return;
        }
        data.push_back(orbit);
        coeffs.push_back(std::move(coeff));
}

template <typename Coeff>
Coeff PauliOrbitSum<Coeff>::contains(const Orbit& orbit) const {
        if (orbit.n_qubits != n_qubits) {
                return pauliengine_detail::coeff_zero<Coeff>();
        }
        const auto it = index_.find(label_key(orbit.p, orbit.q, orbit.r));
        return it == index_.end() ? pauliengine_detail::coeff_zero<Coeff>() : coeffs[it->second];
}

template <typename Coeff>
PauliOrbitSum<Coeff> PauliOrbitSum<Coeff>::commutator(const PauliOrbitSum& other) const {
        // [ sum_i a_i B_i , sum_j b_j B_j ] = sum_{i,j} a_i b_j [B_i, B_j].
        PauliOrbitSum result(n_qubits);
        for (size_t i = 0; i < data.size(); ++i) {
                for (size_t j = 0; j < other.data.size(); ++j) {
                        const auto part = data[i].commutator(other.data[j]);
                        const Coeff scale = coeffs[i] * other.coeffs[j];
                        for (size_t k = 0; k < part.data.size(); ++k) {
                                result.add(part.data[k], scale * pauliengine_detail::coeff_as<Coeff>(part.coeffs[k]));
                        }
                }
        }
        result.prune_zeros();
        return result;
}

template <typename Coeff>
PauliOrbitSum<Coeff> PauliOrbitSum<Coeff>::operator*(double scale) const {
        if (scale == 0.0) {
                return PauliOrbitSum(n_qubits);
        }
        PauliOrbitSum out = *this;
        for (auto& c : out.coeffs) {
                c *= scale;
        }
        return out;
}

template <typename Coeff>
PauliOrbitSum<Coeff> operator*(double scale, const PauliOrbitSum<Coeff>& s) { return s * scale; }
template <typename Coeff>
PauliOrbitSum<Coeff> operator*(int scale, const PauliOrbitSum<Coeff>& s) { return s * scale; }

template <typename Coeff>
PauliOrbitSum<Coeff> PauliOrbitSum<Coeff>::operator*(const PauliOrbitSum& other) const {
        // ( sum_i a_i B_i )( sum_j b_j B_j ) = sum_{i,j} a_i b_j B_i B_j.
        PauliOrbitSum result(n_qubits);
        for (size_t i = 0; i < data.size(); ++i) {
                for (size_t j = 0; j < other.data.size(); ++j) {
                        const auto part = data[i] * other.data[j];
                        const Coeff scale = coeffs[i] * other.coeffs[j];
                        for (size_t k = 0; k < part.data.size(); ++k) {
                                result.add(part.data[k], scale * pauliengine_detail::coeff_as<Coeff>(part.coeffs[k]));
                        }
                }
        }
        result.prune_zeros();
        return result;
}

template <typename Coeff>
Coeff PauliOrbitSum<Coeff>::targeted_commutator(const PauliOrbitSum& other, const Orbit& target) const {
        Coeff total = pauliengine_detail::coeff_zero<Coeff>();
        for (size_t i = 0; i < data.size(); ++i) {
                for (size_t j = 0; j < other.data.size(); ++j) {
                        total = total + coeffs[i] * other.coeffs[j]
                                * pauliengine_detail::coeff_as<Coeff>(data[i].targeted_commutator(other.data[j], target));
                }
        }
        return total;
}

template <typename Coeff>
QubitHamiltonian<Coeff> PauliOrbitSum<Coeff>::to_qubit_hamiltonian() const {
        std::vector<PauliString<Coeff>> terms;
        for (size_t k = 0; k < data.size(); ++k) {
                const auto part = data[k].to_qubit_hamiltonian();
                for (const auto& ps : part.terms()) {
                        terms.emplace_back(ps.x, ps.y, pauliengine_detail::coeff_as<Coeff>(ps.coeff) * coeffs[k]);
                }
        }
        // different labels never share a string; the constructor only drops zeros
        return QubitHamiltonian<Coeff>(std::move(terms));
}

template <typename Coeff>
void PauliOrbitSum<Coeff>::prune_zeros() {
        std::vector<Orbit> kept_data;
        std::vector<Coeff> kept_coeffs;
        for (size_t i = 0; i < data.size(); ++i) {
                if (!pauliengine_detail::coeff_is_zero(coeffs[i])) {
                        kept_data.push_back(data[i]);
                        kept_coeffs.push_back(coeffs[i]);
                }
        }
        data = std::move(kept_data);
        coeffs = std::move(kept_coeffs);

        index_.clear();
        for (size_t i = 0; i < data.size(); ++i) {
                index_.emplace(label_key(data[i].p, data[i].q, data[i].r), i);
        }
}


// The orbits themselves, so that this header can be included on its own.
#include "pauliengine/PauliOrbit.h"
