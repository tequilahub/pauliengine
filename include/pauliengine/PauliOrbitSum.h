#pragma once

#include <complex>
#include <cstddef>
#include <unordered_map>
#include <vector>

class PauliOrbit;

/// A scalar-weighted sum of Pauli orbits, sum_i coeffs[i] * data[i]
/// (arXiv:2604.16701, Section V). See PauliOrbit.h for the orbit data structure.
///
/// Coefficients are complex: commutators of the (anti-Hermitian) orbits have real
/// structure constants, but ordinary products carry a full complex Pauli phase
/// (see PauliOrbit::multiply), so both are representable here.
class PauliOrbitSum {
        public:
                using Coeff = std::complex<double>;

                /// Fixed at construction (exposed read-only to Python).
                int n_qubits;
                std::vector<PauliOrbit> data;   ///< the orbits, each label at most once
                std::vector<Coeff> coeffs;      ///< coefficients, parallel to data

                /// Empty sum on n qubits. Throws std::invalid_argument unless 0 <= n <= PauliOrbit::kMaxQubits.
                explicit PauliOrbitSum(int n);

                /// Defined out-of-line in PauliOrbit.h, where PauliOrbit is complete.
                /// Merges into an existing term with the same label in O(1).
                void add(const PauliOrbit& orbit, Coeff coeff);

                /// Coefficient of `orbit` in the sum, or 0 if it does not occur.
                Coeff contains(const PauliOrbit& orbit) const;

                /// Drop terms whose coefficient is exactly zero.
                void prune_zeros();

                /// Commutator [this, other] = sum_{i,j} a_i b_j [B_i, B_j].
                PauliOrbitSum commutator(const PauliOrbitSum& other) const;


                /// Product this * other = sum_{i,j} a_i b_j B_i B_j.
                PauliOrbitSum multiply(const PauliOrbitSum& other) const;


                /// Coefficient of `target` in [this, other], without building the full commutator.
                Coeff targeted_commutator(const PauliOrbitSum& other, const PauliOrbit& target) const;

                /// Number of terms.
                size_t size() const { return data.size(); }

        private:
                /// label key (p, q, r) -> position in data/coeffs; n_qubits is fixed
                /// per sum, so the label alone identifies an orbit.
                std::unordered_map<int, size_t> index_;

                static int label_key(int p, int q, int r);
};
