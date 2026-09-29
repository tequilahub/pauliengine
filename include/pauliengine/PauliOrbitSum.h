#pragma once

#include <complex>
#include <cstddef>
#include <vector>

// A scalar-weighted sum of Pauli orbits, sum_i coeffs[i] * data[i]
// (arXiv:2604.16701, Section V). See PauliOrbit.h for the orbit data structure.
//
// Coefficients are complex: commutators of the (anti-Hermitian) orbits have real
// structure constants, but ordinary products carry a full complex Pauli phase
// (see PauliOrbit::multiply), so both are representable here.

class PauliOrbit;

class PauliOrbitSum {
        public:
                using Coeff = std::complex<double>;

                int n_qubits;
                std::vector<PauliOrbit> data;
                std::vector<Coeff> coeffs;

                explicit PauliOrbitSum(int n) : n_qubits(n) {}

                // Defined out-of-line in PauliOrbit.h, where PauliOrbit is complete.
                void add(const PauliOrbit& orbit, Coeff coeff);

                void prune_zeros();

                PauliOrbitSum commutator(const PauliOrbitSum& other) const;


                PauliOrbitSum multiply(const PauliOrbitSum& other) const;


                Coeff targeted_commutator(const PauliOrbitSum& other, const PauliOrbit& target) const;

                size_t size() const { return data.size(); }
};
