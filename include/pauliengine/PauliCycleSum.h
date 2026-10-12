#pragma once

#include <complex>
#include <iterator>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"


class PauliCycle;


/// A sum of PauliCycles, sum_i data[i]. Each cycle carries its own coefficient
/// (in its base string); there is no separate global prefactor. The sum is kept
/// merged: no two summands describe the same rotation orbit, and summands whose
/// coefficient cancelled to zero are dropped.
class PauliCycleSum {
        public:
                using Coeff = std::complex<double>;

                /// The summands; no two lie in the same rotation orbit.
                std::vector<PauliCycle> data;

                /// The empty sum (zero operator).
                PauliCycleSum() = default;

                /// Merges summands from the same rotation orbit and drops zeros.
                explicit PauliCycleSum(std::vector<PauliCycle> cycles);

                /// Defined out-of-line in PauliCycle.h
                QubitHamiltonian<Coeff> to_qubit_hamiltonian() const;

                /// Coefficient of the summand in the same rotation orbit as `cycle`
                /// (its own coefficient is ignored), or 0 if there is none.
                Coeff contains(const PauliCycle& cycle) const;

                /// Commutator [this, other] = sum_{i,j} [c_i, d_j], merged.
                PauliCycleSum commutator(const PauliCycleSum& other) const;

                /// Product this * other = sum_{i,j} c_i d_j, merged.
                PauliCycleSum multiply(const PauliCycleSum& other) const;

                /// Coefficient of `target`'s orbit in [this, other], without building
                /// the full commutator.
                Coeff targeted_commutator(const PauliCycleSum& other, const PauliCycle& target) const;

                /// Number of summands.
                size_t size() const { return data.size(); }
};
