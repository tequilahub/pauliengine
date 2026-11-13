#pragma once

#include <complex>
#include <iterator>
#include <vector>

#include "pauliengine/Coefficients.h"
#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"


template <typename Coeff>
class PauliCycle;


/// A sum of PauliCycles, sum_i data[i]. Each cycle carries its own coefficient
/// (in its base string); there is no separate global prefactor. The sum is kept
/// merged: no two summands describe the same rotation orbit, and summands whose
/// coefficient cancelled to zero are dropped. Coeff is std::complex<double> or
/// SymEngine::Expression.
template <typename Coeff>
class PauliCycleSum {
        public:
                /// The summands; no two lie in the same rotation orbit.
                std::vector<PauliCycle<Coeff>> data;

                /// The empty sum (zero operator).
                PauliCycleSum() = default;

                /// Merges summands from the same rotation orbit and drops zeros.
                explicit PauliCycleSum(std::vector<PauliCycle<Coeff>> cycles);

                /// Defined out-of-line in PauliCycle.h
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
};
