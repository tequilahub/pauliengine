#pragma once

#include <complex>
#include <iterator>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"


class PauliCycle;



class PauliCycleSum {
        public:
                using Coeff = std::complex<double>;

                std::vector<PauliCycle> data;
                double coeff;

                PauliCycleSum(double coeff, std::vector<PauliCycle>& data)
                        : data(data), coeff(coeff) {}

                // Defined out-of-line in PauliCycle.h
                QubitHamiltonian<Coeff> to_qubit_hamiltonian() const;
};
