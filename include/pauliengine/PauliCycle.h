#pragma once

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
                PauliString<Coeff> Base;

                PauliCycle(int qubits, const std::unordered_map<int, std::string>& data, Coeff coeff = Coeff(1.0, 0.0))
                        : n_qubits(qubits), Base(coeff, data) {}

                PauliCycle(int qubits, const std::string& pauli_string, Coeff coeff = Coeff(1.0, 0.0))
                        : n_qubits(qubits) {

                        Base = PauliString<Coeff>(coeff, pauli_string);
                }

                PauliCycle(int qubits, const PauliString<Coeff>& base)
                        : n_qubits(qubits), Base(base) {}


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
                //TODO: gleiche cycles erkennen.

        private:
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
