#include "pauliengine/QubitHamiltonian.h"


QubitHamiltonian::QubitHamiltonian(const std::vector<PauliString<>>& data) {
            this->data = data;
    }

QubitHamiltonian::QubitHamiltonian(const Hamiltonian_structure& data) {
                        std::vector<PauliString<>> converted_data;
                        converted_data.reserve(data.size());
                        for (Pauli_structure entry : data) {
                                converted_data.push_back(PauliString<>(entry.second, entry.first));
                        }
                        this->data = converted_data;
                }    
                
#ifdef HAVE_SYMENGINE
QubitHamiltonian::QubitHamiltonian(const Hamiltonian_structure_variable& data) {
        std::vector<PauliString<>> converted_data;
        converted_data.reserve(data.size());
        for (Pauli_structure_variable entry : data) {
                converted_data.push_back(PauliString<>(entry.second, entry.first));
        }
        this->data = converted_data;
}

QubitHamiltonian QubitHamiltonian::compact() {
        std::unordered_map<PauliString<>, SymEngine::Expression, PauliStringHash> merged;
        for (const auto& ps : data) {
                merged[ps] = merged[ps] + ps.coeff;
        }
        data.clear();
        for (const auto& [ps, coeff] : merged) {
                PauliString<> new_ps = ps;
                new_ps.coeff = coeff;
                data.push_back(new_ps);
        }
        return QubitHamiltonian(data);
}
#endif


                

QubitHamiltonian QubitHamiltonian::trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state) {
        std::vector<PauliString<>> reduced(this->data.size());
        for (PauliString<>& current : this-> data) {
                PauliString traced_out_string = current.trace_out_qubits(qubits, state);
                if (traced_out_string.coeff != 0.0) {
                        reduced.push_back(traced_out_string);
                }
        }
        return QubitHamiltonian(reduced);
}


std::string QubitHamiltonian::to_string() {
        std::string result = "";
        for (int i = 0; i < this->data.size(); i++) {
                result.append(data[i].to_string() + "\n");
        }    
        return result;   
}
               




                
