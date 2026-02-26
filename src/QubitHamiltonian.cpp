#include "pauliengine/QubitHamiltonian.h"


QubitHamiltonian::QubitHamiltonian(const std::vector<PauliString<>>& data) {
            this->data = data;
 }

QubitHamiltonian::QubitHamiltonian(const Hamiltonian_structure& data) {
        std::vector<PauliString<>> converted_data;
        converted_data.reserve(data.size());
        for (Pauli_structure entry : data) {
                converted_data.push_back(PauliString<>(entry.first, entry.second));
        }
        this->data = converted_data;
}    
                
#ifdef HAVE_SYMENGINE
QubitHamiltonian::QubitHamiltonian(const Hamiltonian_structure_variable& data) {
        std::vector<PauliString<>> converted_data;
        converted_data.reserve(data.size());
        for (Pauli_structure_variable entry : data) {
                converted_data.push_back(PauliString<>(entry.first, entry.second));
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
               
                
Matrix2D QubitHamiltonian::get_pauli_matrix(const std::string& p){static const std::complex<double> I(0,1);
        if (p == "I") return {{1,0},{0,1}};
        else if (p == "X") return {{0,1},{1,0}};
        else if (p == "Y") return {{0,-I},{I,0}};
        else if (p == "Z") return {{1,0},{0,-1}};
        else throw std::runtime_error("Unbekannter Pauli-Operator");
}


Matrix2D QubitHamiltonian::tensor_product(const Matrix2D& A, const Matrix2D& B)
{
        int rows = A.size() * B.size();
        int cols = A[0].size() * B[0].size();
        Matrix2D result(rows, std::vector<std::complex<double>>(cols, 0));

        for (int i=0; i<(int)A.size(); ++i) {
                for (int j=0; j<(int)A[0].size(); ++j) {
                        for (int k=0; k<(int)B.size(); ++k) {
                                for (int l=0; l<(int)B[0].size(); ++l) {
                                        result[i*B.size()+k][j*B[0].size()+l] = A[i][j] * B[k][l];
                                }
                        }
                }
        }
        return result;
}

Matrix2D QubitHamiltonian::pauli_string_to_matrix(const std::unordered_map<int, std::string>& pauli_map, int num_qubits){
        Matrix2D result = get_pauli_matrix("I");
        for (int q=0; q < num_qubits; ++q) {
                std::string op = "I";
                auto it = pauli_map.find(q);
                if (it != pauli_map.end()) op = it->second;
                Matrix2D pmat = get_pauli_matrix(op);
                result = tensor_product(result, pmat);
        }
        return result;
}