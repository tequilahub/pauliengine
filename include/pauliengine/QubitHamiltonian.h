#pragma once
#include <unordered_map>
#include <string>
#include <complex>
#include <vector>
#include <omp.h>

#include "pauliengine/PauliString.h"

#ifdef HAVE_SYMENGINE
#include <symengine/expression.h>
using SymEngine::Expression;
#else
#endif



#define BITS_IN_INTEGER (sizeof(uint64_t) * 8)


using Pauli_structure = std::pair<std::complex<double>, std::unordered_map<int, std::string>>;
using Pauli_structure_variable = std::pair<SymEngine::Expression, std::unordered_map<int, std::string>>;
using Hamiltonian_structure = std::vector<Pauli_structure>;
using Hamiltonian_structure_variable = std::vector<Pauli_structure_variable>;
using Matrix2D = std::vector<std::vector<std::complex<double>>>;

const std::complex<double> neg_I(0.0 , -1.0);


class QubitHamiltonian{
    public:
    std::vector<PauliString<>> data;

    QubitHamiltonian(const std::vector<PauliString<>>& data); 
    QubitHamiltonian(const Hamiltonian_structure& data);
    QubitHamiltonian(const Hamiltonian_structure_variable& data);
    QubitHamiltonian trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state);
    std::string to_string();

    QubitHamiltonian compact();

    QubitHamiltonian operator*(std::complex<double> scale) const{
                        QubitHamiltonian to_scale = QubitHamiltonian(this->data);
                        for (auto &entry : to_scale.data) {
                                entry.coeff *= scale;
                        }
                        return to_scale;
                }
    QubitHamiltonian operator*(QubitHamiltonian other) const{
                        std::vector<PauliString<>> data;
                        data.reserve(this->data.size() * other.data.size());
                        for (auto &entry_first : this->data) {
                                for (auto &entry_second : other.data) {
                                        data.push_back(entry_first * entry_second);                       
                                }    
                        }
                        return QubitHamiltonian(data).compact();
                }
    QubitHamiltonian operator+(const QubitHamiltonian& other) {
        std::vector<PauliString<>> data;
        QubitHamiltonian first_parsed = QubitHamiltonian(this->data);
        QubitHamiltonian second_parsed = QubitHamiltonian(other.data);
        size_t remaining_in_second = second_parsed.data.size();
        for (auto &entry_first : first_parsed.data) {
                //TODO: Multi-Threading
                for (auto it = second_parsed.data.begin(); it != second_parsed.data.end(); ) {
                        if (entry_first == *it) { 
                                entry_first.coeff += it->coeff; 
                                it = second_parsed.data.erase(it);
                                remaining_in_second--;
                                continue;
                        } else {
                                ++it;
                        }
                }
        }
        data.reserve(first_parsed.data.size() + remaining_in_second);
        for (PauliString first_data : first_parsed.data) {
                data.push_back(first_data);
        }
        for (PauliString remaining : second_parsed.data) {
                data.push_back(remaining);
        }
        return QubitHamiltonian(data);
    }

    QubitHamiltonian diff(std::string symbol) const{
        std::vector<PauliString<>> data;
        for (const PauliString<>& ps : this->data) {
                PauliString<> temp = ps.diff(symbol);
                if (temp.is_zero == false) {
                        data.push_back(temp);
                }
        }
        return QubitHamiltonian(data);
    }

    QubitHamiltonian substitute(const std::unordered_map<std::string, std::complex<double>>& substitution_map) const {
                        std::vector<PauliString<>> temp_data;
                        for (const auto ps : data) {
                                temp_data.push_back(ps.substitute(substitution_map));
                        }
                        return QubitHamiltonian(temp_data);
                }
    
    #ifdef HAVE_SYMENGINE
    Hamiltonian_structure_variable parse_python_format() const{
            Hamiltonian_structure_variable output;
            output.reserve(this->data.size()); 
            for (const auto& entry : this->data) {
                    std::unordered_map<int, std::string> temp;
                    temp.reserve(entry.x.size() * BITS_IN_INTEGER); 
                    for (size_t i = 0; i < entry.x.size(); ++i) {
                            uint64_t x_word = entry.x[i];
                            uint64_t y_word = entry.y[i];
                            for (int j = 0; j < BITS_IN_INTEGER; ++j) {
                                    uint64_t bit = (uint64_t)1 << j;
                                    bool x_bit = x_word & bit;
                                    bool y_bit = y_word & bit;
                                    if (!x_bit && !y_bit) continue;
                                    int index = i * BITS_IN_INTEGER + j;
                                    if (x_bit && y_bit) {
                                            temp[(int)index] = "Z";
                                    } else if (x_bit) {
                                            temp[(int)index] = "X";
                                    } else {
                                            temp[(int)index] = "Y";
                                    }
                            }
                    }

                    output.emplace_back(entry.coeff, std::move(temp));
            }
            return output;
    }
    #else
    Hamiltonian_structure parse_python_format() const{
            Hamiltonian_structure output;
            output.reserve(this->data.size()); 
            for (const auto& entry : this->data) {
                    std::unordered_map<int, std::string> temp;
                    temp.reserve(entry.x.size() * BITS_IN_INTEGER); 
                    for (size_t i = 0; i < entry.x.size(); ++i) {
                            uint64_t x_word = entry.x[i];
                            uint64_t y_word = entry.y[i];
                            for (int j = 0; j < BITS_IN_INTEGER; ++j) {
                                    uint64_t bit = (uint64_t)1 << j;
                                    bool x_bit = x_word & bit;
                                    bool y_bit = y_word & bit;
                                    if (!x_bit && !y_bit) continue;
                                    int index = i * BITS_IN_INTEGER + j;
                                    if (x_bit && y_bit) {
                                            temp[(int)index] = "Z";
                                    } else if (x_bit) {
                                            temp[(int)index] = "X";
                                    } else {
                                            temp[(int)index] = "Y";
                                    }
                            }
                    }

                    output.emplace_back(entry.coeff, std::move(temp));
            }
            return output;
    }
    #endif

    private:
    static Matrix2D get_pauli_matrix(const std::string& p){
        static const std::complex<double> I(0,1);
        if (p == "I") return {{1,0},{0,1}};
        else if (p == "X") return {{0,1},{1,0}};
        else if (p == "Y") return {{0,-I},{I,0}};
        else if (p == "Z") return {{1,0},{0,-1}};
        else throw std::runtime_error("Unbekannter Pauli-Operator");
}
    static Matrix2D tensor_product(const Matrix2D& A, const Matrix2D& B)
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
    static Matrix2D pauli_string_to_matrix(const std::unordered_map<int, std::string>& pauli_map, int num_qubits){
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
};