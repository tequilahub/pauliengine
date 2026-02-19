#pragma once
#include <unordered_map>
#include <string>
#include <complex>
#include <vector>

#include "pauliengine/PauliString.h"

#ifdef HAVE_SYMENGINE
#include <symengine/expression.h>
using SymEngine::Expression;
#else
#endif



#define BITS_IN_INTEGER (sizeof(uint64_t) * 8)


using Pauli_structure = std::pair<std::complex<double>, std::unordered_map<int, std::string>>;
#ifdef HAVE_SYMENGINE
using Pauli_structure_variable = std::pair<SymEngine::Expression, std::unordered_map<int, std::string>>;
using Hamiltonian_structure_variable = std::vector<Pauli_structure_variable>;
#endif

using Hamiltonian_structure = std::vector<Pauli_structure>;
using Matrix2D = std::vector<std::vector<std::complex<double>>>;

const std::complex<double> neg_I(0.0 , -1.0);


class QubitHamiltonian{
    public:
    std::vector<PauliString<>> data;

    QubitHamiltonian(const std::vector<PauliString<>>& data); 
    QubitHamiltonian(const Hamiltonian_structure& data);
    QubitHamiltonian trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state);
    std::string to_string();
    
    
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
                // return QubitHamiltonian(data).compact();
                return QubitHamiltonian(data);
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
            for (const auto &ps : data) {
                    temp_data.push_back(ps.substitute(substitution_map));
                }
                return QubitHamiltonian(temp_data);
        }
        
        #ifdef HAVE_SYMENGINE

        QubitHamiltonian(const Hamiltonian_structure_variable& data);

        QubitHamiltonian compact();

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
    static Matrix2D get_pauli_matrix(const std::string& p);
    static Matrix2D tensor_product(const Matrix2D& A, const Matrix2D& B);
    static Matrix2D pauli_string_to_matrix(const std::unordered_map<int, std::string>& pauli_map, int num_qubits);
};