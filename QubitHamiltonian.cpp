#include "QubitHamiltonian.h"
#include "PauliString.h"

class QubitHamiltonian {

        public:

                std::vector<PauliString<>> data;

                QubitHamiltonian(const std::vector<PauliString<>>& data) {
                        this->data = data;
                }

                QubitHamiltonian(const Hamiltonian_structure& data) {
                        std::vector<PauliString<>> converted_data;
                        converted_data.reserve(data.size());
                        for (Pauli_structure entry : data) {
                                converted_data.push_back(PauliString<>(entry.second, entry.first));
                        }
                        this->data = converted_data;
                }

                QubitHamiltonian(const Hamiltonian_structure_variable& data) {
                        std::vector<PauliString<>> converted_data;
                        converted_data.reserve(data.size());
                        for (Pauli_structure_variable entry : data) {
                                converted_data.push_back(PauliString<>(entry.second, entry.first));
                        }
                        this->data = converted_data;
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
                        return QubitHamiltonian(data);
                }
                // Fehlerbehandlung im Wrapper
                QubitHamiltonian trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state) {
                        std::vector<PauliString<>> reduced(this->data.size());
                        for (PauliString<>& current : this-> data) {
                                PauliString traced_out_string = current.trace_out_qubits(qubits, state);
                                if (traced_out_string.coeff != 0.0) {
                                        reduced.push_back(traced_out_string);
                                }
                        }
                        return QubitHamiltonian(reduced);
                }
                /* TODO: Wenn expression noch variablen haben, nicht ausführen
                Matrix2D hamiltonian_to_matrix(const QubitHamiltonian& H_input, int num_qubits) {
                        int dim = 1 << num_qubits;
                        Hamiltonian_structure H = H_input.parse_python_format();
                        Matrix2D mat(dim, std::vector<std::complex<double>>(dim, 0));

                        for (const auto& [coef, pauli_map] : H) {
                                Matrix2D term_mat = pauli_string_to_matrix(pauli_map, num_qubits);
                                for (int i=0; i<dim; ++i) {
                                        for (int j=0; j<dim; ++j) {
                                                mat[i][j] += coef * term_mat[i][j];
                                        }
                                }
                        }
                        return mat;
                }
                */
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

                std::string to_string() {
                        std::string result = "";
                        for (int i = 0; i < this->data.size(); i++) {
                                result.append(data[i].to_string() + "\n");
                        }    
                        return result;   
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
                /*
                void compact() {
                        std::unordered_map<PauliString<>, std::complex<double>, PauliStringHash> merged;
                        for (const auto& ps : data) {
                                merged[ps] += ps.coeff;
                        }
                        data.clear();
                        for (const auto& [ps, coeff] : merged) {
                                PauliString<> new_ps = ps;
                                new_ps.coeff = coeff;
                                data.push_back(new_ps);
                        }
                }
                */

        private:
                
                

                static Matrix2D get_pauli_matrix(const std::string& p) {
                        static const std::complex<double> I(0,1);
                        if (p == "I") return {{1,0},{0,1}};
                        else if (p == "X") return {{0,1},{1,0}};
                        else if (p == "Y") return {{0,-I},{I,0}};
                        else if (p == "Z") return {{1,0},{0,-1}};
                        else throw std::runtime_error("Unbekannter Pauli-Operator");
                }

                static Matrix2D tensor_product(const Matrix2D& A, const Matrix2D& B) {
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

                static Matrix2D pauli_string_to_matrix(const std::unordered_map<int, std::string>& pauli_map, int num_qubits) {
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

NB_MODULE(PauliEngine, m) {
        nb::class_<PauliString<>>(m, "PauliString", "Represents a Pauli string in binary symplectic form.")
                .def(nb::init<>(), "Default constructor.")
                .def(nb::init<const std::unordered_map<int, std::string>&, std::complex<double>>(),
                "Constructor from a map of qubit indices to Pauli operators and a complex coefficient.")
                .def(nb::init<const std::unordered_map<int, std::string>&, std::string>(),
                "Constructor from a map of qubit indices to Pauli operators and a complex coefficient.")
                .def("to_string", &PauliString<>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("is_all_z", &PauliString<>::is_all_z, "Checks if the Pauli string consists only of Z operators.")
                .def("get_coeff", &PauliString<>::get_coeff, "Returns the complex coefficient of the Pauli string.")
                .def("to_dictionary", &PauliString<>::to_dictionary, "Converts the Pauli string to a dictionary with coefficient and operators.")
                .def("commutator", &PauliString<>::commutator, "Computes the commutator with another Pauli string.")
                .def("map_qubits", &PauliString<>::map_qubits, "Remaps qubit indices according to a given mapping.")
                .def("diff", &PauliString<>::diff, "Differentiation after given Variable")
                .def("__mul__", nb::overload_cast<const PauliString<>&>(&PauliString<>::operator*, nb::const_), "Multiply two PauliStrings")
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<>::operator*), "Scale PauliString by complex scalar")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<>::operator*=), "In-place scale")
                .def("__imul__", nb::overload_cast<PauliString<>>(&PauliString<>::operator*=), "In-place multiply")
                .def("__repr__",  &PauliString<>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__str__",  &PauliString<>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("equals", &PauliString<>::equals, "Checks if 2 PauliStrings have same data and Coefficient")
                .def("diff", &PauliString<>::key_openfermion, "Returns Paulistring in Openfermion format")
                .def("qubits", &PauliString<>::qubits, "Returns list of qubits")
                .def("__eq__", &PauliString<>::operator==, "Checks if 2 PauliStrings have same data")
                .def("set_coeff", nb::overload_cast<SymEngine::Expression>(&PauliString<>::set_coeff), "Sets coefficient using a SymEngine expression")
                .def("set_coeff", nb::overload_cast<std::complex<double>>(&PauliString<>::set_coeff), "Sets coefficient using a complex number")
                .def_static("to_complex", &PauliString<>::to_complex, "Parse SymEngine expression into complex number")
                .def("copy", &PauliString<>::copy, "Create a copy of the PauliString")
                .def("get_pauli_at_index", &PauliString<>::get_pauli_from_index, "Returns list of qubits");

        nb::class_<QubitHamiltonian>(m, "QubitHamiltonian", "Represents a Hamiltonian as a sum of Pauli strings.")
                .def(nb::init<const std::vector<PauliString<>>&>(), "Constructor from a vector of Pauli strings.")
                .def(nb::init<const Hamiltonian_structure&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                .def(nb::init<const Hamiltonian_structure_variable&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                .def("__add__", &QubitHamiltonian::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double> const>(&QubitHamiltonian::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<QubitHamiltonian const>(&QubitHamiltonian::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("trace_out_qubits", &QubitHamiltonian::trace_out_qubits, "Traces out specified qubits in given states.")
                //.def("hamiltonian_to_matrix", &QubitHamiltonian::hamiltonian_to_matrix, "Converts the Hamiltonian to its full matrix representation.")
                .def("to_string", &QubitHamiltonian::to_string, "Converts the Hamiltonian to its full matrix representation.")
                .def("diff", &QubitHamiltonian::diff, "Differentiation after given Variable")
                .def("subs", &QubitHamiltonian::substitute, "Replace Variable with value")
                .def("parse_python_format", &QubitHamiltonian::parse_python_format, "Konvertiert mit SymEngine-Koeffizienten");

                
        //m.def("to_complex", &PauliString::to_complex, "Parse SymEngine expression into complex number");

        nb::class_<SymEngine::Expression>(m, "Expression")
                .def(nb::init<const std::string&>())  // Konstruktor aus String
                .def("__str__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                });
}