#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>

#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/Info.h"

namespace nb = nanobind;
using namespace pauliengine;


NB_MODULE(_core, m) {

        m.doc() = "PauliEngine";
        m.attr("__build_type__") = std::string(build_type());
        m.attr("__compiler_flags__") = compiler_flags();

        // PauliString for complex coefficients
        nb::class_<PauliString<std::complex<double>>>(m, "PauliString", "Represents a Pauli string in binary symplectic form.")
                .def(nb::init<>(), "Default constructor.")
                .def(nb::init<std::complex<double>, const std::unordered_map<int, std::string>&>(),
                "Constructor from a map of qubit indices to Pauli operators and a complex coefficient.")
                .def(nb::init<const std::pair<std::complex<double>, std::vector<std::pair<char, int>>>>(),
                "Constructor OpenFermion style and a complex coefficient.")
                .def("to_string", &PauliString<std::complex<double>>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("is_all_z", &PauliString<std::complex<double>>::is_all_z, "Checks if the Pauli string consists only of Z operators.")
                .def("get_coeff", &PauliString<std::complex<double>>::get_coeff, "Returns the complex coefficient of the Pauli string.")
                .def("to_dictionary", &PauliString<std::complex<double>>::to_dictionary, "Converts the Pauli string to a dictionary with coefficient and operators.")
                .def("commutator", &PauliString<std::complex<double>>::commutator, "Computes the commutator with another Pauli string.")
                .def("map_qubits", &PauliString<std::complex<double>>::map_qubits, "Remaps qubit indices according to a given mapping.")
                .def("__mul__", nb::overload_cast<const PauliString<std::complex<double>>&>(&PauliString<std::complex<double>>::operator*, nb::const_), "Multiply two PauliStrings")
                // TODO: Implement multiplication for PauliString with SymEngine::Expression coeff.
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<std::complex<double>>::operator*), "Scale PauliString by complex scalar")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<std::complex<double>>::operator*=), "In-place scale")
                .def("__imul__", nb::overload_cast<const PauliString<std::complex<double>>&>(&PauliString<std::complex<double>>::operator*=), "In-place multiply")
                .def("__repr__",  &PauliString<std::complex<double>>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__str__",  &PauliString<std::complex<double>>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__eq__", &PauliString<std::complex<double>>::equals, "Checks if 2 PauliStrings have same data and Coefficient")
                .def("diff", &PauliString<std::complex<double>>::key_openfermion, "Returns Paulistring in Openfermion format")
                .def("qubits", &PauliString<std::complex<double>>::qubits, "Returns list of qubits")
                .def("equals", &PauliString<std::complex<double>>::operator==, "Checks if 2 PauliStrings have same data")
                .def("set_coeff", &PauliString<std::complex<double>>::set_coeff, "Sets coefficient using a SymEngine expression")
                .def("copy", &PauliString<std::complex<double>>::copy, "Create a copy of the PauliString")
                .def("get_pauli_at_index", &PauliString<std::complex<double>>::get_pauli_from_index, "Returns list of qubits")
                .def_ro("x", &PauliString<std::complex<double>>::x, "Returns the x vector of the Pauli string.")
                .def_ro("y", &PauliString<std::complex<double>>::y, "Returns the y vector of the Pauli string.")
                .def_ro("coeff", &PauliString<std::complex<double>>::coeff, "Returns the coefficient of the Pauli string.")
                .def_ro("is_zero", &PauliString<std::complex<double>>::is_zero, "Returns whether the Pauli string is zero.");

        // PauliString for SymEngine::Expression coefficients
        nb::class_<PauliString<SymEngine::Expression>>(m, "PauliStringSym", "Represents a Pauli string in binary symplectic form.")
                .def(nb::init<>(), "Default constructor.")
                .def(nb::init<SymEngine::Expression, const std::unordered_map<int, std::string>&>(),
                "Constructor from a map of qubit indices to Pauli operators and a complex coefficient.")
                .def(nb::init<const std::pair<SymEngine::Expression, std::vector<std::pair<char, int>>>>(),
                "Constructor OpenFermion style and a complex coefficient.")
                .def("to_string", &PauliString<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("is_all_z", &PauliString<SymEngine::Expression>::is_all_z, "Checks if the Pauli string consists only of Z operators.")
                .def("get_coeff", &PauliString<SymEngine::Expression>::get_coeff, "Returns the complex coefficient of the Pauli string.")
                .def("to_dictionary", &PauliString<SymEngine::Expression>::to_dictionary, "Converts the Pauli string to a dictionary with coefficient and operators.")
                .def("commutator", &PauliString<SymEngine::Expression>::commutator, "Computes the commutator with another Pauli string.")
                .def("map_qubits", &PauliString<SymEngine::Expression>::map_qubits, "Remaps qubit indices according to a given mapping.")
                .def("__mul__", nb::overload_cast<const PauliString<SymEngine::Expression>&>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "Multiply two PauliStrings")
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<SymEngine::Expression>::operator*), "Scale PauliString by complex scalar")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<SymEngine::Expression>::operator*=), "In-place scale")
                .def("__imul__", nb::overload_cast<const PauliString<SymEngine::Expression>&>(&PauliString<SymEngine::Expression>::operator*=), "In-place multiply")
                .def("__repr__",  &PauliString<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__str__",  &PauliString<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__eq__", &PauliString<SymEngine::Expression>::equals, "Checks if 2 PauliStrings have same data and Coefficient")
                .def("diff", &PauliString<SymEngine::Expression>::key_openfermion, "Returns Paulistring in Openfermion format")
                .def("qubits", &PauliString<SymEngine::Expression>::qubits, "Returns list of qubits")
                .def("equals", &PauliString<SymEngine::Expression>::operator==, "Checks if 2 PauliStrings have same data")
                .def("set_coeff", &PauliString<SymEngine::Expression>::set_coeff, "Sets coefficient using a SymEngine expression")
                .def_static("to_complex", nb::overload_cast<const Expression&>(&PauliString<SymEngine::Expression>::to_complex), "Parse SymEngine expression into complex number")
                .def("copy", &PauliString<SymEngine::Expression>::copy, "Create a copy of the PauliString")
                .def("get_pauli_at_index", &PauliString<SymEngine::Expression>::get_pauli_from_index, "Returns list of qubits")
                .def_ro("x", &PauliString<SymEngine::Expression>::x, "Returns the x vector of the Pauli string.")
                .def_ro("y", &PauliString<SymEngine::Expression>::y, "Returns the y vector of the Pauli string.")
                .def_ro("coeff", &PauliString<SymEngine::Expression>::coeff, "Returns the coefficient of the Pauli string.")
                .def_ro("is_zero", &PauliString<SymEngine::Expression>::is_zero, "Returns whether the Pauli string is zero.");

        // QubitHamiltonian for complex coefficients
        nb::class_<QubitHamiltonian<std::complex<double>>>(m, "QubitHamiltonian", "Represents a Hamiltonian as a sum of Pauli strings.")
                .def(nb::init<const std::vector<PauliString<std::complex<double>>>&>(), "Constructor from a vector of Pauli strings.")
                .def(nb::init<const Hamiltonian_structure<std::complex<double>>&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                .def("__add__", &QubitHamiltonian<std::complex<double>>::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double> const>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<QubitHamiltonian<std::complex<double>> const>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("trace_out_qubits", &QubitHamiltonian<std::complex<double>>::trace_out_qubits, "Traces out specified qubits in given states.")
                .def("to_string", &QubitHamiltonian<std::complex<double>>::to_string, "Converts the Hamiltonian to its full matrix representation.")
                .def("parse_python_format", &QubitHamiltonian<std::complex<double>>::parse_python_format, "Konvertiert mit SymEngine-Koeffizienten")
                .def("__str__",&QubitHamiltonian<std::complex<double>>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.")
                .def("__repr__",&QubitHamiltonian<std::complex<double>>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.");

        // QubitHamiltonian for SymEngine::Expression coefficients
        nb::class_<QubitHamiltonian<SymEngine::Expression>>(m, "QubitHamiltonianSym", "Represents a Hamiltonian as a sum of Pauli strings.")
                .def(nb::init<const std::vector<PauliString<SymEngine::Expression>>&>(), "Constructor from a vector of Pauli strings.")
                .def(nb::init<const Hamiltonian_structure<SymEngine::Expression>&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                .def("__add__", &QubitHamiltonian<SymEngine::Expression>::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double>>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<QubitHamiltonian<SymEngine::Expression> const>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("trace_out_qubits", &QubitHamiltonian<SymEngine::Expression>::trace_out_qubits, "Traces out specified qubits in given states.")
                .def("to_string", &QubitHamiltonian<SymEngine::Expression>::to_string, "Converts the Hamiltonian to its full matrix representation.")
                .def("subs", &QubitHamiltonian<SymEngine::Expression>::substitute, "Replace Variable with value")
                .def("parse_python_format", &QubitHamiltonian<SymEngine::Expression>::parse_python_format, "Konvertiert mit SymEngine-Koeffizienten")
                .def("__str__",&QubitHamiltonian<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.")
                .def("__repr__",&QubitHamiltonian<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.");

        nb::class_<SymEngine::Expression>(m, "Expression")
                .def(nb::init<const std::string&>())  // Konstruktor aus String
                .def("__str__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                })
                .def("__repr__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                });
}
