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
        
        nb::class_<PauliString<>>(m, "PauliString", "Represents a Pauli string in binary symplectic form.")
                .def(nb::init<>(), "Default constructor.")
                .def(nb::init<const std::unordered_map<int, std::string>&, std::complex<double>>(),
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
                .def("__imul__", nb::overload_cast<const PauliString<>&>(&PauliString<>::operator*=), "In-place multiply")
                .def("__repr__",  &PauliString<>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__str__",  &PauliString<>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__eq__", &PauliString<>::equals, "Checks if 2 PauliStrings have same data and Coefficient")
                .def("diff", &PauliString<>::key_openfermion, "Returns Paulistring in Openfermion format")
                .def("qubits", &PauliString<>::qubits, "Returns list of qubits")
                .def("equals", &PauliString<>::operator==, "Checks if 2 PauliStrings have same data")
                .def("set_coeff", &PauliString<>::set_coeff, "Sets coefficient using a SymEngine expression")
                .def_static("to_complex", &PauliString<>::to_complex, "Parse SymEngine expression into complex number")
                .def("copy", &PauliString<>::copy, "Create a copy of the PauliString")
                .def("get_pauli_at_index", &PauliString<>::get_pauli_from_index, "Returns list of qubits")
                .def_ro("x", &PauliString<>::x, "Returns the x vector of the Pauli string.")
                .def_ro("y", &PauliString<>::y, "Returns the y vector of the Pauli string.")
                .def_ro("coeff", &PauliString<>::coeff, "Returns the coefficient of the Pauli string.")
                .def_ro("is_zero", &PauliString<>::is_zero, "Returns whether the Pauli string is zero.");

        nb::class_<QubitHamiltonian>(m, "QubitHamiltonian", "Represents a Hamiltonian as a sum of Pauli strings.")
                .def(nb::init<const std::vector<PauliString<>>&>(), "Constructor from a vector of Pauli strings.")
                .def(nb::init<const Hamiltonian_structure&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                #ifdef HAVE_SYMENGINE
                .def(nb::init<const Hamiltonian_structure_variable&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                #endif
                .def("__add__", &QubitHamiltonian::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double> const>(&QubitHamiltonian::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<QubitHamiltonian const>(&QubitHamiltonian::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("trace_out_qubits", &QubitHamiltonian::trace_out_qubits, "Traces out specified qubits in given states.")
                .def("to_string", &QubitHamiltonian::to_string, "Converts the Hamiltonian to its full matrix representation.")
                .def("diff", &QubitHamiltonian::diff, "Differentiation after given Variable")
                .def("subs", &QubitHamiltonian::substitute, "Replace Variable with value")
                .def("parse_python_format", &QubitHamiltonian::parse_python_format, "Konvertiert mit SymEngine-Koeffizienten")
                .def("__str__",&QubitHamiltonian::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.")
                .def("__repr__",&QubitHamiltonian::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.");


        #ifdef HAVE_SYMENGINE
        nb::class_<SymEngine::Expression>(m, "Expression")
                .def(nb::init<const std::string&>())  // Konstruktor aus String
                .def("__str__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                });
        #endif
}