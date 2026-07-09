#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>

#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/Info.h"

#ifdef PAULIENGINE_HAS_OPENMP
#include <omp.h>
#endif

namespace nb = nanobind;
using namespace pauliengine;


NB_MODULE(_core, m) {

        m.doc() = "PauliEngine";
        m.attr("__build_type__") = std::string(build_type());
        m.attr("__compiler_flags__") = compiler_flags();

        // OpenMP diagnostics — check these at runtime to verify the build
        // actually picked up OpenMP. `__openmp__` reflects the compile-time
        // guard, `__omp_max_threads__` reflects what OpenMP will actually use
        // (respects OMP_NUM_THREADS env var).
#ifdef PAULIENGINE_HAS_OPENMP
        m.attr("__openmp__") = true;
        m.def("__omp_max_threads__", []() { return omp_get_max_threads(); },
              "Number of threads OpenMP would use for a parallel region.");
        m.attr("__omp_version__") = _OPENMP;  // yyyymm date, e.g. 201511 = OpenMP 4.5
#else
        m.attr("__openmp__") = false;
        m.def("__omp_max_threads__", []() { return 1; },
              "Sequential build — always returns 1.");
        m.attr("__omp_version__") = 0;
#endif

        // PauliString for complex coefficients
        nb::class_<PauliString<std::complex<double>>>(m, "PauliStringComplex", "Represents a Pauli string in binary symplectic form.")
                .def(nb::init<>(), "Default constructor.")
                .def(nb::init<std::complex<double>, const std::unordered_map<int, std::string>&>(),
                "Constructor from a map of qubit indices to Pauli operators and a complex coefficient.")
                .def(nb::init<const std::pair<std::complex<double>, std::vector<std::pair<char, int>>>>(), "Constructor OpenFermion style and a complex coefficient.")
                .def(nb::init<std::complex<double>&, const std::string&>(),
                "Constructor from a string representation of the Pauli string and a complex coefficient.")
                .def("to_string", &PauliString<std::complex<double>>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("is_all_z", &PauliString<std::complex<double>>::is_all_z, "Checks if the Pauli string consists only of Z operators.")
                .def("get_coeff", &PauliString<std::complex<double>>::get_coeff, "Returns the complex coefficient of the Pauli string.")
                .def("to_dictionary", &PauliString<std::complex<double>>::to_dictionary, "Converts the Pauli string to a dictionary with coefficient and operators.")
                .def("commutator", &PauliString<std::complex<double>>::commutator, "Computes the commutator with another Pauli string.")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &PauliString<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Project onto a computational-basis state of the listed qubits (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &PauliString<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Project onto general single-qubit states (a, b) representing a|0>+b|1> per traced qubit.")
                .def("map_qubits", &PauliString<std::complex<double>>::map_qubits, "Remaps qubit indices according to a given mapping.")
                .def("__mul__", nb::overload_cast<const PauliString<std::complex<double>>&>(&PauliString<std::complex<double>>::operator*, nb::const_), "Multiply two PauliStrings")
                // Type homogeneous multiplication.
                .def("__mul__", nb::overload_cast<double>(&PauliString<std::complex<double>>::operator*, nb::const_), "Scale PauliString by float")
                .def("__mul__", nb::overload_cast<int>(&PauliString<std::complex<double>>::operator*, nb::const_), "Scale PauliString by int")
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<std::complex<double>>::operator*, nb::const_), "Scale PauliString by complex scalar")
                .def("__rmul__", [](const PauliString<std::complex<double>>& ps, double s){ return ps * s; }, "Scale from left by float")
                .def("__rmul__", [](const PauliString<std::complex<double>>& ps, int s){ return ps * s; }, "Scale from left by int")
                .def("__rmul__", [](const PauliString<std::complex<double>>& ps, const std::complex<double> s){ return ps * s; }, "Scale from left by complex")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<std::complex<double>>::operator*=), "In-place scale")
                .def("__imul__", nb::overload_cast<const PauliString<std::complex<double>>&>(&PauliString<std::complex<double>>::operator*=), "In-place multiply")
                // Type heterogeneous multiplication with SymEngine::Expression.
                .def("__mul__", [](const PauliString<std::complex<double>>& ps1, const PauliString<SymEngine::Expression>& ps2){return ps1 * ps2;}, "Multiply PauliString with complex coefficient by PauliString with SymEngine::Expression coefficient")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& ps1, const PauliString<std::complex<double>>& ps2){return ps2 * ps1;}, "Multiply PauliString with complex coefficient by PauliString with SymEngine::Expression coefficient")
                .def("__repr__",  &PauliString<std::complex<double>>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__str__",  &PauliString<std::complex<double>>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__eq__", &PauliString<std::complex<double>>::operator==, "Checks if 2 PauliStrings have same data and Coefficient")
                .def("__neg__", [](const PauliString<std::complex<double>>& ps){ return ps * -1; }, "Negate PauliString")
                .def("__add__", [](const PauliString<std::complex<double>>& ps1, const PauliString<std::complex<double>>& ps2){
                        return QubitHamiltonian<std::complex<double>>({ps1, ps2});
                }, "Add two PauliStringComplex to get QubitHamiltonianComplex")
                .def("__add__", [](const PauliString<std::complex<double>>& ps1, const PauliString<SymEngine::Expression>& ps2){
                        SymEngine::Expression sym_coeff = ps1.coeff;
                        PauliString<SymEngine::Expression> ps1_sym(ps1.x, ps1.y, sym_coeff);
                        return QubitHamiltonian<SymEngine::Expression>({ps1_sym, ps2});
                }, "Add PauliStringComplex and PauliStringSymbolic to get QubitHamiltonianSymbolic")
                .def("__sub__", [](const PauliString<std::complex<double>>& ps1, const PauliString<std::complex<double>>& ps2){
                        return QubitHamiltonian<std::complex<double>>({ps1, ps2 * -1});
                }, "Subtract two PauliStringComplex to get QubitHamiltonianComplex")
                .def("__sub__", [](const PauliString<std::complex<double>>& ps1, const PauliString<SymEngine::Expression>& ps2){
                        SymEngine::Expression sym_coeff = ps1.coeff;
                        PauliString<SymEngine::Expression> ps1_sym(ps1.x, ps1.y, sym_coeff);
                        return QubitHamiltonian<SymEngine::Expression>({ps1_sym, ps2 * -1});
                }, "Subtract PauliStringComplex and PauliStringSymbolic to get QubitHamiltonianSymbolic")
                .def("qubits", &PauliString<std::complex<double>>::qubits, "Returns list of qubits")
                .def("equals", &PauliString<std::complex<double>>::equals, "Checks if 2 PauliStrings have same data")
                .def("set_coeff", &PauliString<std::complex<double>>::set_coeff, "Sets coefficient using a SymEngine expression")
                .def("copy", &PauliString<std::complex<double>>::copy, "Create a copy of the PauliString")
                .def("get_pauli_at_index", &PauliString<std::complex<double>>::get_pauli_from_index, "Returns Pauli operator at given qubit index ('I', 'X', 'Y', 'Z').")
                .def("naked", &PauliString<std::complex<double>>::naked, "Returns a copy of the PauliString with coefficient set to 1.")
                .def("size", &PauliString<std::complex<double>>::size, "Number of qubits with non-identity operator.")
                .def("__len__", &PauliString<std::complex<double>>::size, "Number of qubits with non-identity operator.")
                .def("count_y", &PauliString<std::complex<double>>::count_y, "Number of Y operators in the Pauli string.")
                .def("key_openfermion", &PauliString<std::complex<double>>::key_openfermion, "OpenFermion-style key as list of (Pauli, qubit) pairs.")
                .def_ro("x", &PauliString<std::complex<double>>::x, "Returns the x vector of the Pauli string.")
                .def_ro("y", &PauliString<std::complex<double>>::y, "Returns the y vector of the Pauli string.")
                .def_ro("coeff", &PauliString<std::complex<double>>::coeff, "Returns the coefficient of the Pauli string.")
                .def_ro("is_zero", &PauliString<std::complex<double>>::is_zero, "Returns whether the Pauli string is zero.");

        // PauliString for SymEngine::Expression coefficients
        nb::class_<PauliString<SymEngine::Expression>>(m, "PauliStringSymbolic", "Represents a Pauli string in binary symplectic form.")
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
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &PauliString<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Project onto a computational-basis state of the listed qubits (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &PauliString<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Project onto general single-qubit states (a, b) representing a|0>+b|1> per traced qubit.")
                .def("map_qubits", &PauliString<SymEngine::Expression>::map_qubits, "Remaps qubit indices according to a given mapping.")
                // Type homogeneous multiplication.
                .def("__mul__", nb::overload_cast<const PauliString<SymEngine::Expression>&>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "Multiply two PauliStrings")
                .def("__mul__", nb::overload_cast<double>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "Scale PauliString by float")
                .def("__mul__", nb::overload_cast<int>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "Scale PauliString by int")
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "Scale PauliString by complex scalar")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& ps, double s){ return ps * s; }, "Scale from left by float")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& ps, int s){ return ps * s; }, "Scale from left by int")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& ps, const std::complex<double> s){ return ps * s; }, "Scale from left by complex")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<SymEngine::Expression>::operator*=), "In-place scale")
                .def("__imul__", nb::overload_cast<const PauliString<SymEngine::Expression>&>(&PauliString<SymEngine::Expression>::operator*=), "In-place multiply")
                // Type heterogeneous multiplication with complex coefficients.
                .def("__mul__", [](const PauliString<SymEngine::Expression>& ps1, const PauliString<std::complex<double>>& ps2){return ps1 * ps2;}, "Multiply PauliString with SymEngine::Expression coefficient by PauliString with complex coefficient")
                .def("__rmul__", [](const PauliString<std::complex<double>>& ps1, const PauliString<SymEngine::Expression>& ps2){return ps2 * ps1;}, "Multiply PauliString with SymEngine::Expression coefficient by PauliString with complex coefficient")
                .def("__repr__",  &PauliString<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__str__",  &PauliString<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the Pauli string.")
                .def("__eq__", &PauliString<SymEngine::Expression>::operator==, "Checks if 2 PauliStrings have same data and Coefficient")
                .def("__neg__", [](const PauliString<SymEngine::Expression>& ps){ return ps * -1; }, "Negate PauliString")
                .def("__add__", [](const PauliString<SymEngine::Expression>& ps1, const PauliString<SymEngine::Expression>& ps2){
                        return QubitHamiltonian<SymEngine::Expression>({ps1, ps2});
                }, "Add two PauliStringSymbolic to get QubitHamiltonianSymbolic")
                .def("__add__", [](const PauliString<SymEngine::Expression>& ps1, const PauliString<std::complex<double>>& ps2){
                        SymEngine::Expression sym_coeff = ps2.coeff;
                        PauliString<SymEngine::Expression> ps2_sym(ps2.x, ps2.y, sym_coeff);
                        return QubitHamiltonian<SymEngine::Expression>({ps1, ps2_sym});
                }, "Add PauliStringSymbolic and PauliStringComplex to get QubitHamiltonianSymbolic")
                .def("__sub__", [](const PauliString<SymEngine::Expression>& ps1, const PauliString<SymEngine::Expression>& ps2){
                        return QubitHamiltonian<SymEngine::Expression>({ps1, ps2 * -1});
                }, "Subtract two PauliStringSymbolic to get QubitHamiltonianSymbolic")
                .def("__sub__", [](const PauliString<SymEngine::Expression>& ps1, const PauliString<std::complex<double>>& ps2){
                        SymEngine::Expression sym_coeff = ps2.coeff;
                        PauliString<SymEngine::Expression> ps2_sym(ps2.x, ps2.y, sym_coeff);
                        return QubitHamiltonian<SymEngine::Expression>({ps1, ps2_sym * -1});
                }, "Subtract PauliStringSymbolic and PauliStringComplex to get QubitHamiltonianSymbolic")
                .def("diff", &PauliString<SymEngine::Expression>::diff, "Symbolic derivative wrt the named symbol.")
                .def("qubits", &PauliString<SymEngine::Expression>::qubits, "Returns list of qubits")
                .def("equals", &PauliString<SymEngine::Expression>::equals, "Checks if 2 PauliStrings have same data")
                .def("set_coeff", &PauliString<SymEngine::Expression>::set_coeff, "Sets coefficient using a SymEngine expression")
                .def_static("to_complex", nb::overload_cast<const Expression&>(&PauliString<SymEngine::Expression>::to_complex), "Parse SymEngine expression into complex number")
                .def("copy", &PauliString<SymEngine::Expression>::copy, "Create a copy of the PauliString")
                .def("get_pauli_at_index", &PauliString<SymEngine::Expression>::get_pauli_from_index, "Returns Pauli operator at given qubit index ('I', 'X', 'Y', 'Z').")
                .def("naked", &PauliString<SymEngine::Expression>::naked, "Returns a copy of the PauliString with coefficient set to 1.")
                .def("size", &PauliString<SymEngine::Expression>::size, "Number of qubits with non-identity operator.")
                .def("__len__", &PauliString<SymEngine::Expression>::size, "Number of qubits with non-identity operator.")
                .def("count_y", &PauliString<SymEngine::Expression>::count_y, "Number of Y operators in the Pauli string.")
                .def("key_openfermion", &PauliString<SymEngine::Expression>::key_openfermion, "OpenFermion-style key as list of (Pauli, qubit) pairs.")
                .def_ro("x", &PauliString<SymEngine::Expression>::x, "Returns the x vector of the Pauli string.")
                .def_ro("y", &PauliString<SymEngine::Expression>::y, "Returns the y vector of the Pauli string.")
                .def_ro("coeff", &PauliString<SymEngine::Expression>::coeff, "Returns the coefficient of the Pauli string.")
                .def_ro("is_zero", &PauliString<SymEngine::Expression>::is_zero, "Returns whether the Pauli string is zero.");

        // QubitHamiltonian for complex coefficients
        nb::class_<QubitHamiltonian<std::complex<double>>>(m, "QubitHamiltonianComplex", "Represents a Hamiltonian as a sum of Pauli strings.")
                .def(nb::init<const std::vector<PauliString<std::complex<double>>>&>(), "Constructor from a vector of Pauli strings.")
                .def(nb::init<const Hamiltonian_structure<std::complex<double>>&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                .def("__eq__", &QubitHamiltonian<std::complex<double>>::operator==, "Check if two QubitHamiltonians are equal")
                .def("__add__", &QubitHamiltonian<std::complex<double>>::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double> const>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<double>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Scales the Hamiltonian by a float.")
                .def("__mul__", nb::overload_cast<int>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Scales the Hamiltonian by an int.")
                .def("__mul__", nb::overload_cast<const QubitHamiltonian<std::complex<double>>&>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("__rmul__", [](const QubitHamiltonian<std::complex<double>>& qh, std::complex<double> s){ return qh * s; }, "Scale from left by complex")
                .def("__rmul__", [](const QubitHamiltonian<std::complex<double>>& qh, double s){ return qh * s; }, "Scale from left by float")
                .def("__rmul__", [](const QubitHamiltonian<std::complex<double>>& qh, int s){ return qh * s; }, "Scale from left by int")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &QubitHamiltonian<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto computational-basis states (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &QubitHamiltonian<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto arbitrary single-qubit states (a, b) per traced qubit.")
                .def("to_string", &QubitHamiltonian<std::complex<double>>::to_string, "Converts the Hamiltonian to its full matrix representation.")
                .def("to_dictionary", &QubitHamiltonian<std::complex<double>>::to_dictionary, "Return list of (coefficient, {qubit: Pauli}) pairs.")
                .def("__str__",&QubitHamiltonian<std::complex<double>>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.")
                .def("__repr__",&QubitHamiltonian<std::complex<double>>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.")
                .def("__sub__", nb::overload_cast<const QubitHamiltonian<std::complex<double>>&>(&QubitHamiltonian<std::complex<double>>::operator-, nb::const_), "Subtracts two Hamiltonians.")
                .def("__neg__", nb::overload_cast<>(&QubitHamiltonian<std::complex<double>>::operator-, nb::const_), "Negates Hamiltonian.")
                .def("__pow__", &QubitHamiltonian<std::complex<double>>::power, "Raise Hamiltonian to integer power.")
                .def("__len__", &QubitHamiltonian<std::complex<double>>::size, "Number of Pauli string terms.")
                .def("size", &QubitHamiltonian<std::complex<double>>::size, "Number of Pauli string terms.")
                .def("commutator", &QubitHamiltonian<std::complex<double>>::commutator, "Returns commutator of two QubitHamiltonians")
                .def("compact", &QubitHamiltonian<std::complex<double>>::compact, "Merges duplicate operator terms and removes zero-coefficient terms.")
                .def("set_all_coeff", &QubitHamiltonian<std::complex<double>>::set_all_coeff, nb::arg("value"), "Return a copy with every term's coefficient replaced by value.")
                .def("simplify", &QubitHamiltonian<std::complex<double>>::simplify, nb::arg("threshold") = 0.0, "Removes terms whose coefficient magnitude is below threshold.")
                .def("qubits", &QubitHamiltonian<std::complex<double>>::qubits, "Sorted list of qubit indices the Hamiltonian acts on non-trivially.")
                .def("n_qubits", &QubitHamiltonian<std::complex<double>>::n_qubits, "Number of qubits the Hamiltonian acts on non-trivially.")
                .def("is_all_z", &QubitHamiltonian<std::complex<double>>::is_all_z, "True if all terms are products of Z operators only.")
                .def("count_measurements", &QubitHamiltonian<std::complex<double>>::count_measurements, "1 if is_all_z, else number of terms.")
                .def("is_hermitian", &QubitHamiltonian<std::complex<double>>::is_hermitian, "True if all coefficients are real.")
                .def("is_antihermitian", &QubitHamiltonian<std::complex<double>>::is_antihermitian, "True if all coefficients are purely imaginary.")
                .def("dagger", &QubitHamiltonian<std::complex<double>>::dagger, "Hermitian conjugate (complex-conjugates each coefficient).")
                .def("conjugate", &QubitHamiltonian<std::complex<double>>::conjugate, "Complex conjugate (flips sign for each Y, conjugates coefficient).")
                .def("transpose", &QubitHamiltonian<std::complex<double>>::transpose, "Transpose (flips sign for each Y).")
                .def("split", &QubitHamiltonian<std::complex<double>>::split, "Split into (hermitian, anti-hermitian) parts.")
                .def("map_qubits", &QubitHamiltonian<std::complex<double>>::map_qubits, "Remap qubit indices.")
                .def("power", &QubitHamiltonian<std::complex<double>>::power, "Raise Hamiltonian to integer power.")
                .def("paulistrings", &QubitHamiltonian<std::complex<double>>::paulistrings, "Return the list of Pauli string terms.")
                .def("to_list", &QubitHamiltonian<std::complex<double>>::to_list, "Return the list of Pauli string terms.")
                .def("to_matrix", &QubitHamiltonian<std::complex<double>>::to_matrix, nb::arg("ignore_unused_qubits") = true, "Return dense 2**N x 2**N matrix representation.")
                .def_static("zero", &QubitHamiltonian<std::complex<double>>::zero, "Return an empty Hamiltonian (the zero operator).")
                .def_static("unit", &QubitHamiltonian<std::complex<double>>::unit, "Return the identity Hamiltonian (single term, coefficient 1, no operators).");

        // QubitHamiltonian for SymEngine::Expression coefficients
        nb::class_<QubitHamiltonian<SymEngine::Expression>>(m, "QubitHamiltonianSymbolic", "Represents a Hamiltonian as a sum of Pauli strings.")
                .def(nb::init<const std::vector<PauliString<SymEngine::Expression>>&>(), "Constructor from a vector of Pauli strings.")
                .def(nb::init<const Hamiltonian_structure<SymEngine::Expression>&>(), "Constructor from a Hamiltonian structure (coefficient and operator map).")
                .def("__eq__", &QubitHamiltonian<SymEngine::Expression>::operator==, "Check if two QubitHamiltonians are equal")
                .def("__add__", &QubitHamiltonian<SymEngine::Expression>::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double>>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<double>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Scales the Hamiltonian by a float.")
                .def("__mul__", nb::overload_cast<int>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Scales the Hamiltonian by an int.")
                .def("__mul__", nb::overload_cast<const QubitHamiltonian<SymEngine::Expression>&>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("__rmul__", [](const QubitHamiltonian<SymEngine::Expression>& qh, std::complex<double> s){ return qh * s; }, "Scale from left by complex")
                .def("__rmul__", [](const QubitHamiltonian<SymEngine::Expression>& qh, double s){ return qh * s; }, "Scale from left by float")
                .def("__rmul__", [](const QubitHamiltonian<SymEngine::Expression>& qh, int s){ return qh * s; }, "Scale from left by int")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &QubitHamiltonian<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto computational-basis states (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &QubitHamiltonian<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto arbitrary single-qubit states (a, b) per traced qubit.")
                .def("to_string", &QubitHamiltonian<SymEngine::Expression>::to_string, "Converts the Hamiltonian to its full matrix representation.")
                .def("subs", &QubitHamiltonian<SymEngine::Expression>::substitute, "Replace Variable with value")
                .def("to_dictionary", &QubitHamiltonian<SymEngine::Expression>::to_dictionary, "Return list of (coefficient, {qubit: Pauli}) pairs.")
                .def("__str__",&QubitHamiltonian<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.")
                .def("__repr__",&QubitHamiltonian<SymEngine::Expression>::to_string, "Returns a human-readable string representation of the QubitHamiltonian operator.")
                .def("__sub__", nb::overload_cast<const QubitHamiltonian<SymEngine::Expression>&>(&QubitHamiltonian<SymEngine::Expression>::operator-, nb::const_), "Subtracts two Hamiltonians.")
                .def("__neg__", nb::overload_cast<>(&QubitHamiltonian<SymEngine::Expression>::operator-, nb::const_), "Negates Hamiltonian.")
                .def("__pow__", &QubitHamiltonian<SymEngine::Expression>::power, "Raise Hamiltonian to integer power.")
                .def("__len__", &QubitHamiltonian<SymEngine::Expression>::size, "Number of Pauli string terms.")
                .def("size", &QubitHamiltonian<SymEngine::Expression>::size, "Number of Pauli string terms.")
                .def("commutator", &QubitHamiltonian<SymEngine::Expression>::commutator, "Returns commutator of two QubitHamiltonians")
                .def("diff", &QubitHamiltonian<SymEngine::Expression>::diff, "Symbolic derivative wrt the named symbol.")
                .def("compact", &QubitHamiltonian<SymEngine::Expression>::compact, "Merges duplicate operator terms and removes zero-coefficient terms.")
                .def("set_all_coeff", &QubitHamiltonian<SymEngine::Expression>::set_all_coeff, nb::arg("value"), "Return a copy with every term's coefficient replaced by value.")
                .def("simplify", &QubitHamiltonian<SymEngine::Expression>::simplify, nb::arg("threshold") = 0.0, "Removes terms whose coefficient magnitude is below threshold.")
                .def("qubits", &QubitHamiltonian<SymEngine::Expression>::qubits, "Sorted list of qubit indices the Hamiltonian acts on non-trivially.")
                .def("n_qubits", &QubitHamiltonian<SymEngine::Expression>::n_qubits, "Number of qubits the Hamiltonian acts on non-trivially.")
                .def("is_all_z", &QubitHamiltonian<SymEngine::Expression>::is_all_z, "True if all terms are products of Z operators only.")
                .def("count_measurements", &QubitHamiltonian<SymEngine::Expression>::count_measurements, "1 if is_all_z, else number of terms.")
                .def("is_hermitian", &QubitHamiltonian<SymEngine::Expression>::is_hermitian, "True if all coefficients are real.")
                .def("is_antihermitian", &QubitHamiltonian<SymEngine::Expression>::is_antihermitian, "True if all coefficients are purely imaginary.")
                .def("dagger", &QubitHamiltonian<SymEngine::Expression>::dagger, "Hermitian conjugate.")
                .def("conjugate", &QubitHamiltonian<SymEngine::Expression>::conjugate, "Complex conjugate (flips sign for each Y, conjugates coefficient).")
                .def("transpose", &QubitHamiltonian<SymEngine::Expression>::transpose, "Transpose (flips sign for each Y).")
                .def("split", &QubitHamiltonian<SymEngine::Expression>::split, "Split into (hermitian, anti-hermitian) parts.")
                .def("map_qubits", &QubitHamiltonian<SymEngine::Expression>::map_qubits, "Remap qubit indices.")
                .def("power", &QubitHamiltonian<SymEngine::Expression>::power, "Raise Hamiltonian to integer power.")
                .def("paulistrings", &QubitHamiltonian<SymEngine::Expression>::paulistrings, "Return the list of Pauli string terms.")
                .def("to_list", &QubitHamiltonian<SymEngine::Expression>::to_list, "Return the list of Pauli string terms.")
                .def("to_matrix", &QubitHamiltonian<SymEngine::Expression>::to_matrix, nb::arg("ignore_unused_qubits") = true, "Return dense 2**N x 2**N matrix representation.")
                .def_static("zero", &QubitHamiltonian<SymEngine::Expression>::zero, "Return an empty Hamiltonian (the zero operator).")
                .def_static("unit", &QubitHamiltonian<SymEngine::Expression>::unit, "Return the identity Hamiltonian (single term, coefficient 1, no operators).");

        nb::class_<SymEngine::Expression>(m, "Expression")
                .def(nb::init<const std::string&>())
                .def("__str__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                })
                .def("__repr__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                });
}
