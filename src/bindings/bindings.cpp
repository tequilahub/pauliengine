#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>

#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/PauliCycle.h"
#include "pauliengine/PauliCycleSum.h"
#include "pauliengine/PauliOrbit.h"
#include "pauliengine/PauliOrbitSum.h"
#include "pauliengine/Info.h"

#ifdef PAULIENGINE_HAS_OPENMP
#include <omp.h>
#endif

namespace nb = nanobind;
using namespace pauliengine;

// PauliString::get_hash() limbs (little-endian 64-bit) -> arbitrary-size Python int.
static nb::int_ hash_limbs_to_int(const std::vector<uint64_t>& limbs) {
        static const char digits[] = "0123456789abcdef";
        std::string hex = "0";
        for (size_t i = limbs.size(); i-- > 0; ) {
                for (int s = 60; s >= 0; s -= 4) {
                        hex.push_back(digits[(limbs[i] >> s) & 0xF]);
                }
        }
        PyObject* value = PyLong_FromString(hex.c_str(), nullptr, 16);
        if (!value) {
                throw nb::python_error();
        }
        return nb::steal<nb::int_>(value);
}

// Python int -> little-endian 64-bit limbs for PauliString::from_hash().
static std::vector<uint64_t> hash_int_to_limbs(const nb::int_& value) {
        PyObject* hex_obj = PyNumber_ToBase(value.ptr(), 16);
        if (!hex_obj) {
                throw nb::python_error();
        }
        const std::string hex = nb::cast<std::string>(nb::steal(hex_obj));  // "0x..." or "-0x..."
        if (hex[0] == '-') {
                throw nb::value_error("from_hash: hash must be non-negative");
        }
        std::vector<uint64_t> limbs;
        for (size_t end = hex.size(); end > 2; ) {
                const size_t begin = end >= 2 + 16 ? end - 16 : 2;
                limbs.push_back(std::stoull(hex.substr(begin, end - begin), nullptr, 16));
                end = begin;
        }
        return limbs;
}


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
        nb::class_<PauliString<std::complex<double>>>(m, "PauliStringComplex", "A Pauli string with a complex coefficient, stored in binary symplectic form (x and y bit vectors).")
                .def(nb::init<>(), "The zero Pauli string (identity operator, coefficient 0).")
                .def(nb::init<std::complex<double>, const std::unordered_map<int, std::string>&>(),
                "From a complex coefficient and a map {qubit: 'X'/'Y'/'Z'/'I'}.")
                .def(nb::init<const std::pair<std::complex<double>, std::vector<std::pair<char, int>>>>(), "From an OpenFermion-style pair (coefficient, [('X', qubit), ...]).")
                .def(nb::init<std::complex<double>&, const std::string&>(),
                "From a complex coefficient and a dense string such as 'XIZ' (character i acts on qubit i).")
                .def("to_string", &PauliString<std::complex<double>>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("is_all_z", &PauliString<std::complex<double>>::is_all_z, "True if every non-identity operator is Z.")
                .def("get_coeff", &PauliString<std::complex<double>>::get_coeff, "The coefficient of the Pauli string.")
                .def("to_dictionary", &PauliString<std::complex<double>>::to_dictionary, "The Pauli string as a (coefficient, {qubit: 'X'/'Y'/'Z'}) pair.")
                .def("commutator", &PauliString<std::complex<double>>::commutator, "Commutator [self, other] = self*other - other*self (a single Pauli string, zero if they commute).")
                .def("targeted_commutator", &PauliString<std::complex<double>>::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's operator in [self, other] (0 if absent); target's coefficient is ignored.")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &PauliString<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Project onto a computational-basis state of the listed qubits (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &PauliString<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Project onto general single-qubit states (a, b) representing a|0>+b|1> per traced qubit.")
                .def("map_qubits", &PauliString<std::complex<double>>::map_qubits, "Copy with qubit indices remapped by {old: new}; an occupied target that is not itself remapped is swapped.")
                .def("__mul__", nb::overload_cast<const PauliString<std::complex<double>>&>(&PauliString<std::complex<double>>::operator*, nb::const_), "Product of two Pauli strings, including the phase of the Pauli algebra.")
                // Type homogeneous multiplication.
                .def("__mul__", nb::overload_cast<double>(&PauliString<std::complex<double>>::operator*, nb::const_), "The Pauli string scaled by a float.")
                .def("__mul__", nb::overload_cast<int>(&PauliString<std::complex<double>>::operator*, nb::const_), "The Pauli string scaled by a int.")
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<std::complex<double>>::operator*, nb::const_), "The Pauli string scaled by a complex number.")
                .def("__rmul__", [](const PauliString<std::complex<double>>& ps, double s){ return ps * s; }, "The Pauli string scaled by a float (scalar on the left).")
                .def("__rmul__", [](const PauliString<std::complex<double>>& ps, int s){ return ps * s; }, "The Pauli string scaled by a int (scalar on the left).")
                .def("__rmul__", [](const PauliString<std::complex<double>>& ps, const std::complex<double> s){ return ps * s; }, "The Pauli string scaled by a complex number (scalar on the left).")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<std::complex<double>>::operator*=), "Scale in place by a complex number.")
                .def("__imul__", nb::overload_cast<const PauliString<std::complex<double>>&>(&PauliString<std::complex<double>>::operator*=), "Multiply in place by another Pauli string.")
                // Type heterogeneous multiplication with SymEngine::Expression.
                .def("__mul__", [](const PauliString<std::complex<double>>& ps1, const PauliString<SymEngine::Expression>& ps2){return ps1 * ps2;}, "Product with a symbolic Pauli string; the result is symbolic.")
                .def("__rmul__", [](const PauliString<std::complex<double>>& self, const PauliString<SymEngine::Expression>& other){return other * self;}, "Product other * self with a symbolic Pauli string; the result is symbolic.")
                .def("__repr__",  &PauliString<std::complex<double>>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__str__",  &PauliString<std::complex<double>>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__eq__", &PauliString<std::complex<double>>::operator==, "True if both the operator and the coefficient are equal (see equals() for the operator only).")
                .def("__hash__", [](const PauliString<std::complex<double>>& ps) {
                        return static_cast<Py_ssize_t>(ps.hash());
                }, "Hash over the operator part (x, y); ignores the coefficient.")
                .def("__neg__", [](const PauliString<std::complex<double>>& ps){ return ps * -1; }, "The Pauli string with negated coefficient.")
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
                .def("qubits", &PauliString<std::complex<double>>::qubits, "Sorted list of the qubits with a non-identity operator.")
                .def("equals", &PauliString<std::complex<double>>::equals, "True if the operators are equal, ignoring the coefficients.")
                .def("set_coeff", &PauliString<std::complex<double>>::set_coeff, "Replace the coefficient (in place) by a complex number.")
                .def("copy", &PauliString<std::complex<double>>::copy, "An independent copy of the Pauli string.")
                .def("get_pauli_at_index", &PauliString<std::complex<double>>::get_pauli_from_index, "Pauli operator ('I', 'X', 'Y' or 'Z') on the given qubit.")
                .def("naked", &PauliString<std::complex<double>>::naked, "Copy of the Pauli string with coefficient 1.")
                .def("size", &PauliString<std::complex<double>>::size, "Number of qubits with a non-identity operator (the weight).")
                .def("__len__", &PauliString<std::complex<double>>::size, "Number of qubits with a non-identity operator (the weight).")
                .def("count_y", &PauliString<std::complex<double>>::count_y, "Number of Y operators in the Pauli string.")
                .def("get_hash", [](const PauliString<std::complex<double>>& ps) { return hash_limbs_to_int(ps.get_hash()); },
                        "Invertible int key of the operator part (coefficient ignored): sum_q p_q * 4**q with I=0, X=1, Y=2, Z=3.")
                .def_static("from_hash", [](const nb::int_& hash, std::complex<double> coeff) {
                                return PauliString<std::complex<double>>::from_hash(hash_int_to_limbs(hash), coeff);
                        }, nb::arg("hash"), nb::arg("coeff") = std::complex<double>(1.0, 0.0),
                        "Build the Pauli string whose get_hash() equals hash, with the given coefficient.")
                .def("key_openfermion", &PauliString<std::complex<double>>::key_openfermion, "OpenFermion-style key: list of (Pauli, qubit) pairs of the non-identity sites, sorted by qubit.")
                .def_prop_ro("x", [](const PauliString<std::complex<double>>& ps) { return ps.x.to_vector(); }, "x bits as 64-bit words (bit q set for X or Z on qubit q).")
                .def_prop_ro("y", [](const PauliString<std::complex<double>>& ps) { return ps.y.to_vector(); }, "y bits as 64-bit words (bit q set for Y or Z on qubit q).")
                .def_ro("coeff", &PauliString<std::complex<double>>::coeff, "The coefficient of the Pauli string (read-only).")
                .def_prop_ro("is_zero", [](const PauliString<std::complex<double>>& ps) { return ps.coeff == std::complex<double>(0.0, 0.0); },
                        "True if the coefficient is zero.");

        // PauliString for SymEngine::Expression coefficients
        nb::class_<PauliString<SymEngine::Expression>>(m, "PauliStringSymbolic", "A Pauli string with a symbolic (SymEngine) coefficient, stored in binary symplectic form (x and y bit vectors).")
                .def(nb::init<>(), "The zero Pauli string (identity operator, coefficient 0).")
                .def(nb::init<SymEngine::Expression, const std::unordered_map<int, std::string>&>(),
                "From a symbolic coefficient (Expression) and a map {qubit: 'X'/'Y'/'Z'/'I'}.")
                .def(nb::init<const std::pair<SymEngine::Expression, std::vector<std::pair<char, int>>>>(),
                "From an OpenFermion-style pair (Expression, [('X', qubit), ...]).")
                .def("to_string", &PauliString<SymEngine::Expression>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("is_all_z", &PauliString<SymEngine::Expression>::is_all_z, "True if every non-identity operator is Z.")
                .def("get_coeff", &PauliString<SymEngine::Expression>::get_coeff, "The coefficient of the Pauli string.")
                .def("to_dictionary", &PauliString<SymEngine::Expression>::to_dictionary, "The Pauli string as a (coefficient, {qubit: 'X'/'Y'/'Z'}) pair.")
                .def("commutator", &PauliString<SymEngine::Expression>::commutator, "Commutator [self, other] = self*other - other*self (a single Pauli string, zero if they commute).")
                .def("targeted_commutator", &PauliString<SymEngine::Expression>::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's operator in [self, other] (0 if absent); target's coefficient is ignored.")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &PauliString<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Project onto a computational-basis state of the listed qubits (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &PauliString<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Project onto general single-qubit states (a, b) representing a|0>+b|1> per traced qubit.")
                .def("map_qubits", &PauliString<SymEngine::Expression>::map_qubits, "Copy with qubit indices remapped by {old: new}; an occupied target that is not itself remapped is swapped.")
                // Type homogeneous multiplication.
                .def("__mul__", nb::overload_cast<const PauliString<SymEngine::Expression>&>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "Product of two Pauli strings, including the phase of the Pauli algebra.")
                .def("__mul__", nb::overload_cast<double>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "The Pauli string scaled by a float.")
                .def("__mul__", nb::overload_cast<int>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "The Pauli string scaled by a int.")
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<SymEngine::Expression>::operator*, nb::const_), "The Pauli string scaled by a complex number.")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& ps, double s){ return ps * s; }, "The Pauli string scaled by a float (scalar on the left).")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& ps, int s){ return ps * s; }, "The Pauli string scaled by a int (scalar on the left).")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& ps, const std::complex<double> s){ return ps * s; }, "The Pauli string scaled by a complex number (scalar on the left).")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<SymEngine::Expression>::operator*=), "Scale in place by a complex number.")
                .def("__imul__", nb::overload_cast<const PauliString<SymEngine::Expression>&>(&PauliString<SymEngine::Expression>::operator*=), "Multiply in place by another Pauli string.")
                // Type heterogeneous multiplication with complex coefficients.
                .def("__mul__", [](const PauliString<SymEngine::Expression>& ps1, const PauliString<std::complex<double>>& ps2){return ps1 * ps2;}, "Product with a complex Pauli string; the result is symbolic.")
                .def("__rmul__", [](const PauliString<SymEngine::Expression>& self, const PauliString<std::complex<double>>& other){return other * self;}, "Product other * self with a complex Pauli string; the result is symbolic.")
                .def("__repr__",  &PauliString<SymEngine::Expression>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__str__",  &PauliString<SymEngine::Expression>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__eq__", &PauliString<SymEngine::Expression>::operator==, "True if both the operator and the coefficient are equal (see equals() for the operator only).")
                .def("__hash__", [](const PauliString<SymEngine::Expression>& ps) {
                        return static_cast<Py_ssize_t>(ps.hash());
                }, "Hash over the operator part (x, y); ignores the coefficient.")
                .def("__neg__", [](const PauliString<SymEngine::Expression>& ps){ return ps * -1; }, "The Pauli string with negated coefficient.")
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
                .def("diff", &PauliString<SymEngine::Expression>::diff, "Derivative of the coefficient with respect to the named symbol.")
                .def("subs", &PauliString<SymEngine::Expression>::substitute, "Substitute values for symbols, given as {name: complex value}.")
                .def("qubits", &PauliString<SymEngine::Expression>::qubits, "Sorted list of the qubits with a non-identity operator.")
                .def("equals", &PauliString<SymEngine::Expression>::equals, "True if the operators are equal, ignoring the coefficients.")
                .def("set_coeff", &PauliString<SymEngine::Expression>::set_coeff, "Replace the coefficient (in place) by a SymEngine Expression.")
                .def_static("to_complex", nb::overload_cast<const Expression&>(&PauliString<SymEngine::Expression>::to_complex), "Evaluate a fully numeric Expression to a complex number; raises RuntimeError if it still contains symbols.")
                .def("copy", &PauliString<SymEngine::Expression>::copy, "An independent copy of the Pauli string.")
                .def("get_pauli_at_index", &PauliString<SymEngine::Expression>::get_pauli_from_index, "Pauli operator ('I', 'X', 'Y' or 'Z') on the given qubit.")
                .def("naked", &PauliString<SymEngine::Expression>::naked, "Copy of the Pauli string with coefficient 1.")
                .def("size", &PauliString<SymEngine::Expression>::size, "Number of qubits with a non-identity operator (the weight).")
                .def("__len__", &PauliString<SymEngine::Expression>::size, "Number of qubits with a non-identity operator (the weight).")
                .def("count_y", &PauliString<SymEngine::Expression>::count_y, "Number of Y operators in the Pauli string.")
                .def("get_hash", [](const PauliString<SymEngine::Expression>& ps) { return hash_limbs_to_int(ps.get_hash()); },
                        "Invertible int key of the operator part (coefficient ignored): sum_q p_q * 4**q with I=0, X=1, Y=2, Z=3.")
                .def_static("from_hash", [](const nb::int_& hash, const SymEngine::Expression& coeff) {
                                return PauliString<SymEngine::Expression>::from_hash(hash_int_to_limbs(hash), coeff);
                        }, nb::arg("hash"), nb::arg("coeff"),
                        "Build the Pauli string whose get_hash() equals hash, with the given coefficient.")
                .def_static("from_hash", [](const nb::int_& hash) {
                                return PauliString<SymEngine::Expression>::from_hash(hash_int_to_limbs(hash), SymEngine::Expression(1));
                        }, nb::arg("hash"),
                        "Build the Pauli string whose get_hash() equals hash, with coefficient 1.")
                .def("key_openfermion", &PauliString<SymEngine::Expression>::key_openfermion, "OpenFermion-style key: list of (Pauli, qubit) pairs of the non-identity sites, sorted by qubit.")
                .def_prop_ro("x", [](const PauliString<SymEngine::Expression>& ps) { return ps.x.to_vector(); }, "x bits as 64-bit words (bit q set for X or Z on qubit q).")
                .def_prop_ro("y", [](const PauliString<SymEngine::Expression>& ps) { return ps.y.to_vector(); }, "y bits as 64-bit words (bit q set for Y or Z on qubit q).")
                .def_ro("coeff", &PauliString<SymEngine::Expression>::coeff, "The coefficient of the Pauli string (read-only).")
                .def_prop_ro("is_zero", [](const PauliString<SymEngine::Expression>& ps) {
                                const SymEngine::Expression e = SymEngine::expand(ps.coeff.get_basic());
                                const auto& b = *e.get_basic();
                                return SymEngine::is_a_Number(b) && SymEngine::down_cast<const SymEngine::Number&>(b).is_zero();
                        }, "True if the coefficient is (symbolically) zero.");

        // QubitHamiltonian for complex coefficients
        nb::class_<QubitHamiltonian<std::complex<double>>>(m, "QubitHamiltonianComplex", "A qubit Hamiltonian: a sum of Pauli strings with complex coefficients. Kept compacted: no two terms share an operator and zero terms are dropped.")
                .def(nb::init<const std::vector<PauliString<std::complex<double>>>&>(), "From a list of Pauli strings (equal operators are merged).")
                .def(nb::init<const Hamiltonian_structure<std::complex<double>>&>(), "From a list of (coefficient, {qubit: 'X'/'Y'/'Z'}) pairs.")
                .def("__eq__", &QubitHamiltonian<std::complex<double>>::operator==, "True if both have the same terms with equal coefficients (order-independent).")
                .def("contains", &QubitHamiltonian<std::complex<double>>::contains, nb::arg("pauli_string"),
                        "True if a term with this Pauli operator is present (coefficient ignored).")
                .def("__contains__", &QubitHamiltonian<std::complex<double>>::contains, nb::arg("pauli_string"),
                        "True if a term with this Pauli operator is present (coefficient ignored).")
                .def("__add__", &QubitHamiltonian<std::complex<double>>::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double> const>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<double>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Scales the Hamiltonian by a float.")
                .def("__mul__", nb::overload_cast<int>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Scales the Hamiltonian by an int.")
                .def("__mul__", nb::overload_cast<const QubitHamiltonian<std::complex<double>>&>(&QubitHamiltonian<std::complex<double>>::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("__rmul__", [](const QubitHamiltonian<std::complex<double>>& qh, std::complex<double> s){ return qh * s; }, "The Hamiltonian scaled by a complex number (scalar on the left).")
                .def("__rmul__", [](const QubitHamiltonian<std::complex<double>>& qh, double s){ return qh * s; }, "The Hamiltonian scaled by a float (scalar on the left).")
                .def("__rmul__", [](const QubitHamiltonian<std::complex<double>>& qh, int s){ return qh * s; }, "The Hamiltonian scaled by an int (scalar on the left).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &QubitHamiltonian<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto computational-basis states (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &QubitHamiltonian<std::complex<double>>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto arbitrary single-qubit states (a, b) per traced qubit.")
                .def("to_string", &QubitHamiltonian<std::complex<double>>::to_string, "Human-readable form, one term per line.")
                .def("to_dictionary", &QubitHamiltonian<std::complex<double>>::to_dictionary, "Return list of (coefficient, {qubit: Pauli}) pairs.")
                .def("__str__",&QubitHamiltonian<std::complex<double>>::to_string, "Human-readable form, one term per line.")
                .def("__repr__",&QubitHamiltonian<std::complex<double>>::to_string, "Human-readable form, one term per line.")
                .def("__sub__", nb::overload_cast<const QubitHamiltonian<std::complex<double>>&>(&QubitHamiltonian<std::complex<double>>::operator-, nb::const_), "Subtracts two Hamiltonians.")
                .def("__neg__", nb::overload_cast<>(&QubitHamiltonian<std::complex<double>>::operator-, nb::const_), "Negates Hamiltonian.")
                .def("__pow__", &QubitHamiltonian<std::complex<double>>::power, "Raise Hamiltonian to integer power.")
                .def("__len__", &QubitHamiltonian<std::complex<double>>::size, "Number of Pauli string terms.")
                .def("size", &QubitHamiltonian<std::complex<double>>::size, "Number of Pauli string terms.")
                .def("commutator", &QubitHamiltonian<std::complex<double>>::commutator, "Commutator [self, other] = self*other - other*self.")
                .def("targeted_commutator", &QubitHamiltonian<std::complex<double>>::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of the target Pauli string in the commutator, via a hash-join, O(min(len(A), len(B))).")
                .def("compact", &QubitHamiltonian<std::complex<double>>::compact, "Merge terms with equal operators and drop zero terms, in place; returns the result. Every operation already does this.")
                .def("set_all_coeff", &QubitHamiltonian<std::complex<double>>::set_all_coeff, nb::arg("value"), "Return a copy with every term's coefficient replaced by value.")
                .def("simplify", &QubitHamiltonian<std::complex<double>>::simplify, nb::arg("threshold") = 0.0, "Copy without the terms whose abs(coefficient) <= threshold.")
                .def("qubits", &QubitHamiltonian<std::complex<double>>::qubits, "Sorted list of qubit indices the Hamiltonian acts on non-trivially.")
                .def("n_qubits", &QubitHamiltonian<std::complex<double>>::n_qubits, "Number of qubits the Hamiltonian acts on non-trivially.")
                .def("is_all_z", &QubitHamiltonian<std::complex<double>>::is_all_z, "True if all terms are products of Z operators only.")
                .def("count_measurements", &QubitHamiltonian<std::complex<double>>::count_measurements, "1 if is_all_z(), else the number of terms.")
                .def("is_hermitian", &QubitHamiltonian<std::complex<double>>::is_hermitian, "True if every coefficient is real (abs(imag) <= 1e-12).")
                .def("is_antihermitian", &QubitHamiltonian<std::complex<double>>::is_antihermitian, "True if every coefficient is purely imaginary (abs(real) <= 1e-12).")
                .def("dagger", &QubitHamiltonian<std::complex<double>>::dagger, "Hermitian conjugate: complex-conjugates each coefficient (Pauli strings are Hermitian).")
                .def("conjugate", &QubitHamiltonian<std::complex<double>>::conjugate, "Complex conjugate: conjugates each coefficient and flips its sign for an odd number of Y.")
                .def("transpose", &QubitHamiltonian<std::complex<double>>::transpose, "Transpose: flips the sign of terms with an odd number of Y.")
                .def("split", &QubitHamiltonian<std::complex<double>>::split, "Split into (hermitian, anti-hermitian) parts; needs numeric coefficients.")
                .def("map_qubits", &QubitHamiltonian<std::complex<double>>::map_qubits, "Copy with qubit indices remapped by {old: new} (see PauliString.map_qubits).")
                .def("power", &QubitHamiltonian<std::complex<double>>::power, "Raise Hamiltonian to integer power.")
                .def("paulistrings", &QubitHamiltonian<std::complex<double>>::paulistrings, "Return the list of Pauli string terms.")
                .def("to_list", &QubitHamiltonian<std::complex<double>>::to_list, "The list of Pauli string terms (same as paulistrings()).")
                .def("to_matrix", &QubitHamiltonian<std::complex<double>>::to_matrix, nb::arg("ignore_unused_qubits") = true, "Dense 2**N x 2**N matrix (list of rows). With ignore_unused_qubits, N counts only the qubits the Hamiltonian acts on; otherwise qubit indices are absolute. Needs numeric coefficients.")
                .def_static("zero", &QubitHamiltonian<std::complex<double>>::zero, "Return an empty Hamiltonian (the zero operator).")
                .def_static("unit", &QubitHamiltonian<std::complex<double>>::unit, "Return the identity Hamiltonian (single term, coefficient 1, no operators).");

        // QubitHamiltonian for SymEngine::Expression coefficients
        nb::class_<QubitHamiltonian<SymEngine::Expression>>(m, "QubitHamiltonianSymbolic", "A qubit Hamiltonian: a sum of Pauli strings with symbolic (SymEngine) coefficients. Kept compacted: no two terms share an operator and zero terms are dropped.")
                .def(nb::init<const std::vector<PauliString<SymEngine::Expression>>&>(), "From a list of Pauli strings (equal operators are merged).")
                .def(nb::init<const Hamiltonian_structure<SymEngine::Expression>&>(), "From a list of (coefficient, {qubit: 'X'/'Y'/'Z'}) pairs.")
                .def("__eq__", &QubitHamiltonian<SymEngine::Expression>::operator==, "True if both have the same terms with equal coefficients (order-independent).")
                .def("contains", &QubitHamiltonian<SymEngine::Expression>::contains, nb::arg("pauli_string"),
                        "True if a term with this Pauli operator is present (coefficient ignored).")
                .def("__contains__", &QubitHamiltonian<SymEngine::Expression>::contains, nb::arg("pauli_string"),
                        "True if a term with this Pauli operator is present (coefficient ignored).")
                .def("__add__", &QubitHamiltonian<SymEngine::Expression>::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double>>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<double>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Scales the Hamiltonian by a float.")
                .def("__mul__", nb::overload_cast<int>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Scales the Hamiltonian by an int.")
                .def("__mul__", nb::overload_cast<const QubitHamiltonian<SymEngine::Expression>&>(&QubitHamiltonian<SymEngine::Expression>::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("__rmul__", [](const QubitHamiltonian<SymEngine::Expression>& qh, std::complex<double> s){ return qh * s; }, "The Hamiltonian scaled by a complex number (scalar on the left).")
                .def("__rmul__", [](const QubitHamiltonian<SymEngine::Expression>& qh, double s){ return qh * s; }, "The Hamiltonian scaled by a float (scalar on the left).")
                .def("__rmul__", [](const QubitHamiltonian<SymEngine::Expression>& qh, int s){ return qh * s; }, "The Hamiltonian scaled by an int (scalar on the left).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &QubitHamiltonian<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto computational-basis states (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &QubitHamiltonian<SymEngine::Expression>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto arbitrary single-qubit states (a, b) per traced qubit.")
                .def("to_string", &QubitHamiltonian<SymEngine::Expression>::to_string, "Human-readable form, one term per line.")
                .def("subs", &QubitHamiltonian<SymEngine::Expression>::substitute, "Substitute values for symbols, given as {name: complex value}.")
                .def("to_dictionary", &QubitHamiltonian<SymEngine::Expression>::to_dictionary, "Return list of (coefficient, {qubit: Pauli}) pairs.")
                .def("__str__",&QubitHamiltonian<SymEngine::Expression>::to_string, "Human-readable form, one term per line.")
                .def("__repr__",&QubitHamiltonian<SymEngine::Expression>::to_string, "Human-readable form, one term per line.")
                .def("__sub__", nb::overload_cast<const QubitHamiltonian<SymEngine::Expression>&>(&QubitHamiltonian<SymEngine::Expression>::operator-, nb::const_), "Subtracts two Hamiltonians.")
                .def("__neg__", nb::overload_cast<>(&QubitHamiltonian<SymEngine::Expression>::operator-, nb::const_), "Negates Hamiltonian.")
                .def("__pow__", &QubitHamiltonian<SymEngine::Expression>::power, "Raise Hamiltonian to integer power.")
                .def("__len__", &QubitHamiltonian<SymEngine::Expression>::size, "Number of Pauli string terms.")
                .def("size", &QubitHamiltonian<SymEngine::Expression>::size, "Number of Pauli string terms.")
                .def("commutator", &QubitHamiltonian<SymEngine::Expression>::commutator, "Commutator [self, other] = self*other - other*self.")
                .def("targeted_commutator", &QubitHamiltonian<SymEngine::Expression>::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of the target Pauli string in the commutator, via a hash-join, O(min(len(A), len(B))).")
                .def("diff", &QubitHamiltonian<SymEngine::Expression>::diff, "Derivative of every coefficient with respect to the named symbol.")
                .def("compact", &QubitHamiltonian<SymEngine::Expression>::compact, "Merge terms with equal operators and drop zero terms, in place; returns the result. Every operation already does this.")
                .def("set_all_coeff", &QubitHamiltonian<SymEngine::Expression>::set_all_coeff, nb::arg("value"), "Return a copy with every term's coefficient replaced by value.")
                .def("simplify", &QubitHamiltonian<SymEngine::Expression>::simplify, nb::arg("threshold") = 0.0, "Copy without the terms whose abs(coefficient) <= threshold.")
                .def("qubits", &QubitHamiltonian<SymEngine::Expression>::qubits, "Sorted list of qubit indices the Hamiltonian acts on non-trivially.")
                .def("n_qubits", &QubitHamiltonian<SymEngine::Expression>::n_qubits, "Number of qubits the Hamiltonian acts on non-trivially.")
                .def("is_all_z", &QubitHamiltonian<SymEngine::Expression>::is_all_z, "True if all terms are products of Z operators only.")
                .def("count_measurements", &QubitHamiltonian<SymEngine::Expression>::count_measurements, "1 if is_all_z(), else the number of terms.")
                .def("is_hermitian", &QubitHamiltonian<SymEngine::Expression>::is_hermitian, "True if every coefficient is real (abs(imag) <= 1e-12).")
                .def("is_antihermitian", &QubitHamiltonian<SymEngine::Expression>::is_antihermitian, "True if every coefficient is purely imaginary (abs(real) <= 1e-12).")
                .def("dagger", &QubitHamiltonian<SymEngine::Expression>::dagger, "Hermitian conjugate: complex-conjugates each coefficient (Pauli strings are Hermitian).")
                .def("conjugate", &QubitHamiltonian<SymEngine::Expression>::conjugate, "Complex conjugate: conjugates each coefficient and flips its sign for an odd number of Y.")
                .def("transpose", &QubitHamiltonian<SymEngine::Expression>::transpose, "Transpose: flips the sign of terms with an odd number of Y.")
                .def("split", &QubitHamiltonian<SymEngine::Expression>::split, "Split into (hermitian, anti-hermitian) parts; needs numeric coefficients.")
                .def("map_qubits", &QubitHamiltonian<SymEngine::Expression>::map_qubits, "Copy with qubit indices remapped by {old: new} (see PauliString.map_qubits).")
                .def("power", &QubitHamiltonian<SymEngine::Expression>::power, "Raise Hamiltonian to integer power.")
                .def("paulistrings", &QubitHamiltonian<SymEngine::Expression>::paulistrings, "Return the list of Pauli string terms.")
                .def("to_list", &QubitHamiltonian<SymEngine::Expression>::to_list, "The list of Pauli string terms (same as paulistrings()).")
                .def("to_matrix", &QubitHamiltonian<SymEngine::Expression>::to_matrix, nb::arg("ignore_unused_qubits") = true, "Dense 2**N x 2**N matrix (list of rows). With ignore_unused_qubits, N counts only the qubits the Hamiltonian acts on; otherwise qubit indices are absolute. Needs numeric coefficients.")
                .def_static("zero", &QubitHamiltonian<SymEngine::Expression>::zero, "Return an empty Hamiltonian (the zero operator).")
                .def_static("unit", &QubitHamiltonian<SymEngine::Expression>::unit, "Return the identity Hamiltonian (single term, coefficient 1, no operators).");

        nb::class_<SymEngine::Expression>(m, "Expression", "A SymEngine expression, used as symbolic coefficient.")
                .def(nb::init<const std::string&>(), "Parse an expression such as 'a*b + 2'.")
                .def("__str__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                }, "The expression as a string.")
                .def("__repr__", [](const SymEngine::Expression &e) {
                        return SymEngine::str(*e.get_basic());
                }, "The expression as a string.");

        // PauliCycle (complex coefficients only)
        nb::class_<PauliCycle>(m, "PauliCycleComplex",
                "A Pauli string together with all of its cyclic rotations on n_qubits qubits, "
                "C(P) = coeff * sum_k R^k P. The sum runs over all n rotations, so a base with a "
                "shorter period counts each distinct rotation several times.")
                .def(nb::init<int, const std::unordered_map<int, std::string>&, std::complex<double>>(),
                        nb::arg("n_qubits"), nb::arg("data"), nb::arg("coeff") = std::complex<double>(1.0, 0.0),
                        "From a qubit count, a base {qubit: 'X'/'Y'/'Z'} and a coefficient. Raises ValueError "
                        "if n_qubits < 1 or the base acts on a qubit >= n_qubits.")
                .def(nb::init<int, const std::string&, std::complex<double>>(),
                        nb::arg("n_qubits"), nb::arg("pauli_string"), nb::arg("coeff") = std::complex<double>(1.0, 0.0),
                        "From a qubit count, a dense base string such as 'XIY' and a coefficient.")
                .def(nb::init<int, const PauliString<std::complex<double>>&>(),
                        nb::arg("n_qubits"), nb::arg("base"),
                        "From a qubit count and a base PauliString (its coefficient becomes the cycle's).")
                .def("unify", &PauliCycle::unify,
                        "Rotate the base to the unified orbit representative. Optional: ==/hash are "
                        "rotation-invariant on their own; this only changes the stored base.")
                .def("rotate", &PauliCycle::rotate, nb::arg("k"),
                        "The stored base shifted by k qubits, new[(i + k) mod n] = old[i] (k may be negative or >= n_qubits).")
                .def("rotations", &PauliCycle::rotations,
                        "All n_qubits rotations rotate(0), ..., rotate(n - 1) of the stored base.")
                .def("to_qubit_hamiltonian", &PauliCycle::to_qubit_hamiltonian,
                        "The cycle expanded into a QubitHamiltonian (the sum of all n rotations).")
                .def("commutator", &PauliCycle::commutator, nb::arg("other"),
                        "Commutator [self, other] as a PauliCycleSum (Eq. 40). Only shifts at which the supports "
                        "overlap contribute, so it is cheap for bounded weight. Raises ValueError on an n_qubits mismatch.")
                .def("targeted_commutator", &PauliCycle::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's orbit in [self, other], without building the full commutator; "
                        "target's coefficient is ignored.")
                .def("multiply", &PauliCycle::multiply, nb::arg("other"),
                        "Product self * other as a PauliCycleSum (at most n_qubits summands after merging).")
                .def("__mul__", &PauliCycle::multiply, nb::arg("other"),
                        "Product self * other as a PauliCycleSum (at most n_qubits summands after merging).")
                .def_ro("n_qubits", &PauliCycle::n_qubits, "Number of qubits (read-only).")
                .def_ro("base", &PauliCycle::Base, "The stored base Pauli string (as constructed, or unified after unify()); carries the coefficient.")
                .def("__repr__", [](const PauliCycle& pc) {
                        return "PauliCycleComplex(n_qubits=" + std::to_string(pc.n_qubits)
                                + ", base=" + pc.Base.to_string() + ")";
                }, "PauliCycleComplex(n_qubits=..., base=...).")
                .def("__eq__", [](const PauliCycle& a, const PauliCycle& b) { return a == b; },
                        "True if both describe the same rotation orbit on the same n_qubits; coefficients are ignored.")
                .def("__ne__", [](const PauliCycle& a, const PauliCycle& b) { return a != b; }, "Negation of ==.")
                .def("__hash__", [](const PauliCycle& pc) {
                        return static_cast<Py_ssize_t>(pc.hash());
                }, "Fixed-size hash, consistent with == (rotation-invariant, coefficient ignored).")
                .def("get_hash", [](const PauliCycle& pc) { return hash_limbs_to_int(pc.get_hash()); },
                        "Invertible int key of the cycle (coefficient ignored, identical for all rotations): "
                        "PauliString.get_hash() of the unified representative * 2**32 + n_qubits.")
                .def_static("from_hash", [](const nb::int_& hash, std::complex<double> coeff) {
                                return PauliCycle::from_hash(hash_int_to_limbs(hash), coeff);
                        }, nb::arg("hash"), nb::arg("coeff") = std::complex<double>(1.0, 0.0),
                        "Build the cycle whose get_hash() equals hash, with the given coefficient "
                        "(its base is the unified representative).");

        // PauliCycleSum (complex coefficients only)
        nb::class_<PauliCycleSum>(m, "PauliCycleSumComplex",
                "A sum of PauliCycles; each cycle carries its own coefficient. "
                "Summands from the same rotation orbit are always merged.")
                .def(nb::init<std::vector<PauliCycle>>(), nb::arg("data"),
                "Construct from a list of PauliCycles (cycles in the same rotation orbit are merged, zeros dropped).")
                .def("to_qubit_hamiltonian", &PauliCycleSum::to_qubit_hamiltonian,
                        "The sum expanded into a QubitHamiltonian.")
                .def("commutator", &PauliCycleSum::commutator, nb::arg("other"),
                        "Commutator [self, other] = sum_ij [c_i, d_j], merged.")
                .def("targeted_commutator", &PauliCycleSum::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's orbit in [self, other], without building the full commutator.")
                .def("multiply", &PauliCycleSum::multiply, nb::arg("other"),
                        "Product self * other = sum_ij c_i d_j, merged.")
                .def("__mul__", &PauliCycleSum::multiply, nb::arg("other"),
                        "Product self * other = sum_ij c_i d_j, merged.")
                .def("contains", &PauliCycleSum::contains, nb::arg("cycle"),
                        "Coefficient of the summand in the same rotation orbit as cycle (0 if absent); "
                        "the coefficient of cycle itself is ignored.")
                .def("__contains__", [](const PauliCycleSum& s, const PauliCycle& c) {
                                return s.contains(c) != PauliCycleSum::Coeff(0.0, 0.0);
                        }, nb::arg("cycle"),
                        "True if a summand lies in the same rotation orbit as cycle.")
                .def("__len__", &PauliCycleSum::size, "Number of summands.")
                .def_ro("data", &PauliCycleSum::data, "The summands; no two lie in the same rotation orbit.");

        // PauliOrbit (permutation-invariant basis, arXiv:2604.16701 Section V)
        nb::class_<PauliOrbit>(m, "PauliOrbit",
                "A permutation-invariant Pauli orbit B_{p,q,r} on n qubits: the normalized sum "
                "(i / n_terms) * sum of all Pauli strings with p X-, q Y- and r Z-sites.")
                .def(nb::init<int, int, int, int>(),
                        nb::arg("n_qubits"), nb::arg("p"), nb::arg("q"), nb::arg("r"),
                        "The orbit with p X-, q Y- and r Z-sites on n_qubits qubits. Raises ValueError unless "
                        "0 <= n_qubits <= max_qubits, p, q, r >= 0 and p + q + r <= n_qubits.")
                .def("identities", &PauliOrbit::identities, "Number of identity sites, n - p - q - r.")
                .def("is_valid", &PauliOrbit::is_valid, "Whether p, q, r >= 0 and p + q + r <= n_qubits; always True, since the constructor rejects anything else.")
                .def("n_terms", &PauliOrbit::n_terms, "Number of distinct Pauli strings in the orbit, n! / (p! q! r! (n-p-q-r)!) (Eq. 52).")
                .def("commutator", &PauliOrbit::commutator, nb::arg("other"),
                        "Commutator [self, other] as a PauliOrbitSum (Theorem 1), without expanding the strings. "
                        "Raises ValueError on an n_qubits mismatch.")
                .def("targeted_commutator", &PauliOrbit::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target in [self, other] (Corollary 5 / Algorithm 2), as a complex number "
                        "with zero imaginary part; O(1) for bounded weights.")
                .def("multiply", &PauliOrbit::multiply, nb::arg("other"),
                        "Product self * other as a PauliOrbitSum (complex coefficients).")
                .def("__mul__", &PauliOrbit::multiply, nb::arg("other"),
                        "Product self * other as a PauliOrbitSum (complex coefficients).")
                .def_ro("n_qubits", &PauliOrbit::n_qubits, "Number of qubits (read-only).")
                .def_ro("p", &PauliOrbit::p, "Number of X sites (read-only).")
                .def_ro("q", &PauliOrbit::q, "Number of Y sites (read-only).")
                .def_ro("r", &PauliOrbit::r, "Number of Z sites (read-only).")
                .def_ro_static("max_qubits", &PauliOrbit::kMaxQubits,
                        "Largest supported n_qubits (float64 factorials overflow beyond it).")
                .def("__eq__", [](const PauliOrbit& a, const PauliOrbit& b) { return a == b; },
                        "True if n_qubits and the label (p, q, r) agree.")
                .def("__ne__", [](const PauliOrbit& a, const PauliOrbit& b) { return a != b; }, "Negation of ==.")
                .def("__hash__", [](const PauliOrbit& o) { return static_cast<Py_ssize_t>(o.hash()); },
                        "Fixed-size hash, consistent with ==.")
                .def("get_hash", &PauliOrbit::get_hash,
                        "Invertible int key of the orbit label: n * 2**24 + p * 2**16 + q * 2**8 + r.")
                .def_static("from_hash", [](long long hash) {
                                if (hash < 0) {
                                        throw nb::value_error("from_hash: hash must be non-negative");
                                }
                                return PauliOrbit::from_hash(static_cast<uint64_t>(hash));
                        }, nb::arg("hash"),
                        "Build the orbit whose get_hash() equals hash.")
                .def("__repr__", [](const PauliOrbit& o) {
                        return "PauliOrbit(n_qubits=" + std::to_string(o.n_qubits)
                                + ", p=" + std::to_string(o.p)
                                + ", q=" + std::to_string(o.q)
                                + ", r=" + std::to_string(o.r) + ")";
                }, "PauliOrbit(n_qubits=..., p=..., q=..., r=...).");

        // PauliOrbitSum
        nb::class_<PauliOrbitSum>(m, "PauliOrbitSum",
                "A weighted sum of Pauli orbits, sum_i coeffs[i] * data[i], on a fixed number of qubits.")
                .def(nb::init<int>(), nb::arg("n_qubits"),
                        "The empty sum on n_qubits qubits. Raises ValueError unless 0 <= n_qubits <= PauliOrbit.max_qubits.")
                .def("add", &PauliOrbitSum::add, nb::arg("orbit"), nb::arg("coeff"),
                        "Add coeff * orbit in place, merging into an existing term with the same label. "
                        "Raises ValueError if orbit.n_qubits differs.")
                .def("prune_zeros", &PauliOrbitSum::prune_zeros, "Drop terms whose coefficient is exactly zero, in place.")
                .def("commutator", &PauliOrbitSum::commutator, nb::arg("other"),
                        "Commutator [self, other] = sum_ij a_i b_j [B_i, B_j].")
                .def("targeted_commutator", &PauliOrbitSum::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target in [self, other], without building the full commutator.")
                .def("multiply", &PauliOrbitSum::multiply, nb::arg("other"),
                        "Product self * other = sum_ij a_i b_j B_i B_j.")
                .def("__mul__", &PauliOrbitSum::multiply, nb::arg("other"),
                        "Product self * other = sum_ij a_i b_j B_i B_j.")
                .def("contains", &PauliOrbitSum::contains, nb::arg("orbit"),
                        "Coefficient of orbit in the sum (0 if absent).")
                .def("__contains__", [](const PauliOrbitSum& s, const PauliOrbit& o) {
                                return s.contains(o) != PauliOrbitSum::Coeff(0.0, 0.0);
                        }, nb::arg("orbit"),
                        "True if orbit occurs in the sum with a non-zero coefficient.")
                .def_ro("n_qubits", &PauliOrbitSum::n_qubits, "Number of qubits (read-only).")
                .def_ro("data", &PauliOrbitSum::data, "The orbits; each label occurs at most once.")
                .def_ro("coeffs", &PauliOrbitSum::coeffs, "The list of coefficients, parallel to data.")
                .def("__len__", [](const PauliOrbitSum& s) { return s.size(); }, "Number of terms.");
}
