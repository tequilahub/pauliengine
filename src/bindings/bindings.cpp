#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/optional.h>
#include <nanobind/ndarray.h>

#include "pauliengine/QubitHamiltonian.h"
#include "pauliengine/UnweightedPauliString.h"
#include "pauliengine/UnweightedQubitHamiltonian.h"
#include "pauliengine/PauliCycle.h"
#include "pauliengine/PauliCycleSum.h"
#include "pauliengine/PauliOrbit.h"
#include "pauliengine/PauliOrbitSum.h"
#include "pauliengine/Info.h"
#include "pauliengine/Arithmetic.h"
#include "pauliengine/PauliDecomposition.h"

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


using PSC = PauliString<std::complex<double>>;
using PSS = PauliString<SymEngine::Expression>;
using QHC = QubitHamiltonian<std::complex<double>>;
using QHS = QubitHamiltonian<SymEngine::Expression>;

// Whether a coefficient returned by contains() is non-zero (for `in`).
static bool coeff_nonzero(const std::complex<double>& c) { return c != std::complex<double>(0.0, 0.0); }
static bool coeff_nonzero(const SymEngine::Expression& c) {
        const SymEngine::Expression e = SymEngine::expand(c.get_basic());
        const auto& b = *e.get_basic();
        return !(SymEngine::is_a_Number(b) && SymEngine::down_cast<const SymEngine::Number&>(b).is_zero());
}

// Register a * s and s * a of class A for an int or float s.
template <typename A, typename Cls>
void bind_scale_ops(Cls& cls, const char* doc) {
        cls.def("__mul__", [](const A& a, int s) { return a * s; }, nb::arg("scale"), doc)
                .def("__mul__", [](const A& a, double s) { return a * s; }, nb::arg("scale"), doc)
                .def("__rmul__", [](const A& a, int s) { return a * s; }, nb::arg("scale"), doc)
                .def("__rmul__", [](const A& a, double s) { return a * s; }, nb::arg("scale"), doc);
}

// Arithmetic of an unweighted string or set U with a weighted operator W
// (Arithmetic.h): U counts with coefficient 1.0 * i**k on every element.
template <typename W, typename U, typename ClsW, typename ClsU>
void bind_unweighted_pair(ClsW& cls_w, ClsU& cls_u) {
        cls_w.def("__add__", [](const W& a, const U& b) { return a + b; },
                  "Sum with an unweighted operator, which counts with coefficient 1.0 * i**phase.")
             .def("__mul__", [](const W& a, const U& b) { return a * b; },
                  "Product with an unweighted operator, which counts with coefficient 1.0 * i**phase.")
             .def("commutator", [](const W& a, const U& b) { return commutator(a, b); }, nb::arg("other"),
                  "Commutator with an unweighted operator, which counts with coefficient 1.0 * i**phase.");
        cls_u.def("__add__", [](const U& a, const W& b) { return a + b; },
                  "Sum with another operator; next to a weighted one this counts with coefficient 1.0 * i**phase.")
             .def("__mul__", [](const U& a, const W& b) { return a * b; },
                  "Product with another operator; next to a weighted one this counts with coefficient 1.0 * i**phase.")
             .def("commutator", [](const U& a, const W& b) { return commutator(a, b); }, nb::arg("other"),
                  "Commutator with another operator; next to a weighted one this counts with coefficient 1.0 * i**phase.");
        cls_w.def("__sub__", [](const W& a, const U& b) { return a - b; },
                  "Difference with an unweighted operator, which counts with coefficient 1.0 * i**phase.")
             .def("inner_product", [](const W& a, const U& b) { return inner_product(a, b); }, nb::arg("other"),
                  "Normalized Hilbert-Schmidt inner product with an unweighted operator, which counts with "
                  "coefficient 1.0 * i**phase.");
        cls_u.def("__sub__", [](const U& a, const W& b) { return a - b; },
                  "Difference with a weighted operator; this counts with coefficient 1.0 * i**phase.")
             .def("inner_product", [](const U& a, const W& b) { return inner_product(a, b); }, nb::arg("other"),
                  "Normalized Hilbert-Schmidt inner product with a weighted operator; this counts with "
                  "coefficient 1.0 * i**phase.");
}

// Coefficient lookup of an unweighted string in a Hamiltonian.
template <typename H, typename ClsH>
void bind_contains_unweighted(ClsH& cls_h) {
        using U = unweighted_pauliString;
        cls_h.def("contains", [](const H& a, const U& b) { return a.contains(as_weighted(b)); }, nb::arg("pauli_string"),
                  "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its phase is ignored.")
             .def("__contains__", [](const H& a, const U& b) { return coeff_nonzero(a.contains(as_weighted(b))); }, nb::arg("pauli_string"),
                  "True if the Hamiltonian has a non-zero term with this operator.");
}


using CYC = PauliCycle<std::complex<double>>;
using CYS = PauliCycle<SymEngine::Expression>;
using CSC = PauliCycleSum<std::complex<double>>;
using CSS = PauliCycleSum<SymEngine::Expression>;
using OSC = PauliOrbitSum<std::complex<double>>;
using OSS = PauliOrbitSum<SymEngine::Expression>;

template <typename... Ts>
struct TypeList {};

/// Every operator representation, for the mixed arithmetic below.
using AllOperators = TypeList<PSC, PSS, QHC, QHS, CYC, CYS, CSC, CSS, PauliOrbit, OSC, OSS>;

// Pairs whose operation is bound by the class itself, so the mixed binding is left out.
template <typename A, typename B>
inline constexpr bool own_add = (pauliengine_detail::is_pauli_string<A> && pauliengine_detail::is_pauli_string<B>)
        || (pauliengine_detail::is_hamiltonian<A> && std::is_same_v<A, B>);
template <typename A, typename B>
inline constexpr bool own_mul = (pauliengine_detail::is_pauli_string<A> && pauliengine_detail::is_pauli_string<B>)
        || std::is_same_v<A, B>;
template <typename A, typename B>
inline constexpr bool own_commutator = std::is_same_v<A, B>;

// Mixed arithmetic (see Arithmetic.h): register +, -, *, commutator and
// inner_product of class A with every B, unless A binds the operation itself.
template <typename A, typename Cls, typename... Bs>
void bind_arithmetic(Cls& cls, TypeList<Bs...>) {
        ([&] {
                using B = Bs;
                if constexpr (!own_add<A, B>) {
                        cls.def("__add__", [](const A& a, const B& b) { return a + b; },
                                "Sum with another operator; a QubitHamiltonian unless both are in the same family.");
                        cls.def("__sub__", [](const A& a, const B& b) { return a - b; },
                                "Difference with another operator; a QubitHamiltonian unless both are in the same family.");
                }
                if constexpr (!own_mul<A, B>) {
                        cls.def("__mul__", [](const A& a, const B& b) { return a * b; },
                                "Product with another operator; a QubitHamiltonian unless both are in the same family.");
                }
                if constexpr (!own_commutator<A, B>) {
                        cls.def("commutator", [](const A& a, const B& b) { return commutator(a, b); }, nb::arg("other"),
                                "Commutator [self, other] with another operator; a QubitHamiltonian unless both are in "
                                "the same family.");
                }
                cls.def("inner_product", [](const A& a, const B& b) { return inner_product(a, b); }, nb::arg("other"),
                        "Normalized Hilbert-Schmidt inner product Tr(self^dagger other) / 2**n, so Pauli strings are "
                        "orthonormal (arXiv:2604.16701, Cor. 1 / Lemma 1). Conjugate-linear in self; multiply by "
                        "2**n for the unnormalized trace on n qubits. Orbits are never expanded.");
        }(), ...);
}

// PauliCycle with coefficients of type C.
template <typename C>
auto bind_cycle(nb::module_& m, const char* name, const char* doc) {
        using Cy = PauliCycle<C>;
        const C one = pauliengine_detail::coeff_one<C>();
        auto cls = nb::class_<Cy>(m, name, doc)
                .def(nb::init<int, const std::unordered_map<int, std::string>&, C>(),
                        nb::arg("n_qubits"), nb::arg("data"), nb::arg("coeff") = one,
                        "From a qubit count, a base {qubit: 'X'/'Y'/'Z'} and a coefficient. Raises ValueError "
                        "if n_qubits < 1 or the base acts on a qubit >= n_qubits.")
                .def(nb::init<int, const std::string&, C>(),
                        nb::arg("n_qubits"), nb::arg("pauli_string"), nb::arg("coeff") = one,
                        "From a qubit count, a dense base string such as 'XIY' and a coefficient.")
                .def(nb::init<int, const PauliString<C>&>(),
                        nb::arg("n_qubits"), nb::arg("base"),
                        "From a qubit count and a base PauliString (its coefficient becomes the cycle's).")
                .def("period", &Cy::period,
                        "Number of distinct rotations d of the base (d divides n_qubits); the cycle counts each "
                        "distinct rotation n_qubits / d times.")
                .def("unify", &Cy::unify,
                        "Rotate the base to the unified orbit representative. Optional: ==/hash are "
                        "rotation-invariant on their own; this only changes the stored base.")
                .def("rotate", &Cy::rotate, nb::arg("k"),
                        "The stored base shifted by k qubits, new[(i + k) mod n] = old[i] (k may be negative or >= n_qubits).")
                .def("rotations", &Cy::rotations,
                        "All n_qubits rotations rotate(0), ..., rotate(n - 1) of the stored base.")
                .def("to_qubit_hamiltonian", &Cy::to_qubit_hamiltonian,
                        "The cycle expanded into a QubitHamiltonian (the sum of all n rotations).")
                .def("commutator", &Cy::commutator, nb::arg("other"),
                        "Commutator [self, other] as a PauliCycleSum (Eq. 40). Only shifts at which the supports "
                        "overlap contribute, so it is cheap for bounded weight. Raises ValueError on an n_qubits mismatch.")
                .def("targeted_commutator", [](const Cy& a, const Cy& b, const CYC& t) { return a.targeted_commutator(b, t); },
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's orbit in [self, other], without building the full commutator; "
                        "target's coefficient is ignored.")
                .def("targeted_commutator", [](const Cy& a, const Cy& b, const CYS& t) { return a.targeted_commutator(b, t); },
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's orbit in [self, other], without building the full commutator; "
                        "target's coefficient is ignored.")
                .def("__mul__", [](const Cy& a, const Cy& b) { return a * b; }, nb::arg("other"),
                        "Product self * other as a PauliCycleSum (at most n_qubits summands after merging).")
                .def_ro("n_qubits", &Cy::n_qubits, "Number of qubits (read-only).")
                .def_ro("base", &Cy::Base, "The stored base Pauli string (as constructed, or unified after unify()); carries the coefficient.")
                .def("__repr__", [type = std::string(name)](const Cy& pc) {
                        return type + "(n_qubits=" + std::to_string(pc.n_qubits) + ", base=" + pc.Base.to_string() + ")";
                }, "Class name(n_qubits=..., base=...).")
                .def("__eq__", [](const Cy& a, const CYC& b) { return a == b; },
                        "True if both describe the same rotation orbit on the same n_qubits; coefficients are ignored.")
                .def("__eq__", [](const Cy& a, const CYS& b) { return a == b; },
                        "True if both describe the same rotation orbit on the same n_qubits; coefficients are ignored.")
                .def("__ne__", [](const Cy& a, const CYC& b) { return a != b; }, "Negation of ==.")
                .def("__ne__", [](const Cy& a, const CYS& b) { return a != b; }, "Negation of ==.")
                .def("__hash__", [](const Cy& pc) { return static_cast<Py_ssize_t>(pc.hash()); },
                        "Fixed-size hash, consistent with == (rotation-invariant, coefficient ignored).")
                .def("get_hash", [](const Cy& pc) { return hash_limbs_to_int(pc.get_hash()); },
                        "Invertible int key of the cycle (coefficient ignored, identical for all rotations): "
                        "PauliString.get_hash() of the unified representative * 2**32 + n_qubits.")
                .def_static("from_hash", [](const nb::int_& hash, C coeff) {
                                return Cy::from_hash(hash_int_to_limbs(hash), coeff);
                        }, nb::arg("hash"), nb::arg("coeff") = one,
                        "Build the cycle whose get_hash() equals hash, with the given coefficient "
                        "(its base is the unified representative).");
        return cls;
}

// PauliCycleSum with coefficients of type C.
template <typename C>
auto bind_cycle_sum(nb::module_& m, const char* name, const char* doc) {
        using Cy = PauliCycle<C>;
        using Sum = PauliCycleSum<C>;
        auto cls = nb::class_<Sum>(m, name, doc)
                .def(nb::init<std::vector<Cy>>(), nb::arg("data"),
                        "Construct from a list of cycles (cycles in the same rotation orbit are merged, zeros dropped).")
                .def("to_qubit_hamiltonian", &Sum::to_qubit_hamiltonian,
                        "The sum expanded into a QubitHamiltonian.")
                .def("commutator", &Sum::commutator, nb::arg("other"),
                        "Commutator [self, other] = sum_ij [c_i, d_j], merged.")
                .def("targeted_commutator", [](const Sum& a, const Sum& b, const CYC& t) { return a.targeted_commutator(b, t); },
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's orbit in [self, other], without building the full commutator.")
                .def("targeted_commutator", [](const Sum& a, const Sum& b, const CYS& t) { return a.targeted_commutator(b, t); },
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's orbit in [self, other], without building the full commutator.")
                .def("__mul__", [](const Sum& a, const Sum& b) { return a * b; }, nb::arg("other"),
                        "Product self * other = sum_ij c_i d_j, merged.")
                .def("contains", [](const Sum& s, const CYC& c) { return s.contains(c); }, nb::arg("cycle"),
                        "Coefficient of the summand in the same rotation orbit as cycle (0 if absent); "
                        "the coefficient of cycle itself is ignored.")
                .def("contains", [](const Sum& s, const CYS& c) { return s.contains(c); }, nb::arg("cycle"),
                        "Coefficient of the summand in the same rotation orbit as cycle (0 if absent); "
                        "the coefficient of cycle itself is ignored.")
                .def("__contains__", [](const Sum& s, const CYC& c) { return coeff_nonzero(s.contains(c)); }, nb::arg("cycle"),
                        "True if a summand lies in the same rotation orbit as cycle.")
                .def("__contains__", [](const Sum& s, const CYS& c) { return coeff_nonzero(s.contains(c)); }, nb::arg("cycle"),
                        "True if a summand lies in the same rotation orbit as cycle.")
                .def("unify", &Sum::unify,
                        "Rewrite every summand's base to the unified representative of its orbit, in place. "
                        "Cheap: the unified forms are already computed when the sum is merged.")
                .def("__len__", &Sum::size, "Number of summands.")
                .def_ro("data", &Sum::data, "The summands; no two lie in the same rotation orbit.");
        return cls;
}

// PauliOrbitSum with coefficients of type C.
template <typename C>
auto bind_orbit_sum(nb::module_& m, const char* name, const char* doc) {
        using Sum = PauliOrbitSum<C>;
        auto cls = nb::class_<Sum>(m, name, doc)
                .def(nb::init<int>(), nb::arg("n_qubits"),
                        "The empty sum on n_qubits qubits. Raises ValueError unless 0 <= n_qubits <= PauliOrbit.max_qubits.")
                .def("add", &Sum::add, nb::arg("orbit"), nb::arg("coeff"),
                        "Add coeff * orbit in place, merging into an existing term with the same label. "
                        "Raises ValueError if orbit.n_qubits differs.");
        if constexpr (std::is_same_v<C, SymEngine::Expression>) {
                cls.def("add", [](Sum& s, const PauliOrbit& o, const std::string& coeff) { s.add(o, SymEngine::Expression(coeff)); },
                        nb::arg("orbit"), nb::arg("coeff"),
                        "Add coeff * orbit in place; the coefficient is parsed like 'a*b'.")
                   .def("add", [](Sum& s, const PauliOrbit& o, std::complex<double> coeff) { s.add(o, SymEngine::Expression(coeff)); },
                        nb::arg("orbit"), nb::arg("coeff"),
                        "Add coeff * orbit in place, with a numeric coefficient.");
        }
        cls.def("to_qubit_hamiltonian", &Sum::to_qubit_hamiltonian,
                        "The sum expanded into a QubitHamiltonian.")
                .def("prune_zeros", &Sum::prune_zeros, "Drop terms whose coefficient is zero, in place.")
                .def("commutator", &Sum::commutator, nb::arg("other"),
                        "Commutator [self, other] = sum_ij a_i b_j [B_i, B_j].")
                .def("targeted_commutator", &Sum::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target in [self, other], without building the full commutator.")
                .def("__mul__", [](const Sum& a, const Sum& b) { return a * b; }, nb::arg("other"),
                        "Product self * other = sum_ij a_i b_j B_i B_j.")
                .def("contains", &Sum::contains, nb::arg("orbit"),
                        "Coefficient of orbit in the sum (0 if absent).")
                .def("__contains__", [](const Sum& s, const PauliOrbit& o) { return coeff_nonzero(s.contains(o)); },
                        nb::arg("orbit"), "True if orbit occurs in the sum with a non-zero coefficient.")
                .def_ro("n_qubits", &Sum::n_qubits, "Number of qubits (read-only).")
                .def_ro("data", &Sum::data, "The orbits; each label occurs at most once.")
                .def_ro("coeffs", &Sum::coeffs, "The list of coefficients, parallel to data.")
                .def("__len__", [](const Sum& s) { return s.size(); }, "Number of terms.");
        return cls;
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
        auto cls_psc = nb::class_<PauliString<std::complex<double>>>(m, "PauliStringComplex", "A Pauli string with a complex coefficient, stored in binary symplectic form (x and y bit vectors).")
                .def(nb::init<>(), "The zero Pauli string (identity operator, coefficient 0).")
                .def(nb::init<std::complex<double>, const std::unordered_map<int, std::string>&>(),
                "From a complex coefficient and a map {qubit: 'X'/'Y'/'Z'/'I'}.")
                .def(nb::init<const std::pair<std::complex<double>, std::vector<std::pair<char, int>>>>(), "From an OpenFermion-style pair (coefficient, [('X', qubit), ...]).")
                .def(nb::init<std::complex<double>&, const std::string&>(),
                "From a complex coefficient and a text: dense such as 'XIZ' (character i acts on qubit i) or sparse such as 'X0 Z2'.")
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
        auto cls_pss = nb::class_<PauliString<SymEngine::Expression>>(m, "PauliStringSymbolic", "A Pauli string with a symbolic (SymEngine) coefficient, stored in binary symplectic form (x and y bit vectors).")
                .def(nb::init<>(), "The zero Pauli string (identity operator, coefficient 0).")
                .def(nb::init<SymEngine::Expression, const std::unordered_map<int, std::string>&>(),
                "From a symbolic coefficient (Expression) and a map {qubit: 'X'/'Y'/'Z'/'I'}.")
                .def(nb::init<SymEngine::Expression&, const std::string&>(),
                "From a symbolic coefficient (Expression) and a text: dense such as 'XIZ' (character i acts on qubit i) or sparse such as 'X0 Z2'.")
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
        auto cls_qhc = nb::class_<QubitHamiltonian<std::complex<double>>>(m, "QubitHamiltonianComplex", "A qubit Hamiltonian: a sum of Pauli strings with complex coefficients. Kept compacted: no two terms share an operator and zero terms are dropped.")
                .def(nb::init<const std::vector<PauliString<std::complex<double>>>&>(), "From a list of Pauli strings (equal operators are merged).")
                .def(nb::init<const Hamiltonian_structure<std::complex<double>>&>(), "From a list of (coefficient, {qubit: 'X'/'Y'/'Z'}) pairs.")
                .def("__eq__", &QubitHamiltonian<std::complex<double>>::operator==, "True if both have the same terms with equal coefficients (order-independent).")
                .def("contains", [](const QubitHamiltonian<std::complex<double>>& qh, const PauliString<std::complex<double>>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
                .def("contains", [](const QubitHamiltonian<std::complex<double>>& qh, const PauliString<SymEngine::Expression>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
                .def("contains", [](const QubitHamiltonian<std::complex<double>>& qh, const std::unordered_map<int, std::string>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string {qubit: 'X'} in the Hamiltonian (0 if absent).")
                .def("contains", [](const QubitHamiltonian<std::complex<double>>& qh, const std::vector<std::pair<char, int>>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string [('X', 0), ...] in the Hamiltonian (0 if absent).")
                .def("contains", [](const QubitHamiltonian<std::complex<double>>& qh, const std::string& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string 'XIZ' or 'X0 Z2' in the Hamiltonian (0 if absent).")
                .def("__contains__", [](const QubitHamiltonian<std::complex<double>>& qh, const PauliString<std::complex<double>>& ps) {
                                return coeff_nonzero(qh.contains(ps));
                        }, nb::arg("pauli_string"), "True if a term with this Pauli operator is present (coefficient ignored).")
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
        auto cls_qhs = nb::class_<QubitHamiltonian<SymEngine::Expression>>(m, "QubitHamiltonianSymbolic", "A qubit Hamiltonian: a sum of Pauli strings with symbolic (SymEngine) coefficients. Kept compacted: no two terms share an operator and zero terms are dropped.")
                .def(nb::init<const std::vector<PauliString<SymEngine::Expression>>&>(), "From a list of Pauli strings (equal operators are merged).")
                .def(nb::init<const Hamiltonian_structure<SymEngine::Expression>&>(), "From a list of (coefficient, {qubit: 'X'/'Y'/'Z'}) pairs.")
                .def("__eq__", &QubitHamiltonian<SymEngine::Expression>::operator==, "True if both have the same terms with equal coefficients (order-independent).")
                .def("contains", [](const QubitHamiltonian<SymEngine::Expression>& qh, const PauliString<SymEngine::Expression>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
                .def("contains", [](const QubitHamiltonian<SymEngine::Expression>& qh, const PauliString<std::complex<double>>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
                .def("contains", [](const QubitHamiltonian<SymEngine::Expression>& qh, const std::unordered_map<int, std::string>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string {qubit: 'X'} in the Hamiltonian (0 if absent).")
                .def("contains", [](const QubitHamiltonian<SymEngine::Expression>& qh, const std::vector<std::pair<char, int>>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string [('X', 0), ...] in the Hamiltonian (0 if absent).")
                .def("contains", [](const QubitHamiltonian<SymEngine::Expression>& qh, const std::string& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string 'XIZ' or 'X0 Z2' in the Hamiltonian (0 if absent).")
                .def("__contains__", [](const QubitHamiltonian<SymEngine::Expression>& qh, const PauliString<SymEngine::Expression>& ps) {
                                return coeff_nonzero(qh.contains(ps));
                        }, nb::arg("pauli_string"), "True if a term with this Pauli operator is present (coefficient ignored).")
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

        // Pauli string without a coefficient
        auto cls_up = nb::class_<unweighted_pauliString>(m, "unweighted_pauliString",
                "A Pauli string without a coefficient: i**phase * P with phase in {0, 1, 2, 3}. Use it only as a basis element "
                "of a Lie algebra; for linear combinations use weighted strings. The product of two such strings has an exact "
                "phase; use to_pauli_string() to attach a coefficient.")
                .def(nb::init<>(), "The identity with phase 0.")
                .def(nb::init<const std::string&, int>(), nb::arg("paulis"), nb::arg("phase") = 0,
                        "From a text, dense such as 'XIZ' or sparse such as 'X0 Z2'; phase is the exponent k of i**k.")
                .def(nb::init<const std::unordered_map<int, std::string>&, int>(), nb::arg("paulis"), nb::arg("phase") = 0,
                        "From a map {qubit: 'X'/'Y'/'Z'/'I'}; phase is the exponent k of i**k.")
                .def(nb::init<const std::vector<std::pair<char, int>>&, int>(), nb::arg("paulis"), nb::arg("phase") = 0,
                        "From a list of (Pauli, qubit) pairs; phase is the exponent k of i**k.")
                .def(nb::init<const PSC&, int>(), nb::arg("string"), nb::arg("phase") = 0,
                        "The operator part of a Pauli string; its coefficient is dropped.")
                .def(nb::init<const PSS&, int>(), nb::arg("string"), nb::arg("phase") = 0,
                        "The operator part of a symbolic Pauli string; its coefficient is dropped.")
                .def_prop_ro("phase", &unweighted_pauliString::phase, "Exponent k of the phase i**k, in {0, 1, 2, 3}.")
                .def_prop_ro("phase_factor", &unweighted_pauliString::phase_factor, "The phase i**k as a complex number.")
                .def_prop_ro("x", [](const unweighted_pauliString& s) { return s.x.to_vector(); }, "x bits as 64-bit words (bit q set for X or Z on qubit q).")
                .def_prop_ro("y", [](const unweighted_pauliString& s) { return s.y.to_vector(); }, "y bits as 64-bit words (bit q set for Y or Z on qubit q).")
                .def("__mul__", &unweighted_pauliString::operator*, "Product of two strings, including the phase of the Pauli algebra.")
                .def("__mul__", [](const unweighted_pauliString& s, std::complex<double> scale) { return s.to_pauli_string(scale); },
                        "The weighted Pauli string with coefficient scale * i**phase.")
                .def("__rmul__", [](const unweighted_pauliString& s, std::complex<double> scale) { return s.to_pauli_string(scale); },
                        "The weighted Pauli string with coefficient scale * i**phase.")
                .def("__imul__", &unweighted_pauliString::operator*=, "Multiply in place by another string.")
                .def("__neg__", [](const unweighted_pauliString& s) { return -s; }, "The string multiplied by -1.")
                .def("__eq__", &unweighted_pauliString::operator==, "True if both the operator and the phase are equal (see equals() for the operator only).")
                .def("__hash__", [](const unweighted_pauliString& s) { return static_cast<Py_ssize_t>(s.hash()); }, "Hash over the operator part; ignores the phase.")
                .def("__len__", &unweighted_pauliString::size, "Number of qubits with a non-identity operator (the weight).")
                .def("__repr__", &unweighted_pauliString::to_string, "Human-readable form, e.g. '+i : X(0)Z(2)'.")
                .def("__str__", &unweighted_pauliString::to_string, "Human-readable form, e.g. '+i : X(0)Z(2)'.")
                .def("to_string", &unweighted_pauliString::to_string, "Human-readable form, e.g. '+i : X(0)Z(2)'.")
                .def("equals", &unweighted_pauliString::equals, nb::arg("other"), "True if the operators are equal, ignoring the phases.")
                .def("commutes", &unweighted_pauliString::commutes, nb::arg("other"), "True if the two strings commute, false if they anti-commute.")
                .def("commutator", &unweighted_pauliString::commutator, nb::arg("other"),
                        "The commutator without its factor 2: self * other if the strings anti-commute, otherwise None.")
                .def("dagger", &unweighted_pauliString::dagger, "Hermitian conjugate: the operator stays and the phase is conjugated.")
                .def("naked", &unweighted_pauliString::naked, "Copy of the string with phase 0.")
                .def("to_pauli_string", [](const unweighted_pauliString& s, std::complex<double> coeff) { return s.to_pauli_string(coeff); },
                        nb::arg("coeff") = std::complex<double>(1.0, 0.0), "The Pauli string with coefficient coeff * i**phase.")
                .def("to_pauli_string", [](const unweighted_pauliString& s, const std::string& coeff) { return s.to_pauli_string(SymEngine::Expression(coeff)); },
                        nb::arg("coeff"), "The symbolic Pauli string with coefficient coeff * i**phase; coeff is parsed like 'a*b'.")
                .def("to_pauli_string", [](const unweighted_pauliString& s, const SymEngine::Expression& coeff) { return s.to_pauli_string(coeff); },
                        nb::arg("coeff"), "The symbolic Pauli string with coefficient coeff * i**phase.")
                .def("to_dictionary", &unweighted_pauliString::to_dictionary, "The string as a (phase factor, {qubit: 'X'/'Y'/'Z'}) pair.")
                .def("qubits", &unweighted_pauliString::qubits, "Sorted list of the qubits with a non-identity operator.")
                .def("size", &unweighted_pauliString::size, "Number of qubits with a non-identity operator (the weight).")
                .def("count_y", &unweighted_pauliString::count_y, "Number of Y operators in the string.")
                .def("get_pauli_at_index", &unweighted_pauliString::get_pauli_from_index, nb::arg("index"), "Pauli operator ('I', 'X', 'Y' or 'Z') on the given qubit.")
                .def("key_openfermion", &unweighted_pauliString::key_openfermion, "OpenFermion-style key: list of (Pauli, qubit) pairs of the non-identity sites, sorted by qubit.")
                .def("get_hash", [](const unweighted_pauliString& s) { return hash_limbs_to_int(s.get_hash()); },
                        "Invertible int key of the operator part (phase ignored): sum_q p_q * 4**q with I=0, X=1, Y=2, Z=3.")
                .def_static("from_hash", [](const nb::int_& hash, int phase) { return unweighted_pauliString::from_hash(hash_int_to_limbs(hash), phase); },
                        nb::arg("hash"), nb::arg("phase") = 0, "Build the string whose get_hash() equals hash, with the given phase.");


        // Set of Pauli strings with one phase per element, for Lie algebras
        auto cls_uq = nb::class_<unweighted_QubitHamiltonian>(m, "unweighted_QubitHamiltonian",
                "An ordered set of Pauli operators with one phase i**k per element. Use it only as the basis of a Lie algebra. "
                "It is not a linear combination: every operator occurs once and the first phase seen for it wins, so A + A == A. "
                "An element keeps its index.")
                .def(nb::init<>(), "The empty set.")
                .def(nb::init<const std::vector<unweighted_pauliString>&>(), nb::arg("terms"),
                        "From a list of strings; later duplicates of an operator are ignored.")
                .def(nb::init<const QHC&>(), nb::arg("hamiltonian"), "The operators of a Hamiltonian with phase 0; its coefficients are dropped.")
                .def(nb::init<const QHS&>(), nb::arg("hamiltonian"), "The operators of a symbolic Hamiltonian with phase 0; its coefficients are dropped.")
                .def("__len__", &unweighted_QubitHamiltonian::size, "Number of elements.")
                .def("size", &unweighted_QubitHamiltonian::size, "Number of elements.")
                .def("insert", &unweighted_QubitHamiltonian::insert, nb::arg("term"),
                        "Add the string unless its operator is present; returns whether it was added.")
                .def("__contains__", &unweighted_QubitHamiltonian::contains, nb::arg("term"),
                        "True if the operator of the string is present, whatever the phases are.")
                .def("index", [](const unweighted_QubitHamiltonian& h, const unweighted_pauliString& term) {
                                const std::optional<size_t> i = h.find(term);
                                if (!i) {
                                        throw nb::value_error("unweighted_QubitHamiltonian: operator is not present");
                                }
                                return *i;
                        }, nb::arg("term"), "Index of the operator of the string; raises ValueError if it is not present.")
                .def("__getitem__", [](const unweighted_QubitHamiltonian& h, Py_ssize_t index) {
                                const Py_ssize_t n = static_cast<Py_ssize_t>(h.size());
                                if (index < -n || index >= n) {
                                        throw nb::index_error("unweighted_QubitHamiltonian: index out of range");
                                }
                                return h.terms()[static_cast<size_t>(index < 0 ? index + n : index)];
                        }, nb::arg("index"), "The element with the given index.")
                .def("__iter__", [](const unweighted_QubitHamiltonian& h) { return nb::iter(nb::cast(h.terms())); }, "Iterate over the elements in insertion order.")
                .def("terms", &unweighted_QubitHamiltonian::terms, "All elements in insertion order.")
                .def("paulistrings", &unweighted_QubitHamiltonian::terms, "All elements in insertion order, as QubitHamiltonian.paulistrings().")
                .def("__add__", &unweighted_QubitHamiltonian::operator+, "Union: the elements of self, then the new ones of other.")
                .def("__iadd__", [](unweighted_QubitHamiltonian& h, const unweighted_QubitHamiltonian& other) -> unweighted_QubitHamiltonian& {
                                h += other;
                                return h;
                        }, nb::rv_policy::none, "Add the new elements of other in place; existing indices stay.")
                .def("__mul__", &unweighted_QubitHamiltonian::operator*, "The set of all products a * b in the order of a, then b.")
                .def("commutator", &unweighted_QubitHamiltonian::commutator, nb::arg("other"),
                        "The set of commutators without their factor 2: a * b for every anti-commuting pair.")
                .def("__eq__", &unweighted_QubitHamiltonian::operator==, "Equal sets with equal phases, whatever the order of the elements is.")
                .def("to_qubit_hamiltonian", [](const unweighted_QubitHamiltonian& h, std::complex<double> coeff) { return h.to_qubit_hamiltonian(coeff); },
                        nb::arg("coeff") = std::complex<double>(1.0, 0.0), "The ordinary Hamiltonian with coefficient coeff * i**k on every element.")
                .def("to_qubit_hamiltonian", [](const unweighted_QubitHamiltonian& h, const std::string& coeff) { return h.to_qubit_hamiltonian(SymEngine::Expression(coeff)); },
                        nb::arg("coeff"), "The symbolic Hamiltonian with coefficient coeff * i**k on every element; coeff is parsed like 'a*b'.");

        // PauliCycle and PauliCycleSum, with complex or symbolic coefficients
        auto cls_cy = bind_cycle<std::complex<double>>(m, "PauliCycleComplex",
                "A Pauli string with a complex coefficient together with all of its cyclic rotations on "
                "n_qubits qubits, C(P) = coeff * sum_k R^k P. The sum runs over all n rotations, so a base with a "
                "shorter period counts each distinct rotation several times.");
        auto cls_cys_sym = bind_cycle<SymEngine::Expression>(m, "PauliCycleSymbolic",
                "A Pauli string with a symbolic (SymEngine) coefficient together with all of its cyclic rotations "
                "on n_qubits qubits, C(P) = coeff * sum_k R^k P.");
        auto cls_cys = bind_cycle_sum<std::complex<double>>(m, "PauliCycleSumComplex",
                "A sum of PauliCycles with complex coefficients; each cycle carries its own coefficient. "
                "Summands from the same rotation orbit are always merged.");
        auto cls_css = bind_cycle_sum<SymEngine::Expression>(m, "PauliCycleSumSymbolic",
                "A sum of PauliCycles with symbolic coefficients; each cycle carries its own coefficient. "
                "Summands from the same rotation orbit are always merged.");

        // PauliOrbit (permutation-invariant basis, arXiv:2604.16701 Section V)
        auto cls_or = nb::class_<PauliOrbit>(m, "PauliOrbit",
                "A permutation-invariant Pauli orbit B_{p,q,r} on n qubits: the normalized sum "
                "(i / n_terms) * sum of all Pauli strings with p X-, q Y- and r Z-sites.")
                .def(nb::init<int, int, int, int>(),
                        nb::arg("n_qubits"), nb::arg("p"), nb::arg("q"), nb::arg("r"),
                        "The orbit with p X-, q Y- and r Z-sites on n_qubits qubits. Raises ValueError unless "
                        "0 <= n_qubits <= max_qubits, p, q, r >= 0 and p + q + r <= n_qubits.")
                .def("identities", &PauliOrbit::identities, "Number of identity sites, n - p - q - r.")
                .def("is_valid", &PauliOrbit::is_valid, "Whether p, q, r >= 0 and p + q + r <= n_qubits; always True, since the constructor rejects anything else.")
                .def("to_qubit_hamiltonian", &PauliOrbit::to_qubit_hamiltonian,
                        "The orbit expanded into a QubitHamiltonian, (i / n_terms) * sum of its Pauli strings. "
                        "Raises an error above 2e7 strings.")
                .def("n_terms", &PauliOrbit::n_terms, "Number of distinct Pauli strings in the orbit, n! / (p! q! r! (n-p-q-r)!) (Eq. 52).")
                .def("commutator", &PauliOrbit::commutator, nb::arg("other"),
                        "Commutator [self, other] as a PauliOrbitSum (Theorem 1), without expanding the strings. "
                        "Raises ValueError on an n_qubits mismatch.")
                .def("targeted_commutator", &PauliOrbit::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target in [self, other] (Corollary 5 / Algorithm 2), as a complex number "
                        "with zero imaginary part; O(1) for bounded weights.")
                .def("__mul__", [](const PauliOrbit& a, const PauliOrbit& b) { return a * b; }, nb::arg("other"),
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

        // PauliOrbitSum, with complex or symbolic coefficients
        auto cls_ors = bind_orbit_sum<std::complex<double>>(m, "PauliOrbitSumComplex",
                "A weighted sum of Pauli orbits with complex coefficients, sum_i coeffs[i] * data[i], on a fixed "
                "number of qubits.");
        auto cls_oss = bind_orbit_sum<SymEngine::Expression>(m, "PauliOrbitSumSymbolic",
                "A weighted sum of Pauli orbits with symbolic coefficients, sum_i coeffs[i] * data[i], on a fixed "
                "number of qubits.");

        // Scaling cycles and orbits by an int or float.
        bind_scale_ops<CYC>(cls_cy, "The cycle with its coefficient scaled by an int or float.");
        bind_scale_ops<CYS>(cls_cys_sym, "The cycle with its coefficient scaled by an int or float.");
        bind_scale_ops<CSC>(cls_cys,
                "The sum with every coefficient scaled by an int or float; scaling by 0 gives the empty sum.");
        bind_scale_ops<CSS>(cls_css,
                "The sum with every coefficient scaled by an int or float; scaling by 0 gives the empty sum.");
        bind_scale_ops<PauliOrbit>(cls_or,
                "scale * self as a PauliOrbitSum, for an int or float; scaling by 0 gives the empty sum.");
        bind_scale_ops<OSC>(cls_ors,
                "The sum with every coefficient scaled by an int or float; scaling by 0 gives the empty sum.");
        bind_scale_ops<OSS>(cls_oss,
                "The sum with every coefficient scaled by an int or float; scaling by 0 gives the empty sum.");

        // Pauli decomposition of a matrix (PauliDecomposition.h).
        cls_qhc.def_static("from_matrix",
                [](nb::ndarray<const std::complex<double>, nb::ndim<2>, nb::c_contig, nb::device::cpu> m,
                   double threshold) {
                        if (m.shape(0) != m.shape(1)) {
                                throw std::invalid_argument("from_matrix: the matrix must be square");
                        }
                        const size_t dim = m.shape(0);
                        std::vector<std::complex<double>> data(m.data(), m.data() + dim * dim);
                        return pauli_decomposition(std::move(data), dim, threshold);
                }, nb::arg("matrix"), nb::arg("threshold") = 0.0,
                "Pauli decomposition of a 2**n x 2**n matrix (Georges et al. 2025, O(4**n n)). Qubit 0 is the "
                "most significant index bit, as in to_matrix(ignore_unused_qubits=False). Terms with "
                "abs(coefficient) <= threshold are dropped. The input is not modified.");
        cls_qhc.def_static("from_matrix",
                [](const Matrix2D& m, double threshold) { return pauli_decomposition(m, threshold); },
                nb::arg("matrix"), nb::arg("threshold") = 0.0,
                "Pauli decomposition of a 2**n x 2**n matrix given as a list of rows.");

        // Arithmetic and inner products between all representations (Arithmetic.h).
        // Pairs that the class binds itself above are left out (see own_add etc.).
        bind_arithmetic<PSC>(cls_psc, AllOperators{});
        bind_arithmetic<PSS>(cls_pss, AllOperators{});
        bind_arithmetic<QHC>(cls_qhc, AllOperators{});
        bind_arithmetic<QHS>(cls_qhs, AllOperators{});
        bind_arithmetic<CYC>(cls_cy, AllOperators{});
        bind_arithmetic<CYS>(cls_cys_sym, AllOperators{});
        bind_arithmetic<CSC>(cls_cys, AllOperators{});
        bind_arithmetic<CSS>(cls_css, AllOperators{});
        bind_arithmetic<PauliOrbit>(cls_or, AllOperators{});
        bind_arithmetic<OSC>(cls_ors, AllOperators{});
        bind_arithmetic<OSS>(cls_oss, AllOperators{});

        // Unweighted strings and sets with every other kind of operator (Arithmetic.h).
        using UPS = unweighted_pauliString;
        using UQH = unweighted_QubitHamiltonian;
        bind_unweighted_pair<PSC, UPS>(cls_psc, cls_up);
        bind_unweighted_pair<PSS, UPS>(cls_pss, cls_up);
        bind_unweighted_pair<QHC, UPS>(cls_qhc, cls_up);
        bind_unweighted_pair<QHS, UPS>(cls_qhs, cls_up);
        bind_unweighted_pair<CYC, UPS>(cls_cy, cls_up);
        bind_unweighted_pair<CYS, UPS>(cls_cys_sym, cls_up);
        bind_unweighted_pair<CSC, UPS>(cls_cys, cls_up);
        bind_unweighted_pair<CSS, UPS>(cls_css, cls_up);
        bind_unweighted_pair<PauliOrbit, UPS>(cls_or, cls_up);
        bind_unweighted_pair<OSC, UPS>(cls_ors, cls_up);
        bind_unweighted_pair<OSS, UPS>(cls_oss, cls_up);
        bind_unweighted_pair<PSC, UQH>(cls_psc, cls_uq);
        bind_unweighted_pair<PSS, UQH>(cls_pss, cls_uq);
        bind_unweighted_pair<QHC, UQH>(cls_qhc, cls_uq);
        bind_unweighted_pair<QHS, UQH>(cls_qhs, cls_uq);
        bind_unweighted_pair<CYC, UQH>(cls_cy, cls_uq);
        bind_unweighted_pair<CYS, UQH>(cls_cys_sym, cls_uq);
        bind_unweighted_pair<CSC, UQH>(cls_cys, cls_uq);
        bind_unweighted_pair<CSS, UQH>(cls_css, cls_uq);
        bind_unweighted_pair<PauliOrbit, UQH>(cls_or, cls_uq);
        bind_unweighted_pair<OSC, UQH>(cls_ors, cls_uq);
        bind_unweighted_pair<OSS, UQH>(cls_oss, cls_uq);
        // A string and a set: the string counts as a set with one element. There is no
        // difference of sets, so only + and * are defined.
        cls_uq.def("__add__", [](const UQH& h, const UPS& s) { return h + s; },
                   "Add the string unless its operator is present, as a new set.")
              .def("__mul__", [](const UQH& h, const UPS& s) { return h * s; }, "The set of all products element * string.")
              .def("commutator", [](const UQH& h, const UPS& s) { return commutator(h, s); }, nb::arg("other"),
                   "The set of commutators without their factor 2 with the string.");
        cls_up.def("__add__", [](const UPS& s, const UQH& h) { return s + h; },
                   "The string followed by the new elements of the set.")
              .def("__mul__", [](const UPS& s, const UQH& h) { return s * h; }, "The set of all products string * element.")
              .def("commutator", [](const UPS& s, const UQH& h) { return commutator(s, h); }, nb::arg("other"),
                   "The set of commutators without their factor 2 with every element.");
        cls_uq.def("__mul__", [](const UQH& h, std::complex<double> s) { return h * s; },
                   "The weighted Hamiltonian with coefficient scale * i**k on every element.")
              .def("__rmul__", [](const UQH& h, std::complex<double> s) { return h * s; },
                   "The weighted Hamiltonian with coefficient scale * i**k on every element.");
        bind_contains_unweighted<QHC>(cls_qhc);
        bind_contains_unweighted<QHS>(cls_qhs);
}
