#include <bit>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/variant.h>
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
#include "pauliengine/SparsePauliStringSum.h"
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
using PSS = PauliString<Expression>;
using QHC = QubitHamiltonian<std::complex<double>>;
using QHS = QubitHamiltonian<Expression>;
using SPC = SparsePauliString<std::complex<double>>;
using SPS = SparsePauliString<Expression>;
using SQC = SparsePauliStringSum<std::complex<double>>;
using SQS = SparsePauliStringSum<Expression>;

// Whether a coefficient returned by contains() is non-zero (for `in`).
static bool coeff_nonzero(const std::complex<double>& c) { return c != std::complex<double>(0.0, 0.0); }
static bool coeff_nonzero(const Expression& c) { return !c.is_zero(); }

// Values for subs() from Python: numbers, Expressions or polynomial texts such as "b + 2*c".
using SubsValue = std::variant<std::complex<double>, std::string, Expression>;
using SubsMap = std::unordered_map<std::string, SubsValue>;

static const char* const kSubsDoc =
        "Substitute values for symbols, given as {name: value}. A value is a number, an Expression or a "
        "polynomial as text such as 'b + 2*c'; all symbols are replaced at once, so {'a': 'b', 'b': 'a'} swaps "
        "them. Symbols without a value stay.";

// x.substitute(...) with the values of `values`: numbers only, or polynomials as
// soon as one value is not a number.
template <typename T>
auto subs_any(const T& x, const SubsMap& values) {
        bool numbers = true;
        for (const auto& [name, v] : values) numbers = numbers && std::holds_alternative<std::complex<double>>(v);
        if (numbers) {
                pauliengine_detail::SymbolValues plain;
                for (const auto& [name, v] : values) plain.emplace(name, std::get<std::complex<double>>(v));
                return x.substitute(plain);
        }
        pauliengine_detail::SymbolPolys polys;
        for (const auto& [name, v] : values) {
                if (const auto* c = std::get_if<std::complex<double>>(&v)) polys.emplace(name, Expression(*c));
                else if (const auto* t = std::get_if<std::string>(&v)) polys.emplace(name, Expression(*t));
                else polys.emplace(name, std::get<Expression>(v));
        }
        return x.substitute(polys);
}

// subs_any for an Expression itself.
struct ExpressionSubs {
        const Expression& e;
        template <typename V>
        Expression substitute(const V& v) const { return e.subs(v); }
};

// Register a * s and s * a of class A for an int, float or complex s, and -a.
template <typename A, typename Cls>
void bind_scale_ops(Cls& cls, const char* doc) {
        cls.def("__mul__", [](const A& a, int s) { return a * s; }, nb::arg("scale"), doc)
                .def("__mul__", [](const A& a, double s) { return a * s; }, nb::arg("scale"), doc)
                .def("__mul__", [](const A& a, std::complex<double> s) { return a * s; }, nb::arg("scale"), doc)
                .def("__rmul__", [](const A& a, int s) { return a * s; }, nb::arg("scale"), doc)
                .def("__rmul__", [](const A& a, double s) { return a * s; }, nb::arg("scale"), doc)
                .def("__rmul__", [](const A& a, std::complex<double> s) { return a * s; }, nb::arg("scale"), doc)
                .def("__neg__", [](const A& a) { return a * -1; }, "-1 * self.");
}

// a / s for an int, float or complex s, as a * (1 / s); weighted(a) turns an
// unweighted operator into its weighted form first.
template <typename A, typename Cls, typename Weighted>
void bind_division(Cls& cls, Weighted weighted) {
        const auto zero_division = [] {
                PyErr_SetString(PyExc_ZeroDivisionError, "division by zero");
                throw nb::python_error();
        };
        cls.def("__truediv__", [weighted, zero_division](const A& a, double s) {
                        if (s == 0.0) zero_division();
                        return weighted(a) * (1.0 / s);
                }, nb::arg("scale"), "self * (1 / scale) for an int or float scale.")
           .def("__truediv__", [weighted, zero_division](const A& a, std::complex<double> s) {
                        if (s == std::complex<double>(0.0)) zero_division();
                        return weighted(a) * (std::complex<double>(1.0) / s);
                }, nb::arg("scale"), "self * (1 / scale) for a complex scale.");
}

// a * e and e * a for a symbolic e (scaled_by in Arithmetic.h).
template <typename A, typename Cls>
void bind_symbolic_scale(Cls& cls) {
        cls.def("__mul__", [](const A& a, const Expression& e) { return scaled_by(a, e); }, nb::arg("scale"),
                "Every coefficient multiplied by a symbolic expression; the result is symbolic.")
           .def("__rmul__", [](const A& a, const Expression& e) { return scaled_by(a, e); }, nb::arg("scale"),
                "Every coefficient multiplied by a symbolic expression; the result is symbolic.");
}

// a == b and a != b for an A and a B, decided by eq(a, b).
template <typename A, typename B, typename Cls, typename Eq>
void bind_eq(Cls& cls, Eq eq, const char* doc) {
        cls.def("__eq__", [eq](const A& a, const B& b) { return static_cast<bool>(eq(a, b)); }, doc)
           .def("__ne__", [eq](const A& a, const B& b) { return !static_cast<bool>(eq(a, b)); }, "Negation of ==.");
}

// == with any object not compared above gives NotImplemented, so Python answers
// False (and != True) instead of raising a TypeError. Registered after all other
// __eq__ overloads, which are tried first.
template <typename A, typename Cls>
void bind_eq_fallback(Cls& cls) {
        cls.def("__eq__", [](const A&, nb::handle) { return nb::borrow<nb::object>(Py_NotImplemented); },
                nb::arg("other").none(), "NotImplemented for objects of other types.")
           .def("__ne__", [](const A&, nb::handle) { return nb::borrow<nb::object>(Py_NotImplemented); },
                nb::arg("other").none(), "NotImplemented for objects of other types.");
}

// A term of a list given to a constructor: a Pauli string with bit vectors or a
// sparse one, in any mix.
template <typename C>
using AnyString = std::variant<PauliString<C>, SparsePauliString<C>>;

template <typename C>
std::vector<PauliString<C>> dense_strings(const std::vector<AnyString<C>>& terms) {
        std::vector<PauliString<C>> out;
        out.reserve(terms.size());
        for (const auto& t : terms) {
                if (const auto* d = std::get_if<PauliString<C>>(&t)) out.push_back(*d);
                else out.push_back(std::get<SparsePauliString<C>>(t).to_dense());
        }
        return out;
}

template <typename C>
std::vector<SparsePauliString<C>> sparse_strings(const std::vector<AnyString<C>>& terms) {
        std::vector<SparsePauliString<C>> out;
        out.reserve(terms.size());
        for (const auto& t : terms) {
                if (const auto* sp = std::get_if<SparsePauliString<C>>(&t)) out.push_back(*sp);
                else out.emplace_back(std::get<PauliString<C>>(t));
        }
        return out;
}

// The result of a sum or difference, as returned to Python: a sparse operator that no
// longer fits the sparse representation (its strings grew too heavy, see
// pauliengine_detail::sparse_stays) is handed back with bit vectors instead.
template <typename T>
T finish(T value) {
        return value;
}
template <typename C>
nb::object finish(SparsePauliStringSum<C> h) {
        if (h.fits_sparse()) return nb::cast(std::move(h));
        return nb::cast(h.to_dense());
}
template <typename C>
nb::object finish(SparsePauliString<C> s) {
        if (s.fits_sparse()) return nb::cast(std::move(s));
        return nb::cast(s.to_dense());
}

// Products and commutators (also square and power) are always handed back with
// bit vectors: they make strings heavier, so their results rarely stay sparse.
// They are still computed in the sparse representation when both operands are sparse.
template <typename T>
auto dense_result(T value) {
        if constexpr (pauliengine_detail::is_sparse<T>) {
                return value.to_dense();
        } else {
                return value;
        }
}

// Arithmetic of an unweighted string or set U with a weighted operator W
// (Arithmetic.h): U counts with coefficient 1.0 * i**k on every element.
template <typename W, typename U, typename ClsW, typename ClsU>
void bind_unweighted_pair(ClsW& cls_w, ClsU& cls_u) {
        cls_w.def("__add__", [](const W& a, const U& b) { return finish(a + b); },
                  "Sum with an unweighted operator, which counts with coefficient 1.0 * i**phase.")
             .def("__mul__", [](const W& a, const U& b) { return dense_result(a * b); },
                  "Product with an unweighted operator, which counts with coefficient 1.0 * i**phase.")
             .def("commutator", [](const W& a, const U& b) { return dense_result(commutator(a, b)); }, nb::arg("other"),
                  "Commutator with an unweighted operator, which counts with coefficient 1.0 * i**phase.");
        cls_u.def("__add__", [](const U& a, const W& b) { return finish(a + b); },
                  "Sum with another operator; next to a weighted one this counts with coefficient 1.0 * i**phase.")
             .def("__mul__", [](const U& a, const W& b) { return dense_result(a * b); },
                  "Product with another operator; next to a weighted one this counts with coefficient 1.0 * i**phase.")
             .def("commutator", [](const U& a, const W& b) { return dense_result(commutator(a, b)); }, nb::arg("other"),
                  "Commutator with another operator; next to a weighted one this counts with coefficient 1.0 * i**phase.");
        cls_w.def("__sub__", [](const W& a, const U& b) { return finish(a - b); },
                  "Difference with an unweighted operator, which counts with coefficient 1.0 * i**phase.")
             .def("inner_product", [](const W& a, const U& b) { return inner_product(a, b); }, nb::arg("other"),
                  "Normalized Hilbert-Schmidt inner product with an unweighted operator, which counts with "
                  "coefficient 1.0 * i**phase.");
        cls_u.def("__sub__", [](const U& a, const W& b) { return finish(a - b); },
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
using CYS = PauliCycle<Expression>;
using CSC = PauliCycleSum<std::complex<double>>;
using CSS = PauliCycleSum<Expression>;
using OSC = PauliOrbitSum<std::complex<double>>;
using OSS = PauliOrbitSum<Expression>;

template <typename... Ts>
struct TypeList {};

/// Every operator representation, for the mixed arithmetic below.
using AllOperators = TypeList<PSC, PSS, QHC, QHS, CYC, CYS, CSC, CSS, PauliOrbit, OSC, OSS, SPC, SPS, SQC, SQS>;

// Pairs whose operation is bound by the class itself, so the mixed binding is left out.
template <typename A, typename B>
inline constexpr bool own_add = (pauliengine_detail::is_pauli_string<A> && pauliengine_detail::is_pauli_string<B>)
        || (pauliengine_detail::is_hamiltonian<A> && std::is_same_v<A, B>);
template <typename A, typename B>
inline constexpr bool own_mul = (pauliengine_detail::is_pauli_string<A> && pauliengine_detail::is_pauli_string<B>)
        || std::is_same_v<A, B>;
template <typename A, typename B>
inline constexpr bool own_commutator = std::is_same_v<A, B>;

// Whether a Pauli string or Hamiltonian with bit vectors would itself be sparse
// (pauliengine_detail::sparse_stays): only then does a sparse operand next to it
// keep the operation sparse.
template <typename C>
bool dense_fits_sparse(const PauliString<C>& ps) {
        size_t weight = 0;
        int top = -1;
        for (size_t w = 0; w < ps.x.size(); ++w) {
                const uint64_t bits = ps.x[w] | ps.y[w];
                weight += static_cast<size_t>(std::popcount(bits));
                if (bits) top = static_cast<int>(64 * w + 63 - std::countl_zero(bits));
        }
        return pauliengine_detail::sparse_stays(top + 1.0, static_cast<double>(weight));
}
template <typename C>
bool dense_fits_sparse(const QubitHamiltonian<C>& h) {
        size_t max_weight = 0;
        int top = -1;
        for (const auto& ps : h.terms()) {
                size_t weight = 0;
                for (size_t w = 0; w < ps.x.size(); ++w) {
                        const uint64_t bits = ps.x[w] | ps.y[w];
                        weight += static_cast<size_t>(std::popcount(bits));
                        if (bits) top = std::max(top, static_cast<int>(64 * w + 63 - std::countl_zero(bits)));
                }
                max_weight = std::max(max_weight, weight);
        }
        return pauliengine_detail::sparse_stays(top + 1.0, static_cast<double>(max_weight));
}

template <typename T>
inline constexpr bool is_dense_string_family = pauliengine_detail::is_pauli_string<T> || pauliengine_detail::is_hamiltonian<T>;

// a + b or a - b for the bindings. A sparse operand next to a Pauli string or Hamiltonian
// with bit vectors that does not fit the sparse representation (for example a
// product that already went back to bit vectors) is turned into bit vectors
// instead, so that chains of operations do not switch back and forth.
template <typename A, typename B, typename Op>
auto mixed(const A& a, const B& b, Op op) {
        using namespace pauliengine_detail;
        if constexpr (is_sparse<A> && is_dense_string_family<B>) {
                if (!dense_fits_sparse(b)) return nb::cast(op(a.to_dense(), b));
                return nb::object(finish(op(a, b)));
        } else if constexpr (is_dense_string_family<A> && is_sparse<B>) {
                if (!dense_fits_sparse(a)) return nb::cast(op(a, b.to_dense()));
                return nb::object(finish(op(a, b)));
        } else {
                return finish(op(a, b));
        }
}

// a * b or [a, b] for the bindings, always with bit vectors (see dense_result).
// A sparse operand next to a Pauli string or Hamiltonian with bit vectors is
// converted first, so the bit-vector algorithms do the work; next to cycles, orbits
// and other sparse operators the sparse algorithms run and only the result is converted.
template <typename A, typename B, typename Op>
auto product_like(const A& a, const B& b, Op op) {
        using namespace pauliengine_detail;
        if constexpr (is_sparse<A> && is_dense_string_family<B>) {
                return dense_result(op(a.to_dense(), b));
        } else if constexpr (is_dense_string_family<A> && is_sparse<B>) {
                return dense_result(op(a, b.to_dense()));
        } else {
                return dense_result(op(a, b));
        }
}

// Mixed arithmetic (see Arithmetic.h): register +, -, *, commutator and
// inner_product of class A with every B, unless A binds the operation itself.
template <typename A, typename Cls, typename... Bs>
void bind_arithmetic(Cls& cls, TypeList<Bs...>) {
        ([&] {
                using B = Bs;
                if constexpr (!own_add<A, B>) {
                        cls.def("__add__", [](const A& a, const B& b) { return mixed(a, b, [](const auto& x, const auto& y) { return x + y; }); },
                                "Sum with another operator; a PauliStringSum unless both are in the same family.");
                        cls.def("__sub__", [](const A& a, const B& b) { return mixed(a, b, [](const auto& x, const auto& y) { return x - y; }); },
                                "Difference with another operator; a PauliStringSum unless both are in the same family.");
                }
                if constexpr (!own_mul<A, B>) {
                        cls.def("__mul__", [](const A& a, const B& b) { return product_like(a, b, [](const auto& x, const auto& y) { return x * y; }); },
                                "Product with another operator; a PauliStringSum unless both are in the same family.");
                }
                if constexpr (!own_commutator<A, B>) {
                        cls.def("commutator", [](const A& a, const B& b) { return product_like(a, b, [](const auto& x, const auto& y) { return commutator(x, y); }); }, nb::arg("other"),
                                "Commutator [self, other] with another operator; a PauliStringSum unless both are in "
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
                .def("__init__", [](Cy* self, int n_qubits, const SparsePauliString<C>& base) {
                                new (self) Cy(n_qubits, base.to_dense());
                        }, nb::arg("n_qubits"), nb::arg("base"),
                        "From a qubit count and a sparse base Pauli string (its coefficient becomes the cycle's).")
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
                .def("to_pauli_string_sum", &Cy::to_qubit_hamiltonian,
                        "The cycle expanded into a PauliStringSum (the sum of all n rotations).")
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
                .def_ro("base", &Cy::Base, "The stored base Pauli string (as constructed, or unified after unify(); the summands of products and "
                        "commutators of cycle sums with sparse bases are unified already); carries the coefficient.")
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
                .def("to_pauli_string_sum", &Sum::to_qubit_hamiltonian,
                        "The sum expanded into a PauliStringSum.")
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
        if constexpr (std::is_same_v<C, Expression>) {
                cls.def("add", [](Sum& s, const PauliOrbit& o, const std::string& coeff) { s.add(o, Expression(coeff)); },
                        nb::arg("orbit"), nb::arg("coeff"),
                        "Add coeff * orbit in place; the coefficient is parsed like 'a*b'.")
                   .def("add", [](Sum& s, const PauliOrbit& o, std::complex<double> coeff) { s.add(o, Expression(coeff)); },
                        nb::arg("orbit"), nb::arg("coeff"),
                        "Add coeff * orbit in place, with a numeric coefficient.");
        }
        cls.def("to_pauli_string_sum", &Sum::to_qubit_hamiltonian,
                        "The sum expanded into a PauliStringSum.")
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


// SparsePauliString with coefficients of type C.
template <typename C>
auto bind_sparse_string(nb::module_& m, const char* name, const char* doc) {
        using S = SparsePauliString<C>;
        using Dense = PauliString<C>;
        auto cls = nb::class_<S>(m, name, doc)
                .def(nb::init<C, const std::unordered_map<int, std::string>&>(), nb::arg("coeff"), nb::arg("data"),
                        "From a coefficient and a map {qubit: 'X'/'Y'/'Z'/'I'}.")
                .def("__init__", [](S* self, C coeff, const std::string& text) { new (self) S(Dense(coeff, text)); },
                        nb::arg("coeff"), nb::arg("pauli_string"),
                        "From a coefficient and a text: dense such as 'XIZ' or sparse such as 'X0 Z2'.")
                .def(nb::init<const Dense&>(), nb::arg("pauli_string"), "The same string as a Pauli string with bit vectors.")
                .def("to_dense", &S::to_dense, "The string as a Pauli string with bit vectors.")
                .def("to_string", &S::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__str__", &S::to_string, "Human-readable form.")
                .def("__repr__", &S::to_string, "Human-readable form.")
                .def("size", &S::size, "Number of qubits with a non-identity operator (the weight).")
                .def("__len__", &S::size, "Number of qubits with a non-identity operator (the weight).")
                .def("qubits", &S::qubits, "Sorted list of the qubits with a non-identity operator.")
                .def("get_pauli_at_index", &S::get_pauli_at_index, "Pauli operator ('I', 'X', 'Y' or 'Z') on the given qubit.")
                .def("count_y", &S::count_y, "Number of Y operators.")
                .def("is_all_z", &S::is_all_z, "True if every non-identity operator is Z.")
                .def("to_dictionary", &S::to_dictionary, "The string as a (coefficient, {qubit: 'X'/'Y'/'Z'}) pair.")
                .def("key_openfermion", &S::key_openfermion, "OpenFermion-style key: list of (Pauli, qubit) pairs, sorted by qubit.")
                .def("equals", &S::equals, "True if the operators are equal, ignoring the coefficients.")
                .def("naked", &S::naked, "Copy with coefficient 1.")
                .def("copy", [](const S& s) { return s; }, "An independent copy.")
                .def("get_coeff", [](const S& s) { return s.coeff; }, "The coefficient.")
                .def("set_coeff", [](S& s, C c) { s.coeff = c; }, "Replace the coefficient (in place).")
                .def_ro("coeff", &S::coeff, "The coefficient (read-only).")
                .def_prop_ro("sites", [](const S& s) {
                        std::vector<std::pair<int, std::string>> out;
                        out.reserve(s.sites.size());
                        for (auto site : s.sites) {
                                out.emplace_back(static_cast<int>(pauliengine_detail::site_qubit(site)),
                                                 std::string(1, pauliengine_detail::pauli_letter(pauliengine_detail::site_pauli(site))));
                        }
                        return out;
                }, "The non-identity sites as (qubit, Pauli) pairs by ascending qubit.")
                .def_prop_ro("x", [](const S& s) { return s.to_dense().x.to_vector(); },
                        "x bits as 64-bit words (as PauliString.x; formed from the sites).")
                .def_prop_ro("y", [](const S& s) { return s.to_dense().y.to_vector(); },
                        "y bits as 64-bit words (as PauliString.y; formed from the sites).")
                .def_prop_ro("is_zero", [](const S& s) { return pauliengine_detail::coeff_is_zero(s.coeff); },
                        "True if the coefficient is zero.")
                .def("__hash__", [](const S& s) { return static_cast<Py_ssize_t>(s.to_dense().hash()); },
                        "Hash over the operator part, the same as for the Pauli string with bit vectors.")
                .def("get_hash", [](const S& s) { return hash_limbs_to_int(s.to_dense().get_hash()); },
                        "Invertible int key of the operator part (as PauliString.get_hash()).")
                .def("__mul__", [](const S& a, const S& b) { return dense_result(a * b); }, "Product with the phase of the Pauli algebra, as a PauliString with bit vectors.")
                .def("commutator", [](const S& a, const S& b) { return dense_result(a.commutator(b)); }, nb::arg("other"),
                        "Commutator [self, other]: 2 * self * other if they anti-commute, else the zero string; "
                        "as a PauliString with bit vectors.")
                .def("fits_sparse", &S::fits_sparse,
                        "Whether the string still fits the sparse representation (products and commutators are "
                        "always returned with bit vectors).")
                .def("anticommutes", &S::anticommutes, nb::arg("other"), "True if the strings anti-commute.")
                .def("targeted_commutator", &S::targeted_commutator, nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's operator in [self, other].")
                .def("map_qubits", &S::map_qubits, "Copy with qubit indices remapped by {old: new}.")
                .def("trace_out_qubits", [](const S& s, const std::vector<int>& qubits, const std::vector<int>& state) {
                                return S(s.to_dense().trace_out_qubits(qubits, state));
                        }, nb::arg("qubits"), nb::arg("state"), "Trace out the given qubits in computational basis states.")
                .def("trace_out_qubits", [](const S& s, const std::vector<int>& qubits,
                                            const std::vector<std::pair<std::complex<double>, std::complex<double>>>& states) {
                                return S(s.to_dense().trace_out_qubits(qubits, states));
                        }, nb::arg("qubits"), nb::arg("states"), "Trace out the given qubits in states a|0> + b|1>.");
        if constexpr (std::is_same_v<C, Expression>) {
                cls.def("diff", &S::diff, nb::arg("symbol"), "Derivative of the coefficient with respect to the named symbol.")
                   .def("subs", [](const S& x, const SubsMap& v) { return subs_any(x, v); }, nb::arg("values"), kSubsDoc);
        }
        return cls;
}

// SparsePauliStringSum with coefficients of type C.
template <typename C>
auto bind_sparse_hamiltonian(nb::module_& m, const char* name, const char* doc) {
        using H = SparsePauliStringSum<C>;
        using S = SparsePauliString<C>;
        using Dense = QubitHamiltonian<C>;
        auto cls = nb::class_<H>(m, name, doc)
                .def(nb::init<const std::vector<S>&>(), nb::arg("data"), "From a list of sparse Pauli strings (equal operators are merged).")
                .def(nb::init<const std::vector<PauliString<C>>&>(), nb::arg("data"), "From a list of Pauli strings with bit vectors.")
                .def("__init__", [](H* self, const std::vector<AnyString<C>>& data) { new (self) H(sparse_strings(data)); },
                        nb::arg("data"), "From a list that mixes sparse Pauli strings and Pauli strings with bit vectors.")
                .def(nb::init<const Hamiltonian_structure<C>&>(), nb::arg("data"), "From a list of (coefficient, {qubit: 'X'/'Y'/'Z'}) pairs.")
                .def(nb::init<const Dense&>(), nb::arg("hamiltonian"), "The same operator as a Hamiltonian with bit vectors.")
                .def("to_dense", &H::to_dense, "The operator as a PauliStringSum with bit vectors.")
                .def("to_string", &H::to_string, "Human-readable form, one term per line.")
                .def("__str__", &H::to_string, "Human-readable form, one term per line.")
                .def("__repr__", &H::to_string, "Human-readable form, one term per line.")
                .def("__len__", &H::size, "Number of terms.")
                .def("size", &H::size, "Number of terms.")
                .def("to_dictionary", &H::to_dictionary, "List of (coefficient, {qubit: Pauli}) pairs.")
                .def("paulistrings", &H::paulistrings, "The terms as sparse Pauli strings.")
                .def("to_list", &H::paulistrings, "The terms as sparse Pauli strings (same as paulistrings()).")
                .def("contains", [](const H& h, const S& s) { return h.contains(s); }, nb::arg("pauli_string"),
                        "Coefficient of the given Pauli string's operator (0 if absent).")
                .def("contains", [](const H& h, const PauliString<C>& s) { return h.contains(S(s)); }, nb::arg("pauli_string"),
                        "Coefficient of the given Pauli string's operator (0 if absent).")
                .def("contains", [](const H& h, const std::unordered_map<int, std::string>& d) {
                                return h.contains(S(pauliengine_detail::coeff_one<C>(), d));
                        }, nb::arg("pauli_string"), "Coefficient of the Pauli string {qubit: 'X'} (0 if absent).")
                .def("__contains__", [](const H& h, const S& s) { return coeff_nonzero(h.contains(s)); }, nb::arg("pauli_string"),
                        "True if the operator has a non-zero term.")
                .def("__contains__", [](const H& h, const PauliString<C>& s) { return coeff_nonzero(h.contains(S(s))); },
                        nb::arg("pauli_string"), "True if the operator has a non-zero term.")
                .def("__add__", [](const H& a, const H& b) { return finish(a + b); }, "Sum.")
                .def("__sub__", [](const H& a, const H& b) { return finish(a - b); }, "Difference.")
                .def("__mul__", [](const H& a, const H& b) { return dense_result(a * b); },
                        "Product, as a PauliStringSum with bit vectors; h * h uses the square formula (only "
                        "commuting pairs, each once).")
                .def("square", [](const H& h) { return dense_result(h.square()); }, "self * self, visiting each commuting pair once.")
                .def("__pow__", [](const H& h, int k) { return dense_result(h.power(k)); }, "Raise to a non-negative integer power.")
                .def("power", [](const H& h, int k) { return dense_result(h.power(k)); }, "Raise to a non-negative integer power.")
                .def("fits_sparse", &H::fits_sparse,
                        "Whether the operator still fits the sparse representation; sums that do not are returned "
                        "with bit vectors (products and commutators always are).")
                .def("commutator", [](const H& a, const H& b) { return dense_result(a.commutator(b)); }, nb::arg("other"),
                        "Commutator [self, other]; only pairs sharing a qubit are visited when they are few.")
                .def("targeted_commutator", [](const H& a, const H& b, const S& t) { return a.targeted_commutator(b, t); },
                        nb::arg("other"), nb::arg("target"), "Coefficient of target's operator in [self, other], via a hash join.")
                .def("targeted_commutator", [](const H& a, const S& b, const S& t) {
                                return a.targeted_commutator(H(std::vector<S>{b}), t);
                        }, nb::arg("other"), nb::arg("target"), "Coefficient of target's operator in [self, other] for a single string other.")
                .def("targeted_commutator", [](const H& a, const H& b, const PauliString<C>& t) { return a.targeted_commutator(b, S(t)); },
                        nb::arg("other"), nb::arg("target"), "Coefficient of target's operator in [self, other], via a hash join.")
                .def("qubits", &H::qubits, "Sorted list of the qubits acted on non-trivially.")
                .def("n_qubits", &H::n_qubits, "Number of qubits acted on non-trivially.")
                .def("is_all_z", &H::is_all_z, "True if all terms are products of Z operators only.")
                .def("count_measurements", &H::count_measurements, "1 if is_all_z(), else the number of terms.")
                .def("simplify", &H::simplify, nb::arg("threshold") = 0.0, "Copy without the terms whose abs(coefficient) <= threshold.")
                .def("compact", [](const H& h) { return h; }, "The sum is always compacted; returns a copy.")
                .def("set_all_coeff", [](const H& h, C v) { return h.set_all_coeff(v); }, nb::arg("value"),
                        "Copy with every coefficient replaced by value.")
                .def("dagger", &H::dagger, "Hermitian conjugate: complex-conjugates each coefficient.")
                .def("conjugate", &H::conjugate, "Complex conjugate: conjugates each coefficient, sign flipped for an odd number of Y.")
                .def("transpose", &H::transpose, "Transpose: sign flipped for an odd number of Y.")
                .def("map_qubits", &H::map_qubits, "Copy with qubit indices remapped by {old: new}.")
                .def("is_hermitian", [](const H& h) { return h.to_dense().is_hermitian(); },
                        "True if every coefficient is real (abs(imag) <= 1e-12).")
                .def("is_antihermitian", [](const H& h) { return h.to_dense().is_antihermitian(); },
                        "True if every coefficient is purely imaginary (abs(real) <= 1e-12).")
                .def("split", [](const H& h) {
                                auto parts = h.to_dense().split();
                                return std::make_pair(H(parts.first), H(parts.second));
                        }, "Split into (hermitian, anti-hermitian) parts; needs numeric coefficients.")
                .def("to_matrix", [](const H& h, bool ignore) { return h.to_dense().to_matrix(ignore); },
                        nb::arg("ignore_unused_qubits") = true, "Dense 2**N x 2**N matrix (see PauliStringSum.to_matrix).")
                .def("trace_out_qubits", [](const H& h, const std::vector<int>& qubits, const std::vector<int>& state) {
                                return H(h.to_dense().trace_out_qubits(qubits, state));
                        }, nb::arg("qubits"), nb::arg("state"), "Trace out the given qubits in computational basis states.")
                .def("trace_out_qubits", [](const H& h, const std::vector<int>& qubits,
                                            const std::vector<std::pair<std::complex<double>, std::complex<double>>>& states) {
                                return H(h.to_dense().trace_out_qubits(qubits, states));
                        }, nb::arg("qubits"), nb::arg("states"), "Trace out the given qubits in states a|0> + b|1>.")
                .def_static("zero", &H::zero, "The zero operator.")
                .def_static("unit", &H::unit, "The identity (one term, coefficient 1).");
        if constexpr (std::is_same_v<C, Expression>) {
                cls.def("diff", &H::diff, nb::arg("symbol"), "Derivative of every coefficient with respect to the named symbol.")
                   .def("subs", [](const H& x, const SubsMap& v) { return subs_any(x, v); }, nb::arg("values"), kSubsDoc);
        }
        return cls;
}

NB_MODULE(_core, m) {

        m.doc() = "PauliEngine";

        // When operators are stored sparse (pauliengine_detail::SparseThreshold).
        m.def("_sparse_threshold", []() {
                const auto& t = pauliengine_detail::sparse_threshold();
                return std::make_pair(t.min_qubits, t.max_density);
        }, "(min_qubits, max_density) of the sparse threshold; see pauliengine.set_sparse_threshold.");
        m.def("_set_sparse_threshold", [](double min_qubits, double max_density) {
                auto& t = pauliengine_detail::sparse_threshold();
                t.min_qubits = min_qubits;
                t.max_density = max_density;
        }, nb::arg("min_qubits"), nb::arg("max_density"), "Set the sparse threshold; see pauliengine.set_sparse_threshold.");
        m.def("_sparse_enters", &pauliengine_detail::sparse_enters, nb::arg("span"), nb::arg("max_weight"),
              "Whether an operator with this span and heaviest string is built sparse.");
        m.def("_fits_sparse", [](const QHC& h) { return dense_fits_sparse(h); }, nb::arg("hamiltonian"),
              "Whether a sum with bit vectors would be sparse by the rule.");
        m.def("_fits_sparse", [](const QHS& h) { return dense_fits_sparse(h); }, nb::arg("hamiltonian"),
              "Whether a sum with bit vectors would be sparse by the rule.");
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
                // Type heterogeneous multiplication with Expression.
                .def("__mul__", [](const PauliString<std::complex<double>>& ps1, const PauliString<Expression>& ps2){return ps1 * ps2;}, "Product with a symbolic Pauli string; the result is symbolic.")
                .def("__rmul__", [](const PauliString<std::complex<double>>& self, const PauliString<Expression>& other){return other * self;}, "Product other * self with a symbolic Pauli string; the result is symbolic.")
                .def("__repr__",  &PauliString<std::complex<double>>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__str__",  &PauliString<std::complex<double>>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__hash__", [](const PauliString<std::complex<double>>& ps) {
                        return static_cast<Py_ssize_t>(ps.hash());
                }, "Hash over the operator part (x, y); ignores the coefficient.")
                .def("__neg__", [](const PauliString<std::complex<double>>& ps){ return ps * -1; }, "The Pauli string with negated coefficient.")
                .def("__add__", [](const PauliString<std::complex<double>>& ps1, const PauliString<std::complex<double>>& ps2){
                        return QubitHamiltonian<std::complex<double>>({ps1, ps2});
                }, "Add two PauliStringComplex to get PauliStringSumComplex")
                .def("__add__", [](const PauliString<std::complex<double>>& ps1, const PauliString<Expression>& ps2){
                        Expression sym_coeff = ps1.coeff;
                        PauliString<Expression> ps1_sym(ps1.x, ps1.y, sym_coeff);
                        return QubitHamiltonian<Expression>({ps1_sym, ps2});
                }, "Add PauliStringComplex and PauliStringSymbolic to get PauliStringSumSymbolic")
                .def("__sub__", [](const PauliString<std::complex<double>>& ps1, const PauliString<std::complex<double>>& ps2){
                        return QubitHamiltonian<std::complex<double>>({ps1, ps2 * -1});
                }, "Subtract two PauliStringComplex to get PauliStringSumComplex")
                .def("__sub__", [](const PauliString<std::complex<double>>& ps1, const PauliString<Expression>& ps2){
                        Expression sym_coeff = ps1.coeff;
                        PauliString<Expression> ps1_sym(ps1.x, ps1.y, sym_coeff);
                        return QubitHamiltonian<Expression>({ps1_sym, ps2 * -1});
                }, "Subtract PauliStringComplex and PauliStringSymbolic to get PauliStringSumSymbolic")
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

        // PauliString for Expression coefficients
        auto cls_pss = nb::class_<PauliString<Expression>>(m, "PauliStringSymbolic", "A Pauli string with a symbolic (polynomial) coefficient, stored in binary symplectic form (x and y bit vectors).")
                .def(nb::init<>(), "The zero Pauli string (identity operator, coefficient 0).")
                .def(nb::init<Expression, const std::unordered_map<int, std::string>&>(),
                "From a symbolic coefficient (Expression) and a map {qubit: 'X'/'Y'/'Z'/'I'}.")
                .def(nb::init<Expression&, const std::string&>(),
                "From a symbolic coefficient (Expression) and a text: dense such as 'XIZ' (character i acts on qubit i) or sparse such as 'X0 Z2'.")
                .def(nb::init<const std::pair<Expression, std::vector<std::pair<char, int>>>>(),
                "From an OpenFermion-style pair (Expression, [('X', qubit), ...]).")
                .def("to_string", &PauliString<Expression>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("is_all_z", &PauliString<Expression>::is_all_z, "True if every non-identity operator is Z.")
                .def("get_coeff", &PauliString<Expression>::get_coeff, "The coefficient of the Pauli string.")
                .def("to_dictionary", &PauliString<Expression>::to_dictionary, "The Pauli string as a (coefficient, {qubit: 'X'/'Y'/'Z'}) pair.")
                .def("commutator", &PauliString<Expression>::commutator, "Commutator [self, other] = self*other - other*self (a single Pauli string, zero if they commute).")
                .def("targeted_commutator", &PauliString<Expression>::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of target's operator in [self, other] (0 if absent); target's coefficient is ignored.")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &PauliString<Expression>::trace_out_qubits, nb::const_),
                        "Project onto a computational-basis state of the listed qubits (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &PauliString<Expression>::trace_out_qubits, nb::const_),
                        "Project onto general single-qubit states (a, b) representing a|0>+b|1> per traced qubit.")
                .def("map_qubits", &PauliString<Expression>::map_qubits, "Copy with qubit indices remapped by {old: new}; an occupied target that is not itself remapped is swapped.")
                // Type homogeneous multiplication.
                .def("__mul__", nb::overload_cast<const PauliString<Expression>&>(&PauliString<Expression>::operator*, nb::const_), "Product of two Pauli strings, including the phase of the Pauli algebra.")
                .def("__mul__", nb::overload_cast<double>(&PauliString<Expression>::operator*, nb::const_), "The Pauli string scaled by a float.")
                .def("__mul__", nb::overload_cast<int>(&PauliString<Expression>::operator*, nb::const_), "The Pauli string scaled by a int.")
                .def("__mul__", nb::overload_cast<const std::complex<double>>(&PauliString<Expression>::operator*, nb::const_), "The Pauli string scaled by a complex number.")
                .def("__rmul__", [](const PauliString<Expression>& ps, double s){ return ps * s; }, "The Pauli string scaled by a float (scalar on the left).")
                .def("__rmul__", [](const PauliString<Expression>& ps, int s){ return ps * s; }, "The Pauli string scaled by a int (scalar on the left).")
                .def("__rmul__", [](const PauliString<Expression>& ps, const std::complex<double> s){ return ps * s; }, "The Pauli string scaled by a complex number (scalar on the left).")
                .def("__imul__", nb::overload_cast<const std::complex<double>>(&PauliString<Expression>::operator*=), "Scale in place by a complex number.")
                .def("__imul__", nb::overload_cast<const PauliString<Expression>&>(&PauliString<Expression>::operator*=), "Multiply in place by another Pauli string.")
                // Type heterogeneous multiplication with complex coefficients.
                .def("__mul__", [](const PauliString<Expression>& ps1, const PauliString<std::complex<double>>& ps2){return ps1 * ps2;}, "Product with a complex Pauli string; the result is symbolic.")
                .def("__rmul__", [](const PauliString<Expression>& self, const PauliString<std::complex<double>>& other){return other * self;}, "Product other * self with a complex Pauli string; the result is symbolic.")
                .def("__repr__",  &PauliString<Expression>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__str__",  &PauliString<Expression>::to_string, "Human-readable form, e.g. '1.000000 + 0.000000i : X(0)Z(2)'.")
                .def("__hash__", [](const PauliString<Expression>& ps) {
                        return static_cast<Py_ssize_t>(ps.hash());
                }, "Hash over the operator part (x, y); ignores the coefficient.")
                .def("__neg__", [](const PauliString<Expression>& ps){ return ps * -1; }, "The Pauli string with negated coefficient.")
                .def("__add__", [](const PauliString<Expression>& ps1, const PauliString<Expression>& ps2){
                        return QubitHamiltonian<Expression>({ps1, ps2});
                }, "Add two PauliStringSymbolic to get PauliStringSumSymbolic")
                .def("__add__", [](const PauliString<Expression>& ps1, const PauliString<std::complex<double>>& ps2){
                        Expression sym_coeff = ps2.coeff;
                        PauliString<Expression> ps2_sym(ps2.x, ps2.y, sym_coeff);
                        return QubitHamiltonian<Expression>({ps1, ps2_sym});
                }, "Add PauliStringSymbolic and PauliStringComplex to get PauliStringSumSymbolic")
                .def("__sub__", [](const PauliString<Expression>& ps1, const PauliString<Expression>& ps2){
                        return QubitHamiltonian<Expression>({ps1, ps2 * -1});
                }, "Subtract two PauliStringSymbolic to get PauliStringSumSymbolic")
                .def("__sub__", [](const PauliString<Expression>& ps1, const PauliString<std::complex<double>>& ps2){
                        Expression sym_coeff = ps2.coeff;
                        PauliString<Expression> ps2_sym(ps2.x, ps2.y, sym_coeff);
                        return QubitHamiltonian<Expression>({ps1, ps2_sym * -1});
                }, "Subtract PauliStringSymbolic and PauliStringComplex to get PauliStringSumSymbolic")
                .def("diff", &PauliString<Expression>::diff, "Derivative of the coefficient with respect to the named symbol.")
                .def("subs", [](const PSS& x, const SubsMap& v) { return subs_any(x, v); }, nb::arg("values"), kSubsDoc)
                .def("qubits", &PauliString<Expression>::qubits, "Sorted list of the qubits with a non-identity operator.")
                .def("equals", &PauliString<Expression>::equals, "True if the operators are equal, ignoring the coefficients.")
                .def("set_coeff", &PauliString<Expression>::set_coeff, "Replace the coefficient (in place) by an Expression.")
                .def_static("to_complex", nb::overload_cast<const Expression&>(&PauliString<Expression>::to_complex), "Evaluate a fully numeric Expression to a complex number; raises RuntimeError if it still contains symbols.")
                .def("copy", &PauliString<Expression>::copy, "An independent copy of the Pauli string.")
                .def("get_pauli_at_index", &PauliString<Expression>::get_pauli_from_index, "Pauli operator ('I', 'X', 'Y' or 'Z') on the given qubit.")
                .def("naked", &PauliString<Expression>::naked, "Copy of the Pauli string with coefficient 1.")
                .def("size", &PauliString<Expression>::size, "Number of qubits with a non-identity operator (the weight).")
                .def("__len__", &PauliString<Expression>::size, "Number of qubits with a non-identity operator (the weight).")
                .def("count_y", &PauliString<Expression>::count_y, "Number of Y operators in the Pauli string.")
                .def("get_hash", [](const PauliString<Expression>& ps) { return hash_limbs_to_int(ps.get_hash()); },
                        "Invertible int key of the operator part (coefficient ignored): sum_q p_q * 4**q with I=0, X=1, Y=2, Z=3.")
                .def_static("from_hash", [](const nb::int_& hash, const Expression& coeff) {
                                return PauliString<Expression>::from_hash(hash_int_to_limbs(hash), coeff);
                        }, nb::arg("hash"), nb::arg("coeff"),
                        "Build the Pauli string whose get_hash() equals hash, with the given coefficient.")
                .def_static("from_hash", [](const nb::int_& hash) {
                                return PauliString<Expression>::from_hash(hash_int_to_limbs(hash), Expression(1));
                        }, nb::arg("hash"),
                        "Build the Pauli string whose get_hash() equals hash, with coefficient 1.")
                .def("key_openfermion", &PauliString<Expression>::key_openfermion, "OpenFermion-style key: list of (Pauli, qubit) pairs of the non-identity sites, sorted by qubit.")
                .def_prop_ro("x", [](const PauliString<Expression>& ps) { return ps.x.to_vector(); }, "x bits as 64-bit words (bit q set for X or Z on qubit q).")
                .def_prop_ro("y", [](const PauliString<Expression>& ps) { return ps.y.to_vector(); }, "y bits as 64-bit words (bit q set for Y or Z on qubit q).")
                .def_ro("coeff", &PauliString<Expression>::coeff, "The coefficient of the Pauli string (read-only).")
                .def_prop_ro("is_zero", [](const PauliString<Expression>& ps) { return ps.coeff.is_zero(); },
                        "True if the coefficient is (symbolically) zero.");

        // QubitHamiltonian for complex coefficients
        auto cls_qhc = nb::class_<QubitHamiltonian<std::complex<double>>>(m, "PauliStringSumComplex", "A sum of Pauli strings (a qubit Hamiltonian) with complex coefficients. Kept compacted: no two terms share an operator and zero terms are dropped.")
                .def(nb::init<const std::vector<PauliString<std::complex<double>>>&>(), "From a list of Pauli strings (equal operators are merged).")
                .def("__init__", [](QubitHamiltonian<std::complex<double>>* self, const std::vector<AnyString<std::complex<double>>>& data) {
                                new (self) QubitHamiltonian<std::complex<double>>(dense_strings(data));
                        }, nb::arg("data"), "From a list of Pauli strings, sparse ones or mixed (equal operators are merged).")
                .def(nb::init<const Hamiltonian_structure<std::complex<double>>&>(), "From a list of (coefficient, {qubit: 'X'/'Y'/'Z'}) pairs.")
                .def("contains", [](const QubitHamiltonian<std::complex<double>>& qh, const PauliString<std::complex<double>>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
                .def("contains", [](const QubitHamiltonian<std::complex<double>>& qh, const PauliString<Expression>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
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
                .def("targeted_commutator", [](const QHC& h, const PSC& other, const PSC& target) {
                                return h.targeted_commutator(QHC(std::vector<PSC>{other}), target);
                        }, nb::arg("other"), nb::arg("target"),
                        "Coefficient of the target Pauli string in [self, other] for a single Pauli string other.")
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

        // QubitHamiltonian for Expression coefficients
        auto cls_qhs = nb::class_<QubitHamiltonian<Expression>>(m, "PauliStringSumSymbolic", "A sum of Pauli strings (a qubit Hamiltonian) with symbolic (polynomial) coefficients. Kept compacted: no two terms share an operator and zero terms are dropped.")
                .def(nb::init<const std::vector<PauliString<Expression>>&>(), "From a list of Pauli strings (equal operators are merged).")
                .def("__init__", [](QubitHamiltonian<Expression>* self, const std::vector<AnyString<Expression>>& data) {
                                new (self) QubitHamiltonian<Expression>(dense_strings(data));
                        }, nb::arg("data"), "From a list of Pauli strings, sparse ones or mixed (equal operators are merged).")
                .def(nb::init<const Hamiltonian_structure<Expression>&>(), "From a list of (coefficient, {qubit: 'X'/'Y'/'Z'}) pairs.")
                .def("contains", [](const QubitHamiltonian<Expression>& qh, const PauliString<Expression>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
                .def("contains", [](const QubitHamiltonian<Expression>& qh, const PauliString<std::complex<double>>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the given Pauli string in the Hamiltonian (0 if absent); its own coefficient is ignored.")
                .def("contains", [](const QubitHamiltonian<Expression>& qh, const std::unordered_map<int, std::string>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string {qubit: 'X'} in the Hamiltonian (0 if absent).")
                .def("contains", [](const QubitHamiltonian<Expression>& qh, const std::vector<std::pair<char, int>>& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string [('X', 0), ...] in the Hamiltonian (0 if absent).")
                .def("contains", [](const QubitHamiltonian<Expression>& qh, const std::string& ps) { return qh.contains(ps); }, nb::arg("pauli_string"), "Coefficient of the Pauli string 'XIZ' or 'X0 Z2' in the Hamiltonian (0 if absent).")
                .def("__contains__", [](const QubitHamiltonian<Expression>& qh, const PauliString<Expression>& ps) {
                                return coeff_nonzero(qh.contains(ps));
                        }, nb::arg("pauli_string"), "True if a term with this Pauli operator is present (coefficient ignored).")
                .def("__add__", &QubitHamiltonian<Expression>::operator+, "Adds two Hamiltonians.")
                .def("__mul__", nb::overload_cast<std::complex<double>>(&QubitHamiltonian<Expression>::operator*, nb::const_), "Scales the Hamiltonian by a complex scalar.")
                .def("__mul__", nb::overload_cast<double>(&QubitHamiltonian<Expression>::operator*, nb::const_), "Scales the Hamiltonian by a float.")
                .def("__mul__", nb::overload_cast<int>(&QubitHamiltonian<Expression>::operator*, nb::const_), "Scales the Hamiltonian by an int.")
                .def("__mul__", nb::overload_cast<const QubitHamiltonian<Expression>&>(&QubitHamiltonian<Expression>::operator*, nb::const_), "Multiplies two Hamiltonians.")
                .def("__rmul__", [](const QubitHamiltonian<Expression>& qh, std::complex<double> s){ return qh * s; }, "The Hamiltonian scaled by a complex number (scalar on the left).")
                .def("__rmul__", [](const QubitHamiltonian<Expression>& qh, double s){ return qh * s; }, "The Hamiltonian scaled by a float (scalar on the left).")
                .def("__rmul__", [](const QubitHamiltonian<Expression>& qh, int s){ return qh * s; }, "The Hamiltonian scaled by an int (scalar on the left).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<int>&>(
                                &QubitHamiltonian<Expression>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto computational-basis states (state[i] in {0, 1}).")
                .def("trace_out_qubits",
                        nb::overload_cast<const std::vector<int>&, const std::vector<std::pair<std::complex<double>, std::complex<double>>>&>(
                                &QubitHamiltonian<Expression>::trace_out_qubits, nb::const_),
                        "Trace out qubits projected onto arbitrary single-qubit states (a, b) per traced qubit.")
                .def("to_string", &QubitHamiltonian<Expression>::to_string, "Human-readable form, one term per line.")
                .def("subs", [](const QHS& x, const SubsMap& v) { return subs_any(x, v); }, nb::arg("values"), kSubsDoc)
                .def("to_dictionary", &QubitHamiltonian<Expression>::to_dictionary, "Return list of (coefficient, {qubit: Pauli}) pairs.")
                .def("__str__",&QubitHamiltonian<Expression>::to_string, "Human-readable form, one term per line.")
                .def("__repr__",&QubitHamiltonian<Expression>::to_string, "Human-readable form, one term per line.")
                .def("__sub__", nb::overload_cast<const QubitHamiltonian<Expression>&>(&QubitHamiltonian<Expression>::operator-, nb::const_), "Subtracts two Hamiltonians.")
                .def("__neg__", nb::overload_cast<>(&QubitHamiltonian<Expression>::operator-, nb::const_), "Negates Hamiltonian.")
                .def("__pow__", &QubitHamiltonian<Expression>::power, "Raise Hamiltonian to integer power.")
                .def("__len__", &QubitHamiltonian<Expression>::size, "Number of Pauli string terms.")
                .def("size", &QubitHamiltonian<Expression>::size, "Number of Pauli string terms.")
                .def("commutator", &QubitHamiltonian<Expression>::commutator, "Commutator [self, other] = self*other - other*self.")
                .def("targeted_commutator", &QubitHamiltonian<Expression>::targeted_commutator,
                        nb::arg("other"), nb::arg("target"),
                        "Coefficient of the target Pauli string in the commutator, via a hash-join, O(min(len(A), len(B))).")
                .def("targeted_commutator", [](const QHS& h, const PSS& other, const PSS& target) {
                                return h.targeted_commutator(QHS(std::vector<PSS>{other}), target);
                        }, nb::arg("other"), nb::arg("target"),
                        "Coefficient of the target Pauli string in [self, other] for a single Pauli string other.")
                .def("diff", &QubitHamiltonian<Expression>::diff, "Derivative of every coefficient with respect to the named symbol.")
                .def("compact", &QubitHamiltonian<Expression>::compact, "Merge terms with equal operators and drop zero terms, in place; returns the result. Every operation already does this.")
                .def("set_all_coeff", &QubitHamiltonian<Expression>::set_all_coeff, nb::arg("value"), "Return a copy with every term's coefficient replaced by value.")
                .def("simplify", &QubitHamiltonian<Expression>::simplify, nb::arg("threshold") = 0.0, "Copy without the terms whose abs(coefficient) <= threshold.")
                .def("qubits", &QubitHamiltonian<Expression>::qubits, "Sorted list of qubit indices the Hamiltonian acts on non-trivially.")
                .def("n_qubits", &QubitHamiltonian<Expression>::n_qubits, "Number of qubits the Hamiltonian acts on non-trivially.")
                .def("is_all_z", &QubitHamiltonian<Expression>::is_all_z, "True if all terms are products of Z operators only.")
                .def("count_measurements", &QubitHamiltonian<Expression>::count_measurements, "1 if is_all_z(), else the number of terms.")
                .def("is_hermitian", &QubitHamiltonian<Expression>::is_hermitian, "True if every coefficient is real (abs(imag) <= 1e-12).")
                .def("is_antihermitian", &QubitHamiltonian<Expression>::is_antihermitian, "True if every coefficient is purely imaginary (abs(real) <= 1e-12).")
                .def("dagger", &QubitHamiltonian<Expression>::dagger, "Hermitian conjugate: complex-conjugates each coefficient (Pauli strings are Hermitian).")
                .def("conjugate", &QubitHamiltonian<Expression>::conjugate, "Complex conjugate: conjugates each coefficient and flips its sign for an odd number of Y.")
                .def("transpose", &QubitHamiltonian<Expression>::transpose, "Transpose: flips the sign of terms with an odd number of Y.")
                .def("split", &QubitHamiltonian<Expression>::split, "Split into (hermitian, anti-hermitian) parts; needs numeric coefficients.")
                .def("map_qubits", &QubitHamiltonian<Expression>::map_qubits, "Copy with qubit indices remapped by {old: new} (see PauliString.map_qubits).")
                .def("power", &QubitHamiltonian<Expression>::power, "Raise Hamiltonian to integer power.")
                .def("paulistrings", &QubitHamiltonian<Expression>::paulistrings, "Return the list of Pauli string terms.")
                .def("to_list", &QubitHamiltonian<Expression>::to_list, "The list of Pauli string terms (same as paulistrings()).")
                .def("to_matrix", &QubitHamiltonian<Expression>::to_matrix, nb::arg("ignore_unused_qubits") = true, "Dense 2**N x 2**N matrix (list of rows). With ignore_unused_qubits, N counts only the qubits the Hamiltonian acts on; otherwise qubit indices are absolute. Needs numeric coefficients.")
                .def_static("zero", &QubitHamiltonian<Expression>::zero, "Return an empty Hamiltonian (the zero operator).")
                .def_static("unit", &QubitHamiltonian<Expression>::unit, "Return the identity Hamiltonian (single term, coefficient 1, no operators).");

        using cplx = std::complex<double>;
        nb::class_<Expression>(m, "Expression",
                "A symbolic coefficient: a polynomial in named symbols with complex coefficients, kept in "
                "canonical form (so equal polynomials compare equal and a - a is 0). conjugate() turns a symbol a "
                "into conjugate(a) unless it was declared real with set_symbol_real().")
                .def(nb::init<const std::string&>(), nb::arg("text"),
                        "Parse a polynomial such as 'a*b + 2', '2a - b**2/4' or '(1 + 2*I)*c'. I is the imaginary "
                        "unit; pi and E are numbers. Division is only allowed by numbers. Raises ValueError.")
                .def(nb::init<cplx>(), nb::arg("value"), "A numeric expression.")
                .def("__str__", &Expression::str, "The expression as a string, which Expression(...) parses back exactly.")
                .def("__repr__", &Expression::str, "The expression as a string, which Expression(...) parses back exactly.")
                .def("__eq__", [](const Expression& a, const Expression& b) { return a == b; }, nb::is_operator(), "Equal polynomials.")
                .def("__eq__", [](const Expression& a, cplx b) { return a == Expression(b); }, nb::is_operator(), "Equal to the number.")
                .def("__hash__", [](const Expression& e) { return static_cast<Py_ssize_t>(e.hash()); })
                .def("__add__", [](const Expression& a, const Expression& b) { return a + b; }, nb::is_operator())
                .def("__add__", [](const Expression& a, cplx b) { return a + b; }, nb::is_operator())
                .def("__radd__", [](const Expression& a, cplx b) { return b + a; }, nb::is_operator())
                .def("__sub__", [](const Expression& a, const Expression& b) { return a - b; }, nb::is_operator())
                .def("__sub__", [](const Expression& a, cplx b) { return a - b; }, nb::is_operator())
                .def("__rsub__", [](const Expression& a, cplx b) { return b - a; }, nb::is_operator())
                .def("__mul__", [](const Expression& a, const Expression& b) { return a * b; }, nb::is_operator())
                .def("__mul__", [](const Expression& a, cplx b) { return a * b; }, nb::is_operator())
                .def("__rmul__", [](const Expression& a, cplx b) { return b * a; }, nb::is_operator())
                .def("__neg__", [](const Expression& a) { return -a; })
                .def("__pow__", [](const Expression& a, unsigned n) { return a.pow(n); }, "Non-negative integer power.")
                .def("__complex__", &Expression::to_complex, "The value; raises RuntimeError if symbols are left.")
                .def("to_complex", &Expression::to_complex, "The value; raises RuntimeError if symbols are left.")
                .def("diff", nb::overload_cast<const std::string&>(&Expression::diff, nb::const_), nb::arg("symbol"),
                        "Derivative with respect to the named symbol; conjugate(symbol) counts as the same variable "
                        "unless the symbol was declared complex.")
                .def("subs", [](const Expression& e, const SubsMap& values) { return subs_any(ExpressionSubs{e}, values); },
                        nb::arg("values"), kSubsDoc)
                .def("conjugate", [](const Expression& e) { return e.conjugate(); },
                        "Complex conjugate: every symbol a becomes conjugate(a), except symbols declared real.")
                .def_prop_ro("is_zero", &Expression::is_zero, "True for the zero polynomial.")
                .def_prop_ro("is_number", &Expression::is_number, "True if no symbols are left.")
                .def("free_symbols", &Expression::free_symbols, "Names of the symbols that occur, sorted.")
                .def("n_terms", &Expression::n_terms, "Number of monomials.");

        m.def("set_symbol_real", &pauliengine::symbolic::set_symbol_real, nb::arg("name"), nb::arg("is_real") = true,
                "Declare a symbol real (with False: undeclared again). By default conjugate() and dagger() turn a symbol "
                "a into conjugate(a), and subs({'a': v}) puts conj(v) in for it, so the result is right for real and "
                "complex values. A real symbol stays a, so a Hermitian sum with real symbols equals its dagger "
                "before anything is substituted.");
        m.def("set_symbol_complex", &pauliengine::symbolic::set_symbol_complex, nb::arg("name"), nb::arg("is_complex") = true,
                "Declare a symbol complex (with False: real). It is conjugated like an undeclared symbol, but diff "
                "treats a and conjugate(a) as independent variables (Wirtinger derivatives); for an undeclared symbol "
                "diff('a') also differentiates conjugate(a).");
        m.def("is_symbol_complex", &pauliengine::symbolic::is_symbol_complex, nb::arg("name"),
                "True if the symbol was declared complex.");
        m.def("is_symbol_real", &pauliengine::symbolic::is_symbol_real, nb::arg("name"),
                "True if the symbol was declared real.");

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
                .def("__init__", [](unweighted_pauliString* self, const SPC& s, int phase) {
                                new (self) unweighted_pauliString(s.to_dense(), phase);
                        }, nb::arg("string"), nb::arg("phase") = 0, "The operator part of a sparse Pauli string; its coefficient is dropped.")
                .def("__init__", [](unweighted_pauliString* self, const SPS& s, int phase) {
                                new (self) unweighted_pauliString(s.to_dense(), phase);
                        }, nb::arg("string"), nb::arg("phase") = 0,
                        "The operator part of a symbolic sparse Pauli string; its coefficient is dropped.")
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
                .def("to_pauli_string", [](const unweighted_pauliString& s, const std::string& coeff) { return s.to_pauli_string(Expression(coeff)); },
                        nb::arg("coeff"), "The symbolic Pauli string with coefficient coeff * i**phase; coeff is parsed like 'a*b'.")
                .def("to_pauli_string", [](const unweighted_pauliString& s, const Expression& coeff) { return s.to_pauli_string(coeff); },
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
        auto cls_uq = nb::class_<unweighted_QubitHamiltonian>(m, "unweighted_PauliStringSum",
                "An ordered set of Pauli operators with one phase i**k per element. Use it only as the basis of a Lie algebra. "
                "It is not a linear combination: every operator occurs once and the first phase seen for it wins, so A + A == A. "
                "An element keeps its index.")
                .def(nb::init<>(), "The empty set.")
                .def(nb::init<const std::vector<unweighted_pauliString>&>(), nb::arg("terms"),
                        "From a list of strings; later duplicates of an operator are ignored.")
                .def(nb::init<const QHC&>(), nb::arg("hamiltonian"), "The operators of a Hamiltonian with phase 0; its coefficients are dropped.")
                .def(nb::init<const QHS&>(), nb::arg("hamiltonian"), "The operators of a symbolic Hamiltonian with phase 0; its coefficients are dropped.")
                .def("__init__", [](unweighted_QubitHamiltonian* self, const SQC& h) { new (self) unweighted_QubitHamiltonian(h.to_dense()); },
                        nb::arg("hamiltonian"), "The operators of a sparse sum with phase 0; its coefficients are dropped.")
                .def("__init__", [](unweighted_QubitHamiltonian* self, const SQS& h) { new (self) unweighted_QubitHamiltonian(h.to_dense()); },
                        nb::arg("hamiltonian"), "The operators of a symbolic sparse sum with phase 0; its coefficients are dropped.")
                .def("__len__", &unweighted_QubitHamiltonian::size, "Number of elements.")
                .def("size", &unweighted_QubitHamiltonian::size, "Number of elements.")
                .def("insert", &unweighted_QubitHamiltonian::insert, nb::arg("term"),
                        "Add the string unless its operator is present; returns whether it was added.")
                .def("__contains__", &unweighted_QubitHamiltonian::contains, nb::arg("term"),
                        "True if the operator of the string is present, whatever the phases are.")
                .def("index", [](const unweighted_QubitHamiltonian& h, const unweighted_pauliString& term) {
                                const std::optional<size_t> i = h.find(term);
                                if (!i) {
                                        throw nb::value_error("unweighted_PauliStringSum: operator is not present");
                                }
                                return *i;
                        }, nb::arg("term"), "Index of the operator of the string; raises ValueError if it is not present.")
                .def("__getitem__", [](const unweighted_QubitHamiltonian& h, Py_ssize_t index) {
                                const Py_ssize_t n = static_cast<Py_ssize_t>(h.size());
                                if (index < -n || index >= n) {
                                        throw nb::index_error("unweighted_PauliStringSum: index out of range");
                                }
                                return h.terms()[static_cast<size_t>(index < 0 ? index + n : index)];
                        }, nb::arg("index"), "The element with the given index.")
                .def("__iter__", [](const unweighted_QubitHamiltonian& h) { return nb::iter(nb::cast(h.terms())); }, "Iterate over the elements in insertion order.")
                .def("terms", &unweighted_QubitHamiltonian::terms, "All elements in insertion order.")
                .def("paulistrings", &unweighted_QubitHamiltonian::terms, "All elements in insertion order, as PauliStringSum.paulistrings().")
                .def("__add__", &unweighted_QubitHamiltonian::operator+, "Union: the elements of self, then the new ones of other.")
                .def("__iadd__", [](unweighted_QubitHamiltonian& h, const unweighted_QubitHamiltonian& other) -> unweighted_QubitHamiltonian& {
                                h += other;
                                return h;
                        }, nb::rv_policy::none, "Add the new elements of other in place; existing indices stay.")
                .def("__mul__", &unweighted_QubitHamiltonian::operator*, "The set of all products a * b in the order of a, then b.")
                .def("commutator", &unweighted_QubitHamiltonian::commutator, nb::arg("other"),
                        "The set of commutators without their factor 2: a * b for every anti-commuting pair.")
                .def("__eq__", &unweighted_QubitHamiltonian::operator==, "Equal sets with equal phases, whatever the order of the elements is.")
                .def("to_pauli_string_sum", [](const unweighted_QubitHamiltonian& h, std::complex<double> coeff) { return h.to_qubit_hamiltonian(coeff); },
                        nb::arg("coeff") = std::complex<double>(1.0, 0.0), "The ordinary Hamiltonian with coefficient coeff * i**k on every element.")
                .def("to_pauli_string_sum", [](const unweighted_QubitHamiltonian& h, const std::string& coeff) { return h.to_qubit_hamiltonian(Expression(coeff)); },
                        nb::arg("coeff"), "The symbolic Hamiltonian with coefficient coeff * i**k on every element; coeff is parsed like 'a*b'.");

        // PauliCycle and PauliCycleSum, with complex or symbolic coefficients
        auto cls_cy = bind_cycle<std::complex<double>>(m, "PauliCycleComplex",
                "A Pauli string with a complex coefficient together with all of its cyclic rotations on "
                "n_qubits qubits, C(P) = coeff * sum_k R^k P. The sum runs over all n rotations, so a base with a "
                "shorter period counts each distinct rotation several times.");
        auto cls_cys_sym = bind_cycle<Expression>(m, "PauliCycleSymbolic",
                "A Pauli string with a symbolic (polynomial) coefficient together with all of its cyclic rotations "
                "on n_qubits qubits, C(P) = coeff * sum_k R^k P.");
        auto cls_cys = bind_cycle_sum<std::complex<double>>(m, "PauliCycleSumComplex",
                "A sum of PauliCycles with complex coefficients; each cycle carries its own coefficient. "
                "Summands from the same rotation orbit are always merged.");
        auto cls_css = bind_cycle_sum<Expression>(m, "PauliCycleSumSymbolic",
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
                .def("to_pauli_string_sum", &PauliOrbit::to_qubit_hamiltonian,
                        "The orbit expanded into a PauliStringSum, (i / n_terms) * sum of its Pauli strings. "
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
        auto cls_oss = bind_orbit_sum<Expression>(m, "PauliOrbitSumSymbolic",
                "A weighted sum of Pauli orbits with symbolic coefficients, sum_i coeffs[i] * data[i], on a fixed "
                "number of qubits.");

        // Sparse Pauli strings and Hamiltonians (SparsePauliString.h, SparsePauliStringSum.h)
        auto cls_spc = bind_sparse_string<std::complex<double>>(m, "SparsePauliStringComplex",
                "A Pauli string with a complex coefficient, stored as the sorted list of its non-identity sites. "
                "Costs grow with the weight instead of the highest qubit; meant for low-weight strings on many qubits.");
        auto cls_sps = bind_sparse_string<Expression>(m, "SparsePauliStringSymbolic",
                "A Pauli string with a symbolic coefficient, stored as the sorted list of its non-identity sites.");
        auto cls_sqc = bind_sparse_hamiltonian<std::complex<double>>(m, "SparsePauliStringSumComplex",
                "A sum of sparse Pauli strings with complex coefficients. Products, commutators and squares cost "
                "O(w) per pair of terms of weight w instead of O(n / 64) words. Kept compacted.");
        auto cls_sqs = bind_sparse_hamiltonian<Expression>(m, "SparsePauliStringSumSymbolic",
                "A sum of sparse Pauli strings with symbolic coefficients. Kept compacted.");

        // Scaling cycles and orbits by an int or float.
        bind_scale_ops<CYC>(cls_cy, "The cycle with its coefficient scaled by an int, float or complex number.");
        bind_scale_ops<CYS>(cls_cys_sym, "The cycle with its coefficient scaled by an int, float or complex number.");
        bind_scale_ops<CSC>(cls_cys,
                "The sum with every coefficient scaled by an int, float or complex number; scaling by 0 gives the empty sum.");
        bind_scale_ops<CSS>(cls_css,
                "The sum with every coefficient scaled by an int, float or complex number; scaling by 0 gives the empty sum.");
        bind_scale_ops<PauliOrbit>(cls_or,
                "scale * self as a PauliOrbitSum, for an int, float or complex number; scaling by 0 gives the empty sum.");
        bind_scale_ops<OSC>(cls_ors,
                "The sum with every coefficient scaled by an int, float or complex number; scaling by 0 gives the empty sum.");
        bind_scale_ops<SPC>(cls_spc, "The string with its coefficient scaled by an int, float or complex number.");
        bind_scale_ops<SPS>(cls_sps, "The string with its coefficient scaled by an int, float or complex number.");
        bind_scale_ops<SQC>(cls_sqc, "Every coefficient scaled by an int, float or complex number; scaling by 0 gives the zero operator.");
        bind_scale_ops<SQS>(cls_sqs, "Every coefficient scaled by an int, float or complex number; scaling by 0 gives the zero operator.");
        bind_scale_ops<OSS>(cls_oss,
                "The sum with every coefficient scaled by an int, float or complex number; scaling by 0 gives the empty sum.");

        // Derivative and substitution of symbolic cycles and orbit sums, as for Pauli strings.
        cls_cys_sym.def("diff", &CYS::diff, nb::arg("symbol"), "Derivative of the coefficient with respect to the named symbol.")
                   .def("subs", [](const CYS& x, const SubsMap& v) { return subs_any(x, v); }, nb::arg("values"), kSubsDoc);
        cls_css.def("diff", &CSS::diff, nb::arg("symbol"),
                    "Derivative of every coefficient with respect to the named symbol; summands that become zero are dropped.")
               .def("subs", [](const CSS& x, const SubsMap& v) { return subs_any(x, v); }, nb::arg("values"),
                    "Substitute values for symbols (see PauliStringSymbolic.subs); summands that become zero are dropped.");
        cls_oss.def("diff", &OSS::diff, nb::arg("symbol"),
                    "Derivative of every coefficient with respect to the named symbol; terms that become zero are dropped.")
               .def("subs", [](const OSS& x, const SubsMap& v) { return subs_any(x, v); }, nb::arg("values"),
                    "Substitute values for symbols (see PauliStringSymbolic.subs); terms that become zero are dropped.");

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

        // Methods that do not depend on the representation exist on both sides, so
        // that code works whether the factories chose bit vectors or sparse sites.
        cls_sqc.def_static("from_matrix",
                [](nb::ndarray<const std::complex<double>, nb::ndim<2>, nb::c_contig, nb::device::cpu> m,
                   double threshold) {
                        if (m.shape(0) != m.shape(1)) {
                                throw std::invalid_argument("from_matrix: the matrix must be square");
                        }
                        const size_t dim = m.shape(0);
                        std::vector<std::complex<double>> data(m.data(), m.data() + dim * dim);
                        return SQC(pauli_decomposition(std::move(data), dim, threshold));
                }, nb::arg("matrix"), nb::arg("threshold") = 0.0,
                "Pauli decomposition of a 2**n x 2**n matrix, as PauliStringSumComplex.from_matrix, stored sparse.");
        cls_sqc.def_static("from_matrix",
                [](const Matrix2D& m, double threshold) { return SQC(pauli_decomposition(m, threshold)); },
                nb::arg("matrix"), nb::arg("threshold") = 0.0,
                "Pauli decomposition of a 2**n x 2**n matrix given as a list of rows, stored sparse.");
        cls_qhc.def("square", [](const QHC& h) { return h.square(); }, "self * self, visiting each commuting pair once.");
        cls_qhs.def("square", [](const QHS& h) { return h.square(); }, "self * self, visiting each commuting pair once.");
        cls_psc.def("anticommutes", [](const PSC& a, const PSC& b) {
                        return pauliengine_detail::anticommute_words(a.x.data(), a.y.data(), a.x.size(),
                                                                     b.x.data(), b.y.data(), b.x.size());
                }, nb::arg("other"), "True if the strings anti-commute.");
        cls_pss.def("anticommutes", [](const PSS& a, const PSS& b) {
                        return pauliengine_detail::anticommute_words(a.x.data(), a.y.data(), a.x.size(),
                                                                     b.x.data(), b.y.data(), b.x.size());
                }, nb::arg("other"), "True if the strings anti-commute.");
        cls_spc.def_static("from_hash", [](const nb::int_& hash, std::complex<double> coeff) {
                        return SPC(PSC::from_hash(hash_int_to_limbs(hash), coeff));
                }, nb::arg("hash"), nb::arg("coeff") = std::complex<double>(1.0, 0.0),
                "Build the sparse Pauli string whose get_hash() equals hash, with the given coefficient.");
        cls_sps.def_static("from_hash", [](const nb::int_& hash, const Expression& coeff) {
                        return SPS(PSS::from_hash(hash_int_to_limbs(hash), coeff));
                }, nb::arg("hash"), nb::arg("coeff"),
                "Build the sparse Pauli string whose get_hash() equals hash, with the given coefficient.");
        cls_sps.def_static("from_hash", [](const nb::int_& hash) {
                        return SPS(PSS::from_hash(hash_int_to_limbs(hash), Expression(1)));
                }, nb::arg("hash"), "Build the sparse Pauli string whose get_hash() equals hash, with coefficient 1.");
        cls_sps.def_static("to_complex", nb::overload_cast<const Expression&>(&PSS::to_complex),
                "Evaluate a fully numeric Expression to a complex number; raises RuntimeError if it still contains symbols.");

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
        bind_arithmetic<SPC>(cls_spc, AllOperators{});
        bind_arithmetic<SPS>(cls_sps, AllOperators{});
        bind_arithmetic<SQC>(cls_sqc, AllOperators{});
        bind_arithmetic<SQS>(cls_sqs, AllOperators{});

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
        bind_unweighted_pair<SPC, UPS>(cls_spc, cls_up);
        bind_unweighted_pair<SPS, UPS>(cls_sps, cls_up);
        bind_unweighted_pair<SQC, UPS>(cls_sqc, cls_up);
        bind_unweighted_pair<SQS, UPS>(cls_sqs, cls_up);
        bind_unweighted_pair<SPC, UQH>(cls_spc, cls_uq);
        bind_unweighted_pair<SPS, UQH>(cls_sps, cls_uq);
        bind_unweighted_pair<SQC, UQH>(cls_sqc, cls_uq);
        bind_unweighted_pair<SQS, UQH>(cls_sqs, cls_uq);
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

        // Division by a number and multiplication by a symbolic expression, for every kind.
        const auto same = [](const auto& a) -> const auto& { return a; };
        const auto weighted = [](const auto& a) { return as_weighted(a); };
        bind_division<PSC>(cls_psc, same);
        bind_division<PSS>(cls_pss, same);
        bind_division<QHC>(cls_qhc, same);
        bind_division<QHS>(cls_qhs, same);
        bind_division<CYC>(cls_cy, same);
        bind_division<CYS>(cls_cys_sym, same);
        bind_division<CSC>(cls_cys, same);
        bind_division<CSS>(cls_css, same);
        bind_division<PauliOrbit>(cls_or, same);
        bind_division<OSC>(cls_ors, same);
        bind_division<OSS>(cls_oss, same);
        bind_division<UPS>(cls_up, weighted);
        bind_division<UQH>(cls_uq, weighted);
        bind_division<SPC>(cls_spc, same);
        bind_division<SPS>(cls_sps, same);
        bind_division<SQC>(cls_sqc, same);
        bind_division<SQS>(cls_sqs, same);
        bind_symbolic_scale<PSC>(cls_psc);
        bind_symbolic_scale<PSS>(cls_pss);
        bind_symbolic_scale<QHC>(cls_qhc);
        bind_symbolic_scale<QHS>(cls_qhs);
        bind_symbolic_scale<CYC>(cls_cy);
        bind_symbolic_scale<CYS>(cls_cys_sym);
        bind_symbolic_scale<CSC>(cls_cys);
        bind_symbolic_scale<CSS>(cls_css);
        bind_symbolic_scale<PauliOrbit>(cls_or);
        bind_symbolic_scale<OSC>(cls_ors);
        bind_symbolic_scale<OSS>(cls_oss);
        bind_symbolic_scale<UPS>(cls_up);
        bind_symbolic_scale<UQH>(cls_uq);
        bind_symbolic_scale<SPC>(cls_spc);
        bind_symbolic_scale<SPS>(cls_sps);
        bind_symbolic_scale<SQC>(cls_sqc);
        bind_symbolic_scale<SQS>(cls_sqs);

        // Two unweighted operators: inner product of their weighted forms; -set.
        const char* unweighted_ip = "Normalized Hilbert-Schmidt inner product; both count with coefficient 1.0 * i**phase.";
        cls_up.def("inner_product", [](const UPS& a, const UPS& b) { return inner_product(a, b); }, nb::arg("other"), unweighted_ip)
              .def("inner_product", [](const UPS& a, const UQH& b) { return inner_product(a, b); }, nb::arg("other"), unweighted_ip);
        cls_uq.def("inner_product", [](const UQH& a, const UPS& b) { return inner_product(a, b); }, nb::arg("other"), unweighted_ip)
              .def("inner_product", [](const UQH& a, const UQH& b) { return inner_product(a, b); }, nb::arg("other"), unweighted_ip)
              .def("__neg__", [](const UQH& h) { return -h; }, "Every element multiplied by -1, as a set.");

        // Equality. Complex and symbolic operators of the same kind compare by value
        // (the complex one counts as symbolic); anything else is not equal.
        const char* ps_eq = "True if both the operator and the coefficient are equal (see equals() for the operator only).";
        const auto sym_ps = [](const PSC& a) { return PSS(a.x, a.y, Expression(a.coeff)); };
        bind_eq<PSC, PSC>(cls_psc, [](const PSC& a, const PSC& b) { return a == b; }, ps_eq);
        bind_eq<PSC, PSS>(cls_psc, [sym_ps](const PSC& a, const PSS& b) { return sym_ps(a) == b; }, ps_eq);
        bind_eq<PSS, PSS>(cls_pss, [](const PSS& a, const PSS& b) { return a == b; }, ps_eq);
        bind_eq<PSS, PSC>(cls_pss, [sym_ps](const PSS& a, const PSC& b) { return a == sym_ps(b); }, ps_eq);
        const char* qh_eq = "True if both have the same terms with equal coefficients (order-independent).";
        const auto sym_qh = [](const QHC& h) { return pauliengine_detail::convert_coeffs<Expression>(h); };
        bind_eq<QHC, QHC>(cls_qhc, [](const QHC& a, const QHC& b) { return a == b; }, qh_eq);
        bind_eq<QHC, QHS>(cls_qhc, [sym_qh](const QHC& a, const QHS& b) { return sym_qh(a) == b; }, qh_eq);
        bind_eq<QHS, QHS>(cls_qhs, [](const QHS& a, const QHS& b) { return a == b; }, qh_eq);
        bind_eq<QHS, QHC>(cls_qhs, [sym_qh](const QHS& a, const QHC& b) { return a == sym_qh(b); }, qh_eq);
        // Sums of cycles or orbits are equal if their difference has no terms left
        // (cancelling coefficients are dropped exactly, symbolic ones after expanding).
        const char* sum_eq = "True if both have the same summands with equal coefficients.";
        const auto same_sum = [](const auto& a, const auto& b) {
                try {
                        return (a - b).size() == 0;
                } catch (const std::invalid_argument&) {
                        return false;  // different numbers of qubits
                }
        };
        bind_eq<CSC, CSC>(cls_cys, same_sum, sum_eq);
        bind_eq<CSC, CSS>(cls_cys, same_sum, sum_eq);
        bind_eq<CSS, CSS>(cls_css, same_sum, sum_eq);
        bind_eq<CSS, CSC>(cls_css, same_sum, sum_eq);
        bind_eq<OSC, OSC>(cls_ors, same_sum, sum_eq);
        bind_eq<OSC, OSS>(cls_ors, same_sum, sum_eq);
        bind_eq<OSS, OSS>(cls_oss, same_sum, sum_eq);
        bind_eq<OSS, OSC>(cls_oss, same_sum, sum_eq);
        // Mutable sums that compare by value cannot be hashed by identity.
        cls_cys.attr("__hash__") = nb::none();
        cls_css.attr("__hash__") = nb::none();
        cls_ors.attr("__hash__") = nb::none();
        cls_oss.attr("__hash__") = nb::none();
        // Sparse operators compare by value with sparse ones and with Pauli strings or
        // Hamiltonians of the same kind (string with string, sum with sum).
        const auto same_string = [](const auto& a, const auto& b) {
                using R = pauliengine_detail::promoted_coeff_t<std::decay_t<decltype(a)>, std::decay_t<decltype(b)>>;
                return pauliengine_detail::as_sparse_string<R>(a) == pauliengine_detail::as_sparse_string<R>(b);
        };
        const auto same_operator = [](const auto& a, const auto& b) {
                using R = pauliengine_detail::promoted_coeff_t<std::decay_t<decltype(a)>, std::decay_t<decltype(b)>>;
                return pauliengine_detail::as_sparse<R>(a) == pauliengine_detail::as_sparse<R>(b);
        };
        bind_eq<SPC, SPC>(cls_spc, same_string, ps_eq);
        bind_eq<SPC, SPS>(cls_spc, same_string, ps_eq);
        bind_eq<SPC, PSC>(cls_spc, same_string, ps_eq);
        bind_eq<SPC, PSS>(cls_spc, same_string, ps_eq);
        bind_eq<SPS, SPS>(cls_sps, same_string, ps_eq);
        bind_eq<SPS, SPC>(cls_sps, same_string, ps_eq);
        bind_eq<SPS, PSC>(cls_sps, same_string, ps_eq);
        bind_eq<SPS, PSS>(cls_sps, same_string, ps_eq);
        bind_eq<PSC, SPC>(cls_psc, same_string, ps_eq);
        bind_eq<PSC, SPS>(cls_psc, same_string, ps_eq);
        bind_eq<PSS, SPC>(cls_pss, same_string, ps_eq);
        bind_eq<PSS, SPS>(cls_pss, same_string, ps_eq);
        bind_eq<SQC, SQC>(cls_sqc, same_operator, qh_eq);
        bind_eq<SQC, SQS>(cls_sqc, same_operator, qh_eq);
        bind_eq<SQC, QHC>(cls_sqc, same_operator, qh_eq);
        bind_eq<SQC, QHS>(cls_sqc, same_operator, qh_eq);
        bind_eq<SQS, SQS>(cls_sqs, same_operator, qh_eq);
        bind_eq<SQS, SQC>(cls_sqs, same_operator, qh_eq);
        bind_eq<SQS, QHC>(cls_sqs, same_operator, qh_eq);
        bind_eq<SQS, QHS>(cls_sqs, same_operator, qh_eq);
        bind_eq<QHC, SQC>(cls_qhc, same_operator, qh_eq);
        bind_eq<QHC, SQS>(cls_qhc, same_operator, qh_eq);
        bind_eq<QHS, SQC>(cls_qhs, same_operator, qh_eq);
        bind_eq<QHS, SQS>(cls_qhs, same_operator, qh_eq);
        cls_sqc.attr("__hash__") = nb::none();
        cls_sqs.attr("__hash__") = nb::none();
        bind_eq_fallback<SPC>(cls_spc);
        bind_eq_fallback<SPS>(cls_sps);
        bind_eq_fallback<SQC>(cls_sqc);
        bind_eq_fallback<SQS>(cls_sqs);
        bind_eq_fallback<PSC>(cls_psc);
        bind_eq_fallback<PSS>(cls_pss);
        bind_eq_fallback<QHC>(cls_qhc);
        bind_eq_fallback<QHS>(cls_qhs);
        bind_eq_fallback<CYC>(cls_cy);
        bind_eq_fallback<CYS>(cls_cys_sym);
        bind_eq_fallback<CSC>(cls_cys);
        bind_eq_fallback<CSS>(cls_css);
        bind_eq_fallback<PauliOrbit>(cls_or);
        bind_eq_fallback<OSC>(cls_ors);
        bind_eq_fallback<OSS>(cls_oss);
        bind_eq_fallback<UPS>(cls_up);
        bind_eq_fallback<UQH>(cls_uq);
        // numpy scalars on the left: let Python call __rmul__ instead of numpy
        // treating the set as a sequence.
        cls_uq.attr("__array_ufunc__") = nb::none();
        // Old Python names (QubitHamiltonian...) stay available as aliases.
        for (const char* name : {"Complex", "Symbolic"}) {
                m.attr((std::string("QubitHamiltonian") + name).c_str()) = m.attr((std::string("PauliStringSum") + name).c_str());
        }
        m.attr("unweighted_QubitHamiltonian") = m.attr("unweighted_PauliStringSum");
        for (nb::handle cls : {nb::handle(cls_cy), nb::handle(cls_cys_sym), nb::handle(cls_cys), nb::handle(cls_css),
                               nb::handle(cls_or), nb::handle(cls_ors), nb::handle(cls_oss), nb::handle(cls_uq)}) {
                nb::setattr(cls, "to_qubit_hamiltonian", nb::getattr(cls, "to_pauli_string_sum"));
        }
}
