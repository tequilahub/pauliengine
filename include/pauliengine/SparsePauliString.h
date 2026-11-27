#pragma once
/// @file
/// A Pauli string stored as the sorted list of its non-identity sites.
///
/// PauliString stores two bit vectors up to the highest qubit it acts on, so a
/// string of weight 2 on qubit 4000 takes 126 words and every product walks all
/// of them. SparsePauliString stores only its w sites (4 bytes each); products and
/// commutators merge the two lists in O(w_a + w_b). This pays off for strings of
/// low weight on many qubits, the typical terms of large local Hamiltonians.

#include <algorithm>
#include <complex>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/detail/Coefficients.h"
#include "pauliengine/detail/SparseSites.h"

template <typename Coeff>
class SparsePauliString {
    public:
        using Site = pauliengine_detail::Site;

        /// The non-identity sites, by ascending qubit (see detail/SparseSites.h).
        std::vector<Site> sites;
        Coeff coeff = pauliengine_detail::coeff_one<Coeff>();

        /// coeff times the identity.
        explicit SparsePauliString(Coeff c = pauliengine_detail::coeff_one<Coeff>()) : coeff(std::move(c)) {}

        /// From a coefficient and a map {qubit: "X"/"Y"/"Z"/"I"}.
        SparsePauliString(Coeff c, const std::unordered_map<int, std::string>& paulis) : coeff(std::move(c)) {
                sites.reserve(paulis.size());
                for (const auto& [qubit, op] : paulis) {
                        if (qubit < 0 || static_cast<uint64_t>(qubit) > pauliengine_detail::kMaxSparseQubit) {
                                throw std::invalid_argument("SparsePauliString: qubit index out of range");
                        }
                        if (op.size() != 1) throw std::invalid_argument("SparsePauliString: invalid Pauli operator '" + op + "'");
                        const uint32_t p = pauliengine_detail::pauli_bits(op[0]);
                        if (p) sites.push_back(pauliengine_detail::make_site(static_cast<uint32_t>(qubit), p));
                }
                std::sort(sites.begin(), sites.end());
        }

        /// From a coefficient and sites (ascending qubits, no identities; checked).
        SparsePauliString(Coeff c, std::vector<Site> s) : sites(std::move(s)), coeff(std::move(c)) {
                pauliengine_detail::check_sites(sites.data(), sites.size());
        }

        /// The same string as a PauliString.
        explicit SparsePauliString(const PauliString<Coeff>& dense) : coeff(dense.coeff) {
                pauliengine_detail::sites_of_words(dense.x.data(), dense.y.data(),
                                                   std::min(dense.x.size(), dense.y.size()), sites);
        }

        /// The string as a PauliString (bit vectors up to its highest qubit).
        PauliString<Coeff> to_dense() const {
                std::vector<uint64_t> x, y;
                pauliengine_detail::words_of_sites(sites.data(), sites.size(), x, y);
                return PauliString<Coeff>(x, y, coeff);
        }

        /// Whether the string still fits the sparse representation (see
        /// pauliengine_detail::sparse_stays); otherwise to_dense() is the better form.
        bool fits_sparse() const {
                const double span = sites.empty() ? 0.0 : pauliengine_detail::site_qubit(sites.back()) + 1.0;
                return pauliengine_detail::sparse_stays(span, static_cast<double>(sites.size()));
        }

        /// Number of non-identity sites.
        size_t size() const { return sites.size(); }

        /// Sorted qubits with a non-identity operator.
        std::vector<int> qubits() const {
                std::vector<int> q;
                q.reserve(sites.size());
                for (Site s : sites) q.push_back(static_cast<int>(pauliengine_detail::site_qubit(s)));
                return q;
        }

        /// Pauli operator ("I", "X", "Y" or "Z") on qubit q.
        std::string get_pauli_at_index(int q) const {
                const auto it = std::lower_bound(sites.begin(), sites.end(),
                        pauliengine_detail::make_site(static_cast<uint32_t>(std::max(q, 0)), 0));
                if (q < 0 || it == sites.end() || static_cast<int>(pauliengine_detail::site_qubit(*it)) != q) return "I";
                return std::string(1, pauliengine_detail::pauli_letter(pauliengine_detail::site_pauli(*it)));
        }

        int count_y() const {
                return static_cast<int>(std::count_if(sites.begin(), sites.end(),
                        [](Site s) { return pauliengine_detail::site_pauli(s) == 2; }));
        }

        bool is_all_z() const {
                return std::all_of(sites.begin(), sites.end(), [](Site s) { return pauliengine_detail::site_pauli(s) == 3; });
        }

        /// (coefficient, {qubit: "X"/"Y"/"Z"}).
        std::pair<Coeff, std::unordered_map<int, std::string>> to_dictionary() const {
                std::unordered_map<int, std::string> ops;
                ops.reserve(sites.size());
                for (Site s : sites) {
                        ops[static_cast<int>(pauliengine_detail::site_qubit(s))] =
                                std::string(1, pauliengine_detail::pauli_letter(pauliengine_detail::site_pauli(s)));
                }
                return {coeff, std::move(ops)};
        }

        /// OpenFermion-style key: [(Pauli, qubit), ...] by ascending qubit.
        std::vector<std::pair<char, int>> key_openfermion() const {
                std::vector<std::pair<char, int>> key;
                key.reserve(sites.size());
                for (Site s : sites) {
                        key.emplace_back(pauliengine_detail::pauli_letter(pauliengine_detail::site_pauli(s)),
                                         static_cast<int>(pauliengine_detail::site_qubit(s)));
                }
                return key;
        }

        std::string to_string() const {
                std::ostringstream oss;
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        oss << std::to_string(coeff.real()) << " + " << std::to_string(coeff.imag()) << "i : ";
                } else {
                        oss << coeff.str() << " : ";
                }
                for (Site s : sites) {
                        oss << pauliengine_detail::pauli_letter(pauliengine_detail::site_pauli(s)) << "("
                            << pauliengine_detail::site_qubit(s) << ")";
                }
                return oss.str();
        }

        /// Hash of the operator part only, consistent with equals() and operator==.
        size_t hash() const { return static_cast<size_t>(pauliengine_detail::hash_sites(sites.data(), sites.size())); }

        /// Whether the operators are equal, ignoring the coefficients.
        bool equals(const SparsePauliString& other) const { return sites == other.sites; }

        bool operator==(const SparsePauliString& other) const {
                return sites == other.sites && pauliengine_detail::coeff_equal_to(coeff, other.coeff);
        }

        /// Product with the phase of the Pauli algebra.
        SparsePauliString operator*(const SparsePauliString& other) const {
                SparsePauliString out;
                out.sites.resize(sites.size() + other.sites.size());
                size_t n = 0;
                const int k = pauliengine_detail::multiply_sites(sites.data(), sites.size(),
                        other.sites.data(), other.sites.size(), out.sites.data(), n);
                out.sites.resize(n);
                out.coeff = pauliengine_detail::phase_product(coeff, other.coeff, k, 1.0);
                return out;
        }

        SparsePauliString operator*(std::complex<double> scale) const { return scaled(Coeff(scale)); }
        SparsePauliString operator*(double scale) const { return scaled(Coeff(scale)); }
        SparsePauliString operator*(int scale) const { return scaled(Coeff(static_cast<double>(scale))); }

        SparsePauliString operator-() const { return scaled(Coeff(-1.0)); }

        /// Whether the strings anti-commute.
        bool anticommutes(const SparsePauliString& other) const {
                return pauliengine_detail::anticommute_sites(sites.data(), sites.size(), other.sites.data(), other.sites.size());
        }

        /// [this, other]: 2 * this * other if they anti-commute, else the zero string
        /// (identity with coefficient 0).
        SparsePauliString commutator(const SparsePauliString& other) const {
                if (!anticommutes(other)) return SparsePauliString(pauliengine_detail::coeff_zero<Coeff>());
                SparsePauliString out = *this * other;
                out.coeff = out.coeff * Coeff(2.0);
                return out;
        }

        /// Coefficient of target's operator in [this, other].
        Coeff targeted_commutator(const SparsePauliString& other, const SparsePauliString& target) const {
                std::vector<Site> product(sites.size() + other.sites.size());
                size_t n = 0;
                const int k = pauliengine_detail::multiply_sites(sites.data(), sites.size(),
                        other.sites.data(), other.sites.size(), product.data(), n);
                if (!(k & 1) || n != target.sites.size() || !std::equal(product.begin(), product.begin() + n, target.sites.begin())) {
                        return pauliengine_detail::coeff_zero<Coeff>();
                }
                return pauliengine_detail::phase_product(coeff, other.coeff, k, 2.0);
        }

        /// Copy with coefficient 1.
        SparsePauliString naked() const {
                SparsePauliString out = *this;
                out.coeff = pauliengine_detail::coeff_one<Coeff>();
                return out;
        }

        /// Copy with qubit indices remapped by {old: new}; an occupied target that is
        /// not itself remapped is swapped (as PauliString::map_qubits).
        SparsePauliString map_qubits(const std::unordered_map<int, int>& qubit_map) const {
                return SparsePauliString(to_dense().map_qubits(qubit_map));
        }

        /// Derivative of the coefficient by the symbol `name` (symbolic coefficients).
        SparsePauliString diff(const std::string& name) const {
                SparsePauliString out = *this;
                out.coeff = pauliengine_detail::coeff_diff(coeff, name);
                return out;
        }

        /// Values or polynomials substituted for symbols of the coefficient.
        template <typename Values>
        SparsePauliString substitute(const Values& values) const {
                SparsePauliString out = *this;
                out.coeff = pauliengine_detail::coeff_subs(coeff, values);
                return out;
        }

    private:
        SparsePauliString scaled(const Coeff& s) const {
                SparsePauliString out = *this;
                out.coeff = out.coeff * s;
                return out;
        }
};
