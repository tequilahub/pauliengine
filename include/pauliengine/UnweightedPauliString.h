#pragma once

#include <algorithm>
#include <bit>
#include <complex>
#include <cstdint>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pauliengine/PauliString.h"

/// A Pauli string without a coefficient: i^phase * P with P a tensor product of
/// I, X, Y and Z and phase in {0, 1, 2, 3}. Use it only as a basis element of a Lie
/// algebra; linear combinations need PauliString and QubitHamiltonian.
///
/// The product of two Pauli operators is again a Pauli operator up to a power of
/// i, so these strings are closed under multiplication: the phase is an exact
/// integer modulo 4 and no floating-point number is involved. P itself is
/// Hermitian, so phase 0 means +P and not a phase relative to an XZ ordering.
class unweighted_pauliString {
        public:
                WordVec x;
                WordVec y;

                /// The identity with phase 0.
                unweighted_pauliString() = default;

                /// From a map {qubit: "X"/"Y"/"Z"/"I"}.
                explicit unweighted_pauliString(const std::unordered_map<int, std::string>& paulis, int phase = 0)
                    : phase_(normalize_phase(phase)) {
                        for (const auto& [qubit, pauli] : paulis) {
                                if (pauli.size() != 1) {
                                        throw std::invalid_argument("invalid Pauli operator '" + pauli + "'");
                                }
                                set_pauli_bits(x, y, qubit, pauli[0]);
                        }
                        trim();
                }

                /// From a list of (Pauli, qubit) pairs.
                explicit unweighted_pauliString(const std::vector<std::pair<char, int>>& paulis, int phase = 0)
                    : phase_(normalize_phase(phase)) {
                        for (const auto& [pauli, qubit] : paulis) {
                                set_pauli_bits(x, y, qubit, pauli);
                        }
                        trim();
                }

                /// From a dense text such as "XIZ" or a sparse one such as "X0 Z2"; empty is the identity.
                explicit unweighted_pauliString(const std::string& paulis, int phase = 0)
                    : phase_(normalize_phase(phase)) {
                        parse_pauli_text(paulis, x, y);
                        trim();
                }

                /// The operator part of a Pauli string; its coefficient is dropped.
                template<typename Coeff>
                explicit unweighted_pauliString(const PauliString<Coeff>& string, int phase = 0)
                    : x(string.x), y(string.y), phase_(normalize_phase(phase)) {}

                /// From symplectic words; trailing zero words are dropped.
                unweighted_pauliString(const uint64_t* x_words, const uint64_t* y_words, size_t n_words, int phase)
                    : phase_(normalize_phase(phase)) {
                        x.assign(x_words, x_words + n_words);
                        y.assign(y_words, y_words + n_words);
                        trim();
                }

                /// Exponent k of the phase i^k, in {0, 1, 2, 3}.
                int phase() const noexcept { return phase_; }

                /// The phase i^k as a complex number, with exact real and imaginary parts.
                std::complex<double> phase_factor() const noexcept {
                        switch (phase_) {
                                case 0: return {1.0, 0.0};
                                case 1: return {0.0, 1.0};
                                case 2: return {-1.0, 0.0};
                                default: return {0.0, -1.0};
                        }
                }

                /// Product including the phase of the Pauli algebra.
                unweighted_pauliString operator*(const unweighted_pauliString& other) const {
                        const size_t n_max = std::max(x.size(), other.x.size());
                        unweighted_pauliString result;
                        result.x.resize(n_max);
                        result.y.resize(n_max);
                        const int power_of_i = multiply_symplectic_words(
                                x.data(), y.data(), x.size(),
                                other.x.data(), other.y.data(), other.x.size(),
                                result.x.data(), result.y.data());
                        result.phase_ = normalize_phase(phase_ + other.phase_ + power_of_i);
                        result.trim();
                        return result;
                }

                unweighted_pauliString& operator*=(const unweighted_pauliString& other) {
                        *this = *this * other;
                        return *this;
                }

                /// The string multiplied by -1.
                unweighted_pauliString operator-() const {
                        unweighted_pauliString result(*this);
                        result.phase_ = normalize_phase(phase_ + 2);
                        return result;
                }

                /// Equal operator and equal phase.
                bool operator==(const unweighted_pauliString& other) const {
                        return phase_ == other.phase_ && equals(other);
                }

                /// Equal operator, ignoring the phase.
                bool equals(const unweighted_pauliString& other) const {
                        return x == other.x && y == other.y;
                }

                /// True if the two strings commute, false if they anti-commute.
                bool commutes(const unweighted_pauliString& other) const {
                        const size_t n_min = std::min(x.size(), other.x.size());
                        uint64_t anti = 0;
                        for (size_t i = 0; i < n_min; ++i) {
                                anti ^= (x[i] & other.y[i]) ^ (other.x[i] & y[i]);
                        }
                        return (std::popcount(anti) & 1) == 0;
                }

                /// The commutator without its factor 2, as used for a Lie algebra basis:
                /// [P, Q] = 2 * P * Q if the strings anti-commute, so this returns P * Q
                /// then and nothing if they commute.
                std::optional<unweighted_pauliString> commutator(const unweighted_pauliString& other) const {
                        if (commutes(other)) {
                                return std::nullopt;
                        }
                        return *this * other;
                }

                /// Hermitian conjugate: the operator stays and the phase is conjugated.
                unweighted_pauliString dagger() const {
                        unweighted_pauliString result(*this);
                        result.phase_ = normalize_phase(-phase_);
                        return result;
                }

                /// Copy with phase 0.
                unweighted_pauliString naked() const {
                        unweighted_pauliString result(*this);
                        result.phase_ = 0;
                        return result;
                }

                /// The Pauli string with coefficient coeff * i^phase.
                template<typename Coeff>
                PauliString<Coeff> to_pauli_string(Coeff coeff) const {
                        if (phase_ & 1) {
                                coeff *= std::complex<double>(0.0, 1.0);
                        }
                        if (phase_ & 2) {
                                coeff = -coeff;
                        }
                        return PauliString<Coeff>(x, y, coeff);
                }

                /// Value hash over the operator part, consistent with operator==; strings
                /// differing only by the phase are allowed to collide.
                size_t hash() const { return operator_part().hash(); }

                /// See PauliString::size().
                size_t size() const { return operator_part().size(); }

                /// See PauliString::count_y().
                int count_y() const { return operator_part().count_y(); }

                /// See PauliString::qubits().
                std::vector<int> qubits() const { return operator_part().qubits(); }

                /// See PauliString::key_openfermion().
                std::vector<std::pair<char, int>> key_openfermion() const {
                        return operator_part().key_openfermion();
                }

                /// See PauliString::get_hash(); the phase is not part of the key.
                std::vector<uint64_t> get_hash() const { return operator_part().get_hash(); }

                /// Inverse of get_hash().
                static unweighted_pauliString from_hash(const std::vector<uint64_t>& limbs, int phase = 0) {
                        return unweighted_pauliString(PauliString<std::complex<double>>::from_hash(limbs, 1.0), phase);
                }

                /// Pauli operator ("I", "X", "Y" or "Z") on the given qubit.
                std::string get_pauli_from_index(int index) const {
                        return operator_part().get_pauli_from_index(index);
                }

                /// The string as a (phase factor, {qubit: "X"/"Y"/"Z"}) pair.
                std::pair<std::complex<double>, std::unordered_map<int, std::string>> to_dictionary() const {
                        return {phase_factor(), operator_part().to_dictionary().second};
                }

                /// Human-readable form, e.g. "+i : X(0)Z(2)".
                std::string to_string() const {
                        static const char* const prefixes[] = {"+1", "+i", "-1", "-i"};
                        std::ostringstream out;
                        out << prefixes[phase_] << " : ";
                        for (const auto& [pauli, qubit] : key_openfermion()) {
                                out << pauli << "(" << qubit << ")";
                        }
                        return out.str();
                }

        private:
                uint8_t phase_{0};

                static uint8_t normalize_phase(int phase) noexcept {
                        return static_cast<uint8_t>(phase & 3);
                }

                PauliString<std::complex<double>> operator_part() const {
                        return PauliString<std::complex<double>>(x, y, std::complex<double>(1.0, 0.0));
                }

                /// Equal operators must have equal word counts, as in PauliString.
                void trim() {
                        size_t n = x.size();
                        while (n > 0 && x[n - 1] == 0 && y[n - 1] == 0) {
                                --n;
                        }
                        x.resize(n);
                        y.resize(n);
                }
};
