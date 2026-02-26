#pragma once

#include <unordered_map>
#include <complex>
#include <vector>
#include <iostream>
#include <bit>
#include <functional>
#include <cstdint>
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <symengine/complex.h>
#include <symengine/visitor.h>

using namespace SymEngine;
using SymEngine::Expression;

#define BITS_IN_INTEGER (sizeof(uint64_t) * 8)

const std::complex<double> Unit_matrix(0.0 , 1.0);

inline int popcount(uint64_t x);
void testSymengine();


template<typename Coeff>
class PauliString {
        public:
                using Pauli_structure = std::pair<Coeff, std::unordered_map<int, std::string>>;
                using Hamiltonian_structure = std::vector<Pauli_structure>;

                std::vector<uint64_t> x;
                std::vector<uint64_t> y;
                Coeff coeff{0.0};
                uint64_t is_zero{true};
                PauliString() = default;

                PauliString(const std::vector<uint64_t>& x, const std::vector<uint64_t>& y, Coeff coeff) {
                        this->x = x;
                        this->y = y;
                        this->coeff = coeff;
                }

                PauliString(Coeff coeff, const std::unordered_map<int, std::string>& data) {
                        this->coeff = coeff;
                        uint64_t mask;
                        for (const auto& [key, value] : data) {
                                get_symplectic_form(key, mask, value);
                        }
                }

                PauliString(const std::pair<Coeff, std::vector<std::pair<char, int>>>& input) {
                        this->coeff = input.first;
                        std::vector<std::pair<char, int>> paulis = input.second;
                        uint64_t mask = 0;
                        for (int i = 0; i < paulis.size(); i++) {
                                get_symplectic_form(paulis[i].second, mask, std::string(1, paulis[i].first));

                        }
                }

                /**
                 * @brief Construct a new Pauli String object
                 *
                 * @param pauli_string : A string representation of the Pauli string.
                 * @param coeff : The coefficient associated with the Pauli string.
                 */
                PauliString(const std::string &pauli_string, Coeff &coeff){
                        this->coeff = coeff;
                        uint64_t mask;
                        for (size_t i = 0; i < pauli_string.size(); i++) {
                                get_symplectic_form(i, mask, pauli_string[i]);
                        }
                }


                uint64_t operator==(const PauliString& other) const {
                        return x == other.x && y == other.y;
                }


                PauliString &operator*=(const PauliString& other) {
                        // Source: https://arxiv.org/pdf/2103.02202 figure 12

                        uint64_t *x1 = x.data();
                        uint64_t *y1 = y.data();
                        const uint64_t *x2 = other.x.data();
                        const uint64_t *y2 = other.y.data();
                        size_t n1 = x.size();
                        size_t n2 = other.x.size();
                        size_t n = std::min(n1, n2);

                        // The 1s and 2s bits of 64 mod4 counters.
                        uint64_t c1 = 0;
                        uint64_t c2 = 0;
                        // Iterate over data in 64 bit chunks.
                        for (size_t k = 0; k < n; k++) {
                            uint64_t old_x1 = x1[k];
                            uint64_t old_y1 = y1[k];
                            // Update the left hand side Paulis.
                            x1[k] ^= x2[k];
                            y1[k] ^= y2[k];
                            // Accumulate anti-commutation counts.
                            uint64_t x1y2 = old_x1 & y2[k];
                            uint64_t anti_commutes = (x2[k] & old_y1) ^ x1y2;
                            c2 ^= (c1 ^ x1[k] ^ y1[k] ^ x1y2) & anti_commutes;
                            c1 ^= anti_commutes;
                        }

                        // Update coefficient accounting for factors of i gained from Pauli multiplications.
                        coeff *= other.coeff;
                        int power_of_i = 2 * std::popcount(c2) - std::popcount(c1);
                        if (power_of_i & 1) {
                            coeff *= std::complex<double>(0.0, 1.0);
                        }
                        if (power_of_i & 2) {
                            coeff = -coeff;
                        }

                        // Automatically extend length if needed.
                        if (n1 < n2) {
                            x.insert(x.end(), x2 + n1, x2 + n2);
                            y.insert(y.end(), y2 + n1, y2 + n2);
                        }

                        // DIDNTDO: update is_zero (because the old code didn't?)

                        return *this;
                }
                PauliString operator*(const PauliString& other) const {
                        PauliString copy = *this;
                        copy *= other;
                        return copy;
                }

                PauliString operator*(const std::complex<double> scalar){
                        return PauliString(this->x, this->y, this->coeff * scalar);
                }


                PauliString& operator*=(const std::complex<double> scalar)  {
                        this->coeff*= scalar;
                        return *this;
                }

                uint64_t is_all_z() {
                        for (size_t i = 0; i < this->x.size(); i++) {
                                if (this->x[i] != this->y[i]) {
                                        return false;
                                }
                        }
                        return true;
                }

                Pauli_structure to_dictionary() {
                        std::unordered_map<int, std::string> temp;
                        temp.reserve(this->x.size() * BITS_IN_INTEGER / 2);
                        for (size_t i = 0; i < this->x.size(); ++i) {
                                uint64_t x_word = this->x[i];
                                uint64_t y_word = this->y[i];
                                for (uint64_t j = 0; j < BITS_IN_INTEGER; ++j) {
                                        uint64_t bit = (uint64_t)1 << j;
                                        uint64_t x_bit = x_word & bit;
                                        uint64_t y_bit = y_word & bit;
                                        if (!x_bit && !y_bit) continue;
                                        uint64_t index = i * BITS_IN_INTEGER + j;
                                        if (x_bit && y_bit) {
                                                temp[(int)index] = "Z";
                                        } else if (x_bit) {
                                                temp[(int)index] = "X";
                                        } else {
                                                temp[(int)index] = "Y";
                                        }
                                }
                        }
                        return {this->coeff, temp};
                }

                PauliString map_qubits(std::unordered_map<int, int> qubit_map) {
                        Pauli_structure result_parsed = this->to_dictionary();
                        std::unordered_map<int, std::string> mapped;
                        for (const auto& [old_idx, pauli_char] : result_parsed.second) {
                                auto it = qubit_map.find(old_idx);
                                if (it != qubit_map.end()) {
                                        int new_idx = it->second;
                                        mapped[new_idx] = pauli_char;
                                } else {
                                        mapped[old_idx] = pauli_char;
                                }
                        }
                        return PauliString(this->coeff, mapped);
                }

                std::string to_string() const {
                        std::ostringstream oss;
                        size_t num_qubits = x.size() * BITS_IN_INTEGER;

                        // if Coeff is complex
                        if constexpr (std::is_same<Coeff, std::complex<double>>::value) {
                                oss << std::to_string(this->coeff.real()) + " + " + std::to_string(this->coeff.imag()) + "i : ";
                        }
                        // if is SymEngine::Expression this should be enough
                        else {
                                oss << str(this->coeff) << " : ";
                        }

                        for (size_t i = 0; i < num_qubits; ++i) {
                                size_t word = i / BITS_IN_INTEGER;
                                size_t bit = i % BITS_IN_INTEGER;
                                uint64_t xi = (x[word] >> bit) & 1;
                                uint64_t yi = (y[word] >> bit) & 1;
                                char p;
                                if (xi && yi)      p = 'Z';
                                else if (xi)       p = 'X';
                                else if (yi)       p = 'Y';
                                else               continue;

                                oss << p << "(" << i << ")";
                        }
                        return oss.str();
                }

                PauliString commutator(const PauliString &other) const{
                        size_t max_length = std::max(this->x.size(), other.x.size());
                        size_t min_length = std::min(this->x.size(), other.x.size());
                        std::vector<uint64_t> new_x(max_length);
                        std::vector<uint64_t> new_y(max_length);
                        Coeff coeff_temp = this->coeff * other.coeff;

                        for (int i = 0; i < min_length; i++) {
                                new_x[i] = this->x[i] ^ other.x[i];
                                new_y[i] = this->y[i] ^ other.y[i];
                        }

                        const std::vector<uint64_t>& bigger_x = (this->x.size() >= other.x.size()) ? this->x : other.x;
                        const std::vector<uint64_t>& bigger_y = (this->y.size() >= other.y.size()) ? this->y : other.y;

                        for (size_t i = min_length; i <max_length; i++){

                                new_x[i]= bigger_x[i];
                                new_y[i]= bigger_y[i];
                        }
                        int first_fac = 0;
                        int second_fac = 0;
                        for (int i = 0; i < min_length; i++) {
                                first_fac += std::popcount((~this->x[i] & other.x[i] & this->y[i] & other.y[i]) ^ (this->x[i] & ~other.x[i] & ~this->y[i] & other.y[i]) ^ (this->x[i] & other.x[i] & this->y[i] & ~other.y[i]));
                                second_fac += std::popcount((~this->x[i] & other.x[i] & this->y[i] & ~other.y[i]) ^ (this->x[i] & ~other.x[i] & this->y[i] & other.y[i]) ^ (this->x[i] & other.x[i] & ~this->y[i] & other.y[i]));
                        }
                        int tau = first_fac - second_fac;
                        std::complex<double> coeff_new{1.0, 0.0};
                        while (tau < 0) {
                                tau += 4;
                        }

                        if (tau % 2 == 0) {
                                return PauliString(new_x, new_y, 0);
                        }
                        for (int i = 0; i < tau % 4; i++) {
                                coeff_new *= Unit_matrix;
                        }
                        return PauliString(new_x, new_y,  coeff_temp * (coeff_new * 2.0));

                }

                Coeff get_coeff() const {
                                return this->coeff;
                }

                void set_coeff(Coeff new_coeff) {
                        this->coeff = new_coeff;
                }

                PauliString trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state) {
                        double factor = 1.0;
                        std::vector<uint64_t> x_y(this->x.size());
                        for (size_t i = 0; i < this->x.size() ; i++) {
                                x_y[i] = this->x[i] ^ this->y[i];
                        }
                        uint64_t mask;
                        for (size_t i = 0; i < qubits.size(); i++) {
                                mask = (uint64_t) 1 << (qubits[i] % BITS_IN_INTEGER);
                                if ((x_y[qubits[i] / BITS_IN_INTEGER] & mask) == mask) {
                                        return PauliString();
                                }
                        }
                        std::vector<uint64_t> z(this->x.size());
                        for (size_t i = 0; i < this->x.size() ; i++) {
                                z[i] = this->x[i] & this->y[i];
                        }
                        for (size_t i = 0; i < qubits.size(); i++) {
                                mask = (uint64_t) 1 << (qubits[i] % BITS_IN_INTEGER);
                                if ((z[qubits[i] / BITS_IN_INTEGER] & mask) == mask) {
                                        if (state[i] == 1) {
                                                factor *= -1.0;
                                        }
                                        this->x[qubits[i] / BITS_IN_INTEGER] -= mask;
                                        this->y[qubits[i] / BITS_IN_INTEGER] -= mask;
                                }
                        }
                        return PauliString(this->x, this->y, this->coeff * factor);
                }

                std::vector<int> qubits() {
                        std::vector<int> qubit_list(this->x.size() * BITS_IN_INTEGER);
                        for (uint64_t i = 0; i < this->x.size(); i++) {
                                uint64_t mask = 1;
                                for (size_t bit = 0; bit < BITS_IN_INTEGER; bit++) {
                                        if (((mask & this->x[i]) == mask) |((mask & this->y[i]) == mask)) {
                                                qubit_list.push_back((int)(i * BITS_IN_INTEGER + bit));
                                        }
                                }
                        }
                        return qubit_list;
                }

                std::vector<std::pair<char, int>> key_openfermion() {
                        std::vector<std::pair<char, int>> result(this->x.size() * BITS_IN_INTEGER);
                        uint64_t mask;
                        for (int i = 0; i < this->x.size(); i++) {
                                mask = 1;
                                for (int j = 0; j < BITS_IN_INTEGER; j++) {
                                        uint64_t x_bit = (mask & this->x[i]) == mask;
                                        uint64_t y_bit = (mask & this->y[i]) == mask;
                                        if (x_bit && y_bit) {
                                                result.push_back({'Z', i * BITS_IN_INTEGER + j});
                                        } else if (x_bit) {
                                                result.push_back({'X', i * BITS_IN_INTEGER + j});
                                        } else if (y_bit) {
                                                result.push_back({'Y', i * BITS_IN_INTEGER + j});
                                        }
                                        mask = mask << 1;
                                }
                        }
                        return result;
                }

                std::string get_pauli_from_index(int index) {
                        uint64_t position = index / BITS_IN_INTEGER;
                        uint64_t mask = (uint64_t) 1 << (index % BITS_IN_INTEGER);
                        uint64_t x_bit = (this->x[index] & mask) == mask;
                        uint64_t y_bit = (this->y[index] & mask) == mask;
                        if (x_bit && y_bit) {
                                return "Z";
                        } else if(x_bit) {
                                return "X";
                        } else if(y_bit) {
                                return "Y";
                        }
                        return "I";
                }

                PauliString diff(const std::string& symbol_name) const{
                        set_basic symbols = get_free_symbols(this->coeff);
                        Coeff coeff_diff;
                        for(const auto& s : symbols) {
                                if (is_a<Symbol>(*s)) {
                                        const Symbol& sym = down_cast<const Symbol&>(*s);
                                        if (sym.get_name() == symbol_name) {
                                                coeff_diff = this->coeff.diff(symbol(symbol_name));
                                                return PauliString(this->x, this->y, coeff_diff);
                                        }
                                }
                        }
                        return PauliString();

                }

                PauliString substitute(const std::unordered_map<std::string, std::complex<double>>& substitution_map) const {
                        map_basic_basic substitutions;
                        for (const auto& [symbol_string, value] : substitution_map) {
                                Expression temp = value;
                                RCP<const Basic> basic_ptr = temp.get_basic();
                                substitutions[symbol(symbol_string)] = basic_ptr;
                        }
                        Coeff temp = this->coeff.subs(substitutions);
                        return PauliString(x, y, temp);
                }

                PauliString copy() const{
                        return PauliString(this->x, this->y, this->coeff);
                }

                bool equals(const PauliString& other) {
                        return this->x == other.x && this->y == other.y && PauliString::to_complex(this->coeff) == PauliString::to_complex(other.coeff);
                }

                static std::complex<double> to_complex(const std::complex<double>& expr) {
                        return expr;
                }

                static std::complex<double> to_complex(const Expression &expr) {
                        const auto &basic = *expr.get_basic();

                        if (is_a<RealDouble>(basic)) {
                                const auto &rd = down_cast<const RealDouble &>(basic);
                                return { rd.as_double(), 0.0 };
                        }
                        else if (is_a<ComplexDouble>(basic)) {
                                const auto &cd = down_cast<const ComplexDouble &>(basic);
                                return {
                                eval_double(*cd.real_part()),
                                eval_double(*cd.imaginary_part())
                                };
                        }
                        else {
                                try {
                                double val = eval_double(basic);
                                return { val, 0.0 };
                                } catch (...) {
                                throw std::runtime_error("Expression cannot be converted to std::complex<double>.");
                                }
                        }
                }

        private:
                static set_basic get_free_symbols(const Expression& expr) {
                        vec_basic args = expr.get_basic()->get_args();
                        set_basic symbols;
                        for (const auto& a : args) {
                                if (is_a<Symbol>(*a)) {
                                        symbols.insert(a);
                                }
                        }
                        return symbols;
                }

                void get_symplectic_form(size_t pauli_index, uint64_t& mask, const std::string& pauli_char) {
                        size_t index = pauli_index / BITS_IN_INTEGER;
                        if ((index) + 1 > x.size()) {
                                x.push_back(0);
                                y.push_back(0);
                        }
                        mask = ((uint64_t) 1 <<  pauli_index % BITS_IN_INTEGER);


                        if (pauli_char == "X" || pauli_char == "Z") {
                                this->x[index] |= mask;
                        }
                        if (pauli_char == "Y" || pauli_char == "Z") {
                                this->y[index] |= mask;
                        }
                }

};


template<typename Coeff>
struct PauliStringHash {
    size_t operator()(const PauliString<Coeff>& ps) const {
        size_t seed = 0;
        for (auto& val : ps.x) {
            seed ^= std::hash<uint64_t>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        for (auto& val : ps.y) {
            seed ^= std::hash<uint64_t>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
