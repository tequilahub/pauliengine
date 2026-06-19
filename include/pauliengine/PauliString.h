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

const std::complex<double> I(0.0 , 1.0);

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
                        size_t n = std::min(x.size(), y.size());
                        while (n > 0 && x[n - 1] == 0 && y[n - 1] == 0) {
                                --n;
                        }
                        this->x.assign(x.begin(), x.begin() + n);
                        this->y.assign(y.begin(), y.begin() + n);
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
                PauliString(Coeff &coeff, const std::string &pauli_string){
                        this->coeff = coeff;
                        uint64_t mask;
                        for (size_t i = 0; i < pauli_string.size(); i++) {
                                get_symplectic_form(i, mask, std::string(1, pauli_string[i]));
                        }
                }


                uint64_t operator==(const PauliString& other) const {
                        if (!((this->x == other.x) && (this->y == other.y))) {
                                return false;
                        }
                        if constexpr (std::is_same<Coeff, std::complex<double>>::value) {
                                return this->coeff == other.coeff;
                        } else {
                                return (this->coeff - other.coeff) == 0;

                        }
                }

                PauliString operator*(const PauliString& other) const {
                        // Source: https://arxiv.org/pdf/2103.02202 figure 12

                        std::vector<uint64_t> x1 = x;
                        std::vector<uint64_t> y1 = y;
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
                        auto new_coeff = this->coeff * other.coeff;
                        int power_of_i = 2 * std::popcount(c2) - std::popcount(c1);
                        if (power_of_i & 1) {
                            new_coeff *= std::complex<double>(0.0, 1.0);
                        }
                        if (power_of_i & 2) {
                            new_coeff = -new_coeff;
                        }

                        if (n1 < n2) {
                            x1.insert(x1.end(), x2 + n1, x2 + n2);
                            y1.insert(y1.end(), y2 + n1, y2 + n2);
                        }

                        return PauliString(x1, y1, new_coeff);
                }

                PauliString &operator*=(const PauliString& other) {
                        *this = *this * other;
                        return *this;
                }

                PauliString operator*(const std::complex<double> scalar) const {
                        return PauliString(this->x, this->y, this->coeff * scalar);
                }

                PauliString operator*(const int scalar) const {
                        return PauliString(this->x, this->y, this->coeff * static_cast<double>(scalar));
                }

                PauliString operator*(const double scalar) const {
                        return PauliString(this->x, this->y, this->coeff * scalar);
                }



                PauliString& operator*=(const std::complex<double> scalar)  {
                        this->coeff*= scalar;
                        return *this;
                }

                uint64_t is_all_z() const {
                        for (size_t i = 0; i < this->x.size(); i++) {
                                if (this->x[i] != this->y[i]) {
                                        return false;
                                }
                        }
                        return true;
                }

                // Number of qubits with non-identity operator.
                size_t size() const {
                        size_t count = 0;
                        for (size_t i = 0; i < this->x.size(); i++) {
                                count += std::popcount(this->x[i] | this->y[i]);
                        }
                        return count;
                }

                // Number of Y operators (x bit = 0, y bit = 1).
                int count_y() const {
                        int y_count = 0;
                        for (size_t i = 0; i < this->x.size(); i++) {
                                y_count += std::popcount(~this->x[i] & this->y[i]);
                        }
                        return y_count;
                }

                // Return a copy with coefficient set to 1.
                PauliString naked() const {
                        Coeff one;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                one = std::complex<double>(1.0, 0.0);
                        } else {
                                one = Coeff(1);
                        }
                        return PauliString(this->x, this->y, one);
                }

                Pauli_structure to_dictionary() const {
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

                PauliString map_qubits(std::unordered_map<int, int> qubit_map) const {
                        const auto& current_ops = this->to_dictionary().second;

                        // For each explicit A→B: if B is occupied and B is not itself a
                        // source in qubit_map, implicitly add B→A (swap instead of overwrite).
                        std::unordered_map<int, int> effective_map = qubit_map;
                        for (const auto& [from, to] : qubit_map) {
                                if (current_ops.count(to) > 0 && qubit_map.count(to) == 0) {
                                        effective_map[to] = from;
                                }
                        }

                        std::unordered_map<int, std::string> mapped;
                        for (const auto& [old_idx, pauli_char] : current_ops) {
                                auto it = effective_map.find(old_idx);
                                if (it != effective_map.end()) {
                                        mapped[it->second] = pauli_char;
                                } else {
                                        mapped[old_idx] = pauli_char;
                                }
                        }
                        return PauliString(this->coeff, mapped);
                }

                std::string to_string() const {
                        std::ostringstream oss;
                        size_t num_qubits = x.size() * BITS_IN_INTEGER;

                        if constexpr (std::is_same<Coeff, std::complex<double>>::value) {
                                oss << std::to_string(this->coeff.real()) + " + " + std::to_string(this->coeff.imag()) + "i : ";
                        } else {
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
                                coeff_new *= std::complex<double>(0.0 , 1.0);
                        }
                        return PauliString(new_x, new_y,  coeff_temp * (coeff_new * 2.0));

                }

                Coeff get_coeff() const {
                                
                        return this->coeff;
                }

                void set_coeff(Coeff new_coeff) {
                        this->coeff = new_coeff;
                }

                // Project this Pauli string onto a computational-basis state
                // for the listed qubits. state[i] = 0 means qubit `qubits[i]` is
                // |0>, state[i] = 1 means |1>. Returns a NEW PauliString with the
                // traced qubits removed;
                //
                // Rules per traced qubit q:
                //   - I:  <s|I|s> = 1                 → keep as-is
                //   - X:  <0|X|0> = <1|X|1> = 0       → result is 0
                //   - Y:  <0|Y|0> = <1|Y|1> = 0       → result is 0
                //   - Z:  <0|Z|0> = +1, <1|Z|1> = -1  → drop the Z, sign-flip if state=1
                PauliString trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state) const {
                        if (qubits.size() != state.size()) {
                                throw std::runtime_error("trace_out_qubits: qubits and state must have the same length");
                        }

                        std::vector<uint64_t> new_x = this->x;
                        std::vector<uint64_t> new_y = this->y;
                        double sign = 1.0;

                        for (size_t i = 0; i < qubits.size(); i++) {
                                size_t word = qubits[i] / BITS_IN_INTEGER;
                                if (word >= new_x.size()) {
                                        // Qubit beyond stored range = identity, nothing to do.
                                        continue;
                                }
                                uint64_t mask = (uint64_t)1 << (qubits[i] % BITS_IN_INTEGER);
                                bool x_bit = (new_x[word] & mask) != 0;
                                bool y_bit = (new_y[word] & mask) != 0;

                                if (x_bit != y_bit) {
                                        // X (x=1,y=0) or Y (x=0,y=1) → expectation value is 0.
                                        return PauliString();
                                }
                                if (x_bit && y_bit) {
                                        // Z on this qubit. Sign-flip if |1>, then remove the operator.
                                        if (state[i] != 0) {
                                                sign *= -1.0;
                                        }
                                        new_x[word] &= ~mask;
                                        new_y[word] &= ~mask;
                                }
                                // else: identity → nothing to do.
                        }

                        return PauliString(new_x, new_y, this->coeff * sign);
                }

                // Project this Pauli string onto an arbitrary single-qubit state
                // |psi_i> = a_i |0> + b_i |1> for each traced qubit. 
                PauliString trace_out_qubits(
                        const std::vector<int>& qubits,
                        const std::vector<std::pair<std::complex<double>, std::complex<double>>>& states
                ) const {
                        if (qubits.size() != states.size()) {
                                throw std::runtime_error("trace_out_qubits: qubits and states must have the same length");
                        }

                        std::vector<uint64_t> new_x = this->x;
                        std::vector<uint64_t> new_y = this->y;
                        std::complex<double> factor(1.0, 0.0);
                        const std::complex<double> I_unit(0.0, 1.0);

                        for (size_t i = 0; i < qubits.size(); i++) {
                                size_t word = qubits[i] / BITS_IN_INTEGER;
                                bool x_bit = false, y_bit = false;
                                if (word < new_x.size()) {
                                        uint64_t mask = (uint64_t)1 << (qubits[i] % BITS_IN_INTEGER);
                                        x_bit = (new_x[word] & mask) != 0;
                                        y_bit = (new_y[word] & mask) != 0;
                                        new_x[word] &= ~mask;
                                        new_y[word] &= ~mask;
                                }

                                const std::complex<double>& a = states[i].first;
                                const std::complex<double>& b = states[i].second;
                                std::complex<double> ac_b = std::conj(a) * b;
                                std::complex<double> a_bc = a * std::conj(b);

                                std::complex<double> expectation;
                                if (!x_bit && !y_bit) {
                                        // I: |a|^2 + |b|^2
                                        expectation = std::norm(a) + std::norm(b);
                                } else if (x_bit && !y_bit) {
                                        // X: a* b + a b*
                                        expectation = ac_b + a_bc;
                                } else if (!x_bit && y_bit) {
                                        // Y: -i a* b + i a b*
                                        expectation = -I_unit * ac_b + I_unit * a_bc;
                                } else {
                                        // Z: |a|^2 - |b|^2
                                        expectation = std::norm(a) - std::norm(b);
                                }

                                factor *= expectation;
                        }

                        return PauliString(new_x, new_y, this->coeff * factor);
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
                        if (index < 0) {
                                return "I";
                        }
                        uint64_t position = index / BITS_IN_INTEGER;
                        if (position >= this->x.size() || position >= this->y.size()) {
                                return "I";
                        }
                        uint64_t mask = (uint64_t) 1 << (index % BITS_IN_INTEGER);
                        uint64_t x_bit = (this->x[position] & mask) == mask;
                        uint64_t y_bit = (this->y[position] & mask) == mask;
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
                        return this->x == other.x && this->y == other.y;
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
                        if (x.size() <= index) {
                                x.resize(index + 1, 0);
                                y.resize(index + 1, 0);
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

template<typename T, typename U, typename V = decltype(T() * U())>
inline auto operator*(const PauliString<T>& left_op, const PauliString<U>& right_op) {
    return PauliString<V>(left_op.x, left_op.y, static_cast<V>(left_op.coeff)) * PauliString<V>(right_op.x, right_op.y, static_cast<V>(right_op.coeff));
}

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

// Equality comparator that only compares the operator part (x and y vectors),
// ignoring the coefficient.  Required so that hash maps keyed by PauliString
// treat two terms with the same operator but different coefficients as the
// same key — which is what compact() and QH::operator== need.
template<typename Coeff>
struct PauliStringOperatorEqual {
    bool operator()(const PauliString<Coeff>& a, const PauliString<Coeff>& b) const {
        return a.x == b.x && a.y == b.y;
    }
};
