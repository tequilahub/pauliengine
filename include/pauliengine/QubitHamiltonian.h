#pragma once
#include <unordered_map>
#include <string>
#include <complex>
#include <vector>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <algorithm>

#include "pauliengine/PauliString.h"

#include <symengine/expression.h>
#include <symengine/functions.h>
using SymEngine::Expression;

#ifdef PAULIENGINE_HAS_OPENMP
#include <omp.h>
#endif

// Threshold for switching on OpenMP. Below this many pair-multiplications the
// thread-spawn overhead dominates.
#ifndef PAULIENGINE_OMP_THRESHOLD
#define PAULIENGINE_OMP_THRESHOLD 4096
#endif


#define BITS_IN_INTEGER (sizeof(uint64_t) * 8)

template<typename Coeff>
using Pauli_structure = std::pair<Coeff, std::unordered_map<int, std::string>>;

template<typename Coeff>
using Hamiltonian_structure = std::vector<Pauli_structure<Coeff>>;

using Matrix2D = std::vector<std::vector<std::complex<double>>>;

const std::complex<double> neg_I(0.0 , -1.0);

// Tag type used by the fast-path constructor below: pass this when the caller
// already guarantees the input vector is compacted (deduplicated + zero-filtered).
// Lets operations like operator* skip a redundant compact() at the end.
struct already_compacted_t {};
inline constexpr already_compacted_t already_compacted{};

template<typename Coeff>
class QubitHamiltonian{
    public:
    std::vector<PauliString<Coeff>> data;

        QubitHamiltonian(const std::vector<PauliString<Coeff>>& data){
                this->data = data;
                this->compact();
        }


        QubitHamiltonian(std::vector<PauliString<Coeff>>&& data){
                this->data = std::move(data);
                this->compact();
        }

        // Fast-path constructor: skips compact(). Only use when you can prove
        // the vector has no duplicates and no zero-coefficient terms.
        QubitHamiltonian(std::vector<PauliString<Coeff>>&& data, already_compacted_t){
                this->data = std::move(data);
        }

        QubitHamiltonian(const Hamiltonian_structure<Coeff>& data){
                std::vector<PauliString<Coeff>> converted_data;
                converted_data.reserve(data.size());
                for (Pauli_structure<Coeff> entry : data) {
                        converted_data.push_back(PauliString<Coeff>(entry.first, entry.second));
                }
                this->data = converted_data;
                this->compact();
        }

        QubitHamiltonian trace_out_qubits(const std::vector<int>& qubits, const std::vector<int>& state) const {
                std::vector<PauliString<Coeff>> reduced;
                reduced.reserve(this->data.size());
                for (const PauliString<Coeff>& current : this->data) {
                        reduced.push_back(current.trace_out_qubits(qubits, state));
                }
                // The constructor runs compact() which drops zero-coefficient terms,
                // including the sentinel default-constructed PauliString returned for
                // tracing X/Y against a computational-basis state.
                return QubitHamiltonian(std::move(reduced));
        }

        // Overload accepting arbitrary single-qubit amplitudes (a, b) for each
        // traced qubit (general superposition states).
        QubitHamiltonian trace_out_qubits(
                const std::vector<int>& qubits,
                const std::vector<std::pair<std::complex<double>, std::complex<double>>>& states
        ) const {
                std::vector<PauliString<Coeff>> reduced;
                reduced.reserve(this->data.size());
                for (const PauliString<Coeff>& current : this->data) {
                        reduced.push_back(current.trace_out_qubits(qubits, states));
                }
                return QubitHamiltonian(std::move(reduced));
        }
        std::string to_string(){
                std::string result = "";
                for (int i = 0; i < this->data.size(); i++) {
                        result.append(data[i].to_string() + "\n");
                }
                return result;
        }


        QubitHamiltonian operator*(std::complex<double> scale) const{
            QubitHamiltonian to_scale = QubitHamiltonian(this->data);
            for (auto &entry : to_scale.data) {
                    entry.coeff *= scale;
                }
                return to_scale;
        }

        
        QubitHamiltonian operator*(int scale) const{
            QubitHamiltonian to_scale = QubitHamiltonian(this->data);
            for (auto &entry : to_scale.data) {
                    entry.coeff *= scale;
                }
                return to_scale;
        }

        QubitHamiltonian operator*(double scale) const{
            QubitHamiltonian to_scale = QubitHamiltonian(this->data);
            for (auto &entry : to_scale.data) {
                    entry.coeff *= scale;
                }
                return to_scale;
        }

        QubitHamiltonian operator*(const QubitHamiltonian& other) const{ 
                const int n1 = static_cast<int>(this->data.size());  // signed for OpenMP 2.0
                const size_t n2 = other.data.size();
                std::vector<PauliString<Coeff>> product(static_cast<size_t>(n1) * n2);

                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
#ifdef PAULIENGINE_HAS_OPENMP
                        #pragma omp parallel for schedule(static) if(static_cast<size_t>(n1) * n2 >= PAULIENGINE_OMP_THRESHOLD)
#endif
                        for (int i = 0; i < n1; ++i) {
                                for (size_t j = 0; j < n2; ++j) {
                                        product[static_cast<size_t>(i) * n2 + j] = this->data[i] * other.data[j];
                                }
                        }
                } else {
                        for (int i = 0; i < n1; ++i) {
                                for (size_t j = 0; j < n2; ++j) {
                                        product[static_cast<size_t>(i) * n2 + j] = this->data[i] * other.data[j];
                                }
                        }
                }
                return QubitHamiltonian(std::move(product));
        }

        QubitHamiltonian operator+(const QubitHamiltonian& other) const {
                std::vector<PauliString<Coeff>> combined;
                combined.reserve(this->data.size() + other.data.size());
                combined.insert(combined.end(), this->data.begin(), this->data.end());
                combined.insert(combined.end(), other.data.begin(), other.data.end());
                return QubitHamiltonian(std::move(combined));  // constructor calls compact()
        }

        QubitHamiltonian operator-(const QubitHamiltonian& other) const {
                return *this + (other * -1.0);
        }

        QubitHamiltonian operator-() const {
                return *this * -1.0;
        }

        bool operator==(const QubitHamiltonian& other) const {
                using Map = std::unordered_map<PauliString<Coeff>, Coeff,
                        PauliStringHash<Coeff>, PauliStringOperatorEqual<Coeff>>;
                Map first_data, second_data;
                for (const auto& ps : this->data) {
                        first_data[ps] = first_data[ps] + ps.coeff;
                }
                for (const auto& ps : other.data) {
                        second_data[ps] = second_data[ps] + ps.coeff;
                }
                return first_data == second_data;
        }

        QubitHamiltonian set_all_coeff(std::complex<double> value) const {
                std::vector<PauliString<Coeff>> temp_data;
                temp_data.reserve(this->data.size());
                Coeff coeff;
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        coeff = value;
                } else {
                        coeff = Expression(value);
                }
                for (const auto& ps : this->data) {
                        temp_data.push_back(PauliString<Coeff>(ps.x, ps.y, coeff));
                }
                return QubitHamiltonian(std::move(temp_data));
        }

        QubitHamiltonian diff(std::string symbol) const{
                std::vector<PauliString<Coeff>> data;
                for (const PauliString<Coeff>& ps : this->data) {
                        PauliString<Coeff> temp = ps.diff(symbol);
                        if (temp.is_zero == false) {
                                data.push_back(temp);
                        }
                }
                return QubitHamiltonian(std::move(data));
        }

        QubitHamiltonian substitute(const std::unordered_map<std::string, std::complex<double>>& substitution_map) const {
                std::vector<PauliString<Coeff>> temp_data;
                temp_data.reserve(this->data.size());
                for (const auto &ps : data) {
                        temp_data.push_back(ps.substitute(substitution_map));
                }
                return QubitHamiltonian(std::move(temp_data));
        }


        QubitHamiltonian compact(){
                if (data.empty()) return *this;

                // Sort-based compact.
                auto op_less = [](const PauliString<Coeff>& a, const PauliString<Coeff>& b) {
                        if (a.x != b.x) return a.x < b.x;
                        return a.y < b.y;
                };
                auto op_equal = [](const PauliString<Coeff>& a, const PauliString<Coeff>& b) {
                        return a.x == b.x && a.y == b.y;
                };

#ifdef PAULIENGINE_HAS_OPENMP
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        if (data.size() >= PAULIENGINE_OMP_THRESHOLD) {
                                parallel_sort_data(op_less);
                        } else {
                                std::sort(data.begin(), data.end(), op_less);
                        }
                } else {
                        std::sort(data.begin(), data.end(), op_less);
                }
#else
                std::sort(data.begin(), data.end(), op_less);
#endif

                size_t write = 0;
                for (size_t read = 1; read < data.size(); ++read) {
                        if (op_equal(data[write], data[read])) {
                                data[write].coeff = data[write].coeff + data[read].coeff;
                        } else {
                                ++write;
                                if (write != read) {
                                        data[write] = std::move(data[read]);
                                }
                        }
                }
                data.resize(write + 1);

                
                data.erase(std::remove_if(data.begin(), data.end(),
                        [](const PauliString<Coeff>& ps) {
                                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                        return ps.coeff == std::complex<double>(0.0, 0.0);
                                } else {

                                        try {
                                                auto c = PauliString<Coeff>::to_complex(ps.coeff);
                                                return c == std::complex<double>(0.0, 0.0);
                                        } catch (...) {
                                                return false;
                                        }
                                }
                        }), data.end());

                return *this;
        }

    private:
     
        template <typename Cmp>
        void parallel_sort_data(Cmp cmp) {
#ifdef PAULIENGINE_HAS_OPENMP
                const int nthreads = omp_get_max_threads();
                if (nthreads <= 1) {
                        std::sort(data.begin(), data.end(), cmp);
                        return;
                }
                const size_t n = data.size();
                const size_t chunk = (n + nthreads - 1) / nthreads;

                #pragma omp parallel for schedule(static)
                for (int t = 0; t < nthreads; ++t) {
                        const size_t begin = static_cast<size_t>(t) * chunk;
                        if (begin >= n) continue;
                        const size_t end = std::min(begin + chunk, n);
                        std::sort(data.begin() + begin, data.begin() + end, cmp);
                }

                for (int step = 1; step < nthreads; step *= 2) {
                        #pragma omp parallel for schedule(dynamic)
                        for (int t = 0; t < nthreads; t += 2 * step) {
                                const size_t begin = static_cast<size_t>(t) * chunk;
                                const size_t mid   = std::min(static_cast<size_t>(t + step) * chunk, n);
                                const size_t end   = std::min(static_cast<size_t>(t + 2 * step) * chunk, n);
                                if (mid >= end || begin >= mid) continue;
                                std::inplace_merge(data.begin() + begin,
                                                   data.begin() + mid,
                                                   data.begin() + end,
                                                   cmp);
                        }
                }
#else
                std::sort(data.begin(), data.end(), cmp);
#endif
        }
    public:

        Hamiltonian_structure<Coeff> to_dictionary() const{
                Hamiltonian_structure<Coeff> output;
                output.reserve(this->data.size());
                for (const auto& entry : this->data) {
                        std::unordered_map<int, std::string> temp;
                        temp.reserve(entry.x.size() * BITS_IN_INTEGER);
                        for (size_t i = 0; i < entry.x.size(); ++i) {
                                uint64_t x_word = entry.x[i];
                                uint64_t y_word = entry.y[i];
                                for (int j = 0; j < BITS_IN_INTEGER; ++j) {
                                        uint64_t bit = (uint64_t)1 << j;
                                        bool x_bit = x_word & bit;
                                        bool y_bit = y_word & bit;
                                        if (!x_bit && !y_bit) continue;
                                        int index = i * BITS_IN_INTEGER + j;
                                        if (x_bit && y_bit) {
                                                temp[(int)index] = "Z";
                                        } else if (x_bit) {
                                                temp[(int)index] = "X";
                                        } else {
                                                temp[(int)index] = "Y";
                                        }
                                }
                        }
                        output.emplace_back(entry.coeff, std::move(temp));
                }
                return output;
        }

        std::vector<PauliString<Coeff>> to_list(){
                return this->data;
        }

        QubitHamiltonian commutator(const QubitHamiltonian& other) const {

                using Map = std::unordered_map<PauliString<Coeff>, Coeff,
                        PauliStringHash<Coeff>, PauliStringOperatorEqual<Coeff>>;

                const int n1 = static_cast<int>(this->data.size());
                const size_t n2 = other.data.size();

                auto insert_into = [](Map& m, PauliString<Coeff>&& temp) {
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                if (temp.coeff == std::complex<double>(0.0, 0.0)) return;
                        }
                        Coeff c = temp.coeff;  
                        auto [it, inserted] = m.try_emplace(std::move(temp), c);
                        if (!inserted) {
                                it->second = it->second + c;
                        }
                };

                Map merged;

                const bool parallel_ok =
#ifdef PAULIENGINE_HAS_OPENMP
                        std::is_same_v<Coeff, std::complex<double>>
                        && static_cast<size_t>(n1) * n2 >= PAULIENGINE_OMP_THRESHOLD;
#else
                        false;
#endif

                if (!parallel_ok) {
                        merged.reserve(static_cast<size_t>(n1) * n2 / 8 + 1);
                        for (int i = 0; i < n1; ++i) {
                                for (size_t j = 0; j < n2; ++j) {
                                        insert_into(merged, this->data[i].commutator(other.data[j]));
                                }
                        }
                } else {
#ifdef PAULIENGINE_HAS_OPENMP
                        const int nthreads = omp_get_max_threads();
                        std::vector<Map> local_maps(nthreads);
                        const size_t reserve_hint = (static_cast<size_t>(n1) * n2) / (8 * nthreads) + 1;
                        for (auto& m : local_maps) m.reserve(reserve_hint);

                        #pragma omp parallel for schedule(static)
                        for (int i = 0; i < n1; ++i) {
                                const int tid = omp_get_thread_num();
                                Map& local = local_maps[tid];
                                for (size_t j = 0; j < n2; ++j) {
                                        insert_into(local, this->data[i].commutator(other.data[j]));
                                }
                        }

                        // Sequential merge of local maps; use the largest as base.
                        size_t largest = 0;
                        for (size_t t = 1; t < local_maps.size(); ++t) {
                                if (local_maps[t].size() > local_maps[largest].size()) largest = t;
                        }
                        merged = std::move(local_maps[largest]);
                        for (size_t t = 0; t < local_maps.size(); ++t) {
                                if (t == largest) continue;
                                for (auto& kv : local_maps[t]) {
                                        auto [it, inserted] = merged.try_emplace(kv.first, kv.second);
                                        if (!inserted) {
                                                it->second = it->second + kv.second;
                                        }
                                }
                                local_maps[t].clear();
                        }
#endif
                }

                std::vector<PauliString<Coeff>> result;
                result.reserve(merged.size());
                while (!merged.empty()) {
                        auto node = merged.extract(merged.begin());
                        if (node.mapped() != Coeff(0.0)) {
                                node.key().coeff = node.mapped();
                                result.push_back(std::move(node.key()));
                        }
                }
                return QubitHamiltonian(std::move(result), already_compacted);
        }

        size_t size() const { return data.size(); }

        std::vector<int> qubits() const {
                std::set<int> qs;
                for (const auto& ps : data) {
                        for (size_t i = 0; i < ps.x.size(); ++i) {
                                uint64_t xy = ps.x[i] | ps.y[i];
                                while (xy) {
                                        int j = std::countr_zero(xy);
                                        qs.insert((int)(i * BITS_IN_INTEGER + j));
                                        xy &= xy - 1;
                                }
                        }
                }
                return std::vector<int>(qs.begin(), qs.end());
        }

        int n_qubits() const { return (int)qubits().size(); }

        bool is_all_z() const {
                for (const auto& ps : data) {
                        if (!ps.is_all_z()) return false;
                }
                return true;
        }

        int count_measurements() const {
                return is_all_z() ? 1 : (int)data.size();
        }

        bool is_hermitian() const {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        for (const auto& ps : data) {
                                if (std::abs(ps.coeff.imag()) > 1e-12) return false;
                        }
                        return true;
                } else {
                        try {
                                for (const auto& ps : data) {
                                        auto c = PauliString<Coeff>::to_complex(ps.coeff);
                                        if (std::abs(c.imag()) > 1e-12) return false;
                                }
                                return true;
                        } catch (...) {
                                return false;
                        }
                }
        }

        bool is_antihermitian() const {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        for (const auto& ps : data) {
                                if (std::abs(ps.coeff.real()) > 1e-12) return false;
                        }
                        return true;
                } else {
                        try {
                                for (const auto& ps : data) {
                                        auto c = PauliString<Coeff>::to_complex(ps.coeff);
                                        if (std::abs(c.real()) > 1e-12) return false;
                                }
                                return true;
                        } catch (...) {
                                return false;
                        }
                }
        }

        std::vector<PauliString<Coeff>> paulistrings() const {
                return this->data;
        }

        static QubitHamiltonian zero() {
                return QubitHamiltonian(std::vector<PauliString<Coeff>>{});
        }

        static QubitHamiltonian unit() {
                Coeff one = make_one();
                PauliString<Coeff> identity_ps(std::vector<uint64_t>{}, std::vector<uint64_t>{}, one);
                return QubitHamiltonian(std::vector<PauliString<Coeff>>{identity_ps});
        }

        QubitHamiltonian dagger() const {
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (const auto& ps : data) {
                        Coeff conj_coeff;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                conj_coeff = std::conj(ps.coeff);
                        } else {
                                conj_coeff = Expression(SymEngine::conjugate(ps.coeff.get_basic()));
                        }
                        result.push_back(PauliString<Coeff>(ps.x, ps.y, conj_coeff));
                }
                return QubitHamiltonian(std::move(result));
        }

        QubitHamiltonian conjugate() const {
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (const auto& ps : data) {
                        int y_count = ps.count_y();
                        Coeff new_coeff;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                new_coeff = std::conj(ps.coeff);
                        } else {
                                new_coeff = Expression(SymEngine::conjugate(ps.coeff.get_basic()));
                        }
                        if (y_count % 2 == 1) {
                                new_coeff = new_coeff * -1.0;
                        }
                        result.push_back(PauliString<Coeff>(ps.x, ps.y, new_coeff));
                }
                return QubitHamiltonian(std::move(result));
        }

        QubitHamiltonian transpose() const {
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (const auto& ps : data) {
                        int y_count = ps.count_y();
                        Coeff new_coeff = ps.coeff;
                        if (y_count % 2 == 1) {
                                new_coeff = new_coeff * -1.0;
                        }
                        result.push_back(PauliString<Coeff>(ps.x, ps.y, new_coeff));
                }
                return QubitHamiltonian(std::move(result));
        }

        QubitHamiltonian simplify(double threshold = 0.0) const {
                std::vector<PauliString<Coeff>> result;
                for (const auto& ps : data) {
                        bool keep;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                keep = std::abs(ps.coeff) > threshold;
                        } else {
                                try {
                                        auto c = PauliString<Coeff>::to_complex(ps.coeff);
                                        keep = std::abs(c) > threshold;
                                } catch (...) {
                                        keep = true;
                                }
                        }
                        if (keep) result.push_back(ps);
                }
                return QubitHamiltonian(std::move(result));
        }

        std::pair<QubitHamiltonian, QubitHamiltonian> split() const {
                std::vector<PauliString<Coeff>> herm_terms;
                std::vector<PauliString<Coeff>> anti_terms;
                for (const auto& ps : data) {
                        std::complex<double> c;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                c = ps.coeff;
                        } else {
                                try {
                                        c = PauliString<Coeff>::to_complex(ps.coeff);
                                } catch (...) {
                                        throw std::runtime_error("split() requires fully evaluated numeric coefficients");
                                }
                        }
                        if (c.real() != 0.0) {
                                Coeff real_coeff;
                                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                        real_coeff = std::complex<double>(c.real(), 0.0);
                                } else {
                                        real_coeff = Expression(c.real());
                                }
                                herm_terms.push_back(PauliString<Coeff>(ps.x, ps.y, real_coeff));
                        }
                        if (c.imag() != 0.0) {
                                Coeff imag_coeff;
                                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                        imag_coeff = std::complex<double>(0.0, c.imag());
                                } else {
                                        imag_coeff = Expression(std::complex<double>(0.0, c.imag()));
                                }
                                anti_terms.push_back(PauliString<Coeff>(ps.x, ps.y, imag_coeff));
                        }
                }
                return {QubitHamiltonian(std::move(herm_terms)), QubitHamiltonian(std::move(anti_terms))};
        }

        QubitHamiltonian map_qubits(std::unordered_map<int, int> qubit_map) const {
                std::vector<PauliString<Coeff>> result;
                result.reserve(data.size());
                for (auto ps : data) {
                        result.push_back(ps.map_qubits(qubit_map));
                }
                return QubitHamiltonian(std::move(result));
        }

        QubitHamiltonian power(int n) const {
                if (n < 0) {
                        throw std::runtime_error("Negative power not supported for QubitHamiltonian");
                }
                if (n == 0) return unit();
                QubitHamiltonian result = *this;
                for (int i = 1; i < n; ++i) {
                        result = result * *this;
                }
                return result;
        }

        Matrix2D to_matrix(bool ignore_unused_qubits = true) const {
                std::vector<int> q_list = qubits();
                int nq;
                if (ignore_unused_qubits) {
                        nq = (int)q_list.size();
                } else if (q_list.empty()) {
                        nq = 0;
                } else {
                        nq = q_list.back() + 1;
                }

                size_t dim = (nq == 0) ? 1 : ((size_t)1 << nq);
                Matrix2D total(dim, std::vector<std::complex<double>>(dim, 0.0));

                std::unordered_map<int, int> q_to_pos;
                if (ignore_unused_qubits) {
                        for (int i = 0; i < (int)q_list.size(); ++i) {
                                q_to_pos[q_list[i]] = i;
                        }
                }

                for (const auto& ps : data) {
                        std::complex<double> c;
                        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                                c = ps.coeff;
                        } else {
                                c = PauliString<Coeff>::to_complex(ps.coeff);
                        }

                        auto dict = ps.to_dictionary().second;
                        std::unordered_map<int, std::string> remapped;
                        for (const auto& [q, op] : dict) {
                                int pos = ignore_unused_qubits ? q_to_pos.at(q) : q;
                                remapped[pos] = op;
                        }

                        Matrix2D mat = {{std::complex<double>(1.0, 0.0)}};
                        for (int q = 0; q < nq; ++q) {
                                auto it = remapped.find(q);
                                std::string op = (it != remapped.end()) ? it->second : "I";
                                mat = tensor_product(mat, get_pauli_matrix(op));
                        }

                        for (size_t i = 0; i < dim; ++i) {
                                for (size_t j = 0; j < dim; ++j) {
                                        total[i][j] += c * mat[i][j];
                                }
                        }
                }
                return total;
        }


    private:
        static Coeff make_one() {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        return std::complex<double>(1.0, 0.0);
                } else {
                        return Coeff(1);
                }
        }
    static Matrix2D get_pauli_matrix(const std::string& p){static const std::complex<double> I(0,1);
        if (p == "I") return {{1,0},{0,1}};
        else if (p == "X") return {{0,1},{1,0}};
        else if (p == "Y") return {{0,-I},{I,0}};
        else if (p == "Z") return {{1,0},{0,-1}};
        else throw std::runtime_error("Unknown Pauli operator");
}
    static Matrix2D tensor_product(const Matrix2D& A, const Matrix2D& B){
        int rows = A.size() * B.size();
        int cols = A[0].size() * B[0].size();
        Matrix2D result(rows, std::vector<std::complex<double>>(cols, 0));

        for (int i=0; i<(int)A.size(); ++i) {
                for (int j=0; j<(int)A[0].size(); ++j) {
                        for (int k=0; k<(int)B.size(); ++k) {
                                for (int l=0; l<(int)B[0].size(); ++l) {
                                        result[i*B.size()+k][j*B[0].size()+l] = A[i][j] * B[k][l];
                                }
                        }
                }
        }
        return result;
}
    static Matrix2D pauli_string_to_matrix(const std::unordered_map<int, std::string>& pauli_map, int num_qubits){
        Matrix2D result = get_pauli_matrix("I");
        for (int q=0; q < num_qubits; ++q) {
                std::string op = "I";
                auto it = pauli_map.find(q);
                if (it != pauli_map.end()) op = it->second;
                Matrix2D pmat = get_pauli_matrix(op);
                result = tensor_product(result, pmat);
        }
        return result;
}

};
