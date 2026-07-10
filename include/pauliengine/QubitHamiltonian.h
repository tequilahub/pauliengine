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

                size_t max_words = 0;
                for (const auto& ps : data) max_words = std::max(max_words, ps.x.size());
                if (try_compact_small_key<>(max_words)) {
                        return *this;
                }

                // sort-based compact. Single-pass ordering: by word
                // count first (constructors trim trailing zero words, so equal
                // operators always have equal sizes), then by content.
                auto op_less = [](const PauliString<Coeff>& a, const PauliString<Coeff>& b) {
                        const size_t na = a.x.size();
                        const size_t nb = b.x.size();
                        if (na != nb) return na < nb;
                        for (size_t i = 0; i < na; ++i) {
                                if (a.x[i] != b.x[i]) return a.x[i] < b.x[i];
                        }
                        for (size_t i = 0; i < na; ++i) {
                                if (a.y[i] != b.y[i]) return a.y[i] < b.y[i];
                        }
                        return false;
                };
                auto op_equal = [](const PauliString<Coeff>& a, const PauliString<Coeff>& b) {
                        return a.x == b.x && a.y == b.y;
                };

                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        parallel_sort_vec(data, op_less);
                } else {
                        // symengine coefficients are not safe to move across threads
                        std::sort(data.begin(), data.end(), op_less);
                }

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
                        [](const PauliString<Coeff>& ps) { return coeff_is_zero(ps.coeff); }),
                        data.end());

                return *this;
        }

    private:
        static bool coeff_is_zero(const Coeff& c) {
                if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                        return c == std::complex<double>(0.0, 0.0);
                } else {
                    const auto& basic = *c.get_basic();
                        return SymEngine::is_a_Number(basic)
                                && SymEngine::down_cast<const SymEngine::Number&>(basic).is_zero();
                }
        }

        // Largest word count handled by the small-key fast path (8 words =
        // 512 qubits). Beyond this, compact() falls back to the generic sort.
        static constexpr size_t kMaxSmallKeyWords = 8;

        // Compile-time dispatch: picks the smallest W in {1, 2, 4, ...,
        // kMaxSmallKeyWords} that fits max_words and runs compact_small_key<W>.
        // Returns false when max_words exceeds the fast-path limit.
        template<size_t W = 1>
        bool try_compact_small_key(size_t max_words) {
                if constexpr (W > kMaxSmallKeyWords) {
                        return false;
                } else {
                        if (max_words <= W) {
                                compact_small_key<W>();
                                return true;
                        }
                        return try_compact_small_key<W * 2>(max_words);
                }
        }

        template<size_t W>
        void compact_small_key() {
                struct Entry {
                        uint64_t key[2 * W];
                        uint32_t idx;
                };
                const size_t n = data.size();
                std::vector<Entry> entries(n);
                for (size_t i = 0; i < n; ++i) {
                        Entry& e = entries[i];
                        const PauliString<Coeff>& ps = data[i];
                        const size_t nw = ps.x.size();
                        for (size_t k = 0; k < W; ++k) {
                                e.key[k]     = (k < nw) ? ps.x[k] : 0;
                                e.key[W + k] = (k < nw) ? ps.y[k] : 0;
                        }
                        e.idx = static_cast<uint32_t>(i);
                }

                auto key_less = [](const Entry& a, const Entry& b) {
                        for (size_t k = 0; k < 2 * W; ++k) {
                                if (a.key[k] != b.key[k]) return a.key[k] < b.key[k];
                        }
                        return false;
                };
                auto key_equal = [](const Entry& a, const Entry& b) {
                        for (size_t k = 0; k < 2 * W; ++k) {
                                if (a.key[k] != b.key[k]) return false;
                        }
                        return true;
                };

                parallel_sort_vec(entries, key_less);

                std::vector<PauliString<Coeff>> out;
                out.reserve(n);
                size_t i = 0;
                while (i < n) {
                        Coeff acc = std::move(data[entries[i].idx].coeff);
                        size_t j = i + 1;
                        while (j < n && key_equal(entries[i], entries[j])) {
                                acc = acc + data[entries[j].idx].coeff;
                                ++j;
                        }
                        if (!coeff_is_zero(acc)) {
                                PauliString<Coeff> ps = std::move(data[entries[i].idx]);
                                ps.coeff = std::move(acc);
                                out.push_back(std::move(ps));
                        }
                        i = j;
                }
                data = std::move(out);
        }

        // Chunked parallel merge sort; falls back to std::sort without OpenMP
        // or for small inputs. Elements must be safe to move across threads.
        template<typename T, typename Cmp>
        static void parallel_sort_vec(std::vector<T>& v, Cmp cmp) {
#ifdef PAULIENGINE_HAS_OPENMP
                const size_t n = v.size();
                const int nthreads = omp_get_max_threads();
                if (n >= PAULIENGINE_OMP_THRESHOLD && nthreads > 1) {
                        const size_t chunk = (n + nthreads - 1) / nthreads;

                        #pragma omp parallel for schedule(static)
                        for (int t = 0; t < nthreads; ++t) {
                                const size_t begin = static_cast<size_t>(t) * chunk;
                                if (begin >= n) continue;
                                const size_t end = std::min(begin + chunk, n);
                                std::sort(v.begin() + begin, v.begin() + end, cmp);
                        }

                        for (int step = 1; step < nthreads; step *= 2) {
                                #pragma omp parallel for schedule(dynamic)
                                for (int t = 0; t < nthreads; t += 2 * step) {
                                        const size_t begin = static_cast<size_t>(t) * chunk;
                                        const size_t mid   = std::min(static_cast<size_t>(t + step) * chunk, n);
                                        const size_t end   = std::min(static_cast<size_t>(t + 2 * step) * chunk, n);
                                        if (mid >= end || begin >= mid) continue;
                                        std::inplace_merge(v.begin() + begin,
                                                           v.begin() + mid,
                                                           v.begin() + end,
                                                           cmp);
                                }
                        }
                        return;
                }
#endif
                std::sort(v.begin(), v.end(), cmp);
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
                const int n1 = static_cast<int>(this->data.size()); 
                const size_t n2 = other.data.size();
                const size_t pairs = static_cast<size_t>(n1) * n2;


                std::vector<PauliString<Coeff>> terms;

                const bool parallel_ok =
#ifdef PAULIENGINE_HAS_OPENMP
                        std::is_same_v<Coeff, std::complex<double>>
                        && pairs >= PAULIENGINE_OMP_THRESHOLD;
#else
                        false;
#endif

                if (!parallel_ok) {
                        terms.reserve(pairs / 8 + 1);
                        for (int i = 0; i < n1; ++i) {
                                for (size_t j = 0; j < n2; ++j) {
                                        PauliString<Coeff> t = this->data[i].commutator(other.data[j]);
                                        if (!t.x.empty()) {
                                                terms.push_back(std::move(t));
                                        }
                                }
                        }
                } else {
#ifdef PAULIENGINE_HAS_OPENMP
                        const int nthreads = omp_get_max_threads();
                        std::vector<std::vector<PauliString<Coeff>>> local(nthreads);
                        for (auto& v : local) v.reserve(pairs / (8 * nthreads) + 1);

                        #pragma omp parallel for schedule(static)
                        for (int i = 0; i < n1; ++i) {
                                auto& mine = local[omp_get_thread_num()];
                                for (size_t j = 0; j < n2; ++j) {
                                        PauliString<Coeff> t = this->data[i].commutator(other.data[j]);
                                        if (!t.x.empty()) {
                                                mine.push_back(std::move(t));
                                        }
                                }
                        }

                        size_t total = 0;
                        for (const auto& v : local) total += v.size();
                        terms.reserve(total);
                        for (auto& v : local) {
                                for (auto& t : v) terms.push_back(std::move(t));
                                v.clear();
                        }
#endif
                }
                return QubitHamiltonian(std::move(terms));  // constructor calls compact()
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
