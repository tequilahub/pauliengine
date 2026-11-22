#pragma once
/// Internal: flat storage of the terms of a QubitHamiltonian with complex
/// coefficients, and small helpers on single terms.

#include <algorithm>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <new>
#include <utility>

#include "pauliengine/PauliString.h"
#include "pauliengine/detail/Parallel.h"

namespace pauliengine_detail {

/// One term with a complex coefficient, read from a string or from flat storage.
/// x and y hold nw words each.
struct TermView {
        const uint64_t* x;
        const uint64_t* y;
        size_t nw;
        std::complex<double> coeff;
};

inline TermView view_of(const PauliString<std::complex<double>>& ps) {
        return {ps.x.data(), ps.y.data(), ps.x.size(), ps.coeff};
}

/// scale * i^power_of_i * a * b.
inline std::complex<double> product_coeff(std::complex<double> a, std::complex<double> b,
                                          int power_of_i, double scale) {
        static constexpr double kRealSign[4] = {1.0, -1.0, -1.0, 1.0};
        static constexpr double kImagSign[4] = {1.0, 1.0, -1.0, -1.0};
        const double parts[2] = {a.real() * b.real() - a.imag() * b.imag(),
                                 a.real() * b.imag() + a.imag() * b.real()};
        const int k = power_of_i & 3;
        return {scale * kRealSign[k] * parts[k & 1], scale * kImagSign[k] * parts[(k & 1) ^ 1]};
}

/// Whether two strings anti-commute: the symplectic product has odd parity. Much
/// cheaper than forming the product, which a commutator needs only for these pairs
/// and a square only for the others.
inline bool anticommute_words(const uint64_t* x1, const uint64_t* y1, size_t n1,
                              const uint64_t* x2, const uint64_t* y2, size_t n2) {
        const size_t n = std::min(n1, n2);
        uint64_t anti = 0;
        for (size_t k = 0; k < n; ++k) anti ^= (x1[k] & y2[k]) ^ (x2[k] & y1[k]);
        return std::popcount(anti) & 1;
}

inline bool anticommute(const TermView& a, const TermView& b) {
        return anticommute_words(a.x, a.y, a.nw, b.x, b.y, b.nw);
}

struct FreeDeleter {
        void operator()(void* p) const noexcept { std::free(p); }
};

/// Uninitialised memory for `count` objects of a trivial type; null if count is zero.
template<typename T>
inline T* allocate_block(size_t count) {
        if (count == 0) return nullptr;
        void* p = std::malloc(count * sizeof(T));
        if (p == nullptr) throw std::bad_alloc();
        return static_cast<T*>(p);
}

/// The terms of a Hamiltonian with complex coefficients, stored flat. Row i holds
/// the x words and then the y words of term i, both padded with zeros to the common
/// width; the coefficients are a second array of (real, imaginary) pairs. A sum of
/// n terms takes two allocations, instead of one string object per term plus heap
/// words for every string beyond 128 qubits.
class FlatTerms {
    public:
        FlatTerms() = default;

        /// n terms of nw words each, contents uninitialised.
        FlatTerms(size_t n, size_t nw)
                : n_(n), nw_(nw),
                  words_(allocate_block<uint64_t>(2 * n * nw)),
                  coeffs_(allocate_block<double>(2 * n)) {}

        FlatTerms(const FlatTerms& o) : FlatTerms(o.n_, o.nw_) {
                std::copy(o.words_.get(), o.words_.get() + 2 * n_ * nw_, words_.get());
                std::copy(o.coeffs_.get(), o.coeffs_.get() + 2 * n_, coeffs_.get());
        }

        FlatTerms(FlatTerms&& o) noexcept
                : n_(std::exchange(o.n_, 0)), nw_(std::exchange(o.nw_, 0)),
                  words_(std::move(o.words_)), coeffs_(std::move(o.coeffs_)) {}

        FlatTerms& operator=(const FlatTerms& o) {
                if (this != &o) *this = FlatTerms(o);
                return *this;
        }

        FlatTerms& operator=(FlatTerms&& o) noexcept {
                n_ = std::exchange(o.n_, 0);
                nw_ = std::exchange(o.nw_, 0);
                words_ = std::move(o.words_);
                coeffs_ = std::move(o.coeffs_);
                return *this;
        }

        size_t size() const noexcept { return n_; }
        /// Words per term, the same for x and y.
        size_t width() const noexcept { return nw_; }

        uint64_t* x(size_t i) noexcept { return words_.get() + 2 * nw_ * i; }
        uint64_t* y(size_t i) noexcept { return x(i) + nw_; }
        const uint64_t* x(size_t i) const noexcept { return words_.get() + 2 * nw_ * i; }
        const uint64_t* y(size_t i) const noexcept { return x(i) + nw_; }

        std::complex<double> coeff(size_t i) const noexcept { return {coeffs_[2 * i], coeffs_[2 * i + 1]}; }
        void set_coeff(size_t i, std::complex<double> c) noexcept {
                coeffs_[2 * i] = c.real();
                coeffs_[2 * i + 1] = c.imag();
        }

        TermView view(size_t i) const noexcept { return {x(i), y(i), nw_, coeff(i)}; }

        /// Words of term i without its trailing zero words, as a PauliString keeps them.
        size_t used_words(size_t i) const noexcept {
                size_t n = nw_;
                while (n > 0 && x(i)[n - 1] == 0 && y(i)[n - 1] == 0) --n;
                return n;
        }

        /// Sets row i to the term `v`, which has at most width() words.
        void set_row(size_t i, const TermView& v) noexcept {
                uint64_t* rx = x(i);
                uint64_t* ry = y(i);
                for (size_t k = 0; k < nw_; ++k) {
                        rx[k] = k < v.nw ? v.x[k] : 0;
                        ry[k] = k < v.nw ? v.y[k] : 0;
                }
                set_coeff(i, v.coeff);
        }

        /// Sets row i to the operator part of a * b and returns the exponent k of the
        /// phase i^k; the coefficient is left to the caller.
        int set_row_product(size_t i, const TermView& a, const TermView& b) noexcept {
                uint64_t* rx = x(i);
                uint64_t* ry = y(i);
                const int power_of_i = multiply_symplectic_words(a.x, a.y, a.nw, b.x, b.y, b.nw, rx, ry);
                for (size_t k = std::max(a.nw, b.nw); k < nw_; ++k) {
                        rx[k] = 0;
                        ry[k] = 0;
                }
                return power_of_i;
        }

        PauliString<std::complex<double>> term(size_t i) const {
                PauliString<std::complex<double>> ps;
                const size_t nw = used_words(i);
                ps.x.assign(x(i), x(i) + nw);
                ps.y.assign(y(i), y(i) + nw);
                ps.coeff = coeff(i);
                return ps;
        }

        std::vector<PauliString<std::complex<double>>> terms() const {
                std::vector<PauliString<std::complex<double>>> out(n_);
                for_each_index(n_, use_threads(n_ * nw_), [&](size_t i) { out[i] = term(i); });
                return out;
        }

        /// Multiplies every coefficient by s; a real s scales both parts alone.
        template<typename S>
        void scale(S s) noexcept {
                for (size_t i = 0; i < n_; ++i) set_coeff(i, coeff(i) * s);
        }

        /// Removes the terms with coefficient zero, keeping the order of the others.
        void drop_zero_terms() noexcept {
                size_t write = 0;
                for (size_t read = 0; read < n_; ++read) {
                        if (coeffs_[2 * read] == 0.0 && coeffs_[2 * read + 1] == 0.0) continue;
                        if (write != read) {
                                std::copy(x(read), x(read) + 2 * nw_, x(write));
                                coeffs_[2 * write] = coeffs_[2 * read];
                                coeffs_[2 * write + 1] = coeffs_[2 * read + 1];
                        }
                        ++write;
                }
                n_ = write;
        }

        bool has_zero_coeff() const noexcept {
                for (size_t i = 0; i < n_; ++i) {
                        if (coeffs_[2 * i] == 0.0 && coeffs_[2 * i + 1] == 0.0) return true;
                }
                return false;
        }

        /// The (real, imaginary) pairs of terms i, i + 1, ...
        double* coeff_data(size_t i) noexcept { return coeffs_.get() + 2 * i; }

        /// Keeps the first n terms. The memory of the others is given back if it
        /// is a relevant share; shrinking in place moves nothing.
        void shrink_to(size_t n) noexcept {
                if (n >= n_) return;
                const bool release = (n_ - n) >= n_ / 8;
                n_ = n;
                if (release) {
                        shrink_block(words_, 2 * n * nw_);
                        shrink_block(coeffs_, 2 * n);
                }
        }

        /// Changes the number of terms to n, keeping the first min(n, size()) of them;
        /// added terms are uninitialised.
        void resize(size_t n) {
                if (n == n_) return;
                grow_block(words_, 2 * n * nw_);
                grow_block(coeffs_, 2 * n);
                n_ = n;
        }

        /// Sets rows first, first + 1, ... to `count` terms given as rows of
        /// 2 * width() words and as (real, imaginary) pairs.
        void assign_rows(size_t first, const uint64_t* words, const double* coeffs, size_t count) noexcept {
                if (count == 0) return;
                std::copy(words, words + 2 * nw_ * count, words_.get() + 2 * nw_ * first);
                std::copy(coeffs, coeffs + 2 * count, coeffs_.get() + 2 * first);
        }

        /// Moves `count` terms from row `from` down to row `to` <= from; the ranges may overlap.
        void move_rows(size_t to, size_t from, size_t count) noexcept {
                if (count == 0) return;
                std::copy(x(from), x(from) + 2 * nw_ * count, x(to));
                std::copy(coeff_data(from), coeff_data(from) + 2 * count, coeff_data(to));
        }

    private:
        size_t n_ = 0;
        size_t nw_ = 0;
        template<typename T>
        static void shrink_block(std::unique_ptr<T[], FreeDeleter>& block, size_t count) noexcept {
                if (count == 0) {
                        block.reset();
                        return;
                }
                void* p = std::realloc(block.get(), count * sizeof(T));
                if (p == nullptr) return;  // the larger block stays valid
                block.release();
                block.reset(static_cast<T*>(p));
        }

        template<typename T>
        static void grow_block(std::unique_ptr<T[], FreeDeleter>& block, size_t count) {
                if (count == 0) {
                        block.reset();
                        return;
                }
                void* p = std::realloc(block.get(), count * sizeof(T));
                if (p == nullptr) throw std::bad_alloc();  // the old block stays valid
                block.release();
                block.reset(static_cast<T*>(p));
        }

        std::unique_ptr<uint64_t[], FreeDeleter> words_;
        std::unique_ptr<double[], FreeDeleter> coeffs_;
};

}  // namespace pauliengine_detail
