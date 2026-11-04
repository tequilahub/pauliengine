#pragma once

/// @file
/// Pauli decomposition: any 2^n x 2^n matrix as a QubitHamiltonian.
///
/// Implements the method of T. N. Georges et al., New J. Phys. (2025),
/// https://doi.org/10.1088/1367-2630/adb44d: a XOR permutation of the matrix
/// entries followed by a Walsh-Hadamard transform of every row. It needs
/// O(4^n n) time and works on one copy of the matrix.
///
/// Conventions match QubitHamiltonian::to_matrix(): qubit 0 is the most
/// significant bit of the row and column index, so
/// pauli_decomposition(H.to_matrix(false)) reproduces H.

#include <bit>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/QubitHamiltonian.h"

namespace pauliengine_detail {

/// In place: turns the row-major dim x dim matrix `data` into its Pauli
/// coefficients. Afterwards data[i * dim + j] is the coefficient of the Pauli
/// string with x bits (i ^ j) and y bits j (see pauli_decomposition).
inline void pauli_coefficients(std::complex<double>* data, size_t dim) {
        // XOR permutation: exchange (j, i ^ j) with (i, i ^ j) for all i < j.
        for (size_t j = 0; j < dim; ++j) {
                for (size_t i = 0; i < j; ++i) {
                        std::swap(data[j * dim + (i ^ j)], data[i * dim + (i ^ j)]);
                }
        }

        const double norm = 1.0 / static_cast<double>(dim);   // exact: dim is a power of two
        for (size_t j = 0; j < dim; ++j) {
                std::complex<double>* row = data + j * dim;
                // Walsh-Hadamard transform of row j.
                for (size_t half = 1; half < dim; half *= 2) {
                        for (size_t start = 0; start < dim; start += 2 * half) {
                                for (size_t i = start; i < start + half; ++i) {
                                        const std::complex<double> a = row[i];
                                        const std::complex<double> b = row[i + half];
                                        row[i] = a + b;
                                        row[i + half] = a - b;
                                }
                        }
                }
                // Phase (-i)^popcount(i & j) and normalisation 1 / dim.
                for (size_t i = 0; i < dim; ++i) {
                        const std::complex<double> v = row[i];
                        switch (std::popcount(static_cast<uint64_t>(i & j)) & 3) {
                                case 0: row[i] = v * norm; break;
                                case 1: row[i] = std::complex<double>(v.imag(), -v.real()) * norm; break;  // * (-i)
                                case 2: row[i] = -v * norm; break;
                                default: row[i] = std::complex<double>(-v.imag(), v.real()) * norm; break; // * i
                        }
                }
        }
}

/// Number of qubits for a matrix of side `dim`; throws unless dim is a power of two.
inline int qubits_for_dimension(size_t dim) {
        if (dim == 0 || (dim & (dim - 1)) != 0) {
                throw std::invalid_argument("pauli_decomposition: the matrix side must be a power of two, got "
                        + std::to_string(dim));
        }
        const int n = std::countr_zero(static_cast<uint64_t>(dim));
        if (n > 31) {
                throw std::invalid_argument("pauli_decomposition: too many qubits (" + std::to_string(n) + ")");
        }
        return n;
}

/// Reverses the lowest n bits of v (qubit 0 is the most significant index bit).
inline uint64_t reverse_bits(uint64_t v, int n) {
        uint64_t r = 0;
        for (int k = 0; k < n; ++k) {
                r = (r << 1) | ((v >> k) & 1ULL);
        }
        return r;
}

/// Pauli decomposition of a row-major matrix that may be overwritten.
inline QubitHamiltonian<std::complex<double>> pauli_decomposition_inplace(
                std::complex<double>* data, size_t dim, double threshold) {
        const int n = qubits_for_dimension(dim);
        pauli_coefficients(data, dim);

        // Coefficient (i, j) belongs to the string with x bits (i ^ j) and y bits j,
        // read with qubit 0 as the most significant bit: X = (1,0), Y = (0,1), Z = (1,1).
        std::vector<uint64_t> rev(dim);
        for (size_t k = 0; k < dim; ++k) {
                rev[k] = reverse_bits(k, n);
        }
        std::vector<PauliString<std::complex<double>>> terms;
        for (size_t i = 0; i < dim; ++i) {
                for (size_t j = 0; j < dim; ++j) {
                        const std::complex<double> c = data[i * dim + j];
                        if (std::abs(c) <= threshold) {
                                continue;
                        }
                        WordVec x(1), y(1);
                        x[0] = rev[i ^ j];
                        y[0] = rev[j];
                        terms.emplace_back(std::move(x), std::move(y), c);
                }
        }
        // every (i, j) gives a different string and zeros are already dropped
        return QubitHamiltonian<std::complex<double>>(std::move(terms), already_compacted);
}

}  // namespace pauliengine_detail

/// Pauli decomposition of a 2^n x 2^n matrix given row-major in `data`
/// (length dim * dim). The input is not modified.
///
/// Terms with |coefficient| <= threshold are dropped; the default 0 drops only
/// exact zeros. Throws std::invalid_argument unless dim is a power of two and
/// data has dim * dim entries.
inline QubitHamiltonian<std::complex<double>> pauli_decomposition(
                std::vector<std::complex<double>> data, size_t dim, double threshold = 0.0) {
        if (data.size() != dim * dim) {
                throw std::invalid_argument("pauli_decomposition: expected " + std::to_string(dim * dim)
                        + " entries for a " + std::to_string(dim) + " x " + std::to_string(dim)
                        + " matrix, got " + std::to_string(data.size()));
        }
        return pauliengine_detail::pauli_decomposition_inplace(data.data(), dim, threshold);
}

/// Pauli decomposition of a square matrix given as rows (the format returned by
/// QubitHamiltonian::to_matrix). The input is not modified.
inline QubitHamiltonian<std::complex<double>> pauli_decomposition(const Matrix2D& matrix, double threshold = 0.0) {
        const size_t dim = matrix.size();
        std::vector<std::complex<double>> data;
        data.reserve(dim * dim);
        for (const auto& row : matrix) {
                if (row.size() != dim) {
                        throw std::invalid_argument("pauli_decomposition: the matrix must be square");
                }
                data.insert(data.end(), row.begin(), row.end());
        }
        return pauliengine_detail::pauli_decomposition_inplace(data.data(), dim, threshold);
}
