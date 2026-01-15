#pragma once

#include <unordered_map>
#include <string>
#include <complex>
#include <vector>
#include <omp.h>

#ifdef HAVE_SYMENGINE
#include <symengine/expression.h>
using SymEngine::Expression;
#else
#endif

#define BITS_IN_INTEGER (sizeof(uint64_t) * 8)

using Pauli_structure = std::pair<std::complex<double>, std::unordered_map<int, std::string>>;
#ifdef HAVE_SYMENGINE
using Pauli_structure_variable = std::pair<SymEngine::Expression, std::unordered_map<int, std::string>>;
using Hamiltonian_structure_variable = std::vector<Pauli_structure_variable>;
#endif
using Hamiltonian_structure = std::vector<Pauli_structure>;
using Matrix2D = std::vector<std::vector<std::complex<double>>>;

const std::complex<double> neg_I(0.0 , -1.0);


