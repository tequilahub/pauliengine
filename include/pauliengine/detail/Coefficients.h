#pragma once

/// @file
/// Internal: small helpers that let the operator classes work with both coefficient
/// types, std::complex<double> and the polynomial Expression.

#include <complex>
#include <string>
#include <type_traits>
#include <unordered_map>

#include "pauliengine/Polynomial.h"

namespace pauliengine_detail {

/// One as a coefficient of type C.
template <typename C>
C coeff_one() {
        if constexpr (std::is_same_v<C, std::complex<double>>) {
                return C(1.0, 0.0);
        } else {
                return C(1);
        }
}

/// Zero as a coefficient of type C.
template <typename C>
C coeff_zero() {
        if constexpr (std::is_same_v<C, std::complex<double>>) {
                return C(0.0, 0.0);
        } else {
                return C(0);
        }
}

/// Whether c is zero. Polynomials are canonical, so this needs no expansion.
inline bool coeff_is_zero(const std::complex<double>& c) { return c == std::complex<double>(0.0, 0.0); }
inline bool coeff_is_zero(const pauliengine::symbolic::Expression& c) { return c.is_zero(); }

/// A coefficient converted to type To (complex to symbolic, or unchanged).
template <typename To, typename From>
To coeff_as(const From& c) {
        if constexpr (std::is_same_v<To, From>) {
                return c;
        } else {
                return To(c);
        }
}

/// Values for symbols by name, as taken by subs().
using SymbolValues = std::unordered_map<std::string, std::complex<double>>;

/// Polynomials for symbols by name, as taken by subs().
using SymbolPolys = std::unordered_map<std::string, pauliengine::symbolic::Expression>;

/// Derivative of a symbolic coefficient by the symbol `name`.
inline pauliengine::symbolic::Expression coeff_diff(const pauliengine::symbolic::Expression& c, const std::string& name) {
        return c.diff(name);
}

/// A symbolic coefficient with the given values or polynomials substituted for its
/// symbols; symbols without one stay. `values` is a SymbolValues or SymbolPolys
/// map, or one prepared with Expression::prepare_subs for many coefficients.
template <typename Values>
pauliengine::symbolic::Expression coeff_subs(const pauliengine::symbolic::Expression& c, const Values& values) {
        return c.subs(values);
}

}  // namespace pauliengine_detail
