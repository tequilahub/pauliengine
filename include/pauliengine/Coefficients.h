#pragma once

/// @file
/// Small helpers that let the cycle and orbit classes work with both coefficient
/// types, std::complex<double> and the polynomial Expression.

#include <complex>
#include <type_traits>

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

}  // namespace pauliengine_detail
