#pragma once

/// @file
/// Small helpers that let the cycle and orbit classes work with both coefficient
/// types, std::complex<double> and SymEngine::Expression.

#include <complex>
#include <type_traits>

#include <symengine/basic.h>
#include <symengine/expression.h>
#include <symengine/number.h>

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

/// Whether c is zero. A symbolic coefficient is expanded first, since SymEngine
/// does not reduce for example (2 + a) - (2 + a) to 0 on its own.
inline bool coeff_is_zero(const std::complex<double>& c) { return c == std::complex<double>(0.0, 0.0); }
inline bool coeff_is_zero(const SymEngine::Expression& c) {
        const auto& direct = *c.get_basic();
        if (SymEngine::is_a_Number(direct)) {
                return SymEngine::down_cast<const SymEngine::Number&>(direct).is_zero();
        }
        const SymEngine::RCP<const SymEngine::Basic> expanded = SymEngine::expand(c.get_basic());
        return SymEngine::is_a_Number(*expanded) && SymEngine::down_cast<const SymEngine::Number&>(*expanded).is_zero();
}

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
