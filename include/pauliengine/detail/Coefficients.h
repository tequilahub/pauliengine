#pragma once

/// @file
/// Internal: small helpers that let the operator classes work with both coefficient
/// types, std::complex<double> and SymEngine::Expression.

#include <complex>
#include <string>
#include <type_traits>
#include <unordered_map>

#include <symengine/basic.h>
#include <symengine/expression.h>
#include <symengine/number.h>
#include <symengine/symbol.h>

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

/// Values for symbols by name, as taken by subs().
using SymbolValues = std::unordered_map<std::string, std::complex<double>>;

/// Derivative of a symbolic coefficient by the symbol `name`.
inline SymEngine::Expression coeff_diff(const SymEngine::Expression& c, const std::string& name) {
        return c.diff(SymEngine::symbol(name));
}

/// A symbolic coefficient with the given values substituted for its symbols;
/// symbols without a value stay.
inline SymEngine::Expression coeff_subs(const SymEngine::Expression& c, const SymbolValues& values) {
        SymEngine::map_basic_basic substitutions;
        for (const auto& [name, value] : values) {
                substitutions[SymEngine::symbol(name)] = SymEngine::Expression(value).get_basic();
        }
        return c.subs(substitutions);
}

}  // namespace pauliengine_detail
