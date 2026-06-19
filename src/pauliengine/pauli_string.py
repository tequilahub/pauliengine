"""Module for PauliString classes and factory."""

from __future__ import annotations

from collections.abc import Iterable
from ._core import (
    Expression,
    PauliStringComplex,
    PauliStringSymbolic,
)


class PauliStringFactory:
    """An implementation of the Factory method pattern."""

    __slots__ = ["_builders", "_name"]

    def __init__(self, name: str, builders: dict[type, type]) -> None:
        self._builders = builders
        self._name = name

    def _convert_terms(self, terms) -> dict[int, str]:
        """Convert various input formats into {index: 'Pauli'} dict."""
        if isinstance(terms, dict):
            return terms

        if isinstance(terms, str):
            result = {}
            for part in terms.split():
                pauli = part[0]
                index = int(part[1:])
                result[index] = pauli
            return result

        if isinstance(terms, Iterable):
            try:
                # Expect [(Pauli, index), ...]
                return {i: p for p, i in terms}
            except Exception as e:
                raise ValueError(
                    f"{self._name}: invalid iterable format for Pauli terms"
                ) from e

        raise ValueError(
            f"{self._name}: unsupported pauli_strings format: {type(terms)}"
        )

    def __call__(
        self,
        coeff: complex | str | tuple = 0,
        pauli_strings: str | dict[int, str] | Iterable = {0: "I"},
    ) -> PauliStringComplex | PauliStringSymbolic:
        """Creates an instance of the specified type.

        Supports:
        - PauliString(coeff, {index: 'X'})
        - PauliString(coeff, [('X', index)])
        - PauliString((coeff, [('X', index)]))   # OpenFermion-like tuple
        """
        if isinstance(coeff, tuple):
            if len(coeff) != 2:
                raise ValueError(
                    f"{self._name}: tuple input must be (coeff, terms)"
                )
            coeff, terms = coeff
            pauli_strings = self._convert_terms(terms)
        else:
            pauli_strings = self._convert_terms(pauli_strings)


        arg_type = type(coeff)
        builder = self._builders.get(arg_type)

        if not builder:
            raise ValueError(
                f"{self._name} of type {arg_type} has no builder registered"
            )

        if builder == PauliStringSymbolic:
            coeff = Expression(coeff)

        return builder(coeff, pauli_strings)

    def to_complex(self, expr: Expression) -> complex:
        return PauliStringSymbolic.to_complex(expr)


PauliString = PauliStringFactory(
    "PauliString",
    {
        complex: PauliStringComplex,
        float: PauliStringComplex,
        int: PauliStringComplex,
        str: PauliStringSymbolic,
    },
)