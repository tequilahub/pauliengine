"""Module for PauliString classes and factory."""

from __future__ import annotations

from ._core import (
    Expression,
    PauliStringComplex,
    PauliStringSymbolic,
)


class PauliStringFactory:
    """An implementation of the Factory method pattern."""

    __slots__ = ["_builders", "_name"]

    def __init__(self, name: str, builders: dict[type, type]) -> None:
        """Initializes the PauliStringFactory."""
        self._builders = builders
        self._name = name

    def __call__(
        self,
        coeff: complex | None = 0,
        pauli_strings: str | dict[int, str] | None = {0: "I"},
    ) -> PauliStringComplex | PauliStringSymbolic:
        """Creates an instance of the specified type."""
        arg_type = type(coeff)
        builder = self._builders.get(arg_type)
        if not builder:
            raise ValueError(
                f"{self._name} of type {arg_type} has no builder registered"
            )
        if builder == PauliStringSymbolic:
            coeff = Expression(coeff)
        return builder(coeff, pauli_strings)


PauliString = PauliStringFactory(
    "PauliString",
    {
        complex: PauliStringComplex,
        float: PauliStringComplex,
        int: PauliStringComplex,
        str: PauliStringSymbolic,
    },
)
