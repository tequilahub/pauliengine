"""Module for wrapping QubitHamiltonian."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ._core import (
    Expression,
    PauliStringComplex,
    PauliStringSymbolic,
    QubitHamiltonianComplex,
    QubitHamiltonianSymbolic,
)

if TYPE_CHECKING:
    from collections.abc import Sequence


class QubitHamiltonianFactory:
    """An implementation of the Factory method pattern."""

    __slots__ = ["_builders", "_name"]

    def __init__(self, name: str, builders: dict[type, type]) -> None:
        """Initializes the QubitHamiltonianFactory."""
        self._builders = builders
        self._name = name

    def __call__(
        self,
        pauli_strings: (
            Sequence[PauliStringComplex | PauliStringSymbolic]
            | Sequence[tuple[complex | str, dict[int, str]]]
        ),
    ) -> QubitHamiltonianComplex | QubitHamiltonianSymbolic:
        """Creates an instance of the specified type."""
        arg_type = _check_pauli_strings_type(pauli_strings)
        builder = self._builders.get(arg_type)
        if not builder:
            raise ValueError(
                f"{self._name} of type {arg_type} has no builder registered"
            )
        if arg_type is str:
            new_ps = [
                (Expression(coeff), pauli_dict) for coeff, pauli_dict in pauli_strings
            ]
            return builder(new_ps)
        return builder(pauli_strings)


def _check_pauli_strings_type(
    pauli_string: (
        Sequence[PauliStringComplex | PauliStringSymbolic]
        | Sequence[tuple[complex | str, dict[int, str]]]
    ),
) -> PauliStringSymbolic | PauliStringComplex:
    first_term = pauli_string[0]
    # TODO: Maybe check consistence of terms in pauli_string, just in case
    if isinstance(first_term, PauliStringComplex):
        return PauliStringComplex
    if isinstance(first_term, PauliStringSymbolic):
        return PauliStringSymbolic
    if isinstance(first_term, tuple):
        # First term in the tuple should be the coeffient
        coeff = first_term[0]
        if isinstance(coeff, complex):
            return complex
        if isinstance(coeff, str):
            return PauliStringSymbolic
    raise ValueError("Invalid input.")


QubitHamiltonian = QubitHamiltonianFactory(
    "QubitHamiltonian",
    {
        PauliStringComplex: QubitHamiltonianComplex,
        complex: QubitHamiltonianComplex,
        PauliStringSymbolic: QubitHamiltonianSymbolic,
        str: QubitHamiltonianSymbolic,
    },
)
