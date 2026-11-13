"""Module for PauliOrbit and PauliOrbitSum classes."""

from __future__ import annotations

from typing import Any

from ._core import (
    PauliOrbit,
    PauliOrbitSumComplex,
    PauliOrbitSumSymbolic,
)
from .pauli_cycle import _as_expression, _CoefficientFactory, _is_symbolic_coeff


class PauliOrbitSumFactory(_CoefficientFactory):
    """Factory behind ``pauliengine.PauliOrbitSum``.

    Builds an empty or filled sum of orbits on ``n_qubits`` qubits::

        PauliOrbitSum(4)                                     # empty, complex coefficients
        PauliOrbitSum(4, symbolic=True)                      # empty, symbolic coefficients
        PauliOrbitSum(4, [(PauliOrbit(4, 1, 0, 0), 0.5)])    # complex
        PauliOrbitSum(4, [(PauliOrbit(4, 1, 0, 0), "a")])    # symbolic

    The sum is symbolic if ``symbolic`` is true or any coefficient is a string such
    as ``"a*b"`` or an ``Expression``. Terms with the same label are merged.
    """

    def __init__(self) -> None:
        super().__init__("PauliOrbitSum", (PauliOrbitSumComplex, PauliOrbitSumSymbolic))

    def __call__(self, n_qubits: int, terms=(), *, symbolic: bool | None = None):
        terms = list(terms)
        if symbolic is None:
            symbolic = any(_is_symbolic_coeff(coeff) for _, coeff in terms)
        if symbolic:
            out = PauliOrbitSumSymbolic(n_qubits)
            for orbit, coeff in terms:
                out.add(orbit, _as_expression(coeff))
        else:
            out = PauliOrbitSumComplex(n_qubits)
            for orbit, coeff in terms:
                out.add(orbit, complex(coeff))
        return out


PauliOrbitSum = PauliOrbitSumFactory()

__all__ = ["PauliOrbit", "PauliOrbitSum"]
