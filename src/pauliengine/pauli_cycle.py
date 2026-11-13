"""Module for PauliCycle and PauliCycleSum classes."""

from __future__ import annotations

from typing import Any

from ._core import (
    Expression,
    PauliCycleComplex,
    PauliCycleSumComplex,
    PauliCycleSumSymbolic,
    PauliCycleSymbolic,
    PauliStringSymbolic,
)
from .qubit_hamiltonian import _complex_coeff_to_expression, _to_symbolic_ps


def _is_symbolic_coeff(coeff: Any) -> bool:
    """A string such as ``"a*b"`` or an ``Expression`` makes a symbolic coefficient."""
    return isinstance(coeff, (str, Expression))


def _as_expression(coeff: Any) -> Expression:
    """A coefficient as a SymEngine Expression."""
    if isinstance(coeff, Expression):
        return coeff
    if isinstance(coeff, str):
        return Expression(coeff)
    return _complex_coeff_to_expression(complex(coeff))


def _as_symbolic_cycle(cycle):
    """A complex cycle as a symbolic one with the same base and coefficient."""
    if isinstance(cycle, PauliCycleSymbolic):
        return cycle
    return PauliCycleSymbolic(cycle.n_qubits, _to_symbolic_ps(cycle.base))


class _CoefficientFactory:
    """Builds the complex or the symbolic class; ``isinstance`` accepts either."""

    def __init__(self, name: str, classes: tuple[type, ...]) -> None:
        self._name = name
        self._classes = classes

    def __instancecheck__(self, obj: Any) -> bool:
        return isinstance(obj, self._classes)

    def __repr__(self) -> str:
        return f"<pauliengine.{self._name} factory>"


class PauliCycleFactory(_CoefficientFactory):
    """Factory behind ``pauliengine.PauliCycle``.

    Builds a ``PauliCycleComplex`` or, for a symbolic coefficient (a string such as
    ``"a*b"`` or an ``Expression``) or a symbolic base string, a ``PauliCycleSymbolic``::

        PauliCycle(4, "XZ")                     # coefficient 1
        PauliCycle(4, {0: "X", 1: "Z"}, 0.5)
        PauliCycle(4, "XZ", "a")                # symbolic coefficient a
        PauliCycle(4, PauliString("a", "XZ"))   # base with its coefficient
    """

    def __init__(self) -> None:
        super().__init__("PauliCycle", (PauliCycleComplex, PauliCycleSymbolic))

    def __call__(self, n_qubits: int, *args: Any, **kwargs: Any):
        base = args[0] if args else kwargs.get("data", kwargs.get("pauli_string", kwargs.get("base")))
        coeff = args[1] if len(args) > 1 else kwargs.get("coeff")
        if isinstance(base, PauliStringSymbolic) or _is_symbolic_coeff(coeff):
            if coeff is not None:
                if len(args) > 1:
                    args = (args[0], _as_expression(coeff), *args[2:])
                else:
                    kwargs["coeff"] = _as_expression(coeff)
            return PauliCycleSymbolic(n_qubits, *args, **kwargs)
        return PauliCycleComplex(n_qubits, *args, **kwargs)

    @staticmethod
    def from_hash(hash: int, coeff: Any = 1.0):
        """The cycle whose ``get_hash()`` equals ``hash``, with the given coefficient."""
        if _is_symbolic_coeff(coeff):
            return PauliCycleSymbolic.from_hash(hash, _as_expression(coeff))
        return PauliCycleComplex.from_hash(hash, complex(coeff))


class PauliCycleSumFactory(_CoefficientFactory):
    """Factory behind ``pauliengine.PauliCycleSum``.

    Builds a ``PauliCycleSumComplex`` from complex cycles, or a ``PauliCycleSumSymbolic``
    as soon as one cycle is symbolic (the complex ones are converted).
    """

    def __init__(self) -> None:
        super().__init__("PauliCycleSum", (PauliCycleSumComplex, PauliCycleSumSymbolic))

    def __call__(self, data=()):
        cycles = list(data)
        if any(isinstance(c, PauliCycleSymbolic) for c in cycles):
            return PauliCycleSumSymbolic([_as_symbolic_cycle(c) for c in cycles])
        return PauliCycleSumComplex(cycles)


PauliCycle = PauliCycleFactory()
PauliCycleSum = PauliCycleSumFactory()

__all__ = ["PauliCycle", "PauliCycleSum"]
