"""Module for PauliString classes and factory."""

from __future__ import annotations

from collections.abc import Iterable
from ._core import (
    Expression,
    PauliStringComplex,
    PauliStringSymbolic,
    unweighted_pauliString,
)


_UNSET = object()


class PauliStringFactory:
    """Factory behind ``pauliengine.PauliString``.

    Builds a ``PauliStringComplex`` for a numeric coefficient, a
    ``PauliStringSymbolic`` for a string coefficient and an
    ``unweighted_pauliString`` without a coefficient, and accepts several
    formats for the Pauli operators.
    """

    __slots__ = ["_builders", "_name"]

    def __init__(self, name: str, builders: dict[type, type]) -> None:
        """Create the factory; ``builders`` maps coefficient types to classes."""
        self._builders = builders
        self._name = name

    def _convert_terms(self, terms) -> dict[int, str] | str:
        """Convert various input formats into a text or an {index: 'Pauli'} dict."""
        if isinstance(terms, (dict, str)):
            return terms

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
        coeff: complex | str | tuple | object = _UNSET,
        pauli_strings: str | dict[int, str] | Iterable | object = _UNSET,
        *,
        phase: int | None = None,
    ) -> PauliStringComplex | PauliStringSymbolic | unweighted_pauliString:
        """Create a Pauli string.

        Supported forms::

            PauliString("XIZ")                        # no coefficient: unweighted
            PauliString("X0 Z2", phase=1)             # unweighted, i**1 * X0 Z2
            PauliString({0: "X", 2: "Z"})             # unweighted
            PauliString([("X", 0), ("Z", 2)])         # unweighted
            PauliString(weighted_string)              # unweighted, coefficient dropped
            PauliString(coeff, "XIZ")                 # dense string
            PauliString(coeff, "X0 Z2")               # space-separated string
            PauliString(coeff, {0: "X", 2: "Z"})      # dict {qubit: Pauli}
            PauliString(coeff, [("X", 0), ("Z", 2)])  # list of (Pauli, qubit)
            PauliString((coeff, [("X", 0)]))          # OpenFermion-style tuple
            PauliString(coeff)                        # coeff times the identity

        A single string is always read as Pauli operators. For a symbolic multiple
        of the identity write ``PauliString("a", {})``.

        Args:
            coeff: The coefficient. ``int``, ``float`` or ``complex`` give a
                ``PauliStringComplex``; a ``str`` such as ``"a*b"`` is parsed as
                a symbolic expression and gives a ``PauliStringSymbolic``.
                May also be a ``(coeff, terms)`` tuple, then ``pauli_strings``
                is ignored. Without a coefficient, the first argument holds the
                Pauli operators and the result is unweighted.
            pauli_strings: The Pauli operators in one of the forms above.
            phase: Only without a coefficient: the exponent ``k`` of the phase
                ``i**k`` of the unweighted string. Defaults to 0.

        Returns:
            The Pauli string: ``PauliStringComplex``, ``PauliStringSymbolic`` or,
            without a coefficient, ``unweighted_pauliString``.

        Raises:
            ValueError: If the coefficient type or the operator format is not
                supported, or ``phase`` is given together with a coefficient.
        """
        if pauli_strings is _UNSET and not isinstance(coeff, (bool, int, float, complex, tuple)):
            terms = {} if coeff is _UNSET else coeff
            return unweighted_pauliString(self._unweighted_terms(terms), phase or 0)
        if coeff is _UNSET:
            return unweighted_pauliString(self._unweighted_terms(pauli_strings), phase or 0)
        if phase is not None:
            raise ValueError(f"{self._name}: phase is only for a string without a coefficient")

        if isinstance(coeff, tuple):
            if len(coeff) != 2:
                raise ValueError(
                    f"{self._name}: tuple input must be (coeff, terms)"
                )
            coeff, terms = coeff
            pauli_strings = self._convert_terms(terms)
        elif pauli_strings is _UNSET:
            pauli_strings = {}
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

    def _unweighted_terms(self, terms):
        """The operators of an unweighted string, as its constructor takes them."""
        if isinstance(terms, (str, dict, PauliStringComplex, PauliStringSymbolic)):
            return terms
        if isinstance(terms, Iterable):
            return list(terms)
        raise ValueError(
            f"{self._name}: unsupported pauli_strings format: {type(terms)}"
        )

    def from_hash(
        self, hash: int, coeff: complex | str = 1
    ) -> PauliStringComplex | PauliStringSymbolic:
        """Inverse of ``get_hash()``: the Pauli string with that operator part.

        The coefficient type selects the class exactly as in ``__call__``.
        """
        builder = self._builders.get(type(coeff))
        if not builder:
            raise ValueError(
                f"{self._name} of type {type(coeff)} has no builder registered"
            )
        if builder == PauliStringSymbolic:
            coeff = Expression(coeff)
        return builder.from_hash(hash, coeff)

    def to_complex(self, expr: Expression) -> complex:
        """Evaluate a fully numeric symbolic coefficient to a complex number.

        Raises:
            RuntimeError: If ``expr`` still contains free symbols.
        """
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