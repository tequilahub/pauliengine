"""Module for wrapping QubitHamiltonian."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ._core import (
    Expression,
    PauliStringComplex,
    PauliStringSymbolic,
    QubitHamiltonianComplex,
    QubitHamiltonianSymbolic,
)

if TYPE_CHECKING:
    from collections.abc import Sequence


def _complex_coeff_to_expression(coeff: complex) -> Expression:
    """Convert a complex<double> coefficient to a SymEngine Expression."""
    real, imag = coeff.real, coeff.imag
    if imag == 0.0:
        return Expression(repr(real))
    return Expression(f"{real!r} + {imag!r}*I")


def _to_symbolic_ps(term) -> PauliStringSymbolic:
    """Convert any single term to PauliStringSymbolic."""
    if isinstance(term, PauliStringSymbolic):
        return term
    if isinstance(term, PauliStringComplex):
        coeff, ops = term.to_dictionary()
        return PauliStringSymbolic(_complex_coeff_to_expression(coeff), ops)
    if isinstance(term, tuple):
        coeff, ops = term
        if isinstance(coeff, str):
            return PauliStringSymbolic(Expression(coeff), ops)
        # complex or float coeff in a tuple
        c = complex(coeff)
        return PauliStringSymbolic(_complex_coeff_to_expression(c), ops)
    raise ValueError(f"Cannot convert {type(term)} to PauliStringSymbolic")


def _is_qubit_operator(obj: Any) -> bool:
    """Duck-typed check for openfermion.QubitOperator without importing it."""
    return (
        hasattr(obj, "terms")
        and isinstance(getattr(obj, "terms", None), dict)
        and type(obj).__name__ == "QubitOperator"
    )


def _qubit_operator_to_ps_list(qop: Any) -> list[PauliStringComplex]:
    """Convert openfermion.QubitOperator.terms to a list of PauliStringComplex."""
    ps_list: list[PauliStringComplex] = []
    for key, coeff in qop.terms.items():
        ops_dict = {int(q): str(p).upper() for (q, p) in key}
        ps_list.append(PauliStringComplex(complex(coeff), ops_dict))
    return ps_list


def _qh_to_openfermion(qh: QubitHamiltonianComplex | QubitHamiltonianSymbolic):
    """Convert a PauliEngine QubitHamiltonian into an openfermion.QubitOperator.

    For symbolic Hamiltonians coefficients are first evaluated to complex
    via PauliStringSymbolic.to_complex, which only succeeds when every
    coefficient is numeric.
    """
    try:
        from openfermion import QubitOperator
    except ImportError as e:
        raise ImportError(
            "openfermion is required for to_openfermion(). "
            "Install it via `pip install openfermion`."
        ) from e

    qop = QubitOperator()
    is_symbolic = isinstance(qh, QubitHamiltonianSymbolic)
    for coeff, ops_dict in qh.to_dictionary():
        if is_symbolic:
            coeff = PauliStringSymbolic.to_complex(coeff)
        else:
            coeff = complex(coeff)
        key = tuple(sorted((int(q), str(p)) for q, p in ops_dict.items()))
        qop += QubitOperator(term=key, coefficient=coeff)
    return qop


def from_openfermion(qop: Any) -> QubitHamiltonianComplex:
    """Build a PauliEngine QubitHamiltonianComplex from an openfermion.QubitOperator.

    Raises:
        TypeError: If ``qop`` is not a ``QubitOperator``.
    """
    if not _is_qubit_operator(qop):
        raise TypeError(
            "from_openfermion expects an openfermion.QubitOperator instance, "
            f"got {type(qop).__name__}"
        )
    ps_list = _qubit_operator_to_ps_list(qop)
    if not ps_list:
        return QubitHamiltonianComplex.zero()
    return QubitHamiltonianComplex(ps_list)


class QubitHamiltonianFactory:
    """Factory behind ``pauliengine.QubitHamiltonian``.

    Builds a ``QubitHamiltonianComplex`` or, as soon as any term is symbolic, a
    ``QubitHamiltonianSymbolic``. Accepted inputs:

    - ``list[PauliString]`` (complex, symbolic or mixed)
    - ``list[tuple[coeff, dict[int, str]]]`` (coeff a number or a symbolic string)
    - ``openfermion.QubitOperator``
    - an existing ``QubitHamiltonianComplex`` / ``QubitHamiltonianSymbolic`` (copied)
    """

    __slots__ = ["_builders", "_name"]

    def __init__(self, name: str, builders: dict[type, type]) -> None:
        """Create the factory; ``builders`` maps term types to Hamiltonian classes."""
        self._builders = builders
        self._name = name

    def __call__(self, pauli_strings: Any) -> QubitHamiltonianComplex | QubitHamiltonianSymbolic:
        """Create a qubit Hamiltonian from one of the accepted inputs.

        Args:
            pauli_strings: The terms, see the class docstring. An empty list
                gives the zero Hamiltonian.

        Returns:
            The compacted Hamiltonian; symbolic if any term is symbolic.

        Raises:
            ValueError: If the input format is not supported.
        """
        # OpenFermion QubitOperator → always Complex
        if _is_qubit_operator(pauli_strings):
            ps_list = _qubit_operator_to_ps_list(pauli_strings)
            if not ps_list:
                return QubitHamiltonianComplex.zero()
            return QubitHamiltonianComplex(ps_list)

        # Already a QubitHamiltonian → pass through (copy via reconstruction)
        if isinstance(pauli_strings, (QubitHamiltonianComplex, QubitHamiltonianSymbolic)):
            builder = type(pauli_strings)
            return builder(pauli_strings.paulistrings())

        if not pauli_strings:
            # Empty list → zero Complex Hamiltonian
            return QubitHamiltonianComplex.zero()

        arg_type = _check_pauli_strings_type(pauli_strings)
        builder = self._builders.get(arg_type)
        if not builder:
            raise ValueError(
                f"{self._name} of type {arg_type} has no builder registered"
            )
        if arg_type is PauliStringSymbolic:
            return builder([_to_symbolic_ps(t) for t in pauli_strings])
        return builder(pauli_strings)

    @staticmethod
    def from_openfermion(qop: Any) -> QubitHamiltonianComplex:
        """Build a ``QubitHamiltonianComplex`` from an ``openfermion.QubitOperator``."""
        return from_openfermion(qop)

    @staticmethod
    def to_openfermion(qh):
        """Convert a Hamiltonian into an ``openfermion.QubitOperator``.

        Symbolic coefficients must be fully numeric. Requires ``openfermion``.
        """
        return _qh_to_openfermion(qh)

    @staticmethod
    def zero() -> QubitHamiltonianComplex:
        """The zero Hamiltonian (no terms)."""
        return QubitHamiltonianComplex.zero()

    @staticmethod
    def unit() -> QubitHamiltonianComplex:
        """The identity Hamiltonian (one term, coefficient 1, no operators)."""
        return QubitHamiltonianComplex.unit()


def _check_pauli_strings_type(
    pauli_string: (
        Sequence[PauliStringComplex | PauliStringSymbolic]
        | Sequence[tuple[complex | str, dict[int, str]]]
    ),
) -> type:
    # If ANY term is symbolic the whole Hamiltonian must be symbolic.
    for term in pauli_string:
        if isinstance(term, PauliStringSymbolic):
            return PauliStringSymbolic
        if isinstance(term, tuple) and isinstance(term[0], str):
            return PauliStringSymbolic

    first_term = pauli_string[0]
    if isinstance(first_term, PauliStringComplex):
        return PauliStringComplex
    if isinstance(first_term, tuple):
        coeff = first_term[0]
        if isinstance(coeff, (complex, float, int)):
            return complex
    raise ValueError("Invalid input.")


QubitHamiltonian = QubitHamiltonianFactory(
    "QubitHamiltonian",
    {
        PauliStringComplex: QubitHamiltonianComplex,
        complex: QubitHamiltonianComplex,
        PauliStringSymbolic: QubitHamiltonianSymbolic,
    },
)
