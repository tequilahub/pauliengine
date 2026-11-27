"""Sums of Pauli strings (qubit Hamiltonians): the ``PauliStringSum`` factory."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ._core import (
    Expression,
    PauliStringComplex,
    PauliStringSymbolic,
    PauliStringSumComplex,
    PauliStringSumSymbolic,
    SparsePauliStringComplex,
    SparsePauliStringSymbolic,
    SparsePauliStringSumComplex,
    SparsePauliStringSumSymbolic,
    unweighted_pauliString,
    _fits_sparse,
    unweighted_PauliStringSum,
)
from .sparse import span_and_weight, use_sparse

_SPARSE_STRINGS = (SparsePauliStringComplex, SparsePauliStringSymbolic)
_SPARSE_HAMILTONIANS = (SparsePauliStringSumComplex, SparsePauliStringSumSymbolic)

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
    if isinstance(term, (SparsePauliStringComplex, SparsePauliStringSymbolic)):
        coeff, ops = term.to_dictionary()
        if isinstance(term, SparsePauliStringComplex):
            coeff = _complex_coeff_to_expression(coeff)
        return PauliStringSymbolic(coeff, ops)
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


def _term_span_and_weight(term) -> tuple[int, int]:
    """(highest qubit + 1, weight) of a term of a list given to ``PauliStringSum``."""
    if isinstance(term, tuple):
        return span_and_weight(term[1])
    qubits = term.qubits()
    return (qubits[-1] + 1 if qubits else 0), len(qubits)


def _sparse_hamiltonian(terms, symbolic: bool):
    """A sparse Hamiltonian of the terms of a list given to ``PauliStringSum``."""
    if symbolic:
        strings = []
        for t in terms:
            if isinstance(t, SparsePauliStringSymbolic):
                strings.append(t)
            elif isinstance(t, SparsePauliStringComplex):
                strings.append(SparsePauliStringSymbolic(_complex_coeff_to_expression(t.coeff), t.to_dictionary()[1]))
            elif isinstance(t, tuple):
                c = Expression(t[0]) if isinstance(t[0], str) else _complex_coeff_to_expression(complex(t[0]))
                strings.append(SparsePauliStringSymbolic(c, t[1]))
            else:
                strings.append(SparsePauliStringSymbolic(_to_symbolic_ps(t)))
        return SparsePauliStringSumSymbolic(strings)
    strings = []
    for t in terms:
        if isinstance(t, SparsePauliStringComplex):
            strings.append(t)
        elif isinstance(t, tuple):
            strings.append(SparsePauliStringComplex(complex(t[0]), t[1]))
        else:
            strings.append(SparsePauliStringComplex(t))
    return SparsePauliStringSumComplex(strings)


def _weighted_term(term):
    """A term of a list that mixes unweighted strings with weighted terms, as a Pauli string."""
    if isinstance(term, unweighted_pauliString):
        return term.to_pauli_string()
    if isinstance(term, tuple) and not isinstance(term[0], str):
        return PauliStringComplex(complex(term[0]), term[1])
    return term


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


def _qh_to_openfermion(qh: PauliStringSumComplex | PauliStringSumSymbolic):
    """Convert a PauliEngine PauliStringSum into an openfermion.QubitOperator.

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
    is_symbolic = isinstance(qh, PauliStringSumSymbolic)
    for coeff, ops_dict in qh.to_dictionary():
        if is_symbolic:
            coeff = PauliStringSymbolic.to_complex(coeff)
        else:
            coeff = complex(coeff)
        key = tuple(sorted((int(q), str(p)) for q, p in ops_dict.items()))
        qop += QubitOperator(term=key, coefficient=coeff)
    return qop


def from_openfermion(qop: Any) -> PauliStringSumComplex:
    """Build a PauliEngine PauliStringSumComplex from an openfermion.QubitOperator.

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
        return PauliStringSumComplex.zero()
    return PauliStringSumComplex(ps_list)


class PauliStringSumFactory:
    """Factory behind ``pauliengine.PauliStringSum``.

    Builds a ``PauliStringSumComplex`` or, as soon as any term is symbolic, a
    ``PauliStringSumSymbolic``. When every string is light on many qubits it gives
    a ``SparsePauliStringSumComplex`` or ``SparsePauliStringSumSymbolic`` instead
    (see ``set_sparse_threshold``), whatever the form of the input; an existing sum
    or an OpenFermion operator is decided by the same rule.
    Accepted inputs:

    - ``list[PauliString]`` (complex, symbolic or mixed)
    - ``list[tuple[coeff, dict[int, str]]]`` (coeff a number or a symbolic string)
    - ``openfermion.QubitOperator``
    - an existing ``PauliStringSumComplex`` / ``PauliStringSumSymbolic`` (copied)

    A list of Pauli strings without a coefficient gives an
    ``unweighted_PauliStringSum``: an ordered set of operators with one phase each,
    meant as the basis of a Lie algebra. It is not a linear combination, so
    ``A + A == A``, and its ``commutator`` omits the factor 2. Next to terms with a
    coefficient, a string without one becomes a term of the Hamiltonian's type with
    the coefficient ``i**phase``.

    ``isinstance(x, PauliStringSum)`` is true for every weighted sum, complex or
    symbolic, with bit vectors or sparse; an ``unweighted_PauliStringSum`` is a set,
    not a sum, and is not one.
    """

    __slots__ = ["_builders", "_name"]

    _CLASSES = (PauliStringSumComplex, PauliStringSumSymbolic) + _SPARSE_HAMILTONIANS

    def __instancecheck__(self, obj: Any) -> bool:
        return isinstance(obj, self._CLASSES)

    def __init__(self, name: str, builders: dict[type, type]) -> None:
        """Create the factory; ``builders`` maps term types to Hamiltonian classes."""
        self._builders = builders
        self._name = name

    def __call__(self, pauli_strings: Any) -> PauliStringSumComplex | PauliStringSumSymbolic:
        """Create a qubit Hamiltonian from one of the accepted inputs.

        Args:
            pauli_strings: The terms, see the class docstring. An empty list
                gives the zero Hamiltonian.

        Returns:
            The compacted Hamiltonian; symbolic if any term is symbolic.

        Raises:
            ValueError: If the input format is not supported.
        """
        # OpenFermion QubitOperator: complex coefficients, then the same rule as a list.
        if _is_qubit_operator(pauli_strings):
            pauli_strings = _qubit_operator_to_ps_list(pauli_strings)
            if not pauli_strings:
                return PauliStringSumComplex.zero()

        # An existing sum: a copy, sparse or with bit vectors by the rule.
        if isinstance(pauli_strings, (PauliStringSumComplex, PauliStringSumSymbolic)):
            if pauli_strings.n_qubits() and _fits_sparse(pauli_strings):
                sparse = SparsePauliStringSumComplex if isinstance(pauli_strings, PauliStringSumComplex) \
                    else SparsePauliStringSumSymbolic
                return sparse(pauli_strings)
            return type(pauli_strings)(pauli_strings.paulistrings())
        if isinstance(pauli_strings, _SPARSE_HAMILTONIANS):
            if pauli_strings.n_qubits() and pauli_strings.fits_sparse():
                return pauli_strings.compact()
            return pauli_strings.to_dense()

        if isinstance(pauli_strings, unweighted_PauliStringSum):
            return unweighted_PauliStringSum(pauli_strings.terms())

        if not pauli_strings:
            # Empty list → zero Complex Hamiltonian
            return PauliStringSumComplex.zero()

        if all(isinstance(t, unweighted_pauliString) for t in pauli_strings):
            return unweighted_PauliStringSum(list(pauli_strings))

        # Next to weighted terms, an unweighted string gets the coefficient i**phase.
        if any(isinstance(t, unweighted_pauliString) for t in pauli_strings):
            pauli_strings = [_weighted_term(t) for t in pauli_strings]

        arg_type = _check_pauli_strings_type(pauli_strings)
        builder = self._builders.get(arg_type)
        if not builder:
            raise ValueError(
                f"{self._name} of type {arg_type} has no builder registered"
            )
        # Every string light on many qubits: store only the sites. The form of the
        # input (tuples, strings with bit vectors or sparse ones) does not matter.
        span, max_weight = 0, 0
        for t in pauli_strings:
            top, weight = _term_span_and_weight(t)
            span = max(span, top)
            max_weight = max(max_weight, weight)
        if use_sparse(span, max_weight):
            return _sparse_hamiltonian(pauli_strings, arg_type is PauliStringSymbolic)
        if arg_type is PauliStringSymbolic:
            return builder([_to_symbolic_ps(t) for t in pauli_strings])
        if any(isinstance(t, _SPARSE_STRINGS) for t in pauli_strings):
            return builder([t.to_dense() if isinstance(t, _SPARSE_STRINGS) else t for t in pauli_strings])
        return builder(pauli_strings)

    @staticmethod
    def from_openfermion(qop: Any) -> PauliStringSumComplex:
        """Build a ``PauliStringSumComplex`` from an ``openfermion.QubitOperator``."""
        return from_openfermion(qop)

    @staticmethod
    def to_openfermion(qh):
        """Convert a Hamiltonian into an ``openfermion.QubitOperator``.

        Symbolic coefficients must be fully numeric. Requires ``openfermion``.
        """
        return _qh_to_openfermion(qh)

    @staticmethod
    def from_matrix(matrix: Any, threshold: float = 0.0) -> PauliStringSumComplex:
        """Decompose a 2**n x 2**n matrix into Pauli strings.

        Uses the method of Georges et al. (2025), O(4**n n). Qubit 0 is the
        most significant bit of the row and column index, matching
        ``to_matrix(ignore_unused_qubits=False)``, so
        ``PauliStringSum.from_matrix(H.to_matrix(False)) == H``.

        Args:
            matrix: A square matrix with side 2**n, as a numpy array or a list
                of rows. It is not modified.
            threshold: Terms with ``abs(coefficient) <= threshold`` are dropped.
                The default 0 drops only exact zeros.

        Returns:
            The matrix as a ``PauliStringSumComplex``.

        Raises:
            ValueError: If the matrix is not square with a power-of-two side.
        """
        try:
            import numpy as np
        except ImportError:
            return PauliStringSumComplex.from_matrix(matrix, threshold)
        arr = np.ascontiguousarray(matrix, dtype=np.complex128)
        if arr.ndim != 2:
            raise ValueError(f"from_matrix: expected a 2-D matrix, got {arr.ndim} dimensions")
        return PauliStringSumComplex.from_matrix(arr, threshold)

    @staticmethod
    def zero() -> PauliStringSumComplex:
        """The zero Hamiltonian (no terms)."""
        return PauliStringSumComplex.zero()

    @staticmethod
    def unit() -> PauliStringSumComplex:
        """The identity Hamiltonian (one term, coefficient 1, no operators)."""
        return PauliStringSumComplex.unit()


def _check_pauli_strings_type(
    pauli_string: (
        Sequence[PauliStringComplex | PauliStringSymbolic]
        | Sequence[tuple[complex | str, dict[int, str]]]
    ),
) -> type:
    # If ANY term is symbolic the whole Hamiltonian must be symbolic.
    for term in pauli_string:
        if isinstance(term, (PauliStringSymbolic, SparsePauliStringSymbolic)):
            return PauliStringSymbolic
        if isinstance(term, tuple) and isinstance(term[0], str):
            return PauliStringSymbolic

    first_term = pauli_string[0]
    if isinstance(first_term, (PauliStringComplex, SparsePauliStringComplex)):
        return PauliStringComplex
    if isinstance(first_term, tuple):
        coeff = first_term[0]
        if isinstance(coeff, (complex, float, int)):
            return complex
    raise ValueError("Invalid input.")


PauliStringSum = PauliStringSumFactory(
    "PauliStringSum",
    {
        PauliStringComplex: PauliStringSumComplex,
        complex: PauliStringSumComplex,
        PauliStringSymbolic: PauliStringSumSymbolic,
    },
)

# Old names, kept so that code written for ``QubitHamiltonian`` keeps working.
QubitHamiltonian = PauliStringSum
QubitHamiltonianFactory = PauliStringSumFactory
