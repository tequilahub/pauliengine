"""When the factories build sparse Pauli strings and Hamiltonians.

``PauliString`` and ``PauliStringSum`` store operators as bit vectors up to the
highest qubit they act on. For low-weight strings on many qubits that wastes memory
and time: a string of weight 2 on qubit 4000 takes 126 words, and every product
walks all of them. The sparse classes (``SparsePauliStringComplex``,
``SparsePauliStringSumComplex`` and their symbolic variants) store only the
non-identity sites, so their cost grows with the weight instead.

The factories pick the sparse classes when the operator spans at least
``min_qubits`` qubits (highest qubit index + 1) and every one of its strings is no
heavier than ``max_density`` times that span; a single heavy string keeps the whole
operator with bit vectors. Only the content decides, not the form of the input
(tuples, strings with bit vectors or sparse ones, an existing sum). Sums, differences and scaling keep a
sparse operand sparse; products and commutators always return the classes with bit
vectors, since they make strings heavier.
"""

from __future__ import annotations

from ._core import _set_sparse_threshold, _sparse_enters, _sparse_threshold


def set_sparse_threshold(min_qubits: float | None = None, max_density: float | None = None) -> None:
    """Change when operators are stored sparse.

    The factories ``PauliString`` and ``PauliStringSum`` build a sparse operator
    when it spans at least ``min_qubits`` qubits and each of its strings has a
    weight of at most ``max_density`` times the span. Products and commutators of sparse operators
    always come back with bit vectors. A sum of a sparse operand with an operator
    with bit vectors that misses this rule is formed with bit vectors as well, so
    chains of operations do not switch back and forth.

    Args:
        min_qubits: Smallest span (highest qubit index + 1) for a sparse operator.
            ``0`` makes every operator within the density limit sparse,
            ``float("inf")`` turns sparse operators off. Default 256.
        max_density: Largest weight of any string, as a fraction of the
            span. Default 1/256, so strings of weight 2 are sparse from 512 qubits on
            and strings of weight 4 from 1024.
    """
    current_min, current_density = _sparse_threshold()
    if min_qubits is not None and min_qubits < 0:
        raise ValueError("set_sparse_threshold: min_qubits must be >= 0")
    if max_density is not None and not 0 <= max_density <= 1:
        raise ValueError("set_sparse_threshold: max_density must be within [0, 1]")
    _set_sparse_threshold(current_min if min_qubits is None else float(min_qubits),
                          current_density if max_density is None else float(max_density))


def get_sparse_threshold() -> dict[str, float]:
    """The current threshold as ``{"min_qubits": ..., "max_density": ...}``."""
    min_qubits, max_density = _sparse_threshold()
    return {"min_qubits": min_qubits, "max_density": max_density}


def use_sparse(span: int, max_weight: int) -> bool:
    """Whether an operator with this span and heaviest string is built sparse."""
    return _sparse_enters(float(span), float(max_weight))


def span_and_weight(ops) -> tuple[int, int]:
    """(highest qubit + 1, number of non-identity operators) of a Pauli string given
    as a dict {qubit: 'X'}, a dense text 'XIZ' or a sparse text 'X0 Z2'."""
    if isinstance(ops, dict):
        top, weight = 0, 0
        for q, p in ops.items():
            if p not in ("I", "i"):
                weight += 1
                if q + 1 > top:
                    top = q + 1
        return top, weight
    if isinstance(ops, str):
        tokens = ops.split()
        if len(tokens) > 1 or (tokens and tokens[0][1:].isdigit()):
            top, weight = 0, 0
            for token in tokens:
                if token[0] not in "Ii" and token[1:].isdigit():
                    weight += 1
                    top = max(top, int(token[1:]) + 1)
            return top, weight
        weight = sum(1 for ch in ops if ch not in "Ii")
        stripped = ops.rstrip("Ii")
        return len(stripped), weight
    raise TypeError(f"unsupported Pauli operator format: {type(ops)}")
