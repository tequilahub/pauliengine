"""Random Hamiltonian generation for benchmarks.

Two knobs matter for scaling:

* ``n_terms``  — how many Pauli strings the Hamiltonian is a sum of.
* ``n_qubits`` — the length (weight) of each individual Pauli string.

We deliberately use a fixed seed by default so successive runs are
comparable across machines and versions.
"""

from __future__ import annotations

import random
from typing import Literal

import pauliengine as pe

CoeffKind = Literal["complex", "symbolic"]

_PAULI_OPS = ("X", "Y", "Z")


def _random_pauli_term_complex(
    rng: random.Random, n_qubits: int,
) -> tuple[complex, dict[int, str]]:
    ops = {q: rng.choice(_PAULI_OPS) for q in range(n_qubits)}
    coeff = complex(rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0))
    return coeff, ops


def _random_pauli_term_symbolic(
    rng: random.Random, n_qubits: int, symbol_pool: int,
) -> tuple[str, dict[int, str]]:
    ops = {q: rng.choice(_PAULI_OPS) for q in range(n_qubits)}
    # Coefficient is a symbolic literal like "a3" — the pool controls how many
    # distinct symbols show up across the Hamiltonian.
    coeff = f"a{rng.randrange(symbol_pool)}"
    return coeff, ops


def random_hamiltonian(
    n_terms: int,
    n_qubits: int,
    coeff_kind: CoeffKind,
    seed: int = 0,
    symbol_pool: int = 8,
):
    """Build a random Hamiltonian with ``n_terms`` strings of length ``n_qubits``.

    Args:
        n_terms: number of Pauli strings summed in the Hamiltonian.
        n_qubits: number of non-identity operators in each string.
        coeff_kind: ``"complex"`` for numeric coefficients, ``"symbolic"``
            for SymEngine-expression coefficients.
        seed: RNG seed. Same seed → same Hamiltonian, on any machine.
        symbol_pool: only used for ``"symbolic"``. Controls how many distinct
            symbols appear across terms.

    Returns:
        A ``pauliengine.QubitHamiltonian`` (Complex or Symbolic dispatch).
    """
    if n_terms <= 0:
        raise ValueError("n_terms must be positive")
    if n_qubits <= 0:
        raise ValueError("n_qubits must be positive")

    rng = random.Random(seed)
    terms = []
    for _ in range(n_terms):
        if coeff_kind == "complex":
            terms.append(_random_pauli_term_complex(rng, n_qubits))
        elif coeff_kind == "symbolic":
            terms.append(_random_pauli_term_symbolic(rng, n_qubits, symbol_pool))
        else:
            raise ValueError(f"unknown coeff_kind: {coeff_kind!r}")
    return pe.QubitHamiltonian(terms)
