"""PauliEngine library."""

from ._core import (
    Expression,
    __build_type__,
    __compiler_flags__,
    __openmp__,
    __omp_max_threads__,
    __omp_version__,
    is_symbol_complex,
    is_symbol_real,
    set_symbol_complex,
    set_symbol_real,
)
from ._version import __version__
from .pauli_cycle import PauliCycle, PauliCycleSum
from .pauli_orbit import PauliOrbit, PauliOrbitSum
from .pauli_string import PauliString
from .pauli_string_sum import PauliStringSum, QubitHamiltonian, from_openfermion
from .sparse import get_sparse_threshold, set_sparse_threshold


def inner_product(a, b):
    """Normalized Hilbert-Schmidt inner product Tr(a^dagger b) / 2**n.

    Works for any two of PauliString, PauliStringSum, PauliCycle, PauliCycleSum,
    PauliOrbit and PauliOrbitSum; Pauli strings are orthonormal. Conjugate-linear
    in ``a``. Multiply by 2**n for the unnormalized trace on n qubits.
    """
    return a.inner_product(b)

__all__ = [
    "Expression",
    "PauliCycle",
    "PauliCycleSum",
    "PauliOrbit",
    "PauliOrbitSum",
    "PauliString",
    "PauliStringSum",
    "QubitHamiltonian",
    "__build_type__",
    "__compiler_flags__",
    "__omp_max_threads__",
    "__omp_version__",
    "__openmp__",
    "__version__",
    "from_openfermion",
    "get_sparse_threshold",
    "inner_product",
    "is_symbol_complex",
    "is_symbol_real",
    "set_sparse_threshold",
    "set_symbol_complex",
    "set_symbol_real",
]
