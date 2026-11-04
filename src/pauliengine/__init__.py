"""PauliEngine library."""

from ._core import (
    __build_type__,
    __compiler_flags__,
    __openmp__,
    __omp_max_threads__,
    __omp_version__,
)
from ._version import __version__
from .pauli_cycle import PauliCycle, PauliCycleSum
from .pauli_orbit import PauliOrbit, PauliOrbitSum
from .pauli_string import PauliString
from .qubit_hamiltonian import QubitHamiltonian, from_openfermion


def inner_product(a, b):
    """Normalized Hilbert-Schmidt inner product Tr(a^dagger b) / 2**n.

    Works for any two of PauliString, QubitHamiltonian, PauliCycle, PauliCycleSum,
    PauliOrbit and PauliOrbitSum; Pauli strings are orthonormal. Conjugate-linear
    in ``a``. Multiply by 2**n for the unnormalized trace on n qubits.
    """
    return a.inner_product(b)

__all__ = [
    "PauliCycle",
    "PauliCycleSum",
    "PauliOrbit",
    "PauliOrbitSum",
    "PauliString",
    "QubitHamiltonian",
    "__build_type__",
    "__compiler_flags__",
    "__omp_max_threads__",
    "__omp_version__",
    "__openmp__",
    "__version__",
    "from_openfermion",
    "inner_product",
]
