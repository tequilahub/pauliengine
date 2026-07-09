"""PauliEngine library."""

from ._core import (
    __build_type__,
    __compiler_flags__,
    __openmp__,
    __omp_max_threads__,
    __omp_version__,
)
from ._version import __version__
from .pauli_string import PauliString
from .qubit_hamiltonian import QubitHamiltonian, from_openfermion

__all__ = [
    "PauliString",
    "QubitHamiltonian",
    "__build_type__",
    "__compiler_flags__",
    "__omp_max_threads__",
    "__omp_version__",
    "__openmp__",
    "__version__",
    "from_openfermion",
]
