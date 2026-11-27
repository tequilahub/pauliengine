Python API
==========

.. currentmodule:: pauliengine

The factories ``PauliString`` and ``PauliStringSum`` pick the numeric or the
symbolic class from the coefficient types, and for low-weight operators on many
qubits the sparse class (``SparsePauliStringComplex``, ``SparsePauliStringSumComplex``
and their symbolic variants; see ``set_sparse_threshold``); the classes below are
what they return. ``isinstance(x, pauliengine.PauliString)`` and
``isinstance(x, pauliengine.PauliStringSum)`` accept all of them. The sparse
classes have the same methods as the ones with bit vectors, plus ``to_dense()``,
``fits_sparse()`` and, for strings, ``sites``.

Pauli strings
-------------

.. autoclass:: pauliengine.pauli_string.PauliStringFactory
   :members: __call__, from_hash, to_complex

.. autoclass:: pauliengine._core.PauliStringComplex

.. autoclass:: pauliengine._core.PauliStringSymbolic

.. autoclass:: pauliengine._core.unweighted_pauliString

Pauli string sums (qubit Hamiltonians)
--------------------------------------

``PauliStringSum`` was called ``QubitHamiltonian`` before; the old names
(``QubitHamiltonian``, ``QubitHamiltonianComplex``, ``QubitHamiltonianSymbolic``,
``unweighted_QubitHamiltonian``, ``to_qubit_hamiltonian()``, the module
``pauliengine.qubit_hamiltonian``) remain as aliases.

.. autoclass:: pauliengine.pauli_string_sum.PauliStringSumFactory
   :members: __call__, from_matrix, from_openfermion, to_openfermion, zero, unit

.. autofunction:: pauliengine.from_openfermion

.. autoclass:: pauliengine._core.PauliStringSumComplex

.. autoclass:: pauliengine._core.PauliStringSumSymbolic

.. autoclass:: pauliengine._core.unweighted_PauliStringSum

Pauli cycles
------------

``pauliengine.PauliCycle`` and ``pauliengine.PauliCycleSum`` build the complex or,
for a symbolic coefficient such as ``"a*b"``, the symbolic class below;
``isinstance(x, pauliengine.PauliCycle)`` accepts both.

.. autoclass:: pauliengine.pauli_cycle.PauliCycleFactory
   :special-members: __call__

.. autoclass:: pauliengine._core.PauliCycleComplex

.. autoclass:: pauliengine._core.PauliCycleSymbolic

.. autoclass:: pauliengine._core.PauliCycleSumComplex

.. autoclass:: pauliengine._core.PauliCycleSumSymbolic

Pauli orbits
------------

An orbit has no coefficient of its own. ``pauliengine.PauliOrbitSum`` builds a sum
with complex or symbolic coefficients; ``isinstance(x, pauliengine.PauliOrbitSum)``
accepts both.

.. autoclass:: pauliengine.PauliOrbit

.. autoclass:: pauliengine.pauli_orbit.PauliOrbitSumFactory
   :special-members: __call__

.. autoclass:: pauliengine._core.PauliOrbitSumComplex

.. autoclass:: pauliengine._core.PauliOrbitSumSymbolic

Functions
---------

.. autofunction:: pauliengine.inner_product
