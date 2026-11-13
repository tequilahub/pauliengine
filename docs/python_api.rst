Python API
==========

.. currentmodule:: pauliengine

The factories ``PauliString`` and ``QubitHamiltonian`` pick the numeric or the
symbolic class from the coefficient types; the classes below are what they
return.

Pauli strings
-------------

.. autoclass:: pauliengine.pauli_string.PauliStringFactory
   :members: __call__, from_hash, to_complex

.. autoclass:: pauliengine._core.PauliStringComplex

.. autoclass:: pauliengine._core.PauliStringSymbolic

.. autoclass:: pauliengine._core.unweighted_pauliString

Qubit Hamiltonians
------------------

.. autoclass:: pauliengine.qubit_hamiltonian.QubitHamiltonianFactory
   :members: __call__, from_matrix, from_openfermion, to_openfermion, zero, unit

.. autofunction:: pauliengine.from_openfermion

.. autoclass:: pauliengine._core.QubitHamiltonianComplex

.. autoclass:: pauliengine._core.QubitHamiltonianSymbolic

.. autoclass:: pauliengine._core.unweighted_QubitHamiltonian

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
