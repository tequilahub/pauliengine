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

Qubit Hamiltonians
------------------

.. autoclass:: pauliengine.qubit_hamiltonian.QubitHamiltonianFactory
   :members: __call__, from_openfermion, to_openfermion, zero, unit

.. autofunction:: pauliengine.from_openfermion

.. autoclass:: pauliengine._core.QubitHamiltonianComplex

.. autoclass:: pauliengine._core.QubitHamiltonianSymbolic

Pauli cycles
------------

``pauliengine.PauliCycle`` and ``pauliengine.PauliCycleSum`` are aliases of the
classes below (complex coefficients).

.. autoclass:: pauliengine._core.PauliCycleComplex

.. autoclass:: pauliengine._core.PauliCycleSumComplex

Pauli orbits
------------

.. autoclass:: pauliengine.PauliOrbit

.. autoclass:: pauliengine.PauliOrbitSum
