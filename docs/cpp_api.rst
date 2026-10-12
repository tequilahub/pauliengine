C++ API
=======

The library is header-only under ``include/pauliengine/``. ``PauliString`` and
``QubitHamiltonian`` are templates over the coefficient type
(``std::complex<double>`` or ``SymEngine::Expression``); the cycle and orbit
classes use ``std::complex<double>``. Every operation exposed in Python is
available in C++ under the same name.

Pauli strings
-------------

.. doxygenclass:: PauliString
   :members:
   :undoc-members:

.. doxygenclass:: WordVec
   :members:
   :undoc-members:

Qubit Hamiltonians
------------------

.. doxygenclass:: QubitHamiltonian
   :members:
   :undoc-members:

Pauli cycles
------------

.. doxygenclass:: PauliCycle
   :members:
   :undoc-members:

.. doxygenclass:: PauliCycleSum
   :members:
   :undoc-members:

Pauli orbits
------------

.. doxygenclass:: PauliOrbit
   :members:
   :undoc-members:

.. doxygenclass:: PauliOrbitSum
   :members:
   :undoc-members:
