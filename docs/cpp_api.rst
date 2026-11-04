C++ API
=======

The library is header-only under ``include/pauliengine/``. ``PauliString`` and
``QubitHamiltonian`` are templates over the coefficient type
(``std::complex<double>`` or ``SymEngine::Expression``); the cycle and orbit
classes use ``std::complex<double>``. Operations exposed in Python are
available in C++ under the same name, with these exceptions:

- ``QubitHamiltonian.from_matrix`` is the free function ``pauli_decomposition``
  (``PauliDecomposition.h``).
- ``inner_product`` and mixed-type ``commutator`` are free functions
  ``inner_product(a, b)`` and ``commutator(a, b)`` (``Arithmetic.h``).
- ``subs`` is ``substitute`` and ``PauliString.get_pauli_at_index`` is
  ``get_pauli_from_index``.

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

Arithmetic between representations
----------------------------------

``include/pauliengine/Arithmetic.h`` overloads ``+``, ``-``, ``*`` and provides
the free functions ``commutator(a, b)`` and ``inner_product(a, b)`` for every pair
of representations.

.. doxygenfile:: Arithmetic.h

Pauli decomposition
-------------------

.. doxygenfile:: PauliDecomposition.h

