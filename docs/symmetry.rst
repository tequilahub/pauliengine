Pauli cycles and Pauli orbits
=============================

For translation- or permutation-invariant operators PauliEngine works directly
in a symmetry-adapted basis (`arXiv:2604.16701
<https://arxiv.org/abs/2604.16701>`_) instead of expanding into exponentially or
polynomially many Pauli strings.

Pauli cycles
------------

A :class:`~pauliengine.PauliCycle` on ``n`` qubits stands for a base Pauli
string together with all of its cyclic rotations,

.. math::

   C(P) = c \sum_{k=0}^{n-1} R^k P ,

where :math:`R` shifts every qubit by one position. The sum always runs over
all ``n`` rotations, so a base with a shorter period counts each distinct
rotation several times (the cycle of ``ZZZ`` on three qubits is ``3 * ZZZ``).

.. code-block:: python

   import pauliengine as pe

   a = pe.PauliCycle(8, "XIZIIIII")
   b = pe.PauliCycle(8, "YYIIIIII", 0.5)

   s = a.commutator(b)            # PauliCycleSum, O(n / 64) for bounded weight
   p = a * b                      # product, also a PauliCycleSum
   s.to_qubit_hamiltonian()       # expand into an explicit QubitHamiltonian

Two cycles are equal if they describe the same rotation orbit; the coefficient
is ignored. ``unify()`` rewrites the stored base to the unified representative
of its orbit (the operator on qubit 0, the smallest gaps first); it is optional,
since ``==`` and ``hash`` are rotation-invariant on their own.

A :class:`~pauliengine.PauliCycleSum` merges summands from the same orbit and
drops zeros; ``contains(cycle)`` returns the coefficient of that orbit.

Pauli orbits
------------

A :class:`~pauliengine.PauliOrbit` :math:`B_{p,q,r}` on ``n`` qubits stands for
the normalized, permutation-invariant sum of *all* Pauli strings with ``p``
X-, ``q`` Y- and ``r`` Z-sites:

.. code-block:: python

   a = pe.PauliOrbit(20, 1, 1, 0)       # all strings with one X and one Y
   b = pe.PauliOrbit(20, 0, 0, 2)

   s = a.commutator(b)                  # PauliOrbitSum, independent of the 20**2 strings
   a.targeted_commutator(b, s.data[0])  # a single structure constant

Orbits are limited to ``n <= PauliOrbit.max_qubits`` (170): the structure
constants use ``n!`` in double precision.

Hash keys
---------

Like Pauli strings, cycles and orbits have invertible integer keys that are
consistent with ``==``:

.. code-block:: python

   key = pe.PauliCycle(8, "IIXIZIII").get_hash()   # identical for every rotation
   pe.PauliCycle.from_hash(key, 2.0)               # base = unified representative

   pe.PauliOrbit.from_hash(pe.PauliOrbit(20, 1, 1, 0).get_hash())
