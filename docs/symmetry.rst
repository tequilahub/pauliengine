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
drops zeros; ``contains(cycle)`` returns the coefficient of that orbit, and
``unify()`` rewrites every summand to its unified representative. This is
nearly free, because merging already computes the unified forms.

Pauli orbits
------------

A :class:`~pauliengine.PauliOrbit` :math:`B_{p,q,r}` on ``n`` qubits stands for
the normalized, permutation-invariant sum of *all* Pauli strings with ``p``
X-, ``q`` Y- and ``r`` Z-sites,

.. math::

   B_{p,q,r} = \frac{i}{N_\text{terms}} \sum_{P \,:\, (\#X, \#Y, \#Z) = (p, q, r)} P ,
   \qquad N_\text{terms} = \frac{n!}{p!\,q!\,r!\,(n-p-q-r)!} ,

as in the paper. ``n_terms()`` returns :math:`N_\text{terms}`, and
``to_qubit_hamiltonian()`` writes the orbit out (only feasible for small orbits):

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

Arithmetic between representations
----------------------------------

All representations (``PauliString``, ``QubitHamiltonian``, ``PauliCycle``,
``PauliCycleSum``, ``PauliOrbit``, ``PauliOrbitSum``) can be combined with
``+``, ``-``, ``*`` and ``commutator``:

.. code-block:: python

   c = pe.PauliCycle(8, "XZIIIIII")
   o = pe.PauliOrbit(8, 1, 1, 0)
   p = pe.PauliString(0.5, {0: "Y"})

   c + c.commutator(c)        # PauliCycleSum: stays in the cycle family
   o * o                      # PauliOrbitSum: stays in the orbit family
   c * o                      # QubitHamiltonian
   p.commutator(o)            # QubitHamiltonian

Within the cycle family and within the orbit family the result stays
compressed. A product or commutator of two Pauli strings stays a Pauli string,
and two Hamiltonians with the same coefficient type use their own operators.
Every other combination expands both operands and returns a
``QubitHamiltonian``; it is symbolic as soon as one operand is symbolic.
Scalar multiples of cycles and orbits are expressed through their coefficients
(``PauliCycle(n, P, coeff)``, ``PauliOrbitSum.add(orbit, coeff)``); ``*`` with a
plain number is not defined for them.
Expanding an orbit materializes all of its ``n_terms()`` Pauli strings, so this
is only feasible for small orbits (the expansion refuses more than ``2e7``
strings).

Inner products
--------------

``inner_product`` is the Hilbert-Schmidt inner product normalized so that Pauli
strings are orthonormal, :math:`\langle A, B\rangle = \mathrm{Tr}(A^\dagger B)/2^n`
(Corollary 1 and Lemma 1 of arXiv:2604.16701). It is conjugate-linear in the
first argument, works between any two representations, and does not depend on
the register size; multiply by :math:`2^n` for the unnormalized trace.

.. code-block:: python

   o = pe.PauliOrbit(6, 1, 1, 0)
   o.inner_product(o)                    # 1 / o.n_terms()   (Lemma 1)

   c = pe.PauliCycle(4, "ZIZI", 1.5)
   c.period()                            # 2 distinct rotations
   c.inner_product(c)                    # |1.5|^2 * 4^2 / 2

   pe.inner_product(c, pe.PauliString(1.0, {0: "Z", 2: "Z"}))

Orbits are never expanded, so inner products with very large orbits are cheap.
For cycles, :math:`\langle C(P), C(Q)\rangle = \bar c_P c_Q\, n^2/d` if both lie
in the same rotation orbit (period :math:`d`), and 0 otherwise. Note that the
paper's cycle basis carries a factor :math:`i/n`
(:math:`O_P = \tfrac{i}{n}\sum_k R^k P`); with PauliEngine this is
``PauliCycle(n, P, 1j / n)``.
