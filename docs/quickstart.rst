Quickstart
==========

Pauli strings
-------------

``pe.PauliString`` accepts several input shapes. The coefficient type selects
the backend: a number gives a numeric string, a string gives a symbolic one.

.. code-block:: python

   import pauliengine as pe

   pe.PauliString(1.0, {0: "Z", 1: "X"})         # dict input
   pe.PauliString(1.0, "X0 Y1 Z2")               # space-separated string input
   pe.PauliString(1.0, "XIZ")                    # dense string input
   pe.PauliString((1.0, [("X", 0), ("Y", 2)]))   # OpenFermion-style (coeff, list)
   pe.PauliString("a", {0: "X"})                 # symbolic coefficient

Products carry the correct phases of the Pauli algebra:

.. code-block:: python

   pe.PauliString(1.0, {0: "X"}) * pe.PauliString(1.0, {0: "Y"})   # 1j * Z(0)

Without a coefficient ``pe.PauliString`` gives an unweighted string. It only keeps the
exponent ``k`` of the phase ``i**k`` as an exact integer modulo 4, so its products stay
unweighted.

.. note::

   Unweighted strings and their sets are meant only as basis elements of a Lie
   algebra, for example to compute the dynamical Lie algebra of a set of generators.
   For linear combinations of Pauli strings use weighted strings and Hamiltonians.
   Give the generators as ``i * P`` (``phase=1``) to keep every element of the closure
   anti-Hermitian, and collect a basis in ``pe.QubitHamiltonian`` rather than in a
   Python ``set``: ``hash`` ignores the phase but ``==`` does not, so a ``set`` keeps
   ``i * Z`` and ``-i * Z`` apart.

``commutator`` is the Lie bracket without the factor 2 of ``[P, Q] = 2 * P * Q``:
antisymmetric, with the Jacobi identity and an exact phase. It returns ``None`` for
commuting strings:

.. code-block:: python

   s = pe.PauliString("X") * pe.PauliString("Y")   # +i : Z(0)
   s.phase                                         # 1
   s.phase_factor                                  # 1j
   s.commutator(pe.PauliString("Z"))               # None
   s * 0.5                                         # 0.5j * Z(0)

A dense string such as ``"XIZ"`` is the fastest input for every coefficient type, because
it is read in C++ without building a Python dict. A single string is always read as Pauli
operators, so a symbolic multiple of the identity is ``pe.PauliString("a", {})``.

Next to weighted terms such a string counts with coefficient ``1.0 * i**phase``: in a list
for ``pe.QubitHamiltonian`` and in ``+``, ``-``, ``*``, ``commutator`` and ``inner_product``
with any weighted operator (Pauli strings, Hamiltonians, cycles and orbits). The result then
follows the rules for that weighted operand, for example a ``QubitHamiltonian`` next to a
cycle. A set of unweighted strings (below) is upgraded in the same way, element by element,
and a number turns it into a weighted Hamiltonian.

``pe.QubitHamiltonian`` of such strings is an ordered set and not a linear combination:
every operator occurs once and the first phase seen for it wins, so ``A + A == A``. With a
single unweighted string, ``+`` adds the string to the set and ``*`` and ``commutator`` act
on every element; two unweighted strings have no sum. An element keeps its index, so the
set can serve as a basis. This computes the Lie closure of
a list of generators:

.. code-block:: python

   generators = [pe.PauliString("ZI", phase=1), pe.PauliString("IZ", phase=1), pe.PauliString("XX", phase=1)]
   basis = pe.QubitHamiltonian(generators)
   while True:
       size = len(basis)
       basis += basis.commutator(basis)
       if len(basis) == size:
           break

Hamiltonians
------------

.. code-block:: python

   H = pe.QubitHamiltonian([(1.0, {0: "X"}), ("a", {1: "Z"})])

   H1 = 0.5 * H + (-H) * 2j
   H3 = H ** 3
   C = H.commutator(H1)

   H.qubits()          # sorted list of qubits with non-identity operators
   H.is_hermitian()
   H.dagger()

Every operation that returns a Hamiltonian merges terms with the same operator
and drops zero coefficients, so no explicit ``simplify()`` is needed.

Looking up and hashing terms
----------------------------

``QubitHamiltonian.contains`` returns the coefficient of a Pauli string (``0``
if it does not occur) and accepts the same formats as the constructors:

.. code-block:: python

   H.contains({0: "X"})          # dict
   H.contains([("X", 0)])        # OpenFermion-style list
   H.contains("XI")              # dense string, character i acts on qubit i
   H.contains("X0")              # sparse string

``PauliString.get_hash()`` is an invertible integer key of the operator part
(the coefficient is ignored), and ``PauliString.from_hash`` rebuilds the string:

.. code-block:: python

   key = pe.PauliString(2.0, {0: "X", 2: "Z"}).get_hash()   # 1 + 3 * 4**2 = 49
   pe.PauliString.from_hash(key, 0.5)

Inner products
--------------

``inner_product`` is the Hilbert-Schmidt inner product normalized so that Pauli
strings are orthonormal, :math:`\mathrm{Tr}(A^\dagger B)/2^n`. It is
conjugate-linear in the first argument:

.. code-block:: python

   p = pe.PauliString(2.0, {0: "X"})
   H2 = pe.QubitHamiltonian([(1.0, {0: "X"}), (3.0, {1: "Z"})])
   p.inner_product(H2)             # conj(2.0) * 1.0 = 2.0
   pe.inner_product(H2, H2)        # 1 + 9 = 10

See :doc:`symmetry` for inner products with cycles and orbits.

From a matrix
-------------

``QubitHamiltonian.from_matrix`` decomposes any :math:`2^n \times 2^n` matrix
into Pauli strings in :math:`O(4^n n)` (method of Georges et al., 2025). Qubit 0
is the most significant bit of the row and column index, as in ``to_matrix``:

.. code-block:: python

   M = [[0, 1], [2, 3]]
   pe.QubitHamiltonian.from_matrix(M)       # 1.5 I + 1.5 X - 0.5i Y - 1.5 Z

   H3 = pe.QubitHamiltonian([(0.5, {0: "X", 1: "Z"})])
   pe.QubitHamiltonian.from_matrix(H3.to_matrix(ignore_unused_qubits=False)) == H3

The input may be a numpy array or a list of rows; it is not modified. Terms with
``abs(coefficient) <= threshold`` can be dropped with the ``threshold`` argument.

Symbolic coefficients
---------------------

.. code-block:: python

   H = pe.QubitHamiltonian([("a", {0: "X"}), ("b", {1: "Z"})])
   dH = H.diff("a")              # symbolic derivative
   H2 = H.subs({"a": 2.0})       # substitute and evaluate
