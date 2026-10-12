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
   pe.PauliString((1.0, [("X", 0), ("Y", 2)]))   # OpenFermion-style (coeff, list)
   pe.PauliString("a", {0: "X"})                 # symbolic coefficient

Products carry the correct phases of the Pauli algebra:

.. code-block:: python

   pe.PauliString(1.0, {0: "X"}) * pe.PauliString(1.0, {0: "Y"})   # 1j * Z(0)

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

Symbolic coefficients
---------------------

.. code-block:: python

   H = pe.QubitHamiltonian([("a", {0: "X"}), ("b", {1: "Z"})])
   dH = H.diff("a")              # symbolic derivative
   H2 = H.subs({"a": 2.0})       # substitute and evaluate
