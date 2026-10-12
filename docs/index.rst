PauliEngine
===========

**Fast arithmetic for quantum operators.** A C++ core with a nanobind Python
frontend for Pauli strings, qubit Hamiltonians and symmetry-adapted operator
bases (Pauli cycles and Pauli orbits). Coefficients can be numeric
(``std::complex<double>``) or fully symbolic (SymEngine).

See `arXiv:2601.02233 <https://arxiv.org/abs/2601.02233>`_ for the algorithmic
background and `arXiv:2604.16701 <https://arxiv.org/abs/2604.16701>`_ for the
cycle and orbit bases.

.. code-block:: python

   import pauliengine as pe

   H = pe.QubitHamiltonian([pe.PauliString(1.0, {0: "Z", 1: "X"}),
                            pe.PauliString("a", {1: "X"})])
   print(H.commutator(H.dagger()))

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   symmetry

.. toctree::
   :maxdepth: 2
   :caption: Reference

   python_api
   cpp_api
