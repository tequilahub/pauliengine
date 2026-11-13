# Changelog

All notable changes to PauliEngine are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Symmetry-adapted operator bases from arXiv:2604.16701 (Pauli cycles and Pauli
orbits), arithmetic and inner products between all representations, and the
Pauli decomposition of matrices.

### Added

- **Symbolic coefficients for Pauli cycles and orbit sums**: `PauliCycleSymbolic`,
  `PauliCycleSumSymbolic` and `PauliOrbitSumSymbolic` next to the complex classes
  (`PauliCycleComplex`, `PauliCycleSumComplex`, `PauliOrbitSumComplex`).
  `pe.PauliCycle(4, "XZ", "a")` builds a symbolic cycle, `pe.PauliCycleSum` turns symbolic as soon
  as one cycle is, and `pe.PauliOrbitSum(n, [(orbit, "a"), ...])` or
  `pe.PauliOrbitSum(n, symbolic=True)` gives a symbolic orbit sum; `add` also takes a text
  such as `"a*b"`. `isinstance(x, pe.PauliCycle)` (and likewise for the sums) accepts both
  kinds. Products, commutators, sums, scaling, `targeted_commutator`, `contains` and inner
  products work as for complex coefficients; mixing complex and symbolic operands gives a
  symbolic result, also with Pauli strings and Hamiltonians. A `PauliOrbit` has no
  coefficient and its structure constants stay numeric. Zero checks of merged symbolic
  coefficients expand them first, so `x - x` is empty.
  C++: `PauliCycle<Coeff>`, `PauliCycleSum<Coeff>` and `PauliOrbitSum<Coeff>` are templates
  like `PauliString<Coeff>`; the Python name of the complex orbit sum class in `_core` is
  now `PauliOrbitSumComplex`.
- **Unweighted Pauli strings for Lie algebras**: `pe.PauliString` without a coefficient
  (`pe.PauliString("XIZ")`, `phase=1` for `i * P`) stores only the exponent `k` of the
  phase `i**k` as an exact integer modulo 4. `pe.QubitHamiltonian` of such strings is an
  ordered set with one phase per element, stable indices and a hash index; it is not a
  linear combination, so `A + A == A`. `commutator` omits the factor 2 for both. A scalar
  or a weighted factor gives a weighted result. Next to any weighted operator (Pauli
  strings, Hamiltonians, cycles and orbits) an unweighted string or set counts with
  coefficient `1.0 * i**phase` per element in `+`, `-`, `*`, `commutator` and
  `inner_product`, and the result follows the rules for the weighted operand. A set with
  a single unweighted string stays unweighted (`+` adds it, `*` and `commutator` act on
  every element); two unweighted strings have no sum.
- `pe.PauliString(coeff, "XIZ")` takes a dense string for complex and symbolic
  coefficients, which is much faster than a dict. Invalid operators in a text and a
  qubit given twice now raise `ValueError`.
- **Pauli cycles** (translation-invariant operators): `PauliCycle(n, base, coeff)`
  stands for `coeff * sum_k R^k P` over all `n` cyclic rotations of `P`; the base
  can be a dict, a dense string such as `"XIY"` or a `PauliString`.
  - `commutator` (Eq. 40) and the product `*` return a `PauliCycleSum`; only shifts at
    which the supports overlap are evaluated, so the commutator of bounded-weight
    cycles costs O(n / 64) per shift.
  - `targeted_commutator(other, target)` returns a single coefficient without
    building the full commutator.
  - `==` and `hash` are rotation-invariant and ignore the coefficient.
    `unify()` rewrites the base to a unique representative of its orbit, computed
    in O(n / 64 + w) for weight `w`.
  - `period()` returns the number of distinct rotations; `rotate(k)`,
    `rotations()` and `to_qubit_hamiltonian()` expand the cycle.
  - `get_hash()` / `from_hash(hash, coeff=1)`: invertible integer key, identical
    for all rotations.
- `PauliCycleSum`: a sum of cycles that always merges summands of the same orbit
  and drops zeros; `commutator`, `*`, `targeted_commutator`, `unify()`,
  `to_qubit_hamiltonian()`, `contains(cycle)` (returns the coefficient), `in` and
  `len`.
- **Pauli orbits** (permutation-invariant operators): `PauliOrbit(n, p, q, r)` is
  `B_{p,q,r} = (i / N_terms) * sum` of all Pauli strings with `p` X-, `q` Y- and
  `r` Z-sites.
  - `commutator` (Theorem 1) and the product `*` work on the contingency tables of the
    labels and never expand the strings; `targeted_commutator` (Corollary 5)
    returns a single structure constant and picks the cheaper of two exact
    routes.
  - `n_terms()`, `to_qubit_hamiltonian()` (refuses more than 2e7 strings),
    `get_hash()` / `from_hash(hash)`, `PauliOrbit.max_qubits` (170).
- `PauliOrbitSum(n)`: weighted sum of orbits with `add`, `contains` (returns the
  coefficient), `prune_zeros`, `commutator`, `*`, `targeted_commutator`,
  `to_qubit_hamiltonian()` and `len`.
- **Arithmetic between all representations**: `+`, `-`, `*` and `commutator`
  work for any pair of `PauliString`, `QubitHamiltonian`, `PauliCycle`,
  `PauliCycleSum`, `PauliOrbit` and `PauliOrbitSum`. Within the cycle and orbit
  families the result stays a `PauliCycleSum` / `PauliOrbitSum`; other
  combinations return a `QubitHamiltonian`, symbolic as soon as one operand is
  symbolic. In C++ these are overloaded operators and the free function
  `commutator(a, b)` (`Arithmetic.h`).
- Cycles and orbits can be scaled by an `int` or `float` on either side. A
  `PauliCycle`, `PauliCycleSum` or `PauliOrbitSum` keeps its type; a `PauliOrbit`,
  which has no coefficient of its own, becomes a `PauliOrbitSum`. Scaling a sum by
  0 gives the empty sum.
- **Inner products**: `a.inner_product(b)` and `pe.inner_product(a, b)` give the
  Hilbert-Schmidt inner product normalized so that Pauli strings are
  orthonormal, `Tr(a^dagger b) / 2^n` (Corollary 1, Lemma 1), for every pair of
  representations. Orbits are never expanded, so inner products with very large
  orbits are cheap.
- **Pauli decomposition**: `QubitHamiltonian.from_matrix(matrix, threshold=0.0)`
  turns any `2^n x 2^n` matrix (numpy array or list of rows) into a
  `QubitHamiltonianComplex` in O(4^n n), following Georges et al., New J. Phys.
  (2025). It is the inverse of `to_matrix(ignore_unused_qubits=False)`; in C++
  `pauli_decomposition` (`PauliDecomposition.h`).
- `targeted_commutator(other, target)` for `PauliString` and `QubitHamiltonian`:
  the coefficient of one target Pauli string in the commutator; for Hamiltonians
  as a hash join instead of the full O(|A| |B|) pair scan.
- Documentation: Sphinx site under `docs/` (user guide, Python and C++
  reference via Breathe/Doxygen); install the tools with the `docs` dependency
  group.
- A smoke test comparing nested commutators with a plain-Python computation
  (`tools/commutator_precision_smoke.py`).

### Changed

- `pe.PauliString()`, `pe.PauliString("a")` and other calls without a coefficient no longer
  give a zero string or a symbolic multiple of the identity, see the Pauli strings above.
  Use `pe.PauliString("a", {})` for the latter.
- `PauliString` takes 80 instead of 104 bytes: `WordVec` stores its size and capacity
  as 32-bit integers, and the unused C++ member `is_zero` is gone (the Python property
  `is_zero` is unchanged). Writing the result of large products gets faster.
- Faster `QubitHamiltonian` product, commutator and compaction for complex
  coefficients on up to 512 qubits: pairs are computed without temporary Pauli
  strings, terms are sorted with a radix sort, and large inputs run through a parallel
  sample sort without a serial merge step. Coefficients of equal operators are added
  in input order, so results no longer depend on the thread count.
- `QubitHamiltonian` with complex coefficients stores its terms flat: one array of
  operator words (every term padded to the same width) and one of coefficients,
  instead of one `PauliString` per term. Products, commutators, sums, differences,
  scaling and compaction work on these arrays directly; beyond 512 qubits they sort
  an index instead of whole records. Commutators skip commuting pairs before forming
  the product. Symbolic Hamiltonians are stored as before.
  C++ only: the public member `data` is gone. Use `terms()` (all terms as Pauli
  strings; built on request for complex coefficients), `term(i)` and `size()`.
  The Python API is unchanged.
- No OpenMP region is opened for work below the threshold; even an inactive one
  costs microseconds with libgomp. `PAULIENGINE_OMP_THRESHOLD` now defaults to 4096
  with MSVC and 32768 with other compilers, counted in records for sorting and in
  pairs times words for forming products.
- Constructors no longer copy all terms through the return value of `compact()`.
- Constructing a `PauliCycle` or `PauliOrbit` validates its input: a base acting
  beyond `n_qubits`, `n_qubits < 1` (cycles), an invalid label or
  `n_qubits > 170` (orbits), and operations on mismatched `n_qubits` raise
  `ValueError`. `n_qubits`, `p`, `q` and `r` are read-only in Python.
- Python docstrings of `PauliString` and `QubitHamiltonian` were reviewed and
  corrected (for example `to_string`, `to_dictionary`, `simplify`, `conjugate`,
  and the constructors of the symbolic classes).

### Fixed

- `PauliString.qubits()` returned 64 zeros (plus all indices 0..63 if qubit 0
  was occupied) instead of the qubits with a non-identity operator.
- `PauliString.key_openfermion()` prepended 64 invalid `('\x00', 0)` entries.
- `PauliString.is_zero` was always true; it now reports whether the coefficient
  is zero.
- The `__rmul__` overloads between complex and symbolic Pauli strings could never
  be called.

## [0.1.4] - 2026-09-30

### Added

- `QubitHamiltonian.contains(pauli_string)` returns the coefficient of a Pauli
  string in the Hamiltonian, or `0` if it does not occur. The coefficient of the
  query itself is ignored. Accepted formats, for both `QubitHamiltonianComplex`
  and `QubitHamiltonianSymbolic`:
  - a `PauliStringComplex` or `PauliStringSymbolic`
  - a dict `{qubit: "X"}`
  - an OpenFermion-style list `[("X", 0), ("Z", 2)]`
  - a dense string `"XIZ"` (character `i` acts on qubit `i`)
  - a sparse string `"X0 Z2"`

  Operator labels are case-insensitive; invalid labels raise `ValueError`.
- `PauliString.get_hash()` returns an invertible, coefficient-independent key of
  the operator part as a Python `int`: `sum_q p_q * 4**q` with `I=0, X=1, Y=2, Z=3`.
  Equal operators always share a key and distinct operators never collide, for
  any number of qubits.
- `PauliString.from_hash(hash, coeff=1)` rebuilds the Pauli string from such a
  key. Through the `pe.PauliString` factory the coefficient type selects the
  class (number → `PauliStringComplex`, string → `PauliStringSymbolic`), as in
  the constructor.
- `PauliStringComplex` and `PauliStringSymbolic` are now hashable (`__hash__`),
  so they can be used in sets and as dict keys. The hash is a fixed-size hash of
  the operator part; the coefficient is ignored. Note that `__eq__` still
  compares coefficients, so a dict lookup only matches a Pauli string with the
  same coefficient. Use `get_hash()` as key to look up by operator alone.

### Fixed

- `pip install .` now works with the Visual Studio (multi-config) CMake
  generator on Windows. Previously CMake also configured a Debug configuration,
  for which Conan provides no GMP/SymEngine libraries, and the build failed with
  `IMPORTED_LOCATION not set for imported target "CONAN_LIB::gmp_gmp_RELEASE"`.

## [0.1.3] - 2026-07-14

Starting point of this changelog; the release on PyPI at the time it was created.

### Fixed

- Deployment workflow for publishing to PyPI.

[Unreleased]: https://github.com/tequilahub/pauliengine/compare/v0.1.4...HEAD
[0.1.4]: https://github.com/tequilahub/pauliengine/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/tequilahub/pauliengine/compare/v0.1.2...v0.1.3
