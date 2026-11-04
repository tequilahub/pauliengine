# Changelog

All notable changes to PauliEngine are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Symmetry-adapted operator bases from arXiv:2604.16701 (Pauli cycles and Pauli
orbits), arithmetic and inner products between all representations, and the
Pauli decomposition of matrices.

### Added

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
- Benchmarks under `tools/benchmarks/` for cycles and orbits against the naive
  `QubitHamiltonian` route and against g-sim, with power-law fits, and a smoke
  test comparing nested commutators with a plain-Python computation
  (`tools/commutator_precision_smoke.py`).

### Changed

- Constructing a `PauliCycle` or `PauliOrbit` validates its input: a base acting
  beyond `n_qubits`, `n_qubits < 1` (cycles), an invalid label or
  `n_qubits > 170` (orbits), and operations on mismatched `n_qubits` raise
  `ValueError`. `n_qubits`, `p`, `q` and `r` are read-only in Python.
- Python docstrings of `PauliString` and `QubitHamiltonian` were reviewed and
  corrected (for example `to_string`, `to_dictionary`, `simplify`, `conjugate`,
  and the constructors of the symbolic classes).
- Benchmark result files store every raw timing sample instead of precomputed
  means (`schema_version` 2); the plotting tools derive the statistics and still
  read older files. The benchmark tools are split into `native/`, `gsim/` and
  `utils/`.

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
