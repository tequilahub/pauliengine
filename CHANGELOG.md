# Changelog

All notable changes to PauliEngine are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
