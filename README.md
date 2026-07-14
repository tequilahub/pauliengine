# PauliEngine — Fast Arithmetic for Quantum Operators

A C++ core with a nanobind Python frontend for working with Pauli strings and
qubit Hamiltonians. Coefficients can be numeric (`std::complex<double>`) or
fully symbolic (SymEngine).

See [arXiv:2601.02233](https://arxiv.org/abs/2601.02233) for the algorithmic
background.

> **Status.** Functional core with a large test suite. The Python
> API surface is stable. Sphinx docs are not published yet.

---

## Highlights

- **Binary symplectic representation** — multiplication, commutators, and
  hashing are O(n / 64) over the qubit count.
- **Two coefficient backends** — numeric (`complex<double>`) and symbolic
  (`SymEngine::Expression`); cross-type arithmetic is handled automatically and
  type-promotes the way you expect (any symbolic term → the whole Hamiltonian
  becomes symbolic).
- **Compaction invariant** — every operation that returns a `QubitHamiltonian`
  merges duplicate-operator terms and drops zero-coefficient terms. You never
  need to call `simplify()` for correctness.
- **Tequila-compatible API** — `qubits`, `n_qubits`, `is_hermitian`,
  `is_antihermitian`, `dagger`, `conjugate`, `transpose`, `simplify(threshold)`,
  `split`, `map_qubits`, `power` / `__pow__`, `to_matrix`, `paulistrings`,
  `count_measurements`, `is_all_z`.
- **OpenFermion bridge** — the `QubitHamiltonian` factory accepts an
  `openfermion.QubitOperator` directly; `to_openfermion` and `from_openfermion`
  helpers are available.

---

## Quickstart

```python
import pauliengine as pe

# A single Pauli string: dict-of-operators style.
p1 = pe.PauliString(1.0, {0: "Z", 1: "X"})

# Symbolic coefficient — anything coercible to a SymEngine Expression.
p2 = pe.PauliString("a", {1: "X"})

# OpenFermion-style: PauliString((coeff, [(Pauli, qubit), ...]))
p3 = pe.PauliString((1.0, [("X", 0), ("Y", 2)]))

# Build a Hamiltonian from a list of PauliStrings (or (coeff, dict) tuples).
H = pe.QubitHamiltonian([p1, p2, p3])

print(H.qubits())        # [0, 1, 2]
print(H.is_hermitian())  # True
print(H.dagger())        # for Hermitian H this is just H
```

Tests are in `tests/`:

```bash
pytest tests/
```

---

## Construction

`pe.PauliString` accepts several input shapes:

```python
pe.PauliString(1.0, {0: "Z", 1: "X"})         # dict input
pe.PauliString(1.0, "X0 Y1 Z2")               # space-separated string input
pe.PauliString((1.0, [("X", 0), ("Y", 2)]))   # OpenFermion-style (coeff, list)
pe.PauliString("a", {0: "X"})                 # symbolic coefficient
```

`pe.QubitHamiltonian` accepts:

```python
pe.QubitHamiltonian([ps1, ps2, ...])                       # list of PauliStrings
pe.QubitHamiltonian([(1.0, {0: "X"}), ("a", {1: "Z"})])    # list of tuples
pe.QubitHamiltonian(openfermion_qubit_operator)            # see OpenFermion bridge
pe.QubitHamiltonian.zero()                                 # empty Hamiltonian
pe.QubitHamiltonian.unit()                                 # identity (single term, coeff 1, no ops)
```

If any term in the list is symbolic, the resulting Hamiltonian is symbolic.

---

## Arithmetic

```python
# PauliString * PauliString, with the right factors of i from Pauli algebra.
p4 = pe.PauliString(1.0, {0: "X"}) * pe.PauliString(1.0, {0: "Y"})  # -> 1j * Z(0)

# Scalar multiplication on both sides; +, -, unary -, and addition between
# PauliStrings (returns a QubitHamiltonian).
H1 = p1 + p2 - p3
H2 = 0.5 * H1 + (-H1) * 2j
H3 = H1 ** 3                  # integer powers
c  = H1.commutator(H2)        # commutator (also available on PauliString)
```

Cross-type multiplication (numeric × symbolic) is supported and promotes the
result to symbolic.

---

## Inspection and properties

```python
H.size()                # number of Pauli-string terms (also len(H))
H.qubits()              # sorted list of qubits with non-identity operators
H.n_qubits()            # len(H.qubits())
H.is_all_z()
H.is_hermitian()        # True iff every coefficient is real
H.is_antihermitian()    # True iff every coefficient is purely imaginary
H.count_measurements()  # 1 if all-Z, else len(H)
```

For a single `PauliString`:

```python
ps.size()               # number of non-identity Pauli ops (also len(ps))
ps.count_y()            # number of Y operators (used by conjugate/transpose)
ps.naked()              # same operator with coefficient 1
ps.key_openfermion()    # OpenFermion-style key
ps.get_pauli_at_index(q)  # "I" / "X" / "Y" / "Z"
```

---

## Transformations

```python
H.dagger()           # complex-conjugate each coefficient
H.conjugate()        # complex conjugation (flips a sign per Y operator)
H.transpose()        # transpose (flips a sign per Y operator, no conjugation)
H.simplify(1e-10)    # drop terms with |coefficient| <= threshold
H.split()            # -> (hermitian, anti_hermitian) pair (numeric coeffs only)
H.map_qubits({0: 5, 1: 2})
H.power(3)           # also via H ** 3
```

> `split()` and `to_matrix()` require coefficients that evaluate to a complex
> number — call `H.subs({...})` first on symbolic Hamiltonians.

### Dense matrix

```python
import numpy as np

M = np.array(H.to_matrix())                          # 2**n x 2**n, ignores unused qubits
M_full = np.array(H.to_matrix(ignore_unused_qubits=False))  # absolute qubit indices
```

---


## Symbolic coefficients

Any string coefficient (or `SymEngine::Expression` from C++) makes the term
symbolic. Symbolic PauliStrings and Hamiltonians support every arithmetic
operation plus:

```python
H = pe.QubitHamiltonian([("a", {0: "X"}), ("b", {1: "Z"})])

dH = H.diff("a")            # symbolic derivative
H2 = H.subs({"a": 2.0})     # substitute and evaluate
```

`diff` is also available on `PauliString`. Mixing symbolic and numeric inputs
is fine: the factory scans every element and uses the symbolic builder if
needed.

---

## OpenFermion bridge

```python
from openfermion import QubitOperator
qop = 1.5 * QubitOperator("X0 Y1") + 0.5j * QubitOperator("Z2")

# Factory accepts QubitOperator directly:
H = pe.QubitHamiltonian(qop)

# Or use the explicit helpers:
H = pe.from_openfermion(qop)
qop_back = pe.QubitHamiltonian.to_openfermion(H)
assert qop == qop_back
```

`openfermion` is an optional dependency — `from_openfermion` /
`to_openfermion` import it lazily and raise `ImportError` with a helpful
message if it is missing.

---

## C++ usage

The library is a header-only template under `include/pauliengine/`. Both
`PauliString<Coeff>` and `QubitHamiltonian<Coeff>` work for
`Coeff = std::complex<double>` and `Coeff = SymEngine::Expression`. Every
operation exposed in Python is available in C++ with the same name.

---

## Installation

### From PyPI

```bash
pip install pauliengine
```

Prebuilt wheels are available for Linux (x86_64/aarch64), macOS (arm64) and
Windows (x86_64/arm64) on Python 3.10–3.13. Installing from the source
distribution requires the build dependencies below.

### Build dependencies

- A C++20 compiler (MSVC 19.3+, GCC 11+, or Clang 14+)
- CMake 3.20+
- Python 3.10–3.13
- [Conan 2](https://conan.io) (to pull in SymEngine)
- [nanobind](https://github.com/wjakob/nanobind) (build-time)

### Install from source

```bash
pip install conan
conan profile detect
conan install . --output-folder=build --build=missing
pip install .
```

The CMake build picks up the Conan toolchain from `build/conan_toolchain.cmake`.
For an editable / development install use `pip install -e .` instead of
`pip install .`.

### Windows notes

If you are on Windows and have not built SymEngine before, the Conan step will
build it from source on first install — that takes a few minutes. Subsequent
builds use the cached artifact.

> **Performance note.** PauliEngine can be built without SymEngine, but
> symbolic coefficients are unavailable in that mode and the runtime cost of
> certain numeric paths increases.

### macOS prerequisite

Conan does not currently ship a prebuilt SymEngine binary for macOS. Build it
from source once before the main install step:

```bash
conan install --build=symengine/0.14.0
```

Subsequent installs pick up the cached artifact, so this only needs to be done
the first time.

---



## Testing

```bash
pip install pytest
pytest tests/
```


---

## Citation

If you use PauliEngine in academic work, please cite
[arXiv:2601.02233](https://arxiv.org/abs/2601.02233).
