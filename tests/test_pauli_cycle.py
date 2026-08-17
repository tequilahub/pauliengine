"""Tests for PauliCycle / PauliCycleSum (complex coefficients only).

The central property under test: a PauliCycle represents a base Pauli string
(with a complex coefficient) together with all of its cyclic rotations.
Expanding two cycles into full Hamiltonians and commuting those must give the
same operator as forming the commutator directly on the cycle level
(``PauliCycle.commutator``, Eq. (40) of arXiv:2604.16701).
"""

import pytest

import pauliengine as pe
from pauliengine._core import PauliCycleComplex, PauliStringComplex


# Helpers


def py_rotate(s: str, k: int) -> str:
    """Pure-Python reference rotation: new[(i + k) mod n] = old[i]."""
    n = len(s)
    out = [""] * n
    for i, ch in enumerate(s):
        out[(i + k) % n] = ch
    return "".join(out)


def ps_to_str(ps, n: int) -> str:
    """Render a PauliString as an n-character string of I/X/Y/Z."""
    return "".join(ps.get_pauli_at_index(i) for i in range(n))


def cycle_strings(base: str) -> list[str]:
    """All n cyclic rotations of ``base`` as strings."""
    return [py_rotate(base, k) for k in range(len(base))]


def anticommute(s1: str, s2: str) -> bool:
    """Do two Pauli strings anticommute? (odd number of anticommuting sites)"""
    count = sum(
        1 for a, b in zip(s1, s2) if a != "I" and b != "I" and a != b
    )
    return count % 2 == 1


def reference_hamiltonian(base: str, coeff: complex = 1.0):
    """Build the full QubitHamiltonian: coeff * sum of all rotations of base."""
    terms = []
    for s in cycle_strings(base):
        ops = {i: ch for i, ch in enumerate(s) if ch != "I"}
        terms.append(pe.PauliString(complex(coeff), ops))
    return pe.QubitHamiltonian(terms)


# Bases chosen so that all n rotations are distinct (no cyclic sub-period).
DISTINCT_BASES = ["XIY", "XYZI", "ZIIX", "XXYZ", "XYIZI"]


# rotate()


class TestRotate:
    @pytest.mark.parametrize("base", DISTINCT_BASES)
    def test_matches_python_rotation(self, base):
        n = len(base)
        pc = PauliCycleComplex(n, base)
        for k in range(-n - 1, 2 * n + 2):
            got = ps_to_str(pc.rotate(k), n)
            expected = py_rotate(base, k % n)
            assert got == expected, f"base={base} k={k}: {got} != {expected}"

    def test_full_turn_is_identity(self):
        base = "XYZI"
        pc = PauliCycleComplex(len(base), base)
        assert ps_to_str(pc.rotate(len(base)), len(base)) == base

    def test_cross_word_boundary(self):
        pc = PauliCycleComplex(70, "X")
        rotated = pc.rotate(65)
        assert rotated.get_pauli_at_index(0) == "I"
        assert rotated.get_pauli_at_index(65) == "X"

    def test_rotations_returns_all(self):
        base = "XYZI"
        pc = PauliCycleComplex(len(base), base)
        got = {ps_to_str(r, len(base)) for r in pc.rotations()}
        assert got == set(cycle_strings(base))


# Coefficients on the base Pauli string


class TestCoefficients:
    def test_rotate_preserves_coefficient(self):
        pc = PauliCycleComplex(4, "XYZI", complex(2.0, -1.0))
        for k in range(4):
            assert pc.rotate(k).coeff == complex(2.0, -1.0)

    def test_coefficient_via_map_constructor(self):
        pc = PauliCycleComplex(3, {0: "X", 2: "Y"}, complex(0.0, 3.0))
        assert pc.base.coeff == complex(0.0, 3.0)

    def test_coefficient_via_pauli_string_constructor(self):
        # The base-PauliString constructor takes the coefficient from the string.
        ps = PauliStringComplex(complex(1.5, -0.5), {0: "X", 1: "Z"})
        pc = PauliCycleComplex(4, ps)
        assert pc.base.coeff == complex(1.5, -0.5)

    def test_to_qubit_hamiltonian_scales_with_coefficient(self):
        base = "XYZI"
        n = len(base)
        pc = PauliCycleComplex(n, base, complex(2.0, 0.0))
        assert pc.to_qubit_hamiltonian() == reference_hamiltonian(base, 2.0)


# to_qubit_hamiltonian()


class TestToQubitHamiltonian:
    @pytest.mark.parametrize("base", DISTINCT_BASES)
    def test_equals_manual_sum_of_rotations(self, base):
        pc = PauliCycleComplex(len(base), base)
        assert pc.to_qubit_hamiltonian() == reference_hamiltonian(base)


# commutator()  -- the main property


class TestCommutator:
    @pytest.mark.parametrize(
        "base_a, base_b",
        [
            ("XYZI", "ZIXI"),
            ("XIY", "ZYX"),
            ("XXYZ", "ZZIX"),
            ("XYIZI", "ZIYXI"),
        ],
    )
    def test_matches_full_hamiltonian_commutator(self, base_a, base_b):
        n = len(base_a)
        assert len(base_b) == n

        a = PauliCycleComplex(n, base_a)
        b = PauliCycleComplex(n, base_b)

        reference = reference_hamiltonian(base_a).commutator(reference_hamiltonian(base_b))
        candidate = a.commutator(b).to_qubit_hamiltonian()

        assert candidate == reference
        assert len(reference - candidate) == 0

    def test_matches_with_coefficients(self):
        base_a, base_b = "XYZI", "ZIXI"
        n = len(base_a)
        alpha, beta = complex(2.0, -1.0), complex(0.5, 1.5)

        a = PauliCycleComplex(n, base_a, alpha)
        b = PauliCycleComplex(n, base_b, beta)

        reference = reference_hamiltonian(base_a, alpha).commutator(
            reference_hamiltonian(base_b, beta)
        )
        candidate = a.commutator(b).to_qubit_hamiltonian()
        assert candidate == reference

    def test_self_commutator_is_zero(self):
        base = "XYZI"
        n = len(base)
        pc = PauliCycleComplex(n, base)
        h = reference_hamiltonian(base)
        assert pc.commutator(pc).to_qubit_hamiltonian() == h.commutator(h)

    def test_commuting_cycles_give_zero(self):
        a = PauliCycleComplex(4, "ZIZI")
        b = PauliCycleComplex(4, "ZZII")
        result = a.commutator(b).to_qubit_hamiltonian()
        assert len(result) == 0


# Corollary 3: only anticommuting shifts appear as cycle terms


class TestCorollary3:
    @pytest.mark.parametrize(
        "base_a, base_b",
        [
            ("XYZI", "ZIXI"),
            ("XIY", "ZYX"),
            ("XXYZ", "ZZIX"),
            ("XYIZI", "ZIYXI"),
        ],
    )
    def test_only_nonzero_shifts_are_kept(self, base_a, base_b):
        n = len(base_a)
        a = PauliCycleComplex(n, base_a)
        b = PauliCycleComplex(n, base_b)
        s = a.commutator(b)

        # The cycle sum should contain exactly one term per shift k for which
        # [P, rotate(P', k)] != 0, i.e. the two operators anticommute.
        expected_shifts = [
            k for k in range(n) if anticommute(base_a, py_rotate(base_b, k))
        ]
        assert len(s.data) == len(expected_shifts)
        # Every retained cycle carries a non-zero (physical factor 2) coefficient.
        for cycle in s.data:
            assert cycle.base.coeff != 0

    def test_commuting_pair_produces_no_terms(self):
        a = PauliCycleComplex(4, "ZIZI")
        b = PauliCycleComplex(4, "ZZII")
        assert len(a.commutator(b).data) == 0


# PauliCycleSum structure


class TestPauliCycleSum:
    def test_prefactor(self):
        n = 4
        a = PauliCycleComplex(n, "XYZI")
        b = PauliCycleComplex(n, "ZIXI")
        s = a.commutator(b)
        assert s.coeff == pytest.approx(1.0 / n)
        assert 0 <= len(s.data) <= n
