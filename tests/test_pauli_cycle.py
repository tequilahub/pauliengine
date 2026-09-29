"""Tests for PauliCycle / PauliCycleSum (complex coefficients only).

The central property under test: a PauliCycle represents a base Pauli string
(with a complex coefficient) together with all of its cyclic rotations.
Expanding two cycles into full Hamiltonians and commuting those must give the
same operator as forming the commutator directly on the cycle level
(``PauliCycle.commutator``, Eq. (40) of arXiv:2604.16701).
"""

import pytest

import pauliengine as pe


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
    # NOTE: PauliCycle stores its Base as the canonical representative of the
    # rotation orbit (not the string passed in), so rotate()/base are relative to
    # that representative. These tests therefore compare against ``pc.base``.
    @pytest.mark.parametrize("base", DISTINCT_BASES)
    def test_matches_python_rotation(self, base):
        n = len(base)
        pc = pe.PauliCycle(n, base)
        stored = ps_to_str(pc.base, n)  # canonical representative
        for k in range(-n - 1, 2 * n + 2):
            got = ps_to_str(pc.rotate(k), n)
            expected = py_rotate(stored, k % n)
            assert got == expected, f"base={base} k={k}: {got} != {expected}"

    def test_full_turn_is_identity(self):
        pc = pe.PauliCycle(4, "XYZI")
        n = 4
        assert ps_to_str(pc.rotate(n), n) == ps_to_str(pc.base, n)

    def test_cross_word_boundary(self):
        n = 70
        pc = pe.PauliCycle(n, "X")
        base_pos = next(i for i in range(n) if pc.base.get_pauli_at_index(i) == "X")
        rotated = pc.rotate(65)
        expected_pos = (base_pos + 65) % n
        for i in range(n):
            want = "X" if i == expected_pos else "I"
            assert rotated.get_pauli_at_index(i) == want

    def test_rotations_returns_all(self):
        base = "XYZI"
        pc = pe.PauliCycle(len(base), base)
        got = {ps_to_str(r, len(base)) for r in pc.rotations()}
        assert got == set(cycle_strings(base))


# Coefficients on the base Pauli string


class TestCoefficients:
    def test_rotate_preserves_coefficient(self):
        pc = pe.PauliCycle(4, "XYZI", complex(2.0, -1.0))
        for k in range(4):
            assert pc.rotate(k).coeff == complex(2.0, -1.0)

    def test_coefficient_via_map_constructor(self):
        pc = pe.PauliCycle(3, {0: "X", 2: "Y"}, complex(0.0, 3.0))
        assert pc.base.coeff == complex(0.0, 3.0)

    def test_coefficient_via_pauli_string_constructor(self):
        # The base-PauliString constructor takes the coefficient from the string.
        ps = pe.PauliString(complex(1.5, -0.5), {0: "X", 1: "Z"})
        pc = pe.PauliCycle(4, ps)
        assert pc.base.coeff == complex(1.5, -0.5)

    def test_to_qubit_hamiltonian_scales_with_coefficient(self):
        base = "XYZI"
        n = len(base)
        pc = pe.PauliCycle(n, base, complex(2.0, 0.0))
        assert pc.to_qubit_hamiltonian() == reference_hamiltonian(base, 2.0)


# to_qubit_hamiltonian()


class TestToQubitHamiltonian:
    @pytest.mark.parametrize("base", DISTINCT_BASES)
    def test_equals_manual_sum_of_rotations(self, base):
        pc = pe.PauliCycle(len(base), base)
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

        a = pe.PauliCycle(n, base_a)
        b = pe.PauliCycle(n, base_b)

        reference = reference_hamiltonian(base_a).commutator(reference_hamiltonian(base_b))
        candidate = a.commutator(b).to_qubit_hamiltonian()

        assert candidate == reference
        assert len(reference - candidate) == 0

    def test_matches_with_coefficients(self):
        base_a, base_b = "XYZI", "ZIXI"
        n = len(base_a)
        alpha, beta = complex(2.0, -1.0), complex(0.5, 1.5)

        a = pe.PauliCycle(n, base_a, alpha)
        b = pe.PauliCycle(n, base_b, beta)

        reference = reference_hamiltonian(base_a, alpha).commutator(
            reference_hamiltonian(base_b, beta)
        )
        candidate = a.commutator(b).to_qubit_hamiltonian()
        assert candidate == reference

    def test_self_commutator_is_zero(self):
        base = "XYZI"
        n = len(base)
        pc = pe.PauliCycle(n, base)
        h = reference_hamiltonian(base)
        assert pc.commutator(pc).to_qubit_hamiltonian() == h.commutator(h)

    def test_commuting_cycles_give_zero(self):
        a = pe.PauliCycle(4, "ZIZI")
        b = pe.PauliCycle(4, "ZZII")
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
        a = pe.PauliCycle(n, base_a)
        b = pe.PauliCycle(n, base_b)
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
        a = pe.PauliCycle(4, "ZIZI")
        b = pe.PauliCycle(4, "ZZII")
        assert len(a.commutator(b).data) == 0


# PauliCycleSum structure


class TestPauliCycleSum:
    def test_prefactor(self):
        n = 4
        a = pe.PauliCycle(n, "XYZI")
        b = pe.PauliCycle(n, "ZIXI")
        s = a.commutator(b)
        assert s.coeff == pytest.approx(1.0 / n)
        assert 0 <= len(s.data) <= n

    def test_contains_matches_rotation_orbit(self):
        n = 4
        c = pe.PauliCycle(n, "XYZI")
        c.canonicalize()
        s = pe.PauliCycleSum(1.0, [c])
        # A rotation, once canonicalized, is recognised as the same orbit.
        rot = pe.PauliCycle(n, py_rotate("XYZI", 1))
        rot.canonicalize()
        assert s.contains(rot)
        assert rot in s
        assert not s.contains(pe.PauliCycle(n, "IIII"))


# multiply()  -- product reproduces the full-Hamiltonian product


class TestMultiply:
    @pytest.mark.parametrize(
        "base_a, base_b",
        [("XYZI", "ZIXI"), ("XIY", "ZYX"), ("XXYZ", "ZZIX"), ("XYIZI", "ZIYXI")],
    )
    def test_matches_full_hamiltonian_product(self, base_a, base_b):
        n = len(base_a)
        a = pe.PauliCycle(n, base_a)
        b = pe.PauliCycle(n, base_b)
        reference = reference_hamiltonian(base_a) * reference_hamiltonian(base_b)
        assert a.multiply(b).to_qubit_hamiltonian() == reference
        assert (a * b).to_qubit_hamiltonian() == reference  # __mul__ agrees

    def test_matches_with_coefficients(self):
        base_a, base_b = "XYZI", "ZIXI"
        n = len(base_a)
        alpha, beta = complex(2.0, -1.0), complex(0.5, 1.5)
        a = pe.PauliCycle(n, base_a, alpha)
        b = pe.PauliCycle(n, base_b, beta)
        reference = reference_hamiltonian(base_a, alpha) * reference_hamiltonian(base_b, beta)
        assert a.multiply(b).to_qubit_hamiltonian() == reference

    def test_result_is_deduplicated_with_merged_coeffs(self):
        # "XII" * "XII" on n=3: two of the three products fall into the same
        # rotation orbit and must be merged into one cycle with summed coeff.
        a = pe.PauliCycle(3, "XII")
        s = a.multiply(a)
        # No two summands share a rotation orbit.
        for i in range(len(s.data)):
            for j in range(i + 1, len(s.data)):
                assert s.data[i] != s.data[j]
        assert len(s.data) == 2  # collapsed from 3
        assert sorted(round(c.base.coeff.real, 3) for c in s.data) == [1.0, 2.0]
        # The operator content is unchanged by the merge.
        assert s.to_qubit_hamiltonian() == a.to_qubit_hamiltonian() * a.to_qubit_hamiltonian()


# Sum-level algebra: quadratic combination of every cycle with every cycle


class TestCycleSumAlgebra:
    def _sum(self, n, bases):
        return pe.PauliCycleSum(1.0, [pe.PauliCycle(n, b) for b in bases])

    def test_commutator_matches_expanded(self):
        n = 4
        s1 = self._sum(n, ["XYZI", "ZIXI"])
        s2 = self._sum(n, ["XXYZ", "IZYX"])
        h1, h2 = s1.to_qubit_hamiltonian(), s2.to_qubit_hamiltonian()
        assert s1.commutator(s2).to_qubit_hamiltonian() == h1.commutator(h2)

    def test_multiply_matches_expanded(self):
        n = 4
        s1 = self._sum(n, ["XYZI", "ZIXI"])
        s2 = self._sum(n, ["XXYZ", "IZYX"])
        h1, h2 = s1.to_qubit_hamiltonian(), s2.to_qubit_hamiltonian()
        assert s1.multiply(s2).to_qubit_hamiltonian() == h1 * h2
        assert (s1 * s2).to_qubit_hamiltonian() == h1 * h2  # __mul__ agrees


# Canonical representative: equality / hashing across rotations


class TestEquality:
    @pytest.mark.parametrize("base", DISTINCT_BASES)
    def test_rotations_are_equal_and_hash_equal(self, base):
        # operator== / hash canonicalize on demand, so no explicit canonicalize().
        n = len(base)
        a = pe.PauliCycle(n, base)
        for k in range(n):
            rot = pe.PauliCycle(n, py_rotate(base, k))
            assert a == rot
            assert hash(a) == hash(rot)
        # A set of all rotations collapses to a single element.
        assert len({pe.PauliCycle(n, py_rotate(base, k)) for k in range(n)}) == 1

    def test_equality_is_orbit_invariant_without_canonicalize(self):
        # Two rotations of the same base are equal straight from construction.
        a = pe.PauliCycle(4, "XYZI")
        rot = pe.PauliCycle(4, py_rotate("XYZI", 1))
        assert a == rot
        assert hash(a) == hash(rot)

    def test_distinct_cycles_differ(self):
        # Not rotations of one another, and different qubit counts, must differ.
        assert pe.PauliCycle(4, "XIIZ") != pe.PauliCycle(4, "XIZI")
        assert pe.PauliCycle(3, "XII") != pe.PauliCycle(4, "XIII")


# Targeted commutators: single-coefficient shortcut vs. the full bracket


class TestTargetedCommutator:
    @staticmethod
    def _coeff_of(sum_result, target):
        # Coefficient of `target` in a PauliCycleSum: sum of base coefficients of
        # the cycles whose rotation orbit equals target.
        return sum(complex(c.base.coeff) for c in sum_result.data if c == target)

    @pytest.mark.parametrize(
        "base_a, base_b",
        [("XYZI", "ZIXI"), ("XIY", "ZYX"), ("XXYZ", "ZZIX"), ("XYIZI", "ZIYXI")],
    )
    def test_cycle_targeted_matches_full(self, base_a, base_b):
        n = len(base_a)
        a, b = pe.PauliCycle(n, base_a), pe.PauliCycle(n, base_b)
        full = a.commutator(b)
        for term in full.data:
            assert a.targeted_commutator(b, term) == pytest.approx(self._coeff_of(full, term))

    def test_cycle_sum_targeted_matches_full(self):
        n = 4
        s1 = pe.PauliCycleSum(1.0, [pe.PauliCycle(n, "XYZI"), pe.PauliCycle(n, "ZIXI")])
        s2 = pe.PauliCycleSum(1.0, [pe.PauliCycle(n, "XXYZ"), pe.PauliCycle(n, "IZYX")])
        full = s1.commutator(s2)
        for term in full.data:
            assert s1.targeted_commutator(s2, term) == pytest.approx(self._coeff_of(full, term))
