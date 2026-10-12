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
    # NOTE: PauliCycle stores its Base exactly as passed in (construction does not
    # unify), and rotate() is relative to that stored base. The tests compare
    # against ``pc.base`` so they also hold after an explicit unify().
    @pytest.mark.parametrize("base", DISTINCT_BASES)
    def test_matches_python_rotation(self, base):
        n = len(base)
        pc = pe.PauliCycle(n, base)
        stored = ps_to_str(pc.base, n)  # the stored base
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
        # Different shifts can land in the same rotation orbit; the sum merges
        # those, so count distinct orbits of the anticommuting products.
        expected_shifts = [
            k for k in range(n) if anticommute(base_a, py_rotate(base_b, k))
        ]
        products = [a.base.commutator(b.rotate(k)) for k in expected_shifts]
        assert len(s.data) == len({pe.PauliCycle(n, ps) for ps in products})
        # Every retained cycle carries a non-zero (physical factor 2) coefficient.
        for cycle in s.data:
            assert cycle.base.coeff != 0

    def test_commuting_pair_produces_no_terms(self):
        a = pe.PauliCycle(4, "ZIZI")
        b = pe.PauliCycle(4, "ZZII")
        assert len(a.commutator(b).data) == 0


# PauliCycleSum structure


class TestPauliCycleSum:
    def test_no_global_prefactor(self):
        # Each cycle carries its own coefficient; there is no separate prefactor.
        s = pe.PauliCycleSum([pe.PauliCycle(3, "XII", 2.0)])
        assert not hasattr(s, "coeff")
        assert s.to_qubit_hamiltonian() == reference_hamiltonian("XII", 2.0)

    def test_contains_matches_rotation_orbit(self):
        n = 4
        c = pe.PauliCycle(n, "XYZI")
        c.unify()
        s = pe.PauliCycleSum([c])
        # A rotation, once unified, is recognised as the same orbit.
        rot = pe.PauliCycle(n, py_rotate("XYZI", 1))
        rot.unify()
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
        return pe.PauliCycleSum([pe.PauliCycle(n, b) for b in bases])

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


# Unified representative: equality / hashing across rotations


class TestEquality:
    @pytest.mark.parametrize("base", DISTINCT_BASES)
    def test_rotations_are_equal_and_hash_equal(self, base):
        # operator== / hash unify on demand, so no explicit unify().
        n = len(base)
        a = pe.PauliCycle(n, base)
        for k in range(n):
            rot = pe.PauliCycle(n, py_rotate(base, k))
            assert a == rot
            assert hash(a) == hash(rot)
        # A set of all rotations collapses to a single element.
        assert len({pe.PauliCycle(n, py_rotate(base, k)) for k in range(n)}) == 1

    def test_equality_is_orbit_invariant_without_unify(self):
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
        s1 = pe.PauliCycleSum([pe.PauliCycle(n, "XYZI"), pe.PauliCycle(n, "ZIXI")])
        s2 = pe.PauliCycleSum([pe.PauliCycle(n, "XXYZ"), pe.PauliCycle(n, "IZYX")])
        full = s1.commutator(s2)
        for term in full.data:
            assert s1.targeted_commutator(s2, term) == pytest.approx(self._coeff_of(full, term))


# Input validation and invariants


class TestValidation:
    def test_base_beyond_n_qubits_raises(self):
        with pytest.raises(ValueError):
            pe.PauliCycle(3, {5: "X"})
        with pytest.raises(ValueError):
            pe.PauliCycle(3, "XIIZ")

    def test_non_positive_n_qubits_raises(self):
        with pytest.raises(ValueError):
            pe.PauliCycle(0, "")

    def test_identity_padding_beyond_n_is_fine(self):
        # Trailing identities carry no support, so they are not out of range.
        assert pe.PauliCycle(3, {0: "X", 7: "I"}) == pe.PauliCycle(3, "XII")

    @pytest.mark.parametrize("op", ["commutator", "multiply"])
    def test_n_qubits_mismatch_raises(self, op):
        a, b = pe.PauliCycle(3, "XYZ"), pe.PauliCycle(4, "ZIXI")
        with pytest.raises(ValueError):
            getattr(a, op)(b)

    def test_targeted_n_qubits_mismatch_raises(self):
        a, b = pe.PauliCycle(4, "XYZI"), pe.PauliCycle(4, "ZIXI")
        with pytest.raises(ValueError):
            a.targeted_commutator(b, pe.PauliCycle(3, "XII"))

    def test_sum_with_mixed_n_qubits_raises(self):
        with pytest.raises(ValueError):
            pe.PauliCycleSum([pe.PauliCycle(3, "XII"), pe.PauliCycle(4, "XIII")])

    def test_n_qubits_is_read_only(self):
        # Mutating it would silently invalidate the cached unified form.
        c = pe.PauliCycle(4, "XIII")
        with pytest.raises(AttributeError):
            c.n_qubits = 3


class TestSumMerging:
    def test_commutator_merges_shifts_in_same_orbit(self):
        # Shifts 0 and 2 of the period-2 base ZIZI give the same product orbit.
        a, b = pe.PauliCycle(4, "XIII"), pe.PauliCycle(4, "ZIZI")
        s = a.commutator(b)
        for i in range(len(s.data)):
            for j in range(i + 1, len(s.data)):
                assert s.data[i] != s.data[j]
        assert s.to_qubit_hamiltonian() == reference_hamiltonian("XIII").commutator(
            reference_hamiltonian("ZIZI")
        )

    def test_constructor_merges_and_drops_zeros(self):
        s = pe.PauliCycleSum([
            pe.PauliCycle(4, "XYZI", 1.0),
            pe.PauliCycle(4, py_rotate("XYZI", 1), 2.0),  # same orbit -> 3.0
            pe.PauliCycle(4, "ZIII", 1.0),
            pe.PauliCycle(4, "IZII", -1.0),  # cancels the previous one
        ])
        assert len(s) == 1
        assert s.data[0].base.coeff == pytest.approx(3.0)

    def test_contains_returns_coefficient(self):
        s = pe.PauliCycleSum([pe.PauliCycle(4, "XYZI", 2.5)])
        # any rotation, with any coefficient, finds the summand
        assert s.contains(pe.PauliCycle(4, py_rotate("XYZI", 3), 7.0)) == pytest.approx(2.5)
        assert s.contains(pe.PauliCycle(4, "ZZZZ")) == 0
        assert pe.PauliCycle(4, "IXYZ") in s
        assert pe.PauliCycle(4, "ZZZZ") not in s


# A cycle is the sum over ALL n rotations, also for periodic bases


class TestPeriodicBases:
    def test_periodic_cycle_counts_multiplicity(self):
        assert pe.PauliCycle(3, "ZZZ").to_qubit_hamiltonian() == pe.QubitHamiltonian(
            [pe.PauliString(3.0, {0: "Z", 1: "Z", 2: "Z"})]
        )
        assert pe.PauliCycle(4, "XIXI").to_qubit_hamiltonian() == reference_hamiltonian("XIXI")

    @pytest.mark.parametrize(
        "base_a, base_b",
        [("XIXI", "ZYII"), ("ZZZZ", "XIII"), ("XYXY", "ZIZI"), ("XXX", "YII")],
    )
    def test_commutator_and_product_match_naive(self, base_a, base_b):
        n = len(base_a)
        a, b = pe.PauliCycle(n, base_a), pe.PauliCycle(n, base_b)
        ha, hb = reference_hamiltonian(base_a), reference_hamiltonian(base_b)
        assert a.commutator(b).to_qubit_hamiltonian() == ha.commutator(hb)
        assert a.multiply(b).to_qubit_hamiltonian() == ha * hb

    def test_targeted_on_periodic(self):
        a, b = pe.PauliCycle(4, "XIII"), pe.PauliCycle(4, "ZIZI")
        full = a.commutator(b)
        for term in full.data:
            assert a.targeted_commutator(b, term) == pytest.approx(complex(term.base.coeff))


# get_hash / from_hash


class TestHash:
    def test_rotation_invariant_and_coefficient_free(self):
        base = "XIZIIYII"
        keys = {pe.PauliCycle(8, py_rotate(base, k), complex(k + 1, -k)).get_hash() for k in range(8)}
        assert len(keys) == 1

    def test_consistent_with_eq(self):
        a, b = pe.PauliCycle(4, "XIIZ"), pe.PauliCycle(4, "XIZI")
        assert a != b and a.get_hash() != b.get_hash()
        # same base string, different ring size -> different cycles, different keys
        assert pe.PauliCycle(3, "XII").get_hash() != pe.PauliCycle(4, "XII").get_hash()

    def test_layout(self):
        # string key of the unified representative * 2**32 + n_qubits
        c = pe.PauliCycle(5, "IIXIZ")
        c.unify()
        assert c.get_hash() == c.base.get_hash() * 2**32 + 5
        assert pe.PauliCycle(7, "IIIIIII").get_hash() == 7  # identity

    @pytest.mark.parametrize(
        "n, base",
        [(1, "X"), (4, "XYZI"), (4, "ZIZI"), (8, "YIIXIZII"), (70, {0: "X", 33: "Y", 69: "Z"}), (300, {2: "Z", 299: "X"})],
    )
    def test_roundtrip(self, n, base):
        c = pe.PauliCycle(n, base, 2.5 - 1j)
        back = pe.PauliCycle.from_hash(c.get_hash(), 2.5 - 1j)
        assert back == c
        assert back.n_qubits == n
        assert back.base.coeff == 2.5 - 1j
        assert back.get_hash() == c.get_hash()
        assert back.to_qubit_hamiltonian() == c.to_qubit_hamiltonian()

    def test_from_hash_base_is_unified(self):
        c = pe.PauliCycle(8, "IIYIIXIZ")
        u = pe.PauliCycle(8, "IIYIIXIZ")
        u.unify()
        assert pe.PauliCycle.from_hash(c.get_hash()).base == u.base

    def test_from_hash_default_coeff(self):
        assert pe.PauliCycle.from_hash(pe.PauliCycle(3, "XYI", 4.0).get_hash()).base.coeff == 1

    def test_from_hash_invalid(self):
        with pytest.raises(ValueError):
            pe.PauliCycle.from_hash(-1)
        with pytest.raises(ValueError):
            pe.PauliCycle.from_hash(0)  # n_qubits = 0
        with pytest.raises(ValueError):
            # string acts on qubit 5 but n_qubits = 3
            pe.PauliCycle.from_hash(pe.PauliString(1.0, {5: "X"}).get_hash() * 2**32 + 3)
