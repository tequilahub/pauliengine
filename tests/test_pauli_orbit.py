"""Tests for PauliOrbit / PauliOrbitSum (permutation-invariant basis).

The commutator is verified against the ground-truth operator algebra: each orbit
B_{p,q,r} is expanded into an explicit QubitHamiltonian (i / N_terms times the sum
of all placements of p X's, q Y's and r Z's), the two Hamiltonians are commuted,
and the result must match the expansion of ``PauliOrbit.commutator``
(arXiv:2604.16701, Theorem 1 / Prop. 7). A tolerance-based comparison is used
because the orbit coefficients are rationals such as 1/12.
"""

from itertools import permutations

import pytest

import pauliengine as pe


# Helpers


def orbit_to_qh(n: int, p: int, q: int, r: int):
    """B_{p,q,r} as an explicit QubitHamiltonian: (i / N_terms) * sum of strings."""
    letters = ["X"] * p + ["Y"] * q + ["Z"] * r + ["I"] * (n - p - q - r)
    perms = set(permutations(letters))
    n_terms = len(perms)
    terms = [
        pe.PauliString(1j / n_terms, {i: c for i, c in enumerate(perm) if c != "I"} or {0: "I"})
        for perm in perms
    ]
    return pe.QubitHamiltonian(terms)


def orbitsum_to_qh(s):
    """Expand a PauliOrbitSum into a QubitHamiltonian."""
    acc = pe.QubitHamiltonian.zero()
    for orbit, coeff in zip(s.data, s.coeffs):
        acc = acc + orbit_to_qh(s.n_qubits, orbit.p, orbit.q, orbit.r) * complex(coeff)
    return acc


def qh_close(a, b, tol: float = 1e-9) -> bool:
    """Coefficient-wise comparison of two QubitHamiltonians within a tolerance."""

    def as_dict(qh):
        return {tuple(sorted(ops.items())): complex(c) for c, ops in qh.to_dictionary()}

    da, db = as_dict(a), as_dict(b)
    return all(abs(da.get(k, 0j) - db.get(k, 0j)) < tol for k in set(da) | set(db))


# Data structures


class TestPauliOrbit:
    def test_n_terms_and_identities(self):
        orb = pe.PauliOrbit(4, 1, 1, 0)  # 4!/(1!1!0!2!) = 12 strings, 2 identities
        assert orb.n_terms() == pytest.approx(12.0)
        assert orb.identities() == 2

    def test_validity_and_equality(self):
        assert pe.PauliOrbit(3, 1, 1, 1).is_valid()
        assert not pe.PauliOrbit(3, 2, 2, 0).is_valid()  # p + q + r > n
        assert pe.PauliOrbit(4, 1, 0, 2) == pe.PauliOrbit(4, 1, 0, 2)
        assert pe.PauliOrbit(4, 1, 0, 2) != pe.PauliOrbit(4, 0, 1, 2)
        assert len({pe.PauliOrbit(4, 1, 0, 2), pe.PauliOrbit(4, 1, 0, 2)}) == 1


class TestPauliOrbitSum:
    def test_add_merges_labels_and_prune(self):
        s = pe.PauliOrbitSum(4)
        s.add(pe.PauliOrbit(4, 1, 0, 0), 1.5)
        s.add(pe.PauliOrbit(4, 1, 0, 0), -0.5)  # same label -> coefficients sum
        s.add(pe.PauliOrbit(4, 0, 1, 0), 2.0)
        assert len(s) == 2
        merged = {(o.p, o.q, o.r): c for o, c in zip(s.data, s.coeffs)}
        assert merged[(1, 0, 0)] == pytest.approx(1.0)
        s.add(pe.PauliOrbit(4, 0, 1, 0), -2.0)  # cancels to zero
        s.prune_zeros()
        assert len(s) == 1


# Commutator against the expanded operator algebra


class TestCommutator:
    @pytest.mark.parametrize(
        "n, a, b",
        [
            (3, (1, 0, 0), (0, 0, 1)),
            (3, (1, 0, 0), (0, 1, 0)),
            (4, (1, 1, 0), (0, 1, 1)),
            (5, (1, 1, 1), (1, 0, 1)),
        ],
    )
    def test_matches_expanded_hamiltonian(self, n, a, b):
        oa, ob = pe.PauliOrbit(n, *a), pe.PauliOrbit(n, *b)
        reference = orbit_to_qh(n, *a).commutator(orbit_to_qh(n, *b))
        candidate = orbitsum_to_qh(oa.commutator(ob))
        assert qh_close(reference, candidate)

    def test_antisymmetry(self):
        n = 4
        a, b = pe.PauliOrbit(n, 1, 1, 0), pe.PauliOrbit(n, 0, 1, 1)
        forward = {(o.p, o.q, o.r): c for o, c in zip(a.commutator(b).data, a.commutator(b).coeffs)}
        backward = {(o.p, o.q, o.r): -c for o, c in zip(b.commutator(a).data, b.commutator(a).coeffs)}
        assert forward.keys() == backward.keys()
        assert all(forward[k] == pytest.approx(backward[k]) for k in forward)


class TestOrbitSumCommutator:
    def test_matches_expanded_hamiltonian(self):
        n = 5
        a = pe.PauliOrbitSum(n)
        a.add(pe.PauliOrbit(n, 1, 0, 0), 2.0)
        a.add(pe.PauliOrbit(n, 0, 1, 0), -1.5)
        b = pe.PauliOrbitSum(n)
        b.add(pe.PauliOrbit(n, 0, 0, 1), 0.5)
        b.add(pe.PauliOrbit(n, 1, 1, 0), 1.0)
        reference = orbitsum_to_qh(a).commutator(orbitsum_to_qh(b))
        candidate = orbitsum_to_qh(a.commutator(b))
        assert qh_close(reference, candidate)


class TestTargetedCommutator:
    @pytest.mark.parametrize(
        "n, a, b",
        [(3, (1, 0, 0), (0, 0, 1)), (4, (1, 1, 0), (0, 1, 1)), (5, (1, 1, 1), (1, 0, 1))],
    )
    def test_orbit_targeted_matches_full(self, n, a, b):
        oa, ob = pe.PauliOrbit(n, *a), pe.PauliOrbit(n, *b)
        full = oa.commutator(ob)
        # Every output orbit's targeted structure constant matches the full result.
        for orbit, c in zip(full.data, full.coeffs):
            target = pe.PauliOrbit(n, orbit.p, orbit.q, orbit.r)
            assert oa.targeted_commutator(ob, target) == pytest.approx(c)
        # The identity orbit is never produced (parity D must be odd) -> 0.
        assert oa.targeted_commutator(ob, pe.PauliOrbit(n, 0, 0, 0)) == pytest.approx(0.0)

    def test_orbit_sum_targeted_matches_full(self):
        n = 5
        a = pe.PauliOrbitSum(n)
        a.add(pe.PauliOrbit(n, 1, 0, 0), 2.0)
        a.add(pe.PauliOrbit(n, 0, 1, 0), -1.5)
        b = pe.PauliOrbitSum(n)
        b.add(pe.PauliOrbit(n, 0, 0, 1), 0.5)
        b.add(pe.PauliOrbit(n, 1, 1, 0), 1.0)
        full = a.commutator(b)
        for orbit, c in zip(full.data, full.coeffs):
            target = pe.PauliOrbit(n, orbit.p, orbit.q, orbit.r)
            assert a.targeted_commutator(b, target) == pytest.approx(c)


class TestMultiply:
    @pytest.mark.parametrize(
        "n, a, b",
        [
            (2, (1, 0, 0), (0, 1, 0)),
            (3, (1, 0, 0), (0, 0, 1)),
            (2, (1, 0, 0), (1, 0, 0)),   # produces the identity orbit B_{0,0,0}
            (4, (1, 1, 0), (0, 1, 1)),
            (5, (1, 1, 1), (1, 0, 1)),
        ],
    )
    def test_orbit_multiply_matches_expanded(self, n, a, b):
        oa, ob = pe.PauliOrbit(n, *a), pe.PauliOrbit(n, *b)
        reference = orbit_to_qh(n, *a) * orbit_to_qh(n, *b)
        assert qh_close(orbitsum_to_qh(oa.multiply(ob)), reference)
        assert qh_close(orbitsum_to_qh(oa * ob), reference)  # __mul__ agrees

    def test_orbit_sum_multiply_matches_expanded(self):
        n = 5
        a = pe.PauliOrbitSum(n)
        a.add(pe.PauliOrbit(n, 1, 0, 0), 2.0)
        a.add(pe.PauliOrbit(n, 0, 1, 0), 1j)
        b = pe.PauliOrbitSum(n)
        b.add(pe.PauliOrbit(n, 0, 0, 1), 0.5)
        b.add(pe.PauliOrbit(n, 1, 1, 0), 1.0)
        reference = orbitsum_to_qh(a) * orbitsum_to_qh(b)
        assert qh_close(orbitsum_to_qh(a.multiply(b)), reference)
