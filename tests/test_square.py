"""h * h and h ** n use the square formula

    H^2 = (sum_i a_i^2) I + 2 sum_{i < j, [P_i, P_j] = 0} a_i a_j P_i P_j

instead of all pairs. The results must agree with the general product of h with a
copy of itself (a different object, so the general path is taken).
"""

import random

import pytest

import pauliengine as pe
import pauliengine._core as _core


def random_hamiltonian(n_terms, n_qubits, weight=None, seed=0, complex_coeffs=False):
    rng = random.Random(seed)
    terms = []
    for _ in range(n_terms):
        sites = range(n_qubits) if weight is None else rng.sample(range(n_qubits), weight)
        ops = {q: rng.choice("IXYZ" if weight is None else "XYZ") for q in sites}
        ops = {q: p for q, p in ops.items() if p != "I"}
        c = complex(rng.uniform(-1, 1), rng.uniform(-1, 1) if complex_coeffs else 0.0)
        terms.append(pe.PauliString(c, ops))
    return pe.QubitHamiltonian(terms)


def as_dict(h, tol=1e-12):
    return {tuple(sorted(ops.items())): complex(c) for c, ops in h.to_dictionary() if abs(complex(c)) > tol}


def assert_same(a, b):
    da, db = as_dict(a), as_dict(b)
    assert da.keys() == db.keys()
    for k in da:
        assert da[k] == pytest.approx(db[k], abs=1e-12, rel=1e-12)


def copy_of(h):
    return pe.QubitHamiltonian(h.paulistrings())


@pytest.mark.parametrize(
    ("n_terms", "n_qubits", "weight", "complex_coeffs"),
    [
        (1, 3, None, False),
        (2, 3, None, True),
        (40, 4, None, False),       # few qubits: many equal products
        (60, 64, None, True),
        (60, 64, 3, False),
        (80, 200, None, False),     # 4 words
        (30, 700, None, True),      # more than 512 qubits: the index-sorted path
    ],
)
def test_square_matches_general_product(n_terms, n_qubits, weight, complex_coeffs):
    h = random_hamiltonian(n_terms, n_qubits, weight, seed=n_terms + n_qubits, complex_coeffs=complex_coeffs)
    assert_same(h * h, h * copy_of(h))


def test_power_matches_repeated_product():
    h = random_hamiltonian(15, 6, seed=1, complex_coeffs=True)
    c = copy_of(h)
    assert_same(h ** 2, h * c)
    assert_same(h.power(3), (h * c) * c)
    assert_same(h ** 4, ((h * c) * c) * c)
    assert_same(h ** 1, h)
    assert as_dict(h ** 0) == {(): 1}


def test_anticommuting_terms_cancel_exactly():
    # (aX + bZ)^2 = (a^2 + b^2) I: the XZ and ZX products cancel
    h = pe.QubitHamiltonian([pe.PauliString(0.5, "X"), pe.PauliString(-2.0, "Z")])
    sq = h * h
    assert len(sq) == 1
    assert as_dict(sq, tol=0) == {(): pytest.approx(4.25)}


def test_commuting_terms_double():
    # (aZ0 + bZ1)^2 = (a^2 + b^2) I + 2ab Z0 Z1
    h = pe.QubitHamiltonian([pe.PauliString(3.0, "Z0"), pe.PauliString(0.5, "Z1")])
    assert as_dict(h * h) == {(): 9.25, ((0, "Z"), (1, "Z")): 3.0}


def test_empty_and_zero_terms():
    assert len(pe.QubitHamiltonian([]) ** 2) == 0
    h = random_hamiltonian(10, 5, seed=3) * 0     # zero coefficients are kept by the scaling
    assert len(h * h) == 0


def test_symbolic_square():
    h = pe.QubitHamiltonian([pe.PauliString("a", "X"), pe.PauliString("b", "Z"), pe.PauliString("c", "Z1")])
    sq = h * h
    assert isinstance(sq, _core.QubitHamiltonianSymbolic)
    general = h * pe.QubitHamiltonian(h.paulistrings())
    values = {"a": 0.3 + 0j, "b": -1.1 + 0j, "c": 2.0 + 0j}

    def evaluated(q):
        return {
            tuple(sorted(ops.items())): _core.PauliStringSymbolic.to_complex(c)
            for c, ops in q.subs(values).to_dictionary()
            if abs(_core.PauliStringSymbolic.to_complex(c)) > 1e-12
        }

    left, right = evaluated(sq), evaluated(general)
    assert left.keys() == right.keys()
    assert all(left[k] == pytest.approx(right[k]) for k in left)
    # no XZ term at all: the anti-commuting pair is never formed
    assert len(sq) == 3       # I, X Z1, Z Z1
