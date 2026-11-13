"""Pauli cycles and orbit sums with symbolic coefficients.

Every symbolic result is checked against the complex classes: substituting numbers
for the symbols must give the same operator.
"""

import pytest

import pauliengine as pe
import pauliengine._core as _core

A, B = 0.7, -1.3
VALUES = {"a": complex(A), "b": complex(B)}


def as_dict(qh, values=None):
    """{operator: coefficient} of a Hamiltonian; symbolic ones are evaluated first."""
    if isinstance(qh, _core.QubitHamiltonianSymbolic):
        qh = qh.subs(values or VALUES)
        items = [(_core.PauliStringSymbolic.to_complex(c), ops) for c, ops in qh.to_dictionary()]
    else:
        items = [(complex(c), ops) for c, ops in qh.to_dictionary()]
    return {tuple(sorted(ops.items())): c for c, ops in items if abs(c) > 1e-12}


def same_operator(sym, num):
    """sym (symbolic, any kind) evaluated at VALUES equals num (complex)."""
    left = as_dict(sym.to_qubit_hamiltonian() if hasattr(sym, "to_qubit_hamiltonian") else sym)
    right = as_dict(num.to_qubit_hamiltonian() if hasattr(num, "to_qubit_hamiltonian") else num)
    return left.keys() == right.keys() and all(abs(left[k] - right[k]) < 1e-9 for k in left)


def value(coeff):
    """A symbolic coefficient evaluated at VALUES."""
    h = pe.QubitHamiltonian([pe.PauliString(str(coeff), {})]).subs(VALUES)
    terms = h.to_dictionary()
    return _core.PauliStringSymbolic.to_complex(terms[0][0]) if terms else 0j


N = 6


@pytest.fixture
def cycles():
    sym = (pe.PauliCycle(N, "XZI", "a"), pe.PauliCycle(N, "YIX", "b"))
    num = (pe.PauliCycle(N, "XZI", A), pe.PauliCycle(N, "YIX", B))
    return sym, num


class TestCycleFactory:
    def test_symbolic_coefficient_gives_symbolic_cycle(self):
        assert isinstance(pe.PauliCycle(4, "XZ", "a"), _core.PauliCycleSymbolic)
        assert isinstance(pe.PauliCycle(4, {0: "X"}, coeff="a*b"), _core.PauliCycleSymbolic)
        assert isinstance(pe.PauliCycle(4, pe.PauliString("a", "XZ")), _core.PauliCycleSymbolic)

    def test_numeric_coefficient_gives_complex_cycle(self):
        assert isinstance(pe.PauliCycle(4, "XZ"), _core.PauliCycleComplex)
        assert isinstance(pe.PauliCycle(4, "XZ", 2), _core.PauliCycleComplex)
        assert isinstance(pe.PauliCycle(4, pe.PauliString(1.0, "XZ")), _core.PauliCycleComplex)

    def test_isinstance_accepts_both(self):
        assert isinstance(pe.PauliCycle(4, "XZ"), pe.PauliCycle)
        assert isinstance(pe.PauliCycle(4, "XZ", "a"), pe.PauliCycle)
        assert not isinstance(pe.PauliOrbit(4, 1, 0, 0), pe.PauliCycle)

    def test_sum_of_mixed_cycles_is_symbolic(self):
        s = pe.PauliCycleSum([pe.PauliCycle(4, "XZ", 2.0), pe.PauliCycle(4, "IXZ", "a")])
        assert isinstance(s, _core.PauliCycleSumSymbolic)
        assert isinstance(s, pe.PauliCycleSum)
        assert len(s) == 1          # same rotation orbit: merged to (2 + a)
        assert value(s.contains(pe.PauliCycle(4, "XZ"))) == pytest.approx(2.0 + A)

    def test_from_hash(self):
        c = pe.PauliCycle(5, "XYI", "a")
        back = pe.PauliCycle.from_hash(c.get_hash(), "a")
        assert isinstance(back, _core.PauliCycleSymbolic)
        assert back == c
        assert isinstance(pe.PauliCycle.from_hash(c.get_hash()), _core.PauliCycleComplex)


class TestCycleAlgebra:
    def test_product(self, cycles):
        (x, y), (xc, yc) = cycles
        assert isinstance(x * y, _core.PauliCycleSumSymbolic)
        assert same_operator(x * y, xc * yc)

    def test_commutator(self, cycles):
        (x, y), (xc, yc) = cycles
        assert same_operator(x.commutator(y), xc.commutator(yc))

    def test_sums(self, cycles):
        (x, y), (xc, yc) = cycles
        s, sc = x + y + x, xc + yc + xc
        assert same_operator(s, sc)
        assert same_operator(s * s, sc * sc)
        assert same_operator(s.commutator(x - y), sc.commutator(xc - yc))
        assert len(x - x) == 0

    def test_mixed_coefficients_promote(self, cycles):
        (x, y), (xc, yc) = cycles
        assert isinstance(x + yc, _core.PauliCycleSumSymbolic)
        assert isinstance(xc * y, _core.PauliCycleSumSymbolic)
        assert same_operator(xc * y, xc * yc)
        assert same_operator(xc.commutator(y), xc.commutator(yc))

    def test_scaling(self, cycles):
        (x, _), (xc, _) = cycles
        assert same_operator(x * 3, xc * 3)
        assert same_operator(2.5 * (x + x), 2.5 * (xc + xc))
        assert len((x + x) * 0) == 0

    def test_targeted_commutator(self, cycles):
        (x, y), (xc, yc) = cycles
        target = pe.PauliCycle(N, "ZZX")
        assert value(x.targeted_commutator(y, target)) == pytest.approx(xc.targeted_commutator(yc, target))

    def test_with_other_operators(self, cycles):
        (x, y), (xc, yc) = cycles
        h = yc.to_qubit_hamiltonian()
        assert isinstance(x * h, _core.QubitHamiltonianSymbolic)
        assert same_operator(x * h, xc * h)
        assert same_operator(x + pe.PauliString(1.0, "XX"), xc + pe.PauliString(1.0, "XX"))

    def test_inner_product(self, cycles):
        (x, y), (xc, yc) = cycles
        # symbolic on the right: no conjugate, so the value can be compared
        assert value(pe.inner_product(xc + yc, x + y)) == pytest.approx(pe.inner_product(xc + yc, xc + yc))
        assert value(pe.inner_product(xc, x)) == pytest.approx(pe.inner_product(xc, xc))
        # symbolic on the left: conjugate-linear
        assert "conjugate(b)" in str(pe.inner_product(x + y, y))


@pytest.fixture
def orbit_sums():
    o1, o2 = pe.PauliOrbit(5, 1, 1, 0), pe.PauliOrbit(5, 0, 1, 1)
    sym = pe.PauliOrbitSum(5, [(o1, "a"), (o2, "b"), (o1, 2)])
    num = pe.PauliOrbitSum(5, [(o1, A + 2), (o2, B)])
    return sym, num


class TestOrbitSum:
    def test_factory(self, orbit_sums):
        sym, num = orbit_sums
        assert isinstance(sym, _core.PauliOrbitSumSymbolic)
        assert isinstance(num, _core.PauliOrbitSumComplex)
        assert isinstance(sym, pe.PauliOrbitSum) and isinstance(num, pe.PauliOrbitSum)
        assert isinstance(pe.PauliOrbitSum(4), _core.PauliOrbitSumComplex)
        assert isinstance(pe.PauliOrbitSum(4, symbolic=True), _core.PauliOrbitSumSymbolic)
        assert len(sym) == 2

    def test_add_accepts_text_and_numbers(self):
        s = pe.PauliOrbitSum(4, symbolic=True)
        o = pe.PauliOrbit(4, 1, 0, 0)
        s.add(o, "a")
        s.add(o, 1.5)
        assert value(s.contains(o)) == pytest.approx(A + 1.5)

    def test_algebra(self, orbit_sums):
        sym, num = orbit_sums
        assert same_operator(sym, num)
        assert same_operator(sym.commutator(sym * 3), num.commutator(num * 3))
        assert same_operator(sym * sym, num * num)
        o3 = pe.PauliOrbit(5, 2, 0, 0)
        assert same_operator(sym + o3, num + o3)
        assert same_operator(o3 * sym, o3 * num)
        assert isinstance(num.commutator(sym), _core.PauliOrbitSumSymbolic)
        assert same_operator(num.commutator(sym), num.commutator(num))
        assert len(sym - sym) == 0

    def test_targeted_commutator(self, orbit_sums):
        sym, num = orbit_sums
        target = pe.PauliOrbit(5, 1, 0, 1)
        assert value(sym.targeted_commutator(sym, target)) == pytest.approx(num.targeted_commutator(num, target))

    def test_inner_product(self, orbit_sums):
        sym, num = orbit_sums
        assert value(pe.inner_product(num, sym)) == pytest.approx(pe.inner_product(num, num))
        h = num.to_qubit_hamiltonian()
        assert value(pe.inner_product(h, sym)) == pytest.approx(pe.inner_product(h, num))
        assert "conjugate(2.0 + a)" in str(pe.inner_product(sym, h))

    def test_with_cycles(self, orbit_sums, cycles):
        sym, num = orbit_sums
        (x, _), (xc, _) = cycles
        assert isinstance(sym * x, _core.QubitHamiltonianSymbolic)
        assert same_operator(sym * x, num * xc)
