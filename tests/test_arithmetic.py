"""Arithmetic between all operator representations (include/pauliengine/Arithmetic.h).

Every combination of PauliString, QubitHamiltonian, PauliCycle, PauliCycleSum,
PauliOrbit and PauliOrbitSum must give the same operator as expanding both
operands into QubitHamiltonians first. Within the cycle and orbit families the
result stays compressed; everything else becomes a QubitHamiltonian.
"""

import itertools
import operator

import pytest

import pauliengine as pe
import pauliengine._core as _core

N = 4


def make_operands():
    cs = pe.PauliCycleSum([pe.PauliCycle(N, "ZIZI", 0.5), pe.PauliCycle(N, "XIII", -1j)])
    os_ = pe.PauliOrbitSum(N)
    os_.add(pe.PauliOrbit(N, 1, 0, 1), 2.0)
    os_.add(pe.PauliOrbit(N, 0, 0, 2), -0.5j)
    return {
        "PauliString": pe.PauliString(0.5 + 0.25j, {0: "X", 2: "Z"}),
        "QubitHamiltonian": pe.QubitHamiltonian(
            [pe.PauliString(1.0, {1: "Y"}), pe.PauliString(-0.3j, {0: "Z", 3: "X"})]
        ),
        "PauliCycle": pe.PauliCycle(N, "XYII", 1.5),
        "PauliCycleSum": cs,
        "PauliOrbit": pe.PauliOrbit(N, 1, 1, 0),
        "PauliOrbitSum": os_,
    }


def to_qh(x):
    """Any operator (or result) as a QubitHamiltonian."""
    if hasattr(x, "to_qubit_hamiltonian"):
        return x.to_qubit_hamiltonian()
    if isinstance(x, (_core.PauliStringComplex, _core.PauliStringSymbolic)):
        return pe.QubitHamiltonian([x])
    return x


def assert_same_operator(got, expected, tol=1e-12):
    diff = to_qh(got) - to_qh(expected)
    assert len(diff.simplify(tol)) == 0, f"difference:\n{diff}"


OPS = {
    "+": operator.add,
    "-": operator.sub,
    "*": operator.mul,
    "commutator": lambda a, b: a.commutator(b),
}
NAMES = list(make_operands())
CYCLES = {"PauliCycle", "PauliCycleSum"}
ORBITS = {"PauliOrbit", "PauliOrbitSum"}


@pytest.mark.parametrize("op", list(OPS))
@pytest.mark.parametrize(("left", "right"), list(itertools.product(NAMES, NAMES)))
def test_matches_expanded_operators(left, right, op):
    ops = make_operands()
    a, b = ops[left], ops[right]
    got = OPS[op](a, b)
    expected = OPS[op](to_qh(a), to_qh(b))
    assert_same_operator(got, expected)


@pytest.mark.parametrize(("left", "right"), list(itertools.product(NAMES, NAMES)))
def test_result_types(left, right):
    ops = make_operands()
    a, b = ops[left], ops[right]
    for op in OPS:
        got = OPS[op](a, b)
        if left in CYCLES and right in CYCLES:
            assert isinstance(got, _core.PauliCycleSumComplex), op
        elif left in ORBITS and right in ORBITS:
            assert isinstance(got, _core.PauliOrbitSumComplex), op
        elif left == right == "PauliString" and op in ("*", "commutator"):
            assert isinstance(got, _core.PauliStringComplex), op
        else:
            assert isinstance(got, _core.QubitHamiltonianComplex), op


class TestSymbolic:
    """A symbolic operand makes the result symbolic; values agree after substitution."""

    @staticmethod
    def as_numeric(qh):
        return {
            tuple(sorted(o.items())): complex(_core.PauliStringSymbolic.to_complex(c))
            for c, o in qh.to_dictionary()
        }

    @pytest.mark.parametrize("other", [n for n in NAMES if n != "PauliString"])
    @pytest.mark.parametrize("op", list(OPS))
    def test_symbolic_string_with_every_type(self, other, op):
        ops = make_operands()
        sym = pe.PauliString("a", {0: "X", 2: "Z"})
        num = pe.PauliString(2.0, {0: "X", 2: "Z"})
        b = ops[other]
        for got, expected in ((OPS[op](sym, b), OPS[op](num, b)), (OPS[op](b, sym), OPS[op](b, num))):
            assert isinstance(got, _core.QubitHamiltonianSymbolic)
            values = self.as_numeric(got.subs({"a": 2.0}))
            ref = {tuple(sorted(o.items())): complex(c) for c, o in to_qh(expected).to_dictionary()}
            assert values.keys() == ref.keys()
            for k in ref:
                assert values[k] == pytest.approx(ref[k], abs=1e-12)

    def test_symbolic_and_complex_hamiltonians(self):
        hs = pe.QubitHamiltonian([("a", {0: "X"})])
        hc = pe.QubitHamiltonian([(1.0, {0: "Y"})])
        assert isinstance(hs + hc, _core.QubitHamiltonianSymbolic)
        assert isinstance(hc * hs, _core.QubitHamiltonianSymbolic)
        assert isinstance(hc.commutator(hs), _core.QubitHamiltonianSymbolic)


class TestFamilies:
    def test_cycle_sum_merges_orbits(self):
        a, b = pe.PauliCycle(N, "XYII", 1.0), pe.PauliCycle(N, "IXYI", 2.0)  # same orbit
        s = a + b
        assert len(s) == 1
        assert s.contains(a) == pytest.approx(3.0)
        assert len(a - pe.PauliCycle(N, "IXYI", 1.0)) == 0

    def test_orbit_sum(self):
        a, b = pe.PauliOrbit(N, 1, 0, 0), pe.PauliOrbit(N, 0, 1, 0)
        s = a + b - a
        assert s.contains(b) == pytest.approx(1.0)
        assert s.contains(a) == 0

    def test_n_qubits_mismatch_in_family_raises(self):
        with pytest.raises(ValueError):
            pe.PauliCycle(3, "XII") + pe.PauliCycle(4, "XIII")
        with pytest.raises(ValueError):
            pe.PauliOrbit(3, 1, 0, 0) + pe.PauliOrbit(4, 1, 0, 0)

    def test_mixed_families_with_different_sizes(self):
        # no size constraint across kinds: both are simply operators
        h = pe.PauliCycle(3, "XII") + pe.PauliOrbit(5, 1, 0, 0)
        assert_same_operator(h, pe.PauliCycle(3, "XII").to_qubit_hamiltonian()
                             + pe.PauliOrbit(5, 1, 0, 0).to_qubit_hamiltonian())


class TestOrbitExpansion:
    def test_normalization(self):
        h = pe.PauliOrbit(4, 1, 1, 0).to_qubit_hamiltonian()
        assert len(h) == 12
        assert all(complex(c) == pytest.approx(1j / 12) for c, _ in h.to_dictionary())

    def test_identity_orbit(self):
        h = pe.PauliOrbit(3, 0, 0, 0).to_qubit_hamiltonian()
        assert len(h) == 1

    def test_too_large_expansion_raises(self):
        big = pe.PauliOrbit(170, 4, 4, 4)
        with pytest.raises((ValueError, RuntimeError)):
            big.to_qubit_hamiltonian()
        with pytest.raises((ValueError, RuntimeError)):
            big + pe.PauliString(1.0, {0: "X"})


class TestConsistentTypes:
    """The same kind of operation gives the same kind of result for every operand."""

    def test_mixed_pauli_string_commutator_is_a_pauli_string(self):
        x = pe.PauliString(1.0, {0: "X"})
        y = _core.PauliStringSymbolic(_core.Expression("a"), {0: "Y"})
        for c in (x.commutator(y), y.commutator(x), y.commutator(pe.PauliString("X"))):
            assert isinstance(c, _core.PauliStringSymbolic)
        assert x.commutator(y) == _core.PauliStringSymbolic(_core.Expression("2*I*a"), {0: "Z"})

    @pytest.mark.parametrize("make", [
        lambda: pe.PauliCycle(4, "XZII", 0.5),
        lambda: pe.PauliCycleSum([pe.PauliCycle(4, "XIII", 0.5)]),
        lambda: pe.PauliOrbit(4, 1, 0, 1),
        lambda: pe.PauliOrbitSum(4, [(pe.PauliOrbit(4, 1, 0, 1), 0.5)]),
    ])
    def test_complex_scaling_and_negation(self, make):
        op = make()
        expected = op.to_qubit_hamiltonian() * 2j
        assert (op * 2j).to_qubit_hamiltonian() == expected
        assert (2j * op).to_qubit_hamiltonian() == expected
        assert (-op).to_qubit_hamiltonian() == op.to_qubit_hamiltonian() * -1

    def test_equality_across_coefficient_types(self):
        x = pe.PauliString(1.0, {0: "X"})
        xs = _core.PauliStringSymbolic(_core.Expression("1"), {0: "X"})
        assert x == xs and xs == x and hash(x) == hash(xs)
        assert x != xs * 2
        assert pe.QubitHamiltonian([x]) == _core.QubitHamiltonianSymbolic([xs])
        assert type(x == x) is bool

    def test_equality_with_other_objects_is_false(self):
        x = pe.PauliString(1.0, {0: "X"})
        for other in (1, "X", None, pe.QubitHamiltonian([x]), pe.PauliOrbit(4, 1, 0, 0)):
            assert not (x == other)
            assert x != other

    def test_sums_compare_by_value(self):
        c = pe.PauliCycle(4, "XIII", 0.5)
        assert pe.PauliCycleSum([c]) == pe.PauliCycleSum([c])
        assert pe.PauliCycleSum([c]) != pe.PauliCycleSum([c * 2])
        o = pe.PauliOrbit(4, 1, 0, 0)
        assert pe.PauliOrbitSum(4, [(o, 1.0)]) == pe.PauliOrbitSum(4, [(o, 1.0)])
        assert pe.PauliOrbitSum(4, [(o, 1.0)]) != pe.PauliOrbitSum(4, [(o, 2.0)])
        with pytest.raises(TypeError):
            hash(pe.PauliCycleSum([c]))


ALL_KINDS = [
    lambda: pe.PauliString(1.5, {0: "X", 1: "Z"}),
    lambda: pe.QubitHamiltonian([pe.PauliString(1.5, {0: "X"}), pe.PauliString(0.5, {1: "Y"})]),
    lambda: pe.PauliCycle(4, "XZII", 0.5),
    lambda: pe.PauliCycleSum([pe.PauliCycle(4, "XIII", 0.5), pe.PauliCycle(4, "YZII", 2.0)]),
    lambda: pe.PauliOrbit(4, 1, 0, 1),
    lambda: pe.PauliOrbitSum(4, [(pe.PauliOrbit(4, 1, 0, 1), 0.5)]),
    lambda: _core.unweighted_pauliString("XZ"),
    lambda: _core.unweighted_QubitHamiltonian([_core.unweighted_pauliString("XZ"), _core.unweighted_pauliString("YI")]),
]


def expanded(op):
    """The operator as a QubitHamiltonian, symbolic ones evaluated at b = 0.7."""
    if isinstance(op, (_core.unweighted_pauliString, _core.unweighted_QubitHamiltonian)):
        op = op * 1.0
    if not isinstance(op, (_core.QubitHamiltonianComplex, _core.QubitHamiltonianSymbolic)):
        op = op.to_qubit_hamiltonian() if hasattr(op, "to_qubit_hamiltonian") else pe.QubitHamiltonian([op])
    if isinstance(op, _core.QubitHamiltonianSymbolic):
        op = op.subs({"b": 0.7})
        op = pe.QubitHamiltonian([pe.PauliString(_core.PauliStringSymbolic.to_complex(c), d) for c, d in op.to_dictionary()])
    return op


class TestScalars:
    @pytest.mark.parametrize("make", ALL_KINDS)
    @pytest.mark.parametrize("scale", [2, 2.5, 2j])
    def test_division_is_multiplication_by_the_inverse(self, make, scale):
        op = make()
        assert expanded(op / scale) == expanded(op) * (1 / scale)

    @pytest.mark.parametrize("make", ALL_KINDS)
    def test_division_by_zero(self, make):
        with pytest.raises(ZeroDivisionError):
            make() / 0

    @pytest.mark.parametrize("make", ALL_KINDS)
    def test_symbolic_scale_on_both_sides(self, make):
        op, b = make(), _core.Expression("b")
        assert expanded(op * b) == expanded(op) * 0.7
        assert expanded(b * op) == expanded(op) * 0.7

    def test_symbolic_scale_by_zero_drops_terms(self):
        h = pe.QubitHamiltonian([pe.PauliString(1.0, {0: "X"})])
        assert len(h * _core.Expression("0")) == 0


class TestUnweightedCompletions:
    def test_inner_product_between_unweighted_operators(self):
        s = _core.unweighted_pauliString("XZ")
        t = _core.unweighted_pauliString("YI")
        h = _core.unweighted_QubitHamiltonian([s, t])
        for a, b in itertools.product((s, t, h), repeat=2):
            assert a.inner_product(b) == pytest.approx(expanded(a).inner_product(expanded(b)))

    def test_negated_set(self):
        h = _core.unweighted_QubitHamiltonian([_core.unweighted_pauliString("XZ"), _core.unweighted_pauliString("YI")])
        assert isinstance(-h, _core.unweighted_QubitHamiltonian)
        assert expanded(-h) == expanded(h) * -1
