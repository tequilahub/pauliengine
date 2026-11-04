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
            assert isinstance(got, _core.PauliOrbitSum), op
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
