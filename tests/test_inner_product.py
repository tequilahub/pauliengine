"""Normalized Hilbert-Schmidt inner product <A, B> = Tr(A^dagger B) / 2^n
(arXiv:2604.16701, Corollary 1 and Lemma 1) between all representations."""

import itertools
import math

import pytest

import pauliengine as pe
import pauliengine._core as _core

np = pytest.importorskip("numpy")

from test_arithmetic import NAMES, make_operands, to_qh  # noqa: E402


def dense(x, n_qubits):
    """Matrix of any operator on n_qubits qubits (qubit 0 = most significant)."""
    h = to_qh(x)
    if len(h) == 0:
        return np.zeros((2**n_qubits, 2**n_qubits), dtype=complex)
    m = np.array(h.to_matrix(ignore_unused_qubits=False))
    used = int(round(math.log2(m.shape[0])))
    return np.kron(m, np.eye(2 ** (n_qubits - used)))


def reference(a, b, n_qubits):
    A, B = dense(a, n_qubits), dense(b, n_qubits)
    return np.trace(A.conj().T @ B) / 2**n_qubits


@pytest.mark.parametrize(("left", "right"), list(itertools.product(NAMES, NAMES)))
def test_matches_dense_trace(left, right):
    ops = make_operands()
    a, b = ops[left], ops[right]
    assert complex(a.inner_product(b)) == pytest.approx(reference(a, b, 4), abs=1e-12)
    assert complex(pe.inner_product(a, b)) == pytest.approx(complex(a.inner_product(b)))


@pytest.mark.parametrize(("left", "right"), list(itertools.product(NAMES, NAMES)))
def test_conjugate_symmetry(left, right):
    ops = make_operands()
    a, b = ops[left], ops[right]
    assert complex(a.inner_product(b)) == pytest.approx(complex(b.inner_product(a)).conjugate(), abs=1e-12)


@pytest.mark.parametrize("name", NAMES)
def test_norm_is_positive(name):
    x = make_operands()[name]
    v = complex(x.inner_product(x))
    assert v.real > 0 and v.imag == pytest.approx(0, abs=1e-14)


class TestPaperFormulas:
    def test_corollary_1_pauli_strings(self):
        p, q = pe.PauliString(2 - 1j, {0: "X", 3: "Z"}), pe.PauliString(0.5j, {0: "X", 3: "Z"})
        assert complex(p.inner_product(q)) == pytest.approx((2 - 1j).conjugate() * 0.5j)
        assert complex(p.inner_product(pe.PauliString(1.0, {0: "Y", 3: "Z"}))) == 0
        # independent of the register the strings live in
        assert complex(pe.PauliString(1.0, {0: "X"}).inner_product(pe.PauliString(1.0, {0: "X"}))) == 1

    @pytest.mark.parametrize("label", [(1, 0, 0), (1, 1, 0), (2, 1, 1), (0, 0, 0)])
    def test_lemma_1_orbits(self, label):
        n = 6
        o = pe.PauliOrbit(n, *label)
        assert complex(o.inner_product(o)) == pytest.approx(1 / o.n_terms())   # 2^n / N_terms, normalized
        other = pe.PauliOrbit(n, label[1], label[0], label[2] + (1 if sum(label) < n else 0))
        if (other.p, other.q, other.r) != label:
            assert complex(o.inner_product(other)) == 0

    def test_orbits_on_different_registers(self):
        a, b = pe.PauliOrbit(3, 1, 0, 0), pe.PauliOrbit(5, 1, 0, 0)
        assert complex(a.inner_product(b)) == pytest.approx(reference(a, b, 5), abs=1e-14)

    def test_large_orbit_is_not_expanded(self):
        big = pe.PauliOrbit(170, 4, 4, 4)                     # ~1e20 strings
        assert complex(big.inner_product(big)) == pytest.approx(1 / big.n_terms(), rel=1e-12)
        s = pe.PauliString(3.0, {0: "X", 1: "X", 2: "X", 3: "X", 10: "Y", 11: "Y", 12: "Y", 13: "Y",
                                 20: "Z", 21: "Z", 22: "Z", 23: "Z"})
        assert complex(big.inner_product(s)) == pytest.approx(-1j / big.n_terms() * 3.0, rel=1e-12)

    @pytest.mark.parametrize(("base", "period"), [("XYZI", 4), ("ZIZI", 2), ("ZZZZ", 1), ("IIII", 1),
                                                  ("XYXYXY", 2), ("XIIXII", 3)])
    def test_cycle_period_and_norm(self, base, period):
        n = len(base)
        c = pe.PauliCycle(n, base, 1.5)
        assert c.period() == period
        # <C(P), C(P)> = |c|^2 n^2 / d
        assert complex(c.inner_product(c)) == pytest.approx(1.5**2 * n * n / period)
        assert complex(c.inner_product(c)) == pytest.approx(reference(c, c, n))

    def test_cycles_of_different_orbits_are_orthogonal(self):
        assert complex(pe.PauliCycle(4, "XZII").inner_product(pe.PauliCycle(4, "XIZI"))) == 0

    def test_cycles_on_different_registers(self):
        a, b = pe.PauliCycle(3, "ZII"), pe.PauliCycle(4, "ZIII")
        assert complex(a.inner_product(b)) == pytest.approx(reference(a, b, 4))


def test_linearity_in_second_argument():
    ops = make_operands()
    a, b, c = ops["PauliCycle"], ops["QubitHamiltonian"], ops["PauliOrbitSum"]
    lhs = complex(a.inner_product(b + c))
    assert lhs == pytest.approx(complex(a.inner_product(b)) + complex(a.inner_product(c)), abs=1e-12)


def test_symbolic():
    s = pe.PauliString("a", {0: "X", 2: "Z"})
    c = pe.PauliCycle(4, "XIZI", 2.0)
    v = c.inner_product(s)
    assert isinstance(v, _core.Expression)
    num = complex(_core.PauliStringSymbolic.to_complex(
        pe.QubitHamiltonian([pe.PauliString(str(v), {0: "X"})]).subs({"a": 3.0}).to_dictionary()[0][0]))
    assert num == pytest.approx(complex(c.inner_product(pe.PauliString(3.0, {0: "X", 2: "Z"}))))
