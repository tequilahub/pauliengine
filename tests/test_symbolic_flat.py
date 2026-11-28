"""Symbolic sums on flat storage: every algorithm against the same operation with numbers.

Each case builds symbolic sums, applies an operation, substitutes values and
compares with the operation on the substituted (complex) sums. The cases cover
the packed product records, grouping by operator (monomials beyond 64 bits),
strings beyond 512 qubits (wide products and merge), squares, sums and maps.
"""

import os
import random
import subprocess
import sys

import pytest

import pauliengine as pe
import pauliengine._core as core

PSS = core.PauliStringSymbolic
QHS = core.PauliStringSumSymbolic
E = core.Expression


def random_sum(rng, n_terms, n_qubits, weight, coeffs):
    terms = []
    for _ in range(n_terms):
        ops = {q: rng.choice("XYZ") for q in rng.sample(range(n_qubits), weight)}
        terms.append(PSS(E(rng.choice(coeffs)), ops))
    return QHS(terms)


SIMPLE = ["a", "2*b", "a*b - 1", "(1 + I)*c", "3"]
# 20 symbols with exponents up to 900: monomials need far more than 64 bits
HEAVY = [f"s{k}**{rng_e}" for k, rng_e in zip(range(20), range(300, 920, 31))] + ["s0*s1", "2"]
VALUES = {"a": 0.7, "b": -1.3 + 0.2j, "c": 0.4j}
HEAVY_VALUES = {f"s{k}": 1.0 + 0.0005 * k for k in range(20)}


def as_dict(h):
    """{operator: complex} of a sum whose coefficients are numbers."""
    out = {}
    for c, d in h.to_dictionary():
        out[tuple(sorted(d.items()))] = PSS.to_complex(c) if isinstance(c, E) else complex(c)
    return out


def assert_close(sym, values, expected):
    got = as_dict(sym.subs(values))
    want = as_dict(expected)
    scale = max([abs(v) for v in want.values()] + [1.0])
    for k in got.keys() | want.keys():
        assert abs(got.get(k, 0) - want.get(k, 0)) <= 1e-9 * scale, k


CASES = [
    # (qubits, weight, coefficients, values): packed records, grouping, wide strings
    (12, 3, SIMPLE, VALUES),
    (12, 3, HEAVY, HEAVY_VALUES),
    (700, 5, SIMPLE, VALUES),
    (700, 5, HEAVY, HEAVY_VALUES),
]


@pytest.fixture(params=range(len(CASES)), ids=["packed", "grouped", "wide", "wide-heavy"])
def case(request):
    n, w, coeffs, values = CASES[request.param]
    rng = random.Random(100 + request.param)
    return random_sum(rng, 40, n, w, coeffs), random_sum(rng, 35, n, w, coeffs), values


def test_product_commutator_square(case):
    a, b, values = case
    an, bn = a.subs(values), b.subs(values)
    assert_close(a * b, values, an * bn)
    assert_close(a.commutator(b), values, an.commutator(bn))
    assert_close(a * a, values, an * an)
    assert_close(a ** 3, values, an ** 3)


def test_sums_and_maps(case):
    a, b, values = case
    an, bn = a.subs(values), b.subs(values)
    assert_close(a + b, values, an + bn)
    assert_close(a - b, values, an - bn)
    assert_close(a - a, values, an - an)
    assert len(a - a) == 0
    assert_close(a * 2.5j, values, an * 2.5j)
    assert_close(a.transpose(), values, an.transpose())
    assert a.dagger().dagger() == a
    t = a.paulistrings()[0]
    assert a.contains(t) == t.coeff
    assert QHS(a.paulistrings()) == a
    assert QHS(a.paulistrings() + b.paulistrings()) == a + b


def test_compaction_merges_and_drops():
    x = {0: "X", 600: "Z"}
    h = QHS([PSS(E("a"), x), PSS(E("b"), {1: "Y"}), PSS(E("-a"), x), PSS(E("c"), {1: "Y"})])
    assert len(h) == 1 and h.contains(PSS(E(1), {1: "Y"})) == E("b + c")
    assert len(QHS([PSS(E("a"), {3: "X"}), PSS(E("-a"), {3: "X"})])) == 0


def sparse_sum(rng, n_terms, n_qubits, coeffs):
    """A symbolic sum that the factory stores sparse (weight 2 on many qubits)."""
    terms = [(rng.choice(coeffs), {q: rng.choice("XYZ") for q in rng.sample(range(n_qubits), 2)})
             for _ in range(n_terms)]
    h = pe.PauliStringSum(terms)
    assert type(h).__name__ == "SparsePauliStringSumSymbolic"
    return h


def cycle_sum(rng, n_cycles, n_qubits, coeffs):
    cycles = [pe.PauliCycle(n_qubits, {q: rng.choice("XYZ") for q in rng.sample(range(6), 3)}, rng.choice(coeffs))
              for _ in range(n_cycles)]
    return pe.PauliCycleSum(cycles)


def cycle_values(cycles):
    """[(base operator, value)] of a cycle sum with numeric coefficients, in its order."""
    out = []
    for c in cycles.data:
        coeff, ops = c.base.to_dictionary()
        out.append((tuple(sorted(ops.items())), PSS.to_complex(coeff) if isinstance(coeff, E) else complex(coeff)))
    return out


def test_sparse_and_cycle_products():
    # large enough for the parallel paths of sparse and cycle products
    rng = random.Random(11)
    a, b = sparse_sum(rng, 120, 1500, SIMPLE), sparse_sum(rng, 100, 1500, SIMPLE)
    an, bn = a.subs(VALUES), b.subs(VALUES)
    assert_close(a * b, VALUES, an * bn)
    assert_close(a.commutator(b), VALUES, an.commutator(bn))
    x, y = cycle_sum(rng, 40, 24, SIMPLE), cycle_sum(rng, 40, 24, SIMPLE)
    xn, yn = x.subs(VALUES), y.subs(VALUES)
    for got, want in ((x * y, xn * yn), (x.commutator(y), xn.commutator(yn))):
        # symbolic sums cancel exactly where numbers leave rounding residues
        got_values = dict(cycle_values(got.subs(VALUES)))
        want_values = dict(cycle_values(want))
        for k in got_values.keys() | want_values.keys():
            assert abs(got_values.get(k, 0) - want_values.get(k, 0)) <= 1e-9, k


SCRIPT = r"""
import random, sys
sys.path.insert(0, {tests!r})
from test_symbolic_flat import random_sum, sparse_sum, cycle_sum, SIMPLE, HEAVY
rng = random.Random(7)
for n, w, c in ((12, 3, SIMPLE), (12, 3, HEAVY), (700, 5, SIMPLE)):
    a, b = random_sum(rng, 150, n, w, c), random_sum(rng, 120, n, w, c)
    print(str(a * b)); print(str(a.commutator(b))); print(str(a * a)); print(str(a + b))
a, b = sparse_sum(rng, 120, 1500, SIMPLE), sparse_sum(rng, 100, 1500, SIMPLE)
print(str(a * b)); print(str(a.commutator(b)))
x, y = cycle_sum(rng, 40, 24, SIMPLE), cycle_sum(rng, 40, 24, SIMPLE)
for z in (x * y, x.commutator(y)):
    print("|".join(str(c.base) for c in z.data))
"""


def test_results_do_not_depend_on_the_thread_count():
    tests = os.path.dirname(os.path.abspath(__file__))
    outputs = []
    for threads in ("1", "4"):
        env = dict(os.environ, OMP_NUM_THREADS=threads)
        outputs.append(subprocess.run([sys.executable, "-c", SCRIPT.format(tests=tests)], env=env,
                                      capture_output=True, text=True, check=True).stdout)
    assert outputs[0] == outputs[1]
