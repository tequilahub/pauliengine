"""Sparse Pauli strings and Hamiltonians (SparsePauliString.h, SparsePauliStringSum.h).

Every operation must give the same operator as the classes with bit vectors.
"""

import random

import pytest

import pauliengine as pe
import pauliengine._core as _core

@pytest.fixture(autouse=True)
def always_sparse():
    """Most tests use few qubits, where operators would not be sparse (see
    TestDenseResults); keep everything sparse unless a test sets its own threshold."""
    saved = pe.get_sparse_threshold()
    pe.set_sparse_threshold(min_qubits=0, max_density=1.0)
    yield
    pe.set_sparse_threshold(**saved)


SPC = _core.SparsePauliStringComplex
SPS = _core.SparsePauliStringSymbolic
SQC = _core.SparsePauliStringSumComplex
SQS = _core.SparsePauliStringSumSymbolic


def random_terms(rng, n_terms, n_qubits, max_weight, complex_coeffs=True):
    terms = []
    for _ in range(n_terms):
        w = rng.randint(0, max_weight)
        qubits = rng.sample(range(n_qubits), w)
        ops = {q: rng.choice("XYZ") for q in qubits}
        c = complex(rng.uniform(-1, 1), rng.uniform(-1, 1) if complex_coeffs else 0.0)
        terms.append((c, ops))
    return terms


def dense_qh(terms):
    return _core.PauliStringSumComplex([_core.PauliStringComplex(c, d) for c, d in terms])


def as_dict(h):
    """{sorted operator items: complex coefficient} of any complex Hamiltonian."""
    out = {}
    for c, d in h.to_dictionary():
        out[tuple(sorted(d.items()))] = complex(c)
    return out


def assert_same(sparse, dense):
    """Equal operators; terms that cancel up to rounding may be kept by one side only."""
    a, b = as_dict(sparse), as_dict(dense)
    scale = max([abs(c) for c in b.values()] + [1.0])
    for k in a.keys() | b.keys():
        assert a.get(k, 0) == pytest.approx(b.get(k, 0), abs=1e-12 * scale), k


CASES = [
    # (terms, qubits, max weight)
    (30, 8, 4),
    (60, 300, 3),
    (40, 2000, 6),
    (25, 130, 130),
]


@pytest.fixture(params=range(len(CASES)), ids=lambda i: f"case{i}")
def pair(request):
    rng = random.Random(1234 + request.param)
    n, q, w = CASES[request.param]
    return random_terms(rng, n, q, w), random_terms(rng, n + 7, q, w)


class TestStrings:
    def test_product_and_commutator_match_dense(self):
        rng = random.Random(1)
        for _ in range(500):
            (ca, da), (cb, db) = random_terms(rng, 2, 200, 8)
            a, b = SPC(ca, da), SPC(cb, db)
            pa, pb = a.to_dense(), b.to_dense()
            prod = a * b
            assert isinstance(prod, _core.PauliStringComplex) and prod == pa * pb
            c = a.commutator(b)
            dc = pa.commutator(pb)
            assert isinstance(c, _core.PauliStringComplex)
            if dc.x:
                assert c == dc
            else:
                assert complex(c.coeff) == 0
            assert a.anticommutes(b) == (len(dc.x) > 0)

    def test_round_trip_and_queries(self):
        s = SPC(0.5j, {700: "Y", 3: "X", 64: "Z"})
        assert s.sites == [(3, "X"), (64, "Z"), (700, "Y")]
        assert SPC(s.to_dense()) == s
        assert s.qubits() == [3, 64, 700]
        assert s.get_pauli_at_index(64) == "Z" and s.get_pauli_at_index(65) == "I"
        assert len(s) == 3 and s.count_y() == 1 and not s.is_all_z()
        assert hash(s) == hash(s.to_dense()) and s == s.to_dense()
        assert SPC(1.0, "XIZ") == SPC(1.0, {0: "X", 2: "Z"})
        assert SPC(1.0, "X0 Z2") == SPC(1.0, {0: "X", 2: "Z"})

    def test_targeted_commutator(self):
        a = SPC(1.0, {0: "X", 3: "Y"})
        b = SPC(2.0, {3: "Z"})
        full = a.commutator(b)
        assert a.targeted_commutator(b, SPC(full.naked())) == pytest.approx(complex(full.coeff))
        assert a.targeted_commutator(b, SPC(1.0, {0: "X"})) == 0


class TestHamiltonians:
    def test_construction_compacts(self):
        h = SQC([SPC(1.0, {5: "X"}), SPC(2.0, {5: "X"}), SPC(1.0, {1: "Z"}), SPC(-1.0, {1: "Z"})])
        assert len(h) == 1 and h.contains(SPC(1.0, {5: "X"})) == 3
        assert_same(h, dense_qh([(3.0, {5: "X"})]))

    def test_conversions(self, pair):
        ta, _ = pair
        d = dense_qh(ta)
        s = SQC(d)
        assert_same(s, d)
        assert s.to_dense() == d
        assert s == d and d == s

    def test_sum_and_difference(self, pair):
        ta, tb = pair
        a, b = SQC(dense_qh(ta)), SQC(dense_qh(tb))
        assert_same(a + b, dense_qh(ta) + dense_qh(tb))
        assert_same(a - b, dense_qh(ta) - dense_qh(tb))
        assert len(a - a) == 0

    def test_product(self, pair):
        ta, tb = pair
        a, b = SQC(dense_qh(ta)), SQC(dense_qh(tb))
        assert_same(a * b, dense_qh(ta) * dense_qh(tb))

    def test_square(self, pair):
        ta, _ = pair
        a = SQC(dense_qh(ta))
        da = dense_qh(ta)
        assert_same(a * a, da * da)
        assert_same(a ** 3, da ** 3)

    def test_commutator(self, pair):
        ta, tb = pair
        a, b = SQC(dense_qh(ta)), SQC(dense_qh(tb))
        assert_same(a.commutator(b), dense_qh(ta).commutator(dense_qh(tb)))

    def test_local_commutator_uses_overlap(self):
        n = 3000
        chain = [(1.0, {j: p, j + 1: p}) for j in range(n - 1) for p in "XYZ"]
        field = [(0.3, {j: "X"}) for j in range(n)]
        a, b = SQC(dense_qh(chain)), SQC(dense_qh(field))
        assert_same(a.commutator(b), dense_qh(chain).commutator(dense_qh(field)))

    def test_targeted_commutator(self, pair):
        ta, tb = pair
        a, b = SQC(dense_qh(ta)), SQC(dense_qh(tb))
        full = as_dict(a.commutator(b))
        for key, c in list(full.items())[:20]:
            target = SPC(1.0, dict(key))
            assert a.targeted_commutator(b, target) == pytest.approx(c)
            assert b.targeted_commutator(a, target) == pytest.approx(-c)

    def test_inner_product_and_contains(self, pair):
        ta, tb = pair
        a, b = SQC(dense_qh(ta)), SQC(dense_qh(tb + ta[:5]))
        assert a.inner_product(b) == pytest.approx(dense_qh(ta).inner_product(dense_qh(tb + ta[:5])))
        for c, d in ta[:5]:
            assert a.contains(SPC(1.0, d)) == pytest.approx(dense_qh(ta).contains(_core.PauliStringComplex(1.0, d)))

    def test_coefficient_maps(self, pair):
        ta, _ = pair
        a, d = SQC(dense_qh(ta)), dense_qh(ta)
        assert_same(a.dagger(), d.dagger())
        assert_same(a.conjugate(), d.conjugate())
        assert_same(a.transpose(), d.transpose())
        assert_same(a.simplify(0.5), d.simplify(0.5))
        assert_same(a * 2j, d * 2j)
        assert_same(-a, -d)
        assert_same(a / 4, d / 4)
        assert a.qubits() == d.qubits() and a.n_qubits() == d.n_qubits()
        assert a.is_hermitian() == d.is_hermitian()

    def test_map_qubits(self):
        a = SQC([SPC(1.0, {0: "X", 900: "Z"}), SPC(0.5, {3: "Y"})])
        mapped = a.map_qubits({900: 1, 3: 2})
        assert_same(mapped, a.to_dense().map_qubits({900: 1, 3: 2}))

    def test_symbolic(self):
        a = SQS([SPS(_core.Expression("a"), {0: "X", 400: "Y"}), SPS(_core.Expression("b"), {400: "Z"})])
        sq = a * a
        assert len(sq.diff("a")) > 0
        values = {"a": 0.7, "b": -0.2}
        dense = a.to_dense()
        expected = (dense * dense).subs(values)
        got = sq.subs(values)
        assert got == expected
        assert a.commutator(a.diff("a")).subs(values) == dense.commutator(dense.diff("a")).subs(values)


class TestMixed:
    def test_sums_stay_sparse_products_do_not(self):
        s = SQC([SPC(1.0, {0: "X", 500: "Z"})])
        d = _core.PauliStringSumComplex([_core.PauliStringComplex(1.0, {500: "X"})])
        p = _core.PauliStringComplex(2.0, {0: "Y"})
        for r in (s + d, d + s, s + p, p - s):
            assert isinstance(r, SQC)
        for r in (s * d, d * s, s.commutator(d), d.commutator(s), s * p, p * s, s * s, s.commutator(s), s ** 2):
            assert isinstance(r, _core.PauliStringSumComplex)
        assert isinstance(SPC(1.0, {0: "X"}) * p, _core.PauliStringComplex)
        assert isinstance(p.commutator(SPC(1.0, {0: "X"})), _core.PauliStringComplex)
        assert isinstance(SPC(1.0, {0: "X"}) * SPC(1.0, {3: "Z"}), _core.PauliStringComplex)
        assert_same(s * d, s.to_dense() * d)

    def test_symbolic_promotion(self):
        s = SQC([SPC(1.0, {0: "X"})])
        t = SQS([SPS(_core.Expression("a"), {0: "Z"})])
        assert isinstance(s * t, _core.PauliStringSumSymbolic) and isinstance(t + s, SQS)
        assert isinstance(s * _core.Expression("b"), SQS)       # scaling keeps the sparse form

    def test_with_cycles_orbits_and_unweighted(self):
        s = SQC([SPC(1.0, {0: "X", 2: "Z"})])
        cycle = pe.PauliCycle(4, "XIII", 0.5)
        r = s * cycle
        assert isinstance(r, _core.PauliStringSumComplex)
        assert_same(r, s.to_dense() * cycle.to_pauli_string_sum())
        assert isinstance(s + cycle, SQC)
        u = pe.PauliString("XIZ")
        assert isinstance(s.commutator(u), _core.PauliStringSumComplex)
        assert s.inner_product(pe.PauliOrbit(4, 1, 0, 1)) == pytest.approx(
            s.to_dense().inner_product(pe.PauliOrbit(4, 1, 0, 1)))


class TestConstructors:
    """Every constructor that takes a Pauli string also takes a sparse one."""

    def test_sums(self):
        a, b = SPC(0.5, {0: "X", 700: "Z"}), _core.PauliStringComplex(2.0, {3: "Y"})
        expected = _core.PauliStringSumComplex([a.to_dense(), b])
        for got in (_core.PauliStringSumComplex([a, b]), _core.PauliStringSumComplex([b, a]),
                    _core.PauliStringSumComplex([a]) + _core.PauliStringSumComplex([b]),
                    SQC([a, b]), SQC([b, a]), pe.PauliStringSum([a, b])):
            assert_same(got, expected)
        assert isinstance(_core.PauliStringSumComplex([a, b]), _core.PauliStringSumComplex)
        assert isinstance(SQC([b, a]), SQC)
        x = _core.Expression("x")
        sym = _core.PauliStringSumSymbolic([SPS(x, {1: "Z"}), _core.PauliStringSymbolic(x, {1: "Z"})])
        assert len(sym) == 1 and sym.subs({"x": 1.5}) == _core.PauliStringSumSymbolic([_core.PauliStringSymbolic(_core.Expression("3"), {1: "Z"})]).subs({})

    def test_cycles(self):
        base = SPC(0.5, {0: "X", 1: "Z"})
        assert pe.PauliCycle(4, base) == pe.PauliCycle(4, base.to_dense())
        assert complex(pe.PauliCycle(4, base).base.coeff) == 0.5
        sym = pe.PauliCycle(4, SPS(_core.Expression("a"), {0: "Y"}))
        assert isinstance(sym, _core.PauliCycleSymbolic)
        assert sym.to_pauli_string_sum() == pe.PauliCycle(4, _core.PauliStringSymbolic(_core.Expression("a"), {0: "Y"})).to_pauli_string_sum()

    def test_unweighted(self):
        s = SPC(3.0, {0: "X", 900: "Y"})
        u = _core.unweighted_pauliString(s, 1)
        assert u == _core.unweighted_pauliString(s.to_dense(), 1)
        assert _core.unweighted_pauliString(SPS(_core.Expression("a"), {2: "Z"})) == pe.PauliString("IIZ")
        assert pe.PauliString(s) == _core.unweighted_pauliString(s.to_dense())
        h = SQC([s, SPC(1.0, {5: "Z"})])
        assert _core.unweighted_PauliStringSum(h) == _core.unweighted_PauliStringSum(h.to_dense())
        hs = SQS([SPS(_core.Expression("a"), {5: "Z"})])
        assert len(_core.unweighted_PauliStringSum(hs)) == 1


class TestSameInterface:
    """Code must not care whether the factories chose bit vectors or sparse sites."""

    def test_isinstance_through_the_factories(self):
        x = _core.Expression("x")
        for s in (_core.PauliStringComplex(1.0, {0: "X"}), _core.PauliStringSymbolic(x, {0: "X"}),
                  SPC(1.0, {0: "X"}), SPS(x, {0: "X"})):
            assert isinstance(s, pe.PauliString) and not isinstance(s, pe.PauliStringSum)
        for h in (_core.PauliStringSumComplex([_core.PauliStringComplex(1.0, {0: "X"})]),
                  _core.PauliStringSumSymbolic([_core.PauliStringSymbolic(x, {0: "X"})]),
                  SQC([SPC(1.0, {0: "X"})]), SQS([SPS(x, {0: "X"})])):
            assert isinstance(h, pe.PauliStringSum) and isinstance(h, pe.QubitHamiltonian)
            assert not isinstance(h, pe.PauliString)
        assert not isinstance(pe.PauliString("XZ"), pe.PauliString)              # unweighted
        assert not isinstance(pe.PauliStringSum([pe.PauliString("XZ")]), pe.PauliStringSum)
        assert not isinstance(3, pe.PauliStringSum)

    def test_same_methods(self):
        def methods(t):
            return {n for n in dir(t) if not n.startswith("_")}

        sparse_only = {"to_dense", "fits_sparse", "sites"}
        for dense, sparse in ((_core.PauliStringComplex, SPC), (_core.PauliStringSymbolic, SPS),
                              (_core.PauliStringSumComplex, SQC), (_core.PauliStringSumSymbolic, SQS)):
            assert methods(dense) <= methods(sparse), (dense, methods(dense) - methods(sparse))
            assert methods(sparse) - methods(dense) <= sparse_only, (sparse, methods(sparse) - methods(dense))

    def test_square_anticommutes_from_hash_from_matrix(self):
        d = _core.PauliStringSumComplex([_core.PauliStringComplex(1.0, {0: "X"}), _core.PauliStringComplex(0.5, {0: "Z", 1: "Y"})])
        assert d.square() == d * d
        assert_same(SQC(d).square(), d.square())
        x = _core.Expression("x")
        ds = _core.PauliStringSumSymbolic([_core.PauliStringSymbolic(x, {0: "X"}), _core.PauliStringSymbolic(x, {1: "Z"})])
        assert ds.square().subs({"x": 2.0}) == (ds * ds).subs({"x": 2.0})
        a, b = _core.PauliStringComplex(1.0, {0: "X", 700: "Z"}), _core.PauliStringComplex(1.0, {0: "Z"})
        assert a.anticommutes(b) and not a.anticommutes(a)
        assert a.anticommutes(b) == SPC(a).anticommutes(SPC(b))
        sa = _core.PauliStringSymbolic(x, {0: "Y"})
        assert sa.anticommutes(_core.PauliStringSymbolic(x, {0: "Z"})) and not sa.anticommutes(sa)
        h = a.get_hash()
        assert SPC.from_hash(h, 0.5) == _core.PauliStringComplex.from_hash(h, 0.5)
        assert SPS.from_hash(h).to_dense() == _core.PauliStringSymbolic.from_hash(h)
        assert SPS.to_complex(_core.Expression("2")) == 2
        m = [[0, 1], [2, 3]]
        assert_same(SQC.from_matrix(m), _core.PauliStringSumComplex.from_matrix(m))
        import numpy as np
        assert isinstance(SQC.from_matrix(np.array(m, dtype=complex)), SQC)


class TestFactory:
    @pytest.fixture(autouse=True)
    def default_threshold(self):
        pe.set_sparse_threshold(min_qubits=256, max_density=1 / 256)
        yield

    def test_strings(self):
        assert isinstance(pe.PauliString(1.0, {2000: "X", 5: "Z"}), SPC)
        assert isinstance(pe.PauliString(1.0, "X5 Z2000"), SPC)
        assert isinstance(pe.PauliString("a", {3000: "Y"}), SPS)
        assert isinstance(pe.PauliString(1.0, {20: "X"}), _core.PauliStringComplex)
        assert isinstance(pe.PauliString(1.0, {600: "X", 1: "Z"}), SPC)          # weight 2 on 601 qubits
        assert isinstance(pe.PauliString(1.0, {300: "X", 1: "Z"}), _core.PauliStringComplex)   # 301 qubits
        assert isinstance(pe.PauliString(1.0, {1023: "X", 1: "Z", 2: "Y", 3: "X"}), SPC)     # weight 4 on 1024
        assert isinstance(pe.PauliString(1.0, {1000: "X", 1: "Z", 2: "Y", 3: "X"}), _core.PauliStringComplex)
        assert pe.PauliString(1.0, {2000: "X"}) == _core.PauliStringComplex(1.0, {2000: "X"})

    def test_hamiltonians(self):
        chain = [(1.0, {j: "Z", j + 1: "Z"}) for j in range(2000)]
        assert isinstance(pe.PauliStringSum(chain), SQC)
        assert isinstance(pe.PauliStringSum(chain[:100]), _core.PauliStringSumComplex)
        assert isinstance(pe.PauliStringSum([("a", d) for _, d in chain]), SQS)
        dense_terms = [_core.PauliStringComplex(c, d) for c, d in chain]
        assert isinstance(pe.PauliStringSum(dense_terms), SQC)
        mixed = [SPC(1.0, {0: "X"}), _core.PauliStringComplex(2.0, {1: "Z"})]
        assert isinstance(pe.PauliStringSum(mixed), _core.PauliStringSumComplex)   # 2 qubits: not sparse
        assert isinstance(pe.PauliStringSum([SPC(1.0, {0: "X"})]), _core.PauliStringSumComplex)
        assert_same(pe.PauliStringSum(chain), dense_qh(chain))

    def test_every_term_must_be_light(self):
        chain = [(1.0, {j: "Z", j + 1: "Z"}) for j in range(2000)]
        heavy = (0.1, {q: "X" for q in range(0, 2000, 200)})              # weight 10 > 2001 / 256
        assert isinstance(pe.PauliStringSum(chain + [heavy]), _core.PauliStringSumComplex)
        light = (0.1, {q: "X" for q in range(0, 2000, 500)})              # weight 4 <= 2001 / 256
        assert isinstance(pe.PauliStringSum(chain + [light]), SQC)

    def test_input_form_does_not_matter(self):
        chain = [(1.0, {j: "Z", j + 1: "Z"}) for j in range(2000)]
        dense_strings = [_core.PauliStringComplex(c, d) for c, d in chain]
        sparse_strings = [SPC(c, d) for c, d in chain]
        mixed = dense_strings[:1000] + sparse_strings[1000:]
        for terms in (chain, dense_strings, sparse_strings, mixed):
            assert isinstance(pe.PauliStringSum(terms), SQC)
        small = [SPC(1.0, {0: "X"}), SPC(1.0, {1: "Z"})]
        assert isinstance(pe.PauliStringSum(small), _core.PauliStringSumComplex)

    def test_existing_sums_follow_the_rule(self):
        chain = [(1.0, {j: "Z", j + 1: "Z"}) for j in range(2000)]
        assert isinstance(pe.PauliStringSum(dense_qh(chain)), SQC)
        small = SQC([SPC(1.0, {0: "X"})])
        assert isinstance(pe.PauliStringSum(small), _core.PauliStringSumComplex)
        assert isinstance(pe.PauliStringSum(dense_qh(chain[:10])), _core.PauliStringSumComplex)
        assert_same(pe.PauliStringSum(dense_qh(chain)), dense_qh(chain))

    def test_openfermion_follows_the_rule(self):
        class QubitOperator:          # duck type of openfermion.QubitOperator
            def __init__(self, terms):
                self.terms = terms

        big = QubitOperator({((j, "Z"), (j + 1, "Z")): 1.0 for j in range(2000)})
        assert isinstance(pe.PauliStringSum(big), SQC)
        assert isinstance(pe.PauliStringSum(QubitOperator({((0, "X"),): 1.0})), _core.PauliStringSumComplex)

    def test_threshold(self):
        pe.set_sparse_threshold(min_qubits=0, max_density=1.0)
        assert isinstance(pe.PauliString(1.0, {2: "X"}), SPC)
        pe.set_sparse_threshold(min_qubits=float("inf"))
        assert isinstance(pe.PauliString(1.0, {5000: "X"}), _core.PauliStringComplex)
        pe.set_sparse_threshold(min_qubits=1024, max_density=1 / 32)
        heavy = {q: "X" for q in range(0, 2000, 10)}     # weight 200 on 2000 qubits
        assert isinstance(pe.PauliString(1.0, heavy), _core.PauliStringComplex)
        assert pe.get_sparse_threshold() == {"min_qubits": 1024, "max_density": 1 / 32}

    def test_default_threshold(self):
        assert pe.get_sparse_threshold() == {"min_qubits": 256, "max_density": 1 / 256}
        with pytest.raises(ValueError):
            pe.set_sparse_threshold(max_density=2)


def _others():
    """One operator of every kind, complex and symbolic."""
    a = _core.Expression("a")
    osc = _core.PauliOrbitSumComplex(4)
    osc.add(pe.PauliOrbit(4, 0, 1, 0), 0.5)
    oss = _core.PauliOrbitSumSymbolic(4)
    oss.add(pe.PauliOrbit(4, 1, 0, 1), a)
    return {
        "PSC": _core.PauliStringComplex(0.5, {0: "X", 1: "Z"}),
        "PSS": _core.PauliStringSymbolic(a, {1: "Y"}),
        "QHC": _core.PauliStringSumComplex([_core.PauliStringComplex(1.5, {0: "Y"}), _core.PauliStringComplex(0.5, {2: "Z"})]),
        "QHS": _core.PauliStringSumSymbolic([_core.PauliStringSymbolic(a, {3: "X"})]),
        "CYC": pe.PauliCycle(4, "YZII", 0.5),
        "CYS": _core.PauliCycleSymbolic(4, "XIZI", a),
        "CSC": pe.PauliCycleSum([pe.PauliCycle(4, "XIII", 0.5), pe.PauliCycle(4, "ZZII", 1.0)]),
        "CSS": _core.PauliCycleSumSymbolic([_core.PauliCycleSymbolic(4, "XIII", a)]),
        "OR": pe.PauliOrbit(4, 1, 0, 1),
        "OSC": osc,
        "OSS": oss,
        "UPS": pe.PauliString("XZ", phase=1),
        "UQH": _core.unweighted_PauliStringSum([pe.PauliString("XZ"), pe.PauliString("IY")]),
    }


def _sparse_operands():
    a = _core.Expression("a")
    return {
        "SPC": SPC(0.7, {0: "Z", 2: "X"}),
        "SPS": SPS(a, {1: "X", 3: "Z"}),
        "SQC": SQC([SPC(1.0, {1: "Y"}), SPC(0.25, {0: "X", 3: "X"})]),
        "SQS": SQS([SPS(a, {1: "Y"}), SPS(_core.Expression("2"), {2: "Z"})]),
    }


def _as_complex_dict(op):
    """{operator: value at a = 0.7} of any weighted operator."""
    if isinstance(op, (SPC, SPS, SQC, SQS)):
        op = op.to_dense()
    if isinstance(op, (_core.PauliStringComplex, _core.PauliStringSymbolic)):
        op = (_core.PauliStringSumSymbolic if isinstance(op, _core.PauliStringSymbolic) else _core.PauliStringSumComplex)([op])
    if isinstance(op, _core.PauliStringSumSymbolic):
        op = op.subs({"a": 0.7})
        return {tuple(sorted(d.items())): _core.PauliStringSymbolic.to_complex(c) for c, d in op.to_dictionary()}
    return {tuple(sorted(d.items())): complex(c) for c, d in op.to_dictionary()}


def _close(x, y):
    keys = x.keys() | y.keys()
    return all(abs(x.get(k, 0) - y.get(k, 0)) < 1e-12 for k in keys)


DENSE = (_core.PauliStringComplex, _core.PauliStringSymbolic, _core.PauliStringSumComplex, _core.PauliStringSumSymbolic)

OPS = {
    "+": lambda x, y: x + y,
    "-": lambda x, y: x - y,
    "*": lambda x, y: x * y,
    "commutator": lambda x, y: x.commutator(y),
}


@pytest.mark.parametrize("other", list(_others()))
@pytest.mark.parametrize("sparse", list(_sparse_operands()))
@pytest.mark.parametrize("op", list(OPS))
def test_sparse_with_every_operator(op, sparse, other):
    """A sparse operand with any other operator, in both orders: the same operator as
    the computation with bit vectors; sums stay sparse, products and commutators
    come back with bit vectors."""
    s, o = _sparse_operands()[sparse], _others()[other]
    f = OPS[op]
    for left, right, dense_left, dense_right in ((s, o, s.to_dense(), o), (o, s, o, s.to_dense())):
        got = f(left, right)
        if op in ("+", "-"):
            assert isinstance(got, (SPC, SPS, SQC, SQS)), (op, sparse, other)
        else:
            assert isinstance(got, DENSE), (op, sparse, other)
        expected = f(dense_left, dense_right)
        assert _close(_as_complex_dict(got), _as_complex_dict(expected)), (op, sparse, other)


@pytest.mark.parametrize("other", list(_others()))
@pytest.mark.parametrize("sparse", list(_sparse_operands()))
def test_sparse_inner_product_with_every_operator(sparse, other):
    s, o = _sparse_operands()[sparse], _others()[other]

    def value(v):
        if isinstance(v, _core.Expression):
            v = _core.PauliStringSymbolic(v, {}).subs({"a": 0.7}).coeff
            return _core.PauliStringSymbolic.to_complex(v)
        return complex(v)

    for got, expected in ((s.inner_product(o), s.to_dense().inner_product(o)),
                          (o.inner_product(s), o.inner_product(s.to_dense()))):
        assert value(got) == pytest.approx(value(expected), abs=1e-12)


class TestDenseResults:
    """Products and commutators of sparse operators come back with bit vectors;
    sums leave the sparse form only next to a heavy operator with bit vectors."""

    @pytest.fixture(autouse=True)
    def default_threshold(self):
        pe.set_sparse_threshold(min_qubits=256, max_density=1 / 256)
        yield

    def test_repeated_products(self):
        h = pe.PauliStringSum([(1.0, {j: "X", j + 1: "Z"}) for j in range(0, 1024, 128)])
        assert isinstance(h, SQC) and h.fits_sparse()
        p = h
        for _ in range(4):
            p = p * h
            assert isinstance(p, _core.PauliStringSumComplex)
        dense = h.to_dense()
        assert_same(p, dense * dense * dense * dense * dense)

    def test_commutator_and_power(self):
        n = 1024
        h = pe.PauliStringSum([(1.0, {j: "X", j + 1: "Y", j + 2: "Z"}) for j in range(0, n - 2, 170)])   # 6 terms
        assert isinstance(h, SQC)
        for r in (h ** 4, h.power(2), h.square(), h.commutator(h * 1.5)):
            assert isinstance(r, _core.PauliStringSumComplex)
        assert_same(h ** 4, h.to_dense() ** 4)
        s = SPC(1.0, {0: "X", 100: "Y", 300: "Z", 511: "X"})
        assert isinstance(s * SPC(1.0, {0: "Z"}), _core.PauliStringComplex)
        assert isinstance(s.commutator(SPC(1.0, {0: "Z"})), _core.PauliStringComplex)

    def test_sum_with_heavy_dense_operand(self):
        h = pe.PauliStringSum([(1.0, {j: "X", j + 1: "Z"}) for j in range(0, 1024, 128)])
        heavy = _core.PauliStringSumComplex([_core.PauliStringComplex(1.0, {q: "Y" for q in range(k, 1024, 50)})
                                               for k in range(3)])       # weight about 20 on 1024 qubits
        for r in (heavy + h, h + heavy, h - heavy):
            assert isinstance(r, _core.PauliStringSumComplex)        # no detour through sparse
        assert_same(heavy + h, heavy + h.to_dense())
        assert isinstance(h + h.to_dense(), SQC)                        # a light dense operand stays sparse

    def test_threshold_is_shared_with_cpp(self):
        pe.set_sparse_threshold(min_qubits=float("inf"))
        assert _core._sparse_threshold()[0] == float("inf")
        assert isinstance(SQC([SPC(1.0, {5000: "X"})]) + SQC([SPC(1.0, {5000: "Z"})]), _core.PauliStringSumComplex)
