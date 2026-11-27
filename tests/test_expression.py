"""The polynomial Expression used as symbolic coefficient."""

import pytest

import pauliengine as pe
from pauliengine import Expression


def test_parse_and_canonical_form():
    assert Expression("a - a").is_zero
    assert Expression("0").is_zero
    assert Expression("(a + b)**2") == Expression("a**2 + 2*a*b + b**2")
    assert Expression("b*a") == Expression("a*b")
    assert Expression("2a") == Expression("2*a")
    assert Expression("a/4") == Expression("0.25*a")
    assert Expression("a^3") == Expression("a**3")
    assert Expression("-a**2") == Expression("-(a**2)")
    assert Expression("2*I") == Expression(2j)


@pytest.mark.parametrize("text", [
    "a", "-a", "2*a*b + 3", "a**2 - b/4 + 7", "(1 + 2*I)*c", "-I*d + 0.1", "1e-300*a", "3.141592653589793*x_1",
])
def test_str_round_trip(text):
    e = Expression(text)
    assert Expression(str(e)) == e


@pytest.mark.parametrize("text", ["a/b", "a**b", "a**-1", "(a", "a +", "2 $ a", ""])
def test_rejects_non_polynomials_and_bad_input(text):
    with pytest.raises(ValueError):
        Expression(text)


def test_numbers():
    assert complex(Expression("2 + 3*I")) == 2 + 3j
    assert complex(Expression("2**10")) == 1024
    assert Expression("pi").to_complex() == pytest.approx(3.141592653589793)
    with pytest.raises(RuntimeError):
        complex(Expression("a"))


def test_arithmetic_in_python():
    a, b = Expression("a"), Expression("b")
    assert a * b + 2 == Expression("a*b + 2")
    assert (a - b) ** 2 == Expression("a**2 - 2*a*b + b**2")
    assert 3 * a - a == 2 * a
    assert -a == Expression("-a")


def test_diff_and_subs():
    e = Expression("3*a**2*b + b + 5")
    assert e.diff("a") == Expression("6*a*b")
    assert e.diff("b") == Expression("3*a**2 + 1")
    assert e.diff("unknown_symbol_xyz").is_zero
    assert complex(e.subs({"a": 2, "b": 1j})) == 3 * 4 * 1j + 1j + 5
    assert e.subs({"a": 2}) == Expression("13*b + 5")


class TestConjugationFlag:
    """conjugate() flags a symbol as conjugate(a); the flag decides what subs and diff do."""

    def test_undeclared_symbols_are_flagged(self):
        e = Expression("(1 + 2*I)*cf_a")
        c = e.conjugate()
        assert c == Expression("(1 - 2*I)*conjugate(cf_a)") and str(c) == "(1 - 2*I)*conjugate(cf_a)"
        assert c.conjugate() == e
        assert not pe.is_symbol_complex("cf_a") and not pe.is_symbol_real("cf_a")
        # whatever is put in later, the conjugated value is used
        assert c.subs({"cf_a": 1 + 1j}) == Expression((1 - 2j) * (1 - 1j))
        assert c.subs({"cf_a": 3}) == Expression((1 - 2j) * 3)
        # diff: conjugate(a) is differentiated like a (right for real parameters)
        assert Expression("cf_a**2 * conjugate(cf_a)").diff("cf_a") == Expression("2*cf_a*conjugate(cf_a) + cf_a**2")

    def test_dagger_of_a_sum(self):
        h = pe.PauliStringSum([("cf_t", {0: "X"}), ("2*I*cf_t", {0: "Y"})])
        d = h.dagger()
        assert d == pe.PauliStringSum([("conjugate(cf_t)", {0: "X"}), ("-2*I*conjugate(cf_t)", {0: "Y"})])
        for v in (0.3, 0.3 + 0.7j):
            assert d.subs({"cf_t": v}) == h.subs({"cf_t": v}).dagger()
        assert h.dagger().diff("cf_t") == h.diff("cf_t").dagger()

    def test_real_symbols_stay(self):
        pe.set_symbol_real("cf_r")
        try:
            assert pe.is_symbol_real("cf_r")
            h = pe.PauliStringSum([("cf_r", {0: "X"}), ("cf_r**2", {1: "Z"})])
            assert h.dagger() == h                      # Hermitian before substituting
            assert Expression("(1 + I)*cf_r").conjugate() == Expression("(1 - I)*cf_r")
        finally:
            pe.set_symbol_real("cf_r", False)
        assert not pe.is_symbol_real("cf_r")
        assert "conjugate(cf_r)" in str(Expression("cf_r").conjugate())

    def test_declared_real_after_conjugating(self):
        c = Expression("cf_late").conjugate()               # conjugate(cf_late) exists
        pe.set_symbol_real("cf_late")
        try:
            assert c.conjugate() == Expression("cf_late")    # goes back to the symbol
            assert c.subs({"cf_late": 2}) == Expression(2)
        finally:
            pe.set_symbol_real("cf_late", False)

    def test_complex_symbols_use_wirtinger_derivatives(self):
        pe.set_symbol_complex("cf_z")
        try:
            e = Expression("cf_z * conjugate(cf_z)")
            assert e.diff("cf_z") == Expression("conjugate(cf_z)")
            assert e.diff("conjugate(cf_z)") == Expression("cf_z")
        finally:
            pe.set_symbol_real("cf_z", False)
        e = Expression("cf_z * conjugate(cf_z)")                 # undeclared again
        assert e.diff("cf_z") == Expression("conjugate(cf_z) + cf_z")


def test_complex_symbol():
    pe.set_symbol_complex("expr_z")
    try:
        assert pe.is_symbol_complex("expr_z")
        z = Expression("expr_z")
        zc = z.conjugate()
        assert str(zc) == "conjugate(expr_z)"
        assert zc.conjugate() == z
        assert Expression(str(zc)) == zc
        # |z|^2 at z = 1 + 2i
        assert complex((z * zc).subs({"expr_z": 1 + 2j})) == pytest.approx(5)
        # Wirtinger: d/dz (z * conj z) = conj z
        assert (z * zc).diff("expr_z") == zc
        # dagger of a Hamiltonian conjugates the symbol
        h = pe.QubitHamiltonian([pe.PauliString("expr_z", {0: "X"})])
        hd = h.dagger()
        assert "conjugate(expr_z)" in str(hd)
        assert complex(hd.subs({"expr_z": 1 + 2j}).to_dictionary()[0][0]) == pytest.approx(1 - 2j)
    finally:
        pe.set_symbol_real("expr_z", False)
    # undeclared again: still conjugated, but diff counts conjugate(z) like z
    assert Expression("expr_z").conjugate() == Expression("conjugate(expr_z)")
    assert Expression("expr_z*conjugate(expr_z)").diff("expr_z") == Expression("conjugate(expr_z) + expr_z")


def test_free_symbols():
    assert Expression("b*a + c**2 + 1").free_symbols() == ["a", "b", "c"]
    assert Expression("3").free_symbols() == []


def test_symbolic_hamiltonian_matches_numeric():
    values = {"a": 0.3, "b": -1.2}
    h = pe.QubitHamiltonian([pe.PauliString("a", {0: "X"}), pe.PauliString("b*a", {1: "Y"}), pe.PauliString("2", {0: "Z"})])
    g = pe.QubitHamiltonian([pe.PauliString("b", {0: "Y"}), pe.PauliString("a**2", {1: "Z", 0: "X"})])
    hn = pe.QubitHamiltonian([pe.PauliString(0.3, {0: "X"}), pe.PauliString(-0.36, {1: "Y"}), pe.PauliString(2.0, {0: "Z"})])
    gn = pe.QubitHamiltonian([pe.PauliString(-1.2, {0: "Y"}), pe.PauliString(0.09, {1: "Z", 0: "X"})])
    for sym, num in [(h * g, hn * gn), (h.commutator(g), hn.commutator(gn)), (h + g, hn + gn), (h ** 3, hn ** 3)]:
        got = {tuple(sorted(o.items())): complex(c) for c, o in sym.subs(values).to_dictionary()}
        want = {tuple(sorted(o.items())): complex(c) for c, o in num.to_dictionary() if abs(c) > 1e-14}
        got = {k: v for k, v in got.items() if abs(v) > 1e-14}
        assert got.keys() == want.keys()
        for k in want:
            assert got[k] == pytest.approx(want[k])


class TestSubstitutePolynomials:
    """subs() also takes polynomials, as Expression or as text."""

    def test_expression(self):
        e = Expression("a**2 + a*b")
        assert e.subs({"a": "b + 1"}) == Expression("(b + 1)**2 + (b + 1)*b")
        assert e.subs({"a": Expression("2*c")}) == Expression("4*c**2 + 2*c*b")
        assert e.subs({"a": "b", "b": "a"}) == Expression("b**2 + a*b")     # all at once
        assert e.subs({"a": "x", "b": 3}) == Expression("x**2 + 3*x")          # mixed values
        assert e.subs({"a": 2, "b": 3}) == Expression("10")                     # numbers only
        assert e.subs({"a": "b - b"}).is_zero
        assert e.subs({"unknown": "x"}) == e

    def test_complex_symbols_get_the_conjugated_polynomial(self):
        pe.set_symbol_complex("zp")
        try:
            e = Expression("zp").conjugate()
            assert e.subs({"zp": "w + 2*I"}) == Expression("conjugate(w) - 2*I")   # w undeclared: flagged too
            pe.set_symbol_real("wr")
            assert e.subs({"zp": "wr + 2*I"}) == Expression("wr - 2*I")
            assert Expression("zp*conjugate(zp)").subs({"zp": "1 + I"}) == Expression("2")
        finally:
            pe.set_symbol_real("zp", False)
            pe.set_symbol_real("wr", False)

    def test_invalid_text(self):
        with pytest.raises(ValueError):
            Expression("a").subs({"a": "1/b"})

    def test_every_operator(self):
        a = Expression("a")
        h = pe.PauliStringSum([("a**2", {0: "X"}), ("a*b", {1: "Z"}), (1.0, {2: "Y"})])
        expected = pe.PauliStringSum([("(c + 1)**2", {0: "X"}), ("(c + 1)*b", {1: "Z"}), (1.0, {2: "Y"})])
        assert h.subs({"a": "c + 1"}) == expected
        assert h.subs({"a": "c + 1"}).subs({"b": 2, "c": 1}) == h.subs({"a": 2, "b": 2})
        assert pe.PauliString("a", {0: "X"}).subs({"a": "2*c"}) == pe.PauliString("2*c", {0: "X"})
        cycle = pe.PauliCycle(4, "XZ", "a")
        assert cycle.subs({"a": "c"}).to_pauli_string_sum() == pe.PauliCycle(4, "XZ", "c").to_pauli_string_sum()
        cycles = pe.PauliCycleSum([cycle, pe.PauliCycle(4, "ZZ", "b")])
        assert cycles.subs({"a": "b - b"}) == pe.PauliCycleSum([pe.PauliCycle(4, "ZZ", "b")])
        orbit = pe.PauliOrbit(4, 1, 0, 0)
        orbits = pe.PauliOrbitSum(4, [(orbit, "a")])
        assert orbits.subs({"a": "x*y"}) == pe.PauliOrbitSum(4, [(orbit, "x*y")])
        assert len(orbits.subs({"a": "x - x"}).data) == 0
        sparse = pe.PauliStringSum([("a", {j: "Z", j + 1: "Z"}) for j in range(2000)])
        assert sparse.subs({"a": "c"}) == pe.PauliStringSum([("c", {j: "Z", j + 1: "Z"}) for j in range(2000)])
        assert sparse.paulistrings()[0].subs({"a": "c"}) == pe.PauliString("c", {0: "Z", 1: "Z"})
        assert (a * 2).subs({"a": "c"}) == Expression("2*c")
