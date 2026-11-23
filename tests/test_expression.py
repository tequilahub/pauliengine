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


def test_symbols_are_real_by_default():
    e = Expression("(1 + 2*I)*a")
    assert e.conjugate() == Expression("(1 - 2*I)*a")
    assert not pe.is_symbol_complex("a")


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
        pe.set_symbol_complex("expr_z", False)
    assert Expression("expr_z").conjugate() == Expression("expr_z")


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
