import itertools
import random

import numpy as np
import pytest

import pauliengine as pe
import pauliengine._core as _core

MATRICES = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def dense(word, n_qubits):
    """The matrix i**phase * P of a word on n_qubits qubits."""
    matrix = np.eye(1, dtype=complex)
    for qubit in range(n_qubits):
        matrix = np.kron(matrix, MATRICES[word.get_pauli_at_index(qubit)])
    return word.phase_factor * matrix


def random_labels(rng, n_qubits):
    return {q: rng.choice("IXYZ") for q in range(n_qubits)}


@pytest.mark.parametrize(
    "build",
    [
        lambda: pe.PauliString("X0 Y1 Z70"),
        lambda: pe.PauliString({0: "X", 1: "Y", 2: "I", 70: "Z"}),
        lambda: pe.PauliString([("X", 0), ("Y", 1), ("Z", 70)]),
        lambda: pe.PauliString("XY" + "I" * 68 + "Z"),
        lambda: pe.PauliString(pe.PauliString(2.5, {0: "X", 1: "Y", 70: "Z"})),
        lambda: pe.PauliString(pe.PauliString("a", {0: "X", 1: "Y", 70: "Z"})),
    ],
)
def test_inputs_give_the_same_word(build):
    word = build()
    assert word == pe.PauliString("X0 Y1 Z70")
    assert word.phase == 0
    assert word.qubits() == [0, 1, 70]
    assert word.key_openfermion() == [("X", 0), ("Y", 1), ("Z", 70)]
    assert len(word) == word.size() == 3
    assert word.count_y() == 1
    assert word.get_pauli_at_index(70) == "Z"
    assert word.get_pauli_at_index(5) == "I"


def test_no_coefficient_gives_an_unweighted_string():
    for s in (pe.PauliString(), pe.PauliString("XIZ"), pe.PauliString("X0"), pe.PauliString({0: "X"}),
              pe.PauliString([("X", 0)]), pe.PauliString(pauli_strings="X0"), pe.PauliString(iter([("X", 0)])),
              pe.PauliString(pe.PauliString(2.5, "X0")), pe.PauliString("X0", phase=1)):
        assert isinstance(s, _core.unweighted_pauliString)


def test_a_coefficient_gives_a_weighted_string():
    assert isinstance(pe.PauliString(1.0, "XIZ"), _core.PauliStringComplex)
    assert isinstance(pe.PauliString(2, "XIZ"), _core.PauliStringComplex)
    assert isinstance(pe.PauliString(1j, "XIZ"), _core.PauliStringComplex)
    assert isinstance(pe.PauliString("a", "XIZ"), _core.PauliStringSymbolic)
    assert isinstance(pe.PauliString(2.0), _core.PauliStringComplex)
    assert isinstance(pe.PauliString("a", {}), _core.PauliStringSymbolic)
    assert isinstance(pe.PauliString((1.0, "X0")), _core.PauliStringComplex)


def test_only_one_public_class():
    assert not hasattr(pe, "unweighted_pauliString")
    assert not hasattr(pe, "unweighted_QubitHamiltonian")


def test_phase_needs_an_unweighted_string():
    with pytest.raises(ValueError):
        pe.PauliString(1.0, "X0", phase=1)
    with pytest.raises(ValueError):
        pe.PauliString(1.0, phase=1)


@pytest.mark.parametrize("coeff", [2.5, 1j, 3, "a*b"])
@pytest.mark.parametrize(
    ("dense", "sparse"),
    [("XIZ", "X0 Z2"), ("IIY", "Y2"), ("III", ""), ("", ""), ("I" * 64 + "X", "X64"), ("I" * 200 + "Z", "Z200")],
)
def test_dense_and_sparse_text_match_the_dict(coeff, dense, sparse):
    ops = {int(t[1:]): t[0] for t in sparse.split()}
    reference = pe.PauliString(coeff, ops)
    for text in (dense, sparse):
        built = pe.PauliString(coeff, text)
        assert built == reference
        assert built.x == reference.x
        assert built.y == reference.y


@pytest.mark.parametrize("coeff", [1.0, "a"])
@pytest.mark.parametrize("text", ["XQZ", "xyz", "X0 Y0", "Q0", "X-1", "X 0"])
def test_invalid_text_is_rejected_for_every_coefficient(coeff, text):
    with pytest.raises(ValueError):
        pe.PauliString(coeff, text)


def test_dense_text_of_many_random_strings():
    rng = random.Random(7)
    for n_qubits in (1, 63, 64, 65, 128, 129, 1000):
        for _ in range(20):
            text = "".join(rng.choice("IXYZ") for _ in range(n_qubits))
            ops = {i: p for i, p in enumerate(text) if p != "I"}
            assert pe.PauliString(text) == pe.PauliString(ops)
            assert pe.PauliString(1.0, text) == pe.PauliString(1.0, ops)


def test_scalars_and_weighted_factors_give_weighted_results():
    s = pe.PauliString("X", phase=1)
    assert s * 0.5 == pe.PauliString(0.5j, "X")
    assert 0.5 * s == pe.PauliString(0.5j, "X")
    assert 2 * pe.PauliString("X") == pe.PauliString(2.0, "X")
    assert pe.PauliString("X") * pe.PauliString(2.0, "Y") == pe.PauliString(2j, "Z")
    assert pe.PauliString(2.0, "X") * pe.PauliString("Y") == pe.PauliString(2j, "Z")
    assert pe.PauliString("X") * pe.PauliString("a", "X") == pe.PauliString("a", {})
    assert pe.PauliString("a", "X") * pe.PauliString("X") == pe.PauliString("a", {})


def test_hamiltonian_factory_recognizes_unweighted_strings():
    h = pe.QubitHamiltonian([pe.PauliString("X"), pe.PauliString("Y")])
    assert isinstance(h, _core.unweighted_QubitHamiltonian)
    assert isinstance(pe.QubitHamiltonian(h), _core.unweighted_QubitHamiltonian)
    assert pe.QubitHamiltonian(h) == h
    assert isinstance(pe.QubitHamiltonian([pe.PauliString(1.0, "X")]), _core.QubitHamiltonianComplex)


def test_unweighted_next_to_complex_terms_becomes_complex():
    h = pe.QubitHamiltonian([pe.PauliString("X"), pe.PauliString(0.5, "Y"), pe.PauliString("Z", phase=1)])
    assert isinstance(h, _core.QubitHamiltonianComplex)
    assert h == pe.QubitHamiltonian([pe.PauliString(1.0, "X"), pe.PauliString(0.5, "Y"), pe.PauliString(1j, "Z")])


def test_unweighted_next_to_symbolic_terms_becomes_symbolic():
    h = pe.QubitHamiltonian([pe.PauliString("X"), pe.PauliString("a", "Y"), pe.PauliString("Z", phase=1)])
    assert isinstance(h, _core.QubitHamiltonianSymbolic)
    assert h == pe.QubitHamiltonian([pe.PauliString(1.0, "X"), pe.PauliString("a", "Y"), pe.PauliString(1j, "Z")])


@pytest.mark.parametrize(
    ("term", "kind"),
    [((0.5, {1: "Z"}), _core.QubitHamiltonianComplex), (("a", {1: "Z"}), _core.QubitHamiltonianSymbolic)],
)
def test_unweighted_next_to_tuples(term, kind):
    assert isinstance(pe.QubitHamiltonian([pe.PauliString("X"), term]), kind)
    assert isinstance(pe.QubitHamiltonian([term, pe.PauliString("X")]), kind)


def test_unweighted_terms_merge_with_weighted_ones_as_a_sum():
    h = pe.QubitHamiltonian([pe.PauliString("X"), pe.PauliString(2.0, "X"), pe.PauliString("Y", phase=2), pe.PauliString(1.0, "Y")])
    assert h == pe.QubitHamiltonian([pe.PauliString(3.0, "X")])


@pytest.mark.parametrize(
    "paulis",
    ["Q0", "XQZ", "xyz", "a", "a*b", "X-1", "X0 Y0", "X 0", "X Y", {0: "XY"}, {0: "Q"}, {-1: "X"}, [("X", 0), ("Z", 0)], [("Q", 0)]],
)
def test_invalid_input_is_rejected(paulis):
    with pytest.raises(ValueError):
        pe.PauliString(paulis)


def test_identity():
    assert pe.PauliString() == pe.PauliString("") == pe.PauliString({0: "I"}) == pe.PauliString("I3")
    assert pe.PauliString().x == []
    assert pe.PauliString().y == []
    assert str(pe.PauliString()) == "+1 : "


@pytest.mark.parametrize(
    ("a", "b", "phase", "result"),
    [
        ("X0", "Y0", 1, "Z0"),
        ("Y0", "Z0", 1, "X0"),
        ("Z0", "X0", 1, "Y0"),
        ("Y0", "X0", 3, "Z0"),
        ("Z0", "Y0", 3, "X0"),
        ("X0", "Z0", 3, "Y0"),
        ("X0", "X0", 0, ""),
        ("X0", "I0", 0, "X0"),
        ("X0 X1", "Y0 Y1", 2, "Z0 Z1"),
        ("X0 X1 X2", "Y0 Y1 Y2", 3, "Z0 Z1 Z2"),
        ("X0 X1 X2 X3", "Y0 Y1 Y2 Y3", 0, "Z0 Z1 Z2 Z3"),
    ],
)
def test_known_products(a, b, phase, result):
    assert pe.PauliString(a) * pe.PauliString(b) == pe.PauliString(result, phase=phase)


def test_product_of_every_two_qubit_pair_matches_the_matrices():
    labels = ["".join(p) for p in itertools.product("IXYZ", repeat=2)]
    for la, lb, pa, pb in itertools.product(labels, labels, range(4), range(4)):
        a = pe.PauliString(dict(enumerate(la)), phase=pa)
        b = pe.PauliString(dict(enumerate(lb)), phase=pb)
        assert np.array_equal(dense(a * b, 2), dense(a, 2) @ dense(b, 2)), (la, lb, pa, pb)


@pytest.mark.parametrize("n_qubits", [8, 64, 65, 128, 129, 512, 1000])
def test_product_matches_weighted_pauli_string(n_qubits):
    rng = random.Random(n_qubits)
    for _ in range(50):
        la, lb = random_labels(rng, n_qubits), random_labels(rng, n_qubits)
        word = pe.PauliString(la) * pe.PauliString(lb)
        string = pe.PauliString(1.0, la) * pe.PauliString(1.0, lb)
        assert word.x == string.x
        assert word.y == string.y
        assert word.phase_factor == string.coeff
        assert word.to_pauli_string() == string


def test_operands_of_different_length():
    short, long = pe.PauliString("X0"), pe.PauliString("Y0 Z200")
    assert short * long == pe.PauliString("Z0 Z200", phase=1)
    assert long * short == pe.PauliString("Z0 Z200", phase=3)


def test_trailing_words_are_trimmed_after_the_product():
    word = pe.PauliString("Z200")
    assert (word * word).x == []
    assert word * word == pe.PauliString()
    assert hash(word * word) == hash(pe.PauliString())


def test_group_properties():
    rng = random.Random(1)
    words = [pe.PauliString(random_labels(rng, 70), phase=rng.randrange(4)) for _ in range(30)]
    identity = pe.PauliString()
    for a, b, c in zip(words, words[1:], words[2:]):
        assert (a * b) * c == a * (b * c)
        assert a * identity == identity * a == a
        assert a * a.dagger() == identity
        assert (a * b).dagger() == b.dagger() * a.dagger()


def test_many_products_keep_an_exact_phase():
    x, y = pe.PauliString("X0"), pe.PauliString("Y0")
    word = pe.PauliString()
    for _ in range(100_000):
        word *= x
        word *= y
    # (XY)^100000 = (iZ)^100000 = i^100000 * I = I
    assert word == pe.PauliString()


def test_commutes_matches_the_phase_and_the_product_order():
    rng = random.Random(2)
    for _ in range(200):
        a = pe.PauliString(random_labels(rng, 130), phase=rng.randrange(4))
        b = pe.PauliString(random_labels(rng, 70), phase=rng.randrange(4))
        ab, ba = a * b, b * a
        assert ab.equals(ba)
        assert a.commutes(b) == (ab == ba)
        assert a.commutes(b) == b.commutes(a)
        if not a.commutes(b):
            assert ab == -ba


@pytest.mark.parametrize(("phase", "factor"), [(0, 1), (1, 1j), (2, -1), (3, -1j), (4, 1), (-1, -1j), (7, -1j)])
def test_phase_is_normalized(phase, factor):
    word = pe.PauliString("X0", phase=phase)
    assert word.phase in range(4)
    assert word.phase_factor == factor
    assert word.to_dictionary() == (factor, {0: "X"})


def test_negation_dagger_and_naked():
    word = pe.PauliString("X0", phase=1)
    assert (-word).phase == 3
    assert word.dagger().phase == 3
    assert word.naked() == pe.PauliString("X0")
    assert word.phase == 1


def test_equality_and_hash():
    assert pe.PauliString("X0") != pe.PauliString("X0", phase=1)
    assert pe.PauliString("X0").equals(pe.PauliString("X0", phase=1))
    assert not pe.PauliString("X0").equals(pe.PauliString("Y0"))
    assert hash(pe.PauliString("X0")) == hash(pe.PauliString({0: "X"}))
    assert len({pe.PauliString("X0"), pe.PauliString({0: "X"}), pe.PauliString("X0", phase=2)}) == 2


@pytest.mark.parametrize(("phase", "text"), [(0, "+1"), (1, "+i"), (2, "-1"), (3, "-i")])
def test_string_form(phase, text):
    word = pe.PauliString("Z2 X0", phase=phase)
    assert str(word) == repr(word) == word.to_string() == f"{text} : X(0)Z(2)"


def test_to_pauli_string():
    word = pe.PauliString("X0", phase=1)
    assert word.to_pauli_string() == pe.PauliString(1j, {0: "X"})
    assert word.to_pauli_string(0.5) == pe.PauliString(0.5j, {0: "X"})
    assert pe.PauliString("X0", phase=2).to_pauli_string(0.5) == pe.PauliString(-0.5, {0: "X"})
    assert pe.PauliString("X0", phase=2).to_pauli_string("a") == pe.PauliString("-a", {0: "X"})
    assert pe.PauliString("X0").to_pauli_string(pe.PauliString("a", {0: "X"}).coeff) == pe.PauliString("a", {0: "X"})


def test_to_pauli_string_enables_sums_and_commutators():
    a, b = pe.PauliString("X0"), pe.PauliString("Y0")
    commutator = a.to_pauli_string().commutator(b.to_pauli_string())
    assert commutator == pe.PauliString("Z0", phase=1).to_pauli_string(2.0)
    hamiltonian = pe.QubitHamiltonian([a.to_pauli_string(0.5), b.to_pauli_string(0.25)])
    assert len(hamiltonian) == 2


@pytest.mark.parametrize("paulis", ["", "X0", "Y63 Z64", "Z0 X700"])
def test_hash_round_trip(paulis):
    word = pe.PauliString(paulis, phase=1)
    key = word.get_hash()
    assert key == pe.PauliString(paulis).to_pauli_string().get_hash()
    assert _core.unweighted_pauliString.from_hash(key, 1) == word
    assert _core.unweighted_pauliString.from_hash(key) == word.naked()


def test_commutator_omits_the_factor_2():
    x, y = pe.PauliString("X0"), pe.PauliString("Y0")
    assert x.commutator(y) == pe.PauliString("Z0", phase=1)
    assert y.commutator(x) == pe.PauliString("Z0", phase=3)
    assert x.commutator(x) is None
    assert x.commutator(pe.PauliString("X1")) is None
    exact = x.to_pauli_string().commutator(y.to_pauli_string())
    assert x.commutator(y).to_pauli_string(2.0) == exact


# The bracket of unweighted strings is a Lie bracket: half the matrix commutator,
# antisymmetric, Jacobi identity, and the identity is central.


def all_words(n_qubits):
    rng = random.Random(7)
    return [pe.PauliString("".join(w), phase=rng.randrange(4)) for w in itertools.product("IXYZ", repeat=n_qubits)]


def bracket_matrix(a, b, n_qubits):
    c = a.commutator(b)
    return np.zeros((2**n_qubits, 2**n_qubits), dtype=complex) if c is None else dense(c, n_qubits)


def test_bracket_is_half_the_matrix_commutator_and_antisymmetric():
    n = 2
    for a, b in itertools.product(all_words(n), repeat=2):
        ma, mb = dense(a, n), dense(b, n)
        assert np.allclose(bracket_matrix(a, b, n), (ma @ mb - mb @ ma) / 2)
        assert np.allclose(bracket_matrix(a, b, n), -bracket_matrix(b, a, n))


def test_jacobi_identity():
    n = 2
    words = all_words(n)

    def nested(a, b, c):
        inner = b.commutator(c)
        return np.zeros((2**n, 2**n), dtype=complex) if inner is None else bracket_matrix(a, inner, n)

    for a, b, c in itertools.product(words, repeat=3):
        assert np.allclose(nested(a, b, c) + nested(b, c, a) + nested(c, a, b), 0)


def test_identity_is_central():
    identity = pe.PauliString("II")
    for a in all_words(2):
        assert identity.commutator(a) is None
        assert a.commutator(identity) is None


def test_closure_of_anti_hermitian_generators_stays_anti_hermitian():
    # i * P has an odd phase; the closure must neither leave the anti-Hermitian
    # elements nor produce the identity, and its size is the real dimension of
    # the algebra the matrices i * P generate (checked independently: 10).
    generators = [pe.PauliString(w, phase=1) for w in ("XI", "ZI", "XX", "IZ")]
    basis = lie_closure(generators)
    assert len(basis) == 10
    assert {t.phase for t in basis} <= {1, 3}
    assert all(t.dagger() == -t for t in basis)
    assert pe.PauliString("II") not in basis


def string_set(*terms):
    return pe.QubitHamiltonian(
        [pe.PauliString(t[0], phase=t[1]) if isinstance(t, tuple) else pe.PauliString(t) for t in terms]
    )


def lie_closure(generators):
    basis = pe.QubitHamiltonian(generators)
    while True:
        size = len(basis)
        basis += basis.commutator(basis)
        if len(basis) == size:
            return basis


def test_set_semantics_first_phase_wins():
    h = string_set("X0", ("X0", 2), "Y0")
    assert len(h) == 2
    assert h[0].phase == 0
    assert h.insert(pe.PauliString("X0", phase=1)) is False
    assert h.insert(pe.PauliString("Z0", phase=1)) is True
    assert h[2].phase == 1
    assert [str(t) for t in h] == ["+1 : X(0)", "+1 : Y(0)", "+i : Z(0)"]
    assert h + h == h
    assert string_set(("X0", 1)) + string_set("X0", "Y0") == string_set(("X0", 1), "Y0")


def test_in_place_union_keeps_the_object_and_its_indices():
    h = string_set("X0", "Y0")
    alias = h
    h += string_set("Y0", "Z0")
    assert h is alias
    assert [str(t) for t in h] == ["+1 : X(0)", "+1 : Y(0)", "+1 : Z(0)"]
    h += h
    assert len(h) == 3


def test_indices_are_stable_across_growth():
    rng = random.Random(3)
    strings, seen = [], set()
    while len(strings) < 2000:
        s = pe.PauliString(random_labels(rng, 8), phase=rng.randrange(4))
        if s.get_hash() not in seen:
            seen.add(s.get_hash())
            strings.append(s)
    h = _core.unweighted_QubitHamiltonian()
    for i, s in enumerate(strings):
        assert h.insert(s)
        assert h.index(s) == i
    assert all(h[i] == s and s in h and h.index(s) == i for i, s in enumerate(strings))
    assert h[-1] == strings[-1]
    assert h.terms() == strings
    assert h.paulistrings() == strings


def test_lookup_and_errors():
    h = string_set("X0", "Z200")
    assert pe.PauliString("X0", phase=3) in h
    assert pe.PauliString("Y0") not in h
    assert pe.PauliString("X900") not in h
    assert pe.PauliString() not in h
    assert pe.PauliString("X0") not in _core.unweighted_QubitHamiltonian()
    with pytest.raises(ValueError):
        h.index(pe.PauliString("Y0"))
    with pytest.raises(IndexError):
        h[2]
    with pytest.raises(IndexError):
        h[-3]


def test_identity_is_an_ordinary_element():
    h = string_set("X0", ("", 1))
    assert len(h) == 2
    assert pe.PauliString() in h
    assert h[1] == pe.PauliString("", phase=1)
    assert (string_set("Z200") * string_set("Z200"))[0] == pe.PauliString()


def test_strings_of_different_length_share_a_set():
    h = string_set("X0")
    h.insert(pe.PauliString("Y300"))
    h.insert(pe.PauliString("Z70"))
    assert [str(t) for t in h] == ["+1 : X(0)", "+1 : Y(300)", "+1 : Z(70)"]
    assert h.insert(pe.PauliString("X0")) is False
    assert h == string_set("Z70", "Y300", "X0")


@pytest.mark.parametrize("n_qubits", [6, 64, 65, 200])
def test_product_and_commutator_match_the_strings(n_qubits):
    rng = random.Random(n_qubits)
    a_terms = [pe.PauliString(random_labels(rng, n_qubits), phase=rng.randrange(4)) for _ in range(30)]
    b_terms = [pe.PauliString(random_labels(rng, n_qubits // 2 + 1), phase=rng.randrange(4)) for _ in range(30)]
    a, b = pe.QubitHamiltonian(a_terms), pe.QubitHamiltonian(b_terms)
    expected_product, expected_commutator = [], []
    for s in a.terms():
        for t in b.terms():
            expected_product.append(s * t)
            if not s.commutes(t):
                expected_commutator.append(s * t)
    assert (a * b).terms() == pe.QubitHamiltonian(expected_product).terms()
    assert a.commutator(b).terms() == pe.QubitHamiltonian(expected_commutator).terms()


def test_conversion_to_and_from_qubit_hamiltonian():
    h = string_set(("X0", 1), ("Y1", 2))
    weighted = h.to_qubit_hamiltonian(0.5)
    assert weighted == pe.QubitHamiltonian([pe.PauliString(0.5j, {0: "X"}), pe.PauliString(-0.5, {1: "Y"})])
    assert _core.unweighted_QubitHamiltonian(weighted) == string_set("X0", "Y1")
    symbolic = h.to_qubit_hamiltonian("a")
    assert _core.unweighted_QubitHamiltonian(symbolic) == string_set("X0", "Y1")


def test_commutator_operators_match_the_weighted_hamiltonian():
    rng = random.Random(5)
    a = pe.QubitHamiltonian([pe.PauliString(random_labels(rng, 6)) for _ in range(40)])
    b = pe.QubitHamiltonian([pe.PauliString(random_labels(rng, 6)) for _ in range(40)])
    unweighted = a.commutator(b)
    weighted = _core.unweighted_QubitHamiltonian(a.to_qubit_hamiltonian().commutator(b.to_qubit_hamiltonian()))
    # The weighted sum may cancel terms, a set of operators never does.
    assert all(t in unweighted for t in weighted)


@pytest.mark.parametrize("n_qubits", [2, 3, 4, 5])
def test_lie_closure_of_the_transverse_field_ising_chain_is_so_2n(n_qubits):
    generators = [pe.PauliString(f"Z{q}", phase=1) for q in range(n_qubits)]
    generators += [pe.PauliString(f"X{q} X{q + 1}", phase=1) for q in range(n_qubits - 1)]
    basis = lie_closure(generators)
    assert len(basis) == n_qubits * (2 * n_qubits - 1)
    assert len(basis + basis.commutator(basis)) == len(basis)
    assert all(t.phase % 2 == 1 for t in basis)  # i*P is anti-Hermitian, and the algebra keeps that


@pytest.mark.parametrize("n_qubits", [2, 3, 4])
def test_lie_closure_of_a_universal_set_is_su_2_to_the_n(n_qubits):
    generators = [pe.PauliString(f"{p}{q}", phase=1) for q in range(n_qubits) for p in "XZ"]
    generators += [pe.PauliString(f"Z{q} Z{q + 1}", phase=1) for q in range(n_qubits - 1)]
    assert len(lie_closure(generators)) == 4**n_qubits - 1


def hamiltonians():
    complex_h = pe.QubitHamiltonian([pe.PauliString(2.0, "XI"), pe.PauliString(0.5j, "ZZ")])
    symbolic_h = pe.QubitHamiltonian([pe.PauliString("a", "XI"), pe.PauliString("b", "ZZ")])
    return [(complex_h, _core.QubitHamiltonianComplex), (symbolic_h, _core.QubitHamiltonianSymbolic)]


@pytest.mark.parametrize("phase", range(4))
@pytest.mark.parametrize("index", range(2))
def test_unweighted_string_next_to_a_hamiltonian_counts_as_weighted(phase, index):
    h, kind = hamiltonians()[index]
    p = pe.PauliString("YI", phase=phase)
    w = pe.PauliString(1j**phase, "YI")
    for got, expected in [
        (h + p, h + w),
        (p + h, w + h),
        (h - p, h - w),
        (p - h, w - h),
        (h * p, h * w),
        (p * h, w * h),
        (h.commutator(p), h.commutator(w)),
        (p.commutator(h), w.commutator(h)),
    ]:
        assert isinstance(got, kind)
        assert got == expected


def test_commutator_of_a_string_and_a_hamiltonian_keeps_the_factor_2():
    h = pe.QubitHamiltonian([pe.PauliString(1.0, "X")])
    assert h.commutator(pe.PauliString("Y")) == pe.QubitHamiltonian([pe.PauliString(2j, "Z")])


def test_inner_product_and_lookup_next_to_a_hamiltonian():
    h = pe.QubitHamiltonian([pe.PauliString(2.0, "X"), pe.PauliString(0.5j, "Z")])
    p = pe.PauliString("Z", phase=1)
    assert h.inner_product(p) == h.inner_product(pe.PauliString(1j, "Z"))
    assert p.inner_product(h) == pe.PauliString(1j, "Z").inner_product(h)
    assert h.contains(p) == 0.5j
    assert p in h
    assert pe.PauliString("Y") not in h
    symbolic = pe.QubitHamiltonian([pe.PauliString("a", "X")])
    assert pe.PauliString("X") in symbolic
    assert pe.PauliString("Y") not in symbolic


def test_sum_cancels_next_to_a_hamiltonian():
    h = pe.QubitHamiltonian([pe.PauliString(1.0, "X")])
    assert len(h - pe.PauliString("X")) == 0
    assert (h + pe.PauliString("X")) == pe.QubitHamiltonian([pe.PauliString(2.0, "X")])


def test_two_unweighted_strings_have_no_sum():
    with pytest.raises(TypeError):
        pe.PauliString("X") + pe.PauliString("Y")


# Unweighted operands with every other kind of operator (Arithmetic.h)


def weighted_operands():
    cycle_sum = pe.PauliCycleSum([pe.PauliCycle(2, "XZ", 0.5), pe.PauliCycle(2, "YI", -1j)])
    orbit_sum = pe.PauliOrbit(2, 1, 0, 1) * 2
    return {
        "complex string": pe.PauliString(0.5 + 0.25j, "XZ"),
        "symbolic string": pe.PauliString("a", "XZ"),
        "complex hamiltonian": pe.QubitHamiltonian([pe.PauliString(2.0, "XI"), pe.PauliString(0.5j, "ZZ")]),
        "symbolic hamiltonian": pe.QubitHamiltonian([pe.PauliString("a", "XI"), pe.PauliString("b", "ZZ")]),
        "cycle": pe.PauliCycle(2, "XY", 1.5),
        "cycle sum": cycle_sum,
        "orbit": pe.PauliOrbit(2, 1, 1, 0),
        "orbit sum": orbit_sum,
    }


def unweighted_operands():
    s = pe.PauliString("YX", phase=1)
    h = pe.QubitHamiltonian([pe.PauliString("ZI", phase=3), pe.PauliString("XY", phase=2)])
    return {
        "string": (s, s.to_pauli_string(1.0)),
        "set": (h, h.to_qubit_hamiltonian(1.0)),
    }


OPERATIONS = {
    "+": lambda a, b: a + b,
    "-": lambda a, b: a - b,
    "*": lambda a, b: a * b,
    "commutator": lambda a, b: a.commutator(b),
    "inner_product": lambda a, b: a.inner_product(b),
}


@pytest.mark.parametrize("op", list(OPERATIONS))
@pytest.mark.parametrize("unweighted", list(unweighted_operands()))
@pytest.mark.parametrize("weighted", list(weighted_operands()))
def test_unweighted_operand_counts_with_coefficient_one_times_its_phase(weighted, unweighted, op):
    w = weighted_operands()[weighted]
    u, u_as_weighted = unweighted_operands()[unweighted]
    f = OPERATIONS[op]
    for got, expected in [(f(w, u), f(w, u_as_weighted)), (f(u, w), f(u_as_weighted, w))]:
        assert type(got) is type(expected)
        if op == "inner_product":  # symbolic numbers have no value comparison in Python
            assert str(got) == str(expected)
        else:
            assert got == expected


def test_upgrade_uses_coefficient_one_times_the_phase():
    h = pe.QubitHamiltonian([pe.PauliString("Z", phase=1)])
    c = pe.QubitHamiltonian([pe.PauliString(2.0, "X")])
    assert c + h == pe.QubitHamiltonian([pe.PauliString(2.0, "X"), pe.PauliString(1j, "Z")])
    assert c * pe.PauliString("Z") == pe.QubitHamiltonian([pe.PauliString(-2j, "Y")])


@pytest.mark.parametrize("scale", [2, -0.5, 1j])
def test_set_times_a_number_is_a_weighted_hamiltonian(scale):
    h = pe.QubitHamiltonian([pe.PauliString("Z", phase=1), pe.PauliString("X")])
    expected = h.to_qubit_hamiltonian(complex(scale))
    for got in (h * scale, scale * h):
        assert isinstance(got, _core.QubitHamiltonianComplex)
        assert got == expected


def test_set_and_string_stay_unweighted():
    x, y, z = pe.PauliString("X"), pe.PauliString("Y"), pe.PauliString("Z")
    h = pe.QubitHamiltonian([x, y])
    assert isinstance(h + z, _core.unweighted_QubitHamiltonian)
    assert [t.to_string() for t in h + z] == [x.to_string(), y.to_string(), z.to_string()]
    assert [t.to_string() for t in z + h] == [z.to_string(), x.to_string(), y.to_string()]
    # the first phase seen for an operator wins
    assert (h + pe.PauliString("X", phase=2))[0].phase == 0
    assert [t.to_string() for t in h * z] == [(x * z).to_string(), (y * z).to_string()]
    assert [t.to_string() for t in z * h] == [(z * x).to_string(), (z * y).to_string()]
    assert [t.to_string() for t in h.commutator(x)] == [(y * x).to_string()]
    assert [t.to_string() for t in x.commutator(h)] == [(x * y).to_string()]
    with pytest.raises(TypeError):
        h - z
