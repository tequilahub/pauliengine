"""Pauli decomposition of matrices (include/pauliengine/PauliDecomposition.h)."""

import itertools
import random

import pytest

import pauliengine as pe

np = pytest.importorskip("numpy")

PAULI = {
    "I": np.eye(2),
    "X": np.array([[0, 1], [1, 0]]),
    "Y": np.array([[0, -1j], [1j, 0]]),
    "Z": np.array([[1, 0], [0, -1]]),
}


def brute_force(M):
    """c_P = Tr(P M) / 2^n with qubit 0 as the leftmost tensor factor."""
    n = int(np.log2(M.shape[0]))
    out = {}
    for labels in itertools.product("IXYZ", repeat=n):
        P = np.array([[1.0]])
        for ch in labels:
            P = np.kron(P, PAULI[ch])
        c = np.trace(P @ M) / 2**n
        if abs(c) > 1e-14:
            out[tuple((q, ch) for q, ch in enumerate(labels) if ch != "I")] = c
    return out


def as_dict(h):
    return {tuple(sorted(ops.items())): complex(c) for c, ops in h.to_dictionary()}


def assert_close(d1, d2, tol=1e-12):
    for k in set(d1) | set(d2):
        assert abs(d1.get(k, 0) - d2.get(k, 0)) <= tol, k


@pytest.mark.parametrize("n", range(0, 6))
def test_matches_trace_formula(n):
    rng = np.random.default_rng(n)
    M = rng.normal(size=(2**n, 2**n)) + 1j * rng.normal(size=(2**n, 2**n))
    assert_close(as_dict(pe.QubitHamiltonian.from_matrix(M)), brute_force(M))


@pytest.mark.parametrize("n", [1, 3, 6, 8])
def test_roundtrip_through_to_matrix(n):
    rnd = random.Random(n)
    terms = [
        pe.PauliString(complex(rnd.uniform(-1, 1), rnd.uniform(-1, 1)),
                       {q: rnd.choice("XYZ") for q in rnd.sample(range(n), rnd.randint(1, n))})
        for _ in range(10)
    ] + [pe.PauliString(0.7, {n - 1: "Z"})]          # make sure the top qubit is used
    H = pe.QubitHamiltonian(terms)
    back = pe.QubitHamiltonian.from_matrix(np.array(H.to_matrix(ignore_unused_qubits=False)))
    assert len((back - H).simplify(1e-12)) == 0


@pytest.mark.parametrize(
    ("matrix", "expected"),
    [
        ([[1, 0], [0, 1]], {(): 1}),
        ([[0, 1], [1, 0]], {((0, "X"),): 1}),
        ([[0, -1j], [1j, 0]], {((0, "Y"),): 1}),
        ([[1, 0], [0, -1]], {((0, "Z"),): 1}),
        ([[1, 1], [1, -1]], {((0, "X"),): 1, ((0, "Z"),): 1}),
        ([[0, 1], [2, 3]], {(): 1.5, ((0, "X"),): 1.5, ((0, "Y"),): -0.5j, ((0, "Z"),): -1.5}),
    ],
)
def test_single_qubit_examples(matrix, expected):
    assert_close(as_dict(pe.QubitHamiltonian.from_matrix(matrix)), expected)


def test_qubit_order_matches_to_matrix():
    # X on qubit 0 = X (x) I: qubit 0 is the most significant index bit
    M = np.kron(PAULI["X"], np.eye(2))
    assert as_dict(pe.QubitHamiltonian.from_matrix(M)) == {((0, "X"),): 1}


def test_hermitian_matrix_gives_real_coefficients():
    rng = np.random.default_rng(1)
    A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    H = pe.QubitHamiltonian.from_matrix(A + A.conj().T, threshold=1e-14)
    assert H.is_hermitian()


def test_list_input_and_real_dtype():
    from_list = as_dict(pe.QubitHamiltonian.from_matrix([[1, 2], [3, 4]]))
    from_int_array = as_dict(pe.QubitHamiltonian.from_matrix(np.array([[1, 2], [3, 4]])))
    assert_close(from_list, from_int_array)


def test_core_list_overload_without_numpy_path():
    h = pe._core.QubitHamiltonianComplex.from_matrix([[0, 1], [1, 0]])
    assert as_dict(h) == {((0, "X"),): 1}


def test_input_is_not_modified():
    M = np.arange(16, dtype=complex).reshape(4, 4)
    before = M.copy()
    pe.QubitHamiltonian.from_matrix(M)
    assert np.array_equal(M, before)


def test_threshold_drops_small_terms():
    M = np.array([[1, 1e-13], [1e-13, -1]], dtype=complex)
    assert len(pe.QubitHamiltonian.from_matrix(M)) == 2          # Z and a tiny X
    assert as_dict(pe.QubitHamiltonian.from_matrix(M, threshold=1e-12)) == {((0, "Z"),): 1}


def test_zero_matrix_is_empty():
    assert len(pe.QubitHamiltonian.from_matrix(np.zeros((4, 4)))) == 0


@pytest.mark.parametrize(
    "bad",
    [np.zeros((3, 3)), np.zeros((2, 4)), np.zeros((0, 0)), np.zeros(4), [[1, 2], [3]]],
)
def test_invalid_shapes_raise(bad):
    with pytest.raises(ValueError):
        pe.QubitHamiltonian.from_matrix(bad)
