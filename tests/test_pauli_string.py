import math

import pytest

import pauliengine as pe
import pauliengine._core as _core


@pytest.mark.parametrize(
    ("pauli_string", "coeff", "expected_x", "expected_y"),
    [
        ({0: "I", 1: "Z", 2: "X"}, 5.0, [6], [2]),  # using dict input
        (
            {10: "Z", 2: "X"},
            6.0,
            [2**10 + 2**2],
            [2**10],
        ),  # using dict input with non-consecutive qubits
    ],
)
def test_pauli_string_valid_pauli_string_input(pauli_string, coeff, expected_x, expected_y):
    ps = pe.PauliString(coeff, pauli_string)
    assert ps.coeff == coeff
    assert ps.x == expected_x
    assert ps.y == expected_y

def test_pauli_string_valid_input_open_fermion():
    coeff = 3.0j
    ps = pe.PauliString((coeff, [("Z", 0), ("X", 2), ("Y", 5)]))
    expected_x = [1 + 2**2]
    expected_y = [1 + 2**5]
    assert ps.coeff == coeff
    assert expected_x == ps.x
    assert expected_y == ps.y


@pytest.mark.xfail
@pytest.mark.parametrize(
    "invalid_pauli_string",
    [
        {0: "A", 2: "B", 3: "Z"},  # invalid character
        {0: "i", 2: "z", 3: "y"},  # lowercase character
    ],
)
def test_pauli_string_invalid_pauli_string_input(invalid_pauli_string):
    with pytest.raises(ValueError):  # noqa: PT011
        pe.PauliString(1.0, invalid_pauli_string)


@pytest.mark.parametrize(
    ("pauli_string", "other", "expected_product"),
    [
        (
            pe.PauliString(1.0, {0: "X", 1: "Z", 2: "Y", 3: "Z"}),
            pe.PauliString(2.0, {0: "Z", 1: "Z", 3: "X"}),
            pe.PauliString(2.0, {0: "Y", 2: "Y", 3: "Y"}),
        ),  # XZYZ * ZZIX = YIYY
        (
            pe.PauliString(2.0j, {0: "Z", 2: "Z"}),
            pe.PauliString(3.0, {0: "Y", 2: "Y"}),
            pe.PauliString(-6.0j, {0: "X", 2: "X"}),
        ),  # XZ * ZX = YY
        (
            pe.PauliString(-4.0, {0: "X", 1: "Y"}),
            1.0 - 2.0j,
            pe.PauliString(-4.0 + 8.0j, {0: "X", 1: "Y"}),
        ),  # -4.0XY * 1.0 - 2.0j
    ],
)

class TestPauliStringMultiplication:
    def test_pauli_string_multiplication(self, pauli_string, other, expected_product):
        product = pauli_string * other
        assert expected_product == product

    def test_pauli_string_multiplication_inplace(
        self, pauli_string, other, expected_product
    ):
        pauli_string *= other
        assert expected_product == pauli_string


def test_SymbolicAndNonSymbolicMultiplication():
    ps1 = pe.PauliString("a", {0: "Z", 1: "Z", 3: "X"})
    ps2 = pe.PauliString(2.0j, {0: "X", 1: "Z", 2: "Y", 3: "Z"})
    product = ps1 * ps2
    result = pe.PauliString("2a", {0: "Y", 2: "Y", 3: "Y"}) * 1.0j
    assert result == product

def test_paulistring_is_all_z():
    ps1 = pe.PauliString(1.0, {0: "Z", 1: "Z", 500: "Z"})
    assert ps1.is_all_z() == True

def test_paulistring_is_all_z_not_all_z():
    ps1 = pe.PauliString(1.0, {5: "Y", 1: "X", 500: "Z"})
    assert ps1.is_all_z() == False

def test_pauli_string_to_dictionary():
    pauli_data = {2: "Z", 5: "Y", 18: "Y"}
    coeff = 2.0j
    ps = pe.PauliString((coeff, pauli_data))
    result = ps.to_dictionary()
    assert result[0] == coeff
    assert result[1] == pauli_data


@pytest.mark.parametrize(
    ("first_pauli_string", "second_pauli_string"),
    [
        (
            pe.PauliString((1.0, {0: "X", 1: "Y", 3:"X"})),
            pe.PauliString((1.0, {0: "X", 1: "Y", 3:"X"})),
           # pe.PauliString((0.0, {})),
        ),
        (
            pe.PauliString((3.0, {0: "X", 1: "Y", 400: "X"})),
            pe.PauliString((1.0j, {0: "X", 1: "X", 3:"X", 1000: "Z"})),
           # pe.PauliString((6.0, {1:"Z", 3: "X", 400: "X", 1000: "Z"})),
        ),
    ],
 )

def test_pauli_string_commutator(first_pauli_string, second_pauli_string):
    ps = first_pauli_string.commutator(second_pauli_string)
    ps_AB = first_pauli_string * second_pauli_string
    ps_BA = second_pauli_string * first_pauli_string
    assert (ps_AB.coeff - ps_BA.coeff) == ps.coeff
    assert ps_AB.x == ps.x
    assert ps_AB.y == ps.y
    assert ps_BA.x == ps.x
    assert ps_BA.y == ps.y


# --- Skalarmultiplikation: ps * scalar und scalar * ps ---

@pytest.mark.parametrize(
    ("scalar", "expected_coeff"),
    [
        (2,       2.0 + 0j),
        (2.0,     2.0 + 0j),
        (2.0+0j,  2.0 + 0j),
        (2.0+1j,  2.0 + 1j),
    ],
)
def test_pauli_string_mul_scalar(scalar, expected_coeff):
    ps = pe.PauliString(1.0, {0: "X", 1: "Z"})
    result = ps * scalar
    assert result.coeff == expected_coeff
    assert result.equals(ps)


@pytest.mark.parametrize(
    ("scalar", "expected_coeff"),
    [
        (2,       2.0 + 0j),
        (2.0,     2.0 + 0j),
        (2.0+0j,  2.0 + 0j),
        (2.0+1j,  2.0 + 1j),
    ],
)
def test_pauli_string_rmul_scalar(scalar, expected_coeff):
    ps = pe.PauliString(1.0, {0: "X", 1: "Z"})
    result = scalar * ps
    assert result.coeff == expected_coeff
    assert result.equals(ps)


def test_pauli_string_mul_rmul_commutative():
    ps = pe.PauliString(3.0, {0: "Y", 2: "Z"})
    assert (ps * 4) == (4 * ps)
    assert (ps * 2.5) == (2.5 * ps)
    assert (ps * (1+2j)) == ((1+2j) * ps)


def test_pauli_string_symbolic_rmul_scalar():
    ps = pe.PauliString("a", {0: "X"})
    result_right = ps * 2.0
    result_left  = 2.0 * ps
    assert result_right == result_left


# --- to_string / __str__ / __repr__ ---

def test_pauli_string_to_string():
    ps = pe.PauliString(1.0, {0: "X", 1: "Z"})
    s = ps.to_string()
    assert isinstance(s, str)
    assert len(s) > 0


def test_pauli_string_str_repr_consistent():
    ps = pe.PauliString(2.0, {0: "Y"})
    assert str(ps) == repr(ps)


# --- get_coeff ---

def test_pauli_string_get_coeff():
    ps = pe.PauliString(3.0, {0: "X"})
    assert ps.get_coeff() == 3.0 + 0j
    assert ps.get_coeff() == ps.coeff


# --- get_pauli_at_index ---

@pytest.mark.parametrize(
    ("pauli_map", "index", "expected"),
    [
        ({0: "X", 1: "Y", 2: "Z"}, 0, "X"),
        ({0: "X", 1: "Y", 2: "Z"}, 1, "Y"),
        ({0: "X", 1: "Y", 2: "Z"}, 2, "Z"),
        ({0: "X"}, 5, "I"),
    ],
)
def test_pauli_string_get_pauli_at_index(pauli_map, index, expected):
    ps = pe.PauliString(1.0, pauli_map)
    assert ps.get_pauli_at_index(index) == expected


# equals() compares only the operator part (x, y).
# == (operator==) compares x, y AND the coefficient.

def test_pauli_string_equals_vs_eq():
    ps1 = pe.PauliString(1.0, {0: "X"})
    ps2 = pe.PauliString(2.0, {0: "X"})
    ps3 = pe.PauliString(1.0, {0: "X"})
    assert ps1.equals(ps2)
    assert not (ps1 == ps2)
    assert ps1 == ps3


# --- copy ---

def test_pauli_string_copy():
    ps = pe.PauliString(5.0, {0: "X", 1: "Y"})
    ps_copy = ps.copy()
    assert ps == ps_copy
    ps_copy *= 2.0
    assert not (ps == ps_copy)


# --- set_coeff ---

def test_pauli_string_set_coeff():
    ps = pe.PauliString(1.0, {0: "X"})
    ps.set_coeff(3.0 + 0j)
    assert ps.get_coeff() == 3.0 + 0j


def test_pauli_string_map_qubits_complete():
    ps = pe.PauliString(1.0, {0: "X", 1: "Y"})
    ps_mapped = ps.map_qubits({0: 2, 1: 3})
    assert ps_mapped.get_pauli_at_index(2) == "X"
    assert ps_mapped.get_pauli_at_index(3) == "Y"
    assert ps_mapped.get_pauli_at_index(0) == "I"
    assert ps_mapped.get_pauli_at_index(1) == "I"


def test_pauli_string_map_qubits_partial():
    """Unmapped qubits keep their original index."""
    ps = pe.PauliString(1.0, {0: "X", 1: "Y", 2: "Z"})
    ps_mapped = ps.map_qubits({0: 5})
    assert ps_mapped.get_pauli_at_index(5) == "X"
    assert ps_mapped.get_pauli_at_index(1) == "Y"
    assert ps_mapped.get_pauli_at_index(2) == "Z"
    assert ps_mapped.get_pauli_at_index(0) == "I"


def test_pauli_string_map_qubits_implicit_swap():
    """2 -> 17, but 17 is already occupied -> implicit swap: 17 -> 2."""
    ps = pe.PauliString(1.0, {2: "X", 17: "Z"})
    ps_mapped = ps.map_qubits({2: 17})
    assert ps_mapped.get_pauli_at_index(17) == "X"
    assert ps_mapped.get_pauli_at_index(2) == "Z"


def test_pauli_string_map_qubits_implicit_swap_leaves_others_unchanged():
    """Swap only touches the affected qubits."""
    ps = pe.PauliString(1.0, {0: "X", 2: "Y", 17: "Z"})
    ps_mapped = ps.map_qubits({2: 17})
    assert ps_mapped.get_pauli_at_index(17) == "Y"
    assert ps_mapped.get_pauli_at_index(2) == "Z"
    assert ps_mapped.get_pauli_at_index(0) == "X"


def test_pauli_string_map_qubits_no_implicit_swap_if_target_explicitly_mapped():
    """Explicit 17 -> 5 disables the implicit swap from the 2 -> 17 mapping."""
    ps = pe.PauliString(1.0, {2: "X", 17: "Z"})
    ps_mapped = ps.map_qubits({2: 17, 17: 5})
    assert ps_mapped.get_pauli_at_index(17) == "X"
    assert ps_mapped.get_pauli_at_index(5) == "Z"
    assert ps_mapped.get_pauli_at_index(2) == "I"


def test_pauli_string_map_qubits_chain():
    """Explicit chain 1->4, 4->5, 11->13, 13->2 with no implicit swap."""
    ps = pe.PauliString(3.0j, {1: "X", 4: "Z", 11: "Y", 13: "X"})
    ps_mapped = ps.map_qubits({1: 4, 4: 5, 11: 13, 13: 2})
    assert ps_mapped.get_pauli_at_index(4) == "X"
    assert ps_mapped.get_pauli_at_index(5) == "Z"
    assert ps_mapped.get_pauli_at_index(13) == "Y"
    assert ps_mapped.get_pauli_at_index(2) == "X"


def test_pauli_string_qubits():
    ps = pe.PauliString(1.0, {0: "X", 2: "Z"})
    qs = ps.qubits()
    assert 0 in qs
    assert 2 in qs


def test_pauli_string_openfermion_format():
    ps = pe.PauliString(1.0, {0: "X", 1: "Z"})
    of_format = ps.key_openfermion()
    assert isinstance(of_format, list)
    qubit_indices = [pair[1] for pair in of_format]
    assert 0 in qubit_indices
    assert 1 in qubit_indices


def test_pauli_string_space_separated_string_constructor():
    ps = pe.PauliString(1.0, "X0 Y1 Z2")
    assert ps.get_pauli_at_index(0) == "X"
    assert ps.get_pauli_at_index(1) == "Y"
    assert ps.get_pauli_at_index(2) == "Z"


def test_pauli_string_default_constructor():
    import pauliengine._core as _core
    ps = _core.PauliStringComplex()
    assert ps.coeff == 0j
    assert ps.is_zero


# --- PauliStringSymbolic ---

class TestPauliStringSymbolic:
    def test_symbolic_construction(self):
        ps = pe.PauliString("a", {0: "X"})
        assert "a" in str(ps.coeff)

    def test_symbolic_to_complex(self):
        import pauliengine._core as _core
        expr = _core.Expression("2.0")
        result = _core.PauliStringSymbolic.to_complex(expr)
        assert result == 2.0 + 0j

    def test_symbolic_get_pauli_at_index(self):
        ps = pe.PauliString("a", {0: "X", 1: "Z"})
        assert ps.get_pauli_at_index(0) == "X"
        assert ps.get_pauli_at_index(1) == "Z"

    def test_symbolic_map_qubits(self):
        ps = pe.PauliString("a", {0: "Y"})
        ps_mapped = ps.map_qubits({0: 3})
        assert ps_mapped.get_pauli_at_index(3) == "Y"
        assert ps_mapped.get_pauli_at_index(0) == "I"

    def test_symbolic_scalar_multiplication(self):
        ps = pe.PauliString("a", {0: "X"})
        result = ps * (2.0 + 0j)
        assert "a" in str(result.coeff)

    def test_symbolic_imul_pauli_string(self):
        ps1 = pe.PauliString("a", {0: "X"})
        ps2 = pe.PauliString("b", {0: "Z"})
        ps1 *= ps2
        assert "a" in str(ps1.coeff) and "b" in str(ps1.coeff)

    def test_symbolic_imul_complex_scalar(self):
        ps = pe.PauliString("a", {0: "X"})
        ps *= (2.0 + 0j)
        assert "a" in str(ps.coeff)

    def test_symbolic_rmul_int(self):
        ps = pe.PauliString("a", {0: "X"})
        result = 3 * ps
        assert result.equals(ps)

    def test_symbolic_rmul_float(self):
        ps = pe.PauliString("a", {0: "X"})
        result = 2.5 * ps
        assert result.equals(ps)


def test_pauli_string_xx_gives_identity_operator():
    """X * X = I: same operator squared gives an empty x/y vector."""
    ps = pe.PauliString(1.0, {0: "X"})
    result = ps * ps
    assert result.get_pauli_at_index(0) == "I"


def test_pauli_string_mul_zero_scalar():
    ps = pe.PauliString(3.0, {0: "X", 1: "Z"})
    result = ps * 0
    assert result.coeff == 0j
    assert result.equals(ps)


def test_pauli_string_mul_negative_scalar():
    ps = pe.PauliString(2.0, {0: "Y"})
    result = ps * (-1)
    assert result.coeff == -2.0 + 0j
    assert result.equals(ps)


def test_pauli_string_mul_zero_from_left():
    ps = pe.PauliString(3.0, {0: "Z"})
    result = 0 * ps
    assert result.coeff == 0j


 
# Comprehensive commutator tests
 

@pytest.mark.parametrize(
    ("op_a", "op_b", "expected_op", "expected_coeff"),
    [
        ("X", "Y", "Z", 2j),
        ("Y", "X", "Z", -2j),
        ("Y", "Z", "X", 2j),
        ("Z", "Y", "X", -2j),
        ("Z", "X", "Y", 2j),
        ("X", "Z", "Y", -2j),
    ],
)
def test_pauli_string_commutator_exact_values(op_a, op_b, expected_op, expected_coeff):
    """[A, B] exact value for all non-trivial single-qubit Pauli pairs."""
    ps_a = pe.PauliString(1.0, {0: op_a})
    ps_b = pe.PauliString(1.0, {0: op_b})
    result = ps_a.commutator(ps_b)
    expected = pe.PauliString(expected_coeff, {0: expected_op})
    assert result == expected


def test_pauli_string_commutator_with_identity_is_zero():
    """[A, I] = 0 and [A, B on a different qubit] = 0 (single check)."""
    ps_x = pe.PauliString(1.0, {0: "X"})
    assert ps_x.commutator(pe.PauliString(1.0, {0: "I"})).coeff == 0j
    assert ps_x.commutator(pe.PauliString(1.0, {1: "Z"})).coeff == 0j


def test_pauli_string_commutator_coefficient_scaling():
    """[2X, 3Z] = -12iY covers both coefficient scaling and antisymmetry sign."""
    result = pe.PauliString(2.0, {0: "X"}).commutator(pe.PauliString(3.0, {0: "Z"}))
    assert result == pe.PauliString(-12j, {0: "Y"})


def test_pauli_string_commutator_imaginary_coefficient():
    """[iX, Y] = i * 2iZ = -2Z."""
    result = pe.PauliString(1j, {0: "X"}).commutator(pe.PauliString(1.0, {0: "Y"}))
    assert result == pe.PauliString(-2.0 + 0j, {0: "Z"})


def test_pauli_string_commutator_multiqubit_anticommuting_pairs():
    """One anticommuting pair -> nonzero, two anticommuting pairs -> zero."""
    # [X(0)Z(1), Z(0)] = -2i Y(0)Z(1)
    one_pair = pe.PauliString(1.0, {0: "X", 1: "Z"}).commutator(pe.PauliString(1.0, {0: "Z"}))
    assert one_pair == pe.PauliString(-2j, {0: "Y", 1: "Z"})
    # [X(0)X(1), Z(0)Z(1)] = 0 (two pairs)
    two_pairs = pe.PauliString(1.0, {0: "X", 1: "X"}).commutator(pe.PauliString(1.0, {0: "Z", 1: "Z"}))
    assert two_pairs.coeff == 0j


def test_pauli_string_commutator_symbolic_contains_both_symbols():
    """[a*X, b*Z] result coefficient contains both symbols."""
    result = pe.PauliString("a", {0: "X"}).commutator(pe.PauliString("b", {0: "Z"}))
    coeff_str = str(result.coeff)
    assert "a" in coeff_str and "b" in coeff_str


def test_pauli_string_commutator_large_qubit_indices():
    """Commutator works correctly with high qubit indices (cross-word boundary)."""
    result = pe.PauliString(1.0, {500: "X"}).commutator(pe.PauliString(1.0, {500: "Z"}))
    assert result == pe.PauliString(-2j, {500: "Y"})


 
# PauliString negation
 

def test_pauli_string_complex_neg_negates_coeff():
    ps = pe.PauliString(3.0, {0: "X"})
    neg = -ps
    assert neg.coeff == pytest.approx(-3.0 + 0j)
    assert neg.equals(ps)


def test_pauli_string_complex_neg_preserves_operator():
    ps = pe.PauliString(1.0, {0: "X", 1: "Z"})
    neg = -ps
    assert neg.get_pauli_at_index(0) == "X"
    assert neg.get_pauli_at_index(1) == "Z"


def test_pauli_string_complex_double_neg_is_identity():
    ps = pe.PauliString(2.0, {0: "Y"})
    assert (- (-ps)) == ps


def test_pauli_string_symbolic_neg():
    ps = pe.PauliString("a", {0: "X"})
    neg = -ps
    assert isinstance(neg, _core.PauliStringSymbolic)
    assert neg.equals(ps)


 
# Inspection helpers: naked, size, len, count_y, key_openfermion
 

class TestPauliStringInspection:
    def test_naked_zeros_coefficient(self):
        ps = pe.PauliString(3.0 + 1.5j, {0: "X", 1: "Z"})
        naked = ps.naked()
        assert naked.coeff == 1.0 + 0j
        assert naked.get_pauli_at_index(0) == "X"
        assert naked.get_pauli_at_index(1) == "Z"

    def test_naked_symbolic(self):
        ps = pe.PauliString("a", {0: "X"})
        naked = ps.naked()
        assert _core.PauliStringSymbolic.to_complex(naked.coeff) == 1.0 + 0j

    def test_size_returns_count_of_non_identity_qubits(self):
        ps = pe.PauliString(1.0, {0: "X", 1: "Z", 3: "Y"})
        assert ps.size() == 3

    def test_len_equals_size(self):
        ps = pe.PauliString(1.0, {0: "X", 5: "Z"})
        assert len(ps) == ps.size()

    def test_size_identity_is_zero(self):
        ps = pe.PauliString(1.0, {0: "I"})
        assert ps.size() == 0

    def test_count_y(self):
        ps = pe.PauliString(1.0, {0: "X", 1: "Y", 2: "Y", 3: "Z"})
        assert ps.count_y() == 2

    def test_count_y_none(self):
        ps = pe.PauliString(1.0, {0: "X", 1: "Z"})
        assert ps.count_y() == 0

    def test_key_openfermion_returns_pairs(self):
        ps = pe.PauliString(1.0, {0: "X", 2: "Z"})
        key = ps.key_openfermion()
        seen = {(p, q) for p, q in key}
        assert ('X', 0) in seen
        assert ('Z', 2) in seen


 
# trace_out_qubits — purity (no in-place mutation, idempotent)
 

class TestPauliStringTraceOutPurity:
    def test_trace_does_not_mutate_receiver(self):
        ps = pe.PauliString(2.0, {0: "Z", 1: "X"})
        x_before = list(ps.x)
        y_before = list(ps.y)
        coeff_before = ps.coeff
        _ = ps.trace_out_qubits([0], [0])
        assert list(ps.x) == x_before
        assert list(ps.y) == y_before
        assert ps.coeff == coeff_before

    def test_trace_is_idempotent(self):
        """Calling trace twice on the same PauliString must give the same result."""
        ps = pe.PauliString(2.0, {0: "Z", 1: "X"})
        r1 = ps.trace_out_qubits([0], [0])
        r2 = ps.trace_out_qubits([0], [0])
        assert r1.coeff == r2.coeff
        assert list(r1.x) == list(r2.x)
        assert list(r1.y) == list(r2.y)

    def test_trace_length_mismatch_raises(self):
        ps = pe.PauliString(1.0, {0: "Z"})
        with pytest.raises(Exception):
            ps.trace_out_qubits([0, 1], [0])


 
# trace_out_qubits general-state overload: matches tequila's
# <psi|P|psi> = |a|^2 <0|P|0> + a*b <0|P|1> + ab* <1|P|0> + |b|^2 <1|P|1>
 

# Flattened (row-major) 2x2 Pauli matrices.
_PAULI_MATRICES = {
    "I": (1 + 0j, 0j, 0j, 1 + 0j),
    "X": (0j, 1 + 0j, 1 + 0j, 0j),
    "Y": (0j, -1j, 1j, 0j),
    "Z": (1 + 0j, 0j, 0j, -1 + 0j),
}


def _tequila_expectation(op: str, a: complex, b: complex) -> complex:
    """Faithful reproduction of tequila's <psi|op|psi> via the flattened
    matrix · density-vector formula."""
    vec = (
        abs(a) ** 2,
        a.conjugate() * b,
        b.conjugate() * a,
        abs(b) ** 2,
    )
    return complex(sum(v * m for v, m in zip(vec, _PAULI_MATRICES[op], strict=True)))


_INV_SQRT2 = 1 / math.sqrt(2)


class TestPauliStringTraceGeneralStates:

    @pytest.mark.parametrize(
        ("op", "a", "b"),
        [
            # Each Pauli paired with a state that exercises a non-trivial branch
            # of the |a|^2/<0|.|1>/|b|^2 expansion.
            ("I", _INV_SQRT2, _INV_SQRT2 * 1j),       # I on |+i>
            ("X", _INV_SQRT2, _INV_SQRT2),            # <+|X|+> = 1
            ("Y", _INV_SQRT2, _INV_SQRT2 * 1j),       # <+i|Y|+i> = 1
            ("Z", 0.6 + 0j, 0.8 + 0j),                # real-amplitude Z
            ("X", 0.5 + 0.5j, 0.5 - 0.5j),            # generic complex on X
            ("Y", 0.5 + 0.5j, 0.5 - 0.5j),            # generic complex on Y
        ],
    )
    def test_matches_tequila_formula(self, op, a, b):
        ps = pe.PauliString(2.5 + 0j, {0: op})
        traced = ps.trace_out_qubits([0], [(complex(a), complex(b))])
        expected = (2.5 + 0j) * _tequila_expectation(op, complex(a), complex(b))
        assert complex(traced.coeff) == pytest.approx(expected, abs=1e-12)

    def test_multi_qubit_product_of_expectations(self):
        ps = pe.PauliString(1.0, {0: "X", 1: "Z"})
        states = [(_INV_SQRT2, _INV_SQRT2), (0.0 + 0j, 1.0 + 0j)]
        assert complex(ps.trace_out_qubits([0, 1], states).coeff) == pytest.approx(-1.0 + 0j)

    def test_basis_state_overload_matches_int_overload(self):
        ps = pe.PauliString(1.5 + 0j, {0: "Z", 1: "X"})
        r_int = ps.trace_out_qubits([0], [0])
        r_amp = ps.trace_out_qubits([0], [(1.0 + 0j, 0.0 + 0j)])
        assert complex(r_int.coeff) == pytest.approx(complex(r_amp.coeff))
        assert list(r_int.x) == list(r_amp.x)
        assert list(r_int.y) == list(r_amp.y)

    def test_general_states_length_mismatch_raises(self):
        ps = pe.PauliString(1.0, {0: "Z"})
        with pytest.raises(Exception):
            ps.trace_out_qubits([0, 1], [(1.0 + 0j, 0.0 + 0j)])

    def test_unnormalized_state_scales_factor(self):
        ps = pe.PauliString(1.0, {0: "Z"})
        traced = ps.trace_out_qubits([0], [(2.0 + 0j, 0.0 + 0j)])
        assert complex(traced.coeff) == pytest.approx(4.0 + 0j)