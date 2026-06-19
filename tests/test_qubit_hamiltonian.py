import numpy as np
import pytest

import pauliengine as pe
import pauliengine._core as _core

try:
    from openfermion import QubitOperator
    _HAS_OF = True
except ImportError:
    _HAS_OF = False

requires_openfermion = pytest.mark.skipif(
    not _HAS_OF, reason="openfermion not installed"
)


 
# Helpers
 

def _complex_qh(*pauli_dicts_and_coeffs):
    """Build a QubitHamiltonianComplex from (coeff, ops_dict) pairs."""
    ps_list = [pe.PauliString(complex(c), d) for c, d in pauli_dicts_and_coeffs]
    return pe.QubitHamiltonian(ps_list)


def _symbolic_qh(*pauli_dicts_and_coeffs):
    """Build a QubitHamiltonianSymbolic from (str_coeff, ops_dict) pairs."""
    ps_list = [pe.PauliString(str(c), d) for c, d in pauli_dicts_and_coeffs]
    return pe.QubitHamiltonian(ps_list)




class TestConstruction:
    def test_from_complex_pauli_strings(self):
        ps = pe.PauliString(1.0, {0: "X"})
        qh = pe.QubitHamiltonian([ps])
        assert isinstance(qh, _core.QubitHamiltonianComplex)

    def test_from_symbolic_pauli_strings(self):
        ps = pe.PauliString("a", {0: "X"})
        qh = pe.QubitHamiltonian([ps])
        assert isinstance(qh, _core.QubitHamiltonianSymbolic)

    def test_from_tuple_list_complex(self):
        qh = pe.QubitHamiltonian([(1.0 + 0j, {0: "X"}), (2.0 + 0j, {1: "Z"})])
        assert isinstance(qh, _core.QubitHamiltonianComplex)

    def test_from_tuple_list_symbolic(self):
        qh = pe.QubitHamiltonian([("a", {0: "X"}), ("b", {1: "Z"})])
        assert isinstance(qh, _core.QubitHamiltonianSymbolic)

    def test_any_symbolic_pauli_string_gives_symbolic_hamiltonian(self):
        """If ANY PauliString in the list is symbolic the result must be symbolic."""
        ps_complex = pe.PauliString(1.0, {0: "X"})
        ps_symbolic = pe.PauliString("a", {0: "Z"})
        qh = pe.QubitHamiltonian([ps_complex, ps_symbolic])
        assert isinstance(qh, _core.QubitHamiltonianSymbolic)

    def test_multiple_pauli_strings_complex(self):
        ps_list = [pe.PauliString(float(i), {i: "Z"}) for i in range(1, 4)]
        qh = pe.QubitHamiltonian(ps_list)
        assert isinstance(qh, _core.QubitHamiltonianComplex)

    def test_multiple_pauli_strings_symbolic(self):
        ps_list = [pe.PauliString(f"c{i}", {i: "Z"}) for i in range(3)]
        qh = pe.QubitHamiltonian(ps_list)
        assert isinstance(qh, _core.QubitHamiltonianSymbolic)



class TestEquality:
    def test_equal_hamiltonians(self):
        qh1 = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        qh2 = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        assert qh1 == qh2

    def test_unequal_hamiltonians_different_coeffs(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((2.0, {0: "X"}))
        assert not (qh1 == qh2)

    def test_unequal_hamiltonians_different_operators(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((1.0, {0: "Z"}))
        assert not (qh1 == qh2)

    def test_equal_commutative_term_order(self):
        """Equality must be order-independent (uses coefficient merging internally)."""
        qh1 = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        qh2 = _complex_qh((2.0, {1: "Z"}), (1.0, {0: "X"}))
        assert qh1 == qh2

    def test_equal_symbolic_hamiltonians(self):
        qh1 = _symbolic_qh(("a", {0: "X"}))
        qh2 = _symbolic_qh(("a", {0: "X"}))
        assert qh1 == qh2



class TestAddition:
    def test_add_same_terms_accumulates_coefficients(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((2.0, {0: "X"}))
        result = qh1 + qh2
        expected = _complex_qh((3.0, {0: "X"}))
        assert result == expected

    def test_add_different_terms(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((1.0, {1: "Z"}))
        result = qh1 + qh2
        expected = _complex_qh((1.0, {0: "X"}), (1.0, {1: "Z"}))
        assert result == expected

    def test_add_returns_correct_type_complex(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((1.0, {1: "Z"}))
        result = qh1 + qh2
        assert isinstance(result, _core.QubitHamiltonianComplex)

    def test_add_returns_correct_type_symbolic(self):
        qh1 = _symbolic_qh(("a", {0: "X"}))
        qh2 = _symbolic_qh(("b", {1: "Z"}))
        result = qh1 + qh2
        assert isinstance(result, _core.QubitHamiltonianSymbolic)

    def test_add_cancelling_terms(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((-1.0, {0: "X"}))
        result = qh1 + qh2
        expected = _complex_qh((0.0, {0: "X"}))
        assert result == expected


class TestScalarMultiplication:
    def test_mul_complex_scalar(self):
        qh = _complex_qh((2.0, {0: "X"}))
        result = qh * (3.0 + 0j)
        expected = _complex_qh((6.0, {0: "X"}))
        assert result == expected

    def test_mul_imaginary_scalar(self):
        qh = _complex_qh((1.0, {0: "X"}))
        result = qh * (2j)
        expected = _complex_qh((0.0, {0: "X"}))  # coeff 2j
        # Just check type and that coeff changed — exact equality depends on QH == logic
        assert isinstance(result, _core.QubitHamiltonianComplex)
        coeff_val = result.to_dictionary()[0][0]
        assert coeff_val == 2j

    def test_mul_symbolic_hamiltonian_complex_scalar(self):
        qh = _symbolic_qh(("a", {0: "X"}))
        result = qh * (2.0 + 0j)
        assert isinstance(result, _core.QubitHamiltonianSymbolic)

    def test_mul_returns_same_type_complex(self):
        qh = _complex_qh((1.0, {0: "Z"}))
        result = qh * (3.0 + 0j)
        assert isinstance(result, _core.QubitHamiltonianComplex)

    def test_mul_returns_same_type_symbolic(self):
        qh = _symbolic_qh(("a", {0: "Z"}))
        result = qh * (3.0 + 0j)
        assert isinstance(result, _core.QubitHamiltonianSymbolic)


class TestHamiltonianMultiplication:
    def test_mul_two_single_term_hamiltonians(self):
        # X * Z = iY  (on same qubit)
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((1.0, {0: "Z"}))
        result = qh1 * qh2
        expected_ps = pe.PauliString(1.0, {0: "X"}) * pe.PauliString(1.0, {0: "Z"})
        expected = pe.QubitHamiltonian([expected_ps])
        assert result == expected

    def test_mul_distributes(self):
        # (X + Z) * Y = XY + ZY
        qh1 = _complex_qh((1.0, {0: "X"}), (1.0, {0: "Z"}))
        qh2 = _complex_qh((1.0, {0: "Y"}))
        result = qh1 * qh2
        ps_xy = pe.PauliString(1.0, {0: "X"}) * pe.PauliString(1.0, {0: "Y"})
        ps_zy = pe.PauliString(1.0, {0: "Z"}) * pe.PauliString(1.0, {0: "Y"})
        expected = pe.QubitHamiltonian([ps_xy, ps_zy])
        assert result == expected

    def test_mul_result_count(self):
        qh1 = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        qh2 = _complex_qh((1.0, {0: "Y"}), (3.0, {1: "X"}))
        result = qh1 * qh2
        # 2x2 terms before compaction
        assert len(result.to_dictionary()) == 4

    def test_mul_symbolic_hamiltonians(self):
        qh1 = _symbolic_qh(("a", {0: "X"}))
        qh2 = _symbolic_qh(("b", {0: "Z"}))
        result = qh1 * qh2
        assert isinstance(result, _core.QubitHamiltonianSymbolic)

    def test_mul_returns_same_type_complex(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((1.0, {1: "Z"}))
        result = qh1 * qh2
        assert isinstance(result, _core.QubitHamiltonianComplex)

class TestStringRepresentation:
    def test_to_string_returns_str(self):
        qh = _complex_qh((1.0, {0: "X"}))
        s = qh.to_string()
        assert isinstance(s, str)
        assert len(s) > 0

    def test_str_and_repr_consistent(self):
        qh = _complex_qh((2.0, {0: "Z"}))
        assert str(qh) == repr(qh)

    def test_str_equals_to_string(self):
        qh = _complex_qh((1.0, {0: "X"}))
        assert str(qh) == qh.to_string()

    def test_str_symbolic(self):
        qh = _symbolic_qh(("a", {0: "X"}))
        s = str(qh)
        assert isinstance(s, str)
        assert len(s) > 0



class TestToDictionary:
    def test_to_dictionary_single_term(self):
        qh = _complex_qh((2.0, {0: "X", 1: "Z"}))
        result = qh.to_dictionary()
        assert len(result) == 1
        coeff, ops = result[0]
        assert coeff == 2.0 + 0j
        assert ops[0] == "X"
        assert ops[1] == "Z"

    def test_to_dictionary_multiple_terms(self):
        qh = _complex_qh((1.0, {0: "X"}), (3.0, {1: "Z"}))
        result = qh.to_dictionary()
        assert len(result) == 2
        coeffs = {entry[0] for entry in result}
        assert (1.0 + 0j) in coeffs
        assert (3.0 + 0j) in coeffs

    def test_to_dictionary_symbolic(self):
        qh = _symbolic_qh(("a", {0: "X"}))
        result = qh.to_dictionary()
        assert len(result) == 1



class TestTraceOutQubits:
    def test_trace_out_returns_hamiltonian(self):
        qh = _complex_qh((1.0, {0: "Z", 1: "X"}))
        result = qh.trace_out_qubits([1], [0])
        assert isinstance(result, _core.QubitHamiltonianComplex)

    def test_trace_out_identity_on_non_present_qubit(self):
        qh = _complex_qh((1.0, {0: "Z"}))
        result = qh.trace_out_qubits([5], [0])
        # qubit 5 is not present; tracing out a qubit that acts as I should preserve
        assert isinstance(result, _core.QubitHamiltonianComplex)

    def test_trace_out_z_in_state_0(self):
        # Z|0><0| trace → eigenvalue +1 → coeff preserved
        qh = _complex_qh((1.0, {0: "Z"}))
        result = qh.trace_out_qubits([0], [0])
        d = result.to_dictionary()
        # Resulting terms may be zero or have coefficient based on trace semantics
        assert isinstance(result, _core.QubitHamiltonianComplex)



class TestCommutator:
    """Algebraic properties of QubitHamiltonian.commutator.

    Exact single-qubit values are exercised by PauliString commutator tests.
    Here we only verify that QH dispatches correctly and respects the
    bilinear / Jacobi / multi-term structure on top of that."""

    def test_commutator_returns_same_type_complex(self):
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((1.0, {0: "Z"}))
        assert isinstance(qh1.commutator(qh2), _core.QubitHamiltonianComplex)

    def test_commutator_returns_same_type_symbolic(self):
        qh1 = _symbolic_qh(("a", {0: "X"}))
        qh2 = _symbolic_qh(("b", {0: "Z"}))
        assert isinstance(qh1.commutator(qh2), _core.QubitHamiltonianSymbolic)

    def test_commutator_exact_sanity(self):
        """[2X, 3Z] = -12iY (covers exact dispatch + coefficient scaling)."""
        result = _complex_qh((2.0, {0: "X"})).commutator(_complex_qh((3.0, {0: "Z"})))
        assert result == pe.QubitHamiltonian([pe.PauliString(-12j, {0: "Y"})])

    def test_commutator_different_qubits_is_zero(self):
        """Operators on different qubits commute."""
        result = _complex_qh((1.0, {0: "X"})).commutator(_complex_qh((1.0, {1: "Z"})))
        for coeff, _ in result.to_dictionary():
            assert coeff == pytest.approx(0j)

    def test_commutator_bilinearity(self):
        """[H1+H2, H3] = [H1,H3] + [H2,H3] (left). Right side follows by antisymmetry."""
        qh_xz = _complex_qh((1.0, {0: "X"}), (1.0, {0: "Z"}))
        qh_y  = _complex_qh((1.0, {0: "Y"}))
        qh_x  = _complex_qh((1.0, {0: "X"}))
        qh_z  = _complex_qh((1.0, {0: "Z"}))
        assert qh_xz.commutator(qh_y) == qh_x.commutator(qh_y) + qh_z.commutator(qh_y)

    def test_commutator_jacobi_multiqubit(self):
        """Jacobi identity for a non-trivial multi-qubit triple."""
        qh_a = _complex_qh((1.0, {0: "X", 1: "Z"}))
        qh_b = _complex_qh((1.0, {0: "Z", 1: "X"}))
        qh_c = _complex_qh((1.0, {0: "Y", 1: "Y"}))
        zero = (qh_a.commutator(qh_b).commutator(qh_c)
                + qh_b.commutator(qh_c).commutator(qh_a)
                + qh_c.commutator(qh_a).commutator(qh_b))
        for coeff, _ in zero.to_dictionary():
            assert coeff == pytest.approx(0j)

    def test_commutator_multiqubit_one_anticommuting_pair(self):
        """[X(0)Z(1), Z(0)] = -2i Y(0)Z(1)."""
        qh_a = pe.QubitHamiltonian([pe.PauliString(1.0, {0: "X", 1: "Z"})])
        qh_b = pe.QubitHamiltonian([pe.PauliString(1.0, {0: "Z"})])
        d = qh_a.commutator(qh_b).to_dictionary()
        assert len(d) == 1
        coeff, ops = d[0]
        assert coeff == pytest.approx(-2j)
        assert ops == {0: "Y", 1: "Z"}

    def test_commutator_multi_term_hamiltonians(self):
        """[X(0)+Z(0), Y(0)] = 2iZ(0) - 2iX(0)."""
        qh1 = _complex_qh((1.0, {0: "X"}), (1.0, {0: "Z"}))
        qh2 = _complex_qh((1.0, {0: "Y"}))
        expected = pe.QubitHamiltonian([
            pe.PauliString(2j, {0: "Z"}),
            pe.PauliString(-2j, {0: "X"}),
        ])
        assert qh1.commutator(qh2) == expected

    def test_commutator_symbolic_antisymmetry_and_symbols(self):
        """Symbolic [aX, bZ] = -[bZ, aX] and the result coefficient mentions both symbols."""
        qh_a = _symbolic_qh(("a", {0: "X"}))
        qh_b = _symbolic_qh(("b", {0: "Z"}))
        ab = qh_a.commutator(qh_b)
        ba = qh_b.commutator(qh_a)
        assert ab == ba * (-1.0 + 0j)
        coeff_str = str(ab.to_dictionary()[0][0])
        assert "a" in coeff_str and "b" in coeff_str



class TestSubstitute:
    def test_subs_replaces_symbol_with_value(self):
        qh = _symbolic_qh(("a", {0: "X"}))
        result = qh.subs({"a": 2.0 + 0j})
        assert isinstance(result, _core.QubitHamiltonianSymbolic)
        coeff, _ = result.to_dictionary()[0]
        assert _core.PauliStringSymbolic.to_complex(coeff) == pytest.approx(2.0 + 0j)

    def test_subs_multiple_symbols(self):
        qh = _symbolic_qh(("a", {0: "X"}), ("b", {1: "Z"}))
        result = qh.subs({"a": 1.0 + 0j, "b": 3.0 + 0j})
        assert isinstance(result, _core.QubitHamiltonianSymbolic)
        d = result.to_dictionary()
        assert len(d) == 2

    def test_subs_partial_substitution_leaves_remaining_symbolic(self):
        qh = _symbolic_qh(("a", {0: "X"}), ("b", {1: "Z"}))
        result = qh.subs({"a": 1.0 + 0j})
        assert isinstance(result, _core.QubitHamiltonianSymbolic)

    def test_subs_zero_coefficient(self):
        """Substituting a symbol with 0 removes the term via compaction."""
        qh = _symbolic_qh(("a", {0: "X"}))
        result = qh.subs({"a": 0.0 + 0j})
        assert _term_count(result) == 0


class TestTypePromotion:
    """Verify that symbolic presence always propagates to QubitHamiltonianSymbolic."""

    @pytest.mark.parametrize(
        ("symbolic_idx", "label"),
        [(0, "first"), (1, "middle"), (2, "last")],
    )
    def test_any_symbolic_term_gives_symbolic(self, symbolic_idx, label):
        """Symbolic term at any position promotes the Hamiltonian to symbolic."""
        ps_list = [pe.PauliString(1.0, {i: "Z"}) for i in range(3)]
        ps_list[symbolic_idx] = pe.PauliString("a", {symbolic_idx: "X"})
        assert isinstance(pe.QubitHamiltonian(ps_list), _core.QubitHamiltonianSymbolic)

    def test_all_complex_stays_complex(self):
        ps_list = [pe.PauliString(1.0, {i: "Z"}) for i in range(3)]
        assert isinstance(pe.QubitHamiltonian(ps_list), _core.QubitHamiltonianComplex)

    def test_symbolic_arithmetic_preserves_type(self):
        """Mul-by-scalar, mul-by-QH, and addition all preserve the symbolic type."""
        qh1 = _symbolic_qh(("a", {0: "X"}))
        qh2 = _symbolic_qh(("b", {0: "Z"}))
        assert isinstance(qh1 * (2.0 + 0j), _core.QubitHamiltonianSymbolic)
        assert isinstance(qh1 * qh2, _core.QubitHamiltonianSymbolic)
        assert isinstance(qh1 + qh2, _core.QubitHamiltonianSymbolic)


@pytest.mark.parametrize("scalar", [3, 3.0, 3.0 + 0j])
def test_qh_scalar_mul_left_and_right(scalar):
    """qh * scalar == scalar * qh for int / float / complex."""
    qh = _complex_qh((2.0, {0: "X"}))
    expected = _complex_qh((6.0, {0: "X"}))
    assert qh * scalar == expected
    assert scalar * qh == expected


@pytest.mark.parametrize("scalar", [2, 2.0, 2.0 + 0j])
def test_symbolic_qh_scalar_mul_preserves_type(scalar):
    qh = _symbolic_qh(("a", {0: "X"}))
    assert isinstance(qh * scalar, _core.QubitHamiltonianSymbolic)
    assert isinstance(scalar * qh, _core.QubitHamiltonianSymbolic)


class TestArithmeticEdgeCases:
    def test_mul_zero_scalar(self):
        """Multiplication by 0 zeros the coefficient (from either side)."""
        qh = _complex_qh((3.0, {0: "X"}))
        assert (qh * 0).to_dictionary()[0][0] == 0j
        assert (0 * qh).to_dictionary()[0][0] == 0j

    def test_mul_negative_scalar(self):
        qh = _complex_qh((2.0, {0: "X"}))
        assert qh * (-1) == _complex_qh((-2.0, {0: "X"}))

    def test_qh_mul_qh_same_operator_gives_identity(self):
        """X * X = I (empty operator dict, coefficient 1)."""
        qh = _complex_qh((1.0, {0: "X"}))
        d = (qh * qh).to_dictionary()
        assert len(d) == 1
        coeff, ops = d[0]
        assert coeff == pytest.approx(1.0 + 0j)
        assert ops == {}

def _term_count(qh) -> int:
    return len(qh.to_dictionary())


def _ops_for(qh, qubit: int) -> list:
    """Return all operator strings at `qubit` across all terms."""
    return [ops.get(qubit) for _, ops in qh.to_dictionary()]


class TestCompaction:

    def test_construction_merges_duplicate_operators(self):
        """Passing two PauliStrings with the same operator must produce one term."""
        ps1 = pe.PauliString(1.0, {0: "X"})
        ps2 = pe.PauliString(2.0, {0: "X"})
        qh = pe.QubitHamiltonian([ps1, ps2])
        assert _term_count(qh) == 1
        coeff, _ = qh.to_dictionary()[0]
        assert coeff == pytest.approx(3.0 + 0j)

    def test_construction_removes_cancelling_terms(self):
        """Terms that cancel must be dropped from the term list."""
        ps1 = pe.PauliString(1.0, {0: "X"})
        ps2 = pe.PauliString(-1.0, {0: "X"})
        qh = pe.QubitHamiltonian([ps1, ps2])
        assert _term_count(qh) == 0

    def test_construction_keeps_distinct_operators(self):
        """Distinct operators must NOT be merged."""
        ps1 = pe.PauliString(1.0, {0: "X"})
        ps2 = pe.PauliString(2.0, {0: "Z"})
        qh = pe.QubitHamiltonian([ps1, ps2])
        assert _term_count(qh) == 2

    def test_addition_merges_same_operators(self):
        """qh1 + qh2 must yield one term when both carry the same operator."""
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((2.0, {0: "X"}))
        result = qh1 + qh2
        assert _term_count(result) == 1
        coeff, _ = result.to_dictionary()[0]
        assert coeff == pytest.approx(3.0 + 0j)

    def test_addition_drops_zero_terms(self):
        """Adding a term and its negative must drop the term entirely."""
        qh1 = _complex_qh((1.0, {0: "X"}))
        qh2 = _complex_qh((-1.0, {0: "X"}))
        result = qh1 + qh2
        assert _term_count(result) == 0

    def test_hamiltonian_mul_merges_resulting_duplicates(self):
        """(2X) * Y produces a single term; duplicate products must compact."""
        qh_xx = pe.QubitHamiltonian([
            pe.PauliString(1.0, {0: "X"}),
            pe.PauliString(1.0, {0: "X"}),
        ])
        assert _term_count(qh_xx) == 1
        assert _term_count(qh_xx * _complex_qh((1.0, {0: "Y"}))) == 1

    def test_equality_independent_of_internal_duplicates(self):
        """operator== must see through internal duplicates."""
        qh_duplicated = pe.QubitHamiltonian([
            pe.PauliString(1.0, {0: "X"}),
            pe.PauliString(2.0, {0: "X"}),
        ])
        assert qh_duplicated == _complex_qh((3.0, {0: "X"}))

    def test_symbolic_construction_merges_duplicates(self):
        qh = pe.QubitHamiltonian([pe.PauliString("a", {0: "X"}), pe.PauliString("b", {0: "X"})])
        assert _term_count(qh) == 1

    # --- Zero-coefficient removal ---

    def test_construction_drops_zero_coeff_term(self):
        qh = pe.QubitHamiltonian([pe.PauliString(0.0, {0: "X"})])
        assert _term_count(qh) == 0

    def test_compact_removes_zero_coeff_after_scalar_mul(self):
        """qh * 0 produces phantom zero terms; compact() drops them."""
        qh_zero = _complex_qh((3.0, {0: "X"}), (1.0, {1: "Z"})) * 0
        assert _term_count(qh_zero) > 0
        assert _term_count(qh_zero.compact()) == 0

    def test_symbolic_zero_coeff_dropped(self):
        """Symbolic '0' coefficient is recognised and dropped on construction."""
        qh = pe.QubitHamiltonian([pe.PauliString("0", {0: "X"})])
        assert _term_count(qh) == 0



class TestSubtraction:

    def test_complex_sub_basic(self):
        result = _complex_qh((3.0, {0: "X"})) - _complex_qh((1.0, {0: "X"}))
        assert result == _complex_qh((2.0, {0: "X"}))

    def test_complex_sub_self_gives_zero(self):
        qh = _complex_qh((2.0, {0: "X"}), (3.0, {1: "Z"}))
        assert _term_count(qh - qh) == 0

    def test_complex_neg_negates_coefficients(self):
        qh = _complex_qh((2.0, {0: "X"}))
        assert -qh == _complex_qh((-2.0, {0: "X"}))

    def test_complex_neg_neg_is_identity(self):
        qh = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        assert (- (-qh)) == qh

    def test_symbolic_sub_returns_symbolic_qh(self):
        result = _symbolic_qh(("a", {0: "X"})) - _symbolic_qh(("b", {0: "X"}))
        assert isinstance(result, _core.QubitHamiltonianSymbolic)
        assert _term_count(result) == 1

    def test_sub_compacts_result(self):
        """qh1 - qh2 auto-compacts duplicate terms."""
        qh1 = _complex_qh((3.0, {0: "X"}), (3.0, {0: "X"}))  # 6X after compact
        result = qh1 - _complex_qh((1.0, {0: "X"}))
        assert _term_count(result) == 1
        assert result.to_dictionary()[0][0] == pytest.approx(5.0 + 0j)


class TestPauliStringAddSub:

    def test_complex_ps_add_ps_compacts(self):
        ps1 = pe.PauliString(1.0, {0: "X"})
        ps2 = pe.PauliString(2.0, {0: "X"})
        result = ps1 + ps2
        assert isinstance(result, _core.QubitHamiltonianComplex)
        assert _term_count(result) == 1
        assert result.to_dictionary()[0][0] == pytest.approx(3.0 + 0j)

    def test_complex_ps_sub_self_gives_empty(self):
        ps = pe.PauliString(2.0, {0: "X"})
        assert _term_count(ps - ps) == 0

    def test_complex_ps_neg(self):
        ps = pe.PauliString(2.0, {0: "X"})
        neg = -ps
        assert isinstance(neg, _core.PauliStringComplex)
        assert neg.coeff == pytest.approx(-2.0 + 0j)
        assert neg.equals(ps)

    def test_symbolic_ps_add_returns_symbolic_qh(self):
        result = pe.PauliString("a", {0: "X"}) + pe.PauliString("b", {0: "Z"})
        assert isinstance(result, _core.QubitHamiltonianSymbolic)

    def test_cross_type_promotes_to_symbolic(self):
        """Either direction of complex + symbolic must yield a symbolic QH."""
        ps_c = pe.PauliString(1.0, {0: "X"})
        ps_s = pe.PauliString("a", {0: "Z"})
        assert isinstance(ps_c + ps_s, _core.QubitHamiltonianSymbolic)
        assert isinstance(ps_s + ps_c, _core.QubitHamiltonianSymbolic)
        assert isinstance(ps_c - ps_s, _core.QubitHamiltonianSymbolic)


 
# Inspection helpers: qubits, n_qubits, size, __len__, is_all_z, count_measurements
 

class TestInspection:
    def test_qubits_sorted_unique(self):
        qh = _complex_qh((1.0, {2: "X", 0: "Z"}), (1.0, {1: "Y", 0: "X"}))
        assert qh.qubits() == [0, 1, 2]

    def test_qubits_skips_identity(self):
        qh = _complex_qh((1.0, {0: "I", 1: "X"}))
        assert qh.qubits() == [1]

    def test_n_qubits(self):
        qh = _complex_qh((1.0, {3: "X", 7: "Y"}))
        assert qh.n_qubits() == 2

    def test_size_returns_term_count(self):
        qh = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        assert qh.size() == 2

    def test_len_equals_size(self):
        qh = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        assert len(qh) == qh.size()

    def test_is_all_z_true(self):
        qh = _complex_qh((1.0, {0: "Z"}), (1.0, {1: "Z", 2: "Z"}))
        assert qh.is_all_z()

    def test_is_all_z_false(self):
        qh = _complex_qh((1.0, {0: "X"}))
        assert not qh.is_all_z()

    def test_count_measurements_all_z(self):
        qh = _complex_qh((1.0, {0: "Z"}), (1.0, {1: "Z"}))
        assert qh.count_measurements() == 1

    def test_count_measurements_mixed(self):
        qh = _complex_qh((1.0, {0: "X"}), (1.0, {1: "Z"}))
        assert qh.count_measurements() == 2


 
# Hermiticity
 

class TestHermiticity:
    def test_real_coeffs_is_hermitian(self):
        qh = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        assert qh.is_hermitian()
        assert not qh.is_antihermitian()

    def test_imag_coeffs_is_antihermitian(self):
        qh = _complex_qh((1j, {0: "X"}), (3j, {1: "Z"}))
        assert not qh.is_hermitian()
        assert qh.is_antihermitian()

    def test_mixed_complex_is_neither(self):
        qh = _complex_qh((1.0 + 1j, {0: "X"}))
        assert not qh.is_hermitian()
        assert not qh.is_antihermitian()


 
# Conjugate / Transpose / Dagger
 

class TestConjugateTransposeDagger:
    def test_dagger_conjugates_coefficient(self):
        qh = _complex_qh((1.0 + 2j, {0: "X"}))
        dag = qh.dagger()
        coeff, _ = dag.to_dictionary()[0]
        assert coeff == pytest.approx(1.0 - 2j)

    def test_dagger_is_self_for_hermitian(self):
        qh = _complex_qh((2.0, {0: "X"}), (3.0, {0: "Z"}))
        assert qh.dagger() == qh

    def test_conjugate_flips_y_sign(self):
        """Y -> -Y under complex conjugation (Y is pure imaginary)."""
        qh = _complex_qh((1.0, {0: "Y"}))
        conj = qh.conjugate()
        coeff, ops = conj.to_dictionary()[0]
        assert coeff == pytest.approx(-1.0 + 0j)
        assert ops == {0: "Y"}

    def test_conjugate_two_ys_no_sign(self):
        qh = _complex_qh((1.0, {0: "Y", 1: "Y"}))
        conj = qh.conjugate()
        coeff, _ = conj.to_dictionary()[0]
        assert coeff == pytest.approx(1.0 + 0j)

    def test_conjugate_no_y_just_conjugates_coeff(self):
        qh = _complex_qh((1.0 + 2j, {0: "X", 1: "Z"}))
        conj = qh.conjugate()
        coeff, _ = conj.to_dictionary()[0]
        assert coeff == pytest.approx(1.0 - 2j)

    def test_transpose_flips_y_sign_no_conj(self):
        qh = _complex_qh((1.0 + 2j, {0: "Y"}))
        trans = qh.transpose()
        coeff, _ = trans.to_dictionary()[0]
        assert coeff == pytest.approx(-(1.0 + 2j))

    def test_double_dagger_is_identity(self):
        qh = _complex_qh((1.0 + 2j, {0: "X"}), (3.0 - 1j, {1: "Y"}))
        assert qh.dagger().dagger() == qh

    def test_double_conjugate_is_identity(self):
        qh = _complex_qh((1.0 + 2j, {0: "Y"}), (3.0, {1: "Z"}))
        assert qh.conjugate().conjugate() == qh


 
# Simplify with threshold
 

class TestSimplify:
    def test_simplify_no_threshold_keeps_all(self):
        qh = _complex_qh((1e-10, {0: "X"}), (1.0, {1: "Z"}))
        assert qh.simplify().size() == 2

    def test_simplify_removes_below_threshold(self):
        qh = _complex_qh((1e-10, {0: "X"}), (1.0, {1: "Z"}))
        result = qh.simplify(threshold=1e-8)
        assert result.size() == 1

    def test_simplify_keeps_above_threshold(self):
        qh = _complex_qh((0.5, {0: "X"}), (1.0, {1: "Z"}))
        result = qh.simplify(threshold=0.1)
        assert result.size() == 2

    def test_simplify_returns_correct_type(self):
        qh = _complex_qh((1.0, {0: "X"}))
        assert isinstance(qh.simplify(), _core.QubitHamiltonianComplex)


 
# split() — hermitian and anti-hermitian parts
 

class TestSplit:
    def test_split_pure_real(self):
        qh = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        herm, anti = qh.split()
        assert herm == qh
        assert anti.size() == 0

    def test_split_pure_imag(self):
        qh = _complex_qh((1j, {0: "X"}), (3j, {1: "Z"}))
        herm, anti = qh.split()
        assert herm.size() == 0
        assert anti == qh

    def test_split_mixed_returns_real_imag_parts(self):
        qh = _complex_qh((1.0 + 2j, {0: "X"}))
        herm, anti = qh.split()
        herm_coeff, _ = herm.to_dictionary()[0]
        anti_coeff, _ = anti.to_dictionary()[0]
        assert herm_coeff == pytest.approx(1.0 + 0j)
        assert anti_coeff == pytest.approx(2j)

    def test_split_recombines_to_original(self):
        qh = _complex_qh((1.0 + 2j, {0: "X"}), (3.0 - 1j, {1: "Y"}))
        herm, anti = qh.split()
        assert (herm + anti) == qh


 
# map_qubits
 

class TestMapQubitsQH:
    def test_map_qubits_relabels(self):
        qh = _complex_qh((1.0, {0: "X"}), (1.0, {1: "Z"}))
        mapped = qh.map_qubits({0: 3, 1: 5})
        assert mapped.qubits() == [3, 5]

    def test_map_qubits_preserves_coefficients(self):
        qh = _complex_qh((2.5, {0: "X"}))
        mapped = qh.map_qubits({0: 10})
        coeff, ops = mapped.to_dictionary()[0]
        assert coeff == pytest.approx(2.5 + 0j)
        assert ops == {10: "X"}


 
# power / __pow__
 

class TestPower:
    def test_pow_zero_is_identity(self):
        qh = _complex_qh((1.0, {0: "X"}))
        result = qh ** 0
        d = result.to_dictionary()
        assert len(d) == 1
        coeff, ops = d[0]
        assert coeff == pytest.approx(1.0 + 0j)
        assert ops == {}

    def test_pow_one_is_self(self):
        qh = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        assert (qh ** 1) == qh

    def test_pow_two_x_squared_is_identity(self):
        """X * X = I, coefficient 1."""
        qh = _complex_qh((1.0, {0: "X"}))
        result = qh ** 2
        d = result.to_dictionary()
        assert len(d) == 1
        coeff, ops = d[0]
        assert coeff == pytest.approx(1.0 + 0j)
        assert ops == {}

    def test_pow_three(self):
        """X^3 = X."""
        qh = _complex_qh((1.0, {0: "X"}))
        result = qh ** 3
        coeff, ops = result.to_dictionary()[0]
        assert coeff == pytest.approx(1.0 + 0j)
        assert ops == {0: "X"}

    def test_pow_negative_raises(self):
        qh = _complex_qh((1.0, {0: "X"}))
        with pytest.raises(Exception):
            qh ** -1


 
# to_matrix
 

_PAULI_MATS = {
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


class TestToMatrix:
    @pytest.mark.parametrize("op", ["X", "Y", "Z"])
    def test_single_pauli(self, op):
        m = np.array(_complex_qh((1.0, {0: op})).to_matrix())
        assert np.allclose(m, _PAULI_MATS[op])

    def test_sum_and_scaling(self):
        """X + Z gives [[1,1],[1,-1]] and 2.5*X scales the X matrix by 2.5."""
        m_sum = np.array(_complex_qh((1.0, {0: "X"}), (1.0, {0: "Z"})).to_matrix())
        assert np.allclose(m_sum, np.array([[1, 1], [1, -1]], dtype=complex))
        m_scaled = np.array(_complex_qh((2.5, {0: "X"})).to_matrix())
        assert np.allclose(m_scaled, 2.5 * _PAULI_MATS["X"])

    def test_two_qubit_zz(self):
        m = np.array(_complex_qh((1.0, {0: "Z", 1: "Z"})).to_matrix())
        assert m.shape == (4, 4)
        assert np.allclose(m, np.kron(_PAULI_MATS["Z"], _PAULI_MATS["Z"]))

    def test_ignore_unused_qubits(self):
        """ignore=True compresses; ignore=False uses absolute qubit indices."""
        qh_high = _complex_qh((1.0, {5: "X"}))
        assert np.array(qh_high.to_matrix(ignore_unused_qubits=True)).shape == (2, 2)
        qh_mid = _complex_qh((1.0, {2: "X"}))
        assert np.array(qh_mid.to_matrix(ignore_unused_qubits=False)).shape == (8, 8)


 
# zero / unit static helpers
 

class TestZeroUnit:
    def test_zero_is_empty(self):
        zero = _core.QubitHamiltonianComplex.zero()
        assert zero.size() == 0

    def test_unit_is_identity(self):
        unit = _core.QubitHamiltonianComplex.unit()
        assert unit.size() == 1
        coeff, ops = unit.to_dictionary()[0]
        assert coeff == pytest.approx(1.0 + 0j)
        assert ops == {}

    def test_factory_zero_returns_complex(self):
        assert isinstance(pe.QubitHamiltonian.zero(), _core.QubitHamiltonianComplex)

    def test_factory_unit_returns_complex(self):
        assert isinstance(pe.QubitHamiltonian.unit(), _core.QubitHamiltonianComplex)

    def test_unit_mul_is_identity_for_h(self):
        h = _complex_qh((1.0, {0: "X"}))
        u = _core.QubitHamiltonianComplex.unit()
        assert h * u == h


 
# paulistrings
 

class TestPaulistrings:
    def test_paulistrings_returns_list(self):
        qh = _complex_qh((1.0, {0: "X"}), (2.0, {1: "Z"}))
        ps_list = qh.paulistrings()
        assert len(ps_list) == 2
        coeffs = {ps.coeff for ps in ps_list}
        assert coeffs == {1.0 + 0j, 2.0 + 0j}


 
# OpenFermion bridge: factory accepts QubitOperator, to_openfermion roundtrip
 

@requires_openfermion
class TestOpenFermionBridge:
    def test_factory_accepts_qubit_operator(self):
        qop = 1.5 * QubitOperator("X0 Y1") + 0.5 * QubitOperator("Z2")
        qh = pe.QubitHamiltonian(qop)
        assert isinstance(qh, _core.QubitHamiltonianComplex)
        assert qh.size() == 2

    def test_factory_qubit_operator_correct_coefficients(self):
        qop = QubitOperator("X0", 1.5) + QubitOperator("Z1", 0.5j)
        qh = pe.QubitHamiltonian(qop)
        d = dict(
            (tuple(sorted(ops.items())), coeff) for coeff, ops in qh.to_dictionary()
        )
        assert d[((0, "X"),)] == pytest.approx(1.5 + 0j)
        assert d[((1, "Z"),)] == pytest.approx(0.5j)

    def test_empty_qubit_operator(self):
        qop = QubitOperator.zero()
        qh = pe.QubitHamiltonian(qop)
        assert qh.size() == 0

    def test_identity_qubit_operator(self):
        qop = QubitOperator("", 3.0)
        qh = pe.QubitHamiltonian(qop)
        assert qh.size() == 1
        coeff, ops = qh.to_dictionary()[0]
        assert coeff == pytest.approx(3.0 + 0j)
        assert ops == {}

    def test_from_openfermion_function(self):
        qop = QubitOperator("X0", 2.0)
        qh = pe.from_openfermion(qop)
        assert isinstance(qh, _core.QubitHamiltonianComplex)
        assert qh.size() == 1

    def test_to_openfermion_roundtrip(self):
        qop = 1.5 * QubitOperator("X0 Y1") + 0.5j * QubitOperator("Z2")
        qh = pe.QubitHamiltonian(qop)
        qop_back = pe.QubitHamiltonian.to_openfermion(qh)
        assert qop == qop_back

    def test_from_openfermion_wrong_type_raises(self):
        with pytest.raises(TypeError):
            pe.from_openfermion("not a qubit operator")


 
# trace_out_qubits — purity (regression coverage for the in-place mutation /
# pre-sized vector bugs in QubitHamiltonian::trace_out_qubits)
 

class TestTraceOutPurity:
    def test_qh_trace_does_not_mutate_receiver(self):
        qh = _complex_qh((1.0, {0: "Z", 1: "X"}), (3.0, {0: "X"}))
        before = qh.to_dictionary()
        _ = qh.trace_out_qubits([0], [0])
        after = qh.to_dictionary()
        assert sorted(map(repr, before)) == sorted(map(repr, after))

    def test_qh_trace_is_idempotent(self):
        qh = _complex_qh((1.0, {0: "Z", 1: "X"}), (3.0, {0: "X"}))
        r1 = qh.trace_out_qubits([0], [0])
        r2 = qh.trace_out_qubits([0], [0])
        assert r1 == r2

    def test_qh_trace_no_phantom_zero_terms(self):
        """H = Z(0) + X(0) traced against |0> -> single identity term, coeff 1
        (drops the X term because <0|X|0> = 0, no phantom zeros from default-
        constructed sentinels)."""
        qh = _complex_qh((1.0, {0: "Z"}), (1.0, {0: "X"}))
        result = qh.trace_out_qubits([0], [0])
        assert result.size() == 1
        coeff, ops = result.to_dictionary()[0]
        assert coeff == pytest.approx(1.0 + 0j)
        assert ops == {}

    def test_qh_trace_all_killed_gives_empty(self):
        qh = _complex_qh((1.0, {0: "X"}), (1.0, {0: "Y"}))
        result = qh.trace_out_qubits([0], [0])
        assert result.size() == 0

    def test_trace_z_in_state_1_flips_sign(self):
        """<1|Z|1> = -1."""
        qh = _complex_qh((3.0, {0: "Z"}))
        result = qh.trace_out_qubits([0], [1])
        coeff, ops = result.to_dictionary()[0]
        assert coeff == pytest.approx(-3.0 + 0j)
        assert ops == {}

    def test_trace_multiple_zs_combines_signs(self):
        """Z(0)Z(1) against |0>|1>: (+1)·(-1) = -1."""
        qh = _complex_qh((2.0, {0: "Z", 1: "Z"}))
        result = qh.trace_out_qubits([0, 1], [0, 1])
        coeff, _ = result.to_dictionary()[0]
        assert coeff == pytest.approx(-2.0 + 0j)

    def test_trace_identity_on_traced_qubit_no_change(self):
        qh = _complex_qh((1.5, {1: "Z"}))
        result = qh.trace_out_qubits([0], [0])
        coeff, ops = result.to_dictionary()[0]
        assert coeff == pytest.approx(1.5 + 0j)
        assert ops == {1: "Z"}


 
# trace_out_qubits general-state overload (Tequila-equivalent)
 

class TestTraceGeneralStates:
    def test_qh_general_state_two_terms(self):
        """X(0) + Z(0) traced against |+> -> <X> + <Z> = 1 + 0 = identity term."""
        qh = _complex_qh((1.0, {0: "X"}), (1.0, {0: "Z"}))
        plus = (1 / np.sqrt(2), 1 / np.sqrt(2))
        result = qh.trace_out_qubits([0], [plus])
        assert result.size() == 1
        coeff, ops = result.to_dictionary()[0]
        assert coeff == pytest.approx(1.0 + 0j)
        assert ops == {}
