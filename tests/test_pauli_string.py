import pytest

import pauliengine as pe


@pytest.mark.parametrize(("pauli_string", "coeff", "expected_x", "expected_z"), [
    ({0:"I", 1:"Z", 2:"X"}, 5.0, [6], [2]), # using dict input
    ({10: 'Z', 2: 'X'}, 6.0, [2**10 + 2**2] , [2**10]), # using dict input with non-consecutive qubits
])
def test_pauli_string_valid_pauli_string_input(pauli_string, coeff, expected_x, expected_z):
    ps = pe.PauliString(pauli_string, coeff)
    # TODO: coeff is an Expression object, need to be exposed to Python somehow. 
    # assert ps.coeff == coeff
    assert ps.x == expected_x
    assert ps.y == expected_z

@pytest.mark.xfail
@pytest.mark.parametrize("invalid_pauli_string", [
    {0: "A", 2: "B", 3: "Z"}, # invalid character
    {0: "i", 2: "z", 3: "y"}, # lowercase character
])
def test_pauli_string_invalid_pauli_string_input(invalid_pauli_string):
    with pytest.raises(ValueError):
        pe.PauliString(invalid_pauli_string, 1.0)


@pytest.mark.parametrize(("pauli_string", "other", "expected_product"), [
    (pe.PauliString({0: "X", 1: "Z", 2: "Y", 3: "Z"}, 1.0),
     pe.PauliString({0: "Z", 1: "Z", 3: "X"}, 2.0),
     pe.PauliString({0: "Y", 2: "Y", 3: "Y"}, 2.0)), # XZYZ * ZZIX = YIYY 

    (pe.PauliString({0: "Z", 2: "Z"}, 2.0j),
     pe.PauliString({0: "Y", 2: "Y"}, 3.0),
     pe.PauliString({0: "X", 2: "X"}, -6.0j)), # XZ * ZX = YY

    (pe.PauliString({0: "X", 1: "Y"}, -4.0),
     1.0 - 2.0j,
     pe.PauliString({0: "X", 1: "Y"}, -4.0 + 8.0j)), # -4.0XY * 1.0 - 2.0j
])

class TestPauliStringMultiplication:

    def test_pauli_string_multiplication(self, pauli_string, other, expected_product):
        product = pauli_string * other
        assert expected_product == product

    def test_pauli_string_multiplication_inplace(self, pauli_string, other, expected_product):
        pauli_string *= other
        assert expected_product == pauli_string

