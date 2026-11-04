"""Smoke test: nested commutators [A, [A, ... [A, B]]] in PauliEngine agree with a
plain-Python computation after many iterations.

The Python reference does its own Pauli algebra (single-qubit product table), so it
is independent of PauliEngine. Both use IEEE doubles; only the summation order
differs, so the coefficients agree up to rounding (checked with a tolerance).

Run from the repo root:

    python tools/commutator_precision_smoke.py
"""

from __future__ import annotations

import random
import sys

import pauliengine as pe

# single-qubit products: (a, b) -> (phase, result), e.g. X*Y = iZ
_MUL = {
    ("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"),
    ("Y", "Z"): (1j, "X"), ("Z", "Y"): (-1j, "X"),
    ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y"),
}


def pauli_product(p: tuple, q: tuple) -> tuple[complex, tuple]:
    """Product of two Pauli strings given as sorted ((qubit, pauli), ...) tuples."""
    a, b = dict(p), dict(q)
    phase, out = 1 + 0j, {}
    for k in a.keys() | b.keys():
        x, y = a.get(k), b.get(k)
        if x is None or y is None:
            out[k] = x or y
        elif x != y:
            ph, r = _MUL[(x, y)]
            phase *= ph
            out[k] = r
        # x == y gives the identity on that qubit
    return phase, tuple(sorted(out.items()))


def py_commutator(A: dict, B: dict) -> dict:
    """[A, B] for {pauli_tuple: coeff} dicts; only anticommuting pairs contribute 2*PQ."""
    out: dict = {}
    for p, ca in A.items():
        for q, cb in B.items():
            ph_pq, r = pauli_product(p, q)
            ph_qp, _ = pauli_product(q, p)
            if ph_pq == ph_qp:              # commuting pair: exactly zero
                continue
            out[r] = out.get(r, 0j) + ca * cb * 2 * ph_pq
    return {k: v for k, v in out.items() if v != 0}


def random_sum(rng: random.Random, n_qubits: int, n_terms: int) -> dict:
    out = {}
    for _ in range(n_terms):
        ops = {q: rng.choice("XYZ") for q in rng.sample(range(n_qubits), rng.randint(1, 3))}
        out[tuple(sorted(ops.items()))] = complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
    return out


def to_pe(d: dict):
    return pe.QubitHamiltonian([pe.PauliString(c, dict(k)) for k, c in d.items()])


def from_pe(h) -> dict:
    return {tuple(sorted(ops.items())): complex(c) for c, ops in h.to_dictionary()}


def main(iterations: int = 8, n_qubits: int = 6, n_terms: int = 10, rtol: float = 1e-12) -> int:
    rng = random.Random(0)
    A, B = random_sum(rng, n_qubits, n_terms), random_sum(rng, n_qubits, n_terms)
    A_pe, cur_pe, cur_py = to_pe(A), to_pe(B), dict(B)

    ok = True
    print(f"{'iter':>4} {'terms':>6} {'bit-identical':>14} {'max rel. diff':>14}")
    for it in range(1, iterations + 1):
        cur_pe = A_pe.commutator(cur_pe)
        cur_py = py_commutator(A, cur_py)
        got = from_pe(cur_pe)

        scale = max(abs(v) for v in cur_py.values())
        keys = got.keys() | cur_py.keys()
        diff = max(abs(got.get(k, 0j) - cur_py.get(k, 0j)) for k in keys) / scale
        same = sum(got.get(k) == cur_py.get(k) for k in keys)
        print(f"{it:>4} {len(keys):>6} {same:>8}/{len(keys):<5} {diff:>14.2e}")

        if diff > rtol:
            ok = False
    print("OK" if ok else f"FAIL: coefficients differ by more than rtol={rtol:g}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
