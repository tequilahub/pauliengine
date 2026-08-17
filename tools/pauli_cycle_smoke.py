"""Smoke test for PauliCycle: prints everything so you can eyeball the results.

Computes the commutator of two Pauli cycles **two ways** and shows them side by
side:

  1. cycle way   — a.commutator(b) stays in the Pauli-cycle basis (a
     PauliCycleSum), then .to_qubit_hamiltonian() materialises it.
  2. naive way   — expand both cycles to full QubitHamiltonians (all rotations)
     and commute those directly.

Both must yield the identical operator. Coefficients are non-trivial on purpose,
to show they are carried through.

Run it directly:

    python tools/pauli_cycle_smoke.py
"""

from __future__ import annotations

from pauliengine._core import (
    PauliCycleComplex,
    QubitHamiltonianComplex,
)


def ps_str(ps, n: int) -> str:
    """Render a PauliString as its n-character I/X/Y/Z pattern."""
    return "".join(ps.get_pauli_at_index(i) for i in range(n))


def fmt_coeff(c: complex) -> str:
    return f"{c.real:+.3f}{c.imag:+.3f}j"


def show_hamiltonian(title: str, h) -> None:
    """Pretty-print a QubitHamiltonian as sorted (coeff, ops) rows."""
    rows = []
    for coeff, ops in h.to_dictionary():
        key = "".join(f"{p}{q} " for q, p in sorted(ops.items())) or "I"
        rows.append((key, complex(coeff)))
    rows.sort()
    print(f"  {title}  ({len(rows)} terms):")
    for key, coeff in rows:
        print(f"      {fmt_coeff(coeff):>16}   {key}")


def banner(text: str) -> None:
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70)


def main() -> None:
    n = 4
    base_a, alpha = "XYII", complex(2.0, -1.0)
    base_b, beta = "ZIII", complex(0.5, 0.0)

    a = PauliCycleComplex(n, base_a, alpha)
    b = PauliCycleComplex(n, base_b, beta)

    banner("INPUT CYCLES")
    print(f"  n_qubits = {n}")
    print(f"  a: base = {base_a}  coeff = {fmt_coeff(alpha)}")
    print(f"     rotations: {[ps_str(r, n) for r in a.rotations()]}")
    print(f"     rotation coeffs: {[fmt_coeff(r.coeff) for r in a.rotations()]}")
    print(f"  b: base = {base_b}  coeff = {fmt_coeff(beta)}")
    print(f"     rotations: {[ps_str(r, n) for r in b.rotations()]}")

    # --- Way 1: commutator in the cycle basis --------------------------------
    banner("WAY 1 — cycle basis:  a.commutator(b)  ->  PauliCycleSum")
    s = a.commutator(b)
    print(f"  PauliCycleSum.coeff (paper 1/n prefactor) = {s.coeff}")
    print(f"  number of cycle terms in the sum          = {len(s.data)}")
    for i, cyc in enumerate(s.data):
        print(f"    term[{i}]: base = {ps_str(cyc.base, n)}  "
              f"coeff = {fmt_coeff(cyc.base.coeff)}  n_qubits = {cyc.n_qubits}")
    cycle_h = s.to_qubit_hamiltonian()
    print()
    show_hamiltonian("expanded to QubitHamiltonian", cycle_h)

    # --- Way 2: naive full-Hamiltonian commutator ----------------------------
    banner("WAY 2 — naive:  expand to full QubitHamiltonians, then commute")
    h_a = QubitHamiltonianComplex(a.rotations())
    h_b = QubitHamiltonianComplex(b.rotations())
    show_hamiltonian("H_a = sum of a's rotations", h_a)
    show_hamiltonian("H_b = sum of b's rotations", h_b)
    naive_h = h_a.commutator(h_b)
    print()
    show_hamiltonian("naive  [H_a, H_b]", naive_h)

    # --- Compare -------------------------------------------------------------
    banner("COMPARISON")
    equal = cycle_h == naive_h
    print(f"  cycle way == naive way : {equal}")
    print(f"  difference has {len(cycle_h - naive_h)} terms "
          f"(0 means identical operators)")
    print(f"  off by factor n? (naive == cycle * n): {naive_h == (cycle_h * float(n))}")
    print()
    print("  RESULT:", "OK — both ways agree" if equal else "MISMATCH!")


if __name__ == "__main__":
    main()
