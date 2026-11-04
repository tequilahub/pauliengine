"""Compatibility shim and loaders for the installed g-sim package.

g-sim's modules ``import PauliEngine as pe`` (capitalized), whereas this project
installs as ``pauliengine``. Its cycle module further constructs Pauli strings as
``PauliString(pauli_dict, coeff)`` (dict first), while pauliengine's factory takes
``PauliString(coeff, pauli_dict)``.

To let g-sim import and run against the locally installed pauliengine — with no
path configuration — we register a thin ``PauliEngine`` alias module in
``sys.modules`` that forwards to ``pauliengine`` and swaps the constructor
argument order. The gorbits (Numba) kernel does not use pauliengine at all; the
alias is only needed so the ``gsim`` package imports (its ``__init__`` pulls in
the cycle module).
"""

from __future__ import annotations

import sys
import types


def install_pauliengine_alias() -> bool:
    """Register a ``PauliEngine`` alias for the installed ``pauliengine``.

    Idempotent. Returns True if ``PauliEngine`` is importable afterwards.
    """
    if "PauliEngine" in sys.modules:
        return True
    try:
        import pauliengine
    except Exception:
        return False

    shim = types.ModuleType("PauliEngine")

    class _PauliStringCompat:
        """g-sim calls PauliString(pauli_dict, coeff); pauliengine wants (coeff, dict)."""

        def __call__(self, pauli_dict, coeff):
            return pauliengine.PauliString(coeff, pauli_dict)

        def to_complex(self, expr):
            if isinstance(expr, complex):
                return expr
            if isinstance(expr, (int, float)):
                return complex(expr)
            return pauliengine.PauliString.to_complex(expr)

    shim.PauliString = _PauliStringCompat()
    sys.modules["PauliEngine"] = shim
    return True


def load_orbits():
    """Return g-sim's (orbit_commutator, _precompute_factorials) or None."""
    install_pauliengine_alias()
    try:
        from gsim.gorbits import orbit_commutator, _precompute_factorials
        return orbit_commutator, _precompute_factorials
    except Exception as e:
        print(f"[gsim] gsim.gorbits not available: {e}", file=sys.stderr)
        return None


def load_cycles():
    """Return g-sim's (cycle_commutator, get_canonical_cycle) or None."""
    install_pauliengine_alias()
    try:
        from gsim.gcycles import cycle_commutator, get_canonical_cycle
        return cycle_commutator, get_canonical_cycle
    except Exception as e:
        print(f"[gsim] gsim.gcycles not available: {e}", file=sys.stderr)
        return None
