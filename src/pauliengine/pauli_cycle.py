"""Module for PauliCycle and PauliCycleSum classes."""

from __future__ import annotations

from ._core import (
    PauliCycleComplex,
    PauliCycleSumComplex,
)


PauliCycle = PauliCycleComplex
PauliCycleSum = PauliCycleSumComplex

__all__ = ["PauliCycle", "PauliCycleSum"]
