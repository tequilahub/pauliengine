"""Timing primitives + sweep logic.

Separated from ``run.py`` so it can be imported and driven from a notebook
or an ad-hoc script.
"""

from __future__ import annotations

import gc
import math
import time
from dataclasses import dataclass
from typing import Any, Callable, Literal

from . import generate

Operation = Literal["multiply", "commutator"]


@dataclass
class Measurement:
    """Result of timing one operation at one point on the sweep axis."""

    op: Operation
    coeff_kind: generate.CoeffKind
    scaling_axis: Literal["n_terms", "n_qubits"]
    n_terms: int
    n_qubits: int
    repeats: int
    time_min: float
    time_mean: float
    time_stdev: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "op": self.op,
            "coeff_kind": self.coeff_kind,
            "scaling_axis": self.scaling_axis,
            "n_terms": self.n_terms,
            "n_qubits": self.n_qubits,
            "repeats": self.repeats,
            "time_min_s": self.time_min,
            "time_mean_s": self.time_mean,
            "time_stdev_s": self.time_stdev,
        }


def _stdev(samples: list[float], mean: float) -> float:
    if len(samples) < 2:
        return 0.0
    var = sum((x - mean) ** 2 for x in samples) / (len(samples) - 1)
    return math.sqrt(var)


def _time_call(fn: Callable[[], Any], repeats: int, warmup: int) -> tuple[float, float, float]:
    """Run ``fn`` ``warmup`` times to warm caches, then ``repeats`` measured runs.

    Returns ``(min, mean, stdev)`` in seconds.
    """
    for _ in range(warmup):
        fn()
    gc.collect()
    gc.disable()
    try:
        samples: list[float] = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            fn()
            samples.append(time.perf_counter() - t0)
    finally:
        gc.enable()
    mean = sum(samples) / len(samples)
    return min(samples), mean, _stdev(samples, mean)


def _make_op(op: Operation, h1, h2) -> Callable[[], Any]:
    if op == "multiply":
        return lambda: h1 * h2
    if op == "commutator":
        return lambda: h1.commutator(h2)
    raise ValueError(f"unknown op: {op!r}")


def measure_point(
    op: Operation,
    coeff_kind: generate.CoeffKind,
    scaling_axis: Literal["n_terms", "n_qubits"],
    n_terms: int,
    n_qubits: int,
    repeats: int,
    warmup: int,
    seed: int,
) -> Measurement:
    """Measure one point of one sweep."""
    h1 = generate.random_hamiltonian(n_terms, n_qubits, coeff_kind, seed=seed)
    h2 = generate.random_hamiltonian(n_terms, n_qubits, coeff_kind, seed=seed + 1)
    fn = _make_op(op, h1, h2)
    t_min, t_mean, t_std = _time_call(fn, repeats=repeats, warmup=warmup)
    return Measurement(
        op=op,
        coeff_kind=coeff_kind,
        scaling_axis=scaling_axis,
        n_terms=n_terms,
        n_qubits=n_qubits,
        repeats=repeats,
        time_min=t_min,
        time_mean=t_mean,
        time_stdev=t_std,
    )


def scan_n_terms(
    ops: list[Operation],
    coeff_kinds: list[generate.CoeffKind],
    n_terms_values: list[int],
    n_qubits_fixed: int,
    repeats: int,
    warmup: int,
    seed: int,
    progress: Callable[[str], None] = lambda _msg: None,
) -> list[Measurement]:
    """Vary ``n_terms`` with ``n_qubits`` held fixed."""
    results: list[Measurement] = []
    for coeff_kind in coeff_kinds:
        for op in ops:
            for n in n_terms_values:
                progress(f"[scan n_terms] {op} {coeff_kind}: n_terms={n} n_qubits={n_qubits_fixed}")
                results.append(measure_point(
                    op=op,
                    coeff_kind=coeff_kind,
                    scaling_axis="n_terms",
                    n_terms=n,
                    n_qubits=n_qubits_fixed,
                    repeats=repeats,
                    warmup=warmup,
                    seed=seed,
                ))
    return results


def scan_n_qubits(
    ops: list[Operation],
    coeff_kinds: list[generate.CoeffKind],
    n_qubits_values: list[int],
    n_terms_fixed: int,
    repeats: int,
    warmup: int,
    seed: int,
    progress: Callable[[str], None] = lambda _msg: None,
) -> list[Measurement]:
    """Vary ``n_qubits`` with ``n_terms`` held fixed."""
    results: list[Measurement] = []
    for coeff_kind in coeff_kinds:
        for op in ops:
            for q in n_qubits_values:
                progress(f"[scan n_qubits] {op} {coeff_kind}: n_terms={n_terms_fixed} n_qubits={q}")
                results.append(measure_point(
                    op=op,
                    coeff_kind=coeff_kind,
                    scaling_axis="n_qubits",
                    n_terms=n_terms_fixed,
                    n_qubits=q,
                    repeats=repeats,
                    warmup=warmup,
                    seed=seed,
                ))
    return results
