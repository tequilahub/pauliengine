"""Timing primitives + sweep logic.

Separated from ``run.py`` so it can be imported and driven from a notebook
or an ad-hoc script.
"""

from __future__ import annotations

import gc
import time
from dataclasses import dataclass
from typing import Any, Callable, Literal

from scipy.stats import gmean, gstd

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
    times: list[float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "op": self.op,
            "coeff_kind": self.coeff_kind,
            "scaling_axis": self.scaling_axis,
            "n_terms": self.n_terms,
            "n_qubits": self.n_qubits,
            "repeats": self.repeats,
            "times_s": self.times,
        }


def _geo_stats(samples: list[float]) -> tuple[float, float]:
    """Geometric mean and an absolute spread from the geometric standard deviation.

    Runtimes are positive and right-skewed (occasional slow samples, never a
    negative one), so the geometric mean is a more robust central value than the
    arithmetic mean and the geometric standard deviation (``scipy.stats.gstd``, a
    dimensionless multiplicative factor >= 1) the natural measure of spread. We
    return the spread as ``gmean * (gstd - 1)`` so it stays an absolute quantity
    in seconds, compatible with the additive error bars the plots draw.
    """
    safe = [s if s > 0 else 1e-12 for s in samples]  # gmean/gstd need positives
    g = float(gmean(safe))
    if len(safe) < 2:
        return g, 0.0
    spread = g * (float(gstd(safe)) - 1.0)
    return g, spread


def summarize(samples: list[float]) -> tuple[float, float, float]:
    """``(min, gmean, spread)`` in seconds of raw timing samples (see ``_geo_stats``)."""
    g, spread = _geo_stats(samples)
    return min(samples), g, spread


def route_stats(m: dict[str, Any], key: str = "") -> tuple[float, float]:
    """``(gmean, spread)`` of one route in a stored measurement dict.

    Result files store the raw samples as ``<key>_times_s`` (or ``times_s`` for
    ``key=""``); the statistics are derived here, at analysis time. Older files
    (schema_version 1) only carry the precomputed ``<key>_time_mean_s`` /
    ``<key>_time_stdev_s`` and are read as-is.
    """
    prefix = f"{key}_" if key else ""
    samples = m.get(f"{prefix}times_s")
    if samples is not None:
        _, g, spread = summarize(samples)
        return g, spread
    return m[f"{prefix}time_mean_s"], m.get(f"{prefix}time_stdev_s", 0.0)


def has_route(measurements: list[dict[str, Any]], key: str) -> bool:
    """Whether any measurement carries timings for route ``key`` (new or old schema)."""
    return any(f"{key}_times_s" in m or f"{key}_time_mean_s" in m for m in measurements)


def _time_call(fn: Callable[[], Any], repeats: int, warmup: int) -> list[float]:
    """Run ``fn`` ``warmup`` times to warm caches, then ``repeats`` measured runs.

    Returns every measured sample in seconds, in run order. Statistics are left to
    the analysis side (``summarize`` / ``route_stats``) so no information is lost.
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
    return samples


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
    return Measurement(
        op=op,
        coeff_kind=coeff_kind,
        scaling_axis=scaling_axis,
        n_terms=n_terms,
        n_qubits=n_qubits,
        repeats=repeats,
        times=_time_call(fn, repeats=repeats, warmup=warmup),
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
